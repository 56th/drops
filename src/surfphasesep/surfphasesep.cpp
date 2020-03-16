/*
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
*/

#include "surfactant/ifacetransp.h"
#include "misc/params.h"
#include "geom/builder.h"
#include "levelset/levelset.h"
#include "levelset/adaptriang.h"
#include "levelset/surfacetension.h"
#include "misc/dynamicload.h"
// #include "out/vtkOut.h"
#include "num/bndData.h"
#include "num/precond.h"
#include "parallel/exchange.h"

#include "num/oseensolver.h"
#include <fstream>

#include "surfphasesep/surfphasesep_funcs.h"
#include "surfnavierstokes/surfnavierstokes_utils.h"
#include "surfnavierstokes/surfnavierstokes_tests.h"

// vtk output
#include "surfnavierstokes/VTKWriter.hpp"
// logger
#include "SingletonLogger.hpp"
// initial data
#include "SurfCahnHilliardData.hpp"
// for chemical potential
#include "hermite_cubic.hpp"

using namespace DROPS;

void InitVecLaplace(const MultiGridCL& MG, LevelsetP2CL& lset, DROPS::VecDescCL& rhs, DROPS::VecDescCL& vSol, DROPS::VecDescCL& pSol,
                    InstatVectorFunction f_rhs, InstatVectorFunction f_vsol, InstatScalarFunction f_psol, double t = 0.0)
{
    if( vSol.RowIdx->NumUnknownsEdge()) {
        DROPS::SetupInterfaceVectorRhsP2(MG, &rhs, lset.Phi, lset.GetBndData(), f_rhs);
    } else {
        DROPS::SetupInterfaceVectorRhsP1(MG, &rhs, lset.Phi, lset.GetBndData(), f_rhs, t);
    }
    InitScalar(MG, pSol, f_psol);
    InitVector(MG, vSol, f_vsol);
}


int main (int argc, char* argv[])
{
    srand(time(NULL));
    auto& logger = SingletonLogger::instance();
    try {
    logger.beg("read input .json");
      auto& inpJSON = P;
      read_parameter_file_from_cmdline( inpJSON, argc, argv, "../../param/surfphasesep/No_Bnd_Condition_ch.json");
      dynamicLoad(inpJSON.get<std::string>("General.DynamicLibsPrefix"), inpJSON.get<std::vector<std::string>>("General.DynamicLibs"));
      auto dirName = inpJSON.get<std::string>("Output.Directory");
      system(("mkdir -p " + dirName).c_str());
      system(("rm -r -f " + dirName + "/*").c_str());
      system((
          "mkdir " + dirName + "/matrices ; "
          "mkdir " + dirName + "/vtk ; "
          "mkdir " + dirName + "/tmp ; "
          "mkdir " + dirName + "/stats"
      ).c_str());
      {
          std::ofstream input(dirName + "/input.json");
          input << inpJSON;
          logger.buf << inpJSON;
          logger.log();
      }
    logger.end();
    auto surfCahnHilliardData = SurfCahnHilliardDataFactory(inpJSON);
    logger.buf << surfCahnHilliardData.description;
    logger.log();
    logger.beg("set up chemical potential");
        auto xi = inpJSON.get<double>("SurfCahnHilliard.ChemicalPotentialScaling");
        auto c0 = inpJSON.get<double>("SurfCahnHilliard.c_0");
        auto c0_l = std::min(c0, 1. - c0) / sqrt(3.);
        auto chemicalPotential = [&](double c) {
            if (c < 0.) return xi * (1. - c0) * c;
            if (c > 1.) return xi * c0 * (c - 1.);
            double x[1], f[1], d[1], s[1], t[1];
            x[0] = c;
            if      (c < c0 - c0_l) hermite_cubic_value(0., 0., 1. - c0, c0 - c0_l, 1. / (12. * sqrt(3.)), 0., 1, x, f, d, s, t);
            else if (c < c0)        hermite_cubic_value(c0 - c0_l, 1. / (12. * sqrt(3.)), 0., c0, 0., -std::max(c0 * c0, (1. - c0) * (1. - c0)), 1, x, f, d, s, t);
            else if (c < c0 + c0_l) hermite_cubic_value(c0, 0., -std::max(c0 * c0, (1. - c0) * (1. - c0)), c0 + c0_l, -1. / (12. * sqrt(3.)), 0., 1, x, f, d, s, t);
            else                    hermite_cubic_value(c0 + c0_l, -1. / (12. * sqrt(3.)), 0., 1., 0., c0, 1, x, f, d, s, t);
            return xi * f[0];
        };
    logger.end();
    // build initial mesh
    std::cout << "Setting up interface-PDE:\n";
    std::auto_ptr<DROPS::MGBuilderCL> builder( DROPS::make_MGBuilder( P));
    DROPS::MultiGridCL mg( *builder);
    const DROPS::ParamCL::ptree_type* ch= 0;
    try {
        ch= &P.get_child( "Mesh.Periodicity");
    }
    catch (DROPS::DROPSParamErrCL) {}
    if (ch)
        read_PeriodicBoundaries( mg, *ch);

    // adaptive mesh refinement based on level set function
    typedef DROPS::DistMarkingStrategyCL InitMarkerT;
    InitMarkerT initmarker( surfCahnHilliardData.surface.phi, P.get<double>("Mesh.AdaptRef.Width"), P.get<int>("Mesh.AdaptRef.CoarsestLevel"), P.get<int>("Mesh.AdaptRef.FinestLevel") );
    //adap.set_marking_strategy( &initmarker );
      DROPS::AdapTriangCL adap( mg, &initmarker );

      adap.MakeInitialTriang();
    adap.set_marking_strategy( 0 );

    // create level set
    InstatScalarFunction sigma (0);
    SurfaceTensionCL sf( sigma, 0);

    BndDataCL<double> lsbnd( 0);
    read_BndData( lsbnd, mg, P.get_child( "Levelset.BndData"));

    DROPS::LevelsetP2CL & lset( * DROPS::LevelsetP2CL::Create( mg, lsbnd, sf) );

    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
//    LinearLSInit( mg, lset.Phi, surfCahnHilliardData.surface.phi);
    lset.Init( surfCahnHilliardData.surface.phi);

    // parse FE types and some other parameters from json file
    std::string FE = P.get<std::string>("SurfCahnHilliard.FE");
    std::string velFE, prFE, LgFE;
    velFE = FE.substr(0,2);
    prFE = FE.substr(2,2);
    if( FE.size()>4) {
        LgFE = FE.substr(4,2);
    } else {
        LgFE = FE.substr(2,2);
    }


    std::string model = P.get<std::string>("SurfCahnHilliard.model");
    std::string testcase = P.get<std::string>("SurfCahnHilliard.testcase");
    double h = P.get<DROPS::Point3DCL>("Mesh.E1")[0]/P.get<double>("Mesh.N1")*std::pow(2., -P.get<double>("Mesh.AdaptRef.FinestLevel"));

    auto T = P.get<double>("Time.FinalTime");
    auto dt = P.get<double>("Time.StepSize");
    auto alpha = inpJSON.get<double>("SurfCahnHilliard.VolumeStab.Factor") * pow(h, inpJSON.get<double>("SurfCahnHilliard.VolumeStab.Power"));
    std::cout << "h is: " << h << std::endl;
    std::cout << "dt is: " << dt << std::endl;
    ParameterNS::h = h;

    auto S = inpJSON.get<double>("SurfCahnHilliard.Beta_s");

    double sigm = P.get<double>("SurfCahnHilliard.MobilityScaling");
    double eps = P.get<double>("SurfCahnHilliard.Epsilon");
    ParameterNS::sigma = sigm;
    ParameterNS::eps = eps;

    // construct FE spaces
    BndDataCL<Point3DCL> vbnd( 0);
    read_BndData( vbnd, mg, P.get_child( "Stokes.VelocityBndData"));

//  BndDataCL<double> vscalarbnd( 0);
//  read_BndData( vscalarbnd, mg, P.get_child( "Stokes.VelocityBndData"));

    BndDataCL<double> pbnd( 0);
    read_BndData( pbnd, mg, P.get_child( "Stokes.PressureBndData"));

    BndDataCL<double> chibnd( 0);
    read_BndData( chibnd, mg, P.get_child( "Stokes.VolumeFractionBndData"));
      BndDataCL<double> omegabnd( 0);
      read_BndData( omegabnd, mg, P.get_child( "Stokes.ChemPotentialBndData"));

    DROPS::IdxDescCL ifaceVecP1idx( vecP1IF_FE, vbnd);
    DROPS::IdxDescCL ifaceP1idx( P1IF_FE, pbnd);
    DROPS::IdxDescCL vecP1idx( vecP1_FE, vbnd);
    DROPS::IdxDescCL P1FEidx( P1_FE, pbnd);

    ifaceVecP1idx.GetXidx().SetBound( P.get<double>("SurfTransp.OmitBound"));
    ifaceVecP1idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
    ifaceP1idx.GetXidx().SetBound( P.get<double>("SurfTransp.OmitBound"));
    ifaceP1idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
    vecP1idx.GetXidx().SetBound( P.get<double>("SurfTransp.OmitBound"));
    vecP1idx.CreateNumbering(mg.GetLastLevel(), mg);
    P1FEidx.GetXidx().SetBound( P.get<double>("SurfTransp.OmitBound"));
    P1FEidx.CreateNumbering(mg.GetLastLevel(), mg);

    // construct FE vectors (initialized with zero)
     DROPS::VecDescCL v, omega, omegaSol, chi, chi_ext, chi_prev, chi_prev_prev, chi_BDF1, chi_BDF2, chi_extrap, chiSol, energy, well_potential, rhs3, instantrhs3, rhs4, instantrhs4;
     if( !FE.compare("P1P1")) {
    	 v.SetIdx( &ifaceVecP1idx);
    	 chi.SetIdx( &ifaceP1idx);
    	 chi_ext.SetIdx(&P1FEidx);
    	 chi_prev.SetIdx(&ifaceP1idx);
         chi_prev_prev.SetIdx(&ifaceP1idx);
         chi_BDF1.SetIdx( &ifaceP1idx);
         chi_BDF2.SetIdx( &ifaceP1idx);
         chi_extrap.SetIdx(&ifaceP1idx);
         well_potential.SetIdx( &ifaceP1idx);
         energy.SetIdx( &ifaceP1idx);
    	 chiSol.SetIdx( &ifaceP1idx);
    	 omega.SetIdx( &ifaceP1idx);
    	 omegaSol.SetIdx( &ifaceP1idx);
         rhs3.SetIdx(&ifaceP1idx);
         instantrhs3.SetIdx(&ifaceP1idx);
         rhs4.SetIdx(&ifaceP1idx);
         instantrhs4.SetIdx(&ifaceP1idx);

     }

    // setup matrices
    DROPS::MatDescCL Laplace, Mass, Normal_stab, Tangent_stab, Volume_stab, LaplaceM, Gprimeprime;
    MatrixCL A,B,C,D, Precond3, Precond4;
    if( !FE.compare("P1P1")) {
        Laplace.SetIdx( &ifaceP1idx, &ifaceP1idx);
        LaplaceM.SetIdx( &ifaceP1idx, &ifaceP1idx);
        Gprimeprime.SetIdx( &ifaceP1idx, &ifaceP1idx);
        Mass.SetIdx(&ifaceP1idx, &ifaceP1idx);
        Normal_stab.SetIdx(&ifaceP1idx, &ifaceP1idx);
        Tangent_stab.SetIdx(&ifaceP1idx, &ifaceP1idx);
        Volume_stab.SetIdx(&ifaceP1idx, &ifaceP1idx);
        SetupCahnHilliardIF_P1P1(mg,  &Mass,&Normal_stab, &Tangent_stab, &Volume_stab, &Laplace, &LaplaceM, &Gprimeprime, lset.Phi, lset.GetBndData(), v, vbnd, chi, chibnd);
        if (model=="CahnHilliard")
            Precond4.LinComb(S, Mass.Data, alpha*eps*eps, Volume_stab.Data, eps*eps, Laplace.Data);
        else if (model=="AllenCahn")
            Precond4.LinComb(1.0    , Mass.Data, alpha, Volume_stab.Data,   eps, Laplace.Data);
    }
      typedef SSORPcCL      SymmPcPcT;
      SymmPcPcT symmPcPc_;

      typedef PCGSolverCL<SymmPcPcT> PCGSolverT;
      std::stringstream Precondstream3;
      std::stringstream Precondstream4;
      PCGSolverT PCGSolver3 (symmPcPc_, P.get<int>("Solver.PcAIter"), P.get<double>("Solver.PcATol"), true, &Precondstream3);
      PCGSolverT PCGSolver4 (symmPcPc_, P.get<int>("Solver.PcBIter"), P.get<double>("Solver.PcBTol"), true, &Precondstream4);

      SchurPreBaseCL *spc3_ = new SurfaceLaplacePreCL<PCGSolverT>( Precond3, PCGSolver3);
      SchurPreBaseCL *spc4_ = new SurfaceLaplacePreCL<PCGSolverT>( Precond4, PCGSolver4);

      std::stringstream Globalstream;
     /* DummyPcCL Ident;
      typedef BlockPreCL<DummyPcCL, DummyPcCL, DiagBlockPreCL>  DiagBlockPcT;
      DiagBlockPcT    *DBlock = new DiagBlockPcT(Ident, Ident);*/


      typedef BlockPreCL<SchurPreBaseCL, SchurPreBaseCL, DiagBlockPreCL>  DiagBlockPcT;
      DiagBlockPcT    *DBlock = new DiagBlockPcT(*spc3_, *spc4_);

/*
      typedef BlockPreCL<ExpensivePreBaseCL, SchurPreBaseCL,  DiagSpdBlockPreCL>  DiagBlockPcT;
      DiagBlockPcT    *DBlock = new DiagBlockPcT(PCGPc_, *spc4_);
*/



      typedef GMResSolverCL<SchurPreBaseCL> GMResT;
      GMResT *AC_solver = new GMResT( *spc3_, P.get<int>("Solver.Iter"), P.get<int>("Solver.Iter"), P.get<double>("Solver.Tol"),
              /*bool relative=*/ false, /*bool calculate2norm=*/ false, /*PreMethGMRES method=*/ LeftPreconditioning,
              /*bool mod =*/ true,     /* bool useModGS = */false, &Globalstream);

      typedef GMResSolverCL<DiagBlockPcT> GMResBlockT;
      GMResBlockT *GMRes_ = new GMResBlockT( *DBlock, P.get<int>("Solver.Iter"), P.get<int>("Solver.Iter"), P.get<double>("Solver.Tol"),
              /*bool relative=*/ false, /*bool calculate2norm=*/ false, /*PreMethGMRES method=*/ LeftPreconditioning,
              /*bool mod = */true, /* bool useModGS =*/ false, &Globalstream);
      BlockMatrixSolverCL<GMResBlockT> *CH_solver= new BlockMatrixSolverCL<GMResBlockT>( *GMRes_);

      SolverBaseCL *solver;
      if (model=="CahnHilliard")
      {
            solver=CH_solver;
      }
      else if (model=="AllenCahn")
      {
          solver=AC_solver;
      }

        VectorCL unityVector(1., Mass.Data.num_rows());

        VTKWriter vtkWriter(dirName + "/vtk/separation", mg, inpJSON.get<bool>("Output.Binary"));
        vtkWriter
                .add(VTKWriter::VTKVar({"level-set", &lset.Phi.Data, VTKWriter::VTKVar::Type::P2}))
                .add(VTKWriter::VTKVar({"c_h", &chi_ext.Data, VTKWriter::VTKVar::Type::P1}));
        auto everyStep = inpJSON.get<int>("Output.EveryStep");
        auto writeVTK = [&](double t) {
            Extend(mg, chi, chi_ext);
            vtkWriter.write(t);
        };

        	BndDataCL<Point3DCL> bndvec = vbnd;
        	BndDataCL<double> bndscalar = pbnd;
        	auto t = 0.;
        	auto e = 0.;
            InitScalar(mg, chi, surfCahnHilliardData.chi, 0.);
            InitScalar(mg, chiSol, surfCahnHilliardData.chi, 0.);
            size_t numbOfTries = 0;
            auto vtkExported = false;
            auto exportStats = [&](size_t i) {
                std::ofstream stats(dirName + "/stats/t_" + std::to_string(i) + ".json");
                ParamCL tJSON;
                tJSON.put("h", h);
                tJSON.put("rho", alpha);
                tJSON.put("t", t);
                tJSON.put("c_h.Max", chi.Data.max());
                tJSON.put("c_h.Min", chi.Data.min());
                tJSON.put("Integral.PerimeterEstimate", eps * dot(Laplace.Data * chi.Data, chi.Data));
                // for (int i = 0; i < energy.Data.size(); i++) energy.Data[i] = potential(chi.Data[i]);
                // tJSON.put("Integral.LyapunovEnergy", (eps / 2.) * dot(Laplace.Data * chi.Data, chi.Data) + (1. / eps) * dot(Mass.Data * energy.Data, unityVector));
                auto surfaceArea = dot(Mass.Data * unityVector, unityVector);
                tJSON.put("Integral.SurfaceArea", surfaceArea);
                tJSON.put("Integral.RaftFraction", dot(Mass.Data * chi.Data, unityVector) / surfaceArea);
                if (i > 0) {
                    tJSON.put("dt", dt);
                    tJSON.put("Integral.Solver.Outer.TotalIters", solver->GetIter());
                    tJSON.put("Integral.Solver.Outer.Residual", solver->GetResid());
                    tJSON.put("Integral.Solver.Outer.Converged", solver->GetResid() <= inpJSON.get<double>("Solver.Tol"));
                }
                if (i > 1) {
                    tJSON.put("ApaptiveTimeStep.NumbOfTries", numbOfTries);
                    tJSON.put("ApaptiveTimeStep.RelativeError", e);
                }
                if (everyStep > 0 && i % everyStep == 0) {
                    vtkExported = true;
                    logger.beg("write vtk");
                        writeVTK(t);
                    auto vtkTime = logger.end();
                    tJSON.put("ElapsedTime.VTK", vtkTime);
                }
                else vtkExported = false;
                stats << tJSON;
                logger.buf << tJSON;
                logger.log();
            };
        exportStats(0);
        logger.beg("do timestep i = 1 w/ BDF1");
            logger.beg("assemble");
                chi_prev = chi;
                SetupCahnHilliardIF_P1P1(mg, &Mass, &Normal_stab, &Tangent_stab, &Volume_stab, &Laplace, &LaplaceM, &Gprimeprime, lset.Phi, lset.GetBndData(), v, vbnd, chi_prev, chibnd);
                // Precond3.LinComb(sigm, LaplaceM.Data, alpha, Volume_stab.Data);
            logger.end();
            logger.beg("update time");
                t = dt;
                logger.buf
                    << "t  = " << t << '\n'
                    << "dt = " << dt;
                logger.log();
                InitScalar(mg, chiSol, surfCahnHilliardData.chi, t);
                InitScalar(mg, omegaSol, surfCahnHilliardData.omega, t);
                SetupInterfaceRhsP1(mg, &rhs3, lset.Phi, lset.GetBndData(), surfCahnHilliardData.rhs3, t);
                SetupInterfaceRhsP1(mg, &rhs4, lset.Phi, lset.GetBndData(), surfCahnHilliardData.rhs4, t);
            logger.end();
            logger.beg("BDF1 step");
                auto doBDF1 = [&](VecDescCL& chi) {
                    Precond3.LinComb(1. / dt, Mass.Data, alpha, Volume_stab.Data, sigm, Laplace.Data);
                    for (int i = 0; i < well_potential.Data.size(); i++)
                        well_potential.Data[i] = chemicalPotential(chi_prev.Data[i]);
                    instantrhs3 = rhs3;
                    instantrhs3.Data += (1. / dt) * (Mass.Data * chi_prev.Data);
                    instantrhs4 = rhs4;
                    instantrhs4.Data += S * (Mass.Data * chi_prev.Data) - Mass.Data * well_potential.Data;
                    if (model == "CahnHilliard") {
                        A.LinComb(sigm, LaplaceM.Data, alpha, Volume_stab.Data);
                        B.LinComb(0., Laplace.Data, 1. / dt, Mass.Data);
                        C.LinComb(0., Laplace.Data, -eps, Mass.Data);
                        D.LinComb(eps * eps, Laplace.Data, S, Mass.Data, alpha * eps * eps, Volume_stab.Data);
                        CH_solver->Solve(A, B, C, D, omega.Data, chi.Data, instantrhs3.Data, instantrhs4.Data, omega.RowIdx->GetEx(), chi.RowIdx->GetEx());
                    } else if (model == "AllenCahn") {
                        instantrhs3.Data += (S / eps) * (Mass.Data * chi_prev.Data) - (1. / eps) * (Mass.Data * well_potential.Data);
                        A.LinComb(sigm * eps, Laplace.Data, (S / eps + 1. / dt), Mass.Data, alpha * eps, Volume_stab.Data);
                        AC_solver->Solve(A, chi.Data, instantrhs3.Data, chi.RowIdx->GetEx());
                    }
                };
                doBDF1(chi);
                exportStats(1);
            logger.end();
            auto F_rho = P.get<double>("Time.Adaptive.rho");
            auto F_tol = P.get<double>("Time.Adaptive.Tol");
            auto F_min = P.get<double>("Time.Adaptive.MinStepSize");
            auto F = [=](double e, double dt) {
                return std::max(F_min, F_rho * std::sqrt(F_tol / e) * dt);
            };
            size_t i = 2; // second time step, apply BDF1/BDF2 apaptive scheme
        	while (t < T) {
        	    logger.beg("do timestep i = " + std::to_string(i) + " w/ BDF1/BDF2 adaptive scheme");
                    logger.beg("assemble");
                        chi_prev_prev = chi_prev;
                        chi_prev = chi;
                        chi_extrap.Data = 2. * chi_prev.Data - chi_prev_prev.Data;
                        SetupCahnHilliardIF_P1P1(mg, &Mass, &Normal_stab, &Tangent_stab, &Volume_stab, &Laplace, &LaplaceM, &Gprimeprime, lset.Phi, lset.GetBndData(), v, vbnd, chi_extrap, chibnd);
                        // Precond3.LinComb(sigm, LaplaceM.Data, alpha, Volume_stab.Data);
                    logger.end();
                    auto t_old = t;
                    do {
                        logger.beg("attempt #" + std::to_string(numbOfTries + 1));
                            logger.beg("update time");
                                if (numbOfTries > 0) dt = F(e, dt);
                                t = t_old + dt;
                                logger.buf
                                    << "t  = " << t << '\n'
                                    << "dt = " << dt;
                                logger.log();
                                InitScalar(mg, chiSol, surfCahnHilliardData.chi, t);
                                InitScalar(mg, omegaSol, surfCahnHilliardData.omega, t);
                                SetupInterfaceRhsP1(mg, &rhs3, lset.Phi, lset.GetBndData(), surfCahnHilliardData.rhs3, t);
                                SetupInterfaceRhsP1(mg, &rhs4, lset.Phi, lset.GetBndData(), surfCahnHilliardData.rhs4, t);
                            logger.end();
                            logger.beg("BDF1 step");
                                doBDF1(chi_BDF1);
                            logger.end();
                            logger.beg("BDF2 step");
                                Precond3.LinComb(1.5 / dt, Mass.Data, alpha, Volume_stab.Data, sigm, Laplace.Data);
                                for (int i = 0; i < well_potential.Data.size(); i++) well_potential.Data[i] = 2. * chemicalPotential(chi_prev.Data[i]) - chemicalPotential(chi_prev_prev.Data[i]);
                                instantrhs3 = rhs3;
                                instantrhs3.Data += (2. / dt) * (Mass.Data * chi_prev.Data) - (.5 / dt) * (Mass.Data * chi_prev_prev.Data);
                                instantrhs4 = rhs4;
                                instantrhs4.Data += S * (Mass.Data * chi_extrap.Data) - Mass.Data * well_potential.Data;
                                if (model == "CahnHilliard") {
                                    B.LinComb(0., Laplace.Data, 1.5 / dt, Mass.Data);
                                    CH_solver->Solve(A, B, C, D, omega.Data, chi_BDF2.Data, instantrhs3.Data, instantrhs4.Data, omega.RowIdx->GetEx(), chi_BDF2.RowIdx->GetEx());
                                } else if (model == "AllenCahn") {
                                    instantrhs3.Data += (S / eps) * (Mass.Data * chi_extrap.Data) - (1. / eps) * (Mass.Data * well_potential.Data);
                                    A.LinComb(sigm * eps, Laplace.Data, (S / eps + 1.5 / dt), Mass.Data, alpha * eps, Volume_stab.Data);
                                    AC_solver->Solve(A, chi_BDF2.Data, instantrhs3.Data, chi_BDF2.RowIdx->GetEx());
                                }
                            logger.end();
                            e = std::sqrt(norm_sq(chi_BDF2.Data - chi_BDF1.Data) / norm_sq(chi_BDF2.Data));
                            logger.buf << "e = " << e;
                            logger.log();
                            ++numbOfTries;
                        logger.end();
                    } while (e > F_tol && dt != F_min);
                    chi = chi_BDF2; // save soln as BDF2
                    // export
                    exportStats(i);
                    // prepare for the next step
                    dt = F(e, dt);
                    ++i;
                    numbOfTries = 0;
                logger.end();
        }
        if (everyStep > 0 && !vtkExported) { // make sure to export final time
            logger.beg("write vtk (last time frame)");
            writeVTK(t);
            auto vtkTime = logger.end();
            // TODOLATER: update JSON
        }
        double average3( 0.), variation3( 0.);
        double average4( 0.), variation4( 0.);
        double averageglobal( 0.), variationglobal( 0.);
        std::stringstream Precondstream3copy( Precondstream3.str());
        std::stringstream Precondstream4copy( Precondstream4.str());
        std::stringstream Globalstreamcopy( Globalstream.str());
        //std::cout << "Globalstream: " << std::endl << Globalstream.str() << std::endl;
        //std::cout << "Precondstream: " << std::endl << Precondstream4.str() << std::endl;

        ComputeAverageIterations(Precondstream3, average3);
        ComputeAverageIterations(Precondstream4, average4);
        RightComputeAverageIterations(Globalstream, averageglobal);

        ComputeVariationFromAverageIterations(Precondstream3copy, average3, variation3);
        ComputeVariationFromAverageIterations(Precondstream4copy, average4, variation4);
        RightComputeVariationFromAverageIterations(Globalstreamcopy, averageglobal, variationglobal);
    delete &lset;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }

}