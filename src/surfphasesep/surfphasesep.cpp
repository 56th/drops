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
// belos (iterative solvers)
#include "BelosSolverFactory.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
// amesos (sparse-direct solvers)
#include "Amesos.h"
#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
// epetra (vectors and operators / matrices)
#include "../surfnavierstokes/Epetra_OperatorApply.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#ifdef HAVE_MPI
    #include "Epetra_MpiComm.h"
#else
    #include "Epetra_SerialComm.h"
#endif

using namespace DROPS;

void InitVecLaplace(const MultiGridCL& MG, LevelsetP2CL& lset, VecDescCL& rhs, VecDescCL& vSol, VecDescCL& pSol,
                    InstatVectorFunction f_rhs, InstatVectorFunction f_vsol, InstatScalarFunction f_psol, double t = 0.0)
{
    if( vSol.RowIdx->NumUnknownsEdge()) {
        SetupInterfaceVectorRhsP2(MG, &rhs, lset.Phi, lset.GetBndData(), f_rhs);
    } else {
        SetupInterfaceVectorRhsP1(MG, &rhs, lset.Phi, lset.GetBndData(), f_rhs, t);
    }
    InitScalar(MG, pSol, f_psol);
    InitVector(MG, vSol, f_vsol);
}


int main (int argc, char* argv[]) {
    srand(time(NULL));
    auto& logger = SingletonLogger::instance();
    try {
        logger.beg("mpi init");
            #ifdef HAVE_MPI
                    MPI_Init(&argc, &argv);
                    Epetra_MpiComm comm(MPI_COMM_WORLD);
                    logger.log("using Epetra_MpiComm");
            #else
                    Epetra_SerialComm comm;
                    logger.log("using Epetra_SerialComm");
            #endif
            auto myRank = comm.MyPID();
            if (myRank == 0) {
                logger.buf << "num proc = " << comm.NumProc();
                logger.log();
            }
        logger.end();
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
                  "mkdir " + dirName + "/stats"
              ).c_str());
              {
                  std::ofstream input(dirName + "/input.json");
                  input << inpJSON;
                  logger.buf << inpJSON;
                  logger.log();
              }
        logger.end();
        logger.beg("set up IC");
            auto surfCahnHilliardData = SurfCahnHilliardDataFactory(inpJSON);
            logger.buf << surfCahnHilliardData.description;
            logger.log();
        logger.end();
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
        std::auto_ptr<MGBuilderCL> builder( make_MGBuilder( P));
        MultiGridCL mg( *builder);
        const ParamCL::ptree_type* ch= 0;
        try {
            ch= &P.get_child( "Mesh.Periodicity");
        }
        catch (DROPSParamErrCL) {}
        if (ch)
            read_PeriodicBoundaries( mg, *ch);

        // adaptive mesh refinement based on level set function
        typedef DistMarkingStrategyCL InitMarkerT;
        auto meshCoarseLevel = std::max(0, P.get<int>("Mesh.AdaptRef.CoarsestLevel"));
        auto meshFineLevel = std::max(0, P.get<int>("Mesh.AdaptRef.FinestLevel"));
        if (meshFineLevel < meshCoarseLevel) meshFineLevel = meshCoarseLevel;
        InitMarkerT initmarker(surfCahnHilliardData.surface.phi, P.get<double>("Mesh.AdaptRef.Width"), meshCoarseLevel, meshFineLevel);
        //adap.set_marking_strategy( &initmarker );
        AdapTriangCL adap( mg, &initmarker );
        adap.MakeInitialTriang();
        adap.set_marking_strategy( 0 );
        // create level set
        InstatScalarFunction sigma (0);
        SurfaceTensionCL sf( sigma, 0);
        BndDataCL<double> lsbnd( 0);
        read_BndData( lsbnd, mg, P.get_child( "Levelset.BndData"));
        std::shared_ptr<LevelsetP2CL> lsetPtr(LevelsetP2CL::Create(mg, lsbnd, sf));
        auto& lset = *lsetPtr;
        lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
        lset.Phi.SetIdx( &lset.idx);
        lset.Init(surfCahnHilliardData.surface.phi);
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
        double h = P.get<Point3DCL>("Mesh.E1")[0]/P.get<double>("Mesh.N1")*std::pow(2., -P.get<double>("Mesh.AdaptRef.FinestLevel"));
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
        BndDataCL<double> pbnd( 0);
        read_BndData( pbnd, mg, P.get_child( "Stokes.PressureBndData"));
        BndDataCL<double> chibnd( 0);
        read_BndData( chibnd, mg, P.get_child( "Stokes.VolumeFractionBndData"));
        BndDataCL<double> omegabnd( 0);
        read_BndData( omegabnd, mg, P.get_child( "Stokes.ChemPotentialBndData"));
        IdxDescCL ifaceVecP1idx( vecP1IF_FE, vbnd);
        IdxDescCL ifaceP1idx( P1IF_FE, pbnd);
        IdxDescCL vecP1idx( vecP1_FE, vbnd);
        IdxDescCL P1FEidx( P1_FE, pbnd);
        ifaceVecP1idx.GetXidx().SetBound( P.get<double>("SurfTransp.OmitBound"));
        ifaceVecP1idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
        ifaceP1idx.GetXidx().SetBound( P.get<double>("SurfTransp.OmitBound"));
        ifaceP1idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
        vecP1idx.GetXidx().SetBound( P.get<double>("SurfTransp.OmitBound"));
        vecP1idx.CreateNumbering(mg.GetLastLevel(), mg);
        P1FEidx.GetXidx().SetBound( P.get<double>("SurfTransp.OmitBound"));
        P1FEidx.CreateNumbering(mg.GetLastLevel(), mg);
        auto coarseLevel = P.get<int>("SurfCahnHilliard.IC.ProlongateFromLevelNo");
        if (coarseLevel < meshCoarseLevel) coarseLevel = meshCoarseLevel;
        if (coarseLevel > meshFineLevel)   coarseLevel = meshFineLevel;
        IdxDescCL P1FEidxCoarse(P1_FE, pbnd);
        P1FEidxCoarse.GetXidx().SetBound(P.get<double>("SurfTransp.OmitBound"));
        P1FEidxCoarse.CreateNumbering(coarseLevel, mg);
        // construct FE vectors (initialized with zero)
        VecDescCL v, omega, omegaSol, chi, chi_ext, chi_coarse_ini_ext, chi_prev, chi_prev_prev, chi_BDF1, chi_BDF2, chi_extrap, chiSol, energy, well_potential, f1, rhs1, f2, rhs2;
        if( !FE.compare("P1P1")) {
             v.SetIdx( &ifaceVecP1idx);
             chi.SetIdx( &ifaceP1idx);
             chi_coarse_ini_ext.SetIdx(&P1FEidxCoarse);
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
             f1.SetIdx(&ifaceP1idx);
             rhs1.SetIdx(&ifaceP1idx);
             f2.SetIdx(&ifaceP1idx);
             rhs2.SetIdx(&ifaceP1idx);
        }
        // setup matrices
        MatDescCL Laplace, Mass, Normal_stab, Tangent_stab, Volume_stab, LaplaceM, Gprimeprime;
        MatrixCL A, B, C, D;
        size_t n_c, n_omega;
        if( !FE.compare("P1P1")) {
            n_c = n_omega = ifaceP1idx.NumUnknowns();
            Laplace.SetIdx( &ifaceP1idx, &ifaceP1idx);
            LaplaceM.SetIdx( &ifaceP1idx, &ifaceP1idx);
            Gprimeprime.SetIdx( &ifaceP1idx, &ifaceP1idx);
            Mass.SetIdx(&ifaceP1idx, &ifaceP1idx);
            Normal_stab.SetIdx(&ifaceP1idx, &ifaceP1idx);
            Tangent_stab.SetIdx(&ifaceP1idx, &ifaceP1idx);
            Volume_stab.SetIdx(&ifaceP1idx, &ifaceP1idx);
            SetupCahnHilliardIF_P1P1(mg,  &Mass,&Normal_stab, &Tangent_stab, &Volume_stab, &Laplace, &LaplaceM, &Gprimeprime, lset.Phi, lset.GetBndData(), v, vbnd, chi, chibnd);
        }
        logger.beg("set up vtk");
            VTKWriter vtkWriter(dirName + "/vtk/separation", mg, inpJSON.get<bool>("Output.Binary"));
            vtkWriter
                .add(VTKWriter::VTKVar({"level-set", &lset.Phi.Data, VTKWriter::VTKVar::Type::P2}))
                .add(VTKWriter::VTKVar({"c_h", &chi_ext.Data, VTKWriter::VTKVar::Type::P1}));
            auto everyStep = inpJSON.get<int>("Output.EveryStep");
            auto writeVTK = [&](double t) {
                Extend(mg, chi, chi_ext);
                vtkWriter.write(t);
            };
        logger.end();
        logger.beg("set up outer solver");
            using MV = Epetra_MultiVector;
            using OP = Epetra_Operator;
            using ST = double;
            using namespace Teuchos;
            using namespace Belos;
            RCP<ParameterList> belosParams = parameterList();
            belosParams->set("Num Blocks", inpJSON.get<int>("Solver.Outer.KrylovSubspaceSize"));
            belosParams->set("Maximum Iterations", inpJSON.get<int>("Solver.Outer.MaxIter"));
            belosParams->set("Convergence Tolerance", inpJSON.get<double>("Solver.Outer.RelResTol"));
            belosParams->set( "Output Frequency", inpJSON.get<int>("Solver.Outer.OutputFrequency"));
            belosParams->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails + Belos::TimingDetails + Belos::FinalSummary + Belos::IterationDetails);
            belosParams->set<int>("Output Style", Belos::Brief);
            SolverFactory<ST, MV, OP> belosFactory;
            logger.buf << "available iterations: " << belosFactory.supportedSolverNames();
            logger.log();
            auto belosSolver = belosFactory.create(inpJSON.get<std::string>("Solver.Outer.Iteration"), belosParams);
            logger.buf << "outer solver: " << belosSolver->description();
            logger.log();
            // matrix, lhs, and rhs
            RCP<OP> belosMTX;
            Epetra_Map belosMap(static_cast<int>(n_c + n_omega), 0, comm);
            Epetra_Vector belosLHS(belosMap), belosRHS(belosMap);
        logger.end();
        logger.beg("do timestep i = 0");
        	auto t = 0.;
        	auto e = 0.;
            VectorCL unityVector(1., n_c);
            {
                IdxDescCL P1FEidxFrom(P1_FE, pbnd), P1FEidxTo(P1_FE, pbnd);
                P1FEidxFrom.GetXidx().SetBound(P.get<double>("SurfTransp.OmitBound"));
                P1FEidxTo.GetXidx().SetBound(P.get<double>("SurfTransp.OmitBound"));
                std::vector<ProlongationCL<double>> P;
                std::vector<ProlongationCL<double>> p;
                for (size_t l = coarseLevel + 1; l <= meshFineLevel; ++l) {
                    P.emplace_back(ProlongationCL<double>(mg));
                    P1FEidxFrom.CreateNumbering(l - 1, mg);
                    P1FEidxTo.CreateNumbering(l, mg);
                    P.back().Create(&P1FEidxFrom, &P1FEidxTo);
                }
                double raftFraction, raftFractionError;
                do {
                    InitScalar(mg, chi_coarse_ini_ext, surfCahnHilliardData.chi, 0.);
                    chi_ext.Data = chi_coarse_ini_ext.Data;
                    for (auto const &p : P)
                        chi_ext.Data = p * chi_ext.Data;
                    Restrict(mg, chi_ext, chi);
                    if (surfCahnHilliardData.raftRatio > 0.) {
                        raftFraction = dot(Mass.Data * chi.Data, unityVector) / dot(Mass.Data * unityVector, unityVector);
                        raftFractionError = std::fabs(surfCahnHilliardData.raftRatio - raftFraction) / surfCahnHilliardData.raftRatio;
                        logger.buf << "raft ratio error (%) = " << raftFractionError;
                        logger.log();
                    }
                } while (surfCahnHilliardData.raftRatio > 0. && raftFractionError > .001);
            }
            chiSol = chi;
            size_t numbOfTries = 0;
            auto vtkExported = false;
            auto solveTime = 0.;
            auto factorizationTime = 0.;
            auto belosSolverResult = Belos::Converged;
            auto exportStats = [&](size_t i) {
                std::ofstream stats(dirName + "/stats/t_" + std::to_string(i) + ".json");
                ParamCL tJSON;
                tJSON.put("DOF.c", n_c);
                tJSON.put("DOF.mu", n_omega);
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
                    tJSON.put("Integral.Solver.Outer.TotalIters", belosSolver->getNumIters());
                    tJSON.put("Integral.Solver.Outer.ResidualNorm.SolverRelative", belosSolver->achievedTol());
                    tJSON.put("Integral.Solver.Outer.Converged", belosSolverResult == Belos::Converged);
                    tJSON.put("ElapsedTime.Factorization", factorizationTime);
                    factorizationTime = 0.;
                    tJSON.put("ElapsedTime.LinearSolve", solveTime);
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
        logger.end();
        if (T > 0.) {
            logger.beg("do timestep i = 1 w/ BDF1");
                logger.beg("assemble");
                    chi_prev = chi;
                    SetupCahnHilliardIF_P1P1(mg, &Mass, &Normal_stab, &Tangent_stab, &Volume_stab, &Laplace, &LaplaceM, &Gprimeprime, lset.Phi, lset.GetBndData(), v, vbnd, chi_prev, chibnd);
                logger.end();
                logger.beg("update time");
                    t = dt;
                    logger.buf
                        << "t  = " << t << '\n'
                        << "dt = " << dt;
                    logger.log();
                    InitScalar(mg, chiSol, surfCahnHilliardData.chi, t);
                    InitScalar(mg, omegaSol, surfCahnHilliardData.omega, t);
                    SetupInterfaceRhsP1(mg, &f1, lset.Phi, lset.GetBndData(), surfCahnHilliardData.rhs3, t);
                    SetupInterfaceRhsP1(mg, &f2, lset.Phi, lset.GetBndData(), surfCahnHilliardData.rhs4, t);
                logger.end();
                logger.beg("BDF1 step");
                    logger.beg("set up blocks");
                        for (int i = 0; i < well_potential.Data.size(); i++)
                            well_potential.Data[i] = chemicalPotential(chi_prev.Data[i]);
                        rhs1 = f1;
                        rhs1.Data += (1. / dt) * (Mass.Data * chi_prev.Data);
                        rhs2 = f2;
                        rhs2.Data += S * (Mass.Data * chi_prev.Data) - Mass.Data * well_potential.Data;
                        A.LinComb(0., Laplace.Data, 1. / dt, Mass.Data);
                        B.LinComb(sigm, LaplaceM.Data, alpha, Volume_stab.Data);
                        C.LinComb(eps * eps, Laplace.Data, S, Mass.Data, alpha * eps * eps, Volume_stab.Data);
                        D.LinComb(0., Laplace.Data, -1., Mass.Data);
                        MatrixCL ABCD(A, B, C, D);
                        /* A.exportMAT(dirName + "/matrices/A.mat");
                        B.exportMAT(dirName + "/matrices/B.mat");
                        C.exportMAT(dirName + "/matrices/C.mat");
                        D.exportMAT(dirName + "/matrices/D.mat");
                        ABCD.exportMAT(dirName + "/matrices/ABCD.mat"); */
                    logger.end();
                    logger.beg("convert to Epetra");
                        auto ABCD_Epetra = static_cast<Epetra_CrsMatrix>(ABCD);
                        auto printStat = [&](std::string const & name, Epetra_CrsMatrix const & A) {
                            logger.buf << name << ": " << A.NumGlobalRows() << 'x' << A.NumGlobalCols() << ", " << A.NumGlobalNonzeros() << " nonzeros (" << (100. * A.NumGlobalNonzeros()) / (static_cast<double>(A.NumGlobalRows()) * A.NumGlobalCols()) << "%)";
                            logger.log();
                        };
                        printStat("{A, B; C, D} block mtx", ABCD_Epetra);
                        belosMTX = rcpFromRef(ABCD_Epetra);
                        for (size_t i = 0; i < n_c; ++i) belosRHS[i] = rhs1.Data[i];
                        for (size_t i = 0; i < n_omega; ++i) belosRHS[i + n_c] = rhs2.Data[i];
                    logger.end();
                    RCP<OP> belosPRE;
                    Epetra_LinearProblem amesosProblem;
                    RCP<Amesos_BaseSolver> amesosSolver;
                    Amesos amesosFactory;
                    std::function<void()> runFactorization = [](){};
                    if (inpJSON.get<bool>("Solver.Inner.Use")) {
                        logger.beg("set up preconditioner");
                            amesosProblem.SetOperator(&ABCD_Epetra);
                            amesosSolver = rcp(amesosFactory.Create("Amesos_Klu", amesosProblem));
                            belosPRE = rcp(new Epetra_OperatorApply([&](MV const &X, MV &Y) {
                                amesosProblem.SetLHS(&Y);
                                amesosProblem.SetRHS(const_cast<MV*>(&X));
                                amesosSolver->Solve();
                            }));
                            runFactorization = [&]() {
                                logger.beg("factorization");
                                    logger.beg("symbolic factorization");
                                        amesosSolver->SymbolicFactorization();
                                    logger.end();
                                    logger.beg("numeric factorization");
                                        amesosSolver->NumericFactorization();
                                    logger.end();
                                factorizationTime = logger.end();
                            };
                        logger.end();
                    }
                    runFactorization();
                    logger.beg("linear solve");
                        belosLHS.PutScalar(0.);
                        LinearProblem<ST, MV, OP> belosProblem(belosMTX, rcpFromRef(belosLHS), rcpFromRef(belosRHS));
                        if (inpJSON.get<bool>("Solver.Inner.Use")) belosProblem.setRightPrec(belosPRE);
                        belosProblem.setProblem();
                        belosSolver->setProblem(rcpFromRef(belosProblem));
                        std::cout << std::scientific;
                        belosSolverResult = belosSolver->solve();
                        if (myRank == 0) {
                            if (belosSolverResult == Belos::Converged) logger.log("belos converged");
                            else logger.wrn("belos did not converge");
                        }
                    solveTime = logger.end();
                    logger.beg("convert from Epetra");
                        for (size_t i = 0; i < n_c; ++i) chi.Data[i] = belosLHS[i];
                        for (size_t i = 0; i < n_omega; ++i) omega.Data[i] = belosLHS[i + n_c];
                    logger.end();
                    exportStats(1);
                logger.end();
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
                                SetupInterfaceRhsP1(mg, &f1, lset.Phi, lset.GetBndData(), surfCahnHilliardData.rhs3, t);
                                SetupInterfaceRhsP1(mg, &f2, lset.Phi, lset.GetBndData(), surfCahnHilliardData.rhs4, t);
                            logger.end();
                            logger.beg("BDF1 step");
                                logger.beg("set up blocks");
                                    for (int i = 0; i < well_potential.Data.size(); i++)
                                        well_potential.Data[i] = chemicalPotential(chi_prev.Data[i]);
                                    rhs1 = f1;
                                    rhs1.Data += (1. / dt) * (Mass.Data * chi_prev.Data);
                                    rhs2 = f2;
                                    rhs2.Data += S * (Mass.Data * chi_prev.Data) - Mass.Data * well_potential.Data;
                                    A.LinComb(0., Laplace.Data, 1. / dt, Mass.Data);
                                    B.LinComb(sigm, LaplaceM.Data, alpha, Volume_stab.Data);
                                    C.LinComb(eps * eps, Laplace.Data, S, Mass.Data, alpha * eps * eps, Volume_stab.Data);
                                    D.LinComb(0., Laplace.Data, -1., Mass.Data);
                                    ABCD = MatrixCL(A, B, C, D);
                                logger.end();
                                logger.beg("convert to Epetra");
                                    ABCD_Epetra = static_cast<Epetra_CrsMatrix>(ABCD);
                                    printStat("{A, B; C, D} block mtx", ABCD_Epetra);
                                    belosMTX = rcpFromRef(ABCD_Epetra);
                                    for (size_t i = 0; i < n_c; ++i) belosRHS[i] = rhs1.Data[i];
                                    for (size_t i = 0; i < n_omega; ++i) belosRHS[i + n_c] = rhs2.Data[i];
                                logger.end();
                                logger.beg("linear solve");
                                    belosLHS.PutScalar(0.);
                                    LinearProblem<ST, MV, OP> belosProblemBDF1(belosMTX, rcpFromRef(belosLHS), rcpFromRef(belosRHS));
                                    if (inpJSON.get<bool>("Solver.Inner.Use")) belosProblemBDF1.setRightPrec(belosPRE);
                                    belosProblemBDF1.setProblem();
                                    belosSolver->setProblem(rcpFromRef(belosProblemBDF1));
                                    std::cout << std::scientific;
                                    belosSolverResult = belosSolver->solve();
                                    if (myRank == 0) {
                                        if (belosSolverResult == Belos::Converged) logger.log("belos converged");
                                        else logger.wrn("belos did not converge");
                                    }
                                solveTime = logger.end();
                                if (belosSolverResult != Belos::Converged) {
                                    runFactorization();
                                    logger.beg("linear solve w/ new factorization");
                                        belosLHS.PutScalar(0.);
                                        belosSolverResult = belosSolver->solve();
                                        if (belosSolverResult == Belos::Converged) logger.log("belos converged");
                                        else logger.wrn("belos did not converge");
                                    solveTime += logger.end();
                                }
                                logger.beg("convert from Epetra");
                                    for (size_t i = 0; i < n_c; ++i) chi_BDF1.Data[i] = belosLHS[i];
                                logger.end();
                            logger.end();
                            logger.beg("BDF2 step");
                                logger.beg("set up blocks");
                                    for (int i = 0; i < well_potential.Data.size(); i++)
                                        well_potential.Data[i] = 2. * chemicalPotential(chi_prev.Data[i]) - chemicalPotential(chi_prev_prev.Data[i]);
                                    rhs1 = f1;
                                    rhs1.Data += (2. / dt) * (Mass.Data * chi_prev.Data) - (.5 / dt) * (Mass.Data * chi_prev_prev.Data);
                                    rhs2 = f2;
                                    rhs2.Data += S * (Mass.Data * chi_extrap.Data) - Mass.Data * well_potential.Data;
                                    A.LinComb(0., Laplace.Data, 1.5 / dt, Mass.Data);
                                    ABCD = MatrixCL(A, B, C, D);
                                logger.end();
                                logger.beg("convert to Epetra");
                                    ABCD_Epetra = static_cast<Epetra_CrsMatrix>(ABCD);
                                    printStat("{A, B; C, D} block mtx", ABCD_Epetra);
                                    belosMTX = rcpFromRef(ABCD_Epetra);
                                    for (size_t i = 0; i < n_c; ++i) belosRHS[i] = rhs1.Data[i];
                                    for (size_t i = 0; i < n_omega; ++i) belosRHS[i + n_c] = rhs2.Data[i];
                                logger.end();
                                logger.beg("linear solve");
                                    belosLHS.PutScalar(0.);
                                    LinearProblem<ST, MV, OP> belosProblemBDF2(belosMTX, rcpFromRef(belosLHS), rcpFromRef(belosRHS));
                                    if (inpJSON.get<bool>("Solver.Inner.Use")) belosProblemBDF2.setRightPrec(belosPRE);
                                    belosProblemBDF2.setProblem();
                                    belosSolver->setProblem(rcpFromRef(belosProblemBDF2));
                                    std::cout << std::scientific;
                                    belosSolverResult = belosSolver->solve();
                                    if (myRank == 0) {
                                        if (belosSolverResult == Belos::Converged) logger.log("belos converged");
                                        else logger.wrn("belos did not converge");
                                    }
                                solveTime += logger.end();
                                if (belosSolverResult != Belos::Converged) {
                                    runFactorization();
                                    logger.beg("linear solve w/ new factorization");
                                        belosLHS.PutScalar(0.);
                                        belosSolverResult = belosSolver->solve();
                                        if (belosSolverResult == Belos::Converged) logger.log("belos converged");
                                        else logger.wrn("belos did not converge");
                                    solveTime += logger.end();
                                }
                                logger.beg("convert from Epetra");
                                    for (size_t i = 0; i < n_c; ++i) chi_BDF2.Data[i] = belosLHS[i];
                                logger.end();
                            logger.end();
                            e = std::sqrt(norm_sq(chi_BDF2.Data - chi_BDF1.Data) / norm_sq(chi_BDF2.Data));
                            logger.buf << "e = " << e;
                            logger.log();
                            ++numbOfTries;
                        logger.end();
                    } while (e > F_tol && dt != F_min);
                    logger.beg("convert $\\omega$ from Epetra");
                        for (size_t i = 0; i < n_omega; ++i) omega.Data[i] = belosLHS[i + n_c];
                    logger.end();
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
        }
        #ifdef HAVE_MPI
                MPI_Finalize() ;
        #endif
        return 0;
    } catch (std::exception const & e) {
        logger.err(e.what());
    } catch (DROPSErrCL const & e) {
        e.what(logger.buf);
        logger.err(logger.buf.str());
    } catch (...) {
        logger.err("unknown error");
    }
    return 1;
}