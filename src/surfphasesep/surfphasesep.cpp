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
#include "out/vtkOut.h"
#include "num/bndData.h"
#include "num/precond.h"
#include "parallel/exchange.h"


#include "num/oseensolver.h"
#include <fstream>

#include "surfphasesep/surfphasesep_funcs.h"
#include "surfnavierstokes/surfnavierstokes_utils.h"
#include "surfnavierstokes/surfnavierstokes_tests.h"

using namespace DROPS;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Model Test Cases ////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void InitVecLaplace(const MultiGridCL& MG, LevelsetP2CL& lset, DROPS::VecDescCL& rhs, DROPS::VecDescCL& vSol, DROPS::VecDescCL& pSol,
                   instat_vector_fun_ptr f_rhs, instat_vector_fun_ptr f_vsol, instat_scalar_fun_ptr f_psol, double t = 0.0)
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

  try {

	DROPS::read_parameter_file_from_cmdline( P, argc, argv, "../../param/surfphasesep/No_Bnd_Condition_ch.json");
    std::cout << P << std::endl;
    DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"), P.get<std::vector<std::string> >("General.DynamicLibs") );

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


    // choose level set
    instat_scalar_fun_ptr levelset_fun;
    std::string levelset_fun_str = P.get<std::string>("Levelset.case");

    if( !levelset_fun_str.compare("sphere_2")) {
        levelset_fun = &sphere_2;
        std::cout << "The levelset is the unit sphere." << std::endl;
    } else  if( !levelset_fun_str.compare("tamarind")) {
        levelset_fun = &tamarind;
        std::cout << "The levelset is the tamarind." << std::endl;
    }
    else  if( !levelset_fun_str.compare("spindle")) {
        levelset_fun = &spindle;
        std::cout << "The levelset is the spindle." << std::endl;
    }

    // adaptive mesh refinement based on level set function
    typedef DROPS::DistMarkingStrategyCL InitMarkerT;
    InitMarkerT initmarker( levelset_fun, P.get<double>("Mesh.AdaptRef.Width"), P.get<int>("Mesh.AdaptRef.CoarsestLevel"), P.get<int>("Mesh.AdaptRef.FinestLevel") );
    //adap.set_marking_strategy( &initmarker );
      DROPS::AdapTriangCL adap( mg, &initmarker );

      adap.MakeInitialTriang();
    adap.set_marking_strategy( 0 );

    // create level set
    instat_scalar_fun_ptr sigma (0);
    SurfaceTensionCL sf( sigma, 0);

    BndDataCL<double> lsbnd( 0);
    read_BndData( lsbnd, mg, P.get_child( "Levelset.BndData"));

    DROPS::LevelsetP2CL & lset( * DROPS::LevelsetP2CL::Create( mg, lsbnd, sf) );

    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
//    LinearLSInit( mg, lset.Phi, levelset_fun);
    lset.Init( levelset_fun);

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

    double tau = P.get<double>("Time.StepSize");
    if (P.get<std::string>("SurfCahnHilliard.instationary") == "none") tau=1;

    double eta_order=-2.0;
    double epsilon_order=1.0;
    //double alpha_order=1.0;
    double alpha_order=1.0;
    double eta 		   = 1.e0  * pow(h, eta_order);//std::pow(2.e0,eta_index);; //constant for tangential penalty
    double epsilon     = 1.e0  * pow(h, epsilon_order); //constant for velocity stabilisation
    double alpha       = 1.e0  * pow(h, alpha_order); //constant for volume stabilisation
    double rho         = 1.e0 * pow(h, 1); //constant for Schur complement preconditioner

    std::cout << "h is: " << h << std::endl;
    std::cout << "tau is: " << tau << std::endl;
    ParameterNS::h = h;

      double S=0.29;

      double sigm = P.get<double>("SurfCahnHilliard.mobility");
    double eps = P.get<double>("SurfCahnHilliard.epsilon");
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
     DROPS::VecDescCL v, omega, omegaSol, chi,chi_old, chi_new, chiSol, energy, well_potential, well_potential_concave,  well_potential_convex, chiInit, rhs3, instantrhs3, rhs4, instantrhs4;
     if( !FE.compare("P1P1")) {
    	 v.SetIdx( &ifaceVecP1idx);
    	 chi.SetIdx( &ifaceP1idx);
    	 chi_old.SetIdx( &ifaceP1idx);
         chi_new.SetIdx( &ifaceP1idx);
         well_potential.SetIdx( &ifaceP1idx);
         energy.SetIdx( &ifaceP1idx);
         well_potential_convex.SetIdx( &ifaceP1idx);
         well_potential_concave.SetIdx( &ifaceP1idx);
         chiInit.SetIdx( &ifaceP1idx);
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
        {
            Precond3.LinComb(1.0/tau, Mass.Data, rho, Volume_stab.Data, sigm, Laplace.Data);
            Precond4.LinComb(1.0, Mass.Data, rho*eps*eps, Volume_stab.Data, eps*eps, Laplace.Data);
        }
        else if (model=="AllenCahn")
        {
            Precond3.LinComb(1.0/tau, Mass.Data, rho, Volume_stab.Data, sigm, Laplace.Data);
            Precond4.LinComb(1.0    , Mass.Data, rho, Volume_stab.Data,   eps, Laplace.Data);        }


    }

    // construct preconditioners

    /*typedef GMResSolverCL<DummyPcCL> GMResSolverT;

    std::stringstream PCstream;

    DummyPcCL Ident;

    GMResSolverT GMResSolver_( Ident, P.get<int>("Solver.PcAIter"), P.get<int>("Solver.PcAIter"), P.get<double>("Solver.PcATol"),
    		 *//*bool relative=*//* true, *//*bool calculate2norm=*//* false, *//*PreMethGMRES method=*//* LeftPreconditioning,
			 *//*bool mod =*//* true,      *//*bool useModGS =*//* false, &PCstream);*/

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

      /*typedef SolverAsPreCL<PCGSolverT> PCGPcT;
      PCGPcT PCGPc_( PCGSolver3);
      typedef BlockPreCL<ExpensivePreBaseCL, SchurPreBaseCL, DiagSpdBlockPreCL>  DiagBlockEPcT;
      typedef PLanczosONBCL<VectorCL, DiagBlockEPcT> LanczosT;
      typedef PMResSolverCL<LanczosT> MinResT;
      DiagBlockEPcT    *DBlock_ = new DiagBlockEPcT( PCGPc_, *spc4_);
      LanczosT *lanczos_ = new LanczosT( *DBlock_);
      MinResT *MinRes_ = new MinResT( *lanczos_,  P.get<int>("Solver.Iter"), P.get<double>("Solver.Tol"), *//*relative*//* false);
      BlockMatrixSolverCL<MinResT> *stokessolver= new BlockMatrixSolverCL<MinResT>( *MinRes_);*/

      SolverBaseCL *solver;
      if (model=="CahnHilliard")
      {
            solver=CH_solver;
      }
      else if (model=="AllenCahn")
      {
          solver=AC_solver;
      }

    // set function pointers and rhs vectors for different test cases
    DROPS::instat_scalar_fun_ptr extchisol = &ZeroScalarFun,
    							 extrhs3=&ZeroScalarFun,
								 extrhs4=&ZeroScalarFun,
								 extomegasol=&ZeroScalarFun;

    if( !levelset_fun_str.compare("sphere_2")) {
        if( !testcase.compare("1")) {
            std::cout << "Test case 1: mobility=const, stationary omega=harmonic, chi=harmonic/2" << std::endl;
            extchisol = &Test1_chiSol;
            extomegasol = &Test1_omegaSol;
            extrhs3 = &Test1_rhs3;
        }
        if( !testcase.compare("2")) {
                    std::cout << "Test case 2: mobility=const, instationary omega=(1-exp(-4t))*harmonic, chi=(1-exp(-4t))*harmonic/2" << std::endl;
                    extchisol = &Test2_chiSol;
                    extomegasol = &Test2_omegaSol;
                    extrhs3 = &Test2_rhs3;

        }
        if( !testcase.compare("3")) {
                           std::cout << "Test case 3: mobility=nonlinear, instationary omega=(1-exp(-4t))*harmonic, chi=(1-exp(-4t))*harmonic/2 " << std::endl;
                           extchisol = &Test3_chiSol;
                           extomegasol = &Test3_omegaSol;
                           extrhs3 = &Test3_rhs3;
                           extrhs4 = &Test3_rhs4;

        }
        if( !testcase.compare("4")) {
            std::cout << "Test case 4: mobility=nonlinear and strictly positive, instationary omega=(1-exp(-4t))*harmonic, chi=(1-exp(-4t))*harmonic/2  " << std::endl;
            extchisol = &Test3_chiSol;
            extomegasol = &Test3_omegaSol;
            extrhs3 = &Test4_rhs3;
            extrhs4 = &Test3_rhs4;

        }
        if( !testcase.compare("5")) {
            std::cout << "Test case 5: Allen-Cahn with space-uniform solution" << std::endl;
            extchisol = &Test5_chiSol;
        }
        if( !testcase.compare("6")) {
            std::cout << "Test case 6: Cahn-Hilliard with decay*harmonic initial condition" << std::endl;
            extchisol = &Test6_chiSol;
            extomegasol = &Test6_omegaSol;
            extrhs3 = &Test6_rhs3;
        }
        if( !testcase.compare("7")) {
            std::cout << "Test case 7: Allen-Cahn with decay*harmonic initial condition" << std::endl;
            extchisol = &Test7_chiSol;
            extrhs3 = &Test7_rhs3;

        }
        if( !testcase.compare("8")) {
            std::cout << "Test case 8: Allen-Cahn with decay*harmonic initial condition but linear well potential" << std::endl;
            extchisol = &Test7_chiSol;
            extrhs3 = &Test8_rhs3;

        }
        if( !testcase.compare("9")) {
            std::cout << "Test case 9: Allen-Cahn with random initial condition " << std::endl;
            extchisol = &Test9_chiSol;

        }
        if( !testcase.compare("10")) {
            std::cout << "Test case 10: Cahn-Hilliard with decay*second-order initial condition" << std::endl;
            extchisol = &Test10_chiSol;
            extomegasol = &Test10_omegaSol;
            extrhs3 = &Test10_rhs3;

        }
        if( !testcase.compare("11")) {
            std::cout << "Test case 11: Allen-Cahn with decay*second-order initial condition" << std::endl;
            extchisol = &Test11_chiSol;
            extrhs3 = &Test11_rhs3;

        }
        if( !testcase.compare("12")) {
            std::cout << "Test case 12: Cahn-Hilliard with  harmonic initial condition" << std::endl;
            extchisol = &Test12_chiSol;
        }

    } else if( !levelset_fun_str.compare("tamarind")) {
        if (!testcase.compare("9")) {
            std::cout << "Test case 9: Allen-Cahn with random initial condition " << std::endl;
            extchisol = &Test9_chiSol;
        }
    }

    else if( !levelset_fun_str.compare("spindle")) {
            if (!testcase.compare("9")) {
                std::cout << "Test case 9: Allen-Cahn with random initial condition " << std::endl;
                extchisol = &Test9_chiSol;

            }
        }

/////////////////////////////////////// Cahn-Hilliard ///////////////////////////////////////

	if( !model.compare("CahnHilliard") || !model.compare("AllenCahn")) {
        VectorCL id;
        id.resize(rhs3.Data.size(), 1);
		//set up discrete vectors and matrices
        DROPS::VecDescCL omegaxtent, chixtent;

        //set up output
        std::string filename = "test" + testcase + "_";
       	std::string dirname  = P.get<std::string>("Output.Directory") + "/" + model + "_" + levelset_fun_str + "/" + "test" + testcase + "_h=" + std::to_string(float(h));

       	std::cout << "dirname: " << dirname << std::endl;

        std::ofstream log_global;

       	//set up VTK output
       	VTKOutCL * vtkwriter = NULL;
       	vtkwriter = new VTKOutCL(mg, "DROPS data", (int)P.get<double>("Time.NumSteps")/P.get<int>("Output.every timestep"), dirname , filename, "none", 1 );
        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "level-set") );
        vtkwriter->Register( make_VTKIfaceScalar(mg, chi, "volume fraction", pbnd));
        vtkwriter->Register( make_VTKIfaceScalar(mg, omega, "GL-potential", pbnd));
        vtkwriter->Register( make_VTKIfaceScalar(mg, chiSol, "volume fraction exact", pbnd));
        vtkwriter->Register( make_VTKIfaceScalar(mg, omegaSol, "GL-potential exact", pbnd));
        vtkwriter->Register( make_VTKIfaceVector(mg, v, "velocity", velFE, vbnd));
        vtkwriter->Register( make_VTKIfaceScalar(mg, rhs3, "rhs3", pbnd));
        vtkwriter->Register( make_VTKIfaceScalar(mg, rhs4, "rhs4", pbnd));


        InitScalar(mg, chiInit, extchisol, 0.);
        InitScalar(mg, chiSol, extchisol, 0.);

        //NAVIER-STOKES starts here
        if ( P.get<std::string>("SurfCahnHilliard.instationary") != "none" )
        {
        	BndDataCL<Point3DCL> bndvec = vbnd;
        	BndDataCL<double> bndscalar = pbnd;

        	double T0=P.get<double>("Time.Read");
        	if (T0!= 0)
            {
                std::ifstream dump(dirname + "/chi_dump"+"_t=" + std::to_string((float)(T0)) + ".dat");
                std::ifstream dump2(dirname + "/omega_dump"+"_t=" + std::to_string((float)(T0)) + ".dat");
                if (dump.is_open()) chiInit.Read(dump);
                if (dump2.is_open()) omega.Read(dump2);
                else std::cout<< "Couldn't find " <<  dirname + "/chi_dump"+"_t=" + std::to_string((float)(T0)) + ".dat"<<std::endl;

            }
        	//initial velocity for backward Euler formula
        	chi=chiInit;
        	chi_old=chiInit;

        	//output of initial data and exact solutions to vtk and custom
        	vtkwriter->Write(T0);//send v to vtk

        	//right constant for temporal BDF
        	/*double c;
        	if ( P.get<std::string>("SurfCahnHilliard.instationary") == "BDF1")
        	{
        		//leading time term mass coeff
        		c=1.0;
        	}
        	else if ( P.get<std::string>("SurfCahnHilliard.instationary") == "BDF2")
        	{
        		//leading time term mass coeff
        		c=3.0/2.0;
        	}
        	else
        	{
        		std::cout << "problem is instationary, pick BDF1 or BDF2" << std::endl;
        		return 0;
        	}
*/
            //set up a txt file for error time output
            std::ofstream log_plot( dirname +"/"
                                     + "Plot_" + filename
                                     + "l="  + P.get<std::string>("Mesh.AdaptRef.FinestLevel") + "_"
                                     + P.get<std::string>("SurfCahnHilliard.instationary") + "="
                                     + std::to_string(float(tau))
                                     + ".txt");

            log_global.open( dirname +"/"
                                     + "Global_" + filename
                                     + "l="  + P.get<std::string>("Mesh.AdaptRef.FinestLevel") + "_"
                                     + P.get<std::string>("SurfCahnHilliard.instationary") + "="
                                     + std::to_string(float(tau)) +
                                     ".txt");

            double Energy_L2_omega=0.0;
            double Energy_L2_chi=0.0;
            double C_norm_chi=0.0;
            double C_norm_omega=0.0;
            log_plot << "Time\t" << "L_2(omega)\t" << "L_2(chi)" << "L_2(omega-ext_omega)\t" << "L_2(chi-ext_chi)\t" << "Lyapunov_energy" << std::endl;

            if ( !prFE.compare("P1")) {
                log_plot <<  std::to_string((float)T0) << "\t";
                omegaxtent.SetIdx(&P1FEidx);
                Extend(mg, omega, omegaxtent);
                chixtent.SetIdx(&P1FEidx);
                Extend(mg, chi, chixtent);
                double L2_omega_error = L2_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, omegabnd, omegaxtent),
                                                 extomegasol, T0);
                double L2_chi_error = L2_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, chibnd, chixtent),
                                               extchisol, T0);
                double L2_omega = L2_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, omegabnd, omegaxtent),
                                           ZeroScalarFun, T0);
                double L2_chi = L2_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, chibnd, chixtent),
                                         ZeroScalarFun, T0);
                for (int i = 0; i < energy.Data.size(); i++) {
                    energy.Data[i] = Potential_function(chi.Data[i]);
                }
                double Lyapunov_energy = (eps / 2.) * dot(Laplace.Data * chi.Data, chi.Data) +
                                         (1. / eps) * dot(Mass.Data * energy.Data, id);
                log_plot << std::to_string((float) L2_omega) << "\t" << std::to_string((float) L2_chi) << "\t" <<
                         std::to_string((float) L2_omega_error) << "\t" << std::to_string((float) L2_chi_error) << "\t"
                         <<
                         std::to_string((float) Lyapunov_energy) <<
                         std::endl;
            }
        	for (int i = 0; i < P.get<double>("Time.NumSteps"); i++)
        	{
                double t=T0+(i+1)*tau;
                //right constant for temporal BDF
                double c;
                if ( P.get<std::string>("SurfCahnHilliard.instationary") == "BDF1")
                {
                    //leading time term mass coeff
                    c=1.0;
                    chi_new.Data=chi.Data;
                    for (int i=0; i<well_potential.Data.size();i++)
                    {
                        well_potential.Data[i]=Potential_prime_function(chi.Data[i]);
                    }
                }
                else if ( P.get<std::string>("SurfCahnHilliard.instationary") == "BDF2")
                {
                    //leading time term mass coeff
                    if (i!=0) {
                        c = 3.0 / 2.0;
                    }
                    else c=1.0;
                    chi_new.Data = 2.0*chi.Data + (-1.)*chi_old.Data;
                    for (int i=0; i<well_potential.Data.size();i++)
                    {
                        well_potential.Data[i]=2*Potential_prime_function(chi.Data[i]) - Potential_prime_function(chi_old.Data[i]);
                    }
                }
                else {return 0;}
                InitScalar(mg,   chiSol,   extchisol,t);
                InitScalar(mg, omegaSol, extomegasol,t);

        		//current timestep logfile
        		std::ofstream log( dirname +"/"+ filename   + "_time="+std::to_string(t) + ".txt");
        		SetupInterfaceRhsP1(mg, &rhs3, lset.Phi, lset.GetBndData(), extrhs3, t);
        		SetupInterfaceRhsP1(mg, &rhs4, lset.Phi, lset.GetBndData(), extrhs4,t);


        		for (int i=0; i<well_potential_concave.Data.size();i++)
                {
                    well_potential_concave.Data[i]=Potential_prime_concave_function(chi.Data[i]);
                }

                for (int i=0; i<well_potential_convex.Data.size();i++)
                {
                    well_potential_convex.Data[i]=Potential_prime_convex_function(chi.Data[i]);
                }
        		//well_potential.Data = chi_new.Data.apply(&Potential_prime_function);
        		//well_potential.Data += (-1)*chi_old.Data.apply(&Potential_prime_function);

                SetupCahnHilliardIF_P1P1(mg,  &Mass, &Normal_stab, &Tangent_stab, &Volume_stab, &Laplace, &LaplaceM,
                                         &Gprimeprime, lset.Phi, lset.GetBndData(), v, vbnd,chi_new, chibnd);
        		//current timestep logfile
                std::ofstream log_slice( dirname +"/"+ filename   + "_time="+std::to_string(t) + ".txt");

        		//set actual external force to instant rhs
        		instantrhs3 = rhs3;
        		instantrhs4 = rhs4;



        		// pick inertial term and reinitialise unknowns
        		if ((P.get<std::string>("SurfCahnHilliard.instationary") == "BDF1")||(i==0))
        		{
        			std::cout << "inertial term: BDF1 " << std::endl;
        		    instantrhs3.Data += ( 1.0/tau ) * ( Mass.Data * chi.Data) ;
                    chi_old=chi;
                }
        		else if ( P.get<std::string>("SurfCahnHilliard.instationary") == "BDF2")
        		{
        			std::cout << "inertial term: BDF2 " << std::endl;
        			instantrhs3.Data += ( 2.0/tau ) * ( Mass.Data * chi.Data) - ( 1.0/(2.0*tau)) * ( Mass.Data * chi_old.Data) ;
        			chi_old=chi;
        		}
        		if (model == "CahnHilliard") {
                    //add well-potential
                    //instantrhs3.Data *= (tau/c);
                    instantrhs4.Data -= (1.)*(Mass.Data*well_potential.Data);
                    //instantrhs4.Data += (eps*eps/2.)*(Laplace.Data*chi_new.Data);
                    instantrhs4.Data += (S) * (Mass.Data *chi_new.Data);
                    /* instantrhs4.Data += (Mass.Data*well_potential_convex.Data);
                    instantrhs4.Data += (Mass.Data*well_potential_concave.Data);*/
                    //instantrhs4.Data -= (1.) * (Gprimeprime.Data * chi_new.Data);

                    A.LinComb(sigm, LaplaceM.Data, 0, Mass.Data, alpha, Volume_stab.Data);
                    B.LinComb(0.,   Laplace.Data, (c/tau), Mass.Data);
                    C.LinComb(0.,   Laplace.Data, -1.0*eps, Mass.Data);
                    D.LinComb(eps*eps,  Laplace.Data,
                           //   -1., Gprimeprime.Data,
                              S, Mass.Data,
                              alpha*eps*eps, Volume_stab.Data);
                    //chi.Clear(0);
                    //omega.Clear(0);

                    CH_solver->Solve(A, B, C, D, omega.Data, chi.Data, instantrhs3.Data, instantrhs4.Data,
                                     omega.RowIdx->GetEx(), chi.RowIdx->GetEx());
                    //omega.Data *= 1./eps;
                    //chi.Data *= -1.;
                }
                else if (model=="AllenCahn") {
//                  instantrhs3.Data += (-sigm / eps) * (Mass.Data * well_potential_convex.Data);
//                  instantrhs3.Data += (-sigm / eps) * (Mass.Data * well_potential_concave.Data);
                    //instantrhs3.Data *= eps;
                    instantrhs3.Data -= (1./eps)*(Mass.Data*well_potential.Data);
                    instantrhs3.Data += (S/eps) * (Mass.Data *chi_new.Data);
                    //instantrhs3.Data += 1*(sigm / eps) * (Gprimeprime.Data * chi_new.Data);
                    A.LinComb(sigm * eps, Laplace.Data,
                              //1*sigm / eps, Gprimeprime.Data,
                              (S/eps + c/tau), Mass.Data, alpha*eps,
                              Volume_stab.Data);
                    AC_solver->Solve(A, chi.Data, instantrhs3.Data, chi.RowIdx->GetEx());

                }


        		//skip vtk output of needed
        		if ((i+1) % P.get<int>("Output.every timestep")  == 0)
        		{
        			std::cout << "output # : " << i << std::endl;
        			vtkwriter->Write(t);
        		}

                log_plot <<  std::to_string((float)t) << "\t";
                if ( !prFE.compare("P1")) {
                    omegaxtent.SetIdx(&P1FEidx);
                    Extend(mg, omega, omegaxtent);
                    chixtent.SetIdx(&P1FEidx);
                    Extend(mg, chi, chixtent);
                    double L2_omega_error = L2_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, omegabnd, omegaxtent),
                                                  extomegasol, t);
                    double L2_chi_error = L2_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, chibnd, chixtent),
                                               extchisol, t);
                    double L2_omega = L2_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, omegabnd, omegaxtent),
                                               ZeroScalarFun, t);
                    double L2_chi = L2_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, chibnd, chixtent),
                                             ZeroScalarFun, t);
                    for (int i=0; i<energy.Data.size();i++)
                    {
                        energy.Data[i]=Potential_function(chi.Data[i]);
                    }
                    double Lyapunov_energy = (eps/2.)*dot(Laplace.Data*chi.Data,chi.Data) + (1./eps)*dot(Mass.Data*energy.Data,id);
                    log_plot <<  std::to_string((float)L2_omega) << "\t"<<  std::to_string((float)L2_chi) <<"\t"<<
                                 std::to_string((float)L2_omega_error) << "\t"<<  std::to_string((float)L2_chi_error) << "\t"<<
                                                                                                                      std::to_string((float)Lyapunov_energy) <<
                                 std::endl;
                    Energy_L2_omega += L2_omega_error*L2_omega_error*tau;
                    Energy_L2_chi += L2_chi_error*L2_chi_error*tau;
                    if (L2_chi_error > C_norm_chi)
                    {
                        C_norm_chi = L2_chi_error;
                    }
                    if (L2_omega_error > C_norm_omega)
                    {
                        C_norm_omega = L2_omega_error;
                    }
                }
                std::cout << "Total iterations of inner 3-solver, t=" <<  t << ": " << PCGSolver3.GetIter() << '\n';
                std::cout << "Total iterations of inner 4-solver, t=" <<  t << ": " << PCGSolver4.GetIter() << '\n';
                std::cout << "Total iterations of outer solver, t=" <<  t << ": " << solver->GetIter() << '\n';

                log_slice << "Total iterations of inner 3-solver, t=" <<  t << ": " << PCGSolver3.GetIter() << '\n';
                log_slice << "Total iterations of inner 4-solver, t=" <<  t << ": " << PCGSolver4.GetIter() << '\n';
                log_slice << "Total iterations of outer solver, t=" << t << ": " << solver->GetIter() << '\n';

/*
                prev_iter=CH_solver->GetIter();
*/
                std::cout	<< "Final residual, t=" <<  t << ": " << solver->GetResid() << '\n';
                log_slice	<< "Final residual, t=" <<  t << ": " << solver->GetResid() << '\n';

            }

            std::cout << "L_2 energy error norm of 3- and 4-variable:\t" << std::sqrt(Energy_L2_omega) <<"\t" << std::sqrt(Energy_L2_chi) << std::endl;
            log_global << "L_2 energy error norm of 3- and 4-variable:\t" << std::sqrt(Energy_L2_omega) <<"\t" << std::sqrt(Energy_L2_chi) << std::endl;
            log_global << "C error norm of 3- and 4-variable:\t" << C_norm_omega <<"\t" << C_norm_chi << "\t" <<  std::endl;

            if ( P.get<std::string>("Time.Write") != "none" )
            {
                std::ofstream dump(dirname  + "/chi_dump"+"_t=" + std::to_string((float)(T0+tau*P.get<double>("Time.NumSteps"))) + ".dat");
                std::ofstream dump2(dirname  + "/omega_dump"+"_t=" + std::to_string((float)(T0+tau*P.get<double>("Time.NumSteps"))) + ".dat");

                chi.Write(dump);
                omega.Write(dump2);
            }


        }
        else if ( P.get<std::string>("SurfCahnHilliard.instationary") == "none" )
        {
            std::ofstream log_plot( dirname +"/"
                                    + "Plot_" + filename
                                    + "l="  + P.get<std::string>("Mesh.AdaptRef.FinestLevel")
                                    + ".txt");
            log_global.open( dirname +"/"
                                      + "Global_" + filename
                                      + "l="  + P.get<std::string>("Mesh.AdaptRef.FinestLevel")
                                      + ".txt");
            SetupInterfaceRhsP1(mg, &rhs3, lset.Phi, lset.GetBndData(), extrhs3, 0);
            SetupInterfaceRhsP1(mg, &rhs4, lset.Phi, lset.GetBndData(), extrhs4, 0);


            SetupCahnHilliardIF_P1P1(mg,  &Mass, &Normal_stab, &Tangent_stab, &Volume_stab, &Laplace, &LaplaceM,&Gprimeprime, lset.Phi, lset.GetBndData(), v, vbnd,chi, chibnd);


            A.LinComb(sigm,    Laplace.Data, 0 , Mass.Data, alpha, Volume_stab.Data);
            B.LinComb(0.   ,    Laplace.Data, 0 , Mass.Data);
            C.LinComb(0.   ,    Laplace.Data, 1., Mass.Data);
            D.LinComb(eps  ,    Laplace.Data, 0 , Mass.Data, alpha, Volume_stab.Data);

//            CH_solver->Solve(A, B, C, D, omega.Data, chi.Data, rhs3.Data, rhs4.Data,omega.RowIdx->GetEx(), chi.RowIdx->GetEx() );

            chi.Data *= -1.0;

            InitScalar(mg,   chiSol,   extchisol,1);
            InitScalar(mg, omegaSol, extomegasol,1);

            //skip vtk output of needed
            if ( P.get<int>("Output.every timestep")  > 0)
            {
                vtkwriter->Write((1));
            }

            if ( !prFE.compare("P1")) {
                omegaxtent.SetIdx(&P1FEidx);
                Extend(mg, omega, omegaxtent);
                chixtent.SetIdx(&P1FEidx);
                Extend(mg, chi, chixtent);
                double L2_omega = L2_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, omegabnd, omegaxtent),
                                                 &ZeroScalarFun, 0);
                double L2_chi = L2_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, chibnd, chixtent),
                                               &ZeroScalarFun, 0);
                log_global <<  "L2 3-variable: " << std::to_string((float)L2_omega) << std::endl;
                log_global <<  "L2 4-variable: " << std::to_string((float)L2_chi) << std::endl;
                double L2_omega_error = L2_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, omegabnd, omegaxtent),
                                           extomegasol, 0);
                double L2_chi_error = L2_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, chibnd, chixtent),
                                         extchisol, 0);
                log_global <<  "L2 3-variable error: " << std::to_string((float)L2_omega_error) << std::endl;
                log_global <<  "L2 4-variable error: " << std::to_string((float)L2_chi_error) << std::endl;

            }


            std::cout << "Total iterations of inner 3-solver"  << ": " << PCGSolver3.GetIter() << '\n';
            std::cout << "Total iterations of inner 4-solver" << ": " << PCGSolver4.GetIter() << '\n';
            std::cout << "Total iterations of outer solver" << ": " << solver->GetIter() << '\n';

            log_global << "Total iterations of inner 3-solver" <<   ": " << PCGSolver3.GetIter() << '\n';
            log_global << "Total iterations of inner 4-solver" <<   ": " << PCGSolver4.GetIter() << '\n';
            log_global << "Total iterations of outer solver" <<    ": " <<solver->GetIter() << '\n';

            std::cout	<< "Final residual" << ": " << solver->GetResid() << '\n';

        }

        //
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
        log_global << "The average iterationsnumber of the 3-preconditioner is: " <<     average3     << '\n' << " ...with a variation of: " << variation3 << std::endl;
        log_global << "The average iterationsnumber of the 4-preconditioner is: " << average4 << '\n' << " ...with a variation of: " << variation4 << std::endl;
        log_global << "The average iterationsnumber of the global solver is: " << averageglobal << '\n' << " ...with a standart deviation of: " << variationglobal  << std::endl;

        std::cout << "Output is located: " << dirname << std::endl;
    }

/////////////////////////////////////// Error ///////////////////////////////////////

    delete &lset;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }

}