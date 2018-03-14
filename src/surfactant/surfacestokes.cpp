/// \file surfacestokes.cpp
/// \brief Solve Vector Laplace on manifold
/// \author LNM RWTH Aachen: Thomas Jankuhn

/*
 * This file is part of DROPS.
 *
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
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
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

#include "num/oseensolver.h"
#include <fstream>

#include "surfactant/surfacestokes_funcs.h"
#include "surfactant/surfacestokes_utils.h"
#include "surfactant/surfacestokes_tests.h"

using namespace DROPS;


/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Model Test Cases ////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void InitVecLaplace(const MultiGridCL& MG, LevelsetP2CL& lset, DROPS::VecDescCL& rhs, DROPS::VecDescCL& vSol, DROPS::VecDescCL& pSol,
                   instat_vector_fun_ptr f_rhs, instat_vector_fun_ptr f_vsol, instat_scalar_fun_ptr f_psol)
{
    if( vSol.RowIdx->NumUnknownsEdge()) {
        DROPS::SetupInterfaceVectorRhsP2(MG, &rhs, lset.Phi, lset.GetBndData(), f_rhs);
    } else {
        DROPS::SetupInterfaceVectorRhsP1(MG, &rhs, lset.Phi, lset.GetBndData(), f_rhs);
    }
    InitScalar(MG, pSol, f_psol);
    InitVector(MG, vSol, f_vsol);
}


int main (int argc, char* argv[])
{
  try {
    DROPS::read_parameter_file_from_cmdline( P, argc, argv, "../../param/surfactant/surfactant/default.json");
    std::cout << P << std::endl;

    DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"), P.get<std::vector<std::string> >("General.DynamicLibs") );

    // build initial mesh
    std::cout << "Setting up interface-PDE:\n";
//    DROPS::BrickBuilderCL brick( DROPS::MakePoint3D( -2., -2., -53./48),
//                                 DROPS::MakePoint3D(4., 0., -1.),
//                                 4.*DROPS::std_basis<3>( 2),
//                                 4.*DROPS::std_basis<3>( 3),
//                                 P.get<int>("InitialDivisions"), P.get<int>("InitialDivisions"), P.get<int>("InitialDivisions"));
//    DROPS::MultiGridCL mg( brick);
    std::auto_ptr<DROPS::MGBuilderCL> builder( DROPS::make_MGBuilder( P));
    DROPS::MultiGridCL mg( *builder);
    const DROPS::ParamCL::ptree_type* ch= 0;
    try {
        ch= &P.get_child( "Domain.Periodicity");
    }
    catch (DROPS::DROPSParamErrCL) {}
    if (ch)
        read_PeriodicBoundaries( mg, *ch);

    DROPS::AdapTriangCL adap( mg );

    // choose level set
    instat_scalar_fun_ptr levelset_fun;
    std::string levelset_fun_str = P.get<std::string>("Levelset.case");
    if( !levelset_fun_str.compare("sphere_2")) {
        levelset_fun = &sphere_2;
        std::cout << "The levelset is the unit sphere." << std::endl;
    } else if( !levelset_fun_str.compare("xy_plane")) {
        levelset_fun = &xy_plane;
        std::cout << "The levelset is the xy-plane." << std::endl;
    } else if( !levelset_fun_str.compare("tilted_plane")) {
        levelset_fun = &tilted_plane;
        std::cout << "The levelset is the tilted plane." << std::endl;
    } else if( !levelset_fun_str.compare("tilted_plane_xy")) {
        levelset_fun = &tilted_plane_xy;
        std::cout << "The levelset is the tilted xy-plane." << std::endl;
    }


    // adaptive mesh refinement based on level set function
    typedef DROPS::DistMarkingStrategyCL InitMarkerT;
    InitMarkerT initmarker( levelset_fun, P.get<double>("Mesh.AdaptRef.Width"), 0, P.get<int>("Mesh.AdaptRef.FinestLevel") );
    adap.set_marking_strategy( &initmarker );
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
    std::string FE = P.get<std::string>("SurfStokes.FE");
    std::string velFE, prFE, LgFE;
    velFE = FE.substr(0,2);
    prFE = FE.substr(2,2);
    if( FE.size()>4) {
        LgFE = FE.substr(4,2);
    } else {
        LgFE = FE.substr(2,2);
    }
    bool fullgrad = P.get<bool>("SurfStokes.fullgrad");
    std::string model = P.get<std::string>("SurfStokes.model");
    std::string testcase = P.get<std::string>("SurfStokes.testcase");
    double h = P.get<DROPS::Point3DCL>("Mesh.E1")[0]/P.get<double>("Mesh.N1")*std::pow(2., -P.get<double>("Mesh.AdaptRef.FinestLevel"));
    double eta = 0.; //Constant for penalty
    double epsilon = 1.0*std::pow(h,1.0); //Constant for A_stab
    double hat_epsilon = epsilon; //Constant for L_stab
    double rho = hat_epsilon; //Constant for Schur complement preconditioner
    std::cout << "h is: " << h << std::endl;


    // construct FE spaces
    BndDataCL<Point3DCL> vbnd( 0);
    read_BndData( vbnd, mg, P.get_child( "Stokes.VelocityBndData"));

//    BndDataCL<double> vscalarbnd( 0);
//    read_BndData( vscalarbnd, mg, P.get_child( "Stokes.VelocityBndData"));

    BndDataCL<double> pbnd( 0);
    read_BndData( pbnd, mg, P.get_child( "Stokes.PressureBndData"));

    DROPS::IdxDescCL ifaceVecP2idx( vecP2IF_FE, vbnd);
    DROPS::IdxDescCL ifaceVecP1idx( vecP1IF_FE, vbnd);
    DROPS::IdxDescCL ifaceP1idx( P1IF_FE, pbnd);
    DROPS::IdxDescCL ifaceP2idx( P2IF_FE, pbnd);
    DROPS::IdxDescCL vecP2idx( vecP2_FE, vbnd);
    DROPS::IdxDescCL vecP1idx( vecP1_FE, vbnd);
    DROPS::IdxDescCL P1FEidx( P1_FE, pbnd);
    DROPS::IdxDescCL P2FEidx( P2_FE, pbnd);

    ifaceVecP2idx.GetXidx().SetBound( P.get<double>("SurfTransp.XFEMReduced"));
    ifaceVecP2idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
    ifaceVecP1idx.GetXidx().SetBound( P.get<double>("SurfTransp.XFEMReduced"));
    ifaceVecP1idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
    ifaceP1idx.GetXidx().SetBound( P.get<double>("SurfTransp.XFEMReduced"));
    ifaceP1idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
    ifaceP2idx.GetXidx().SetBound( P.get<double>("SurfTransp.XFEMReduced"));
    ifaceP2idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
    vecP2idx.GetXidx().SetBound( P.get<double>("SurfTransp.XFEMReduced"));
    vecP2idx.CreateNumbering(mg.GetLastLevel(), mg);
    vecP1idx.GetXidx().SetBound( P.get<double>("SurfTransp.XFEMReduced"));
    vecP1idx.CreateNumbering(mg.GetLastLevel(), mg);
    P1FEidx.GetXidx().SetBound( P.get<double>("SurfTransp.XFEMReduced"));
    P1FEidx.CreateNumbering(mg.GetLastLevel(), mg);
    P2FEidx.GetXidx().SetBound( P.get<double>("SurfTransp.XFEMReduced"));
    P2FEidx.CreateNumbering(mg.GetLastLevel(), mg);
    std::cout << "NumUnknowns Vector IFP2: " << ifaceVecP2idx.NumUnknowns() << std::endl;
    std::cout << "NumUnknowns Vector IFP1: " << ifaceVecP1idx.NumUnknowns() << std::endl;
    std::cout << "NumUnknowns Scalar IFP2: " << ifaceP2idx.NumUnknowns() << std::endl;
    std::cout << "NumUnknowns Scalar IFP1: " << ifaceP1idx.NumUnknowns() << std::endl;
    std::cout << "NumUnknowns Vector P2: " << vecP2idx.NumUnknowns() << std::endl;
    std::cout << "NumUnknowns Vector P1: " << vecP1idx.NumUnknowns() << std::endl;
    std::cout << "NumUnknowns Scalar P2: " << P2FEidx.NumUnknowns() << std::endl;
    std::cout << "NumUnknowns Scalar P1: " << P1FEidx.NumUnknowns() << std::endl;

    //TestAllP2MatricesWithP1Matrices(mg, lset, ifaceVecP2idx, ifaceVecP1idx, ifaceP1idx, fullgrad);

    // setup matrices
    DROPS::MatDescCL A, A_stab, B, M, S, L, L_stab, Schur, Schur_stab;
    MatrixCL Schur_hat;

    if( !FE.compare("P2P1")) {
        A.SetIdx( &ifaceVecP2idx, &ifaceVecP2idx);
        A_stab.SetIdx( &ifaceVecP2idx, &ifaceVecP2idx);
        B.SetIdx( &ifaceP1idx, &ifaceVecP2idx);
        M.SetIdx( &ifaceVecP2idx, &ifaceVecP2idx);
        S.SetIdx( &ifaceVecP2idx, &ifaceVecP2idx);
        L.SetIdx( &ifaceP1idx, &ifaceVecP2idx);
        L_stab.SetIdx( &ifaceP1idx, &ifaceVecP2idx);
        Schur.SetIdx(&ifaceP1idx, &ifaceP1idx);
        Schur_stab.SetIdx(&ifaceP1idx, &ifaceP1idx);

        SetupStokesIF_P2P1(mg, &A, &A_stab, &B, &M, &S, &L, &L_stab, &Schur, &Schur_stab, lset.Phi, lset.GetBndData(), fullgrad);

        Schur_hat.LinComb(1., Schur.Data, rho, Schur_stab.Data);
    } else if( !FE.compare("P1P1")) {
        A.SetIdx( &ifaceVecP1idx, &ifaceVecP1idx);
        A_stab.SetIdx( &ifaceVecP1idx, &ifaceVecP1idx);
        B.SetIdx( &ifaceP1idx, &ifaceVecP1idx);
        M.SetIdx( &ifaceVecP1idx, &ifaceVecP1idx);
        S.SetIdx( &ifaceVecP1idx, &ifaceVecP1idx);
        L.SetIdx( &ifaceP1idx, &ifaceVecP1idx);
        L_stab.SetIdx( &ifaceP1idx, &ifaceVecP1idx);
        Schur.SetIdx(&ifaceP1idx, &ifaceP1idx);
        Schur_stab.SetIdx(&ifaceP1idx, &ifaceP1idx);

        SetupStokesIF_P1P1(mg, &A, &A_stab, &B, &M, &S, &L, &L_stab, &Schur, &Schur_stab, lset.Phi, lset.GetBndData(), fullgrad);

        Schur_hat.LinComb(1., Schur.Data, rho, Schur_stab.Data);
    } else if( !FE.compare("P2P2")) {
        A.SetIdx( &ifaceVecP2idx, &ifaceVecP2idx);
        A_stab.SetIdx( &ifaceVecP2idx, &ifaceVecP2idx);
        B.SetIdx( &ifaceP2idx, &ifaceVecP2idx);
        M.SetIdx( &ifaceVecP2idx, &ifaceVecP2idx);
        S.SetIdx( &ifaceVecP2idx, &ifaceVecP2idx);
        L.SetIdx( &ifaceP2idx, &ifaceVecP2idx);
        L_stab.SetIdx( &ifaceP2idx, &ifaceVecP2idx);
        Schur.SetIdx(&ifaceP2idx, &ifaceP2idx);
        Schur_stab.SetIdx(&ifaceP2idx, &ifaceP2idx);

        SetupStokesIF_P2P2(mg, &A, &A_stab, &B, &M, &S, &L, &L_stab, &Schur, &Schur_stab, lset.Phi, lset.GetBndData(), fullgrad);

        Schur_hat.LinComb(1., Schur.Data, rho, Schur_stab.Data);
    } else if( !FE.compare("P1P2")) {
        A.SetIdx( &ifaceVecP1idx, &ifaceVecP1idx);
        A_stab.SetIdx( &ifaceVecP1idx, &ifaceVecP1idx);
        B.SetIdx( &ifaceP2idx, &ifaceVecP1idx);
        M.SetIdx( &ifaceVecP1idx, &ifaceVecP1idx);
        S.SetIdx( &ifaceVecP1idx, &ifaceVecP1idx);
        L.SetIdx( &ifaceP2idx, &ifaceVecP1idx);
        L_stab.SetIdx( &ifaceP2idx, &ifaceVecP1idx);
        Schur.SetIdx(&ifaceP2idx, &ifaceP2idx);
        Schur_stab.SetIdx(&ifaceP2idx, &ifaceP2idx);

        SetupStokesIF_P1P2(mg, &A, &A_stab, &B, &M, &S, &L, &L_stab, &Schur, &Schur_stab, lset.Phi, lset.GetBndData(), fullgrad);

        Schur_hat.LinComb(1., Schur.Data, rho, Schur_stab.Data);
    }

    //TestP2Matrices(mg, lset, ifaceVecP2idx, vecP2idx, vbnd, A.Data, M.Data);

    // construct iterative solver and preconditioners
    typedef SSORPcCL      SymmPcPcT;
    SymmPcPcT symmPcPc_;
    

    std::stringstream Astream;
    typedef PCGSolverCL<SymmPcPcT> PCGSolverT;
    PCGSolverT PCGSolver_( symmPcPc_, P.get<int>("Solver.PcAIter"), P.get<double>("Solver.PcATol"), true, &Astream);
    typedef SolverAsPreCL<PCGSolverT> PCGPcT;
    PCGPcT PCGPc_( PCGSolver_);


    std::stringstream Schurstream;
    PCGSolverT SchurPCGSolver (symmPcPc_, P.get<int>("Solver.PcBIter"), P.get<double>("Solver.PcBTol"), true, &Schurstream);
    SchurPreBaseCL *spc_ = new SurfaceLaplacePreCL<PCGSolverT>( Schur_hat, SchurPCGSolver);
    
    //ExpensivePreBaseCL *apc_;
//    SchurPreBaseCL  *spc_ = new DummyPreCL(1,1);  // no preconditioning for Schur
    typedef BlockPreCL<ExpensivePreBaseCL, SchurPreBaseCL, DiagSpdBlockPreCL>  DiagBlockPcT;
    typedef PLanczosONBCL<VectorCL, DiagBlockPcT> LanczosT;
    typedef PMResSolverCL<LanczosT> MinResT;
    MinResT *MinRes_;
    DiagBlockPcT    *DBlock_;
    LanczosT *lanczos_;

    DBlock_= new DiagBlockPcT( PCGPc_, *spc_);
    lanczos_= new LanczosT( *DBlock_);
    MinRes_= new MinResT( *lanczos_,  P.get<int>("Solver.Iter"), P.get<double>("Solver.Tol"), /*relative*/ false);
    BlockMatrixSolverCL<MinResT> *stokessolver= new BlockMatrixSolverCL<MinResT>( *MinRes_);


    // construct FE vectors (initialized with zero)
    DROPS::VecDescCL v, vSol, p, pSol, rhs, ZeroVec;

    if( !FE.compare("P2P1")) {
        v.SetIdx( &ifaceVecP2idx);
        vSol.SetIdx( &ifaceVecP2idx);
        p.SetIdx( &ifaceP1idx);
        pSol.SetIdx( &ifaceP1idx);
        rhs.SetIdx( &ifaceVecP2idx);
        ZeroVec.SetIdx( &ifaceVecP2idx);
    } else if( !FE.compare("P1P1")) {
        v.SetIdx( &ifaceVecP1idx);
        vSol.SetIdx( &ifaceVecP1idx);
        p.SetIdx( &ifaceP1idx);
        pSol.SetIdx( &ifaceP1idx);
        rhs.SetIdx( &ifaceVecP1idx);
        ZeroVec.SetIdx( &ifaceVecP1idx);
    } else if( !FE.compare("P2P2")) {
        v.SetIdx( &ifaceVecP2idx);
        vSol.SetIdx( &ifaceVecP2idx);
        p.SetIdx( &ifaceP2idx);
        pSol.SetIdx( &ifaceP2idx);
        rhs.SetIdx( &ifaceVecP2idx);
        ZeroVec.SetIdx( &ifaceVecP2idx);
    } else if( !FE.compare("P1P2")) {
        v.SetIdx( &ifaceVecP1idx);
        vSol.SetIdx( &ifaceVecP1idx);
        p.SetIdx( &ifaceP2idx);
        pSol.SetIdx( &ifaceP2idx);
        rhs.SetIdx( &ifaceVecP1idx);
        ZeroVec.SetIdx( &ifaceVecP1idx);
    }

    // set function pointers and rhs vectors for different test cases
    DROPS::instat_vector_fun_ptr extvsol, extsol_grad1, extsol_grad2, extsol_grad3;
    DROPS::instat_scalar_fun_ptr extpsol = &ZeroScalarFun;

    if( !levelset_fun_str.compare("sphere_2")) {
        if( !testcase.compare("1")) {
            std::cout << "Test case 1 with vSol = P[-z^2, y, x]^T" << std::endl;
            extvsol = &Test_A_plus_M_vSolVectorFun1;
            extsol_grad1 = &Test_A_plus_M_vSolVectorFun1_Gradient1;
            extsol_grad2 = &Test_A_plus_M_vSolVectorFun1_Gradient2;
            extsol_grad3 = &Test_A_plus_M_vSolVectorFun1_Gradient3;
            extpsol = &Test_A_plus_M_pSolScalarFun1;
            InitVecLaplace(mg, lset, rhs, vSol, pSol, Test_A_plus_M_RhsVectorFun1, extvsol, extpsol);
        } else if( !testcase.compare("2")) {
            std::cout << "Test case 2 with vSol = [-y, x, 0]^T" << std::endl;
            extvsol = &Test_A_plus_M_vSolVectorFun2;
            extsol_grad1 = &Test_A_plus_M_vSolVectorFun2_Gradient1;
            extsol_grad2 = &Test_A_plus_M_vSolVectorFun2_Gradient2;
            extsol_grad3 = &Test_A_plus_M_vSolVectorFun2_Gradient3;
            extpsol = &Test_A_plus_M_pSolScalarFun2;
            InitVecLaplace(mg, lset, rhs, vSol, pSol, Test_A_plus_M_RhsVectorFun2, extvsol, extpsol);
        } else if( !testcase.compare("3")) {
            std::cout << "Test case 3 with vSol = [-x*y-x*z+y^2+z^2, x^2-x*y-y*z+z^2, x^2-x*z+y^2-y*z]^T" << std::endl;
            extvsol = &Test_A_plus_M_vSolVectorFun3;
            extsol_grad1 = &Test_A_plus_M_vSolVectorFun3_Gradient1;
            extsol_grad2 = &Test_A_plus_M_vSolVectorFun3_Gradient2;
            extsol_grad3 = &Test_A_plus_M_vSolVectorFun3_Gradient3;
            extpsol = &Test_A_plus_M_pSolScalarFun3;
            InitVecLaplace(mg, lset, rhs, vSol, pSol, Test_A_plus_M_RhsVectorFun3, extvsol, extpsol);
        } else if( !testcase.compare("4")) {
            std::cout << "Test case 4 (only mass matrix) with vSol = [-x*y-x*z+y^2+z^2, x^2-x*y-y*z+z^2, x^2-x*z+y^2-y*z]^T" << std::endl;
            extvsol = &Test_A_plus_M_vSolVectorFun3;
            extsol_grad1 = &Test_A_plus_M_vSolVectorFun3_Gradient1;
            extsol_grad2 = &Test_A_plus_M_vSolVectorFun3_Gradient2;
            extsol_grad3 = &Test_A_plus_M_vSolVectorFun3_Gradient3;
            InitVecLaplace(mg, lset, rhs, vSol, pSol, Test_A_plus_M_vSolVectorFun3, extvsol, extpsol);
        } else if( !testcase.compare("5")) {
            std::cout << "Test case 5 (1_in_spherical) with vSol = P[-z^2, y, x]^T" << std::endl;
            extvsol = &Test_A_plus_M_vSolVectorFun5;
            extsol_grad1 = &Test_A_plus_M_vSolVectorFun5_Gradient1;
            extsol_grad2 = &Test_A_plus_M_vSolVectorFun5_Gradient2;
            extsol_grad3 = &Test_A_plus_M_vSolVectorFun5_Gradient3;
            extpsol = &Test_A_plus_M_pSolScalarFun5;
            InitVecLaplace(mg, lset, rhs, vSol, pSol, Test_A_plus_M_RhsVectorFun5, extvsol, extpsol);
        }
    } else if( !levelset_fun_str.compare("xy_plane")) {
        if( !testcase.compare("1")) {
            std::cout << "Test case 1 with vSol = [sin(Pi*x)sin(Pi*y), sin(Pi*x)sin(Pi*y), 0]^T" << std::endl;
            extvsol = &Test_A_plus_M_xy_plane_vSolVectorFun1;
            extsol_grad1 = &Test_A_plus_M_xy_plane_vSolVectorFun1_Gradient1;
            extsol_grad2 = &Test_A_plus_M_xy_plane_vSolVectorFun1_Gradient2;
            extsol_grad3 = &Test_A_plus_M_xy_plane_vSolVectorFun1_Gradient3;
            InitVecLaplace(mg, lset, rhs, vSol, pSol, Test_A_plus_M_xy_plane_RhsVectorFun1, extvsol, extpsol);
        } else if( !testcase.compare("2")) {
            std::cout << "Test case 2 with vSol = [sin(Pi/2*x)sin(Pi/2*y), sin(Pi/2*x)sin(Pi/2*y), 0]^T" << std::endl;
            extvsol = &Test_A_plus_M_xy_plane_vSolVectorFun2;
            extsol_grad1 = &Test_A_plus_M_xy_plane_vSolVectorFun2_Gradient1;
            extsol_grad2 = &Test_A_plus_M_xy_plane_vSolVectorFun2_Gradient2;
            extsol_grad3 = &Test_A_plus_M_xy_plane_vSolVectorFun2_Gradient3;
            InitVecLaplace(mg, lset, rhs, vSol, pSol, Test_A_plus_M_xy_plane_RhsVectorFun2, extvsol, extpsol);
        } else if( !testcase.compare("3")) {
            std::cout << "Test case 3 with vSol = [sin(Pi/4*x)sin(Pi/4*y), sin(Pi/4*x)sin(Pi/4*y), 0]^T" << std::endl;
            extvsol = &Test_A_plus_M_xy_plane_vSolVectorFun3;
            extsol_grad1 = &Test_A_plus_M_xy_plane_vSolVectorFun3_Gradient1;
            extsol_grad2 = &Test_A_plus_M_xy_plane_vSolVectorFun3_Gradient2;
            extsol_grad3 = &Test_A_plus_M_xy_plane_vSolVectorFun3_Gradient3;
            InitVecLaplace(mg, lset, rhs, vSol, pSol, Test_A_plus_M_xy_plane_RhsVectorFun3, extvsol, extpsol);
        } else if( !testcase.compare("Zero")) {
            extvsol = &ZeroVectorFun;
            extsol_grad1 = &ZeroVectorFun;
            extsol_grad2 = &ZeroVectorFun;
            extsol_grad3 = &ZeroVectorFun;
        } else if( !testcase.compare("4")) {
            std::cout << "Test case 4 with vSol = [(x^2-4)*(y^2-4), (x^2-4)*(y^2-4), 0]^T" << std::endl;
            extvsol = &Test_A_plus_M_xy_plane_vSolVectorFun4;
            extsol_grad1 = &Test_A_plus_M_xy_plane_vSolVectorFun4_Gradient1;
            extsol_grad2 = &Test_A_plus_M_xy_plane_vSolVectorFun4_Gradient2;
            extsol_grad3 = &Test_A_plus_M_xy_plane_vSolVectorFun4_Gradient3;
            InitVecLaplace(mg, lset, rhs, vSol, pSol, Test_A_plus_M_xy_plane_RhsVectorFun4, extvsol, extpsol);
        }
    } else if( !levelset_fun_str.compare("tilted_plane")) {
        if( !testcase.compare("1")) {
            std::cout << "Test case 1 with vSol = [4/5*sin(Pi*x)sin(Pi*z), 2/5*sin(Pi*x)sin(Pi*z), sin(Pi*x)sin(Pi*z)]^T" << std::endl;
            extvsol = &Test_A_plus_M_tilted_plane_vSolVectorFun1;
            extsol_grad1 = &Test_A_plus_M_tilted_plane_vSolVectorFun1_Gradient1;
            extsol_grad2 = &Test_A_plus_M_tilted_plane_vSolVectorFun1_Gradient2;
            extsol_grad3 = &Test_A_plus_M_tilted_plane_vSolVectorFun1_Gradient3;
            InitVecLaplace(mg, lset, rhs, vSol, pSol, Test_A_plus_M_tilted_plane_RhsVectorFun1, extvsol, extpsol);
        } else if( !testcase.compare("2")) {
            std::cout << "Test case 2 with vSol = [4/5*(x^2-4)*(z^2-4), 2/5*(x^2-4)*(z^2-4), (x^2-4)*(z^2-4)]^T" << std::endl;
            extvsol = &Test_A_plus_M_tilted_plane_vSolVectorFun2;
            extsol_grad1 = &Test_A_plus_M_tilted_plane_vSolVectorFun2_Gradient1;
            extsol_grad2 = &Test_A_plus_M_tilted_plane_vSolVectorFun2_Gradient2;
            extsol_grad3 = &Test_A_plus_M_tilted_plane_vSolVectorFun2_Gradient3;
            InitVecLaplace(mg, lset, rhs, vSol, pSol, Test_A_plus_M_tilted_plane_RhsVectorFun2, extvsol, extpsol);
        }
    } else if( !levelset_fun_str.compare("tilted_plane_xy")) {
        if( !testcase.compare("1")) {
            std::cout << "Test case 1 with vSol = [4/17*sin(Pi*y)*sin(Pi*z), sin(Pi*y)*sin(Pi*z), 16/17*sin(Pi*y)*sin(Pi*z)]^T" << std::endl;
            extvsol = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun1;
            extsol_grad1 = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun1_Gradient1;
            extsol_grad2 = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun1_Gradient2;
            extsol_grad3 = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun1_Gradient3;
            InitVecLaplace(mg, lset, rhs, vSol, pSol, Test_A_plus_M_tilted_plane_xy_RhsVectorFun1, extvsol, extpsol);
        } else if( !testcase.compare("Zero")) {
            extvsol = &ZeroVectorFun;
            extsol_grad1 = &ZeroVectorFun;
            extsol_grad2 = &ZeroVectorFun;
            extsol_grad3 = &ZeroVectorFun;
        } else if( !testcase.compare("2")) {
            std::cout << "Test case 2 with vSol = [4/17*(z^2-4)*(y^2-4), (z^2-4)*(y^2-4), 16/17*(z^2-4)*(y^2-4)]^T" << std::endl;
            extvsol = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun2;
            extsol_grad1 = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun2_Gradient1;
            extsol_grad2 = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun2_Gradient2;
            extsol_grad3 = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun2_Gradient3;
            InitVecLaplace(mg, lset, rhs, vSol, pSol, Test_A_plus_M_tilted_plane_xy_RhsVectorFun2, extvsol, extpsol);
        } else if( !testcase.compare("3")) {
            std::cout << "Test case 3 with vSol = [4/17*(z^2-4)^2*(y^2-4)^2, (z^2-4)^2*(y^2-4)^2, 16/17*(z^2-4)^2*(y^2-4)^2]^T" << std::endl;
            extvsol = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun3;
            extsol_grad1 = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun3_Gradient1;
            extsol_grad2 = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun3_Gradient2;
            extsol_grad3 = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun3_Gradient3;
            InitVecLaplace(mg, lset, rhs, vSol, pSol, Test_A_plus_M_tilted_plane_xy_RhsVectorFun3, extvsol, extpsol);
        }
    }

/////////////////////////////////////// A + M with Lagrange-Multiplier (and Penalty) ///////////////////////////////////////

    if( !model.compare("VectorLaplaceLagrange")) {

        MatrixCL Ahat, Bhat;
        Ahat.LinComb(1., A.Data, 1., M.Data, eta, S.Data, epsilon, A_stab.Data);
        Bhat.LinComb(1., L.Data, hat_epsilon, L_stab.Data);


        WriteToFile(Ahat, "MatrixAhat2", "Test");
        WriteToFile(Bhat, "MatrixBhat2", "Test");
        WriteToFile(Schur_hat, "MatrixSchur_hat2", "Test");

        MatrixCL BhatTranspose;
        transpose(Bhat, BhatTranspose);

        VectorCL calcrhs1 = Ahat*vSol.Data;
        calcrhs1 += BhatTranspose*pSol.Data;
        //VectorCL calcrhs2 = Bhat*vSol.Data;
        //sfsdf
        VectorCL calcrhs2 = Bhat*ZeroVec.Data;
        MatrixCL Zero(calcrhs2);

        stokessolver->Solve(Ahat, Bhat, Zero, v.Data, p.Data, rhs.Data, calcrhs2, v.RowIdx->GetEx(), p.RowIdx->GetEx() ); // Does not work in parallel!!!

        std::cout << "Iter: " << stokessolver->GetIter() << "\tres: " << stokessolver->GetResid() << '\n';

        std::cout << "DER WERT norm(calcrhs1 - rhs) ist: " << norm (calcrhs1 -rhs.Data) << std::endl;
        std::cout << "DER WERT norm(Ahat*v + BhatTranspose*p - rhs) ist: " << norm(Ahat*v.Data + BhatTranspose*p.Data - rhs.Data) << std::endl;
        std::cout << "DER WERT norm(Bhat*v - calcrhs2) ist: " << norm(Bhat*v.Data - calcrhs2) << std::endl;

        std::cout << "DER WERT norm(Ahat*vSol + BhatTranspose*pSol - rhs) ist: " << norm(Ahat*vSol.Data + BhatTranspose*pSol.Data - rhs.Data) << std::endl;
        std::cout << "DER WERT norm(Bhat*vSol - calcrhs2) ist: " << norm(Bhat*vSol.Data - calcrhs2) << std::endl;

        DROPS::VecDescCL resxtent, res;
        resxtent.SetIdx( &vecP2idx);
        res.SetIdx( &ifaceVecP2idx);
        res.Data = Ahat*vSol.Data + BhatTranspose*pSol.Data - rhs.Data;
        Extend(mg, res, resxtent);
        std::cout << "The L2-Norm of Ahat*vSol + BhatTranspose*pSol - rhs is: " << L2_Vector_error(mg, lset.Phi, lset.GetBndData(), make_P2Eval(mg, vbnd, resxtent), ZeroVectorFun) << std::endl;

        std::cout << "DER WERT norm(A*v + M*v - rhs) ist: " << norm(A.Data*v.Data + M.Data*v.Data - rhs.Data) << std::endl;
        std::cout << "DER WERT norm(A*vSol + M*vSol - rhs) ist: " << norm(A.Data*vSol.Data + M.Data*vSol.Data - rhs.Data) << std::endl;

        std::cout << "Es sollte (vSol) " << dot(vSol.Data, A.Data*vSol.Data) + dot(vSol.Data, M.Data*vSol.Data) << "=" << dot(rhs.Data, vSol.Data) << " gelten." << std::endl;
        std::cout << "Es sollte (v) " << dot(v.Data, A.Data*v.Data) + dot(v.Data, M.Data*v.Data) << "=" << dot(rhs.Data, v.Data) << " gelten." << std::endl;

        std::cout << "Es sollte (v) res " << std::abs(dot(v.Data, A.Data*v.Data) + dot(v.Data, M.Data*v.Data) - dot(rhs.Data, v.Data)) << " gelten." << std::endl;
    }

//////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////// A + M with Penalty ///////////////////////////////////////

    if( !model.compare("VectorLaplacePenalty")) {
        PCGSolverT *solver = new PCGSolverT (symmPcPc_, P.get<int>("Solver.Iter"), P.get<double>("Solver.Tol"),false, &std::cout);

        MatrixCL Lefthandside;
        Lefthandside.LinComb(1., A.Data, 1., M.Data, eta*std::pow(h,-2), S.Data, epsilon, A_stab.Data);
        VectorCL righthandside = Lefthandside*vSol.Data;

//        WriteToFile(A.Data, "MatrixA", "Test");
//        WriteToFile(M.Data, "MatrixM", "Test");
//        WriteToFile(S.Data, "MatrixS", "Test");
//        WriteToFile(A_stab.Data, "MatrixA_stab", "Test");
//        WriteToFile(Lefthandside, "MatrixLefthandside", "Test");

        solver->Solve(Lefthandside, v.Data, rhs.Data, v.RowIdx->GetEx());

        std::cout << "Es sollte (vSol) " << dot(vSol.Data, A.Data*vSol.Data) + dot(vSol.Data, M.Data*vSol.Data) << "=" << dot(rhs.Data, vSol.Data) << " gelten." << std::endl;
        std::cout << "Es sollte (v) " << dot(v.Data, A.Data*v.Data) + dot(v.Data, M.Data*v.Data) << "=" << dot(rhs.Data, v.Data) << " gelten." << std::endl;
    }

//////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////// Stokes (not working yet) ///////////////////////////////////////


//    Test_Stokes_1(mg, lset, &rhs, &vSol, pSol);

//    MatrixCL Ahat, Bhat;

//    Ahat.LinComb(1., A.Data, 1., M.Data, 500., S.Data, 0.1, N.Data);
//    Bhat.concat_under(B.Data, L.Data);

//    MatrixCL Null(K.Data.num_rows(), K.Data.num_cols(), 0.);

//    BlockMatrixBaseCL<MatrixCL> C( &K.Data, MUL, &Null, MUL, &Null, MUL, &K.Data, MUL );
//    MatrixCL Chat = BuildMatrix(C);
//    Chat *= 0.1;

//    VectorCL calcrhs1 = Ahat*vSol.Data;


////    //WriteToFile(Ahat, "MatrixAhat", "Test");
////    //WriteToFile(K.Data, "MatrixK", "Test");
////    //WriteToFile(B.Data, "MatrixB", "Test");


//    VectorCL w(p.Data.size()+q.Data.size());
//    w[std::slice(0, p.Data.size(),1)] = p.Data;
//    w[std::slice(p.Data.size(), q.Data.size(), 1)] = q.Data;
////    VectorCL b(v.Data.size());
//    VectorCL c1(w.size());
//    c1[std::slice(0, pSol.Data.size(),1)] = K.Data*VectorCL(0.1*pSol.Data);
//    VectorCL c2(w.size());
//    c2[std::slice(0, pSol.Data.size(),1)] = pSol.Data;

//    MatrixCL BhatTranspose;
//    transpose(Bhat, BhatTranspose);
////    calcrhs1+= transp_mul(Bhat,c);

//    MatrixCL BTranspose;
//    transpose(B.Data, BTranspose);
//    calcrhs1 += transp_mul(B.Data, pSol.Data);

//    MatrixCL BTranspose;
//    transpose(B.Data, BTranspose);
//    MatrixCL LTranspose;
//    transpose(L.Data, LTranspose);
//    std::cout << "B*vSol = " << dot(k.Data, B.Data*vSol.Data) << std::endl;
//    std::cout << "L*vSol = " << dot(k.Data, L.Data*vSol.Data) << std::endl;
//    std::cout << "B^T*pSol = " << dot(vSol3.Data, BTranspose*pSol.Data) << std::endl;
//    std::cout << "L^T*q = " << dot(vSol3.Data, LTranspose*q.Data) << std::endl;
//    std::cout << "Ahat*vSol + B^T*pSol + L^T*q = " << dot(vSol3.Data, Ahat*vSol.Data) + dot(vSol3.Data, BTranspose*pSol.Data) + dot(vSol3.Data, LTranspose*q.Data) << std::endl;
//    std::cout << "rhs*vSol3 = " << dot(rhs.Data, vSol3.Data) << std::endl;

//    std::cout << "A: " << dot(vSol3.Data, A.Data*vSol.Data) << " M: " << dot(vSol3.Data, M.Data*vSol.Data) << " S: " << dot(vSol3.Data, S.Data*vSol.Data) << " N: " << dot(vSol3.Data, N.Data*vSol.Data) << " B: " << dot(pSol.Data, B.Data*vSol3.Data) << " L: " << dot(q.Data, L.Data*vSol3.Data) << " rhs: " << dot(rhs.Data, vSol3.Data) << std::endl;
//    std::cout << "Ahat: " << dot(vSol.Data, Ahat*vSol.Data) << std::endl;

//    VectorCL calcrhs2 = Bhat*vSol.Data;
//    calcrhs2 += c1;

//    VectorCL calcrhs2 = B.Data*vSol.Data;
//    calcrhs2 += K.Data*VectorCL(pSol.Data);

//    MatrixCL Kscale = K.Data;
    //Kscale *= 0.1;

//    std::cout << "DER WERT norm(calcrhs1 - rhs) ist: " << norm(calcrhs1 - rhs.Data) << std::endl;

//    stokessolver->Solve(Ahat, Bhat, Chat, v.Data, w, rhs.Data, calcrhs2, v.RowIdx->GetEx(), v.RowIdx->GetEx() ); // Does not work in parallel!!!

//    stokessolver->Solve(Ahat, B.Data, Kscale, v.Data, p.Data, rhs.Data, calcrhs2, v.RowIdx->GetEx(), v.RowIdx->GetEx() ); // Does not work in parallel!!!



//    std::cout << "Iter: " << stokessolver->GetIter() << "\tres: " << stokessolver->GetResid() << '\n';

//    std::cout << "DER WERT norm(Ahat*v +BhatTranspose*w - rhs) ist: " << norm(Ahat*v.Data + BhatTranspose*w - rhs.Data) << std::endl;
//    std::cout << "DER WERT norm(Bhat*v - calcrhs2) ist: " << norm(Bhat*v.Data - calcrhs2) << std::endl;

//    std::cout << "DER WERT norm(Ahat*vSol +BhatTranspose*c - rhs) ist: " << norm(Ahat*vSol.Data + BhatTranspose*c2 - rhs.Data) << std::endl;
//    std::cout << "DER WERT norm(Bhat*vSol + c1 - calcrhs2) ist: " << norm(Bhat*vSol.Data + c1 - calcrhs2) << std::endl;

//    std::cout << "DER WERT norm(Ahat*v +BTranspose*p - rhs) ist: " << norm(Ahat*v.Data + BTranspose*p.Data - rhs.Data) << std::endl;
//    std::cout << "DER WERT norm(B*v + K*p - calcrhs2) ist: " << norm(B.Data*v.Data + K.Data*VectorCL(p.Data)  - calcrhs2) << std::endl;

//    std::cout << "DER WERT norm(Ahat*vSol +BTranspose*pSol - rhs) ist: " << norm(Ahat*vSol.Data + BTranspose*pSol.Data - rhs.Data) << std::endl;
//    std::cout << "DER WERT norm(B*vSol + K*pSol - calcrhs2) ist: " << norm(B.Data*vSol.Data + K.Data*VectorCL(pSol.Data) - calcrhs2) << std::endl;

//    p.Data = w[std::slice(0, p.Data.size(),1)];

//////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////// Error ///////////////////////////////////////

    double Aaverage( 0.), Avariation( 0.), Schuraverage( 0.), Schurvariation( 0.);

    std::stringstream Astreamcopy( Astream.str());
    std::stringstream Schurstreamcopy( Schurstream.str());

    ComputeAverageIterations(Astream, Aaverage);
    ComputeAverageIterations(Schurstream, Schuraverage);
    ComputeVariationFromAverageIterations(Astreamcopy, Aaverage, Avariation);
    ComputeVariationFromAverageIterations(Schurstreamcopy, Schuraverage, Schurvariation);

    std::cout << "The average iterationsnumber of the A-preconditioner is " << Aaverage << " with a variation of " << Avariation << std::endl;
    std::cout << "The average iterationsnumber of the Schur-preconditioner is: " << Schuraverage << " with a variation of " << Schurvariation << std::endl;

    //if( !testcase.compare("Zero")) {
//        if( !levelset_fun_str.compare("xy_plane")) {
//            InitVector(mg, &v, &Test_A_plus_M_xy_plane_vSolVectorFun2);
//        } else if( !levelset_fun_str.compare("tilted_plane_xy")) {
//            InitVector(mg, &v, &Test_A_plus_M_tilted_plane_xy_vSolVectorFun2);
//        }
    //}
//    TestRhsP2(mg, lset, &ifaceVecP2idx);
//extvsol = &ZeroVectorFun;
//InitVector(mg, v, &Test_A_plus_M_vSolVectorFun5);
//extvsol = &Test_A_plus_M_vSolVectorFun5;
//DROPS::VecDescCL vSol3(&ifaceVecP1idx);
//int s = vSol3.Data.size();
//vSol3.Data = v.Data[std::slice(0,s,1)];

std::cout << "The value of v^T*A_stab*v is: " << dot(A_stab.Data*v.Data, v.Data) << std::endl;
std::cout << "The value of pSol^T*L_stab*v is: " << dot(L_stab.Data*v.Data, pSol.Data) << std::endl;

    DROPS::VecDescCL vxtent, pxtent;

    if( !velFE.compare("P2")) {
        vxtent.SetIdx( &vecP2idx);
        Extend(mg, v, vxtent);
        BndDataCL<Point3DCL> bndvec = vbnd;
        std::cout << "The L2-Norm of v - vSol is: " << L2_Vector_error(mg, lset.Phi, lset.GetBndData(), make_P2Eval(mg, bndvec, vxtent), extvsol) << std::endl;
        double H1( 0.), surfH1( 0.), advanced_surfH1( 0.), normal_velocity;
        H1_Vector_error_P2(mg, lset.Phi, lset.GetBndData(), v, vbnd, extvsol, extsol_grad1, extsol_grad2, extsol_grad3, eta, H1, surfH1, advanced_surfH1, normal_velocity);
//        std::cout << "The H1-Norm of v - vSol is: " << H1 << std::endl;
//        std::cout << "The surfH1-Norm of v - vSol is: " << surfH1 << std::endl;
        std::cout << "The advanced surfH1-Norm of v - vSol is: " << advanced_surfH1 << std::endl;
        std::cout << "The U-Norm of v - vSol is: " << std::sqrt(advanced_surfH1*advanced_surfH1 + epsilon*dot(A_stab.Data*v.Data, v.Data)) << std::endl;
        std::cout << "The L2-Norm of v * n is: " << normal_velocity << std::endl;
        if ( !prFE.compare("P1")) {
            pxtent.SetIdx( &P1FEidx);
            Extend(mg, p, pxtent);
            BndDataCL<double> bndscalar = pbnd;
            double L2_Lagrange =L2_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, bndscalar, pxtent), extpsol);
            std::cout << "The L2-Norm of p - pSol is: " << L2_Lagrange << std::endl;
            std::cout << "The M-Norm of p - pSol is: " << std::sqrt(L2_Lagrange*L2_Lagrange + hat_epsilon*dot(Schur_stab.Data*p.Data, p.Data)) << std::endl;
        }
    } else if( !velFE.compare("P1")) {
        vxtent.SetIdx( &vecP1idx);
        Extend(mg, v, vxtent);
        BndDataCL<Point3DCL> bndvec = vbnd;
        std::cout << "The L2-Norm of v - vSol is: " << L2_Vector_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, bndvec, vxtent), extvsol) << std::endl;
        double H1( 0.), surfH1( 0.), advanced_surfH1( 0.), normal_velocity( 0.);
        H1_Vector_error_P1(mg, lset.Phi, lset.GetBndData(), v, vbnd, extvsol, extsol_grad1, extsol_grad2, extsol_grad3, eta, H1, surfH1, advanced_surfH1, normal_velocity);
//        std::cout << "The H1-Norm of v - vSol is: " << H1 << std::endl;
//        std::cout << "The surfH1-Norm of v - vSol is: " << surfH1 << std::endl;
        std::cout << "The advanced surfH1-Norm of v - vSol is: " << advanced_surfH1 << std::endl;
        std::cout << "The U-Norm of v - vSol is: " << std::sqrt(advanced_surfH1*advanced_surfH1 + epsilon*dot(A_stab.Data*v.Data, v.Data)) << std::endl;
        std::cout << "The L2-Norm of v * n is: " << normal_velocity << std::endl;
        if ( !prFE.compare("P1")) {
            pxtent.SetIdx( &P1FEidx);
            Extend(mg, p, pxtent);
            BndDataCL<double> bndscalar = pbnd;
            double L2_Lagrange =L2_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, bndscalar, pxtent), extpsol);
            std::cout << "The L2-Norm of p - pSol is: " << L2_Lagrange << std::endl;
            std::cout << "The M-Norm of p - pSol is: " << std::sqrt(L2_Lagrange*L2_Lagrange + hat_epsilon*dot(Schur_stab.Data*p.Data, p.Data)) << std::endl;
        }
    }

//////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////// VTK-output///////////////////////////////////////

    VTKOutCL * vtkwriter = NULL;

    vtkwriter = new VTKOutCL(mg, "DROPS data", 1, "vtk", "IFStokes", "none", 0);
    vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "level-set") );
    vtkwriter->Register( make_VTKIfaceVector(mg, vSol, "velSol", velFE, vbnd));
    vtkwriter->Register( make_VTKIfaceVector(mg, rhs, "rhs", velFE, vbnd));
//    vtkwriter->Register( make_VTKIfaceVector(mg, vSol3, "velSol3", "P1", vbnd));
    vtkwriter->Register( make_VTKIfaceVector(mg, v, "velocity", velFE, vbnd));
    vtkwriter->Register( make_VTKIfaceScalar(mg, p, "pressure", pbnd));
    vtkwriter->Register( make_VTKIfaceScalar(mg, pSol, "pSol", pbnd));
    vtkwriter->Write(0);

//////////////////////////////////////////////////////////////////////////////

    delete &lset;
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
