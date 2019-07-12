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
 * GNU Lesser General Public License for more dtau_uils.
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
#include <surfactant/ifacetransp.h>

#include "surfnavierstokes/surfnavierstokes_funcs.h"
#include "surfnavierstokes/surfnavierstokes_utils.h"
#include "surfnavierstokes/surfnavierstokes_tests.h"

using namespace DROPS;

//void InitVecLaplace(const MultiGridCL& MG, LevelsetP2CL& lset, DROPS::VecDescCL& rhs, DROPS::VecDescCL& vSol, DROPS::VecDescCL& pSol,
//                    instat_vector_fun_ptr f_rhs, instat_vector_fun_ptr f_vsol, instat_scalar_fun_ptr f_psol, double t = 0.0) {
//    if( vSol.RowIdx->NumUnknownsEdge()) {
//        DROPS::SetupInterfaceVectorRhsP2(MG, &rhs, lset.Phi, lset.GetBndData(), f_rhs);
//    } else {
//        DROPS::SetupInterfaceVectorRhsP1(MG, &rhs, lset.Phi, lset.GetBndData(), f_rhs, t);
//    }
//    InitScalar(MG, pSol, f_psol);
//    InitVector(MG, vSol, f_vsol);
//}

DROPS::Point3DCL shift;
DROPS::Point3DCL shiftTransform(const Point3DCL& p) {
    return p - shift;
}


int main(int argc, char* argv[]) {
    try {
        DROPS::read_parameter_file_from_cmdline( P, argc, argv, "../../param/surfnavierstokes/No_Bnd_Condition.json");
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
        std::auto_ptr<DROPS::MGBuilderCL> builder(DROPS::make_MGBuilder(P));
        DROPS::MultiGridCL mg(*builder);
        double h = P.get<DROPS::Point3DCL>("Mesh.E1")[0] / P.get<double>("Mesh.N1") * std::pow(2., -P.get<double>("Mesh.AdaptRef.FinestLevel"));
        std::cout << "h is: " << std::to_string(float(h)) << '\n';
        // levelset shift
        shift = P.get<DROPS::Point3DCL>("Levelset.ShiftDir", DROPS::Point3DCL(0., 0., 0.));
        shift /= shift.norm();
        auto shiftNorm = fabs(P.get<double>("Levelset.ShiftNorm", 0.));
        shift *= shiftNorm;
        std::cout << "surface shift is: " << shift << '\n';
        mg.Transform(shiftTransform);
        const DROPS::ParamCL::ptree_type* ch= 0;
        try {
            ch= &P.get_child("Mesh.Periodicity");
        } catch (DROPS::DROPSParamErrCL) {}
        if (ch) read_PeriodicBoundaries(mg, *ch);
        DROPS::AdapTriangCL adap(mg);
        // choose level set
        instat_scalar_fun_ptr levelset_fun;
        instat_vector_fun_ptr exact_normal;
        instat_matrix_fun_ptr exact_shape;
        std::string levelset_fun_str = P.get<std::string>("Levelset.case");
        if( !levelset_fun_str.compare("sphere_2")) {
            levelset_fun = &sphere_2;
            exact_normal = &sphere_2_normal;
            exact_shape  = &sphere_2_shape;
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
        } else if( !levelset_fun_str.compare("cube_madeof_edges")) {
            levelset_fun = &cube_madeof_edges;
            std::cout << "The levelset is the cube_madeof_edges." << std::endl;
        }
        else if( !levelset_fun_str.compare("torus")) {
            levelset_fun = &torus;
            std::cout << "The levelset is the torus." << std::endl;
        }
        else if( !levelset_fun_str.compare("torus_flower")) {
            levelset_fun = &torus_flower;
            std::cout << "The levelset is the torus_flower." << std::endl;
        }
        LocalStokesParam param;
        param.input.exactNormal = P.get<bool>("SurfNavStokes.ComputeNormalErr") ? exact_normal : nullptr;
        param.input.exactShape  = P.get<bool>("SurfNavStokes.ComputeShapeErr")  ? exact_shape  : nullptr;
        param.input.computeMatrices = P.get<bool>("SurfNavStokes.ComputeMatrices");
        // adaptive mesh refinement based on level set function
        typedef DROPS::DistMarkingStrategyCL InitMarkerT;
        InitMarkerT initmarker(levelset_fun, P.get<double>("Mesh.AdaptRef.Width"), 0, P.get<int>("Mesh.AdaptRef.FinestLevel") );
        adap.set_marking_strategy(&initmarker );
        adap.MakeInitialTriang();
        adap.set_marking_strategy(0);
        // create level set
        instat_scalar_fun_ptr sigma (0);
        SurfaceTensionCL sf( sigma, 0);
        BndDataCL<double> lsbnd( 0);
        read_BndData( lsbnd, mg, P.get_child( "Levelset.BndData"));
        DROPS::LevelsetP2CL & lset( * DROPS::LevelsetP2CL::Create( mg, lsbnd, sf) );
        param.input.numbOfVirtualSubEdges = P.get<size_t >("SurfNavStokes.NumbOfVirtualSubEdges");
        param.input.usePatchNormal = P.get<bool>("SurfNavStokes.UsePatchNormals");
        if (param.input.usePatchNormal) std::cout << "using patch normals in surf integrands\n";
        lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
        lset.Phi.SetIdx( &lset.idx);
        // LinearLSInit( mg, lset.Phi, levelset_fun);
        lset.Init(levelset_fun);
        // parse FE types and some other parameters from json file
        std::string FE = P.get<std::string>("SurfNavStokes.FE");
        std::string velFE, prFE, LgFE;
        velFE = FE.substr(0,2);
        prFE = FE.substr(2,2);
        LgFE = FE.size() > 4 ? FE.substr(4,2) : FE.substr(2,2);
        std::string model = P.get<std::string>("SurfNavStokes.model");
        std::string testcase = P.get<std::string>("SurfNavStokes.testcase");
        double tau = P.get<double>("Time.StepSize");
        ParameterNS::nu = P.get<double>("SurfNavStokes.kinematic_viscosity");
        auto formulation    = P.get<std::string>("SurfNavStokes.formulation");
        if (formulation == "consistent") {
            param.input.formulation = LocalStokesParam::Formulation::consistent;
            std::cout << "using consistent penalty formulation\n";
        }
        else {
            param.input.formulation = LocalStokesParam::Formulation::inconsistent;
            std::cout << "using inconsistent penalty formulation\n";
        }
        auto tau_u_order  = P.get<double>("SurfNavStokes.normal_penalty_pow");
        auto tau_u_factor = P.get<double>("SurfNavStokes.normal_penalty_fac");
        auto tau_u 		  = tau_u_factor * pow(h, tau_u_order); // constant for normal penalty
        auto rho_u_order  = P.get<double>("SurfNavStokes.vel_volumestab_pow");
        auto rho_u_factor = P.get<double>("SurfNavStokes.vel_volumestab_fac");
        auto rho_u        = rho_u_factor  * pow(h, rho_u_order); // constant for velocity stabilisation
        auto rho_p_order  = 1.0;
        auto rho_p        = 1.e0  * pow(h, rho_p_order); // constant for pressure stabilisation
        auto rho          = 1.e0  * pow(h, 1); // constant for Schur complement preconditioner
        double hat_rho_u  = rho_u; //Constant for L_stab
        std::cout << "tau is: " << tau << std::endl;
        ParameterNS::h = h;
        // construct FE spaces
        BndDataCL<Point3DCL> vbnd( 0);
        read_BndData( vbnd, mg, P.get_child( "Stokes.VelocityBndData"));
    //  BndDataCL<double> vscalarbnd( 0);
    //  read_BndData( vscalarbnd, mg, P.get_child( "Stokes.VelocityBndData"));
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
        ifaceVecP2idx.GetXidx().SetBound( P.get<double>("SurfTransp.OmitBound"));
        ifaceVecP2idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
        ifaceVecP1idx.GetXidx().SetBound( P.get<double>("SurfTransp.OmitBound"));
        ifaceVecP1idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
        ifaceP1idx.GetXidx().SetBound( P.get<double>("SurfTransp.OmitBound"));
        ifaceP1idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
        ifaceP2idx.GetXidx().SetBound( P.get<double>("SurfTransp.OmitBound"));
        ifaceP2idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
        vecP2idx.GetXidx().SetBound( P.get<double>("SurfTransp.OmitBound"));
        vecP2idx.CreateNumbering(mg.GetLastLevel(), mg);
        vecP1idx.GetXidx().SetBound( P.get<double>("SurfTransp.OmitBound"));
        vecP1idx.CreateNumbering(mg.GetLastLevel(), mg);
        P1FEidx.GetXidx().SetBound( P.get<double>("SurfTransp.OmitBound"));
        P1FEidx.CreateNumbering(mg.GetLastLevel(), mg);
        P2FEidx.GetXidx().SetBound( P.get<double>("SurfTransp.OmitBound"));
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
        // construct FE vectors (initialized with zero)
        DROPS::VecDescCL v, vSol, vInit, p, curl,curlSol, pSol, curl_proj, ZeroVec;
        DROPS::VecDescCL v_iter, v_temp,v_aux,v_old,v_oldold,p_aux,p_temp,p_iter,p_old;//for steady nonlinear procedure or BDF2
        StokesSystem system;
        auto& fRHS = system.fRHS;
        auto& gRHS = system.gRHS;
        if( !FE.compare("P2P1")) {
             v.SetIdx(&ifaceVecP2idx);
             vSol.SetIdx(&ifaceVecP2idx);
             vInit.SetIdx( &ifaceVecP2idx);
             p.SetIdx(&ifaceP1idx);
             pSol.SetIdx(&ifaceP1idx);
             fRHS.SetIdx(&ifaceVecP2idx);
             gRHS.SetIdx(&ifaceP1idx);
             ZeroVec.SetIdx(&ifaceVecP2idx);
             //
             curl.SetIdx( &ifaceP1idx);
             curlSol.SetIdx( &ifaceP1idx);
        } else if( !FE.compare("P1P1")) {
             v.SetIdx( &ifaceVecP1idx);
             vSol.SetIdx( &ifaceVecP1idx);
             vInit.SetIdx( &ifaceVecP1idx);
             v_old.SetIdx( &ifaceVecP1idx);
             v_oldold.SetIdx( &ifaceVecP1idx);

             v_iter.SetIdx( &ifaceVecP1idx);
             v_temp.SetIdx( &ifaceVecP1idx);
             v_aux.SetIdx( &ifaceVecP1idx);
             p.SetIdx( &ifaceP1idx);
             curl.SetIdx( &ifaceP1idx);
             curlSol.SetIdx( &ifaceP1idx);
             p_old.SetIdx( &ifaceP1idx);
             p_aux.SetIdx( &ifaceP1idx);
             p_temp.SetIdx( &ifaceP1idx);
             p_iter.SetIdx( &ifaceP1idx);
             pSol.SetIdx( &ifaceP1idx);
             gRHS.SetIdx(&ifaceP1idx);
             curl_proj.SetIdx(&ifaceP1idx);
             fRHS.SetIdx( &ifaceVecP1idx);
             ZeroVec.SetIdx( &ifaceVecP1idx);
        } else if( !FE.compare("P2P2")) {
             v.SetIdx( &ifaceVecP2idx);
             vSol.SetIdx( &ifaceVecP2idx);
             p.SetIdx( &ifaceP2idx);
             pSol.SetIdx( &ifaceP2idx);
             fRHS.SetIdx( &ifaceVecP2idx);
             ZeroVec.SetIdx( &ifaceVecP2idx);
        } else if( !FE.compare("P1P2")) {
             v.SetIdx( &ifaceVecP1idx);
             vSol.SetIdx( &ifaceVecP1idx);
             p.SetIdx( &ifaceP2idx);
             pSol.SetIdx( &ifaceP2idx);
             fRHS.SetIdx( &ifaceVecP1idx);
             ZeroVec.SetIdx( &ifaceVecP1idx);
        }
        // set function pointers and rhs vectors for different test cases
        DROPS::instat_vector_fun_ptr extvsol, extsol_grad1, extsol_grad2, extsol_grad3, extvinit, extrhs;
        DROPS::instat_scalar_fun_ptr extpsol = &ZeroScalarFun, extrhs2, extcurlsol=&ZeroScalarFun;
        if(!levelset_fun_str.compare("sphere_2")) {
              if(!testcase.compare("1")) {
                  std::cout << "Test case 1 with vSol = P[-z^2, y, x]^T" << std::endl;
                  extvsol = &Test_A_plus_M_vSolVectorFun1;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun1_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun1_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun1_Gradient3;
                  extpsol = &Test_A_plus_M_pSolScalarFun1;
                  extrhs = &Test_A_plus_M_RhsVectorFun1;
                  extrhs2 = &Test_A_plus_M_rhs2Fun1;
              } else if(!testcase.compare("2")) {
                  std::cout << "Test case 2 with vSol = [-y, x, 0]^T" << std::endl;
                  extvsol = &Test_A_plus_M_vSolVectorFun2;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun2_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun2_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun2_Gradient3;
                  extpsol = &Test_A_plus_M_pSolScalarFun2;
                  extrhs = &Test_A_plus_M_RhsVectorFun2;
              } else if(!testcase.compare("3")) {
                  std::cout << "Test case 3 with vSol = [-x*y-x*z+y^2+z^2, x^2-x*y-y*z+z^2, x^2-x*z+y^2-y*z]^T" << std::endl;
                  extvsol = &Test_A_plus_M_vSolVectorFun3;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun3_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun3_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun3_Gradient3;
                  extpsol = &Test_A_plus_M_pSolScalarFun3;
                  extrhs = &Test_A_plus_M_RhsVectorFun3;
              } else if(!testcase.compare("4")) {
                  std::cout << "Test case 4 (only mass matrix) with vSol = [-x*y-x*z+y^2+z^2, x^2-x*y-y*z+z^2, x^2-x*z+y^2-y*z]^T" << std::endl;
                  extvsol = &Test_A_plus_M_vSolVectorFun3;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun3_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun3_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun3_Gradient3;
                  extrhs = &Test_A_plus_M_RhsVectorFun3;
              } else if(!testcase.compare("5")) {
                  std::cout << "Test case 5 (1_in_spherical) with vSol = P[-z^2, y, x]^T" << std::endl;
                  extvsol = &Test_A_plus_M_vSolVectorFun5;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun5_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun5_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun5_Gradient3;
                  extpsol = &Test_A_plus_M_pSolScalarFun5;//&TestL_LagrangeFun24;
                  //SetupInterfaceRhsP1(mg, &rhs2, lset.Phi, lset.GetBndData(), Test_A_plus_M_rhs2Fun1);
                  extrhs = &Test_A_plus_M_RhsVectorFun5;
              } else if(!testcase.compare("0")) {
                  std::cout << "Test case 0 zero solution" << std::endl;
                  extvsol = &ZeroVectorFun;
                  extsol_grad1 = &ZeroVectorFun;
                  extsol_grad2 = &ZeroVectorFun;
                  extsol_grad3 = &ZeroVectorFun;
                  extpsol = &ZeroScalarFun    ;
                  extrhs = &ZeroVectorFun;
                  extrhs2 = &ZeroScalarFun;
                  InitVector(mg, vInit, &ZeroVectorFun);
              }
              else if(!testcase.compare("6")) {
                  std::cout << "Test case 6 with vSol = P[1, 0, 0]^T, pSol = y/(x^2+y^2+z^2)^(1/2)" << std::endl;
                  extvsol = &Test_A_plus_M_vSolVectorFun6;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun6_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun6_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun6_Gradient3;
                  extpsol = &Test_A_plus_M_pSolScalarFun6;
                  extrhs = &Test_A_plus_M_RhsVectorFun6;
                  extrhs2 = &Test_A_plus_M_rhs2Fun6;
                  InitVector(mg, vInit, &Test_A_plus_M_vSolVectorFun6);
              } else if(!testcase.compare("7")) {
                  std::cout << "Test case 7 with pSol = x*y^3/(x^2+y^2+z^2)^2+z/(x^2+y^2+z^2)^(1/2); " << std::endl;
                  extvsol = &Test_A_plus_M_vSolVectorFun7;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun7_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun7_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun7_Gradient3;
                  extpsol = &Test_A_plus_M_pSolScalarFun7;
                  extrhs = &Test_A_plus_M_RhsVectorFun7;
                  extrhs2 = &Test_A_plus_M_rhs2Fun7;
                  InitVector(mg, vInit, &Test_A_plus_M_vSolVectorFun7);
              }
              else if(!testcase.compare("8")) {
                  std::cout << "Test case 8 with vSol = P[-z^2, y, x]^T, pSol = (x*y^3/(x^2+y^2+z^2)^2+z/(x^2+y^2+z^2)^(1/2)); " << std::endl;
                  extvsol = &Test_A_plus_M_vSolVectorFun8;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun8_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun8_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun8_Gradient3;
                  extpsol = &Test_A_plus_M_pSolScalarFun8;
                  extrhs = &Test_A_plus_M_RhsVectorFun8;
                  extrhs2 = &Test_A_plus_M_rhs2Fun8;
                  InitVector(mg, vInit, &Test_A_plus_M_vSolVectorFun8);
              }
              else if(!testcase.compare("85")) {
                  std::cout << "Test case 85 with vSol = P normal_ext ( [-z^2, y, x]^T ), pSol = (x*y^3/(x^2+y^2+z^2)^2+z/(x^2+y^2+z^2)^(1/2)); " << std::endl;
                  extvsol = &Test_A_plus_M_vSolVectorFun85;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun85_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun85_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun85_Gradient3;
                  extpsol = &Test_A_plus_M_pSolScalarFun85;
                  extrhs = &Test_A_plus_M_RhsVectorFun85;
                  extrhs2 = &Test_A_plus_M_rhs2Fun85;
                  InitVector(mg, vInit, &Test_A_plus_M_vSolVectorFun85);
              }
              else if(!testcase.compare("9")) {
                  std::cout << "Test case 9 with vSol = n x grad(y*z), pSol = (x*y^3/(x^2+y^2+z^2)^2+z/(x^2+y^2+z^2)^(1/2)); " << std::endl;
                  extvsol = &Test_A_plus_M_vSolVectorFun9;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun9_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun9_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun9_Gradient3;
                  extpsol = &Test_A_plus_M_pSolScalarFun9;
                  extrhs = &Test_A_plus_M_RhsVectorFun9;
                  extrhs2 = &Test_A_plus_M_rhs2Fun9;
                  InitVector(mg, vInit, &Test_A_plus_M_vSolVectorFun9);
              }
              else if(!testcase.compare("95")) {
                  std::cout << "Test case 95 with vSol = P[1,0,0]^T, pSol = (x*y^3/(x^2+y^2+z^2)^2+z/(x^2+y^2+z^2)^(1/2)); " << std::endl;
                  extvsol = &Test_A_plus_M_vSolVectorFun95;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun95_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun95_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun95_Gradient3;
                  extpsol = &Test_A_plus_M_pSolScalarFun9;
                  extrhs = &Test_A_plus_M_RhsVectorFun95;
                  extrhs2 = &Test_A_plus_M_rhs2Fun95;
                  InitVector(mg, vInit, &Test_A_plus_M_vSolVectorFun95);
              }
              else if(!testcase.compare("10")) {
                  std::cout << "Test case 10 initial velocity = (2+3)spherical harmonics; " << std::endl;
                  extvsol = &Test_A_plus_M_vSolVectorFun10;
                  extvinit = &Test_A_plus_M_vInitVectorFun10;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun10_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun10_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun10_Gradient3;
                  extpsol = &Test_A_plus_M_pSolScalarFun10;
                  extrhs = &Test_A_plus_M_RhsVectorFun10;
                  extrhs2 = &Test_A_plus_M_rhs2Fun10;
                  InitVector(mg, vInit, extvinit);
              }
              else if(!testcase.compare("11")) {
                  std::cout << "Test case 11 = (6 + convection) with vSol = P[1, 0, 0]^T, pSol = y/(x^2+y^2+z^2)^(1/2)" << std::endl;
                  extvsol = &Test_A_plus_M_vSolVectorFun6;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun6_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun6_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun6_Gradient3;
                  extpsol = &Test_A_plus_M_pSolScalarFun6;
                  extrhs2 = &Test_A_plus_M_rhs2Fun6;
              }
              else if(!testcase.compare("12")) {
                  std::cout << "Test case 12 = (11 + Bernouli pressure) with vSol = P[1, 0, 0]^T, pSol = vSol^2/2 + y/(x^2+y^2+z^2)^(1/2)" << std::endl;
                  extvsol = &Test_A_plus_M_vSolVectorFun6;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun6_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun6_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun6_Gradient3;
                  extpsol = &Test_A_plus_M_pBerSolScalarFun12;
                  extrhs2 = &Test_A_plus_M_rhs2Fun6;
              }
              else if(!testcase.compare("13")) {
                  std::cout << "Test case 13 = (8 + convection) with vSol = P[-z^2, y, x]^T, pSol = (x*y^3/(x^2+y^2+z^2)^2+z/(x^2+y^2+z^2)^(1/2)); " << std::endl;
                  extvsol = &Test_A_plus_M_vSolVectorFun8;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun8_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun8_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun8_Gradient3;
                  extrhs = &Test_A_plus_M_RhsVectorFun13;
                  extrhs2 = &Test_A_plus_M_rhs2Fun8;
                  extpsol = &Test_A_plus_M_pSolScalarFun8;
                  InitVector(mg, vInit, &Test_A_plus_M_vSolVectorFun8);

              }
              else if(!testcase.compare("14")) {
                  std::cout << "Test case 14 = (13 + Bernouli pressure) with vSol = P[-z^2, y, x]^T, pSol = vSol^2/2 +  (x*y^3/(x^2+y^2+z^2)^2+z/(x^2+y^2+z^2)^(1/2)); " << std::endl;
                  extvsol = &Test_A_plus_M_vSolVectorFun8;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun8_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun8_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun8_Gradient3;
                  extpsol = &Test_A_plus_M_pBerSolScalarFun14;
                  extrhs2 = &Test_A_plus_M_rhs2Fun8;
              }
              else if(!testcase.compare("15")) {
                  std::cout << "Test case 15 = (13 - Bernouli pressure) with vSol = P[-z^2, y, x]^T, pSol = - vSol^2/2 +  (x*y^3/(x^2+y^2+z^2)^2+z/(x^2+y^2+z^2)^(1/2)); " << std::endl;
                  extvsol = &Test_A_plus_M_vSolVectorFun8;
                  extrhs  = &Test_A_plus_M_RhsVectorFun13;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun8_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun8_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun8_Gradient3;
                  extpsol = &Test_A_plus_M_pBerSolScalarFun15;
                  extrhs2 = &Test_A_plus_M_rhs2Fun8;
              }
              else if(!testcase.compare("16")) {
                  std::cout << "Test case 16 = from Killing to Killing; " << std::endl;
                  extvsol = &Test_A_plus_M_vSolVectorFun16;
                  extpsol = &Test_A_plus_M_pSolScalarFun10;
                  extrhs  = //&ZeroVectorFun;
                          &Test_A_plus_M_rhs_convFun16;
                  extrhs2  = &ZeroScalarFun;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun10_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun10_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun10_Gradient3;
                  InitVector(mg, vInit, &Test_A_plus_M_vSolVectorFun16);
              }
              else if(!testcase.compare("17")) {
                  std::cout << "Test case 17 = swirled Killing; " << std::endl;
                  ParameterNS::nu = P.get<double>("SurfNavStokes.kinematic_viscosity");
                  extvsol = &Test_A_plus_M_vSolVectorFun17;
                  extpsol = &Test_A_plus_M_pSolScalarFun10;
                  extrhs  = &Test_A_plus_M_RhsVectorFun17  ;//viscosity is here!
                  extrhs2  = &ZeroScalarFun;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun17_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun17_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun17_Gradient3;
                  InitVector(mg, vInit, &Test_A_plus_M_vSolVectorFun17);
              }
              else if(!testcase.compare("18")) {
                  std::cout << "Test case 18 = Kelvin-Helmholtz; " << std::endl;
                  ParameterNS::nu = P.get<double>("SurfNavStokes.kinematic_viscosity");
                  extvsol = &Test_A_plus_M_vSolVectorFun18;
                  extpsol = &Test_A_plus_M_pSolScalarFun10;
                  extrhs  = &ZeroVectorFun  ;
                  extrhs2  = &ZeroScalarFun;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun10_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun10_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun10_Gradient3;
                  InitVector(mg, vInit, &Test_A_plus_M_vSolVectorFun18);
              }
              else if(!testcase.compare("19")) {
                  std::cout << "Test case 19 = vorteces; " << std::endl;
                  ParameterNS::nu = P.get<double>("SurfNavStokes.kinematic_viscosity");
                  extvsol = &ZeroVectorFun;
                  extpsol = &ZeroScalarFun;
                  extrhs  = &ZeroVectorFun  ;
                  extrhs2  = &ZeroScalarFun;
                  extsol_grad1 = &ZeroVectorFun;
                  extsol_grad2 = &ZeroVectorFun;
                  extsol_grad3 = &ZeroVectorFun;
                  InitVector(mg, vInit, &Test_A_plus_M_disturbed_test18);
              }
              else if(!testcase.compare("20")) {
                  std::cout << "Test case 20 = test 17 + swirling depends on x; " << std::endl;
                  ParameterNS::nu = P.get<double>("SurfNavStokes.kinematic_viscosity");
                  extvsol = &Test_A_plus_M_vSolVectorFun20;
                  extpsol = &ZeroScalarFun;
                  extcurlsol = &Test_A_plus_M_curlFun20;
                  extrhs  = &Test_A_plus_M_RhsVectorFun20  ;//viscosity is here!
                  extrhs2  = &Test_A_plus_M_rhs2Fun20;
                  extsol_grad1 = &Test_A_plus_M_vSolVectorFun20_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_vSolVectorFun20_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_vSolVectorFun20_Gradient3;
                  InitVector(mg, vInit, &Test_A_plus_M_vSolVectorFun20);
              }
        } else if (!levelset_fun_str.compare("xy_plane")) {
              if( !testcase.compare("1")) {
                  std::cout << "Test case 1 with vSol = [sin(Pi*x)sin(Pi*y), sin(Pi*x)sin(Pi*y), 0]^T" << std::endl;
                  extvsol = &Test_A_plus_M_xy_plane_vSolVectorFun1;
                  extsol_grad1 = &Test_A_plus_M_xy_plane_vSolVectorFun1_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_xy_plane_vSolVectorFun1_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_xy_plane_vSolVectorFun1_Gradient3;
              } else if(!testcase.compare("2")) {
                  std::cout << "Test case 2 with vSol = [sin(Pi/2*x)sin(Pi/2*y), sin(Pi/2*x)sin(Pi/2*y), 0]^T" << std::endl;
                  extvsol = &Test_A_plus_M_xy_plane_vSolVectorFun2;
                  extsol_grad1 = &Test_A_plus_M_xy_plane_vSolVectorFun2_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_xy_plane_vSolVectorFun2_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_xy_plane_vSolVectorFun2_Gradient3;
              } else if(!testcase.compare("3")) {
                  std::cout << "Test case 3 with vSol = [sin(Pi/4*x)sin(Pi/4*y), sin(Pi/4*x)sin(Pi/4*y), 0]^T" << std::endl;
                  extvsol = &Test_A_plus_M_xy_plane_vSolVectorFun3;
                  extsol_grad1 = &Test_A_plus_M_xy_plane_vSolVectorFun3_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_xy_plane_vSolVectorFun3_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_xy_plane_vSolVectorFun3_Gradient3;
              } else if(!testcase.compare("Zero")) {
                  extvsol = &ZeroVectorFun;
                  extsol_grad1 = &ZeroVectorFun;
                  extsol_grad2 = &ZeroVectorFun;
                  extsol_grad3 = &ZeroVectorFun;
              } else if(!testcase.compare("4")) {
                  std::cout << "Test case 4 with vSol = [(x^2-4)*(y^2-4), (x^2-4)*(y^2-4), 0]^T" << std::endl;
                  extvsol = &Test_A_plus_M_xy_plane_vSolVectorFun4;
                  extsol_grad1 = &Test_A_plus_M_xy_plane_vSolVectorFun4_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_xy_plane_vSolVectorFun4_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_xy_plane_vSolVectorFun4_Gradient3;
              }
        } else if( !levelset_fun_str.compare("tilted_plane")) {
              if( !testcase.compare("1")) {
                  std::cout << "Test case 1 with vSol = [4/5*sin(Pi*x)sin(Pi*z), 2/5*sin(Pi*x)sin(Pi*z), sin(Pi*x)sin(Pi*z)]^T" << std::endl;
                  extvsol = &Test_A_plus_M_tilted_plane_vSolVectorFun1;
                  extsol_grad1 = &Test_A_plus_M_tilted_plane_vSolVectorFun1_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_tilted_plane_vSolVectorFun1_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_tilted_plane_vSolVectorFun1_Gradient3;
              } else if( !testcase.compare("2")) {
                  std::cout << "Test case 2 with vSol = [4/5*(x^2-4)*(z^2-4), 2/5*(x^2-4)*(z^2-4), (x^2-4)*(z^2-4)]^T" << std::endl;
                  extvsol = &Test_A_plus_M_tilted_plane_vSolVectorFun2;
                  extsol_grad1 = &Test_A_plus_M_tilted_plane_vSolVectorFun2_Gradient1;
                  extsol_grad2 = &Test_A_plus_M_tilted_plane_vSolVectorFun2_Gradient2;
                  extsol_grad3 = &Test_A_plus_M_tilted_plane_vSolVectorFun2_Gradient3;
              }
        } else if( !levelset_fun_str.compare("tilted_plane_xy")) {
            if (!testcase.compare("1")) {
                std::cout
                        << "Test case 1 with vSol = [4/17*sin(Pi*y)*sin(Pi*z), sin(Pi*y)*sin(Pi*z), 16/17*sin(Pi*y)*sin(Pi*z)]^T"
                        << std::endl;
                extvsol = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun1;
                extsol_grad1 = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun1_Gradient1;
                extsol_grad2 = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun1_Gradient2;
                extsol_grad3 = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun1_Gradient3;
            } else if (!testcase.compare("Zero")) {
                extvsol = &ZeroVectorFun;
                extsol_grad1 = &ZeroVectorFun;
                extsol_grad2 = &ZeroVectorFun;
                extsol_grad3 = &ZeroVectorFun;
            } else if (!testcase.compare("2")) {
                std::cout << "Test case 2 with vSol = [4/17*(z^2-4)*(y^2-4), (z^2-4)*(y^2-4), 16/17*(z^2-4)*(y^2-4)]^T"
                          << std::endl;
                extvsol = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun2;
                extsol_grad1 = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun2_Gradient1;
                extsol_grad2 = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun2_Gradient2;
                extsol_grad3 = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun2_Gradient3;
            } else if (!testcase.compare("3")) {
                std::cout
                        << "Test case 3 with vSol = [4/17*(z^2-4)^2*(y^2-4)^2, (z^2-4)^2*(y^2-4)^2, 16/17*(z^2-4)^2*(y^2-4)^2]^T"
                        << std::endl;
                extvsol = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun3;
                extsol_grad1 = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun3_Gradient1;
                extsol_grad2 = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun3_Gradient2;
                extsol_grad3 = &Test_A_plus_M_tilted_plane_xy_vSolVectorFun3_Gradient3;
            }
        } else if( !levelset_fun_str.compare("cube_madeof_edges")) {
            if (!testcase.compare("1")) {
                std::cout << "Test case 1 with nonzero divergence" << std::endl;
                extrhs = &Test_cube_madeof_edges_RhsVectorFun1;
                extrhs2 = &Test_cube_madeof_edges_rhs2Fun1;
                extvsol = &Test_cube_madeof_edges_vSolVectorFun1;
                extsol_grad1 = &Test_cube_madeof_edges_vSolVectorFun1_Gradient1;
                extsol_grad2 = &Test_cube_madeof_edges_vSolVectorFun1_Gradient2;
                extsol_grad3 = &Test_cube_madeof_edges_vSolVectorFun1_Gradient3;
                extpsol = &Test_cube_madeof_edges_pSolScalarFun1;
            }
        } else if( !levelset_fun_str.compare("torus")) {
            if (!testcase.compare("1")) {
                std::cout << "Test case 1 = Arnold example; " << std::endl;
                ParameterNS::nu = P.get<double>("SurfNavStokes.kinematic_viscosity");
                extvsol = &Torus_vSolVectorArnold;
                extpsol = &ZeroScalarFun;
                extrhs = &ZeroVectorFun;
                extrhs2 = &ZeroScalarFun;
                extsol_grad1 = &ZeroVectorFun;
                extsol_grad2 = &ZeroVectorFun;
                extsol_grad3 = &ZeroVectorFun;
                InitVector(mg, vInit, extvsol);
            }
            if (!testcase.compare("2")) {
                std::cout << "Test case 2 = linear combinations; " << std::endl;
                ParameterNS::nu = P.get<double>("SurfNavStokes.kinematic_viscosity");
                extvsol = &Torus_vSolVectorHarmonic;
                extpsol = &ZeroScalarFun;
                extrhs = &ZeroVectorFun;
                extrhs2 = &ZeroScalarFun;
                extsol_grad1 = &ZeroVectorFun;
                extsol_grad2 = &ZeroVectorFun;
                extsol_grad3 = &ZeroVectorFun;
                InitVector(mg, vInit, extvsol);
            }
        } else if (!levelset_fun_str.compare("torus_flower")) {
            if (!testcase.compare("1")) {
                std::cout << "Test case 1 = harmonic; " << std::endl;
                ParameterNS::nu = P.get<double>("SurfNavStokes.kinematic_viscosity");
                extvsol = &Torus_vSolVectorHarmonicPhi;
                extpsol = &ZeroScalarFun;
                extrhs = &ZeroVectorFun;
                extrhs2 = &ZeroScalarFun;
                extsol_grad1 = &ZeroVectorFun;
                extsol_grad2 = &ZeroVectorFun;
                extsol_grad3 = &ZeroVectorFun;
                InitVector(mg, vInit, extvsol);
            }
        }
        // setup matrices
        DROPS::MatDescCL Omega, N, NT, D, L, L_stab;
        auto& A = system.A;
        auto& A_stab = system.A_stab;
        auto& B = system.B;
        auto& M = system.M;
        auto& S = system.S;
        auto& Schur = system.Schur;
        auto& Schur_stab = system.Schur_stab;
        auto& Schur_normal_stab = system.Schur_normal_stab;
        param.input.f = extrhs;
        param.input.g = extrhs2;
        MatrixCL Schur_hat;
        if( !FE.compare("P2P1")) {
            A.SetIdx( &ifaceVecP2idx, &ifaceVecP2idx);
            A_stab.SetIdx( &ifaceVecP2idx, &ifaceVecP2idx);
            N.SetIdx( &ifaceVecP2idx, &ifaceVecP2idx);
            B.SetIdx( &ifaceP1idx, &ifaceVecP2idx);
            M.SetIdx( &ifaceVecP2idx, &ifaceVecP2idx);
            D.SetIdx( &ifaceVecP2idx, &ifaceVecP2idx);
            S.SetIdx( &ifaceVecP2idx, &ifaceVecP2idx);
            L.SetIdx( &ifaceP1idx, &ifaceVecP2idx);
            L_stab.SetIdx( &ifaceP1idx, &ifaceVecP2idx);
            Schur.SetIdx(&ifaceP1idx, &ifaceP1idx);
            Schur_stab.SetIdx(&ifaceP1idx, &ifaceP1idx);
            Schur_normal_stab.SetIdx(&ifaceP1idx, &ifaceP1idx);
            SetupStokesIF_P2P1(mg, lset, &system, &param);
        } else if( !FE.compare("P1P1")) {
            A.SetIdx( &ifaceVecP1idx, &ifaceVecP1idx);
            A_stab.SetIdx( &ifaceVecP1idx, &ifaceVecP1idx);
            B.SetIdx( &ifaceP1idx, &ifaceVecP1idx);
            Omega.SetIdx( &ifaceP1idx, &ifaceVecP1idx);
            N.SetIdx( &ifaceVecP1idx, &ifaceVecP1idx);
            NT.SetIdx( &ifaceVecP1idx, &ifaceVecP1idx);
            M.SetIdx( &ifaceVecP1idx, &ifaceVecP1idx);
            D.SetIdx( &ifaceVecP1idx, &ifaceVecP1idx);
            S.SetIdx( &ifaceVecP1idx, &ifaceVecP1idx);
            L.SetIdx( &ifaceP1idx, &ifaceVecP1idx);
            L_stab.SetIdx( &ifaceP1idx, &ifaceVecP1idx);
            Schur.SetIdx(&ifaceP1idx, &ifaceP1idx);
            Schur_stab.SetIdx(&ifaceP1idx, &ifaceP1idx);
            Schur_normal_stab.SetIdx(&ifaceP1idx, &ifaceP1idx);
            SetupNavierStokesIF_P1P1(mg, &A, &A_stab, &B, &Omega, &N, &NT, &M, &D, &S, &L, &L_stab, &Schur, &Schur_stab, &Schur_normal_stab, lset, v, vbnd, &param);
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
            SetupStokesIF_P2P2(mg, &A, &A_stab, &B, &M, &S, &L, &L_stab, &Schur, &Schur_stab, lset.Phi, lset.GetBndData(), &param);
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
            SetupStokesIF_P1P2(mg, &A, &A_stab, &B, &M, &S, &L, &L_stab, &Schur, &Schur_stab, lset.Phi, lset.GetBndData(), &param);
        }
        std::cout << "stiffness mtx is              " << A.Data.num_rows() << " * " << A.Data.num_cols() << '\n'
                  << "velocity soln size is         " << v.Data.size() << '\n'
                  << "exact velocity interp size is " << vSol.Data.size() << '\n'
                  << "pre mass mtx is               " << Schur.Data.num_rows() << " * " << Schur.Data.num_cols() << '\n'
                  << "pressure soln size is         " << p.Data.size() << '\n'
                  << "exact pressure interp size is " << pSol.Data.size() << '\n'
                  << "f size is                     " << fRHS.Data.size() << '\n'
                  << "g size is                     " << gRHS.Data.size() << '\n';
        std::string dirname  = P.get<std::string>("Output.Directory") + "/" + model + "_" + levelset_fun_str + "/" + "test" + testcase + "_h=" + std::to_string(float(h));
        // export matrices to test inf-sup constant for P1-P1 / P2-P1
        if (param.input.computeMatrices && P.get<bool>("SurfNavStokes.ExportMatrices")) {
            std::cout << "test inf-sup constant\n";
            MatrixCL A_final;
            A_final.LinComb(1., A.Data, 1., M.Data, tau_u, S.Data, rho_u, A_stab.Data);
            auto C_full = Schur_stab.Data;
            C_full *= rho_p;
            auto C_n = Schur_normal_stab.Data;
            C_n *= rho_p;
            auto& B_final = B.Data;
            auto& M_final = Schur.Data;
            auto surfName = P.get<std::string>("Levelset.case") +"_shift=" + P.get<std::string>("Levelset.ShiftNorm", "0");
            auto format = P.get<std::string>("SurfNavStokes.ExportMatricesFormat") == "mtx" ? ".mtx" : ".mat";
            auto expFunc = format == ".mtx" ? &MatrixCL::exportMTX : &MatrixCL::exportMAT;
            std::cout << "exporting " << format << " matrices to " + dirname + "*\n";
            (A_final.*expFunc)(dirname + "/A" + format);
            (B_final.*expFunc)(dirname + "/B" + format);
            (M_final.*expFunc)(dirname + "/M" + format);
            (C_full.*expFunc)(dirname + "/C_full" + format);
            (C_n.*expFunc)(dirname + "/C_full" + format);
        }
        if (P.get<bool>("SurfNavStokes.ComputeNormalErr") || P.get<bool>("SurfNavStokes.ComputeShapeErr")) {
            std::ofstream log(dirname + "/normal_and_shape_errs_m=" + std::to_string(param.input.numbOfVirtualSubEdges) + ".txt");
            if (P.get<bool>("SurfNavStokes.ComputeNormalErr")) {
                log << "levelset normal L2 error is: " << sqrt(param.output.normalErrSq.lvset) << '\n';
                log << "patch    normal L2 error is: " << sqrt(param.output.normalErrSq.patch) << '\n';
            }
            if (P.get<bool>("SurfNavStokes.ComputeShapeErr"))
                log << "shape operator L2 error is: " << sqrt(param.output.shapeErrSq) << '\n';
        }
        if (!P.get<bool>("SurfNavStokes.Solve")) return 0;
        // Schur precond
        if (P.get<std::string>("SurfNavStokes.stab") == "full")
            Schur_hat.LinComb(1., Schur.Data, rho, Schur_stab.Data);
        else if (P.get<std::string>("SurfNavStokes.stab") == "normal")
            Schur_hat.LinComb(1., Schur.Data, rho, Schur_normal_stab.Data);
        else {
            std::cout << "WRN: no stabilization is used! Schur precond = I\n";
            Schur_hat = std::valarray<double>(1., Schur.Data.num_rows());
        }
        // construct preconditioners
        typedef SSORPcCL      SymmPcPcT;
        typedef PCGSolverCL<SymmPcPcT> PCGSolverT;
        typedef SolverAsPreCL<PCGSolverT> PCGPcT;
        SymmPcPcT symmPcPc_;
        std::stringstream Astream;
        PCGSolverT PCGSolver_( symmPcPc_, P.get<int>("Solver.PcAIter"), P.get<double>("Solver.PcATol"), true, &Astream);
        PCGPcT PCGPc_( PCGSolver_);//velocity preconditioner
        std::stringstream Schurstream;
        PCGSolverT SchurPCGSolver (symmPcPc_, P.get<int>("Solver.PcBIter"), P.get<double>("Solver.PcBTol"), true, &Schurstream);
        SchurPreBaseCL *spc_ = new SurfaceLaplacePreCL<PCGSolverT>( Schur_hat, SchurPCGSolver);//pressure preconditioner
        //construct symmteric block iterative solver
        typedef BlockPreCL<ExpensivePreBaseCL, SchurPreBaseCL, DiagSpdBlockPreCL>  DiagBlockPcT;
        typedef PLanczosONBCL<VectorCL, DiagBlockPcT> LanczosT;
        typedef PMResSolverCL<LanczosT> MinResT;
        DiagBlockPcT    *DBlock_ = new DiagBlockPcT( PCGPc_, *spc_);
        LanczosT *lanczos_ = new LanczosT( *DBlock_);
        MinResT *MinRes_ = new MinResT( *lanczos_,  P.get<int>("Solver.Iter"), P.get<double>("Solver.Tol"), /*relative*/ false);
        BlockMatrixSolverCL<MinResT> *stokessolver= new BlockMatrixSolverCL<MinResT>( *MinRes_);
        // construct preconditioners for nonsymmetric case
        typedef /*JACPcCL*/ SSORPcCL NonSymmPcPcT;
        typedef GMResSolverCL<NonSymmPcPcT> GMResSolverT;
        typedef SolverAsPreCL<GMResSolverT> GMResPcT;
        NonSymmPcPcT nonsymmPcPc_;
        std::stringstream PCstream;
        GMResSolverT GMResSolver_( nonsymmPcPc_, P.get<int>("Solver.PcAIter"), P.get<int>("Solver.PcAIter"), P.get<double>("Solver.PcATol"),
                 /*bool relative=*/ true, /*bool calculate2norm=*/ false, /*PreMethGMRES method=*/ LeftPreconditioning,
                     /*bool mod =*/ true,      /*bool useModGS =*/ false, &PCstream);
        GMResPcT GMResPc_nonsym( GMResSolver_);
        //setup iterative block solver for nonsymmetric case
        typedef GMResSolverCL<DiagBlockPcT> GMResT;
        DiagBlockPcT    *DBlock_nonsym = new DiagBlockPcT( GMResPc_nonsym, *spc_);
        GMResT *GMRes_ = new GMResT( *DBlock_nonsym, P.get<int>("Solver.Iter"), P.get<int>("Solver.Iter"), P.get<double>("Solver.Tol"),
                 /*bool relative=*/ false, /*bool calculate2norm=*/ false, /*PreMethGMRES method=*/ LeftPreconditioning,
                      /*bool mod =*/ true,      /*bool useModGS =*/ false);
        BlockMatrixSolverCL<GMResT> *nonsym_stokessolver= new BlockMatrixSolverCL<GMResT>( *GMRes_);
        //pick correct solver
        StokesSolverBaseCL *Solver;
        if (  P.get<std::string>("SurfNavStokes.nonlinear_term") == "convective"
           || P.get<std::string>("SurfNavStokes.nonlinear_term") == "rotational"
           || P.get<std::string>("SurfNavStokes.nonlinear_term") == "strain"
           || P.get<std::string>("SurfNavStokes.nonlinear_term") == "skew-symmetric"
           || P.get<std::string>("SurfNavStokes.nonlinear_term") == "conservative"
           || P.get<std::string>("SurfNavStokes.nonlinear_term") == "EMAC"
            ) Solver = nonsym_stokessolver;
        else  Solver = stokessolver;
        // set up interpolants for exact solution
        InitScalar(mg, pSol, extpsol);
        InitVector(mg, vSol, extvsol);
        if( !model.compare("NavierStokes")) {
            //set up discrete vectors and matrices
            VectorCL id, id2, random ;
            id.resize(fRHS.Data.size(), 1);
            random.resize(fRHS.Data.size(), 0);
            for (int i=0; i<fRHS.Data.size(); i++)
                random[i] = - 5 + (rand() % 10);
            id2.resize(gRHS.Data.size(), 1);
            MatrixCL Adyn, Ahat, Bhat, Chat, Mhat;
            MatrixCL BTranspose, NTranspose;
            DROPS::VecDescCL instantrhs;
            DROPS::VecDescCL vxtent, pxtent;
            //set up output
            std::string filename = "test" + testcase + "_";
            std::cout << "dirname: " << dirname << std::endl;
            //set up VTK output
            VTKOutCL* vtkwriter = nullptr;
            vtkwriter = new VTKOutCL(mg, "DROPS data", (int)P.get<double>("Time.NumSteps")/P.get<int>("Output.every timestep"), dirname , filename, "none", 0);
            vtkwriter->Register(make_VTKScalar(lset.GetSolution(), "level-set"));
            vtkwriter->Register(make_VTKIfaceVector(mg, vSol, "velSol", velFE, vbnd));
            vtkwriter->Register(make_VTKIfaceVector(mg, fRHS, "rhs", velFE, vbnd));
            vtkwriter->Register(make_VTKIfaceVector(mg, v, "velocity", velFE, vbnd));
            vtkwriter->Register(make_VTKIfaceScalar(mg, p, "pressure", /*prFE,*/ pbnd));
            vtkwriter->Register(make_VTKIfaceScalar(mg, curl, "vorticity", /*prFE,*/ pbnd));
            vtkwriter->Register(make_VTKIfaceScalar(mg, curlSol, "vortSol", /*prFE,*/ pbnd));
            vtkwriter->Register(make_VTKIfaceScalar(mg, gRHS, "gRHS", /*prFE,*/ pbnd));
            vtkwriter->Register(make_VTKIfaceScalar(mg, pSol, "pSol", /*prFE,*/ pbnd));

            //set up a txt file for custom time output
            std::ofstream log_solo( dirname +"/"
                                  + "tau_u="+ std::to_string(float(tau_u))+ "_"
                                  + "l="  + P.get<std::string>("Mesh.AdaptRef.FinestLevel")
                                  + "_nu=" + P.get<std::string>("SurfNavStokes.kinematic_viscosity")
                                  + "_Plot_" + filename
                                  + P.get<std::string>("SurfNavStokes.nonlinear_term") + "_"
                                  + P.get<std::string>("SurfNavStokes.instationary") + "="
                                  + std::to_string(float(tau)) + ".txt");
            log_solo << "Time" <<  "\tKinetic\t" << "\tMomentum\t" << "\tWork" <<  std::endl;

            //set up a txt file for error time output
            std::ofstream log_error( dirname +"/"
                    + "tau_u="+ std::to_string(float(tau_u))+ "_"
                    + "l="  + P.get<std::string>("Mesh.AdaptRef.FinestLevel")
                    + "_nu=" + P.get<std::string>("SurfNavStokes.kinematic_viscosity")
                    + "_Error_" + filename
                    + P.get<std::string>("SurfNavStokes.nonlinear_term") + "_"
                    + P.get<std::string>("SurfNavStokes.instationary") + "="
                    + std::to_string(float(tau)) + ".txt");
            log_error << "Time" << "\tL_2(uT-ext_uT)\t" <<  "\tadvH_1(u-ext_u)\t" <<  "\tsurfH_1(u-ext_u)\t" <<  "\tH_1(u-ext_u)\t" << "\tL_2(uN)\t" << "\tL_2(p-ext_p)" << std::endl;
            double mu = P.get<double>("SurfNavStokes.kinematic_viscosity");
            std::cout << "viscosity: " << mu << std::endl;
            std::cout << "test is instationary: " << P.get<std::string>("SurfNavStokes.instationary") << '\n';

            //NAVIER-STOKES starts here
            if ( P.get<std::string>("SurfNavStokes.instationary") != "none" ) {
                //initial velocity for backward Euler formula
                v=vInit;
                v_old=vInit;

                Mhat.LinComb(1., Schur.Data, rho_p, Schur_normal_stab.Data);
                //renormalization to fullfill \int rhs = 0 in discrete sence

                curl_proj.Data = Omega.Data*v.Data  - ( dot(Omega.Data*v.Data,id2) / dot(id2,id2) ) * id2;
                //GMResSolver_.Solve(Mhat, curl.Data,curl_proj.Data , gRHS.Data);
                InitScalar(mg, p, extpsol);
                //output of initial data and exact solutions to vtk and custom
                vtkwriter->Write(0);//send v to vtk
                vxtent.SetIdx( &vecP1idx);
                Extend(mg, vInit, vxtent);
                BndDataCL<Point3DCL> bndvec = vbnd;
                BndDataCL<double> bndscalar = pbnd;
                double l2norm = L2_Vector_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, bndvec, vxtent), &ZeroVectorFun);
                double mom = dot(M.Data * vInit.Data,id);
                double work = dot(M.Data * vInit.Data,fRHS.Data);
                log_solo << std::to_string(0.0) << "\t" << std::to_string((float)(l2norm*l2norm*0.5)) << "\t" << std::to_string((float)(mom))<< "\t" << std::to_string((float)(work)) << std::endl;

                bool experimental=false;

                std::ofstream log( dirname +"/"+ filename   + "_time="+std::to_string(0) + ".txt");
                log_error <<  std::to_string((float)0) << "\t";

                if (!velFE.compare("P1")) {
                    vxtent.SetIdx(&vecP1idx);
                    Extend(mg, v, vxtent);
                    BndDataCL<Point3DCL> bndvec = vbnd;
                    double error_L2 = L2_Vector_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, bndvec, vxtent), extvsol, 0);
                    log << "The L2-Norm of v - vSol is: " << error_L2 << std::endl;
                    log_error << std::to_string((float) error_L2) << "\t";
                    double l2norm = L2_Vector_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, bndvec, vxtent),
                                                    &ZeroVectorFun, 0);
                    log << "The L2-Norm of v  is: " << l2norm << std::endl;

                    double H1(0.), surfH1(0.), advanced_surfH1(0.), normal_velocity(0.);
                    H1_Vector_error_P1(mg, lset.Phi, lset.GetBndData(), v, vbnd, extvsol, extsol_grad1, extsol_grad2,
                                       extsol_grad3, tau_u, H1, surfH1, advanced_surfH1, normal_velocity, 0);
                    //        std::cout << "The H1-Norm of v - vSol is: " << H1 << std::endl;
                    //        std::cout << "The surfH1-Norm of v - vSol is: " << surfH1 << std::endl;
                    log << "The advanced surfH1-Norm of v - vSol is: " << advanced_surfH1 << std::endl;
                    log_error << std::to_string((float) advanced_surfH1) << "\t";
                    log_error << std::to_string((float) surfH1) << "\t";
                    log_error << std::to_string((float) H1) << "\t";
                    //log << "The U-Norm of v - vSol is: " << std::sqrt(advanced_surfH1*advanced_surfH1 + rho_u*dot(A_stab.Data*v.Data, v.Data)) << std::endl;
                    log << "The L2-Norm of v * n is: " << normal_velocity << std::endl;
                    log_error << std::to_string((float) normal_velocity) << "\t";
                }
                if (!prFE.compare("P1")) {
                    pxtent.SetIdx( &P1FEidx);
                    Extend(mg, p, pxtent);
                    BndDataCL<double> bndscalar = pbnd;
                    double L2_Lagrange =L2_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, bndscalar, pxtent), extpsol, 0);
                    log << "The L2-Norm of p - pSol is: " << L2_Lagrange << std::endl;
                    log_error <<  std::to_string((float)L2_Lagrange) << std::endl;
                    log << "The L2-Norm of p  is: " << L2_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, bndscalar, pxtent), &ZeroScalarFun);
                    //log << "The M-Norm of p - pSol is: " << std::sqrt(L2_Lagrange*L2_Lagrange + rho_p*dot(Schur_stab.Data*p.Data, p.Data)) << std::endl;
                }

                for (int i = 0; i < P.get<double>("Time.NumSteps"); i++)
                {
                    InitVector(mg, vSol, extvsol,(i+1)*tau);
                    InitScalar(mg, curlSol, extcurlsol,(i+1)*tau);
                    InitScalar(mg, pSol, extpsol,(i+1)*tau);

                    //current timestep logfile
                    std::ofstream log( dirname +"/"+ filename   + "_time="+std::to_string((i+1)*tau) + ".txt");

                    //right constant for temporal BDF
                    double c;
                    if ( P.get<std::string>("SurfNavStokes.instationary") == "BDF1")
                    {
                        //leading time term mass coeff
                        c=1.0;
                        //interpolated from previous timesteps, "wind"
                        v_aux.Data = 1.0*v.Data;
                    }
                    else if ( P.get<std::string>("SurfNavStokes.instationary") == "BDF2")
                    {
                        //leading time term mass coeff
                        if (i!=0) {
                            c = 3.0 / 2.0;
                        }
                        else c=1.0;
                        //interpolated from previous timesteps, "wind"
                        v_aux.Data = 2.0*v.Data-v_old.Data;
                    }
                    else
                    {
                        std::cout << "problem is instationary, pick BDF1 or BDF2" << std::endl;
                        return 0;
                    }

                    //set up current matrices based in previous timestep velocity and levelset
                    SetupNavierStokesIF_P1P1(mg, &A, &A_stab, &B, &Omega, &N, &NT, &M, &D, &S,  &L, &L_stab, &Schur, &Schur_stab, &Schur_normal_stab, lset, v_aux, vbnd, &param);

                    //construct final matrices for a linear solver
                    Ahat.LinComb(mu, A.Data, c/tau, M.Data, tau_u, S.Data, rho_u, A_stab.Data);
                    Bhat.LinComb(1., B.Data, 0., B.Data);
                    if (P.get<std::string>("SurfNavStokes.stab") == "full")
                        Chat.LinComb(0, Schur.Data, -rho_p, Schur_stab.Data);
                    else if (P.get<std::string>("SurfNavStokes.stab") == "normal")
                        Chat.LinComb(0, Schur.Data, -rho_p, Schur_normal_stab.Data);
                    else {
                        std::cout << "WRN: no stabilization is used!\n";
                        Chat.LinComb(0, Schur.Data, 0, Schur.Data);
                    }
                    transpose(B.Data, BTranspose);

                    //pick the form of nonlinear term
                    if ( P.get<std::string>("SurfNavStokes.nonlinear_term") == "convective" )
                    {
                        std::cout << "nonlinear_term: convective " << std::endl;
                        Adyn.LinComb(1, Ahat, 1, N.Data);
                    }
                    else if (P.get<std::string>("SurfNavStokes.nonlinear_term") == "rotational")
                    {
                        std::cout << "nonlinear_term: rotational " << std::endl;

                        Adyn.LinComb(1, Ahat, 1, N.Data, -1, NT.Data);
                    }
                    else if (P.get<std::string>("SurfNavStokes.nonlinear_term") == "strain")
                    {
                        std::cout << "nonlinear_term: strain " << std::endl;

                        Adyn.LinComb(1, Ahat, 1, N.Data, 1, NT.Data);
                    }
                    else if (P.get<std::string>("SurfNavStokes.nonlinear_term") == "skew-symmetric")
                    {
                        std::cout << "nonlinear_term: skew-symmetric " << std::endl;

                        Adyn.LinComb(1, Ahat, 1, N.Data, 0.5, D.Data);
                    }
                    else if (P.get<std::string>("SurfNavStokes.nonlinear_term") == "conservative")
                    {
                        std::cout << "nonlinear_term: conservative " << std::endl;

                        Adyn.LinComb(1, Ahat, 1, N.Data, 1, D.Data);
                    }
                    else if (P.get<std::string>("SurfNavStokes.nonlinear_term") == "EMAC")
                    {
                        std::cout << "nonlinear_term: EMAC " << std::endl;

                        Adyn.LinComb(1, Ahat, 1, N.Data, 1, NT.Data, 1, D.Data);
                    }
                    else
                    {
                        std::cout << "nonlinear_term: none" << std::endl;
                        Adyn=Ahat;
                    }
                    //set up rhs from external forces
                    InitVector(mg, fRHS,  extrhs, (i+1)*tau);
                    InitScalar(mg, gRHS, extrhs2, (i+1)*tau);
                    fRHS.Data = M.Data * fRHS.Data;
                    gRHS.Data = Schur.Data * gRHS.Data;
                    //renormalization to fullfill \int rhs = 0 in discrete sence
                    gRHS.Data   -= ( dot(gRHS.Data,id2) / dot(id2,id2) ) * id2;
                    //set actual external force to instant rhs
                    instantrhs = fRHS;
                    // pick inertial term and reinitialise unknowns
                    if ( P.get<std::string>("SurfNavStokes.instationary") == "BDF1")
                    {
                        std::cout << "inertial term: BDF1 " << std::endl;
                        instantrhs.Data += ( 1.0/tau ) * ( M.Data * v.Data) ;
                    }
                    else if ( P.get<std::string>("SurfNavStokes.instationary") == "BDF2")
                    {

                        if (i!=0) {
                            std::cout << "inertial term: BDF2 " << std::endl;
                            instantrhs.Data += ( 2.0/tau ) * ( M.Data * v.Data) - ( 1.0/(2.0*tau)) * ( M.Data * v_old.Data) ;
                        } else {
                            std::cout << "inertial term: BDF1 " << std::endl;
                            instantrhs.Data += ( 1.0/tau ) * ( M.Data * v.Data) ;
                        }
                        //v_oldold=v_old;
                        v_old=v;

                    }
                    else
                    {
                        std::cout << "problem is instationary, pick _BDF1_ or _BDF2_ or _none_" << std::endl;
                        return 0;
                    }

                    /*std::cout << "REMOVE THIS!" << std::endl;
                    MatrixCL Adyncopy=Adyn;
                    Adyn.LinComb(1.0, Adyncopy, 1.0, M.Data);*/

                    //solve linear system
                    Solver->Solve(Adyn, Bhat, Chat, v.Data, p.Data, instantrhs.Data, gRHS.Data,v.RowIdx->GetEx(), p.RowIdx->GetEx() );

                    //std::cout <<  (3./(2.*tau)) * v.Data  - (2./tau) * v_old.Data + (1./(2.*tau)) * v_oldold.Data;

                    //rhs.Data    -= ( dot(rhs.Data,id)   / dot(id,id)   ) * id;
                    std::cout << "norm(Omega * v ) ist: " << norm(Omega.Data*v.Data) << std::endl;
                    std::cout << "norm(M * v ) ist: " << norm(M.Data*v.Data) << std::endl;

                    std::cout << "dot(Omega*v, M*id2)" << dot(Omega.Data*id,Schur.Data*id2) << std::endl;

                    Mhat.LinComb(1., Schur.Data, rho_p, Schur_normal_stab.Data);
                    //renormalization to fullfill \int rhs = 0 in discrete sence

                    curl_proj.Data = Omega.Data*v.Data  - ( dot(Omega.Data*v.Data,id2) / dot(id2,id2) ) * id2;
                    //GMResSolver_.Solve(Mhat, curl.Data,curl_proj.Data , gRHS.Data);
                    //curl.Data = gRHS.Data;

                    //postprocess output to normalize pressure
                    p.Data -=  dot(Schur.Data * p.Data,id2) / dot(Schur.Data*id2,id2) * id2;

                    //check residuals
                    std::cout << "norm( Adyn * v + B^T  * p  -  instantrhs) ist: " << norm(Adyn*v.Data + BTranspose*p.Data - instantrhs.Data) << std::endl;
                    std::cout << "norm( B    * v + Chat * p - gRHS) ist: " << norm(B.Data*v.Data + Chat * p.Data - gRHS.Data) << std::endl;

                    //skip vtk output of needed
                    if ((i+1) % P.get<int>("Output.every timestep")  == 0)
                    {
                        std::cout << "output # : " << i << std::endl;
                        vtkwriter->Write((i+1)*tau);
                    }

                    //send to logfiles
                    log_error <<  std::to_string((float)(i+1)*tau) << "\t";
                    log_solo <<  std::to_string((float)(i+1)*tau) << "\t";
                    log << "Time is: "    << std::to_string((float)(i+1)*tau) << std::endl;
                    log << "h is: "       << h << std::endl;
                    log << "rho_p is: "   << rho_p/**pow(h, -rho_p_order) )*/<< std::endl;
                    log << "rho_u is: " << rho_u/**pow(h, -rho_u_order) )*/<< std::endl;
                    log << "tau_u is: "     << tau_u/**pow(h, -tau_u_order) )<< std::endl*/;
                    log << "Total iterations: " << nonsym_stokessolver->GetIter() << '\n';
                    log	<< "Final residual: " << nonsym_stokessolver->GetResid() << '\n';
                    if( !velFE.compare("P2")) {
                        vxtent.SetIdx( &vecP2idx);
                        Extend(mg, v, vxtent);
                        BndDataCL<Point3DCL> bndvec = vbnd;
                        log << "The L2-Norm of v - vSol is: " << L2_Vector_error(mg, lset.Phi, lset.GetBndData(), make_P2Eval(mg, bndvec, vxtent), extvsol) << std::endl;
                        double H1( 0.), surfH1( 0.), advanced_surfH1( 0.), normal_velocity;
                        H1_Vector_error_P2(mg, lset.Phi, lset.GetBndData(), v, vbnd, extvsol, extsol_grad1, extsol_grad2, extsol_grad3, tau_u, H1, surfH1, advanced_surfH1, normal_velocity);
                        //        std::cout << "The H1-Norm of v - vSol is: " << H1 << std::endl;
                        //        std::cout << "The surfH1-Norm of v - vSol is: " << surfH1 << std::endl;
                        log << "The advanced surfH1-Norm of v - vSol is: " << advanced_surfH1 << std::endl;
                        log << "The U-Norm of v - vSol is: " << std::sqrt(advanced_surfH1*advanced_surfH1 + rho_u*dot(A_stab.Data*v.Data, v.Data)) << std::endl;
                        log << "The L2-Norm of v * n is: " << normal_velocity << std::endl;
                        if ( !prFE.compare("P1")) {
                            pxtent.SetIdx( &P1FEidx);
                            Extend(mg, p, pxtent);
                            BndDataCL<double> bndscalar = pbnd;
                            double L2_Lagrange =L2_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, bndscalar, pxtent), extpsol);
                            log << "The L2-Norm of p - pSol is: " << L2_Lagrange << std::endl;
                            log << "The M-Norm of p - pSol is: " << std::sqrt(L2_Lagrange*L2_Lagrange + hat_rho_u*dot(Schur_stab.Data*p.Data, p.Data)) << std::endl;
                        }
                    } else if( !velFE.compare("P1")) {
                        vxtent.SetIdx( &vecP1idx);
                        Extend(mg, v, vxtent);
                        BndDataCL<Point3DCL> bndvec = vbnd;
                        double error_L2 = L2_Vector_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, bndvec, vxtent), extvsol, (i+1)*tau);
                        log << "The L2-Norm of v - vSol is: " << error_L2 << std::endl;
                        log_error <<  std::to_string((float)error_L2) << "\t";
                        double l2norm = L2_Vector_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, bndvec, vxtent), &ZeroVectorFun, (i+1)*tau);
                        log << "The L2-Norm of v  is: " << l2norm<< std::endl;
                        log_solo << std::to_string((float)(l2norm*l2norm*0.5)) << "\t" << std::to_string((float) dot(M.Data * v.Data,id) ) << "\t" << dot(M.Data * v.Data,fRHS.Data) << std::endl;
                        double H1( 0.), surfH1( 0.), advanced_surfH1( 0.), normal_velocity( 0.);
                        H1_Vector_error_P1(mg, lset.Phi, lset.GetBndData(), v, vbnd, extvsol, extsol_grad1, extsol_grad2, extsol_grad3, tau_u, H1, surfH1, advanced_surfH1, normal_velocity, (i+1)*tau);
                        //        std::cout << "The H1-Norm of v - vSol is: " << H1 << std::endl;
                        //        std::cout << "The surfH1-Norm of v - vSol is: " << surfH1 << std::endl;
                        log << "The advanced surfH1-Norm of v - vSol is: " << advanced_surfH1 << std::endl;
                        log_error <<  std::to_string((float)advanced_surfH1) << "\t";
                        log_error <<  std::to_string((float)surfH1) << "\t";
                        log_error <<  std::to_string((float)H1) << "\t";
                        //log << "The U-Norm of v - vSol is: " << std::sqrt(advanced_surfH1*advanced_surfH1 + rho_u*dot(A_stab.Data*v.Data, v.Data)) << std::endl;
                        log << "The L2-Norm of v * n is: " << normal_velocity << std::endl;
                        log_error <<  std::to_string((float)normal_velocity) <<  "\t";

                        if ( !prFE.compare("P1")) {
                            pxtent.SetIdx( &P1FEidx);
                            Extend(mg, p, pxtent);
                            BndDataCL<double> bndscalar = pbnd;
                            double L2_Lagrange =L2_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, bndscalar, pxtent), extpsol, (i+1)*tau);
                            log << "The L2-Norm of p - pSol is: " << L2_Lagrange << std::endl;
                            log_error <<  std::to_string((float)L2_Lagrange) << std::endl;
                            log << "The L2-Norm of p  is: " << L2_error(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, bndscalar, pxtent), &ZeroScalarFun);
                            //log << "The M-Norm of p - pSol is: " << std::sqrt(L2_Lagrange*L2_Lagrange + rho_p*dot(Schur_stab.Data*p.Data, p.Data)) << std::endl;
                        }
                    }

                    log.close();
                }//eng of time iteration

            }
            else if ( P.get<std::string>("SurfNavStokes.instationary") == "none" ) {
                //renormalization to fullfill \int rhs = 0 in discrete sence
                gRHS.Data -= (dot(gRHS.Data, id2) / dot(id2, id2)) * id2;
                std::ofstream log(
                     dirname + "/" + filename +
                    "time=1_m=" + std::to_string(param.input.numbOfVirtualSubEdges) +
                    "_" + FE +
                    "_patchnormals=" + std::to_string(param.input.usePatchNormal) +
                    "_shift=" + P.get<std::string>("Levelset.ShiftNorm", "0") +
                    "_form=" + formulation +
                    "_rhofac=" + P.get<std::string>("SurfNavStokes.vel_volumestab_fac") +
                    "_rhopow=" + P.get<std::string>("SurfNavStokes.vel_volumestab_pow") + ".txt"
                );

                if (  ( P.get<std::string>("SurfNavStokes.nonlinear_term") == "convective" )
                   || ( P.get<std::string>("SurfNavStokes.nonlinear_term") == "rotational" )
                   || ( P.get<std::string>("SurfNavStokes.nonlinear_term") == "strain" )
                   || ( P.get<std::string>("SurfNavStokes.nonlinear_term") == "skew-symmetric" )
                   || ( P.get<std::string>("SurfNavStokes.nonlinear_term") == "conservative" )
                   || ( P.get<std::string>("SurfNavStokes.nonlinear_term") == "EMAC" )
                   )
                {
                    log << "nonlinear_term term is: " << P.get<std::string>("SurfNavStokes.nonlinear_term") << std::endl;
                    //nonlinear procedure parameters
                    double omega=1.0;
                    double max_iter=30;
                    int iter = 0;

                    std::cout<<"REMOVE IT!" << std::endl;
                    //v=vInit;

                    //fixed point iterations
                    double res;
                    do
                    {
                        iter++;
                        v_old = v;
                        p_old = p;

                        //setup matrices based on v_old
                        SetupNavierStokesIF_P1P1(mg, &A, &A_stab, &B, &Omega, &N, &NT, &M, &D, &S, &L, &L_stab, &Schur, &Schur_stab, &Schur_normal_stab, lset, v_old, vbnd, &param);

                        Ahat.LinComb(mu, A.Data, 1.0, M.Data, tau_u, S.Data, rho_u, A_stab.Data);
                        Bhat.LinComb(1., B.Data, 0.,B.Data);
                        transpose(B.Data, BTranspose);
                        if (P.get<std::string>("SurfNavStokes.stab") == "full")
                            Chat.LinComb(0, Schur.Data, -rho_p, Schur_stab.Data);
                        else if (P.get<std::string>("SurfNavStokes.stab") == "normal")
                            Chat.LinComb(0, Schur.Data, -rho_p, Schur_normal_stab.Data);
                        else {
                            std::cout << "WRN: no stabilization is used!\n";
                            Chat.LinComb(0, Schur.Data, 0, Schur.Data);
                        }
                        if ( P.get<std::string>("SurfNavStokes.nonlinear_term") == "convective" )
                        {
                            std::cout << "nonlinear_term: convective " << std::endl;
                            Adyn.LinComb(1, Ahat, 1, N.Data);

                        }
                        else if (P.get<std::string>("SurfNavStokes.nonlinear_term") == "rotational")
                        {
                            std::cout << "nonlinear_term: rotational " << std::endl;
                            //transpose(N.Data, NTranspose);
                            Adyn.LinComb(1, Ahat, 1, N.Data, -1, NT.Data);
                        }
                        else if (P.get<std::string>("SurfNavStokes.nonlinear_term") == "strain")
                        {
                            std::cout << "nonlinear_term: strain " << std::endl;
                            //transpose(N.Data, NTranspose);
                            Adyn.LinComb(1, Ahat, 1, N.Data, 1, NT.Data);
                        }else if (P.get<std::string>("SurfNavStokes.nonlinear_term") == "skew-symmetric")
                        {
                            std::cout << "nonlinear_term: skew-symmetric " << std::endl;

                            Adyn.LinComb(1, Ahat, 1, N.Data, 0.5, D.Data);
                        }
                        else if (P.get<std::string>("SurfNavStokes.nonlinear_term") == "conservative")
                        {
                            std::cout << "nonlinear_term: conservative " << std::endl;

                            Adyn.LinComb(1, Ahat, 1, N.Data, 1, D.Data);
                        }
                        else if (P.get<std::string>("SurfNavStokes.nonlinear_term") == "EMAC")
                        {
                            std::cout << "nonlinear_term: EMAC " << std::endl;

                            Adyn.LinComb(1, Ahat, 1, N.Data, 1, NT.Data, 1, D.Data);
                        }


                        //solve for v_aux, p_aux
                        nonsym_stokessolver->Solve(Adyn, Bhat, Chat, v_aux.Data, p_aux.Data, fRHS.Data, gRHS.Data,v.RowIdx->GetEx(), p.RowIdx->GetEx() );

                        p_aux.Data -=  dot(Schur.Data * p_aux.Data,id2) / dot(Schur.Data*id2,id2) * id2;

                        //preconditioning
                        v.Data = omega*v_aux.Data + (1-omega)*v_old.Data;
                        p.Data = omega*p_aux.Data + (1-omega)*p_old.Data;

                        std::cout << "norm( Adyn*v + BTranspose*p - rhs ): " << norm(Adyn*v.Data + BTranspose*p.Data - fRHS.Data)<<std::endl;
                        std::cout << "norm( B * v + Chat*p - gRHS       ): " << norm( B.Data*v.Data + Chat*p.Data -gRHS.Data ) << std::endl;

                        res = pow( norm(Adyn*v.Data + BTranspose*p.Data - fRHS.Data), 2.0) + pow( norm(B.Data*v.Data + Chat * p.Data - gRHS.Data), 2.0);
                        std::cout << "residual^2 after nonlinear iteration #" << iter << ": " << res << std::endl;
                    }
                    while (  ( iter < max_iter)
                          && ( res > pow( P.get<double>("Solver.Tol"), 2.0) ) );

                    std::cout << "# nonlinear iterations: " << iter << std::endl;
                }
                else //standart Stokes' case
                {
                    Adyn.LinComb(1., A.Data, 1., M.Data, tau_u, S.Data, rho_u, A_stab.Data);
                    Bhat.LinComb(1., B.Data, 0., B.Data);
                    if (P.get<std::string>("SurfNavStokes.stab") == "full")
                        Chat.LinComb(0, Schur.Data, -rho_p, Schur_stab.Data);
                    else if (P.get<std::string>("SurfNavStokes.stab") == "normal")
                        Chat.LinComb(0, Schur.Data, -rho_p, Schur_normal_stab.Data);
                    else {
                        std::cout << "WRN: no stabilization is used!\n";
                        Chat.LinComb(0, Schur.Data, 0, Schur.Data);
                    }
                    transpose(B.Data, BTranspose);
                    stokessolver->Solve(Adyn, Bhat, Chat, v.Data, p.Data, fRHS.Data, gRHS.Data, v.RowIdx->GetEx(), p.RowIdx->GetEx());
                }
                if (P.get<int>("Output.every timestep") > 0) {
                    std::cout << "writing vtk...\n";
                    vtkwriter->Write(1); // skip vtk output if needed
                }
                else std::cout << "skipping vtk output\n";
                auto velResSq = norm_sq(Adyn * v.Data + BTranspose * p.Data - fRHS.Data);
                auto preResSq = norm_sq(B.Data * v.Data + Chat * p.Data - gRHS.Data);
                auto residual = sqrt(velResSq + preResSq);
                std::cout << "residual: " << residual << '\n';
                // postprocess output pressure up to constant
                p.Data -=  dot(Schur.Data * p.Data, id2) / dot(Schur.Data*id2, id2) * id2;
                log << "norm(Adyn * v + B^T  * p -  rhs): " << sqrt(velResSq) << '\n'
                    << "norm(B    * v + Chat * p - gRHS): " << sqrt(preResSq) << '\n';
                // send to logfiles
                log_solo <<  std::to_string((float)1) << "\t";
                log << "Time is: " << std::to_string((float)1) << std::endl;
                log << "h is: " << h << std::endl;
                log << "rho_p is: "   << rho_p << std::endl;
                log << "rho_u is: " << rho_u << std::endl;
                log << "tau_u is: "     << tau_u << std::endl;
                log << "Total iterations: " << Solver->GetIter() << '\n';
                log	<< "Final MINRES residual: " << Solver->GetResid() << '\n';
                VectorCL vSolMinusV = vSol.Data - v.Data, pSolMinusP = pSol.Data - p.Data;
                auto velL2          = sqrt(dot(v.Data, M.Data * v.Data));
                auto velL2err       = dot(vSolMinusV, M.Data * vSolMinusV);
                auto velNormalL2    = dot(v.Data, S.Data * v.Data);
                auto velTangenL2    = sqrt(velL2err - velNormalL2);
                velL2err = sqrt(velL2err);
                velNormalL2 = sqrt(velNormalL2);
                auto velH1err       = sqrt(dot(vSolMinusV, A.Data * vSolMinusV));
                auto preL2          = sqrt(dot(p.Data, Schur.Data * p.Data));
                auto preL2err       = sqrt(dot(pSolMinusP, Schur.Data * pSolMinusP));
                log << "The L2-Norm of v - vSol is: " << velL2err << '\n';
                log << "The L2-Norm of v  is: " << velL2 << '\n';
                log_solo << std::to_string((float)(velL2*velL2*0.5)) << '\n';
                log << "The H1-Norm of v - vSol is: " << velH1err << '\n';
                log << "The L2-Norm of v * n is: " << velNormalL2 << '\n';
                log << "The L2-Norm of p - pSol is: " << preL2err << '\n';
                log << "The L2-Norm of p is: " << preL2 << '\n';
                log << "The L2-Norm of v_T - vSol is: " << velTangenL2 << '\n';
                log << "Actual residual is: " << residual << '\n';
            }
            log_solo.close();
            log_error.close();
            std::cout << "Output is located: " << dirname << std::endl;
        }
        // error
        double Aaverage( 0.), Avariation( 0.), Schuraverage( 0.), Schurvariation( 0.);
        double PCaverage( 0.), PCvariation( 0.);
        std::stringstream Astreamcopy( Astream.str());
        std::stringstream PCstreamcopy( PCstream.str());
        std::stringstream Schurstreamcopy( Schurstream.str());
        RightComputeAverageIterations(Astream, Aaverage);
        ComputeAverageIterations(Schurstream, Schuraverage);
        RightComputeAverageIterations(PCstream, PCaverage);
        RightComputeVariationFromAverageIterations(Astreamcopy, Aaverage, Avariation);
        RightComputeVariationFromAverageIterations(PCstreamcopy, PCaverage, PCvariation);
        ComputeVariationFromAverageIterations(Schurstreamcopy, Schuraverage, Schurvariation);
        std::cout << "The average iterationsnumber of the A-preconditioner is: " <<     Aaverage     << '\n' << " ...with a variation of: " << Avariation << std::endl;
        std::cout << "The average iterationsnumber of the nonsymmetric A-preconditioner is: " <<     PCaverage     << '\n' << " ...with a variation of: " << PCvariation << std::endl;
        std::cout << "The average iterationsnumber of the Schur-preconditioner is: " << Schuraverage << '\n' << " ...with a variation of: " << Schurvariation << std::endl;
        delete &lset;
        return 0;
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
}