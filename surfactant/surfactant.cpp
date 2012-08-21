/// \file
/// \brief Solve a non-stationary convection-diffusion-equation on a moving interface
/// \author LNM RWTH Aachen: Joerg Grande

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
#include "out/output.h"
#include "out/vtkOut.h"
#include "misc/funcmap.h"

#include <fstream>
#include <string>

using namespace DROPS;

DROPS::ParamCL P;

std::auto_ptr<VTKOutCL> vtkwriter;

DROPS::InVecMap& invecmap= DROPS::InVecMap::getInstance();
DROPS::InScaMap& inscamap= DROPS::InScaMap::getInstance();

instat_vector_fun_ptr the_wind_fun;
instat_scalar_fun_ptr the_lset_fun;
instat_scalar_fun_ptr the_rhs_fun;
instat_scalar_fun_ptr the_sol_fun;

typedef DROPS::Point3DCL (*bnd_val_fun) (const DROPS::Point3DCL&, double);

DROPS::BndCondT bc_wind[6]= { DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC };
bnd_val_fun bf_wind[6];

instat_scalar_fun_ptr sigma( 0);
SurfaceTensionCL sf( sigma, 0);
DROPS::LsetBndDataCL lsbnd( 6);


DROPS::MGBuilderCL* make_MGBuilder (const DROPS::ParamCL& P)
{
    const std::string type= P.get<std::string>( "Type");
    if (type != std::string("BrickBuilder"))
        throw DROPS::DROPSErrCL(std::string( "make_MGBuilder: Unknown Domain: ") + type + std::string("\n"));

    const Point3DCL orig( P.get<Point3DCL>( "Origin")),
                      e1( P.get<Point3DCL>( "E1")),
                      e2( P.get<Point3DCL>( "E2")),
                      e3( P.get<Point3DCL>( "E3"));
    const Uint n1= P.get<Uint>( "N1"),
               n2= P.get<Uint>( "N2"),
               n3= P.get<Uint>( "N3");
    return new BrickBuilderCL( orig, e1, e2, e3, n1, n2, n3);
}

// Surface divergence of a vector field w
inline double div_gamma_wind (const Point3DCL& n, const SMatrixCL<3,3>& dn,
                              const Point3DCL& w, const SMatrixCL<3,3>& dw)
{
    const double tr_Pdw= trace( dw) - inner_prod( n, dw*n),
                 tr_dn_nw= trace( dn)*inner_prod( n, w),
                 wT_dn_n= inner_prod( w, dn*n);
    return tr_Pdw - tr_dn_nw - wT_dn_n;
}

// laplace-beltrami of a function u
inline double laplace_beltrami_u (const Point3DCL& n,      const SMatrixCL<3,3>& dn,
                                  const Point3DCL& grad_u, const SMatrixCL<3,3>& Hess_u)
{
    return div_gamma_wind( n, dn, grad_u, Hess_u);
//     const double tr_PHessu= trace( Hess_u) - inner_prod( n, Hess_u*n),
//                  tr_dn_ngradu= trace( dn)*inner_prod( n, grad_u),
//                  graduT_dn_n= iner_prod( grad_u, dn*n);
//     return tr_PHessu - tr_dn_ngradu - graduT_dn_n;
}


DROPS::Point3DCL WindVelocity;
DROPS::Point3DCL constant_wind (const DROPS::Point3DCL&, double)
{
    return WindVelocity;
}
static RegisterVectorFunction regvec_constant_wind( "ConstantWind", constant_wind);

DROPS::Point3DCL RadDrop;
DROPS::Point3DCL PosDrop;
double ellipsoid (const DROPS::Point3DCL& p, double)
{
    const DROPS::Point3DCL x= (p - PosDrop)/RadDrop;
    return x.norm_sq() - 1.;
}
static RegisterScalarFunction regsca_ellipsoid_lset( "Ellipsoid", ellipsoid);


DROPS::Point2DCL RadTorus; // R= RadTorus[0], r= RadTorus[1]; R > r for tori with a hole.
double torus (const Point3DCL& p, double)
{
    return std::sqrt( std::pow( RadTorus[0] - std::sqrt(p[0]*p[0] + p[1]*p[1]), 2) + std::pow( p[2], 2)) - RadTorus[1];
}
static RegisterScalarFunction regsca_torus_lset( "Torus", torus);

// ==non-stationary test case "HeatConduction"
// "ConstantWind" 0
// "Ellipsoid" pos:[0,0,0], rad:[1,1,1]

double heat_conduction_u0 (const Point3DCL& p)
{
    return p[0]*p[1];
}

double heat_conduction_rhs (const Point3DCL&, double)
{
    return 0.;
}
static RegisterScalarFunction regsca_heat_conduction_rhs( "HeatConductionRhs", heat_conduction_rhs);

double heat_conduction_sol (const Point3DCL& p, double t)
{
    return std::exp( -6.*t)*heat_conduction_u0( p);
}
static RegisterScalarFunction regsca_heat_conduction_sol( "HeatConductionSol", heat_conduction_sol);


// ==non-stationary test-case "ToroidalFlow"==
// Level set: "torus" with RadTorus

double angular_velocity (const Point3DCL& p, double)
{
    return 1. + p[2];
}

DROPS::SMatrixCL<3,3> rotation_matrix (const DROPS::Point3DCL& p, double t)
{
    const double omega= angular_velocity( p, 0.);
    SMatrixCL<3,3> m;
    m(0,0)= std::cos(omega*t); m(0,1)= -std::sin(omega*t);
    m(1,0)= std::sin(omega*t); m(1,1)=  std::cos(omega*t);
    m(2,2)= 1.;
    return m;
}

DROPS::Point3DCL toroidal_flow (const DROPS::Point3DCL& p, double t)
{
    return rotation_matrix( p, t)*p;
}

DROPS::Point3DCL toroidal_flow_wind (const DROPS::Point3DCL& p, double)
{
    const double omega= angular_velocity( p, 0.);
    return MakePoint3D( -p[1]*omega, p[0]*omega, 0.);
}
static RegisterVectorFunction regvec_toroidal_flow_wind( "ToroidalFlowWind", toroidal_flow_wind);

double toroidal_flow_sol (const Point3DCL& p, double t)
{
    const Point3DCL q( toroidal_flow( p, -t));
    return std::exp( -t)*q[0]*q[1];
}
static RegisterScalarFunction regsca_toroidal_flow_sol( "ToroidalFlowSol", toroidal_flow_sol);

double toroidal_flow_rhs (const Point3DCL& p, double t)
{
    const Point3DCL w( toroidal_flow_wind( p, t));
    const double omega= angular_velocity( p, 0.);
    SMatrixCL<3,3> dw;
    dw(0,1)= -omega; dw(0,2)= -p[1];
    dw(1,0)=  omega; dw(1,2)=  p[0];

    const Point2DCL xhat( MakePoint2D( p[0], p[1]));
    const double norm_xhat= xhat.norm();
    const double l= torus( p, t) + RadTorus[1];
    const Point2DCL tt= (norm_xhat - RadTorus[0])/(l*norm_xhat)*xhat;
    const Point3DCL n( MakePoint3D( tt[0], tt[1], p[2]/l));
    SMatrixCL<3,3> dn;
    dn= eye<3,3>() - outer_product( n, n);
    SMatrixCL<2,2> dnhat= RadTorus[0]/norm_xhat*(eye<2,2>() - outer_product( xhat/norm_xhat, xhat/norm_xhat));
    dn(0,0)-= dnhat(0,0); dn(0,1)-= dnhat(0,1);
    dn(1,0)-= dnhat(1,0); dn(1,1)-= dnhat(1,1);
    dn*= 1./l;

    const double c= std::cos( omega*t),
                 s= std::sin( omega*t);
    const Point3DCL z( toroidal_flow( p, -t));
    const Point3DCL dz0( MakePoint3D(  c, s, t*(-s*p[0] + c*p[1]))),
                    dz1( MakePoint3D( -s, c, t*(-c*p[0] - s*p[1])));
    const Point3DCL grad_u= std::exp( -t)*(z[1]*dz0 + z[0]*dz1);
    SMatrixCL<3,3> Hess_u;
    Hess_u(0,2)= -c*z[0] - s*z[1];
    Hess_u(1,2)= -s*z[0] + c*z[1];
    Hess_u(2,0)= -c*z[0] - s*z[1];
    Hess_u(2,1)= -s*z[0] + c*z[1];
    Hess_u(2,2)= t*(z[0]*(s*p[0] - c*p[1]) + z[1]*(-c*p[0] - s*p[1]));
    Hess_u*= t;
    Hess_u+= outer_product( dz0, dz1) + outer_product( dz1, dz0);
    Hess_u*= std::exp( -t);

    const double u= toroidal_flow_sol( p, t),
                 mat_der=  -u,
                 reaction= div_gamma_wind( n, dn, w, dw)*u,
                 diffusion= -laplace_beltrami_u( n, dn, grad_u, Hess_u);

    return mat_der + reaction + diffusion;
}
static RegisterScalarFunction regsca_toroidal_flow_rhs( "ToroidalFlowRhs", toroidal_flow_rhs);


// ==non-stationary test case "AxisScaling"==

double axis_scaling (double t)
{
    return 1. + 0.25*std::sin( t);
}

DROPS::Point3DCL axis_scaling_wind (const DROPS::Point3DCL& p, double t)
{
    return MakePoint3D( 0.25*std::cos( t)/axis_scaling( t)*p[0], 0., 0.);
}
static RegisterVectorFunction regvec_axis_scaling_wind( "AxisScalingWind", axis_scaling_wind);

double axis_scaling_lset (const Point3DCL& p, double t)
{
    return std::pow( p[0]/axis_scaling( t), 2) + std::pow( p[1], 2) + std::pow( p[2], 2) - 1.;
}
static RegisterScalarFunction regsca_axis_scaling_lset( "AxisScalingLset", axis_scaling_lset);

double axis_scaling_sol (const Point3DCL& p, double t)
{
    return std::exp( -0.5*t)*heat_conduction_u0( p); // The same as for heat conduction test case.
}
static RegisterScalarFunction regsca_axis_scaling_sol( "AxisScalingSol", axis_scaling_sol);

double axis_scaling_rhs (const Point3DCL& p, double t)
{
    const double a= axis_scaling( t);
    const double bf4= (p/MakePoint3D( std::pow( a, 2), 1., 1.)).norm_sq();
    const double bf6= (p/MakePoint3D( std::pow( a, 3), 1., 1.)).norm_sq();

    const double mat_der=  0.25*std::cos( t)/a - 0.5;
    const double reaction= 0.25*std::cos( t)/a/bf4*( std::pow( p[1], 2) + std::pow( p[2], 2)
        - 2.*std::pow( p[0]/a, 2)*(1. + 1./std::pow( a, 2) - bf6/bf4));
    const double diffusion= (1. + 1./std::pow( a, 2))/bf4*(-3. - 2./std::pow( a, 2) + 2.*bf6/bf4);

//     const Point3DCL tt( p/MakePoint3D( std::pow( a, 2), 1., 1.));
//     const double l= tt.norm();
//     const Point3DCL n( tt/l);
//     SMatrixCL<3,3> dn( (eye<3,3>() - outer_product( n, n))/l );
//     dn(0,0)/= std::pow( a, 2); dn(1,0)/= std::pow( a, 2); dn(2,0)/= std::pow( a, 2);
// 
//     const Point3DCL w( MakePoint3D( 0.25*p[0]*std::cos( t)/a, 0., 0.));
//     SMatrixCL<3,3> dw;
//     dw(0,0)= 0.25*std::cos( t)/a;
// 
//     const Point3DCL grad_u( std::exp( -0.5*t)*MakePoint3D( p[1], p[0], 0.));
//     SMatrixCL<3,3> Hess_u;
//     Hess_u(0,1)= 1.;
//     Hess_u(1,0)= 1.;
//     Hess_u*= std::exp( -0.5*t);
// 
//     const double err= div_gamma_wind( n, dn, w, dw) - reaction,
//                errlb= laplace_beltrami_u( n, dn, grad_u, Hess_u) - diffusion*axis_scaling_sol( p, t);
//     if (std::fabs( err) > 1e-12 || std::fabs( errlb) > 1e-12) {
//         std::cerr << err   << " " << div_gamma_wind( n, dn, w, dw) << " " << reaction << "\n"
//                   << errlb << " " << laplace_beltrami_u( n, dn, grad_u, Hess_u) << " " << diffusion*axis_scaling_sol( p, t) << "\n";
//         exit( 1);
//     }

    return (mat_der + reaction - diffusion)*axis_scaling_sol( p, t);
}
static RegisterScalarFunction regsca_axis_scaling_rhs( "AxisScalingRhs", axis_scaling_rhs);


double sphere_dist (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL x( p - PosDrop);
    return x.norm() - RadDrop[0];
}
static RegisterScalarFunction regsca_sphere_dist_lset( "SphereDist", sphere_dist);


typedef double (*dist_funT) (const DROPS::Point3DCL&, double);

double sphere_2move (const DROPS::Point3DCL& p, double t)
{
    DROPS::Point3DCL x( p - (PosDrop + t*constant_wind(p, t)));
    return x.norm() - RadDrop[0];
}

// ==stationary test case "LaplaceBeltrami0"==
// Sphere around 0, RadDrop 1, wind == 0
// A right hand side from C.J. Heine...
const double a( -13./8.*std::sqrt( 35./M_PI));
double laplace_beltrami_0_rhs (const DROPS::Point3DCL& p, double)
{
    return a*(3.*p[0]*p[0]*p[1] - p[1]*p[1]*p[1]);
}
static RegisterScalarFunction regsca_laplace_beltrami_0_rhs( "LaplaceBeltrami0Rhs", laplace_beltrami_0_rhs);

// ...and the corresponding solution (extended)
double laplace_beltrami_0_sol (const DROPS::Point3DCL& p, double)
{
    return p.norm_sq()/(12. + p.norm_sq())*laplace_beltrami_0_rhs( p, 0.);
//    return 1. + p.norm_sq()/(12. + p.norm_sq())*laplace_beltrami_0_rhs( p, 0.);
}
static RegisterScalarFunction regsca_laplace_beltrami_0_sol( "LaplaceBeltrami0Sol", laplace_beltrami_0_sol);


double sol0t (const DROPS::Point3DCL& p, double t)
{
    const Point3DCL q( p - (PosDrop + t*constant_wind(p, t)));
    const double val( a*(3.*q[0]*q[0]*q[1] - q[1]*q[1]*q[1]));

//    return q.norm_sq()/(12. + q.norm_sq())*val;
    return 1. + q.norm_sq()/(12. + q.norm_sq())*val;
}


template<class DiscP1FunType>
double L2_error (const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
    const DiscP1FunType& discsol, DROPS::instat_scalar_fun_ptr extsol)
{
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( 2);
    const double t= discsol.GetTime();
    QuadDomain2DCL qdom;
    std::valarray<double> qsol,
                          qdiscsol;

    double d( 0.);
    DROPS_FOR_TRIANG_CONST_TETRA( discsol.GetMG(), ls.GetLevel(), it) {
        make_CompositeQuad5Domain2D (qdom, *it, lat, ls, lsbnd);
        resize_and_evaluate_on_vertexes ( discsol, *it, qdom, qdiscsol);
        resize_and_evaluate_on_vertexes ( extsol,  *it, qdom, t, qsol);
        d+= quad_2D( std::pow( qdiscsol - qsol, 2), qdom);
    }
    return std::sqrt( d);
}

/// The nodal interpolant of extsol on the interface-FE-space ist computed first.
/// The H1-error is then computed between the interpolant and the numerical solution.
template<class DiscP1FunType>
double H1_error (const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
    const DiscP1FunType& discsol, DROPS::instat_scalar_fun_ptr extsol)
{
    IdxDescCL* idx= const_cast<IdxDescCL*>( discsol.GetSolution()->RowIdx);
    MatDescCL A( idx, idx);
    SetupLBP1( discsol.GetMG(), &A, ls, lsbnd, 1.);
    VecDescCL sol_vec( idx);
    P1Init (extsol, sol_vec, discsol.GetMG(), discsol.GetTime());
    // sol_vec.t= discsol.GetTime();
    // SetupInterfaceRhsP1( discsol.GetMG(), &sol_vec, ls, lsbnd, extsol);
    const VectorCL diff( discsol.GetSolution()->Data - sol_vec.Data);
    return std::sqrt( dot(A.Data*diff, diff));
}

double L2_norm (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
    DROPS::instat_scalar_fun_ptr extsol)
{
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( 2);
    const double t= ls.t;
    QuadDomain2DCL qdom;
    std::valarray<double> qsol;

    double d( 0.);
    DROPS_FOR_TRIANG_CONST_TETRA( mg, ls.GetLevel(), it) {
        make_CompositeQuad5Domain2D (qdom, *it, lat, ls, lsbnd);
        resize_and_evaluate_on_vertexes ( extsol,  *it, qdom, t, qsol);
        d+= quad_2D( qsol*qsol, qdom);
    }
    return std::sqrt( d);
}

void LinearLSInit (const DROPS::MultiGridCL& mg, DROPS::VecDescCL& ls, DROPS::instat_scalar_fun_ptr d, double t= 0.)
{
    const DROPS::Uint lvl= ls.GetLevel(),
                      idx= ls.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= d( it->GetCoord(), t);

    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= ls.Data[it->Unknowns( idx)]=
            0.5*(ls.Data[it->GetVertex( 0)->Unknowns( idx)] + ls.Data[it->GetVertex( 1)->Unknowns( idx)]);
    ls.t= t;
}

void LSInit (const DROPS::MultiGridCL& mg, DROPS::VecDescCL& ls, dist_funT d, double t= 0.)
{
    const DROPS::Uint lvl= ls.GetLevel(),
                      idx= ls.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= d( it->GetCoord(), t);

    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= d( 0.5*(it->GetVertex( 0)->GetCoord() + it->GetVertex( 1)->GetCoord()), t);
    ls.t= t;
}

void InitVel ( const MultiGridCL& mg, VecDescCL* vec, BndDataCL<Point3DCL>& Bnd, instat_vector_fun_ptr LsgVel, double t= 0.)
{
    VectorCL& lsgvel= vec->Data;
    const Uint lvl  = vec->GetLevel(),
               vidx = vec->RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, sit) {
        if (!Bnd.IsOnDirBnd( *sit))
            DoFHelperCL<Point3DCL, VectorCL>::set( lsgvel, sit->Unknowns( vidx),
                LsgVel(sit->GetCoord(), t));
    }
    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, sit) {
        if (!Bnd.IsOnDirBnd( *sit))
            DoFHelperCL<Point3DCL, VectorCL>::set( lsgvel, sit->Unknowns( vidx),
                LsgVel( (sit->GetVertex(0)->GetCoord() + sit->GetVertex(1)->GetCoord())/2., t));
    }
    vec->t= t;
}

void Strategy (DROPS::MultiGridCL& mg, DROPS::AdapTriangCL& adap, DROPS::LevelsetP2CL& lset)
{
    using namespace DROPS;

    adap.MakeInitialTriang( the_lset_fun);

    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
    // LinearLSInit( mg, lset.Phi, the_lset_fun, 0.);
    LSInit( mg, lset.Phi, the_lset_fun, 0.);

    //DROPS::LevelsetP2CL lset2( mg, lsbnd, sf, P.get<double>("Levelset.Theta"), P.get<double>("Levelset.SD")); // Only for output
    //lset2.idx.CreateNumbering( mg.GetLastLevel(), mg);
    //lset2.Phi.SetIdx( &lset2.idx);
    //LSInit( mg, lset2.Phi, the_lset_fun, 0.);

    const double Vol= lset.GetVolume();
    std::cout << "droplet volume: " << Vol << std::endl;

    BndDataCL<Point3DCL> Bnd_v( 6, bc_wind, bf_wind);
    IdxDescCL vidx( vecP2_FE);
    vidx.CreateNumbering( mg.GetLastLevel(), mg, Bnd_v);
    VecDescCL v( &vidx);
    InitVel( mg, &v, Bnd_v, the_wind_fun, 0.);

    //lset2.SetupSystem( make_P2Eval( mg, Bnd_v, v), P.get<double>("Time.StepSize"));

    SurfactantcGP1CL timedisc( mg,
        P.get<double>("SurfTransp.Theta"), P.get<double>("SurfTransp.Visc"),
        &v, Bnd_v, lset.Phi, lset.GetBndData(),
        P.get<int>("SurfTransp.Iter"), P.get<double>("SurfTransp.Tol"),
        P.get<double>("SurfTransp.OmitBound"));
    timedisc.SetRhs( the_rhs_fun);

    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);
    InterfaceP1RepairCL ic_repair( mg, lset.Phi, lset.GetBndData(), timedisc.ic);
    adap.push_back( &ic_repair);
    //LevelsetRepairCL lset2repair( lset2);
    //adap.push_back( &lset2repair);

    // Init Interface-Sol
    timedisc.idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
    std::cout << "NumUnknowns: " << timedisc.idx.NumUnknowns() << std::endl;
    timedisc.ic.SetIdx( &timedisc.idx);
    timedisc.SetInitialValue( the_sol_fun, 0.);

    ScalarFunAsP2EvalCL the_sol_eval( the_sol_fun, 0., &mg);
    if (vtkwriter.get() != 0) {
        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "Levelset") );
        vtkwriter->Register( make_VTKIfaceScalar( mg, timedisc.ic,  "InterfaceSol"));
        vtkwriter->Register( make_VTKVector(      make_P2Eval( mg, Bnd_v, v),      "Velocity"));
        //vtkwriter->Register( make_VTKScalar(      lset2.GetSolution(),             "Levelset2"));
        vtkwriter->Register( make_VTKScalar( the_sol_eval, "TrueSol"));
        //vtkwriter->Register( make_VTKScalar( ScalarFunAsP2EvalCL( the_sol_fun, 0., &mg), "TrueSol"));
        vtkwriter->Write( 0.);
    }
    //if (P.get<int>( "SolutionOutput.Freq") > 0)
    //    DROPS::WriteFEToFile( timedisc.ic, mg, P.get<std::string>( "SolutionOutput.Path"), P.get<bool>( "SolutionOutput.Binary"));

    const double dt= P.get<double>("Time.StepSize");
    double L_2x_err= L2_error( lset.Phi, lset.GetBndData(), timedisc.GetSolution(), the_sol_fun);
    std::cout << "L_2x-error: " << L_2x_err
              << "\nnorm of true solution: " << L2_norm( mg, lset.Phi, lset.GetBndData(), the_sol_fun)
              << std::endl;
    double L_inftL_2x_err= L_2x_err;
    std::cout << "L_inftL_2x-error: " <<  L_inftL_2x_err << std::endl;
    double H_1x_err= H1_error( lset.Phi, lset.GetBndData(), timedisc.GetSolution(), the_sol_fun);
    std::cout << "H_1x-error: " << H_1x_err << std::endl;
    double L_2tH_1x_err_sq= 0.5*dt*std::pow( H_1x_err, 2);
    BndDataCL<> ifbnd( 0);
    std::cerr << "initial surfactant on \\Gamma: " << Integral_Gamma( mg, lset.Phi, lset.GetBndData(), make_P1Eval(  mg, ifbnd, timedisc.ic)) << '\n';

    for (int step= 1; step <= P.get<int>("Time.NumSteps"); ++step) {
        std::cout << "======================================================== step " << step << ":\n";
        const double cur_time= step*dt;
        // Assumes (as the rest of Drops), that the current triangulation is acceptable to perform the time-step.
        // If dt is large and AdapRef.Width is small, this may not be true.
        // Watch for large differences in numbers of old and new dof.
        timedisc.InitTimeStep();
        LSInit( mg, lset.Phi, the_lset_fun, cur_time);
        InitVel( mg, &v, Bnd_v, the_wind_fun, cur_time);
        timedisc.DoStep( cur_time);
        std::cout << "surfactant on \\Gamma: " << Integral_Gamma( mg, lset.Phi, lset.GetBndData(), make_P1Eval(  mg, ifbnd, timedisc.ic)) << '\n';
        L_2x_err= L2_error( lset.Phi, lset.GetBndData(), timedisc.GetSolution(), the_sol_fun);
        std::cout << "L_2x-error: " << L_2x_err
                  << "\nnorm of true solution: " << L2_norm( mg, lset.Phi, lset.GetBndData(), the_sol_fun)
                  << std::endl;
        L_inftL_2x_err= std::max( L_inftL_2x_err, L_2x_err);
        std::cout << "L_inftL_2x-error: " << L_inftL_2x_err << std::endl;
        L_2tH_1x_err_sq+= (step > 1 ? 0.5 : 0.)*dt*std::pow( H_1x_err, 2);
        H_1x_err= H1_error( lset.Phi, lset.GetBndData(), timedisc.GetSolution(), the_sol_fun);
        std::cout << "H_1x-error: " << H_1x_err << std::endl;
        L_2tH_1x_err_sq+= 0.5*dt*std::pow( H_1x_err, 2);
        std::cout << "L_2tH_1x-error: " << std::sqrt( L_2tH_1x_err_sq) << std::endl;
        if (vtkwriter.get() != 0 && step % P.get<int>( "VTK.VTKOut") == 0) {
            the_sol_eval.SetTime( cur_time);
            vtkwriter->Write( cur_time);
        }
        if (P.get<int>( "SolutionOutput.Freq") > 0 && step % P.get<int>( "SolutionOutput.Freq") == 0)
            DROPS::WriteFEToFile( timedisc.ic, mg, P.get<std::string>( "SolutionOutput.Path"), P.get<bool>( "SolutionOutput.Binary"));

//        lset2.DoStep();
//        VectorCL rhs( lset2.Phi.Data.size());
//        lset2.ComputeRhs( rhs);
//        lset2.SetupSystem( make_P2Eval( mg, Bnd_v, v, cur_time));
//        lset2.SetTimeStep( P.get<double>("Time.StepSize"));
//        lset2.DoStep( rhs);

//         std::cout << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
//         if (P.get("Levelset.VolCorr", 0)) {
//             double dphi= lset.AdjustVolume( Vol, 1e-9);
//             std::cout << "volume correction is " << dphi << std::endl;
//             lset.Phi.Data+= dphi;
//             std::cout << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
//         }
        //if (C.rpm_Freq && step%C.rpm_Freq==0) { // reparam levelset function
            // lset.ReparamFastMarching( C.rpm_Method);
        const bool doGridMod= P.get<int>("AdaptRef.Freq") && step%P.get<int>("AdaptRef.Freq") == 0;
        const bool gridChanged= doGridMod ? adap.UpdateTriang( lset) : false;
        if (gridChanged) {
            std::cout << "Triangulation changed.\n";
            vidx.DeleteNumbering( mg);
            vidx.CreateNumbering( mg.GetLastLevel(), mg, Bnd_v);
            v.SetIdx( &vidx);
            InitVel( mg, &v, Bnd_v, the_wind_fun, cur_time);
            LSInit( mg, lset.Phi, the_lset_fun, cur_time);
            // timedisc.Update(); // Called unconditionally in DoStep.

            //lset2.SetupSystem( make_P2Eval( mg, Bnd_v, v), P.get<double>("Time.StepSize"));

            std::cout << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            if (P.get<int>( "Levelset.VolCorr")) {
                double dphi= lset.AdjustVolume( Vol, 1e-9);
                std::cout << "volume correction is " << dphi << std::endl;
                lset.Phi.Data+= dphi;
                std::cout << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            }
        }
    }
    std::cout << std::endl;
}

void StationaryStrategy (DROPS::MultiGridCL& mg, DROPS::AdapTriangCL& adap, DROPS::LevelsetP2CL& lset)
{
    adap.MakeInitialTriang( sphere_dist);

    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
    // LinearLSInit( mg, lset.Phi, &sphere_dist);
    LSInit( mg, lset.Phi, sphere_dist, 0.);

    DROPS::IdxDescCL ifaceidx( P1IF_FE);
    ifaceidx.GetXidx().SetBound( P.get<double>("SurfTransp.OmitBound"));
    ifaceidx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
    std::cout << "NumUnknowns: " << ifaceidx.NumUnknowns() << std::endl;

    DROPS::MatDescCL M( &ifaceidx, &ifaceidx);
    DROPS::SetupInterfaceMassP1( mg, &M, lset.Phi, lset.GetBndData());
    std::cout << "M is set up.\n";
    DROPS::MatDescCL A( &ifaceidx, &ifaceidx);
    DROPS::SetupLBP1( mg, &A, lset.Phi, lset.GetBndData(), P.get<double>("SurfTransp.Visc"));
    DROPS::MatrixCL L;
    L.LinComb( 1.0, A.Data, 1.0, M.Data);
    DROPS::VecDescCL b( &ifaceidx);
    DROPS::SetupInterfaceRhsP1( mg, &b, lset.Phi, lset.GetBndData(), laplace_beltrami_0_rhs);

    //DROPS::WriteToFile( M.Data, "m_iface.txt", "M");
    //DROPS::WriteToFile( A.Data, "a_iface.txt", "A");
    //DROPS::WriteFEToFile( b, mg, "rhs_iface.txt", /*binary=*/ true);

    typedef DROPS::SSORPcCL SurfPcT;
    SurfPcT surfpc;
    typedef DROPS::PCGSolverCL<SurfPcT> SurfSolverT;
    SurfSolverT surfsolver( surfpc, P.get<int>("SurfTransp.Iter"), P.get<double>("SurfTransp.Tol"), true);

    DROPS::VecDescCL x( &ifaceidx);
    surfsolver.Solve( L, x.Data, b.Data);
    std::cout << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';

    if (P.get<int>( "SolutionOutput.Freq") > 0)
        DROPS::WriteFEToFile( x, mg, P.get<std::string>( "SolutionOutput.Path"), P.get<bool>( "SolutionOutput.Binary"));

    DROPS::IdxDescCL ifacefullidx( DROPS::P1_FE);
    ifacefullidx.CreateNumbering( mg.GetLastLevel(), mg);
    DROPS::VecDescCL xext( &ifacefullidx);
    DROPS::Extend( mg, x, xext);
    DROPS::NoBndDataCL<> nobnd;
    if (vtkwriter.get() != 0) {
        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "Levelset") );
        vtkwriter->Register( make_VTKIfaceScalar( mg, x, "InterfaceSol"));
        vtkwriter->Write( 0.);
    }

    double L2_err( L2_error( lset.Phi, lset.GetBndData(), make_P1Eval( mg, nobnd, xext), &laplace_beltrami_0_sol));
    std::cout << "L_2-error: " << L2_err << std::endl;
}

int main (int argc, char* argv[])
{
  try {
    std::ifstream param;
    if (argc != 2) {
        std::cout << "Using default parameter file: surfactant.json\n";
        param.open( "surfactant.json");
    }
    else
        param.open( argv[1]);
    if (!param)
        throw DROPS::DROPSErrCL( "main: error while opening parameter file\n");
    param >> P;
    param.close();
    std::cout << P << std::endl;

    std::cout << "Setting up interface-PDE.\n";
    WindVelocity= P.get<DROPS::Point3DCL>("Exp.Velocity");
    RadDrop=      P.get<DROPS::Point3DCL>("Exp.RadDrop");
    PosDrop=      P.get<DROPS::Point3DCL>("Exp.PosDrop");
    RadTorus=     P.get<DROPS::Point2DCL>("Exp.RadTorus");
    the_wind_fun= invecmap[P.get<std::string>("Exp.Wind")];
    the_lset_fun= inscamap[P.get<std::string>("Exp.Levelset")];
    the_rhs_fun=  inscamap[P.get<std::string>("Exp.Rhs")];
    the_sol_fun=  inscamap[P.get<std::string>("Exp.Solution")];
    for (Uint i= 0; i < 6; ++i)
        bf_wind[i]= the_wind_fun;

    std::cout << "Setting up domain:\n";
    std::auto_ptr<MGBuilderCL> builder( make_MGBuilder( P.get_child( "Domain")));
    DROPS::MultiGridCL mg( *builder);
    DROPS::AdapTriangCL adap( mg, P.get<double>("AdaptRef.Width"),
                              P.get<int>("AdaptRef.CoarsestLevel"), P.get<int>("AdaptRef.FinestLevel"));

    DROPS::LevelsetP2CL lset( mg, lsbnd, sf);

    if (P.get<int>("VTK.VTKOut",0))
        vtkwriter= std::auto_ptr<VTKOutCL>( new VTKOutCL(
            adap.GetMG(),
            "DROPS data",
            P.get<int>("Time.NumSteps")/P.get<int>("VTK.VTKOut") + 1,
            P.get<std::string>("VTK.VTKDir"),
            P.get<std::string>("VTK.VTKName"),
            P.get<std::string>("VTK.TimeFileName"),
            P.get<int>("VTK.Binary"), 
            P.get<bool>("VTK.UseOnlyP1"),
            -1,  /* <- level */
            P.get<bool>("VTK.ReUseTimeFile")));

    if (P.get<bool>( "Exp.StationaryPDE"))
        StationaryStrategy( mg, adap, lset);
    else
        Strategy( mg, adap, lset);

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
