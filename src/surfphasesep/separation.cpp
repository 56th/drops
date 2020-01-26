/// \file
/// \brief Solve a non-stationary Cahn-Hilliard equation on a moving interface
/// \author UH: Vladimir Yushutin

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
#include "levelset/levelsetmapper.h"
#include "out/output.h"
#include "misc/dynamicload.h"
#include "misc/funcmap.h"
#include "misc/omp_variable.h"
#include "geom/subtriangulation.h"
#include "num/gradient_recovery.h"
#include <cmath>
#include <fstream>
#include <string>
#include <tr1/unordered_map>
#include <tr1/unordered_set>
//#include <stokes/instatstokes2phase.h>

// vtk output
#include "surfnavierstokes/VTKWriter.hpp"
// logger
#include "SingletonLogger.hpp"

using namespace DROPS;

DROPS::ParamCL P;

// std::unique_ptr<VTKOutCL> vtkwriter;

DROPS::InVecMap& invecmap= DROPS::InVecMap::getInstance();
DROPS::InScaMap& inscamap= DROPS::InScaMap::getInstance();

InstatVectorFunction the_wind_fun;
InstatVectorFunction the_normal_fun;
InstatScalarFunction the_lset_fun;
InstatScalarFunction the_rhs_fun;
InstatScalarFunction the_conc_sol_fun;
InstatScalarFunction the_poten_sol_fun;
InstatScalarFunction the_species_sol_fun;

InstatScalarFunction the_zero_fun;

InstatVectorFunction the_sol_grad_fun;

using bnd_val_fun = InstatVectorFunction;

DROPS::BndCondT bc_wind[6]= { DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC };
bnd_val_fun bf_wind[6];

InstatScalarFunction sigma( 0);
SurfaceTensionCL sf( sigma, 0);
DROPS::LsetBndDataCL lsbnd( 6);

namespace ParameterNS {
    double h=0.1;
    double nu;
    double sigma;
    double eps;
}

// Surface divergence of a vector field w
inline double div_gamma_wind (const Point3DCL& n, const SMatrixCL<3,3>& dw)
{
   return trace( dw) - inner_prod( n, dw*n);
}

// laplace-beltrami of a function u
inline double laplace_beltrami_u (const Point3DCL& n,      const SMatrixCL<3,3>& dn,
                                  const Point3DCL& grad_u, const SMatrixCL<3,3>& Hess_u)
{
     const double tr_PHessu= trace( Hess_u) - inner_prod( n, Hess_u*n),
                  tr_Pdn= trace( dn) - inner_prod( n, dn*n),
                  ngradu= inner_prod( n, grad_u);
     return tr_PHessu - tr_Pdn*ngradu;
}


DROPS::Point3DCL WindVelocity, AngularVelocity;
std::string dirname;
DROPS::Point3DCL constant_wind (const DROPS::Point3DCL&, double)
{
    return ((0.0)*WindVelocity);
}
static RegisterVectorFunction regvec_constant_wind( "ConstantWind", constant_wind);

DROPS::Point3DCL RadDrop;
DROPS::Point3DCL PosDrop;

DROPS::Point3DCL constant_rotation (const DROPS::Point3DCL& pp, double)
{
    DROPS::Point3DCL velocity;
    cross_product(velocity, AngularVelocity, pp);
    return velocity;
}
static RegisterVectorFunction regvec_constant_rotation( "ConstantRotation", constant_rotation);

DROPS::Point3DCL isometric(const DROPS::Point3DCL& pp, double t)
{
     return (constant_rotation(pp,t)+constant_wind(pp,t));
}
static RegisterVectorFunction regvec_isometric( "Isometric", isometric);

DROPS::Point3DCL Rotate(const DROPS::Point3DCL& pp,
        double t)
{
    double u=AngularVelocity[0];
    double v=AngularVelocity[1];
    double w=AngularVelocity[2];
    double L = (u*u + v * v + w * w);
    double angle = std::sqrt(L)*(t);
    double u2 = u * u;
    double v2 = v * v;
    double w2 = w * w;

    SMatrixCL<3,3> rotationMatrix;

    if (!(L>0)) return pp;

    rotationMatrix(0,0) = (u2 + (v2 + w2) * cos(angle)) / L;
    rotationMatrix(0,1) = (u * v * (1 - cos(angle)) - w * sqrt(L) * sin(angle)) / L;
    rotationMatrix(0,2) = (u * w * (1 - cos(angle)) + v * sqrt(L) * sin(angle)) / L;

    rotationMatrix(1,0) = (u * v * (1 - cos(angle)) + w * sqrt(L) * sin(angle)) / L;
    rotationMatrix(1,1) = (v2 + (u2 + w2) * cos(angle)) / L;
    rotationMatrix(1,2) = (v * w * (1 - cos(angle)) - u * sqrt(L) * sin(angle)) / L;

    rotationMatrix(2,0) = (u * w * (1 - cos(angle)) - v * sqrt(L) * sin(angle)) / L;
    rotationMatrix(2,1) = (v * w * (1 - cos(angle)) + u * sqrt(L) * sin(angle)) / L;
    rotationMatrix(2,2) = (w2 + (u2 + v2) * cos(angle)) / L;
    return(rotationMatrix*pp);

}

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

// ==non-stationary test case "Constant function on rigid moving surface"
// "ConstantWind" 0.2
// "MovingEllipsoid" pos:[0,0,0], rad:[1,1,1]


double constant_sol (const Point3DCL& pp, double )
{
    return (0.1);
}
static RegisterScalarFunction regsca_constant_sol( "ConstantSolution", constant_sol);

double coscoscos_sol (const Point3DCL& pp, double )
{
    double a=0.05;
    Point3DCL p=(1./pp.norm())*pp;
    return (0.5+a*std::cos(2*3.141592*p[0])*std::cos(2*3.141592*p[1])*std::cos(2*3.141592*p[2]));
}
static RegisterScalarFunction regsca_coscoscos_sol( "CosxCosyCoszSolution", coscoscos_sol);



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
                 reaction= div_gamma_wind( n, dw)*u,
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

double axis_scaling_lset_ini (const Point3DCL& p, double)
{
    static const double t_end= P.get<double>( "Time.FinalTime");
    const double tout= t_end <= M_PI/2. ? t_end : M_PI/2.,
                 tin= t_end >= M_PI*3./2. ? M_PI*3./2. : (t_end >= M_PI ? t_end : 0.),
                 lout= axis_scaling_lset( p, tout),
                 lin= axis_scaling_lset( p, tin);
    return lout >= 0. ? lout :
           lin  <= 0. ? lin  :
                        0.;
}

double axis_scaling_sol (const Point3DCL& p, double t)
{
    return std::exp( -0.5*t)*heat_conduction_u0( p);
}
static RegisterScalarFunction regsca_axis_scaling_sol( "AxisScalingSol", axis_scaling_sol);

double axis_scaling_rhs (const Point3DCL& p, double t)
{
    const double a= axis_scaling( t);
    const double bf4= (p/MakePoint3D( std::pow( a, 2), 1., 1.)).norm_sq();
    const double bf6= (p/MakePoint3D( std::pow( a, 3), 1., 1.)).norm_sq();

    const double mat_der=  0.25*std::cos( t)/a - 0.5;
    const double reaction= 0.25*std::cos( t)/a/bf4*( std::pow( p[1], 2) + std::pow( p[2], 2));
    const double diffusion= (-2./std::pow( a, 2) - (1. + 1./std::pow( a, 2))*( 2. + 1./std::pow( a, 2) - bf6/bf4))/bf4;

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
//     const double err= div_gamma_wind( n, dw) - reaction,
//                errlb= laplace_beltrami_u( n, dn, grad_u, Hess_u) - diffusion*axis_scaling_sol( p, t);
//     if (std::fabs( err) > 1e-12 || std::fabs( errlb) > 1e-12) {
//         std::cerr << err   << " " << div_gamma_wind( n, dw) << " " << reaction << "\n"
//                   << errlb << " " << laplace_beltrami_u( n, dn, grad_u, Hess_u) << " " << diffusion*axis_scaling_sol( p, t) << "\n";
//         exit( 1);
//     }
//     return (mat_der + div_gamma_wind( n, dw))*axis_scaling_sol( p, t) - laplace_beltrami_u( n, dn, grad_u, Hess_u);

    return (mat_der + reaction - diffusion)*axis_scaling_sol( p, t);
}
static RegisterScalarFunction regsca_axis_scaling_rhs( "AxisScalingRhs", axis_scaling_rhs);


// ==non-stationary test case "Collision"==

const double collision_p= 3.0;

// sphere 1: c1= (-1.5 0 0) + v*t
Point3DCL collision_center_1 (double t)
{
    return MakePoint3D( -1.5, 0., 0.) + WindVelocity*t;
}
// sphere 2: c2= ( 1.5 0 0) - v*t = -c1
Point3DCL collision_center_2 (double t)
{
    return -collision_center_1( t);
}

double collision_Dt_lset (const DROPS::Point3DCL& x, double t)
{
    Point3DCL x1= x - collision_center_1( t),
              x2= x - collision_center_2( t);
    const double n1= x1.norm() < 1e-3 ? 1e-3 : std::pow( x1.norm(), collision_p + 2.),
                 n2= x2.norm() < 1e-3 ? 1e-3 : std::pow( x2.norm(), collision_p + 2.);
    return collision_p*(DROPS::inner_prod( WindVelocity, x1)/n1
                      - DROPS::inner_prod( WindVelocity, x2)/n2);
}

DROPS::Point3DCL collision_Dx_lset (const DROPS::Point3DCL& x, double t)
{
    Point3DCL x1= x - collision_center_1( t),
              x2= x - collision_center_2( t);
    const double n1= x1.norm() < 1e-3 ? 1e-3 : std::pow( x1.norm(), collision_p + 2.),
                 n2= x2.norm() < 1e-3 ? 1e-3 : std::pow( x2.norm(), collision_p + 2.);
    return collision_p*(x1/n1 + x2/n2);
}

DROPS::Point3DCL collision_wind (const DROPS::Point3DCL& x, double t)
{
    const Point3DCL Dphi= collision_Dx_lset( x, t);
    const double Dtphi=   collision_Dt_lset( x, t);
    const double n= Dphi.norm_sq();
    return n < 1e-6 ? DROPS::Point3DCL() : (Dtphi/n)*Dphi;
}
static RegisterVectorFunction regvec_collision_wind( "CollisionWind", collision_wind);

double collision_lset (const Point3DCL& x, double t)
{
    const Point3DCL c1= collision_center_1( t);
    const Point3DCL x1= x - c1;
    const Point3DCL c2= collision_center_2( t);
    const Point3DCL x2= x - c2;

    const double n1= x1.norm() < 1e-3 ? 1e-3 : std::pow( x1.norm(), collision_p),
                 n2= x2.norm() < 1e-3 ? 1e-3 : std::pow( x2.norm(), collision_p);
    return RadDrop[0] - 1./n1  - 1./n2;
}
static RegisterScalarFunction regsca_collision_lset( "CollisionLset", collision_lset);

double collision_sol (const Point3DCL& x, double)
{
    return 2.*std::cos( x[0])*std::cos(M_PI*x[1]);
}
static RegisterScalarFunction regsca_collision_sol( "CollisionSol", collision_sol);

double collision_sol2 (const Point3DCL& x, double)
{
    return x[0] < 0. ? 0 : 3. - x[0];
}
static RegisterScalarFunction regsca_collision_sol2( "CollisionSol2", collision_sol2);

double zero_function (const Point3DCL& , double)
{
    return 0;
}
static RegisterScalarFunction regsca_zero_sca_fun( "ZeroScalarFun", zero_function);


double collision_rhs (const Point3DCL&, double)
{
    return 0;
}
static RegisterScalarFunction regsca_collision_rhs( "CollisionRhs", collision_rhs);


DROPS::Point3DCL mergedrop_wind (const DROPS::Point3DCL& p, double t)
{
    Point3DCL wd;

    DROPS::Point3DCL  c1,c2;
    double velocity=0.1;
    c1=MakePoint3D( -1.5, 0., 0.) + MakePoint3D( 1.0, 0., 0.)*t;
    c2=-c1;
   /* c1[0]=1.5*(t-1),c1[1]=0,c1[2]=0;
    c2[0]=-1.5*(t-1),c2[1]=0,c2[2]=0;*/

    DROPS::Point3DCL x1( p - c1);
    DROPS::Point3DCL x2( p - c2);

    //Point3DCL wd;
    /* if(x1.norm()!=0&&x2.norm()!=0)
     {
         double a=-4.5*x1[0]/std::pow(x1.norm(),5) + 4.5*x2[0]/std::pow(x2.norm(),5);
         wd=3/std::pow(x1.norm(),5)*x1+3/std::pow(x2.norm(),5)*x2;
         if(wd.norm_sq()!=0)
             wd=-a/wd.norm_sq()*wd;
         else
            wd=-a*wd;
     }*/
    double a=-4.5*x1[0]*std::pow(x2.norm(),5) + 4.5*x2[0]*std::pow(x1.norm(),5);
    wd=3*std::pow(x2.norm(),5)*x1+3*std::pow(x1.norm(),5)*x2;
    if(wd.norm_sq()!=0)
        wd=-a/wd.norm_sq()*wd ;
    else
        wd=-a*wd;
    return wd;
}
static RegisterVectorFunction regvec_mergedrop_wind( "MergeDropWind", mergedrop_wind);

DROPS::Point3DCL mergedrop_compositewind (const DROPS::Point3DCL& p, double t)
{
    double time_rate=1.0;
    double time_shift=0.2385;
    t=time_shift+time_rate*t;

    double scaler;
    scaler=1.;


    //if (t>1.5) return Point3DCL(0,0,0);
    double width=0.1;
   double a=    1.-  std::tanh((1./(width))*(std::abs(p[0])));

   double b= std::tanh((1./(width))*(-p[0]));

   double norm=1;//std::sqrt(a*a+b*b);

   return (time_rate)*((a/norm)*mergedrop_wind(p,t)+(b/norm)*MakePoint3D( 1, 0., 0.));//constant_wind(p,t));
}
static RegisterVectorFunction regvec_mergedrop_compositewind( "MergeDropCompositeWind", mergedrop_compositewind);


double mergedrop_lset (const Point3DCL& p, double t)
{
    double time_rate=1.0;
    double time_shift=0.2385;
    t=time_shift+t*time_rate;

    //if (t>1.5) t=1.5;

    DROPS::Point3DCL  c1,c2;
/*    c1[0]=1.5*(t-5),c1[1]=0,c1[2]=0;
    c2[0]=-1.5*(t-5),c2[1]=0,c2[2]=0;*/
    c1=MakePoint3D( -1.5, 0., 0.) + MakePoint3D( 1.0, 0., 0.)*t;
    c2=-c1;
    DROPS::Point3DCL x1( p - c1);
    DROPS::Point3DCL x2( p - c2);
    // std::cout<<1.0-1./std::pow(x1.norm(),3)<<std::endl;
    if(x1.norm()!=0&&x2.norm()!=0)
        return  1- 1./std::pow(x1.norm(),3)-1.0/std::pow(x2.norm(),3);
    else
        return   -10*(std::pow(x1.norm(),3)*std::pow(x2.norm(),3)-std::pow(x1.norm(),3)-std::pow(x2.norm(),3));//-1./std::pow(x2.norm(),3);

}
static RegisterScalarFunction regsca_mergedrop_lset( "MergeDropLset", mergedrop_lset);

DROPS::Point3DCL normal_mergedrop (const Point3DCL& p, double t)
{
    double time_rate=1.0;
    double time_shift=0.2385;
    t=time_shift+t*time_rate;

    DROPS::Point3DCL  c1,c2;
  /*  c1[0]=1.5*(t-5),c1[1]=0,c1[2]=0;
    c2[0]=-1.5*(t-5),c2[1]=0,c2[2]=0;*/
    c1=MakePoint3D( -1.5, 0., 0.) + MakePoint3D( 1.0, 0., 0.)*t;
    c2=-c1;

    DROPS::Point3DCL x1( p - c1);
    DROPS::Point3DCL x2( p - c2);

    DROPS::Point3DCL grad;
    if(x1.norm()!=0&&x2.norm()!=0){
        grad=-3*(std::pow(x1.norm(),5)*x1+std::pow(x2.norm(),5)*x2);
        return  ((time_rate/std::abs(time_rate))*(grad.norm()>0.000001?grad/grad.norm():grad));
    }
    else
        return   0*grad;//-1./std::pow(x2.norm(),3);

}
static RegisterVectorFunction regsca_normal_mergedrop( "NormalMergeDrop", normal_mergedrop);


double mergedrop_sol (const Point3DCL& p, double t)
{
    double random = ((double)rand() / RAND_MAX);
    double k=0.5;
    double ampl=0.1;
    //return 0.5*(1.+p[0]);
    //if(p[0]>0)
    {
        if (random<0.5) return(0.05 + ampl*(2*k*random - 0.5));
        else return(0.05+ampl*(2*(1-k)*random + 2*k-1- 0.5));
    }
//    if(p[0]>=0)
     //   return(2./3. + ampl*);
//
//    else     return(0);

   /* if(p[0]>=0)
            return((std::tanh((1./(2*std::sqrt(2)*ParameterNS::eps))*(p[0]-1.5+0.235)/pow(pow(p[0]-1.5, 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.5))  + 0.1e1) / 0.2e1);

    else     return((std::tanh((1./(2*std::sqrt(2)*ParameterNS::eps))*(-p[1])/pow(pow(p[0]+1.5, 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.5))  + 0.1e1) / 0.2e1);
*/
}
static RegisterScalarFunction regsca_mergedrop_sol( "MergeDropSol", mergedrop_sol);

double mergedrop_randomsol (const Point3DCL& p, double t)
{
    double random = ((double)rand() / RAND_MAX);
    double k=0.5;
    //return 0.5*(1.+p[0]);
    //if(p[0]>0)
    {
        if (random<0.5) return(2*k*random);
        else return(2*(1-k)*random + 2*k-1);
    }
    //else return 0;
}
static RegisterScalarFunction regsca_mergedrop_randomsol( "MergeDropRandomSol", mergedrop_randomsol);

DROPS::Point3DCL mergedrop_surfgradsol (const DROPS::Point3DCL& p, double t)
{
    Point3DCL surfgrad;
    return surfgrad;
}
static RegisterVectorFunction regvec_mergedrop_surfgradsol( "MergeDropSurfGradSol", mergedrop_surfgradsol);


double mergedrop_rhs (const Point3DCL& p, double t)
{
    return 0; //Ex3
}
static RegisterScalarFunction regsca_mergedrop_rhs( "MergeDropRhs", mergedrop_rhs);
//////////////////////////////////////////////////////////////////////////////////////////////

// ==Some spheres==


double sphere_dist (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL x( p - PosDrop);
    return x.norm() - RadDrop[0];
}
static RegisterScalarFunction regsca_sphere_dist_lset( "Sphere", sphere_dist);

DROPS::Point3DCL d_sphere_dist (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL x( p - PosDrop);
    return x/x.norm();
}
static RegisterVectorFunction regsca_sphere_dist_normal( "NormalSphere", d_sphere_dist);


using dist_funT = InstatScalarFunction;

double sphere_2move (const DROPS::Point3DCL& p, double t)
{
    DROPS::Point3DCL x( p - (PosDrop + t*constant_wind(p, t)));
    return x.norm() - RadDrop[0];
}

static RegisterScalarFunction regsca_moving_sphere_lset( "MovingEllipsoid", sphere_2move);

DROPS::Point3DCL normal_sphere_2move (const DROPS::Point3DCL& p, double t)
{
    DROPS::Point3DCL x= (p - PosDrop - t*constant_wind(p, t));//--constant wind
    if(x.norm()>0.000001)
        return x/x.norm();
    else
        return x;
}
static RegisterVectorFunction regvec_moving_sphere_normal( "NormalMovingEllipsoid", normal_sphere_2move);

double dilating_radius(double t)
{
    double A=0.2;
    double T=5;
    return(1.+A*std::sin(T*2.*3.141592*t));
}
double dilating_radius_prime(double t)
{
    double A=0.2;
    double T=5;
    return(T*2.*3.141592*A*std::cos(T*2.*3.141592*t));
}

double sphere_2dilate (const DROPS::Point3DCL& p, double t)
{
    DROPS::Point3DCL x((1.0/dilating_radius(t))*(p - PosDrop));
    //x[0]=p[0]/dilating_radius(t);
    return x.norm() - RadDrop[0];
}

static RegisterScalarFunction regsca_dilatating_sphere_lset( "DilatatingEllipsoid", sphere_2dilate);

DROPS::Point3DCL normal_sphere_2dilate (const DROPS::Point3DCL& p, double t)
{
    DROPS::Point3DCL x= (1.0/(dilating_radius(t)*dilating_radius(t)))*(p - PosDrop );
    //x[0]=p[0]/(dilating_radius(t)*dilating_radius(t));

    if(x.norm()>0.000001)
        return x/x.norm();
    else
        return x;
}
static RegisterVectorFunction regvec_dilate_sphere_normal( "NormalDilatatingEllipsoid", normal_sphere_2dilate);

DROPS::Point3DCL sin_dilatation (const DROPS::Point3DCL& pp, double t)
{
    DROPS::Point3DCL x((dilating_radius_prime(t)/dilating_radius(t))*(pp - PosDrop));
    return x;

}
static RegisterVectorFunction regvec_sin_dilatation( "UniformDilatation", sin_dilatation);


double sphere_2squeez (const DROPS::Point3DCL& p, double t)
{
    double a=1.0*dilating_radius(t);
    double b=1.0;
    double c=1.0;
    DROPS::Point3DCL x;
    x[0]=(p[0]-PosDrop[0])/a;
    x[1]=(p[1]-PosDrop[1])/b;
    x[2]=(p[2]-PosDrop[2])/c;
    return x.norm() - 1.;
}

static RegisterScalarFunction regsca_squeezing_sphere_lset( "SqueezingEllipsoid", sphere_2squeez);

DROPS::Point3DCL normal_sphere_2squeez (const DROPS::Point3DCL& p, double t)
{
    double a=1.0*dilating_radius(t);
    double b=1.0;
    double c=1.0;

    DROPS::Point3DCL x;
    x[0]=(p[0]-PosDrop[0])/a;
    x[1]=(p[1]-PosDrop[1])/b;
    x[2]=(p[2]-PosDrop[2])/c;

    if(x.norm()>0.000001)
        return x/x.norm();
    else
        return x;
}
static RegisterVectorFunction regvec_squeezing_sphere_normal( "NormalSqueezingEllipsoid", normal_sphere_2squeez);

DROPS::Point3DCL sin_squeezing (const DROPS::Point3DCL& p, double t)
{
    double a=1.0*dilating_radius(t);
    double b=1.0;
    double c=1.0;

    DROPS::Point3DCL primes;
    primes[0]=1.0*dilating_radius_prime(t);
    primes[1]=0.0;
    primes[2]=0.0;

    DROPS::Point3DCL x;
    x[0]=std::pow( (p[0]-PosDrop[0])/a , 2);
    x[1]=std::pow( (p[1]-PosDrop[1])/b , 2);
    x[2]=std::pow( (p[2]-PosDrop[2])/c , 2);

    DROPS::Point3DCL gradf;
    gradf[0]=(p[0]-PosDrop[0])/a;
    gradf[1]=(p[1]-PosDrop[1])/b;
    gradf[2]=(p[2]-PosDrop[2])/c;

    double norm_comp=inner_prod(x,primes)/x.norm_sq();
    return norm_comp*gradf;

}
static RegisterVectorFunction regvec_squeezing( "SqueezingWind", sin_squeezing);


SMatrixCL<3,3> dp_sphere (const DROPS::Point3DCL& x, double)
{
    const double normx= x.norm();
    return normx == 0. ? SMatrixCL<3,3>() : RadDrop[0]/normx*(eye<3,3>() - outer_product( x/normx, x/normx));
}

double abs_det_sphere (const TetraCL& tet, const BaryCoordCL& xb, const SurfacePatchCL& p)
{
    if (p.empty())
        return 0.;

    // Compute the jacobian of p.
    SMatrixCL<3,3> dp= dp_sphere( GetWorldCoord( tet, xb), 0.);

    const Bary2WorldCoordCL b2w( tet);
    QRDecompCL<3,2> qr;
    SMatrixCL<3,2>& M= qr.GetMatrix();
    M.col(0, b2w( p.vertex_begin()[1]) - b2w( p.vertex_begin()[0]));
    M.col(1, b2w( p.vertex_begin()[2]) - b2w( p.vertex_begin()[0]));
    qr.prepare_solve();
    SMatrixCL<3,2> U;
    Point3DCL tmp;
    for (Uint i= 0; i < 2; ++i) {
        tmp= std_basis<3>( i + 1);
        qr.apply_Q( tmp);
        U.col( i, tmp);
    }
    const SMatrixCL<2,2> Gram= GramMatrix( dp*U);
    return std::sqrt( Gram(0,0)*Gram(1,1) - Gram(0,1)*Gram(1,0));
}

// ==stationary test case "LaplaceBeltrami0"==
// Sphere around 0, RadDrop 1, wind == 0
// "Levelset": "SphereDist"
// A right hand side from C.J. Heine...
// const double a( -13./8.*std::sqrt( 35./M_PI));
const double a( 12.);
double laplace_beltrami_0_rhs (const DROPS::Point3DCL& p, double)
{
    return a/std::pow( p.norm(), 3.)*(3.*p[0]*p[0]*p[1] - p[1]*p[1]*p[1]);
}
static RegisterScalarFunction regsca_laplace_beltrami_0_rhs( "LaplaceBeltrami0Rhs", laplace_beltrami_0_rhs);

// ...and the corresponding solution (extended)
double laplace_beltrami_0_sol (const DROPS::Point3DCL& p, double)
{
    return p.norm_sq()/(12. + p.norm_sq())*laplace_beltrami_0_rhs( p, 0.);
//     return 1./12.*laplace_beltrami_0_rhs( p, 0.);
//    return 1. + p.norm_sq()/(12. + p.norm_sq())*laplace_beltrami_0_rhs( p, 0.);
}
static RegisterScalarFunction regsca_laplace_beltrami_0_sol( "LaplaceBeltrami0Sol", laplace_beltrami_0_sol);

// The tangential gradient of laplace_beltrami_0_sol with respect to the exact sphere.
DROPS::Point3DCL laplace_beltrami_0_sol_grad (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL tmp= 3./std::pow( p.norm(), 3)*( MakePoint3D(2.*p[0]*p[1], p[0]*p[0] - p[1]*p[1], 0.) - (3.*p[0]*p[0]*p[1] - std::pow(p[1], 3))/p.norm_sq()*p);
    return tmp; // This equals tmp - inner_prod( p/p.norm(), tmp)*p/p.norm().
}

double sol0t (const DROPS::Point3DCL& p, double t)
{
    const Point3DCL q( p - (PosDrop + t*constant_wind(p, t)));
    const double val( a*(3.*q[0]*q[0]*q[1] - q[1]*q[1]*q[1]));

//    return q.norm_sq()/(12. + q.norm_sq())*val;
    return 1. + q.norm_sq()/(12. + q.norm_sq())*val;
}

// test on a unit sphere:  mobility=const, stationary omega=harmonic, chi=harmonic/2
double Test1_chiSol (const DROPS::Point3DCL& pp, double t)
{
    DROPS::Point3DCL p(pp  - t*constant_wind(pp,t));
    return (1./2.+(1./2.)*(sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) / 0.2e1));

}
static RegisterScalarFunction regsca_spherical_harmonic( "1/2+1/2*First spherical harmonic", Test1_chiSol);

double Test1_omegaSol (const DROPS::Point3DCL& pp, double t)
{
    DROPS::Point3DCL p(pp  - t*constant_wind(pp,t));

    return ((1.)*(sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) / 0.2e1));

}
static RegisterScalarFunction regsca_half_spherical_harmonic( "First spherical harmonic", Test1_omegaSol);

double Test1_rhs3 (const DROPS::Point3DCL& pp, double t)
{
    DROPS::Point3DCL p(pp  - t*constant_wind(pp,t));

    if (p.norm()<0.01) return (0);
    return ((2.)*(sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) / 0.2e1));

}
static RegisterScalarFunction regsca_rhs1_harmonic( "First spherical harmonic rhs", Test1_rhs3);


double TestStationary_chiSol (const DROPS::Point3DCL& pp, double t)
{
    //DROPS::Point3DCL p(pp  - t*constant_wind(pp,t));
    DROPS::Point3DCL p=Rotate(pp  - t*constant_wind(pp,t), t) - PosDrop;
    double shift=0;

   return((std::tanh(shift+(1./(2*std::sqrt(2)*ParameterNS::eps))*p[1]/pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.5))  + 0.1e1) / 0.2e1);
}

static RegisterScalarFunction regsca_stationary_harmonic( "SteadyPhases", TestStationary_chiSol);



double Test10_chiSol (const DROPS::Point3DCL& pp, double t)
{
    //DROPS::Point3DCL p(pp  - t*constant_wind(pp,t));
    DROPS::Point3DCL p=Rotate(pp  - t*constant_wind(pp,t), -t) - PosDrop;

    double a=0.0;
    double T=0.25;
    return((0.1e1 - a * exp(-t / T)) * (p[0] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[1] + 0.1e1) / 0.2e1);
}
static RegisterScalarFunction regsca_secondspherical_harmonic( "Second spherical harmonic", Test10_chiSol);

double Test10_omegaSol (const DROPS::Point3DCL& pp, double t)
{
    //DROPS::Point3DCL p(pp  - t*constant_wind(pp,t));
    //DROPS::Point3DCL p=Rotate(pp, -t);
 DROPS::Point3DCL p=Rotate(pp  - t*constant_wind(pp,t), -t)- PosDrop;
    double a=0.0;
    double T=0.25;
    if (p.norm()<0.01) return (0);

    double diffusion= -0.3e1 * ParameterNS::eps * ParameterNS::eps * (a * exp(-t / T) - 0.1e1) * p[1] * p[0] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    double well_potential= -(a * exp(-t / T) - 0.1e1) * (pow(p[0], 0.2e1) + p[0] * p[1] + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * (exp(-t / T) * a * pow(p[0], 0.2e1) + exp(-t / T) * a * p[0] * p[1] + exp(-t / T) * a * pow(p[1], 0.2e1) + exp(-t / T) * a * pow(p[2], 0.2e1) + pow(p[0], 0.2e1) - p[0] * p[1] + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * (exp(-t / T) * a * pow(p[0], 0.2e1) + exp(-t / T) * a * p[0] * p[1] + exp(-t / T) * a * pow(p[1], 0.2e1) + exp(-t / T) * a * pow(p[2], 0.2e1) - p[0] * p[1]) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1) / 0.8e1;
    return((diffusion+well_potential)/ParameterNS::eps);
}
static RegisterScalarFunction regsca_omega_secondspherical_harmonic( "Omega from Second spherical harmonic", Test10_omegaSol);

double Test10_rhs3 (const DROPS::Point3DCL& pp, double  t) {
    //DROPS::Point3DCL p(pp  - t*constant_wind(pp,t));
    //DROPS::Point3DCL p=Rotate(pp, -t);
 DROPS::Point3DCL p=Rotate(pp - t*constant_wind(pp,t), -t)- PosDrop;
    double a=0.0;
    double T=0.25;
    double inertia=a / T * exp(-t / T) * (p[0] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[1] + 0.1e1) / 0.2e1;
    double bidiffusion =  -0.3e1 / 0.16e2 * ParameterNS::sigma * ((pow(p[1], 0.4e1) - 0.3e1 / 0.2e1 * pow(p[1], 0.3e1) * p[0] + (-0.7e1 / 0.2e1 * pow(p[0], 0.2e1) + pow(p[2], 0.2e1)) * pow(p[1], 0.2e1) + (-0.3e1 / 0.2e1 * pow(p[0], 0.3e1) - 0.3e1 / 0.2e1 * p[0] * pow(p[2], 0.2e1)) * p[1] + pow(p[0], 0.4e1) + pow(p[0], 0.2e1) * pow(p[2], 0.2e1)) * pow(a, 0.3e1) * pow(pow(p[0], 0.2e1) + p[0] * p[1] + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.3e1) * pow(exp(-t / T), 0.3e1) - 0.9e1 / 0.2e1 * p[0] * (pow(p[1], 0.4e1) - 0.2e1 / 0.3e1 * pow(p[1], 0.3e1) * p[0] + (-0.5e1 / 0.3e1 * pow(p[0], 0.2e1) + 0.4e1 / 0.3e1 * pow(p[2], 0.2e1)) * pow(p[1], 0.2e1) + (-0.2e1 / 0.3e1 * pow(p[0], 0.3e1) - 0.2e1 / 0.3e1 * p[0] * pow(p[2], 0.2e1)) * p[1] + pow(p[0], 0.4e1) + 0.4e1 / 0.3e1 * pow(p[0], 0.2e1) * pow(p[2], 0.2e1) + pow(p[2], 0.4e1) / 0.3e1) * a * a * p[1] * pow(pow(p[0], 0.2e1) + p[0] * p[1] + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.2e1) * pow(exp(-t / T), 0.2e1) - 0.2e1 / 0.3e1 * a * (pow(p[0], 0.2e1) + p[0] * p[1] + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * (pow(p[1], 0.8e1) - 0.3e1 / 0.4e1 * pow(p[1], 0.7e1) * p[0] + (0.3e1 * pow(p[2], 0.2e1) - (double) (6 * ParameterNS::eps * ParameterNS::eps) - 0.39e2 / 0.4e1 * pow(p[0], 0.2e1)) * pow(p[1], 0.6e1) - 0.9e1 / 0.4e1 * p[0] * (-(double) (8 * ParameterNS::eps * ParameterNS::eps) + pow(p[2], 0.2e1)) * pow(p[1], 0.5e1) + (0.3e1 * pow(p[2], 0.4e1) + (-(double) (12 * ParameterNS::eps * ParameterNS::eps) - 0.14e2 * pow(p[0], 0.2e1)) * pow(p[2], 0.2e1) + 0.24e2 * pow(p[0], 0.2e1) * (double) (ParameterNS::eps * ParameterNS::eps) + 0.13e2 / 0.4e1 * pow(p[0], 0.4e1)) * pow(p[1], 0.4e1) - 0.9e1 / 0.4e1 * p[0] * (p[2] - (double) (4 * ParameterNS::eps)) * (p[2] + (double) (4 * ParameterNS::eps)) * (pow(p[0], 0.2e1) + pow(p[2], 0.2e1)) * pow(p[1], 0.3e1) - 0.39e2 / 0.4e1 * (-0.4e1 / 0.39e2 * pow(p[2], 0.4e1) + (0.17e2 / 0.39e2 * pow(p[0], 0.2e1) + 0.8e1 / 0.13e2 * (double) ParameterNS::eps * (double) ParameterNS::eps) * pow(p[2], 0.2e1) + pow(p[0], 0.4e1) - 0.32e2 / 0.13e2 * pow(p[0], 0.2e1) * (double) (ParameterNS::eps * ParameterNS::eps)) * (pow(p[0], 0.2e1) + pow(p[2], 0.2e1)) * pow(p[1], 0.2e1) - 0.3e1 / 0.4e1 * p[0] * pow(pow(p[0], 0.2e1) + pow(p[2], 0.2e1), 0.2e1) * (-(double) (24 * ParameterNS::eps * ParameterNS::eps) + pow(p[0], 0.2e1) + pow(p[2], 0.2e1)) * p[1] + pow(p[0], 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[2], 0.2e1), 0.2e1) * (-(double) (6 * ParameterNS::eps * ParameterNS::eps) + pow(p[0], 0.2e1) + pow(p[2], 0.2e1))) * exp(-t / T) + 0.7e1 / 0.6e1 * (pow(p[1], 0.8e1) + (0.24e2 / 0.7e1 * pow(p[2], 0.2e1) - 0.96e2 / 0.7e1 * (double) ParameterNS::eps * (double) ParameterNS::eps - 0.6e1 / 0.7e1 * pow(p[0], 0.2e1)) * pow(p[1], 0.6e1) + (0.30e2 / 0.7e1 * pow(p[2], 0.4e1) + (-0.264e3 / 0.7e1 * (double) ParameterNS::eps * (double) ParameterNS::eps + 0.10e2 / 0.7e1 * pow(p[0], 0.2e1)) * pow(p[2], 0.2e1) - 0.120e3 / 0.7e1 * pow(p[0], 0.2e1) * (double) (ParameterNS::eps * ParameterNS::eps) + pow(p[0], 0.4e1)) * pow(p[1], 0.4e1) - 0.6e1 / 0.7e1 * (pow(p[0], 0.2e1) + pow(p[2], 0.2e1)) * (-0.8e1 / 0.3e1 * pow(p[2], 0.4e1) + (-0.8e1 / 0.3e1 * pow(p[0], 0.2e1) + (double) (40 * ParameterNS::eps * ParameterNS::eps)) * pow(p[2], 0.2e1) + pow(p[0], 0.4e1) + 0.20e2 * pow(p[0], 0.2e1) * (double) (ParameterNS::eps * ParameterNS::eps)) * pow(p[1], 0.2e1) + (0.3e1 / 0.7e1 * pow(p[2], 0.4e1) + (0.10e2 / 0.7e1 * pow(p[0], 0.2e1) - 0.72e2 / 0.7e1 * (double) ParameterNS::eps * (double) ParameterNS::eps) * pow(p[2], 0.2e1) + pow(p[0], 0.4e1) - 0.96e2 / 0.7e1 * pow(p[0], 0.2e1) * (double) (ParameterNS::eps * ParameterNS::eps)) * pow(pow(p[0], 0.2e1) + pow(p[2], 0.2e1), 0.2e1)) * p[0] * p[1]) * (a * (pow(p[0], 0.2e1) + p[0] * p[1] + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * exp(-t / T) + pow(p[0], 0.2e1) - p[0] * p[1] + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(a * exp(-t / T) - 0.1e1, 0.3e1) * (pow(p[0], 0.2e1) + p[0] * p[1] + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(pow(0.1e1 - a * exp(-t / T), 0.2e1) * pow(p[0] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[1] + 0.1e1, 0.2e1) * pow(0.1e1 - (0.1e1 - a * exp(-t / T)) * (p[0] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[1] + 0.1e1) / 0.2e1, 0.2e1), -0.1e1 / 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.8e1);
    return(inertia - bidiffusion/ParameterNS::eps);
}
static RegisterScalarFunction regsca_rhs10_harmonic( "Rhs from Second spherical harmonic", Test10_rhs3);



double Test6_chiSol (const DROPS::Point3DCL& pp, double t)
{
        DROPS::Point3DCL p(pp  - t*constant_wind(pp,t));

    double a=0.8;
    double T=5;
    return ((0.1e1 - a * exp(-t / T)) * (sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] *
                                         pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) / 0.2e1 +
                                         0.1e1 / 0.2e1));
}
static RegisterScalarFunction regsca_conc_decay( "DecayFirstHarmonicConc", Test6_chiSol);

double Test6_omegaSol (const DROPS::Point3DCL& pp, double t)
{
        DROPS::Point3DCL p(pp  - t*constant_wind(pp,t));

    double a=0.8;
    double T=5;
    double diffusion= - ParameterNS::eps * (a * exp(-t / T) - 0.1e1) * sqrt(0.3e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1);
    double well_potential= ( -1. / ParameterNS::eps) * (a * exp(-t / T) - 0.1e1) * (sqrt(0.3e1) * p[2] +
                                                                                    sqrt(0.3141592654e1) * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1))) *
                           (exp(-t / T) * sqrt(0.3e1) * a * p[2] + exp(-t / T) * sqrt(0.3141592654e1) * sqrt(pow(p[0], 0.2e1) +
                                                                                                             pow(p[1], 0.2e1) +
                                                                                                             pow(p[2], 0.2e1)) * a -
                            sqrt(0.3e1) * p[2] + sqrt(0.3141592654e1) * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)))
                           * (exp(-t / T) * sqrt(0.3e1) * a * p[2] + exp(-t / T) * sqrt(0.3141592654e1) * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * a - sqrt(0.3e1) * p[2])
                           * pow(0.3141592654e1, -0.3e1 / 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1) / 0.8e1;
    return(diffusion+well_potential);
}
static RegisterScalarFunction regsca_poten_decay( "DecayFirstHarmonicPoten", Test6_omegaSol);


double Test6_rhs3 (const DROPS::Point3DCL& pp, double  t) {
        DROPS::Point3DCL p(pp  - t*constant_wind(pp,t));

    double a=0.8;
    double T=5;
    double inertia=a / T * exp(-t / T) * (sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] *
                                          pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1)
                                          / 0.2e1 + 0.1e1 / 0.2e1);
    double bidiffusion = -0.9e1 / 0.32e2 * ParameterNS::sigma * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.15e2)
                         * pow(pow(a * exp(-t / T) - 0.1e1, 0.2e1) * pow(sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] *
                                                                         (a * exp(-t / T) - 0.1e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) +
                                                                         a * exp(-t / T) + 0.1e1, 0.2e1) * pow(sqrt(0.3e1) *pow(0.3141592654e1,-0.1e1 / 0.2e1) * p[2] *
                                                                                                               pow(pow(p[0], 0.2e1)+ pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) + 0.1e1, 0.2e1), -0.1e1 / 0.2e1)
                         * pow(a * exp(-t / T) - 0.1e1, 0.3e1) * (-pow(0.3141592654e1, 0.6e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.21e2 / 0.2e1) * sqrt(0.3e1) * a * (pow(a, 0.3e1) * ((pow(0.3141592654e1, 0.3e1) + 0.21e2 * 0.3141592654e1 * 0.3141592654e1 + 0.63e2 * 0.3141592654e1 + 0.27e2) * pow(p[2], 0.6e1) + 0.3e1 * (0.3141592654e1 - 0.3e1) * (0.3141592654e1 * 0.3141592654e1 + 0.11e2 * 0.3141592654e1 + 0.6e1) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * pow(p[2], 0.4e1) + 0.3e1 * 0.3141592654e1 * (0.3141592654e1 * 0.3141592654e1 - 0.5e1 * 0.3141592654e1 - 0.48e2) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1), 0.2e1) * pow(p[2], 0.2e1) + 0.3141592654e1 * 0.3141592654e1 * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1), 0.3e1) * (0.3141592654e1 - 0.18e2)) * pow(exp(-t / T), 0.3e1) + 0.2e1 * a * a * ((pow(0.3141592654e1, 0.3e1) - 0.63e2 * 0.3141592654e1 - 0.54e2) * pow(p[2], 0.6e1) + 0.3e1 * (pow(0.3141592654e1, 0.3e1) + 0.27e2 * 0.3141592654e1 + 0.36e2) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * pow(p[2], 0.4e1) + 0.3e1 * 0.3141592654e1 * (0.3141592654e1 * 0.3141592654e1 + 0.48e2) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1), 0.2e1) * pow(p[2], 0.2e1) + pow(0.3141592654e1, 0.3e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1), 0.3e1)) * pow(exp(-t / T), 0.2e1) + 0.2e1 / 0.3e1 * a * ((pow(0.3141592654e1, 0.3e1) - 0.27e2 * 0.3141592654e1 * 0.3141592654e1 + 0.63e2 * 0.3141592654e1 + 0.243e3) * pow(p[2], 0.6e1) + ((0.3e1 * pow(0.3141592654e1, 0.3e1) - 0.27e2 * 0.3141592654e1 * 0.3141592654e1 - 0.108e3 * 0.3141592654e1 - 0.486e3) * pow(p[0], 0.2e1) + (0.3e1 * pow(0.3141592654e1, 0.3e1) - 0.27e2 * 0.3141592654e1 * 0.3141592654e1 - 0.108e3 * 0.3141592654e1 - 0.486e3) * pow(p[1], 0.2e1) + 0.4e1 * 0.3141592654e1 * ParameterNS::eps * ParameterNS::eps * (0.3141592654e1 * 0.3141592654e1 + 0.18e2 * 0.3141592654e1 + 0.9e1)) * pow(p[2], 0.4e1) + 0.3e1 * 0.3141592654e1 * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * ((0.3141592654e1 * 0.3141592654e1 + 0.9e1 * 0.3141592654e1 - 0.57e2) * pow(p[0], 0.2e1) + (0.3141592654e1 * 0.3141592654e1 + 0.9e1 * 0.3141592654e1 - 0.57e2) * pow(p[1], 0.2e1) + 0.8e1 / 0.3e1 * (0.3141592654e1 * 0.3141592654e1 + 0.9e1 / 0.2e1 * 0.3141592654e1 - 0.9e1 / 0.2e1) * ParameterNS::eps * ParameterNS::eps) * pow(p[2], 0.2e1) + 0.3141592654e1 * 0.3141592654e1 * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1), 0.2e1) * ((0.3141592654e1 + 0.27e2) * pow(p[0], 0.2e1) + (0.3141592654e1 + 0.27e2) * pow(p[1], 0.2e1) + 0.4e1 * ParameterNS::eps * ParameterNS::eps * (0.3141592654e1 - 0.9e1))) * exp(-t / T) + (-0.2e1 / 0.3e1 * pow(0.3141592654e1, 0.3e1) + 0.42e2 * 0.3141592654e1 - 0.108e3) * pow(p[2], 0.6e1) + ((-0.2e1 * pow(0.3141592654e1, 0.3e1) - 0.18e2 * 0.3141592654e1 + 0.216e3) * pow(p[0], 0.2e1) + (-0.2e1 * pow(0.3141592654e1, 0.3e1) - 0.18e2 * 0.3141592654e1 + 0.216e3) * pow(p[1], 0.2e1) + 0.16e2 / 0.3e1 * pow(0.3141592654e1, 0.3e1) * ParameterNS::eps * ParameterNS::eps - 0.48e2 * 0.3141592654e1 * ParameterNS::eps * ParameterNS::eps) * pow(p[2], 0.4e1) - 0.2e1 * 0.3141592654e1 * ((0.3141592654e1 * 0.3141592654e1 + 0.30e2) * pow(p[0], 0.2e1) + (0.3141592654e1 * 0.3141592654e1 + 0.30e2) * pow(p[1], 0.2e1) - 0.16e2 / 0.3e1 * 0.3141592654e1 * 0.3141592654e1 * ParameterNS::eps * ParameterNS::eps - 0.24e2 * ParameterNS::eps * ParameterNS::eps) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * pow(p[2], 0.2e1) - 0.2e1 / 0.3e1 * pow(0.3141592654e1, 0.3e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1), 0.2e1) * (-0.8e1 * ParameterNS::eps * ParameterNS::eps + pow(p[0], 0.2e1) + pow(p[1], 0.2e1))) * p[2] * exp(-t / T) / 0.6e1 + 0.2e1 * pow(exp(-t / T), 0.2e1) * a * a * p[2] * pow(0.3141592654e1, 0.8e1) * sqrt(0.3e1) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - 0.2e1 * pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.25e2 / 0.2e1) * (a * exp(-t / T) - 0.1e1) * (a * exp(-t / T) + 0.1e1) + pow(0.3141592654e1, 0.6e1) * ((0.3141592654e1 - 0.3e1) * pow(p[2], 0.2e1) + 0.3141592654e1 * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1))) * ((0.3141592654e1 * 0.3141592654e1 - 0.12e2 * 0.3141592654e1 + 0.27e2) * pow(p[2], 0.4e1) + ((0.2e1 * 0.3141592654e1 * 0.3141592654e1 - 0.54e2) * pow(p[0], 0.2e1) + (0.2e1 * 0.3141592654e1 * 0.3141592654e1 - 0.54e2) * pow(p[1], 0.2e1) - 0.8e1 * 0.3141592654e1 * ParameterNS::eps * ParameterNS::eps * (0.3141592654e1 - 0.3e1)) * pow(p[2], 0.2e1) + 0.3141592654e1 * ((0.3141592654e1 + 0.12e2) * pow(p[0], 0.2e1) + (0.3141592654e1 + 0.12e2) * pow(p[1], 0.2e1) - 0.8e1 * ParameterNS::eps * ParameterNS::eps * (0.3141592654e1 + 0.3e1)) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1))) * sqrt(0.3e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.21e2 / 0.2e1) / 0.18e2 + pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.11e2) * (0.6e1 * exp(-t / T) * a * pow(p[2], 0.3e1) * pow(0.3141592654e1, 0.7e1) * sqrt(0.3e1) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - 0.2e1 * pow(p[2], 0.2e1)) * pow(a * exp(-t / T) - 0.1e1, 0.2e1) * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) + ((-0.30e2 * pow(p[2], 0.6e1) + 0.30e2 * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1), 0.2e1) * pow(p[2], 0.2e1)) * pow(0.3141592654e1, 0.15e2 / 0.2e1) + 0.45e2 * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - 0.3e1 / 0.5e1 * pow(p[2], 0.2e1)) * pow(p[2], 0.4e1) * pow(0.3141592654e1, 0.13e2 / 0.2e1) + pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.2e1) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - 0.3e1 * pow(p[2], 0.2e1)) * pow(0.3141592654e1, 0.17e2 / 0.2e1)) * pow(a, 0.3e1) * pow(exp(-t / T), 0.3e1) + ((0.30e2 * pow(p[2], 0.6e1) - 0.30e2 * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1), 0.2e1) * pow(p[2], 0.2e1)) * pow(0.3141592654e1, 0.15e2 / 0.2e1) - 0.135e3 * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - 0.3e1 / 0.5e1 * pow(p[2], 0.2e1)) * pow(p[2], 0.4e1) * pow(0.3141592654e1, 0.13e2 / 0.2e1) + pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.2e1) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - 0.3e1 * pow(p[2], 0.2e1)) * pow(0.3141592654e1, 0.17e2 / 0.2e1)) * a * a * pow(exp(-t / T), 0.2e1) - 0.2e1 / 0.3e1 * (0.45e2 / 0.2e1 * (-0.14e2 / 0.15e2 * pow(p[2], 0.4e1) + (pow(p[0], 0.2e1) / 0.15e2 + pow(p[1], 0.2e1) / 0.15e2 + 0.16e2 / 0.15e2 * ParameterNS::eps * ParameterNS::eps) * pow(p[2], 0.2e1) + (pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - 0.4e1 / 0.5e1 * ParameterNS::eps * ParameterNS::eps)) * pow(p[2], 0.2e1) * pow(0.3141592654e1, 0.15e2 / 0.2e1) - 0.405e3 / 0.2e1 * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - 0.3e1 / 0.5e1 * pow(p[2], 0.2e1)) * pow(p[2], 0.4e1) * pow(0.3141592654e1, 0.13e2 / 0.2e1) + (-0.5e1 / 0.2e1 * pow(p[2], 0.4e1) + (-0.3e1 / 0.2e1 * pow(p[0], 0.2e1) - 0.3e1 / 0.2e1 * pow(p[1], 0.2e1) + 0.8e1 * ParameterNS::eps * ParameterNS::eps) * pow(p[2], 0.2e1) + (pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * (-0.2e1 * ParameterNS::eps * ParameterNS::eps + pow(p[0], 0.2e1) + pow(p[1], 0.2e1))) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(0.3141592654e1, 0.17e2 / 0.2e1)) * a * exp(-t / T) + 0.15e2 * (-0.14e2 / 0.15e2 * pow(p[2], 0.4e1) + (pow(p[0], 0.2e1) / 0.15e2 + pow(p[1], 0.2e1) / 0.15e2 + 0.16e2 / 0.15e2 * ParameterNS::eps * ParameterNS::eps) * pow(p[2], 0.2e1) + (pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - 0.4e1 / 0.5e1 * ParameterNS::eps * ParameterNS::eps)) * pow(p[2], 0.2e1) * pow(0.3141592654e1, 0.15e2 / 0.2e1) - 0.45e2 * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - 0.3e1 / 0.5e1 * pow(p[2], 0.2e1)) * pow(p[2], 0.4e1) * pow(0.3141592654e1, 0.13e2 / 0.2e1) - 0.2e1 / 0.3e1 * (-0.5e1 / 0.2e1 * pow(p[2], 0.4e1) + (-0.3e1 / 0.2e1 * pow(p[0], 0.2e1) - 0.3e1 / 0.2e1 * pow(p[1], 0.2e1) + 0.8e1 * ParameterNS::eps * ParameterNS::eps) * pow(p[2], 0.2e1) + (pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * (-0.2e1 * ParameterNS::eps * ParameterNS::eps + pow(p[0], 0.2e1) + pow(p[1], 0.2e1))) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(0.3141592654e1, 0.17e2 / 0.2e1)) * a * exp(-t / T)) * sqrt(0.16e2) * pow(0.3141592654e1, -0.19e2 / 0.2e1) / ParameterNS::eps;
    return(inertia - bidiffusion);
}
static RegisterScalarFunction regsca_rhs_decay( "DecayFirstHarmonicRhs", Test6_rhs3);

// ==stationary test case "LaplaceBeltrami1"==
// Torus with R= RadTorus[0]= 1., r= RadTorus[1]= 0.6, wind == 0
// "Levelset": "Torus"

// angle from positive x-axis to (x,y,0)
double t_angle (const Point3DCL& p, double)
{
    return std::atan2( p[1], p[0]);
}

// distance from the circle in the x-y-plane around 0 with radius R
double rho (const Point3DCL& p, double)
{
    return std::sqrt( p[2]*p[2] + std::pow( std::sqrt( p[0]*p[0] + p[1]*p[1]) - RadTorus[0], 2));
}

// angle from positive (x,y,0)-direction to p.
double p_angle (const Point3DCL& p, double)
{
    return std::atan2( p[2], std::sqrt( p[0]*p[0] + p[1]*p[1]) - RadTorus[0]);
}

double laplace_beltrami_1_rhs (const Point3DCL& p, double)
{
    const double pa= p_angle( p, 0.);
    const double ta= t_angle( p, 0.);

    using std::sin;
    using std::cos;
    using std::pow;
    const double t0= (9.*sin( 3.*ta)*cos( 3.*pa + ta))/(RadTorus[1]*RadTorus[1]);
    const double t1= -(-10.*sin( 3.*ta)*cos(3.*pa + ta) - 6.*cos( 3.*ta)*sin( 3.*pa + ta))
                     /pow(RadTorus[0] + RadTorus[1]*cos( pa), 2);
    const double t2= -(3.*sin( pa)*sin( 3.*ta)*sin( 3.*pa + ta))/(RadTorus[1]*(RadTorus[0] + RadTorus[1]*cos( pa)));
    return t0 + t1 + t2;
}
static RegisterScalarFunction regsca_laplace_beltrami_1_rhs( "LaplaceBeltrami1Rhs", laplace_beltrami_1_rhs);

double laplace_beltrami_1_sol (const Point3DCL& p, double)
{
    return std::sin(3.*t_angle( p, 0.))*std::cos( 3.*p_angle( p, 0.) + t_angle( p, 0.));
}
static RegisterScalarFunction regsca_laplace_beltrami_1_sol( "LaplaceBeltrami1Sol", laplace_beltrami_1_sol);

template<class DiscP1FunType>
double L2_error (const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
    const DiscP1FunType& discsol, DROPS::InstatScalarFunction extsol)
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
double H1_error_sq (const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
    const DiscP1FunType& discsol, DROPS::InstatScalarFunction extsol)
{
    IdxDescCL* idx= const_cast<IdxDescCL*>( discsol.GetSolution()->RowIdx);
    MatDescCL A( idx, idx);
    SetupLBP1( discsol.GetMG(), &A, ls, lsbnd, 1.);
    VecDescCL sol_vec( idx);
    P1Init (extsol, sol_vec, discsol.GetMG(), discsol.GetTime());
    // sol_vec.t= discsol.GetTime();
    // SetupInterfaceRhsP1( discsol.GetMG(), &sol_vec, ls, lsbnd, extsol);
    const VectorCL diff( discsol.GetSolution()->Data - sol_vec.Data);
    return dot(A.Data*diff, diff);
}

template<class DiscP1FunType>
double H1_error (const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
                 const DiscP1FunType& discsol, DROPS::InstatScalarFunction extsol) {
    return sqrt(H1_error_sq(ls, lsbnd, discsol, extsol));
}

double L2_norm (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
    DROPS::InstatScalarFunction extsol)
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

void LinearLSInit (const DROPS::MultiGridCL& mg, DROPS::VecDescCL& ls, DROPS::InstatScalarFunction d, double t= 0.)
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

void InitVel ( const MultiGridCL& mg, VecDescCL* vec, BndDataCL<Point3DCL>& Bnd, InstatVectorFunction LsgVel, double t= 0.)
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


CahnHilliardP1BaseCL* make_cahnhilliard_timedisc( MultiGridCL& mg, LevelsetP2CL& lset,
                                              VecDescCL& v, const BndDataCL<Point3DCL>& Bnd_v,
                                              const ParamCL& P, const double & dist=0)
{
    CahnHilliardP1BaseCL* ret= 0;
    const std::string method= P.get<std::string>( "SurfSeparation.Method");

    //std::valarray<double> id(1.,v.Data.size()/3);
    //MatrixCL Ident;

    if (method == std::string( "cGcG"))
        ret= new CahnHilliardcGP1CL( mg,
            P.get<double>("SurfSeparation.Theta"), P.get<double>("SurfSeparation.Mobility"),
                    P.get<double>("SurfSeparation.Epsilon"),
            &v, Bnd_v, lset.Phi, lset.GetBndData(), the_normal_fun, dist,
            P.get<DROPS::Point3DCL>("Mesh.E1")[0]/(P.get<double>("Mesh.N1")), 1.,//S=beta_s,
            P.get<int>("SurfSeparation.Solver.Iter"), P.get<double>("SurfSeparation.Solver.Tol"),
            P.get<int>("SurfSeparation.Solver.PcAIter"), P.get<double>("SurfSeparation.Solver.PcATol"),
            P.get<int>("SurfSeparation.Solver.PcBIter"), P.get<double>("SurfSeparation.Solver.PcBTol"),
            P.get<double>("SurfSeparation.XFEMReduced"));
    else if (method == std::string( "NarrowBandStabilization"))
        ret= new CahnHilliardNarrowBandStblP1CL( mg,
                                     P.get<double>("SurfSeparation.Theta"), P.get<double>("SurfSeparation.Mobility"),
                                     P.get<double>("SurfSeparation.Epsilon"),
                                     &v, Bnd_v, lset.Phi, lset.GetBndData(), the_normal_fun, dist,
                                     1./(P.get<DROPS::Point3DCL>("Mesh.E1")[0]/(P.get<double>("Mesh.N1"))),
                                     1.,//S=beta_s,
                                     P.get<int>("SurfSeparation.Solver.Iter"), P.get<double>("SurfSeparation.Solver.Tol"),
                                     P.get<int>("SurfSeparation.Solver.PcAIter"), P.get<double>("SurfSeparation.Solver.PcATol"),
                                     P.get<int>("SurfSeparation.Solver.PcBIter"), P.get<double>("SurfSeparation.Solver.PcBTol"),
                                     P.get<double>("SurfSeparation.XFEMReduced"));
    else
        throw DROPSErrCL( std::string( "make_cahnhilliard_timedisc: Unknown method '") + method + std::string( "'.\n"));

    return ret;
}

//phase separation, Cahn-Hilliard
void Strategy (DROPS::MultiGridCL& mg, DROPS::AdapTriangCL& adap, DROPS::LevelsetP2CL& lset)
{
    using namespace DROPS;
    auto& logger = SingletonLogger::instance();
    logger.beg("read input .json");
        auto& inpJSON = P;
        // read_parameter_file_from_cmdline(inpJSON, argc, argv, "../../param/surfnavierstokes/No_Bnd_Condition.json");
        // dynamicLoad(inpJSON.get<std::string>("General.DynamicLibsPrefix"), inpJSON.get<std::vector<std::string>>("General.DynamicLibs"));
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

    //set up output
    //dirname  = P.get<std::string>("VTK.VTKDir")  + "CahnHilliard" + "_" + P.get<std::string>("SurfSeparation.Exp.Levelset")  + "/"  "l=" + std::to_string(int(P.get<int>("Mesh.AdaptRef.FinestLevel")))  + "_N=" + std::to_string(int(P.get<double>("Time.NumSteps")));

    std::cout << "dirname: " << dirname << std::endl;

    std::ofstream log_errors;
    log_errors.open( dirname +"/"
                     + "Errors.txt");
    log_errors.close();

    if (!log_errors.is_open()) {
        log_errors.open( dirname +"/"
                         + "Errors.txt");}

    std::ofstream log_global;
    log_global.open( dirname +"/"
                         + "Global.txt");
    log_global.close();

    if (!log_global.is_open()) {
            log_global.open( dirname +"/"
                             + "Global.txt");

        //std::cout<< "can't open file with computed errors in: " << dirname << std::endl;
        //return;
    }

    log_global << "Time " << "\t" << "Lyapunov_energy " << "H1_energy " << "\t" << "Number_iterations " << std::endl;

    if (P.get<std::string>("SurfSeparation.Exp.Levelset") == std::string( "AxisScalingLset"))
        dynamic_cast<DistMarkingStrategyCL*>( adap.get_marking_strategy())->SetDistFct( axis_scaling_lset_ini);
    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
    // LinearLSInit( mg, lset.Phi, the_lset_fun, 0.);
    LSInit( mg, lset.Phi, the_lset_fun, 0.);

    //DROPS::LevelsetP2CL& lset2( *LevelsetP2CL::Create( mg, lsbnd, sf, P.get_child("Levelset")) );
    //lset2.idx.CreateNumbering( mg.GetLastLevel(), mg);
    //lset2.Phi.SetIdx( &lset2.idx);
    //LSInit( mg, lset2.Phi, the_lset_fun, 0.);

    const double Vol= lset.GetVolume();
    lset.InitVolume( Vol);
    std::cout << "droplet volume: " << Vol << std::endl;

    BndDataCL<Point3DCL> Bnd_v( 6, bc_wind, bf_wind);
    IdxDescCL vidx( vecP2_FE);
    vidx.CreateNumbering( mg.GetLastLevel(), mg, Bnd_v);
    VecDescCL v( &vidx);
    InitVel( mg, &v, Bnd_v, the_wind_fun, 0.);

    //lset2.SetupSystem( make_P2Eval( mg, Bnd_v, v), P.get<double>("Time.FinalTime")/P.get<double>("Time.NumSteps"));
    double dist= 1*(2.0*P.get<DROPS::Point3DCL>("SurfSeparation.Exp.Velocity").norm()*P.get<double>("Time.FinalTime")/P.get<double>("Time.NumSteps")
                +P.get<DROPS::Point3DCL>("Mesh.E1")[0]/P.get<double>("Mesh.N1")/pow(2,P.get<int>("Mesh.AdaptRef.FinestLevel")));

    std::unique_ptr<CahnHilliardP1BaseCL> timediscp( make_cahnhilliard_timedisc( mg, lset, v, Bnd_v, P, dist));
    CahnHilliardP1BaseCL& timedisc= *timediscp;
    timedisc.SetRhs( the_rhs_fun, the_zero_fun);

    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);

    InterfaceP1RepairCL ic_repair( mg, lset.Phi, lset.GetBndData(), timedisc.ic, dist);//dist should be nonzero only for NarrowBand!
    adap.push_back( &ic_repair);

    InterfaceP1RepairCL imu_repair( mg, lset.Phi, lset.GetBndData(), timedisc.imu, dist);
    adap.push_back( &imu_repair);


    InterfaceP1RepairCL is_repair( mg, lset.Phi, lset.GetBndData(), timedisc.is, dist);
    adap.push_back( &is_repair);

    //LevelsetRepairCL lset2repair( lset2);
    //adap.push_back( &lset2repair);

    // Init concentration     // Init chemical potential
    timedisc.idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData(), dist);//WARNING
    std::cout << "Concentration NumUnknowns: " << timedisc.idx.NumUnknowns() << std::endl;
    timedisc.ic.SetIdx( &timedisc.idx);
    timedisc.imu.SetIdx( &timedisc.idx);
    timedisc.is.SetIdx( &timedisc.idx);
    timedisc.iface.SetIdx( &timedisc.idx);

    timedisc.SetInitialValue( the_poten_sol_fun, the_conc_sol_fun, the_species_sol_fun, 0.);

    BndDataCL<> nobnd( 0);
    VecDescCL the_conc_sol_vd( &lset.idx);
    LSInit( mg, the_conc_sol_vd, the_conc_sol_fun, /*t*/ 0.);
    VecDescCL the_poten_sol_vd( &lset.idx);
    LSInit( mg, the_poten_sol_vd, the_poten_sol_fun, /*t*/ 0.);
    VecDescCL the_species_sol_vd( &lset.idx);
    LSInit( mg, the_species_sol_vd, the_species_sol_fun, /*t*/ 0.);

    auto& cDOF = *(timedisc.GetConcentr().GetSolution());
    VecDescCL cDOF_ext;
    IdxDescCL P1FEidx(P1_FE);
    P1FEidx.CreateNumbering(mg.GetLastLevel(), mg);
    cDOF_ext.SetIdx(&P1FEidx);
    logger.buf << "num unknowns for bulk c: " << P1FEidx.NumUnknowns() << '\n';

    logger.log();

    VTKWriter vtkWriter(dirName + "/vtk/separation", mg, inpJSON.get<bool>("Output.Binary"));
    auto everyStep = inpJSON.get<int>("Output.EveryStep");
    auto writeVTK = [&](double t) {
        Extend(mg, cDOF, cDOF_ext);
        vtkWriter.write(t);
    };
    if (everyStep > 0) {
        vtkWriter
            .add(VTKWriter::VTKVar({"level-set", &lset.Phi.Data, VTKWriter::VTKVar::Type::P2}))
            .add(VTKWriter::VTKVar({"c_h", &cDOF_ext.Data, VTKWriter::VTKVar::Type::P1}));
        writeVTK(0.);
    }
//    if (vtkwriter.get() != 0) {
//        vtkwriter->Register( make_VTKScalar(      lset.GetSolution(),              "Levelset") );
//        vtkwriter->Register( make_VTKIfaceScalar( mg, timedisc.iface,                 "interface_mesh"));
//        vtkwriter->Register( make_VTKIfaceScalar( mg, timedisc.ic,                 "concentration"));
//        vtkwriter->Register( make_VTKIfaceScalar( mg, timedisc.imu,                 "chem_potential"));
//        vtkwriter->Register( make_VTKIfaceScalar( mg, timedisc.is,                 "species"));
//        vtkwriter->Register( make_VTKVector(      make_P2Eval( mg, Bnd_v, v),      "Velocity"));
//        //vtkwriter->Register( make_VTKScalar(      lset2.GetSolution(),             "Levelset2"));
//        vtkwriter->Register( make_VTKScalar(      make_P2Eval( mg, nobnd, the_conc_sol_vd),  "concentration_exact"));
//        vtkwriter->Register( make_VTKScalar(      make_P2Eval( mg, nobnd, the_poten_sol_vd),  "chem_potential_exact"));
//
//        vtkwriter->Write( 0.);
//    }
    //if (P.get<int>( "SurfSeparation.SolutionOutput.Freq") > 0)
    //    DROPS::WriteFEToFile( timedisc.ic, mg, P.get<std::string>( "SurfSeparation.SolutionOutput.Path"), P.get<bool>( "SolutionOutput.Binary"));

    const double dt= P.get<double>("Time.FinalTime")/P.get<double>("Time.NumSteps");
    double L_2x_chi_err= L2_error( lset.Phi, lset.GetBndData(), timedisc.GetConcentr(), the_conc_sol_fun);
    double L_2x_omega_err= L2_error( lset.Phi, lset.GetBndData(), timedisc.GetPotential(), the_poten_sol_fun);

    std::cout << "L_2x-3-error: " << L_2x_chi_err
              << "\nnorm of true solution: " << L2_norm( mg, lset.Phi, lset.GetBndData(), the_conc_sol_fun)
              << std::endl;
    double L_inftL_2x_chi_err= L_2x_chi_err;
    double L_inftL_2x_omega_err= L_2x_omega_err;

    std::cout << "L_inftL_2x-chi_error: " <<  L_inftL_2x_chi_err << std::endl;
    std::cout << "L_inftL_2x-omega_error: " <<  L_inftL_2x_omega_err << std::endl;

    double H_1x_chi_err= H1_error( lset.Phi, lset.GetBndData(), timedisc.GetConcentr(), the_conc_sol_fun);
    std::cout << "H_1x-chi_error: " << H_1x_chi_err << std::endl;
    double H_1x_omega_err= H1_error( lset.Phi, lset.GetBndData(), timedisc.GetPotential(), the_poten_sol_fun);
    std::cout << "H_1x-omega_error: " << H_1x_omega_err << std::endl;
    double L_2tH_1x_chi_err_sq= 0.5*dt*std::pow( H_1x_chi_err, 2);
    double L_2tH_1x_omega_err_sq= 0.5*dt*std::pow( H_1x_omega_err, 2);
    double L_2tL_2x_chi_err_sq= 0.5*dt*std::pow( L_2x_chi_err, 2);
    double L_2tL_2x_omega_err_sq= 0.5*dt*std::pow( L_2x_omega_err, 2);

    BndDataCL<> ifbnd( 0);
    std::cout << "initial concentration on \\Gamma: " << Integral_Gamma( mg, lset.Phi, lset.GetBndData(), make_P1Eval(  mg, ifbnd, timedisc.ic)) << '\n';
    std::cout << "initial chemical potential on \\Gamma: " << Integral_Gamma( mg, lset.Phi, lset.GetBndData(), make_P1Eval(  mg, ifbnd, timedisc.imu)) << '\n';

    dynamic_cast<DistMarkingStrategyCL*>( adap.get_marking_strategy())->SetDistFct( lset);

    //subdivide first timesteps into substeps
    int total_subs=1;
    int sub_ratio=100;
    for (int step= 1, STEP = 1; step <= P.get<int>("Time.NumSteps"); ++step) {
        if (step>total_subs) sub_ratio=1;

        for (int substep= 1; substep <= sub_ratio; ++substep) {
            std::cout << "======================================================== step " << STEP  << ":\n";
            ScopeTimerCL timer("Strategy: Time-loop");
            const double cur_time = (step-1 + (double)substep/(double)sub_ratio) * dt;
            double cur_dt=(1./(double)sub_ratio) * dt;
            // Assumes (as the rest of Drops), that the current triangulation is acceptable to perform the time-step.
            // If dt is large and AdapRef.Width is small, this may not be true.
            // Watch for large differences in numbers of old and new dof.
            timedisc.InitTimeStep();
            LSInit(mg, lset.Phi, the_lset_fun, cur_time);
            InitVel(mg, &v, Bnd_v, the_wind_fun, cur_time);

            timedisc.DoStep0(cur_time);

            VectorCL unityVector(1., timedisc.Mass.Data.num_rows());
            auto surfaceArea = dot(timedisc.Mass.Data * unityVector, unityVector);
            auto raftFraction = dot(timedisc.Mass.Data * timedisc.ic.Data, unityVector) / surfaceArea;
            // auto surfaceArea = Integral_Gamma(mg, lset.Phi, lset.GetBndData(), unityFunc);
            // auto raftFraction = Integral_Gamma(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, ifbnd, timedisc.ic)) / surfaceArea;

            std::cout << "% of raft area: " << 100. * raftFraction << '\n';
            std::cout << "chemical potential on \\Gamma: "
                      << Integral_Gamma(mg, lset.Phi, lset.GetBndData(), make_P1Eval(mg, ifbnd, timedisc.imu)) << '\n';

            L_2x_chi_err = L2_error(lset.Phi, lset.GetBndData(), timedisc.GetConcentr(), the_conc_sol_fun);
            std::cout << "L_2x-chi_error: " << L_2x_chi_err
                      << "\nnorm of concentration solution: "
                      << L2_norm(mg, lset.Phi, lset.GetBndData(), the_conc_sol_fun)
                      << std::endl;

            L_2x_omega_err = L2_error(lset.Phi, lset.GetBndData(), timedisc.GetPotential(), the_poten_sol_fun);
            std::cout << "L_2x-omega_error: " << L_2x_omega_err
                      << "\nnorm of omega solution: " << L2_norm(mg, lset.Phi, lset.GetBndData(), the_poten_sol_fun)
                      << std::endl;

            L_inftL_2x_chi_err = std::max(L_inftL_2x_chi_err, L_2x_chi_err);
            std::cout << "L_inftL_2x-chi_error: " << L_inftL_2x_chi_err << std::endl;

            L_inftL_2x_omega_err = std::max(L_inftL_2x_omega_err, L_2x_omega_err);
            std::cout << "L_inftL_2x-omega_error: " << L_inftL_2x_omega_err << std::endl;

            //L_2tH_1x_chi_err_sq+= (step > 1 ? 0.5 : 0.)*dt*std::pow( H_1x_chi_err, 2);
            H_1x_chi_err = H1_error(lset.Phi, lset.GetBndData(), timedisc.GetConcentr(), the_conc_sol_fun);
            std::cout << "H_1x-chi_error: " << H_1x_chi_err << std::endl;

            L_2tH_1x_chi_err_sq += 0.5 * cur_dt * std::pow(H_1x_chi_err, 2);
            std::cout << "L_2tH_1x-chi_error: " << std::sqrt(L_2tH_1x_chi_err_sq) << std::endl;

            //L_2tH_1x_omega_err_sq+= (step > 1 ? 0.5 : 0.)*dt*std::pow( H_1x_omega_err, 2);
            H_1x_omega_err = H1_error(lset.Phi, lset.GetBndData(), timedisc.GetPotential(), the_poten_sol_fun);
            std::cout << "H_1x-omega_error: " << H_1x_omega_err << std::endl;

            L_2tH_1x_omega_err_sq += 0.5 * cur_dt * std::pow(H_1x_omega_err, 2);
            std::cout << "L_2tH_1x-omega_error: " << std::sqrt(L_2tH_1x_omega_err_sq) << std::endl;

            L_2tL_2x_omega_err_sq += 0.5 * cur_dt * std::pow(L_2x_omega_err, 2);
            std::cout << "L_2tL_2x-omega_error: " << std::sqrt(L_2tL_2x_omega_err_sq) << std::endl;

            L_2tL_2x_chi_err_sq += 0.5 * cur_dt * std::pow(L_2x_chi_err, 2);
            std::cout << "L_2tL2x-chi_error: " << std::sqrt(L_2tL_2x_chi_err_sq) << std::endl;

            double Lyapunov_energy = (ParameterNS::eps / 2.) *
                                     H1_error(lset.Phi, lset.GetBndData(), timedisc.GetConcentr(),
                                              the_zero_fun)//dot(Laplace.Data*chi.Data,chi.Data)
                                     + (1. / ParameterNS::eps) * Integral_Gamma(mg, lset.Phi, lset.GetBndData(),
                                                                                make_P1Eval(mg, ifbnd,
                                                                                            timedisc.ienergy));

            double perimeter_estimator = ParameterNS::eps * dot(timedisc.Laplace.Data * cDOF.Data, cDOF.Data);
            log_global << cur_time << "\t" << Lyapunov_energy << "\t" << perimeter_estimator << "\t" << timedisc.GetSolver().GetIter() << "\n" << std::endl;

            if (everyStep > 0 && STEP % everyStep == 0) {
                logger.beg("write vtk");
                    writeVTK(cur_time);
                logger.end();
//                if (
//                        (((substep) % 20) == 0)
//                        //(step<100)
//                        ||
//                        (vtkwriter.get() != 0 && ((step) % P.get<int>("VTK.Freq")) == 0)) {
//                    LSInit(mg, the_conc_sol_vd, the_conc_sol_fun, /*t*/ cur_time);
//                    LSInit(mg, the_poten_sol_vd, the_poten_sol_fun, /*t*/ cur_time);
//
//                    vtkwriter->Write(cur_time);
//                }

            }
//        lset2.DoStep();
//        VectorCL rhs( lset2.Phi.Data.size());
//        lset2.ComputeRhs( rhs);
//        lset2.SetupSystem( make_P2Eval( mg, Bnd_v, v, cur_time));
//        lset2.SetTimeStep( dt);
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
            /*int freq = (int) (0.5 * (dist / P.get<DROPS::Point3DCL>("SurfSeparation.Exp.Velocity").norm()) /
                    ((1./(double)sub_ratio)*P.get<double>("Time.FinalTime") / P.get<double>("Time.NumSteps")));*/
            int freq=P.get<int>("Mesh.AdaptRef.Freq");
            const bool doGridMod = freq && step % freq == 0;
            const bool gridChanged = doGridMod ? adap.UpdateTriang() : false;

            // output
            std::ofstream stats(dirName + "/stats/t_" + std::to_string(STEP) + ".json");
            ParamCL tJSON;
            tJSON.put("t", cur_time);
            tJSON.put("dt", cur_dt);
            tJSON.put("Integral.PerimeterEstimate", perimeter_estimator);
            tJSON.put("Integral.LyapunovEnergy", Lyapunov_energy);
            tJSON.put("Integral.SurfaceArea", surfaceArea);
            tJSON.put("Integral.RaftFraction", raftFraction);
            stats << tJSON;
            logger.buf << tJSON;
            logger.log();

            if (gridChanged) {

                std::cout << "Triangulation changed.\n";
                vidx.DeleteNumbering(mg);
                vidx.CreateNumbering(mg.GetLastLevel(), mg, Bnd_v);
                v.SetIdx(&vidx);
                InitVel(mg, &v, Bnd_v, the_wind_fun, cur_time);
                LSInit(mg, lset.Phi, the_lset_fun, cur_time);
                the_conc_sol_vd.SetIdx(&lset.idx);
                LSInit(mg, the_conc_sol_vd, the_conc_sol_fun, /*t*/ cur_time);
                the_poten_sol_vd.SetIdx(&lset.idx);
                LSInit(mg, the_poten_sol_vd, the_poten_sol_fun, /*t*/ cur_time);
                the_species_sol_vd.SetIdx(&lset.idx);
                LSInit(mg, the_species_sol_vd, the_species_sol_fun, /*t*/ cur_time);
                // timedisc.Update(); // Called unconditionally in DoStep.

                //lset2.SetupSystem( make_P2Eval( mg, Bnd_v, v), dt);

                std::cout << "rel. Volume: " << lset.GetVolume() / Vol << std::endl;
                lset.AdjustVolume();
                lset.GetVolumeAdjuster()->DebugOutput(std::cout);
                //vtkwriter->Write( cur_time+dt/2.);

            }
            ++STEP;
        }
    }
    std::cout << std::endl;
    std::cout << "L2t_L2x_omega L2t_L2x_chi L2t_H1x_omega L2t_H1x_chi Linft_L2x_omega Linft_L2x_chi  H1x_omega(T) H1x_chi(T) L2x_omega(T) L2x_chi(T)\t" <<std::endl
        << std::sqrt(L_2tL_2x_omega_err_sq) <<"\t" << std::sqrt( L_2tL_2x_chi_err_sq)<<"\t"
        << std::sqrt(L_2tH_1x_omega_err_sq) <<"\t" << std::sqrt(L_2tH_1x_chi_err_sq)<<"\t"
        << L_inftL_2x_omega_err<<"\t" << L_inftL_2x_chi_err << "\t"
        << H_1x_omega_err <<"\t" <<  H_1x_chi_err <<"\t"
        << L_2x_omega_err <<"\t" << L_2x_chi_err << "\t"
        <<std::endl;
    log_errors << "L2t_L2x_omega L2t_L2x_chi L2t_H1x_omega L2t_H1x_chi Linft_L2x_omega Linft_L2x_chi  H1x_omega(T) H1x_chi(T) L2x_omega(T) L2x_chi(T)\t" <<std::endl
        << std::sqrt(L_2tL_2x_omega_err_sq) <<"\t" << std::sqrt( L_2tL_2x_chi_err_sq)<<"\t"
        << std::sqrt(L_2tH_1x_omega_err_sq) <<"\t" << std::sqrt(L_2tH_1x_chi_err_sq)<<"\t"
        << L_inftL_2x_omega_err<<"\t" << L_inftL_2x_chi_err << "\t"
        << H_1x_omega_err <<"\t" <<  H_1x_chi_err <<"\t"
        << L_2x_omega_err <<"\t" << L_2x_chi_err << "\t"
        <<std::endl;
    //delete &lset2;
    log_global.close();
        log_errors.close();

    }

void StationaryStrategyP1 (DROPS::MultiGridCL& mg, DROPS::AdapTriangCL& adap, DROPS::LevelsetP2CL& lset)
{
    adap.MakeInitialTriang();

    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
    // LinearLSInit( mg, lset.Phi, the_lset_fun);
    LSInit( mg, lset.Phi, the_lset_fun, 0.);

    DROPS::IdxDescCL ifaceidx( P1IF_FE);
    ifaceidx.GetXidx().SetBound( P.get<double>("SurfSeparation.XFEMReduced"));
    ifaceidx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
    std::cout << "NumUnknowns: " << ifaceidx.NumUnknowns() << std::endl;

    DROPS::MatDescCL M( &ifaceidx, &ifaceidx);
    DROPS::SetupInterfaceMassP1( mg, &M, lset.Phi, lset.GetBndData());
    std::cout << "M is set up.\n";
    DROPS::MatDescCL A( &ifaceidx, &ifaceidx);
    DROPS::SetupLBP1( mg, &A, lset.Phi, lset.GetBndData(), P.get<double>("SurfSeparation.Visc"));
//     DROPS::MatrixCL L;
//     L.LinComb( 1.0, A.Data, 1.0, M.Data);
    DROPS::MatrixCL& L= A.Data;
    DROPS::VecDescCL b( &ifaceidx);
    DROPS::SetupInterfaceRhsP1( mg, &b, lset.Phi, lset.GetBndData(), the_rhs_fun);

    //DROPS::WriteToFile( M.Data, "m_iface.txt", "M");
    //DROPS::WriteToFile( A.Data, "a_iface.txt", "A");
    //DROPS::WriteFEToFile( b, mg, "rhs_iface.txt", /*binary=*/ true);

    typedef DROPS::SSORPcCL SurfPcT;
    SurfPcT surfpc;
    typedef DROPS::PCGSolverCL<SurfPcT> SurfSolverT;
    SurfSolverT surfsolver( surfpc, P.get<int>("SurfSeparation.Solver.Iter"), P.get<double>("SurfSeparation.Solver.Tol"), true);

    DROPS::VecDescCL x( &ifaceidx);
    surfsolver.Solve( L, x.Data, b.Data, x.RowIdx->GetEx());
    std::cout << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';

    if (P.get<int>( "SurfSeparation.SolutionOutput.Freq") > 0)
        DROPS::WriteFEToFile( x, mg, P.get<std::string>( "SurfSeparation.SolutionOutput.Path"), P.get<bool>( "SurfSeparation.SolutionOutput.Binary"));

    DROPS::IdxDescCL ifacefullidx( DROPS::P1_FE);
    ifacefullidx.CreateNumbering( mg.GetLastLevel(), mg);
    DROPS::VecDescCL xext( &ifacefullidx);
    DROPS::Extend( mg, x, xext);
    DROPS::NoBndDataCL<> nobnd;
//    if (vtkwriter.get() != 0) {
//        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "Levelset") );
//        vtkwriter->Register( make_VTKIfaceScalar( mg, x, "InterfaceSol"));
//        vtkwriter->Write( 0.);
//    }

    double L2_err( L2_error( lset.Phi, lset.GetBndData(), make_P1Eval( mg, nobnd, xext), the_conc_sol_fun));
    std::cout << "L_2-error: " << L2_err << std::endl;
}

/// \brief Accumulate L2-norms and errors on the higher order zero level.
/// Works for P1IF_FE, P2IF_FE, and C-functions. All functions are evaluated on the P2-levelset.
class InterfaceL2AccuP2CL : public TetraAccumulatorCL
{
  private:
    const InterfaceCommonDataP2CL& cdata_;
    const MultiGridCL& mg;
    std::string name_;

    NoBndDataCL<> nobnddata;
    const VecDescCL* fvd;

    InstatScalarFunction f;
    double f_time;

    LocalLaplaceBeltramiP2CL loc_lb;

    InstatVectorFunction f_grad;
    double f_grad_time;

    InterfaceL2AccuP2CL* tid0p; // The object in OpenMP-thread 0, in which the following variables are updated.
    std::vector<double> f_grid_norm,
                        f_grid_int,
                        f_norm,
                        f_int,
                        err,
                        area,
                        f_grid_grad_norm,
                        f_grad_norm,
                        grad_err;

  public:
    double f_grid_norm_acc,
           f_grid_int_acc,
           f_norm_acc,
           f_int_acc,
           err_acc,
           area_acc,
           f_grid_grad_norm_acc,
           f_grad_norm_acc,
           grad_err_acc;

    InterfaceL2AccuP2CL (const InterfaceCommonDataP2CL& cdata, const MultiGridCL& mg_arg, std::string name= std::string())
        : cdata_( cdata), mg( mg_arg), name_( name), fvd( 0), f( 0), f_time( 0.),  loc_lb( 1.), f_grad( 0), f_grad_time( 0.){}
    virtual ~InterfaceL2AccuP2CL () {}

    void set_name (const std::string& n) { name_= n; }
    void set_grid_function (const VecDescCL& fvdarg) { fvd= &fvdarg; }
    void set_function (const InstatScalarFunction farg, double f_time_arg= 0.) {
        f= farg;
        f_time= f_time_arg;
    }
    void set_grad_function (const InstatVectorFunction farg, double f_time_arg= 0.) {
        f_grad= farg;
        f_grad_time= f_time_arg;
    }

    virtual void begin_accumulation () {
        std::cout << "InterfaceL2AccuP2CL::begin_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\".\n";

        tid0p= this;
        f_grid_norm.clear();
        f_grid_norm.resize( omp_get_max_threads(), 0.);
        f_grid_int.clear();
        f_grid_int.resize( omp_get_max_threads(), 0.);

        f_norm.clear();
        f_norm.resize( omp_get_max_threads(), 0.);
        f_int.clear();
        f_int.resize( omp_get_max_threads(), 0.);

        err.clear();
        err.resize( omp_get_max_threads(), 0.);
        area.clear();
        area.resize( omp_get_max_threads(), 0.);

        f_grid_grad_norm.clear();
        f_grid_grad_norm.resize( omp_get_max_threads(), 0.);
        f_grad_norm.clear();
        f_grad_norm.resize( omp_get_max_threads(), 0.);
        grad_err.clear();
        grad_err.resize( omp_get_max_threads(), 0.);

        f_grid_norm_acc= f_grid_int_acc= f_norm_acc= f_int_acc= err_acc= area_acc= 0.;
        f_grid_grad_norm_acc= f_grad_norm_acc= grad_err_acc= 0.;
    }

    virtual void finalize_accumulation() {
        std::cout << "InterfaceL2AccuP2CL::finalize_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\":";
        area_acc= std::accumulate( area.begin(), area.end(), 0.);
        std::cout << "\n\tarea: " << area_acc;
        if (fvd != 0) {
            f_grid_norm_acc=  std::sqrt( std::accumulate( f_grid_norm.begin(), f_grid_norm.end(), 0.));
            f_grid_int_acc= std::accumulate( f_grid_int.begin(), f_grid_int.end(), 0.);
            std::cout << "\n\t|| f_grid ||_L2: " << f_grid_norm_acc
                      << "\tintegral: " << f_grid_int_acc;
        }
        if (f != 0) {
            f_norm_acc=  std::sqrt( std::accumulate( f_norm.begin(), f_norm.end(), 0.));
            f_int_acc= std::accumulate( f_int.begin(), f_int.end(), 0.);
            std::cout << "\n\t|| f ||_L2: " << f_norm_acc
                      << "\t integral: " << f_int_acc;
        }
        if (fvd != 0 && f != 0) {
            err_acc=  std::sqrt( std::accumulate( err.begin(), err.end(), 0.));
            std::cout << "\n\t|| f - f_grid ||_L2: " << err_acc;

            const double mvf_err= std::sqrt( std::pow( err_acc, 2) - std::pow( f_grid_int_acc - f_int_acc, 2)/area_acc);
            std:: cout << "\t|| f - c_f - (f_grid -c_{f_grid}) ||_L2: " << mvf_err;
        }

        if (fvd != 0) {
            f_grid_grad_norm_acc=  std::sqrt( std::accumulate( f_grid_grad_norm.begin(), f_grid_grad_norm.end(), 0.));
            std::cout << "\n\t|| f_grid_grad ||_L2: " << f_grid_grad_norm_acc;
        }
        if (f_grad != 0) {
            f_grad_norm_acc=  std::sqrt( std::accumulate( f_grad_norm.begin(), f_grad_norm.end(), 0.));
            std::cout << "\n\t|| f_grad ||_L2: " << f_grad_norm_acc;
        }
        if (fvd != 0 && f_grad != 0) {
            grad_err_acc=  std::sqrt( std::accumulate( grad_err.begin(), grad_err.end(), 0.));
            std::cout << "\n\t|| f_grad - f_grid_grad ||_L2: " << grad_err_acc;
        }
        std::cout << std::endl;
    }

    virtual void visit (const TetraCL& t) {
        const InterfaceCommonDataP2CL& cdata= cdata_.get_clone();
        if (cdata.empty())
            return;

        const int tid= omp_get_thread_num();

        tid0p->area[tid]+= quad_2D( cdata.qdom_projected.absdets(), cdata.qdom);

        std::valarray<double> qfgrid,
                              qf;
        if (fvd != 0) {
            // XXX: Check, whether incomplete P2-Data exists locally (which is allowed for P2IF_FE, but not handled correctly by this class --> Extend fvd). Likewise for P1IF_FE...
            if (fvd->RowIdx->GetFE() == P2IF_FE)
                resize_and_evaluate_on_vertexes( make_P2Eval( mg, nobnddata, *fvd), t, cdata.qdom, qfgrid);
            else if (fvd->RowIdx->GetFE() == P1IF_FE)
                resize_and_evaluate_on_vertexes( make_P1Eval( mg, nobnddata, *fvd), t, cdata.qdom, qfgrid);
//             resize_and_evaluate_on_vertexes( make_P2Eval( mg, nobnddata, *fvd), cdata.qdom_projected, qfgrid);
            tid0p->f_grid_int[tid]+=  quad_2D( cdata.qdom_projected.absdets()*qfgrid,        cdata.qdom);
            tid0p->f_grid_norm[tid]+= quad_2D( cdata.qdom_projected.absdets()*qfgrid*qfgrid, cdata.qdom);
        }
        if (f != 0) {
            resize_and_evaluate_on_vertexes( f, cdata.qdom_projected.vertexes(), f_time, qf);
            tid0p->f_int[tid]+= quad_2D( cdata.qdom_projected.absdets()*qf, cdata.qdom);
            tid0p->f_norm[tid]+= quad_2D( cdata.qdom_projected.absdets()*qf*qf, cdata.qdom);
        }
        if (fvd != 0 && f != 0) {
            std::valarray<double> qerr= qfgrid - qf;
            tid0p->err[tid]+= quad_2D( cdata.qdom_projected.absdets()*qerr*qerr, cdata.qdom);
        }

        GridFunctionCL<Point3DCL> qfgradgrid,
                                  qfgrad;
        if (fvd != 0) {
            qfgradgrid.resize( cdata.qdom.vertex_size());
            if (fvd->RowIdx->GetFE() == P2IF_FE) {
                LocalP2CL<> lp2( t, *fvd, nobnddata);
                loc_lb.setup( t, cdata);
                for (Uint i= 0; i < 10; ++i)
                    qfgradgrid+= lp2[i]*loc_lb.get_qgradp2( i);
            }
            else if (fvd->RowIdx->GetFE() == P1IF_FE) {
// // XXX Implement this case.
//                 LocalP1CL<> lp1( t, *fvd, nobnddata);
//                 loc_lb_p1.setup( t, cdata);
//                 for (Uint i= 0; i < 4; ++i)
//                     qfgradgrid+= lp1[i]*loc_lb_p1.get_qgradp1( i);
            }
            tid0p->f_grid_grad_norm[tid]+= quad_2D( cdata.qdom_projected.absdets()*dot( qfgradgrid, qfgradgrid), cdata.qdom);
        }
        if (f_grad != 0) {
            resize_and_evaluate_on_vertexes( f_grad, cdata.qdom_projected.vertexes(), f_grad_time, qfgrad);
            for (Uint i= 0; i < cdata.qdom_projected.vertexes().size(); ++i) {
                Point3DCL n= cdata.quaqua.local_ls_grad( *cdata.qdom_projected.vertexes()[i].first, cdata.qdom_projected.vertexes()[i].second);
                n/= n.norm();
                qfgrad[i]-= inner_prod( n, qfgrad[i])*n;
            }
            tid0p->f_grad_norm[tid]+= quad_2D( cdata.qdom_projected.absdets()*dot( qfgrad, qfgrad), cdata.qdom);
        }
        if (fvd != 0 && f_grad != 0) {
            GridFunctionCL<Point3DCL> qerr( qfgradgrid - qfgrad);
            tid0p->grad_err[tid]+= quad_2D( cdata.qdom_projected.absdets()*dot( qerr, qerr), cdata.qdom);
        }
    }

    virtual InterfaceL2AccuP2CL* clone (int /*clone_id*/) { return new InterfaceL2AccuP2CL( *this); }
};


/// \brief Accumulate different error measures for the approximation of the level
/// set function $\varphi$ by the piecewise linear approximation $\varphi_h$.
class InterfaceApproxErrorDeformAccuCL : public TetraAccumulatorCL
{
  private:
    VecDescCL* yG_, ///< error-term for H^1 smooth \varphi
             * ydist_; ///< distance per QuaQuamapper

    InterfaceCommonDataDeformP2CL* cdata_;

    IdxT numry[4];
    double vec[4];

    InstatScalarFunction d_; // exact distance
    InstatVectorFunction Dd_; // gradient of exact distance
    GridFunctionCL<> qd;
    GridFunctionCL<Point3DCL> qDderr;

    OpenMPVar_MinInit_Max_CL<double> max_h,
                                     max_d,
                                     max_Dderr;

  public:
    InterfaceApproxErrorDeformAccuCL (InterfaceCommonDataDeformP2CL* cdataarg, VecDescCL* yg, VecDescCL* ydist)
        : yG_( yg), ydist_( ydist), cdata_ (cdataarg), d_ (0), Dd_ (0) {}
    virtual ~InterfaceApproxErrorDeformAccuCL () {}

    InterfaceApproxErrorDeformAccuCL& set_d  (InstatScalarFunction darg)  {  d_= darg;  return *this; }
    InterfaceApproxErrorDeformAccuCL& set_Dd (InstatVectorFunction Ddarg) { Dd_= Ddarg; return *this; }

    virtual void begin_accumulation () {
        std::cout << "#InterfaceApproxErrorDeformAccuCL::begin_accumulation"
                     ": " << ydist_->RowIdx->NumUnknowns() << " rows.\n";
        max_h.scatter ();
        max_d.scatter ();
        max_Dderr.scatter ();
    }

    virtual void finalize_accumulation() {
        std::cout << "#InterfaceApproxErrorDeformAccuCL::finalize_accumulation: ";
        max_h.reduce();
        std::cout << "\n\tmax_h: " << max_h.value() << "\n";
        max_d.reduce();
        std::cout << "\tmax_d: " << max_d.value() << "\n";
        max_Dderr.reduce();
        std::cout << "\tmax_dDerr: " << max_Dderr.value() << "\n";
    }

    virtual void visit (const TetraCL& t) {
        InterfaceCommonDataDeformP2CL& cdata= cdata_->get_clone ();
        const int tid= omp_get_thread_num();
        if (cdata.empty ())
            return;

        resize_and_evaluate_on_vertexes( d_, t, cdata.qdom2d_full, 0., qd);
        qd= std::abs (qd);
        max_d.value (tid)= std::max (max_d.value (tid), *std::max_element (&qd[0], &qd[0] + cdata.qdom2d_full.vertex_size ()));
        resize_and_evaluate_on_vertexes( Dd_, t, cdata.qdom2d_full, 0., qDderr);
        const SurfacePatchCL::FacetT& facet= cdata.surf.facet_begin()[0];
        const BaryCoordCL verts[3]= { cdata.surf.vertex_begin()[facet[0]],
                                      cdata.surf.vertex_begin()[facet[1]],
                                      cdata.surf.vertex_begin()[facet[2]] };
        cdata.Phi.set_surface_patch (verts, cdata.pos_pt);
        for (Uint i= 0; i < cdata.qdom2d_full.vertex_size (); ++i) {
            cdata.Phi.set_point (cdata.qdom2d_only_weights.vertex_begin ()[i], true);
            Point3DCL n= cdata.Phi.dPhi(cdata.qdom2d_only_weights.vertex_begin ()[i])*cdata.Phi.w;
            n/=n.norm();
            qDderr[i]-= n;
            max_Dderr.value (tid)= std::max (max_Dderr.value (tid), qDderr[i].norm ());
        }

//         for (Uint d= 0; d < 10; ++d) {
//             G+= cdata.gradp2[d]( bc)*cdata.locp2_ls[d];
//         }
        max_h.value (tid)= std::max (max_h.value (tid), ::cbrt( std::abs( cdata.det_T)));
        GetLocalNumbP1NoBnd( numry, t, *ydist_->RowIdx);
        for (int i= 0; i < 4; ++i) {
            ydist_->Data[numry[i]]= std::max (ydist_->Data[numry[i]], max_d.value (tid));
//             yG_->Data[numry[i]]= std::max( yG_->Data[numry[i]], G.norm ());
        }
    }

    virtual InterfaceApproxErrorDeformAccuCL* clone (int /*clone_id*/) {
        InterfaceApproxErrorDeformAccuCL* p= new InterfaceApproxErrorDeformAccuCL( *this);
        p->max_h.make_reference_to (max_h);
        p->max_d.make_reference_to (max_d);
        p->max_Dderr.make_reference_to (max_Dderr);
        return p;
    }
};

/// \brief Accumulate L2-norms and errors on the deformed zero level.
/// Works for P1IF_FE, P2IF_FE, and C-functions. All functions are evaluated on the deformed levelset.
class InterfaceL2AccuDeformP2CL : public TetraAccumulatorCL
{
  private:
    const InterfaceCommonDataDeformP2CL& cdata_;
    const MultiGridCL& mg;
    std::string name_;

    NoBndDataCL<> nobnddata;
    const VecDescCL* fvd;

    InstatScalarFunction f;
    double f_time;

    LocalLaplaceBeltramiDeformP2CL loc_lb;
    LocalNormalLaplaceDeformP2CL loc_ngrad;

    InstatVectorFunction f_grad;
    double f_grad_time;

    InterfaceL2AccuDeformP2CL* tid0p; // The object in OpenMP-thread 0, in which the following variables are updated.
    std::vector<double> f_grid_norm,
                        f_grid_int,
                        f_norm,
                        f_int,
                        err,
                        area,
                        f_grid_grad_norm,
                        f_grad_norm,
                        grad_err,
                        normal_grad;

  public:
    double f_grid_norm_acc,
           f_grid_int_acc,
           f_norm_acc,
           f_int_acc,
           err_acc,
           area_acc,
           f_grid_grad_norm_acc,
           f_grad_norm_acc,
           grad_err_acc,
           normal_grad_acc;

    InterfaceL2AccuDeformP2CL (const InterfaceCommonDataDeformP2CL& cdata, const MultiGridCL& mg_arg, std::string name= std::string())
        : cdata_( cdata), mg( mg_arg), name_( name), fvd( 0), f( 0), f_time( 0.),  loc_lb( 1.), loc_ngrad (1.), f_grad( 0), f_grad_time( 0.){}
    virtual ~InterfaceL2AccuDeformP2CL () {}

    void set_name (const std::string& n) { name_= n; }
    void set_grid_function (const VecDescCL& fvdarg) { fvd= &fvdarg; }
    void set_function (const InstatScalarFunction farg, double f_time_arg= 0.) {
        f= farg;
        f_time= f_time_arg;
    }
    void set_grad_function (const InstatVectorFunction farg, double f_time_arg= 0.) {
        f_grad= farg;
        f_grad_time= f_time_arg;
    }

    virtual void begin_accumulation () {
        std::cout << "InterfaceL2AccuDeformP2CL::begin_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\".\n";

        tid0p= this;
        f_grid_norm.clear();
        f_grid_norm.resize( omp_get_max_threads(), 0.);
        f_grid_int.clear();
        f_grid_int.resize( omp_get_max_threads(), 0.);

        f_norm.clear();
        f_norm.resize( omp_get_max_threads(), 0.);
        f_int.clear();
        f_int.resize( omp_get_max_threads(), 0.);

        err.clear();
        err.resize( omp_get_max_threads(), 0.);
        area.clear();
        area.resize( omp_get_max_threads(), 0.);

        f_grid_grad_norm.clear();
        f_grid_grad_norm.resize( omp_get_max_threads(), 0.);
        f_grad_norm.clear();
        f_grad_norm.resize( omp_get_max_threads(), 0.);
        grad_err.clear();
        grad_err.resize( omp_get_max_threads(), 0.);
        normal_grad.clear();
        normal_grad.resize( omp_get_max_threads(), 0.);

        f_grid_norm_acc= f_grid_int_acc= f_norm_acc= f_int_acc= err_acc= area_acc= 0.;
        f_grid_grad_norm_acc= f_grad_norm_acc= grad_err_acc= 0.;
    }

    virtual void finalize_accumulation() {
        std::cout << "InterfaceL2AccuDeformP2CL::finalize_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\":";
        area_acc= std::accumulate( area.begin(), area.end(), 0.);
        std::cout << "\n\tarea: " << area_acc;
        if (fvd != 0) {
            f_grid_norm_acc=  std::sqrt( std::accumulate( f_grid_norm.begin(), f_grid_norm.end(), 0.));
            f_grid_int_acc= std::accumulate( f_grid_int.begin(), f_grid_int.end(), 0.);
            std::cout << "\n\t|| f_grid ||_L2: " << f_grid_norm_acc
                      << "\tintegral: " << f_grid_int_acc;
        }
        if (f != 0) {
            f_norm_acc=  std::sqrt( std::accumulate( f_norm.begin(), f_norm.end(), 0.));
            f_int_acc= std::accumulate( f_int.begin(), f_int.end(), 0.);
            std::cout << "\n\t|| f ||_L2: " << f_norm_acc
                      << "\t integral: " << f_int_acc;
        }
        if (fvd != 0 && f != 0) {
            err_acc=  std::sqrt( std::accumulate( err.begin(), err.end(), 0.));
            std::cout << "\n\t|| f - f_grid ||_L2: " << err_acc;

            const double mvf_err= std::sqrt( std::pow( err_acc, 2) - std::pow( f_grid_int_acc - f_int_acc, 2)/area_acc);
            std:: cout << "\t|| f - c_f - (f_grid -c_{f_grid}) ||_L2: " << mvf_err;
        }

        if (fvd != 0) {
            f_grid_grad_norm_acc=  std::sqrt( std::accumulate( f_grid_grad_norm.begin(), f_grid_grad_norm.end(), 0.));
            std::cout << "\n\t|| f_grid_grad ||_L2: " << f_grid_grad_norm_acc;
        }
        if (f_grad != 0) {
            f_grad_norm_acc=  std::sqrt( std::accumulate( f_grad_norm.begin(), f_grad_norm.end(), 0.));
            std::cout << "\n\t|| f_grad ||_L2: " << f_grad_norm_acc;
        }
        if (fvd != 0 && f_grad != 0) {
            grad_err_acc=  std::sqrt( std::accumulate( grad_err.begin(), grad_err.end(), 0.));
            std::cout << "\n\t|| f_grad - f_grid_grad ||_L2: " << grad_err_acc;
        }
        if (fvd != 0) {
            normal_grad_acc=  std::sqrt( std::accumulate( normal_grad.begin(), normal_grad.end(), 0.));
            std::cout << "\n\t|| n^Tf_grid_grad ||_L2(vol): " << normal_grad_acc;
        }
        std::cout << std::endl;
    }

    virtual void visit (const TetraCL& t) {
        const InterfaceCommonDataDeformP2CL& cdata= cdata_.get_clone();
        if (cdata.empty())
            return;

        const int tid= omp_get_thread_num();

        std::valarray<double> ones (1., cdata.qdom2d_only_weights.vertex_size ());
        tid0p->area[tid]+= quad_2D( ones, cdata.qdom2d_only_weights);

        std::valarray<double> qfgrid,
                              qf;
        if (fvd != 0) {
            // XXX: Check, whether incomplete P2-Data exists locally (which is allowed for P2IF_FE, but not handled correctly by this class --> Extend fvd). Likewise for P1IF_FE...
            if (fvd->RowIdx->GetFE() == P2IF_FE)
                resize_and_evaluate_on_vertexes( make_P2Eval( mg, nobnddata, *fvd), t, cdata.qdom2d_only_weights, qfgrid);
            else if (fvd->RowIdx->GetFE() == P1IF_FE)
                resize_and_evaluate_on_vertexes( make_P1Eval( mg, nobnddata, *fvd), t, cdata.qdom2d_only_weights, qfgrid);
//             resize_and_evaluate_on_vertexes( make_P2Eval( mg, nobnddata, *fvd), cdata.qdom_projected, qfgrid);
            tid0p->f_grid_int[tid]+=  quad_2D( qfgrid,        cdata.qdom2d_only_weights);
            tid0p->f_grid_norm[tid]+= quad_2D( qfgrid*qfgrid, cdata.qdom2d_only_weights);
        }
        if (f != 0) {
            resize_and_evaluate_on_vertexes( f, t, cdata.qdom2d_full, f_time, qf);
            tid0p->f_int[tid]+= quad_2D( qf, cdata.qdom2d_full);
            tid0p->f_norm[tid]+= quad_2D( qf*qf, cdata.qdom2d_full);
        }
        if (fvd != 0 && f != 0) {
            std::valarray<double> qerr= qfgrid - qf;
            tid0p->err[tid]+= quad_2D( qerr*qerr, cdata.qdom2d_full);
        }

        if (fvd != 0) {
            if (fvd->RowIdx->GetFE() == P2IF_FE) {
                LocalP2CL<> lp2( t, *fvd, nobnddata);
                loc_ngrad.setup( t, cdata);
                for (Uint i= 0; i < 10; ++i)
                    for (Uint j= 0; j < 10; ++j)
                        tid0p->normal_grad[tid]+= lp2[i]*lp2[j]*loc_ngrad.coup[i][j];
            }
        }
//         GridFunctionCL<Point3DCL> qfgradgrid,
//                                   qfgrad;
//         if (fvd != 0) {
//             qfgradgrid.resize( cdata.qdom.vertex_size());
//             if (fvd->RowIdx->GetFE() == P2IF_FE) {
//                 LocalP2CL<> lp2( t, *fvd, nobnddata);
//                 loc_lb.setup( t, cdata);
//                 for (Uint i= 0; i < 10; ++i)
//                     qfgradgrid+= lp2[i]*loc_lb.get_qgradp2( i);
//             }
//             else if (fvd->RowIdx->GetFE() == P1IF_FE) {
// // // XXX Implement this case.
// //                 LocalP1CL<> lp1( t, *fvd, nobnddata);
// //                 loc_lb_p1.setup( t, cdata);
// //                 for (Uint i= 0; i < 4; ++i)
// //                     qfgradgrid+= lp1[i]*loc_lb_p1.get_qgradp1( i);
//             }
//             tid0p->f_grid_grad_norm[tid]+= quad_2D( cdata.qdom_projected.absdets()*dot( qfgradgrid, qfgradgrid), cdata.qdom);
//         }
//         if (f_grad != 0) {
//             resize_and_evaluate_on_vertexes( f_grad, cdata.qdom_projected.vertexes(), f_grad_time, qfgrad);
//             for (Uint i= 0; i < cdata.qdom_projected.vertexes().size(); ++i) {
//                 Point3DCL n= cdata.quaqua.local_ls_grad( *cdata.qdom_projected.vertexes()[i].first, cdata.qdom_projected.vertexes()[i].second);
//                 n/= n.norm();
//                 qfgrad[i]-= inner_prod( n, qfgrad[i])*n;
//             }
//             tid0p->f_grad_norm[tid]+= quad_2D( cdata.qdom_projected.absdets()*dot( qfgrad, qfgrad), cdata.qdom);
//         }
//         if (fvd != 0 && f_grad != 0) {
//             GridFunctionCL<Point3DCL> qerr( qfgradgrid - qfgrad);
//             tid0p->grad_err[tid]+= quad_2D( cdata.qdom_projected.absdets()*dot( qerr, qerr), cdata.qdom);
//         }
    }

    virtual InterfaceL2AccuDeformP2CL* clone (int /*clone_id*/) { return new InterfaceL2AccuDeformP2CL( *this); }
};

void StationaryStrategyP2 (DROPS::MultiGridCL& mg, DROPS::AdapTriangCL& adap, DROPS::LevelsetP2CL& lset)
{
    // Initialize level set and triangulation
    adap.MakeInitialTriang();
    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
    // LinearLSInit( mg, lset.Phi, &the_lset_fun);
    LSInit( mg, lset.Phi, the_lset_fun, 0.);

    // Setup an interface-P2 numbering
    DROPS::IdxDescCL ifacep2idx( P2IF_FE);
    ifacep2idx.GetXidx().SetBound( P.get<double>("SurfSeparation.XFEMReduced"));
    ifacep2idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
    std::cout << "P2-NumUnknowns: " << ifacep2idx.NumUnknowns() << std::endl;

    // Recover the gradient of the level set function
    IdxDescCL vecp2idx( vecP2_FE);
    vecp2idx.CreateNumbering( mg.GetLastLevel(), mg);
    VecDescCL lsgradrec( &vecp2idx);
    averaging_P2_gradient_recovery( mg, lset.Phi, lset.GetBndData(), lsgradrec);

    // Compute neighborhoods of the tetras at the interface
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( 1);
    TetraToTetrasT tetra_neighborhoods;
    compute_tetra_neighborhoods( mg, lset.Phi, lset.GetBndData(), lat, tetra_neighborhoods);

    QuaQuaMapperCL quaqua( mg, lset.Phi, lsgradrec, &tetra_neighborhoods,
                           P.get<int>( "LevelsetMapper.Iter"),
                           P.get<double>( "LevelsetMapper.Tol"),
                           P.get<std::string>( "LevelsetMapper.Method") == "FixedPointWithLineSearch",
                           P.get<double>( "LevelsetMapper.ArmijoConstant"));

    VecDescCL to_iface( &vecp2idx);
//     {
//         TetraAccumulatorTupleCL accus;
//         InterfaceCommonDataP2CL cdatap2( lset.Phi, lset.GetBndData(), quaqua, lat);
//         accus.push_back( &cdatap2);
//         InterfaceDebugP2CL p2debugaccu( cdatap2);
// //         p2debugaccu.store_offsets( to_iface);
//         p2debugaccu.set_true_area( 4.*M_PI*RadDrop[0]*RadDrop[0]);
//         p2debugaccu.set_ref_dp( &dp_sphere);
//         p2debugaccu.set_ref_abs_det( &abs_det_sphere);
//         accus.push_back( &p2debugaccu);
//         accumulate( accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetMatchingFunction(), ifacep2idx.GetBndInfo());
//     }

    TetraAccumulatorTupleCL accus;
    InterfaceCommonDataP2CL cdatap2( lset.Phi, lset.GetBndData(), quaqua, lat);
    accus.push_back( &cdatap2);
    DROPS::MatDescCL Mp2( &ifacep2idx, &ifacep2idx);
    InterfaceMatrixAccuCL<LocalMassP2CL, InterfaceCommonDataP2CL> accuMp2( &Mp2, LocalMassP2CL(), cdatap2, "Mp2");
    accus.push_back( &accuMp2);
    DROPS::MatDescCL Ap2( &ifacep2idx, &ifacep2idx);
    InterfaceMatrixAccuCL<LocalLaplaceBeltramiP2CL, InterfaceCommonDataP2CL> accuAp2( &Ap2, LocalLaplaceBeltramiP2CL( P.get<double>("SurfSeparation.Visc")), cdatap2, "Ap2");
    accus.push_back( &accuAp2);
    DROPS::VecDescCL bp2( &ifacep2idx);
    InterfaceVectorAccuCL<LocalVectorP2CL, InterfaceCommonDataP2CL> acculoadp2( &bp2, LocalVectorP2CL( the_rhs_fun, bp2.t), cdatap2);
    accus.push_back( &acculoadp2);

    accumulate( accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetBndInfo());

//     TetraAccumulatorTupleCL mean_accus;
//     mean_accus.push_back( &cdatap2);
//     InterfaceL2AccuP2CL L2_mean_accu( cdatap2, mg, "P2-mean");
//     L2_mean_accu.set_grid_function( bp2);
//     mean_accus.push_back( &L2_mean_accu);
//     accumulate( mean_accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetMatchingFunction(), ifacep2idx.GetBndInfo());
//     bp2.Data-= L2_mean_accu.f_grid_int_acc/L2_mean_accu.area_acc;
//     accumulate( mean_accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetMatchingFunction(), ifacep2idx.GetBndInfo());

//     VectorCL e( 1., bp2.Data.size());
//     VectorCL Ldiag( Ap2.Data.GetDiag());
//     bp2.Data-= dot( VectorCL( e/Ldiag), bp2.Data)/std::sqrt( dot( VectorCL( e/Ldiag), e));

//     DROPS::MatrixCL Lp2;
//     Lp2.LinComb( 1.0, Ap2.Data, 1.0, Mp2.Data);
    MatrixCL& Lp2= Ap2.Data;

    DROPS::WriteToFile( Ap2.Data, "ap2_iface.txt", "Ap2");
    DROPS::WriteToFile( Mp2.Data, "mp2_iface.txt", "Mp2");
    DROPS::WriteFEToFile( bp2, mg, "rhsp2_iface.txt", /*binary=*/ false);

    typedef DROPS::SSORPcCL SurfPcT;
//     typedef DROPS::JACPcCL SurfPcT;
    SurfPcT surfpc;
    typedef DROPS::PCGSolverCL<SurfPcT> SurfSolverT;
    SurfSolverT surfsolver( surfpc, P.get<int>("SurfSeparation.Solver.Iter"), P.get<double>("SurfSeparation.Solver.Tol"), true);

    DROPS::VecDescCL xp2( &ifacep2idx);
    surfsolver.Solve( Lp2, xp2.Data, bp2.Data, xp2.RowIdx->GetEx());
    std::cout << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';

    TetraAccumulatorTupleCL err_accus;
    err_accus.push_back( &cdatap2);
    InterfaceL2AccuP2CL L2_accu( cdatap2, mg, "P2-solution");
    L2_accu.set_grid_function( xp2);
    L2_accu.set_function( the_conc_sol_fun, 0.);
    L2_accu.set_grad_function( the_sol_grad_fun, 0.);
    err_accus.push_back( &L2_accu);
    accumulate( err_accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetBndInfo());

    {
        std::ofstream os( "quaqua_num_outer_iter.txt");
        for (Uint i= 0; i != quaqua.num_outer_iter.size(); ++i)
            os << i << '\t' << quaqua.num_outer_iter[i] << '\n';
        os << '\n';
        for (Uint i= 0; i != quaqua.num_inner_iter.size(); ++i)
            os << i << '\t' << quaqua.num_inner_iter[i] << '\n';
    }

    if (P.get<int>( "SurfSeparation.SolutionOutput.Freq") > 0)
        DROPS::WriteFEToFile( xp2, mg, P.get<std::string>( "SurfSeparation.SolutionOutput.Path") + "_p2", P.get<bool>( "SurfSeparation.SolutionOutput.Binary"));

    DROPS::NoBndDataCL<> nobnd;
    DROPS::NoBndDataCL<Point3DCL> nobnd_vec;
    VecDescCL the_sol_vd( &lset.idx);
    LSInit( mg, the_sol_vd, the_conc_sol_fun, /*t*/ 0.);
//    if (vtkwriter.get() != 0) {
//        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "Levelset") );
//        vtkwriter->Register( make_VTKIfaceScalar( mg, xp2, "InterfaceSolP2"));
//        vtkwriter->Register( make_VTKScalar(      make_P2Eval( mg, nobnd, the_sol_vd),  "TrueSol"));
//        vtkwriter->Register( make_VTKVector( make_P2Eval( mg, nobnd_vec, lsgradrec), "LSGradRec") );
//        vtkwriter->Register( make_VTKVector( make_P2Eval( mg, nobnd_vec, to_iface), "to_iface") );
//        vtkwriter->Write( 0.);
//    }
}

void StationaryStrategyDeformationP2 (DROPS::MultiGridCL& mg, DROPS::AdapTriangCL& adap, DROPS::LevelsetP2CL& lset)
{
    // Initialize level set and triangulation
    adap.MakeInitialTriang();
    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
    // LinearLSInit( mg, lset.Phi, &the_lset_fun);
    LSInit( mg, lset.Phi, the_lset_fun, 0.);

    // Setup an interface-P2 numbering
    DROPS::IdxDescCL ifacep2idx( P2IF_FE);
    ifacep2idx.GetXidx().SetBound( P.get<double>("SurfSeparation.XFEMReduced"));
    ifacep2idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
    std::cout << "P2-NumUnknowns: " << ifacep2idx.NumUnknowns() << std::endl;

    // Compute the mesh deformation
    IdxDescCL vecp2idx( vecP2_FE);
    vecp2idx.CreateNumbering( mg.GetLastLevel(), mg);
    IdxDescCL p2idx( P2_FE);
    p2idx.CreateNumbering( mg.GetLastLevel(), mg);
    VecDescCL deformation( &vecp2idx);
    LocalQuaMapperCL locqua (mg, lset.Phi,
        /*maxiter*/ P.get<int>( "LevelsetMapper.Iter"),
        /*tol*/ P.get<double>( "LevelsetMapper.Tol"),
        /*armijo_c*/ P.get<double>( "LevelsetMapper.ArmijoConstant"),
        /*max_damping_steps*/ P.get<Uint>( "LevelsetMapper.MaxDampingSteps"));
    locqua.set_trust_region (P.get<double>( "LevelsetMapper.TrustRegion"))
          .set_deformation_method (P.get<std::string>( "LevelsetMapper.DeformationMethod") == "map_local_level_sets" ? LocalQuaMapperCL::MAP_LOCAL_LEVEL_SETS : LocalQuaMapperCL::MAP_ZERO_LEVEL_SETS);
    LocalQuaMapperDistanceP2CL locquap2(locqua); // Provides the interface for the Oswald-projection class.
    VecDescCL locdist_vd ( &p2idx);
    OswaldProjectionP2AccuCL<LocalQuaMapperDistanceP2CL> loc_dist_accu(locquap2, locdist_vd);
        loc_dist_accu.set_level_set_function (&lset.Phi, &lset.GetBndData(), &PrincipalLatticeCL::instance (1))
                     .set_check_averaging (true);
    TetraAccumulatorTupleCL accus2;
        accus2.push_back( &loc_dist_accu);
    LocalQuaMapperDeformationP2CL locquadefp2(locqua); // Provides the interface for the Oswald-projection class.
    OswaldProjectionP2AccuCL<LocalQuaMapperDeformationP2CL> loc_def_accu(locquadefp2, deformation);
    loc_def_accu.set_level_set_function (&lset.Phi, &lset.GetBndData(), &PrincipalLatticeCL::instance (1))
                .set_check_averaging (true);
    accus2.push_back( &loc_def_accu);
    accumulate( accus2, mg, p2idx.TriangLevel(), p2idx.GetBndInfo());

    // Compute neighborhoods of the tetras at the interface
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( 1);

    VecDescCL to_iface( &vecp2idx);
//     {
//         TetraAccumulatorTupleCL accus;
//         InterfaceCommonDataP2CL cdatap2( lset.Phi, lset.GetBndData(), quaqua, lat);
//         accus.push_back( &cdatap2);
//         InterfaceDebugP2CL p2debugaccu( cdatap2);
// //         p2debugaccu.store_offsets( to_iface);
//         p2debugaccu.set_true_area( 4.*M_PI*RadDrop[0]*RadDrop[0]);
//         p2debugaccu.set_ref_dp( &dp_sphere);
//         p2debugaccu.set_ref_abs_det( &abs_det_sphere);
//         accus.push_back( &p2debugaccu);
//         accumulate( accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetMatchingFunction(), ifacep2idx.GetBndInfo());
//     }

    TetraAccumulatorTupleCL accus;
    InterfaceCommonDataDeformP2CL cdatap2( lset.Phi, lset.GetBndData(), deformation, lat);
    accus.push_back( &cdatap2);
    // Setup a P1 numbering
    DROPS::IdxDescCL p1idx( P1_FE);
    p1idx.CreateNumbering( mg.GetLastLevel(), mg);
    std::cout << "P1-NumUnknowns: " << p1idx.NumUnknowns() << std::endl;
    VecDescCL d_iface_vd (&p1idx);
    InterfaceApproxErrorDeformAccuCL ifaceerroraccu (&cdatap2, /*yg*/ 0, &d_iface_vd);
    ifaceerroraccu.set_d  (&sphere_dist)
                  .set_Dd (&d_sphere_dist);
    accus.push_back( &ifaceerroraccu);

//     accumulate( accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetMatchingFunction(), ifacep2idx.GetBndInfo());
//     exit (0);

    DROPS::MatDescCL Mp2( &ifacep2idx, &ifacep2idx);
    InterfaceMatrixAccuCL<LocalMassDeformP2CL, InterfaceCommonDataDeformP2CL> accuMp2( &Mp2, LocalMassDeformP2CL(), cdatap2, "Mp2");
    accus.push_back( &accuMp2);
    DROPS::MatDescCL Ap2( &ifacep2idx, &ifacep2idx);
    InterfaceMatrixAccuCL<LocalLaplaceBeltramiDeformP2CL, InterfaceCommonDataDeformP2CL> accuAp2( &Ap2, LocalLaplaceBeltramiDeformP2CL( P.get<double>("SurfSeparation.Visc")), cdatap2, "Ap2");
    accus.push_back( &accuAp2);
    DROPS::MatDescCL Anp2( &ifacep2idx, &ifacep2idx);
    InterfaceMatrixAccuCL<LocalNormalLaplaceDeformP2CL, InterfaceCommonDataDeformP2CL> accuAnp2( &Anp2, LocalNormalLaplaceDeformP2CL (P.get<double>("SurfSeparation.NormalLaplaceCoefficient")), cdatap2, "Anp2");
    accus.push_back( &accuAnp2);
    DROPS::VecDescCL bp2( &ifacep2idx);
    InterfaceVectorAccuCL<LocalVectorDeformP2CL, InterfaceCommonDataDeformP2CL> acculoadp2( &bp2, LocalVectorDeformP2CL( the_rhs_fun, bp2.t), cdatap2);
    accus.push_back( &acculoadp2);

    accumulate( accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetBndInfo());

    DROPS::MatrixCL Lp2;
    Lp2.LinComb (1.0, Ap2.Data, 1.0, Mp2.Data, 1.0, Anp2.Data);
//     MatrixCL& Lp2= Mp2.Data;
// 
//     DROPS::WriteToFile( Ap2.Data, "ap2_iface.txt", "Ap2");
//     DROPS::WriteToFile( Mp2.Data, "mp2_iface.txt", "Mp2");
//     DROPS::WriteToFile( Anp2.Data, "anp2_vol.txt", "Anp2");
//     DROPS::WriteFEToFile( bp2, mg, "rhsp2_iface.txt", /*binary=*/ false);

    typedef DROPS::SSORPcCL SurfPcT;
// //     typedef DROPS::JACPcCL SurfPcT;
    SurfPcT surfpc;
    typedef DROPS::PCGSolverCL<SurfPcT> SurfSolverT;
    SurfSolverT surfsolver( surfpc, P.get<int>("SurfSeparation.Solver.Iter"), P.get<double>("SurfSeparation.Solver.Tol"), true);

    DROPS::VecDescCL xp2( &ifacep2idx);
    surfsolver.Solve( Lp2, xp2.Data, bp2.Data, xp2.RowIdx->GetEx());
    std::cout << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';

    TetraAccumulatorTupleCL err_accus;
    err_accus.push_back( &cdatap2);
    InterfaceL2AccuDeformP2CL L2_accu( cdatap2, mg, "deformed P2-solution");
    L2_accu.set_grid_function( xp2);
    L2_accu.set_function( the_conc_sol_fun, 0.);
    L2_accu.set_grad_function( the_sol_grad_fun, 0.);
    err_accus.push_back( &L2_accu);
    accumulate( err_accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetBndInfo());

    if (P.get<int>( "SurfSeparation.SolutionOutput.Freq") > 0)
        DROPS::WriteFEToFile( xp2, mg, P.get<std::string>( "SurfSeparation.SolutionOutput.Path") + "_p2", P.get<bool>( "SurfSeparation.SolutionOutput.Binary"));

    DROPS::NoBndDataCL<> nobnd;
    DROPS::NoBndDataCL<Point3DCL> nobnd_vec;
    VecDescCL the_sol_vd( &lset.idx);
    LSInit( mg, the_sol_vd, the_conc_sol_fun, /*t*/ 0.);
//    if (vtkwriter.get() != 0) {
//        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "Levelset") );
//        vtkwriter->Register( make_VTKIfaceScalar( mg, xp2, "InterfaceSolP2"));
//        vtkwriter->Register( make_VTKScalar(      make_P2Eval( mg, nobnd, the_sol_vd),  "TrueSol"));
//        vtkwriter->Register( make_VTKVector( make_P2Eval( mg, nobnd_vec, deformation), "deformation") );
////         vtkwriter->Register( make_VTKVector( make_P2Eval( mg, nobnd_vec, to_iface), "to_iface") );
//        vtkwriter->Register( make_VTKScalar( make_P1Eval( mg, nobnd, d_iface_vd), "d_iface") );
//        vtkwriter->Register( make_VTKScalar( make_P2Eval( mg, nobnd, locdist_vd), "locdist") );
//        vtkwriter->Write( 0.);
//    }
}

int  main (int argc, char* argv[]) {
    srand(time(NULL));
    auto& logger = SingletonLogger::instance();
    try {
        ScopeTimerCL timer( "main");

        DROPS::read_parameter_file_from_cmdline( P, argc, argv, "../../param/surfphasesep/separation.json");
        std::cout << P << std::endl;

        DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"), P.get<std::vector<std::string> >("General.DynamicLibs") );

        // dirname  = P.get<std::string>("VTK.VTKDir")  + "CahnHilliard" + "_" + P.get<std::string>("SurfSeparation.Exp.Levelset")  + "/"   "l=" + std::to_string(int(P.get<int>("Mesh.AdaptRef.FinestLevel")))  + "_N=" + std::to_string(int(P.get<double>("Time.NumSteps")));
        dirname = P.get<std::string>("Output.Directory") + "/tmp";
        std::cout << "Setting up interface-PDE.\n";
        WindVelocity= P.get<DROPS::Point3DCL>("SurfSeparation.Exp.Velocity");
        AngularVelocity= P.get<DROPS::Point3DCL>("SurfSeparation.Exp.Angular");
        RadDrop=      P.get<DROPS::Point3DCL>("SurfSeparation.Exp.RadDrop");
        PosDrop=      P.get<DROPS::Point3DCL>("SurfSeparation.Exp.PosDrop");
        RadTorus=     P.get<DROPS::Point2DCL>("SurfSeparation.Exp.RadTorus");
        the_wind_fun= invecmap[P.get<std::string>("SurfSeparation.Exp.Wind")];
        the_normal_fun= invecmap[P.get<std::string>("SurfSeparation.Exp.Normal")];
        the_lset_fun= inscamap[P.get<std::string>("SurfSeparation.Exp.Levelset")];
        the_rhs_fun=  inscamap[P.get<std::string>("SurfSeparation.Exp.Rhs")];
        the_conc_sol_fun=  inscamap[P.get<std::string>("SurfSeparation.Exp.ConcentrationSolution")];
        auto raftRatio = P.get<double>("SurfSeparation.Exp.RaftRatio");
        the_conc_sol_fun = [=](Point3DCL const &, double) {
            auto random = (double) rand() / RAND_MAX;
            auto k = .5;
            auto ampl = .1;
            if (random < .5) return raftRatio + ampl*(2*k*random - 0.5);
            return raftRatio + ampl*(2*(1-k)*random + 2*k-1- 0.5);
        };
        logger.wrn("concentration soln set to RaftRatio");
        the_poten_sol_fun=  inscamap[P.get<std::string>("SurfSeparation.Exp.ChemicalPotentialSolution")];
            the_species_sol_fun=  inscamap[P.get<std::string>("SurfSeparation.Exp.SpeciesSolution")];

            the_zero_fun= inscamap["ZeroScalarFun"];
        double sigm = P.get<double>("SurfSeparation.Mobility");
        double eps = P.get<double>("SurfSeparation.Epsilon");
        ParameterNS::sigma = sigm;
        ParameterNS::eps = eps;
        if (P.get<std::string>("SurfSeparation.Exp.ConcentrationSolution") == "LaplaceBeltrami0Sol")
            the_sol_grad_fun=  &laplace_beltrami_0_sol_grad;
        for (Uint i= 0; i < 6; ++i)
            bf_wind[i]= the_wind_fun;

        std::cout << "Setting up domain:\n";
        std::unique_ptr<MGBuilderCL> builder( make_MGBuilder( P));
        DROPS::MultiGridCL mg( *builder);
        typedef DistMarkingStrategyCL MarkerT;
        double dist=3*(2.0*P.get<DROPS::Point3DCL>("SurfSeparation.Exp.Velocity").norm()*P.get<double>("Time.FinalTime")/P.get<double>("Time.NumSteps")
                    +P.get<DROPS::Point3DCL>("Mesh.E1")[0]/P.get<double>("Mesh.N1")/pow(2,P.get<int>("Mesh.AdaptRef.FinestLevel")));

        MarkerT marker( the_lset_fun,
                dist,
               // P.get<double>( "Mesh.AdaptRef.Width"),
                        P.get<int>( "Mesh.AdaptRef.CoarsestLevel"), P.get<int>( "Mesh.AdaptRef.FinestLevel"));

        DROPS::AdapTriangCL adap( mg, &marker);
        //Yushutin
        adap.MakeInitialTriang();
        //

        // DROPS::LevelsetP2CL lset( mg, lsbnd, sf);
        DROPS::LevelsetP2CL& lset( *LevelsetP2CL::Create( mg, lsbnd, sf, P.get_child("Levelset")) );

//        if (P.get<int>("VTK.Freq",0))
//            vtkwriter= std::unique_ptr<VTKOutCL>( new VTKOutCL(
//                adap.GetMG(),
//                "DROPS data",
//                P.get<int>("Time.NumSteps")/P.get<int>("VTK.Freq") + 1,
//                dirname,
//                P.get<std::string>("VTK.VTKName"),
//                P.get<std::string>("VTK.TimeFileName"),
//                P.get<int>("VTK.Binary"),
//                P.get<bool>("VTK.UseOnlyP1"),
//                false, /* <- P2DG */
//                -1,    /* <- level */
//                P.get<bool>("VTK.ReUseTimeFile")));

        if (P.get<bool>( "SurfSeparation.Exp.StationaryPDE")) {
    //        if (P.get<std::string>("LevelsetMapper.DeformationMethod") != "")
    //            StationaryStrategyDeformationP2( mg, adap, lset);
    //        else {
    //            if (P.get<int>( "SurfSeparation.FEDegree") == 1)
    //                StationaryStrategyP1( mg, adap, lset);
    //            else
    //                StationaryStrategyP2( mg, adap, lset);
    //        }
        }
        else
            Strategy( mg, adap, lset);

        delete &lset;
        rusage usage;
        getrusage( RUSAGE_SELF, &usage);
        std::cout << "ru_maxrss: " << usage.ru_maxrss << " kB.\n";
        std::cout<<"output: " << dirname ;

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
