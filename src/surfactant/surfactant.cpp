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
#include "levelset/levelsetmapper.h"
#include "out/output.h"
#include "out/vtkOut.h"
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
#include <misc/auto_diff.h>

using namespace DROPS;

DROPS::ParamCL P;

std::unique_ptr<VTKOutCL> vtkwriter;

DROPS::InVecMap& invecmap= DROPS::InVecMap::getInstance();
DROPS::InScaMap& inscamap= DROPS::InScaMap::getInstance();

instat_vector_fun_ptr the_wind_fun;
instat_scalar_fun_ptr the_lset_fun;
instat_scalar_fun_ptr the_rhs_fun;
instat_scalar_fun_ptr the_sol_fun;
instat_vector_fun_ptr the_sol_grad_fun;

typedef DROPS::Point3DCL (*bnd_val_fun) (const DROPS::Point3DCL&, double);

DROPS::BndCondT bc_wind[6]= { DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC };
bnd_val_fun bf_wind[6];

instat_scalar_fun_ptr sigma( 0);
SurfaceTensionCL sf( sigma, 0);
DROPS::LsetBndDataCL lsbnd( 6);


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

// ==non-stationary test case "HeatConduction"( MG_, Bnd_v_, oldv_, *v_)
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
// ==instationary test case constant velocity

DROPS::Point3DCL const_vel_wind (const DROPS::Point3DCL& , double )
{
    return MakePoint3D( 1., 0., 0.);
}
static RegisterVectorFunction regvec_const_vel_wind( "ConstVelWind", const_vel_wind);

double const_vel_lset (const Point3DCL& p, double t)
{
    return std::pow( p[0]-t, 2.) + std::pow( p[1], 2.) + std::pow( p[2], 2.) - 0.89;
}
static RegisterScalarFunction regsca_const_vel_lset( "ConstVelLset", const_vel_lset);
//double const_vel_lset_ini (const Point3DCL& p, double)
//{
//    // static const double t_end= P.get<double>( "Time.FinalTime");
//    const double //tout= t_end,
//                 tin= 0,
//                // lout= const_vel_lset( p, tout),
//                 lin= const_vel_lset( p, tin);
//    return lin;
//}

double const_vel_sol (const Point3DCL& , double )
{
    return 1.;
}
static RegisterScalarFunction regsca_const_vel_sol( "ConstVelSol", const_vel_sol);

double const_vel_rhs (const Point3DCL& , double )
{
    return 0.;
}

static RegisterScalarFunction regsca_const_vel_rhs( "ConstVelRhs", const_vel_rhs);
double const_vel_sol2 (const Point3DCL& p, double t )
{
    return std::exp(-t)*p[0]*p[1];
}
static RegisterScalarFunction regsca_const_vel_sol2( "ConstVelSol2", const_vel_sol2);

double const_vel_rhs2 (const Point3DCL& p, double t)
{
    double x1=p[0]; double x2=p[1]; double x3=p[2];
    return -(double) x2 * (std::pow(x1, 0.3e1) + (double) (-2 * t - 1) * x1 * x1 + (double) (t * t + x2 * x2 + x3 * x3 + 2 * t - 6) * x1 - (double) (t * t) - (double) (x2 * x2) - (double) (x3 * x3) + (double) (4 * t)) * std::exp(-(double) t) / ((double) (t * t) - 0.2e1 * (double) t * x1 + x1 * x1 + (double) (x2 * x2) + (double) (x3 * x3));
}

static RegisterScalarFunction regsca_const_vel_rhs2( "ConstVelRhs2", const_vel_rhs2);
double const_vel_sol3 (const Point3DCL& , double t )
{
    return t;
}
static RegisterScalarFunction regsca_const_vel_sol3( "ConstVelSol3", const_vel_sol3);

double const_vel_rhs3 (const Point3DCL& , double )
{
    return 1.;
}

static RegisterScalarFunction regsca_const_vel_rhs3( "ConstVelRhs3", const_vel_rhs3);
double const_vel_sol4 (const Point3DCL& p, double t )
{
    return p[0]*p[1]*t;
}
static RegisterScalarFunction regsca_const_vel_sol4( "ConstVelSol4", const_vel_sol4);

double const_vel_rhs4 (const Point3DCL& p, double t)
{
    double x1=p[0]; double x2=p[1]; double x3=p[2];
    return (double) x2 * (pow(t, 0.3e1) + (double) (-x1 - 4) * t * t + (double) (-x1 * x1 + x2 * x2 + x3 * x3 + 6 * x1) * t + (double) (x1 * (x1 * x1 + x2 * x2 + x3 * x3))) / (t * t - 0.2e1 * t * (double) x1 + (double) (x1 * x1) + (double) (x2 * x2) + (double) (x3 * x3));
}

static RegisterScalarFunction regsca_const_vel_rhs4( "ConstVelRhs4", const_vel_rhs4);

// instationary test case quartervelociy
DROPS::Point3DCL quarter_const_vel_wind (const DROPS::Point3DCL& , double )
{
    return MakePoint3D( 0.25, 0., 0.);
}
static RegisterVectorFunction regvec_quarter_const_vel_wind( "QuarterConstVelWind", quarter_const_vel_wind);

double quarter_const_vel_lset (const Point3DCL& p, double t)
{
    return std::pow( p[0]-0.25*t, 2.) + std::pow( p[1], 2.) + std::pow( p[2], 2.) - 0.89;
}
static RegisterScalarFunction regsca_quarter_const_vel_lset( "QuarterConstVelLset", quarter_const_vel_lset);


double quarter_const_vel_sol (const Point3DCL& p, double t)
{
    return exp(-t)*p[0]*p[1];
}
static RegisterScalarFunction regsca_quarter_const_vel_sol( "QuarterConstVelSol", quarter_const_vel_sol);

double quarter_const_vel_rhs (const Point3DCL& p, double t)
{
    double x1=p[0], x2=p[1], x3=p[2];
    return (-0.625e-1 * t * t * x1 + 0.5e0 * t * x1 * x1 - pow(x1, 0.3e1) - x2 * x2 * x1 - x1 * x3 * x3 + 0.15625e-1 * t * t - 0.125e0 * x1 * t + 0.25e0 * x1 * x1 + 0.25e0 * x2 * x2 + 0.25e0 * x3 * x3 - t + 0.6e1 * x1) * exp(-t) * x2 / (x1 * x1 - 0.5e0 * x1 * t + 0.625e-1 * t * t + x2 * x2 + x3 * x3);
}
static RegisterScalarFunction regsca_quarter_const_vel_rhs( "QuarterConstVelRhs", quarter_const_vel_rhs);


// ==instationary test case shrinking_sphere

DROPS::Point3DCL shrinking_sphere_wind (const DROPS::Point3DCL& p, double t)
{
    double x1=p[0]; double x2=p[1]; double x3=p[2];
    if (x1==0 && x2==0 && x3==0){
        return MakePoint3D(0,0,0);
    }
    else{
        return MakePoint3D(-0.3e1 / 0.4e1 * sqrt(exp(-t)) * pow(x1 * x1 + x2 * x2 + x3 * x3, -0.1e1 / 0.2e1) * x1, -0.3e1 / 0.4e1 * sqrt(exp(-t)) * pow(x1 * x1 + x2 * x2 + x3 * x3, -0.1e1 / 0.2e1) * x2, -0.3e1 / 0.4e1 * sqrt(exp(-t)) * pow(x1 * x1 + x2 * x2 + x3 * x3, -0.1e1 / 0.2e1) * x3 );
    }
}
static RegisterVectorFunction regvec_shrinking_sphere_wind( "shrinking_sphereWind", shrinking_sphere_wind);

double shrinking_sphere_lset (const Point3DCL& p, double t)
{
    double x1=p[0]; double x2=p[1]; double x3=p[2];
    return x1 * x1 + x2 * x2 + x3 * x3 - 0.9e1 / 0.4e1 / exp(t);
}
static RegisterScalarFunction regsca_shrinking_sphere_lset( "shrinking_sphereLset", shrinking_sphere_lset);

double shrinking_sphere_sol (const Point3DCL& , double )
{
    return 1.;
}
static RegisterScalarFunction regsca_shrinking_sphere_sol( "shrinking_sphereSol", shrinking_sphere_sol);

double shrinking_sphere_rhs (const Point3DCL& , double )
{
    return -1.;
}

static RegisterScalarFunction regsca_shrinking_sphere_rhs( "shrinking_sphereRhs", shrinking_sphere_rhs);


double shrinking_sphere_solPaper1 (const Point3DCL& p, double t)
{
    return (1+p[0]*p[1]*p[2])*exp(t);
}
static RegisterScalarFunction regsca_shrinking_sphere_solPaper1( "shrinking_sphereSolPaper1", shrinking_sphere_solPaper1);

double shrinking_sphere_rhsPaper1 (const Point3DCL& p, double t)
{
    return (-1.5*exp(t)+16./3.*exp(2.*t))*p[0]*p[1]*p[2];
}
static RegisterScalarFunction regsca_shrinking_sphere_rhsPaper1( "shrinking_sphereRhsPaper1", shrinking_sphere_rhsPaper1);


double shrinking_sphere_solPaper2 (const Point3DCL& , double t)
{
    return exp(t);
}
static RegisterScalarFunction regsca_shrinking_sphere_solPaper2( "shrinking_sphereSolPaper2", shrinking_sphere_solPaper2);

double shrinking_sphere_rhsPaper2 (const Point3DCL& , double )
{
    return 0;
}

static RegisterScalarFunction regsca_shrinking_sphere_rhsPaper2( "shrinking_sphereRhsPaper2", shrinking_sphere_rhsPaper2);

// ==instationary test case shrinking_sphere with constant velocity

DROPS::Point3DCL shrinking_sphere_wind2 (const DROPS::Point3DCL& p, double )
{
    double x1=p[0]; double x2=p[1]; double x3=p[2];
    if (x1==0 && x2==0 && x3==0){
        return MakePoint3D(0,0,0);
    }
    else{
        return MakePoint3D(-pow(x1 * x1 + x2 * x2 + x3 * x3, -0.1e1 / 0.2e1) * x1,-pow(x1 * x1 + x2 * x2 + x3 * x3, -0.1e1 / 0.2e1) * x2,-pow(x1 * x1 + x2 * x2 + x3 * x3, -0.1e1 / 0.2e1) * x3);
    }
}
static RegisterVectorFunction regvec_shrinking_sphere_wind2( "shrinking_sphereWind2", shrinking_sphere_wind2);

double shrinking_sphere_lset2 (const Point3DCL& p, double t)
{
    double x1=p[0]; double x2=p[1]; double x3=p[2];
    return x1 * x1 + x2 * x2 + x3 * x3 - pow(0.9e1 / 0.4e1 - t, 0.2e1);
}
static RegisterScalarFunction regsca_shrinking_sphere_lset2( "shrinking_sphereLset2", shrinking_sphere_lset2);

double shrinking_sphere_sol2 (const Point3DCL& , double )
{
    return 1.;
}
static RegisterScalarFunction regsca_shrinking_sphere_sol2( "shrinking_sphereSol2", shrinking_sphere_sol2);

double shrinking_sphere_rhs2 (const Point3DCL& , double t)
{
    return 8/(4*t-9);
}

static RegisterScalarFunction regsca_shrinking_sphere_rhs2( "shrinking_sphereRhs2", shrinking_sphere_rhs2);

// ==instationary test case rotating sphere

DROPS::Point3DCL rotating_sphere_wind (const DROPS::Point3DCL& p, double )
{
    double x1=p[0]; double x2=p[1];
        return MakePoint3D(-x2,x1,0);
}
static RegisterVectorFunction regvec_rotating_sphere_wind( "rotating_sphereWind", rotating_sphere_wind);

double rotating_sphere_lset (const Point3DCL& p, double )
{
        return std::pow( p[0], 2.) + std::pow( p[1], 2.) + std::pow( p[2], 2.) - 0.89;
}
static RegisterScalarFunction regsca_rotating_sphere_lset( "rotating_sphereLset", rotating_sphere_lset);

double rotating_sphere_sol (const Point3DCL& , double )
{
    return 1.;
}
static RegisterScalarFunction regsca_rotating_sphere_sol( "rotating_sphereSol", rotating_sphere_sol);

double rotating_sphere_rhs (const Point3DCL& , double )
{
    return 0;
}

static RegisterScalarFunction regsca_rotating_sphere_rhs( "rotating_sphereRhs", rotating_sphere_rhs);

// ==instationary test case thin_ellipsoid

DROPS::Point3DCL thin_ellipsoid_wind (const DROPS::Point3DCL& p, double t)
{
    double x1=p[0]; double x2=p[1];double x3=p[2];
        return MakePoint3D((-0.20000e5 * x2 * x2 - 0.20000e5 * x3 * x3 + 0.178e3) * x1 / (0.320000e6 * pow(t, 0.3e1) * x2 * x2 + 0.320000e6 * pow(t, 0.3e1) * x3 * x3 + 0.2400000e7 * t * t * x2 * x2 + 0.2400000e7 * t * t * x3 * x3 + 0.5980000e7 * t * x2 * x2 + 0.5980000e7 * t * x3 * x3 + 0.4950000e7 * x2 * x2 + 0.4950000e7 * x3 * x3 + 0.178e3 * t + 0.445e3)
,-0.8e1 / (0.160000e6 * t * t * x2 * x2 + 0.160000e6 * t * t * x3 * x3 + 0.800000e6 * t * x2 * x2 + 0.800000e6 * t * x3 * x3 + 0.990000e6 * x2 * x2 + 0.990000e6 * x3 * x3 + 0.89e2) * (0.20000e5 * t * x2 * x2 + 0.20000e5 * t * x3 * x3 + 0.50000e5 * x2 * x2 + 0.50000e5 * x3 * x3 - 0.178e3 * t - 0.445e3) * x2
,-0.8e1 / (0.160000e6 * t * t * x2 * x2 + 0.160000e6 * t * t * x3 * x3 + 0.800000e6 * t * x2 * x2 + 0.800000e6 * t * x3 * x3 + 0.990000e6 * x2 * x2 + 0.990000e6 * x3 * x3 + 0.89e2) * (0.20000e5 * t * x2 * x2 + 0.20000e5 * t * x3 * x3 + 0.50000e5 * x2 * x2 + 0.50000e5 * x3 * x3 - 0.178e3 * t - 0.445e3) * x3
);
}
static RegisterVectorFunction regvec_thin_ellipsoid_wind( "thin_ellipsoidWind", thin_ellipsoid_wind);

double thin_ellipsoid_lset (const Point3DCL& p, double t)
{
    double x1=p[0]; double x2=p[1];double x3=p[2];
        return x1 * x1 * pow(0.10e1 + 0.4e0 * t, -0.2e1) + 0.1e3 * x2 * x2 + 0.1e3 * x3 * x3 - 0.89e0;;
}
static RegisterScalarFunction regsca_thin_ellipsoid_lset( "thin_ellipsoidLset", thin_ellipsoid_lset);

double thin_ellipsoid_sol (const Point3DCL& , double )
{
    return 1.;
}
static RegisterScalarFunction regsca_thin_ellipsoid_sol( "thin_ellipsoidSol", thin_ellipsoid_sol);

double thin_ellipsoid_rhs (const Point3DCL& p, double t)
{
    double x1=p[0]; double x2=p[1];double x3=p[2];
    return ((-0.9266896436e28 * pow(x3, 0.10e2) + (-0.7549475917e24 * x1 * x1 + 0.8156312605e26 - 0.4633448217e29 * x2 * x2) * pow(x3, 0.8e1) + (-0.9266896434e29 * pow(x2, 0.4e1) + (-0.3019790366e25 * x1 * x1 + 0.3262525043e27) * x2 * x2 + x1 * x1 * 0.666174561400000086e22 - 0.7137290240e19 * pow(x1, 0.4e1) + 0.809888462300000236e22 + pow(x1, 0.3e1) * (-0.253567486652173102e5) + x1 * (-0.307869141455447716e5)) * pow(x3, 0.6e1) + (-0.9266896434e29 * pow(x2, 0.6e1) + (-0.4529685549e25 * x1 * x1 + 0.4893787561e27) * pow(x2, 0.4e1) + (-0.2141187072e20 * pow(x1, 0.4e1) + x1 * x1 * 0.199852368399999959e23 + 0.242966538799999940e23) * x2 * x2 + 0.6324865741e17 * pow(x1, 0.4e1) + x1 * x1 * 0.509010506700000000e18 + 0.179425326000000032e18) * pow(x3, 0.4e1) + (-0.4633448217e29 * pow(x2, 0.8e1) + (-0.3019790367e25 * x1 * x1 + 0.3262525043e27) * pow(x2, 0.6e1) + (-0.2141187072e20 * pow(x1, 0.4e1) + x1 * x1 * 0.199852368399999959e23 + 0.242966538799999940e23) * pow(x2, 0.4e1) + (0.1264973148e18 * pow(x1, 0.4e1) + x1 * x1 * 0.1018021013000000000e19 + 0.358850651900000000e18) * x2 * x2 + 0.7585240859e13 * x1 * x1 + 0.2430897869e13 * pow(x1, 0.4e1)) * x3 * x3 - 0.9266896438e28 * pow(x2, 0.10e2) + (-0.7549475917e24 * x1 * x1 + 0.8156312605e26) * pow(x2, 0.8e1) + (x1 * x1 * 0.666174561400000086e22 - 0.7137290240e19 * pow(x1, 0.4e1) + 0.809888462300000236e22 + pow(x1, 0.3e1) * (-0.253567486652173102e5) + x1 * (-0.307869141455447716e5)) * pow(x2, 0.6e1) + (x1 * x1 * 0.509010506600000000e18 + 0.6324865741e17 * pow(x1, 0.4e1) + 0.179425325900000000e18) * pow(x2, 0.4e1) + (0.7585240859e13 * x1 * x1 + 0.2430897869e13 * pow(x1, 0.4e1)) * x2 * x2 + 0.721888256e7 * pow(x1, 0.4e1)) * pow(t, 0.4e1) + (-0.1677298689e27 * pow(x3, 0.10e2) + (-0.4151312383e21 * x1 * x1 + 0.1489059170e25 - 0.8386493445e27 * x2 * x2) * pow(x3, 0.8e1) + (-0.1677298688e28 * pow(x2, 0.4e1) + (-0.1660524954e22 * x1 * x1 + 0.5956236680e25) * x2 * x2 + 0.3693548142e19 * x1 * x1 + 0.3324533647e20) * pow(x3, 0.6e1) + (-0.1677298687e28 * pow(x2, 0.6e1) + (-0.2490787430e22 * x1 * x1 + 0.8934355018e25) * pow(x2, 0.4e1) + (0.1108064442e20 * x1 * x1 + 0.9973600938e20) * x2 * x2 + 0.9966924595e13 * x1 * x1 + 0.9759097608e14) * pow(x3, 0.4e1) + (-0.8386493452e27 * pow(x2, 0.8e1) + (-0.1660524954e22 * x1 * x1 + 0.5956236678e25) * pow(x2, 0.6e1) + (0.1108064442e20 * x1 * x1 + 0.9973600938e20) * pow(x2, 0.4e1) + (0.1993384919e14 * x1 * x1 + 0.1951819522e15) * x2 * x2) * x3 * x3 - 0.1677298689e27 * pow(x2, 0.10e2) + (-0.4151312383e21 * x1 * x1 + 0.1489059170e25) * pow(x2, 0.8e1) + (0.3693548142e19 * x1 * x1 + 0.3324533646e20) * pow(x2, 0.6e1) + (0.9966924595e13 * x1 * x1 + 0.9759097606e14) * pow(x2, 0.4e1)) * pow(t, 0.10e2) + (-0.3810670240e28 * pow(x3, 0.10e2) + (-0.6264090080e24 * x1 * x1 + 0.3339665478e26 - 0.1905335120e29 * x2 * x2) * pow(x3, 0.8e1) + (-0.3810670240e29 * pow(x2, 0.4e1) + (-0.2505636032e25 * x1 * x1 + 0.1335866191e27) * x2 * x2 + x1 * x1 * 0.549846400199999934e22 - 0.1772838400e20 * pow(x1, 0.4e1) + 0.459611906200000004e22) * pow(x3, 0.6e1) + (-0.3810670240e29 * pow(x2, 0.6e1) + (-0.3758454048e25 * x1 * x1 + 0.2003799288e27) * pow(x2, 0.4e1) + (-0.5318515200e20 * pow(x1, 0.4e1) + x1 * x1 * 0.164953920299999953e23 + pow(x1, 0.3e1) * (-0.188951617019483820e6) + x1 * (-0.157519645478514052e6) + 0.137883571799999990e23) * x2 * x2 + 0.1560806528e18 * pow(x1, 0.4e1) + x1 * x1 * 0.679394942599999872e18 + 0.149904608199999968e18) * pow(x3, 0.4e1) + (-0.1905335120e29 * pow(x2, 0.8e1) + (-0.2505636032e25 * x1 * x1 + 0.1335866191e27) * pow(x2, 0.6e1) + (-0.5318515200e20 * pow(x1, 0.4e1) + x1 * x1 * 0.164953920299999953e23 + pow(x1, 0.3e1) * (-0.188951617019483820e6) + x1 * (-0.157519645478514052e6) + 0.137883571799999990e23) * pow(x2, 0.4e1) + (0.3121613056e18 * pow(x1, 0.4e1) + x1 * x1 * 0.1358789885000000512e19 + 0.299809216300000192e18) * x2 * x2 + 0.1898340523e14 * x1 * x1 + 0.1511707008e14 * pow(x1, 0.4e1)) * x3 * x3 - 0.3810670240e28 * pow(x2, 0.10e2) + (-0.6264090080e24 * x1 * x1 + 0.3339665478e26) * pow(x2, 0.8e1) + (x1 * x1 * 0.549846400199999934e22 - 0.1772838400e20 * pow(x1, 0.4e1) + 0.459611906200000004e22) * pow(x2, 0.6e1) + (0.1560806528e18 * pow(x1, 0.4e1) + x1 * x1 * 0.679394942599999872e18 + 0.149904608199999968e18) * pow(x2, 0.4e1) + (0.1898340523e14 * x1 * x1 + 0.1511707008e14 * pow(x1, 0.4e1)) * x2 * x2 + 0.270708096e9 * pow(x1, 0.4e1)) * t * t + (-0.6103852834e25 * pow(x3, 0.10e2) + (-0.1006632960e19 * x1 * x1 + 0.5426993876e23 - 0.3051926417e26 * x2 * x2) * pow(x3, 0.8e1) + (-0.6103852836e26 * pow(x2, 0.4e1) + (-0.4026531840e19 * x1 * x1 + 0.2170797552e24) * x2 * x2 + 0.8959033344e16 * x1 * x1 + 0.4837014951e18) * pow(x3, 0.6e1) + (0.2365483439e12 - 0.6103852836e26 * pow(x2, 0.6e1) + (-0.6039797760e19 * x1 * x1 + 0.3256196328e24) * pow(x2, 0.4e1) + (0.2687710003e17 * x1 * x1 + 0.1451104484e19) * x2 * x2) * pow(x3, 0.4e1) + (-0.3051926418e26 * pow(x2, 0.8e1) + (-0.4026531840e19 * x1 * x1 + 0.2170797552e24) * pow(x2, 0.6e1) + (0.2687710003e17 * x1 * x1 + 0.1451104484e19) * pow(x2, 0.4e1) + 0.4730966878e12 * x2 * x2) * x3 * x3 - 0.6103852834e25 * pow(x2, 0.10e2) + (-0.1006632960e19 * x1 * x1 + 0.5426993878e23) * pow(x2, 0.8e1) + (0.8959033344e16 * x1 * x1 + 0.4837014952e18) * pow(x2, 0.6e1) + 0.2365483439e12 * pow(x2, 0.4e1)) * pow(t, 0.12e2) + (-0.1682954066e28 * pow(x3, 0.10e2) + (-0.1943470211e23 * x1 * x1 + 0.1490822886e26 - 0.8414770331e28 * x2 * x2) * pow(x3, 0.8e1) + (-0.1682954065e29 * pow(x2, 0.4e1) + (-0.7773880843e23 * x1 * x1 + 0.5963291543e26) * x2 * x2 + x1 * x1 * 0.172653929499999994e21 - 0.2621440000e16 * pow(x1, 0.4e1) + 0.623040398500000039e21 + pow(x1, 0.3e1) * (-0.149011611938476563e3) + x1 * (-0.537695076907534258e3)) * pow(x3, 0.6e1) + (-0.1682954067e29 * pow(x2, 0.6e1) + (-0.1166082126e24 * x1 * x1 + 0.8944937315e26) * pow(x2, 0.4e1) + (-0.7864320000e16 * pow(x1, 0.4e1) + x1 * x1 * 0.517961788200000029e21 + 0.186912119799999994e22) * x2 * x2 + 0.2333081600e14 * pow(x1, 0.4e1) + x1 * x1 * 0.280247078700000050e16 + 0.4578042065000000e16) * pow(x3, 0.4e1) + (-0.8414770323e28 * pow(x2, 0.8e1) + (-0.7773880842e23 * x1 * x1 + 0.5963291543e26) * pow(x2, 0.6e1) + (-0.7864320000e16 * pow(x1, 0.4e1) + x1 * x1 * 0.517961788200000029e21 + 0.186912119799999994e22) * pow(x2, 0.4e1) + (0.4666163200e14 * pow(x1, 0.4e1) + x1 * x1 * 0.5604941576000000e16 + 0.9156084129000002e16) * x2 * x2 + 0.2772050903e10 * x1 * x1) * x3 * x3 - 0.1682954065e28 * pow(x2, 0.10e2) + (-0.1943470211e23 * x1 * x1 + 0.1490822887e26) * pow(x2, 0.8e1) + (x1 * x1 * 0.172653929499999961e21 - 0.2621440000e16 * pow(x1, 0.4e1) + 0.623040398299999896e21) * pow(x2, 0.6e1) + (0.2333081600e14 * pow(x1, 0.4e1) + x1 * x1 * 0.280247078700000050e16 + 0.4578042065000000e16) * pow(x2, 0.4e1) + 0.2772050903e10 * x2 * x2 * x1 * x1) * pow(t, 0.8e1) + (-0.8904501355e28 * pow(x3, 0.10e2) + (-0.4840254681e24 * x1 * x1 + 0.7852025415e26 - 0.4452250679e29 * x2 * x2) * pow(x3, 0.8e1) + (-0.8904501356e29 * pow(x2, 0.4e1) + (-0.1936101873e25 * x1 * x1 + 0.3140810167e27) * x2 * x2 + x1 * x1 * 0.428030799200000043e22 - 0.2288844800e19 * pow(x1, 0.4e1) + 0.648240145400000204e22) * pow(x3, 0.6e1) + (-0.8904501356e29 * pow(x2, 0.6e1) + (-0.2904152810e25 * x1 * x1 + 0.4711215251e27) * pow(x2, 0.4e1) + (-0.6866534400e19 * pow(x1, 0.4e1) + x1 * x1 * 0.128409239800000010e23 + pow(x1, 0.3e1) * (-0.487896613776683807e5) + x1 * (-0.738307741518386174e5) + 0.194472043700000045e23) * x2 * x2 + 0.2032697344e17 * pow(x1, 0.4e1) + x1 * x1 * 0.244643582499999936e18 + 0.114710931399999984e18) * pow(x3, 0.4e1) + (-0.4452250679e29 * pow(x2, 0.8e1) + (-0.1936101872e25 * x1 * x1 + 0.3140810168e27) * pow(x2, 0.6e1) + (-0.6866534400e19 * pow(x1, 0.4e1) + x1 * x1 * 0.128409239800000010e23 + pow(x1, 0.3e1) * (-0.487896613776683807e5) + x1 * (-0.738307741518386174e5) + 0.194472043700000045e23) * pow(x2, 0.4e1) + (0.4065394688e17 * pow(x1, 0.4e1) + x1 * x1 * 0.489287164999999936e18 + 0.229421863000000000e18) * x2 * x2 + 0.2426410807e13 * x1 * x1 + 0.3893329920e12 * pow(x1, 0.4e1)) * x3 * x3 - 0.8904501361e28 * pow(x2, 0.10e2) + (-0.4840254681e24 * x1 * x1 + 0.7852025415e26) * pow(x2, 0.8e1) + (x1 * x1 * 0.428030799200000043e22 - 0.2288844800e19 * pow(x1, 0.4e1) + 0.648240145400000204e22) * pow(x2, 0.6e1) + (0.2032697344e17 * pow(x1, 0.4e1) + x1 * x1 * 0.244643582499999936e18 + 0.114710931399999984e18) * pow(x2, 0.4e1) + (0.2426410807e13 * x1 * x1 + 0.3893329920e12 * pow(x1, 0.4e1)) * x2 * x2) * pow(t, 0.5e1) + (-0.3661089514e26 * pow(x3, 0.10e2) + (-0.3019898880e20 * x1 * x1 + 0.3252934523e24 - 0.1830544758e27 * x2 * x2) * pow(x3, 0.8e1) + (-0.3661089516e27 * pow(x2, 0.4e1) + (-0.1207959552e21 * x1 * x1 + 0.1301173809e25) * x2 * x2 + 0.2687710003e18 * x1 * x1 + 0.4836483383e19) * pow(x3, 0.6e1) + (0.7096450306e13 - 0.3661089516e27 * pow(x2, 0.6e1) + (-0.1811939328e21 * x1 * x1 + 0.1951760714e25) * pow(x2, 0.4e1) + (0.8063130008e18 * x1 * x1 + 0.1450945014e20) * x2 * x2) * pow(x3, 0.4e1) + (-0.1830544759e27 * pow(x2, 0.8e1) + (-0.1207959552e21 * x1 * x1 + 0.1301173809e25) * pow(x2, 0.6e1) + (0.8063130008e18 * x1 * x1 + 0.1450945014e20) * pow(x2, 0.4e1) + 0.1419290063e14 * x2 * x2) * x3 * x3 - 0.3661089514e26 * pow(x2, 0.10e2) + (-0.3019898880e20 * x1 * x1 + 0.3252934523e24) * pow(x2, 0.8e1) + (0.2687710003e18 * x1 * x1 + 0.4836483384e19) * pow(x2, 0.6e1) + 0.7096450306e13 * pow(x2, 0.4e1)) * pow(t, 0.11e2) + (-0.3435973836e22 * pow(x3, 0.10e2) + (-0.1717986918e23 * x2 * x2 + 0.3058016715e20) * pow(x3, 0.8e1) + (0.1223206686e21 * x2 * x2 - 0.3435973836e23 * pow(x2, 0.4e1)) * pow(x3, 0.6e1) + (0.1834810029e21 * pow(x2, 0.4e1) - 0.3435973836e23 * pow(x2, 0.6e1)) * pow(x3, 0.4e1) + (-0.1717986918e23 * pow(x2, 0.8e1) + 0.1223206686e21 * pow(x2, 0.6e1)) * x3 * x3 + 0.3058016715e20 * pow(x2, 0.8e1) - 0.3435973836e22 * pow(x2, 0.10e2)) * pow(t, 0.15e2) + (-0.1268704800e28 * pow(x3, 0.10e2) + (-0.2840064480e24 * x1 * x1 + 0.1109212696e26 - 0.6343524000e28 * x2 * x2) * pow(x3, 0.8e1) + (-0.1268704800e29 * pow(x2, 0.4e1) + (-0.1136025792e25 * x1 * x1 + 0.4436850784e26) * x2 * x2 + x1 * x1 * 0.248516608200000039e22 - 0.1260864000e20 * pow(x1, 0.4e1) + 0.176650973599999957e22) * pow(x3, 0.6e1) + (-0.1268704800e29 * pow(x2, 0.6e1) + (-0.1704038688e25 * x1 * x1 + 0.6655276176e26) * pow(x2, 0.4e1) + (-0.3782592000e20 * pow(x1, 0.4e1) + x1 * x1 * 0.745549824500000044e22 + 0.529952920800000161e22) * x2 * x2 + 0.1105194880e18 * pow(x1, 0.4e1) + x1 * x1 * 0.376647985800000000e18 + 0.68240999199999992e17) * pow(x3, 0.4e1) + (-0.6343524000e28 * pow(x2, 0.8e1) + (-0.1136025792e25 * x1 * x1 + 0.4436850784e26) * pow(x2, 0.6e1) + (-0.3782592000e20 * pow(x1, 0.4e1) + x1 * x1 * 0.745549824500000044e22 + 0.529952920800000161e22) * pow(x2, 0.4e1) + (0.2210389760e18 * pow(x1, 0.4e1) + x1 * x1 * 0.753295971500000128e18 + 0.136481998400000000e18) * x2 * x2 + 0.1356924331e14 * x1 * x1 + 0.1505623680e14 * pow(x1, 0.4e1)) * x3 * x3 - 0.1268704800e28 * pow(x2, 0.10e2) + (-0.2840064480e24 * x1 * x1 + 0.1109212696e26) * pow(x2, 0.8e1) + (x1 * x1 * 0.248516608200000039e22 - 0.1260864000e20 * pow(x1, 0.4e1) + 0.176650973599999957e22) * pow(x2, 0.6e1) + (0.1105194880e18 * pow(x1, 0.4e1) + x1 * x1 * 0.376647985800000000e18 + 0.68240999199999992e17) * pow(x2, 0.4e1) + (0.1356924331e14 * x1 * x1 + 0.1505623680e14 * pow(x1, 0.4e1)) * x2 * x2 + 0.451180160e9 * pow(x1, 0.4e1)) * t + (-0.7514313718e24 * pow(x3, 0.10e2) + (-0.3757156860e25 * x2 * x2 + 0.668439450299999768e22) * pow(x3, 0.8e1) + (-0.7514313719e25 * pow(x2, 0.4e1) + x2 * x2 * 0.267375780300000015e23 + 0.29767881449999992e17) * pow(x3, 0.6e1) + (x2 * x2 * 0.89303644330000000e17 - 0.7514313719e25 * pow(x2, 0.6e1) + pow(x2, 0.4e1) * 0.401063670600000057e23) * pow(x3, 0.4e1) + (-0.3757156860e25 * pow(x2, 0.8e1) + pow(x2, 0.4e1) * 0.89303644330000016e17 + pow(x2, 0.6e1) * 0.267375780300000057e23) * x3 * x3 + 0.6684394505e22 * pow(x2, 0.8e1) - 0.7514313718e24 * pow(x2, 0.10e2) + pow(x2, 0.6e1) * 0.29767881450000004e17) * pow(t, 0.13e2) + (-0.6535477504e28 * pow(x3, 0.10e2) + (-0.2262223093e24 * x1 * x1 + 0.5772786258e26 - 0.3267738752e29 * x2 * x2) * pow(x3, 0.8e1) + (-0.6535477503e29 * pow(x2, 0.4e1) + (-0.9048892377e24 * x1 * x1 + 0.2309114504e27) * x2 * x2 + x1 * x1 * 0.200420021400000030e22 - 0.4584243200e18 * pow(x1, 0.4e1) + 0.389118664000000072e22 + pow(x1, 0.3e1) * (-0.651460140943527222e4) + x1 * (-0.126425904944297981e5)) * pow(x3, 0.6e1) + (-0.6535477502e29 * pow(x2, 0.6e1) + (-0.1357333857e25 * x1 * x1 + 0.3463671755e27) * pow(x2, 0.4e1) + (-0.1375272960e19 * pow(x1, 0.4e1) + x1 * x1 * 0.601260064099999941e22 + 0.116735599299999983e23) * x2 * x2 + 0.4077060096e16 * pow(x1, 0.4e1) + x1 * x1 * 0.81632726719999952e17 + 0.53483257099999976e17) * pow(x3, 0.4e1) + (-0.3267738752e29 * pow(x2, 0.8e1) + (-0.9048892377e24 * x1 * x1 + 0.2309114504e27) * pow(x2, 0.6e1) + (-0.1375272960e19 * pow(x1, 0.4e1) + x1 * x1 * 0.601260064099999941e22 + 0.116735599299999983e23) * pow(x2, 0.4e1) + (0.8154120192e16 * pow(x1, 0.4e1) + x1 * x1 * 0.163265453500000000e18 + 0.106966514200000016e18) * x2 * x2 + 0.4851666591e12 * x1 * x1 + 0.2595553280e11 * pow(x1, 0.4e1)) * x3 * x3 - 0.6535477512e28 * pow(x2, 0.10e2) + (-0.2262223094e24 * x1 * x1 + 0.5772786257e26) * pow(x2, 0.8e1) + (x1 * x1 * 0.200420021400000057e22 - 0.4584243200e18 * pow(x1, 0.4e1) + 0.389118664100000144e22) * pow(x2, 0.6e1) + (0.4077060096e16 * pow(x1, 0.4e1) + x1 * x1 * 0.81632726719999952e17 + 0.53483257099999976e17) * pow(x2, 0.4e1) + (0.4851666591e12 * x1 * x1 + 0.2595553280e11 * pow(x1, 0.4e1)) * x2 * x2) * pow(t, 0.6e1) + (-0.5900202000e23 * x1 * x1 + 0.1726427340e25 - 0.9899010000e27 * x2 * x2) * pow(x3, 0.8e1) + (-0.9899010000e27 * pow(x2, 0.8e1) + (-0.2360080800e24 * x1 * x1 + 0.6905709360e25) * pow(x2, 0.6e1) + (-0.1176120000e20 * pow(x1, 0.4e1) + x1 * x1 * 0.154352860200000029e22 + pow(x1, 0.3e1) * (-0.208920880595542258e5) + x1 * (-0.127400973168038654e5) + 0.945624821999999910e21) * pow(x2, 0.4e1) + (0.6837336000e17 * pow(x1, 0.4e1) + x1 * x1 * 0.187876614800000032e18 + 0.28480747599999996e17) * x2 * x2 + 0.4243913380e13 * x1 * x1 + 0.6241748000e13 * pow(x1, 0.4e1)) * x3 * x3 + (-0.5900202000e23 * x1 * x1 + 0.1726427340e25) * pow(x2, 0.8e1) + (x1 * x1 * 0.514509534000000074e21 - 0.3920400000e19 * pow(x1, 0.4e1) + 0.315208273999999992e21 + pow(x1, 0.3e1) * (-0.696402935318474192e4) + x1 * (-0.424669910560128847e4)) * pow(x2, 0.6e1) + (0.3418668000e17 * pow(x1, 0.4e1) + x1 * x1 * 0.93938307400000016e17 + 0.14240373799999998e17) * pow(x2, 0.4e1) + (0.4243913380e13 * x1 * x1 + 0.6241748000e13 * pow(x1, 0.4e1)) * x2 * x2 + (-0.1979802000e28 * pow(x2, 0.6e1) + (-0.3540121200e24 * x1 * x1 + 0.1035856404e26) * pow(x2, 0.4e1) + (-0.1176120000e20 * pow(x1, 0.4e1) + x1 * x1 * 0.154352860200000029e22 + pow(x1, 0.3e1) * (-0.208920880595542258e5) + x1 * (-0.127400973168038654e5) + 0.945624821999999910e21) * x2 * x2 + 0.3418668000e17 * pow(x1, 0.4e1) + x1 * x1 * 0.93938307400000016e17 + 0.14240373799999998e17) * pow(x3, 0.4e1) + (-0.3737387312e28 * pow(x3, 0.10e2) + (-0.7766018620e23 * x1 * x1 + 0.3306256935e26 - 0.1868693656e29 * x2 * x2) * pow(x3, 0.8e1) + (-0.3737387310e29 * pow(x2, 0.4e1) + (-0.3106407451e24 * x1 * x1 + 0.1322502777e27) * x2 * x2 + x1 * x1 * 0.689076817099999937e21 - 0.5242880000e17 * pow(x1, 0.4e1) + 0.177952274800000002e22) * pow(x3, 0.6e1) + (-0.3737387310e29 * pow(x2, 0.6e1) + (-0.4659611173e24 * x1 * x1 + 0.1983754165e27) * pow(x2, 0.4e1) + (-0.1572864000e18 * pow(x1, 0.4e1) + x1 * x1 * 0.206723045100000012e22 + pow(x1, 0.3e1) * (-0.223517417907714844e4) + x1 * (-0.577114476225855378e4) + 0.533856824699999853e22) * x2 * x2 + 0.4666163200e15 * pow(x1, 0.4e1) + x1 * x1 * 0.18673448520000004e17 + 0.18323256470000000e17) * pow(x3, 0.4e1) + (-0.1868693657e29 * pow(x2, 0.8e1) + (-0.3106407449e24 * x1 * x1 + 0.1322502777e27) * pow(x2, 0.6e1) + (-0.1572864000e18 * pow(x1, 0.4e1) + x1 * x1 * 0.206723045100000038e22 + pow(x1, 0.3e1) * (-0.223517417907714844e4) + x1 * (-0.577114476333936909e4) + 0.533856824799999925e22) * pow(x2, 0.4e1) + (0.9332326400e15 * pow(x1, 0.4e1) + x1 * x1 * 0.37346897020000000e17 + 0.36646512930000008e17) * x2 * x2 + 0.5544101804e11 * x1 * x1) * x3 * x3 - 0.3737387312e28 * pow(x2, 0.10e2) + (-0.7766018620e23 * x1 * x1 + 0.3306256936e26) * pow(x2, 0.8e1) + (x1 * x1 * 0.689076817099999937e21 - 0.5242880000e17 * pow(x1, 0.4e1) + 0.177952274800000002e22) * pow(x2, 0.6e1) + (x1 * x1 * 0.18673448510000000e17 + 0.4666163200e15 * pow(x1, 0.4e1) + 0.18323256460000004e17) * pow(x2, 0.4e1) + 0.5544101804e11 * x2 * x2 * x1 * x1) * pow(t, 0.7e1) + (-0.7121127936e28 * pow(x3, 0.10e2) + (-0.8371271936e24 * x1 * x1 + 0.6254872838e26 - 0.3560563967e29 * x2 * x2) * pow(x3, 0.8e1) + (-0.7121127935e29 * pow(x2, 0.4e1) + (-0.3348508774e25 * x1 * x1 + 0.2501949134e27) * x2 * x2 + x1 * x1 * 0.736866526099999818e22 - 0.1423370240e20 * pow(x1, 0.4e1) + 0.735843410499999891e22) * pow(x3, 0.6e1) + (-0.7121127935e29 * pow(x2, 0.6e1) + (-0.5022763162e25 * x1 * x1 + 0.3752923702e27) * pow(x2, 0.4e1) + (-0.4270110720e20 * pow(x1, 0.4e1) + x1 * x1 * 0.221059957600000032e23 + pow(x1, 0.3e1) * (-0.151704807649366558e6) + x1 * (-0.151203064742564718e6) + 0.220753023200000035e23) * x2 * x2 + 0.1257704141e18 * pow(x1, 0.4e1) + x1 * x1 * 0.726018835900000000e18 + 0.199602102699999968e18) * pow(x3, 0.4e1) + (-0.3560563967e29 * pow(x2, 0.8e1) + (-0.3348508774e25 * x1 * x1 + 0.2501949134e27) * pow(x2, 0.6e1) + (-0.4270110720e20 * pow(x1, 0.4e1) + x1 * x1 * 0.221059957600000032e23 + pow(x1, 0.3e1) * (-0.151704807649366558e6) + x1 * (-0.151203064742564718e6) + 0.220753023200000035e23) * pow(x2, 0.4e1) + (0.2515408283e18 * pow(x1, 0.4e1) + x1 * x1 * 0.1452037672000000000e19 + 0.399204205400000000e18) * x2 * x2 + 0.1517770058e14 * x1 * x1 + 0.8086770688e13 * pow(x1, 0.4e1)) * x3 * x3 - 0.7121127940e28 * pow(x2, 0.10e2) + (-0.8371271936e24 * x1 * x1 + 0.6254872838e26) * pow(x2, 0.8e1) + (x1 * x1 * 0.736866526099999818e22 - 0.1423370240e20 * pow(x1, 0.4e1) + 0.735843410499999891e22) * pow(x2, 0.6e1) + (x1 * x1 * 0.726018835800000000e18 + 0.1257704141e18 * pow(x1, 0.4e1) + 0.199602102800000000e18) * pow(x2, 0.4e1) + (0.1517770058e14 * x1 * x1 + 0.8086770688e13 * pow(x1, 0.4e1)) * x2 * x2 + 0.721888256e8 * pow(x1, 0.4e1)) * pow(t, 0.3e1) + (-0.6441914071e23 * pow(x3, 0.10e2) + (-0.3220957036e24 * x2 * x2 + 0.573234789400000004e21) * pow(x3, 0.8e1) + (-0.6441914072e24 * pow(x2, 0.4e1) + x2 * x2 * 0.229293915799999978e22 + 0.850510898900000e15) * pow(x3, 0.6e1) + (x2 * x2 * 0.255153269599999950e16 - 0.6441914072e24 * pow(x2, 0.6e1) + pow(x2, 0.4e1) * 0.343940873899999966e22) * pow(x3, 0.4e1) + (-0.3220957036e24 * pow(x2, 0.8e1) + pow(x2, 0.4e1) * 0.255153269600000050e16 + pow(x2, 0.6e1) * 0.229293915800000030e22) * x3 * x3 + 0.5732347895e21 * pow(x2, 0.8e1) - 0.6441914071e23 * pow(x2, 0.10e2) + pow(x2, 0.6e1) * 0.850510898900000125e15) * pow(t, 0.14e2) + (-0.8589934588e20 * pow(x3, 0.10e2) + (-0.4294967296e21 * x2 * x2 + 0.7645041787e18) * pow(x3, 0.8e1) + (0.3058016714e19 * x2 * x2 - 0.8589934592e21 * pow(x2, 0.4e1)) * pow(x3, 0.6e1) + (0.4587025074e19 * pow(x2, 0.4e1) - 0.8589934592e21 * pow(x2, 0.6e1)) * pow(x3, 0.4e1) + (-0.4294967296e21 * pow(x2, 0.8e1) + 0.3058016715e19 * pow(x2, 0.6e1)) * x3 * x3 + 0.7645041787e18 * pow(x2, 0.8e1) - 0.8589934588e20 * pow(x2, 0.10e2)) * pow(t, 0.16e2) + (-0.5987346947e27 * pow(x3, 0.10e2) + (-0.3457679360e22 * x1 * x1 + 0.5310055471e25 - 0.2993673473e28 * x2 * x2) * pow(x3, 0.8e1) + (-0.5987346944e28 * pow(x2, 0.4e1) + (-0.1383071744e23 * x1 * x1 + 0.2124022190e26) * x2 * x2 + 0.3074534932e20 * x1 * x1 + 0.1661900372e21) * pow(x3, 0.6e1) + (-0.5987346944e28 * pow(x2, 0.6e1) + (-0.2074607617e23 * x1 * x1 + 0.3186033286e26) * pow(x2, 0.4e1) + (0.9223604800e20 * x1 * x1 + 0.4985701116e21) * x2 * x2 + 0.2491731149e15 * x1 * x1 + 0.8135045384e15) * pow(x3, 0.4e1) + (-0.2993673473e28 * pow(x2, 0.8e1) + (-0.1383071744e23 * x1 * x1 + 0.2124022190e26) * pow(x2, 0.6e1) + (0.9223604800e20 * x1 * x1 + 0.4985701116e21) * pow(x2, 0.4e1) + (0.4983462297e15 * x1 * x1 + 0.1627009077e16) * x2 * x2) * x3 * x3 - 0.5987346946e27 * pow(x2, 0.10e2) + (-0.3457679360e22 * x1 * x1 + 0.5310055472e25) * pow(x2, 0.8e1) + (0.3074534932e20 * x1 * x1 + 0.1661900372e21) * pow(x2, 0.6e1) + (0.2491731149e15 * x1 * x1 + 0.8135045377e15) * pow(x2, 0.4e1)) * pow(t, 0.9e1) + (-0.1979802000e28 * pow(x2, 0.4e1) + (-0.2360080800e24 * x1 * x1 + 0.6905709360e25) * x2 * x2 + x1 * x1 * 0.514509534000000074e21 - 0.3920400000e19 * pow(x1, 0.4e1) + 0.315208273999999992e21 + pow(x1, 0.3e1) * (-0.696402935318474192e4) + x1 * (-0.424669910560128847e4)) * pow(x3, 0.6e1) + 0.281987600e9 * pow(x1, 0.4e1) - 0.1979802000e27 * pow(x2, 0.10e2) - 0.1979802000e27 * pow(x3, 0.10e2)) * pow(0.1e1 + 0.4e0 * t, -0.2e1) * pow(0.256e3 * pow(t, 0.4e1) * x2 * x2 + 0.256e3 * pow(t, 0.4e1) * x3 * x3 + 0.2560e4 * pow(t, 0.3e1) * x2 * x2 + 0.2560e4 * pow(t, 0.3e1) * x3 * x3 + 0.9600e4 * t * t * x2 * x2 + 0.9600e4 * t * t * x3 * x3 + 0.16000e5 * t * x2 * x2 + 0.16000e5 * t * x3 * x3 + x1 * x1 + 0.10000e5 * x2 * x2 + 0.10000e5 * x3 * x3, -0.2e1) * pow(0.160000e6 * t * t * x2 * x2 + 0.160000e6 * t * t * x3 * x3 + 0.800000e6 * t * x2 * x2 + 0.800000e6 * t * x3 * x3 + 0.990000e6 * x2 * x2 + 0.990000e6 * x3 * x3 + 0.89e2, -0.2e1) / (0.320000e6 * pow(t, 0.3e1) * x2 * x2 + 0.320000e6 * pow(t, 0.3e1) * x3 * x3 + 0.2400000e7 * t * t * x2 * x2 + 0.2400000e7 * t * t * x3 * x3 + 0.5980000e7 * t * x2 * x2 + 0.5980000e7 * t * x3 * x3 + 0.4950000e7 * x2 * x2 + 0.4950000e7 * x3 * x3 + 0.178e3 * t + 0.445e3);
}

static RegisterScalarFunction regsca_thin_ellipsoid_rhs( "thin_ellipsoidRhs", thin_ellipsoid_rhs);









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
double ellipsoid2_lset (const DROPS::Point3DCL& p, double)
{
 return p[0]*p[0]+p[1]*p[1]+p[2]*p[2]/2.25-1.0;
}
static RegisterScalarFunction regsca_ellipsoid2_lset( "Ellipsoid2_Lset", ellipsoid2_lset);
double Testu_scalar_fun (const DROPS::Point3DCL& p, double)
{
     if(p[0]==0 && p[1]==0 && p[2]==0){
         return p[0] * 0.0;
     }
     else {
    //return 3 * p[0] * p[0] - p[1] * p[2];  // Fall 1
     //return (0.3e1 * pow(p[0], 0.9e1) - pow(p[1], 0.6e1) * p[2] * p[2]) / (p[0] * p[0] + p[1] * p[1] + 0.4e1 / 0.9e1 * p[2] * p[2]);
         return p[0] * p[0] * p[1] - p[0] * p[1] * p[2] - 0.1e-3 * p[0];
    }
}
static RegisterScalarFunction regsca_philipskram_sol( "PhilipsKramSol", Testu_scalar_fun);
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

double toroidal_flow_sol2 (const Point3DCL& , double )
{
    return 1.;
}
static RegisterScalarFunction regsca_toroidal_flow_sol2( "ToroidalFlowSol2", toroidal_flow_sol2);

double toroidal_flow_rhs2 (const Point3DCL& , double )
{
    return 0.;
}
static RegisterScalarFunction regsca_toroidal_flow_rhs2( "ToroidalFlowRhs2", toroidal_flow_rhs2);
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
    const double bf4= (p/MakePoint3D( std::pow( a, 2.), 1., 1.)).norm_sq();
    const double bf6= (p/MakePoint3D( std::pow( a, 3.), 1., 1.)).norm_sq();

    const double mat_der=  0.25*std::cos( t)/a - 0.5;
    const double reaction= 0.25*std::cos( t)/a/bf4*( std::pow( p[1], 2.) + std::pow( p[2], 2.));
    const double diffusion= (-2./std::pow( a, 2.) - (1. + 1./std::pow( a, 2.))*( 2. + 1./std::pow( a, 2.) - bf6/bf4))/bf4;

//     const Point3DCL tt( p/MakePoint3D( std::pow( a, 2), 1., 1.));
//     const double l= tt.norm();
//     const Point3DCL n( tt/l);
//     SMatrixCL<3,3> dn( (eye<3,3>() - outer_product( n, n))/l );
//     dn(0,0)/= std::pow( a, 2); dn(1,0)/= std::pow( a, 2); dn(2,0)/= std::pow( a, 2);

//     const Point3DCL w( MakePoint3D( 0.25*p[0]*std::cos( t)/a, 0., 0.));
//     SMatrixCL<3,3> dw;
//     dw(0,0)= 0.25*std::cos( t)/a;

//     const Point3DCL grad_u( std::exp( -0.5*t)*MakePoint3D( p[1], p[0], 0.));
//     SMatrixCL<3,3> Hess_u;
//     Hess_u(0,1)= 1.;
//     Hess_u(1,0)= 1.;
//     Hess_u*= std::exp( -0.5*t);

//     const double err= div_gamma_wind( n, dw) - reaction,
//                errlb= laplace_beltrami_u( n, dn, grad_u, Hess_u) - diffusion*axis_scaling_sol( p, t);
//     if (std::fabs( err) > 1e-12 || std::fabs( errlb) > 1e-12) {
//         std::cerr << err   << " " << div_gamma_wind( n, dw) << " " << reaction << "\n"
//                   << errlb << " " << laplace_beltrami_u( n, dn, grad_u, Hess_u) << " " << diffusion*axis_scaling_sol( p, t) << "\n";
//         exit( 1);
//     }
//     return (mat_der + div_gamma_wind( n, dw))*axis_scaling_sol( p, t) - laplace_beltrami_u( n, dn, grad_u, Hess_u);
    return (mat_der + reaction - diffusion)*axis_scaling_sol( p, t);
//    double x1=p[0],x2=p[1],x3=p[2];
//    return (0.3125e-1 * x1 * x2 * (0.2e1 * pow(x2, 0.4e1) + 0.2e1 * pow(x3, 0.4e1) + 0.4e1 * x2 * x2 * x3 * x3) * pow(cos(t), 0.5e1) + (0.3125e-1 * x1 * x2 * (x2 * x2 + x3 * x3 - 0.5e0 * pow(x2, 0.4e1) - 0.5e0 * pow(x3, 0.4e1) - 0.1e1 * x2 * x2 * x3 * x3) * sin(t) + 0.3125e-1 * x1 * x2 * (0.20e2 * x2 * x2 + 0.20e2 * x3 * x3 - 0.10e2 * pow(x2, 0.4e1) - 0.10e2 * pow(x3, 0.4e1) - 0.20e2 * x2 * x2 * x3 * x3)) * pow(cos(t), 0.4e1) + (0.3125e-1 * x1 * x2 * (-0.32e2 * pow(x2, 0.4e1) - 0.32e2 * pow(x3, 0.4e1) - 0.64e2 * x2 * x2 * x3 * x3) * sin(t) + 0.3125e-1 * x1 * x2 * (-0.148e3 * pow(x2, 0.4e1) - 0.48e2 * x2 * x2 - 0.148e3 * pow(x3, 0.4e1) - 0.296e3 * x2 * x2 * x3 * x3 - 0.48e2 * x3 * x3)) * pow(cos(t), 0.3e1) + (0.3125e-1 * x1 * x2 * (0.65e2 * pow(x2, 0.4e1) + 0.65e2 * pow(x3, 0.4e1) + 0.130e3 * x2 * x2 * x3 * x3 - 0.32e2 - 0.178e3 * x2 * x2 - 0.178e3 * x3 * x3) * sin(t) + 0.3125e-1 * x1 * x2 * (0.148e3 * pow(x2, 0.4e1) + 0.148e3 * pow(x3, 0.4e1) + 0.296e3 * x2 * x2 * x3 * x3 - 0.384e3 - 0.872e3 * x2 * x2 - 0.872e3 * x3 * x3)) * pow(cos(t), 0.2e1) + (0.3125e-1 * x1 * x2 * (0.160e3 * pow(x2, 0.4e1) + 0.160e3 * pow(x3, 0.4e1) + 0.320e3 * x2 * x2 * x3 * x3 + 0.384e3 * x2 * x2 + 0.384e3 * x3 * x3) * sin(t) + 0.3125e-1 * x1 * x2 * (0.146e3 * pow(x2, 0.4e1) + 0.146e3 * pow(x3, 0.4e1) + 0.292e3 * x2 * x2 * x3 * x3 + 0.304e3 * x2 * x2 + 0.304e3 * x3 * x3 + 0.256e3)) * cos(t) + 0.3125e-1 * x1 * x2 * (0.2464e4 - 0.645e2 * pow(x2, 0.4e1) - 0.645e2 * pow(x3, 0.4e1) - 0.129e3 * x2 * x2 * x3 * x3 + 0.1713e4 * x2 * x2 + 0.1713e4 * x3 * x3) * sin(t) + 0.3125e-1 * x1 * x2 * (0.6016e4 - 0.138e3 * pow(x2, 0.4e1) - 0.138e3 * pow(x3, 0.4e1) - 0.276e3 * x2 * x2 * x3 * x3 + 0.852e3 * x2 * x2 + 0.852e3 * x3 * x3)) * pow(exp(t), -0.1e1 / 0.2e1) / (((0.625e-1 * x2 * x2 * x3 * x3 + 0.3125e-1 * pow(x2, 0.4e1) + 0.3125e-1 * pow(x3, 0.4e1)) * sin(t) + 0.625e0 * pow(x2, 0.4e1) + 0.625e0 * pow(x3, 0.4e1) + 0.125e1 * x2 * x2 * x3 * x3) * pow(cos(t), 0.4e1) + ((-0.40625e1 * pow(x2, 0.4e1) + (-0.1e1 - 0.8125e1 * x3 * x3) * x2 * x2 - 0.1e1 * x3 * x3 - 0.40625e1 * pow(x3, 0.4e1)) * sin(t) - 0.925e1 * pow(x2, 0.4e1) + (-0.12e2 - 0.185e2 * x3 * x3) * x2 * x2 - 0.925e1 * pow(x3, 0.4e1) - 0.12e2 * x3 * x3) * pow(cos(t), 0.2e1) + (0.403125e1 * pow(x2, 0.4e1) + (0.80625e1 * x3 * x3 + 0.33e2) * x2 * x2 + 0.403125e1 * pow(x3, 0.4e1) + 0.33e2 * x3 * x3 + 0.8e1) * sin(t) + 0.8625e1 * pow(x2, 0.4e1) + (0.1725e2 * x3 * x3 + 0.12e2) * x2 * x2 + 0.12e2 * x3 * x3 + 0.32e2 + 0.8625e1 * pow(x3, 0.4e1));
}
static RegisterScalarFunction regsca_axis_scaling_rhs( "AxisScalingRhs", axis_scaling_rhs);

double axis_scaling_sol2 (const Point3DCL& , double )
{
    return 1.;
}
static RegisterScalarFunction regsca_axis_scaling_sol2( "AxisScalingSol2", axis_scaling_sol2);

double axis_scaling_rhs2 (const Point3DCL& p, double t)
{
  //double x1=p[0];
  double x2=p[1];
  double x3=p[2];
  return -cos(t) * ((x2 * x2 + x3 * x3) * sin(t) + 0.4e1 * x2 * x2 + 0.4e1 * x3 * x3) / (-0.8e1 * sin(t) * x2 * x2 - 0.8e1 * sin(t) * x3 * x3 + pow(cos(t), 0.2e1) * x2 * x2 + pow(cos(t), 0.2e1) * x3 * x3 - x2 * x2 - x3 * x3 - 0.16e2);
}
static RegisterScalarFunction regsca_axis_scaling_rhs2( "AxisScalingRhs2", axis_scaling_rhs2);


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

AutoDiff::ADFwdCL<double, 4> lset_AD(const DROPS::Point3DCL& x, double t)
{
    AutoDiff::ADFwdCL<double, 4> p[4];
    for (int i= 0; i < 3; ++i)
            p[i].seed (x[i], i);
    p[3].seed(t,3);
    AutoDiff::ADFwdCL<double, 4> n1= pow((pow((p[0].value()-1.5*p[3].value()+1.5),2)+pow(p[1].value(),2)+pow(p[2].value(),2)),0.5)
                                  < 1e-3 ? 1e-3 : pow((pow((p[0]-1.5*p[3]+1.5),2)+pow(p[1],2)+pow(p[2],2)),1.5);
    AutoDiff::ADFwdCL<double, 4> n2= pow((pow((p[0].value()+1.5*p[3].value()-1.5),2)+pow(p[1].value(),2)+pow(p[2].value(),2)),0.5)
                                  < 1e-3 ? 1e-3 : pow((pow((p[0]+1.5*p[3]-1.5),2)+pow(p[1],2)+pow(p[2],2)),1.5);
    return 1-(1/n1)-(1/n2);
}
AutoDiff::ADFwdCL<double, 4> norm_phi_sq_AD(const DROPS::Point3DCL& x, double t)
{
    AutoDiff::ADFwdCL<double, 4> p[4];
    for (int i= 0; i < 3; ++i)
            p[i].seed (x[i], i);
    p[3].seed(t,3);
	return pow(0.3e1 / 0.2e1 * pow(pow(p[0] - 0.3e1 / 0.2e1 * p[3] + 0.3e1 / 0.2e1, 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (0.2e1 * p[0] - 0.3e1 * p[3] + 0.3e1) + 0.3e1 / 0.2e1 * pow(pow(p[0] + 0.3e1 / 0.2e1 * p[3] - 0.3e1 / 0.2e1, 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (0.2e1 * p[0] + 0.3e1 * p[3] - 0.3e1), 0.2e1) + pow(0.3e1 * pow(pow(p[0] - 0.3e1 / 0.2e1 * p[3] + 0.3e1 / 0.2e1, 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[1] + 0.3e1 * pow(pow(p[0] + 0.3e1 / 0.2e1 * p[3] - 0.3e1 / 0.2e1, 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[1], 0.2e1) + pow(0.3e1 * pow(pow(p[0] - 0.3e1 / 0.2e1 * p[3] + 0.3e1 / 0.2e1, 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[2] + 0.3e1 * pow(pow(p[0] + 0.3e1 / 0.2e1 * p[3] - 0.3e1 / 0.2e1, 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[2], 0.2e1);
}
AutoDiff::ADFwdCL<double, 4> wind_AD1(DROPS::Point3DCL x, double t)
{
    AutoDiff::ADFwdCL<double, 4> p[4];
    for (int i= 0; i < 3; ++i)
            p[i].seed (x[i], i);
    p[3].seed(t,3);
	return norm_phi_sq_AD(x,t).value() < 1e-6 ? 0 : (1/norm_phi_sq_AD(x,t))*(-(0.3e1 / 0.2e1 * pow(pow(p[0] - 0.3e1 / 0.2e1 * p[3] + 0.3e1 / 0.2e1, 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (-0.3e1 * p[0] + 0.9e1 / 0.2e1 * p[3] - 0.9e1 / 0.2e1) + 0.3e1 / 0.2e1 * pow(pow(p[0] + 0.3e1 / 0.2e1 * p[3] - 0.3e1 / 0.2e1, 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (0.3e1 * p[0] + 0.9e1 / 0.2e1 * p[3] - 0.9e1 / 0.2e1)) * (0.3e1 / 0.2e1 * pow(pow(p[0] - 0.3e1 / 0.2e1 * p[3] + 0.3e1 / 0.2e1, 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (0.2e1 * p[0] - 0.3e1 * p[3] + 0.3e1) + 0.3e1 / 0.2e1 * pow(pow(p[0] + 0.3e1 / 0.2e1 * p[3] - 0.3e1 / 0.2e1, 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (0.2e1 * p[0] + 0.3e1 * p[3] - 0.3e1)));
}
AutoDiff::ADFwdCL<double, 4> wind_AD2(DROPS::Point3DCL x, double t)
{
    AutoDiff::ADFwdCL<double, 4> p[4];
    for (int i= 0; i < 3; ++i)
            p[i].seed (x[i], i);
    p[3].seed(t,3);
return norm_phi_sq_AD(x,t).value() < 1e-6 ? 0 : (1/norm_phi_sq_AD(x,t))*(-(0.3e1 / 0.2e1 * pow(pow(p[0] - 0.3e1 / 0.2e1 * p[3] + 0.3e1 / 0.2e1, 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (-0.3e1 * p[0] + 0.9e1 / 0.2e1 * p[3] - 0.9e1 / 0.2e1) + 0.3e1 / 0.2e1 * pow(pow(p[0] + 0.3e1 / 0.2e1 * p[3] - 0.3e1 / 0.2e1, 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (0.3e1 * p[0] + 0.9e1 / 0.2e1 * p[3] - 0.9e1 / 0.2e1)) * (0.3e1 * pow(pow(p[0] - 0.3e1 / 0.2e1 * p[3] + 0.3e1 / 0.2e1, 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[1] + 0.3e1 * pow(pow(p[0] + 0.3e1 / 0.2e1 * p[3] - 0.3e1 / 0.2e1, 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[1]));
}
AutoDiff::ADFwdCL<double, 4> wind_AD3(DROPS::Point3DCL x, double t)
{
    AutoDiff::ADFwdCL<double, 4> p[4];
    for (int i= 0; i < 3; ++i)
            p[i].seed (x[i], i);
    p[3].seed(t,3);
return norm_phi_sq_AD(x,t).value() < 1e-6 ? 0 : (1/norm_phi_sq_AD(x,t))*(-(0.3e1 / 0.2e1 * pow(pow(p[0] - 0.3e1 / 0.2e1 * p[3] + 0.3e1 / 0.2e1, 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (-0.3e1 * p[0] + 0.9e1 / 0.2e1 * p[3] - 0.9e1 / 0.2e1) + 0.3e1 / 0.2e1 * pow(pow(p[0] + 0.3e1 / 0.2e1 * p[3] - 0.3e1 / 0.2e1, 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (0.3e1 * p[0] + 0.9e1 / 0.2e1 * p[3] - 0.9e1 / 0.2e1)) * (0.3e1 * pow(pow(p[0] - 0.3e1 / 0.2e1 * p[3] + 0.3e1 / 0.2e1, 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[2] + 0.3e1 * pow(pow(p[0] + 0.3e1 / 0.2e1 * p[3] - 0.3e1 / 0.2e1, 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[2]));
}
double tang_div(AutoDiff::ADFwdCL<double, 4> w1, AutoDiff::ADFwdCL<double, 4> w2,
                          AutoDiff::ADFwdCL<double, 4> w3, AutoDiff::ADFwdCL<double, 4> phi)
{
    DROPS::Point3DCL grad_phi;
    grad_phi[0]=phi.derivative()[0];
    grad_phi[1]=phi.derivative()[1];
    grad_phi[2]=phi.derivative()[2];
    DROPS::Point3DCL n=grad_phi.norm()<1e-3? DROPS::Point3DCL() :   (1/grad_phi.norm())*grad_phi;
    SMatrixCL<3,3> dw;
    dw(0,0)=w1.derivative()[0];
    dw(0,1)=w1.derivative()[1];
    dw(0,2)=w1.derivative()[2];
    dw(1,0)=w2.derivative()[0];
    dw(1,1)=w2.derivative()[1];
    dw(1,2)=w2.derivative()[2];
    dw(2,0)=w3.derivative()[0];
    dw(2,1)=w3.derivative()[1];
    dw(2,2)=w3.derivative()[2];
  return trace(dw)-inner_prod(n, dw*n);
}

DROPS::Point3DCL canonical_wind(const DROPS::Point3DCL& x, double t)
{
    AutoDiff::ADFwdCL<double,4> phi=lset_AD(x,t);
    double phi_t=phi.derivative()[3];
    DROPS::Point3DCL grad_phi;
    grad_phi[0]=phi.derivative()[0];
    grad_phi[1]=phi.derivative()[1];
    grad_phi[2]=phi.derivative()[2];
    return (-phi_t/grad_phi.norm_sq())*grad_phi;
}

double collision_sol (const Point3DCL& , double )
{
    return 1;
}
static RegisterScalarFunction regsca_collision_sol( "CollisionSol", collision_sol);
double collision_sol2 (const Point3DCL& , double )
{
     return  0;
}
static RegisterScalarFunction regsca_collision_sol2( "CollisionSol2", collision_sol2);

double collision_rhs (const Point3DCL& x, double t)
{
    double rr=tang_div(wind_AD1(x,t),wind_AD2(x,t),wind_AD3(x,t),lset_AD(x,t));
    return rr;
}
static RegisterScalarFunction regsca_collision_rhs( "CollisionRhs", collision_rhs);

double collision_rhs2 (const Point3DCL& , double )
{
    return 0;
}
static RegisterScalarFunction regsca_collision_rhs2( "CollisionRhs2", collision_rhs2);

//======================================================================================
//Non-stationay test case Split

const double split_p= 3.0;


Point3DCL split_center_1 (double t)
{
    return -WindVelocity*t;
}

Point3DCL split_center_2 (double t)
{
    return -split_center_1( t);
}
double split_Dt_lset (const DROPS::Point3DCL& x, double t)
{
    Point3DCL x1= x - split_center_1( t),
              x2= x - split_center_2( t);
    const double n1= x1.norm() < 1e-3 ? 1e-3 : std::pow( x1.norm(), split_p + 2.),
                 n2= x2.norm() < 1e-3 ? 1e-3 : std::pow( x2.norm(), split_p + 2.);
    return split_p*(DROPS::inner_prod(WindVelocity, x1)/n1
                      - DROPS::inner_prod( WindVelocity, x2)/n2);
}
DROPS::Point3DCL split_Dx_lset (const DROPS::Point3DCL& x, double t)
{
    Point3DCL x1= x - split_center_1( t),
              x2= x - split_center_2( t);
    const double n1= x1.norm() < 1e-3 ? 1e-3 : std::pow( x1.norm(), split_p + 2.),
                 n2= x2.norm() < 1e-3 ? 1e-3 : std::pow( x2.norm(), split_p + 2.);
    return split_p*(x1/n1 + x2/n2);
}
DROPS::Point3DCL split_wind (const DROPS::Point3DCL& x, double t)
{
    const Point3DCL Dphi= split_Dx_lset( x, t);
    const double Dtphi=   split_Dt_lset( x, t);
    const double n= Dphi.norm_sq();
    return n < 1e-6 ? DROPS::Point3DCL() : (-Dtphi/n)*Dphi;
}
static RegisterVectorFunction regvec_split_wind( "SplitWind", split_wind);

double split_lset (const Point3DCL& x, double t)
{
    const Point3DCL c1= split_center_1( t);
    const Point3DCL x1= x - c1;
    const Point3DCL c2= split_center_2( t);
    const Point3DCL x2= x - c2;

    const double n1= x1.norm() < 1e-3 ? 1e-3 : std::pow( x1.norm(), collision_p),
                 n2= x2.norm() < 1e-3 ? 1e-3 : std::pow( x2.norm(), collision_p);
    return RadDrop[0] - 1./n1  - 1./n2;
}
static RegisterScalarFunction regsca_split_lset( "SplitLset", split_lset);

AutoDiff::ADFwdCL<double, 4> lset_AD_split(const DROPS::Point3DCL& x, double t)
{
    AutoDiff::ADFwdCL<double, 4> p[4];
    for (int i= 0; i < 3; ++i)
            p[i].seed (x[i], i);
    p[3].seed(t,3);
    return 1-(1/pow((pow((p[0]-1.5*p[3]),2)+pow(p[1],2)+pow(p[2],2)),1.5))-(1/pow((pow((p[0]+1.5*p[3]),2)+pow(p[1],2)+pow(p[2],2)),1.5));
}
AutoDiff::ADFwdCL<double, 4> wind_AD1_split(DROPS::Point3DCL x, double t)
{
    AutoDiff::ADFwdCL<double, 4> p[4];
    for (int i= 0; i < 3; ++i)
            p[i].seed (x[i], i);
    p[3].seed(t,3);
    return -(0.3e1 / 0.2e1 * pow(pow(p[0] + 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (0.30e1 * p[0] + 0.450e1 * p[3]) + 0.3e1 / 0.2e1 * pow(pow(p[0] - 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (-0.30e1 * p[0] + 0.450e1 * p[3])) / (pow(0.3e1 / 0.2e1 * pow(pow(p[0] + 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (0.2e1 * p[0] + 0.30e1 * p[3]) + 0.3e1 / 0.2e1 * pow(pow(p[0] - 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (0.2e1 * p[0] - 0.30e1 * p[3]), 0.2e1) + pow(0.3e1 * pow(pow(p[0] + 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[1] + 0.3e1 * pow(pow(p[0] - 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[1], 0.2e1) + pow(0.3e1 * pow(pow(p[0] + 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[2] + 0.3e1 * pow(pow(p[0] - 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[2], 0.2e1)) * (0.3e1 / 0.2e1 * pow(pow(p[0] + 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (0.2e1 * p[0] + 0.30e1 * p[3]) + 0.3e1 / 0.2e1 * pow(pow(p[0] - 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (0.2e1 * p[0] - 0.30e1 * p[3]));
}
AutoDiff::ADFwdCL<double, 4> wind_AD2_split(DROPS::Point3DCL x, double t)
{
    AutoDiff::ADFwdCL<double, 4> p[4];
    for (int i= 0; i < 3; ++i)
            p[i].seed (x[i], i);
    p[3].seed(t,3);
    return -(0.3e1 / 0.2e1 * pow(pow(p[0] + 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (0.30e1 * p[0] + 0.450e1 * p[3]) + 0.3e1 / 0.2e1 * pow(pow(p[0] - 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (-0.30e1 * p[0] + 0.450e1 * p[3])) / (pow(0.3e1 / 0.2e1 * pow(pow(p[0] + 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (0.2e1 * p[0] + 0.30e1 * p[3]) + 0.3e1 / 0.2e1 * pow(pow(p[0] - 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (0.2e1 * p[0] - 0.30e1 * p[3]), 0.2e1) + pow(0.3e1 * pow(pow(p[0] + 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[1] + 0.3e1 * pow(pow(p[0] - 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[1], 0.2e1) + pow(0.3e1 * pow(pow(p[0] + 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[2] + 0.3e1 * pow(pow(p[0] - 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[2], 0.2e1)) * (0.3e1 * pow(pow(p[0] + 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[1] + 0.3e1 * pow(pow(p[0] - 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[1]);
}
AutoDiff::ADFwdCL<double, 4> wind_AD3_split(DROPS::Point3DCL x, double t)
{
    AutoDiff::ADFwdCL<double, 4> p[4];
    for (int i= 0; i < 3; ++i)
            p[i].seed (x[i], i);
    p[3].seed(t,3);
    return -(0.3e1 / 0.2e1 * pow(pow(p[0] + 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (0.30e1 * p[0] + 0.450e1 * p[3]) + 0.3e1 / 0.2e1 * pow(pow(p[0] - 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (-0.30e1 * p[0] + 0.450e1 * p[3])) / (pow(0.3e1 / 0.2e1 * pow(pow(p[0] + 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (0.2e1 * p[0] + 0.30e1 * p[3]) + 0.3e1 / 0.2e1 * pow(pow(p[0] - 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * (0.2e1 * p[0] - 0.30e1 * p[3]), 0.2e1) + pow(0.3e1 * pow(pow(p[0] + 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[1] + 0.3e1 * pow(pow(p[0] - 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[1], 0.2e1) + pow(0.3e1 * pow(pow(p[0] + 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[2] + 0.3e1 * pow(pow(p[0] - 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[2], 0.2e1)) * (0.3e1 * pow(pow(p[0] + 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[2] + 0.3e1 * pow(pow(p[0] - 0.15e1 * p[3], 0.2e1) + p[1] * p[1] + p[2] * p[2], -0.5e1 / 0.2e1) * p[2]);
}
double split_sol (const Point3DCL&, double)
{
    return 1;
}
static RegisterScalarFunction regsca_split_sol( "SplitSol", split_sol);

double split_rhs (const Point3DCL& x, double t)
{
    double rr=tang_div(wind_AD1_split(x,t),wind_AD2_split(x,t),wind_AD3_split(x,t),lset_AD_split(x,t));
    return rr;
}
static RegisterScalarFunction regsca_split_rhs( "SplitRhs", split_rhs);



// ==Some spheres==

double sphere_dist (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL x( p - PosDrop);
    return x.norm() - RadDrop[0];
}
static RegisterScalarFunction regsca_sphere_dist_lset( "SphereDist", sphere_dist);

DROPS::Point3DCL d_sphere_dist (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL x( p - PosDrop);
    return x/x.norm();
}


typedef double (*dist_funT) (const DROPS::Point3DCL&, double);

double sphere_2move (const DROPS::Point3DCL& p, double t)
{
    DROPS::Point3DCL x( p - (PosDrop + t*constant_wind(p, t)));
    return x.norm() - RadDrop[0];
}

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


SurfactantP1BaseCL* make_surfactant_timedisc( MultiGridCL& mg, LevelsetP2CL& lset,
                                              VecDescCL& v, const BndDataCL<Point3DCL>& Bnd_v,
                                              const ParamCL& P)
{
    SurfactantP1BaseCL* ret= 0;
    const std::string method= P.get<std::string>( "SurfTransp.Method");

    if (method == std::string( "cGcG"))
        ret= new SurfactantcGP1CL( mg,
            P.get<double>("SurfTransp.Theta"), P.get<double>("SurfTransp.Visc"),
            &v, Bnd_v, lset.Phi, lset.GetBndData(),
            P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"),
            P.get<double>("SurfTransp.XFEMReduced"));
    else if (method == std::string( "spacetime-cGdG"))
        ret= new SurfactantSTP1CL( mg,
            P.get<double>("SurfTransp.Theta"), P.get<double>("SurfTransp.Visc"),
            &v, Bnd_v, lset.Phi, lset.GetBndData(),
            /* cG_in_t_ */ false, /* use_mass_div */ P.get<bool>( "SurfTransp.UseMassDiv"),
            P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"),
            P.get<double>("SurfTransp.XFEMReduced"));
    else if (method == std::string( "spacetime-cGcG"))
        ret= new SurfactantSTP1CL( mg,
            P.get<double>("SurfTransp.Theta"), P.get<double>("SurfTransp.Visc"),
            &v, Bnd_v, lset.Phi, lset.GetBndData(),
            /* cG_in_t_ */ true, /* use_mass_div */ P.get<bool>( "SurfTransp.UseMassDiv"),
            P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"),
            P.get<double>("SurfTransp.XFEMReduced"));
    else if (method == std::string( "characteristic-transport"))
        ret= new SurfactantCharTransportP1CL ( mg,
            P.get<double>("SurfTransp.Theta"), P.get<double>("SurfTransp.Visc"),
            &v, Bnd_v, lset.Phi, lset.GetBndData(),
            P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"),
            P.get<double>("SurfTransp.XFEMReduced"));
    else if (method == std::string( "spacetime-cGdG-Mixed"))
        ret= new SurfactantSTP1CL( mg,
            P.get<double>("SurfTransp.Theta"), P.get<double>("SurfTransp.Visc"),
            &v, Bnd_v, lset.Phi, lset.GetBndData(),
            /* cG_in_t_ */ false, /* use_mass_div */ P.get<bool>( "SurfTransp.UseMassDiv"),
            P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"),
            P.get<double>("SurfTransp.XFEMReduced"),/*use MixedFormulation*/ true);
    else
        throw DROPSErrCL( std::string( "make_surfactant_timedisc: Unknown method '") + method + std::string( "'.\n"));

    return ret;
}


void Strategy (DROPS::MultiGridCL& mg, DROPS::AdapTriangCL& adap, DROPS::LevelsetP2CL& lset)
{
    using namespace DROPS;
    adap.MakeInitialTriang();
    if (P.get<std::string>("SurfTransp.Exp.Levelset") == std::string( "AxisScalingLset"))
        dynamic_cast<DistMarkingStrategyCL*>( adap.get_marking_strategy())->SetDistFct( axis_scaling_lset_ini);
   // if (P.get<std::string>("SurfTransp.Exp.Levelset") == std::string( "ConstVelLset"))
   //     dynamic_cast<DistMarkingStrategyCL*>( adap.get_marking_strategy())->SetDistFct( const_vel_lset_ini);
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

    std::unique_ptr<SurfactantP1BaseCL> timediscp( make_surfactant_timedisc( mg, lset, v, Bnd_v, P));
    SurfactantP1BaseCL& timedisc= *timediscp;
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

    BndDataCL<> nobnd( 0);
    VecDescCL the_sol_vd( &lset.idx);
    LSInit( mg, the_sol_vd, the_sol_fun, /*t*/ 0.);
    if (vtkwriter.get() != 0) {
        vtkwriter->Register( make_VTKScalar(      lset.GetSolution(),              "Levelset") );
        vtkwriter->Register( make_VTKIfaceScalar( mg, timedisc.ic,                 "InterfaceSol"));
        vtkwriter->Register( make_VTKVector(      make_P2Eval( mg, Bnd_v, v),      "Velocity"));
        //vtkwriter->Register( make_VTKScalar(      lset2.GetSolution(),             "Levelset2"));
        vtkwriter->Register( make_VTKScalar(      make_P2Eval( mg, nobnd, the_sol_vd),  "TrueSol"));
        vtkwriter->Write( 0.);
    }
    //if (P.get<int>( "SurfTransp.SolutionOutput.Freq") > 0)
    //    DROPS::WriteFEToFile( timedisc.ic, mg, P.get<std::string>( "SurfTransp.SolutionOutput.Path"), P.get<bool>( "SolutionOutput.Binary"));

    const double dt= P.get<double>("Time.FinalTime")/P.get<double>("Time.NumSteps");
    double L_2x_err= L2_error( lset.Phi, lset.GetBndData(), timedisc.GetSolution(), the_sol_fun);
    std::cout << "L_2x-error: " << L_2x_err
              << "\nnorm of true solution: " << L2_norm( mg, lset.Phi, lset.GetBndData(), the_sol_fun)
              << std::endl;
    double L_inftL_2x_err= L_2x_err;
    std::cout << "L_inftL_2x-error: " <<  L_inftL_2x_err << std::endl;
    double H_1x_err= H1_error( lset.Phi, lset.GetBndData(), timedisc.GetSolution(), the_sol_fun);
    std::cout << "H_1x-error: " << H_1x_err << std::endl;
    double L_2tH_1x_err_sq= 0.5*dt*std::pow( H_1x_err, 2);
    double L_2tL_2x_err_sq= 0.5*dt*std::pow( L_2x_err, 2);
    BndDataCL<> ifbnd( 0);
    std::cout << "initial surfactant on \\Gamma: " << Integral_Gamma( mg, lset.Phi, lset.GetBndData(), make_P1Eval(  mg, ifbnd, timedisc.ic)) << '\n';

    dynamic_cast<DistMarkingStrategyCL*>( adap.get_marking_strategy())->SetDistFct( lset);
    for (int step= 1; step <= P.get<int>("Time.NumSteps"); ++step) {
        std::cout << "======================================================== step " << step << ":\n";
        ScopeTimerCL timer( "Strategy: Time-loop");
        const double cur_time= step*dt;
        // Assumes (as the rest of Drops), that the current triangulation is acceptable to perform the time-step.
        // If dt is large and AdapRef.Width is small, this may not be true.
        // Watch for large differences in numbers of old and new dof.
        timedisc.InitTimeStep();
        LSInit( mg, lset.Phi, the_lset_fun, cur_time);
        InitVel( mg, &v, Bnd_v, the_wind_fun, cur_time);
        timedisc.DoStep( cur_time);
//        P1EvalCL discsol;
//        discsol=make_P1Eval(  mg, ifbnd, timedisc.ic);
        std::cout << "surfactant on \\Gamma: " << Integral_Gamma( mg, lset.Phi, lset.GetBndData(), make_P1Eval(  mg, ifbnd, timedisc.ic)) << '\n';
        L_2tL_2x_err_sq+= (step > 1 ? 0.5 : 0.)*dt*std::pow( L_2x_err, 2);
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



        L_2tL_2x_err_sq+= 0.5*dt*std::pow( L_2x_err, 2);
        std::cout << "L_2tL_2x-error: " << std::sqrt( L_2tL_2x_err_sq) << std::endl;


        if (vtkwriter.get() != 0 && step % P.get<int>( "VTK.Freq") == 0) {
            LSInit( mg, the_sol_vd, the_sol_fun, /*t*/ cur_time);
            vtkwriter->Write( cur_time);
        }
        if (P.get<int>( "SurfTransp.SolutionOutput.Freq") > 0 && step % P.get<int>( "SurfTransp.SolutionOutput.Freq") == 0) {
            std::ostringstream os1,
                               os2;
            os1 << P.get<int>( "Time.NumSteps");
            os2 << P.get<std::string>( "SurfTransp.SolutionOutput.Path") << std::setw( os1.str().size()) << step;
            DROPS::WriteFEToFile( timedisc.ic, mg, os2.str(), P.get<bool>( "SurfTransp.SolutionOutput.Binary"));
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
        const bool doGridMod= P.get<int>("Mesh.AdaptRef.Freq") && step%P.get<int>("Mesh.AdaptRef.Freq") == 0;
        const bool gridChanged= doGridMod ? adap.UpdateTriang() : false;
         std::cout << "dogridmod " << doGridMod << std::endl;
         std::cout << "gridchanged " << gridChanged << std::endl;
        if (gridChanged) {
            std::cout << "Triangulation changed.\n";
            vidx.DeleteNumbering( mg);
            vidx.CreateNumbering( mg.GetLastLevel(), mg, Bnd_v);
            v.SetIdx( &vidx);
            InitVel( mg, &v, Bnd_v, the_wind_fun, cur_time);
            LSInit( mg, lset.Phi, the_lset_fun, cur_time);
            the_sol_vd.SetIdx( &lset.idx);
            LSInit( mg, the_sol_vd, the_sol_fun, /*t*/ cur_time);
            // timedisc.Update(); // Called unconditionally in DoStep.

            //lset2.SetupSystem( make_P2Eval( mg, Bnd_v, v),

            std::cout << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            lset.AdjustVolume();
            lset.GetVolumeAdjuster()->DebugOutput( std::cout);
        }
    }
    std::cout << std::endl;
    //delete &lset2;
}

void StationaryStrategyP1 (DROPS::MultiGridCL& mg, DROPS::AdapTriangCL& adap, DROPS::LevelsetP2CL& lset)
{
    adap.MakeInitialTriang();

    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
    // LinearLSInit( mg, lset.Phi, the_lset_fun);
    LSInit( mg, lset.Phi, the_lset_fun, 0.);

    DROPS::IdxDescCL ifaceidx( P1IF_FE);
    ifaceidx.GetXidx().SetBound( P.get<double>("SurfTransp.XFEMReduced"));
    ifaceidx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
    std::cout << "NumUnknowns: " << ifaceidx.NumUnknowns() << std::endl;

    DROPS::MatDescCL M( &ifaceidx, &ifaceidx);
    DROPS::SetupInterfaceMassP1( mg, &M, lset.Phi, lset.GetBndData());
    std::cout << "M is set up.\n";
    DROPS::MatDescCL A( &ifaceidx, &ifaceidx);
    DROPS::SetupLBP1( mg, &A, lset.Phi, lset.GetBndData(), P.get<double>("SurfTransp.Visc"));
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
    SurfSolverT surfsolver( surfpc, P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"), true);

    DROPS::VecDescCL x( &ifaceidx);
    surfsolver.Solve( L, x.Data, b.Data, x.RowIdx->GetEx());
    std::cout << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';

    if (P.get<int>( "SurfTransp.SolutionOutput.Freq") > 0)
        DROPS::WriteFEToFile( x, mg, P.get<std::string>( "SurfTransp.SolutionOutput.Path"), P.get<bool>( "SurfTransp.SolutionOutput.Binary"));

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

    double L2_err( L2_error( lset.Phi, lset.GetBndData(), make_P1Eval( mg, nobnd, xext), the_sol_fun));
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

    instat_scalar_fun_ptr f;
    double f_time;

    LocalLaplaceBeltramiP2CL loc_lb;

    instat_vector_fun_ptr f_grad;
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
    void set_function (const instat_scalar_fun_ptr farg, double f_time_arg= 0.) {
        f= farg;
        f_time= f_time_arg;
    }
    void set_grad_function (const instat_vector_fun_ptr farg, double f_time_arg= 0.) {
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

    instat_scalar_fun_ptr d_; // exact distance
    instat_vector_fun_ptr Dd_; // gradient of exact distance
    GridFunctionCL<> qd;
    GridFunctionCL<Point3DCL> qDderr;

    OpenMPVar_MinInit_Max_CL<double> max_h,
                                     max_d,
                                     max_Dderr;

  public:
    InterfaceApproxErrorDeformAccuCL (InterfaceCommonDataDeformP2CL* cdataarg, VecDescCL* yg, VecDescCL* ydist)
        : yG_( yg), ydist_( ydist), cdata_ (cdataarg), d_ (0), Dd_ (0) {}
    virtual ~InterfaceApproxErrorDeformAccuCL () {}

    InterfaceApproxErrorDeformAccuCL& set_d  (instat_scalar_fun_ptr darg)  {  d_= darg;  return *this; }
    InterfaceApproxErrorDeformAccuCL& set_Dd (instat_vector_fun_ptr Ddarg) { Dd_= Ddarg; return *this; }

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

    instat_scalar_fun_ptr f;
    double f_time;

    LocalLaplaceBeltramiDeformP2CL loc_lb;
    LocalNormalLaplaceDeformP2CL loc_ngrad;

    instat_vector_fun_ptr f_grad;
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
    void set_function (const instat_scalar_fun_ptr farg, double f_time_arg= 0.) {
        f= farg;
        f_time= f_time_arg;
    }
    void set_grad_function (const instat_vector_fun_ptr farg, double f_time_arg= 0.) {
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
    ifacep2idx.GetXidx().SetBound( P.get<double>("SurfTransp.XFEMReduced"));
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
    InterfaceMatrixAccuCL<LocalLaplaceBeltramiP2CL, InterfaceCommonDataP2CL> accuAp2( &Ap2, LocalLaplaceBeltramiP2CL( P.get<double>("SurfTransp.Visc")), cdatap2, "Ap2");
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
    SurfSolverT surfsolver( surfpc, P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"), true);

    DROPS::VecDescCL xp2( &ifacep2idx);
    surfsolver.Solve( Lp2, xp2.Data, bp2.Data, xp2.RowIdx->GetEx());
    std::cout << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';

    TetraAccumulatorTupleCL err_accus;
    err_accus.push_back( &cdatap2);
    InterfaceL2AccuP2CL L2_accu( cdatap2, mg, "P2-solution");
    L2_accu.set_grid_function( xp2);
    L2_accu.set_function( the_sol_fun, 0.);
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

    if (P.get<int>( "SurfTransp.SolutionOutput.Freq") > 0)
        DROPS::WriteFEToFile( xp2, mg, P.get<std::string>( "SurfTransp.SolutionOutput.Path") + "_p2", P.get<bool>( "SurfTransp.SolutionOutput.Binary"));

    DROPS::NoBndDataCL<> nobnd;
    DROPS::NoBndDataCL<Point3DCL> nobnd_vec;
    VecDescCL the_sol_vd( &lset.idx);
    LSInit( mg, the_sol_vd, the_sol_fun, /*t*/ 0.);
    if (vtkwriter.get() != 0) {
        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "Levelset") );
        vtkwriter->Register( make_VTKIfaceScalar( mg, xp2, "InterfaceSolP2"));
        vtkwriter->Register( make_VTKScalar(      make_P2Eval( mg, nobnd, the_sol_vd),  "TrueSol"));
        vtkwriter->Register( make_VTKVector( make_P2Eval( mg, nobnd_vec, lsgradrec), "LSGradRec") );
        vtkwriter->Register( make_VTKVector( make_P2Eval( mg, nobnd_vec, to_iface), "to_iface") );
        vtkwriter->Write( 0.);
    }
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
    ifacep2idx.GetXidx().SetBound( P.get<double>("SurfTransp.XFEMReduced"));
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
    InterfaceMatrixAccuCL<LocalLaplaceBeltramiDeformP2CL, InterfaceCommonDataDeformP2CL> accuAp2( &Ap2, LocalLaplaceBeltramiDeformP2CL( P.get<double>("SurfTransp.Visc")), cdatap2, "Ap2");
    accus.push_back( &accuAp2);
    DROPS::MatDescCL Anp2( &ifacep2idx, &ifacep2idx);
    InterfaceMatrixAccuCL<LocalNormalLaplaceDeformP2CL, InterfaceCommonDataDeformP2CL> accuAnp2( &Anp2, LocalNormalLaplaceDeformP2CL (P.get<double>("SurfTransp.NormalLaplaceCoefficent")), cdatap2, "Anp2");
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
    SurfSolverT surfsolver( surfpc, P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"), true);

    DROPS::VecDescCL xp2( &ifacep2idx);
    surfsolver.Solve( Lp2, xp2.Data, bp2.Data, xp2.RowIdx->GetEx());
    std::cout << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';

    TetraAccumulatorTupleCL err_accus;
    err_accus.push_back( &cdatap2);
    InterfaceL2AccuDeformP2CL L2_accu( cdatap2, mg, "deformed P2-solution");
    L2_accu.set_grid_function( xp2);
    L2_accu.set_function( the_sol_fun, 0.);
    L2_accu.set_grad_function( the_sol_grad_fun, 0.);
    err_accus.push_back( &L2_accu);
    accumulate( err_accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetBndInfo());

    if (P.get<int>( "SurfTransp.SolutionOutput.Freq") > 0)
        DROPS::WriteFEToFile( xp2, mg, P.get<std::string>( "SurfTransp.SolutionOutput.Path") + "_p2", P.get<bool>( "SurfTransp.SolutionOutput.Binary"));

    DROPS::NoBndDataCL<> nobnd;
    DROPS::NoBndDataCL<Point3DCL> nobnd_vec;
    VecDescCL the_sol_vd( &lset.idx);
    LSInit( mg, the_sol_vd, the_sol_fun, /*t*/ 0.);
    if (vtkwriter.get() != 0) {
        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "Levelset") );
        vtkwriter->Register( make_VTKIfaceScalar( mg, xp2, "InterfaceSolP2"));
        vtkwriter->Register( make_VTKScalar(      make_P2Eval( mg, nobnd, the_sol_vd),  "TrueSol"));
        vtkwriter->Register( make_VTKVector( make_P2Eval( mg, nobnd_vec, deformation), "deformation") );
//         vtkwriter->Register( make_VTKVector( make_P2Eval( mg, nobnd_vec, to_iface), "to_iface") );
        vtkwriter->Register( make_VTKScalar( make_P1Eval( mg, nobnd, d_iface_vd), "d_iface") );
        vtkwriter->Register( make_VTKScalar( make_P2Eval( mg, nobnd, locdist_vd), "locdist") );
        vtkwriter->Write( 0.);
    }
}

int main (int argc, char* argv[])
{
  try {
    ScopeTimerCL timer( "main");

    DROPS::read_parameter_file_from_cmdline( P, argc, argv, "../../param/surfactant/surfactant/surfactant.json");
    std::cout << P << std::endl;

    DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"), P.get<std::vector<std::string> >("General.DynamicLibs") );
    std::cout << "Setting up interface-PDE.\n";
    WindVelocity= P.get<DROPS::Point3DCL>("SurfTransp.Exp.Velocity");
    RadDrop=      P.get<DROPS::Point3DCL>("SurfTransp.Exp.RadDrop");
    PosDrop=      P.get<DROPS::Point3DCL>("SurfTransp.Exp.PosDrop");
    RadTorus=     P.get<DROPS::Point2DCL>("SurfTransp.Exp.RadTorus");
    the_wind_fun= invecmap[P.get<std::string>("SurfTransp.Exp.Wind")];
    the_lset_fun= inscamap[P.get<std::string>("SurfTransp.Exp.Levelset")];
    the_rhs_fun=  inscamap[P.get<std::string>("SurfTransp.Exp.Rhs")];
    the_sol_fun=  inscamap[P.get<std::string>("SurfTransp.Exp.Solution")];
    if (P.get<std::string>("SurfTransp.Exp.Solution") == "LaplaceBeltrami0Sol")
        the_sol_grad_fun=  &laplace_beltrami_0_sol_grad;
    for (Uint i= 0; i < 6; ++i)
        bf_wind[i]= the_wind_fun;

    std::cout << "Setting up domain:\n";
    std::unique_ptr<MGBuilderCL> builder( make_MGBuilder( P));
    DROPS::MultiGridCL mg( *builder);
    typedef DistMarkingStrategyCL MarkerT;
    MarkerT marker( the_lset_fun, P.get<double>( "Mesh.AdaptRef.Width"),
                    P.get<int>( "Mesh.AdaptRef.CoarsestLevel"), P.get<int>( "Mesh.AdaptRef.FinestLevel"));

    DROPS::AdapTriangCL adap( mg, &marker);

    // DROPS::LevelsetP2CL lset( mg, lsbnd, sf);
    DROPS::LevelsetP2CL& lset( *LevelsetP2CL::Create( mg, lsbnd, sf, P.get_child("Levelset")) );

    if (P.get<int>("VTK.Freq",0))
        vtkwriter= std::unique_ptr<VTKOutCL>( new VTKOutCL(
            adap.GetMG(),
            "DROPS data",
            P.get<int>("Time.NumSteps")/P.get<int>("VTK.Freq") + 1,
            P.get<std::string>("VTK.VTKDir"),
            P.get<std::string>("VTK.VTKName"),
            P.get<std::string>("VTK.TimeFileName"),
            P.get<int>("VTK.Binary"), 
            P.get<bool>("VTK.UseOnlyP1"),
            false, /* <- P2DG */
            -1,    /* <- level */
            P.get<bool>("VTK.ReUseTimeFile")));
    if (P.get<bool>( "SurfTransp.Exp.StationaryPDE")) {
        if (P.get<std::string>("LevelsetMapper.DeformationMethod") != "")
            StationaryStrategyDeformationP2( mg, adap, lset);
        else {
            if (P.get<int>( "SurfTransp.FEDegree") == 1)
                StationaryStrategyP1( mg, adap, lset);
            else
                StationaryStrategyP2( mg, adap, lset);
        }
    }
    else
        Strategy( mg, adap, lset);

    delete &lset;
    rusage usage;
    getrusage( RUSAGE_SELF, &usage);
    std::cout << "ru_maxrss: " << usage.ru_maxrss << " kB.\n";
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
