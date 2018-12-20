/// \file surfacestokes_funcs.h
/// \brief Functions for surfacestokes.cpp
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

#ifndef DROPS_SURFACESTOKES_FUNCS_H
#define DROPS_SURFACESTOKES_FUNCS_H

#include "misc/params.h"

using namespace DROPS;

DROPS::ParamCL P;

namespace ParameterNS {
	double h=0.1;
	double nu;
	double sigma;
	double eps;
}

DROPS::Point3DCL u_func (const DROPS::Point3DCL&, double)
{
    return P.get<DROPS::Point3DCL>("Exp.Velocity");
}

typedef DROPS::Point3DCL (*bnd_val_fun) (const DROPS::Point3DCL&, double);

DROPS::BndCondT bc[6]= {
    DROPS::DirBC, DROPS::DirBC,
    DROPS::DirBC, DROPS::DirBC,
    DROPS::DirBC, DROPS::DirBC
};

bnd_val_fun bf[6]= {
    &u_func, &u_func, &u_func, &u_func, &u_func, &u_func
};

///////////////////////// Levelset-functions //////////////////////////////////////

double sphere_2 (const DROPS::Point3DCL& p, double)
{
//    DROPS::Point3DCL x( p - P.get<DROPS::Point3DCL>("Exp.PosDrop"));

//    return x.norm() - P.get<DROPS::Point3DCL>("Exp.RadDrop")[0];
//    return p[0]*p[0]+p[1]*p[1]+p[2]*p[2]-1;
    return std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]) - 1;
}

double xy_plane (const DROPS::Point3DCL& p, double)
{
    const double PI = 3.141592653589793;
    return p[2] - PI/30.;
}

double levelset1 (const DROPS::Point3DCL& p, double)
{
    return 3*p[1] + p[0] + p[0]*p[0] + p[1]*p[1] + p[2]*p[2] - 1 + p[2]*p[2] + ((std::sqrt((p[0]-0.25)*(p[0]-0.25)+(p[1]-0.25)*(p[1]-0.25)))*(std::sqrt((p[0]-0.25)*(p[0]-0.25)+(p[1]-0.25)*(p[1]-0.25))) - 0.75)*((std::sqrt((p[0]-0.25)*(p[0]-0.25)+(p[1]-0.25)*(p[1]-0.25)))*(std::sqrt((p[0]-0.25)*(p[0]-0.25)+(p[1]-0.25)*(p[1]-0.25))) - 0.75) + p[2]*p[2] + ((std::sqrt((p[0])*(p[0])+(p[1])*(p[1])))*(std::sqrt((p[0]+0.25)*(p[0]+0.25)+(p[1]+0.25)*(p[1]+0.25))) - 0.75)*((std::sqrt((p[0]+0.25)*(p[0]+0.25)+(p[1]+0.25)*(p[1]+0.25)))*(std::sqrt((p[0]+0.25)*(p[0]+0.25)+(p[1]+0.25)*(p[1]+0.25))) - 0.5);
}

double tilted_plane (const DROPS::Point3DCL& p, double)
{
    double h = 2.0;
//    double InitialDivision = P.get<double>("Domain.N1");
//    double FinestLevel = P.get<double>("AdaptRef.FinestLevel");
//    double dx = P.get<DROPS::Point3DCL>("Domain.E1")[0];
//    if (P.get<int>("AdaptRef.FinestLevel") == 0) {
//        h = dx*(1/(InitialDivision));
//    } else {
//        h = dx*(1/InitialDivision)*(pow(2.0, -FinestLevel));
//    }
    /*if(-h*0.5*p[0]+2*p[1]>0 && -h*0.5*p[0]+2*p[1]<2.*h){
        return (-h*0.5*p[0]+2*p[1]-h);
    }else{
        return -(-h*0.5*p[0]+2*p[1]+h);
    }*/
   // return p[1]-0.05*p[0] - 0.1;

     return -h*0.5*p[0]+2*p[1]-4./3;
}

double tilted_plane_xy (const DROPS::Point3DCL& p, double)
{
    return -p[2]+4*p[0]-13./3;
}

typedef double (*dist_funT) (const DROPS::Point3DCL&, double);

double cube_madeof_edges (const DROPS::Point3DCL& p, double)
{
//    DROPS::Point3DCL x( p - P.get<DROPS::Point3DCL>("Exp.PosDrop"));

//    return x.norm() - P.get<DROPS::Point3DCL>("Exp.RadDrop")[0];
    return pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - 0.4e1, 0.2e1) + pow(pow(p[1], 0.2e1) - 0.1e1, 0.2e1) + pow(pow(p[1], 0.2e1) + pow(p[2], 0.2e1) - 0.4e1, 0.2e1) + pow(pow(p[0], 0.2e1) - 0.1e1, 0.2e1) + pow(pow(p[0], 0.2e1) + pow(p[2], 0.2e1) - 0.4e1, 0.2e1) + pow(pow(p[2], 0.2e1) - 0.1e1, 0.2e1) - 0.13e2;

;
}

namespace Torus {
    double R=1;
    double r=0.2;
       double k=0.0;
       double f=0.0;
//    double k=0.3;
//    double f=4.0;
}

double torus (const DROPS::Point3DCL& p, double)
{
    //return pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - 0.4e1, 0.2e1) + pow(pow(p[1], 0.2e1) - 0.1e1, 0.2e1) + pow(pow(p[1], 0.2e1) + pow(p[2], 0.2e1) - 0.4e1, 0.2e1) + pow(pow(p[0], 0.2e1) - 0.1e1, 0.2e1) + pow(pow(p[0], 0.2e1) + pow(p[2], 0.2e1) - 0.4e1, 0.2e1) + pow(pow(p[2], 0.2e1) - 0.1e1, 0.2e1) - 0.13e2;

   double R=Torus::R;
//    double r=0.4;
    //return pow(p[2],2.0) + pow(std::sqrt(pow(p[0],2.0) + pow(p[1], 2.0)) - R, 2.0) - r*r;
    return pow( pow(p[0],2.0) + pow(p[1], 2.0) + pow(p[2],2.0) + Torus::R*Torus::R - Torus::r*Torus::r, 2.0) - 4*Torus::R*Torus::R*(pow(p[0],2.0) + pow(p[1], 2.0));
    //return (p[0]*p[0]+p[1]*p[1]+p[2]*p[2]-1)*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]-0.33);
}

double flower_shape(double x)
{
// double   f=4.;
//double    k=0.;
    return(1+Torus::k*cos(Torus::f*(x+3.141582/Torus::f)));
}


double torus_flower (const DROPS::Point3DCL& p, double)
{
//      Torus::R=1;
//    Torus::r=0.2;
//    Torus::k=0.3;
//    Torus::f=4.0;
    //return pow( pow(p[0],2.0) + pow(p[1], 2.0) + pow(p[2],2.0) + R*R - r*r, 2.0) - 4*R*R*(pow(p[0],2.0) + pow(p[1], 2.0));


    double theta=atan(p[2]/(std::sqrt(pow(p[0],2.0) + pow(p[1], 2.0)) - Torus::R));
    return pow(p[2],2.0) + pow(std::sqrt(pow(p[0],2.0) + pow(p[1], 2.0)) - Torus::R, 2.0) - Torus::r*Torus::r*pow(flower_shape(theta),2.0);
}

/////////////////////////////////////////////////////////////////////////////////////
////////////////////////// General functions for test cases /////////////////////////////////

double ConstantOneScalarFun (const DROPS::Point3DCL&, double)
{
    return 1;
}

double ZeroScalarFun (const DROPS::Point3DCL&, double)
{
    return 0;
}

DROPS::Point3DCL ZeroVectorFun (const DROPS::Point3DCL&, double)
{
    DROPS::Point3DCL v(0,0,0);
    return v;
}

DROPS::Point3DCL ConstantOneVectorFun (const DROPS::Point3DCL&, double)
{
    DROPS::Point3DCL v(1,1,1);
    return v;
}

DROPS::Point3DCL Constant_e1_VectorFun (const DROPS::Point3DCL&, double)
{
    DROPS::Point3DCL v(1,0,0);
    return v;
}

DROPS::Point3DCL Constant_e2_VectorFun (const DROPS::Point3DCL&, double)
{
    DROPS::Point3DCL v(0,1,0);
    return v;
}

DROPS::Point3DCL Constant_e3_VectorFun (const DROPS::Point3DCL&, double)
{
    DROPS::Point3DCL v(0,0,1);
    return v;
}



DROPS::Point3DCL ProjectedConstantOneVectorFun (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL e(1,1,1);
    double norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    DROPS::Point3DCL n(p[0]/norm, p[1]/norm, p[2]/norm);
    return e - inner_prod(n,e)*n;
}

DROPS::Point3DCL IdentityVectorFun (const DROPS::Point3DCL& p, double)
{
    return p;
}

DROPS::Point3DCL NormalsVectorFun (const DROPS::Point3DCL& p, double)
{
    double norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    //DROPS::Point3DCL n(p[0]/norm, p[1]/norm, p[2]/norm);
    return p/norm;
}

DROPS::Point3DCL MotionVectorFun (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v((pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)), -p[0] * p[1] * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)), -p[0] * p[2] * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)));
    //DROPS::Point3DCL v(0, 0, p[2]);
    //double norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    //double norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    //DROPS::Point3DCL n(p[0]/norm, p[1]/norm, p[2]/norm);
    //return v - inner_prod(n,v)*n;
    return v;
}

DROPS::Point3DCL TrigVectorFun (const DROPS::Point3DCL& p, double)
{
    double rad = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    const double PI = 3.141592653589793;
    double theta = acos(p[2]/rad);
    double phi;
    if(p[0] == 0){
        if(p[1] > 0) phi = PI/2;
        if(p[1] < 0) phi = 3*PI/2;
        if(p[1] == 0) phi = 0;
    }
    if(p[0] > 0) {
        if(p[1] > 0) phi = atan(p[1]/p[0]);
        if(p[1] < 0) phi = -atan(std::abs(p[1])/p[0]);
        if(p[1] == 0) phi = 0;
    }
    if(p[0] < 0) {
        if(p[1] > 0) phi = atan(std::abs(p[0])/p[1]) + PI/2;
        if(p[1] < 0) phi = atan(std::abs(p[1])/std::abs(p[0])) + PI;
        if(p[1] == 0) phi = PI;
    }
    if(std::abs(rad*sin(phi)*sin(theta) - p[1]) > 1e-14) { std::cout << "Y IST NICHT KORREKT: " << std::abs(rad*sin(phi)*sin(theta) - p[1]) << std::endl;

    std::cout << p[0] << "   " << p[1] << "   " << p[2] << "   " << phi << "   " << theta << "   "  << -rad*sin(phi)*sin(theta) << "   " << rad*cos(phi)*sin(theta) << std::endl;
    }
    DROPS::Point3DCL v(-rad*sin(phi)*sin(theta)*cos(3*phi+theta)*cos(phi), rad*cos(phi)*sin(theta), 0);
//    DROPS::Point3DCL v(-p[2]*sin(3*phi +theta)*cos(phi)*cos(phi) - cos(phi +3*theta)*sin(3*phi)*p[1], cos(phi +3*theta)*sin(3*phi)*p[0] - sin(3*phi+theta)*cos(phi)*sin(phi)*p[2], sin(3*phi + theta)*cos(phi)*cos(theta));
//    DROPS::Point3DCL v(-p[2] - p[1], p[0] - p[2], 0);
    //DROPS::Point3DCL l(-p[1], p[0], 0);
    //double norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    //double norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    //DROPS::Point3DCL n(p[0]/norm, p[1]/norm, p[2]/norm);
    //return v - inner_prod(n,v)*n;
    //return l;
    return v;
}

DROPS::Point3DCL TrigVectorFun2 (const DROPS::Point3DCL& p, double)
{
    double rad = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    const double PI = M_PI;
    double theta = acos(p[2]/rad);
    double phi;
    if(p[0] == 0){
        if(p[1] > 0) phi = PI/2;
        if(p[1] < 0) phi = 3*PI/2;
        if(p[1] == 0) phi = 0;
    }
    if(p[0] > 0) {
        if(p[1] > 0) phi = atan(p[1]/p[0]);
        if(p[1] < 0) phi = -atan(std::abs(p[1])/p[0]);
        if(p[1] == 0) phi = 0;
    }
    if(p[0] < 0) {
        if(p[1] > 0) phi = atan(std::abs(p[0])/p[1]) + PI/2;
        if(p[1] < 0) phi = atan(std::abs(p[1])/std::abs(p[0])) + PI;
        if(p[1] == 0) phi = PI;
    }
    if(std::abs(rad*cos(phi)*sin(theta) - p[0]) > 1e-11) { std::cout << "X IST NICHT KORREKT: " << std::abs(rad*cos(phi)*sin(theta) - p[0]) << std::endl; }
    if(std::abs(rad*sin(phi)*sin(theta) - p[1]) > 1e-11) { std::cout << "Y IST NICHT KORREKT: " << std::abs(rad*sin(phi)*sin(theta) - p[1]) << std::endl; }
    if(std::abs(rad*cos(theta) - p[2]) > 1e-11) { std::cout << "Z IST NICHT KORREKT: " << std::abs(rad*cos(theta) - p[2]) << std::endl; }

    DROPS::Point3DCL v(0, 0, 0);
    return v;
}

DROPS::Point3DCL MotionVectorFun2 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(-(pow(p[1], 0.2e1) * pow(p[2], 0.2e1) + pow(p[2], 0.4e1) + pow(p[0], 0.2e1) * p[2] + p[0] * pow(p[1], 0.2e1)) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1))
            , p[1] * (p[0] * pow(p[2], 0.2e1) + pow(p[0], 0.2e1) - p[0] * p[2] + pow(p[2], 0.2e1)) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1))
, (p[0] * pow(p[2], 0.3e1) + pow(p[0], 0.3e1) + p[0] * pow(p[1], 0.2e1) - pow(p[1], 0.2e1) * p[2]) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)));
//    DROPS::Point3DCL v(p[0], p[2], 0);
    //double norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
//    double norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
//    DROPS::Point3DCL n(p[0]/norm, p[1]/norm, p[2]/norm);
//    return v - inner_prod(n,v)*n;
    return v;
}

DROPS::Point3DCL MotionVectorFun4 (const DROPS::Point3DCL& p, double)
{
   // const double PI = 3.141592653589793;
   // DROPS::Point3DCL v (pow(p[2], 0.2e1), pow(p[0], 0.2e1) - p[0] * p[1], p[2]*p[2]);
    DROPS::Point3DCL v (0, 0, p[0]*p[0]);
    return v;
}

DROPS::Point3DCL MotionVectorFun3 (const DROPS::Point3DCL& p, double)
{
    const double PI = 3.141592653589793;
    DROPS::Point3DCL v ((4./5.)*sin(PI*p[0])*sin(PI*p[2]), (2./5.)*sin(PI*p[0])*sin(PI*p[2]), sin(PI*p[0])*sin(PI*p[2]));
    return v;
}

DROPS::Point3DCL MotionVector_xy_plane1 (const DROPS::Point3DCL& p, double)
{
    const double PI = 3.141592653589793;
    DROPS::Point3DCL v (sin(PI*p[0])*sin(PI*p[1]), sin(PI*p[0])*sin(PI*p[1]), 0);
    return v;
}

DROPS::Point3DCL RhsVector_xy_plane1 (const DROPS::Point3DCL& p, double)
{
    double theta = 1;
    DROPS::Point3DCL v (-cos(0.3141592654e1 * p[0]) * cos(0.3141592654e1 * p[1]) * 0.3141592654e1 * 0.3141592654e1 * theta + 0.3e1 * 0.3141592654e1 * 0.3141592654e1 * sin(0.3141592654e1 * p[0]) * sin(0.3141592654e1 * p[1]) * theta + sin(0.3141592654e1 * p[0]) * sin(0.3141592654e1 * p[1])
    , -cos(0.3141592654e1 * p[0]) * cos(0.3141592654e1 * p[1]) * 0.3141592654e1 * 0.3141592654e1 * theta + 0.3e1 * 0.3141592654e1 * 0.3141592654e1 * sin(0.3141592654e1 * p[0]) * sin(0.3141592654e1 * p[1]) * theta + sin(0.3141592654e1 * p[0]) * sin(0.3141592654e1 * p[1])
    , 0);
    return v;
}

double PressureFun1 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    return (p[0]-2)*(p[0]-2)*(p[0]+2)*(p[0]+2);
}

double PressureFun2 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    return p[1];
}

DROPS::Point3DCL RhsVectorFun (const DROPS::Point3DCL& p, double)
{
    double theta = 1;
    DROPS::Point3DCL v((-pow(p[0], 0.2e1) * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - pow(p[0], 0.2e1) * pow(p[2], 0.4e1) - pow(p[1], 0.4e1) * pow(p[2], 0.2e1) - 0.2e1 * pow(p[1], 0.2e1) * pow(p[2], 0.4e1) - pow(p[2], 0.6e1) + 0.2e1 * theta * pow(p[0], 0.2e1) * pow(p[1], 0.2e1) + 0.14e2 * theta * pow(p[0], 0.2e1) * pow(p[2], 0.2e1) + 0.2e1 * theta * pow(p[1], 0.4e1) - 0.6e1 * theta * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - 0.8e1 * theta * pow(p[2], 0.4e1) - pow(p[0], 0.4e1) * p[2] - pow(p[0], 0.3e1) * pow(p[1], 0.2e1) - pow(p[0], 0.2e1) * pow(p[1], 0.2e1) * p[2] - pow(p[0], 0.2e1) * pow(p[2], 0.3e1) - p[0] * pow(p[1], 0.4e1) - p[0] * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - 0.5e1 * theta * pow(p[0], 0.2e1) * p[2] - 0.10e2 * theta * p[0] * pow(p[1], 0.2e1) + 0.5e1 * theta * pow(p[1], 0.2e1) * p[2] + 0.5e1 * theta * pow(p[2], 0.3e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1)
, -p[1] * (-pow(p[0], 0.3e1) * pow(p[2], 0.2e1) - p[0] * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - p[0] * pow(p[2], 0.4e1) + 0.2e1 * theta * pow(p[0], 0.3e1) + 0.2e1 * theta * p[0] * pow(p[1], 0.2e1) - 0.20e2 * theta * p[0] * pow(p[2], 0.2e1) - pow(p[0], 0.4e1) + pow(p[0], 0.3e1) * p[2] - pow(p[0], 0.2e1) * pow(p[1], 0.2e1) - 0.2e1 * pow(p[0], 0.2e1) * pow(p[2], 0.2e1) + p[0] * pow(p[1], 0.2e1) * p[2] + p[0] * pow(p[2], 0.3e1) - pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - pow(p[2], 0.4e1) - 0.10e2 * theta * pow(p[0], 0.2e1) + 0.10e2 * theta * p[0] * p[2] - 0.10e2 * theta * pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1)
, -(-pow(p[0], 0.3e1) * pow(p[2], 0.3e1) - p[0] * pow(p[1], 0.2e1) * pow(p[2], 0.3e1) - p[0] * pow(p[2], 0.5e1) + 0.14e2 * theta * pow(p[0], 0.3e1) * p[2] + 0.14e2 * theta * p[0] * pow(p[1], 0.2e1) * p[2] - 0.8e1 * theta * p[0] * pow(p[2], 0.3e1) - pow(p[0], 0.5e1) - 0.2e1 * pow(p[0], 0.3e1) * pow(p[1], 0.2e1) - pow(p[0], 0.3e1) * pow(p[2], 0.2e1) + pow(p[0], 0.2e1) * pow(p[1], 0.2e1) * p[2] - p[0] * pow(p[1], 0.4e1) - p[0] * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) + pow(p[1], 0.4e1) * p[2] + pow(p[1], 0.2e1) * pow(p[2], 0.3e1) - 0.5e1 * theta * pow(p[0], 0.3e1) - 0.5e1 * theta * p[0] * pow(p[1], 0.2e1) + 0.5e1 * theta * p[0] * pow(p[2], 0.2e1) + 0.10e2 * theta * pow(p[1], 0.2e1) * p[2]) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1)
);
    //DROPS::Point3DCL v(0, 0, p[2]);
    //double norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    //double norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    //DROPS::Point3DCL n(p[0]/norm, p[1]/norm, p[2]/norm);
    //return v - inner_prod(n,v)*n;
    return v;
}

double RhsFun (const DROPS::Point3DCL& p, double)
{
    return p[0]*p[1];
}

DROPS::Point3DCL RhsVectorFun2 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(0.5e1 * (pow(p[0], 0.2e1) * p[1] - pow(p[0], 0.2e1) * p[2] - pow(p[1], 0.3e1) + pow(p[1], 0.2e1) * p[2] - p[1] * pow(p[2], 0.2e1) + pow(p[2], 0.3e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1)
            , -0.5e1 * p[0] * (pow(p[0], 0.2e1) - pow(p[1], 0.2e1) + 0.2e1 * p[1] * p[2] + pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1)
            , 0.5e1 * p[0] * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + 0.2e1 * p[1] * p[2] - pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1));
    //DROPS::Point3DCL v(0, 0, p[2]);
    //double norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    //double norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    //DROPS::Point3DCL n(p[0]/norm, p[1]/norm, p[2]/norm);
    //return v - inner_prod(n,v)*n;
    return v;
}

/////////////////////////////////////////////////////////////////////////////
//////////////////////// Functions for TestRhs //////////////////////////////

DROPS::Point3DCL TestRhsVectorFun1 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(-p[1], p[0], 0);
    return v;
}

DROPS::Point3DCL TestRhsVectorFun2 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(-p[2], p[0], 0);
    double norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    DROPS::Point3DCL n(p[0]/norm, p[1]/norm, p[2]/norm);
    return v - inner_prod(n,v)*n;
}

DROPS::Point3DCL TestRhsVectorFun3 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(p[0]*p[2], -p[1], p[2]);
    return v;
}

DROPS::Point3DCL TestRhsVectorFun4 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(0, p[1]*p[2], p[0]);
    double norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    DROPS::Point3DCL n(p[0]/norm, p[1]/norm, p[2]/norm);
    return v - inner_prod(n,v)*n;
}

DROPS::Point3DCL TestRhsVectorFun_v1 (const DROPS::Point3DCL&, double)
{
    DROPS::Point3DCL v(1, 0, 0);
    return v;
}

DROPS::Point3DCL TestRhsVectorFun_v2 (const DROPS::Point3DCL&, double)
{
    DROPS::Point3DCL v(0, 1, 0);
    return v;
}

DROPS::Point3DCL TestRhsVectorFun_v3 (const DROPS::Point3DCL&, double)
{
    DROPS::Point3DCL v(0, 0, 1);
    return v;
}

DROPS::Point3DCL TestRhsVectorFun_v4 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(p[0], p[2], 0);
    return v;
}

DROPS::Point3DCL TestRhsVectorFun_v5 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(0, p[0]*p[1], p[0]);
    return v;
}

DROPS::Point3DCL TestRhsVectorFun_v6 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(p[1], p[2]*p[2], p[0]);
    return v;
}

/////////////////////////////////////////////////////////////////////////////
///////////////////// Functions to test matrices ////////////////////////////

DROPS::Point3DCL TestL_stabVectorFun_v1 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(0, 0, p[2]);
    return v;
}

DROPS::Point3DCL TestL_stabVectorFun_v2 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(0, 0, p[0]*p[2]);
    return v;
}

DROPS::Point3DCL TestL_stabVectorFun_v3 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(0, 0, p[1]*p[2]);
    return v;
}

DROPS::Point3DCL TestL_stabVectorFun_v4 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(0, 0, p[2]*p[2]);
    return v;
}

DROPS::Point3DCL TestA_VectorFun1 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(0.2e1 * pow(p[0], 0.2e1) + 0.3e1 * pow(p[1], 0.2e1) + p[0] + 0.2e1 * p[1]
    , 0.13e2 * pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + 0.7e1 * p[0] + 0.5e1 * p[1]
    , 0);
    return v;
}

DROPS::Point3DCL TestA_VectorFun_v1 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(0.8e1 * pow(p[0], 0.2e1) + 0.9e1 * pow(p[1], 0.2e1) + 0.3e1 * p[0] + 0.4e1 * p[1]
    , 0.5e1 * pow(p[0], 0.2e1) + 0.3e1 * pow(p[1], 0.2e1) + 0.4e1 * p[0] + 0.2e1 * p[1]
    , 0);
    return v;
}

double TestB_PressureFun1 (const DROPS::Point3DCL& p, double)
{
    return p[0];
}

double TestB_PressureFun2 (const DROPS::Point3DCL& p, double)
{
    return p[1]*p[2];
}

double TestB_PressureFun3 (const DROPS::Point3DCL& p, double)
{
    return p[0]*p[0];
}

double TestB_PressureFun4 (const DROPS::Point3DCL& p, double)
{
    return p[0]*(p[2]-p[1]);
}

double TestL_LagrangeFun1 (const DROPS::Point3DCL& p, double)
{
    return p[0];
}

double TestL_LagrangeFun2 (const DROPS::Point3DCL& p, double)
{
    return p[1]*p[2];
}

double TestL_LagrangeFun3 (const DROPS::Point3DCL& p, double)
{
    return p[0]*p[0];
}

double TestL_LagrangeFun4 (const DROPS::Point3DCL& p, double)
{
    return p[0]*(p[2]-p[1]);
}

double TestL_LagrangeFun21 (const DROPS::Point3DCL& p, double)
{
    return p[0];
}

double TestL_LagrangeFun22 (const DROPS::Point3DCL& p, double)
{
    return p[1];
}

double TestL_LagrangeFun23 (const DROPS::Point3DCL& p, double)
{
    return p[2];
}

double TestL_LagrangeFun24 (const DROPS::Point3DCL&, double)
{
    return 1;
}

double TestL_stab_LagrangeFun1 (const DROPS::Point3DCL& p, double)
{
    return p[2];
}

double TestL_stab_LagrangeFun2 (const DROPS::Point3DCL& p, double)
{
    return p[0]*p[2];
}

double TestL_stab_LagrangeFun3 (const DROPS::Point3DCL& p, double)
{
    return p[1]*p[2];
}

double TestL_stab_LagrangeFun4 (const DROPS::Point3DCL& p, double)
{
    return p[2]*p[2];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////// Functions for test cases with rhs computed in Maple //////////////////////////////////////

DROPS::Point3DCL Test_A_plus_M_RhsVectorFun1 (const DROPS::Point3DCL& p, double)
{
    double theta = 1;
    DROPS::Point3DCL v((-pow(p[0], 0.2e1) * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - pow(p[0], 0.2e1) * pow(p[2], 0.4e1) - pow(p[1], 0.4e1) * pow(p[2], 0.2e1) - 0.2e1 * pow(p[1], 0.2e1) * pow(p[2], 0.4e1) - pow(p[2], 0.6e1) + 0.2e1 * theta * pow(p[0], 0.2e1) * pow(p[1], 0.2e1) + 0.14e2 * theta * pow(p[0], 0.2e1) * pow(p[2], 0.2e1) + 0.2e1 * theta * pow(p[1], 0.4e1) - 0.6e1 * theta * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - 0.8e1 * theta * pow(p[2], 0.4e1) - pow(p[0], 0.4e1) * p[2] - pow(p[0], 0.3e1) * pow(p[1], 0.2e1) - pow(p[0], 0.2e1) * pow(p[1], 0.2e1) * p[2] - pow(p[0], 0.2e1) * pow(p[2], 0.3e1) - p[0] * pow(p[1], 0.4e1) - p[0] * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - 0.5e1 * theta * pow(p[0], 0.2e1) * p[2] - 0.10e2 * theta * p[0] * pow(p[1], 0.2e1) + 0.5e1 * theta * pow(p[1], 0.2e1) * p[2] + 0.5e1 * theta * pow(p[2], 0.3e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1)
    , -p[1] * (-pow(p[0], 0.3e1) * pow(p[2], 0.2e1) - p[0] * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - p[0] * pow(p[2], 0.4e1) + 0.2e1 * theta * pow(p[0], 0.3e1) + 0.2e1 * theta * p[0] * pow(p[1], 0.2e1) - 0.20e2 * theta * p[0] * pow(p[2], 0.2e1) - pow(p[0], 0.4e1) + pow(p[0], 0.3e1) * p[2] - pow(p[0], 0.2e1) * pow(p[1], 0.2e1) - 0.2e1 * pow(p[0], 0.2e1) * pow(p[2], 0.2e1) + p[0] * pow(p[1], 0.2e1) * p[2] + p[0] * pow(p[2], 0.3e1) - pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - pow(p[2], 0.4e1) - 0.10e2 * theta * pow(p[0], 0.2e1) + 0.10e2 * theta * p[0] * p[2] - 0.10e2 * theta * pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1)
    , -(-pow(p[0], 0.3e1) * pow(p[2], 0.3e1) - p[0] * pow(p[1], 0.2e1) * pow(p[2], 0.3e1) - p[0] * pow(p[2], 0.5e1) + 0.14e2 * theta * pow(p[0], 0.3e1) * p[2] + 0.14e2 * theta * p[0] * pow(p[1], 0.2e1) * p[2] - 0.8e1 * theta * p[0] * pow(p[2], 0.3e1) - pow(p[0], 0.5e1) - 0.2e1 * pow(p[0], 0.3e1) * pow(p[1], 0.2e1) - pow(p[0], 0.3e1) * pow(p[2], 0.2e1) + pow(p[0], 0.2e1) * pow(p[1], 0.2e1) * p[2] - p[0] * pow(p[1], 0.4e1) - p[0] * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) + pow(p[1], 0.4e1) * p[2] + pow(p[1], 0.2e1) * pow(p[2], 0.3e1) - 0.5e1 * theta * pow(p[0], 0.3e1) - 0.5e1 * theta * p[0] * pow(p[1], 0.2e1) + 0.5e1 * theta * p[0] * pow(p[2], 0.2e1) + 0.10e2 * theta * pow(p[1], 0.2e1) * p[2]) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1)
    );
    return v;
}

double Test_A_plus_M_pSolScalarFun1 (const DROPS::Point3DCL& p, double)
{
    return 2*((-0.4e1 * p[0] * pow(p[2], 0.2e1) - pow(p[0], 0.2e1) + 0.3e1 * p[0] * p[2] + 0.2e1 * pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1));
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun1 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL w (0,0,0);
    DROPS::Point3DCL v(-(pow(p[1], 0.2e1) * pow(p[2], 0.2e1) + pow(p[2], 0.4e1) + pow(p[0], 0.2e1) * p[2] + p[0] * pow(p[1], 0.2e1)) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1))
    , p[1] * (p[0] * pow(p[2], 0.2e1) + pow(p[0], 0.2e1) - p[0] * p[2] + pow(p[2], 0.2e1)) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1))
    , (p[0] * pow(p[2], 0.3e1) + pow(p[0], 0.3e1) + p[0] * pow(p[1], 0.2e1) - pow(p[1], 0.2e1) * p[2]) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1))
    );
    if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
        v=w;
    }
    return v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun1_Gradient1 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL w (0,0,0);
    DROPS::Point3DCL v((0.2e1 * p[0] * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) + 0.2e1 * p[0] * pow(p[2], 0.4e1) + pow(p[0], 0.2e1) * pow(p[1], 0.2e1) - 0.2e1 * p[0] * pow(p[1], 0.2e1) * p[2] - 0.2e1 * p[0] * pow(p[2], 0.3e1) - pow(p[1], 0.4e1) - pow(p[1], 0.2e1) * pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1)
    , -0.2e1 * p[1] * (p[0] * pow(p[2], 0.2e1) + pow(p[0], 0.2e1) - p[0] * p[2] + pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * p[0]
    , -(0.2e1 * pow(p[0], 0.2e1) * pow(p[1], 0.2e1) * p[2] + 0.4e1 * pow(p[0], 0.2e1) * pow(p[2], 0.3e1) + 0.2e1 * pow(p[1], 0.4e1) * p[2] + 0.4e1 * pow(p[1], 0.2e1) * pow(p[2], 0.3e1) + 0.2e1 * pow(p[2], 0.5e1) + pow(p[0], 0.4e1) + pow(p[0], 0.2e1) * pow(p[1], 0.2e1) - pow(p[0], 0.2e1) * pow(p[2], 0.2e1) - 0.2e1 * p[0] * pow(p[1], 0.2e1) * p[2]) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1)
    );
    if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
        v=w;
    }
    return v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun1_Gradient2 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL w (0,0,0);
    DROPS::Point3DCL v(-p[1] * (pow(p[0], 0.2e1) * pow(p[2], 0.2e1) - pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - pow(p[2], 0.4e1) - pow(p[0], 0.2e1) * p[2] - 0.2e1 * p[0] * pow(p[1], 0.2e1) + pow(p[1], 0.2e1) * p[2] + pow(p[2], 0.3e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1)
    , (p[0] * pow(p[2], 0.2e1) + pow(p[0], 0.2e1) - p[0] * p[2] + pow(p[2], 0.2e1)) * (pow(p[0], 0.2e1) - pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1)
    , p[1] * (0.2e1 * pow(p[0], 0.3e1) * p[2] + 0.2e1 * p[0] * pow(p[1], 0.2e1) * p[2] - pow(p[0], 0.3e1) - p[0] * pow(p[1], 0.2e1) + p[0] * pow(p[2], 0.2e1) + 0.2e1 * pow(p[1], 0.2e1) * p[2]) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1)
    );
    if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
        v=w;
    }
    return v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun1_Gradient3 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL w (0,0,0);
    DROPS::Point3DCL v((-pow(p[0], 0.2e1) * pow(p[2], 0.3e1) + pow(p[1], 0.2e1) * pow(p[2], 0.3e1) + pow(p[2], 0.5e1) + pow(p[0], 0.4e1) + 0.2e1 * pow(p[0], 0.2e1) * pow(p[1], 0.2e1) + 0.3e1 * pow(p[0], 0.2e1) * pow(p[2], 0.2e1) + 0.2e1 * p[0] * pow(p[1], 0.2e1) * p[2] + pow(p[1], 0.4e1) + pow(p[1], 0.2e1) * pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1)
    , -0.2e1 * p[1] * (p[0] * pow(p[2], 0.2e1) + pow(p[0], 0.2e1) - p[0] * p[2] + pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * p[2]
    , (0.3e1 * pow(p[0], 0.3e1) * pow(p[2], 0.2e1) + 0.3e1 * p[0] * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) + p[0] * pow(p[2], 0.4e1) - 0.2e1 * pow(p[0], 0.3e1) * p[2] - pow(p[0], 0.2e1) * pow(p[1], 0.2e1) - 0.2e1 * p[0] * pow(p[1], 0.2e1) * p[2] - pow(p[1], 0.4e1) + pow(p[1], 0.2e1) * pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1)
    );
    if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
        v=w;
    }
    return v;
}

DROPS::Point3DCL Test_A_plus_M_RhsVectorFun2 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(-p[1], p[0], 0);
    return v;
}

double Test_A_plus_M_pSolScalarFun2 (const DROPS::Point3DCL&, double)
{
    return 0.;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun2 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(-p[1], p[0], 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun2_Gradient1 (const DROPS::Point3DCL&, double)
{
    DROPS::Point3DCL v(0, -1, 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun2_Gradient2 (const DROPS::Point3DCL&, double)
{
    DROPS::Point3DCL v(1, 0, 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun2_Gradient3 (const DROPS::Point3DCL&, double)
{
    DROPS::Point3DCL v(0, 0, 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_RhsVectorFun3 (const DROPS::Point3DCL& x, double)
{
    double theta = 1;
    DROPS::Point3DCL v(-(pow(x[0], 0.3e1) * x[1] + pow(x[0], 0.3e1) * x[2] - pow(x[0], 0.2e1) * pow(x[1], 0.2e1) - pow(x[0], 0.2e1) * pow(x[2], 0.2e1) + x[0] * pow(x[1], 0.3e1) + x[0] * pow(x[1], 0.2e1) * x[2] + x[0] * x[1] * pow(x[2], 0.2e1) + x[0] * pow(x[2], 0.3e1) - pow(x[1], 0.4e1) - 0.2e1 * pow(x[1], 0.2e1) * pow(x[2], 0.2e1) - pow(x[2], 0.4e1) + 0.2e1 * theta * x[0] * x[1] + 0.2e1 * theta * x[0] * x[2] - 0.2e1 * theta * pow(x[1], 0.2e1) - 0.2e1 * theta * pow(x[2], 0.2e1)) / (pow(x[0], 0.2e1) + pow(x[1], 0.2e1) + pow(x[2], 0.2e1))
    , (pow(x[0], 0.4e1) - pow(x[0], 0.3e1) * x[1] + pow(x[0], 0.2e1) * pow(x[1], 0.2e1) - pow(x[0], 0.2e1) * x[1] * x[2] + 0.2e1 * pow(x[0], 0.2e1) * pow(x[2], 0.2e1) - x[0] * pow(x[1], 0.3e1) - x[0] * x[1] * pow(x[2], 0.2e1) - pow(x[1], 0.3e1) * x[2] + pow(x[1], 0.2e1) * pow(x[2], 0.2e1) - x[1] * pow(x[2], 0.3e1) + pow(x[2], 0.4e1) + 0.2e1 * theta * pow(x[0], 0.2e1) - 0.2e1 * theta * x[0] * x[1] - 0.2e1 * theta * x[1] * x[2] + 0.2e1 * theta * pow(x[2], 0.2e1)) / (pow(x[0], 0.2e1) + pow(x[1], 0.2e1) + pow(x[2], 0.2e1))
    , (pow(x[0], 0.4e1) - pow(x[0], 0.3e1) * x[2] + 0.2e1 * pow(x[0], 0.2e1) * pow(x[1], 0.2e1) - pow(x[0], 0.2e1) * x[1] * x[2] + pow(x[0], 0.2e1) * pow(x[2], 0.2e1) - x[0] * pow(x[1], 0.2e1) * x[2] - x[0] * pow(x[2], 0.3e1) + pow(x[1], 0.4e1) - pow(x[1], 0.3e1) * x[2] + pow(x[1], 0.2e1) * pow(x[2], 0.2e1) - x[1] * pow(x[2], 0.3e1) + 0.2e1 * theta * pow(x[0], 0.2e1) - 0.2e1 * theta * x[0] * x[2] + 0.2e1 * theta * pow(x[1], 0.2e1) - 0.2e1 * theta * x[1] * x[2]) / (pow(x[0], 0.2e1) + pow(x[1], 0.2e1) + pow(x[2], 0.2e1))
    );
    return v;
}

double Test_A_plus_M_pSolScalarFun3 (const DROPS::Point3DCL& p, double)
{
    return 0.4e1 * (p[0] + p[1] + p[2]) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1);
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun3 (const DROPS::Point3DCL& x, double)
{
    DROPS::Point3DCL v(-x[0] * x[1] - x[0] * x[2] + pow(x[1], 0.2e1) + pow(x[2], 0.2e1)
    , pow(x[0], 0.2e1) - x[0] * x[1] - x[1] * x[2] + pow(x[2], 0.2e1)
    , pow(x[0], 0.2e1) - x[0] * x[2] + pow(x[1], 0.2e1) - x[1] * x[2]);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun3_Gradient1 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(-(double) p[1] - (double) p[2], -p[0] + 2 * p[1], -p[0] + 2 * p[2]);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun3_Gradient2 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(2 * p[0] - p[1], -p[0] - p[2], -p[1] + 2 * p[2]);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun3_Gradient3 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(2 * p[0] - p[2], 2 * p[1] - p[2], -p[0] - p[1]);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_RhsVectorFun5 (const DROPS::Point3DCL& p, double)
{
    double rad = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    const double PI = M_PI;
    double theta = acos(p[2]/rad);
    double phi;
    if(p[0] == 0){
        if(p[1] > 0) phi = PI/2;
        if(p[1] < 0) phi = 3*PI/2;
        if(p[1] == 0) phi = 0;
    }
    if(p[0] > 0) {
        if(p[1] > 0) phi = atan(p[1]/p[0]);
        if(p[1] < 0) phi = -atan(std::abs(p[1])/p[0]);
        if(p[1] == 0) phi = 0;
    }
    if(p[0] < 0) {
        if(p[1] > 0) phi = atan(std::abs(p[0])/p[1]) + PI/2;
        if(p[1] < 0) phi = atan(std::abs(p[1])/std::abs(p[0])) + PI;
        if(p[1] == 0) phi = PI;
    }

    DROPS::Point3DCL v(-(pow(cos(phi), 0.2e1) * pow(sin(theta), 0.4e1) * pow(sin(phi), 0.2e1) * pow(cos(theta), 0.2e1) + pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.4e1) + pow(sin(phi), 0.4e1) * pow(sin(theta), 0.4e1) * pow(cos(theta), 0.2e1) + 0.2e1 * pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.4e1) + pow(cos(theta), 0.6e1) + pow(cos(phi), 0.4e1) * pow(sin(theta), 0.4e1) * cos(theta) + pow(cos(phi), 0.3e1) * pow(sin(theta), 0.5e1) * pow(sin(phi), 0.2e1) + pow(cos(phi), 0.2e1) * pow(sin(theta), 0.4e1) * pow(sin(phi), 0.2e1) * cos(theta) + pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.3e1) + cos(phi) * pow(sin(theta), 0.5e1) * pow(sin(phi), 0.4e1) + cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1) * pow(cos(theta), 0.2e1) - 0.2e1 * pow(cos(phi), 0.2e1) * pow(sin(theta), 0.4e1) * pow(sin(phi), 0.2e1) - 0.14e2 * pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.2e1) - 0.2e1 * pow(sin(phi), 0.4e1) * pow(sin(theta), 0.4e1) + 0.6e1 * pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.2e1) + 0.8e1 * pow(cos(theta), 0.4e1) + 0.5e1 * pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) * cos(theta) + 0.10e2 * cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1) - 0.5e1 * pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) * cos(theta) - 0.5e1 * pow(cos(theta), 0.3e1)) * pow(pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(cos(theta), 0.2e1), -0.2e1)
    , sin(phi) * sin(theta) * (pow(cos(phi), 0.3e1) * pow(sin(theta), 0.3e1) * pow(cos(theta), 0.2e1) + cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1) * pow(cos(theta), 0.2e1) + cos(phi) * sin(theta) * pow(cos(theta), 0.4e1) + pow(cos(phi), 0.4e1) * pow(sin(theta), 0.4e1) - pow(cos(phi), 0.3e1) * pow(sin(theta), 0.3e1) * cos(theta) + pow(cos(phi), 0.2e1) * pow(sin(theta), 0.4e1) * pow(sin(phi), 0.2e1) + 0.2e1 * pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.2e1) - cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1) * cos(theta) - cos(phi) * sin(theta) * pow(cos(theta), 0.3e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.2e1) + pow(cos(theta), 0.4e1) - 0.2e1 * pow(cos(phi), 0.3e1) * pow(sin(theta), 0.3e1) - 0.2e1 * cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1) + 0.20e2 * cos(phi) * sin(theta) * pow(cos(theta), 0.2e1) + 0.10e2 * pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) - 0.10e2 * cos(theta) * cos(phi) * sin(theta) + 0.10e2 * pow(cos(theta), 0.2e1)) * pow(pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(cos(theta), 0.2e1), -0.2e1)
    , (pow(cos(phi), 0.3e1) * pow(sin(theta), 0.3e1) * pow(cos(theta), 0.3e1) + cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1) * pow(cos(theta), 0.3e1) + cos(phi) * sin(theta) * pow(cos(theta), 0.5e1) + pow(cos(phi), 0.5e1) * pow(sin(theta), 0.5e1) + 0.2e1 * pow(cos(phi), 0.3e1) * pow(sin(theta), 0.5e1) * pow(sin(phi), 0.2e1) + pow(cos(phi), 0.3e1) * pow(sin(theta), 0.3e1) * pow(cos(theta), 0.2e1) - pow(cos(phi), 0.2e1) * pow(sin(theta), 0.4e1) * pow(sin(phi), 0.2e1) * cos(theta) + cos(phi) * pow(sin(theta), 0.5e1) * pow(sin(phi), 0.4e1) + cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1) * pow(cos(theta), 0.2e1) - pow(sin(phi), 0.4e1) * pow(sin(theta), 0.4e1) * cos(theta) - pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.3e1) - 0.14e2 * pow(cos(phi), 0.3e1) * pow(sin(theta), 0.3e1) * cos(theta) - 0.14e2 * cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1) * cos(theta) + 0.8e1 * cos(phi) * sin(theta) * pow(cos(theta), 0.3e1) + 0.5e1 * pow(cos(phi), 0.3e1) * pow(sin(theta), 0.3e1) + 0.5e1 * cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1) - 0.5e1 * cos(phi) * sin(theta) * pow(cos(theta), 0.2e1) - 0.10e2 * pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) * cos(theta)) * pow(pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(cos(theta), 0.2e1), -0.2e1)
    );
    return v;
}

double Test_A_plus_M_pSolScalarFun5 (const DROPS::Point3DCL& p, double)
{
    double rad = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    const double PI = M_PI;
    double theta = acos(p[2]/rad);
    double phi;
    if(p[0] == 0){
        if(p[1] > 0) phi = PI/2;
        if(p[1] < 0) phi = 3*PI/2;
        if(p[1] == 0) phi = 0;
    }
    if(p[0] > 0) {
        if(p[1] > 0) phi = atan(p[1]/p[0]);
        if(p[1] < 0) phi = -atan(std::abs(p[1])/p[0]);
        if(p[1] == 0) phi = 0;
    }
    if(p[0] < 0) {
        if(p[1] > 0) phi = atan(std::abs(p[0])/p[1]) + PI/2;
        if(p[1] < 0) phi = atan(std::abs(p[1])/std::abs(p[0])) + PI;
        if(p[1] == 0) phi = PI;
    }
    return -1*(0.2e1 * (0.4e1 * cos(phi) * sin(theta) * pow(cos(theta), 0.2e1) + pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) - 0.3e1 * cos(theta) * cos(phi) * sin(theta) - 0.2e1 * pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(cos(theta), 0.2e1)) * pow(pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(cos(theta), 0.2e1), -0.3e1 / 0.2e1));
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun5 (const DROPS::Point3DCL& p, double)
{
    double rad = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    const double PI = M_PI;
    double theta = acos(p[2]/rad);
    double phi;
    if(p[0] == 0){
        if(p[1] > 0) phi = PI/2;
        if(p[1] < 0) phi = 3*PI/2;
        if(p[1] == 0) phi = 0;
    }
    if(p[0] > 0) {
        if(p[1] > 0) phi = atan(p[1]/p[0]);
        if(p[1] < 0) phi = -atan(std::abs(p[1])/p[0]);
        if(p[1] == 0) phi = 0;
    }
    if(p[0] < 0) {
        if(p[1] > 0) phi = atan(std::abs(p[0])/p[1]) + PI/2;
        if(p[1] < 0) phi = atan(std::abs(p[1])/std::abs(p[0])) + PI;
        if(p[1] == 0) phi = PI;
    }

    DROPS::Point3DCL v(-(pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.2e1) + pow(cos(theta), 0.4e1) + pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) * cos(theta) + cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1)) / (pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(cos(theta), 0.2e1))
    , sin(phi) * sin(theta) * (cos(phi) * sin(theta) * pow(cos(theta), 0.2e1) + pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) - cos(theta) * cos(phi) * sin(theta) + pow(cos(theta), 0.2e1)) / (pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(cos(theta), 0.2e1))
    , (cos(phi) * sin(theta) * pow(cos(theta), 0.3e1) + pow(cos(phi), 0.3e1) * pow(sin(theta), 0.3e1) + cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1) - pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) * cos(theta)) / (pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(cos(theta), 0.2e1))
    );
    return v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun5_Gradient1 (const DROPS::Point3DCL& p, double)
{
    double rad = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    const double PI = M_PI;
    double theta = acos(p[2]/rad);
    double phi;
    if(p[0] == 0){
        if(p[1] > 0) phi = PI/2;
        if(p[1] < 0) phi = 3*PI/2;
        if(p[1] == 0) phi = 0;
    }
    if(p[0] > 0) {
        if(p[1] > 0) phi = atan(p[1]/p[0]);
        if(p[1] < 0) phi = -atan(std::abs(p[1])/p[0]);
        if(p[1] == 0) phi = 0;
    }
    if(p[0] < 0) {
        if(p[1] > 0) phi = atan(std::abs(p[0])/p[1]) + PI/2;
        if(p[1] < 0) phi = atan(std::abs(p[1])/std::abs(p[0])) + PI;
        if(p[1] == 0) phi = PI;
    }
    DROPS::Point3DCL v((0.2e1 * cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1) * pow(cos(theta), 0.2e1) + 0.2e1 * cos(phi) * sin(theta) * pow(cos(theta), 0.4e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.4e1) * pow(cos(phi), 0.2e1) - 0.2e1 * cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1) * cos(theta) - 0.2e1 * cos(phi) * sin(theta) * pow(cos(theta), 0.3e1) - pow(sin(phi), 0.4e1) * pow(sin(theta), 0.4e1) - pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.2e1)) * pow(pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(cos(theta), 0.2e1), -0.2e1)
    , -0.2e1 * sin(phi) * pow(sin(theta), 0.2e1) * (cos(phi) * sin(theta) * pow(cos(theta), 0.2e1) + pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) - cos(phi) * sin(theta) * cos(theta) + pow(cos(theta), 0.2e1)) * pow(pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(cos(theta), 0.2e1), -0.2e1) * cos(phi)
    , -(0.2e1 * pow(cos(phi), 0.2e1) * pow(sin(theta), 0.4e1) * pow(sin(phi), 0.2e1) * cos(theta) + 0.4e1 * pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.3e1) + 0.2e1 * pow(sin(phi), 0.4e1) * pow(sin(theta), 0.4e1) * cos(theta) + 0.4e1 * pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.3e1) + 0.2e1 * pow(cos(theta), 0.5e1) + pow(cos(phi), 0.4e1) * pow(sin(theta), 0.4e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.4e1) * pow(cos(phi), 0.2e1) - pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.2e1) - 0.2e1 * cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1) * cos(theta)) * pow(pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(cos(theta), 0.2e1), -0.2e1)
    );
    return v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun5_Gradient2 (const DROPS::Point3DCL& p, double)
{
    double rad = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    const double PI = M_PI;
    double theta = acos(p[2]/rad);
    double phi;
    if(p[0] == 0){
        if(p[1] > 0) phi = PI/2;
        if(p[1] < 0) phi = 3*PI/2;
        if(p[1] == 0) phi = 0;
    }
    if(p[0] > 0) {
        if(p[1] > 0) phi = atan(p[1]/p[0]);
        if(p[1] < 0) phi = -atan(std::abs(p[1])/p[0]);
        if(p[1] == 0) phi = 0;
    }
    if(p[0] < 0) {
        if(p[1] > 0) phi = atan(std::abs(p[0])/p[1]) + PI/2;
        if(p[1] < 0) phi = atan(std::abs(p[1])/std::abs(p[0])) + PI;
        if(p[1] == 0) phi = PI;
    }
    DROPS::Point3DCL v(-sin(phi) * sin(theta) * (pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.2e1) - pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.2e1) - pow(cos(theta), 0.4e1) - pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) * cos(theta) - 0.2e1 * cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) * cos(theta) + pow(cos(theta), 0.3e1)) * pow(pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(cos(theta), 0.2e1), -0.2e1)
    , (cos(phi) * sin(theta) * pow(cos(theta), 0.2e1) + pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) - cos(phi) * sin(theta) * cos(theta) + pow(cos(theta), 0.2e1)) * (pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) - pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(cos(theta), 0.2e1)) * pow(pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(cos(theta), 0.2e1), -0.2e1)
    , sin(phi) * sin(theta) * (0.2e1 * pow(cos(phi), 0.3e1) * pow(sin(theta), 0.3e1) * cos(theta) + 0.2e1 * cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1) * cos(theta) - pow(cos(phi), 0.3e1) * pow(sin(theta), 0.3e1) - cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1) + cos(phi) * sin(theta) * pow(cos(theta), 0.2e1) + 0.2e1 * pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) * cos(theta)) * pow(pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(cos(theta), 0.2e1), -0.2e1)
    );
    return v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun5_Gradient3 (const DROPS::Point3DCL& p, double)
{
    double rad = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    const double PI = M_PI;
    double theta = acos(p[2]/rad);
    double phi;
    if(p[0] == 0){
        if(p[1] > 0) phi = PI/2;
        if(p[1] < 0) phi = 3*PI/2;
        if(p[1] == 0) phi = 0;
    }
    if(p[0] > 0) {
        if(p[1] > 0) phi = atan(p[1]/p[0]);
        if(p[1] < 0) phi = -atan(std::abs(p[1])/p[0]);
        if(p[1] == 0) phi = 0;
    }
    if(p[0] < 0) {
        if(p[1] > 0) phi = atan(std::abs(p[0])/p[1]) + PI/2;
        if(p[1] < 0) phi = atan(std::abs(p[1])/std::abs(p[0])) + PI;
        if(p[1] == 0) phi = PI;
    }
    DROPS::Point3DCL v((-pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.3e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.3e1) + pow(cos(theta), 0.5e1) + pow(cos(phi), 0.4e1) * pow(sin(theta), 0.4e1) + 0.2e1 * pow(sin(phi), 0.2e1) * pow(sin(theta), 0.4e1) * pow(cos(phi), 0.2e1) + 0.3e1 * pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.2e1) + 0.2e1 * cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1) * cos(theta) + pow(sin(phi), 0.4e1) * pow(sin(theta), 0.4e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.2e1)) * pow(pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(cos(theta), 0.2e1), -0.2e1)
    , -0.2e1 * sin(phi) * sin(theta) * (cos(phi) * sin(theta) * pow(cos(theta), 0.2e1) + pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) - cos(phi) * sin(theta) * cos(theta) + pow(cos(theta), 0.2e1)) * pow(pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(cos(theta), 0.2e1), -0.2e1) * cos(theta)
    , (0.3e1 * pow(cos(phi), 0.3e1) * pow(sin(theta), 0.3e1) * pow(cos(theta), 0.2e1) + 0.3e1 * cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1) * pow(cos(theta), 0.2e1) + cos(phi) * sin(theta) * pow(cos(theta), 0.4e1) - 0.2e1 * pow(cos(phi), 0.3e1) * pow(sin(theta), 0.3e1) * cos(theta) - pow(sin(phi), 0.2e1) * pow(sin(theta), 0.4e1) * pow(cos(phi), 0.2e1) - 0.2e1 * cos(phi) * pow(sin(theta), 0.3e1) * pow(sin(phi), 0.2e1) * cos(theta) - pow(sin(phi), 0.4e1) * pow(sin(theta), 0.4e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) * pow(cos(theta), 0.2e1)) * pow(pow(cos(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(sin(phi), 0.2e1) * pow(sin(theta), 0.2e1) + pow(cos(theta), 0.2e1), -0.2e1)
    );
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_RhsVectorFun1 (const DROPS::Point3DCL& p, double)
{
    double theta = 1;
    DROPS::Point3DCL v (-cos(0.3141592654e1 * p[0]) * cos(0.3141592654e1 * p[1]) * 0.3141592654e1 * 0.3141592654e1 * theta + 0.3e1 * 0.3141592654e1 * 0.3141592654e1 * sin(0.3141592654e1 * p[0]) * sin(0.3141592654e1 * p[1]) * theta + sin(0.3141592654e1 * p[0]) * sin(0.3141592654e1 * p[1])
    , -cos(0.3141592654e1 * p[0]) * cos(0.3141592654e1 * p[1]) * 0.3141592654e1 * 0.3141592654e1 * theta + 0.3e1 * 0.3141592654e1 * 0.3141592654e1 * sin(0.3141592654e1 * p[0]) * sin(0.3141592654e1 * p[1]) * theta + sin(0.3141592654e1 * p[0]) * sin(0.3141592654e1 * p[1])
    , 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_vSolVectorFun1 (const DROPS::Point3DCL& p, double)
{
    const double PI = 3.141592653589793;
    DROPS::Point3DCL v (sin(PI*p[0])*sin(PI*p[1]), sin(PI*p[0])*sin(PI*p[1]), 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_vSolVectorFun1_Gradient1 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (cos(0.3141592654e1 * p[0]) * 0.3141592654e1 * sin(0.3141592654e1 * p[1]), sin(0.3141592654e1 * p[0]) * cos(0.3141592654e1 * p[1]) * 0.3141592654e1, 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_vSolVectorFun1_Gradient2 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (cos(0.3141592654e1 * p[0]) * 0.3141592654e1 * sin(0.3141592654e1 * p[1]), sin(0.3141592654e1 * p[0]) * cos(0.3141592654e1 * p[1]) * 0.3141592654e1, 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_vSolVectorFun1_Gradient3 (const DROPS::Point3DCL&, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0, 0, 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_RhsVectorFun2 (const DROPS::Point3DCL& p, double)
{
    double theta = 1;
    DROPS::Point3DCL v (-cos(0.3141592654e1 * p[1] / 0.2e1) * cos(0.3141592654e1 * p[0] / 0.2e1) * 0.3141592654e1 * 0.3141592654e1 * theta / 0.4e1 + 0.3e1 / 0.4e1 * sin(0.3141592654e1 * p[0] / 0.2e1) * sin(0.3141592654e1 * p[1] / 0.2e1) * 0.3141592654e1 * 0.3141592654e1 * theta + sin(0.3141592654e1 * p[0] / 0.2e1) * sin(0.3141592654e1 * p[1] / 0.2e1)
    , -cos(0.3141592654e1 * p[1] / 0.2e1) * cos(0.3141592654e1 * p[0] / 0.2e1) * 0.3141592654e1 * 0.3141592654e1 * theta / 0.4e1 + 0.3e1 / 0.4e1 * sin(0.3141592654e1 * p[0] / 0.2e1) * sin(0.3141592654e1 * p[1] / 0.2e1) * 0.3141592654e1 * 0.3141592654e1 * theta + sin(0.3141592654e1 * p[0] / 0.2e1) * sin(0.3141592654e1 * p[1] / 0.2e1)
    , 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_vSolVectorFun2 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (sin(0.3141592654e1 * p[0] / 0.2e1) * sin(0.3141592654e1 * p[1] / 0.2e1), sin(0.3141592654e1 * p[0] / 0.2e1) * sin(0.3141592654e1 * p[1] / 0.2e1), 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_vSolVectorFun2_Gradient1 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (cos(0.3141592654e1 * p[0] / 0.2e1) * 0.3141592654e1 * sin(0.3141592654e1 * p[1] / 0.2e1) / 0.2e1, sin(0.3141592654e1 * p[0] / 0.2e1) * cos(0.3141592654e1 * p[1] / 0.2e1) * 0.3141592654e1 / 0.2e1, 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_vSolVectorFun2_Gradient2 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (cos(0.3141592654e1 * p[0] / 0.2e1) * 0.3141592654e1 * sin(0.3141592654e1 * p[1] / 0.2e1) / 0.2e1, sin(0.3141592654e1 * p[0] / 0.2e1) * cos(0.3141592654e1 * p[1] / 0.2e1) * 0.3141592654e1 / 0.2e1, 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_vSolVectorFun2_Gradient3 (const DROPS::Point3DCL&, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0, 0, 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_RhsVectorFun3 (const DROPS::Point3DCL& p, double)
{
    double theta = 1;
    DROPS::Point3DCL v (0.3e1 / 0.16e2 * sin(0.3141592654e1 * p[0] / 0.4e1) * sin(0.3141592654e1 * p[1] / 0.4e1) * 0.3141592654e1 * 0.3141592654e1 * theta - cos(0.3141592654e1 * p[0] / 0.4e1) * cos(0.3141592654e1 * p[1] / 0.4e1) * 0.3141592654e1 * 0.3141592654e1 * theta / 0.16e2 + sin(0.3141592654e1 * p[0] / 0.4e1) * sin(0.3141592654e1 * p[1] / 0.4e1)
    , 0.3e1 / 0.16e2 * sin(0.3141592654e1 * p[0] / 0.4e1) * sin(0.3141592654e1 * p[1] / 0.4e1) * 0.3141592654e1 * 0.3141592654e1 * theta - cos(0.3141592654e1 * p[0] / 0.4e1) * cos(0.3141592654e1 * p[1] / 0.4e1) * 0.3141592654e1 * 0.3141592654e1 * theta / 0.16e2 + sin(0.3141592654e1 * p[0] / 0.4e1) * sin(0.3141592654e1 * p[1] / 0.4e1)
    , 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_vSolVectorFun3 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (sin(0.3141592654e1 * p[0] / 0.4e1) * sin(0.3141592654e1 * p[1] / 0.4e1), sin(0.3141592654e1 * p[0] / 0.4e1) * sin(0.3141592654e1 * p[1] / 0.4e1), 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_vSolVectorFun3_Gradient1 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (cos(0.3141592654e1 * p[0] / 0.4e1) * 0.3141592654e1 * sin(0.3141592654e1 * p[1] / 0.4e1) / 0.4e1, sin(0.3141592654e1 * p[0] / 0.4e1) * cos(0.3141592654e1 * p[1] / 0.4e1) * 0.3141592654e1 / 0.4e1, 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_vSolVectorFun3_Gradient2 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (cos(0.3141592654e1 * p[0] / 0.4e1) * 0.3141592654e1 * sin(0.3141592654e1 * p[1] / 0.4e1) / 0.4e1, sin(0.3141592654e1 * p[0] / 0.4e1) * cos(0.3141592654e1 * p[1] / 0.4e1) * 0.3141592654e1 / 0.4e1, 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_vSolVectorFun3_Gradient3 (const DROPS::Point3DCL&, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0, 0, 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_RhsVectorFun4 (const DROPS::Point3DCL& p, double)
{
    double theta = 1;
    DROPS::Point3DCL v (pow(p[0], 0.2e1) * pow(p[1], 0.2e1) - 0.2e1 * theta * pow(p[0], 0.2e1) - 0.4e1 * theta * p[0] * p[1] - 0.4e1 * theta * pow(p[1], 0.2e1) - 0.4e1 * pow(p[0], 0.2e1) - 0.4e1 * pow(p[1], 0.2e1) + 0.24e2 * theta + 0.16e2
    , pow(p[0], 0.2e1) * pow(p[1], 0.2e1) - 0.4e1 * theta * pow(p[0], 0.2e1) - 0.4e1 * theta * p[0] * p[1] - 0.2e1 * theta * pow(p[1], 0.2e1) - 0.4e1 * pow(p[0], 0.2e1) - 0.4e1 * pow(p[1], 0.2e1) + 0.24e2 * theta + 0.16e2
    , 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_vSolVectorFun4 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v ((pow(p[0], 0.2e1) - 0.4e1) * (pow(p[1], 0.2e1) - 0.4e1), (pow(p[0], 0.2e1) - 0.4e1) * (pow(p[1], 0.2e1) - 0.4e1), 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_vSolVectorFun4_Gradient1 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0.2e1 * p[0] * (pow(p[1], 0.2e1) - 0.4e1), 0.2e1 * (pow(p[0], 0.2e1) - 0.4e1) * p[1], 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_vSolVectorFun4_Gradient2 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0.2e1 * p[0] * (pow(p[1], 0.2e1) - 0.4e1), 0.2e1 * (pow(p[0], 0.2e1) - 0.4e1) * p[1], 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_xy_plane_vSolVectorFun4_Gradient3 (const DROPS::Point3DCL&, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0, 0, 0);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_RhsVectorFun1 (const DROPS::Point3DCL& p, double)
{
    double theta = 1;
    DROPS::Point3DCL v (-0.4e1 / 0.5e1 * cos(0.3141592654e1 * p[0]) * cos(0.3141592654e1 * p[2]) * 0.3141592654e1 * 0.3141592654e1 * theta + 0.52e2 / 0.25e2 * sin(0.3141592654e1 * p[2]) * sin(0.3141592654e1 * p[0]) * 0.3141592654e1 * 0.3141592654e1 * theta + 0.4e1 / 0.5e1 * sin(0.3141592654e1 * p[0]) * sin(0.3141592654e1 * p[2])
    , -0.2e1 / 0.5e1 * cos(0.3141592654e1 * p[0]) * cos(0.3141592654e1 * p[2]) * 0.3141592654e1 * 0.3141592654e1 * theta + 0.26e2 / 0.25e2 * sin(0.3141592654e1 * p[2]) * sin(0.3141592654e1 * p[0]) * 0.3141592654e1 * 0.3141592654e1 * theta + 0.2e1 / 0.5e1 * sin(0.3141592654e1 * p[0]) * sin(0.3141592654e1 * p[2])
    , -0.4e1 / 0.5e1 * cos(0.3141592654e1 * p[0]) * cos(0.3141592654e1 * p[2]) * 0.3141592654e1 * 0.3141592654e1 * theta + 0.14e2 / 0.5e1 * sin(0.3141592654e1 * p[2]) * sin(0.3141592654e1 * p[0]) * 0.3141592654e1 * 0.3141592654e1 * theta + sin(0.3141592654e1 * p[0]) * sin(0.3141592654e1 * p[2])
    );
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_vSolVectorFun1 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v (0.4e1 / 0.5e1 * sin(0.3141592654e1 * p[0]) * sin(0.3141592654e1 * p[2]), 0.2e1 / 0.5e1 * sin(0.3141592654e1 * p[0]) * sin(0.3141592654e1 * p[2]), sin(0.3141592654e1 * p[0]) * sin(0.3141592654e1 * p[2]));
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_vSolVectorFun1_Gradient1 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0.4e1 / 0.5e1 * cos(0.3141592654e1 * p[0]) * 0.3141592654e1 * sin(0.3141592654e1 * p[2]), 0, 0.4e1 / 0.5e1 * sin(0.3141592654e1 * p[0]) * cos(0.3141592654e1 * p[2]) * 0.3141592654e1);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_vSolVectorFun1_Gradient2 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0.2e1 / 0.5e1 * cos(0.3141592654e1 * p[0]) * 0.3141592654e1 * sin(0.3141592654e1 * p[2]), 0, 0.2e1 / 0.5e1 * sin(0.3141592654e1 * p[0]) * cos(0.3141592654e1 * p[2]) * 0.3141592654e1);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_vSolVectorFun1_Gradient3 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (cos(0.3141592654e1 * p[0]) * 0.3141592654e1 * sin(0.3141592654e1 * p[2]), 0, sin(0.3141592654e1 * p[0]) * cos(0.3141592654e1 * p[2]) * 0.3141592654e1);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_RhsVectorFun2 (const DROPS::Point3DCL& p, double)
{
    double theta = 1;
    DROPS::Point3DCL v (-0.64e2 / 0.25e2 * theta * pow(p[2], 0.2e1) + 0.416e3 / 0.25e2 * theta - 0.8e1 / 0.5e1 * theta * pow(p[0], 0.2e1) - 0.16e2 / 0.5e1 * theta * p[0] * p[2] + 0.4e1 / 0.5e1 * pow(p[0], 0.2e1) * pow(p[2], 0.2e1) - 0.16e2 / 0.5e1 * pow(p[0], 0.2e1) - 0.16e2 / 0.5e1 * pow(p[2], 0.2e1) + 0.64e2 / 0.5e1
    , -0.32e2 / 0.25e2 * theta * pow(p[2], 0.2e1) + 0.208e3 / 0.25e2 * theta - 0.4e1 / 0.5e1 * theta * pow(p[0], 0.2e1) - 0.8e1 / 0.5e1 * theta * p[0] * p[2] + 0.2e1 / 0.5e1 * pow(p[0], 0.2e1) * pow(p[2], 0.2e1) - 0.8e1 / 0.5e1 * pow(p[0], 0.2e1) - 0.8e1 / 0.5e1 * pow(p[2], 0.2e1) + 0.32e2 / 0.5e1
    , -0.16e2 / 0.5e1 * theta * p[0] * p[2] - 0.8e1 / 0.5e1 * theta * pow(p[2], 0.2e1) + 0.112e3 / 0.5e1 * theta - 0.4e1 * theta * pow(p[0], 0.2e1) + pow(p[0], 0.2e1) * pow(p[2], 0.2e1) - 0.4e1 * pow(p[0], 0.2e1) - 0.4e1 * pow(p[2], 0.2e1) + 0.16e2
    );
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_vSolVectorFun2 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0.4e1 / 0.5e1 * (pow(p[0], 0.2e1) - 0.4e1) * (pow(p[2], 0.2e1) - 0.4e1), 0.2e1 / 0.5e1 * (pow(p[0], 0.2e1) - 0.4e1) * (pow(p[2], 0.2e1) - 0.4e1), (pow(p[0], 0.2e1) - 0.4e1) * (pow(p[2], 0.2e1) - 0.4e1));
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_vSolVectorFun2_Gradient1 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0.8e1 / 0.5e1 * p[0] * (pow(p[2], 0.2e1) - 0.4e1), 0, 0.8e1 / 0.5e1 * (pow(p[0], 0.2e1) - 0.4e1) * p[2]);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_vSolVectorFun2_Gradient2 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0.4e1 / 0.5e1 * p[0] * (pow(p[2], 0.2e1) - 0.4e1), 0, 0.4e1 / 0.5e1 * p[0] * (pow(p[2], 0.2e1) - 0.4e1));
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_vSolVectorFun2_Gradient3 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0.2e1 * p[0] * (pow(p[2], 0.2e1) - 0.4e1), 0, 0.2e1 * (pow(p[0], 0.2e1) - 0.4e1) * p[2]);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_xy_RhsVectorFun1 (const DROPS::Point3DCL& p, double)
{
    double theta = 1;
    DROPS::Point3DCL v (0.196e3 / 0.289e3 * sin(0.3141592654e1 * p[1]) * 0.3141592654e1 * 0.3141592654e1 * sin(0.3141592654e1 * p[2]) * theta - 0.4e1 / 0.17e2 * theta * cos(0.3141592654e1 * p[1]) * 0.3141592654e1 * 0.3141592654e1 * cos(0.3141592654e1 * p[2]) + 0.4e1 / 0.17e2 * sin(0.3141592654e1 * p[1]) * sin(0.3141592654e1 * p[2])
    , 0.50e2 / 0.17e2 * sin(0.3141592654e1 * p[1]) * 0.3141592654e1 * 0.3141592654e1 * sin(0.3141592654e1 * p[2]) * theta - 0.16e2 / 0.17e2 * theta * cos(0.3141592654e1 * p[1]) * 0.3141592654e1 * 0.3141592654e1 * cos(0.3141592654e1 * p[2]) + sin(0.3141592654e1 * p[1]) * sin(0.3141592654e1 * p[2])
    , 0.784e3 / 0.289e3 * sin(0.3141592654e1 * p[1]) * 0.3141592654e1 * 0.3141592654e1 * sin(0.3141592654e1 * p[2]) * theta - 0.16e2 / 0.17e2 * theta * cos(0.3141592654e1 * p[1]) * 0.3141592654e1 * 0.3141592654e1 * cos(0.3141592654e1 * p[2]) + 0.16e2 / 0.17e2 * sin(0.3141592654e1 * p[1]) * sin(0.3141592654e1 * p[2])
    );
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_xy_vSolVectorFun1 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v (0.4e1 / 0.17e2 * sin(0.3141592654e1 * p[1]) * sin(0.3141592654e1 * p[2]), sin(0.3141592654e1 * p[1]) * sin(0.3141592654e1 * p[2]), 0.16e2 / 0.17e2 * sin(0.3141592654e1 * p[1]) * sin(0.3141592654e1 * p[2]));
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_xy_vSolVectorFun1_Gradient1 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0, 0.4e1 / 0.17e2 * cos(0.3141592654e1 * p[1]) * 0.3141592654e1 * sin(0.3141592654e1 * p[2]), 0.4e1 / 0.17e2 * sin(0.3141592654e1 * p[1]) * cos(0.3141592654e1 * p[2]) * 0.3141592654e1);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_xy_vSolVectorFun1_Gradient2 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0, cos(0.3141592654e1 * p[1]) * 0.3141592654e1 * sin(0.3141592654e1 * p[2]), sin(0.3141592654e1 * p[1]) * cos(0.3141592654e1 * p[2]) * 0.3141592654e1);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_xy_vSolVectorFun1_Gradient3 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0, 0.16e2 / 0.17e2 * cos(0.3141592654e1 * p[1]) * 0.3141592654e1 * sin(0.3141592654e1 * p[2]), 0.16e2 / 0.17e2 * sin(0.3141592654e1 * p[1]) * cos(0.3141592654e1 * p[2]) * 0.3141592654e1);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_xy_RhsVectorFun2 (const DROPS::Point3DCL& p, double)
{
    double theta = 1;
    DROPS::Point3DCL v (-0.256e3 / 0.289e3 * theta * pow(p[1], 0.2e1) + 0.1568e4 / 0.289e3 * theta - 0.8e1 / 0.17e2 * theta * pow(p[2], 0.2e1) - 0.16e2 / 0.17e2 * theta * p[2] * p[1] + 0.4e1 / 0.17e2 * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - 0.16e2 / 0.17e2 * pow(p[1], 0.2e1) - 0.16e2 / 0.17e2 * pow(p[2], 0.2e1) + 0.64e2 / 0.17e2
    , -0.64e2 / 0.17e2 * theta * p[2] * p[1] - 0.32e2 / 0.17e2 * theta * pow(p[1], 0.2e1) + 0.400e3 / 0.17e2 * theta - 0.4e1 * theta * pow(p[2], 0.2e1) + pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - 0.4e1 * pow(p[1], 0.2e1) - 0.4e1 * pow(p[2], 0.2e1) + 0.16e2
    , -0.1024e4 / 0.289e3 * theta * pow(p[1], 0.2e1) + 0.6272e4 / 0.289e3 * theta - 0.32e2 / 0.17e2 * theta * pow(p[2], 0.2e1) - 0.64e2 / 0.17e2 * theta * p[2] * p[1] + 0.16e2 / 0.17e2 * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - 0.64e2 / 0.17e2 * pow(p[1], 0.2e1) - 0.64e2 / 0.17e2 * pow(p[2], 0.2e1) + 0.256e3 / 0.17e2
    );
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_xy_vSolVectorFun2 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0.4e1 / 0.17e2 * (pow(p[2], 0.2e1) - 0.4e1) * (pow(p[1], 0.2e1) - 0.4e1), (pow(p[2], 0.2e1) - 0.4e1) * (pow(p[1], 0.2e1) - 0.4e1), 0.16e2 / 0.17e2 * (pow(p[2], 0.2e1) - 0.4e1) * (pow(p[1], 0.2e1) - 0.4e1));
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_xy_vSolVectorFun2_Gradient1 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0, 0.8e1 / 0.17e2 * (pow(p[2], 0.2e1) - 0.4e1) * p[1], 0.8e1 / 0.17e2 * p[2] * (pow(p[1], 0.2e1) - 0.4e1));
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_xy_vSolVectorFun2_Gradient2 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0, 0.2e1 * (pow(p[2], 0.2e1) - 0.4e1) * p[1], 0.2e1 * p[2] * (pow(p[1], 0.2e1) - 0.4e1));
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_xy_vSolVectorFun2_Gradient3 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0, 0.32e2 / 0.17e2 * (pow(p[2], 0.2e1) - 0.4e1) * p[1], 0.32e2 / 0.17e2 * p[2] * (pow(p[1], 0.2e1) - 0.4e1));
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_xy_RhsVectorFun3 (const DROPS::Point3DCL& p, double)
{
    double theta = 1;
    DROPS::Point3DCL v (-0.1536e4 / 0.289e3 * theta * pow(p[1], 0.4e1) * pow(p[2], 0.2e1) + 0.18816e5 / 0.289e3 * theta * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - 0.33280e5 / 0.289e3 * theta * pow(p[2], 0.2e1) + 0.2048e4 / 0.289e3 * theta * pow(p[1], 0.4e1) - 0.29440e5 / 0.289e3 * theta * pow(p[1], 0.2e1) + 0.50176e5 / 0.289e3 * theta - 0.48e2 / 0.17e2 * theta * pow(p[1], 0.2e1) * pow(p[2], 0.4e1) + 0.64e2 / 0.17e2 * theta * pow(p[2], 0.4e1) - 0.64e2 / 0.17e2 * theta * pow(p[1], 0.3e1) * pow(p[2], 0.3e1) + 0.256e3 / 0.17e2 * theta * pow(p[1], 0.3e1) * p[2] + 0.256e3 / 0.17e2 * theta * p[1] * pow(p[2], 0.3e1) - 0.1024e4 / 0.17e2 * theta * p[2] * p[1] + 0.4e1 / 0.17e2 * pow(p[1], 0.4e1) * pow(p[2], 0.4e1) - 0.32e2 / 0.17e2 * pow(p[1], 0.4e1) * pow(p[2], 0.2e1) - 0.32e2 / 0.17e2 * pow(p[1], 0.2e1) * pow(p[2], 0.4e1) + 0.64e2 / 0.17e2 * pow(p[1], 0.4e1) + 0.256e3 / 0.17e2 * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) + 0.64e2 / 0.17e2 * pow(p[2], 0.4e1) - 0.512e3 / 0.17e2 * pow(p[1], 0.2e1) - 0.512e3 / 0.17e2 * pow(p[2], 0.2e1) + 0.1024e4 / 0.17e2
    , -0.256e3 / 0.17e2 * theta * pow(p[1], 0.3e1) * pow(p[2], 0.3e1) + 0.1024e4 / 0.17e2 * theta * pow(p[1], 0.3e1) * p[2] + 0.1024e4 / 0.17e2 * theta * p[1] * pow(p[2], 0.3e1) - 0.4096e4 / 0.17e2 * theta * p[2] * p[1] - 0.192e3 / 0.17e2 * theta * pow(p[1], 0.4e1) * pow(p[2], 0.2e1) + 0.256e3 / 0.17e2 * theta * pow(p[1], 0.4e1) + 0.4800e4 / 0.17e2 * theta * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - 0.8576e4 / 0.17e2 * theta * pow(p[1], 0.2e1) - 0.7424e4 / 0.17e2 * theta * pow(p[2], 0.2e1) + 0.12800e5 / 0.17e2 * theta - 0.24e2 * theta * pow(p[1], 0.2e1) * pow(p[2], 0.4e1) + 0.32e2 * theta * pow(p[2], 0.4e1) + pow(p[1], 0.4e1) * pow(p[2], 0.4e1) - 0.8e1 * pow(p[1], 0.4e1) * pow(p[2], 0.2e1) - 0.8e1 * pow(p[1], 0.2e1) * pow(p[2], 0.4e1) + 0.16e2 * pow(p[1], 0.4e1) + 0.64e2 * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) + 0.16e2 * pow(p[2], 0.4e1) - 0.128e3 * pow(p[1], 0.2e1) - 0.128e3 * pow(p[2], 0.2e1) + 0.256e3
    , -0.6144e4 / 0.289e3 * theta * pow(p[1], 0.4e1) * pow(p[2], 0.2e1) + 0.75264e5 / 0.289e3 * theta * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - 0.133120e6 / 0.289e3 * theta * pow(p[2], 0.2e1) + 0.8192e4 / 0.289e3 * theta * pow(p[1], 0.4e1) - 0.117760e6 / 0.289e3 * theta * pow(p[1], 0.2e1) + 0.200704e6 / 0.289e3 * theta - 0.192e3 / 0.17e2 * theta * pow(p[1], 0.2e1) * pow(p[2], 0.4e1) + 0.256e3 / 0.17e2 * theta * pow(p[2], 0.4e1) - 0.256e3 / 0.17e2 * theta * pow(p[1], 0.3e1) * pow(p[2], 0.3e1) + 0.1024e4 / 0.17e2 * theta * pow(p[1], 0.3e1) * p[2] + 0.1024e4 / 0.17e2 * theta * p[1] * pow(p[2], 0.3e1) - 0.4096e4 / 0.17e2 * theta * p[2] * p[1] + 0.16e2 / 0.17e2 * pow(p[1], 0.4e1) * pow(p[2], 0.4e1) - 0.128e3 / 0.17e2 * pow(p[1], 0.4e1) * pow(p[2], 0.2e1) - 0.128e3 / 0.17e2 * pow(p[1], 0.2e1) * pow(p[2], 0.4e1) + 0.256e3 / 0.17e2 * pow(p[1], 0.4e1) + 0.1024e4 / 0.17e2 * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) + 0.256e3 / 0.17e2 * pow(p[2], 0.4e1) - 0.2048e4 / 0.17e2 * pow(p[1], 0.2e1) - 0.2048e4 / 0.17e2 * pow(p[2], 0.2e1) + 0.4096e4 / 0.17e2
    );
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_xy_vSolVectorFun3 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0.4e1 / 0.17e2 * pow(pow(p[1], 0.2e1) - 0.4e1, 0.2e1) * pow(pow(p[2], 0.2e1) - 0.4e1, 0.2e1), pow(pow(p[1], 0.2e1) - 0.4e1, 0.2e1) * pow(pow(p[2], 0.2e1) - 0.4e1, 0.2e1), 0.16e2 / 0.17e2 * pow(pow(p[1], 0.2e1) - 0.4e1, 0.2e1) * pow(pow(p[2], 0.2e1) - 0.4e1, 0.2e1));
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_xy_vSolVectorFun3_Gradient1 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0, 0.16e2 / 0.17e2 * (pow(p[1], 0.2e1) - 0.4e1) * pow(pow(p[2], 0.2e1) - 0.4e1, 0.2e1) * p[1], 0.16e2 / 0.17e2 * pow(pow(p[1], 0.2e1) - 0.4e1, 0.2e1) * (pow(p[2], 0.2e1) - 0.4e1) * p[2]);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_xy_vSolVectorFun3_Gradient2 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0, 0.4e1 * (pow(p[1], 0.2e1) - 0.4e1) * pow(pow(p[2], 0.2e1) - 0.4e1, 0.2e1) * p[1], 0.4e1 * pow(pow(p[1], 0.2e1) - 0.4e1, 0.2e1) * (pow(p[2], 0.2e1) - 0.4e1) * p[2]);
    return v;
}

DROPS::Point3DCL Test_A_plus_M_tilted_plane_xy_vSolVectorFun3_Gradient3 (const DROPS::Point3DCL& p, double)
{
    //const double PI = 3.141592653589793;
    DROPS::Point3DCL v (0, 0.64e2 / 0.17e2 * (pow(p[1], 0.2e1) - 0.4e1) * pow(pow(p[2], 0.2e1) - 0.4e1, 0.2e1) * p[1], 0.64e2 / 0.17e2 * pow(pow(p[1], 0.2e1) - 0.4e1, 0.2e1) * (pow(p[2], 0.2e1) - 0.4e1) * p[2]);
    return v;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////Functions for TestP2Matrices////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DROPS::Point3DCL TestP2Matrices_A_Vector_Fun_1 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v ((-p[1] - p[2]) * p[0] + p[1] * p[1] + p[2] * p[2], (-p[0] - p[2]) * p[1] + p[0] * p[0] + p[2] * p[2], (-p[0] - p[1]) * p[2] + p[0] * p[0] + p[1] * p[1]);
    return v;
}

DROPS::Point3DCL TestP2Matrices_A_Rhs_1 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v (0.2e1 * pow(p[0], 0.4e1) + (0.5e1 * p[1] + 0.5e1 * p[2]) * pow(p[0], 0.3e1) + (-p[1] * p[1] - p[2] * p[2]) * p[0] * p[0] + 0.5e1 * (p[1] * p[1] + p[2] * p[2] - 0.7e1 / 0.5e1) * (p[1] + p[2]) * p[0] - 0.3e1 * (p[1] * p[1] + p[2] * p[2] - 0.1e1 / 0.3e1) * (p[1] * p[1] + p[2] * p[2] - 0.2e1)
    , 0.2e1 * pow(p[1], 0.4e1) + (0.5e1 * p[0] + 0.5e1 * p[2]) * pow(p[1], 0.3e1) + (-p[0] * p[0] - p[2] * p[2]) * p[1] * p[1] + 0.5e1 * (p[0] * p[0] + p[2] * p[2] - 0.7e1 / 0.5e1) * (p[0] + p[2]) * p[1] - 0.3e1 * (p[0] * p[0] + p[2] * p[2] - 0.1e1 / 0.3e1) * (p[0] * p[0] + p[2] * p[2] - 0.2e1)
    , 0.2e1 * pow(p[2], 0.4e1) + (0.5e1 * p[0] + 0.5e1 * p[1]) * pow(p[2], 0.3e1) + (-p[0] * p[0] - p[1] * p[1]) * p[2] * p[2] + 0.5e1 * (p[0] * p[0] + p[1] * p[1] - 0.7e1 / 0.5e1) * (p[0] + p[1]) * p[2] - 0.3e1 * (p[0] * p[0] + p[1] * p[1] - 0.2e1) * (p[0] * p[0] + p[1] * p[1] - 0.1e1 / 0.3e1)
    );
    return v;
}

#endif
