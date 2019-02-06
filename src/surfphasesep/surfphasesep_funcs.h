/// \file surfnavierstokes_funcs.h
/// \brief Functions for surfnavierstokes.cpp
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

using namespace DROPS;

DROPS::ParamCL P;

namespace ParameterNS {
	double h=0.1;
	double nu;
	double sigma;
	double eps;
}

namespace Torus {
    double R=1;
    double r=0.2;
       double k=0.0;
       double f=0.0;
//    double k=0.3;
//    double f=4.0;
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
    return p[0]*p[0]+p[1]*p[1]+p[2]*p[2]-1;
}

double tamarind (const DROPS::Point3DCL& p, double)
{
//    DROPS::Point3DCL x( p - P.get<DROPS::Point3DCL>("Exp.PosDrop"));
    double PI = 3.141592653589793;
//    return x.norm() - P.get<DROPS::Point3DCL>("Exp.RadDrop")[0];
    return p[0]*p[0]/4.+p[1]*p[1]+4*p[2]*p[2]*pow(1+0.5*sin(PI*p[0]),-2.0)-1;

}

double spindle(const DROPS::Point3DCL& p, double)
{
//    DROPS::Point3DCL x( p - P.get<DROPS::Point3DCL>("Exp.PosDrop"));
    double PI = 3.141592653589793;
//    return x.norm() - P.get<DROPS::Point3DCL>("Exp.RadDrop")[0];
        if (std::abs(p[0])<5.) {
            return (pow(cos((p[0]/5.)*(PI/2.0)),2.0)
                    //(1.-pow(p[0],2.0)
     -pow(p[1]/0.25,2.0)-pow(p[2]/0.25,2.0));
        } else return -pow(5.,2.0);

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
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
			return DROPS::Point3DCL(0,0,0);
		}
	DROPS::Point3DCL v;
	v[0] = 0.2e1 * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * (pow(p[1], 0.4e1) + (pow(p[0], 0.2e1) - 0.3e1 * pow(p[2], 0.2e1) - 0.5e1 * p[0] + 0.5e1 / 0.2e1 * p[2]) * pow(p[1], 0.2e1) + 0.7e1 * pow(p[0], 0.2e1) * pow(p[2], 0.2e1) - 0.4e1 * pow(p[2], 0.4e1) - 0.5e1 / 0.2e1 * pow(p[0], 0.2e1) * p[2] + 0.5e1 / 0.2e1 * pow(p[2], 0.3e1)) + (-pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - pow(p[2], 0.4e1) - pow(p[0], 0.2e1) * p[2] - p[0] * pow(p[1], 0.2e1)) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));
	v[1] = -0.2e1 * p[1] * ((-0.10e2 * p[0] - 0.5e1) * pow(p[2], 0.2e1) + 0.5e1 * p[0] * p[2] + p[0] * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - 0.5e1 * p[0])) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) + (pow(p[0], 0.2e1) + (pow(p[2], 0.2e1) - p[2]) * p[0] + pow(p[2], 0.2e1)) * p[1] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));
	v[2] = (0.8e1 * p[0] * pow(p[2], 0.3e1) - 0.5e1 * p[0] * pow(p[2], 0.2e1) + ((-0.14e2 * p[0] - 0.10e2) * pow(p[1], 0.2e1) - 0.14e2 * pow(p[0], 0.3e1)) * p[2] + 0.5e1 * pow(p[0], 0.3e1) + 0.5e1 * p[0] * pow(p[1], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) + (pow(p[1], 0.2e1) * (p[0] - p[2]) + p[0] * pow(p[2], 0.3e1) + pow(p[0], 0.3e1)) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));
    return v;
}

double Test_A_plus_M_rhs2Fun1 (const DROPS::Point3DCL& p, double)
{

	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
		return(0.0);
	}

    double v;
    v = (pow(p[0], 0.2e1) + (0.4e1 * pow(p[2], 0.2e1) - 0.3e1 * p[2]) * p[0] - 0.2e1 * pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));

    return  (-1)*v;
}



double Test_A_plus_M_pSolScalarFun1 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
			return(0.0);
		}
    return 0;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun1 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
			return DROPS::Point3DCL(0,0,0);
		}
    DROPS::Point3DCL v;
    v[0] = (-pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - pow(p[2], 0.4e1) - pow(p[0], 0.2e1) * p[2] - p[0] * pow(p[1], 0.2e1)) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));
    v[1] = (pow(p[0], 0.2e1) + (pow(p[2], 0.2e1) - p[2]) * p[0] + pow(p[2], 0.2e1)) * p[1] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));
    v[2] = (pow(p[0], 0.3e1) + (pow(p[2], 0.3e1) + pow(p[1], 0.2e1)) * p[0] - pow(p[1], 0.2e1) * p[2]) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));

    return v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun1_Gradient1 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
			return DROPS::Point3DCL(0,0,0);
		}
    DROPS::Point3DCL v;
    v[0] = (pow(p[0], 0.2e1) * pow(p[1], 0.2e1) + 0.2e1 * p[2] * (p[2] - 0.1e1) * (pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[0] - pow(p[1], 0.4e1) - pow(p[1], 0.2e1) * pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    v[1] = -0.2e1 * (pow(p[0], 0.2e1) + (pow(p[2], 0.2e1) - p[2]) * p[0] + pow(p[2], 0.2e1)) * p[1] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * p[0];
    v[2] = (-0.2e1 * pow(p[2], 0.5e1) + (-0.4e1 * pow(p[0], 0.2e1) - 0.4e1 * pow(p[1], 0.2e1)) * pow(p[2], 0.3e1) + pow(p[0], 0.2e1) * pow(p[2], 0.2e1) - 0.2e1 * pow(p[1], 0.2e1) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - p[0]) * p[2] - pow(p[0], 0.4e1) - pow(p[0], 0.2e1) * pow(p[1], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);

    return v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun1_Gradient2 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
			return DROPS::Point3DCL(0,0,0);
		}

	DROPS::Point3DCL v;
    v[0] = -p[1] * (-pow(p[2], 0.4e1) + pow(p[2], 0.3e1) + (pow(p[0], 0.2e1) - pow(p[1], 0.2e1)) * pow(p[2], 0.2e1) + (-pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * p[2] - 0.2e1 * p[0] * pow(p[1], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    v[1] = (pow(p[0], 0.2e1) + (pow(p[2], 0.2e1) - p[2]) * p[0] + pow(p[2], 0.2e1)) * (pow(p[0], 0.2e1) - pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    v[2] = 0.2e1 * p[1] * ((p[2] - 0.1e1 / 0.2e1) * pow(p[0], 0.3e1) + ((p[2] - 0.1e1 / 0.2e1) * pow(p[1], 0.2e1) + pow(p[2], 0.2e1) / 0.2e1) * p[0] + pow(p[1], 0.2e1) * p[2]) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);

    return v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun1_Gradient3 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
			return DROPS::Point3DCL(0,0,0);
		}
	DROPS::Point3DCL v;
    v[0] = (pow(p[2], 0.5e1) + (-pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * pow(p[2], 0.3e1) + (0.3e1 * pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * pow(p[2], 0.2e1) + 0.2e1 * p[0] * pow(p[1], 0.2e1) * p[2] + pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1), 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    v[1] = -0.2e1 * ((p[0] + 0.1e1) * pow(p[2], 0.2e1) - p[0] * p[2] + pow(p[0], 0.2e1)) * p[1] * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    v[2] = ((0.3e1 * pow(p[2], 0.2e1) - 0.2e1 * p[2]) * pow(p[0], 0.3e1) - pow(p[0], 0.2e1) * pow(p[1], 0.2e1) + (0.3e1 * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) + pow(p[2], 0.4e1) - 0.2e1 * pow(p[1], 0.2e1) * p[2]) * p[0] - pow(p[1], 0.4e1) + pow(p[1], 0.2e1) * pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);

    return v;
}





DROPS::Point3DCL Test_A_plus_M_RhsVectorFun2 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(-p[1], p[0], 0);
    return v;
}

double Test_A_plus_M_pSolScalarFun2 (const DROPS::Point3DCL&, double)
{
    return 0. ;
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


DROPS::Point3DCL Test_A_plus_M_vSolVectorFun6 (const DROPS::Point3DCL& p, double)
{
    if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
    	return DROPS::Point3DCL(0,0,0);
    }
    DROPS::Point3DCL v;
    v[0] = 0.1e1 - 0.1e1 *  pow(p[0], 0.2e1)
    				     /( pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1) );
    v[1] =       - 0.1e1 *      p[0] * p[1]
						 /( pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1) );
    v[2] =       - 0.1e1 *      p[0] * p[2]
					     /( pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1) );
    return v;
}

double Test_A_plus_M_pSolScalarFun6 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
	    	return 0;
	    }
    return pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) * p[1];;
}



DROPS::Point3DCL Test_A_plus_M_RhsVectorFun6 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
		return DROPS::Point3DCL(0,0,0);
	}

	DROPS::Point3DCL v;
	v[0] = (0.2e1 * pow(p[1], 0.2e1) + 0.2e1 * pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) + 0.1e1 - 0.1e1 / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(p[0], 0.2e1) - pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1) * p[1] * p[0];
	v[1] = -0.2e1 * p[0] * p[1] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) - 0.1e1 / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[0] * p[1] + (pow(p[0], 0.2e1) + pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1);
	v[2] = -0.2e1 * p[0] * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) - 0.1e1 / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[0] * p[2] - pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1) * p[1] * p[2];
    return  v ;
}

double Test_A_plus_M_rhs2Fun6 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
		return(0.0);
	}
    double v;
    v = -0.2e1 / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[0];

    return  (-1)*v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun6_Gradient1 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
		return DROPS::Point3DCL(0,0,0);
	}

	DROPS::Point3DCL v;
	v[0] = -0.2e1 * p[0] * (pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
	v[1] = 0.2e1 * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * pow(p[0], 0.2e1) * p[1];
	v[2] = 0.2e1 * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * pow(p[0], 0.2e1) * p[2];
	return v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun6_Gradient2 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
		return DROPS::Point3DCL(0,0,0);
	}

    DROPS::Point3DCL v;
    v[0] = p[1] * (pow(p[0], 0.2e1) - pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    v[1] = -p[0] * (pow(p[0], 0.2e1) - pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    v[2] = 0.2e1 * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * p[0] * p[1] * p[2];
    return v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun6_Gradient3 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
		return DROPS::Point3DCL(0,0,0);
	}

	DROPS::Point3DCL v;
	v[0] = p[2] * (pow(p[0], 0.2e1) - pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
	v[1] = 0.2e1 * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * p[0] * p[1] * p[2];
	v[2] = -p[0] * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
	return v;
}



DROPS::Point3DCL Test_A_plus_M_vSolVectorFun7 (const DROPS::Point3DCL& p, double)
{
    return DROPS::Point3DCL(0,0,0);

}

double Test_A_plus_M_pSolScalarFun7 (const DROPS::Point3DCL& p, double)
{
    return p[0] * pow(p[1], 0.3e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) + pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) * p[2];
}

DROPS::Point3DCL Test_A_plus_M_RhsVectorFun7 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
		return DROPS::Point3DCL(0,0,0);
	}


    DROPS::Point3DCL v;
    v[0] = (-p[0] * pow(p[2], 0.5e1) + (-0.2e1 * pow(p[0], 0.3e1) - 0.2e1 * p[0] * pow(p[1], 0.2e1)) * pow(p[2], 0.3e1) + sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(p[1], 0.3e1) * pow(p[2], 0.2e1) - p[0] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1), 0.2e1) * p[2] - sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(p[1], 0.3e1) * (0.3e1 * pow(p[0], 0.2e1) - pow(p[1], 0.2e1))) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.7e1 / 0.2e1);
    v[1] = p[1] * (-pow(p[2], 0.5e1) + (-0.2e1 * pow(p[0], 0.2e1) - 0.2e1 * pow(p[1], 0.2e1)) * pow(p[2], 0.3e1) + 0.3e1 * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[0] * p[1] * pow(p[2], 0.2e1) - pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1), 0.2e1) * p[2] + sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[0] * p[1] * (0.3e1 * pow(p[0], 0.2e1) - pow(p[1], 0.2e1))) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.7e1 / 0.2e1);
    v[2] = (pow(p[0], 0.6e1) + (0.3e1 * pow(p[1], 0.2e1) + 0.2e1 * pow(p[2], 0.2e1)) * pow(p[0], 0.4e1) + (0.3e1 * pow(p[1], 0.4e1) + 0.4e1 * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) + pow(p[2], 0.4e1)) * pow(p[0], 0.2e1) - 0.4e1 * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[0] * pow(p[1], 0.3e1) * p[2] + pow(p[1], 0.2e1) * pow(pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.7e1 / 0.2e1);
    return v;
}

double Test_A_plus_M_rhs2Fun7 (const DROPS::Point3DCL& p, double)
{
	return(0.0);
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun7_Gradient1 (const DROPS::Point3DCL& p, double)
{
	return DROPS::Point3DCL(0,0,0);

}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun7_Gradient2 (const DROPS::Point3DCL& p, double)
{
	return DROPS::Point3DCL(0,0,0);

}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun7_Gradient3 (const DROPS::Point3DCL& p, double)
{
	return DROPS::Point3DCL(0,0,0);
}



DROPS::Point3DCL Test_A_plus_M_RhsVectorFun8 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
			return DROPS::Point3DCL(0,0,0);
		}
	DROPS::Point3DCL v,w;
	v[0] = 0.2e1 * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * (pow(p[1], 0.4e1) + (pow(p[0], 0.2e1) - 0.3e1 * pow(p[2], 0.2e1) - 0.5e1 * p[0] + 0.5e1 / 0.2e1 * p[2]) * pow(p[1], 0.2e1) + 0.7e1 * pow(p[0], 0.2e1) * pow(p[2], 0.2e1) - 0.4e1 * pow(p[2], 0.4e1) - 0.5e1 / 0.2e1 * pow(p[0], 0.2e1) * p[2] + 0.5e1 / 0.2e1 * pow(p[2], 0.3e1)) + (-pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - pow(p[2], 0.4e1) - pow(p[0], 0.2e1) * p[2] - p[0] * pow(p[1], 0.2e1)) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));
	v[1] = -0.2e1 * p[1] * ((-0.10e2 * p[0] - 0.5e1) * pow(p[2], 0.2e1) + 0.5e1 * p[0] * p[2] + p[0] * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - 0.5e1 * p[0])) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) + (pow(p[0], 0.2e1) + (pow(p[2], 0.2e1) - p[2]) * p[0] + pow(p[2], 0.2e1)) * p[1] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));
	v[2] = (0.8e1 * p[0] * pow(p[2], 0.3e1) - 0.5e1 * p[0] * pow(p[2], 0.2e1) + ((-0.14e2 * p[0] - 0.10e2) * pow(p[1], 0.2e1) - 0.14e2 * pow(p[0], 0.3e1)) * p[2] + 0.5e1 * pow(p[0], 0.3e1) + 0.5e1 * p[0] * pow(p[1], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) + (pow(p[1], 0.2e1) * (p[0] - p[2]) + p[0] * pow(p[2], 0.3e1) + pow(p[0], 0.3e1)) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));

	    w[0] = (-p[0] * pow(p[2], 0.5e1) + (-0.2e1 * pow(p[0], 0.3e1) - 0.2e1 * p[0] * pow(p[1], 0.2e1)) * pow(p[2], 0.3e1) + sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(p[1], 0.3e1) * pow(p[2], 0.2e1) - p[0] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1), 0.2e1) * p[2] - sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(p[1], 0.3e1) * (0.3e1 * pow(p[0], 0.2e1) - pow(p[1], 0.2e1))) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.7e1 / 0.2e1);
	    w[1] = p[1] * (-pow(p[2], 0.5e1) + (-0.2e1 * pow(p[0], 0.2e1) - 0.2e1 * pow(p[1], 0.2e1)) * pow(p[2], 0.3e1) + 0.3e1 * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[0] * p[1] * pow(p[2], 0.2e1) - pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1), 0.2e1) * p[2] + sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[0] * p[1] * (0.3e1 * pow(p[0], 0.2e1) - pow(p[1], 0.2e1))) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.7e1 / 0.2e1);
	    w[2] = (pow(p[0], 0.6e1) + (0.3e1 * pow(p[1], 0.2e1) + 0.2e1 * pow(p[2], 0.2e1)) * pow(p[0], 0.4e1) + (0.3e1 * pow(p[1], 0.4e1) + 0.4e1 * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) + pow(p[2], 0.4e1)) * pow(p[0], 0.2e1) - 0.4e1 * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[0] * pow(p[1], 0.3e1) * p[2] + pow(p[1], 0.2e1) * pow(pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.7e1 / 0.2e1);

	return v+w;
}

double Test_A_plus_M_rhs2Fun8 (const DROPS::Point3DCL& p, double)
{

	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
		return(0.0);
	}

    double v;
    v = (pow(p[0], 0.2e1) + (0.4e1 * pow(p[2], 0.2e1) - 0.3e1 * p[2]) * p[0] - 0.2e1 * pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));

    return  (-1)*v;
}

double Test_A_plus_M_pSolScalarFun8 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
			return(0.0);
		}
	return (p[0] * pow(p[1], 0.3e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) + pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) * p[2]);
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun8 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
			return DROPS::Point3DCL(0,0,0);
		}
    DROPS::Point3DCL v;
    v[0] = (-pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - pow(p[2], 0.4e1) - pow(p[0], 0.2e1) * p[2] - p[0] * pow(p[1], 0.2e1)) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));
    v[1] = (pow(p[0], 0.2e1) + (pow(p[2], 0.2e1) - p[2]) * p[0] + pow(p[2], 0.2e1)) * p[1] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));
    v[2] = (pow(p[0], 0.3e1) + (pow(p[2], 0.3e1) + pow(p[1], 0.2e1)) * p[0] - pow(p[1], 0.2e1) * p[2]) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));

    return v;
}


DROPS::Point3DCL Test_A_plus_M_vSolVectorFun8_Gradient1 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
			return DROPS::Point3DCL(0,0,0);
		}
    DROPS::Point3DCL v;
    v[0] = (pow(p[0], 0.2e1) * pow(p[1], 0.2e1) + 0.2e1 * p[2] * (p[2] - 0.1e1) * (pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[0] - pow(p[1], 0.4e1) - pow(p[1], 0.2e1) * pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    v[1] = -0.2e1 * (pow(p[0], 0.2e1) + (pow(p[2], 0.2e1) - p[2]) * p[0] + pow(p[2], 0.2e1)) * p[1] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * p[0];
    v[2] = (-0.2e1 * pow(p[2], 0.5e1) + (-0.4e1 * pow(p[0], 0.2e1) - 0.4e1 * pow(p[1], 0.2e1)) * pow(p[2], 0.3e1) + pow(p[0], 0.2e1) * pow(p[2], 0.2e1) - 0.2e1 * pow(p[1], 0.2e1) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - p[0]) * p[2] - pow(p[0], 0.4e1) - pow(p[0], 0.2e1) * pow(p[1], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);

    return v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun8_Gradient2 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
			return DROPS::Point3DCL(0,0,0);
		}

	DROPS::Point3DCL v;
    v[0] = -p[1] * (-pow(p[2], 0.4e1) + pow(p[2], 0.3e1) + (pow(p[0], 0.2e1) - pow(p[1], 0.2e1)) * pow(p[2], 0.2e1) + (-pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * p[2] - 0.2e1 * p[0] * pow(p[1], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    v[1] = (pow(p[0], 0.2e1) + (pow(p[2], 0.2e1) - p[2]) * p[0] + pow(p[2], 0.2e1)) * (pow(p[0], 0.2e1) - pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    v[2] = 0.2e1 * p[1] * ((p[2] - 0.1e1 / 0.2e1) * pow(p[0], 0.3e1) + ((p[2] - 0.1e1 / 0.2e1) * pow(p[1], 0.2e1) + pow(p[2], 0.2e1) / 0.2e1) * p[0] + pow(p[1], 0.2e1) * p[2]) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);

    return v;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun8_Gradient3 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
			return DROPS::Point3DCL(0,0,0);
		}
	DROPS::Point3DCL v;
    v[0] = (pow(p[2], 0.5e1) + (-pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * pow(p[2], 0.3e1) + (0.3e1 * pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * pow(p[2], 0.2e1) + 0.2e1 * p[0] * pow(p[1], 0.2e1) * p[2] + pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1), 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    v[1] = -0.2e1 * ((p[0] + 0.1e1) * pow(p[2], 0.2e1) - p[0] * p[2] + pow(p[0], 0.2e1)) * p[1] * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    v[2] = ((0.3e1 * pow(p[2], 0.2e1) - 0.2e1 * p[2]) * pow(p[0], 0.3e1) - pow(p[0], 0.2e1) * pow(p[1], 0.2e1) + (0.3e1 * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) + pow(p[2], 0.4e1) - 0.2e1 * pow(p[1], 0.2e1) * p[2]) * p[0] - pow(p[1], 0.4e1) + pow(p[1], 0.2e1) * pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);

    return v;
}



DROPS::Point3DCL Test_A_plus_M_vSolVectorFun9 (const DROPS::Point3DCL& p, double t )
{
    DROPS::Point3DCL gen_c;
    double a=-1.;
    gen_c[0] = (0.1e1 + a * exp(-(double) (6 * t))) * (pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1);
    gen_c[1] = -(0.1e1 + a * exp(-(double) (6 * t))) * p[0] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) * p[1];
    gen_c[2] = (0.1e1 + a * exp(-(double) (6 * t))) * p[0] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) * p[2];
    return gen_c;

}


DROPS::Point3DCL Test_A_plus_M_vSolVectorFun95 (const DROPS::Point3DCL& p, double t )
{
    DROPS::Point3DCL gen_c;
    double a=-1;
    gen_c[0] = (0.1e1 + a * exp(-(double) (6 * t))) * (0.1e1 - pow(p[0], 0.2e1) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)));
    gen_c[1] = -(0.1e1 + a * exp(-(double) (6 * t))) * p[0] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[1];
    gen_c[2] = -(0.1e1 + a * exp(-(double) (6 * t))) * p[0] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[2];
    return gen_c;

}


DROPS::Point3DCL Test_A_plus_M_RhsLaplacianVectorFun9 (const DROPS::Point3DCL& p, double t )
{
    if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
        return DROPS::Point3DCL(0,0,0);
    }
    double a=-1;
    DROPS::Point3DCL gen_c;
    gen_c[0] = 0.4e1 * (0.1e1 + a * exp(-(double) (6 * t))) * (p[1] - p[2]) * (p[1] + p[2]) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1);
    gen_c[1] = -0.4e1 * (0.1e1 + a * exp(-(double) (6 * t))) * p[1] * p[0] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1);
    gen_c[2] = 0.4e1 * (0.1e1 + a * exp(-(double) (6 * t))) * p[2] * p[0] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1);
    return gen_c;
}

DROPS::Point3DCL Test_A_plus_M_RhsPressureVectorFun9 (const DROPS::Point3DCL& p, double t )
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
			return DROPS::Point3DCL(0,0,0);
		}
	DROPS::Point3DCL gen_c;
    gen_c[0] = (-p[0] * pow(p[2], 0.5e1) + (-0.2e1 * pow(p[0], 0.3e1) - 0.2e1 * p[0] * pow(p[1], 0.2e1)) * pow(p[2], 0.3e1) + sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(p[1], 0.3e1) * pow(p[2], 0.2e1) - p[0] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1), 0.2e1) * p[2] - sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(p[1], 0.3e1) * (0.3e1 * pow(p[0], 0.2e1) - pow(p[1], 0.2e1))) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.7e1 / 0.2e1);
    gen_c[1] = p[1] * (-pow(p[2], 0.5e1) + (-0.2e1 * pow(p[0], 0.2e1) - 0.2e1 * pow(p[1], 0.2e1)) * pow(p[2], 0.3e1) + 0.3e1 * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[0] * p[1] * pow(p[2], 0.2e1) - pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1), 0.2e1) * p[2] + sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[0] * p[1] * (0.3e1 * pow(p[0], 0.2e1) - pow(p[1], 0.2e1))) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.7e1 / 0.2e1);
    gen_c[2] = (pow(p[0], 0.6e1) + (0.3e1 * pow(p[1], 0.2e1) + 0.2e1 * pow(p[2], 0.2e1)) * pow(p[0], 0.4e1) + (0.3e1 * pow(p[1], 0.4e1) + 0.4e1 * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) + pow(p[2], 0.4e1)) * pow(p[0], 0.2e1) - 0.4e1 * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[0] * pow(p[1], 0.3e1) * p[2] + pow(p[1], 0.2e1) * pow(pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.7e1 / 0.2e1);
	return gen_c;
}

DROPS::Point3DCL Test_A_plus_M_rhs_inertial_Fun9(const DROPS::Point3DCL& p, double t )
{
    if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
        return DROPS::Point3DCL(0,0,0);
    }
    double a=-1;
    DROPS::Point3DCL gen_c;
    gen_c[0] = -0.6e1 * a * exp(-(double) (6 * t)) * (pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1);
    gen_c[1] = 0.6e1 * a * exp(-(double) (6 * t)) * p[0] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) * p[1];
    gen_c[2] = -0.6e1 * a * exp(-(double) (6 * t)) * p[0] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) * p[2];
    return  gen_c;
}

DROPS::Point3DCL Test_A_plus_M_rhs_inertial_Fun95(const DROPS::Point3DCL& p, double t )
{
    if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
        return DROPS::Point3DCL(0,0,0);
    }
    double a=-1;
    DROPS::Point3DCL gen_c;
    gen_c[0] = -0.6e1 * a * exp(-(double) (6 * t)) * (0.1e1 - pow(p[0], 0.2e1) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)));
    gen_c[1] = 0.6e1 * a * exp(-(double) (6 * t)) * p[0] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[1];
    gen_c[2] = 0.6e1 * a * exp(-(double) (6 * t)) * p[0] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[2];
    return  gen_c;
}

DROPS::Point3DCL Test_A_plus_M_rhs_convFun9(const DROPS::Point3DCL& p, double t )
{
    if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
        return DROPS::Point3DCL(0,0,0);
    }
    double a=-1;
    DROPS::Point3DCL gen_c;
    gen_c[0] = -pow(0.1e1 + a * exp(-(double) (6 * t)), 0.2e1) * p[0] * (pow(p[0], 0.2e1) * pow(p[1], 0.2e1) + pow(p[0], 0.2e1) * pow(p[2], 0.2e1) + pow(p[1], 0.4e1) + 0.6e1 * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) + pow(p[2], 0.4e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    gen_c[1] = pow(0.1e1 + a * exp(-(double) (6 * t)), 0.2e1) * p[1] * (pow(p[0], 0.4e1) + pow(p[0], 0.2e1) * pow(p[1], 0.2e1) + 0.3e1 * pow(p[0], 0.2e1) * pow(p[2], 0.2e1) - 0.2e1 * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) + 0.2e1 * pow(p[2], 0.4e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    gen_c[2] = pow(0.1e1 + a * exp(-(double) (6 * t)), 0.2e1) * p[2] * (pow(p[0], 0.4e1) + 0.3e1 * pow(p[0], 0.2e1) * pow(p[1], 0.2e1) + pow(p[0], 0.2e1) * pow(p[2], 0.2e1) + 0.2e1 * pow(p[1], 0.4e1) - 0.2e1 * pow(p[1], 0.2e1) * pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    return  gen_c;
}


DROPS::Point3DCL Test_A_plus_M_rhs_convFun95(const DROPS::Point3DCL& p, double t )
{
    if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
        return DROPS::Point3DCL(0,0,0);
    }
    double a=-1;
    DROPS::Point3DCL gen_c;
    gen_c[0] = -pow(0.1e1 + a * exp(-(double) (6 * t)), 0.2e1) * p[0] * (pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    gen_c[1] = pow(0.1e1 + a * exp(-(double) (6 * t)), 0.2e1) * pow(p[0], 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * p[1];
    gen_c[2] = pow(0.1e1 + a * exp(-(double) (6 * t)), 0.2e1) * pow(p[0], 0.2e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    return  gen_c;
}

DROPS::Point3DCL Test_A_plus_M_RhsVectorFun9 (const DROPS::Point3DCL& p, double t)
{
    DROPS::Point3DCL rhs =
            Test_A_plus_M_rhs_inertial_Fun9(p,t)+
            Test_A_plus_M_rhs_convFun9(p,t)+
            ParameterNS::nu*Test_A_plus_M_RhsLaplacianVectorFun9(p,t) +
            Test_A_plus_M_RhsPressureVectorFun9(p,t);
    return rhs;
}

double Test_A_plus_M_rhs2Fun9 (const DROPS::Point3DCL& p, double)
{
	return 0;

}


DROPS::Point3DCL Test_A_plus_M_RhsLaplacianVectorFun95 (const DROPS::Point3DCL& p, double t )
{
    if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
        return DROPS::Point3DCL(0,0,0);
    }
    double a=-1.;
    DROPS::Point3DCL gen_c;
    gen_c[0] = 0.2e1 * (0.1e1 + a * exp(-(double) (6 * t))) * (pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    gen_c[1] = -0.2e1 * (0.1e1 + a * exp(-(double) (6 * t))) * p[0] * p[1] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    gen_c[2] = -0.2e1 * (0.1e1 + a * exp(-(double) (6 * t))) * p[0] * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    return gen_c;
}

DROPS::Point3DCL Test_A_plus_M_RhsVectorFun95 (const DROPS::Point3DCL& p, double t)
{
    DROPS::Point3DCL rhs =
            Test_A_plus_M_rhs_inertial_Fun95(p,t)+
            Test_A_plus_M_rhs_convFun95(p,t)+
            ParameterNS::nu*Test_A_plus_M_RhsLaplacianVectorFun95(p,t) +
            Test_A_plus_M_RhsPressureVectorFun9(p,t);
    return rhs;
}

double Test_A_plus_M_rhs2Fun95 (const DROPS::Point3DCL& p, double t )
{
    double a=-1.;
    return (0.2e1 * (0.1e1 + a * exp(-(double) (6 * t))) * p[0] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)));

}

double Test_A_plus_M_pSolScalarFun9 (const DROPS::Point3DCL& p, double)
{

	return ( p[0] * pow(p[1], 0.3e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) + p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1));
}



DROPS::Point3DCL Test_A_plus_M_vSolVectorFun9_Gradient1 (const DROPS::Point3DCL& p, double t )
{
    DROPS::Point3DCL gen_c;
    double a=-1.;
    gen_c[0] = -(0.1e1 + a * exp(-(double) (6 * t))) * (pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1) * p[0];
    gen_c[1] = (0.1e1 + a * exp(-(double) (6 * t))) * p[1] * (0.2e1 * pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + 0.3e1 * pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1);
    gen_c[2] = -(0.1e1 + a * exp(-(double) (6 * t))) * p[2] * (0.2e1 * pow(p[0], 0.2e1) + 0.3e1 * pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1);
    return gen_c;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun9_Gradient2 (const DROPS::Point3DCL& p, double t )
{
    DROPS::Point3DCL gen_c;
    double a=-1.;
    gen_c[0] = -(0.1e1 + a * exp(-(double) (6 * t))) * p[1] * (pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1);
    gen_c[1] = -(0.1e1 + a * exp(-(double) (6 * t))) * p[0] * (pow(p[0], 0.2e1) + pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1);
    gen_c[2] = (0.1e1 + a * exp(-(double) (6 * t))) * p[0] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1) * p[1] * p[2];
    return gen_c;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun9_Gradient3 (const DROPS::Point3DCL& p, double t )
{
    DROPS::Point3DCL gen_c;
    double a=-1.;
    gen_c[0] = (0.1e1 + a * exp(-(double) (6 * t))) * p[2] * (pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1);
    gen_c[1] = -(0.1e1 + a * exp(-(double) (6 * t))) * p[0] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1) * p[1] * p[2];
    gen_c[2] = (0.1e1 + a * exp(-(double) (6 * t))) * p[0] * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1);
    return gen_c;
}


DROPS::Point3DCL Test_A_plus_M_vSolVectorFun95_Gradient1 (const DROPS::Point3DCL& p, double t )
{
    DROPS::Point3DCL gen_c;
    double a=-1.;

    gen_c[0] = -0.2e1 * (0.1e1 + a * exp(-(double) (6 * t))) * p[0] * (pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    gen_c[1] = 0.2e1 * (0.1e1 + a * exp(-(double) (6 * t))) * pow(p[0], 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * p[1];
    gen_c[2] = 0.2e1 * (0.1e1 + a * exp(-(double) (6 * t))) * pow(p[0], 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * p[2];
   return gen_c;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun95_Gradient2 (const DROPS::Point3DCL& p, double t )
{
    DROPS::Point3DCL gen_c;
    double a=-1.;

    gen_c[0] = (0.1e1 + a * exp(-(double) (6 * t))) * p[1] * (pow(p[0], 0.2e1) - pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    gen_c[1] = -(0.1e1 + a * exp(-(double) (6 * t))) * p[0] * (pow(p[0], 0.2e1) - pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    gen_c[2] = 0.2e1 * (0.1e1 + a * exp(-(double) (6 * t))) * p[0] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * p[1] * p[2];
   return gen_c;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun95_Gradient3 (const DROPS::Point3DCL& p, double t )
{
    DROPS::Point3DCL gen_c;
        double a=-1.;

    gen_c[0] = (0.1e1 + a * exp(-(double) (6 * t))) * p[2] * (pow(p[0], 0.2e1) - pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    gen_c[1] = 0.2e1 * (0.1e1 + a * exp(-(double) (6 * t))) * p[0] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * p[1] * p[2];
    gen_c[2] = -(0.1e1 + a * exp(-(double) (6 * t))) * p[0] * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
     return gen_c;
}






DROPS::Point3DCL Test_A_plus_M_RhsVectorFun10 (const DROPS::Point3DCL& p, double)
{
	return DROPS::Point3DCL(0,0,0);

}


double Test_A_plus_M_rhs2Fun10 (const DROPS::Point3DCL& p, double)
{
	return 0.;

}

double Test_A_plus_M_pSolScalarFun10 (const DROPS::Point3DCL& p, double)
{

	return 0.;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun10 (const DROPS::Point3DCL& p, double)
{
	DROPS::Point3DCL omega(0, -sqrt(3/0.3141592654e1)/2, -sqrt(3/0.3141592654e1)/2);
	DROPS::Point3DCL r( p[0], p[1], p[2]);
	r *= pow( pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1) , -5.e-1);

	return DROPS::Point3DCL(omega[1]*r[2] - omega[2]*r[1] , omega[2]*r[0] - omega[0]*r[2]  , omega[0]*r[1] - omega[1]*r[0] );

}


DROPS::Point3DCL Test_A_plus_M_vInitVectorFun10 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
				return DROPS::Point3DCL(0,0,0);
			}
	DROPS::Point3DCL v,w;

   //p(-Z^2,Y,X)
	/*v[0] = (-pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - pow(p[2], 0.4e1) - pow(p[0], 0.2e1) * p[2] - p[0] * pow(p[1], 0.2e1)) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));
	    v[1] = (pow(p[0], 0.2e1) + (pow(p[2], 0.2e1) - p[2]) * p[0] + pow(p[2], 0.2e1)) * p[1] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));
	    v[2] = (pow(p[0], 0.3e1) + (pow(p[2], 0.3e1) + pow(p[1], 0.2e1)) * p[0] - pow(p[1], 0.2e1) * p[2]) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));
*/

	//1Yz+1Yy
	v[0] = sqrt(0.3e1) * (p[1] - p[2]) * pow(0.3141592654e1, -0.1e1 / 0.2e1) / (0.2e1 * pow(p[0], 0.2e1) + 0.2e1 * pow(p[1], 0.2e1) + 0.2e1 * pow(p[2], 0.2e1));
	v[1] = -sqrt(0.3e1) * p[0] * pow(0.3141592654e1, -0.1e1 / 0.2e1) / (0.2e1 * pow(p[0], 0.2e1) + 0.2e1 * pow(p[1], 0.2e1) + 0.2e1 * pow(p[2], 0.2e1));
	v[2] = sqrt(0.3e1) * p[0] * pow(0.3141592654e1, -0.1e1 / 0.2e1) / (0.2e1 * pow(p[0], 0.2e1) + 0.2e1 * pow(p[1], 0.2e1) + 0.2e1 * pow(p[2], 0.2e1));
	//2Yz+3Yz
	w[0] = -0.3e1 / 0.4e1 * pow(0.3141592654e1, -0.1e1 / 0.2e1) * (-0.2e1 * sqrt(0.5e1) * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[2] + sqrt(0.7e1) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - 0.4e1 * pow(p[2], 0.2e1))) * p[1] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
	w[1] = 0.3e1 / 0.4e1 * p[0] * pow(0.3141592654e1, -0.1e1 / 0.2e1) * (-0.2e1 * sqrt(0.5e1) * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[2] + sqrt(0.7e1) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - 0.4e1 * pow(p[2], 0.2e1))) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
	w[2] = 0;

	return 1*(v+w);

}


DROPS::Point3DCL Test_A_plus_M_vSolVectorFun10_Gradient1 (const DROPS::Point3DCL& p, double)
{
	return DROPS::Point3DCL(0,0,0);

}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun10_Gradient2 (const DROPS::Point3DCL& p, double)
{
	return DROPS::Point3DCL(0,0,0);

}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun10_Gradient3 (const DROPS::Point3DCL& p, double)
{
	return DROPS::Point3DCL(0,0,0);

}

DROPS::Point3DCL Test_A_plus_M_rhs_convFun11(const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
		return DROPS::Point3DCL(0,0,0);
	}

	DROPS::Point3DCL conv_term;

	conv_term[0] = -p[0] * (  pow(p[1], 0.2e1) + pow(p[2], 0.2e1) ) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
	conv_term[1] =  p[1] * (  pow(p[0], 0.2e1)                    ) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
	conv_term[2] =  p[2] * (  pow(p[0], 0.2e1)                    ) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);

	return  conv_term;
}

DROPS::Point3DCL Test_A_plus_M_RhsVectorFun11 (const DROPS::Point3DCL& p, double t)
{
	DROPS::Point3DCL rhs = Test_A_plus_M_RhsVectorFun6(p,t) +  Test_A_plus_M_rhs_convFun11(p,t);
	return rhs;
}



double Test_A_plus_M_pDynSolScalarFun12 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
	    	return 0;
	    }
	double normalization = 1.0/3.0;
    return -normalization + ( pow(p[1], 0.2e1) + pow(p[2], 0.2e1) ) / (0.2e1 * pow(p[0], 0.2e1) + 0.2e1 * pow(p[1], 0.2e1) + 0.2e1 * pow(p[2], 0.2e1));
}

double Test_A_plus_M_pBerSolScalarFun12 (const DROPS::Point3DCL& p, double t=0)
{
	return Test_A_plus_M_pSolScalarFun6(p, t) + Test_A_plus_M_pDynSolScalarFun12(p, t);
}


DROPS::Point3DCL Test_A_plus_M_rhs_convFun13(const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
		return DROPS::Point3DCL(0,0,0);
	}

	DROPS::Point3DCL conv_term;
	conv_term[0] = ((-0.2e1 * pow(p[1], 0.2e1) * p[2] - 0.3e1 * pow(p[2], 0.3e1) - pow(p[1], 0.2e1) + 0.2e1 * pow(p[2], 0.2e1)) * pow(p[0], 0.3e1) - 0.2e1 * pow(p[1], 0.2e1) * p[2] * (p[2] - 0.2e1) * pow(p[0], 0.2e1) + (-0.3e1 * pow(p[1], 0.2e1) * pow(p[2], 0.4e1) - 0.3e1 * pow(p[2], 0.6e1) - 0.2e1 * pow(p[1], 0.4e1) * p[2] + 0.2e1 * pow(p[2], 0.5e1) + pow(p[1], 0.4e1) - pow(p[1], 0.2e1) * pow(p[2], 0.2e1)) * p[0] + 0.3e1 * pow(p[1], 0.4e1) * pow(p[2], 0.2e1) + 0.3e1 * pow(p[1], 0.2e1) * pow(p[2], 0.4e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
	conv_term[1] = p[1] * (0.3e1 * pow(p[0], 0.2e1) * pow(p[2], 0.4e1) + 0.2e1 * pow(p[0], 0.4e1) * p[2] + 0.2e1 * pow(p[0], 0.3e1) * pow(p[2], 0.2e1) + 0.2e1 * pow(p[0], 0.2e1) * pow(p[1], 0.2e1) * p[2] - 0.2e1 * pow(p[0], 0.2e1) * pow(p[2], 0.3e1) - 0.3e1 * p[0] * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) + 0.2e1 * p[0] * pow(p[2], 0.4e1) + pow(p[1], 0.2e1) * pow(p[2], 0.3e1) + pow(p[2], 0.5e1) + pow(p[0], 0.4e1) - 0.2e1 * pow(p[0], 0.3e1) * p[2] - pow(p[0], 0.2e1) * pow(p[1], 0.2e1) + 0.4e1 * pow(p[0], 0.2e1) * pow(p[2], 0.2e1) + 0.2e1 * p[0] * pow(p[1], 0.2e1) * p[2] - 0.2e1 * p[0] * pow(p[2], 0.3e1) - pow(p[1], 0.2e1) * pow(p[2], 0.2e1) + pow(p[2], 0.4e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
	conv_term[2] = (0.3e1 * pow(p[0], 0.2e1) * pow(p[2], 0.5e1) + (-0.2e1 * pow(p[0], 0.2e1) - pow(p[1], 0.2e1)) * pow(p[2], 0.4e1) - 0.5e1 * (p[0] + 0.1e1 / 0.5e1) * pow(p[1], 0.2e1) * pow(p[2], 0.3e1) + (0.3e1 * pow(p[0], 0.4e1) + 0.2e1 * pow(p[0], 0.2e1) * pow(p[1], 0.2e1) - pow(p[1], 0.4e1) + 0.2e1 * p[0] * pow(p[1], 0.2e1)) * pow(p[2], 0.2e1) + (-0.2e1 * pow(p[0], 0.4e1) - 0.3e1 * pow(p[0], 0.2e1) * pow(p[1], 0.2e1) + pow(p[1], 0.4e1)) * p[2] - 0.2e1 * pow(p[0], 0.3e1) * pow(p[1], 0.2e1) - 0.2e1 * p[0] * pow(p[1], 0.4e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);

	return  conv_term;
}

DROPS::Point3DCL Test_A_plus_M_RhsVectorFun13 (const DROPS::Point3DCL& p, double t)
{
	DROPS::Point3DCL rhs = Test_A_plus_M_RhsVectorFun8(p,t) +  Test_A_plus_M_rhs_convFun13(p,t);
	return rhs;
}


double Test_A_plus_M_pDynSolScalarFun14 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
	    	return 0;
	    }
	double normalization = 2.0/7.0;
    return -normalization +  (pow(p[2], 0.6e1) + pow(p[1], 0.2e1) * pow(p[2], 0.4e1) + 0.2e1 * pow(p[0], 0.2e1) * pow(p[2], 0.3e1) + 0.2e1 * (p[0] + 0.1e1 / 0.2e1) * pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - 0.2e1 * p[0] * pow(p[1], 0.2e1) * p[2] + pow(p[0], 0.4e1) + 0.2e1 * pow(p[0], 0.2e1) * pow(p[1], 0.2e1)) / (0.2e1 * pow(p[0], 0.2e1) + 0.2e1 * pow(p[1], 0.2e1) + 0.2e1 * pow(p[2], 0.2e1));


}

double Test_A_plus_M_pBerSolScalarFun14 (const DROPS::Point3DCL& p, double t=0)
{
	return Test_A_plus_M_pSolScalarFun8(p, t) + Test_A_plus_M_pDynSolScalarFun14(p, t);
}

double Test_A_plus_M_pBerSolScalarFun15 (const DROPS::Point3DCL& p, double t=0)
{
	return Test_A_plus_M_pSolScalarFun8(p, t) - Test_A_plus_M_pDynSolScalarFun14(p, t);
}


double Test_A_plus_M_decay_test16 (const DROPS::Point3DCL& p, double t)
{
	double a = 1.0;
	double b = 0.0;
	double T = 1.0;
	return(a+b*exp(-t/T));
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun16 (const DROPS::Point3DCL& p, double t)
{
	DROPS::Point3DCL omega(0, 0, -sqrt(3/0.3141592654e1)/2);
	DROPS::Point3DCL r( p[0], p[1], p[2]);
	r *= Test_A_plus_M_decay_test16(p, t)  *  pow( pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1) , -5.e-1);

	return DROPS::Point3DCL(omega[1]*r[2] - omega[2]*r[1] , omega[2]*r[0] - omega[0]*r[2]  , omega[0]*r[1] - omega[1]*r[0] );
}

DROPS::Point3DCL Test_A_plus_M_rhs_convFun16(const DROPS::Point3DCL& p, double t)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
		return DROPS::Point3DCL(0,0,0);
	}
	double a = 1.0;
		double b = 0.0;
		double T = 1.0;
	DROPS::Point3DCL conv_term;
	conv_term[0] = -0.3e1 / 0.4e1 * pow(Test_A_plus_M_decay_test16(p, t), 0.2e1) * pow(p[2], 0.2e1) * p[0] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1) / 0.3141592654e1;
	conv_term[1] = -0.3e1 / 0.4e1 * pow(Test_A_plus_M_decay_test16(p, t), 0.2e1) * pow(p[2], 0.2e1) * p[1] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1) / 0.3141592654e1;
	conv_term[2] =  0.3e1 / 0.4e1 * pow(Test_A_plus_M_decay_test16(p, t), 0.2e1) * p[2] * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1) / 0.3141592654e1;
	return  conv_term + (-b/T*exp(-t/T))*Test_A_plus_M_vSolVectorFun16(p,t)/Test_A_plus_M_decay_test16(p, t);
}



double Test_A_plus_M_decay_test17 (const DROPS::Point3DCL& p, double t)
{
	double a = 1.0;
	double b = 1.0;
	double c = 1.0;
	double d = -3.0;
	double T = 30;
	return (a+b*p[2]*(c + d*exp(-t/T)));
}

double Test_A_plus_M_time_deriv_decay_test17 (const DROPS::Point3DCL& p, double t)
{
	double a = 1.0;
		double b = 1.0;
		double c = 1.0;
		double d = -3.0;
		double T = 30;
	return (-b*p[2]*d*exp(-t/T)/T );
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun17 (const DROPS::Point3DCL& p, double t)
{
	DROPS::Point3DCL omega(0, 0, -sqrt(3/0.3141592654e1)/2);
	DROPS::Point3DCL r( p[0], p[1], p[2]);
	r *= ( Test_A_plus_M_decay_test17(p, t)  *  pow( pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1) , -5.e-1));

	return DROPS::Point3DCL(omega[1]*r[2] - omega[2]*r[1] , omega[2]*r[0] - omega[0]*r[2]  , omega[0]*r[1] - omega[1]*r[0] );
}



DROPS::Point3DCL Test_A_plus_M_rhs_inertial_VectorFun17(const DROPS::Point3DCL& p, double t)
{
	double a = 1.0;
		double b = 1.0;
		double c = 1.0;
		double d = -3.0;
		double T = 30;
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
			return DROPS::Point3DCL(0,0,0);
		}
	DROPS::Point3DCL in;// = Test_A_plus_M_time_deriv_decay_test17(p,t)*Test_A_plus_M_vSolVectorFun17(p,t)/Test_A_plus_M_decay_test17(p,t);
	in[0] = -b * p[2] * d / T * exp(-t / T) * p[1] * sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) / (0.2e1 * pow(p[0], 0.2e1) + 0.2e1 * pow(p[1], 0.2e1) + 0.2e1 * pow(p[2], 0.2e1));
	in[1] = b * p[2] * d / T * exp(-t / T) * p[0] * sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) / (0.2e1 * pow(p[0], 0.2e1) + 0.2e1 * pow(p[1], 0.2e1) + 0.2e1 * pow(p[2], 0.2e1));
	in[2] = 0;
	return in;
}

DROPS::Point3DCL Test_A_plus_M_rhs_laplacian_VectorFun17(const DROPS::Point3DCL& p, double t)
{
	double a = 1.0;
		double b = 1.0;
		double c = 1.0;
		double d = -3.0;
		double T = 30;
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
			return DROPS::Point3DCL(0,0,0);
		}
	DROPS::Point3DCL viscous;
	viscous[0] = 0.2e1 *  sqrt(0.3e1) * p[1] * b * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * (c + d * exp(-t / T)) * pow(0.3141592654e1, -0.1e1 / 0.2e1);
	viscous[1] = -0.2e1 * sqrt(0.3e1) * p[0] * b * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * (c + d * exp(-t / T)) * pow(0.3141592654e1, -0.1e1 / 0.2e1);
	viscous[2] = 0;
	return viscous;
}

DROPS::Point3DCL Test_A_plus_M_rhs_convFun17(const DROPS::Point3DCL& p, double t)
{
	double a = 1.0;
		double b = 1.0;
		double c = 1.0;
		double d = -3.0;
		double T = 30;
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
		return DROPS::Point3DCL(0,0,0);
	}

	DROPS::Point3DCL conv_term;
	conv_term[0] = -0.3e1 / 0.4e1 * p[0] * pow(p[2], 0.2e1) * pow(exp(-t / T) * b * d * p[2] + b * c * p[2] + a, 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1) / 0.3141592654e1;
	conv_term[1] = -0.3e1 / 0.4e1 * pow(p[2], 0.2e1) * pow(exp(-t / T) * b * d * p[2] + b * c * p[2] + a, 0.2e1) * p[1] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1) / 0.3141592654e1;
	conv_term[2] = 0.3e1 / 0.4e1 * p[2] * pow(exp(-t / T) * b * d * p[2] + b * c * p[2] + a, 0.2e1) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1) / 0.3141592654e1;
	return  conv_term ;
}


DROPS::Point3DCL Test_A_plus_M_RhsVectorFun17 (const DROPS::Point3DCL& p, double t)
{
	DROPS::Point3DCL rhs;
	//if (t<4)
	{
		rhs = Test_A_plus_M_rhs_inertial_VectorFun17(p,t) + ParameterNS::nu*Test_A_plus_M_rhs_laplacian_VectorFun17(p,t) +  Test_A_plus_M_rhs_convFun17(p,t);
	}
	return rhs;
}


DROPS::Point3DCL Test_A_plus_M_vSolVectorFun17_Gradient1 (const DROPS::Point3DCL& p, double t)
{
	double a = 1.0;
		double b = 1.0;
		double c = 1.0;
		double d = -3.0;
		double T = 30;
	DROPS::Point3DCL gen_c;
	gen_c[0] = -0.4e1 * (a + b * p[2] * (c + d * exp(-t / T))) * p[1] * sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * pow(0.2e1 * pow(p[0], 0.2e1) + 0.2e1 * pow(p[1], 0.2e1) + 0.2e1 * pow(p[2], 0.2e1), -0.2e1) * p[0];
	gen_c[1] = (exp(-t / T) * b * d * p[2] + b * c * p[2] + a) * sqrt(0.3e1) * (pow(p[0], 0.2e1) - pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) / 0.2e1;
	gen_c[2] = sqrt(0.3e1) * p[1] * pow(0.3141592654e1, -0.1e1 / 0.2e1) * (b * d * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * exp(-t / T) + c * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * b - 0.2e1 * a * p[2]) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) / 0.2e1;
	return gen_c;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun17_Gradient2 (const DROPS::Point3DCL& p, double t)
{
	double a = 1.0;
		double b = 1.0;
		double c = 1.0;
		double d = -3.0;
		double T = 30;
	DROPS::Point3DCL gen_c;
	gen_c[0] = (exp(-t / T) * b * d * p[2] + b * c * p[2] + a) * sqrt(0.3e1) * (pow(p[0], 0.2e1) - pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) / 0.2e1;
	gen_c[1] = (exp(-t / T) * b * d * p[2] + b * c * p[2] + a) * p[1] * sqrt(0.3e1) * p[0] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1);
	gen_c[2] = -p[0] * sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * (b * d * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * exp(-t / T) + c * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * b - 0.2e1 * a * p[2]) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) / 0.2e1;
	return gen_c;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun17_Gradient3 (const DROPS::Point3DCL& p, double)
{
	return DROPS::Point3DCL(0,0,0);

}


double Test_A_plus_M_smooth_step_test18 (const DROPS::Point3DCL& p, double t)
{
	double V=0.375;
	double delta = 0.0005;
	return (V*std::tanh((p[2])/delta));
}




DROPS::Point3DCL Test_A_plus_M_point_vortex_test18 (const DROPS::Point3DCL& x0,const DROPS::Point3DCL& p, double t)
{
	DROPS::Point3DCL omega(x0[0]*sqrt(3/0.3141592654e1)/2, x0[1]*sqrt(3/0.3141592654e1)/2, x0[2]*sqrt(3/0.3141592654e1)/2);
	DROPS::Point3DCL r( p[0], p[1], p[2]);
	double delta=0.0005;
	double inner_prod = (p[0]*x0[0]+p[1]*x0[1]+p[2]*x0[2]);
	double gamma = std::acos( inner_prod
							  * pow( ( pow( p[0], 0.2e1) + pow( p[1], 0.2e1) + pow( p[2], 0.2e1))
							       * ( pow(x0[0], 0.2e1) + pow(x0[1], 0.2e1) + pow(x0[2], 0.2e1))
								   , -5.e-1
								   )
						    );

	r *= std::exp(- pow(gamma/delta, 2.0));

	r *=  pow( pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1) , -5.e-1);

	return DROPS::Point3DCL(omega[1]*r[2] - omega[2]*r[1] , omega[2]*r[0] - omega[0]*r[2]  , omega[0]*r[1] - omega[1]*r[0] );
}


DROPS::Point3DCL Test_A_plus_M_disturbed_test18 (const DROPS::Point3DCL& p, double t)
{
	double c=0.001;
	DROPS::Point3DCL v = - Test_A_plus_M_point_vortex_test18(Point3DCL(1.0,0.0,0.0),p,t)
					  // + Test_A_plus_M_point_vortex_test18(Point3DCL(0.0,1.0,0.0),p,t)
					   //+ Test_A_plus_M_point_vortex_test18(Point3DCL(-1.0,0.0,0.0),p,t)
					  // + Test_A_plus_M_point_vortex_test18(Point3DCL(0.0,-1.0,0.0),p,t)
					  // - Test_A_plus_M_point_vortex_test18(Point3DCL(  0.5,  0.866025403784,0.0),p,t)
					   // - Test_A_plus_M_point_vortex_test18(Point3DCL(  0.5, -0.866025403784,0.0),p,t)
					//   + Test_A_plus_M_point_vortex_test18(Point3DCL(-0.70710678118, 0.70710678118,0.0),p,t)
					  // + Test_A_plus_M_point_vortex_test18(Point3DCL( 0.70710678118,-0.70710678118,0.0),p,t)

					   ;

	v*=c;
	return (v);
}

DROPS::Point3DCL Test_A_plus_M_base_test18 (const DROPS::Point3DCL& p, double t)
{
	DROPS::Point3DCL omega(0, 0, sqrt(3/0.3141592654e1)/2);
	DROPS::Point3DCL r( p[0], p[1], p[2]);
	r *= ( Test_A_plus_M_smooth_step_test18(p, t)  *  pow( pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1) , -5.e-1));

	return DROPS::Point3DCL(omega[1]*r[2] - omega[2]*r[1] , omega[2]*r[0] - omega[0]*r[2]  , omega[0]*r[1] - omega[1]*r[0] );
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun18 (const DROPS::Point3DCL& p, double t)
{

	return (
			Test_A_plus_M_base_test18(p,t)
			+ Test_A_plus_M_disturbed_test18(p,t)
			)
			;
}


double Test_A_plus_M_decay_test20 (const DROPS::Point3DCL& p, double t)
{
	double a = 2.0;
			double b = 4.0;
			double c = -10.0;
			double T = 1;
	return (a+p[2]*(b + c*exp(-t/T)));
}

double Test_A_plus_M_time_deriv_decay_test20 (const DROPS::Point3DCL& p, double t)
{
	double a = 2.0;
			double b = 4.0;
			double c = -10.0;
			double T = 1;
	return (-p[2] * c * exp(-t / T) / T );
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun20 (const DROPS::Point3DCL& p, double t)
{
	DROPS::Point3DCL omega(0, 0, -sqrt(3/0.3141592654e1)/2);
	DROPS::Point3DCL r( p[0], p[1], p[2]);
	r *= ( p[0]*Test_A_plus_M_decay_test20(p, t)  *  pow( pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1) , -5.e-1));

	return DROPS::Point3DCL(omega[1]*r[2] - omega[2]*r[1] , omega[2]*r[0] - omega[0]*r[2]  , omega[0]*r[1] - omega[1]*r[0] );
}



DROPS::Point3DCL Test_A_plus_M_rhs_inertial_VectorFun20(const DROPS::Point3DCL& p, double t)
{
	double a = 2.0;
			double b = 4.0;
			double c = -10.0;
			double T = 1;
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
			return DROPS::Point3DCL(0,0,0);
		}
	DROPS::Point3DCL in;
	in[0] = -p[2] * c / T * exp(-t / T) * p[0] * p[1] * sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) / (0.2e1 * pow(p[0], 0.2e1) + 0.2e1 * pow(p[1], 0.2e1) + 0.2e1 * pow(p[2], 0.2e1));
	in[1] = p[2] * c / T * exp(-t / T) * pow(p[0], 0.2e1) * sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) / (0.2e1 * pow(p[0], 0.2e1) + 0.2e1 * pow(p[1], 0.2e1) + 0.2e1 * pow(p[2], 0.2e1));
	in[2] = 0;
	return in;
}

DROPS::Point3DCL Test_A_plus_M_rhs_laplacian_VectorFun20(const DROPS::Point3DCL& p, double t)
{	double a = 2.0;
double b = 4.0;
double c = -10.0;
double T = 1;
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
			return DROPS::Point3DCL(0,0,0);
		}
	DROPS::Point3DCL viscous;
	viscous[0] = 0.5e1 *  pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * p[0] * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[1] * (exp(-t / T) * c * p[2] + b * p[2] + 0.3e1 / 0.10e2 * a) * sqrt(0.3e1);
	viscous[1] = -0.9e1 / 0.2e1 *  pow(0.3141592654e1, -0.1e1 / 0.2e1) * (c * p[2] * (pow(p[0], 0.2e1) - pow(p[1], 0.2e1) / 0.9e1 - pow(p[2], 0.2e1) / 0.9e1) * exp(-t / T) - b * pow(p[2], 0.3e1) / 0.9e1 - a * pow(p[2], 0.2e1) / 0.9e1 + b * (pow(p[0], 0.2e1) - pow(p[1], 0.2e1) / 0.9e1) * p[2] + a * pow(p[0], 0.2e1) / 0.3e1) * sqrt(0.3e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
	viscous[2] = - pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[1] * (c * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * exp(-t / T) + b * pow(p[0], 0.2e1) + b * pow(p[1], 0.2e1) + b * pow(p[2], 0.2e1) + a * p[2]) * sqrt(0.3e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) / 0.2e1;
	return viscous;
}

DROPS::Point3DCL Test_A_plus_M_rhs_convFun20(const DROPS::Point3DCL& p, double t)
{
	double a = 2.0;
			double b = 4.0;
			double c = -10.0;
			double T = 1;
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
		return DROPS::Point3DCL(0,0,0);
	}

	DROPS::Point3DCL conv_term;
	conv_term[0] = 0.3e1 / 0.4e1 * (exp(-t / T) * c * p[2] + b * p[2] + a) * p[0] * (exp(-t / T) * c * pow(p[0], 0.2e1) * pow(p[1], 0.2e1) * p[2] - exp(-t / T) * c * pow(p[0], 0.2e1) * pow(p[2], 0.3e1) + exp(-t / T) * c * pow(p[1], 0.4e1) * p[2] + exp(-t / T) * c * pow(p[1], 0.2e1) * pow(p[2], 0.3e1) + b * pow(p[0], 0.2e1) * pow(p[1], 0.2e1) * p[2] - b * pow(p[0], 0.2e1) * pow(p[2], 0.3e1) + b * pow(p[1], 0.4e1) * p[2] + b * pow(p[1], 0.2e1) * pow(p[2], 0.3e1) + a * pow(p[0], 0.2e1) * pow(p[1], 0.2e1) - a * pow(p[0], 0.2e1) * pow(p[2], 0.2e1) + a * pow(p[1], 0.4e1) + a * pow(p[1], 0.2e1) * pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1) / 0.3141592654e1;
	conv_term[1] = -0.3e1 / 0.4e1 * pow(p[0], 0.2e1) * (exp(-t / T) * c * p[2] + b * p[2] + a) * p[1] * (c * p[2] * 	(-t / T) * pow(p[0], 0.2e1) + c * p[2] * exp(-t / T) * pow(p[1], 0.2e1) + 0.2e1 * c * pow(p[2], 0.3e1) * exp(-t / T) + b * pow(p[0], 0.2e1) * p[2] + b * pow(p[1], 0.2e1) * p[2] + 0.2e1 * b * pow(p[2], 0.3e1) + a * pow(p[0], 0.2e1) + a * pow(p[1], 0.2e1) + 0.2e1 * a * pow(p[2], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1) / 0.3141592654e1;
	conv_term[2] = 0.3e1 / 0.4e1 * pow(p[0], 0.2e1) * p[2] * pow(exp(-t / T) * c * p[2] + b * p[2] + a, 0.2e1) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1) / 0.3141592654e1;
	return  conv_term ;
}


DROPS::Point3DCL Test_A_plus_M_RhsVectorFun20 (const DROPS::Point3DCL& p, double t)
{
	//return DROPS::Point3DCL(t,t,t);
	DROPS::Point3DCL rhs;
	//if (t<4)
	{
		rhs = Test_A_plus_M_rhs_inertial_VectorFun20(p,t)
		        + ParameterNS::nu*Test_A_plus_M_rhs_laplacian_VectorFun20(p,t)
              +  Test_A_plus_M_rhs_convFun20(p,t);
	}
	return rhs;
}


DROPS::Point3DCL Test_A_plus_M_vSolVectorFun20_Gradient1 (const DROPS::Point3DCL& p, double t)
{
	double a = 2.0;
			double b = 4.0;
			double c = -10.0;
			double T = 1;
	DROPS::Point3DCL gen_c;
	gen_c[0] = -(exp(-t / T) * c * p[2] + b * p[2] + a) * p[1] * sqrt(0.3e1) * (pow(p[0], 0.2e1) - pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) / 0.2e1;
	gen_c[1] = (exp(-t / T) * c * p[2] + b * p[2] + a) * p[0] * sqrt(0.3e1) * (pow(p[0], 0.2e1) - pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) / 0.2e1;
	gen_c[2] = p[0] * pow(0.3141592654e1, -0.1e1 / 0.2e1) * (c * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * exp(-t / T) + b * pow(p[0], 0.2e1) + b * pow(p[1], 0.2e1) - b * pow(p[2], 0.2e1) - 0.2e1 * a * p[2]) * p[1] * sqrt(0.3e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) / 0.2e1;
	return gen_c;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun20_Gradient2 (const DROPS::Point3DCL& p, double t)
{
	double a = 2.0;
			double b = 4.0;
			double c = -10.0;
			double T = 1;
	DROPS::Point3DCL gen_c;
	gen_c[0] = -(exp(-t / T) * c * p[2] + b * p[2] + a) * p[0] * sqrt(0.3e1) * (pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
	gen_c[1] = (exp(-t / T) * c * p[2] + b * p[2] + a) * pow(p[0], 0.2e1) * p[1] * sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
	gen_c[2] = -pow(p[0], 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * (c * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * exp(-t / T) + b * pow(p[0], 0.2e1) + b * pow(p[1], 0.2e1) - b * pow(p[2], 0.2e1) - 0.2e1 * a * p[2]) * sqrt(0.3e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) / 0.2e1;
	return gen_c;
}

DROPS::Point3DCL Test_A_plus_M_vSolVectorFun20_Gradient3 (const DROPS::Point3DCL& p, double)
{
	return DROPS::Point3DCL(0,0,0);

}

double Test_A_plus_M_rhs2Fun20 (const DROPS::Point3DCL& p, double t)
{
	//return t;
	double a = 2.0;
		double b = 4.0;
		double c = -10.0;
		double T = 1;
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
		return(0.0);
	}

    double div;
    div =  p[1] * sqrt(0.3e1) * (exp(-t / T) * c * p[2] + b * p[2] + a) * pow(0.3141592654e1, -0.1e1 / 0.2e1) / (0.2e1 * pow(p[0], 0.2e1) + 0.2e1 * pow(p[1], 0.2e1) + 0.2e1 * pow(p[2], 0.2e1));
    return  (-1)*div;
}

double Test_A_plus_M_curlFun20 (const DROPS::Point3DCL& p, double t)
{
	//return t;
	double a = 2.0;
		double b = 4.0;
		double c = -10.0;
		double T = 1;
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
		return(0.0);
	}

    double curl;
    curl =  ((0.4e1 * sqrt(0.3e1) * exp(t / T) * b + 0.4e1 * sqrt(0.3e1) * c) * pow(p[0], 0.3e1) + (0.4e1 * sqrt(0.3e1) * exp(t / T) * b + 0.4e1 * sqrt(0.3e1) * c) * p[0] * pow(p[1], 0.2e1) - 0.3e1 * sqrt(0.3e1) * p[0] * a * p[2] * exp(t / T) + (-0.3e1 * sqrt(0.3e1) * exp(t / T) * b - 0.3e1 * sqrt(0.3e1) * c) * p[0]) / exp(t / T) * pow(0.3141592654e1, -0.1e1 / 0.2e1) / 0.2e1;
    return  curl;
}



DROPS::Point3DCL Torus_vSolVectorArnold (const DROPS::Point3DCL& p, double t)
{
    double r=0.6;
    double R=1;
    double sign =0;
    if (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) > R*R) sign=1.0;
    else sign = -1.0;

    double norm = 0.1*p[2]/r; //sin()

    DROPS::Point3DCL v (-p[1]*norm/R, p[0]*norm/R, 0);
    return v;

}


DROPS::Point3DCL Torus_vSolVectorHarmonicPhi (const DROPS::Point3DCL& p, double t)
{
   // double r=0.2;
    //double R=1;
    double Lambda=1;//pow(flower_shape(atan(p[2]/(std::sqrt(pow(p[0],2.0) + pow(p[1], 2.0)) - Torus::R))),-1.0);
    double norm = 1*(p[0]*p[0] + p[1]*p[1])/(Torus::R*Torus::R);

    DROPS::Point3DCL v (-Lambda*p[1]/(Torus::R*norm), Lambda*p[0]/(norm*Torus::R), 0);
    return v;

}

DROPS::Point3DCL Torus_vSolVectorHarmonicTheta (const DROPS::Point3DCL& p, double t)
{
    double r=0.4;
    double R=1;
    double sign =0;
    if (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) > R*R) sign=1.0;
    else sign = -1.0;

    double norm = pow(p[0]*p[0] + p[1]*p[1],0.5);

    DROPS::Point3DCL v (-p[0]*p[2]/(r*norm*norm), -p[1]*p[2]/(r*norm*norm), (norm-R)/(r*norm));
    return v;

}

DROPS::Point3DCL Torus_vSolVectorHarmonic(const DROPS::Point3DCL& p, double t)
{
    double c1=1;
    double c2=1;
    return(c1*Torus_vSolVectorHarmonicPhi(p,t) + c2*Torus_vSolVectorHarmonicTheta(p,t));
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


DROPS::Point3DCL Test_cube_madeof_edges_RhsVectorFun1 (const DROPS::Point3DCL& p, double)
{
	//if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
			return DROPS::Point3DCL(0,0,0);
	//	}
	//DROPS::Point3DCL v;
	/*v[0] = 0.2e1 * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) * (pow(p[1], 0.4e1) + (pow(p[0], 0.2e1) - 0.3e1 * pow(p[2], 0.2e1) - 0.5e1 * p[0] + 0.5e1 / 0.2e1 * p[2]) * pow(p[1], 0.2e1) + 0.7e1 * pow(p[0], 0.2e1) * pow(p[2], 0.2e1) - 0.4e1 * pow(p[2], 0.4e1) - 0.5e1 / 0.2e1 * pow(p[0], 0.2e1) * p[2] + 0.5e1 / 0.2e1 * pow(p[2], 0.3e1)) + (-pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - pow(p[2], 0.4e1) - pow(p[0], 0.2e1) * p[2] - p[0] * pow(p[1], 0.2e1)) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));
	v[1] = -0.2e1 * p[1] * ((-0.10e2 * p[0] - 0.5e1) * pow(p[2], 0.2e1) + 0.5e1 * p[0] * p[2] + p[0] * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - 0.5e1 * p[0])) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) + (pow(p[0], 0.2e1) + (pow(p[2], 0.2e1) - p[2]) * p[0] + pow(p[2], 0.2e1)) * p[1] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));
	v[2] = (0.8e1 * p[0] * pow(p[2], 0.3e1) - 0.5e1 * p[0] * pow(p[2], 0.2e1) + ((-0.14e2 * p[0] - 0.10e2) * pow(p[1], 0.2e1) - 0.14e2 * pow(p[0], 0.3e1)) * p[2] + 0.5e1 * pow(p[0], 0.3e1) + 0.5e1 * p[0] * pow(p[1], 0.2e1)) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1) + (pow(p[1], 0.2e1) * (p[0] - p[2]) + p[0] * pow(p[2], 0.3e1) + pow(p[0], 0.3e1)) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));
    */
	//return v;
}


double Test_cube_madeof_edges_rhs2Fun1 (const DROPS::Point3DCL& p, double )
{

	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
		return(0.0);
	}
	double z = std::sqrt((7.0 + std::sqrt(19.0))/3.0);
	DROPS::Point3DCL P1 (-1.0,  1.0 ,  z);
	DROPS::Point3DCL P2 ( 1.0, -1.0 , -z );
	DROPS::Point3DCL P3 (2.04, 0.0, 1.0);
	DROPS::Point3DCL P4 (0.0, -1.0, -2.04);
    double v;
    v = pow(ParameterNS::h,-2.0) *
    	(  exp( - pow(ParameterNS::h,-2.0) * pow( (P1-p).norm(), 2.0) )
    	 - exp( - pow(ParameterNS::h,-2.0) * pow( (P2-p).norm(), 2.0) )
		);
    	//+ exp(-pow( (P3-p).norm(), 2.0)) + exp(-pow( (P4-p).norm(), 2.0)));
    return  (-1)*v;
}



double Test_cube_madeof_edges_pSolScalarFun1 (const DROPS::Point3DCL& p, double)
{
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
			return(0.0);
		}
    return 0;
}

DROPS::Point3DCL Test_cube_madeof_edges_vSolVectorFun1 (const DROPS::Point3DCL& p, double)
{
	/*
	if(p[0] == 0 && p[1] == 0 && p[2] == 0) {
	*/
			return DROPS::Point3DCL(0,0,0);
	/*	}
    DROPS::Point3DCL v;
    v[0] = (-pow(p[1], 0.2e1) * pow(p[2], 0.2e1) - pow(p[2], 0.4e1) - pow(p[0], 0.2e1) * p[2] - p[0] * pow(p[1], 0.2e1)) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));
    v[1] = (pow(p[0], 0.2e1) + (pow(p[2], 0.2e1) - p[2]) * p[0] + pow(p[2], 0.2e1)) * p[1] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));
    v[2] = (pow(p[0], 0.3e1) + (pow(p[2], 0.3e1) + pow(p[1], 0.2e1)) * p[0] - pow(p[1], 0.2e1) * p[2]) / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1));

    return v;*/
}

DROPS::Point3DCL Test_cube_madeof_edges_vSolVectorFun1_Gradient1 (const DROPS::Point3DCL& p, double)
{

	return DROPS::Point3DCL(0,0,0);

}

DROPS::Point3DCL Test_cube_madeof_edges_vSolVectorFun1_Gradient2 (const DROPS::Point3DCL& p, double)
{
	return DROPS::Point3DCL(0,0,0);

}

DROPS::Point3DCL Test_cube_madeof_edges_vSolVectorFun1_Gradient3 (const DROPS::Point3DCL& p, double)
{

	return DROPS::Point3DCL(0,0,0);

}

double Test1_chiSol (const DROPS::Point3DCL& p, double)
{
	return ((1./2.)*(sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) / 0.2e1));

}

double Test1_omegaSol (const DROPS::Point3DCL& p, double)
{
	return ((sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) / 0.2e1));

}

double Test1_rhs3 (const DROPS::Point3DCL& p, double)
{

	return (sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1));

}

double Test2_chiSol (const DROPS::Point3DCL& p, double t)
{
	return (((1-std::exp(-4*t))/2.)*(sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) / 0.2e1));

}

double Test2_omegaSol (const DROPS::Point3DCL& p, double t)
{
	return ((1-std::exp(-4*t))*(sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) / 0.2e1));

}

double Test2_rhs3 (const DROPS::Point3DCL& p, double )
{

	return (sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1));

}


double Test3_chiSol (const DROPS::Point3DCL& p, double t)
{
	return (((1-std::exp(-4*t))/2.)*(sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) / 0.2e1));

}

double Test3_omegaSol (const DROPS::Point3DCL& p, double t)
{
	return ((1-std::exp(-4*t))*(sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) / 0.2e1));

}

double Test3_rhs3 (const DROPS::Point3DCL& p, double  t)
{
 /*double sigma=1.;
 double ParameterNS::eps =1.;*/


 return (exp(-(double) (4 * t)) * sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) / ParameterNS::eps - 0.9e1 / 0.512e3 * ((pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - 0.4e1 / 0.3e1 * pow(p[2], 0.2e1)) * sqrt(0.3e1) * (-0.1e1 + exp(-(double) (4 * t))) * p[2] * ParameterNS::eps * pow(0.3141592654e1, 0.3e1 / 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.7e1 / 0.2e1) + pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.3e1) * (pow(p[2], 0.2e1) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * pow(exp(-(double) (4 * t)), 0.2e1) - 0.2e1 * pow(p[2], 0.2e1) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * exp(-(double) (4 * t)) + (-0.16e2 / 0.3e1 * ParameterNS::eps * ParameterNS::eps * 0.3141592654e1 - 0.1e1) * pow(p[2], 0.4e1) - 0.8e1 / 0.3e1 * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * (ParameterNS::eps * ParameterNS::eps * 0.3141592654e1 - 0.3e1 / 0.8e1) * pow(p[2], 0.2e1) + 0.8e1 / 0.3e1 * ParameterNS::eps * ParameterNS::eps * 0.3141592654e1 * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1), 0.2e1)) * 0.3141592654e1 / 0.2e1) * sqrt(pow(p[2], 0.2e1) * pow(-0.1e1 + exp(-(double) (4 * t)), 0.2e1) * pow((-0.1e1 + exp(-(double) (4 * t))) * sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) * p[2] + 0.4e1 * ParameterNS::eps, 0.2e1) / 0.3141592654e1 / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(ParameterNS::eps, -0.4e1)) * ParameterNS::sigma * pow(0.3141592654e1, -0.3e1 / 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.9e1 / 0.2e1) * sqrt(0.16e2) * pow(ParameterNS::eps * sqrt(0.3141592654e1) * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) + sqrt(0.3e1) * p[2] * (-0.1e1 + exp(-(double) (4 * t))) / 0.4e1, -0.2e1) / p[2]);
}

double Test3_rhs4 (const DROPS::Point3DCL& p, double  t)
{
 /*   double ParameterNS::sigma=1.;
    double ParameterNS::eps =1.;*/

    return ((0.1e1 - exp(-(double) (4 * t))) * sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) / 0.2e1 + (-0.1e1 + exp(-(double) (4 * t))) * sqrt(0.3e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) / 0.2e1);
}

double Test4_rhs3 (const DROPS::Point3DCL& p, double  t)
{
    /*double ParameterNS::sigma=0.01;
    double ParameterNS::eps =1.;*/
    double delta=1;
    return (exp(-(double) (4 * t)) * sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) / ParameterNS::eps + sqrt(pow((-0.1e1 + exp(-(double) (4 * t))) * sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) * p[2] + 0.4e1 * ParameterNS::eps, 0.2e1) * pow(-0.1e1 + exp(-(double) (4 * t)), 0.2e1) * pow(p[2], 0.2e1) / 0.3141592654e1 / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(ParameterNS::eps, -0.4e1)) * ParameterNS::sigma * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.11e2 / 0.2e1) * sqrt(0.16e2) * (sqrt(0.16e2) * sqrt(0.3e1) * 0.3141592654e1 * 0.3141592654e1 * delta * pow(ParameterNS::eps, 0.4e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.6e1) * sqrt(pow((-0.1e1 + exp(-(double) (4 * t))) * sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) * p[2] + 0.4e1 * ParameterNS::eps, 0.2e1) * pow(-0.1e1 + exp(-(double) (4 * t)), 0.2e1) * pow(p[2], 0.2e1) / 0.3141592654e1 / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(ParameterNS::eps, -0.4e1)) - 0.9e1 / 0.2e1 * (sqrt(0.3141592654e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.7e1 / 0.2e1) * (-0.1e1 + exp(-(double) (4 * t))) * sqrt(0.3e1) * p[2] * (-0.2e1 / 0.3e1 * pow(p[2], 0.4e1) + (pow(p[0], 0.2e1) / 0.3e1 - 0.4e1 / 0.3e1 * pow(p[1], 0.2e1)) * pow(p[2], 0.2e1) + pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1), 0.2e1)) * ParameterNS::eps - 0.2e1 / 0.3e1 * (pow(p[2], 0.2e1) + pow(p[0], 0.2e1) - 0.3e1 / 0.2e1 * pow(p[1], 0.2e1)) * sqrt(0.3141592654e1) * (-0.1e1 + exp(-(double) (4 * t))) * sqrt(0.3e1) * pow(p[2], 0.3e1) * ParameterNS::eps * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.7e1 / 0.2e1) + pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.4e1) * (pow(p[2], 0.2e1) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * pow(exp(-(double) (4 * t)), 0.2e1) - 0.2e1 * pow(p[2], 0.2e1) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) - pow(p[2], 0.2e1)) * exp(-(double) (4 * t)) + (-0.16e2 / 0.3e1 * 0.3141592654e1 * ParameterNS::eps * ParameterNS::eps - 0.1e1) * pow(p[2], 0.4e1) - 0.8e1 / 0.3e1 * (0.3141592654e1 * ParameterNS::eps * ParameterNS::eps - 0.3e1 / 0.8e1) * (pow(p[0], 0.2e1) + pow(p[1], 0.2e1)) * pow(p[2], 0.2e1) + 0.8e1 / 0.3e1 * 0.3141592654e1 * ParameterNS::eps * ParameterNS::eps * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1), 0.2e1)) / 0.2e1) * pow(-0.1e1 + exp(-(double) (4 * t)), 0.2e1)) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * pow(-0.1e1 + exp(-(double) (4 * t)), -0.2e1) * pow(ParameterNS::eps * sqrt(0.3141592654e1) * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) + p[2] * sqrt(0.3e1) * (-0.1e1 + exp(-(double) (4 * t))) / 0.4e1, -0.2e1) / p[2] / 0.256e3);
}



double Test5_chiSol (const DROPS::Point3DCL& p, double tt) {
    double a = 0.4;
    double E = std::exp(1.0);
    double t=tt*ParameterNS::sigma/ParameterNS::eps;
    if (a < 0.5)
        return (4 * a - 4 * std::pow(a, 2) + std::pow(E, t / 2.) - 4 * a * std::pow(E, t / 2.) + \
4 * std::pow(a, 2) * std::pow(E, t / 2.) - std::sqrt((-4 * (-1 + \
a) * a * std::pow(E, t / 2.)) / std::pow(1 - 2 * a, 2) + std::pow(E, t)) + \
4 * a * std::sqrt((-4 * (-1 + a) * a * std::pow(E, t / 2.)) / std::pow(1 - 2 * a, 2) + \
std::pow(E, t)) - 4 * std::pow(a, 2) * std::sqrt((-4 * (-1 + \
a) * a * std::pow(E, t / 2.)) / std::pow(1 - 2 * a, 2) + \
std::pow(E, t))) / (2. * (std::pow(E, t / 2.) - 4 * a * (-1 + std::pow(E, t / 2.)) + \
4 * std::pow(a, 2) * (-1 + std::pow(E, t / 2.))));
    else if (a >= 0.5)
        return (4 * a - 4 * std::pow(a, 2) + std::pow(E, t / 2.) - 4 * a * std::pow(E, t / 2.) + \
4 * std::pow(a, 2) * std::pow(E, t / 2.) + std::sqrt((-4 * (-1 + \
a) * a * std::pow(E, t / 2.)) / std::pow(1 - 2 * a, 2) + std::pow(E, t)) - \
4 * a * std::sqrt((-4 * (-1 + a) * a * std::pow(E, t / 2.)) / std::pow(1 - 2 * a, 2) + \
std::pow(E, t)) + 4 * std::pow(a, 2) * std::sqrt((-4 * (-1 + \
a) * a * std::pow(E, t / 2.)) / std::pow(1 - 2 * a, 2) + \
std::pow(E, t))) / (2. * (std::pow(E, t / 2.) - 4 * a * (-1 + std::pow(E, t / 2.)) + \
4 * std::pow(a, 2) * (-1 + std::pow(E, t / 2.))));
}

double Test6_chiSol (const DROPS::Point3DCL& p, double t)
{
    double a=0.8;
    double T=0.25;
    return ((0.1e1 - a * exp(-t / T)) * (sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] *
                                         pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) / 0.2e1 +
                                         0.1e1 / 0.2e1));
}

double Test6_omegaSol (const DROPS::Point3DCL& p, double t)
{
    double a=0.8;
    double T=0.25;
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

double Test6_rhs3 (const DROPS::Point3DCL& p, double  t) {
    double a=0.8;
    double T=0.25;
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

double Test7_chiSol (const DROPS::Point3DCL& p, double t) {
    double a=0.8;
    double T=0.25;
    return ((0.1e1 - a * exp(-t / T)) * (sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] *
                                 pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) / 0.2e1 +
                                 0.1e1 / 0.2e1));
}

double Test7_rhs3 (const DROPS::Point3DCL& p, double  t) {
    double a=0.8;
    double T=0.25;
    double inertia=a / T * exp(-t / T) * (sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] *
                            pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1)
                                          / 0.2e1 + 0.1e1 / 0.2e1);
    double diffusion=-ParameterNS::sigma * ParameterNS::eps * (a * exp(-t / T) - 0.1e1) * sqrt(0.3e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1);
    double well_potential=-ParameterNS::sigma / ParameterNS::eps * (a * exp(-t / T) - 0.1e1) * (sqrt(0.3e1) * p[2] +
                            sqrt(0.3141592654e1) * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1))) *
                          (exp(-t / T) * sqrt(0.3e1) * a * p[2] + exp(-t / T) * sqrt(0.3141592654e1) * sqrt(pow(p[0], 0.2e1) +
                                                                                                            pow(p[1], 0.2e1) +
                                                                                                            pow(p[2], 0.2e1)) * a -
                           sqrt(0.3e1) * p[2] + sqrt(0.3141592654e1) * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)))
                          * (exp(-t / T) * sqrt(0.3e1) * a * p[2] + exp(-t / T) * sqrt(0.3141592654e1) * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * a - sqrt(0.3e1) * p[2])
                          * pow(0.3141592654e1, -0.3e1 / 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1) / 0.8e1;
    return(inertia+diffusion+well_potential);
}

double Test8_rhs3 (const DROPS::Point3DCL& p, double  t) {
    double a=0.8;
    double T=0.25;
    double inertia=a / T * exp(-t / T) * (sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) / 0.2e1 + 0.1e1 / 0.2e1);
    double diffusion=-(ParameterNS::sigma * ParameterNS::eps) * (a * exp(-t / T) - 0.1e1) * sqrt(0.3e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1);
    double well_potential_linear=-(ParameterNS::sigma / ParameterNS::eps) * (a * exp(-t / T) - 0.1e1) * (sqrt(0.3e1) * p[2] + sqrt(0.3141592654e1) * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1))) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) / 0.2e1;
    double well_potential_qubic=-(ParameterNS::sigma / ParameterNS::eps) *pow(a * exp(-t / T) - 0.1e1, 0.3e1) * pow(sqrt(0.3e1) * p[2] + sqrt(0.3141592654e1) * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)), 0.3e1) * pow(0.3141592654e1, -0.3e1 / 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1) / 0.8e1;
    double well_potential_parabolic = (ParameterNS::sigma / ParameterNS::eps) * pow(a * exp(-t / T) - 0.1e1, 0.2e1) * pow(sqrt(0.3e1) * p[2] + sqrt(0.3141592654e1) * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)), 0.2e1) / 0.3141592654e1 / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) / 0.4e1;
    double well_potential =  -(ParameterNS::sigma / ParameterNS::eps) * (a * exp(-t / T) - 0.1e1) * (sqrt(0.3e1) * p[2] + sqrt(0.3141592654e1) * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1))) * (exp(-t / T) * sqrt(0.3e1) * a * p[2] + exp(-t / T) * sqrt(0.3141592654e1) * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * a - sqrt(0.3e1) * p[2] + sqrt(0.3141592654e1) * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1))) * (exp(-t / T) * sqrt(0.3e1) * a * p[2] + exp(-t / T) * sqrt(0.3141592654e1) * sqrt(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * a - sqrt(0.3e1) * p[2]) * pow(0.3141592654e1, -0.3e1 / 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1 / 0.2e1) / 0.8e1;

    return(inertia+diffusion -(0.)*(1.*3./2.)*well_potential_parabolic+(1.)*0.5*well_potential_linear+0.*well_potential_qubic);//
}

double  Test9_chiSol (const DROPS::Point3DCL& p, double  t) {

    //return(0.5+2.*Test1_chiSol(p,t));
    //return(std::cos(p[0]*p[1]*pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1 / 0.2e1)));
    return((double)rand() / RAND_MAX);

    //return(std::sin(p[0])*std::sin(p[1])*std::sin(p[0])*std::sin(p[1]));
    //return(0.55);
}

double Test10_chiSol (const DROPS::Point3DCL& p, double t)
{
    double a=0.8;
    double T=0.25;
    return ( (0.1e1 - a * exp(-t / T)) * (p[0] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[1] + 0.1e1) / 0.2e1);
}

double Test10_omegaSol (const DROPS::Point3DCL& p, double t)
{
    double a=0.8;
    double T=0.25;
    double diffusion=-0.3e1 * ParameterNS::eps * (a * exp(-t / T) - 0.1e1) * p[1] * p[0] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    double well_potential= -0.1e1 / ParameterNS::eps * (a * exp(-t / T) - 0.1e1) * (pow(p[0], 0.2e1) + p[0] * p[1] + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * (exp(-t / T) * a * pow(p[0], 0.2e1) + exp(-t / T) * a * p[0] * p[1] + exp(-t / T) * a * pow(p[1], 0.2e1) + exp(-t / T) * a * pow(p[2], 0.2e1) + pow(p[0], 0.2e1) - p[0] * p[1] + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * (exp(-t / T) * a * pow(p[0], 0.2e1) + exp(-t / T) * a * p[0] * p[1] + exp(-t / T) * a * pow(p[1], 0.2e1) + exp(-t / T) * a * pow(p[2], 0.2e1) - p[0] * p[1]) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1) / 0.8e1; ;
    return(diffusion+well_potential);
}

double Test10_rhs3 (const DROPS::Point3DCL& p, double  t) {
    double a=0.8;
    double T=0.25;
    double inertia=a / T * exp(-t / T) * (p[0] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[1] + 0.1e1) / 0.2e1;
    double bidiffusion = -0.3e1 / 0.16e2 * ParameterNS::sigma * (a * (pow(p[0], 0.2e1) + p[0] * p[1] + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * exp(-t / T) + pow(p[0], 0.2e1) - p[0] * p[1] + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * pow(a * exp(-t / T) - 0.1e1, 0.3e1) * (pow(p[0], 0.2e1) + p[0] * p[1] + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * (pow(a, 0.3e1) * pow(pow(p[0], 0.2e1) + p[0] * p[1] + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.3e1) * (pow(p[1], 0.4e1) - 0.3e1 / 0.2e1 * pow(p[1], 0.3e1) * p[0] + (-0.7e1 / 0.2e1 * pow(p[0], 0.2e1) + pow(p[2], 0.2e1)) * pow(p[1], 0.2e1) + (-0.3e1 / 0.2e1 * pow(p[0], 0.3e1) - 0.3e1 / 0.2e1 * p[0] * pow(p[2], 0.2e1)) * p[1] + pow(p[0], 0.4e1) + pow(p[0], 0.2e1) * pow(p[2], 0.2e1)) * pow(exp(-t / T), 0.3e1) - 0.9e1 / 0.2e1 * a * a * pow(pow(p[0], 0.2e1) + p[0] * p[1] + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), 0.2e1) * p[0] * (pow(p[1], 0.4e1) - 0.2e1 / 0.3e1 * pow(p[1], 0.3e1) * p[0] + (-0.5e1 / 0.3e1 * pow(p[0], 0.2e1) + 0.4e1 / 0.3e1 * pow(p[2], 0.2e1)) * pow(p[1], 0.2e1) + (-0.2e1 / 0.3e1 * pow(p[0], 0.3e1) - 0.2e1 / 0.3e1 * p[0] * pow(p[2], 0.2e1)) * p[1] + pow(p[0], 0.4e1) + 0.4e1 / 0.3e1 * pow(p[0], 0.2e1) * pow(p[2], 0.2e1) + pow(p[2], 0.4e1) / 0.3e1) * p[1] * pow(exp(-t / T), 0.2e1) - 0.2e1 / 0.3e1 * a * (pow(p[1], 0.8e1) - 0.3e1 / 0.4e1 * pow(p[1], 0.7e1) * p[0] + (-(double) (6 * ParameterNS::eps * ParameterNS::eps) + 0.3e1 * pow(p[2], 0.2e1) - 0.39e2 / 0.4e1 * pow(p[0], 0.2e1)) * pow(p[1], 0.6e1) - 0.9e1 / 0.4e1 * p[0] * (-(double) (8 * ParameterNS::eps * ParameterNS::eps) + pow(p[2], 0.2e1)) * pow(p[1], 0.5e1) + (0.3e1 * pow(p[2], 0.4e1) + (-(double) (12 * ParameterNS::eps * ParameterNS::eps) - 0.14e2 * pow(p[0], 0.2e1)) * pow(p[2], 0.2e1) + 0.24e2 * pow(p[0], 0.2e1) * (double) (ParameterNS::eps * ParameterNS::eps) + 0.13e2 / 0.4e1 * pow(p[0], 0.4e1)) * pow(p[1], 0.4e1) - 0.9e1 / 0.4e1 * p[0] * (p[2] - (double) (4 * ParameterNS::eps)) * (p[2] + (double) (4 * ParameterNS::eps)) * (pow(p[0], 0.2e1) + pow(p[2], 0.2e1)) * pow(p[1], 0.3e1) - 0.39e2 / 0.4e1 * (-0.4e1 / 0.39e2 * pow(p[2], 0.4e1) + (0.17e2 / 0.39e2 * pow(p[0], 0.2e1) + 0.8e1 / 0.13e2 * (double) ParameterNS::eps * (double) ParameterNS::eps) * pow(p[2], 0.2e1) + pow(p[0], 0.4e1) - 0.32e2 / 0.13e2 * pow(p[0], 0.2e1) * (double) (ParameterNS::eps * ParameterNS::eps)) * (pow(p[0], 0.2e1) + pow(p[2], 0.2e1)) * pow(p[1], 0.2e1) - 0.3e1 / 0.4e1 * p[0] * pow(pow(p[0], 0.2e1) + pow(p[2], 0.2e1), 0.2e1) * (-(double) (24 * ParameterNS::eps * ParameterNS::eps) + pow(p[0], 0.2e1) + pow(p[2], 0.2e1)) * p[1] + pow(p[0], 0.2e1) * pow(pow(p[0], 0.2e1) + pow(p[2], 0.2e1), 0.2e1) * (-(double) (6 * ParameterNS::eps * ParameterNS::eps) + pow(p[0], 0.2e1) + pow(p[2], 0.2e1))) * (pow(p[0], 0.2e1) + p[0] * p[1] + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * exp(-t / T) + 0.7e1 / 0.6e1 * p[0] * (pow(p[1], 0.8e1) + (-0.96e2 / 0.7e1 * (double) ParameterNS::eps * (double) ParameterNS::eps + 0.24e2 / 0.7e1 * pow(p[2], 0.2e1) - 0.6e1 / 0.7e1 * pow(p[0], 0.2e1)) * pow(p[1], 0.6e1) + (0.30e2 / 0.7e1 * pow(p[2], 0.4e1) + (-0.264e3 / 0.7e1 * (double) ParameterNS::eps * (double) ParameterNS::eps + 0.10e2 / 0.7e1 * pow(p[0], 0.2e1)) * pow(p[2], 0.2e1) - 0.120e3 / 0.7e1 * pow(p[0], 0.2e1) * (double) (ParameterNS::eps * ParameterNS::eps) + pow(p[0], 0.4e1)) * pow(p[1], 0.4e1) - 0.6e1 / 0.7e1 * (-0.8e1 / 0.3e1 * pow(p[2], 0.4e1) + (-0.8e1 / 0.3e1 * pow(p[0], 0.2e1) + (double) (40 * ParameterNS::eps * ParameterNS::eps)) * pow(p[2], 0.2e1) + pow(p[0], 0.4e1) + 0.20e2 * pow(p[0], 0.2e1) * (double) (ParameterNS::eps * ParameterNS::eps)) * (pow(p[0], 0.2e1) + pow(p[2], 0.2e1)) * pow(p[1], 0.2e1) + (0.3e1 / 0.7e1 * pow(p[2], 0.4e1) + (0.10e2 / 0.7e1 * pow(p[0], 0.2e1) - 0.72e2 / 0.7e1 * (double) ParameterNS::eps * (double) ParameterNS::eps) * pow(p[2], 0.2e1) + pow(p[0], 0.4e1) - 0.96e2 / 0.7e1 * pow(p[0], 0.2e1) * (double) (ParameterNS::eps * ParameterNS::eps)) * pow(pow(p[0], 0.2e1) + pow(p[2], 0.2e1), 0.2e1)) * p[1]) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.8e1) * pow(pow(0.1e1 - a * exp(-t / T), 0.2e1) * pow(p[0] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[1] + 0.1e1, 0.2e1) * pow(0.1e1 - (0.1e1 - a * exp(-t / T)) * (p[0] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[1] + 0.1e1) / 0.2e1, 0.2e1), -0.1e1 / 0.2e1) / (double) ParameterNS::eps ;
    return(inertia - bidiffusion);
}

double Test11_chiSol (const DROPS::Point3DCL& p, double t) {
    double a=0.8;
    double T=0.25;
    return ( (0.1e1 - a * exp(-t / T)) * (p[0] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[1] + 0.1e1) / 0.2e1);
}

double Test12_chiSol (const DROPS::Point3DCL& pp, double t)
{
    DROPS::Point3DCL p(pp );
    return (1./2.+(1./2.)*(sqrt(0.3e1) * pow(0.3141592654e1, -0.1e1 / 0.2e1) * p[2] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.1e1 / 0.2e1) / 0.2e1));

}

double Test11_rhs3 (const DROPS::Point3DCL& p, double  t) {
    double a=0.8;
    double T=0.25;
    double inertia=a / T * exp(-t / T) * (p[0] / (pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * p[1] + 0.1e1) / 0.2e1;
    double diffusion=-0.3e1 * ParameterNS::eps * (a * exp(-t / T) - 0.1e1) * p[1] * p[0] * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.2e1);
    double well_potential= -0.1e1 / ParameterNS::eps * (a * exp(-t / T) - 0.1e1) * (pow(p[0], 0.2e1) + p[0] * p[1] + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * (exp(-t / T) * a * pow(p[0], 0.2e1) + exp(-t / T) * a * p[0] * p[1] + exp(-t / T) * a * pow(p[1], 0.2e1) + exp(-t / T) * a * pow(p[2], 0.2e1) + pow(p[0], 0.2e1) - p[0] * p[1] + pow(p[1], 0.2e1) + pow(p[2], 0.2e1)) * (exp(-t / T) * a * pow(p[0], 0.2e1) + exp(-t / T) * a * p[0] * p[1] + exp(-t / T) * a * pow(p[1], 0.2e1) + exp(-t / T) * a * pow(p[2], 0.2e1) - p[0] * p[1]) * pow(pow(p[0], 0.2e1) + pow(p[1], 0.2e1) + pow(p[2], 0.2e1), -0.3e1) / 0.8e1; ;
    return(inertia+diffusion+well_potential);
}