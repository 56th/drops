/// \file filmCoeff.cpp
/// \brief boundary and source functions for film-type problems
/// \author LNM RWTH Aachen: Christoph Lehrenfeld, Thorolf Schulte

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

#include "misc/funcmap.h"
#include "misc/params.h"

extern DROPS::ParamCL P;

//========================================================================
//                        Functions for the film problem
//========================================================================
//A remark : The functions declared here still use MeshSize from json files. But the json style is different from general one.
namespace filminflow{

    DROPS::Point3DCL FilmInflow( const DROPS::Point3DCL& p, double t)
    {
        static DROPS::Point3DCL dx;
        static bool first = true;
        //dirty hack
        if (first){
            std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
            size_t idx_;
            while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
                mesh[idx_]= ' ';
            std::istringstream brick_info( mesh);
            brick_info >> dx[0] >> dx[1] >> dx[2] ;
            first = false;
        }
	const double Ly= dx[1];
        static double PumpFreq = P.get<double>("Exp.PumpFreq");
        static double PumpAmpl = P.get<double>("Exp.PumpAmpl");
        static double Thickness= P.get<double>("Exp.Thickness");
        static double DensFluid= P.get<double>("Mat.DensFluid");
        static double ViscFluid= P.get<double>("Mat.ViscFluid");
        static double GravityX = P.get<DROPS::Point3DCL>("Exp.Gravity")[0];
        DROPS::Point3DCL ret(0.);
        const double d= p[1]/Thickness;
        static const double u= DensFluid*GravityX*Thickness*Thickness/ViscFluid/2;
        ret[0]= d<=1 ? (2*d-d*d)*u * (1 + PumpAmpl*std::sin(2*M_PI*t*PumpFreq))
                     : (Ly-p[1])/(Ly-Thickness)*u;
        return ret;
    }

    //========================================================================
    //        Registration of the function(s) in the func-container
    //========================================================================
    static DROPS::RegisterVectorFunction regvelfilminflow("FilmInflow", FilmInflow);
}

//========================================================================
//                        Functions for the levelset function
//========================================================================
namespace filmdistance{
    double WavyDistanceFct( const DROPS::Point3DCL& p, double)
    {	
        // wave length = 100 x film width
        static DROPS::Point3DCL MeshSize;
        static bool first = true;
        //dirty hack
        if (first){
            std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
            size_t idx_;
            while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
                mesh[idx_]= ' ';
            std::istringstream brick_info( mesh);
            brick_info >> MeshSize[0] >> MeshSize[1] >> MeshSize[2] ;
            first = false;
        }
        static double Ampl_zDir= P.get<double>("Exp.Ampl_zDir");
        static double PumpAmpl = P.get<double>("Exp.PumpAmpl");
        static double Thickness= P.get<double>("Exp.Thickness");
        const double wave_x= std::cos(2*M_PI*p[0]/MeshSize[0]),
    //   const double wave= std::sin(2*M_PI*p[0]/MeshSize[0]),
            wave_z= std::cos(2*M_PI*p[2]/MeshSize[2]); // z \in [-1,1]
    //    return p[1] - P.get<double>("Exp.Thickness") * (1 + P.get<double>("Exp.PumpAmpl")*wave);
    //    return p[1] - P.get<double>("Exp.Thickness") * (1 + P.get<double>("Exp.PumpAmpl")*(wave + P.get<double>("Exp.Ampl_zDir")*std::cos(z*M_PI)));
    //    const double z_fac=  (1 + Ampl_zDir/2*std::cos(z*M_PI));  // (z=+-1) 1-P.get<double>("Exp.Ampl_zDir") <= z_fac <= 1+P.get<double>("Exp.Ampl_zDir") (z=0)
    //    return p[1] - Thickness * (1 + PumpAmpl*wave) * z_fac;
        return p[1] - Thickness * ( 1 + PumpAmpl*wave_x + Ampl_zDir * wave_z);
    }
    DROPS::Point3DCL Nusselt_film( const DROPS::Point3DCL& p, double)
    {
        static DROPS::Point3DCL MeshSize;
        static bool first = true;
        //dirty hack
        if (first){
            std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
            size_t idx_;
            while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
                mesh[idx_]= ' ';
            std::istringstream brick_info( mesh);
            brick_info >> MeshSize[0] >> MeshSize[1] >> MeshSize[2] ;
            first = false;
        }
        static double Ampl_zDir= P.get<double>("Exp.Ampl_zDir");
        static double PumpAmpl = P.get<double>("Exp.PumpAmpl");
        static double Thickness= P.get<double>("Exp.Thickness");
        const double wave_x= std::cos(2*M_PI*p[0]/MeshSize[0]),
            wave_z= std::cos(2*M_PI*p[2]/MeshSize[2]);
        double delta= Thickness * ( 1 + PumpAmpl*wave_x + Ampl_zDir * wave_z);

        static double DensFluid= P.get<double>("Mat.DensFluid");
        static double ViscFluid= P.get<double>("Mat.ViscFluid");
        static double GravityX = P.get<DROPS::Point3DCL>("Exp.Gravity")[0];
        DROPS::Point3DCL ret(0.);
        const double d= p[1]/delta;
        static const double u= DensFluid * GravityX * delta * delta /ViscFluid/2;
        ret[0]= d<=1 ? (2*d-d*d)*u 
                     : 0.;
        return ret;
    }
    //========================================================================
    //        Registration of the function(s) in the func-container
    //========================================================================
    static DROPS::RegisterVectorFunction regvelfilmNusselt("NusseltFilm", Nusselt_film);
    static DROPS::RegisterScalarFunction regscafilmlset("WavyFilm", WavyDistanceFct);
}



//========================================================================
//                        Functions for matching function
//========================================================================
namespace filmperiodic{
    template<int A, int B>
    bool periodic_2sides( const DROPS::Point3DCL& p, const DROPS::Point3DCL& q)
    {
        static DROPS::Point3DCL MeshSize;
        static bool first = true;
        //dirty hack
        if (first){
            std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
            size_t idx_;
            while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
                mesh[idx_]= ' ';
            std::istringstream brick_info( mesh);
            brick_info >> MeshSize[0] >> MeshSize[1] >> MeshSize[2] ;
            first = false;
        }
        const DROPS::Point3DCL d= fabs(p-q);
        const DROPS::Point3DCL L= fabs(MeshSize);

        const int D = 3 - A - B;
        return (d[B] + d[D] < 1e-12 && std::abs( d[A] - L[A]) < 1e-12)  // dB=dD=0 and dA=LA
          ||   (d[A] + d[D] < 1e-12 && std::abs( d[B] - L[B]) < 1e-12)  // dA=dD=0 and dB=LB
          ||   (d[D] < 1e-12 && std::abs( d[A] - L[A]) < 1e-12 && std::abs( d[B] - L[B]) < 1e-12);  // dD=0 and dA=LA and dB=LB
    }

    template<int A>
    bool periodic_1side( const DROPS::Point3DCL& p, const DROPS::Point3DCL& q)
    {
        static DROPS::Point3DCL MeshSize;
        static bool first = true;
        //dirty hack
        if (first){
            std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
            size_t idx_;
            while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
                mesh[idx_]= ' ';
            std::istringstream brick_info( mesh);
            brick_info >> MeshSize[0] >> MeshSize[1] >> MeshSize[2] ;
            first = false;
        }
        const int B = (A+1)%2;
        const int D = (B+1)%2;
        const DROPS::Point3DCL d= fabs(p-q);
        const DROPS::Point3DCL L= fabs(MeshSize); //fabs(P.get<DROPS::Point3DCL>("MeshSize"));
        return (d[B] + d[D] < 1e-12 && std::abs( d[A] - L[A]) < 1e-12);
    }

    //========================================================================
    //        Registration of the function(s) in the func-container
    //========================================================================
    static DROPS::RegisterMatchingFunction regmatch2_xy("periodicxy", periodic_2sides<0,1>);
    static DROPS::RegisterMatchingFunction regmatch2_xz("periodicxz", periodic_2sides<0,2>);
    static DROPS::RegisterMatchingFunction regmatch2_yz("periodicyz", periodic_2sides<1,2>);
    static DROPS::RegisterMatchingFunction regmatch1_x("periodicx", periodic_1side<0>);
    static DROPS::RegisterMatchingFunction regmatch1_y("periodicy", periodic_1side<1>);
    static DROPS::RegisterMatchingFunction regmatch1_z("periodicz", periodic_1side<2>);
}
