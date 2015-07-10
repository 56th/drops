/// \file twophaseCoeff.cpp
/// \brief boundary and source functions for the twophasedrops-type problems
/// \author LNM RWTH Aachen: Christoph Lehrenfeld

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
//          Functions for twophasedrops-executable (inflow)
//========================================================================
namespace tpd_inflow{
    /// \name inflow condition
    DROPS::SVectorCL<3> InflowBrick( const DROPS::Point3DCL& p, double t)
    {
        DROPS::SVectorCL<3> ret(0.);
        const double x = p[0]*(2* P.get<double>("Exp.RadInlet") -p[0]) / (P.get<double>("Exp.RadInlet")*P.get<double>("Exp.RadInlet")),
                     z = p[2]*(2*P.get<double>("Exp.RadInlet")-p[2]) / (P.get<double>("Exp.RadInlet")*P.get<double>("Exp.RadInlet"));

        ret[1]= x * z * P.get<double>("Exp.InflowVel") * (1-P.get<double>("Exp.InflowAmpl")*std::cos(2*M_PI*P.get<double>("Exp.InflowFreq")*t));

        return ret;
    }


    ///microchannel (eindhoven)
    DROPS::SVectorCL<3> InflowChannel( const DROPS::Point3DCL& p, double t)
    {
        DROPS::SVectorCL<3> ret(0.);
        const double y = p[1]*(2*25e-6-p[1]) / (25e-6*25e-6),
                     z = p[2]*(2*50e-6-p[2]) / (50e-6*50e-6);

        ret[0]= y * z * P.get<double>("Exp.InflowVel") * (1-P.get<double>("Exp.InflowAmpl")*std::cos(2*M_PI*P.get<double>("Exp.InflowFreq")*t));

        return ret;
    }

    ///mzelle_ns_adap.cpp + mzelle_instat.cpp
    DROPS::SVectorCL<3> InflowCell( const DROPS::Point3DCL& p, double)
    {
        DROPS::SVectorCL<3> ret(0.);
        const double s2= P.get<double>("Exp.RadInlet")*P.get<double>("Exp.RadInlet"),
                     r2= p.norm_sq() - p[P.get<int>("Exp.FlowDir")]*p[P.get<int>("Exp.FlowDir")];
        ret[P.get<int>("Exp.FlowDir")]= -(r2-s2)/s2*P.get<double>("Exp.InflowVel");

        return ret;
    }


    //========================================================================
    //                       Functions for brick_transp.cpp
    //========================================================================

    DROPS::SVectorCL<3> InflowBrickTransp (const DROPS::Point3DCL& p, double)
    {
        DROPS::SVectorCL<3> ret(0.);
        const double x = p[0]*(2*P.get<double>("Exp.RadInlet")-p[0]) / (P.get<double>("Exp.RadInlet")*P.get<double>("Exp.RadInlet")),
                     z = p[2]*(2*P.get<double>("Exp.RadInlet")-p[2]) / (P.get<double>("Exp.RadInlet")*P.get<double>("Exp.RadInlet"));
        ret[1]= x * z * P.get<double>("Exp.InflowVel");
        return ret;
    }

    //========================================================================
    //            Registration of functions in the func-container
    //========================================================================
    static DROPS::RegisterVectorFunction regvelbrick("InflowBrick", InflowBrick);
    static DROPS::RegisterVectorFunction regvelcell("InflowCell", InflowCell);
    static DROPS::RegisterVectorFunction regvelchannel("InflowChannel", InflowChannel);
    static DROPS::RegisterVectorFunction regvelbricktransp("InflowBrickTransp", InflowBrickTransp);
}


//========================================================================
//          Functions for twophasedrops-executable (Volume Force)
//========================================================================
namespace tpd_volforce{
    /// \name inflow condition
    template<int D>
    DROPS::SVectorCL<3> PeriodicDropletPressure( const DROPS::Point3DCL& , double )
    {
        DROPS::SVectorCL<3> ret(0.);
        ret[D] = -P.get<DROPS::Point3DCL>("Exp.Gravity")[D];

        static bool first = true;
        static DROPS::Point3DCL dx;
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

        static double voldrop = 4./3.*M_PI* P.get<DROPS::Point3DCL>("Exp.RadDrop")[0]*P.get<DROPS::Point3DCL>("Exp.RadDrop")[1]*P.get<DROPS::Point3DCL>("Exp.RadDrop")[2] ;
        static double brickvol = dx[0]* dx[1]* dx[2];
        static double volforce = P.get<double>("Mat.DensFluid") * brickvol - (P.get<double>("Mat.DensFluid") - P.get<double>("Mat.DensDrop")) * voldrop;
        ret[D] *= volforce/brickvol;
        return ret;
    }

    //========================================================================
    //            Registration of functions in the func-container
    //========================================================================
    static DROPS::RegisterVectorFunction regvelppx("PeriodicDropletPressurex", PeriodicDropletPressure<0>);
    static DROPS::RegisterVectorFunction regvelppy("PeriodicDropletPressurey", PeriodicDropletPressure<1>);
    static DROPS::RegisterVectorFunction regvelppz("PeriodicDropletPressurez", PeriodicDropletPressure<2>);
}


//========================================================================
//                       Functions for LevelSet Distance
//========================================================================
namespace levelsetdistance{

    ///mzelle_ns_adap.cpp + mzelle_instat.cpp
    double CubeDistance( const DROPS::Point3DCL& p, double)
    {
        double maxd = - P.get<DROPS::Point3DCL>("Exp.RadDrop")[0] - P.get<DROPS::Point3DCL>("Exp.RadDrop")[1]- P.get<DROPS::Point3DCL>("Exp.RadDrop")[2];
        for (int i=0;i<3;i++){
          double x = std::abs(p[i] - P.get<DROPS::Point3DCL>("Exp.PosDrop")[i]) - P.get<DROPS::Point3DCL>("Exp.RadDrop")[i];
          if (x>maxd) maxd=x;
        }
        return maxd;
    }


    ///mzelle_ns_adap.cpp + mzelle_instat.cpp
    template<int i>
    double PeriodicEllipsoidDistance( const DROPS::Point3DCL& p, double)
    {

        static bool first = true;
        static DROPS::Point3DCL dx;
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

        DROPS::Point3DCL dp;
        dp[i] = dx[i];
        DROPS::Point3DCL ExpRadDrop = P.get<DROPS::Point3DCL>("Exp.RadDrop");
        DROPS::Point3DCL ExpPosDrop = P.get<DROPS::Point3DCL>("Exp.PosDrop");
        DROPS::Point3DCL d= p - ExpPosDrop;
        DROPS::Point3DCL d1= p + dp - ExpPosDrop;
        DROPS::Point3DCL d2= p - dp - ExpPosDrop;
        const double avgRad= cbrt(ExpRadDrop[0]*ExpRadDrop[1]*ExpRadDrop[2]);
        d/= ExpRadDrop;
        d1/= ExpRadDrop;
        d2/= ExpRadDrop;
        double dd = std::min(std::min(d.norm(),d1.norm()),d2.norm());
        return std::abs( avgRad)*dd - avgRad;
    }

    template<int i>
    double planedistance( const DROPS::Point3DCL& p, double)
    {
        double x=p[i]-P.get<DROPS::Point3DCL>("Exp.PosDrop")[i];
        return x;
    }

    static DROPS::RegisterScalarFunction regsca_pdistx("planedistancex", planedistance<0>);
    static DROPS::RegisterScalarFunction regsca_pdisty("planedistancey", planedistance<1>);
    static DROPS::RegisterScalarFunction regsca_pdistz("planedistancez", planedistance<2>);
    static DROPS::RegisterScalarFunction regscacube("CubeDistance", CubeDistance);
    static DROPS::RegisterScalarFunction regscaperellx("perEllipsoidx", PeriodicEllipsoidDistance<0>);
    static DROPS::RegisterScalarFunction regscaperelly("perEllipsoidy", PeriodicEllipsoidDistance<1>);
    static DROPS::RegisterScalarFunction regscaperellz("perEllipsoidz", PeriodicEllipsoidDistance<2>);
}

//========================================================================
//                       Functions for transport
//========================================================================
namespace transpfunctions{
    double tInitialcneg (const DROPS::Point3DCL& , double)
    {
        return P.get<double>("Transp.IniCNeg");
    }

    double tInitialcpos (const DROPS::Point3DCL& , double)
    {
        return P.get<double>("Transp.IniCPos");
    }
    static DROPS::RegisterScalarFunction regscainineg("Initialcneg", tInitialcneg);
    static DROPS::RegisterScalarFunction regscainipos("Initialcpos", tInitialcpos);
}

//========================================================================
//                       Functions for surfactants
//========================================================================
namespace surffunctions{
    /// \name Initial data and rhs for surfactant transport
    //@{
    const double a( -13./8.*std::sqrt( 35./M_PI));
    double surf_rhs (const DROPS::Point3DCL& p, double)
    {
        return a*(3.*p[0]*p[0]*p[1] - p[1]*p[1]*p[1]);
    }
    double surf_sol (const DROPS::Point3DCL& p, double)
    {
        return 1. + std::sin( atan2( p[0] - P.get<DROPS::Point3DCL>("Exp.PosDrop")[0], p[2] - P.get<DROPS::Point3DCL>("Exp.PosDrop")[2]));
    }
    //@}
    static DROPS::RegisterScalarFunction regscasurfrhs("surf_rhs", surf_rhs);
    static DROPS::RegisterScalarFunction regscasurfsol("surf_sol", surf_sol);
}

//TODO: unification with filmCoeff stoff
//========================================================================
//                        Functions for matching function
//========================================================================
namespace filmperiodic{
    template<int A, int B>
    bool periodic_2sides( const DROPS::Point3DCL& p, const DROPS::Point3DCL& q)
    {

        static bool first = true;
        static DROPS::Point3DCL dx;
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

        const DROPS::Point3DCL d= fabs(p-q),
                               L= fabs(dx);

        const int D = 3 - A - B;
        return (d[B] + d[D] < 1e-12 && std::abs( d[A] - L[A]) < 1e-12)  // dB=dD=0 and dA=LA
          ||   (d[A] + d[D] < 1e-12 && std::abs( d[B] - L[B]) < 1e-12)  // dA=dD=0 and dB=LB
          ||   (d[D] < 1e-12 && std::abs( d[A] - L[A]) < 1e-12 && std::abs( d[B] - L[B]) < 1e-12);  // dD=0 and dA=LA and dB=LB
    }

    template<int A>
    bool periodic_1side( const DROPS::Point3DCL& p, const DROPS::Point3DCL& q)
    {

        static bool first = true;
        static DROPS::Point3DCL dx;
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
        const int B = (A+1)%3;
        const int D = (B+1)%3;
        const DROPS::Point3DCL d= fabs(p-q), L= fabs(dx);

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
double TaylorFlowDistance( const DROPS::Point3DCL& p, double)
{
    static const double taylor_len = P.get<double>("Taylor.Length");
    static const double taylor_width = P.get<double>("Taylor.Width");
    static const double rel_filmthickness = P.get<double>("Taylor.RelFilmThickness");
    static const DROPS::Point3DCL taylor_bubble_center = P.get<DROPS::Point3DCL>("Taylor.Center");
    //const double taylor_top_z = taylor_bubble_center[2] + 0.5 * taylor_bubble_len;
    //const double taylor_bottom_z = taylor_bubble_center[2] - 0.5 * taylor_bubble_len;
    const double taylor_radius = (0.5-rel_filmthickness) * taylor_width;
    const double taylor_top_z = taylor_bubble_center[2] + 0.5 * (taylor_len-2*taylor_radius);
    const double taylor_bottom_z = taylor_bubble_center[2] - 0.5 * (taylor_len-2*taylor_radius);
    DROPS::Point3DCL diff = p;
    diff[0] -= taylor_bubble_center[0];
    diff[1] -= taylor_bubble_center[1];
    if (p[2] > taylor_top_z){
        diff[2] -= taylor_top_z;
        return diff.norm() - taylor_radius;
    }
    if (p[2] < taylor_bottom_z){
        diff[2] -= taylor_bottom_z;
        return diff.norm() - taylor_radius;
    }
    diff[2] = 0.0;
    return diff.norm() - taylor_radius;
}

static DROPS::RegisterScalarFunction regscataylor("TaylorFlowDistance", TaylorFlowDistance);


namespace slipBnd{
    /// \name inflow condition

    DROPS::SVectorCL<3> InflowPoiseuille (const DROPS::Point3DCL& p, double)
    {
        DROPS::SVectorCL<3> ret(0.);
        static bool first = true;
        static DROPS::Point3DCL dx;
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

        static double half_z= dx[2];
        ret[0] = 0.0125 - 5.*(p[2]-half_z)*(p[2]-half_z);
        return ret;
    }
	
    DROPS::SVectorCL<3> InflowSine (const DROPS::Point3DCL& p, double)
    {
        DROPS::SVectorCL<3> ret(0.);
        static bool first = true;
        static DROPS::Point3DCL dx;
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

        static double half_z= dx[2];
        ret[0] = sin( p[2]/half_z * M_PI_2);
        return ret;
    }
	
	DROPS::SVectorCL<3> SineF (const DROPS::Point3DCL& p, double)
    {
        DROPS::SVectorCL<3> ret(0.);
        static bool first = true;
        static DROPS::Point3DCL dx;
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

        static double half_z= dx[2];
        ret[0] = (1.0/half_z * M_PI_2)* (1.0/half_z * M_PI_2) * sin( p[2]/half_z * M_PI_2);
        return ret;
    }
	
    DROPS::SVectorCL<3> WallVel( const DROPS::Point3DCL& p, double t)
    {
        DROPS::SVectorCL<3> V_w(0.);
        static bool first = true;
		static double nu, l_s;
		if(first){
			nu=P.get<double>("Mat.ViscFluid")/P.get<double>("Mat.DensFluid");
			l_s=1./(P.get<double>("SpeBnd.beta2"));
			first = false;
		}
		double		    F=std::sin(2*nu*t);
		V_w[0] = std::sin(p[0]/l_s) * F * ( std::cos(M_PI_2/l_s) - std::sin(M_PI_2/l_s) );
        return V_w;
    }
    DROPS::SVectorCL<3> TopVel( const DROPS::Point3DCL&, double)
    {
        DROPS::SVectorCL<3> ret(0.);
		ret[0] = 0.4;
		//ret[0] = 2.4 * p[0] * (1.0-p[0]);
        return ret;
    }
    //========================================================================
    //            Registration of functions in the func-container
    //========================================================================
    static DROPS::RegisterVectorFunction regvelwall("WallVel", WallVel);
	static DROPS::RegisterVectorFunction regveltop("TopVel", TopVel);
    static DROPS::RegisterVectorFunction regvelpoiseuille("InflowPoiseuille", InflowPoiseuille);
    static DROPS::RegisterVectorFunction regvelsine("InflowSine", InflowSine);
	static DROPS::RegisterVectorFunction regvelsineF("SineF", SineF);
}
namespace InstatSlip{

	DROPS::SVectorCL<3> Velocity( const DROPS::Point3DCL& p,double t)
	{
		DROPS::SVectorCL<3> v(0.);
        static bool first = true;
		static double nu, l_s;
		if(first){
			nu=P.get<double>("Mat.ViscFluid")/P.get<double>("Mat.DensFluid");
			l_s=1./(P.get<double>("SpeBnd.beta2"));
			first = false;
		}
		double		        F=std::sin(2*nu*t);

		v[0]=  std::sin(p[0]/l_s)*std::cos(p[1]/l_s)*F;
		v[1]= -std::cos(p[0]/l_s)*std::sin(p[1]/l_s)*F;
		v[2]= 0;

		return v;
	}

	DROPS::SVectorCL<3> PressureGr(const DROPS::Point3DCL& p, double t)
	{
		DROPS::SVectorCL<3> delp(0.);
        static bool first = true;
		static double nu, l_s;
		if(first){
			nu=P.get<double>("Mat.ViscFluid")/P.get<double>("Mat.DensFluid");
			l_s=1./(P.get<double>("SpeBnd.beta2"));
			first = false;
		}
		double F=std::sin(2*nu*t);
		delp[0]=-1./(2.*l_s)*std::sin(2*p[0]/l_s)*F*F;
		delp[1]=-1./(2.*l_s)*std::sin(2*p[1]/l_s)*F*F;
		delp[2]= 0;

		return delp;

	}
	
	double Pressure (const DROPS::Point3DCL& p, double t)
	{
		double norm2=0.;
		double ret=0;		
		static bool first = true;
		static DROPS::Point3DCL origin;
		static double radius, surftension;
		if(first){
			DROPS::Point3DCL org0 = P.get<DROPS::Point3DCL>("Exp.PosDrop");
			DROPS::Point3DCL VelR = P.get<DROPS::Point3DCL>("Exp.RadDrop");
			double angle=P.get<double>("SpeBnd.contactangle")*M_PI/180;
			//radius=VelR[0];
			radius=VelR[0]*std::pow(2/(2+std::cos(angle))/std::pow(1-std::cos(angle),2),1.0/3);
			//assume the initial droplet is semi-spherical;
			org0[1]=org0[1]-radius*std::cos(angle);//The drop is located in the plain normal to the y direction
			origin=org0;
			surftension = P.get<double>("SurfTens.SurfTension");
			first = false;
		}

		for (int i=0; i< 3; i++)
		  norm2 += (p[i]-origin[i]) * (p[i]-origin[i]);
         ret = (std::sqrt(norm2) > radius) ? 0: 2.*surftension/radius; 

		if( t>3. && t<5 )
			ret = 0;
		else if(t >5.)
			ret = 2.*surftension/radius;
		return ret;
	}
	

	DROPS::SVectorCL<3> VolForce( const DROPS::Point3DCL& p, double t)
    {
        DROPS::SVectorCL<3> f(0.);
        static bool first = true;
		 static double nu, l_s;
		 if(first){
			 nu=P.get<double>("Mat.ViscFluid")/P.get<double>("Mat.DensFluid");
			 l_s=1./(P.get<double>("SpeBnd.beta2"));
			 first = false;
		 }
		 f[0] =  2.* nu * std::sin(p[0]/l_s) * std::cos(p[1]/l_s) * ( std::cos(2*nu*t) + 1./(l_s*l_s) * std::sin(2*nu*t) ) ;
		 f[1] = -2.* nu * std::cos(p[0]/l_s) * std::sin(p[1]/l_s) * ( std::cos(2*nu*t) + 1./(l_s*l_s) * std::sin(2*nu*t) ) ;
		 f[2] = 0;

        return f+PressureGr(p,t);
    }
    static DROPS::RegisterVectorFunction regvelVel("TestSlipVel", Velocity);
    static DROPS::RegisterVectorFunction regvelf("TestSlipF", VolForce);
	static DROPS::RegisterVectorFunction regvelgpr("TestSlipPrGrad",PressureGr);
	static DROPS::RegisterScalarFunction regscatestpr("TestSlipPressure",Pressure);

}

namespace AnalPressure{

	double Pressure (const DROPS::Point3DCL& , double t)
	{
		double ret=0;		
		static bool first = true;
		static DROPS::Point3DCL origin, VelR;
		static double radius, surftension;
		if(first){
			origin = P.get<DROPS::Point3DCL>("Exp.PosDrop");
			VelR = P.get<DROPS::Point3DCL>("Exp.RadDrop");
			radius = P.get<double>("Exp.EquiRadius");
			surftension = P.get<double>("SurfTens.SurfTension");
			first = false;
		}
		
		if( t>3. && t<5 )
			ret = 0;
		else if(t >5.)
			ret = 2.*surftension/radius;
		return ret;
	}
    static DROPS::RegisterScalarFunction regscatestpr("AnalPressure",Pressure);
}


namespace StatSlip{

	DROPS::SVectorCL<3> Velocity( const DROPS::Point3DCL& p,double)
	{
		DROPS::SVectorCL<3> v(0.);
		double l_s=2.;
		v[0]=  std::sin(p[0]/l_s)*std::cos(p[1]/l_s);
		v[1]= -std::cos(p[0]/l_s)*std::sin(p[1]/l_s);
		v[2]= 0;
		return v;
	}

	DROPS::SVectorCL<3> PressureGr(const DROPS::Point3DCL& p, double)
	{
		DROPS::SVectorCL<3> delp(0.);
		double l_s=2.;
		delp[0]=-1./(2.*l_s)*std::sin(2*p[0]/l_s);
		delp[1]=-1./(2.*l_s)*std::sin(2*p[1]/l_s);
		delp[2]= 0;

		return delp;

	}

	DROPS::SVectorCL<3> VolForce( const DROPS::Point3DCL& p, double)
    {
        DROPS::SVectorCL<3> f(0.);
		double nu=1., l_s=2.;
		f[0] =  2.* nu * std::sin(p[0]/l_s) * std::cos(p[1]/l_s) *  1./(l_s*l_s) ;
		f[1] = -2.* nu * std::cos(p[0]/l_s) * std::sin(p[1]/l_s) *  1./(l_s*l_s) ;
		f[2] = 0;

        return f+PressureGr(p,0);
    }
    static DROPS::RegisterVectorFunction regvelVel("StatSlipVel", Velocity);
    static DROPS::RegisterVectorFunction regvelf("StatSlipF", VolForce);
	static DROPS::RegisterVectorFunction regvelgpr("StatSlipPrGrad",PressureGr);

}


namespace StatSlip2{

	DROPS::SVectorCL<3> Velocity( const DROPS::Point3DCL& p,double)
	{
		DROPS::SVectorCL<3> v(0.);
		v[0]= p[0]*p[1]-2.*p[0];
		v[1]= -0.5*p[1]*p[1]+2.*p[1]-1.5;
		v[2]= 0;
		return v;
	}

	DROPS::SVectorCL<3> PressureGr(const DROPS::Point3DCL& p, double)
	{
		DROPS::SVectorCL<3> delp(0.);
		delp[0]= std::sin(p[0]);
		delp[1]= std::cos(p[1]);
		delp[2]= 0;

		return delp;

	}

	DROPS::SVectorCL<3> VolForce( const DROPS::Point3DCL& p, double)
    {
        DROPS::SVectorCL<3> f(0.);
		f[0] = 0;
		f[1] = 1;
		f[2] = 0;

        return f+PressureGr(p,0);
    }
    static DROPS::RegisterVectorFunction regvelVel2("StatSlip2Vel", Velocity);
    static DROPS::RegisterVectorFunction regvelf2("StatSlip2F", VolForce);
    static DROPS::RegisterVectorFunction regvelgpr2("StatSlip2PrGrad",PressureGr);

}

namespace StatSlip3{

	DROPS::SVectorCL<3> Velocity( const DROPS::Point3DCL& p,double)
	{
		DROPS::SVectorCL<3> v(0.);
		v[0]= p[0]*p[0]*p[1]-2.*p[0]*p[0];  //x^2y-2x^2
		v[1]= -p[0]*p[1]*p[1]+3.*p[0]*p[1]-2.*p[0];
		v[2]=  p[0]*p[2];
		return v;
	}
	
	double Pressure(const DROPS::Point3DCL& p, double)
	{
		double pr = - std::cos(p[0]) + std::sin(p[1]);

		return pr;
	}

	DROPS::SVectorCL<3> PressureGr(const DROPS::Point3DCL& p, double)
	{
		DROPS::SVectorCL<3> delp(0.);
		delp[0]= std::sin(p[0]);
		delp[1]= std::cos(p[1]);
		delp[2]= 0;

		return delp;

	}

	DROPS::SVectorCL<3> VolForce( const DROPS::Point3DCL& p, double)
    {
        DROPS::SVectorCL<3> f(0.);
		f[0] = -2.*p[1]+4.;
		f[1] = 2.*p[0];
		f[2] = 0;

        return f+PressureGr(p,0);
    }
	
	DROPS::SVectorCL<3> flux( const DROPS::Point3DCL& p,double)
	{
		DROPS::SVectorCL<3> v(0.);
		v[0]= 0;
		v[1]= 0;
		v[2]= p[0];
		return v;
	}
	
	DROPS::SVectorCL<3> Wall3( const DROPS::Point3DCL& p,double)
	{
		DROPS::SVectorCL<3> v(0.);
		v[0]= 0;
		v[1]= 0;
		v[2]= p[0]*p[2];
		return v;
	}
    static DROPS::RegisterVectorFunction regvelVel3("StatSlip3Vel", Velocity);
    static DROPS::RegisterVectorFunction regvelWall3("Wall3", Wall3);
    static DROPS::RegisterVectorFunction regvelf3("StatSlip3F", VolForce);
	static DROPS::RegisterVectorFunction regvelgpr3("StatSlip3PrGrad",PressureGr);
	static DROPS::RegisterVectorFunction regvelflux("StatSlip3flux",flux);
	static DROPS::RegisterScalarFunction regscapr3("StatSlip3Pr",Pressure);
}

namespace contactangle{
	
	double ConstantAngle(const DROPS::Point3DCL&,double)
	{
       double angle = P.get<double>("SpeBnd.contactangle");
		return angle/180.0*M_PI;
	}

	double PeriodicAngle(const DROPS::Point3DCL& pt,double)
	{
		DROPS::Point3DCL pt1(pt-P.get<DROPS::Point3DCL>("SpeBnd.posDrop"));
		double r=std::sqrt(pt1[0]*pt1[0]+pt1[2]*pt1[2]);
		double theta=r<0.001? 0: (pt1[2]>0?std::acos(pt1[0]/r): 2*M_PI- std::acos(pt1[0]/r));
		return (P.get<double>("SpeBnd.contactangle"))*(1+0.5*std::sin(30*theta))/180.0*M_PI;
	}

	double PatternAngle(const DROPS::Point3DCL& pt,double)
	{
		DROPS::Point3DCL pt1(pt-P.get<DROPS::Point3DCL>("SpeBnd.posDrop"));
		double r=std::sqrt(pt1[0]*pt1[0]+pt1[2]*pt1[2]);
		double theta= int(r/P.get<DROPS::Point3DCL>("Exp.RadDrop")[0]/5.0)%2==0?1:-1;
		return P.get<double>("SpeBnd.contactangle")/180.0*M_PI*(1+0.5*theta);
	}
	double PatternAngle1(const DROPS::Point3DCL& pt,double)
		{
			DROPS::Point3DCL pt1(pt-P.get<DROPS::Point3DCL>("SpeBnd.posDrop"));
			double angl=P.get<double>("SpeBnd.contactangle")/180.0*M_PI;
			return pt1[0]*pt1[2]<0?angl:M_PI-angl;
		}
	double PatternAngle2(const DROPS::Point3DCL& pt,double)
		{
			DROPS::Point3DCL pt1(pt-P.get<DROPS::Point3DCL>("SpeBnd.posDrop"));
			double r=std::sqrt(pt1[0]*pt1[0]+pt1[2]*pt1[2]);
			double theta=r<0.001? 0: (pt1[2]>0? std::acos(pt1[0]/r)/M_PI*180 : 360 - std::acos(pt1[0]/r)/M_PI*180);
			double angl=P.get<double>("SpeBnd.contactangle")/180.0*M_PI;
			return (int(theta)/60)%2==0?angl:M_PI-angl;
		}
	DROPS::Point3DCL OutNormalBottomPlane(const DROPS::Point3DCL&,double)
	{
		DROPS::Point3DCL outnormal(0.0);
	//	outnormal[0]=0;
		outnormal[1]=-1.0;
	//	outnormal[2]=0;
		return outnormal;
	}

	DROPS::Point3DCL OutNormalBrick(const DROPS::Point3DCL& pt,double)
	{
		if(P.get<int>("DomainCond.GeomType")!=1)
			 throw DROPS::DROPSErrCL("Error: compute out normal of brick, please use other functions");

       double dx, dy, dz;

        std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
        size_t idx;
        while ((idx= mesh.find_first_of( delim)) != std::string::npos )
               mesh[idx]= ' ';
        std::istringstream brick_info( mesh);
        brick_info >> dx >> dy >> dz;

        if (!brick_info)
            throw DROPS::DROPSErrCL("error while reading geometry information: " + mesh);

		double EPS=1e-10;
		DROPS::Point3DCL outnormal(0.0);
		if(std::fabs(pt[0])<EPS)
		{	outnormal[0]=-1.0; }
		else if(std::fabs(pt[0]-dx)<EPS)
		{	outnormal[0]=1.0; }
		else if(std::fabs(pt[1])<EPS)
		{	outnormal[1]=-1.0; }
		else if(std::fabs(pt[1]-dy)<EPS)
		{	outnormal[1]=1.0; }
		else if(std::fabs(pt[2])<EPS)
		{	outnormal[2]=-1.0; }
		else if(std::fabs(pt[2]-dz)<EPS)
		{	outnormal[2]=1.0; }
		else
			 throw DROPS::DROPSErrCL("error while computing outnormal");
		return outnormal;
	}
	
	static DROPS::RegisterScalarFunction regconstangle("ConstantAngle", ConstantAngle);
	static DROPS::RegisterScalarFunction regperangle("PeriodicAngle", PeriodicAngle);
	static DROPS::RegisterScalarFunction regpatangle("PatternAngle", PatternAngle);
	static DROPS::RegisterScalarFunction regpatangle1("PatternAngle1", PatternAngle1);
	static DROPS::RegisterScalarFunction regpatangle2("PatternAngle2", PatternAngle2);
	static DROPS::RegisterVectorFunction regunitbottomoutnomal("OutNormalBottomPlane", OutNormalBottomPlane);
	static DROPS::RegisterVectorFunction regunitcubicoutnomal("OutNormalBrick", OutNormalBrick);
}

namespace curvebndDomain{

	DROPS::Point3DCL OutNormalSphere(const DROPS::Point3DCL& pt,double)
	{
		DROPS::Point3DCL outnormal(0.0);
		outnormal=pt/pt.norm();
		return outnormal;
	}

	static DROPS::RegisterVectorFunction regunitoutnomalsphere("OutNormalSphere", OutNormalSphere);
    DROPS::Point3DCL OutNormalCylinder(const DROPS::Point3DCL& pt,double)
	{
		DROPS::Point3DCL outnormal(0.0);
		outnormal=pt/std::sqrt(pt[1]*pt[1]+pt[2]*pt[2]);
       outnormal[0]=0;
		return outnormal;
	}
	static DROPS::RegisterVectorFunction regunitoutnomalcylinder("OutNormalCylinder", OutNormalCylinder);

}
namespace CouetteFlow{

	DROPS::Point3DCL UpwallVel(const DROPS::Point3DCL& ,double)
	{ 
        double WallVel = P.get<double>("Exp.WallVelocity");
        DROPS::Point3DCL vel(0.0);
        vel[0]=WallVel; 
        return vel;
	}
	DROPS::Point3DCL DownwallVel(const DROPS::Point3DCL& ,double)
	{ 
        double WallVel = P.get<double>("Exp.WallVelocity");
        DROPS::Point3DCL vel(0.0);
        vel[0]=-WallVel; 
        return vel;
	}
	DROPS::Point3DCL LeftwallVel(const DROPS::Point3DCL& pt,double)
	{

		 double dx, dy, dz;
		 const std::string meshfile_name=P.get<std::string>("DomainCond.MeshFile");
		 std::string mesh( meshfile_name), delim("x@");
		 size_t idx;
		 while ((idx= mesh.find_first_of( delim)) != std::string::npos )
		       mesh[idx]= ' ';
		 std::istringstream brick_info( mesh);
		 brick_info >> dx >> dy >> dz;
		 DROPS::Point3DCL vel(0.0);
		 vel[0]=P.get<double>("Exp.WallVelocity")*((pt[2]-dz/2)*2/dz)*(dz/(dz+2*P.get<double>("Mat.ViscDrop")/P.get<double>("SpeBnd.beta1")));
		 return vel;
	}
	DROPS::Point3DCL RightwallVel(const DROPS::Point3DCL& pt,double)
		{

			 double dx, dy, dz;
			 const std::string meshfile_name=P.get<std::string>("DomainCond.MeshFile");
			 std::string mesh( meshfile_name), delim("x@");
			 size_t idx;
			 while ((idx= mesh.find_first_of( delim)) != std::string::npos )
			       mesh[idx]= ' ';
			 std::istringstream brick_info( mesh);
			 brick_info >> dx >> dy >> dz;
			 DROPS::Point3DCL vel(0.0);
			 vel[0]=P.get<double>("Exp.WallVelocity")*((pt[2]-dz/2)*2/dz)*(dz/(dz+2*P.get<double>("Mat.ViscFluid")/P.get<double>("SpeBnd.beta2")));
			 return vel;
		}
	static DROPS::RegisterVectorFunction regunitupwallvel("UpwallVel", UpwallVel);
	static DROPS::RegisterVectorFunction regunitdownwallvel("DownwallVel", DownwallVel);
	static DROPS::RegisterVectorFunction regunitleftwallvel("LeftwallVel", LeftwallVel);
	static DROPS::RegisterVectorFunction regunitrightwallvel("RightwallVel", RightwallVel);
}

