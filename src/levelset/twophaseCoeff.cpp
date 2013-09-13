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

        static double half_z= dx[2]/2.0;
        ret[0] = 0.0125 - 5.*(p[2]-half_z)*(p[2]-half_z);
        return ret;
    }
    DROPS::SVectorCL<3> WallVel( const DROPS::Point3DCL& p, double t)
    {
        DROPS::SVectorCL<3> V_w(0.);
        static bool first = true;
		static double nu, l_s;
		if(first){
			nu=P.get<double>("Mat.ViscFluid")/P.get<double>("Mat.DensFluid");
			l_s=1./(P.get<double>("SpeBnd.SlipLength2"));
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
}
namespace InstatSlip{

	DROPS::SVectorCL<3> Velocity( const DROPS::Point3DCL& p,double t)
	{
		DROPS::SVectorCL<3> v(0.);
        static bool first = true;
		static double nu, l_s;
		if(first){
			nu=P.get<double>("Mat.ViscFluid")/P.get<double>("Mat.DensFluid");
			l_s=1./(P.get<double>("SpeBnd.SlipLength2"));
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
			l_s=1./(P.get<double>("SpeBnd.SlipLength2"));
			first = false;
		}
		double F=std::sin(2*nu*t);
		delp[0]=-1./(2.*l_s)*std::sin(2*p[0]/l_s)*F*F;
		delp[1]=-1./(2.*l_s)*std::sin(2*p[1]/l_s)*F*F;
		delp[2]= 0;

		return delp;

	}

	DROPS::SVectorCL<3> VolForce( const DROPS::Point3DCL& p, double t)
    {
        DROPS::SVectorCL<3> f(0.);
        static bool first = true;
		static double nu, l_s;
		if(first){
			nu=P.get<double>("Mat.ViscFluid")/P.get<double>("Mat.DensFluid");
			l_s=1./(P.get<double>("SpeBnd.SlipLength2"));
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

}

namespace StatSlip{

	DROPS::SVectorCL<3> Velocity( const DROPS::Point3DCL& p,double)
	{
		DROPS::SVectorCL<3> v(0.);
        static bool first = true;
		static double l_s;
		if(first){
			l_s=1./(P.get<double>("SpeBnd.SlipLength2"));
			first = false;
		}
		v[0]=  std::sin(p[0]/l_s)*std::cos(p[1]/l_s);
		v[1]= -std::cos(p[0]/l_s)*std::sin(p[1]/l_s);
		v[2]= 0;
		return v;
	}

	DROPS::SVectorCL<3> PressureGr(const DROPS::Point3DCL& p, double)
	{
		DROPS::SVectorCL<3> delp(0.);
        static bool first = true;
		static double l_s;
		if(first){
			l_s=1./(P.get<double>("SpeBnd.SlipLength2"));
			first = false;
		}
		delp[0]=-1./(2.*l_s)*std::sin(2*p[0]/l_s);
		delp[1]=-1./(2.*l_s)*std::sin(2*p[1]/l_s);
		delp[2]= 0;

		return delp;

	}

	DROPS::SVectorCL<3> VolForce( const DROPS::Point3DCL& p, double)
    {
        DROPS::SVectorCL<3> f(0.);
        static bool first = true;
		static double nu, l_s;
		if(first){
			nu=P.get<double>("Mat.ViscFluid")/P.get<double>("Mat.DensFluid");
			l_s=1./(P.get<double>("SpeBnd.SlipLength2"));
			first = false;
		}
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

}