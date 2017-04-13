/// \file osmosisCoeff.cpp
/// \brief boundary and source functions for the osmosis Problem
/// \author LNM RWTH Aachen: Thorolf Schulte, Christoph Lehrenfeld

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
 * Copyright 2012 LNM/SC RWTH Aachen, Germany
*/

#include "misc/funcmap.h"
#include "misc/params.h"
#include "misc/utils.h"

double Initialcneg (const DROPS::Point3DCL& p, double )
{
	/*extern DROPS::ParamCL P;
	double t = 0.0;
	static double r0 = 1.0;
	static double dn = P.get<double>("Osmosis.Diffusivity");
	static double vn = P.get<double>("Osmosis.GrowVelocity");
	static DROPS::Point3DCL drop = P.get<DROPS::Point3DCL>("Levelset.PosDrop");
	double r2 = (p[0]-drop[0])*(p[0]-drop[0])
        		+ (p[1]-drop[1])*(p[1]-drop[1])
        		+ (p[2]-drop[2])*(p[2]-drop[2]);
	//EXPANSION / CONTRACTION
	return vn/dn * r2 - 2*(r0+vn*t) - vn/dn * (r0+vn*t)*(r0+vn*t);
	//STATIC
	//return cos(r2 - 1);*/
	return p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
}

double Initialcpos (const DROPS::Point3DCL& p, double )
{
	/*extern DROPS::ParamCL P;
	double t = 0.0;
	static double r0 = 1.0;
	static double dn = P.get<double>("Osmosis.Diffusivity");
	static double vn = P.get<double>("Osmosis.GrowVelocity");
	static DROPS::Point3DCL drop = P.get<DROPS::Point3DCL>("Levelset.PosDrop");
	double r2 = (p[0]-drop[0])*(p[0]-drop[0])
        		+ (p[1]-drop[1])*(p[1]-drop[1])
        		+ (p[2]-drop[2])*(p[2]-drop[2]);
	//EXPANSION / CONTRACTION
	return vn/dn * r2 - 2*(r0+vn*t) - vn/dn * (r0+vn*t)*(r0+vn*t);
	//STATIC
	//return cos(r2 - 1);*/
	return p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
}

double Reaction (const DROPS::Point3DCL& , double )
{
  return 0.0;
}

double linearConc (const DROPS::Point3DCL& p , double )
{
	return p[0];
}

double Sinx (const DROPS::Point3DCL& p , double )
{
	extern DROPS::ParamCL P;
	static double scal = P.get("Levelset.SinxScale", 1.0);
	static double angl = P.get("Levelset.SinxAngle", 1.0);
    const double a = sin(scal*(p[0]+angl*p[1]));
    return a*a;
}

double Heaviside (const DROPS::Point3DCL& p , double )
{
	return 1.0 - p[0]/4.0;
}

double TransverseLinear (const DROPS::Point3DCL& p , double )
{
    extern DROPS::ParamCL P;
	static DROPS::Point3DCL center = P.get<DROPS::Point3DCL>("Levelset.PosDrop");
    const double axval = (p[0]-center[0])+2.0*(p[1]-center[1]);
	if (0.125*(4.0+axval) > 0.5)
		return 0.1;
	else
		return 1.;
    // return 0.5*(1+std::tanh(5.0*(axval-0.5)));
}

double ZeroFct (const DROPS::Point3DCL& ,  double )
{  
    return 0.;
}

double Torus( const DROPS::Point3DCL& pin, double)
{
    extern DROPS::ParamCL P;
    const double R_ = 1.2;
    const double r_ = 0.4;
	static DROPS::Point3DCL center = P.get<DROPS::Point3DCL>("Levelset.PosDrop");
    DROPS::Point3DCL p = pin - center;
    return std::sqrt( p[2]*p[2] + std::pow( std::sqrt( p[0]*p[0] + p[1]*p[1]) - R_, 2) ) - r_;
}


DROPS::Point3DCL ZeroFct_vec (const DROPS::Point3DCL& ,  double )
{  
    return DROPS::Point3DCL(0.);
}

double Rhs (const DROPS::Point3DCL& ,  double )
{
	return 0.;
}

double Dirichlet (const DROPS::Point3DCL& p, double )
{
/*  static double x0 = P.get<DROPS::Point3DCL>("Levelset.PosDrop")[0];
  static double y0 = P.get<DROPS::Point3DCL>("Levelset.PosDrop")[1];
  static double R = P.get<DROPS::Point3DCL>("Levelset.RadDrop")[0];  */
//  double x = p[0];
  double y = p[1];  
  
  int c= (int) (y * 2.0);
  switch(c)
  {
    case 0: return 0; break;
    case 1: return (y-0.5)*(1-y); break;
    case 2: return -(y-1)*(1.5-y); break;
/*    case 1:
    case 2:
            return (y-0.5)*(1.5-y);
*/    case 3: 
    default:
        return 0;
        break; 
  }
}

double DirichletConst (const DROPS::Point3DCL& p, double )
{
  extern DROPS::ParamCL P;
  double y=p[1]-P.get<DROPS::Point3DCL>("Levelset.PosDrop")[1];
  if (y>0)
    return P.get<double>("Transp.IniCPos");
  else
    return P.get<double>("Transp.IniCNeg");
}

double DirichletPos (const DROPS::Point3DCL&, double )
{
  extern DROPS::ParamCL P;
  return P.get<double>("Transp.IniCPos");
}

double DirichletConstt (const DROPS::Point3DCL& p, double )
{
  extern DROPS::ParamCL P;
  double y=p[1]-P.get<DROPS::Point3DCL>("Levelset.PosDrop")[1];
  if (y>0)
    return P.get<double>("Transp.IniCPos");
  else
    return P.get<double>("Transp.IniCNeg") * P.get<double>("Transp.HNeg") / P.get<double>("Transp.HPos");
}

DROPS::SVectorCL<3> PotentialFlowfield (const DROPS::Point3DCL& p, double )
{  
    extern DROPS::ParamCL P;
    DROPS::SVectorCL<3> ret(0.);
    double x=p[0]-P.get<DROPS::Point3DCL>("Levelset.PosDrop")[0]; double y=p[1]-P.get<DROPS::Point3DCL>("Levelset.PosDrop")[1];
    double r2 = x*x+y*y;
    double r4 = r2*r2;
    double R=P.get<DROPS::Point3DCL>("Levelset.RadDrop")[0];
    static bool test=true;
    if(test) { std::cout << "R = " << R << std::endl; test=false;}
    double R2 = R*R;
    if (R2>=r2){
      ;//nothing, v=0
    }else{
      ret[0] = 1 + R2*(y*y-x*x)/r4;
      ret[1] = -2*R2*(y*x)/r4;
    }
//    std::cout << ret[0] << " " << ret[1] << std::endl;
    return ret;
}

template<int i>
DROPS::SVectorCL<3> StraightFlowfield (const DROPS::Point3DCL&, double )
{  
    DROPS::SVectorCL<3> ret(0.);
    ret[i] = 1.0;
    return ret;
}

template<int i>
double cylinderdistance( const DROPS::Point3DCL& p, double)
{
    extern DROPS::ParamCL P;
    double x=p[0]-P.get<DROPS::Point3DCL>("Levelset.PosDrop")[0];
    double y=p[1]-P.get<DROPS::Point3DCL>("Levelset.PosDrop")[1];
    double z=p[2]-P.get<DROPS::Point3DCL>("Levelset.PosDrop")[2];
    double R=P.get<DROPS::Point3DCL>("Levelset.RadDrop")[i];
    if (i == 0)
        return std::sqrt(std::abs(y*y+z*z)) - R;
    if (i == 1)
    	return std::sqrt(std::abs(x*x+z*z)) - R;
    if (i == 2)
        return std::sqrt(std::abs(y*y+x*x)) - R;
}

template<int i>
double planedistance( const DROPS::Point3DCL& p, double)
{
    extern DROPS::ParamCL P;
    double x=p[i]-P.get<DROPS::Point3DCL>("Levelset.PosDrop")[i];
    return x;
}

double ellipsoid_left( const DROPS::Point3DCL& p)
{
	extern DROPS::ParamCL P;
	DROPS::Point3DCL pos = P.get<DROPS::Point3DCL>("Exp.PosLeftDrop");
	DROPS::Point3DCL rad = P.get<DROPS::Point3DCL>("Exp.RadLeftDrop");
	DROPS::Point3DCL d = p - pos;
	const double avgRad= cbrt(rad[0]*rad[1]*rad[2]);
	d/= rad;
	return std::abs( avgRad)*d.norm() - avgRad;
}

double ellipsoid_right( const DROPS::Point3DCL& p)
{
	extern DROPS::ParamCL P;
	DROPS::Point3DCL pos = P.get<DROPS::Point3DCL>("Exp.PosRightDrop");
	DROPS::Point3DCL rad = P.get<DROPS::Point3DCL>("Exp.RadRightDrop");
	DROPS::Point3DCL d = p - pos;
	const double avgRad= cbrt(rad[0]*rad[1]*rad[2]);
	d/= rad;
	return std::abs( avgRad)*d.norm() - avgRad;
}

double TwoEllipsoid (const DROPS::Point3DCL& p, double)
{
	double EllRight = ellipsoid_right(p);
	double EllLeft = ellipsoid_left(p);
	if (EllRight > EllLeft)
		return EllLeft;
	return EllRight;

}


double GrowingCylinder (const DROPS::Point3DCL& p, double)
{
	extern DROPS::ParamCL P;
	DROPS::Point3DCL pos = P.get<DROPS::Point3DCL>("Levelset.PosDrop");
	DROPS::Point3DCL rad = P.get<double>("Exp.GrowDrop")*P.get<DROPS::Point3DCL>("Levelset.RadDrop");
	DROPS::Point3DCL d = p - pos;
	const double avgRad= sqrt(rad[0]*rad[1]);
	//d/= rad;
    const double dnorm2D = std::sqrt(d[0]*d[0]+d[1]*d[1]);

	return dnorm2D - avgRad;

}


double GrowingEllipsoid (const DROPS::Point3DCL& p, double)
{
	extern DROPS::ParamCL P;
	DROPS::Point3DCL pos = P.get<DROPS::Point3DCL>("Levelset.PosDrop");
	DROPS::Point3DCL rad = P.get<double>("Exp.GrowDrop")*P.get<DROPS::Point3DCL>("Levelset.RadDrop");
	DROPS::Point3DCL d = p - pos;
	const double avgRad= cbrt(rad[0]*rad[1]*rad[2]);
	d/= rad;
    const double dnorm = std::sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);

	return avgRad*dnorm - avgRad;

}

double SimpleEllipsoid (const DROPS::Point3DCL& p, double)
{
	extern DROPS::ParamCL P;
	DROPS::Point3DCL pos = P.get<DROPS::Point3DCL>("Levelset.PosDrop");
	DROPS::Point3DCL rad = P.get<DROPS::Point3DCL>("Levelset.RadDrop");
	DROPS::Point3DCL d = p - pos;
	const double avgRad= cbrt(rad[0]*rad[1]*rad[2]);
	d/= rad;
    const double dnorm = std::sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);

	return avgRad*dnorm - avgRad;
}


double TubeEll (const DROPS::Point3DCL& p, double)
{
	extern DROPS::ParamCL P;
	static DROPS::Point3DCL pos = P.get<DROPS::Point3DCL>("Levelset.PosDrop");
	static DROPS::Point3DCL rad = P.get<DROPS::Point3DCL>("Levelset.RadDrop");
	DROPS::Point3DCL d = p - pos;
	const double avgRad= cbrt(rad[0]*rad[1]*rad[2]);
	d/= rad;
	const double dnorm = std::sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);

	return avgRad*dnorm - avgRad;
}

template<int i>
double Tube (const DROPS::Point3DCL& p, double)
{
	extern DROPS::ParamCL P;
	DROPS::Point3DCL ellp = p;
	DROPS::Point3DCL ellp2 = p;
	static DROPS::Point3DCL pos = P.get<DROPS::Point3DCL>("Levelset.PosDrop");
	ellp -= pos;
	if (std::abs(ellp[i]) <= 0.5)
		return cylinderdistance<i>(p, 0.0);
	if (ellp[i] >= 0.5)
		ellp2[i] -= 0.5;
	else
		ellp2[i] += 0.5;
	return TubeEll(ellp2, 0.0);
}


double Prokert_InitConc_II (const DROPS::Point3DCL& p, double)
{
	extern DROPS::ParamCL P;
	static DROPS::Point3DCL pos = P.get<DROPS::Point3DCL>("Exp.PosConcCenter");
    DROPS::Point3DCL d = p - pos;
	static double a = P.get<double>("Exp.ExpFactorA",1.0);
	// static double b = P.get<double>("Exp.FactorB",1.0);
    const double exparg = - a * d.norm() * d.norm();
	return exp(exparg);
}





static DROPS::RegisterScalarFunction regsca_ellipsoid("GrowingEllipsoid", GrowingEllipsoid);
static DROPS::RegisterScalarFunction regsca_simpleellipsoid("SimpleEllipsoid", SimpleEllipsoid);
static DROPS::RegisterScalarFunction regsca_heavi("Heaviside", Heaviside);
static DROPS::RegisterScalarFunction regsca_linearconc("linearConc", linearConc);
static DROPS::RegisterScalarFunction regsca_tdistx("Tubex", Tube<0>);
static DROPS::RegisterScalarFunction regsca_tdisty("Tubey", Tube<1>);
static DROPS::RegisterScalarFunction regsca_tdistz("Tubez", Tube<2>);
static DROPS::RegisterScalarFunction regsca_twoellipsoid("TwoEllipsoids", TwoEllipsoid);
static DROPS::RegisterScalarFunction regsca_inicneg("IniCnegFct", Initialcneg);
static DROPS::RegisterScalarFunction regsca_iniprokert("IniProkertII", Prokert_InitConc_II);
static DROPS::RegisterScalarFunction regsca_inicpos("IniCposFct", Initialcpos);
static DROPS::RegisterScalarFunction regsca_reaction("ReactionFct", Reaction);
static DROPS::RegisterScalarFunction regsca_zero("ZeroFct", ZeroFct);
static DROPS::RegisterScalarFunction regsca_Rhs("Rhs", Rhs);
static DROPS::RegisterScalarFunction regsca_dir("Dirichlet", DirichletPos);
static DROPS::RegisterScalarFunction regsca_sinx("Sinx", Sinx);
static DROPS::RegisterScalarFunction regsca_tl("TransverseLinear", TransverseLinear);
static DROPS::RegisterScalarFunction regsca_dirt("Dirichlett", DirichletPos);
//static DROPS::RegisterScalarFunction regsca_dir("Dirichlet", DirichletConst);
//static DROPS::RegisterScalarFunction regsca_dirt("Dirichlett", DirichletConstt);
static DROPS::RegisterScalarFunction regsca_cdistx("cylinderdistancex", cylinderdistance<0>);
static DROPS::RegisterScalarFunction regsca_cdisty("cylinderdistancey", cylinderdistance<1>);
static DROPS::RegisterScalarFunction regsca_cdistz("cylinderdistancez", cylinderdistance<2>);
static DROPS::RegisterScalarFunction regsca_pdistx("planedistancex", planedistance<0>);
static DROPS::RegisterScalarFunction regsca_pdisty("planedistancey", planedistance<1>);
static DROPS::RegisterScalarFunction regsca_pdistz("planedistancez", planedistance<2>);
static DROPS::RegisterScalarFunction regsca_torus("torus", Torus);
static DROPS::RegisterVectorFunction regvec_potflow("potflow", PotentialFlowfield);
static DROPS::RegisterVectorFunction regvec_strflowx("straightflowx", StraightFlowfield<0>);
static DROPS::RegisterVectorFunction regvec_strflowy("straightflowy", StraightFlowfield<1>);
static DROPS::RegisterVectorFunction regvec_strflowz("straightflowz", StraightFlowfield<2>);
static DROPS::RegisterVectorFunction regvec_zero("ZeroFct", ZeroFct_vec);


