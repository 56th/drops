/// \file transportCoeff.cpp
/// \brief boundary and source functions for the transport Problem
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

double Special (const DROPS::Point3DCL& p, double )
{
    extern DROPS::ParamCL P;
    static DROPS::Point3DCL x0 = P.get<DROPS::Point3DCL>("Levelset.PosDrop");
    static DROPS::Point3DCL shift = P.get<DROPS::Point3DCL>("Exp.Shift");
    static double R = P.get<DROPS::Point3DCL>("Levelset.RadDrop")[0];
    const double r = 0.05;
    DROPS::Point3DCL X = p-x0-shift;

    DROPS::Point3DCL Y(std::abs(X[0]),std::abs(X[1]),std::abs(X[2]));
    
    int order[3];
    for (int i = 0; i < 3; ++i)
        order[i] = i;

    if (Y[order[0]] < Y[order[1]])
        std::swap(order[0],order[1]);

    if (Y[order[1]] < Y[order[2]])
        std::swap(order[1],order[2]);

    if (Y[order[0]] < Y[order[1]])
        std::swap(order[0],order[1]);

    //cnt: how many are outsider inner(st) box
    int cnt = 0;
    for (int i = 0; i < 3; ++i)
    {
        if (Y[i] > R - r) cnt++;;
    }

    if (cnt <= 1)
    {
        return Y[order[0]] - R;
    }
    else if (cnt == 2)
    {
        double rx = Y[order[0]] - (R - r);
        double ry = Y[order[1]] - (R - r);
        return std::sqrt(rx*rx+ry*ry) - r;
    }
    else
    {
        double rx = Y[order[0]] - (R - r);
        double ry = Y[order[1]] - (R - r);
        double rz = Y[order[2]] - (R - r);
        return std::sqrt(rx*rx+ry*ry+rz*rz) - r;
    }

}

static DROPS::RegisterScalarFunction regsca_special("Special", Special);

double Initialcneg (const DROPS::Point3DCL& , double )
{
    extern DROPS::ParamCL P;
    return P.get<double>("Transp.IniCNeg");
}

double Initialcpos (const DROPS::Point3DCL& , double )
{
    extern DROPS::ParamCL P;
    return P.get<double>("Transp.IniCPos");
}

double Reaction (const DROPS::Point3DCL& , double )
{
  return 0.0;
}

double ZeroFct (const DROPS::Point3DCL& ,  double )
{
    return 0.;
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


double cylinderdistance( const DROPS::Point3DCL& p, double)
{
    extern DROPS::ParamCL P;
    double x=p[0]-P.get<DROPS::Point3DCL>("Levelset.PosDrop")[0]; double y=p[1]-P.get<DROPS::Point3DCL>("Levelset.PosDrop")[1];
    double R=P.get<DROPS::Point3DCL>("Levelset.RadDrop")[0];
    static bool test=true;
    if(test) { std::cout << "R = " << R << std::endl; test=false;}
    return std::sqrt(std::abs(x*x+y*y)) - R;
}

template<int i>
double planedistance( const DROPS::Point3DCL& p, double)
{
    extern DROPS::ParamCL P;
    double x=p[i]-P.get<DROPS::Point3DCL>("Levelset.PosDrop")[i];
    return x;
}


static DROPS::RegisterScalarFunction regsca_inicneg("IniCnegFct", Initialcneg);
static DROPS::RegisterScalarFunction regsca_inicpos("IniCposFct", Initialcpos);
static DROPS::RegisterScalarFunction regsca_reaction("ReactionFct", Reaction);
static DROPS::RegisterScalarFunction regsca_zero("ZeroFct", ZeroFct);
static DROPS::RegisterScalarFunction regsca_Rhs("Rhs", Rhs);
static DROPS::RegisterScalarFunction regsca_dir("Dirichlet", DirichletPos);
static DROPS::RegisterScalarFunction regsca_dirt("Dirichlett", DirichletPos);
//static DROPS::RegisterScalarFunction regsca_dir("Dirichlet", DirichletConst);
//static DROPS::RegisterScalarFunction regsca_dirt("Dirichlett", DirichletConstt);
static DROPS::RegisterScalarFunction regsca_cdist("cylinderdistance", cylinderdistance);
static DROPS::RegisterScalarFunction regsca_pdistx("planedistancex", planedistance<0>);
static DROPS::RegisterScalarFunction regsca_pdisty("planedistancey", planedistance<1>);
static DROPS::RegisterScalarFunction regsca_pdistz("planedistancez", planedistance<2>);
static DROPS::RegisterVectorFunction regvec_potflow("potflow", PotentialFlowfield);
static DROPS::RegisterVectorFunction regvec_strflowx("straightflowx", StraightFlowfield<0>);
static DROPS::RegisterVectorFunction regvec_strflowy("straightflowy", StraightFlowfield<1>);
static DROPS::RegisterVectorFunction regvec_strflowz("straightflowz", StraightFlowfield<2>);


