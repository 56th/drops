/// \file transportCoeff.cpp
/// \brief boundary and source functions for the transport Problem | copy of transportCoeff.cpp ... so far
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

DROPS::Point3DCL MeshSize()
{
    extern DROPS::ParamCL P;
    static DROPS::Point3DCL meshsize( norm(P.get<DROPS::Point3DCL>("Mesh.E1")),
                                      norm(P.get<DROPS::Point3DCL>("Mesh.E2")),
                                      norm(P.get<DROPS::Point3DCL>("Mesh.E3")));
    return meshsize;
}


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

double Inineg (const DROPS::Point3DCL& p, double )
{
    extern DROPS::ParamCL P;
    return P.get<double>("Transp.IniCNeg") * (p[0]*p[0]+p[1]*p[1]);
}

double Inipos (const DROPS::Point3DCL& p, double )
{
    extern DROPS::ParamCL P;
    return P.get<double>("Transp.IniCPos") * (1.0-p[0]*p[0]-p[1]*p[1]);
}

double Reaction (const DROPS::Point3DCL& , double )
{
  return 0.0;
}

double ZeroFct (const DROPS::Point3DCL& ,  double )
{  
    return 0.;
}

double tid_ZeroFct (const DROPS::Point3DCL& ,  double )
{  
    return 0.;
}

double Rhs (const DROPS::Point3DCL& ,  double )
{  
    return 0.;
}


double linear_in_t (const DROPS::Point3DCL&,  double t )
{  
    return t;
}

double linear_in_x (const DROPS::Point3DCL& p,  double  )
{  
    return p[0];
}

double linear_in_y (const DROPS::Point3DCL& p,  double  )
{  
    return p[1];
}

double linear_in_z (const DROPS::Point3DCL& p,  double  )
{  
    return p[2];
}


double Ten (const DROPS::Point3DCL&,  double )
{  
    return 10;
}

double LinT (const DROPS::Point3DCL& p,  double t)
{  
    return 10*t*p[1];
}

double LinmT (const DROPS::Point3DCL& p,  double t)
{  
    return 10.0*(1.0-t)*p[1];
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
DROPS::SVectorCL<3> CounterFlowfield (const DROPS::Point3DCL&, double )
{  
    extern DROPS::ParamCL P;
    double c=P.get<double>("Exp.CounterFlowVel");
    DROPS::SVectorCL<3> ret(0.);
    ret[i] = -c;
    return ret;
}


template<int i>
DROPS::SVectorCL<3> StraightSinus (const DROPS::Point3DCL&, double t)
{  
    DROPS::SVectorCL<3> ret(0.);
    ret[i] = cos(t);
    return ret;
}


double cylinderdistance( const DROPS::Point3DCL& p ,  double )
{
    extern DROPS::ParamCL P;
    double x=p[0]-P.get<DROPS::Point3DCL>("Levelset.PosDrop")[0]; double y=p[1]-P.get<DROPS::Point3DCL>("Levelset.PosDrop")[1];
    double R=P.get<DROPS::Point3DCL>("Levelset.RadDrop")[0];
    static bool test=true;
    if(test) { std::cout << "R = " << R << std::endl; test=false;}
    return std::sqrt(std::abs(x*x+y*y)) - R;
}

template<int i>
double planedistance( const DROPS::Point3DCL& p ,  double )
{
    extern DROPS::ParamCL P;
    double x=p[i]-P.get<DROPS::Point3DCL>("Levelset.PosDrop")[i];
    return x;
}

double MovingEllipsoid( const DROPS::Point3DCL& p ,  double t )
{
    extern DROPS::ParamCL P;
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    static DROPS::Point3DCL rad (P.get<DROPS::Point3DCL>("Levelset.RadDrop"));
    DROPS::Point3DCL p2(DROPS::MakePoint3D(pos[0]+sin(t),pos[1],pos[2]));
    DROPS::Point3DCL p3(DROPS::MakePoint3D(pos[0]+1.0+sin(t),pos[1],pos[2]));
    DROPS::Point3DCL p4(DROPS::MakePoint3D(pos[0]-1.0+sin(t),pos[1],pos[2]));
    DROPS::Point3DCL tmp ((p-p2)/rad);
    DROPS::Point3DCL tmp2 ((p-p3)/rad);
    DROPS::Point3DCL tmp3 ((p-p4)/rad);
    return std::min(std::min(tmp.norm(),tmp2.norm()),tmp3.norm()) - 1.0;
}

static DROPS::RegisterScalarFunction regsca_lint("linear_in_t", linear_in_t);
static DROPS::RegisterScalarFunction regsca_linx("linear_in_x", linear_in_x);
static DROPS::RegisterScalarFunction regsca_liny("linear_in_y", linear_in_y);
static DROPS::RegisterScalarFunction regsca_linz("linear_in_z", linear_in_z);
static DROPS::RegisterScalarFunction regsca_inicneg("IniCnegFct", Initialcneg);
static DROPS::RegisterScalarFunction regsca_inicpos("IniCposFct", Initialcpos);
static DROPS::RegisterScalarFunction regsca_inineg("inineg", Inineg);
static DROPS::RegisterScalarFunction regsca_inipos("inipos", Inipos);
static DROPS::RegisterScalarFunction regsca_reaction("ReactionFct", Reaction);
static DROPS::RegisterScalarFunction regsca_zero("ZeroFct", ZeroFct);
static DROPS::RegisterScalarFunction regsca_tidzero("ZeroFct", tid_ZeroFct);
static DROPS::RegisterScalarFunction regsca_Rhs("Rhs", Rhs);
static DROPS::RegisterScalarFunction regsca_Ten("Ten", Ten);
static DROPS::RegisterScalarFunction regsca_LinT("LinT", LinT);
static DROPS::RegisterScalarFunction regsca_LinmT("LinmT", LinmT);
static DROPS::RegisterScalarFunction regsca_dir("Dirichlet", DirichletPos);
static DROPS::RegisterScalarFunction regsca_dirt("Dirichlett", DirichletPos);
//static DROPS::RegisterScalarFunction regsca_dir("Dirichlet", DirichletConst);
//static DROPS::RegisterScalarFunction regsca_dirt("Dirichlett", DirichletConstt);
static DROPS::RegisterScalarFunction regsca_cdist("cylinderdistance", cylinderdistance);
static DROPS::RegisterScalarFunction regsca_pdistx("planedistancex", planedistance<0>);
static DROPS::RegisterScalarFunction regsca_pdisty("planedistancey", planedistance<1>);
static DROPS::RegisterScalarFunction regsca_pdistz("planedistancez", planedistance<2>);
static DROPS::RegisterScalarFunction regsca_ell2("MovingEllipsoid", MovingEllipsoid);
static DROPS::RegisterVectorFunction regvec_potflow("potflow", PotentialFlowfield);
static DROPS::RegisterVectorFunction regvec_strsinus("straightsinus", StraightSinus<0>);
static DROPS::RegisterVectorFunction regvec_strflowx("straightflowx", StraightFlowfield<0>);
static DROPS::RegisterVectorFunction regvec_strflowy("straightflowy", StraightFlowfield<1>);
static DROPS::RegisterVectorFunction regvec_strflowz("straightflowz", StraightFlowfield<2>);

static DROPS::RegisterVectorFunction regvec_ctrflowx("counterflowx", CounterFlowfield<0>);
static DROPS::RegisterVectorFunction regvec_ctrflowy("counterflowy", CounterFlowfield<1>);
static DROPS::RegisterVectorFunction regvec_ctrflowz("counterflowz", CounterFlowfield<2>);

double gtc1( double t)
{
    extern DROPS::ParamCL P;
    static bool bcc (P.get<int>("TestCase1.bcc",0));
    if (bcc)
    {
        return t;
    }
    else
    {
        static double k (P.get<double>("TestCase1.k"));
        return sin(k*M_PI*t);
    }
}

double ddtgtc1( double t)
{
    extern DROPS::ParamCL P;
    static bool bcc (P.get<int>("TestCase1.bcc",0));
    if (bcc)
    {
        return 1;
    }
    else
    {
        static double k (P.get<double>("TestCase1.k"));
        return k*M_PI*cos(k*M_PI*t);
    }
}


double testcase1_shift( double v, double t)
{
    extern DROPS::ParamCL P;
    static double r (P.get<double>("TestCase1.r"));
    static bool sinw (P.get<int>("TestCase1.sinmove",0));
    if (sinw)
    {
        return 0.25/M_PI * v * sin(2.0*M_PI*t);
    }
    else
    {
        if (t<r)
            return v*t;
        else
            return v*(2*r-t);
    }
}

DROPS::Point3DCL TestCase1_vel (const DROPS::Point3DCL& ,  double t )
{  
    extern DROPS::ParamCL P;
    static double v (P.get<double>("TestCase1.v"));
    static DROPS::Point3DCL vel( DROPS::MakePoint3D(v,0,0));
    static double r (P.get<double>("TestCase1.r"));
    static bool sinw (P.get<int>("TestCase1.sinmove",0));

    if (sinw)
    {
        return DROPS::MakePoint3D(0.5 * v * cos(2*M_PI*t),0,0);
    }
    else
    {
        if (t > r) 
            vel = DROPS::MakePoint3D(-v,0,0);
    }
    return vel;
}

void TestCase1_push()
{
    static bool pushed = false;
    if (!pushed)
    {
        extern DROPS::ParamCL P;
        static double hpos (P.get<double>("Transp.HPos"));
        static double apos (P.get<double>("Transp.DiffPos"));
        static double hneg (P.get<double>("Transp.HNeg"));
        static double aneg (P.get<double>("Transp.DiffNeg"));

        static double R (P.get<DROPS::Point3DCL>("Levelset.RadDrop")[0]);

        const double b = -M_PI*sin(M_PI*R)/(2.0*R)*apos/aneg;
        P.put_if_unset<double>("TestCase1.b",b);

        const double a = cos(M_PI*R)*hpos/hneg-b*R*R;
        P.put_if_unset<double>("TestCase1.a",a);

        std::cout << " a = " << a << std::endl;
        std::cout << " b = " << b << std::endl;
        pushed = true;
    }
}

double TestCase1_lset( const DROPS::Point3DCL& p ,  double t )
{
    extern DROPS::ParamCL P;
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    static DROPS::Point3DCL rad (P.get<DROPS::Point3DCL>("Levelset.RadDrop"));
    static double v (P.get<double>("TestCase1.v"));
    DROPS::Point3DCL p2(DROPS::MakePoint3D(pos[0]+testcase1_shift(v,t),pos[1],pos[2]));
    // DROPS::Point3DCL p3(DROPS::MakePoint3D(pos[0]+2.0+testcase1_shift(v,t),pos[1],pos[2]));
    // DROPS::Point3DCL p4(DROPS::MakePoint3D(pos[0]-2.0+testcase1_shift(v,t),pos[1],pos[2]));
    DROPS::Point3DCL tmp ((p-p2)/rad);
    // DROPS::Point3DCL tmp2 ((p-p3)/rad);
    // DROPS::Point3DCL tmp3 ((p-p4)/rad);
    // return std::min(std::min(tmp.norm(),tmp2.norm()),tmp3.norm()) - 1.0;
    return tmp.norm() - 1.0;
}



double TestCase1_rhs_positive (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    // static double hpos (P.get<double>("Transp.HPos"));
    static double apos (P.get<double>("Transp.DiffPos"));
    // static double hneg (P.get<double>("Transp.HNeg"));
    // static double aneg (P.get<double>("Transp.DiffNeg"));
    // static double k (P.get<double>("TestCase1.k"));
    static double C (P.get<double>("TestCase1.C"));
    static double v (P.get<double>("TestCase1.v"));

    DROPS::Point3DCL shiftedpos;
    shiftedpos = pos;
    shiftedpos[0] += testcase1_shift(v,t);
    DROPS::Point3DCL diff = p - shiftedpos;
    const double absx = diff.norm();
    const double timefactor = gtc1(t);

    const double diffpart1 = M_PI*M_PI*cos(M_PI*absx) + sin(M_PI*absx)*M_PI*2.0/absx;
    const double diffpart2 = apos * diffpart1 * timefactor;
    const double timeder = cos(M_PI*absx) * ddtgtc1(t);
    //convection missing
    return C*(timeder + diffpart2);
}


double TestCase1_rhs_negative (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    // static double hpos (P.get<double>("Transp.HPos"));
    // static double apos (P.get<double>("Transp.DiffPos"));
    // static double hneg (P.get<double>("Transp.HNeg"));
    static double aneg (P.get<double>("Transp.DiffNeg"));
    // static double k (P.get<double>("TestCase1.k"));
    TestCase1_push();
    static double a (P.get<double>("TestCase1.a"));
    static double b (P.get<double>("TestCase1.b"));
    static double C (P.get<double>("TestCase1.C"));
    static double v (P.get<double>("TestCase1.v"));
    DROPS::Point3DCL shiftedpos;
    shiftedpos = pos;
    shiftedpos[0] += testcase1_shift(v,t);
    DROPS::Point3DCL diff = p - shiftedpos;
    const double absx = diff.norm();
    const double timefactor = gtc1(t);

    const double diffpart1 = -6*b;
    const double diffpart2 = aneg * diffpart1 * timefactor;

    const double timeder = (a+b*absx*absx) * ddtgtc1(t);
    //convection missing
    return C*(timeder + diffpart2);
}

double TestCase1_sol_positive (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    // static double k (P.get<double>("TestCase1.k"));
    static double C (P.get<double>("TestCase1.C"));
    static double v (P.get<double>("TestCase1.v"));
    DROPS::Point3DCL shiftedpos;
    shiftedpos = pos;
    shiftedpos[0] += testcase1_shift(v,t);
    DROPS::Point3DCL diff = p - shiftedpos;
    const double absx = diff.norm();
    const double timefactor = gtc1(t);

    return C*timefactor*cos(M_PI*absx);
}

double TestCase1_sol_negative (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    // static double k (P.get<double>("TestCase1.k"));
    TestCase1_push();
    static double a (P.get<double>("TestCase1.a"));
    static double b (P.get<double>("TestCase1.b"));
    static double C (P.get<double>("TestCase1.C"));
    static double v (P.get<double>("TestCase1.v"));
    DROPS::Point3DCL shiftedpos;
    shiftedpos = pos;
    shiftedpos[0] += testcase1_shift(v,t);
    DROPS::Point3DCL diff = p - shiftedpos;
    const double absx = diff.norm();
    const double timefactor = gtc1(t);

    // static double lasttime = -1;
    // if (t != lasttime){
    //     std::cout << " timefactor = " << timefactor << std::endl;
    //     lasttime = t;
    // }
    return C*timefactor*(a+b*absx*absx);
}

static DROPS::RegisterScalarFunction regsca_testcase1_RhsPos("testcase1_rhs_pos", TestCase1_rhs_positive);
static DROPS::RegisterScalarFunction regsca_testcase1_RhsNeg("testcase1_rhs_neg", TestCase1_rhs_negative);
static DROPS::RegisterScalarFunction regsca_testcase1_SolPos("testcase1_sol_pos", TestCase1_sol_positive);
static DROPS::RegisterScalarFunction regsca_testcase1_SolNeg("testcase1_sol_neg", TestCase1_sol_negative);
static DROPS::RegisterScalarFunction regsca_testcase1_lset("testcase1_lset", TestCase1_lset);
static DROPS::RegisterVectorFunction regvec_testcase1_velocity("testcase1_vel", TestCase1_vel);


double testcase2_shift( double v, double t)
{
    extern DROPS::ParamCL P;
    static double r (P.get<double>("TestCase2.r"));
    if (t<r)
        return v*t;
    else
        return v*(2*r-t);
}

DROPS::Point3DCL TestCase2_vel (const DROPS::Point3DCL& ,  double t )
{  
    extern DROPS::ParamCL P;
    static double v (P.get<double>("TestCase2.v"));
    static DROPS::Point3DCL vel( DROPS::MakePoint3D(v,0,0));
    static double r (P.get<double>("TestCase2.r"));

    if (t > r) 
        vel = DROPS::MakePoint3D(-v,0,0);
    return vel;
}

double TestCase2_lset( const DROPS::Point3DCL& p ,  double t )
{
    extern DROPS::ParamCL P;
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    static DROPS::Point3DCL rad (P.get<DROPS::Point3DCL>("Levelset.RadDrop"));
    static double v (P.get<double>("TestCase2.v"));
    DROPS::Point3DCL p2(DROPS::MakePoint3D(pos[0]+testcase2_shift(v,t),pos[1],pos[2]));
    DROPS::Point3DCL p3(DROPS::MakePoint3D(pos[0]+2.0+testcase2_shift(v,t),pos[1],pos[2]));
    DROPS::Point3DCL p4(DROPS::MakePoint3D(pos[0]-2.0+testcase2_shift(v,t),pos[1],pos[2]));
    DROPS::Point3DCL tmp ((p-p2)/rad);
    DROPS::Point3DCL tmp2 ((p-p3)/rad);
    DROPS::Point3DCL tmp3 ((p-p4)/rad);
    return std::min(std::min(tmp.norm(),tmp2.norm()),tmp3.norm()) - 1.0;
}



double TestCase2_rhs_positive (const DROPS::Point3DCL& ,  double )
{  
    // extern DROPS::ParamCL P;
    // static double v (P.get<double>("TestCase2.v"));
    return 1.0;
}


double TestCase2_rhs_negative (const DROPS::Point3DCL&  ,  double )
{  
    // extern DROPS::ParamCL P;
    // static double v (P.get<double>("TestCase2.v"));
    return 2;
}

double TestCase2_sol_positive (const DROPS::Point3DCL&  ,  double t)
{  
    return t;
}

double TestCase2_sol_negative (const DROPS::Point3DCL&  ,  double t)
{  
    return 2*t;
}

static DROPS::RegisterScalarFunction regsca_testcase2_RhsPos("testcase2_rhs_pos", TestCase2_rhs_positive);
static DROPS::RegisterScalarFunction regsca_testcase2_RhsNeg("testcase2_rhs_neg", TestCase2_rhs_negative);
static DROPS::RegisterScalarFunction regsca_testcase2_SolPos("testcase2_sol_pos", TestCase2_sol_positive);
static DROPS::RegisterScalarFunction regsca_testcase2_SolNeg("testcase2_sol_neg", TestCase2_sol_negative);
static DROPS::RegisterScalarFunction regsca_testcase2_lset("testcase2_lset", TestCase2_lset);
static DROPS::RegisterVectorFunction regvec_testcase2_velocity("testcase2_vel", TestCase2_vel);


double testcase3_shift( double v, double t)
{
    extern DROPS::ParamCL P;
    static double r (P.get<double>("TestCase3.r"));
    if (t<r)
        return v*t;
    else
        return v*(2*r-t);
}

DROPS::Point3DCL TestCase3_vel (const DROPS::Point3DCL& ,  double t )
{  
    extern DROPS::ParamCL P;
    static double v (P.get<double>("TestCase3.v"));
    static DROPS::Point3DCL vel( DROPS::MakePoint3D(v,0,0));
    static double r (P.get<double>("TestCase3.r"));

    if (t > r) 
        vel = DROPS::MakePoint3D(-v,0,0);
    return vel;
}

double TestCase3_lset( const DROPS::Point3DCL& p ,  double t )
{
    extern DROPS::ParamCL P;
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    static DROPS::Point3DCL rad (P.get<DROPS::Point3DCL>("Levelset.RadDrop"));
    static double v (P.get<double>("TestCase3.v"));
    DROPS::Point3DCL p2(DROPS::MakePoint3D(pos[0]+testcase3_shift(v,t),pos[1],pos[2]));
    DROPS::Point3DCL p3(DROPS::MakePoint3D(pos[0]+2.0+testcase3_shift(v,t),pos[1],pos[2]));
    DROPS::Point3DCL p4(DROPS::MakePoint3D(pos[0]-2.0+testcase3_shift(v,t),pos[1],pos[2]));
    DROPS::Point3DCL tmp ((p-p2)/rad);
    DROPS::Point3DCL tmp2 ((p-p3)/rad);
    DROPS::Point3DCL tmp3 ((p-p4)/rad);
    return std::min(std::min(tmp.norm(),tmp2.norm()),tmp3.norm()) - 1.0;
}



double TestCase3_rhs_positive (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    static double apos (P.get<double>("Transp.DiffPos"));
    static double v (P.get<double>("TestCase3.v"));
    static double C (P.get<double>("TestCase3.C"));
    static bool noconv (P.get<int>("Transp.NoConvection")==1);

    if (!noconv)
        return C*apos*M_PI*M_PI*sin(M_PI*(p[0]-v*t));
    else
    {
        const double factor1 = C*M_PI;
        const double sinterm = apos*M_PI*sin(M_PI*(p[0]-v*t));
        const double costerm = -v * cos(M_PI*(p[0]-v*t));
        return factor1*(sinterm+costerm);
    }
}


double TestCase3_rhs_negative (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    static double hneg (P.get<double>("Transp.HNeg"));
    static double hpos (P.get<double>("Transp.HPos"));
    static double aneg (P.get<double>("Transp.DiffNeg"));
    static double C (P.get<double>("TestCase3.C"));
    static double v (P.get<double>("TestCase3.v"));
    static bool noconv (P.get<int>("Transp.NoConvection")==1);
    if (!noconv)
        return hpos/hneg*C*aneg*M_PI*M_PI*sin(M_PI*(p[0]-v*t));
    else
    {
        const double factor1 = hpos/hneg*C*M_PI;
        const double sinterm = aneg*M_PI*sin(M_PI*(p[0]-v*t));
        const double costerm = - v * cos(M_PI*(p[0]-v*t));
        return factor1*(sinterm+costerm);
    }
}

double TestCase3_sol_positive (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    static double v (P.get<double>("TestCase3.v"));
    static double C (P.get<double>("TestCase3.C"));
    return C*sin(M_PI*(p[0]-v*t));
}

double TestCase3_sol_negative (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    static double v (P.get<double>("TestCase3.v"));
    static double C (P.get<double>("TestCase3.C"));
    static double hneg (P.get<double>("Transp.HNeg"));
    static double hpos (P.get<double>("Transp.HPos"));
    return hpos/hneg*C* sin(M_PI*(p[0]-v*t));
}

static DROPS::RegisterScalarFunction regsca_testcase3_RhsPos("testcase3_rhs_pos", TestCase3_rhs_positive);
static DROPS::RegisterScalarFunction regsca_testcase3_RhsNeg("testcase3_rhs_neg", TestCase3_rhs_negative);
static DROPS::RegisterScalarFunction regsca_testcase3_SolPos("testcase3_sol_pos", TestCase3_sol_positive);
static DROPS::RegisterScalarFunction regsca_testcase3_SolNeg("testcase3_sol_neg", TestCase3_sol_negative);
static DROPS::RegisterScalarFunction regsca_testcase3_lset("testcase3_lset", TestCase3_lset);
static DROPS::RegisterVectorFunction regvec_testcase3_velocity("testcase3_vel", TestCase3_vel);


DROPS::Point3DCL Testcase4_vel (const DROPS::Point3DCL& ,  double )
{  
    return DROPS::MakePoint3D(1.0,0,0);
}

double Testcase4_rhs_positive (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    static double a (P.get<double>("Transp.DiffPos"));
    return a*M_PI*M_PI*sin(M_PI*(p[0]-t));
}

// double Testcase4_neumann_left (const DROPS::Point3DCL& p ,  double t)
// {  
//     return 2.0*t;
// }

// double Testcase4_neumann_right (const DROPS::Point3DCL& p ,  double t)
// {  
//     return 2.0*(2.0-t);
// }


double Testcase4_rhs_negative (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    static double a (P.get<double>("Transp.DiffNeg"));
    return a*M_PI*M_PI*sin(M_PI*(p[0]-t));
}

double Testcase4_sol_positive (const DROPS::Point3DCL& p ,  double t)
{  
    return sin(M_PI*(p[0]-t));
    // return sin(t);
}

double Testcase4_sol_negative (const DROPS::Point3DCL& p ,  double t)
{  
    return sin(M_PI*(p[0]-t));
}

static DROPS::RegisterScalarFunction regsca_testcase4_RhsPos("testcase4_rhs_pos", Testcase4_rhs_positive);
static DROPS::RegisterScalarFunction regsca_testcase4_RhsNeg("testcase4_rhs_neg", Testcase4_rhs_negative);
static DROPS::RegisterScalarFunction regsca_testcase4_SolPos("testcase4_sol_pos", Testcase4_sol_positive);
static DROPS::RegisterScalarFunction regsca_testcase4_SolNeg("testcase4_sol_neg", Testcase4_sol_negative);
static DROPS::RegisterVectorFunction regvec_testcase4_velocity("testcase4_vel", Testcase4_vel);

void TestCase5_push()
{
    static bool pushed = false;
    if (!pushed)
    {
        extern DROPS::ParamCL P;
        static double hpos (P.get<double>("Transp.HPos"));
        static double apos (P.get<double>("Transp.DiffPos"));
        static double hneg (P.get<double>("Transp.HNeg"));
        static double aneg (P.get<double>("Transp.DiffNeg"));

        static double xhat (P.get<DROPS::Point3DCL>("Levelset.RadDrop")[0]);

        const double b = 0.5 / (xhat*xhat) * (apos/aneg * M_PI *cos(M_PI*xhat)-hpos/hneg*sin(M_PI*xhat)/xhat);
        const double a = hpos/hneg * sin(M_PI*xhat)/xhat - b * xhat*xhat;

        P.put_if_unset<double>("TestCase5.b",b);
        P.put_if_unset<double>("TestCase5.a",a);

        std::cout << " a = " << a << std::endl;
        std::cout << " b = " << b << std::endl;

        pushed = true;
    }
}

double TestCase5_scalar_vel (double t )
{  
    extern DROPS::ParamCL P;
    static std::string v_case (P.get<std::string>("TestCase5.velocity"));
    
    if (v_case == "no")
        return 0;

    if (v_case == "planar")
        return 0.25;

    if (v_case == "curved")
        return 0.5*cos(M_PI*2.0*t);

    return 0.0;
}

double testcase5_shift(double t)
{
    extern DROPS::ParamCL P;
    static std::string v_case (P.get<std::string>("TestCase5.velocity"));
    
    double shift;
    if (v_case == "no")
        shift = 0.0;

    if (v_case == "planar")
        shift = 0.25*t;

    if (v_case == "curved")
        shift = 0.25/M_PI*sin(M_PI*2.0*t);
    
    return shift;
}

DROPS::Point3DCL Testcase5_vel (const DROPS::Point3DCL& ,  double t )
{  
    return DROPS::MakePoint3D( TestCase5_scalar_vel(t),0,0);
}


double periodic_correction(const double & x, const double & s)
{
    const double dist0 = std::abs(x-s);
    const double dist1 = std::abs(x-s-2.0);
    const double dist2 = std::abs(x-s+2.0);
    
    if (dist0 < dist1 && dist0 < dist2)
        return x-s;
    else if (dist1 < dist2)
        return x-s-2.0;
    else
        return x-s+2.0;
}

// g(y,z)
double testcase5_deformation (const DROPS::Point3DCL& p )
{
    const double & y(p[1]);    // const double & z(p[2]);

    extern DROPS::ParamCL P;
    static int defcase (P.get<int>("TestCase5.DeformationCase"));

    if (defcase == 1)
        return 0.25*y*y*(2-y)*(2-y)-1.0/8.0;
    else
        return 0.0;
}

// dg/dy(y,z)
double testcase5_deformation_dy (const DROPS::Point3DCL& p )
{
    const double & y(p[1]);    // const double & z(p[2]);

    extern DROPS::ParamCL P;
    static int defcase (P.get<int>("TestCase5.DeformationCase"));

    if (defcase == 1)
        return 2.0*y*(1-y);
    else
        return 0.0;
}


// ||nabla g||^2 + 1
double testcase5_relsurf (const DROPS::Point3DCL& p)
{
    const double & y(p[1]);    // const double & z(p[2]);
    extern DROPS::ParamCL P;
    static int defcase (P.get<int>("TestCase5.DeformationCase"));

    if (defcase == 1)
        return 1.0+y*y*(1-y)*(1-y)*(2-y)*(2-y);
    else
        return 1.0;
}

// lap g
double testcase5_curv (const DROPS::Point3DCL& p)
{
    const double & y(p[1]);    // const double & z(p[2]);
    extern DROPS::ParamCL P;
    static int defcase (P.get<int>("TestCase5.DeformationCase"));

    if (defcase == 1)
        return 2.0-6.0*y+3.0*y*y;
    else
        return 0.0;
}


// phi(x,y,z,t))
double Testcase5_lset( const DROPS::Point3DCL& p ,  double t )
{
    extern DROPS::ParamCL P;
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    static double R (P.get<DROPS::Point3DCL>("Levelset.RadDrop")[0]);
    return std::abs(periodic_correction(p[0],pos[0]+testcase5_shift(t)+testcase5_deformation(p)))-R;
    // const double val1 = std::abs(p[0]-pos[0]-testcase5_shift(t)-testcase5_deformation(p))-R;
    // const double val2 = std::abs(p[0]+2.0-pos[0]-testcase5_shift(t)-testcase5_deformation(p))-R;
    // const double val3 = std::abs(p[0]-2.0-pos[0]-testcase5_shift(t)-testcase5_deformation(p))-R;
    // return std::min(val1,std::min(val2,val3));
}

// v (x,t)
double testcase5_scalarfield_pos (const double& x ,  double t)
{  
    extern DROPS::ParamCL P;
    TestCase5_push();
    static double k (P.get<double>("TestCase5.k"));
    return sin(k*M_PI*t) * sin(M_PI*x);
}

// v_x (x,t)
double testcase5_scalarfield_dx_pos (const double& x ,  double t)
{  
    extern DROPS::ParamCL P;
    TestCase5_push();
    static double k (P.get<double>("TestCase5.k"));
    return sin(k*M_PI*t) * M_PI * cos(M_PI*x);
}

// v_t (x,t)
double testcase5_scalarfield_dt_pos (const double & x, const double & t)
{
    extern DROPS::ParamCL P;
    TestCase5_push();
    static double k (P.get<double>("TestCase5.k"));
    return k*M_PI*cos(k*M_PI*t)*sin(M_PI*x);
}

// v_xx (x,t)
double testcase5_scalarfield_dxx_pos (const double & x, const double & t)
{
    extern DROPS::ParamCL P;
    TestCase5_push();
    static double k (P.get<double>("TestCase5.k"));
    return -M_PI*M_PI*sin(k*M_PI*t)*sin(M_PI*x);
}


double Testcase5_rhs_positive (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    const double x = periodic_correction(p[0],pos[0]+testcase5_shift(t)+testcase5_deformation(p));
    static double apos (P.get<double>("Transp.DiffPos"));
    return testcase5_scalarfield_dt_pos(x,t) 
        - apos * ( testcase5_relsurf(p)*testcase5_scalarfield_dxx_pos(x,t)
                   -testcase5_curv(p)*testcase5_scalarfield_dx_pos(x,t)   );
}

double Testcase5_sol_positive (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    const double x = periodic_correction(p[0],pos[0]+testcase5_shift(t)+testcase5_deformation(p));

    return testcase5_scalarfield_pos(x,t); 
}

double Testcase5_sol_dt_positive (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    const double x = periodic_correction(p[0],pos[0]+testcase5_shift(t)+testcase5_deformation(p));

    return testcase5_scalarfield_dt_pos(x,t); 
}

DROPS::Point3DCL Testcase5_sol_grad_positive (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    const double x = periodic_correction(p[0],pos[0]+testcase5_shift(t)+testcase5_deformation(p));
    return testcase5_scalarfield_dx_pos(x,t) * DROPS::MakePoint3D(1.0, -testcase5_deformation_dy(p),0);
}


// v (x,t)
double testcase5_scalarfield_neg (const double& x ,  double t)
{  
    extern DROPS::ParamCL P;
    TestCase5_push();
    static double a (P.get<double>("TestCase5.a"));
    static double b (P.get<double>("TestCase5.b"));
    static double k (P.get<double>("TestCase5.k"));
    return sin(k*M_PI*t) * (a*x+b*x*x*x);
}

// v_x (x,t)
double testcase5_scalarfield_dx_neg (const double& x ,  double t)
{  
    extern DROPS::ParamCL P;
    TestCase5_push();
    static double a (P.get<double>("TestCase5.a"));
    static double b (P.get<double>("TestCase5.b"));
    static double k (P.get<double>("TestCase5.k"));
    return sin(k*M_PI*t) * (a+3*b*x*x);
}

// v_t (x,t)
double testcase5_scalarfield_dt_neg (const double & x, const double & t)
{
    extern DROPS::ParamCL P;
    TestCase5_push();
    static double a (P.get<double>("TestCase5.a"));
    static double b (P.get<double>("TestCase5.b"));
    static double k (P.get<double>("TestCase5.k"));
    return k*M_PI*cos(k*M_PI*t)*(a*x+b*x*x*x);
}

// v_xx (x,t)
double testcase5_scalarfield_dxx_neg (const double & x, const double & t)
{
    extern DROPS::ParamCL P;
    TestCase5_push();
    static double b (P.get<double>("TestCase5.b"));
    static double k (P.get<double>("TestCase5.k"));
    return 6*b*x*sin(k*M_PI*t);
}


double Testcase5_rhs_negative (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;

    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    const double x = periodic_correction(p[0],pos[0]+testcase5_shift(t)+testcase5_deformation(p));
    static double aneg (P.get<double>("Transp.DiffNeg"));

    return testcase5_scalarfield_dt_neg(x,t) 
        - aneg * ( testcase5_relsurf(p)*testcase5_scalarfield_dxx_neg(x,t)
                   -testcase5_curv(p)*testcase5_scalarfield_dx_neg(x,t)   );
}

double Testcase5_sol_negative (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    const double x = periodic_correction(p[0],pos[0]+testcase5_shift(t)+testcase5_deformation(p));
    return testcase5_scalarfield_neg(x,t);
}

double Testcase5_sol_dt_negative (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    const double x = periodic_correction(p[0],pos[0]+testcase5_shift(t)+testcase5_deformation(p));
    return testcase5_scalarfield_dt_neg(x,t);
}

DROPS::Point3DCL Testcase5_sol_grad_negative (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    const double x = periodic_correction(p[0],pos[0]+testcase5_shift(t)+testcase5_deformation(p));
    return testcase5_scalarfield_dx_neg(x,t) * DROPS::MakePoint3D(1.0, -testcase5_deformation_dy(p),0);
}


static DROPS::RegisterScalarFunction regsca_testcase5_RhsPos("testcase5_rhs_pos", Testcase5_rhs_positive);
static DROPS::RegisterScalarFunction regsca_testcase5_RhsNeg("testcase5_rhs_neg", Testcase5_rhs_negative);
static DROPS::RegisterScalarFunction regsca_testcase5_SolPos("testcase5_sol_pos", Testcase5_sol_positive);
static DROPS::RegisterScalarFunction regsca_testcase5_SolNeg("testcase5_sol_neg", Testcase5_sol_negative);
static DROPS::RegisterScalarFunction regsca_testcase5_SoldtPos("testcase5_sol_dt_pos", Testcase5_sol_dt_positive);
static DROPS::RegisterScalarFunction regsca_testcase5_SoldtNeg("testcase5_sol_dt_neg", Testcase5_sol_dt_negative);
static DROPS::RegisterVectorFunction regsca_testcase5_SolgradNeg("testcase5_sol_grad_neg", Testcase5_sol_grad_negative);
static DROPS::RegisterVectorFunction regsca_testcase5_SolgradPos("testcase5_sol_grad_pos", Testcase5_sol_grad_positive);
static DROPS::RegisterScalarFunction regsca_testcase5_lset("testcase5_lset", Testcase5_lset);
static DROPS::RegisterVectorFunction regvec_testcase5_velocity("testcase5_vel", Testcase5_vel);




double gtc6( double t)
{
    extern DROPS::ParamCL P;
    static bool bcc (P.get<int>("TestCase6.bcc",0));
    if (bcc)
    {
        return t;
    }
    else
    {
        static double k (P.get<double>("TestCase6.k"));
        return sin(k*M_PI*t);
    }
}

double ddtgtc6( double t)
{
    extern DROPS::ParamCL P;
    static bool bcc (P.get<int>("TestCase6.bcc",0));
    if (bcc)
    {
        return 1;
    }
    else
    {
        static double k (P.get<double>("TestCase6.k"));
        return k*M_PI*cos(k*M_PI*t);
    }
}


double blendf( double x, double y, double z)
{
    return x*(2-x) * y*(2-y) * z*(2-z);
}

double blendf( const DROPS::Point3DCL & p){
    return blendf(p[0],p[1],p[2]);
}

double lapblendf( double x, double y, double z)
{
    return -2 * y*(2-y) * z*(2-z)
        -   2 * x*(2-x) * z*(2-z)
        -   2 * x*(2-x) * y*(2-y);
}

double lapblendf( const DROPS::Point3DCL & p){
    return lapblendf(p[0],p[1],p[2]);
}

DROPS::Point3DCL gradblendf( double x, double y, double z)
{
    return DROPS::MakePoint3D(2 *(1-x) * y*(2-y) * z*(2-z), 
                       2 *(1-y) * x*(2-x) * z*(2-z), 
                       2 *(1-z) * x*(2-x) * y*(2-y));
}

DROPS::Point3DCL gradblendf( const DROPS::Point3DCL & p){
    return gradblendf(p[0],p[1],p[2]);
}



double testcase6_shift( double v, double t)
{
    extern DROPS::ParamCL P;
    static double r (P.get<double>("TestCase6.r"));
    static bool sinw (P.get<int>("TestCase6.sinmove",0));
    if (sinw)
    {
        return 0.25/M_PI * v * sin(2.0*M_PI*t);
    }
    else
    {
        if (t<r)
            return v*t;
        else
            return v*(2*r-t);
    }
}

DROPS::Point3DCL TestCase6_vel (const DROPS::Point3DCL& ,  double t )
{  
    extern DROPS::ParamCL P;
    static double v (P.get<double>("TestCase6.v"));
    static DROPS::Point3DCL vel( DROPS::MakePoint3D(v,0,0));
    static double r (P.get<double>("TestCase6.r"));
    static bool sinw (P.get<int>("TestCase6.sinmove",0));

    if (sinw)
    {
        return DROPS::MakePoint3D(0.5 * v * cos(2*M_PI*t),0,0);
    }
    else
    {
        if (t > r) 
            vel = DROPS::MakePoint3D(-v,0,0);
    }
    return vel;
}

void TestCase6_push()
{
    static bool pushed = false;
    if (!pushed)
    {
        extern DROPS::ParamCL P;
        static double hpos (P.get<double>("Transp.HPos"));
        static double apos (P.get<double>("Transp.DiffPos"));
        static double hneg (P.get<double>("Transp.HNeg"));
        static double aneg (P.get<double>("Transp.DiffNeg"));

        static double R (P.get<DROPS::Point3DCL>("Levelset.RadDrop")[0]);

        const double b = -M_PI*sin(M_PI*R)/(2.0*R)*apos/aneg;
        P.put_if_unset<double>("TestCase6.b",b);

        const double a = cos(M_PI*R)*hpos/hneg-b*R*R;
        P.put_if_unset<double>("TestCase6.a",a);

        std::cout << " a = " << a << std::endl;
        std::cout << " b = " << b << std::endl;
        pushed = true;
    }
}

double TestCase6_lset( const DROPS::Point3DCL& p ,  double t )
{
    extern DROPS::ParamCL P;
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    static DROPS::Point3DCL rad (P.get<DROPS::Point3DCL>("Levelset.RadDrop"));
    static double v (P.get<double>("TestCase6.v"));
    DROPS::Point3DCL p2(DROPS::MakePoint3D(pos[0]+testcase6_shift(v,t),pos[1],pos[2]));
    // DROPS::Point3DCL p3(DROPS::MakePoint3D(pos[0]+2.0+testcase6_shift(v,t),pos[1],pos[2]));
    // DROPS::Point3DCL p4(DROPS::MakePoint3D(pos[0]-2.0+testcase6_shift(v,t),pos[1],pos[2]));
    DROPS::Point3DCL tmp ((p-p2)/rad);
    // DROPS::Point3DCL tmp2 ((p-p3)/rad);
    // DROPS::Point3DCL tmp3 ((p-p4)/rad);
    // return std::min(std::min(tmp.norm(),tmp2.norm()),tmp3.norm()) - 1.0;
    return tmp.norm() - 1.0;
}



double TestCase6_rhs_positive (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    // static double hpos (P.get<double>("Transp.HPos"));
    static double apos (P.get<double>("Transp.DiffPos"));
    // static double hneg (P.get<double>("Transp.HNeg"));
    // static double aneg (P.get<double>("Transp.DiffNeg"));
    // static double k (P.get<double>("TestCase6.k"));
    static double C (P.get<double>("TestCase6.C"));
    static double v (P.get<double>("TestCase6.v"));

    DROPS::Point3DCL shiftedpos;
    shiftedpos = pos;
    shiftedpos[0] += testcase6_shift(v,t);
    DROPS::Point3DCL diff = p - shiftedpos;
    const double absx = diff.norm();
    const double timefactor = gtc6(t);

    const double diffpart1 = M_PI*M_PI*cos(M_PI*absx) + sin(M_PI*absx)*M_PI*2.0/absx;
    const double diffpart2 = apos * diffpart1 * timefactor;
    const double blendval = blendf(p);
    const double u = cos(M_PI*absx);
    const double timeder = u * ddtgtc6(t);

    DROPS::Point3DCL gradu = diff;
    gradu *= -M_PI * sin(M_PI * absx) / absx;
    DROPS::Point3DCL gradblend = gradblendf(p);
    const double blendpart = - apos * (2*inner_prod(gradblend,gradu) + lapblendf(p)*u) + u*v*gradblend[0];
        
    //convection missing
    return C*((timeder + diffpart2)* blendval + timefactor*blendpart);
}


double TestCase6_rhs_negative (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    // static double hpos (P.get<double>("Transp.HPos"));
    // static double apos (P.get<double>("Transp.DiffPos"));
    // static double hneg (P.get<double>("Transp.HNeg"));
    static double aneg (P.get<double>("Transp.DiffNeg"));
    // static double k (P.get<double>("TestCase6.k"));
    TestCase6_push();
    static double a (P.get<double>("TestCase6.a"));
    static double b (P.get<double>("TestCase6.b"));
    static double C (P.get<double>("TestCase6.C"));
    static double v (P.get<double>("TestCase6.v"));
    DROPS::Point3DCL shiftedpos;
    shiftedpos = pos;
    shiftedpos[0] += testcase6_shift(v,t);
    DROPS::Point3DCL diff = p - shiftedpos;
    const double absx = diff.norm();
    const double timefactor = gtc6(t);

    const double diffpart = -6*b * aneg * timefactor;
    const double blendval = blendf(p);
    const double u = a+b*absx*absx;
    const double timeder = u * ddtgtc6(t);

    DROPS::Point3DCL gradu = diff;
    gradu *= 2*b;
    DROPS::Point3DCL gradblend = gradblendf(p);
    const double blendpart = - aneg * (2*inner_prod(gradblend,gradu) + lapblendf(p)*u) + u*v*gradblend[0];

    //convection missing
    return C*((timeder + diffpart)*blendval + timefactor*blendpart);
}

double TestCase6_sol_positive (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    // static double k (P.get<double>("TestCase6.k"));
    static double C (P.get<double>("TestCase6.C"));
    static double v (P.get<double>("TestCase6.v"));
    DROPS::Point3DCL shiftedpos;
    shiftedpos = pos;
    shiftedpos[0] += testcase6_shift(v,t);
    DROPS::Point3DCL diff = p - shiftedpos;
    const double absx = diff.norm();
    const double timefactor = gtc6(t);

    return C*timefactor*cos(M_PI*absx)*blendf(p);
}

double TestCase6_sol_negative (const DROPS::Point3DCL& p ,  double t)
{  
    extern DROPS::ParamCL P;
    
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    // static double k (P.get<double>("TestCase6.k"));
    TestCase6_push();
    static double a (P.get<double>("TestCase6.a"));
    static double b (P.get<double>("TestCase6.b"));
    static double C (P.get<double>("TestCase6.C"));
    static double v (P.get<double>("TestCase6.v"));
    DROPS::Point3DCL shiftedpos;
    shiftedpos = pos;
    shiftedpos[0] += testcase6_shift(v,t);
    DROPS::Point3DCL diff = p - shiftedpos;
    const double absx = diff.norm();
    const double timefactor = gtc6(t);

    // static double lasttime = -1;
    // if (t != lasttime){
    //     std::cout << " timefactor = " << timefactor << std::endl;
    //     lasttime = t;
    // }
    return C*timefactor*(a+b*absx*absx)*blendf(p);
}

static DROPS::RegisterScalarFunction regsca_testcase6_RhsPos("testcase6_rhs_pos", TestCase6_rhs_positive);
static DROPS::RegisterScalarFunction regsca_testcase6_RhsNeg("testcase6_rhs_neg", TestCase6_rhs_negative);
static DROPS::RegisterScalarFunction regsca_testcase6_SolPos("testcase6_sol_pos", TestCase6_sol_positive);
static DROPS::RegisterScalarFunction regsca_testcase6_SolNeg("testcase6_sol_neg", TestCase6_sol_negative);
static DROPS::RegisterScalarFunction regsca_testcase6_lset("testcase6_lset", TestCase6_lset);
static DROPS::RegisterVectorFunction regvec_testcase6_velocity("testcase6_vel", TestCase6_vel);

double qr_func(double r)
{
    const double fac = pow(5.0,10)/pow(3.0,13) * 4.0/3.0;
    //const double fac = pow(5.0,4)/pow(2.0,7);

    if (r>1)
        return 0.0;
    else if (r>0.9)
        return 0.1*M_PI * 1e4 * ((0.8-r)*(0.8-r)*(1-r)*(1-r));
    else
        return 0.1*M_PI*(1.0+  fac * (r*r*(0.9-r)*(0.9-r)*(0.9-r)));

    // the snake:
    // else
    //     return 0.1*M_PI*((1-r*r)*(1-r*r));

    // else if (r>0.8)
    //     return 0.1*M_PI * 1e4/16.0 * ((0.6-r)*(0.6-r)*(1-r)*(1-r));
    // else
    //     return 0.1*M_PI * (1.0+  fac * (r*r*(0.8-r)*(0.8-r)));

    // if (r>1)
    //     return 0.1*M_PI*0.0;
    // else
    // {
    //     return 0.1*M_PI*(0.0+3125.0/108.0 * 1.0/8.0 * (r*r*(1-r)*(1-r)*(1-r)));
    // }

}

DROPS::Point3DCL TestCase7_vel (const DROPS::Point3DCL& p ,  double )
{  
    static DROPS::Point3DCL rotcent (DROPS::MakePoint3D(1, 1, 0 ));
    DROPS::Point3DCL xdiff (p-rotcent);
    const double r = sqrt(xdiff[0]*xdiff[0]+xdiff[1]*xdiff[1]);
    const double qr = qr_func(r);
    // const double qr = r > 1 ? 0 : r * (1-r);

    if (false){
        std::cout << " p = " << p << std::endl;
        std::cout << " r = " << r << std::endl;
        std::cout << " qr = " << qr << std::endl;
        std::cout << " - qr *  DROPS::MakePoint3D(p[1]-rotcent[1],p[0]-rotcent[0],0) = " << - qr *  DROPS::MakePoint3D(p[1]-rotcent[1],p[0]-rotcent[0],0) << std::endl;
        getchar();
    }
    if (r>1.0) 
        return DROPS::MakePoint3D(0,0,0);
    else
        return qr *  DROPS::MakePoint3D(p[1]-rotcent[1],rotcent[0]-p[0],0);
}

// void TestCase7_push()
// {
//     static bool pushed = false;
//     if (!pushed)
//     {
//         extern DROPS::ParamCL P;
//         static double hpos (P.get<double>("Transp.HPos"));
//         static double apos (P.get<double>("Transp.DiffPos"));
//         static double hneg (P.get<double>("Transp.HNeg"));
//         static double aneg (P.get<double>("Transp.DiffNeg"));

//         static double R (P.get<DROPS::Point3DCL>("Levelset.RadDrop")[0]);

//         const double b = -M_PI*sin(M_PI*R)/(2.0*R)*apos/aneg;
//         P.put_if_unset<double>("TestCase7.b",b);

//         const double a = cos(M_PI*R)*hpos/hneg-b*R*R;
//         P.put_if_unset<double>("TestCase7.a",a);

//         std::cout << " a = " << a << std::endl;
//         std::cout << " b = " << b << std::endl;
//         pushed = true;
//     }
// }

double TestCase7_lset( const DROPS::Point3DCL& p ,  double t )
{
    extern DROPS::ParamCL P;
    static DROPS::Point3DCL pos (P.get<DROPS::Point3DCL>("Levelset.PosDrop"));
    static DROPS::Point3DCL rad (P.get<DROPS::Point3DCL>("Levelset.RadDrop"));

    static DROPS::Point3DCL rotcent (DROPS::MakePoint3D(1, 1, 0 ));

    DROPS::Point3DCL xdiff (p-rotcent);
    const double r = sqrt(xdiff[0]*xdiff[0]+xdiff[1]*xdiff[1]);
    const double qr = qr_func(r);
    // const double qr = r > 1 ? 0 : r * (1-r);

    double alphapdt = atan2(xdiff[1],xdiff[0]);

    double alpha = alphapdt + t * qr;

    DROPS::Point3DCL oldpos ( DROPS::MakePoint3D(rotcent[0] + r * cos(alpha), rotcent[1] + r * sin(alpha), p[2] ));

    DROPS::Point3DCL olddiff ((oldpos - pos)/rad);
    if (false && olddiff.norm() < 1)
    {
        std::cout << " t = " << t << std::endl;
        std::cout << " alphapdt = " << alphapdt << std::endl;
        std::cout << " alpha = " << alpha << std::endl;
        std::cout << " p = " << p << std::endl;
        std::cout << " pos = " << pos << std::endl;
        std::cout << " rad = " << rad << std::endl;
        std::cout << " r = " << r << std::endl;
        std::cout << " oldpos = " << oldpos << std::endl;
        std::cout << " olddiff = " << olddiff << std::endl;
        std::cout << " olddiff.norm() - 1.0 = " << olddiff.norm() - 1.0 << std::endl;

        getchar();
    }
    return olddiff.norm() - 1.0;
}



double TestCase7_rhs_positive (const DROPS::Point3DCL& ,  double )
{
    // extern DROPS::ParamCL P;

    return 0;
}


double TestCase7_rhs_negative (const DROPS::Point3DCL& ,  double )
{
    //extern DROPS::ParamCL P;

    return 0;
}

static DROPS::RegisterScalarFunction regsca_testcase7_RhsPos("testcase7_rhs_pos", TestCase7_rhs_positive);
static DROPS::RegisterScalarFunction regsca_testcase7_RhsNeg("testcase7_rhs_neg", TestCase7_rhs_negative);
// static DROPS::RegisterScalarFunction regsca_testcase7_SolPos("testcase7_sol_pos", TestCase7_sol_positive);
// static DROPS::RegisterScalarFunction regsca_testcase7_SolNeg("testcase7_sol_neg", TestCase7_sol_negative);
static DROPS::RegisterScalarFunction regsca_testcase7_lset("testcase7_lset", TestCase7_lset);
static DROPS::RegisterVectorFunction regvec_testcase7_velocity("testcase7_vel", TestCase7_vel);





    bool periodic_3sides( const DROPS::Point3DCL& p, const DROPS::Point3DCL& q)
    {
        static DROPS::Point3DCL dx= MeshSize();
        const DROPS::Point3DCL d= fabs(p-q), L= fabs(dx);
        
        bool matches = true;
        for (int i = 0; i < 3; ++i)
            if (!(d[i] < 1e-12 || std::abs(d[i]-L[i]) < 1e-12))
                matches = false;
        return matches;
    }


    template<int A, int B>
    bool periodic_2sides( const DROPS::Point3DCL& p, const DROPS::Point3DCL& q)
    {
        static DROPS::Point3DCL dx= MeshSize();
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
        static DROPS::Point3DCL dx= MeshSize();

        const int B = (A+1)%3;
        const int D = (B+1)%3;
        const DROPS::Point3DCL d= fabs(p-q), L= fabs(dx);

        return (d[B] + d[D] < 1e-12 && std::abs( d[A] - L[A]) < 1e-12);
    }

    //========================================================================
    //        Registration of the function(s) in the func-container
    //========================================================================
    static DROPS::RegisterMatchingFunction regmatch3("periodicxyz", periodic_3sides);
    static DROPS::RegisterMatchingFunction regmatch2_xy("periodicxy", periodic_2sides<0,1>);
    static DROPS::RegisterMatchingFunction regmatch2_xz("periodicxz", periodic_2sides<0,2>);
    static DROPS::RegisterMatchingFunction regmatch2_yz("periodicyz", periodic_2sides<1,2>);
    static DROPS::RegisterMatchingFunction regmatch1_x("periodicx", periodic_1side<0>);
    static DROPS::RegisterMatchingFunction regmatch1_y("periodicy", periodic_1side<1>);
    static DROPS::RegisterMatchingFunction regmatch1_z("periodicz", periodic_1side<2>);

