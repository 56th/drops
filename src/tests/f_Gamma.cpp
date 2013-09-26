/// \file f_Gamma.cpp
/// \brief test surface force term
/// \author LNM RWTH Aachen: Sven Gross, Joerg Peters, Volker Reichelt; SC RWTH Aachen:

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

#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "stokes/instatstokes2phase.h"
#include "num/krylovsolver.h"
#include "num/precond.h"
#include "levelset/coupling.h"
#include "misc/params.h"
#include "levelset/surfacetension.h"
#include "misc/dynamicload.h"
#include <fstream>
#include <sstream>


DROPS::ParamCL P;

int FctCode=9;

class ZeroFlowCL
{
// \Omega_1 = Tropfen,    \Omega_2 = umgebendes Fluid
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
        { DROPS::Point3DCL ret(0.0); return ret; }
    const DROPS::SmoothedJumpCL rho, mu;
    const DROPS::Point3DCL g;

    ZeroFlowCL( const DROPS::ParamCL& P)
      : rho( DROPS::JumpCL( 1, 1), DROPS::H_sm, P.get<double>("Mat.SmoothZone")),
         mu( DROPS::JumpCL( 1, 1),   DROPS::H_sm, P.get<double>("Mat.SmoothZone")),
        g()    {}
};

double sigmaf( const DROPS::Point3DCL&, double) { return -1.0; }

/*
double DistanceFct( const DROPS::Point3DCL& p)
{ // cube
    const DROPS::Point3DCL d= C.Mitte-p;
    double max=-1;
    for (int i=0; i<3; ++i)
        if (fabs(d[i])>max) max= fabs(d[i]);
    return max-C.Radius;
}

double DistanceFct( const DROPS::Point3DCL& p)
{ // plane perpendicular to n=PosDrop with distance Radius from origin.
    return inner_prod( C.Mitte/norm(C.Mitte), p) - C.Radius;
}

double DistanceFct( const DROPS::Point3DCL& p)
{ // ball
    const DROPS::Point3DCL d= C.Mitte-p;
    return d.norm()-C.Radius;
}
*/

double DistanceFct( const DROPS::Point3DCL& p, double)
{ // ball
    const DROPS::Point3DCL d= P.get<DROPS::Point3DCL>("Exp.PosDrop")-p;
//    return d.norm_sq()-C.Radius*C.Radius; // exakte Darstellung mit P2-FE, aber keine Abstandsfunktion
    return d.norm()-P.get<DROPS::Point3DCL>("Exp.RadDrop")[0];
}

DROPS::Point3DCL Null( const DROPS::Point3DCL&, double)
{ return DROPS::Point3DCL(0.); }

void WriteFct( std::ostream& os)
{
    const int idx= FctCode%10,
        fct= FctCode/10;
    os << "(";
    std::string term;
    switch (fct%10)
    {
        case 1:
            term= "x"; break;
        case 2:
            term= "y"; break;
        case 3:
            term= "z"; break;
        case 0: // 1
            term= "1"; break;
        default:
            std::cout << "Fehler in WriteFct, fct = " << fct << std::endl;
    }
    if (fct>=9)
        switch (fct/10)
        {
            case 1:
                term+= "x"; break;
            case 2:
                term+= "y"; break;
            case 3:
                term+= "z"; break;
            case 0: // 1
                term+= "1"; break;
            default:
                std::cout << "Fehler in WriteFct (quadratisch), fct = " << fct << std::endl;
        }
    for (int i=0; i<3; ++i)
        os << ' ' << (i==idx-1 ? term : "0");
    os << " )";
}


DROPS::Point3DCL TestFct( const DROPS::Point3DCL& p, double)
{
    const int idx= FctCode%10,
        fct= FctCode/10,
        fct1= fct%10, fct2= fct/10;
    double val= 0;
    switch (fct1)
    {
        case 1: case 2: case 3: // x/y/z
            val= p[fct1-1]; break;
        case 0: // 1
            val= 1; break;
        default:
            std::cout << "Fehler in TestFct, fct = " << fct << std::endl;
    }
    if (fct>=9)
        switch (fct2)
        {
            case 1: case 2: case 3: // x/y/z
                val*= p[fct2-1]; break;
            case 0: // 1
                val*= 1; break;
            default:
                std::cout << "Fehler in TestFct (quadratisch), fct = " << fct << std::endl;
        }

    DROPS::Point3DCL ret;
    if (idx<1 || idx>3) std::cout << "Fehler in TestFct, idx = " << idx << std::endl;
    ret[idx-1]= val;
    return ret;
}

namespace DROPS // for Strategy
{

void ApplyToTestFct( InstatStokes2PhaseP2P1CL& Stokes, const LsetBndDataCL& lsbnd)
// program for testing the approximation order of
// the discretization of the surface force term f_Gamma
// using simple trial functions:
//     | f_Gamma(v) - f_{Gamma_h}(v) |  =   O(h^p)
{
    MultiGridCL& MG= Stokes.GetMG();
    const double curv= 2/P.get<DROPS::Point3DCL>("Exp.RadDrop")[0];
    SurfaceTensionCL sf( sigmaf, 0);

    LevelsetP2CL & lset( * LevelsetP2CL::Create( MG, lsbnd, sf, P.get_child("Levelset")) );

//    lset.SetSurfaceForce( SF_Const);

    MLIdxDescCL* lidx= &lset.idx;
    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;

    lset.CreateNumbering( MG.GetLastLevel(), lidx);

    lset.Phi.SetIdx( lidx);
    lset.Init( DistanceFct);

    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, NULL, &lset);

    VecDescCL f_Gamma( vidx), v( vidx);
    MG.SizeInfo( std::cout);
    Stokes.b.SetIdx( vidx);
    Stokes.c.SetIdx( pidx);
    Stokes.p.SetIdx( pidx);
    Stokes.v.SetIdx( vidx);

    std::cout << Stokes.p.Data.size() << " pressure unknowns,\n";
    std::cout << Stokes.v.Data.size() << " velocity unknowns,\n";
    std::cout << lset.Phi.Data.size() << " levelset unknowns.\n";

    const double Vol= 4./3*M_PI*std::pow( P.get<DROPS::Point3DCL>("Exp.RadDrop")[0], 3);
    std::cout << "Volumen = " << Vol << "\tKruemmung = " << curv << "\n\n";
    typedef std::map<int,double> LsgMapT;
    LsgMapT reflsg;
    reflsg[1]= 0;
//    reflsg[2]= 0;
//    reflsg[3]= 0;
    reflsg[11]= curv*Vol;
    reflsg[12]= 0;
//    reflsg[13]= 0;
//    reflsg[21]= 0;
    reflsg[22]= curv*Vol;
//    reflsg[23]= 0;
//    reflsg[31]= 0;
//    reflsg[32]= 0;
    reflsg[33]= curv*Vol;
    reflsg[111]= curv*2*Vol*P.get<DROPS::Point3DCL>("Exp.PosDrop")[0];
    reflsg[121]= curv*Vol*P.get<DROPS::Point3DCL>("Exp.PosDrop")[1];
    reflsg[131]= curv*Vol*P.get<DROPS::Point3DCL>("Exp.PosDrop")[2];
//    reflsg[231]= 0;
//    reflsg[221]= 0;
//    reflsg[331]= 0;

    f_Gamma.Clear( Stokes.v.t);
    lset.AccumulateBndIntegral( f_Gamma);

    std::vector<double> errVec, refVec;
    for (LsgMapT::iterator it= reflsg.begin(), end= reflsg.end(); it!=end; ++it)
    {
        FctCode= it->first;
        const double Lsg= it->second;
        Stokes.InitVel( &v, TestFct);
        const double fv= dot( f_Gamma.Data, v.Data);
        const double err= std::abs(fv - Lsg);
        WriteFct( std::cout);
        std::cout << "\nfct " << FctCode << ":\tval = " << fv << "\tref = " << Lsg << "\nerror = " << err << "\n\n";
        errVec.push_back( err);
        refVec.push_back( Lsg);
    }

    std::cout << "\n\n\"f\",";
    for (LsgMapT::iterator it= reflsg.begin(), end= reflsg.end(); it!=end; ++it)
    { std::cout << "\""; FctCode= it->first; WriteFct( std::cout); std::cout << "\","; }

    std::cout << "\n\n\"Referenz\",";
    for (size_t i=0; i<refVec.size(); ++i)
        std::cout << refVec[i] << ",\t";

    std::cout << "\n\n" << P.get<int>("AdaptRef.FinestLevel") << ",\t";
    for (size_t i=0; i<errVec.size(); ++i)
        std::cout << errVec[i] << ",\t";
    std::cout << "\n\n";
    delete &lset;
}

void Compare_LaplBeltramiSF_ConstSF( InstatStokes2PhaseP2P1CL& Stokes, const LsetBndDataCL& lsbnd)
// computation of order of LB discretization, cf. paper
// S. Gross, A. Reusken: Finite element discretization error analysis of a surface tension force in two-phase incompressible flows,
// SIAM J. Numer. Anal. 45, 1679--1700 (2007)

{
    MultiGridCL& MG= Stokes.GetMG();
    // Levelset-Disc.: Crank-Nicholson
    const double curv= 2/P.get<DROPS::Point3DCL>("Exp.RadDrop")[0];
    SurfaceTensionCL sf( sigmaf, 0);
    LevelsetP2CL & lset( * LevelsetP2CL::Create( MG, lsbnd, sf, P.get_child("Levelset")) );


    MLIdxDescCL* lidx= &lset.idx;
    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;

    lset.CreateNumbering( MG.GetLastLevel(), lidx);

    lset.Phi.SetIdx( lidx);
    lset.Init( DistanceFct);

    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, NULL, &lset);

    VecDescCL f_Const( vidx), f_LaplBeltrami( vidx), v( vidx);
    MG.SizeInfo( std::cout);
    Stokes.b.SetIdx( vidx);
    Stokes.c.SetIdx( pidx);
    Stokes.p.SetIdx( pidx);
    Stokes.v.SetIdx( vidx);
    Stokes.A.SetIdx( vidx, vidx);
    Stokes.M.SetIdx( vidx, vidx);

    std::cout << Stokes.p.Data.size() << " pressure unknowns,\n";
    std::cout << Stokes.v.Data.size() << " velocity unknowns,\n";
    std::cout << lset.Phi.Data.size() << " levelset unknowns.\n";

    const double Vol= 4./3*M_PI*std::pow( P.get<DROPS::Point3DCL>("Exp.RadDrop")[0], 3);
    std::cout << "Volumen = " << Vol << "\tKruemmung = " << curv << "\n\n";

    Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &Stokes.b, lset, 0.);

    f_LaplBeltrami.Clear( Stokes.v.t);
    lset.SetSurfaceForce( SF_ImprovedLB);
//     lset.SetSurfaceForce( SF_LB);
    lset.AccumulateBndIntegral( f_LaplBeltrami);

    f_Const.Clear( Stokes.v.t);
    lset.SetSurfaceForce( SF_Const);
    lset.AccumulateBndIntegral( f_Const);

    VectorCL d( curv*f_Const.Data - f_LaplBeltrami.Data);
    std::cout << "|d| = \t\t" << norm(d) << std::endl;
    VectorCL A_inv_d( d.size());
    SSORPcCL pc;
    typedef PCGSolverCL<SSORPcCL>     PCG_SsorCL;
    PCG_SsorCL cg( pc, 1000, 1e-18);
    std::cout << "Solving system with stiffness matrix:\t";
    cg.Solve( Stokes.A.Data, A_inv_d, d, vidx->GetEx());
    std::cout << cg.GetIter() << " iter,\tresid = " << cg.GetResid();
    const double sup= std::sqrt(dot( A_inv_d, d));

    std::cout << "\n\nsup |f1(v)-f2(v)|/|v|_1 = \t\t" << sup
              << "\n|A^-1 d| = \t\t" << norm( A_inv_d) << std::endl;
    MLMatrixCL MA;
    MA.LinComb( 1, Stokes.M.Data, 1, Stokes.A.Data);
    VectorCL MA_inv_d( A_inv_d);
    std::cout << "Solving system with MA matrix:\t";
    cg.Solve( MA, MA_inv_d, d, vidx->GetEx());
    std::cout << cg.GetIter() << " iter,\tresid = " << cg.GetResid();
    const double sup2= std::sqrt(dot( MA_inv_d, d));
    std::cout << "\n\nsup |f1(v)-f2(v)|/||v||_1 = \t\t" << sup2
              << "\n|(MA)^-1 d| = " << norm( MA_inv_d) << std::endl;
    delete &lset;
}

} // end of namespace DROPS


void MarkDrop (DROPS::MultiGridCL& mg, int maxLevel= -1)
{
    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
/*            for (int i=0; i<4; ++i)
            {
                const double val= DistanceFct( It->GetVertex(i)->GetCoord() );
                if ( val<C.ref_width && val > -C.ref_width)
                    It->SetRegRefMark();
            }
            const double val= DistanceFct( GetBaryCenter(*It));
            if ( val<C.ref_width && val > -C.ref_width)
                It->SetRegRefMark();
*/
            int neg= 0, zero= 0;
            for (int i=0; i<4; ++i)
            {
                const double val= DistanceFct( It->GetVertex(i)->GetCoord(), 0.);
                neg+= val<0;
                zero+= fabs(val)<1e-8;
            }
            const double val= DistanceFct( GetBaryCenter(*It), 0.);
            neg+= val<0;
            zero+= fabs(val)<1e-8;

            if ( (neg!=0 && neg!=5) || zero) // change of sign or zero in tetra
               It->SetRegRefMark();
    }
}

/// \brief Set Default parameters here s.t. they are initialized.
/// The result can be checked when Param-list is written to the output.
void SetMissingParameters(DROPS::ParamCL& P){
    P.put_if_unset<std::string>("Exp.VolForce", "ZeroVel");
    P.put_if_unset<double>("SurfTens.ShearVisco", 0.0);
    P.put_if_unset<double>("SurfTens.DilatationalVisco", 0.0);
    P.put_if_unset<double>("Mat.DensDrop", 1000);
    P.put_if_unset<double>("Mat.ViscDrop", 0.001);
    P.put_if_unset<double>("Mat.DensFluid", 1000);
    P.put_if_unset<double>("Mat.ViscFluid", 0.001);
    P.put_if_unset<double>("Mat.SmoothZone", 1e-05);
    P.put_if_unset<DROPS::Point3DCL>("Exp.Gravity", DROPS::Point3DCL());
}

int main (int argc, char** argv)
{
  try
  {
    if (argc>2)
    {
        std::cout << "You have to specify at most one parameter:\n\t"
                  << argv[0] << " [<param_file>]" << std::endl;
        return 1;
    }
    std::ifstream param;
    if (argc>1)
        param.open( argv[1]);
    else
        param.open( "f_Gamma.json");
    if (!param)
    {
        std::cout << "error while opening parameter file\n";
        return 1;
    }
    param >> P;
    param.close();

    SetMissingParameters(P);

    std::cout << P << std::endl;

    DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"), P.get<std::vector<std::string> >("General.DynamicLibs") );

    typedef DROPS::InstatStokes2PhaseP2P1CL    MyStokesCL;

    int nx, ny, nz;
    double dx, dy, dz;
    std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
    size_t idx;
    while ((idx= mesh.find_first_of( delim)) != std::string::npos )
        mesh[idx]= ' ';
    std::istringstream brick_info( mesh);
    brick_info >> dx >> dy >> dz >> nx >> ny >> nz;
    if (!brick_info)
        DROPS::DROPSErrCL("error while reading geometry information: " + mesh);
    P.put("Exp.RadInlet", dx/2);
    DROPS::Point3DCL orig, px, py, pz;
    px[0]= dx; py[1]= dy; pz[2]= dz;
    orig= -0.5*(px+py+pz);
    DROPS::BrickBuilderCL builder ( orig, px, py, pz, nx, ny, nz);

    const DROPS::BndCondT bc[6]=
        { DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &Null, &Null, &Null, &Null, &Null, &Null};

    const DROPS::BndCondT bcls[6]= { DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC };
    const DROPS::LsetBndDataCL::bnd_val_fun bfunls[6]= { 0,0,0,0,0,0};
    DROPS::LsetBndDataCL lsbnd( 6, bcls, bfunls);

    MyStokesCL prob(builder, DROPS::TwoPhaseFlowCoeffCL(P), DROPS::StokesBndDataCL( 6, bc, bnd_fun));

    DROPS::MultiGridCL& mg = prob.GetMG();

    for (int i=0; i<P.get<int>("AdaptRef.FinestLevel"); ++i)
    {
        MarkDrop( mg);
        mg.Refine();
    }
    std::cout << DROPS::SanityMGOutCL(mg) << std::endl;
    DROPS::GeomMGOutCL out( mg, -1, false, 0);
    std::ofstream fil("cube.off");
    fil << out;
    fil.close();
    Compare_LaplBeltramiSF_ConstSF( prob, lsbnd);
    return 0;
  }
  catch (DROPS::DROPSErrCL& err) { err.handle(); }
}

