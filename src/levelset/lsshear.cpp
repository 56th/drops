/// \file lsshear.cpp
/// \brief drop in shear flow
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt;  SC RWTH Aachen:

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
#include "navstokes/instatnavstokes2phase.h"
#include "num/krylovsolver.h"
#include "num/precond.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "levelset/coupling.h"
#include "levelset/surfacetension.h"
#include <fstream>

#include "misc/scalarFunctions.cpp"
#include "misc/vectorFunctions.cpp"

const double      delta_t= 0.01;
const DROPS::Uint num_steps= 50;
const int         FPsteps= -1;
DROPS::ParamCL P;

// du/dt - q*u - nu*laplace u + Dp = f - okn
//                          -div u = 0
//                               u = u0, t=t0

// Randdaten: x=0, x=1, y=0, y=1:  Dirichlet 0
//            z=0 und x<0.5        Neumann   0   (aus Impl.gruenden: Dir.)
//            z=0 und x>0.5        Inflow Dirichlet  parabol.
//            z=1 und x<0.5        Inflow Dirichlet  parabol.
//            z=1 und x>0.5        Neumann   0   (aus Impl.gruenden: Dir.)

//DROPS::Point3DCL ZeroVel( const DROPS::Point3DCL&, double) { return DROPS::Point3DCL(0.); }

// Tropfendaten:
DROPS::Point3DCL Mitte(0.5);
double           Radius= 0.2;

DROPS::SVectorCL<3> Parabol( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.);
    if (p[0]<0.5)
        ret[2]= 4*p[0]*(p[0]-0.5);
    else
        ret[2]= 4*(1-p[0])*(p[0]-0.5);
    return ret;
}

double DistanceFct( const DROPS::Point3DCL& p, double)
{
    return (Mitte-p).norm()-Radius;
}

double sigma;
double sigmaf (const DROPS::Point3DCL&, double) { return sigma; }

namespace DROPS // for Strategy
{

template<class StokesProblemT>
void Strategy( StokesProblemT& Stokes, const BndDataCL<>& lsbnd)
// flow control
{
    MultiGridCL& MG= Stokes.GetMG();
    SurfaceTensionCL sf( sigmaf);
    LevelsetP2CL & lset( * LevelsetP2CL::Create( MG, lsbnd, sf, false, 0.1) );

    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;
    VelVecDescCL* v= &Stokes.v;
    VecDescCL*    p= &Stokes.p;
    VelVecDescCL* b= &Stokes.b;
    VecDescCL* c= &Stokes.c;
    VelVecDescCL cpl_M;
    MLMatDescCL* A= &Stokes.A;
    MLMatDescCL* B= &Stokes.B;
    MLMatDescCL* C= &Stokes.C;
    MLMatDescCL* M= &Stokes.M;
    MLMatDescCL* prM = &Stokes.prM;

    TimerCL time;
    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr ( MG.GetLastLevel(), pidx);
    lset.CreateNumbering     ( MG.GetLastLevel());
    lset.Init( DistanceFct);

    MG.SizeInfo( std::cout);
    b->SetIdx( vidx);
    c->SetIdx( pidx);
    cpl_M.SetIdx( vidx);
    p->SetIdx( pidx);
    v->SetIdx( vidx);
    std::cout << "Anzahl der Druck-Unbekannten: " << p->Data.size() << std::endl;
    std::cout << "Anzahl der Geschwindigkeitsunbekannten: " << v->Data.size() << std::endl;
    A->Reset();
    B->Reset();
    M->Reset();
    A->SetIdx(vidx, vidx);
    B->SetIdx(pidx, vidx);
    M->SetIdx(vidx, vidx);
    Stokes.N.SetIdx(vidx, vidx);
    prM->SetIdx( pidx, pidx);
    time.Reset();
    time.Start();
    Stokes.SetupSystem1( A, M, b, b, &cpl_M, lset, Stokes.v.t);
    Stokes.SetupSystem2( B, C, c, lset, Stokes.v.t);
    Stokes.SetupPrMass( prM, lset);
    time.Stop();
    std::cout << time.GetTime() << " seconds for setting up all systems!" << std::endl;

    Stokes.InitVel( v, ZeroVel);
    lset.SetupSystem( Stokes.GetVelSolution(), delta_t);

    time.Reset();

    double outer_tol;
    std::cout << "tol = "; std::cin >> outer_tol;
    GSPcCL pc;
    PCGSolverCL<GSPcCL> PCGsolver( pc, 200, 1e-2, true);
    typedef SolverAsPreCL<PCGSolverCL<GSPcCL> > PCGPcT;
    PCGPcT apc( PCGsolver);
    ISBBTPreCL bbtispc( &Stokes.B.Data.GetFinest(), &Stokes.prM.Data.GetFinest(), &Stokes.M.Data.GetFinest(), Stokes.pr_idx.GetFinest(), 0.0, 1.0, 1e-4, 1e-4);
    InexactUzawaCL<PCGPcT, ISBBTPreCL, APC_OTHER> inexactuzawasolver( apc, bbtispc, 200, outer_tol, 0.6, 50);

    {
        std::cout << "Computing initial velocity..." << std::endl;

        inexactuzawasolver.Solve( A->Data, B->Data, C->Data, v->Data, p->Data, b->Data, c->Data, v->RowIdx->GetEx(),  p->RowIdx->GetEx());
    }

    // Initialize Ensight6 output
    std::string ensf( "ensight/shear");
    Ensight6OutCL ensight( "shear.case", num_steps + 1);
    ensight.Register( make_Ensight6Geom      ( MG, MG.GetLastLevel(),   "shear flow field", ensf + ".geo"));
    ensight.Register( make_Ensight6Scalar    ( lset.GetSolution(),      "Levelset",         ensf + ".scl", true));
    ensight.Register( make_Ensight6Scalar    ( Stokes.GetPrSolution(),  "Pressure",         ensf + ".pr",  true));
    ensight.Register( make_Ensight6Vector    ( Stokes.GetVelSolution(), "Velocity",         ensf + ".vel", true));
    ensight.Write();
    typedef GMResSolverCL<SSORPcCL> LSetSolver;
    SSORPcCL ssorpc;
    LSetSolver gm( ssorpc, 100, 1000, 1e-7);
    LevelsetModifyCL lsetmod( 0, 0, 0, 0);
    typedef NSSolverBaseCL<StokesProblemT> SolverT;
    SolverT dummyFP( Stokes, inexactuzawasolver);
    LinThetaScheme2PhaseCL<LSetSolver>
        cpl( Stokes, lset, dummyFP, gm, lsetmod, /*theta*/ 0.5, 0.5, /*nonlinear*/ 0.);
    cpl.SetTimeStep( delta_t);

    bbtispc.SetWeights(1.0/delta_t, 0.5);

    for (Uint step= 1; step<=num_steps; ++step)
    {
        std::cout << "============= Schritt " << step << ":\n";
        cpl.DoStep( FPsteps);
        ensight.Write( step*delta_t);
    }
    std::cout << "Iterationen: " << inexactuzawasolver.GetIter()
              << "\tNorm des Res.: " << inexactuzawasolver.GetResid() << std::endl;

    std::cout << std::endl;

    delete &lset;
}

} // end of namespace DROPS


int main (int argc, char** argv)
{
  try
  {
    if (argc!=3)
    {
        std::cout << "You have to specify three parameters:\n\tlsshear <num_subdiv> <surf.tension>" << std::endl;
        return 1;
    }

    int sub_div= std::atoi(argv[1]);
    sigma= std::atof(argv[2]);
    std::cout << "sub divisions:   " << sub_div << std::endl;
    std::cout << "surface tension: " << sigma << std::endl;
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.;

    typedef DROPS::InstatNavierStokes2PhaseP2P1CL MyStokesCL;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, sub_div, sub_div, sub_div);

    const bool IsNeumann[6]=
        {false, false, false, false, false, false};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &ZeroVel, &ZeroVel, &ZeroVel, &ZeroVel,  &Parabol, &Parabol };
    // parabol. Einstroembedingungen bei z=0 und z=1

    const DROPS::BndCondT bcls[6]= { DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC };
    const DROPS::LsetBndDataCL::bnd_val_fun bfunls[6]= { 0,0,0,0,0,0};
    DROPS::LsetBndDataCL lsbnd( 6, bcls, bfunls);

    DROPS::TwoPhaseFlowCoeffCL coeff( 1, 1, 1, 1, 0, DROPS::Point3DCL(0.));

    MyStokesCL prob(brick, coeff, DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
    DROPS::MultiGridCL& mg = prob.GetMG();
    Strategy(prob, lsbnd);
    std::cout << DROPS::SanityMGOutCL(mg) << std::endl;
    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    std::cout << "pressure min/max: "<<min<<", "<<max<<std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
