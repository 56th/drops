/// \file st_transp.cpp
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
 * Copyright 2012 LNM/SC RWTH Aachen, Germany
*/

//multigrid
#include "geom/multigrid.h"
#include "geom/builder.h"
//time integration
#include "navstokes/instatnavstokes2phase.h"
#include "stokes/integrTime.h"
//output
#include "out/output.h"
#include "out/ensightOut.h"
#include "out/vtkOut.h"
#include "levelset/coupling.h"
#include "misc/params.h"
#include "levelset/adaptriang.h"
#include "levelset/marking_strategy.h"
#include "levelset/mzelle_hdr.h"
#include "levelset/surfacetension.h"
#include "surfactant/ifacetransp.h"
#include "levelset/twophaseutils.h"
#include "misc/funcmap.h"
#include "misc/progressaccu.h"
#include "misc/dynamicload.h"

#include "num/spacetime_geom.h"
#include "spacetimetransp/spacetime_sol.h"
#include "spacetimetransp/spacetime_setup.h"
#include "spacetimetransp/spacetime_error.h"
#include "spacetimetransp/stxfem.h"
#include "spacetimetransp/indicator.h"

#include "num/stokessolverfactory.h"
#ifdef _PAR
#include "parallel/loadbal.h"
#include "parallel/parmultigrid.h"
#endif
#include <fstream>
#include <sstream>
//#include "num/directsolver.h"
#ifndef _PAR
#include "num/stokespardiso.h"
#ifdef DROPS_PARDISO
#include "num/pardisosolver.h"
#endif
#endif
#include <sys/resource.h>

DROPS::ParamCL P;

typedef DROPS::BndDataCL<> cBndDataCL;
typedef cBndDataCL::bnd_val_fun  c_bnd_val_fun;

namespace DROPS // for Strategy
{

double GetTimeOffset(){
    double timeoffset = 0.0;
    const std::string restartfilename = P.get<std::string>("DomainCond.InitialFile");
    if (P.get<int>("DomainCond.InitialCond") == -1){
        const std::string timefilename = restartfilename + "time";
        std::ifstream f_(timefilename.c_str());
        f_ >> timeoffset;
        std::cout << "used time offset file is " << timefilename << std::endl;
        std::cout << "time offset is " << timeoffset << std::endl;
    }
    return timeoffset;
}

void Strategy( InstatNavierStokes2PhaseP2P1CL& Stokes, LsetBndDataCL& lsetbnddata, AdapTriangCL& adap)
// flow control
{
    // rusage usage;

    // InVecMap & tdvectormap = InVecMap::getInstance();
    InScaMap & scalarmap = InScaMap::getInstance();
    InScaMap & inscamap = DROPS::InScaMap::getInstance();
    MatchMap & matchmap = MatchMap::getInstance();

    // instat_scalar_fun_ptr Reaction = scalarmap["ReactionFct"];
    instat_scalar_fun_ptr Initialcneg = scalarmap[P.get<std::string>("Transp.InitialConcNeg")];
    instat_scalar_fun_ptr Initialcpos = scalarmap[P.get<std::string>("Transp.InitialConcPos")];
    instat_scalar_fun_ptr RhsNeg = scalarmap[P.get<std::string>("Transp.RhsNeg")];
    instat_scalar_fun_ptr RhsPos = scalarmap[P.get<std::string>("Transp.RhsPos")];
    // instat_scalar_fun_ptr SolNeg = scalarmap[P.get<std::string>("Transp.SolNeg")];
    // instat_scalar_fun_ptr SolPos = scalarmap[P.get<std::string>("Transp.SolPos")];

    // instat_scalar_fun_ptr distance = scalarmap[P.get<std::string>("Transp.Levelset")];

    // cBndDataCL *pBnd_pos, *pBnd_neg;

    bool is_periodic = P.get<std::string>("DomainCond.PeriodicMatching", "none") != "none";
    match_fun periodic_match = is_periodic ? matchmap[P.get("DomainCond.PeriodicMatching", std::string("periodicx"))] : 0;

    MultiGridCL& MG= Stokes.GetMG();

    // initialization of surface tension
    // choose a proper model for surface tension coefficient, see levelset/surfacetension.h
    SurfaceTensionCL * sf;
    sf = new SurfaceTensionCL( inscamap[P.get<std::string>("NavStokes.Coeff.SurfTens.VarTensionFunc")]);
    sf->SetInputMethod( Sigma_X);


    LevelsetP2CL & lset( * LevelsetP2CL::Create( MG, lsetbnddata, *sf, P.get_child("Levelset")) );
    LevelsetP2CL & oldlset( * LevelsetP2CL::Create( MG, lsetbnddata, *sf, P.get_child("Levelset")) );

    if (is_periodic) //CL: Anyone a better idea? perDirection from ParameterFile?
    {
        DROPS::Point3DCL dx;
        //hack:
        std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
        size_t idx_;
        while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
            mesh[idx_]= ' ';
        std::istringstream brick_info( mesh);
        brick_info >> dx[0] >> dx[1] >> dx[2] ;
        int n = 0;
        if (P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicx" || P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicy" || P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicz")
            n = 1;
        if (P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicxy" || P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicxz" || P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicyz")
            n = 2;
        LevelsetP2CL::perDirSetT pdir(n);
        if (P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicx") pdir[0][0] = dx[0];
        if (P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicy") pdir[0][1] = dx[1];
        if (P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicz") pdir[0][2] = dx[2];
        if (P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicxy") {pdir[0][0] = dx[0]; pdir[1][1] = dx[1];}
        if (P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicxz") {pdir[0][0] = dx[0]; pdir[1][2] = dx[2];}
        if (P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicyz") {pdir[0][1] = dx[1]; pdir[1][2] = dx[2];}
        if (P.get("DomainCond.PeriodicMatching", std::string("periodicx")) != "periodicx" && P.get("DomainCond.PeriodicMatching", std::string("periodicx")) != "periodicy" && P.get("DomainCond.PeriodicMatching", std::string("periodicx")) != "periodicz" &&
          P.get("DomainCond.PeriodicMatching", std::string("periodicx")) != "periodicxy" && P.get("DomainCond.PeriodicMatching", std::string("periodicx")) != "periodicxz" && P.get("DomainCond.PeriodicMatching", std::string("periodicx")) != "periodicyz"){
            std::cout << "WARNING: could not set periodic directions! Reparametrization can not work correctly now!" << std::endl;
            std::cout << "Press any key to continue" << std::endl; getchar();
        }
        lset.SetPeriodicDirections(&pdir);
        oldlset.SetPeriodicDirections(&pdir);
    }

    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);
    // LevelsetRepairCL oldlsetrepair( oldlset);
    // adap.push_back( &oldlsetrepair);
    VelocityRepairCL velrepair( Stokes);
    adap.push_back( &velrepair);
    PressureRepairCL prrepair( Stokes, lset);
    adap.push_back( &prrepair);

    MLIdxDescCL* lidx= &lset.idx;
    // MLIdxDescCL* oldlidx= &oldlset.idx;
    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;

    if ( StokesSolverFactoryHelperCL().VelMGUsed(P)){
        Stokes.SetNumVelLvl ( Stokes.GetMG().GetNumLevel());
        lset.SetNumLvl(Stokes.GetMG().GetNumLevel());
    }
    if ( StokesSolverFactoryHelperCL().PrMGUsed(P)){
        Stokes.SetNumPrLvl  ( Stokes.GetMG().GetNumLevel());
        lset.SetNumLvl(Stokes.GetMG().GetNumLevel());
    }

    lset.CreateNumbering( MG.GetLastLevel(), lidx, periodic_match);
    lset.Phi.SetIdx( lidx);
    // oldlset.CreateNumbering( MG.GetLastLevel(), lidx, periodic_match);
    oldlset.Phi.SetIdx( lidx);
    PermutationT lset_downwind;
    lset.SetSurfaceForce( SF_ImprovedLBVar); // see levelset/levelset.h
    
    SetInitialLevelsetConditions( lset, MG, P);
    // SetInitialLevelsetConditions( oldlset, MG, P);
    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx, periodic_match);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, periodic_match, &lset);
    PermutationT vel_downwind;
    // For a two-level MG-solver: P2P1 -- P2P1X; comment out the preceding CreateNumberings
//     Stokes.SetNumVelLvl ( 2);
//     Stokes.SetNumPrLvl  ( 2);
//     Stokes.vel_idx.GetCoarsest().CreateNumbering( MG.GetLastLevel(), MG, Stokes.GetBndData().Vel);
//     Stokes.vel_idx.GetFinest().  CreateNumbering( MG.GetLastLevel(), MG, Stokes.GetBndData().Vel);
//     Stokes.pr_idx.GetCoarsest(). GetXidx().SetBound( 1e99);
//     Stokes.pr_idx.GetCoarsest(). CreateNumbering( MG.GetLastLevel(), MG, Stokes.GetBndData().Pr, 0, &lset.Phi);
//     Stokes.pr_idx.GetFinest().   CreateNumbering( MG.GetLastLevel(), MG, Stokes.GetBndData().Pr, 0, &lset.Phi);

    StokesVelBndDataCL::bnd_val_fun ZeroVel = InVecMap::getInstance().find("ZeroVel")->second;
    Stokes.SetIdx();
    Stokes.v.SetIdx  ( vidx);
    Stokes.p.SetIdx  ( pidx);
    Stokes.InitVel( &Stokes.v, ZeroVel);
    SetInitialConditions( Stokes, lset, MG, P);
    SetInitialLevelsetConditions( oldlset, MG, P);
    lset.Phi.t = 0;
    oldlset.Phi.t = 0;
    IteratedDownwindCL navstokes_downwind( P.get_child( "NavStokes.Downwind"));
    if (P.get<int>( "NavStokes.Downwind.Frequency") > 0) {
        if (StokesSolverFactoryHelperCL().VelMGUsed( P))
            throw DROPSErrCL( "Strategy: Multigrid-solver and downwind-numbering cannot be used together. Sorry.\n");
        vel_downwind= Stokes.downwind_numbering( lset, navstokes_downwind);
    }
    IteratedDownwindCL levelset_downwind( P.get_child( "Levelset.Downwind"));
    if (P.get<int>( "Levelset.Downwind.Frequency") > 0)
        lset_downwind= lset.downwind_numbering( Stokes.GetVelSolution(), levelset_downwind);

    DisplayDetailedGeom( MG);
    DisplayUnks(Stokes, lset, MG);

    double Vol = 0;

    if (P.get("Exp.InitialLSet", std::string("Ellipsoid")) == "Ellipsoid"){
        Vol = EllipsoidCL::GetVolume();
        std::cout << "initial volume: " << lset.GetVolume()/Vol << std::endl;
        double dphi= lset.AdjustVolume( Vol, 1e-9);
        std::cout << "initial volume correction is " << dphi << std::endl;
        lset.Phi.Data+= dphi;
        std::cout << "new initial volume: " << lset.GetVolume()/Vol << std::endl;
    }else{
        Vol = lset.GetVolume();
    }

    // TransportP1CL * massTransp = NULL;
    // TransportRepairCL *  transprepair = NULL;

    cBndDataCL *transp_pBnd_pos, *transp_pBnd_neg;

    DROPS::BuildBoundaryData( &MG, transp_pBnd_pos,  P.get<std::string>("Transp.BoundaryType"), P.get<std::string>("Transp.BoundaryFncsPos"));
    DROPS::BuildBoundaryData( &MG, transp_pBnd_neg, P.get<std::string>("Transp.BoundaryType"), P.get<std::string>("Transp.BoundaryFncsNeg"));

    cBndDataCL & transp_Bnd_neg(*transp_pBnd_neg);
    cBndDataCL & transp_Bnd_pos(*transp_pBnd_pos); 



    // if (P.get<int>("Transp.DoTransp"))
    // {
    //     // CL: the following could be moved outside of strategy to some function like
    //     //" InitializeMassTransport(P,MG,Stokes,lset,adap, TransportP1CL * & massTransp,TransportRepairCL * & transprepair)"
    //     static DROPS::BndCondT c_bc[6]= {
    //         DROPS::OutflowBC, DROPS::OutflowBC, DROPS::OutflowBC,
    //         DROPS::OutflowBC, DROPS::OutflowBC, DROPS::OutflowBC
    //     };
    //     static DROPS::BndDataCL<>::bnd_val_fun c_bfun[6]= {0, 0, 0, 0, 0, 0};
    //     static DROPS::BndDataCL<> Bnd_c( 6, c_bc, c_bfun);
    //     double D[2] = {P.get<double>("Transp.DiffPos"), P.get<double>("Transp.DiffNeg")};

    //     massTransp = new TransportP1CL( MG, Bnd_c, Stokes.GetBndData().Vel, P.get<double>("Transp.Theta"),
    //                               D, P.get<double>("Transp.HNeg")/P.get<double>("Transp.HPos"), &Stokes.v, lset,
    //                               P.get<double>("Time.StepSize"), P.get<int>("Transp.Iter"), P.get<double>("Transp.Tol"));

    //     transprepair = new TransportRepairCL(*massTransp, MG);
    //     adap.push_back(transprepair);

    //     MLIdxDescCL* cidx= &massTransp->idx;
    //     massTransp->CreateNumbering( MG.GetLastLevel(), cidx);
    //     massTransp->ct.SetIdx( cidx);
    //     if (P.get<int>("DomainCond.InitialCond") != -1)
    //         massTransp->Init( inscamap["Initialcneg"], inscamap["Initialcpos"]);
    //     else
    //         ReadFEFromFile( massTransp->ct, MG, P.get<std::string>("DomainCond.InitialFile")+"concentrationTransf");

    //     massTransp->Update();
    //     std::cout << massTransp->c.Data.size() << " concentration unknowns,\n";
    // }

    /// \todo rhs beruecksichtigen
    SurfactantcGP1CL surfTransp( MG, Stokes.GetBndData().Vel, P.get<double>("SurfTransp.Theta"), P.get<double>("SurfTransp.Visc"), &Stokes.v, lset.Phi, lset.GetBndData(),
                                 P.get<double>("Time.StepSize"), P.get<int>("SurfTransp.Iter"), P.get<double>("SurfTransp.Tol"), P.get<double>("SurfTransp.OmitBound"));
    InterfaceP1RepairCL surf_repair( MG, lset.Phi, lset.GetBndData(), surfTransp.ic);
    if (P.get("SurfTransp.DoTransp", 0))
    {
        adap.push_back( &surf_repair);
        surfTransp.idx.CreateNumbering( MG.GetLastLevel(), MG, &lset.Phi, &lset.GetBndData());
        std::cout << "Surfactant transport: NumUnknowns: " << surfTransp.idx.NumUnknowns() << std::endl;
        surfTransp.ic.SetIdx( &surfTransp.idx);
        surfTransp.Init( inscamap["surf_sol"]);
    }

    // Stokes-Solver
    StokesSolverFactoryCL<InstatNavierStokes2PhaseP2P1CL> stokessolverfactory(Stokes, P);


    StokesSolverBaseCL* stokessolver;

    if (! P.get("Stokes.DirectSolve", 0))
        stokessolver = stokessolverfactory.CreateStokesSolver();
#ifndef _PAR
    else
        stokessolver = new StokesPardisoSolverCL(); 
#else
    else
        throw DROPSErrCL("no direct solver in parallel");
#endif
//     StokesSolverAsPreCL pc (*stokessolver1, 1);
//     GCRSolverCL<StokesSolverAsPreCL> gcr(pc, C.stk_OuterIter, C.stk_OuterIter, C.stk_OuterTol, /*rel*/ false);
//     BlockMatrixSolverCL<GCRSolverCL<StokesSolverAsPreCL> >* stokessolver =
//             new BlockMatrixSolverCL<GCRSolverCL<StokesSolverAsPreCL> > (gcr);

    // Navier-Stokes-Solver
    NSSolverBaseCL<InstatNavierStokes2PhaseP2P1CL>* navstokessolver = 0;
    if (P.get<double>("NavStokes.Nonlinear")==0.0)
        navstokessolver = new NSSolverBaseCL<InstatNavierStokes2PhaseP2P1CL>(Stokes, *stokessolver);
    else
        navstokessolver = new AdaptFixedPtDefectCorrCL<InstatNavierStokes2PhaseP2P1CL>(Stokes, *stokessolver, P.get<int>("NavStokes.Iter"), P.get<double>("NavStokes.Tol"), P.get<double>("NavStokes.Reduction"));

    // Level-Set-Solver
#ifndef _PAR
    // SSORPcCL lset_pc;
    typedef GSPcCL LsetPcT;
#else
    typedef JACPcCL LsetPcT;
#endif
    LsetPcT lset_pc;
    GMResSolverCL<LsetPcT>* gm = new GMResSolverCL<LsetPcT>( lset_pc, 200, P.get<int>("Levelset.Iter"), P.get<double>("Levelset.Tol"));

    LevelsetModifyCL lsetmod( P.get<int>("Reparam.Freq"), P.get<int>("Reparam.Method"), P.get<double>("Reparam.MaxGrad"), P.get<double>("Reparam.MinGrad"), P.get<int>("Levelset.VolCorrection"), Vol, is_periodic);

    UpdateProlongationCL<Point3DCL> PVel( Stokes.GetMG(), stokessolverfactory.GetPVel(), &Stokes.vel_idx, &Stokes.vel_idx);
    adap.push_back( &PVel);
    UpdateProlongationCL<double> PPr ( Stokes.GetMG(), stokessolverfactory.GetPPr(), &Stokes.pr_idx, &Stokes.pr_idx);
    adap.push_back( &PPr);
    UpdateProlongationCL<double> PLset( lset.GetMG(), lset.GetProlongation(), &lset.idx, &lset.idx);
    adap.push_back( &PLset);
    Stokes.P_ = stokessolverfactory.GetPVel();

    // For a two-level MG-solver: P2P1 -- P2P1X;
//     MakeP1P1XProlongation ( Stokes.vel_idx.NumUnknowns(), Stokes.pr_idx.NumUnknowns(),
//         Stokes.pr_idx.GetFinest().GetXidx().GetNumUnknownsStdFE(),
//         stokessolverfactory.GetPVel()->GetFinest(), stokessolverfactory.GetPPr()->GetFinest());

    // SetInitialConditions( Stokes, lset, MG, P);

    // Time discretisation + coupling
    TimeDisc2PhaseCL* timedisc= CreateTimeDisc(Stokes, lset, navstokessolver, gm, P, lsetmod);
    if (P.get<int>("Time.NumSteps") != 0){
        timedisc->SetSchurPrePtr( stokessolverfactory.GetSchurPrePtr() );
    }
    if (P.get<double>("NavStokes.Nonlinear")!=0.0 || P.get<int>("Time.NumSteps") == 0) {
        stokessolverfactory.SetMatrixA( &navstokessolver->GetAN()->GetFinest());
            //for Stokes-MGM
        stokessolverfactory.SetMatrices( navstokessolver->GetAN(), &Stokes.B.Data,
                                         &Stokes.M.Data, &Stokes.prM.Data, &Stokes.pr_idx);
    }
    else {
        stokessolverfactory.SetMatrixA( &timedisc->GetUpperLeftBlock()->GetFinest());
            //for Stokes-MGM
        stokessolverfactory.SetMatrices( timedisc->GetUpperLeftBlock(), &Stokes.B.Data,
                                         &Stokes.M.Data, &Stokes.prM.Data, &Stokes.pr_idx);
    }

    std::ofstream* infofile = 0;
    IF_MASTER {
        infofile = new std::ofstream ((P.get<std::string>("VTK.VTKName","twophasedrops")+".info").c_str());
    }
    IFInfo.Init(infofile);
    IFInfo.WriteHeader();

    if (P.get<int>("Time.NumSteps") == 0)
        SolveStatProblem( Stokes, lset, *navstokessolver);


    const double vmax = P.get<double>("Transp.MaxVelocity");
    SpaceTimeXSolutionCL sol(MG, transp_Bnd_neg, transp_Bnd_pos, P, periodic_match);
    sol.UpdateTimeSlab(oldlset,lset,vmax,NULL);
    sol.EvalTraces();

    // for serialization of geometry and numerical data
    TwoPhaseStoreCL<InstatNavierStokes2PhaseP2P1CL> ser(MG, Stokes, lset, sol,
                                                        P.get<std::string>("Restart.Outputfile"),
                                                        P.get<int>("Restart.Overwrite"),
                                                        P.get<int>("Restart.Binary"),
                                                        vel_downwind, lset_downwind);
    Stokes.v.t += GetTimeOffset();
    lset.Phi.t += GetTimeOffset();
    oldlset.Phi.t += GetTimeOffset();



    MassTranspRepairCL massrepair1( MG, sol.GetFutureTrace_Neg(), transp_Bnd_neg, lset);
    MassTranspRepairCL massrepair2( MG, sol.GetFutureTrace_Pos(), transp_Bnd_pos, lset);

    adap.push_back( &massrepair1);
    adap.push_back( &massrepair2);

    bool massTransp = P.get("Transp.DoTransp", 1);

    // for serialization of geometry and numerical data
    // if (massTransp)
    //   std::cout << "WARNING: mass transport data is not serialized, yet!" << std::endl;


    {
        VecDescCL & vpos(sol.GetFutureTrace_Pos());
        VecDescCL & vneg(sol.GetFutureTrace_Neg());

        
        if (P.get<std::string>("Restart.Inputfile") == "none"){
            for (MultiGridCL::const_TriangVertexIteratorCL sit= const_cast<const MultiGridCL&>(MG).GetTriangVertexBegin(-1), send= const_cast<const MultiGridCL&>(MG).GetTriangVertexEnd(-1);
                 sit != send; ++sit)
            {
                if (sit->Unknowns.Exist(vpos.RowIdx->GetIdx()))
                {
                    vneg.Data[sit->Unknowns(vneg.RowIdx->GetIdx())]= Initialcneg( sit->GetCoord(), 0.0);
                    vpos.Data[sit->Unknowns(vpos.RowIdx->GetIdx())]= Initialcpos( sit->GetCoord(), 0.0);
                }
            }
        }
        else
        {
            std::cout << " reading in concentration " << std::endl;
            std::string filename(P.get<std::string>("Restart.Inputfile"));
            bool binary(P.get<int>("Restart.Binary")!=0);
            ReadFEFromFile(sol.GetFutureTrace_Pos(), MG, filename + "concentration_pos", binary);
            ReadFEFromFile(sol.GetFutureTrace_Neg(), MG, filename + "concentration_neg", binary);
        }

    }

    // mass transport:
    // IdxDescCL err_idx(P0_FE);
    // VecDescCL err_vec(&err_idx);
    // err_idx.CreateNumbering( -1, MG);
    // err_vec.SetIdx(&err_idx);

    // MassTranspErrorIndicatorCL errorindicator ( MG, 
    //                                             lset,
    //                                             sol.GetFutureTrace_Neg(),
    //                                             sol.GetFutureTrace_Pos(),
    //                                             transp_Bnd_neg,
    //                                             transp_Bnd_pos,
    //                                             P.get_child("Transp"),
    //                                             err_vec);
    // if (P.get<int>("AdaptRef.AdaptiveMassTransport",0)!=0)
    // {
    //     adap.push_back( &errorindicator);
    //     // adap.SetErrorField(err_vec);
    //     // adap.SetErrorUpperThreshold(P.get<double>("AdaptRef.ErrorThresholdUpper",2e99));
    //     // adap.SetErrorLowerThreshold(P.get<double>("AdaptRef.ErrorThresholdLower",0));
    // }

    // Output-Registrations:
#ifndef _PAR
    Ensight6OutCL* ensight = NULL;
    if (P.get<int>("Ensight.EnsightOut",0)){
        // Initialize Ensight6 output
        std::string ensf( P.get<std::string>("Ensight.EnsDir") + "/" + P.get<std::string>("Ensight.EnsCase"));
        ensight = new Ensight6OutCL( P.get<std::string>("Ensight.EnsCase") + ".case",
                                     P.get<int>("Time.NumSteps")/P.get("Ensight.EnsightOut", 0)+1,
                                     P.get<int>("Ensight.MasterOut"));
        ensight->Register( make_Ensight6Geom      ( MG, MG.GetLastLevel(), P.get<std::string>("Ensight.GeomName"),
                                                    ensf + ".geo", true));
        ensight->Register( make_Ensight6Scalar    ( lset.GetSolution(),      "Levelset",      ensf + ".scl", true));
        ensight->Register( make_Ensight6Scalar    ( Stokes.GetPrSolution(),  "Pressure",      ensf + ".pr",  true));
        ensight->Register( make_Ensight6Vector    ( Stokes.GetVelSolution(), "Velocity",      ensf + ".vel", true));
        ensight->Register( make_Ensight6Scalar    ( ScalarFunAsP2EvalCL( inscamap[P.get<std::string>("NavStokes.Coeff.SurfTens.VarTensionFunc")], 0., &MG, MG.GetLastLevel()),
                                                    "Surfaceforce",  ensf + ".sf",  true));

        // if (massTransp) {
        //     ensight->Register( make_Ensight6Scalar( massTransp->GetSolution(),"Concentration", ensf + ".c",   true));
        //     ensight->Register( make_Ensight6Scalar( massTransp->GetSolution( massTransp->ct),
        //                                             "TransConc",     ensf + ".ct",  true));
        // }
        if (P.get("SurfTransp.DoTransp", 0)) {
            ensight->Register( make_Ensight6IfaceScalar( MG, surfTransp.ic,  "InterfaceSol",  ensf + ".sur", true));
        }
        if (Stokes.UsesXFEM())
            ensight->Register( make_Ensight6P1XScalar( MG, lset.Phi, Stokes.p, "XPressure",   ensf + ".pr", true));

        ensight->Write( Stokes.v.t);
    }
#endif

    // writer for vtk-format
    VTKOutCL * vtkwriter = NULL;

    // VecDescCL & vpos(sol.GetFutureTrace_Pos());
    // VecDescCL & vneg(sol.GetFutureTrace_Neg());
    P1EvalCL<double,cBndDataCL,VecDescCL > solpos(&sol.GetFutureTrace_Pos(),&transp_Bnd_pos,&MG);
    P1EvalCL<double,cBndDataCL,VecDescCL > solneg(&sol.GetFutureTrace_Neg(),&transp_Bnd_neg,&MG);

    if (P.get<int>("VTK.VTKOut",0)){
        vtkwriter = new VTKOutCL(adap.GetMG(), "DROPS data",
                                 P.get<int>("Time.NumSteps")/P.get("VTK.VTKOut", 0)+1,
                                 P.get<std::string>("VTK.VTKDir"), P.get<std::string>("VTK.VTKName"),
                                 P.get<std::string>("VTK.TimeFileName"),
                                 P.get<int>("VTK.Binary"),
                                 P.get<int>("VTK.UseOnlyP1"),
                                 false,
                                 -1,  /* <- level */
                                 P.get<int>("VTK.ReUseTimeFile"),
                                 P.get<int>("VTK.UseDeformation"));
        vtkwriter->Register( make_VTKVector( Stokes.GetVelSolution(), "velocity") );
        vtkwriter->Register( make_VTKScalar( Stokes.GetPrSolution(), "pressure") );
        if (P.get<int>("VTK.AddP1XPressure",0) && Stokes.UsesXFEM())
            vtkwriter->Register( make_VTKP1XScalar( MG, lset.Phi, Stokes.p, "xpressure"));
        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "level-set") );

        if (massTransp) {

            // vtkwriter->Register( make_VTKP1XScalar(MG, lset.Phi, sol.GetFutureTrace(), transp_Bnd_neg, transp_Bnd_pos, "conc+"));
            vtkwriter->Register( make_VTKScalar( solneg, "conc+_neg") );
            vtkwriter->Register( make_VTKScalar( solpos, "conc+_pos") );
        }

        if (P.get("SurfTransp.DoTransp", 0)) {
            vtkwriter->Register( make_VTKIfaceScalar( MG, surfTransp.ic,  "InterfaceSol"));
        }
        vtkwriter->Write(Stokes.v.t);
    }

    StrategyCombinerCL combinedmarker;

    DistMarkingStrategyCL markerlset( lset,
                                      P.get<double>("AdaptRef.Width"),
                                      P.get<double>("AdaptRef.CoarsestLevel"), 
                                      std::max(P.get<int>("AdaptRef.FinestLevel"), 
                                               P.get<int>("AdaptRef.CoarsestLevel")));

    // thresholdlist for classification of refinement lvl. in heuristically hacked
    // adaptivity ( don't develop this further, replace it!)
    double thresholdlist [10];
    for (Uint i = 0; i < 10; ++i)
    {
        std::stringstream sstr;
        sstr << "AdaptRef.Threshold";
        sstr << i;
        std::cout << sstr.str() << std::endl;
        thresholdlist[i] = P.get<double>(sstr.str(), std::numeric_limits<double>::max());
        std::cout << " thresholdlist[" << i << "] = " << thresholdlist[i] << std::endl;
    }

    std::ofstream * concout = NULL;
    ConcentrationMarkingStrategyCL markerconc (lset, solneg, solpos, 
                                               thresholdlist,
                                               P.get<int>("AdaptRef.CoarsestLevel"), 
                                               P.get<double>("AdaptRef.FinestLevel"), 
                                               P.get<int>("AdaptRef.Hacked",0) == 1, 
                                               P.get<double>("AdaptRef.HackedWidth",0.0), 
                                               concout );
    combinedmarker.push_back(markerlset);

    if (P.get<int>("AdaptRef.AddConcMarker",0)>0)
        combinedmarker.push_back(markerconc);

    adap.set_marking_strategy(&combinedmarker);
// #endif
    // if (P.get("Restart.Serialization", 0))
    //     ser.Write();

    const int nsteps = P.get<int>("Time.NumSteps");
    const double dt = P.get<double>("Time.StepSize");
    for (int step= 1; step<=nsteps; ++step)
    {
        std::cout << "============================================================ step " << step << std::endl;

        const double time_old = Stokes.v.t;
        const double time_new = Stokes.v.t + dt;

        const double told= time_old;
        const double tnew= time_new;

        IFInfo.Update( lset, Stokes.GetVelSolution());
        IFInfo.Write(time_old);

        if (P.get("SurfTransp.DoTransp", 0)) surfTransp.InitOld();
        oldlset.Phi = lset.Phi;
        oldlset.Phi.t = told;
        VelVecDescCL vold(Stokes.v);
        // vold.Data = Stokes.v.Data;
        timedisc->DoStep( P.get<int>("Coupling.Iter"));
        lset.Phi.t = tnew;

        if (massTransp){


            STVelocityContainer Flowfield (vold, Stokes.v, Stokes.GetBndData().Vel, MG);


            lset.Phi.t = tnew;
            oldlset.Phi.t = told;
            sol.UpdateTimeSlab(oldlset,lset,vmax,NULL);


            IdxDescCL & Idx(sol.GetIdx());

            TetraAccumulatorTupleCL accus ;
            MatrixCL A;
            VecDescCL b(&Idx);      
            ProgressBarTetraAccumulatorCL accup(MG,"STTranspAcc");
            accus.push_back(&accup);
            STTransportVolumeAccumulator_P1SP1TXCL accu(MG, transp_pBnd_neg, transp_pBnd_pos, &oldlset, &lset,
                                                        NULL,
                                                        RhsNeg, RhsPos, Flowfield, &A, &b, 
                                                        Idx, Idx, told, tnew, 
                                                        sol.GetFutureTrace_Neg(), 
                                                        sol.GetFutureTrace_Pos(), 
                                                        P.get_child("Transp"));
            accus.push_back(&accu);
            std::cout << " accumulate " << std::endl;
            accumulate( accus, MG, Idx.TriangLevel(), Idx.GetMatchingFunction(), Idx.GetBndInfo());

            std::cout << " solve " << std::endl;
            
            if ( P.get("Transp.DirectSolve", 0))
            {
#ifndef _PAR
#ifdef DROPS_PARDISO
                DROPS::PardisoSolverCL SolveA( A);
                SolveA.Solve(A, sol.GetSolution().Data, b.Data);
                DROPS::VectorCL r(A * sol.GetSolution().Data - b.Data);
                std::cout << "PARDISO 1 Residual:  "<< norm(r) << std::endl;
#else
    throw DROPSErrCL("PardisoSolverCL called, but MKL_HOME is not specified in CMake");
#endif
#else
                throw DROPSErrCL("no direct solver in parallel");
#endif
            }
            else
            {
                ScopeTimer scopetiming("ItSolver");


#ifndef _PAR
    typedef GSPcCL STConcPcT;
#else
    typedef JACPcCL STConcPcT;
#endif
                STConcPcT                  pc_;
                GMResSolverCL<STConcPcT> gm_( pc_, 20, P.get<int>("Transp.Iter"), P.get<double>("Transp.Tol"), 
                                              P.get<int>("Transp.RelTol",0), false, RightPreconditioning);
                gm_.Solve( A, sol.GetSolution().Data, b.Data,Idx.GetEx());
                std::cout << "res = " << gm_.GetResid() << ", iter = " << gm_.GetIter()<<"\n";
            }


            sol.EvalTraces();
        }



        if (P.get("SurfTransp.DoTransp", 0)) {
            surfTransp.DoStep( time_new);
            BndDataCL<> ifbnd( 0);
            std::cout << "surfactant on \\Gamma: " << Integral_Gamma( MG, lset.Phi, lset.GetBndData(), make_P1Eval(  MG, ifbnd, surfTransp.ic)) << '\n';
        }


        // grid modification
        const bool doGridMod= P.get<int>("AdaptRef.Freq") && step%P.get<int>("AdaptRef.Freq") == 0;
        bool gridChanged= false;
        if (doGridMod) {
            adap.UpdateTriang();
            gridChanged= adap.WasModified();
        }
        // downwind-numbering for Navier-Stokes
        const bool doNSDownwindNumbering= P.get<int>("NavStokes.Downwind.Frequency")
            && step%P.get<int>("NavStokes.Downwind.Frequency") == 0;
        if (doNSDownwindNumbering) {
            if (!gridChanged) { // We must ensure that the permutation maps the original numbering of CreateNumbering to the downwind-numbering and that the renumbering starts in a known state, i.e. the original numbering.
                vidx->DeleteNumbering( MG);
                Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx, periodic_match);
                permute_Vector( Stokes.v.Data, invert_permutation( vel_downwind), 3);
            }
            vel_downwind= Stokes.downwind_numbering( lset, navstokes_downwind);
        }
        // downwind-numbering for Levelset
        const bool doLsetDownwindNumbering= P.get<int>("Levelset.Downwind.Frequency")
            && step%P.get<int>("Levelset.Downwind.Frequency") == 0;
        if (doLsetDownwindNumbering) {
            if (!gridChanged) { // We must ensure that the permutation maps the original numbering of CreateNumbering to the downwind-numbering and that the renumbering starts in a known state, i.e. the original numbering.
                lset.DeleteNumbering( lidx);
                lset.CreateNumbering( MG.GetLastLevel(), lidx, periodic_match);
                lset.Phi.SetIdx( lidx);
                permute_Vector( lset.Phi.Data, invert_permutation( lset_downwind));
            }
            lset_downwind= lset.downwind_numbering( Stokes.GetVelSolution(), levelset_downwind);
        }
        if (gridChanged || doNSDownwindNumbering || doLsetDownwindNumbering) {
                timedisc->Update();
                // if (massTransp) massTransp->Update();
        }

#ifndef _PAR
        if (ensight && step%P.get("Ensight.EnsightOut", 0)==0)
            ensight->Write( time_new);
#endif
        if (vtkwriter && step%P.get("VTK.VTKOut", 0)==0)
            vtkwriter->Write( time_new);
        if (P.get("Restart.Serialization", 0) && step%P.get("Restart.Serialization", 0)==0)
            ser.Write();

        markerconc.ResetOutput(step);

    }

    adap.set_marking_strategy(0);
    IFInfo.Update( lset, Stokes.GetVelSolution());
    IFInfo.Write(Stokes.v.t);
    std::cout << std::endl;
    delete &lset;
    delete &oldlset;
    delete timedisc;
    delete navstokessolver;
    delete stokessolver;
    delete gm;
    // if (massTransp) delete massTransp;
    // if (transprepair) delete transprepair;
#ifndef _PAR
    if (ensight) delete ensight;
#endif
    if (vtkwriter) delete vtkwriter;
    if (infofile) delete infofile;
    if (sf) delete sf;
//     delete stokessolver1;
}


 
void  OnlyTransportStrategy( MultiGridCL& MG, LsetBndDataCL& lsetbnddata, AdapTriangCL& adap)    // do just the transport stuff
{
    rusage usage;

    InVecMap & tdvectormap = InVecMap::getInstance();
    InScaMap & scalarmap = InScaMap::getInstance();
    MatchMap & matchmap = MatchMap::getInstance();

    //instat_vector_fun_ptr Flowfield = tdvectormap[P.get<std::string>("Transp.Flow")];
    STVelocityContainer Flowfield (tdvectormap[P.get<std::string>("Transp.Flow")]);

    // instat_scalar_fun_ptr Reaction = scalarmap["ReactionFct"];
    instat_scalar_fun_ptr Initialcneg = scalarmap[P.get<std::string>("Transp.InitialConcNeg")];
    instat_scalar_fun_ptr Initialcpos = scalarmap[P.get<std::string>("Transp.InitialConcPos")];
    instat_scalar_fun_ptr RhsNeg = scalarmap[P.get<std::string>("Transp.RhsNeg")];
    instat_scalar_fun_ptr RhsPos = scalarmap[P.get<std::string>("Transp.RhsPos")];
    instat_scalar_fun_ptr SolNeg = scalarmap[P.get<std::string>("Transp.SolNeg")];
    instat_scalar_fun_ptr SolPos = scalarmap[P.get<std::string>("Transp.SolPos")];

    instat_scalar_fun_ptr distance = scalarmap[P.get<std::string>("Transp.Levelset")];

    cBndDataCL *pBnd_pos, *pBnd_neg;

    DROPS::BuildBoundaryData( &MG, pBnd_pos,  P.get<std::string>("Transp.BoundaryType"), P.get<std::string>("Transp.BoundaryFncsPos"));
    DROPS::BuildBoundaryData( &MG, pBnd_neg, P.get<std::string>("Transp.BoundaryType"), P.get<std::string>("Transp.BoundaryFncsNeg"));

    bool is_periodic = P.get<std::string>("DomainCond.PeriodicMatching", "none") != "none";
    match_fun periodic_match = is_periodic ? matchmap[P.get("DomainCond.PeriodicMatching", std::string("periodicx"))] : 0;

    cBndDataCL & Bnd_neg(*pBnd_neg);
    cBndDataCL & Bnd_pos(*pBnd_pos); 
   
    DROPS::instat_scalar_fun_ptr sigmap = 0;
    SurfaceTensionCL sf( sigmap, Bnd_neg);    
    LevelsetP2CL & lset( * LevelsetP2CL::Create( MG, lsetbnddata, sf, P.get_child("Levelset")) );
    LevelsetP2CL & oldlset( * LevelsetP2CL::Create( MG, lsetbnddata, sf, P.get_child("Levelset")) );

    //Prolongate and Restrict solution vector levelset from old mesh to new mesh after mesh adaptation:
    //always act on the same grid with possibly different interface position
    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);
    LevelsetRepairCL oldlsetrepair( oldlset);
    adap.push_back( &oldlsetrepair);

    MLIdxDescCL* lidx= &lset.idx;
    // index wrt the interface at previous time step
    MLIdxDescCL* oldlidx= &oldlset.idx; 

    lset.CreateNumbering( MG.GetLastLevel(), lidx);
    lset.Phi.SetIdx( lidx);

    oldlset.CreateNumbering( MG.GetLastLevel(), oldlidx);
    oldlset.Phi.SetIdx( oldlidx);

    SetInitialLevelsetConditions( lset, MG, P);
    SetInitialLevelsetConditions( oldlset, MG, P);

    
    const double dt = P.get<double>("Time.EndTime") / P.get<double>("Time.NumSteps");

    lset.Init( distance, dt);
    lset.Phi.t = 0;
    oldlset.Init( distance, 0);
    oldlset.Phi.t = 0;

    DisplayDetailedGeom( MG);

    //const double Vol= lset.GetVolume(); //0.5 * 0.125 * M_PI; //EllipsoidCL::GetVolume();
    std::cout << "initial volume(abs value): " << lset.GetVolume() << std::endl;
    
    //VelocityContainer vel(Stokes.v,Stokes.GetBndData().Vel,MG);
    //-> VelocityContainer vel(Flowfield);
    
    // TransportP1XCL massTransp( MG, Bnd_c, Bnd_ct, vel, lsetbnddata, lset.Phi, oldlset.Phi,P,0,Reaction,Rhs);
    // TransportXRepairCL transprepair(massTransp, MG.GetLastLevel());
    
    // index of the concentration wrt the interface at actual time step:
    // MLIdxDescCL* cidx= &massTransp.idx;
    
    // index of the concentration wrt the interface at previous time step:
    // MLIdxDescCL* cidx_old= &massTransp.oldidx; 
    
    const double vmax = P.get<double>("Transp.MaxVelocity");

    SpaceTimeXSolutionCL sol(MG, Bnd_neg, Bnd_pos, P, periodic_match);
    sol.UpdateTimeSlab(oldlset,lset,vmax,P.get<int>("Transp.Quadrature.LevelsetLinearInTime")==1?
                       NULL:distance);
    sol.EvalTraces();

    // for serialization of geometry and numerical data
    if (P.get("Transp.DoTransp", 0))
      std::cout << "WARNING: mass transport data is not serialized, yet!" << std::endl;

    // writer for vtk-format
    VTKOutCL * vtkwriter = NULL;
    if (P.get<int>("VTK.VTKOut")){
        vtkwriter = new VTKOutCL(adap.GetMG(), "DROPS data", 
                                 P.get<int>("Time.NumSteps")/P.get<int>("VTK.VTKOut")+1,
                                 P.get<std::string>("VTK.VTKDir"), P.get<std::string>("VTK.VTKName"), 
                                 P.get<std::string>("VTK.TimeFileName"),
                                 P.get<int>("VTK.Binary"), 
                                 P.get<int>("VTK.UseOnlyP1"),
                                 false,
                                 -1,  /* <- level */
                                 P.get<int>("VTK.ReUseTimeFile") );
        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "level-set+") );
        vtkwriter->Register( make_VTKScalar( oldlset.GetSolution(), "level-set-") );

        vtkwriter->Register( make_VTKP1XScalar(MG, lset.Phi, sol.GetFutureTrace(), Bnd_neg, Bnd_pos, "conc+"));
        vtkwriter->Register( make_VTKP1XScalar(MG, oldlset.Phi, sol.GetPastTrace(), Bnd_neg, Bnd_pos, "conc-"));

        vtkwriter->Write(0);
    }

    // massTransp.CheckSolution(Solutioncneg,Solutioncpos,0);
    // double cmean_old = massTransp.MeanDropConcentration();
    // double cmean_old = 0.0;

    {
        VecDescCL & vpos(sol.GetFutureTrace_Pos());
        VecDescCL & vneg(sol.GetFutureTrace_Neg());
        for (MultiGridCL::const_TriangVertexIteratorCL sit= const_cast<const MultiGridCL&>(MG).GetTriangVertexBegin(-1), send= const_cast<const MultiGridCL&>(MG).GetTriangVertexEnd(-1);
             sit != send; ++sit)
        {
            if (sit->Unknowns.Exist(vpos.RowIdx->GetIdx()))
            {
                vneg.Data[sit->Unknowns(vneg.RowIdx->GetIdx())]= Initialcneg( sit->GetCoord(), 0.0);
                vpos.Data[sit->Unknowns(vpos.RowIdx->GetIdx())]= Initialcpos( sit->GetCoord(), 0.0);
            }
        }
    }

    int step= 1;

    if (P.get<int>("Transp.CompareResults")==1)
        step = P.get<int>("Time.NumSteps");

    for (; step<=P.get<int>("Time.NumSteps"); ++step)
    {
        std::cout << "============================================================ step " << std::setw(8) << step << "  /  " << std::setw(8) << P.get<int>("Time.NumSteps") << std::endl;

        const double told= dt * (step-1);
        const double tnew= dt * step;

        lset.Init( distance, tnew);
        lset.Phi.t = tnew;
        oldlset.Init( distance, told);
        oldlset.Phi.t = told;
        sol.UpdateTimeSlab(oldlset,lset,vmax,P.get<int>("Transp.Quadrature.LevelsetLinearInTime")==1?
                           NULL:distance);

        // sol.EvalTraces();

        if (P.get<int>("Transp.CompareResults")==1)
        {
            sol.EvalTraces();
            continue;
        }



        IdxDescCL & Idx(sol.GetIdx());

        TetraAccumulatorTupleCL accus ;
        MatrixCL A;
        VecDescCL b(&Idx);

        // MatrixCL A2;
        // VecDescCL b2(&Idx);

        // MassTestAccumulator_P1SP1TXCL accu(MG, pBnd_neg, pBnd_pos, &oldlset, &lset, 
        //                                    RhsNeg, RhsPos, &A1, &b1, Idx, Idx, told, tnew);

        // SpatialLaplaceAccumulator_P1SP1TXCL accu2(MG, pBnd_neg, pBnd_pos, &oldlset, &lset, 
        //                                          RhsNeg, RhsPos, &A2, &b2, Idx, Idx, told, tnew);


        ProgressBarTetraAccumulatorCL accup(MG,"STTranspAcc");
        accus.push_back(&accup);


        MassTestAccumulator_P1SP1TXCL accumass(MG, pBnd_neg, pBnd_pos, &oldlset, &lset, 
                                               SolNeg, SolPos, &A, &b, Idx, Idx, told, tnew, P.get_child("Transp"));

        STTransportVolumeAccumulator_P1SP1TXCL accu(MG, pBnd_neg, pBnd_pos, &oldlset, &lset, 
                                                    P.get<int>("Transp.Quadrature.LevelsetLinearInTime")==1?
                                                    NULL:distance,
                                                    RhsNeg, RhsPos, Flowfield, &A, &b, 
                                                    Idx, Idx, told, tnew, 
                                                    sol.GetFutureTrace_Neg(), 
                                                    sol.GetFutureTrace_Pos(), 
                                                    P.get_child("Transp"));
        if (P.get<int>("Transp.Mass")==1)
        {
            accus.push_back(&accumass);
        }
        else
        {
            accus.push_back(&accu);
        }


        std::cout << " accumulate " << std::endl;
        accumulate( accus, MG, Idx.TriangLevel(), Idx.GetMatchingFunction(), Idx.GetBndInfo());

        std::cout << " solve " << std::endl;
        // former hacked stuff: if ( P.get("Transp.SpecialSolve", 0))
        if (false)
        {
        }
        else
        {
            if ( P.get("Transp.DirectSolve", 0))
            {
#ifndef _PAR
#ifdef DROPS_PARDISO
                DROPS::PardisoSolverCL SolveA( A);
                SolveA.Solve(A, sol.GetSolution().Data, b.Data);
                DROPS::VectorCL r(A * sol.GetSolution().Data - b.Data);
                std::cout << "PARDISO 1 Residual:  "<< norm(r) << std::endl;
#else
    throw DROPSErrCL("PardisoSolverCL called, but MKL_HOME is not specified in CMake");
#endif
#else
                throw DROPSErrCL("no direct solver in parallel");
#endif
            }
            else
            {
                ScopeTimer scopetiming("ItSolver");

#ifndef _PAR
                typedef GSPcCL STConcPcT;
#else
                typedef JACPcCL STConcPcT;
#endif
                STConcPcT                  pc_;
                GMResSolverCL<STConcPcT> gm_( pc_, 20, P.get<int>("Transp.Iter"), P.get<double>("Transp.Tol"), 
                                              P.get<int>("Transp.RelTol",0), false, RightPreconditioning);
                gm_.Solve( A, sol.GetSolution().Data, b.Data,Idx.GetEx());
                std::cout << "res = " << gm_.GetResid() << ", iter = " << gm_.GetIter()<<"\n";
            }
        }

        std::cout << " visualize " << std::endl;
        sol.EvalTraces();
        sol.CheckFutureSolution(SolNeg,SolPos,tnew);
        {

            ProgressBarTetraAccumulatorCL accup(MG,"Error calculation");

            TetraAccumulatorTupleCL accus ;
            accus.push_back(&accup);

            STGeomApproxTestAccumulatorCL accu3(MG, &oldlset, &lset,  
                                                P.get<int>("Transp.Quadrature.LevelsetLinearInTime")==1?
                                                NULL:distance,
                                                told, tnew, P.get_child("Transp"));
            if (P.get<int>("Transp.Debug")==1)
            {
                accus.push_back(&accu3);
            }



            InterfaceJumpAccumulatorCL accu(MG, &oldlset, &lset, 
                                            P.get<int>("Transp.Quadrature.LevelsetLinearInTime")==1?
                                            NULL:distance,
                                            told, tnew, 
                                            sol.GetPastTrace_Neg(), 
                                            sol.GetPastTrace_Pos(), 
                                            sol.GetFutureTrace_Neg(), 
                                            sol.GetFutureTrace_Pos(), 
                                            *pBnd_neg, *pBnd_pos,
                                            P.get_child("Transp"));

            accus.push_back(&accu);
            EnergyNormErrorAccumulatorCL accuen(MG, &oldlset, &lset, 
                                            P.get<int>("Transp.Quadrature.LevelsetLinearInTime")==1?
                                            NULL:distance,
                                            told, tnew, 
                                            sol.GetPastTrace_Neg(), 
                                            sol.GetPastTrace_Pos(), 
                                            sol.GetFutureTrace_Neg(), 
                                            sol.GetFutureTrace_Pos(), 
                                            *pBnd_neg, *pBnd_pos,
                                            P.get_child("Transp"));
            
            accus.push_back(&accuen);
            std::cout << " accumulate " << std::endl;
            accumulate( accus, MG, Idx.TriangLevel(), Idx.GetMatchingFunction(), Idx.GetBndInfo());
        }


        getrusage( RUSAGE_SELF, &usage);
        printf( "Memory usage | ru_maxrss: %li kB.\n", usage.ru_maxrss);

        bool vtkoutnow = vtkwriter && (step%P.get<int>("VTK.VTKOut")==0 || step < 20);
        if (vtkoutnow)
            vtkwriter->Write(tnew);
    }
    
    // std::cout << " before destructors ... "<< std::endl;
    std::cout << std::endl;


    if ( P.get("Transp.SaveResult", 0))
    {
        std::cout << " outputting " << std::endl;

        WriteFEToFile( sol.GetFutureTrace_Neg(), MG, "result_src_neg_out", true);
        WriteFEToFile( sol.GetFutureTrace_Pos(), MG, "result_src_pos_out", true);

        std::cout << " outputting end... " << std::endl;
    }
    // VecDescCL 


    if ( P.get("Transp.CompareResult", 0))
    {
        std::cout << " inputting " << std::endl;

        VectorCL tmp_neg = sol.GetFutureTrace_Neg().Data;
        VectorCL tmp_pos = sol.GetFutureTrace_Pos().Data;

        ReadFEFromFile( sol.GetFutureTrace_Neg(), MG, "result_src_neg_in", true);
        ReadFEFromFile( sol.GetFutureTrace_Pos(), MG, "result_src_pos_in", true);

        sol.GetFutureTrace_Neg().Data -= tmp_neg;
        sol.GetFutureTrace_Pos().Data -= tmp_pos;

        std::cout << " inputting end... " << std::endl;

        // std::cout << " sol.GetFutureTrace_Neg().Data = " << sol.GetFutureTrace_Neg().Data << std::endl;
        // std::cout << " sol.GetFutureTrace_Pos().Data = " << sol.GetFutureTrace_Pos().Data << std::endl;

        sol.CalcL2NormOfFuture();

        sol.GetFutureTrace_Neg().Data = tmp_neg;
        sol.GetFutureTrace_Pos().Data = tmp_pos;

        std::cout << " comparing end... " << std::endl;
    }

    if ( P.get("Transp.CompareResults", 0))
    {
        std::cout << " inputting " << std::endl;

        ReadFEFromFile( sol.GetFutureTrace_Neg(), MG, P.get<std::string>("Transp.Results.Neg1","result_neg1"), true);
        ReadFEFromFile( sol.GetFutureTrace_Pos(), MG, P.get<std::string>("Transp.Results.Pos1","result_pos1"), true);

        VectorCL tmp_neg = sol.GetFutureTrace_Neg().Data;
        VectorCL tmp_pos = sol.GetFutureTrace_Pos().Data;

        ReadFEFromFile( sol.GetFutureTrace_Neg(), MG, P.get<std::string>("Transp.Results.Neg2","result_neg2"), true);
        ReadFEFromFile( sol.GetFutureTrace_Pos(), MG, P.get<std::string>("Transp.Results.Pos2","result_pos2"), true);

        sol.GetFutureTrace_Neg().Data -= tmp_neg;
        sol.GetFutureTrace_Pos().Data -= tmp_pos;

        std::cout << " inputting end... " << std::endl;

        // std::cout << " sol.GetFutureTrace_Neg().Data = " << sol.GetFutureTrace_Neg().Data << std::endl;
        // std::cout << " sol.GetFutureTrace_Pos().Data = " << sol.GetFutureTrace_Pos().Data << std::endl;

        sol.CalcL2NormOfFuture();

        std::cout << " comparing end... " << std::endl;
    }


    delete &lset;
    delete &oldlset;
    delete pBnd_neg;
    delete pBnd_pos;
    delete vtkwriter;
}


} // end of namespace DROPS

/// \brief Set Default parameters here s.t. they are initialized.
/// The result can be checked when Param-list is written to the output.
void SetMissingParameters(DROPS::ParamCL& P){

    P.put_if_unset<int>("Transp.DoTransp",0);
    P.put_if_unset<std::string>("Transp.Levelset","Ellipsoid");
    P.put_if_unset<std::string>("Transp.Flow","ZeroVel");
    P.put_if_unset<std::string>("Transp.RhsNeg","Zero");
    P.put_if_unset<std::string>("Transp.RhsPos","Zero");
    P.put_if_unset<std::string>("Transp.InitialConcNeg","IniCnegFct");
    P.put_if_unset<std::string>("Transp.InitialConcPos","IniCposFct");
    P.put_if_unset<std::string>("Transp.BoundaryType","21!21!21!21!21!21");
    P.put_if_unset<std::string>("Transp.BoundaryFncs","Zero!Zero!Zero!Zero!Zero!Zero");
    P.put_if_unset<std::string>("Transp.BoundaryFncst","Zero!Zero!Zero!Zero!Zero!Zero");
    P.put_if_unset<int>("Transp.KappaRule",0);
    P.put_if_unset<int>("Transp.Debug",0);
    P.put_if_unset<int>("Transp.Mass",0);
    P.put_if_unset<int>("Transp.NoConvection",0);
    P.put_if_unset<int>("Transp.MaxVelocity",1.0);
    P.put_if_unset<int>("Transp.AntiSymmetricConvection",0);
    P.put_if_unset<int>("Transp.RhsApproximationOrderInTime",2);
    P.put_if_unset<int>("Transp.ConvectionApproximationOrderInTime",1);
    P.put_if_unset<int>("TestCase5.DeformationCase",0);

    P.put_if_unset<int>("SurfTransp.DoTransp",0);

    P.put_if_unset<std::string>("VTK.TimeFileName",P.get<std::string>("VTK.VTKName"));
    P.put_if_unset<int>("VTK.ReUseTimeFile",0);
    P.put_if_unset<int>("VTK.UseOnlyP1",0);
    P.put_if_unset<int>("VTK.VTKOut",0);
    P.put_if_unset<std::string>("VTK.VTKName","ns_transp");

    P.put_if_unset<int>("Ensight.EnsightOut",0);
    P.put_if_unset<int>("Levelset.Discontinuous", 0);
    P.put_if_unset<int>("Levelset.SD", 1);
    P.put_if_unset<int>("Levelset.CurvDiff", 0);

    P.put_if_unset<double>("NavStokes.Nonlinear", 0.0);

    P.put_if_unset<std::string>("Restart.Inputfile","none");
    P.put_if_unset<int>("NavStokes.Downwind.Frequency", 0);
    P.put_if_unset<double>("NavStokes.Downwind.MaxRelComponentSize", 0.05);
    P.put_if_unset<double>("NavStokes.Downwind.WeakEdgeRatio", 0.2);
    P.put_if_unset<double>("NavStokes.Downwind.CrosswindLimit", std::cos( M_PI/6.));
    P.put_if_unset<int>("Levelset.Downwind.Frequency", 0);
    P.put_if_unset<double>("Levelset.Downwind.MaxRelComponentSize", 0.05);
    P.put_if_unset<double>("Levelset.Downwind.WeakEdgeRatio", 0.2);
    P.put_if_unset<double>("Levelset.Downwind.CrosswindLimit", std::cos( M_PI/6.));

    P.put_if_unset<std::string>("Exp.VolForce", "ZeroVel");
    P.put_if_unset<double>("Mat.DensDrop", 0.0);
    P.put_if_unset<double>("Mat.ShearVisco", 0.0);
    P.put_if_unset<double>("Mat.DilatationalVisco", 0.0);
    P.put_if_unset<double>("SurfTens.ShearVisco", 0.0);
    P.put_if_unset<double>("SurfTens.DilatationalVisco", 0.0);

    P.put_if_unset<int>("General.ProgressBar", 0);
    P.put_if_unset<std::string>("General.DynamicLibsPrefix", "../");
}

int main (int argc, char** argv)
{
#ifdef _PAR
    DROPS::ProcCL::Instance(&argc, &argv);
#endif
  try
  {
    std::cout << "Boost version: " << BOOST_LIB_VERSION << std::endl;

    DROPS::read_parameter_file_from_cmdline( P, argc, argv, "../../param/spacetimetransp/st_transp/risingdroplet.json");
    SetMissingParameters(P);
    std::cout << P << std::endl;

    DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"), P.get<std::vector<std::string> >("General.DynamicLibs") );

    if (P.get<int>("General.ProgressBar"))
        DROPS::ProgressBarTetraAccumulatorCL::Activate();

    // check parameter file
    if (P.get<double>("SurfTens.DilatationalVisco")< P.get<double>("SurfTens.ShearVisco"))
    {
        throw DROPS::DROPSErrCL("Parameter error : Dilatational viscosity must be larger than surface shear viscosity");
    }

    DROPS::MatchMap & matchmap = DROPS::MatchMap::getInstance();
    bool is_periodic = P.get<std::string>("DomainCond.PeriodicMatching", "none") != "none";
    DROPS::match_fun periodic_match = is_periodic ? matchmap[P.get<std::string>("DomainCond.PeriodicMatching", "periodicx")] : 0;

    DROPS::MultiGridCL* mg= 0;
    typedef DROPS::BndDataCL<DROPS::Point3DCL> VelBndDataCL;
    typedef DROPS::BndDataCL<double>    PrBndDataCL;
    DROPS::LsetBndDataCL* lsetbnddata= 0;
    VelBndDataCL *velbnddata = 0;
    PrBndDataCL *prbnddata = 0;

    //you cannot pass a double& per P.get, so you need to use this indirect way
    double ExpRadInlet = P.get<double>("Exp.RadInlet");

    try
    {
        std::auto_ptr<DROPS::MGBuilderCL> builder( DROPS::make_MGBuilder( P));
        mg = new DROPS::MultiGridCL( *builder);
    }
    catch (DROPS::DROPSParamErrCL& e)
    {
        std::cout << "\n"
                  << "  /----------------------------------------------------------------\\ \n"
                  << "  | WARNING: It seems you are using the old domain descriptions    | \n"
                  << "  |          or your \"Domain\" section is not correct.              | \n"
                  << "  |          Please adapt your json-file to the new description.   | \n"
                  <<"  \\----------------------------------------------------------------/ \n"
                  << std::endl;
        DROPS::BuildDomain( mg, P.get<std::string>("DomainCond.MeshFile"), P.get<int>("DomainCond.GeomType"), P.get<std::string>("Restart.Inputfile"), ExpRadInlet);
    }




    P.put("Exp.RadInlet", ExpRadInlet);

    std::cout << "Generated MG of " << mg->GetLastLevel() << " levels." << std::endl;

    std::string perbndtypestr;
    std::string zerobndfun;
    for( size_t i= 1; i<=mg->GetBnd().GetNumBndSeg(); ++i) {
        zerobndfun += "Zero";
        if (i!=mg->GetBnd().GetNumBndSeg())
          zerobndfun += "!";
    }
    DROPS::BuildBoundaryData( mg, velbnddata, P.get<std::string>("DomainCond.BoundaryType"), P.get<std::string>("DomainCond.BoundaryFncs"), periodic_match, &perbndtypestr);
    std::cout << "Generated boundary conditions for velocity, ";
    DROPS::BuildBoundaryData( mg, prbnddata, perbndtypestr, zerobndfun, periodic_match);
    std::cout << "pressure, ";
    DROPS::BuildBoundaryData( mg, lsetbnddata, perbndtypestr, zerobndfun, periodic_match);
    std::cout << "and levelset." << std::endl;
    DROPS::StokesBndDataCL bnddata(*velbnddata,*prbnddata);

    std::cout << "Generated MG of " << mg->GetLastLevel() << " levels." << std::endl;

    std::string InitialLSet= P.get("Exp.InitialLSet", std::string("Ellipsoid"));
    if (InitialLSet == "Ellipsoid")
        DROPS::EllipsoidCL::Init( P.get<DROPS::Point3DCL>("Exp.PosDrop"), P.get<DROPS::Point3DCL>("Exp.RadDrop"));
    if  (InitialLSet == "TwoEllipsoid")
        DROPS::TwoEllipsoidCL::Init( P.get<DROPS::Point3DCL>("Exp.PosDrop"), P.get<DROPS::Point3DCL>("Exp.RadDrop"), P.get<DROPS::Point3DCL>("Exp.PosDrop2"), P.get<DROPS::Point3DCL>("Exp.RadDrop2"));
    if (InitialLSet.find("Cylinder")==0) {
        DROPS::CylinderCL::Init( P.get<DROPS::Point3DCL>("Exp.PosDrop"), P.get<DROPS::Point3DCL>("Exp.RadDrop"), InitialLSet[8]-'X');
        P.put("Exp.InitialLSet", InitialLSet= "Cylinder");
    }

    using DROPS::EllipsoidCL;
    using DROPS::AdapTriangCL;
    using DROPS::DistMarkingStrategyCL;

    DROPS::instat_scalar_fun_ptr distance = NULL;
    if (!P.get<int>("Transp.UseNSSol",1)){
        DROPS::InScaMap & scalarmap = DROPS::InScaMap::getInstance();
        distance = scalarmap[P.get("Transp.Levelset", std::string("Ellipsoid"))];
    }
    else if (P.get<std::string>("Restart.Inputfile") == "none"){
        DROPS::InScaMap & scalarmap = DROPS::InScaMap::getInstance();
        distance = scalarmap[P.get("Exp.InitialLSet", std::string("Ellipsoid"))];
    }

    DROPS::DistMarkingStrategyCL InitialMarker( distance,
                           P.get<double>("AdaptRef.Width"),
                           P.get<int>("AdaptRef.CoarsestLevel"), P.get<int>("AdaptRef.FinestLevel"));

    DROPS::AdapTriangCL adap( *mg, &InitialMarker,
                              ((P.get<std::string>("Restart.Inputfile") == "none") ? P.get<int>("AdaptRef.LoadBalStrategy") : -P.get<int>("AdaptRef.LoadBalStrategy")));
    // If we read the Multigrid, it shouldn't be modified;
    // otherwise the pde-solutions from the ensight files might not fit.
    if (P.get("Restart.Inputfile", std::string("none")) == "none")
    {
        adap.MakeInitialTriang();
    }

    std::cout << DROPS::SanityMGOutCL(*mg) << std::endl;
#ifdef _PAR
    if ( DROPS::ProcCL::Check( DROPS::DiST::InfoCL::Instance().IsSane( std::cerr)))
        std::cout << " DiST-module seems to be alright!" << std::endl;
    else
        std::cout << " DiST-module seems to be broken!" << std::endl;
    if ( DROPS::CheckParMultiGrid())
        std::cout << "As far as I can tell the ParMultigridCL is sane\n";
#endif

    
    if (P.get<int>("Transp.UseNSSol") == 1) 
    {
        DROPS::InstatNavierStokes2PhaseP2P1CL prob( *mg, DROPS::TwoPhaseFlowCoeffCL(P), bnddata, P.get<double>("Stokes.XFEMStab")<0 ? DROPS::P1_FE : DROPS::P1X_FE, P.get<double>("Stokes.XFEMStab"));

        Strategy( prob, *lsetbnddata, adap);    // do all the stuff
    }
    else
        OnlyTransportStrategy( *mg, *lsetbnddata, adap);    // do just the transport stuff

    delete mg;
    delete velbnddata;
    delete prbnddata;
    delete lsetbnddata;

    rusage usage;
    getrusage( RUSAGE_SELF, &usage);

#ifdef _PAR
    printf( "[%i]: ru_maxrss: %li kB.\n", DROPS::ProcCL::MyRank(), usage.ru_maxrss);
#else
    printf( "ru_maxrss: %li kB.\n", usage.ru_maxrss);
#endif
    std::cout << " twophasedrops finished regularly" << std::endl;
    return 0;
  }
  catch (DROPS::DROPSErrCL& err) { err.handle(); }
}

