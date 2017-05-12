/// \file twophasedrops.cpp
/// \brief flow in measurement cell or brick
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Christoph Lehrenfeld; SC RWTH Aachen: Oliver Fortmeier

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
#ifndef _PAR
#include "out/ensightOut.h"
#endif
#include "out/vtkOut.h"
//levelset
#include "levelset/coupling.h"
#include "levelset/marking_strategy.h"
#include "levelset/adaptriang.h"
#include "levelset/mzelle_hdr.h"
#include "levelset/twophaseutils.h"
//surfactants
#include "surfactant/ifacetransp.h"
//function map
#include "misc/funcmap.h"
//solver factory for stokes
#include "num/stokessolverfactory.h"
#include "num/oseensolver.h"
#include "num/prolongation.h"
#ifdef _PAR
#include "parallel/loadbal.h"
#include "parallel/parmultigrid.h"
#endif
//general: streams
#include <fstream>
#include <sstream>

#ifndef _PAR
#include "num/stokespardiso.h"
#endif
#include "misc/progressaccu.h"
#include "misc/dynamicload.h"

#include <sys/resource.h>

DROPS::ParamCL P;

// rho*du/dt - mu*laplace u + Dp = f + rho*g - okn
//                        -div u = 0
//                             u = u0, t=t0

namespace DROPS // for Strategy
{

double GetTimeOffset(){
    double timeoffset = 0.0;
    const std::string restartfilename = P.get<std::string>("Restart.InputData");
    if( ReadInitialConditionFromFile(P) )
    {
        const std::string timefilename = restartfilename + "time";
        std::ifstream f_(timefilename.c_str());
        f_ >> timeoffset;
        std::cout << "used time offset file is " << timefilename << std::endl;
        std::cout << "time offset is " << timeoffset << std::endl;
    }
    return timeoffset;
}

void Strategy( InstatNavierStokes2PhaseP2P1CL& Stokes, LsetBndDataCL& lsetbnddata, AdapTriangCL& adap, const bool is_periodic, const std::string& perMatchName)
// flow control
{
    DROPS::InScaMap & inscamap = DROPS::InScaMap::getInstance();
    //DROPS::ScaMap & scamap = DROPS::ScaMap::getInstance();
    DROPS::InVecMap & invecmap = DROPS::InVecMap::getInstance();
    //DROPS::MatchMap & matchmap = DROPS::MatchMap::getInstance();

    MultiGridCL& MG= Stokes.GetMG();

    // initialization of surface tension
    // choose a proper model for surface tension coefficient, see levelset/surfacetension.h
    instat_scalar_fun_ptr sigmap = inscamap[P.get<std::string>("NavStokes.Coeff.SurfTens.VarTensionFunc", "ConstTau")];
    SurfaceTensionCL * sf;
    sf = new SurfaceTensionCL( sigmap);
    sf->SetInputMethod( Sigma_X);

    // Creates new Levelset-Object, has to be cleaned manually
    LevelsetP2CL & lset( * LevelsetP2CL::Create( MG, lsetbnddata, *sf, P.get_child("Levelset")) );

    //required to simulate flows with moving contact line
    instat_scalar_fun_ptr Young_angle = inscamap[P.get<std::string>("NavStokes.BoundaryData.SlipBnd.ContactAngleFunc")];
    instat_vector_fun_ptr bnd_outnormal = invecmap[P.get<std::string>("NavStokes.BoundaryData.SlipBnd.BndOuterNormal")];
    Stokes.SetYoungAngle(Young_angle);
    Stokes.SetBndOutNormal(bnd_outnormal);
    Stokes.SetSurfTension(sf);
    
    if (is_periodic) /// \todo Periodic directions (used for reparam) should be set based on Mesh.PeriodicBnd. Export to function!
    {
        int n = 0;
        if (perMatchName == "periodicx" || perMatchName == "periodicy" || perMatchName == "periodicz")
            n = 1;
        if (perMatchName == "periodicxy" || perMatchName == "periodicxz" || perMatchName == "periodicyz")
            n = 2;
        LevelsetP2CL::perDirSetT pdir(n);
        if (perMatchName == "periodicx") pdir[0] = P.get<Point3DCL>("Mesh.E1");
        if (perMatchName == "periodicy") pdir[0] = P.get<Point3DCL>("Mesh.E2");
        if (perMatchName == "periodicz") pdir[0] = P.get<Point3DCL>("Mesh.E3");
        if (perMatchName == "periodicxy") {pdir[0] = P.get<Point3DCL>("Mesh.E1"); pdir[1] = P.get<Point3DCL>("Mesh.E2");}
        if (perMatchName == "periodicxz") {pdir[0] = P.get<Point3DCL>("Mesh.E1"); pdir[1] = P.get<Point3DCL>("Mesh.E3");}
        if (perMatchName == "periodicyz") {pdir[0] = P.get<Point3DCL>("Mesh.E2"); pdir[1] = P.get<Point3DCL>("Mesh.E3");}
        if (perMatchName != "periodicx" && perMatchName != "periodicy" && perMatchName != "periodicz" &&
          perMatchName != "periodicxy" && perMatchName != "periodicxz" && perMatchName != "periodicyz"){
            std::cout << "WARNING: could not set periodic directions! Reparametrization can not work correctly now!" << std::endl;
            std::cout << "Press any key to continue" << std::endl; getchar();
        }
        lset.SetPeriodicDirections(&pdir);
    }

    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);
    VelocityRepairCL velrepair( Stokes);
    adap.push_back( &velrepair);
    PressureRepairCL prrepair( Stokes, lset);
    adap.push_back( &prrepair);

    MLIdxDescCL* lidx= &lset.idx;
    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;

    ParamCL PSolver( P.get_child("CouplingSolver.NavStokesSolver.OseenSolver") );
    if ( StokesSolverFactoryHelperCL().VelMGUsed(PSolver)){
        Stokes.SetNumVelLvl ( Stokes.GetMG().GetNumLevel());
        lset.SetNumLvl(Stokes.GetMG().GetNumLevel());
    }
    if ( StokesSolverFactoryHelperCL().PrMGUsed(PSolver)){
        Stokes.SetNumPrLvl  ( Stokes.GetMG().GetNumLevel());
        lset.SetNumLvl(Stokes.GetMG().GetNumLevel());
    }
    lset.CreateNumbering( MG.GetLastLevel());

    if (lset.IsDiscontinuous())
    {
        LevelsetP2DiscontCL& lsetD (dynamic_cast<LevelsetP2DiscontCL&>(lset));
        MLIdxDescCL* lidxc = lsetD.idxC;
        lsetD.CreateNumbering( MG.GetLastLevel(), lidxc);
        lsetD.PhiContinuous.SetIdx( lidxc);
    }

    PermutationT lset_downwind;

    ///\todo Is this really necessary? Already done ~20 lines above!
    if ( StokesSolverFactoryHelperCL().VelMGUsed(PSolver))
        Stokes.SetNumVelLvl ( Stokes.GetMG().GetNumLevel());
    if ( StokesSolverFactoryHelperCL().PrMGUsed(PSolver))
        Stokes.SetNumPrLvl  ( Stokes.GetMG().GetNumLevel());

    SetInitialLevelsetConditions( lset, MG, P);

    double Vol = 0;
    std::string InitialLSet= P.get("Levelset.InitialValue", std::string("Ellipsoid"));
    if ( (InitialLSet == "Ellipsoid"     || InitialLSet == "Cylinder" || InitialLSet == "ContactDroplet"
        || InitialLSet == "HalfEllipsoid" || InitialLSet == "TaylorFlowDistance"))
    {  
        if (P.get<double>("Levelset.InitialVolume",-1.0) > 0 )
            Vol = P.get<double>("Levelset.InitialVolume");
        if (InitialLSet == "Ellipsoid")
            Vol = EllipsoidCL::GetVolume();
        if (InitialLSet == "HalfEllipsoid")
            Vol = HalfEllipsoidCL::GetVolume();
        if (InitialLSet == "ContactDroplet")
            Vol = ContactDropletCL::GetVolume();
        if (InitialLSet.find("Cylinder")==0)
            Vol = CylinderCL::GetVolume();
        std::cout << "initial rel. volume: " << lset.GetVolume()/Vol << std::endl;
        lset.InitVolume( Vol);
        lset.AdjustVolume();
        std::cout << "initial lset volume adjustment:\n";
        lset.GetVolumeAdjuster()->DebugOutput( std::cout);
    } else {
        Vol = lset.GetVolume();
        lset.InitVolume( Vol);
    }

    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, &lset);
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
    if (P.get<int>("NavStokes.ShiftFrame") == 1)
        Stokes.InitVel( &Stokes.v, InVecMap::getInstance().find("InflowShiftFrame")->second);  // shifted zero velocity initial condition
    else
        Stokes.InitVel( &Stokes.v, ZeroVel);

    IteratedDownwindCL navstokes_downwind( P.get_child( "CouplingSolver.NavStokesSolver.Downwind"));
    if (P.get<int>( "CouplingSolver.NavStokesSolver.Downwind.Frequency") > 0) {
        if (StokesSolverFactoryHelperCL().VelMGUsed( PSolver))
            throw DROPSErrCL( "Strategy: Multigrid-solver and downwind-numbering cannot be used together. Sorry.\n");
        vel_downwind= Stokes.downwind_numbering( lset, navstokes_downwind);
    }
    IteratedDownwindCL levelset_downwind( P.get_child( "CouplingSolver.LevelsetSolver.Downwind"));
    if (P.get<int>( "CouplingSolver.LevelsetSolver.Downwind.Frequency") > 0)
        lset_downwind= lset.downwind_numbering( Stokes.GetVelSolution(), levelset_downwind);

    DisplayDetailedGeom( MG);
    DisplayUnks(Stokes, lset, MG);

    const int nsteps = P.get<int>("Time.NumSteps");
    const double tEnd = P.get<double>("Time.FinalTime");
    const double dt = tEnd / nsteps;

    TransportP1CL * massTransp = nullptr;
    TransportRepairCL *  transprepair = nullptr;

    if( P.get<bool>("Transp.Enable") )
    {
        // CL: the following could be moved outside of strategy to some function like
        //" InitializeMassTransport(P,MG,Stokes,lset,adap, TransportP1CL * & massTransp,TransportRepairCL * & transprepair)"
        static DROPS::BndCondT c_bc[6]= {
            DROPS::OutflowBC, DROPS::OutflowBC, DROPS::OutflowBC,
            DROPS::OutflowBC, DROPS::OutflowBC, DROPS::OutflowBC
        };
        static DROPS::BndDataCL<>::bnd_val_fun c_bfun[6]= {0, 0, 0, 0, 0, 0};
        static DROPS::BndDataCL<> Bnd_c( 6, c_bc, c_bfun);
        double D[2] = {P.get<double>("Transp.DiffPos"), P.get<double>("Transp.DiffNeg")};

        massTransp = new TransportP1CL( MG, Bnd_c, Stokes.GetBndData().Vel, P.get<double>("Time.Theta"),
                                  D, P.get<double>("Transp.HenryNeg")/P.get<double>("Transp.HenryPos"), &Stokes.v, lset,
                                  dt, P.get<int>("Transp.Solver.Iter"), P.get<double>("Transp.Solver.Tol"));

        transprepair = new TransportRepairCL(*massTransp, MG);
        adap.push_back(transprepair);

        MLIdxDescCL* cidx= &massTransp->idx;
        massTransp->CreateNumbering( MG.GetLastLevel(), cidx);
        massTransp->ct.SetIdx( cidx);
        //if (P.get<int>("DomainCond.InitialCond") != -1)
        if( !ReadInitialConditionFromFile(P) )
            massTransp->Init( inscamap["Initialcneg"], inscamap["Initialcpos"]);
        else
            ReadFEFromFile( massTransp->ct, MG, P.get<std::string>("Restart.InputData")+"concentrationTransf");

        massTransp->Update();
        std::cout << massTransp->c.Data.size() << " concentration unknowns,\n";
    }

    // TL: can we make a pointer out of this? like massTransp
    /// \todo rhs beruecksichtigen
    SurfactantcGP1CL *surfTransp = nullptr;
    InterfaceP1RepairCL *surf_repair = nullptr;

    if( P.get<bool>("SurfTransp.Enable") )
    {
        surfTransp = new SurfactantcGP1CL( MG, Stokes.GetBndData().Vel, P.get<double>("Time.Theta"), P.get<double>("SurfTransp.Visc"), &Stokes.v, *lset.PhiC, lset.GetBndData(),
                                     dt, P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"), P.get<double>("SurfTransp.XFEMReduced"));
        surf_repair = new InterfaceP1RepairCL ( MG, *lset.PhiC, lset.GetBndData(), surfTransp->ic);
        adap.push_back( surf_repair);
        surfTransp->idx.CreateNumbering( MG.GetLastLevel(), MG, lset.PhiC, &lset.GetBndData());
        std::cout << "Surfactant transport: NumUnknowns: " << surfTransp->idx.NumUnknowns() << std::endl;
        surfTransp->ic.SetIdx( &surfTransp->idx);
        surfTransp->Init( inscamap["surf_sol"]);
    }

    // Stokes-Solver    
    ParamCL PTime( P.get_child("Time") );
    StokesSolverFactoryCL<InstatNavierStokes2PhaseP2P1CL> stokessolverfactory(Stokes, PSolver, PTime );
    StokesSolverBaseCL* stokessolver;

    if (! P.get<int>("CouplingSolver.NavStokesSolver.OseenSolver.DirectSolve"))
        stokessolver = stokessolverfactory.CreateStokesSolver();
#ifndef _PAR
    else
        stokessolver = new StokesPardisoSolverCL();
#else
    else
        throw DROPSErrCL("no direct solver in parallel");
#endif

//  comment: construction of a oseen solver, preconditioned by another oseen solver,
//           e.g. GCR preconditioned by Vanka-MG, do not forget to delete stokessolver1 at the end of strategy
//
//    StokesSolverBaseCL* stokessolver1 = stokessolverfactory.CreateStokesSolver();
//    StokesSolverAsPreCL pc (*stokessolver1, 1);
//    GCRSolverCL<StokesSolverAsPreCL> gcr(pc, C.stk_OuterIter, C.stk_OuterIter, C.stk_OuterTol, /*rel*/ false);
//    BlockMatrixSolverCL<GCRSolverCL<StokesSolverAsPreCL> >* stokessolver =
//            new BlockMatrixSolverCL<GCRSolverCL<StokesSolverAsPreCL> > (gcr);

    // Navier-Stokes-Solver
    NSSolverBaseCL<InstatNavierStokes2PhaseP2P1CL>* navstokessolver = 0;
    if (P.get<double>("CouplingSolver.NavStokesSolver.Nonlinear")==0.0)
        navstokessolver = new NSSolverBaseCL<InstatNavierStokes2PhaseP2P1CL>(Stokes, *stokessolver);
    else
        navstokessolver = new AdaptFixedPtDefectCorrCL<InstatNavierStokes2PhaseP2P1CL>(Stokes, *stokessolver,
                                                                                       P.get<int>("CouplingSolver.NavStokesSolver.Iter"),
                                                                                       P.get<double>("CouplingSolver.NavStokesSolver.Tol"),
                                                                                       P.get<double>("CouplingSolver.NavStokesSolver.Reduction"));

    // Level-Set-Solver
#ifndef _PAR
    typedef GSPcCL  LsetPcT;
#else
    typedef JACPcCL LsetPcT;
#endif
    LsetPcT lset_pc;
    GMResSolverCL<LsetPcT>* gm = new GMResSolverCL<LsetPcT>( lset_pc, 200, P.get<int>("CouplingSolver.LevelsetSolver.Iter"), P.get<double>("CouplingSolver.LevelsetSolver.Tol"));

    LevelsetModifyCL lsetmod( P.get<int>("Levelset.Reparam.Freq"), P.get<int>("Levelset.Reparam.Method"), P.get<double>("Levelset.Reparam.MaxGrad"), P.get<double>("Levelset.Reparam.MinGrad"), is_periodic);

    UpdateProlongationCL<Point3DCL> PVel( Stokes.GetMG(), stokessolverfactory.GetPVel(), &Stokes.vel_idx, &Stokes.vel_idx);
    adap.push_back( &PVel);
    UpdateProlongationCL<double> PPr ( Stokes.GetMG(), stokessolverfactory.GetPPr(), &Stokes.pr_idx, &Stokes.pr_idx);
    adap.push_back( &PPr);
    UpdateProlongationCL<double> PLset( lset.GetMG(), lset.GetProlongation(), lset.idxC, lset.idxC);
    adap.push_back( &PLset);
    Stokes.P_ = stokessolverfactory.GetPVel();

    // For a two-level MG-solver: P2P1 -- P2P1X;
//     MakeP1P1XProlongation ( Stokes.vel_idx.NumUnknowns(), Stokes.pr_idx.NumUnknowns(),
//         Stokes.pr_idx.GetFinest().GetXidx().GetNumUnknownsStdFE(),
//         stokessolverfactory.GetPVel()->GetFinest(), stokessolverfactory.GetPPr()->GetFinest());

    SetInitialConditions( Stokes, lset, MG, P);

    // Time discretisation + coupling
    TimeDisc2PhaseCL* timedisc= CreateTimeDisc(Stokes, lset, navstokessolver, gm, P, lsetmod);
    if ( nsteps != 0){
        timedisc->SetSchurPrePtr( stokessolverfactory.GetSchurPrePtr() );
    }
    if (P.get<double>("CouplingSolver.NavStokesSolver.Nonlinear")!=0.0 || nsteps == 0) {
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

    if ( nsteps == 0)
        SolveStatProblem( Stokes, lset, *navstokessolver);

    // for serialization of geometry and numerical data
    TwoPhaseStoreCL<InstatNavierStokes2PhaseP2P1CL> ser(MG, Stokes, lset, massTransp,
                                                        P.get<std::string>("Restart.OutputData"),
                                                        P.get<std::string>("Restart.OutputGrid"),
                                                        P.get<int>("Restart.OutputOverwrite"),
                                                        P.get<int>("Restart.Binary"),
                                                        vel_downwind, lset_downwind);
    Stokes.v.t += GetTimeOffset();

    // Output-Registrations:
#ifndef _PAR
    Ensight6OutCL* ensight = NULL;
    if (P.get<int>("Ensight.Freq",0)){
        // Initialize Ensight6 output
        std::string ensf( P.get<std::string>("Ensight.EnsDir") + "/" + P.get<std::string>("Ensight.EnsCase"));
        ensight = new Ensight6OutCL( P.get<std::string>("Ensight.EnsCase") + ".case",
                                     nsteps/P.get("Ensight.Freq", 0)+1,
                                     P.get<int>("Ensight.Binary"));
        ensight->Register( make_Ensight6Geom      ( MG, MG.GetLastLevel(), P.get<std::string>("Ensight.GeomName"),
                                                    ensf + ".geo", true));
        ensight->Register( make_Ensight6Scalar    ( lset.GetSolution(),      "Levelset",      ensf + ".scl", true));
        ensight->Register( make_Ensight6Scalar    ( Stokes.GetPrSolution(),  "Pressure",      ensf + ".pr",  true));
        ensight->Register( make_Ensight6Vector    ( Stokes.GetVelSolution(), "Velocity",      ensf + ".vel", true));
        ensight->Register( make_Ensight6Scalar    ( ScalarFunAsP2EvalCL( sigmap, 0., &MG, MG.GetLastLevel()),
                                                    "Surfaceforce",  ensf + ".sf",  true));

        if (massTransp) {
            ensight->Register( make_Ensight6Scalar( massTransp->GetSolution(),"Concentration", ensf + ".c",   true));
            ensight->Register( make_Ensight6Scalar( massTransp->GetSolution( massTransp->ct),
                                                    "TransConc",     ensf + ".ct",  true));
        }
        if (P.get("SurfTransp.DoTransp", 0)) {
            ensight->Register( make_Ensight6IfaceScalar( MG, surfTransp->ic,  "InterfaceSol",  ensf + ".sur", true));
        }
        if (Stokes.UsesXFEM())
            ensight->Register( make_Ensight6P1XScalar( MG, lset.Phi, Stokes.p, "XPressure",   ensf + ".pr", true));

        ensight->Write( Stokes.v.t);
    }
#endif

    // writer for vtk-format
    VecDescCL * sigma_vtk = NULL;
    sigma_vtk = new VecDescCL;

    VTKOutCL * vtkwriter = NULL;
    if (P.get<int>("VTK.Freq",0)){
        vtkwriter = new VTKOutCL(adap.GetMG(), "DROPS data",
                                 nsteps/P.get("VTK.Freq", 0)+1,
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
            vtkwriter->Register( make_VTKP1XScalar( MG, *lset.PhiC, Stokes.p, "xpressure"));
        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "level-set") );
        vtkwriter->Register( lset.GetVolumeAdjuster()->make_VTKComponentMap("ComponentMap") );

        if (massTransp) {
            vtkwriter->Register( make_VTKScalar( massTransp->GetSolution(), "massTransport") );
        }

        if (P.get("SurfTransp.DoTransp", 0)) {
            vtkwriter->Register( make_VTKIfaceScalar( MG, surfTransp->ic,  "InterfaceSol"));
        }
        vtkwriter->Write(Stokes.v.t);
        vtkwriter->Register( make_VTKScalar( P1EvalCL<double, const StokesPrBndDataCL, const VecDescCL>( sigma_vtk, &Stokes.GetBndData().Pr, &MG), "tau"));
        sf->SetVtkOutput( sigma_vtk);
    }

    VTKOutCL * dgvtkwriter = NULL;
    if ((P.get<int>("Levelset.Discontinuous")&&(P.get<int>("VTK.Freq",0))&&(P.get<int>("VTK.AddDGOutput",0))))
    {
        dgvtkwriter = new VTKOutCL(adap.GetMG(), "DROPS data",
                                   nsteps/P.get("VTK.Freq", 0)+1,
                                   P.get<std::string>("VTK.VTKDir"), P.get<std::string>("VTK.VTKName")+"_dg",
                                   P.get<std::string>("VTK.TimeFileName"),
                                   P.get<int>("VTK.Binary"),
                                   false, /*onlyP1*/
                                   true, /*P2DG*/
                                   -1, /*,-level*/
                                   P.get<int>("VTK.ReUseTimeFile") + "_dg");
        dgvtkwriter->Register( make_VTKScalar( dynamic_cast<LevelsetP2DiscontCL&>(lset).GetDSolution(), "dg-level-set") );
        dgvtkwriter->Write(Stokes.v.t);
    }


    double time = 0.0;

    typedef DistMarkingStrategyCL MarkerT;
    MarkerT marker( lset,
                    P.get<double>("Mesh.AdaptRef.Width"),
                    P.get<double>("Mesh.AdaptRef.CoarsestLevel"), P.get<double>("Mesh.AdaptRef.FinestLevel") );
    adap.set_marking_strategy(&marker);

    IdxDescCL p1idx;

    for (int step= 1; step<=nsteps; ++step)
    {
        std::cout << "============================================================ step " << step << std::endl;
        time += dt;
        const double time_old = Stokes.v.t;
        const double time_new = Stokes.v.t + dt;
        IFInfo.Update( lset, Stokes.GetVelSolution());
        IFInfo.Write(time_old);

        if (P.get<int>("VTK.Freq")) {
            p1idx.CreateNumbering( Stokes.p.RowIdx->TriangLevel(), MG);
            sigma_vtk->SetIdx( &p1idx);
            sigma_vtk->Data = 0.;
        }

        if ( surfTransp ) surfTransp->InitOld();
        timedisc->DoStep( P.get<int>("CouplingSolver.Iter"));
        if (massTransp) massTransp->DoStep( time_new);
        if ( surfTransp ) {
            surfTransp->DoStep( time_new);
            BndDataCL<> ifbnd( 0);
            std::cout << "surfactant on \\Gamma: " << Integral_Gamma( MG, *lset.PhiC, lset.GetBndData(), make_P1Eval(  MG, ifbnd, surfTransp->ic)) << '\n';
        }

        // WriteMatrices( Stokes, step);

        // grid modification
        const bool doGridMod= P.get<int>("Mesh.AdaptRef.Freq") && step%P.get<int>("Mesh.AdaptRef.Freq") == 0;
        bool gridChanged= false;
        if (doGridMod)
        {
            gridChanged = adap.UpdateTriang();

        }
        // downwind-numbering for Navier-Stokes
        const bool doNSDownwindNumbering= P.get<int>("CouplingSolver.NavStokesSolver.Downwind.Frequency")
            && step%P.get<int>("CouplingSolver.NavStokesSolver.Downwind.Frequency") == 0;
        if (doNSDownwindNumbering) {
            if (!gridChanged) { // We must ensure that the permutation maps the original numbering of CreateNumbering to the downwind-numbering and that the renumbering starts in a known state, i.e. the original numbering.
                vidx->DeleteNumbering( MG);
                Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
                permute_Vector( Stokes.v.Data, invert_permutation( vel_downwind), 3);
            }
            vel_downwind= Stokes.downwind_numbering( lset, navstokes_downwind);
        }
        // downwind-numbering for Levelset
        const bool doLsetDownwindNumbering= P.get<int>("CouplingSolver.LevelsetSolver.Downwind.Frequency")
            && step%P.get<int>("CouplingSolver.LevelsetSolver.Downwind.Frequency") == 0;
        if (doLsetDownwindNumbering) {
            if (!gridChanged) { // We must ensure that the permutation maps the original numbering of CreateNumbering to the downwind-numbering and that the renumbering starts in a known state, i.e. the original numbering.
                lset.DeleteNumbering( lidx);
                lset.CreateNumbering( MG.GetLastLevel(), lidx);
                lset.Phi.SetIdx( lidx);
                permute_Vector( lset.Phi.Data, invert_permutation( lset_downwind));
            }
            lset_downwind= lset.downwind_numbering( Stokes.GetVelSolution(), levelset_downwind);
        }
        if (gridChanged || doNSDownwindNumbering || doLsetDownwindNumbering) {
                timedisc->Update();
                if (massTransp) massTransp->Update();
        }

#ifndef _PAR
        if (ensight && step%P.get("Ensight.Freq", 0)==0)
            ensight->Write( time_new);
#endif
        if (dgvtkwriter && step%P.get("VTK.Freq", 0)==0)
            dgvtkwriter->Write( time_new);
        if (vtkwriter && step%P.get("VTK.Freq", 0)==0)
            vtkwriter->Write( time_new);
        if (P.get("Restart.OutputFreq", 0) && step%P.get("Restart.OutputFreq", 0)==0)
            ser.Write();
    }
    IFInfo.Update( lset, Stokes.GetVelSolution());
    IFInfo.Write(Stokes.v.t);
    std::cout << std::endl;
    if(sf) delete sf;
    delete timedisc;
    delete navstokessolver;
    delete stokessolver;
    delete gm;
    delete &lset;
    if (massTransp) delete massTransp;
    if (transprepair) delete transprepair;
#ifndef _PAR
    if (ensight) delete ensight;
#endif
    if (vtkwriter) delete vtkwriter;
    if (dgvtkwriter) delete dgvtkwriter;
    if (sigma_vtk) delete sigma_vtk;
    if (infofile) delete infofile;

//     delete stokessolver1;
}

} // end of namespace DROPS

int main (int argc, char** argv)
{
#ifdef _PAR
    DROPS::ProcCL::Instance(&argc, &argv);
#endif
  try
  {
    std::cout << "Boost version: " << BOOST_LIB_VERSION << std::endl;

    DROPS::read_parameter_file_from_cmdline( P, argc, argv, "../../param/levelset/twophasedrops/risingdroplet.json");
    P.put_if_unset<std::string>("VTK.TimeFileName",P.get<std::string>("VTK.VTKName"));

    std::cout << P << std::endl;

    DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"), P.get<std::vector<std::string> >("General.DynamicLibs") );

    if (P.get<int>("General.ProgressBar"))
        DROPS::ProgressBarTetraAccumulatorCL::Activate();

    // check parameter file
    if ( P.get<double>("NavStokes.Coeff.SurfTens.DilatationalVisco") <
         P.get<double>("NavStokes.Coeff.SurfTens.ShearVisco") )
    {
        throw DROPS::DROPSErrCL("Parameter error : Dilatational viscosity must be larger than surface shear viscosity");
    }

    const std::string perMatchName= P.get( "Mesh.PeriodicBnd.PeriodicMatching", std::string());
    const bool is_periodic = !perMatchName.empty();
    DROPS::match_fun periodic_match = is_periodic ? DROPS::MatchMap::getInstance()[perMatchName] : nullptr;

    DROPS::MultiGridCL* mg= 0;
    typedef DROPS::BndDataCL<DROPS::Point3DCL> VelBndDataCL;
    typedef DROPS::BndDataCL<double>    PrBndDataCL;
    VelBndDataCL *velbnddata = 0;
    PrBndDataCL *prbnddata = 0;
    DROPS::LsetBndDataCL* lsetbnddata= 0;

    try
    {
        std::unique_ptr<DROPS::MGBuilderCL> builder( DROPS::make_MGBuilder( P));
        mg = new DROPS::MultiGridCL( *builder);
    }
    catch (DROPS::DROPSParamErrCL& e)
    {
        std::cout << "\n"
                  << "  /----------------------------------------------------------------\\ \n"
                  << "  | WARNING: It seems you are using the old domain descriptions    | \n"
                  << "  |          or your \"Mesh\" section is not correct.                | \n"
                  << "  |          Please adapt your json-file to the new description.   | \n"
                  <<"  \\----------------------------------------------------------------/ \n"
                  << std::endl;
        DROPS::BuildDomain( mg, P.get<std::string>("DomainCond.MeshFile"), P.get<int>("DomainCond.GeomType"), P.get<std::string>("Mesh.RestartFile",""));
    }
    if (P.exists("Mesh.PeriodicBnd"))
        DROPS::read_PeriodicBoundaries ( *mg, P.get_child("Mesh.PeriodicBnd"));
    const bool noRestart= (P.get<std::string>("Mesh.RestartFile","")).empty() || P.get<std::string>("Mesh.RestartFile","") == "none";

    std::cout << "Generated MG of " << mg->GetLastLevel() << " levels." << std::endl;

    try
    {
        velbnddata=  new VelBndDataCL(0);
        prbnddata=   new PrBndDataCL(0);
        lsetbnddata= new DROPS::LsetBndDataCL(0);
        read_BndData( *velbnddata, *mg, P.get_child( "NavStokes.BoundaryData.Velocity"));
        std::cout << "Generated boundary conditions for velocity, ";
        read_BndData( *prbnddata,  *mg, P.get_child( "NavStokes.BoundaryData.Pressure"));
        std::cout << "pressure, ";
        read_BndData( *lsetbnddata,*mg, P.get_child( "Levelset.BoundaryData"));
        std::cout << "and levelset." << std::endl;
    }
    catch (DROPS::DROPSParamErrCL& e)
    {
        e.what( std::cout);
        delete velbnddata;
        delete prbnddata;
        delete lsetbnddata;
        std::cout << "\n"
                  << "  /----------------------------------------------------------------\\ \n"
                  << "  | WARNING: It seems you are using the old bnd data descriptions  | \n"
                  << "  |          or your \"BoundaryData\" section is not correct.        | \n"
                  << "  |          Please adapt your json-file to the new description.   | \n"
                  <<"  \\----------------------------------------------------------------/ \n"
                  << std::endl;

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
    }
    DROPS::StokesBndDataCL bnddata(*velbnddata,*prbnddata);

    std::string InitialLSet= P.get("Levelset.InitialValue", std::string("Ellipsoid"));
    if (InitialLSet == "Ellipsoid")
        DROPS::EllipsoidCL::Init( P.get<DROPS::Point3DCL>("Levelset.PosDrop"), P.get<DROPS::Point3DCL>("Levelset.RadDrop"));
    if (InitialLSet == "HalfEllipsoid")
        DROPS::HalfEllipsoidCL::Init( P.get<DROPS::Point3DCL>("Levelset.PosDrop"), P.get<DROPS::Point3DCL>("Levelset.RadDrop"));
    if (InitialLSet == "ContactDroplet")
        DROPS::ContactDropletCL::Init( P.get<DROPS::Point3DCL>("Levelset.PosDrop"), P.get<DROPS::Point3DCL>("Levelset.RadDrop"), P.get<double>("Levelset.AngleDrop"));
    if  (InitialLSet == "TwoEllipsoid")
        DROPS::TwoEllipsoidCL::Init( P.get<DROPS::Point3DCL>("Levelset.PosDrop"), P.get<DROPS::Point3DCL>("Levelset.RadDrop"), P.get<DROPS::Point3DCL>("Levelset.PosDrop2"), P.get<DROPS::Point3DCL>("Levelset.RadDrop2"));
    if  (InitialLSet.find("Layer")==0){
        DROPS::LayerCL::Init( P.get<DROPS::Point3DCL>("Levelset.PosDrop"), P.get<DROPS::Point3DCL>("Levelset.RadDrop"), InitialLSet[5]-'X');
        P.put("Exp.InitialLSet", InitialLSet= "Layer");
    }
    if (InitialLSet.find("Cylinder")==0){
        DROPS::CylinderCL::Init( P.get<DROPS::Point3DCL>("Levelset.PosDrop"), P.get<DROPS::Point3DCL>("Levelset.RadDrop"), InitialLSet[8]-'X');
        P.put("Exp.InitialLSet", InitialLSet= "Cylinder");
    }
    typedef DROPS::DistMarkingStrategyCL MarkerT;
    MarkerT InitialMarker( DROPS::InScaMap::getInstance()[InitialLSet],
                           P.get<double>("Mesh.AdaptRef.Width"),
                           P.get<double>("Mesh.AdaptRef.CoarsestLevel"), P.get<double>("Mesh.AdaptRef.FinestLevel") );

    DROPS::AdapTriangCL adap( *mg, &InitialMarker,
                              (noRestart ? P.get<int>("Mesh.AdaptRef.LoadBalStrategy") : -P.get<int>("Mesh.AdaptRef.LoadBalStrategy")));
    // If we read the Multigrid, it shouldn't be modified;
    // otherwise the pde-solutions from the ensight files might not fit.
    if (noRestart)
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

    DROPS::InstatNavierStokes2PhaseP2P1CL prob( *mg, DROPS::TwoPhaseFlowCoeffCL(P), bnddata,
                                                P.get<double>("NavStokes.XFEMReduced")<0 ? DROPS::P1_FE : DROPS::P1X_FE,
                                                P.get<double>("NavStokes.XFEMReduced"), DROPS::vecP2_FE,
                                                P.get<double>("NavStokes.GhostPenalty",0.0));

    Strategy( prob, *lsetbnddata, adap, is_periodic, perMatchName);    // do all the stuff

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

