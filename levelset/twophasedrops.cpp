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
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
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
//levelset
#include "levelset/coupling.h"
#include "levelset/adaptriang.h"
#include "levelset/mzelle_hdr.h"
#include "levelset/twophaseutils.h"
//surfactants
#include "surfactant/ifacetransp.h"
//function map
#include "misc/bndmap.h"
//solver factory for stokes
#include "num/stokessolverfactory.h"
#ifndef _PAR
#include "num/stokessolver.h"
#else
#include "num/parstokessolver.h"
#include "parallel/loadbal.h"
#include "parallel/parmultigrid.h"
#endif
//general: streams
#include <fstream>
#include <sstream>

DROPS::ParamCL P;

// rho*du/dt - mu*laplace u + Dp = f + rho*g - okn
//                        -div u = 0
//                             u = u0, t=t0

namespace DROPS // for Strategy
{

double GetTimeOffset(){
    double timeoffset = 0.0;
    const std::string restartfilename = P.get<std::string>("Restart.Inputfile");
    if (restartfilename != "none"){
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
    DROPS::InScaMap & inscamap = DROPS::InScaMap::getInstance();
    //DROPS::ScaMap & scamap = DROPS::ScaMap::getInstance();
    //DROPS::InVecMap & vecmap = DROPS::InVecMap::getInstance();
    DROPS::MatchMap & matchmap = DROPS::MatchMap::getInstance();


    bool is_periodic = P.get<std::string>("DomainCond.PeriodicMatching", "none") != "none";
    match_fun periodic_match = is_periodic ? matchmap[P.get("DomainCond.PeriodicMatching", std::string("periodicx"))] : 0;

    MultiGridCL& MG= Stokes.GetMG();

    // initialization of surface tension
    sigma= Stokes.GetCoeff().SurfTens;
    eps= P.get<double>("SurfTens.JumpWidth");    lambda= P.get<double>("SurfTens.RelPos");    sigma_dirt_fac= P.get<double>("SurfTens.DirtFactor");
    instat_scalar_fun_ptr sigmap  = 0;
    if (P.get<double>("SurfTens.VarTension"))
    {
        sigmap  = &sigma_step;
    }
    else
    {
        sigmap  = &sigmaf;
    }
    SurfaceTensionCL sf( sigmap);

    LevelsetP2CL lset( MG, lsetbnddata, sf, P.get<double>("Levelset.SD"), P.get<double>("Levelset.CurvDiff"));

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
    }

    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);
    VelocityRepairCL velrepair( Stokes);
    adap.push_back( &velrepair);
    PressureRepairCL prrepair( Stokes, lset);
    adap.push_back( &prrepair);

    IdxDescCL* lidx= &lset.idx;
    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;

    lset.CreateNumbering( MG.GetLastLevel(), lidx, periodic_match);
    lset.Phi.SetIdx( lidx);
    PermutationT lset_downwind;
    if (P.get<double>("SurfTens.VarTension"))
        lset.SetSurfaceForce( SF_ImprovedLBVar);
    else
        lset.SetSurfaceForce( SF_ImprovedLB);

    if ( StokesSolverFactoryHelperCL().VelMGUsed(P))
        Stokes.SetNumVelLvl ( Stokes.GetMG().GetNumLevel());
    if ( StokesSolverFactoryHelperCL().PrMGUsed(P))
        Stokes.SetNumPrLvl  ( Stokes.GetMG().GetNumLevel());

    SetInitialLevelsetConditions( lset, MG, P);
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

    TransportP1CL * massTransp = NULL;
    TransportRepairCL *  transprepair = NULL;

    if (P.get<int>("Transp.DoTransp"))
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

        massTransp = new TransportP1CL( MG, Bnd_c, Stokes.GetBndData().Vel, P.get<double>("Transp.Theta"),
                                  D, P.get<double>("Transp.HNeg")/P.get<double>("Transp.HPos"), &Stokes.v, lset,
                                  P.get<double>("Time.StepSize"), P.get<int>("Transp.Iter"), P.get<double>("Transp.Tol"));

        transprepair = new TransportRepairCL(*massTransp, MG);
        adap.push_back(transprepair);

        MLIdxDescCL* cidx= &massTransp->idx;
        massTransp->CreateNumbering( MG.GetLastLevel(), cidx);
        massTransp->ct.SetIdx( cidx);
        if (P.get<int>("DomainCond.InitialCond") != -1)
            massTransp->Init( inscamap["Initialcneg"], inscamap["Initialcpos"]);
        else
            ReadFEFromFile( massTransp->ct, MG, P.get<std::string>("DomainCond.InitialFile")+"concentrationTransf");

        massTransp->Update();
        std::cout << massTransp->c.Data.size() << " concentration unknowns,\n";
    }

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
    StokesSolverBaseCL* stokessolver = stokessolverfactory.CreateStokesSolver();
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
    GSPcCL lset_pc;
    // GMResSolverCL<SSORPcCL>* gm = new GMResSolverCL<SSORPcCL>( lset_pc, 100, P.get<int>("Levelset.Iter"), P.get<double>("Levelset.Tol"));
    GMResSolverCL<GSPcCL>* gm = new GMResSolverCL<GSPcCL>( lset_pc, 200, P.get<int>("Levelset.Iter"), P.get<double>("Levelset.Tol"));
#else
    ParJac0CL jacparpc( *lidx);
    ParPreGMResSolverCL<ParJac0CL>* gm = new ParPreGMResSolverCL<ParJac0CL>
           (/*restart*/100, P.get<int>("Levelset.Iter"), P.get<double>("Levelset.Tol"), *lidx, jacparpc,/*rel*/true, /*acc*/ true, /*modGS*/false, LeftPreconditioning, /*parmod*/true);
#endif

    LevelsetModifyCL lsetmod( P.get<int>("Reparam.Freq"), P.get<int>("Reparam.Method"), P.get<double>("Reparam.MaxGrad"), P.get<double>("Reparam.MinGrad"), P.get<int>("Levelset.VolCorrection"), Vol, is_periodic);

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

    UpdateProlongationCL PVel( Stokes.GetMG(), stokessolverfactory.GetPVel(), &Stokes.vel_idx, &Stokes.vel_idx);
    adap.push_back( &PVel);
    UpdateProlongationCL PPr ( Stokes.GetMG(), stokessolverfactory.GetPPr(), &Stokes.pr_idx, &Stokes.pr_idx);
    adap.push_back( &PPr);
    // For a two-level MG-solver: P2P1 -- P2P1X;
//     MakeP1P1XProlongation ( Stokes.vel_idx.NumUnknowns(), Stokes.pr_idx.NumUnknowns(),
//         Stokes.pr_idx.GetFinest().GetXidx().GetNumUnknownsStdFE(),
//         stokessolverfactory.GetPVel()->GetFinest(), stokessolverfactory.GetPPr()->GetFinest());

    std::ofstream* infofile = 0;
    IF_MASTER {
        infofile = new std::ofstream ((P.get<std::string>("VTK.VTKName","twophasedrops")+".info").c_str());
    }
    IFInfo.Init(infofile);
    IFInfo.WriteHeader();

    if (P.get<int>("Time.NumSteps") == 0)
        SolveStatProblem( Stokes, lset, *navstokessolver);

    // for serialization of geometry and numerical data
    TwoPhaseStoreCL<InstatNavierStokes2PhaseP2P1CL> ser(MG, Stokes, lset, massTransp,
                                                        P.get<std::string>("Restart.Outputfile"),
                                                        P.get<int>("Restart.Overwrite"),
                                                        P.get<int>("Restart.Binary"),
                                                        vel_downwind, lset_downwind);
    Stokes.v.t += GetTimeOffset();
    // Output-Registrations:
    Ensight6OutCL* ensight = NULL;
    if (P.get<int>("Ensight.EnsightOut",0)){
        // Initialize Ensight6 output
        std::string ensf( P.get<std::string>("Ensight.EnsDir") + "/" + P.get<std::string>("Ensight.EnsCase"));
        ensight = new Ensight6OutCL( P.get<std::string>("Ensight.EnsCase") + ".case",
                                     P.get<int>("Time.NumSteps")/P.get("Ensight.EnsightOut", 0)+1,
                                     P.get<int>("Ensight.Binary"), P.get<int>("Ensight.MasterOut"));
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
            ensight->Register( make_Ensight6IfaceScalar( MG, surfTransp.ic,  "InterfaceSol",  ensf + ".sur", true));
        }

#ifndef _PAR
        if (Stokes.UsesXFEM())
            ensight->Register( make_Ensight6P1XScalar( MG, lset.Phi, Stokes.p, "XPressure",   ensf + ".pr", true));
#endif
        ensight->Write( Stokes.v.t);
    }

    // writer for vtk-format
    VTKOutCL * vtkwriter = NULL;
    if (P.get<int>("VTK.VTKOut",0)){
        vtkwriter = new VTKOutCL(adap.GetMG(), "DROPS data",
                                 P.get<int>("Time.NumSteps")/P.get("VTK.VTKOut", 0)+1,
                                 P.get<std::string>("VTK.VTKDir"), P.get<std::string>("VTK.VTKName"),
                                 P.get<int>("VTK.Binary"));
        vtkwriter->Register( make_VTKVector( Stokes.GetVelSolution(), "velocity") );
        vtkwriter->Register( make_VTKScalar( Stokes.GetPrSolution(), "pressure") );
        if (P.get<int>("VTK.AddP1XPressure",0))
            vtkwriter->Register( make_VTKP1XScalar( MG, lset.Phi, Stokes.p, "xpressure"));
        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "level-set") );

        if (massTransp) {
            vtkwriter->Register( make_VTKScalar( massTransp->GetSolution(), "massTransport") );
        }

        if (P.get("SurfTransp.DoTransp", 0)) {
            vtkwriter->Register( make_VTKIfaceScalar( MG, surfTransp.ic,  "InterfaceSol"));
        }
        vtkwriter->Write(Stokes.v.t);
    }

    // if (P.get("Restart.Serialization", 0))
    //     ser.Write();

    const int nsteps = P.get<int>("Time.NumSteps");
    const double dt = P.get<double>("Time.StepSize");
    for (int step= 1; step<=nsteps; ++step)
    {
        std::cout << "============================================================ step " << step << std::endl;
        const double time_old = Stokes.v.t;
        const double time_new = Stokes.v.t + dt;
        IFInfo.Update( lset, Stokes.GetVelSolution());
        IFInfo.Write(time_old);

        if (P.get("SurfTransp.DoTransp", 0)) surfTransp.InitOld();
        timedisc->DoStep( P.get<int>("Coupling.Iter"));
        if (massTransp) massTransp->DoStep( time_new);
        if (P.get("SurfTransp.DoTransp", 0)) {
            surfTransp.DoStep( time_new);
            BndDataCL<> ifbnd( 0);
            std::cout << "surfactant on \\Gamma: " << Integral_Gamma( MG, lset.Phi, lset.GetBndData(), make_P1Eval(  MG, ifbnd, surfTransp.ic)) << '\n';
        }

        // WriteMatrices( Stokes, step);

        // grid modification
        const bool doGridMod= P.get<int>("AdaptRef.Freq") && step%P.get<int>("AdaptRef.Freq") == 0;
        bool gridChanged= false;
        if (doGridMod) {
            adap.UpdateTriang( lset);
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
                if (massTransp) massTransp->Update();
        }

        if (ensight && step%P.get("Ensight.EnsightOut", 0)==0)
            ensight->Write( time_new);
        if (vtkwriter && step%P.get("VTK.VTKOut", 0)==0)
            vtkwriter->Write( time_new);
        if (P.get("Restart.Serialization", 0) && step%P.get("Restart.Serialization", 0)==0)
            ser.Write();
    }
    IFInfo.Update( lset, Stokes.GetVelSolution());
    IFInfo.Write(Stokes.v.t);
    std::cout << std::endl;
    delete timedisc;
    delete navstokessolver;
    delete stokessolver;
    delete gm;
    if (massTransp) delete massTransp;
    if (transprepair) delete transprepair;
    if (ensight) delete ensight;
    if (vtkwriter) delete vtkwriter;
    if (infofile) delete infofile;
//     delete stokessolver1;
}

} // end of namespace DROPS


/// \brief Set Default parameters here s.t. they are initialized.
/// The result can be checked when Param-list is written to the output.
void SetMissingParameters(DROPS::ParamCL& P){
    P.put_if_unset<int>("Transp.DoTransp",0);
    P.put_if_unset<std::string>("Restart.Inputfile","none");
    P.put_if_unset<int>("NavStokes.Downwind.Frequency", 0);
    P.put_if_unset<double>("NavStokes.Downwind.MaxRelComponentSize", 0.05);
    P.put_if_unset<double>("NavStokes.Downwind.WeakEdgeRatio", 0.2);
    P.put_if_unset<double>("NavStokes.Downwind.CrosswindLimit", std::cos( M_PI/6.));
    P.put_if_unset<int>("Levelset.Downwind.Frequency", 0);
    P.put_if_unset<double>("Levelset.Downwind.MaxRelComponentSize", 0.05);
    P.put_if_unset<double>("Levelset.Downwind.WeakEdgeRatio", 0.2);
    P.put_if_unset<double>("Levelset.Downwind.CrosswindLimit", std::cos( M_PI/6.));
}

int main (int argc, char** argv)
{
#ifdef _PAR
    DROPS::ProcInitCL procinit(&argc, &argv);
#endif
  try
  {
#ifdef _PAR
    DROPS::ParMultiGridInitCL pmginit;
#endif
    std::ifstream param;

    std::cout << "Boost version: " << BOOST_LIB_VERSION << std::endl;

    if (argc!=2)
    {
        std::cout << "Using default parameter file: risingdroplet.json\n";
        param.open( "risingdroplet.json");
    }
    else
        param.open( argv[1]);
    if (!param)
    {
        std::cerr << "error while opening parameter file\n";
        return 1;
    }
    param >> P;
    param.close();

    SetMissingParameters(P);

    std::cout << P << std::endl;

    DROPS::MatchMap & matchmap = DROPS::MatchMap::getInstance();
    bool is_periodic = P.get<std::string>("DomainCond.PeriodicMatching", "none") != "none";
    DROPS::match_fun periodic_match = is_periodic ? matchmap[P.get<std::string>("DomainCond.PeriodicMatching", "periodicx")] : 0;

    DROPS::MultiGridCL* mg= 0;
    typedef DROPS::BndDataCL<DROPS::Point3DCL> VelBndDataCL;
    typedef DROPS::BndDataCL<double>    PrBndDataCL;
    VelBndDataCL *velbnddata = 0;
    PrBndDataCL *prbnddata = 0;
    DROPS::LsetBndDataCL* lsetbnddata= 0;

    //you cannot pass a double& per P.get, so you need to use this indirect way
    double ExpRadInlet = P.get<double>("Exp.RadInlet");
    DROPS::BuildDomain( mg, P.get<std::string>("DomainCond.MeshFile"), P.get<int>("DomainCond.GeomType"), P.get<std::string>("Restart.Inputfile"), ExpRadInlet);
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

    std::string InitialLSet= P.get("Exp.InitialLSet", std::string("Ellipsoid"));
    if (InitialLSet == "Ellipsoid")
        DROPS::EllipsoidCL::Init( P.get<DROPS::Point3DCL>("Exp.PosDrop"), P.get<DROPS::Point3DCL>("Exp.RadDrop"));
    if  (InitialLSet == "TwoEllipsoid")
        DROPS::TwoEllipsoidCL::Init( P.get<DROPS::Point3DCL>("Exp.PosDrop"), P.get<DROPS::Point3DCL>("Exp.RadDrop"), P.get<DROPS::Point3DCL>("Exp.PosDrop2"), P.get<DROPS::Point3DCL>("Exp.RadDrop2"));
    if (InitialLSet.find("Cylinder")==0) {
        DROPS::CylinderCL::Init( P.get<DROPS::Point3DCL>("Exp.PosDrop"), P.get<DROPS::Point3DCL>("Exp.RadDrop"), InitialLSet[8]-'X');
        P.put("Exp.InitialLSet", InitialLSet= "Cylinder");
    }
    DROPS::AdapTriangCL adap( *mg, P.get<double>("AdaptRef.Width"), P.get<int>("AdaptRef.CoarsestLevel"), P.get<int>("AdaptRef.FinestLevel"),
                              ((P.get<std::string>("Restart.Inputfile") == "none") ? P.get<int>("AdaptRef.LoadBalStrategy") : -P.get<int>("AdaptRef.LoadBalStrategy")),
                              P.get<int>("AdaptRef.Partitioner"));
    // If we read the Multigrid, it shouldn't be modified;
    // otherwise the pde-solutions from the ensight files might not fit.
    if (P.get("Restart.Inputfile", std::string("none")) == "none")
        adap.MakeInitialTriang( * DROPS::ScaMap::getInstance()[InitialLSet]);

    std::cout << DROPS::SanityMGOutCL(*mg) << std::endl;
#ifdef _PAR
    adap.GetLb().GetLB().SetWeightFnct(1);
    if (DROPS::ProcCL::Check( CheckParMultiGrid( adap.GetPMG())))
        std::cout << "As far as I can tell the ParMultigridCl is sane\n";
#endif

    DROPS::InstatNavierStokes2PhaseP2P1CL prob( *mg, DROPS::TwoPhaseFlowCoeffCL(P), bnddata, P.get<double>("Stokes.XFEMStab")<0 ? DROPS::P1_FE : DROPS::P1X_FE, P.get<double>("Stokes.XFEMStab"));

    Strategy( prob, *lsetbnddata, adap);    // do all the stuff

    delete mg;
    delete velbnddata;
    delete prbnddata;
    delete lsetbnddata;
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}

