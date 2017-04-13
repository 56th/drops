/// \file  scalar.cpp
/// \brief Solver for scalar problems with P1 functions or P2 functions
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Eva Loch, Volker Reichelt, Thorolf Schulte, Liang Zhang, Yuanjun Zhang; SC RWTH Aachen: Oliver Fortmeier

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


 // include geometric computing
#include "geom/multigrid.h"             // multigrid on each processor
#include "geom/builder.h"               // construct the initial multigrid
#include "out/output.h"
#include "geom/geomselect.h"
#include "misc/funcmap.h"

#include "geom/deformation.h"

// include numeric computing!
#include "num/fe.h"
#include "num/krylovsolver.h"
#include "num/MGsolver.h"
#include "poisson/integrTime.h"
#include "num/prolongation.h"

 // include problem class
#include "misc/params.h"
#include "poisson/poissonParam.h"      // poissonCoeffCL
#include "poisson/poisson.h"           // setting up the Poisson problem
#include "num/bndData.h"

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>

//include output
#include "out/ensightOut.h"
#include "out/vtkOut.h"

// include parallel computing!
#ifdef _PAR
#include "parallel/parallel.h"          // proc handling, reduce operations, ...
#include "parallel/partime.h"           // parallel time-messurement
#include "parallel/parmultigrid.h"      // handle multigrid over different procs
#include "parallel/loadbal.h"           // distribute multigrid
#endif

// include function container
#include "misc/funcmap.h"
#include "num/poissonsolverfactory.h"
#include "poisson/ale.h"

#include "misc/progressaccu.h"
#include "misc/dynamicload.h"

const char line[] ="----------------------------------------------------------------------------------\n";

DROPS::ParamCL P;   //Parameter class, read in json file in main function

namespace DROPS
{

template<class PoissonCL, class SolverT>
void SolveStatProblem( PoissonCL& Poisson, SolverT& solver, ParamCL& Param)
{
    // time measurements
#ifndef _PAR
    TimerCL timer;
    const bool doErrorEstimate= P.get<int>("Error.DoErrorEstimate");
#else
    const bool doErrorEstimate= false;
    if (Param.get<int>("Error.DoErrorEstimate"))
        std::cout << "Skipping Error-Estimation ..." << std::endl;
    ParTimerCL timer;
#endif

    if ( !doErrorEstimate) {
        // discretize (setup linear equation system)
        std::cout << line << "Discretize (setup linear equation system) in stationary problem...\n";
        timer.Reset();
        Poisson.SetupSystem( Poisson.A, Poisson.b);
        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s" << std::endl;

        //If we need to add convection
        if(Param.get<int>("Poisson.Coeff.withConvection"))
        {
            std::cout << line << "Setup convection...\n";
            timer.Reset();
            Poisson.vU.SetIdx( &Poisson.idx);
            Poisson.SetupConvection(Poisson.U, Poisson.vU, 0.0);
            Poisson.A.Data.LinComb(1., Poisson.A.Data, 1., Poisson.U.Data);
            Poisson.b.Data+=Poisson.vU.Data;
            timer.Stop();
            std::cout << " o time " << timer.GetTime() << " s" << std::endl;
        }
        timer.Reset();
        solver.Solve( Poisson.A.Data, Poisson.x.Data, Poisson.b.Data, Poisson.x.RowIdx->GetEx());
        timer.Stop();

#ifndef _PAR
        double realresid = norm( VectorCL(Poisson.A.Data*Poisson.x.Data-Poisson.b.Data));
#else
        double realresid = Poisson.idx.GetEx().Norm( VectorCL(Poisson.A.Data*Poisson.x.Data-Poisson.b.Data), false);
#endif
        std::cout << " o Solved system with:\n"
                  << "   - time          " << timer.GetTime()   << " s\n"
                  << "   - iterations    " << solver.GetIter()  << '\n'
                  << "   - residuum      " << solver.GetResid() << '\n'
                  << "   - real residuum " << realresid         << std::endl;

        if (P.get<int>("Poisson.SolutionIsKnown")) {
            std::cout << line << "Check result against known solution ...\n";
            timer.Reset();
            Poisson.CheckSolution( Poisson.x, Poisson.Coeff_.Solution);
            timer.Stop();
            std::cout << " o time " << timer.GetTime() << " s" << std::endl;
        }
    }
    else{//Error estimation, only implemented in P1 FEM not P2 FEM
        MultiGridCL& MG= Poisson.GetMG();
        const typename PoissonCL::BndDataCL& BndData= Poisson.GetBndData();

        MLIdxDescCL  loc_idx;
        VecDescCL  loc_x;
        MLIdxDescCL* new_idx= &Poisson.idx;
        MLIdxDescCL* old_idx= &loc_idx;
        VecDescCL* new_x= &Poisson.x;
        VecDescCL* old_x= &loc_x;
        MeshDeformationCL& md= MeshDeformationCL::getInstance();

        DoerflerMarkCL<typename PoissonP1CL<PoissonCoeffCL>::est_fun, typename PoissonP1CL<PoissonCoeffCL>::base_>
            Estimator( P.get<double>("Error.RelReduction"), P.get<double>("Error.MinRatio"), P.get<double>("Error.Threshold"), P.get<double>("Error.Meas"), P.get<int>("Error.DoMark"),
                       &PoissonP1CL<PoissonCoeffCL>::ResidualErrEstimator, *static_cast<typename PoissonP1CL<PoissonCoeffCL>::base_*>(&Poisson) );

        int step= 0;
        bool new_marks;

        new_idx->SetFE( P1_FE);
        old_idx->SetFE( P1_FE);

        do{
            timer.Reset();
            std::cout << DROPS::SanityMGOutCL(MG) << std::endl;
            MG.Refine();

            Poisson.CreateNumbering( MG.GetLastLevel(), new_idx);    // create numbering for this idx
            std::cout << "new triangLevel: " << Poisson.idx.TriangLevel() << std::endl;
            Poisson.b.SetIdx( new_idx);                              // tell b about numbering
            new_x->SetIdx( new_idx);                    			 // second vector with the same idx

            std::cout << line << "Problem size\no number of unknowns             " << new_x->Data.size() << std::endl;

            MG.SizeInfo(std::cout);
            if ( step == 0)
                Estimator.Init( typename PoissonP1CL<PoissonCoeffCL>::DiscSolCL( new_x, &BndData, &MG));

            if ( old_x->RowIdx)
            {
                P1EvalCL<double, const PoissonBndDataCL, const VecDescCL>  oldx( old_x, &BndData, &MG);
                P1EvalCL<double, const PoissonBndDataCL, VecDescCL>        newx( new_x, &BndData, &MG);
                Interpolate( newx, oldx);
          //            CheckSolution(*new_x, &::Lsg);
                old_x->Reset();
             }

            if( Poisson.usesALE())
            {
                md.SetMeshTransformation(PoissonCoeffCL::ALEDeform, -1, P.get<int>("ALE.OnlyBndCurved"), P.get<int>("ALE.P1")==0);
                //ALE.InitGrid();
            }
            else
                md.SetMeshIdentity();

            Poisson.A.SetIdx( new_idx, new_idx);             // tell A about numbering
            Poisson.SetupSystem( Poisson.A, Poisson.b);
            timer.Stop();
            timer.Reset();
            solver.Solve( Poisson.A.Data, new_x->Data, Poisson.b.Data, new_x->RowIdx->GetEx());
            timer.Stop();
            double realresid = norm( VectorCL(Poisson.A.Data*new_x->Data-Poisson.b.Data));
            std::cout << " o Solved system with:\n"
                      << "   - time          " << timer.GetTime()   << " s\n"
                      << "   - iterations    " << solver.GetIter()  << '\n'
                      << "   - residuum      " << solver.GetResid() << '\n'
                      << "   - real residuum " << realresid         << std::endl;
            Poisson.A.Reset();
            Poisson.b.Reset();
            if (P.get<int>("Poisson.SolutionIsKnown")) {
                std::cout << line << "Check result against known solution ...\n";
                Poisson.CheckSolution( *new_x, Poisson.Coeff_.Solution);
            }
            new_marks = Estimator.Estimate( typename PoissonP1CL<PoissonCoeffCL>::const_DiscSolCL( new_x, &BndData, &MG) );

            std::swap( old_x, new_x);
            std::swap( old_idx, new_idx);
        } while ( new_marks && step++ < P.get<int>("Error.NumRef"));
        // I want the solution to be in Poisson.x
        if ( old_x == &loc_x)
        {
            Poisson.idx.swap( loc_idx);
            Poisson.x.SetIdx( &Poisson.idx);

            Poisson.x.Data.resize( loc_x.Data.size());
            Poisson.x.Data = loc_x.Data;
        }
    }
}

/// \brief Strategy to solve the Poisson problem on a given triangulation
template< class PoissonCL>
void Strategy(PoissonCL& Poisson)
{
    // time measurements
#ifndef _PAR
    TimerCL timer;
#else
    ParTimerCL timer;
#endif

    // the triangulation
    MultiGridCL& mg= Poisson.GetMG();

    MeshDeformationCL& md = MeshDeformationCL::getInstance();
    md.Initialize(&mg);

    //ALECL ALE(P, mg);
    // connection triangulation and vectors
    // -------------------------------------------------------------------------
    std::cout << line << "Connecting triangulation and matrices/vectors ...\n";
    timer.Reset();

    if(P.get<int>("Poisson.P1"))
        Poisson.idx.SetFE( P1_FE);
    else
        Poisson.idx.SetFE( P2_FE);
    // set quadratic finite elements
    //see class for explanation: template didnt work
    if ( PoissonSolverFactoryHelperCL().MGUsed(P))
        Poisson.SetNumLvl ( mg.GetNumLevel());
    Poisson.CreateNumbering( mg.GetLastLevel(), &Poisson.idx);  // number vertices and edges
    Poisson.b.SetIdx( &Poisson.idx);                            // tell b about numbering
    Poisson.x.SetIdx( &Poisson.idx);                            // tell x about numbering
    Poisson.A.SetIdx( &Poisson.idx, &Poisson.idx);              // tell A about numbering
    Poisson.M.SetIdx( &Poisson.idx, &Poisson.idx);              // tell M about numbering
    Poisson.U.SetIdx( &Poisson.idx, &Poisson.idx);              // tell U about numbering

    timer.Stop();
    std::cout << " o time " << timer.GetTime() << " s" << std::endl;

    // display problem size
    // -------------------------------------------------------------------------
    std::cout << line << "Problem size\n";
#ifdef _PAR
    std::vector<size_t> UnkOnProc= ProcCL::Gather( Poisson.x.Data.size(), 0);
    const IdxT numUnk   = Poisson.idx.GetGlobalNumUnknowns(),
               numAccUnk= std::accumulate(UnkOnProc.begin(), UnkOnProc.end(), 0);
#else
    std::vector<size_t> UnkOnProc( 1);
    UnkOnProc[0]  = Poisson.x.Data.size();
    IdxT numUnk   = Poisson.x.Data.size(),
         numAccUnk= Poisson.x.Data.size();
#endif
    std::cout << " o number of unknowns             " << numUnk    << '\n'
              << " o number of accumulated unknowns " << numAccUnk << '\n'
              << " o number of unknowns on proc\n";
    for (size_t i=0; i<UnkOnProc.size(); ++i)
        std::cout << " - Proc " << i << ": " << UnkOnProc[i]<< '\n';

    std::cout << line << "Choose the poisson solver...\n";
    timer.Reset();
    // type of preconditioner and solver
    PoissonSolverFactoryCL<> factory( P.get_child("Poisson.Solver"), Poisson.idx);
    PoissonSolverBaseCL* solver = factory.CreatePoissonSolver();

    if ( factory.GetProlongation() != 0)
        SetupProlongationMatrix( mg, *(factory.GetProlongation()), &Poisson.idx, &Poisson.idx);

    timer.Stop();
    std::cout << " o time " << timer.GetTime() << " s" << std::endl;

    //If it is NOT a stationary problem, set up the system and the initial condition
    if (P.get<int>("Time.NumSteps") != 0)
    {
        // discretize (setup linear equation system)
        std::cout << line << "Discretize (setup linear equation system) for instationary problem...\n";
        timer.Reset();
        if (Poisson.usesALE())
        {
            md.SetMeshTransformation(PoissonCoeffCL::ALEDeform, -1, P.get<int>("ALE.OnlyBndCurved"), P.get<int>("ALE.P1")==0);
            //ALE.InitGrid();
        }
        Poisson.SetupInstatSystem( Poisson.A, Poisson.M, Poisson.x.t);
        Poisson.Init( Poisson.x, Poisson.Coeff_.InitialCondition, 0.0);
        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s" << std::endl;
    }
    else
    {//if it is a stationary problem, call SolveStatProblem
        std::cout << line << "Solve the linear equation system ...\n";
        SolveStatProblem( Poisson, *solver, P);
    }

    // Output-Registrations:
#ifndef _PAR
    //Ensight format
    Ensight6OutCL* ensight = NULL;
    if (P.get<int>("Ensight.Freq",0)){
        // Initialize Ensight6 output
        const std::string filename= P.get<std::string>("Ensight.EnsDir") + "/" + P.get<std::string>("Ensight.EnsCase");
        ensight = new Ensight6OutCL(P.get<std::string>("Ensight.EnsCase")+".case", P.get<int>("Time.NumSteps")+1,
                                    P.get<int>("Ensight.Binary"));
        ensight->Register( make_Ensight6Geom  ( mg, mg.GetLastLevel(),
                                                P.get<std::string>("Ensight.GeomName"), filename + ".geo"));
        ensight->Register( make_Ensight6Scalar( Poisson.GetSolution(), "Temperatur", filename + ".tp", true));
        ensight->Write();
    }
#endif
    //VTK format
    VTKOutCL * vtkwriter = NULL;
    if (P.get<int>("VTK.Freq",0)){
        vtkwriter = new VTKOutCL(mg, "DROPS data",
                                 P.get<int>("Time.NumSteps")+1,
                                 P.get<std::string>("VTK.VTKDir"), P.get<std::string>("VTK.VTKName"),
                                 P.get<std::string>("VTK.TimeFileName"),
                                 P.get<int>("VTK.Binary"),
                                 P.get<int>("VTK.UseOnlyP1"),
                                 -1,  /* <- level */
                                 P.get<int>("VTK.ReUseTimeFile"),
                                 P.get<int>("VTK.UseDeformation"));
        vtkwriter->Register( make_VTKScalar( Poisson.GetSolution(), "ConcenT"));
        vtkwriter->Write( Poisson.x.t);
    }
    //Do we have an instationary problem?
    if(P.get<int>("Time.NumSteps")!=0)
    {
        //Creat instationary ThetaschemeCL to handle time integration for instationary problem and set time steps
        InstatPoissonThetaSchemeCL<PoissonCL, PoissonSolverBaseCL>
        ThetaScheme( Poisson, *solver, P);
        ThetaScheme.SetTimeStep(P.get<double>("Time.FinalTime")/P.get<int>("Time.NumSteps") );
        //Solve linear systerm in each time step
        for ( int step = 1; step <= P.get<int>("Time.NumSteps") ; ++step) {
            timer.Reset();
            std::cout << line << "Step: " << step << std::endl;
            if (Poisson.usesALE())
            {
                md.SetMeshTransformation(PoissonCoeffCL::ALEDeform, Poisson.x.t, P.get<int>("ALE.OnlyBndCurved"), P.get<int>("ALE.P1")==0 );
                //ALE.MovGrid(Poisson.x.t);
            }
            else
                md.SetMeshIdentity();
            ThetaScheme.DoStep( Poisson.x);

            timer.Stop();
            std::cout << " o Solved system with:\n"
                      << "   - time          " << timer.GetTime()    << " s\n"
                      << "   - iterations    " << solver->GetIter()  << '\n'
                      << "   - residuum      " << solver->GetResid() << '\n';

            if (P.get("Poisson.SolutionIsKnown", 0)) {
                std::cout << " o Check result against known solution ...\n";
                timer.Reset();
                Poisson.CheckSolution( Poisson.x, Poisson.Coeff_.Solution, Poisson.x.t);
                timer.Stop();
                std::cout << " o -time " << timer.GetTime() << " s" << std::endl;
            }

#ifndef _PAR
            if (ensight && step%P.get<int>("Ensight.Freq", 0)==0){
                std::cout << " o Ensight output ...\n";
                timer.Reset();
                ensight->Write( Poisson.x.t);
                timer.Stop();
                std::cout << " o -time " << timer.GetTime() << " s" << std::endl;
            }
#endif
            if (vtkwriter && step%P.get<int>("VTK.Freq", 0)==0){
                std::cout << " o VTK output ...\n";
                timer.Reset();
                vtkwriter->Write( Poisson.x.t);
                timer.Stop();
                std::cout << " o -time " << timer.GetTime() << " s" << std::endl;
            }
        }
    }

    if (vtkwriter) delete vtkwriter;
#ifndef _PAR
    if (ensight) delete ensight;
#endif
    delete solver;
}

} // end of namespace DROPS


int main (int argc, char** argv)
{
#ifdef _PAR
    DROPS::ProcCL::Instance(&argc, &argv);
#endif
    try
    {
        // time measurements
#ifndef _PAR
        DROPS::TimerCL timer;
#else
        DROPS::ParTimerCL timer;
#endif

        DROPS::read_parameter_file_from_cmdline( P, argc, argv, "../../param/poisson/cdrdrops/statpoissonEx.json");
        P.put_if_unset<std::string>("VTK.TimeFileName",P.get<std::string>("VTK.VTKName"));
        //output all the parameters
        std::cout << P << std::endl;

        DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"),
                           P.get<std::vector<std::string> >("General.DynamicLibs") );

        if (P.get<int>("General.ProgressBar"))
            DROPS::ProgressBarTetraAccumulatorCL::Activate();
        // set up data structure to represent a poisson problem
        // ---------------------------------------------------------------------
        std::cout << line << "Set up data structure to represent a Poisson problem ...\n";
        timer.Reset();

        //create geometry
        DROPS::MultiGridCL* mg= 0;
        DROPS::PoissonBndDataCL* bdata = new DROPS::PoissonBndDataCL(0);

        //build computational domain
        std::unique_ptr<DROPS::MGBuilderCL> builder( DROPS::make_MGBuilder( P));
        mg = new DROPS::MultiGridCL( *builder);
        //Setup boundary conditions
        read_BndData( *bdata, *mg, P.get_child( "Poisson.BoundaryData"));
        for (int i=0; i<mg->GetBnd().GetNumBndSeg(); ++i)
            std::cout << i << ": BC = " << bdata->GetBndSeg(i).GetBC() << std::endl;
        //Initialize SUPGCL class
        DROPS::SUPGCL supg;
        //SUPG stabilization, ALE method  and error estimation are not yet implemented for P2 case!
        if(!P.get<int>("Poisson.P1"))
        {
              P.put<int>("Stabilization.SUPG",0);
              P.put<int>("Error.DoErrorEstimate",0);
        }
        if(P.get<int>("Stabilization.SUPG"))
        {
            supg.init(P);
            std::cout << line << "The SUPG stabilization will be added ...\n"<<line;
        }
        // Setup the problem
        DROPS::PoissonCoeffCL tmp = DROPS::PoissonCoeffCL( P);

        DROPS::PoissonP1CL<DROPS::PoissonCoeffCL> *probP1 = 0;
        DROPS::PoissonP2CL<DROPS::PoissonCoeffCL> *probP2 = 0;
        if(P.get<int>("Poisson.P1"))
            probP1 = new DROPS::PoissonP1CL<DROPS::PoissonCoeffCL>( *mg, tmp, *bdata, supg, P.get<int>("ALE.wavy"));
        else
        {
            probP2 = new DROPS::PoissonP2CL<DROPS::PoissonCoeffCL>( *mg, tmp, *bdata, P.get<int>("ALE.wavy"));
        }

#ifdef _PAR
        // Set parallel data structures
        DROPS::LoadBalCL lb( *mg);                    // loadbalancing
        lb.DoMigration();    // distribute initial grid
#endif

        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s" << std::endl;

        // Refine the grid
        std::cout << "Refine the grid " << P.get<int>("Mesh.AdaptRef.FinestLevel") << " times regulary ...\n";
        timer.Reset();
        // Create new tetrahedra
        for ( int ref=1; ref <= P.get<int>("Mesh.AdaptRef.FinestLevel"); ++ref){
            std::cout << " refine (" << ref << ")\n";
            DROPS::MarkAll( *mg);
            mg->Refine();
            // do loadbalancing
    #ifdef _PAR
            lb.DoMigration();
    #endif
        }

        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s" << std::endl;
        mg->SizeInfo( std::cout);

        // Solve the problem
        if(P.get<int>("Poisson.P1"))
           DROPS::Strategy<DROPS::PoissonP1CL<DROPS::PoissonCoeffCL> >(*probP1);
        else
           DROPS::Strategy<DROPS::PoissonP2CL<DROPS::PoissonCoeffCL> >(*probP2);
        //Check if Multigrid is sane
        std::cout << line << "Check if multigrid works properly...\n";
        timer.Reset();
        if(P.get<int>("ALE.wavy"))
            std::cout << "Because of ALE method, we don't check the sanity of multigrid here!" << std::endl;
        else
            std::cout << DROPS::SanityMGOutCL(*mg) << std::endl;
        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s" << std::endl;
        // delete dynamically allocated objects
        delete mg;
        delete bdata;
        delete probP1;
        delete probP2;
        std::cout << "cdrdrops finished regularly" << std::endl;
        return 0;
    }
    catch (DROPS::DROPSErrCL& err) { err.handle(); }
}
