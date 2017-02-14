/// \file sdropsP2.cpp
/// \brief stokes problem
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Eva Loch, Volker Reichelt, Yuanjun Zhang; SC RWTH Aachen: Oliver Fortmeier

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
 * Copyright 2010 LNM/SC RWTH Aachen, Germany
*/

 // include geometric computing
#include "geom/multigrid.h"             // multigrid on each processor
#include "geom/builder.h"               // construct the initial multigrid
#include "geom/geomselect.h"

//output
#include "out/output.h"
#include "out/ensightOut.h"
#include "out/vtkOut.h"

 // include numeric computing!
#include "num/fe.h"
#include "num/krylovsolver.h"
#include "num/MGsolver.h"
#include "stokes/integrTime.h"
#include "num/stokessolverfactory.h"

 // include problem class
#include "stokes/stokes.h"
#include "num/bndData.h"

//include coefficient class
#include "stokes/stokesCoeff.h"

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>

#include "misc/params.h"
#include "out/ensightOut.h"
#include "misc/funcmap.h"

#include "misc/progressaccu.h"
#include "misc/dynamicload.h"

using namespace std;

const char line[] ="----------------------------------------------------------------------------------\n";

void MarkLower( DROPS::MultiGridCL& mg, double tresh)
{
    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin()),
             ItEnd(mg.GetTriangTetraEnd()); It!=ItEnd; ++It)
    {
        if (GetBaryCenter(*It)[2]<=tresh )
            It->SetRegRefMark();
    }
}

const double radiusorbit= 0.3; // Radius of the drops' orbit.
const double radiusdrop= 0.15; // Initial radius of the drop.

// positive outside the drop, negative inside the drop.
double
SignedDistToInterface(const DROPS::Point3DCL& p, double t)
{
   DROPS::Point3DCL c;
   c[0]= 0.5 + radiusorbit*std::cos( 2.*M_PI*t);
   c[1]= 0.5 + radiusorbit*std::sin( 2.*M_PI*t);
   c[2]= 0.5;
   return (p-c).norm() - radiusdrop;
}

typedef double (*signed_dist_fun)(const DROPS::Point3DCL& p, double t);

bool
ModifyGridStep(DROPS::MultiGridCL& mg,
               const signed_dist_fun Dist,
               const double width,         // Thickness of refined shell on each side of the interface
               const DROPS::Uint c_level,  // Outside the shell, use this level
               const DROPS::Uint f_level,  // Inside the shell, use this level
               const double t)             // Time of evaluation
// One step of grid change; returns true if modifications were necessary,
// false, if nothing changed.
{
    using namespace DROPS;
    bool shell_not_ready= false;
        for (MultiGridCL::TriangTetraIteratorCL it= mg.GetTriangTetraBegin(),
             theend= mg.GetTriangTetraEnd(); it!=theend; ++it) {
            double d= 1.;
            for (Uint j=0; j<4; ++j) {
                d= std::min( d, std::abs( Dist( it->GetVertex( j)->GetCoord(), t)));
            }
            const Uint l= it->GetLevel();
            if (d<=width) { // In the shell; level should be f_level.
                if (l < f_level) {
                    shell_not_ready= true;
                    it->SetRegRefMark();
                }
                else
                    if (l > f_level) { it->SetRemoveMark(); }
                    else {} // nothing
            }
            else { // Outside the shell; level should be c_level;
                if (l < c_level) { it->SetRegRefMark(); }
                else
                    if (l> c_level) { it->SetRemoveMark(); }
                    else {} // nothing
            }
        }
        mg.Refine();
    return shell_not_ready;
}


void
MakeInitialTriangulation(DROPS::MultiGridCL& mg,
                         const signed_dist_fun Dist,
                         const double width,         // Thickness of refined shell on eache side of the interface
                         const DROPS::Uint c_level,  // Outside the shell, use this level
                         const DROPS::Uint f_level)  // Inside the shell, use this level
{
    using namespace DROPS;
    Assert( c_level<=f_level, "MakeInitialTriangulation: Levels are cheesy.\n", ~0);
    TimerCL time;

    time.Reset();
    time.Start();
    bool shell_not_ready= true;
    const Uint min_ref_num= f_level - c_level;
    Uint i;
    for(i=0; shell_not_ready || i<min_ref_num; ++i)
        shell_not_ready=  ModifyGridStep( mg, Dist, width, c_level, f_level, 0.);
    time.Stop();
    std::cout << "MakeTriangulation: " << i
              << " refinements in " << time.GetTime() << " seconds\n"
              << "last level: " << mg.GetLastLevel() << '\n';
    mg.SizeInfo( std::cout);
}

typedef DROPS::StokesP2P1CL<DROPS::StokesFlowCoeffCL>
        StokesOnBrickCL;
typedef StokesOnBrickCL MyStokesCL;

DROPS::ParamCL P;

namespace DROPS // for Strategy
{

using ::MyStokesCL;

template <class StokesProblemT>
void ComputeProjectedDivergence(VecDescCL& divu, const StokesProblemT& Stokes)
{
    VectorCL rhs( divu.Data.size() );
    rhs = Stokes.B.Data.GetFinest() * Stokes.v.Data - Stokes.c.Data;

    JACPcCL jacpc;
    PCGSolverCL<JACPcCL> pcgsolver( jacpc, 100, 1e-6, true );

    int numiter;
    double resid;
    pcgsolver.Solve( Stokes.prM.Data.GetFinest(), divu.Data, rhs, DummyExchangeCL(), numiter, resid );

    std::cout << "divergence projection: numiter: " << numiter << std::endl;
    std::cout << "divergence projection:   reisd: " << resid << std::endl;

}

template <class StokesProblemT>
void SolveStatProblem( StokesProblemT& Stokes, StokesSolverBaseCL& solver)
{
    TimerCL timer;
    timer.Reset();

    timer.Start();
    solver.Solve( Stokes.A.Data, Stokes.B.Data, Stokes.C.Data, Stokes.v.Data, Stokes.p.Data, Stokes.b.Data, Stokes.c.Data, Stokes.v.RowIdx->GetEx(), Stokes.p.RowIdx->GetEx());
    timer.Stop();
    const double duration = timer.GetTime();
    std::cout << "Solving Stokes took "<<  duration << " sec.\n";
    std::cout << "iter: " << solver.GetIter() << "\tresid: " << solver.GetResid() << std::endl;

    if( P.get<std::string>("NavStokes.Coeff.Solution_Vel").compare("None")!=0)  // check whether solution is given
        Stokes.CheckSolution( &Stokes.v, &Stokes.p, StokesFlowCoeffCL::LsgVel, StokesFlowCoeffCL::DLsgVel, StokesFlowCoeffCL::LsgPr, true);

}

template< class StokesProblemT>
void Strategy( StokesProblemT& Stokes)
// flow control
{
    //Timer function
    TimerCL timer;

    //the triangulation
    MultiGridCL& MG= Stokes.GetMG();

    // connection triangulation and vectors
    // -------------------------------------------------------------------------
    std::cout << line << "Connecting triangulation and matrices/vectors ...\n";
    timer.Reset();

    Stokes.vel_idx.SetFE( vecP2_FE);
    Stokes.pr_idx.SetFE( P1_FE);

    //Modify Triangulation
    if( P.get("Mesh.AdaptRef.ModifyGrid", 0) == 1)
        MakeInitialTriangulation( MG, &SignedDistToInterface, P.get<double>("Mesh.AdaptRef.Width"), P.get<int>("Mesh.AdaptRef.CoarsestLevel"), P.get<int>("Mesh.AdaptRef.FinestLevel"));

    ParamCL PSolver( P.get_child("NavStokesSolver.OseenSolver") );
    if( StokesSolverFactoryHelperCL().VelMGUsed(PSolver))
    	Stokes.SetNumVelLvl( MG.GetNumLevel());
    if( StokesSolverFactoryHelperCL().PrMGUsed(PSolver))
        Stokes.SetNumPrLvl( MG.GetNumLevel());

    Stokes.CreateNumberingVel( MG.GetLastLevel(), &Stokes.vel_idx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), &Stokes.pr_idx);

    Stokes.SetIdx();
    Stokes.v.SetIdx( &Stokes.vel_idx);
    Stokes.p.SetIdx( &Stokes.pr_idx);

    timer.Stop();
    std::cout << " o time " << timer.GetTime() << " s" << std::endl;

    // display problem size
    // -------------------------------------------------------------------------
    std::cout << line << "Problem size\no number of velocity unknowns             " << Stokes.v.Data.size() << std::endl;
    std::cout << "o number of pressure unknowns             " << Stokes.p.Data.size() << std::endl;

    // discretize (setup linear equation system)
    // -------------------------------------------------------------------------
    std::cout << line << "Discretize (setup linear equation system) ...\n";

    timer.Reset();
    VelVecDescCL  cplM( &Stokes.vel_idx);
    Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &cplM, 0.0);
    Stokes.SetupSystem2( &Stokes.B, &Stokes.c, 0.0);
    Stokes.SetupPrMass ( &Stokes.prM);
    Stokes.SetupPrStiff( &Stokes.prA);
    timer.Stop();
    std::cout << " o time " << timer.GetTime() << " s" << std::endl;

    // solve the linear equation system
    // -------------------------------------------------------------------------
    std::cout << line << "Solve the linear equation system ...\n";

    // type of preconditioner and solver    
    ParamCL PTime( P.get_child("Time") );
    StokesSolverFactoryCL< StokesProblemT>         factory( Stokes, PSolver, PTime );
    StokesSolverBaseCL* stokessolver = factory.CreateStokesSolver();

    if( StokesSolverFactoryHelperCL().VelMGUsed(PSolver))
        SetupProlongationMatrix( MG, *factory.GetPVel(), &Stokes.vel_idx, &Stokes.vel_idx);

    if( StokesSolverFactoryHelperCL().PrMGUsed(PSolver))
        SetupProlongationMatrix( MG, *factory.GetPPr(), &Stokes.pr_idx, &Stokes.pr_idx);

    // choose time discretization scheme
    TimeDiscStokesCL< StokesProblemT,  StokesSolverBaseCL>* TimeScheme;
    switch ( P.get<int>("Time.Scheme"))
    {
        case 1 :
            TimeScheme = new InstatStokesThetaSchemeCL<StokesProblemT, StokesSolverBaseCL>( Stokes, *stokessolver, P.get<double>("Time.Theta"));
            break;
        case 2 :
            TimeScheme = new StokesFracStepSchemeCL<InstatStokesThetaSchemeCL, StokesProblemT, StokesSolverBaseCL> ( Stokes, *stokessolver);
            break;
        default : throw DROPSErrCL("Unknown TimeDiscMethod");
    }

    StokesVelBndDataCL::bnd_val_fun ZeroVel = InVecMap::getInstance().find("ZeroVel")->second;

    if (P.get<int>("Time.NumSteps") == 0) {
        if ( P.get<int>("NavStokesSolver.OseenSolver.Solver") < 500000) {
            factory.SetMatrixA( &Stokes.A.Data.GetFinest());
            //for Stokes-MGM: coarse level solver uses bbt
            factory.SetMatrices( &Stokes.A.Data, &Stokes.B.Data,
                                 &Stokes.M.Data, &Stokes.prM.Data, &Stokes.pr_idx);
        }
        SolveStatProblem<MyStokesCL>( Stokes, *stokessolver);

    }
    else {
        TimeScheme->SetTimeStep( P.get<double>("Time.FinalTime") / P.get<double>("Time.NumSteps") );
        if ( P.get<std::string>("NavStokes.Coeff.Solution_Vel").compare("None")!=0)
            Stokes.InitVel( &Stokes.v, StokesFlowCoeffCL::LsgVel);
        else
            Stokes.InitVel( &Stokes.v, ZeroVel);
        if (P.get<int>("NavStokesSolver.OseenSolver.Solver") < 500000) {
            factory.SetMatrixA ( &TimeScheme->GetUpperLeftBlock()->GetFinest());
            factory.SetMatrices( TimeScheme->GetUpperLeftBlock(), &Stokes.B.Data,
                                 &Stokes.M.Data, &Stokes.prM.Data, &Stokes.pr_idx);
        }
    }

    // Output-Registrations:
    Ensight6OutCL* ensight = nullptr;

    if (P.get<int>("Ensight.EnsightOut",0)){
        // Initialize Ensight6 output
        const std::string ensf = P.get<string>("Ensight.EnsDir") + "/" + P.get<string>("Ensight.EnsCase");
        ensight = new Ensight6OutCL (P.get<string>("Ensight.EnsCase")+".case",
                                     P.get<int>("Time.NumSteps")/P.get("Ensight.EnsightOut", 0)+1,
                                     P.get<int>("Ensight.Binary"));

        ensight->Register( make_Ensight6Geom  ( MG, MG.GetLastLevel(), P.get<string>("Ensight.GeomName"),
                                                ensf + ".geo", /*time_dependent*/ false));
        ensight->Register( make_Ensight6Scalar( Stokes.GetPrSolution(),  "Pressure", ensf + ".pr",  true));
        ensight->Register( make_Ensight6Vector( Stokes.GetVelSolution(), "Velocity", ensf + ".vel", true));

        ensight->Write();
    }

    //create output for divergence of velocity
    VecDescCL divu( Stokes.p.RowIdx );
    ComputeProjectedDivergence( divu, Stokes );

    VTKOutCL * vtkwriter = nullptr;
    if (P.get<int>("VTK.VTKOut",0)){
        vtkwriter = new VTKOutCL(MG, "DROPS data",
                                 P.get<int>("Time.NumSteps")/P.get("VTK.VTKOut", 0)+1,
                                 P.get<std::string>("VTK.VTKDir"), P.get<std::string>("VTK.VTKName"),
                                 P.get<std::string>("VTK.TimeFileName"),
                                 P.get<int>("VTK.Binary"),
                                 P.get<int>("VTK.UseOnlyP1"),
                                 false, /*P2DG*/
                                 -1,  /* <- level */
                                 P.get<int>("VTK.ReUseTimeFile"),
                                 P.get<int>("VTK.UseDeformation"));
        vtkwriter->Register( make_VTKVector( Stokes.GetVelSolution(), "velocity") );
        vtkwriter->Register( make_VTKScalar( Stokes.GetPrSolution(), "pressure") );
        vtkwriter->Register( make_VTKScalar( Stokes.GetPrSolution(divu), "divu"));
        vtkwriter->Write(Stokes.v.t);
    }



    for ( int step = 1; step <= P.get<int>("Time.NumSteps"); ++step) {
        timer.Reset();

        std::cout << line << "Step: " << step << std::endl;

        TimeScheme->DoStep( Stokes.v.Data, Stokes.p.Data);
        ComputeProjectedDivergence( divu, Stokes );

        timer.Stop();
        std::cout << " o Solved system with:\n"
                  << "   - time          " << timer.GetTime()    << " s\n";

        // check the result (<= what is meant here?)
        if (ensight && step%P.get("Ensight.EnsightOut", 0)==0)
            ensight->Write( Stokes.v.t);
        if (vtkwriter && step%P.get("VTK.VTKOut", 0)==0)
            vtkwriter->Write(Stokes.v.t);

    }

    if( P.get<string>("NavStokes.Coeff.Solution_Vel").compare("None")!=0 && P.get<int>("Time.NumSteps") != 0)  // check whether solution is given
        Stokes.CheckSolution( &Stokes.v, &Stokes.p, StokesFlowCoeffCL::LsgVel, StokesFlowCoeffCL::DLsgVel, StokesFlowCoeffCL::LsgPr, false);

    if (vtkwriter) delete vtkwriter;
    if (ensight) delete ensight;
    delete stokessolver;
}

} // end of namespace DROPS

int main ( int argc, char** argv)
{
    try
    {

        DROPS::read_parameter_file_from_cmdline( P, argc, argv, "../../param/stokes/sdropsP2/stokes.json");
        P.put_if_unset<std::string>("VTK.TimeFileName",P.get<std::string>("VTK.VTKName"));
        std::cout << P << std::endl;

        DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"), P.get<std::vector<std::string> >("General.DynamicLibs") );
        if (P.get<int>("General.ProgressBar"))
            DROPS::ProgressBarTetraAccumulatorCL::Activate();

        // Check MarkLower value
        if( P.get<std::string>("Mesh.Type").compare("ReadMeshBuilder") == 0 )
            P.put("Misc.MarkLower", 0);
        else {
            double dy = norm(P.get<DROPS::Point3DCL>("Mesh.E2"));
            if (P.get("Misc.MarkLower", 0)<0 || P.get("Misc.MarkLower", 0) > dy)
            {
              std::cerr << "Wrong value of MarkLower\n";
              return 1;
            }
        }

        // time measurement
        DROPS::TimerCL timer;

        // set up data structure to represent a Stokes problem
        // ---------------------------------------------------------------------
        std::cout << line << "Set up data structure to represent a Stokes problem ...\n";
        timer.Reset();

        //create geometry
        DROPS::MultiGridCL* mg= 0;

        std::unique_ptr<DROPS::MGBuilderCL> builder( DROPS::make_MGBuilder( P));
        mg = new DROPS::MultiGridCL( *builder );

        DROPS::BndDataCL<DROPS::Point3DCL> velbnddata(0);
        DROPS::BndDataCL<double>           prbnddata(0);
        DROPS::read_BndData( velbnddata, *mg, P.get_child( "NavStokes.BoundaryData.Velocity"));
        std::cout << "Generated boundary conditions for velocity ";
        DROPS::read_BndData( prbnddata,  *mg, P.get_child( "NavStokes.BoundaryData.Pressure"));
        std::cout << "and pressure." << std::endl;

        // Setup the problem
        DROPS::StokesFlowCoeffCL tmp = DROPS::StokesFlowCoeffCL( P);
        StokesOnBrickCL prob( *mg, tmp, DROPS::StokesBndDataCL(velbnddata,prbnddata) );
        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s" << std::endl;

        // Refine the grid
        // ---------------------------------------------------------------------
        std::cout << "Refine the grid " << P.get<int>("Mesh.InitialRef") << " times regulary ...\n";
        timer.Reset();

        // Create new tetrahedra
        for ( int ref=1; ref<=P.get<int>("Mesh.InitialRef"); ++ref){
            std::cout << " refine (" << ref << ")\n";
            DROPS::MarkAll( *mg);
            mg->Refine();
        }

        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s" << std::endl;
        mg->SizeInfo(std::cout);

        // Solve the problem
        DROPS::Strategy( prob); //do all the stuff
        std::cout << DROPS::SanityMGOutCL(*mg) << std::endl;

        double min= prob.p.Data.min(),
               max= prob.p.Data.max();
        std::cout << "pressure min/max: "<<min<<", "<<max<<std::endl;

        // maple/geomview-output
//      DROPS::RBColorMapperCL colormap;
//      std::ofstream maple("maple.txt");
//      DROPS::Point3DCL e3(0.0); e3[2]= 1.0;
//      maple << DROPS::MapleMGOutCL(*mg, -1, false, true, DROPS::PlaneCL(e3, 0.6)) << std::endl;
//      std::ofstream fil("geom.off");
//      fil << DROPS::GeomSolOutCL<DROPS::PoissonP1CL<PoissonCoeffCL<DROPS::Params> >::DiscSolCL>( *mg, prob.GetSolution(), &colormap, -1, false, 0.0, prob.x.Data.min(), prob.x.Data.max()) << std::endl;
//      std::cout << DROPS::GeomMGOutCL(*mg, -1, true) << std::endl;
        delete mg;
        return 0;
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
}
