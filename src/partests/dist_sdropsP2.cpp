//**************************************************************************
// File:    drops_statStokes.cpp                                           *
// Content: Solver for Stokes problem with Taylor-Hood FE                  *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file drops_statStokes.cpp
/// \brief Solver for Stokes problem with Taylor-Hood FE
/** This main problem solves
    \f{eqnarray*}{
      -\Delta u + \nabla p &=& f \text{ in } \Omega:=[0,1]^3 \\
      -\operatorname{div} p &=& 0 \text{ in } \Omega:=[0,1]^3
    \f}
    for a given solution \f$ (u,p) \f$.
*/

 // include parallel computing!
#ifdef _PAR
#include "parallel/parallel.h"          // proc handling, reduce operations, ...
#include "parallel/partime.h"           // parallel time-messurement
#include "parallel/parmultigrid.h"      // handle multigrid over different procs
#include "parallel/loadbal.h"           // distribute multigrid
#include "num/krylovsolver.h"           // various parallel solvers
#include "num/precond.h"                // various parallel preconditioners
#include "num/oseensolver.h"            // various parallel stokes solvers
#endif

 // include geometric computing
#include "geom/multigrid.h"             // multigrid on each processor
#include "geom/builder.h"               // construuct the initial multigrid

 // include numeric computing!
#include "num/fe.h"

 // include problem class
#include "stokes/stokes.h"      // setting up the Stokes problem
#include "num/bndData.h"

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

using namespace std;

const char line[] ="----------------------------------------------------------------------------------\n";

namespace DROPS
{

/// \brief Coefficients of the Stokes problem
struct StatStokesCL
{
    /// \brief Velocity solution
    static Point3DCL SolVel(const Point3DCL& p, double)
    {
        SVectorCL<3> ret;
        ret[0]=    std::sin(p[0])*std::sin(p[1])*std::sin(p[2]);
        ret[1]=  - std::cos(p[0])*std::cos(p[1])*std::sin(p[2]);
        ret[2]= 2.*std::cos(p[0])*std::sin(p[1])*std::cos(p[2]);
        return ret/3.;
    }

    /// \brief Pressure solution
    static double SolPr(const Point3DCL& p, double)
    {
        return std::cos(p[0])*std::sin(p[1])*std::sin(p[2]);
    }

    /// \brief Coefficients of the Stokes problem
    /// du/dt + q*u - nu*laplace u + Dp = f
    ///                          -div u = 0
    class StokesCoeffCL
    {
      public:
        static double q(const Point3DCL&) { return 0.0; }
        /// \brief right hand side
        static SVectorCL<3> f(const Point3DCL& p, double)
            { SVectorCL<3> ret(0.0); ret[2]= 3*std::cos(p[0])*std::sin(p[1])*std::cos(p[2]); return ret; }
        /// \brief diffusion
        const double nu;
        /// \brief Constructor
        StokesCoeffCL() : nu(1.0) {}
    };
    static StokesCoeffCL Coeff;
};
StatStokesCL::StokesCoeffCL StatStokesCL::Coeff;

/// \brief Solve an Oseen problem, i.e. discretized Stokes problem
/** Discretizing a Stokes problem results in a linear equation system of the
    form:
    \f[
      \left[ \begin{array}{cc} A&B^T\\B&0\end{array}\right]
      \cdot \left( \begin{array}{c} u\\p \end{array} \right)
      =     \left( \begin{array}{c} b\\c \end{array} \right)
    \f]
    Although this is a symmetric matrix, we are using GMRes-type solvers instead
    of CG-type solvers, because these are "more" representativ of solving
    two-phase flow problems.
    \param A coefficient matrix
    \param B coefficient matrix
    \param u velocity solution
    \param p pressure solution
    \param b rhs
    \param c rhs
    \param vel_idx description of velocity DoF (only in parallel version)
    \param pr_idx description of pressure DoF (only in parallel version)
*/
template <typename Mat, typename Vec>
void Solve(const Mat &A, const Mat &B, const Mat &C, Vec &u, Vec &p, const Vec &b, const Vec &c,
           __UNUSED__ const IdxDescCL& vel_idx, __UNUSED__ const IdxDescCL& pr_idx,
           __UNUSED__ Mat& prM)
{
    // parameter for solver:
    __UNUSED__ const int    PCrestart = 100;
    const int    PCmaxiter = 1000;
    const double PCtol     = 0.01;
    const bool   PCrelative= true;
    const int    OuterIter = 200;
    const double OuterTol  = 1e-12;
    const double InnerRed  = 0.01;
    const int    InnerIter = 1000;
    __UNUSED__ std::ostream *output   = &std::cout;    // set pointer to 0 to make solver quiet (only in parallel version)

    // time measurement
    TimerCL timer;

    // type of preconditioner and solver
    typedef JACPcCL                               PCT;
    typedef GMResSolverCL<PCT>                    GMResSolverT;
    typedef SolverAsPreCL<GMResSolverT>           APcT;         // preconditioner of A
    typedef DummyPcCL                             SPcT;         // preconditioner of Schur-complement
    typedef InexactUzawaCL<APcT, SPcT, APC_OTHER> OseenSolverT; // solver of above system

    PCT          pc;
    GMResSolverT gmres( pc, PCrestart, PCmaxiter, PCtol, PCrelative);
    APcT         Apc (gmres);
    SPcT         Spc;
    OseenSolverT oseen( Apc, Spc, OuterIter, OuterTol, InnerRed, InnerIter);


    std::cout << " o Solving system with Inexact-Uzawa-Method: ... \n";

    // Solve the linear equation system
    timer.Reset();
    oseen.Solve( A, B, C, u, p, b, c, vel_idx.GetEx(), pr_idx.GetEx());
    timer.Stop();

    double realresid1, realresid2;
    Vec resid1( A*u + transp_mul(B, p) - b);
    Vec resid2( B*u - c);
#ifndef _PAR
    realresid1= norm( resid1);
    realresid2= norm( resid2);
#else
    realresid1= vel_idx.GetEx().Norm( resid1, false);
    realresid2=  pr_idx.GetEx().Norm( resid2, false);
#endif

    std::cout << " o Solved system with:\n"
              << "   - time                   " << timer.GetTime() << " s\n"
              << "   - iterations             " << oseen.GetIter() << '\n'
              << "   - real residuum velocity " << realresid1      << '\n'
              << "   - real residuum pressure " << realresid2      << std::endl;
}

/// \brief Strategy to solve the Stokes problem on a given triangulation
template <typename CoeffCL>
void Strategy( StokesP2P1CL<CoeffCL>& Stokes)
{
    // time measurement
#ifndef _PAR
    TimerCL timer;
#else
    ParTimerCL timer;
#endif

    // the triangulation
    MultiGridCL& mg= Stokes.GetMG();

    // connection triangulation and vectors
    // -------------------------------------------------------------------------
    std::cout << line << "Connecting triangulation and matrices/vectors ...\n";
    timer.Reset();

    // Set FE types:
    Stokes.vel_idx.SetFE( vecP2_FE);
    Stokes.pr_idx.SetFE(  P1_FE);
    // Number velocity and pressure DoF
    Stokes.CreateNumberingVel( mg.GetLastLevel(), &Stokes.vel_idx);
    Stokes.CreateNumberingPr ( mg.GetLastLevel(), &Stokes.pr_idx);
    // Tell matrices and vectors about DoF-numbering
    Stokes.v.SetIdx( &Stokes.vel_idx);
    Stokes.p.SetIdx( &Stokes.pr_idx);
    Stokes.b.SetIdx( &Stokes.vel_idx);
    Stokes.c.SetIdx( &Stokes.pr_idx);
    Stokes.A.SetIdx( &Stokes.vel_idx, &Stokes.vel_idx);
    Stokes.B.SetIdx( &Stokes.pr_idx, &Stokes.vel_idx);
    Stokes.C.SetIdx( &Stokes.pr_idx, &Stokes.pr_idx);

    timer.Stop();
    std::cout << " o time " << timer.GetTime() << " s" << std::endl;


    // display problem size
    // -------------------------------------------------------------------------
    std::cout << line << "Problem size\n";
#ifdef _PAR
    std::vector<size_t> VelUnkOnProc= ProcCL::Gather( Stokes.v.Data.size(), 0);
    std::vector<size_t> PrUnkOnProc = ProcCL::Gather( Stokes.p.Data.size(), 0);
    const IdxT VelNumUnk   = Stokes.vel_idx.GetGlobalNumUnknowns(),
               VelNumAccUnk= std::accumulate( VelUnkOnProc.begin(), VelUnkOnProc.end(), 0),
               PrNumUnk    = Stokes.pr_idx.GetGlobalNumUnknowns(),
               PrNumAccUnk = std::accumulate( PrUnkOnProc.begin(), PrUnkOnProc.end(), 0);
#else
    std::vector<size_t> VelUnkOnProc( 1);
    std::vector<size_t> PrUnkOnProc( 1);
    VelUnkOnProc[0]  = Stokes.v.Data.size();
    PrUnkOnProc[0]   = Stokes.p.Data.size();
    IdxT VelNumUnk   = Stokes.v.Data.size(),
         VelNumAccUnk= Stokes.v.Data.size(),
         PrNumUnk    = Stokes.p.Data.size(),
         PrNumAccUnk = Stokes.p.Data.size();
#endif
    std::cout << " o number of velocity unknowns, accumulated "
              << VelNumAccUnk << ", global " << VelNumUnk << '\n'
              << " o number of pressure unknowns, accumulated "
              << PrNumAccUnk << ", global " << PrNumUnk << '\n'
              << " o number of unknowns on proc: velocity, pressure\n";
    for (size_t i=0; i<VelUnkOnProc.size(); ++i){
        std::cout << " - Proc " << i << ": "
                  << VelUnkOnProc[i] << ", " << PrUnkOnProc[i] << '\n';
    }


    // discretize (setup linear equation system)
    // -------------------------------------------------------------------------
    std::cout << line << "Discretize (setup linear equation system) ...\n";

    timer.Reset();
    VelVecDescCL  cplM( &Stokes.vel_idx);
    Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &cplM, 0.0);
    Stokes.SetupSystem2( &Stokes.B, &Stokes.c, 0.0);
    MLMatDescCL prM( &Stokes.pr_idx, &Stokes.pr_idx);
    Stokes.SetupPrMass( &prM);
    timer.Stop();
    std::cout << " o time " << timer.GetTime() << " s" << std::endl;


    // solve the linear equation system
    // -------------------------------------------------------------------------
    std::cout << line << "Solve the linear equation system ...\n";
    Solve( Stokes.A.Data.GetFinest(), Stokes.B.Data.GetFinest(), Stokes.C.Data.GetFinest(),
           Stokes.v.Data, Stokes.p.Data,
           Stokes.b.Data, Stokes.c.Data,
           Stokes.vel_idx.GetFinest(), Stokes.pr_idx.GetFinest(), prM.Data.GetFinest());


    // check the result
    // -------------------------------------------------------------------------
/*    std::cout << line << "Check result against known solution ...\n";

    timer.Reset();
    Stokes.CheckSolution( &Stokes.v, &Stokes.p, &StatStokesCL::SolVel, &StatStokesCL::SolPr, 0);
    timer.Stop();
    std::cout << " o time " << timer.GetTime() << " s" << std::endl;*/
}

} // end of namespace DROPS

int main (int argc, char** argv)
{
#ifdef _PAR
    DROPS::ProcCL::Instance(&argc, &argv);
#endif
    try
    {
        if (argc!=5){
            std::cout << "Usage " << argv[0] << " <ref x> <ref y> <ref z> <refinement>\n"
                      << " with\n"
                      << " o ref x: spatial resolution in x direction\n"
                      << " o ref y: spatial resolution in y direction\n"
                      << " o ref z: spatial resolution in z direction\n"
                      << " o refinement: number of regular refinements of each tetrahedron\n\n"
                      << "Ending program" << std::endl;
            return 0;
        }

        // time measurement
#ifndef _PAR
        DROPS::TimerCL timer;
#else
        DROPS::ParTimerCL timer;
#endif

        // parameter for geometry
        const int refX  = atoi( argv[1]),
                  refY  = atoi( argv[2]),
                  refZ  = atoi( argv[3]);
        __UNUSED__ const int refAll= atoi( argv[4]);
        DROPS::Point3DCL orig, e1, e2, e3;  // origin and orientation of unit cube
        e1[0]= e2[1]= e3[2]= 1.0;

        // set up data structure to represent a poisson problem
        // ---------------------------------------------------------------------
        std::cout << line << "Set up data structure to represent a Poisson problem ...\n";
        timer.Reset();

        // create builder for geometry
        DROPS::BrickBuilderCL mgb( orig, e1, e2, e3, refX, refY, refZ);

        // boundary conditions
        DROPS::BndCondT bndcond[6] = { DROPS::DirBC, DROPS::DirBC, DROPS::DirBC,
                                       DROPS::DirBC, DROPS::DirBC, DROPS::DirBC };
        // boundary function
        DROPS::BndDataCL<DROPS::Point3DCL>::bnd_val_fun bndfunc[6] = {
            &DROPS::StatStokesCL::SolVel, &DROPS::StatStokesCL::SolVel, &DROPS::StatStokesCL::SolVel,
            &DROPS::StatStokesCL::SolVel, &DROPS::StatStokesCL::SolVel, &DROPS::StatStokesCL::SolVel };

        // boundary data ( = condition & function)
        DROPS::StokesBndDataCL bdata( 6, bndcond, bndfunc);

        // Setup the problem
        DROPS::StokesP2P1CL<DROPS::StatStokesCL::StokesCoeffCL> prob(
            mgb, DROPS::StatStokesCL::Coeff, bdata);
        DROPS::MultiGridCL& mg= prob.GetMG();

#ifdef _PAR
        // Set parallel data structures
        DROPS::LoadBalCL lb( mg);   // loadbalancing
        lb.DoMigration();           // distribute initial grid
#endif

        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s" << std::endl;

        // Refine the grid
        // ---------------------------------------------------------------------
/*
        std::cout << "Refine the grid " << refAll << " times regulary ...\n";
        timer.Reset();

        // Create new tetrahedra
        for ( int ref=1; ref<=refAll; ++ref){
            std::cout << " refine (" << ref << ")\n";
            DROPS::MarkAll( mg);
            mg.Refine();
        }

        // do loadbalancing
#ifdef _PAR
        lb.DoMigration();
#endif
        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s\n"
                  << " o distribution of elements" << '\n';
*/
        mg.SizeInfo(cout);

        // Solve the problem
        // ---------------------------------------------------------------------
        DROPS::Strategy( prob);

        return 0;
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
}

