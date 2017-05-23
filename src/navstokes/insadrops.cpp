/// \file
/// \brief Solve the non-stationary Navier-Stokes-equations on a adaptive grid.
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt

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
#include "out/ensightOut.h"
#include "geom/builder.h"
#include "stokes/stokes.h"
#include "num/nssolver.h"
#include "navstokes/navstokes.h"
#include "navstokes/integrTime.h"
#include "num/fe_repair.h"
#include <fstream>
#include <sstream>


struct InstatNSCL
{
    static DROPS::SVectorCL<3> LsgVel(const DROPS::Point3DCL& p, double t)
    {
        DROPS::SVectorCL<3> ret;
        ret[0]= (2.*t - 1.)*p[0];
        ret[1]= (2.*t - 1.)*p[1];
        ret[2]= (2. - 4.*t)*p[2];
        return ret;
    }

    // int_{x=0..1, y=0..1,z=0..1} p(x,y,z,t) dx dy dz = 0 for all t.
    static double LsgPr(const DROPS::Point3DCL& p, double t)
    {
        return (0.5 - t)*p.norm_sq() + t - 0.5;
    }

    // Jacobi-matrix of exact solution (only in the spatial variables)
    static inline DROPS::SMatrixCL<3, 3> DxLsgVel(const DROPS::Point3DCL&, double t)
    {
        DROPS::SMatrixCL<3, 3> ret(0.0);
        ret(0,0)= 2.*t - 1.;
        ret(1,1)= 2.*t - 1.;
        ret(2,2)= 2. - 4.*t;
        return ret;
    }

    // Time-derivative of exact solution
    static inline DROPS::SVectorCL<3> DtLsgVel(const DROPS::Point3DCL& p, double)
    {
        DROPS::SVectorCL<3> ret(0.0);
        ret[0]= 2.*p[0];
        ret[1]= 2.*p[1];
        ret[2]= -4.*p[2];
        return ret;
    }
    // u_t + q*u - nu*laplace u + (u*D)u + Dp = f
    //                                 -div u = 0
    class StokesCoeffCL
    {
      public:
        static double q(const DROPS::Point3DCL&, double) { return 0.0; }
        static DROPS::SVectorCL<3> f(const DROPS::Point3DCL& p, double t)
        {
            DROPS::SVectorCL<3> ret;
            ret[0]= (4.*t*t - 6.*t + 4.)*p[0];
            ret[1]= (4.*t*t - 6.*t + 4.)*p[1];
            ret[2]= (16.*t*t -18.*t + 1.)*p[2];
            return ret;
        }
        const double nu;

        StokesCoeffCL() : nu(1.0) {}
    };

    static StokesCoeffCL Coeff;
};


InstatNSCL::StokesCoeffCL InstatNSCL::Coeff;

typedef InstatNSCL MyPdeCL;

typedef DROPS::SVectorCL<3> (*fun_ptr)(const DROPS::SVectorCL<3>&, double);



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
               const double width,  // Thickness of refined shell on each side of the interface
               const int c_level,   // Outside the shell, use this level
               const int f_level,   // Inside the shell, use this level
               const double t)      // Time of evaluation
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
            const int l= it->GetLevel();
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

template<class Coeff>
void
UpdateTriangulation(DROPS::NavierStokesP2P1CL<Coeff>& NS,
                    const signed_dist_fun Dist,
                    const double t,
                    const double width,   // Thickness of refined shell on each side of the interface
                    const int c_level,    // Outside the shell, use this level
                    const int f_level,    // Inside the shell, use this level
                    DROPS::VelVecDescCL* v1,
                    DROPS::VecDescCL* p1)
{
    using namespace DROPS;
    Assert( 0<=c_level && c_level<=f_level, "UpdateTriangulation: Levels are cheesy.\n", ~0);
    TimerCL time;

    time.Reset();
    time.Start();
    MultiGridCL& mg= NS.GetMG();
    IdxDescCL  loc_vidx, loc_pidx;
    IdxDescCL* vidx1= v1->RowIdx;
    IdxDescCL* pidx1= p1->RowIdx;
    VelVecDescCL  loc_v;
    VecDescCL     loc_p;

    loc_vidx.SetFE( vecP2_FE);
    loc_pidx.SetFE( P1_FE);

    bool shell_not_ready= true;
    const Uint min_ref_num= f_level - c_level;
    const StokesBndDataCL& BndData= NS.GetBndData();
    Uint i;
    for(i=0; shell_not_ready || i<min_ref_num; ++i) {
        RepairP2CL<Point3DCL>::type p2repair(mg, *v1, BndData.Vel);
        RepairP1CL<double>::type p1repair(mg, *p1, BndData.Pr);

        shell_not_ready= ModifyGridStep( mg, Dist, width, c_level, f_level, t);

        // Repair velocity
        loc_vidx.CreateNumbering( mg.GetLastLevel(), mg, NS.GetBndData().Vel);
        if ( mg.GetLastLevel() != vidx1->TriangLevel()) {
            std::cout << "LastLevel: " << mg.GetLastLevel()
                      << " loc_vidx->TriangLevel: " << loc_vidx.TriangLevel() << std::endl;
            throw DROPSErrCL( "Strategy: Sorry, not yet implemented.");
        }
        loc_v.SetIdx( &loc_vidx);
        p2repair.repair( loc_v);
        vidx1->DeleteNumbering( mg);
        vidx1->swap(loc_vidx);
        v1->Data.swap(loc_v.Data);
        loc_v.Clear( t);

        // Repair pressure
        loc_pidx.CreateNumbering( mg.GetLastLevel(), mg, NS.GetBndData().Pr);
        loc_p.SetIdx( &loc_pidx);
        p1repair.repair( loc_p);
        pidx1->DeleteNumbering( mg);
        pidx1->swap(loc_pidx);
        p1->Data.swap(loc_p.Data);
        loc_p.Clear( t);
    }

    time.Stop();
    std::cout << "UpdateTriangulation: " << i
              << " refinements in " << time.GetTime() << " seconds\n"
              << "last level: " << mg.GetLastLevel() << '\n';
    mg.SizeInfo( std::cout);
}

void
MakeInitialTriangulation(DROPS::MultiGridCL& mg,
                         const signed_dist_fun Dist,
                         const double width,  // Thickness of refined shell on eache side of the interface
                         const int c_level,   // Outside the shell, use this level
                         const int f_level)   // Inside the shell, use this level
{
    using namespace DROPS;
    Assert( 0<=c_level && c_level<=f_level, "MakeInitialTriangulation: Levels are cheesy.\n", ~0);
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

template<class Coeff>
void
SetMatVecIndices(DROPS::NavierStokesP2P1CL<Coeff>& NS,
                 DROPS::MLIdxDescCL* const vidx,
                 DROPS::MLIdxDescCL* const pidx)
{
    std::cout << "#Druck-Unbekannte: " << pidx->NumUnknowns() << std::endl;
    std::cout << "#Geschwindigkeitsunbekannte: " << vidx->NumUnknowns() << std::endl;
    NS.b.SetIdx( vidx);
    NS.c.SetIdx( pidx);
    NS.cplM.SetIdx( vidx);
    NS.cplN.SetIdx( vidx);
    NS.A.SetIdx( vidx, vidx);
    NS.B.SetIdx( pidx, vidx);
    NS.M.SetIdx( vidx, vidx);
    NS.N.SetIdx( vidx, vidx);
}

template<class Coeff>
void
ResetSystem(DROPS::NavierStokesP2P1CL<Coeff>& NS)
{
    NS.A.Reset(); NS.B.Reset();
    NS.M.Reset(); NS.N.Reset();
    NS.b.Reset(); NS.c.Reset();
    NS.cplM.Reset(); NS.cplN.Reset();
}

namespace DROPS{

typedef PCGSolverCL<SSORPcCL>     PCG_SsorCL;

template <typename PoissonSolverT>
class UzawaSolverCL : public StokesSolverBaseCL
{
  private:
    PoissonSolverT& _poissonSolver;
    MatrixCL&       _M;
    double          _tau;

    template <typename Mat, typename Vec>
    void doSolve( const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c);

  public:
    UzawaSolverCL (PoissonSolverT& solver, MatrixCL& M, int maxiter, double tol, double tau= 1.)
        : StokesSolverBaseCL(maxiter,tol), _poissonSolver(solver), _M(M), _tau(tau) {}
    UzawaSolverCL (PoissonSolverT& solver, MLMatrixCL& M, int maxiter, double tol, double tau= 1.)
    : StokesSolverBaseCL(maxiter,tol), _poissonSolver(solver), _M(M.GetFinest()), _tau(tau) {}


    double GetTau()            const { return _tau; }
    void   SetTau( double tau)       { _tau= tau; }

    void Solve( const MatrixCL& A, const MatrixCL& B, const MatrixCL&, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c, const DummyExchangeCL&, const DummyExchangeCL&);
    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, const MLMatrixCL&, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c, const DummyExchangeCL&, const DummyExchangeCL&);
};

template <class PoissonSolverT>
template <typename Mat, typename Vec>
void UzawaSolverCL<PoissonSolverT>::doSolve(
    const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c)
{
    Vec v_corr(v.size()),
        p_corr(p.size()),
        res1(v.size()),
        res2(p.size());

    double tol= tol_;
    tol*= tol;
    Uint output= 50;//max_iter/20;  // nur 20 Ausgaben pro Lauf

    double res1_norm= 0., res2_norm= 0.;
    for( iter_=0; iter_<maxiter_; ++iter_) {
        z_xpay(res2, B*v, -1.0, c);
        res2_norm= norm_sq( res2);
        _poissonSolver.SetTol( std::sqrt( res2_norm)/20.0);
        _poissonSolver.Solve(_M, p_corr, res2, DummyExchangeCL());
//        p+= _tau * p_corr;
        axpy(_tau, p_corr, p);
//        res1= A*v + transp_mul(B,p) - b;
        z_xpaypby2(res1, A*v, 1.0, transp_mul(B,p), -1.0, b);
        res1_norm= norm_sq( res1);
        if (res1_norm + res2_norm < tol) {
            res_= std::sqrt( res1_norm + res2_norm );
            return;
        }
        if( (iter_%output)==0)
            std::cout << "step " << iter_ << ": norm of 1st eq= " << std::sqrt( res1_norm)
                      << ", norm of 2nd eq= " << std::sqrt( res2_norm) << std::endl;

        _poissonSolver.SetTol( std::sqrt( res1_norm)/20.0);
        _poissonSolver.Solve( A, v_corr, res1, DummyExchangeCL());
        v-= v_corr;
    }
    res_= std::sqrt( res1_norm + res2_norm );
}

template <class PoissonSolverT>
void UzawaSolverCL<PoissonSolverT>::Solve(
    const MatrixCL& A, const MatrixCL& B, const MatrixCL&, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c, const DummyExchangeCL&, const DummyExchangeCL&)
{
    doSolve( A, B, v, p, b, c);
}

template <class PoissonSolverT>
void UzawaSolverCL<PoissonSolverT>::Solve(
    const MLMatrixCL& A, const MLMatrixCL& B, const MLMatrixCL&, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c, const DummyExchangeCL&, const DummyExchangeCL&)
{
    doSolve( A, B, v, p, b, c);
}

class Uzawa_PCG_CL : public UzawaSolverCL<PCG_SsorCL>
{
  private:
    SSORPcCL   _ssor;
    PCG_SsorCL _PCGsolver;
  public:
    Uzawa_PCG_CL( MatrixCL& M, int outer_iter, double outer_tol, int inner_iter, double inner_tol, double tau= 1., double omega=1.)
        : UzawaSolverCL<PCG_SsorCL>( _PCGsolver, M, outer_iter, outer_tol, tau),
          _ssor( omega), _PCGsolver(_ssor, inner_iter, inner_tol)
        {}
};

template <class NavStokesT>
class FPDeCo_Uzawa_PCG_CL: public AdaptFixedPtDefectCorrCL<NavStokesT>
{
  private:
    Uzawa_PCG_CL _uzawaSolver;

  public:
    FPDeCo_Uzawa_PCG_CL( NavStokesT& NS, MatrixCL& M, int fp_maxiter, double fp_tol, int stokes_maxiter,
                         int poiss_maxiter, double poiss_tol, double reduction= 0.1)
        : AdaptFixedPtDefectCorrCL<NavStokesT>( NS, _uzawaSolver, fp_maxiter, fp_tol, reduction, false),
          _uzawaSolver( M, stokes_maxiter, fp_tol, poiss_maxiter, poiss_tol) // outer_tol will be set by the AFPDeCo-solver!
        {}
};
}//end of namespace DROPS

template<class Coeff>
void
Strategy(DROPS::NavierStokesP2P1CL<Coeff>& NS,
         double fp_tol, int fp_maxiter,
         double deco_red, int stokes_maxiter,
         double poi_tol, int poi_maxiter,
         double theta,
         DROPS::Uint num_timestep,
         double shell_width, int c_level, int f_level)
{
    using namespace DROPS;
    typedef NavierStokesP2P1CL<Coeff> NavStokesCL;

    MultiGridCL& MG= NS.GetMG();
    MLIdxDescCL* vidx1= &NS.vel_idx;
    MLIdxDescCL* pidx1= &NS.pr_idx;
    VelVecDescCL* v1= &NS.v;
    VecDescCL*    p1= &NS.p;
    MLMatDescCL  M_pr;
    vidx1->SetFE( vecP2_FE);
    pidx1->SetFE( P1_FE);
    TimerCL time;
    double t= 0.;
    const double dt= 1./num_timestep;
    NS.v.t= 0;
    Uint timestep= 0;

    FPDeCo_Uzawa_PCG_CL<NavStokesCL>* statsolver= 0;
//    FPDeCo_Uzawa_CG_CL<NavStokesCL>* statsolver= 0;
//    FPDeCo_Uzawa_SGSPCG_CL<NavStokesCL>* statsolver= 0;
//    AFPDeCo_Schur_PCG_CL<NavStokesCL>* statsolver= 0;
//    FPDeCo_Schur_PCG_CL<NavStokesCL>* statsolver= 0;
//    InstatNavStokesThetaSchemeCL<NavStokesCL, FPDeCo_Schur_PCG_CL<NavStokesCL> >* instatsolver= 0;
    InstatNavStokesThetaSchemeCL<NavStokesCL, FPDeCo_Uzawa_PCG_CL<NavStokesCL> >* instatsolver= 0;

    MakeInitialTriangulation( MG, &SignedDistToInterface, shell_width, c_level, f_level);
    NS.CreateNumberingVel( MG.GetLastLevel(), vidx1);
    v1->SetIdx( vidx1);
    NS.InitVel( v1, &MyPdeCL::LsgVel);
    NS.CreateNumberingPr( MG.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);

    // Initialize Ensight6 output
    std::string ensf( "insa");
    Ensight6OutCL ensight( ensf + ".case", num_timestep + 1);
    ensight.Register( make_Ensight6Geom      ( MG, MG.GetLastLevel(),   "insa geometry", ensf + ".geo", true));
    ensight.Register( make_Ensight6Scalar    ( NS.GetPrSolution(),  "p",      ensf + ".pr",  true));
    ensight.Register( make_Ensight6Vector    ( NS.GetVelSolution(), "v",      ensf + ".vel", true));
    ensight.Write( 0.);

    for (; timestep<num_timestep; ++timestep, t+= dt, NS.v.t+= dt) {
        std::cout << "----------------------------------------------------------------------------"
                  << std::endl << "t: " << t << std::endl;
        if (timestep%(num_timestep/10) == 0) { // modify the grid
            if (timestep>0) { // perform cleanup, which is not neccessary for t==0.
                delete statsolver; statsolver= 0;
                delete instatsolver; instatsolver= 0;
                M_pr.Reset();
                ResetSystem( NS);
                UpdateTriangulation( NS, &SignedDistToInterface, t,
                                     shell_width, c_level, f_level, v1, p1);
            }
            SetMatVecIndices( NS, vidx1, pidx1);
            time.Reset(); time.Start();
        	NS.SetupSystem1( &NS.A, &NS.M, &NS.b, &NS.b, &NS.b, 0.0);
            NS.SetupSystem2( &NS.B, &NS.c, 0.0);
            NS.b.Clear(0.0);
            NS.c.Clear(0.0);
            time.Stop();
            std::cout << "SetupInstatSystem: " << time.GetTime() << " seconds" << std::endl;
            time.Reset();
            M_pr.SetIdx( pidx1, pidx1);
            NS.SetupPrMass( &M_pr);
//            AFPDeCo_Uzawa_PCG_CL<NavStokesCL> statsolver( NS, M_pr.Data, fp_maxiter, fp_tol,
//                                                          stokes_maxiter, poi_maxiter, poi_tol, deco_red);
            statsolver= new FPDeCo_Uzawa_PCG_CL<NavStokesCL> ( NS, M_pr.Data.GetFinest(), fp_maxiter, fp_tol,
                           stokes_maxiter, poi_maxiter, poi_tol, deco_red);
//            statsolver= new FPDeCo_Uzawa_CG_CL<NavStokesCL>( NS, M_pr.Data, fp_maxiter, fp_tol,
//                           stokes_maxiter, poi_maxiter, poi_tol, deco_red);
//            statsolver= new FPDeCo_Uzawa_SGSPCG_CL<NavStokesCL>( NS, M_pr.Data, fp_maxiter, fp_tol,
//                           stokes_maxiter, poi_maxiter, poi_tol, deco_red);
//            statsolver= new AFPDeCo_Schur_PCG_CL<NavStokesCL>( NS, fp_maxiter, fp_tol,
//                           stokes_maxiter, poi_maxiter, poi_tol, deco_red);
//            statsolver= new FPDeCo_Schur_PCG_CL<NavStokesCL>( NS, fp_maxiter, fp_tol,
//                           stokes_maxiter, poi_maxiter, poi_tol, deco_red);
            // If the saddlepoint-problem is solved via an Uzawa-method, the mass-matrix alone is
            // not an appropriate preconditioner for the Schur-Complement-Matrix. M has to be scaled
            // by 1/(theta*dt).
            static_cast<Uzawa_PCG_CL&>(statsolver->GetStokesSolver()).SetTau( theta*dt); // Betrachte den Code in num/stokessolver.h: M ist bei zeitabhaqengigen Problemen kein geeigneter Vorkonditionierer.
            time.Reset(); time.Start();
            NS.SetupNonlinear( &NS.N, v1, &NS.cplN, t);
            time.Stop();
            std::cout << "SetupNonlinear: " << time.GetTime() << " seconds" << std::endl;
            NS.SetupInstatRhs( &NS.b, &NS.c, &NS.cplM, t, &NS.b, t);
//            instatsolver= new InstatNavStokesThetaSchemeCL<NavStokesCL,
//                              FPDeCo_Schur_PCG_CL<NavStokesCL> >(NS, *statsolver, theta);
            instatsolver= new InstatNavStokesThetaSchemeCL<NavStokesCL,
                             FPDeCo_Uzawa_PCG_CL<NavStokesCL> >( NS, *statsolver, theta, t);
            if (timestep == 0) // check initial velocities
                NS.CheckSolution( v1, vidx1, p1, &MyPdeCL::LsgVel, &MyPdeCL::LsgPr);
        }
        NS.SetTime( t+dt); // We have to set the new time!
        instatsolver->SetTimeStep( dt);
        std::cout << "Before timestep." << std::endl;
        instatsolver->DoStep( *v1, p1->Data);
        std::cout << "After timestep." << std::endl;
        NS.CheckSolution( v1, vidx1, p1, &MyPdeCL::LsgVel, &MyPdeCL::LsgPr);
        ensight.Write( t+dt);
    }
    delete statsolver; statsolver= 0;
    delete instatsolver; instatsolver= 0;
    ResetSystem( NS);
    M_pr.Reset();
}


int main (int argc, char** argv)
{
  try
  {
    if (argc!=12)
    {
        std::cout << "Usage (insadrops):  <fp_tol> <fp_maxiter> "
                  << "<deco_red> <stokes_maxiter> <poi_tol> <poi_maxiter> "
                  << "<theta> <num_timestep> <shell_width> <c_level> <f_level>" << std::endl;
        return 1;
    }
    DROPS::BrickBuilderCL brick(DROPS::std_basis<3>(0),
                                DROPS::std_basis<3>(1),
                                DROPS::std_basis<3>(2),
                                DROPS::std_basis<3>(3),
                                2,2,2);
    const bool IsNeumann[6]= {false, false, false, false, false, false};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel,
          &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel };
    DROPS::RBColorMapperCL colormap;

    double fp_tol= std::atof(argv[1]);
    int fp_maxiter= std::atoi(argv[2]);
    double deco_red= std::atof(argv[3]);
    int stokes_maxiter= std::atoi(argv[4]);
    double poi_tol= std::atof(argv[5]);
    int poi_maxiter= std::atoi(argv[6]);
    double theta= std::atof(argv[7]);
    int num_timestep= std::atoi(argv[8]);
    double shell_width= std::atof(argv[9]);
    int c_level= std::atoi(argv[10]);
    int f_level= std::atoi(argv[11]);

    std::cout << "fp_tol: " << fp_tol<< ", ";
    std::cout << "fp_maxiter: " << fp_maxiter << ", ";
    std::cout << "deco_red: " << deco_red << ", ";
    std::cout << "stokes_maxiter: " << stokes_maxiter << ", ";
    std::cout << "poi_tol: " << poi_tol << ", ";
    std::cout << "poi_maxiter: " << poi_maxiter << ", ";
    std::cout << "theta: " << theta << ", ";
    std::cout << "num_timestep: " << num_timestep <<  ", ";
    std::cout << "shell_width: " << shell_width <<  ", ";
    std::cout << "c_level: " << c_level << ", ";
    std::cout << "f_level: " << f_level << std::endl;

    typedef DROPS::NavierStokesP2P1CL<MyPdeCL::StokesCoeffCL> NSOnBrickCL;
    typedef NSOnBrickCL MyNavierStokesCL;
    MyNavierStokesCL prob(brick, MyPdeCL::StokesCoeffCL(),
                          DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
    DROPS::MultiGridCL& mg = prob.GetMG();

    Strategy(prob, fp_tol, fp_maxiter, deco_red, stokes_maxiter, poi_tol, poi_maxiter,
             theta, num_timestep, shell_width, c_level, f_level);

    std::cout << "hallo" << std::endl;
    std::cout << DROPS::SanityMGOutCL(mg) << std::endl;
    std::ofstream fil("navstokespr.off");
    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    fil << DROPS::GeomSolOutCL<MyNavierStokesCL::const_DiscPrSolCL>(mg, prob.GetPrSolution(), &colormap, -1, false, 0.0, min, max) << std::endl;
    fil.close();

    DROPS::MLIdxDescCL tecIdx( DROPS::P1_FE);
    prob.CreateNumberingPr( mg.GetLastLevel(), &tecIdx);
    std::ofstream v2d("navstokestec2D.dat");
    DROPS::TecPlot2DSolOutCL< MyNavierStokesCL::const_DiscVelSolCL, MyNavierStokesCL::const_DiscPrSolCL>
        tecplot2d( mg, prob.GetVelSolution(), prob.GetPrSolution(), tecIdx.GetFinest(), -1, 2, 0.5); // cutplane is z=0.5
    v2d << tecplot2d;
    v2d.close();
    v2d.open("navstokestec2D2.dat");
    DROPS::TecPlot2DSolOutCL< MyNavierStokesCL::const_DiscVelSolCL, MyNavierStokesCL::const_DiscPrSolCL>
        tecplot2d2( mg, prob.GetVelSolution(), prob.GetPrSolution(), tecIdx.GetFinest(), -1, 1, 0.5); // cutplane is z=0.5
    v2d << tecplot2d2;
    v2d.close();
    return 0;
  }
  catch (DROPS::DROPSErrCL& err) { err.handle(); }
}
