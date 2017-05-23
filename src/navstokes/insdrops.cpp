/// \file
/// \brief Solve the non-stationary Navier-Stokes-equations an a uniform grid.
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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
#include "stokes/stokes.h"
#include "num/nssolver.h"
#include "navstokes/navstokes.h"
#include "navstokes/integrTime.h"
#include <fstream>


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

namespace DROPS // for Strategy
{

using ::MyPdeCL;
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
void UzawaSolverCL<PoissonSolverT>::Solve(const MLMatrixCL& A, const MLMatrixCL& B, const MLMatrixCL &, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c, const DummyExchangeCL&, const DummyExchangeCL &)
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
          _ssor( omega), _PCGsolver( _ssor, inner_iter, inner_tol)
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

template<class Coeff>
void Strategy(NavierStokesP2P1CL<Coeff>& NS, int num_ref, double fp_tol, int fp_maxiter,
              double deco_red, int stokes_maxiter, double poi_tol, int poi_maxiter,
              double theta, double dt)
// flow control
{
    typedef NavierStokesP2P1CL<Coeff> NavStokesCL;

    MultiGridCL& MG= NS.GetMG();

    MLIdxDescCL  loc_vidx, loc_pidx;
    MLIdxDescCL* vidx1= &NS.vel_idx;
    MLIdxDescCL* pidx1= &NS.pr_idx;
    MLIdxDescCL* vidx2= &loc_vidx;
    MLIdxDescCL* pidx2= &loc_pidx;

    VecDescCL     loc_p;
    VelVecDescCL  loc_v;
    VelVecDescCL* v1= &NS.v;
    VelVecDescCL* v2= &loc_v;
    VecDescCL*    p1= &NS.p;
    VecDescCL*    p2= &loc_p;
    VelVecDescCL* b= &NS.b;
    VelVecDescCL* cplN= &NS.cplN;
    VelVecDescCL* cplM= &NS.cplM;
    VelVecDescCL* c= &NS.c;

    MLMatDescCL* A= &NS.A;
    MLMatDescCL* B= &NS.B;
    MLMatDescCL* N= &NS.N;
    MLMatDescCL* M= &NS.M;
    int refstep= 0;

    vidx1->SetFE( vecP2_FE);
    vidx2->SetFE( vecP2_FE);
    pidx1->SetFE( P1_FE);
    pidx2->SetFE( P1_FE);

    TimerCL time;
    do
    {
        MG.Refine();
        NS.CreateNumberingVel( MG.GetLastLevel(), vidx1);
        NS.CreateNumberingPr ( MG.GetLastLevel(), pidx1);
        std::cout << "altes und neues TriangLevel: " << (refstep!=0 ?
        int(vidx2->TriangLevel()) : -1) << ", " << vidx1->TriangLevel() << std::endl;
        MG.SizeInfo(std::cout);
        b->SetIdx(vidx1);
        c->SetIdx(pidx1);
        p1->SetIdx(pidx1);
        v1->SetIdx(vidx1);
        cplN->SetIdx(vidx1);
        cplM->SetIdx(vidx1);
        std::cout << "#Druck-Unbekannte: " << p2->Data.size() << ", "
                  << p1->Data.size() << std::endl;
        std::cout << "#Geschwindigkeitsunbekannte: " << v2->Data.size() << ", "
                  << v1->Data.size() << std::endl;
        if (v2->RowIdx)
        {
            const StokesBndDataCL& BndData= NS.GetBndData();
            P1EvalCL<double, const StokesPrBndDataCL, const VecDescCL>  pr2(p2, &BndData.Pr, &MG);
            P1EvalCL<double, const StokesPrBndDataCL, VecDescCL>        pr1(p1, &BndData.Pr, &MG);
            Interpolate(pr1, pr2);
            v2->Reset();
            p2->Reset();
        }
        A->SetIdx(vidx1, vidx1);
        B->SetIdx(pidx1, vidx1);
        N->SetIdx(vidx1, vidx1);
        NS.InitVel(v1, &MyPdeCL::LsgVel);
        time.Reset();
        time.Start();
    	NS.SetupSystem1( A, M, b, b, b, 0.0);
        NS.SetupSystem2( B, c, 0.0);
        NS.b.Clear(0.0);
        NS.c.Clear(0.0);

        time.Stop();
        std::cout << "SetupInstatSystem: " << time.GetTime() << " seconds" << std::endl;
        time.Reset();
        time.Start();
        A->Data * v1->Data;
        time.Stop();
        std::cout << " A*x: " << time.GetTime() << " seconds" << std::endl;
        time.Reset();

        MLMatDescCL M_pr;
        M_pr.SetIdx( pidx1, pidx1);
        NS.SetupPrMass( &M_pr);
//        AFPDeCo_Uzawa_PCG_CL<NavStokesCL> statsolver(NS, M_pr.Data, fp_maxiter, fp_tol,
//                                                     stokes_maxiter, poi_maxiter, poi_tol, deco_red);
        FPDeCo_Uzawa_PCG_CL<NavStokesCL> statsolver(NS, M_pr.Data.GetFinest(), fp_maxiter, fp_tol,
                                                    stokes_maxiter, poi_maxiter, poi_tol, deco_red);
//        FPDeCo_Uzawa_CG_CL<NavStokesCL> statsolver(NS, M_pr.Data, fp_maxiter, fp_tol,
//                                                   stokes_maxiter, poi_maxiter, poi_tol, deco_red);
//        FPDeCo_Uzawa_SGSPCG_CL<NavStokesCL> statsolver(NS, M_pr.Data, fp_maxiter, fp_tol,
//                                                   stokes_maxiter, poi_maxiter, poi_tol, deco_red);
//        AFPDeCo_Schur_PCG_CL<NavStokesCL> statsolver(NS, fp_maxiter, fp_tol,
//                                                   stokes_maxiter, poi_maxiter, poi_tol, deco_red);
//        FPDeCo_Schur_PCG_CL<NavStokesCL> statsolver(NS, fp_maxiter, fp_tol,
//                                                    stokes_maxiter, poi_maxiter, poi_tol, deco_red);
        // If the saddlepoint-problem is solved via an Uzawa-method, the mass-matrix alone is
        // not an appropriate preconditioner for the Schur-Complement-Matrix. M has to be scaled
        // by 1/(theta*dt).
        static_cast<Uzawa_PCG_CL&>(statsolver.GetStokesSolver()).SetTau(theta*dt); // Betrachte den Code in num/stokessolver.h: M ist bei zeitabhaqengigen Problemen kein geeigneter Vorkonditionierer.
        time.Reset();
        time.Start();
        NS.SetupNonlinear(N, v1, cplN, 0.);
        time.Stop();
        std::cout << "SetupNonlinear: " << time.GetTime() << " seconds" << std::endl;
        NS.SetupInstatRhs( b, c, cplM, 0., b, 0.);
//        InstatNavStokesThetaSchemeCL<NavStokesCL, FPDeCo_Schur_PCG_CL<NavStokesCL> >
//          instatsolver(NS, statsolver, theta);
        InstatNavStokesThetaSchemeCL<NavStokesCL, FPDeCo_Uzawa_PCG_CL<NavStokesCL> >
            instatsolver(NS, statsolver, theta);
        std::cout << "After constructor." << std::endl;
        double t= 0.;
        NS.v.t= 0;
        for (int timestep=0; t<1.; ++timestep, t+= dt, NS.v.t+= dt)
        {
//if (timestep==25) return;
            std::cout << "------------------------------------------------------------------"
                      << std::endl << "t: " << t << std::endl;
            NS.SetTime(t+dt); // We have to set the new time!
            if (timestep==0) // Check the initial solution, at least velocities.
                NS.CheckSolution(v1, vidx1, p1, &MyPdeCL::LsgVel, &MyPdeCL::LsgPr);
            instatsolver.SetTimeStep(dt);
            std::cout << "Before timestep." << std::endl;
            instatsolver.DoStep(*v1, p1->Data);
            std::cout << "After timestep." << std::endl;
            NS.CheckSolution(v1, vidx1, p1, &MyPdeCL::LsgVel, &MyPdeCL::LsgPr);
        }
        MarkAll(MG);

        A->Reset();
        B->Reset();
        M->Reset();
        N->Reset();
//      M_pr.Reset();
        b->Reset();
        c->Reset();
        cplN->Reset();
        cplM->Reset();
        std::swap(v2, v1);
        std::swap(p2, p1);
        std::swap(vidx2, vidx1);
        std::swap(pidx2, pidx1);
        std::cout << std::endl;
    }
    while (refstep++ < num_ref);
    // we want the solution to be in _v
    if (v2 == &loc_v)
    {
        NS.vel_idx.swap( loc_vidx);
        NS.pr_idx.swap( loc_pidx);
        NS.v.SetIdx(&NS.vel_idx);
        NS.p.SetIdx(&NS.pr_idx);

        NS.v.Data= loc_v.Data;
        NS.p.Data= loc_p.Data;
    }
}

} // end of namespace DROPS


int main (int argc, char** argv)
{
  try
  {
    if (argc!=10)
    {
        std::cout << "Usage (insdrops): <num_refinement> <fp_tol> <fp_maxiter> "
                  << "<deco_red> <stokes_maxiter> <poi_tol> <poi_maxiter> "
                  << "<theta> <dt>" << std::endl;
        return 1;
    }
    DROPS::BrickBuilderCL brick(DROPS::std_basis<3>(0),
                                DROPS::std_basis<3>(1),
                                DROPS::std_basis<3>(2),
                                DROPS::std_basis<3>(3),
                                4,4,4);
    const bool IsNeumann[6]= {false, false, false, false, false, false};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel,
          &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel };
    DROPS::RBColorMapperCL colormap;

    int num_ref= std::atoi(argv[1]);
    double fp_tol= std::atof(argv[2]);
    int fp_maxiter= std::atoi(argv[3]);
    double deco_red= std::atof(argv[4]);
    int stokes_maxiter= std::atoi(argv[5]);
    double poi_tol= std::atof(argv[6]);
    int poi_maxiter= std::atoi(argv[7]);
    double theta= std::atof(argv[8]);
    double dt= std::atof(argv[9]);
    std::cout << "num_ref: " << num_ref << ", ";
    std::cout << "fp_tol: " << fp_tol<< ", ";
    std::cout << "fp_maxiter: " << fp_maxiter << ", ";
    std::cout << "deco_red: " << deco_red << ", ";
    std::cout << "stokes_maxiter: " << stokes_maxiter << ", ";
    std::cout << "poi_tol: " << poi_tol << ", ";
    std::cout << "poi_maxiter: " << poi_maxiter << ", ";
    std::cout << "theta: " << theta << ", ";
    std::cout << "dt: " << dt << std::endl;

    typedef DROPS::NavierStokesP2P1CL<MyPdeCL::StokesCoeffCL>
            NSOnBrickCL;
    typedef NSOnBrickCL MyNavierStokesCL;
    MyNavierStokesCL prob(brick, MyPdeCL::StokesCoeffCL(),
                          DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
    DROPS::MultiGridCL& mg = prob.GetMG();

    Strategy(prob, num_ref, fp_tol, fp_maxiter, deco_red, stokes_maxiter, poi_tol, poi_maxiter,
             theta, dt);

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
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
