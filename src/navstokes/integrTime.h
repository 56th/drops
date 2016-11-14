/// \file
/// \brief classes that perform time-integration steps
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen:

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

/// History: begin - Nov, 19 2001

#ifndef DROPS_NS_INTEGRTIME_H
#define DROPS_NS_INTEGRTIME_H

#include "navstokes/navstokes.h"

// TODO: FracStepScheme fuer instat. NavierStokes

namespace DROPS
{

template <class NavStokesT, class SolverT>
class InstatNavStokesThetaSchemeCL
/*****************************************************************************
*   for solving the instat. Navier-Stokes equation of type NavStokesT with a
*   1-step-theta-scheme: theta=1   -> impl. Euler (BDF1, backward Euler)
*                        theta=1/2 -> Crank-Nicholson (Trapezregel)
*
*   Inner stat. Navier-Stokes-type problems are solved with a SolverT-solver.
*   The matrices A, B, M, N and the rhs b, c, cplN of the NavStokes class have
*   to be set properly! After construction, SetTimeStep has to be called once.
*   Then every DoStep performs one step in time. Changing time steps require
*   further calls to SetTimeStep.
******************************************************************************/
{
  private:
//    typedef typename NavStokesT::VelVecDescCL VelVecDescCL;

    NavStokesT& NS_;
    SolverT&    solver_;

    VelVecDescCL *b_, *old_b_;        // rhs + couplings with poisson matrix A
    VelVecDescCL *cplM_, *old_cplM_;  // couplings with mass matrix M
    VelVecDescCL *cplN_;              // couplings with nonlinearity N
    VectorCL      rhs_;
    MLMatrixCL    L_;                 // 1./dt*M + theta*A  = linear part

    double theta_, dt_;

  public:
    InstatNavStokesThetaSchemeCL(NavStokesT& NS, SolverT& solver,
                                 double theta= 0.5, double t= 0.0)
        : NS_( NS), solver_( solver), b_( &NS.b), old_b_( new VelVecDescCL),
          cplM_( &NS.cplM), old_cplM_( new VelVecDescCL),
          cplN_( &NS.cplN), rhs_( NS.b.RowIdx->NumUnknowns()), theta_( theta)
    {
        old_b_->SetIdx( b_->RowIdx);
        old_cplM_->SetIdx( b_->RowIdx);
        // Redundant for NS_.c but does not change its value
        NS_.SetupInstatRhs( old_b_, &NS_.c, old_cplM_, t, old_b_, t);
    }

    ~InstatNavStokesThetaSchemeCL()
    {
        if (old_b_ == &NS_.b)
            delete b_;
        else
            delete old_b_;
        if (old_cplM_ == &NS_.cplM)
            delete cplM_;
        else
            delete old_cplM_;
    }

    double GetTheta()    const { return theta_; }
    double GetTimeStep() const { return dt_; }

    void SetTimeStep( double dt)
    {
        dt_= dt;
        L_.LinComb( 1./dt_, NS_.M.Data, theta_, NS_.A.Data);
    }

    void DoStep( VecDescCL& v, VectorCL& p);
};


//=================================
//     template definitions
//=================================

template <class NavStokesT, class SolverT>
void InstatNavStokesThetaSchemeCL<NavStokesT,SolverT>::DoStep( VecDescCL& v, VectorCL& p)
{
    // NS_.t contains the new time!
    NS_.SetupInstatRhs( b_, &NS_.c, cplM_, NS_.v.t, b_, NS_.v.t);
    const double alpha= theta_;
    const double beta= (theta_ - 1.);
    rhs_=  alpha*b_->Data
         - beta*(old_b_->Data + cplN_->Data)
         + beta*( NS_.A.Data*v.Data + NS_.N.Data*v.Data )
         +(1./dt_)*(cplM_->Data - old_cplM_->Data + NS_.M.Data*v.Data);
    solver_.Solve( L_, NS_.B.Data, NS_.C.Data, v, p, rhs_, *cplN_, NS_.c.Data, NS_.vel_idx.GetEx(), NS_.pr_idx.GetEx(), alpha);
    std::swap( old_b_, b_);
    std::swap( old_cplM_, cplM_);
}


} // end of namespace DROPS

#endif
