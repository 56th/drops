/// \file parsolver.h
/// \brief Parallel basic iterative solvers
/// \author LNM RWTH Aachen: Patrick Esser, Sven Gross, Joerg Peters, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef _DROPS_PARSOLVER_H_
#define _DROPS_PARSOLVER_H_

#include "num/spmat.h"
#include "num/parlanczos.h"
#include "num/solver.h"
#include "parallel/parallel.h"
#include "misc/problem.h"
#include "parallel/exchange.h"
#include <algorithm>

// for file i/o
//#define FILE_OUT
#ifdef FILE_OUT
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#endif

// *****************************************************************************
// *                         T O C                                             *
// *****************************************************************************
// * Definition of parallel solver functions                                    *
// * Definition of parallel solver base classes                                 *
// *   - ParSolverBaseCL                                                        *
// *   - ParPreSolverBaseCL                                                     *
// * Definition of parallel solver classes                                      *
// *   - Conjugate Gradients (ParCGSolverCL)                                    *
// *   - preconditioned Conjugate Gradients (ParPCGSolverCL)                    *
// *   - preconditioned Generalized Minimal Residual  (ParPreGMResSolverCL)     *
// *   - preconditioned Bi-Conjugate Gradient Stabilized (ParBiCGSTABSolverCL)  *
// *   - preconditioned Generalized Conjugate Residuals (ParPreGCRSolverCL)     *
// *   - (preconditioned) Quasi Minimal Residual (ParQMRSolverCL)               *
// *   - (preconditioned) CG for the normal equations (ParCGNeSolverCL)         *
// * Declaration of parallel solver functions                                   *
// ******************************************************************************

namespace DROPS
{

// const double BreakDownC=1e-35;

//***************************************************************************
// implemented parallel methods for solving a linear system
//***************************************************************************
template <typename Mat, typename Vec, typename ExCL>
bool ParCG(const Mat& A,Vec& x_acc,const Vec& b, const ExCL& ExX, int& max_iter, double& tol, 
           bool measure_relative_tol, std::ostream* output=0);

template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParPCG(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M, int& max_iter, 
            double& tol, bool measure_relative_tol=false, std::ostream* output=0);

// Preconditioned GMRES
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParGMRES(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M,
              int m, int& max_iter, double& tol, bool measure_relative_tol=true,
              PreMethGMRES method=LeftPreconditioning,
              std::ostream* output=0);

// Preconditioned GMRES with modifications for better scalability.
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParModGMRES(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M,
                 int m, int& max_iter, double& tol, bool measure_relative_tol=true,
                 bool useMGS=false, PreMethGMRES method=LeftPreconditioning,
                 std::ostream* output=0);

// BiCGSTAB
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParBiCGSTAB(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M, int& max_iter, double& tol, bool measure_relative_tol=true);

// QMR-Solver
template <typename Mat, typename Vec, typename Lanczos, typename ExCL>
bool ParQMR(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, Lanczos lan,
            int& max_iter, double& tol, bool measure_relative_tol);

// Preconditioned GCR
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParPGCR(const Mat& A, Vec& x, const Vec& b, const ExCL& ExX, PreCon& M,
    int m, int& max_iter, double& tol, bool measure_relative_tol= true, std::ostream* output=0);

template <typename Mat, typename Vec, typename PreCon, typename ExACL, typename ExATranspCL>
bool ParPCGNE(const Mat& A, Vec& u, const Vec& b, const ExACL& ExAX, const ExATranspCL& ExATranspX, const PreCon& M,
    int& max_iter, double& tol, bool measure_relative_tol= false, std::ostream* output=0);


//***************************************************************************
//                      S O L V E R   B A S E  C L A S S E S
//***************************************************************************

// ***************************************************************************
/// \brief Parallel base solver class
// ***************************************************************************
class ParSolverBaseCL : public SolverBaseCL
/** All parallel solvers needs additionally to SolverBaseCL a class for
    exchanging numerical values. Since these information is stored as a part
    of the IdxDescCL, the parallel solvers stores a reference to the IdxDescCL
    as well.
*/
{
  private:
    const IdxDescCL* idx_;      // for getting the ExchangeCL

  public:
    /// \brief Constructor
    ParSolverBaseCL(int maxiter, double tol, const IdxDescCL& idx, bool rel= false, std::ostream* output=0)
      : SolverBaseCL(maxiter, tol, rel, output), idx_(&idx) {}
    /// \brief Constructor, that does not initialize the index description
    ParSolverBaseCL(int maxiter, double tol, bool rel= false, std::ostream* output=0)
      : SolverBaseCL(maxiter, tol, rel, output), idx_(0) {}

    /// \brief Ask for ExchangeCL
    const ExchangeCL& GetEx() const {
        Assert(idx_, DROPSErrCL("ParSolverBaseCL::GetEx: Index not set, do you want to use an ExchangeBlockCL?"), DebugParallelNumC);
        return idx_->GetEx();
    }
};

// ***************************************************************************
/// \brief Parallel preconditioned base solver class
// ***************************************************************************
template <typename PC>
class ParPreSolverBaseCL : public ParSolverBaseCL
/** See ParSolverBaseCL with an preconditioner and a flag if the residual should be measured relative.*/
{
  private:
    typedef ParSolverBaseCL base;

  protected:
    PC   &pc_;                                      // preconditioner

  public:
    typedef PC PrecondT;                            ///< Preconditioner

    /// \brief Constructor
    ParPreSolverBaseCL(int maxiter, double tol, const IdxDescCL& idx, PC& pc, bool rel=true, std::ostream* output=0)
      : base(maxiter, tol, idx, rel, output), pc_(pc) {}
    /// \brief Constructor, that does not initialize the index description
    ParPreSolverBaseCL(int maxiter, double tol, PC& pc, bool rel=true, std::ostream* output=0)
      : base(maxiter, tol, rel, output), pc_(pc) {}

    PC& GetPC()             {return pc_;}          ///< return reference on preconditioner
    const PC& GetPC() const {return pc_;}          ///< return constant reference on preconditioner

    void SetPC(PC& pc) { pc_= pc; }                ///< Set preconditioner
};

//***************************************************************************
//                      S O L V E R  C L A S S E S
//***************************************************************************

// ***************************************************************************
/// \brief Parallel CG-Solver class
// ***************************************************************************
class ParCGSolverCL : public ParSolverBaseCL
{
  private:
    typedef ParSolverBaseCL base;

  public:
    /// \brief Constructor for parallel CG-Solver
    /** Tries to solve a linear equation system within \a maxiter steps with
        accuracy \a tol. The ExCL \a ex is used to do parallel inner products.*/
    ParCGSolverCL(int maxiter, double tol, const IdxDescCL& idx, bool rel= false, std::ostream* output=0)
      : base(maxiter, tol, idx, rel, output) {}

    /// \brief Solve a linear equation system with Conjugate-Gradients-Method
    template <typename Mat, typename Vec>
      void Solve(const Mat& A, Vec& x, const Vec& b)
    /// Solve the linear equation system with coefficient matrix \a A and rhs \a b iterative with
    /// CG algorithm, uses \a x as start-vector and result vector.
    /// \post x has accumulated form
    {
        base::_res=  base::_tol;
        base::_iter= base::_maxiter;
        ParCG(A, x, b, base::GetEx(), base::_iter, base::_res, base::rel_, base::output_);
    }
};


// ***************************************************************************
/// \brief Parallel preconditioned CG-Solver class
// ***************************************************************************
template <typename PC>
class ParPCGSolverCL : public ParPreSolverBaseCL<PC>
{
  private:
    typedef ParPreSolverBaseCL<PC> base;

  public:
    /// \brief Constructor for the parallel preconditioned CG Solver
    /** Tries to solve a linear equation system within \a maxiter steps with
        accuracy \a tol. The ExCL \a ex is used to do parallel inner products. \a pc is
        the given preconditioner. If \a rel is given, the residual is computed relative and
        with \a acc the inner products are determined with accure variant (see ExchangeCL). */
    ParPCGSolverCL(int maxiter, double tol, const IdxDescCL &idx, PC& pc, bool rel=false, std::ostream* output=0)
      : base(maxiter, tol, idx, pc, rel, output) {}
    ParPCGSolverCL(int maxiter, double tol, PC& pc, bool rel=false, std::ostream* output=0)
      : base(maxiter, tol, pc, rel, output) {}

    /// \brief Solve a linear equation system with Conjugate Gradients-Method
    template <typename Mat, typename Vec>
      void Solve(const Mat& A, Vec& x, const Vec& b)
    /// Solve the linear equation system with coefficient matrix \a A and rhs \a b iterative with
    /// preconditioned CG algorithm, uses \a x as start-vector and result vector.
    /// \post x has accumulated form
    {
        base::_res=  base::_tol;
        base::_iter= base::_maxiter;
        ParPCG(A, x, b, base::GetEx(), base::GetPC(), base::_iter, base::_res, base::rel_, base::output_);
    }
    template <typename Mat, typename Vec>
      void Solve(const Mat& A, Vec& x, const Vec &b, const ExchangeBlockCL& ex)
    /// Solve the linear equation system with coefficient matrix \a A and rhs \a b iterative with
    /// preconditioned GCR algorithm, uses \a x as start-vector and result vector.
    /// \post x has accumulated form
    {
        base::_res  = base::_tol;
        base::_iter = base::_maxiter;
        ParPCG(A, x, b, ex, base::GetPC(), base::_iter, base::_res, base::rel_, base::output_);
    }

};

// ***************************************************************************
/// \brief Parallel preconditioned GMRES-Solver class
// ***************************************************************************
template <typename PC>
class ParPreGMResSolverCL : public ParPreSolverBaseCL<PC>
{
  private:
    typedef ParPreSolverBaseCL<PC> base;
    int          restart_;                  // number of iterations before restart
    bool         useModGS_;                 // which Gram-Schmidt method should be used to compute Krylov basis
    PreMethGMRES method_;                   // left or right preconditioning
    bool         mod_;                      // use modified variant for better scalability


  public:
    /// \brief Constructor of the parallel preconditioned GMRES-Solver
    /** Tries to solve a linear equation system within \a maxiter steps with
        accuracy \a tol. After \a restart steps, a restart is performed. The ExCL
        \a ex is used to do parallel inner products. \a pc is
        the given preconditioner. If \a rel is given, the residual is computed relative.
        (this configuration needs less memory!). By setting \a ModGS the modified Gram-Schmidt
        algorithm is used for the Arnoldi method.*/
    ParPreGMResSolverCL(int restart, int maxiter, double tol, const IdxDescCL& idx, PC &pc,
                        bool rel=true, bool ModGS=false,
                        PreMethGMRES method=LeftPreconditioning, bool mod=true,
                        std::ostream* output=0)
      : base(maxiter, tol, idx, pc, rel, output),
        restart_(restart), useModGS_(ModGS), method_(method), mod_(mod) {}
    ParPreGMResSolverCL(int restart, int maxiter, double tol, PC &pc,
                        bool rel=true, bool ModGS=false,
                        PreMethGMRES method=LeftPreconditioning, bool mod=true,
                        std::ostream* output=0)
      : base(maxiter, tol, pc, rel, output),
        restart_(restart), useModGS_(ModGS), method_(method), mod_(mod) {}

    int  GetRestart()           const { return restart_; }  ///< number of iterations before restart
    void SetRestart(int restart)      { restart_=restart; } ///< set number of iterations before restart

    /// \brief Solve a linear equation system with a preconditioned Generalized Minimal Residuals-Method
    template <typename Mat, typename Vec>
      void Solve(const Mat& A, Vec& x, const Vec& b)
    /// Solve the linear equation system with coefficient matrix \a A and rhs \a b iterative with
    /// preconditioned GMRES algorithm, uses \a x as start-vector and result vector.
    /// \post x has accumulated form
    {
        base::_res=  base::_tol;
        base::_iter= base::_maxiter;

        if (mod_)
            ParModGMRES(A, x, b, base::GetEx(), base::GetPC(), restart_,
                        base::_iter, base::_res, base::GetRelError(),
                        useModGS_, method_, base::output_);
        else
            ParGMRES(A, x, b, base::GetEx(),  base::GetPC(), restart_,
                     base::_iter, base::_res, base::GetRelError(),
                     method_, base::output_);
    }
    /// \brief Solve a linear equation system with a preconditioned Generalized Minimal Residuals-Method
    template <typename Mat, typename Vec>
      void Solve(const Mat& A, Vec& x, const Vec& b, const ExchangeBlockCL& ex)
    /// Solve the linear equation system with coefficient matrix \a A and rhs \a b iterative with
    /// preconditioned GMRES algorithm, uses \a x as start-vector and result vector.
    /// \post x has accumulated form
    {
        base::_res=  base::_tol;
        base::_iter= base::_maxiter;

        if (mod_)
            ParModGMRES(A, x, b, ex, base::GetPC(), restart_,
                        base::_iter, base::_res, base::GetRelError(),
                        useModGS_, method_, base::output_);
        else
            ParGMRES(A, x, b, ex, base::GetPC(), restart_,
                     base::_iter, base::_res, base::GetRelError(),
                     method_, base::output_);
    }
};

// ***************************************************************************
/// \brief Parallel BiCGSTAB-Solver class
// ***************************************************************************
template<typename PC>
class ParBiCGSTABSolverCL : public ParPreSolverBaseCL<PC>
{
  private:
    typedef ParPreSolverBaseCL<PC> base;

  public:
    /// \brief Constructor of the parallel preconditioned BiCGStab-Solver
    /** Tries to solve a linear equation system within \a maxiter steps with
        accuracy \a tol. The ExCL \a ex is used to do parallel inner products. \a pc is
        the given preconditioner. If \a rel is given, the residual is computed relative. */
    ParBiCGSTABSolverCL(int maxiter, double tol, const IdxDescCL &idx, PC& pc, bool rel=true,
                        std::ostream* output=0)
      : base(maxiter, tol, idx, pc, rel, output) {}

    /// \brief Solve a linear equation system with preconditioned Bi-Conjugate Gradient Stabilized-Method
    template <typename Mat, typename Vec>
      void Solve(const Mat& A, Vec& x, const Vec& b)
    /// Solve the linear equation system with coefficient matrix \a A and rhs \a b iterative with
    /// preconditioned BiCGStab algorithm, uses \a x as start-vector and result vector.
    /// \post x has accumulated form
    {
        base::_res=  base::_tol;
        base::_iter= base::_maxiter;
        ParBiCGSTAB(A, x, b, base::GetEx(), base::GetPC(), base::_iter, base::_res, base::GetRelError());
    }
};

// ***************************************************************************
/// \brief Parallel GCR-Solver class with preconditioning
// ***************************************************************************
template <typename PC>
class ParPreGCRSolverCL : public ParPreSolverBaseCL<PC>
{
  private:
    typedef ParPreSolverBaseCL<PC> base;
    int  trunc_;
    bool mod_;

  public:
    /// \brief Constructor of the parallel preconditioned GCR-Solver
    /** Tries to solve a linear equation system within \a maxiter steps with
        accuracy \a tol. \a trunc vectors are used to span the Krylov subspace. The ExCL
        \a ex is used to do parallel inner products. \a pc is
        the given preconditioner. If \a rel is given, the residual is computed relative and
        with \a acc the inner products are determined with accure variant (see ExchnageCL).
        If \a mod is set than a modified variant for computing the Krylov subspace
        is used, to reduce sync-points.
        \todo (of) <b>truncation strategy with modified GCR do not work!</b>
    */
    ParPreGCRSolverCL(int trunc, int maxiter, double tol, const IdxDescCL& idx, PC &pc, bool mod=false,
                      bool rel=true,std::ostream* output=0)
      : base(maxiter, tol, idx, pc, rel, output), trunc_(trunc), mod_(mod) {}
    /// \brief Constructor, that does not initialize the index description
    ParPreGCRSolverCL(int trunc, int maxiter, double tol, PC &pc, bool mod=false,
                      bool rel=true, std::ostream* output=0)
      : base(maxiter, tol, pc, rel, output), trunc_(trunc), mod_(mod) {}

    /// \brief Solve a linear equation system with preconditioned Generalized Conjugate Residuals-Method
    template <typename Mat, typename Vec>
      void Solve(const Mat& A, Vec& x, const Vec &b)
    /// Solve the linear equation system with coefficient matrix \a A and rhs \a b iterative with
    /// preconditioned GCR algorithm, uses \a x as start-vector and result vector.
    /// \post x has accumulated form
    {
        base::_res  = base::_tol;
        base::_iter = base::_maxiter;
        ParPGCR(A, x, b, base::GetEx(), base::GetPC(), trunc_, base::_iter, base::_res, base::GetRelError(), base::output_);
    }

    /// \brief Solve a linear equation system with preconditioned Generalized Conjugate Residuals-Method
    template <typename Mat, typename Vec>
      void Solve(const Mat& A, Vec& x, const Vec &b, const ExchangeBlockCL& ex)
    /// Solve the linear equation system with coefficient matrix \a A and rhs \a b iterative with
    /// preconditioned GCR algorithm, uses \a x as start-vector and result vector.
    /// \post x has accumulated form
    {
        base::_res  = base::_tol;
        base::_iter = base::_maxiter;
        ParPGCR(A, x, b, ex, base::GetPC(), trunc_, base::_iter, base::_res, base::GetRelError(), base::output_);
    }

};

// ***************************************************************************
/// \brief Parallel QMR-Solver class
// ***************************************************************************
template<typename Lanczos>
class ParQMRSolverCL : public ParSolverBaseCL
{
  private:
    Lanczos *lan_;
    typedef ParSolverBaseCL base;

  public:
    /// \brief Constructor of the parallel preconditioned QMR-Solver
    /** Tries to solve a linear equation system within \a maxiter steps with
        accuracy \a tol. The ExCL \a ex is used to do parallel inner products. The
        Lanczos-Algorithm to compute the bi-orthogonal-basis is given by a Lanczos class \a lan.
        If \a measure_relative_tol is given, the residual is computed relative.*/
    ParQMRSolverCL(int maxiter, double tol, const IdxDescCL& idx, Lanczos &lan, bool rel=true, std::ostream* output=0) :
        base(maxiter, tol, idx, rel, output), lan_(&lan) {}

    /// \brief Solve a linear equation system with Quasi Minimal Residuals-Method
    template <typename Mat, typename Vec>
      void Solve(const Mat& A, Vec& x, const Vec &b)
    /// Solve the linear equation system with coefficient matrix \a A and rhs \a b iterative with
    /// QMR algorithm, uses \a x as start-vector and result vector.
    /// \post x has accumulated form
    {
        base::_res  = base::_tol;
        base::_iter = base::_maxiter;
        ParQMR(A, x, b, base::GetEx(), *lan_, base::_iter, base::_res, base::GetRelError());
    }
};

///\brief Solver for A*A^Tx=b with Craig's method and left preconditioning.
///
/// A preconditioned CG version for matrices of the form A*A^T. Note that *A* must be
/// supplied, not A*A^T, to the Solve-method.
template <typename PC>
class ParPCGNESolverCL : public ParPreSolverBaseCL<PC>
{
  private:
    typedef ParPreSolverBaseCL<PC> base;
    const IdxDescCL* idxAtransp_;

  public:
    ParPCGNESolverCL( int maxiter, double tol, const IdxDescCL& idxA, const IdxDescCL& idxAtransp, PC& pc, bool rel= false, std::ostream* output=0)
        : base( maxiter, tol, idxA, pc, rel, output), idxAtransp_( &idxAtransp) {}

    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b)
    {
        base::_res=  base::_tol;
        base::_iter= base::_maxiter;
        ParPCGNE(A, x, b, base::GetEx(), idxAtransp_->GetEx(), base::GetPC(), base::_iter, base::_res, base::rel_, base::output_);
    }

    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b, int& numIter, double& resid) const
    {
        resid=   base::_tol;
        numIter= base::_maxiter;
        ParPCGNE(A, x, b, base::GetEx(), idxAtransp_->GetEx(), base::GetPC(), numIter, resid, base::rel_, base::output_);
    }
};

//***************************************************************************
// Implementations of the methods
//***************************************************************************

/// \brief Parallel CG-Algorithm
template <typename Mat, typename Vec, typename ExCL>
bool ParCG(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, int& max_iter, 
        double& tol, bool measure_relative_tol, std::ostream* output)
    /// \param[in]     A        local distributed coefficients-matrix of the linear equation system
    /// \param[in,out] x_acc    start vector and the solution in accumulated form
    /// \param[in]     b        rhs of the linear equation system (distributed form)
    /// \param[in]     ExX      ExchangeCL corresponding to the RowIdx of x and the ColIdx of A
    /// \param[in,out] max_iter IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol      IN: tolerance for the residual, OUT: residual
    /// \return                 convergence within max_iter iterations
{
    Vec r( b - A*x_acc ), r_acc(r);
    double normb= ExX.Norm( b, false), 
           res,   
           resid=ExX.Norm_sq( r, false, &r_acc);
    Vec d_acc( -r_acc);

    if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;

    if ((res= std::sqrt( resid)/normb) <= tol)
    {
        tol= res;
        max_iter= 0;
        return true;
    }

    for (int i= 1; i <= max_iter; ++i)
    {
        const Vec    Ad= A*d_acc;
        const double delta= ExX.ParDot( Ad, false, d_acc, true);
        const double alpha= resid/delta;
        double       beta= resid;

        axpy(alpha, d_acc, x_acc);  // x+= alpha*d;
        axpy(alpha, Ad, r); // r+= alpha*Ad;

        resid= ExX.Norm_sq( r, false, &r_acc);

        if ( output){
            (*output) << "ParCG: " << i << " resid " << std::sqrt(resid) << std::endl;
        }

        if ((res= std::sqrt( resid)/normb) <= tol)
        {
            tol= res;
            max_iter= i;
            return true;
        }
        beta= resid / beta;
        d_acc= beta*d_acc-r_acc;
    }
    tol= res;
    return false;
}


/// \brief Parallel preconditioned CG-Algorithm
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParPCG(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX,
            PreCon& M, int& max_iter, double& tol, bool measure_relative_tol,
            std::ostream* output)
    /// \param[in]     A                    local distributed coefficients-matrix of the linear equation system
    /// \param[in,out] x_acc                start vector and the solution in accumulated form
    /// \param[in]     b                    rhs of the linear equation system (distributed form)
    /// \param[in]     ExX                  ExchangeCL corresponding to the RowIdx of x and the ColIdx of A
    /// \param[in,out] M                    Preconditioner
    /// \param[in,out] max_iter             IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol                  IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol measure resid relative
    /// \return                             convergence within max_iter iterations
{
    // Check if preconditioner needs diagonal of matrix. The preconditioner
    // only computes the diagonal new, if the matrix has changed
    if (M.NeedDiag())
        M.SetDiag(A);

    const size_t n= b.size();
    Vec p_acc(n), z(n), z_acc(n), q(n), q_acc(n), r( b - A*x_acc), r_acc(r);

    double rho,
           rho_1= 0.0,
           resid= ExX.Norm( r, false, &r_acc),
           normb= ExX.Norm( b, false);

    if (normb == 0.0 || measure_relative_tol == false)
        normb= 1.0;
    resid = resid/normb;


    if (resid<=tol){
        tol= resid;
        max_iter= 0;
        return true;
    }

    for (int i=1; i<=max_iter; ++i)
    {
        M.Apply(A, z, r);
        rho = ExX.ParDot( z, M.RetAcc(), r_acc, true, &z_acc);

        if (i == 1)
            p_acc= z_acc;
        else
            p_acc = z_acc + (rho/rho_1)*p_acc;

        q= A*p_acc;

        const double lambda = ExX.ParDot( q, false, p_acc, true, &q_acc);
        const double alpha  = rho/lambda;

        x_acc += alpha * p_acc;
        r     -= alpha * q;
        r_acc -= alpha * q_acc;

        resid= ExX.Norm( r_acc, true) / normb;

        if ( output){
            (*output) << "ParPCG: " << i << " resid " << resid << std::endl;
        }

        if (resid<=tol){
            tol= resid;
            max_iter= i;
            return true;
        }
        rho_1= rho;
    }
    tol= resid;
    return false;
}


/// \brief Computes an orthogonal vector on i vectors by the standard Gram-Schmidt method
template <typename Vec, typename ExCL>
void StandardGramSchmidt(DMatrixCL<double>& H,
                          Vec& w, bool acc_w, const std::vector<Vec>& v, bool acc_v,
                          int i, const ExCL& ex, Vec& tmpHCol)
/// This routine uses only one synchronization point.
/// \param H       Hessenberg matrix
/// \param w       orthogonal vector
/// \param acc_w   is w accumulated
/// \param v       vectors used to orthogonalize w
/// \param acc_v   is v accumulated
/// \param i       number of vectors v
/// \param ex      class to accumulate a vector
/// \param tmpHCol temporary vector to store a column of H (as parameter so no mem is allocated in this routine)
{
    if (acc_v&&acc_w){                          // both vectors are accumulated
        for (int k=0; k<=i; ++k)
            tmpHCol[k] = ex.LocalDot( w, true, v[k], true);
    }
    else if (!acc_v&&!acc_w){                   // one of the both vectors have to be accumulated: really bad!
        for (int k=0; k<=i; ++k)
            tmpHCol[k] = ex.LocalDot( w, false, v[k], false);
    }
    else        // update of w do only works on same types
        throw DROPSErrCL("StandardGramSchmidt: Cannot do Gram Schmidt on that kind of vectors!");

    // Syncpoint!
    ProcCL::GlobalSum(Addr(tmpHCol), H.GetCol(i), i+1);
    for (int k=0; k<=i; ++k)
        w -= H(k,i)*v[k];
}


/// \brief Computes an orthogonal vector on i vectors by the modified Gram-Schmidt method
template <typename Vec, typename ExCL>
void ModifiedGramSchmidt(DMatrixCL<double>& H, Vec& w, bool acc_w, const std::vector<Vec>& v, bool acc_v, int i, const ExCL& ex)
/// \param H     Hessenberg matrix
/// \param w     orthogonal vector
/// \param acc_w is w accumulated
/// \param v     vectors used to orthogonalize w
/// \param acc_v is v accumulated
/// \param i     number of vectors v
/// \param ex    class to accumulate a vector
{
    if (acc_v&&acc_w){
        for (int k=0; k<=i; ++k){
            H( k, i)= ex.ParDot( w, true, v[k], true);
            w-= H( k, i)*v[k];
        }
    }
    else if (!acc_v&&!acc_w){
        for (int k=0; k<=i; ++k){
            H( k, i)= ex.ParDot( w, false, v[k], false);
            w-= H( k, i)*v[k];
        }
    }
    else
        throw DROPSErrCL("StandardGramSchmidt: Cannot do Gram Schmidt on that kind of vectors!");
}

/// \brief Parallel GMRES-method.
///
/// This method is the same algorithm as the serial algorithm. For performance issues take the ParModGMRES
/// procedure.
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
  bool ParGMRES(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M,
                int m, int& max_iter, double& tol,
                bool measure_relative_tol, PreMethGMRES method, std::ostream* output)
    /// \param[in]     A                    local distributed coefficients-matrix of the linear equation system
    /// \param[in,out] x_acc                start vector and the solution in accumulated form
    /// \param[in]     b                    rhs of the linear equation system (distributed form)
    /// \param[in]     ExX                  ExchangeCL corresponding to the RowIdx of x and the ColIdx of A
    /// \param[in,out] M                    Preconditioner
    /// \param[in]     m                    number of steps after a restart is performed
    /// \param[in,out] max_iter             IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol                  IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol if true stop if |M^(-1)(b-Ax)|/|M^(-1)b| <= tol, else stop if |M^(-1)(b-Ax)|<=tol
    /// \param[in]     method               left or right preconditioning (see solver.h for definition and declaration)
    /// \return  convergence within max_iter iterations
    /// \pre     the preconditioner should be able to handle a accumulated b
{
    // Check if preconditioner needs diagonal of matrix. The preconditioner
    // only computes the diagonal new, if the matrix has changed
    if (M.NeedDiag())
        M.SetDiag(A);

    m= (m <= max_iter) ? m : max_iter; // m > max_iter only wastes memory.

    // The preconditioner can have one of the following properties by calling M.Apply(A,x,b) with b distributed:
    // 1. x has distributed form    (like Dummy, Jacobi, BlockSSOR)
    // 2. x has accumulated form    (like SolverAsPreCL)
    // Unfortunately there are a lot of differences, so we have to implement two strategies

    if (!M.RetAcc()){                                   // case 1
        if (method==LeftPreconditioning)
            throw DROPSErrCL("ParGMRES: Left preconditioning does not work!");
        /// \todo (of) Left preconditioning funktioniert nicht!
        DMatrixCL<double> H( m, m);
        Vec               s( m), cs( m), sn( m),
                          w( b.size()), r( b.size());
        std::vector<Vec>  v( m);
        double            beta, normb, resid;
        Vec z(x_acc.size()), t(x_acc.size());
        for (int i= 0; i < m; ++i)
            v[i].resize( b.size());

        if (method == RightPreconditioning)
        {
            r= b - A*x_acc;
            beta = ExX.Norm(r, false);
            normb= ExX.Norm(b, false);
        }
        else
        {
            M.Apply( A, r, Vec( b - A*x_acc));
            beta= ExX.Norm(r, false);
            M.Apply( A, w, b);
            normb= ExX.Norm(w, false);
        }
        if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;

        resid = beta/normb;
        if (resid <= tol) {
            tol= resid;
            max_iter= 0;
            return true;
        }

        int j= 1;
        while (j <= max_iter) {
            v[0]= r*(1.0/beta);
            s= 0.0;
            s[0]= beta;

            int i;
            for (i= 0; i<m-1 && j<=max_iter; ++i, ++j) {
                if (method == RightPreconditioning)
                {
                    M.Apply( A, w, v[i]);
                    w=A*ExX.GetAccumulate(w);
                }
                else
                    M.Apply( A, w, A*ExX.GetAccumulate(v[i]));
                for (int k= 0; k <= i; ++k ) {
                    H( k, i)= ExX.ParDot( w, false, v[k], false);
                    w-= H( k, i)*v[k];
                }

                if (i == m - 1) break;

                H( i + 1, i)= ExX.Norm(w, false);
                v[i + 1]= w*(1.0/H( i + 1, i));

                for (int k= 0; k < i; ++k)
                    GMRES_ApplyPlaneRotation( H(k,i), H(k + 1, i), cs[k], sn[k]);

                GMRES_GeneratePlaneRotation( H(i,i), H(i+1,i), cs[i], sn[i]);
                GMRES_ApplyPlaneRotation( H(i,i), H(i+1,i), cs[i], sn[i]);
                GMRES_ApplyPlaneRotation( s[i], s[i+1], cs[i], sn[i]);

                resid= std::abs( s[i+1])/normb;
                if (output && j%10 == 0)
                    (*output) << "ParGMRES: " << j << " resid " << resid << std::endl;
                if (resid <= tol) {
                    if (method == RightPreconditioning)
                    {
                        z=0.;
                        GMRES_Update( z, i, H, s, v);
                        M.Apply( A, t, z);
                        x_acc+=ExX.GetAccumulate(t);
                    }
                    else{
                        z=0.;
                        GMRES_Update( z, i, H, s, ExX.GetAccumulate(v));
                        x_acc += z;
                    }
                    tol= resid;
                    max_iter= j;
                    return true;
                }
            }

            if (method == RightPreconditioning)
            {
                z=0.;
                GMRES_Update( z, i-1, H, s, v);
                M.Apply( A, t, z);
                x_acc+=ExX.GetAccumulate(t);
                r= b - A*x_acc;
            }
            else
            {
                GMRES_Update( x_acc, i-1, H, s, ExX.GetAccumulate(v));
                M.Apply( A, r, Vec( b - A*x_acc));
            }
            beta=ExX.Norm(r, false);
            resid= beta/normb;
            if (resid <= tol) {
                tol= resid;
                max_iter= j;
                return true;
            }
        }
        tol= resid;
        return false;
    }
    else                                                   // case 2: Apply returns accumulated vector
    {
        DMatrixCL<double> H( m, m);
        Vec               s( m), cs( m), sn( m),
                          w( b.size()), r( b.size());
        std::vector<Vec>  v( m);
        double            beta, normb, resid;
        Vec z(x_acc.size()), t(x_acc.size());
        for (int i= 0; i < m; ++i)
            v[i].resize( b.size());

        if (method == RightPreconditioning)
        {
            Vec r_tmp(b - A*x_acc);
            beta = ExX.Norm(r_tmp, false, &r);
            normb= ExX.Norm(b, false);
        }
        else
        {
            M.Apply( A, r, Vec( b - A*x_acc));
            beta = ExX.Norm(r, true);
            M.Apply( A, w, b);
            normb= ExX.Norm(w, true);
        }
        if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;

        resid = beta/normb;
        if (resid <= tol) {
            tol= resid;
            max_iter= 0;
            return true;
        }

        int j= 1;
        while (j <= max_iter) {
            v[0]= r*(1.0/beta);                             // v has accumulated form because r is accumulated
            s= 0.0;
            s[0]= beta;

            int i;
            for (i= 0; i<m-1 && j<=max_iter; ++i, ++j) {
                if (method == RightPreconditioning)
                {
                    M.Apply( A, w, v[i]);                   // hopefully, preconditioner do right things with accumulated v[i]
                    w=A*w;
                    ExX.Accumulate(w);
                }
                else
                    M.Apply( A, w, A*v[i]);
                for (int k= 0; k <= i; ++k ) {
                    H( k, i)= ExX.ParDot(w, true, v[k], true);
                    w-= H( k, i)*v[k];
                }

                if (i == m - 1) break;

                H( i + 1, i)= ExX.Norm(w, true);
                v[i + 1]= w*(1.0/H( i + 1, i));

                for (int k= 0; k < i; ++k)
                    GMRES_ApplyPlaneRotation( H(k,i), H(k + 1, i), cs[k], sn[k]);

                GMRES_GeneratePlaneRotation( H(i,i), H(i+1,i), cs[i], sn[i]);
                GMRES_ApplyPlaneRotation( H(i,i), H(i+1,i), cs[i], sn[i]);
                GMRES_ApplyPlaneRotation( s[i], s[i+1], cs[i], sn[i]);

                resid= std::abs( s[i+1])/normb;
                if (output)
                    (*output) << "ParGMRES: " << j << " resid " << resid << std::endl;

                if (resid <= tol) {
                    if (method == RightPreconditioning)
                    {
                        z=0.;
                        GMRES_Update( z, i, H, s, v);
                        M.Apply( A, t, z);
                        x_acc+=t;
                    }
                    else
                        GMRES_Update( x_acc, i, H, s, v);
                    tol= resid;
                    max_iter= j;
                    return true;
                }
            }

            if (method == RightPreconditioning)
            {
                z=0.;
                GMRES_Update( z, i-1, H, s, v);
                M.Apply( A, t, z);
                x_acc+=t;
                r= ExX.GetAccumulate(Vec(b - A*x_acc));
            }
            else
            {
                GMRES_Update( x_acc, i-1, H, s, v);
                M.Apply( A, r, Vec( b - A*x_acc));
            }
            beta=ExX.Norm(r, true);
            resid= beta/normb;
            if (resid <= tol) {
                tol= resid;
                max_iter= j;
                return true;
            }
        }
        tol= resid;
        return false;
    }
}


/// \brief Parallel GMRES-method with modifications for better scalability.
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
  bool ParModGMRES(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M,
                   int m, int& max_iter, double& tol,
                   bool measure_relative_tol, bool useMGS, PreMethGMRES method, std::ostream* output)
    /// \param[in]     A                    local distributed coefficients-matrix of the linear equation system
    /// \param[in,out] x_acc                start vector and the solution in accumulated form
    /// \param[in]     b                    rhs of the linear equation system (distributed form)
    /// \param[in]     ExX                  ExchangeCL corresponding to the RowIdx of x and the ColIdx of A
    /// \param[in,out] M                    Preconditioner
    /// \param[in]     m                    number of steps after a restart is performed
    /// \param[in,out] max_iter             IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol                  IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol if true stop if |M^(-1)(b-Ax)|/|M^(-1)b| <= tol, else stop if |M^(-1)(b-Ax)|<=tol
    /// \param[in]     useMGS               use modified Gram-Schmidt ortogonalization (many more sync-points exists!)
    /// \param[in]     method               left or right preconditioning (see solver.h for definition and declaration)
    /// \return  convergence within max_iter iterations
    /// \pre     the preconditioner should be able to handle a accumulated b
{
    Assert(x_acc.size()==b.size() && x_acc.size()==ExX.GetNum(), DROPSErrCL("ParModGMRES: Incompatible dimension"), DebugParallelNumC);

    // Check if preconditioner needs diagonal of matrix. The preconditioner
    // only computes the diagonal new, if the matrix has changed
    if (M.NeedDiag())
        M.SetDiag(A);

    const size_t n = x_acc.size();      // dimension

    DMatrixCL<double> H(m,m);           // upper Hessenberg-matrix
    VectorCL tmpHCol(m);
    double beta, normb, resid;

    Vec r(n), w(n), w_acc(n), r_acc(n), z_acc(n), t_acc(n);
    Vec c(m), s(m), gamma(m);

    std::vector<Vec> v_acc(m);    // basis of the krylov-subspaces
    for (int i=0; i<m; ++i)
        v_acc[i].resize(n);

    if (method == RightPreconditioning){
        r    = b - A*x_acc;
        beta = ExX.Norm(r, false, &r_acc);
        normb= ExX.Norm(b, false);
    }
    else{
        M.Apply(A, r, VectorCL( b-A*x_acc));
        beta = ExX.Norm(r, M.RetAcc(), &r_acc);
        M.Apply(A, w, b);
        normb = ExX.Norm(w, M.RetAcc(), &w_acc);
    }

    if (normb == 0. || measure_relative_tol==false) normb=1.0;

    resid = beta/normb;
    if (resid<=tol){                        // finished
        tol = resid; max_iter = 0; return true;
    }

    int j=1;                                // number of steps
    while (j<=max_iter)
    {
        v_acc[0] = r_acc * (1./beta);
        gamma    = 0.;
        gamma[0] = beta;

        int i;
        for (i=0; i<m-1 && j<=max_iter; ++i, ++j)
        {
            if (method == RightPreconditioning){
                M.Apply(A, w_acc, v_acc[i]);                // hopefully M does the right thing
                w = A*w_acc;
                w_acc = ExX.GetAccumulate(w);
            }
            else{
                M.Apply(A, w, A*v_acc[i]);
                if (!M.RetAcc())
                    w_acc = ExX.GetAccumulate(w);
                else
                    w_acc = w;
            }

            // Gram-Schmidt ortogonalization without  update of w!
            if (!useMGS)
                StandardGramSchmidt(H, w_acc, true, v_acc, true, i, ExX, tmpHCol);
            else
                ModifiedGramSchmidt(H, w_acc, true, v_acc, true, i, ExX);

            H(i+1,i) = ExX.Norm(w_acc, true);
            v_acc[i+1] = w_acc * (1.0 / H(i+1,i));

            for (int k=0; k<i; ++k)
                GMRES_ApplyPlaneRotation(H(k,i),H(k+1,i), c[k], s[k]);

            GMRES_GeneratePlaneRotation(H(i,i), H(i+1,i), c[i], s[i]);
            GMRES_ApplyPlaneRotation(H(i,i), H(i+1,i), c[i], s[i]);
            GMRES_ApplyPlaneRotation(gamma[i], gamma[i+1], c[i], s[i]);

            resid = std::abs(gamma[i+1])/normb;
            if (output)
                (*output) << "ParModGMRES: " << j << " resid " << resid << std::endl;

            if (resid<=tol){            // finished
                if (method == RightPreconditioning){
                    z_acc=0.;
                    GMRES_Update( z_acc, i, H, gamma, v_acc);
                    M.Apply( A, t_acc, z_acc);              // hopefully M does the right thing
                    x_acc+=t_acc;
                }
                else
                    GMRES_Update(x_acc, i, H, gamma, v_acc);
                tol = resid; max_iter = j; return true;
            }
        }
        if (method == RightPreconditioning){
            z_acc=0.;
            GMRES_Update( z_acc, i-1, H, gamma, v_acc);
            M.Apply( A, t_acc, z_acc);                      // hopefully M does the right thing
            x_acc += t_acc;
            r      = b-A*x_acc;
            beta = ExX.Norm(r, false, &r_acc);
        }
        else{
            GMRES_Update(x_acc, i-1, H, gamma, v_acc);
            M.Apply(A, r, static_cast<Vec>( b-A*x_acc));
            beta = ExX.Norm(r, M.RetAcc(), &r_acc);
        }


        resid = beta/normb;
        if (resid<=tol){                // finished
            tol = resid; max_iter=j; return true;
        }
    }
    tol = std::fabs(gamma[0]);
    return false;
}


/// \brief Preconditioned BiCGStab-Method with accure inner products
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
  bool ParBiCGSTAB(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX,
                   PreCon& M, int& max_iter, double& tol, bool measure_relative_tol)
    /// \param[in]     A                    local distributed coefficients-matrix of the linear equation system
    /// \param[in,out] x_acc                start vector and the solution in accumulated form
    /// \param[in]     b                    rhs of the linear equation system (distributed form)
    /// \param[in]     ExX                  ExchangeCL corresponding to the RowIdx of x and the ColIdx of A
    /// \param[in]     M                    Preconditioner
    /// \param[in,out] max_iter             IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol                  IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol if true stop if |M^(-1)(b-Ax)|/|M^(-1)b| <= tol, else stop if |M^(-1)(b-Ax)|<=tol
    /// \pre                                the preconditioner must fulfill the following condition: M.Apply(A,x,b): b is accumulated => x is accumulated
    /// \return                             convergence within max_iter iterations
{
    // Check if preconditioner needs diagonal of matrix. The preconditioner
    // only computes the diagonal new, if the matrix has changed
    if (M.NeedDiag())
        M.SetDiag(A);

    const size_t n = x_acc.size();
    Vec r(b-A*x_acc), r_acc(n), r0hat_acc(n),
        p_acc(n), phat_acc(n),
        v(n), v_acc(n),
        s_acc(n), shat_acc(n),
        t(n), t_acc(n), that(n), that_acc(n);

    double resid, rho=1, alpha=1, beta, omega=1, sigma;

    VectorCL dots(4), glob_dots(4);

    double normb = ExX.Norm(b, false);
    if ( normb==0. || !measure_relative_tol) normb= 1.;

    sigma = ExX.Norm(r, false, &r_acc);
    resid = std::sqrt(sigma) / normb;

    M.Apply(A, r0hat_acc, r);
    if (!M.RetAcc())
        r0hat_acc= ExX.GetAccumulate(r0hat_acc);

    for (int i=0; i<max_iter; ++i)
    {
        beta = (sigma*alpha)/(rho*omega);
        rho  = sigma;
        if (i>0)
            p_acc    = r_acc + beta*p_acc - (beta*omega)*v_acc;
        else
            p_acc    = r_acc;

        M.Apply(A, phat_acc, p_acc);

        v= A*phat_acc;

        dots[0]= ExX.LocalDot( v, false, r0hat_acc, true, &v_acc);
        dots[1]= ExX.LocalNorm_sq( r_acc, true);

        ProcCL::GlobalSum(Addr(dots), Addr(glob_dots), 2);

        resid = std::sqrt(glob_dots[1])/normb;
        if (resid<tol){
            tol= resid;
            max_iter=i;
            return true;
        }

        sigma = glob_dots[0];

        if (glob_dots[0]==0.)
        {
            if (ProcCL::MyRank()==0)
                std::cout << ">>>>> BREAKDOWN of BiCGStab!" <<std::endl;
            tol = resid;
            max_iter=i;
            return false;
        }

        alpha = rho/sigma;
        s_acc = r_acc -alpha*v_acc;

        M.Apply(A, shat_acc, s_acc);
        t = A*shat_acc;

        if (!M.RetAcc()){
            M.Apply(A,that,t);
            dots[0]= ExX.LocalDot( that, false, shat_acc, true, &that_acc);
        }
        else{
            M.Apply(A,that_acc,t);
            dots[0]= ExX.LocAccDot(that_acc, shat_acc);
        }

        dots[1]= ExX.LocAccDot(that_acc, that_acc);
        dots[2]= ExX.LocAccDot(r0hat_acc, s_acc);
        dots[3]= ExX.LocalDot( t, false, r0hat_acc, true, &t_acc);

        ProcCL::GlobalSum(Addr(dots), Addr(glob_dots), 4);

        if (glob_dots[1]==0.)
        {
            if (ProcCL::MyRank()==0)
                std::cout << ">>>>> BREAKDOWN of BiCGStab!" <<std::endl;
            tol = resid;
            max_iter=i;
            return false;
        }

        omega = glob_dots[0]/glob_dots[1];
        sigma = glob_dots[2] - omega * glob_dots[3];

        r_acc =  s_acc - omega * t_acc;
        x_acc += alpha * phat_acc + omega * shat_acc;
    }

    tol = resid;
    return false;
}


/// \brief Preconditioned QMR-Method
template <typename Mat, typename Vec, typename Lanczos, typename ExCL>
bool ParQMR(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, Lanczos lan,
            int& max_iter, double& tol, bool measure_relative_tol)
    /// \param[in]     A                    local distributed coefficients-matrix of the linear equation system
    /// \param[in,out] x_acc                start vector and the solution in accumulated form
    /// \param[in]     b                    rhs of the linear equation system (distributed form)
    /// \param[in]     ExX                  ExchangeCL corresponding to the RowIdx of x and the ColIdx of A
    /// \param[in]     lan                  Lanczos-Class for computing a bi-orthogonal basis of Krylov subspaces
    /// \param[in,out] max_iter             IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol                  IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol if true stop if |M^(-1)(b-Ax)|/|M^(-1)b| <= tol, else stop if |M^(-1)(b-Ax)|<=tol
    /// \return                             convergence within max_iter iterations
{
    tol*=tol;

    Vec d_acc(x_acc.size()), s(d_acc), r(b-A*x_acc), r_acc(r);

    double normb = ExX.Norm_sq(b, false);
    if (normb==0.0 || measure_relative_tol==false) normb=1.0;
    double norm_r = ExX.Norm( r, false, &r_acc);

    if (norm_r/normb<std::sqrt(tol))
    {
        tol = std::fabs(norm_r)/std::sqrt(normb);
        max_iter=1;
        return true;
    }

    lan.Init(A,static_cast<Vec>(r*(1./norm_r)),static_cast<Vec>(r*(1./norm_r)),ExX,1);
    double tetha,
    kappa=-1,
    lambda=1,
    beta,
    rho = norm_r,
    rho_last = norm_r;
//  rho_last=lan.GetRho();
    double quot;

    for (int j=1; j<=max_iter; ++j)
    {
        if (!lan.Step()){
            tol = std::sqrt(norm_r/normb);
            max_iter=j-1;
            return false;
        }

        beta    = lan.GetBeta();
        rho_last= rho;
        rho     = lan.GetRho();
        quot    = (lambda*beta*beta+rho*rho);
        tetha   = beta*beta * (1-lambda)/quot;
        kappa   = -rho_last*beta*kappa/quot;
        lambda  = lambda*beta*beta/quot;

        d_acc   = tetha * d_acc + kappa*lan.GetAccP();
        s       = tetha * s     + kappa*lan.GetAp();
        x_acc  += d_acc;
        r      -= s;

        norm_r  = ExX.Norm(r, false, &r_acc);

        if (norm_r<0)
            std::cout << "["<<ProcCL::MyRank()<<"]==> negative squared norm of resid in QMR because of accumulation!" << std::endl;
        if (norm_r/normb<tol)
        {
            tol = std::sqrt(std::fabs(norm_r)/normb);
            max_iter=j;
            return true;
        }
    }
    tol = std::sqrt(norm_r);
    return false;
}


template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParPGCR(const Mat& A, Vec& x, const Vec& b, const ExCL& ExX, PreCon& M,
    int m, int& max_iter, double& tol, bool measure_relative_tol, std::ostream* output)
{
    if (M.NeedDiag())
        M.SetDiag(A);

    m= (m <= max_iter) ? m : max_iter; // m > max_iter only wastes memory.

    Vec r( b - A*x);
    Vec racc( r.size());
    Vec sn( b.size()), vn( b.size());
    Vec vnacc( b.size());
    std::vector<Vec> s, v;
    std::vector<Vec> vacc;         // memory usage vs communication overhead
    std::vector<double> a( m);

    double normb= ExX.Norm( b, false);
    if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;
    double resid= ExX.Norm( r, false, &racc)/normb;
    for (int k= 0; k < max_iter; ++k) {
        if (k%1==0 && output)
            (*output) << "GCR: k: " << k << "\tresidual: " << resid << std::endl;
        if (resid < tol) {
            tol= resid;
            max_iter= k;
            return true;
        }
        M.Apply( A, sn, r);
        if (!M.RetAcc())
            vn=A*ExX.GetAccumulate(sn);
        else
            vn= A*sn;
        vnacc = ExX.GetAccumulate(vn);
        for (int i= 0; i < k && i < m; ++i) {
//            const double alpha= ExX.ParDot( vn, false, v[i], false);
            const double alpha= ExX.ParDot( vnacc, true, vacc[i], true);
            a[i]= alpha;
            vn-= alpha*v[i];
            vnacc-= alpha*vacc[i];  // memory usage vs communication overhead
            sn-= alpha*s[i];
        }
//        const double beta= ExX.Norm( vn, false, &vnacc);
        const double beta= ExX.Norm( vnacc, true);
        vn/= beta;
        vnacc/= beta;
        sn/= beta;
        const double gamma= ExX.ParDot( racc, true, vnacc, true);
        if (!M.RetAcc())
            x+= gamma*ExX.GetAccumulate(sn);
        else
            x+= gamma*sn;
        r-= gamma*vn;
        racc -= gamma*vnacc;
        resid= ExX.Norm( racc, true)/normb;
        if (k < m) {
            s.push_back( sn);
            v.push_back( vn);
            vacc.push_back( vnacc);   // memory usage vs communication overhead
        }
        else {
            throw DROPSErrCL("ParPGCR: Sorry, truncation not implemented");
            int min_idx= 0;
            double a_min= std::fabs( a[0]); // m >= 1, thus this access is valid.
            for (int i= 1; i < k && i < m; ++i)
                if ( std::fabs( a[i]) < a_min) {
                    min_idx= i;
                    a_min= std::fabs( a[i]);
                }
            s[min_idx]= sn;
            v[min_idx]= vn;
        }
    }
    tol= resid;
    return false;
}

/// \brief PCGNE: Preconditioned CG for the normal equations (error-minimization)
///
/// Solve A*A^T x = b with left preconditioner M. This is more stable than PCG with
/// a CompositeMatrixCL.
///
/// The return value indicates convergence within max_iter (input)
/// iterations (true), or no convergence within max_iter iterations (false).
/// Upon successful return, output arguments have the following values:
///
/// \param A - matrix (not necessarily quadratic)
/// \param b - right hand side
/// \param u - approximate solution to A*A^T u = b
/// \param M - preconditioner
/// \param max_iter - number of iterations performed before tolerance was reached
/// \param tol - 2-norm of the (relative, see below) residual after the final iteration
/// \param measure_relative_tol - If true, stop if |b - A*A^T u|/|b| <= tol,
///        if false, stop if |b - A*A^T u| <= tol.
template <typename Mat, typename Vec, typename PreCon, typename ExACL, typename ExATranspCL>
bool
ParPCGNE(const Mat& A, Vec& u, const Vec& b, const ExACL& ExAX,  const ExATranspCL& ExATranspX, const PreCon& M,
    int& max_iter, double& tol, bool measure_relative_tol= false, std::ostream* output=0)
{
    Vec Atranspu( transp_mul( A, u));
    ExAX.Accumulate( Atranspu);
    Vec r( b - A*Atranspu);

    double normb= ExATranspX.Norm( b, false);
    if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;

    double resid= ExATranspX.Norm( r, false)/normb;
    if (output)
        (*output) << "PCGNE: iter: 0 resid: " << resid <<'\n';
    if (resid <= tol) {
        tol= resid;
        max_iter= 0;
        return true;
    }

    const size_t n= A.num_rows();
    const size_t num_cols= A.num_cols();

    Vec z( n);
    M.Apply( A, z, r);
    Vec pt( num_cols);
    Vec qtacc( z), ptacc( pt), racc(r), zacc(z);
    
    double rho= ExATranspX.ParDot( z, false, r, false, &qtacc), rho_1; // M returns non accumulated vector ???

    for (int i= 1; i <= max_iter; ++i) {
        pt= transp_mul( A, qtacc);
        const double alpha= rho/ExAX.Norm_sq( pt, false, &ptacc);
        u+= alpha*qtacc;
        r-= alpha*(A*ptacc);
        M.Apply( A, z, r); // M returns non accumulated vector ???

        resid= ExATranspX.Norm( r, false, &racc)/normb;
        if ( output && i%10 == 0) (*output) << "PCGNE: iter: " << i << " resid: " << resid <<'\n';
        if (resid <= tol) {
            tol= resid;
            max_iter= i;
            return true;
        }
        rho_1= rho;
        rho= ExATranspX.ParDot( z, false, racc, true, &zacc);
        qtacc= zacc + (rho/rho_1)*qtacc;
    }
    tol= resid;
    return false;
}

/// \brief (Symmetric) Gauss-Seidel preconditioner (i.e. start-vector 0) for A*A^T (Normal Equations).
class ParNEGSPcCL
{
  private:
    const IdxDescCL& RowIdx_;
    const IdxDescCL& ColIdx_;
    bool symmetric_;             ///< If true, SGS is performed, else GS.
    mutable const void*  Aaddr_; ///< only used to validate, that the diagonal is for the correct matrix.
    mutable size_t Aversion_;
    mutable VectorCL D_;         ///< diagonal of AA^T
    mutable VectorCL y_;         ///< temp-variable: A^T*x.

    template <typename Mat>
    void Update (const Mat& A) const;

    ///\brief Forward Gauss-Seidel-step with start-vector 0 for A*A^T
    template <typename Mat, typename Vec>
    void ForwardGS(const Mat& A, Vec& x, const Vec& b) const;
    ///\brief Backward Gauss-Seidel-step with start-vector 0 for A*A^T
    template <typename Mat, typename Vec>
    void BackwardGS(const Mat& A, Vec& x, const Vec& b) const;

    ///\brief Inverse of ForwardGS
    template <typename Mat, typename Vec>
    void ForwardMulGS(const Mat& A, Vec& x, const Vec& b) const;
    ///\brief Inverse of BackwardGS
    template <typename Mat, typename Vec>
    void BackwardMulGS(const Mat& A, Vec& x, const Vec& b) const;

  public:
    ParNEGSPcCL (const IdxDescCL& Row, const IdxDescCL& Col, bool symmetric= true) : RowIdx_(Row), ColIdx_(Col), symmetric_( symmetric), Aaddr_( 0), Aversion_( 0) {}

    ///@{ Note, that A and not A*A^T is the first argument.
    ///\brief Execute a (symmetric) Gauss-Seidel preconditioning step.
    template <typename Mat, typename Vec>
    void Apply(const Mat& A, Vec& x, const Vec& b) const;
    ///\brief If symmetric == false, this  performs a backward Gauss-Seidel step, else it is identical to Apply.
    template <typename Mat, typename Vec>
    void ApplyTranspose(const Mat& A, Vec& x, const Vec& b) const;

    ///\brief Multiply with the preconditioning matrix -- needed for right preconditioning.
    template <typename Mat, typename Vec>
    Vec mul (const Mat& A, const Vec& b) const;
    ///\brief Multiply with the transpose of the preconditioning matrix -- needed for right preconditioning.
    template <typename Mat, typename Vec>
    Vec transp_mul(const Mat& A, const Vec& b) const;
    ///@}

    ///\brief Apply if A*A^T is given as CompositeMatrixCL
    void Apply(const CompositeMatrixCL& AAT, VectorCL& x, const VectorCL& b) const
    { Apply( *AAT.GetBlock1(), x, b); }
};

template <typename Mat>
void ParNEGSPcCL::Update (const Mat& A) const
{
    if (&A == Aaddr_ && Aversion_ == A.Version()) return;
    Aaddr_= &A;
    Aversion_= A.Version();

    D_.resize( A.num_rows());
    // Accumulate matrix
    ExchangeMatrixCL exMat_;
    exMat_.BuildCommPattern(A, RowIdx_.GetEx(), ColIdx_.GetEx());
    MatrixCL Aacc(exMat_.Accumulate(A));
    // Determine diagonal of AA^T
    D_.resize(Aacc.num_rows());
    for (size_t i = 0; i < Aacc.num_rows(); ++i)
        for (size_t nz = Aacc.row_beg(i); nz < Aacc.row_beg(i + 1); ++nz)
            D_[i] += Aacc.val(nz) * A.val(nz);
    RowIdx_.GetEx().Accumulate(D_);
    y_.resize( A.num_cols());
}

template <typename Mat, typename Vec>
void ParNEGSPcCL::ForwardGS(const Mat& A, Vec& x, const Vec& b) const
{
    // x= 0.; // implied, but superfluous, because we can assign the x-values below, not update.
    y_= 0.;
    double t;
    const ExchangeCL& RowEx_(RowIdx_.GetEx());

    for (size_t i= 0; i < RowEx_.LocalIndex.size(); ++i) {
        t= (b[RowEx_.LocalIndex[i]] - mul_row( A, y_, RowEx_.LocalIndex[i]))/D_[RowEx_.LocalIndex[i]];
        x[RowEx_.LocalIndex[i]]= t;
        add_row_to_vec( A, t, y_, RowEx_.LocalIndex[i]); // y+= t* (i-th row of A)
    }
    for (size_t i= 0; i < RowEx_.DistrIndex.size(); ++i) {
        t= (b[RowEx_.DistrIndex[i]])/D_[RowEx_.DistrIndex[i]];
        x[RowEx_.DistrIndex[i]]= t;
        add_row_to_vec( A, t, y_, RowEx_.DistrIndex[i]); // y+= t* (i-th row of A)
    }
}

template <typename Mat, typename Vec>
void ParNEGSPcCL::BackwardGS(const Mat& A, Vec& x, const Vec& b) const
{
    // x= 0.; // implied, but superfluous, because we can assign the x-values below, not update.
    y_= 0.;
    const ExchangeCL& RowEx_(RowIdx_.GetEx());

    for (size_t i= RowEx_.LocalIndex.size() - 1; i < RowEx_.LocalIndex.size(); --i) {
        x[RowEx_.LocalIndex[i]]= (b[RowEx_.LocalIndex[i]] - mul_row( A, y_, RowEx_.LocalIndex[i]))/D_[RowEx_.LocalIndex[i]];
        add_row_to_vec( A, x[RowEx_.LocalIndex[i]], y_, RowEx_.LocalIndex[i]); // y+= t* (i-th row of A)
    }
    for (size_t i= RowEx_.DistrIndex.size() - 1; i < RowEx_.DistrIndex.size(); --i) {
        x[RowEx_.DistrIndex[i]]= (b[RowEx_.DistrIndex[i]] )/D_[RowEx_.DistrIndex[i]];
        add_row_to_vec( A, x[RowEx_.DistrIndex[i]], y_, RowEx_.DistrIndex[i]); // y+= t* (i-th row of A)
    }
}

template <typename Mat, typename Vec>
void ParNEGSPcCL::Apply(const Mat& A, Vec& x, const Vec& b) const
{
    Update( A);

    ForwardGS( A, x, b);
    if (!symmetric_) return;

    BackwardGS( A, x, VectorCL( D_*x));
}

template <typename Mat, typename Vec>
void ParNEGSPcCL::ApplyTranspose(const Mat& A, Vec& x, const Vec& b) const
{
    Update( A);

    BackwardGS( A, x, b);
    if (!symmetric_) return;
    ForwardGS( A, x, VectorCL( D_*x));
}

template <typename Mat, typename Vec>
void ParNEGSPcCL::ForwardMulGS(const Mat& A, Vec& x, const Vec& b) const
{
    // x= 0.; // implied, but superfluous, because we can assign the x-values below, not update.
    y_= 0.;
    for (size_t i= 0 ; i < b.size(); ++i) {
        add_row_to_vec( A, b[i], y_, i); // y+= b[i]* (i-th row of A)
        x[i]= mul_row( A, y_, i);
    }
}

template <typename Mat, typename Vec>
void ParNEGSPcCL::BackwardMulGS(const Mat& A, Vec& x, const Vec& b) const
{
    // x= 0.; // implied, but superfluous, because we can assign the x-values below, not update.
    y_= 0.;
    for (size_t i= b.size() - 1 ; i < b.size(); --i) {
        add_row_to_vec( A, b[i], y_, i); // y+= b[i]* (i-th row of A)
        x[i]= mul_row( A, y_, i);
    }
}

template <typename Mat, typename Vec>
Vec ParNEGSPcCL::mul (const Mat& A, const Vec& b) const
{
    Update( A);

    Vec x( A.num_rows());
    VectorCL b2( b);
    if (symmetric_) {
        BackwardMulGS( A, x, b);
        b2= x/D_;
    }
    ForwardMulGS( A, x, b2);
    return x;
}

template <typename Mat, typename Vec>
Vec ParNEGSPcCL::transp_mul(const Mat& A, const Vec& b) const
{
    Update( A);

    Vec x( A.num_rows());
    VectorCL b2( b);
    if (symmetric_) {
        ForwardMulGS( A, x, b);
        b2= x/D_;
    }
    BackwardMulGS( A, x, b2);
    return x;
}

} // end of namespace DROPS

#endif
