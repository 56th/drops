/// \file oseensolver.h
/// \brief Solvers for the Oseen problem
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen:

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

#ifndef DROPS_STOKESSOLVER_H
#define DROPS_STOKESSOLVER_H

#include "num/krylovsolver.h"
#include "num/precond.h"
#include "num/MGsolver.h"
#include "num/oseenprecond.h"
#ifdef _PAR
#  include "parallel/exchange.h"
#endif
namespace DROPS
{

//=============================================================================
//  The Stokes solvers solve systems of the form
//    A v + BT p = b
//    B v + C p  = c
//  where C may be an empty matrix
//=============================================================================

// What every iterative stokes-solver should have
class StokesSolverBaseCL: public SolverBaseCL
{
  public:
    StokesSolverBaseCL (int maxiter, double tol, bool rel= false, std::ostream* output= 0)
        : SolverBaseCL(maxiter, tol, rel, output){}

#ifdef _PAR
    virtual void Solve( const MatrixCL& A, const MatrixCL& B, const MatrixCL& C, VectorCL& v, VectorCL& p,
                        const VectorCL& b, const VectorCL& c, const ExchangeCL& vel_ex, const ExchangeCL& pr_ex) = 0;
    virtual void Solve( const MLMatrixCL& A, const MLMatrixCL& B, const MLMatrixCL& C, VectorCL& v, VectorCL& p,
                        const VectorCL& b, const VectorCL& c, const ExchangeCL& vel_ex, const ExchangeCL& pr_ex) = 0;
#endif    
    virtual void Solve( const MatrixCL&, const MatrixCL&, const MatrixCL&, VectorCL&, VectorCL&,
                        const VectorCL&, const VectorCL&, const DummyExchangeCL&, const DummyExchangeCL&) = 0;

    virtual void Solve( const MLMatrixCL&, const MLMatrixCL&, const MLMatrixCL&, VectorCL&, VectorCL&,
                        const VectorCL&, const VectorCL&, const DummyExchangeCL&, const DummyExchangeCL&) = 0;

};

//=============================================================================
// Ready-to-use solver-class with inexact Uzawa method InexactUzawa.
//=============================================================================

// Characteristics of the preconditioner for the A-block
enum InexactUzawaApcMethodT { APC_OTHER, APC_SYM, APC_SYM_LINEAR };

//=============================================================================
// Inexact Uzawa-method from "Fast Iterative Solvers for Discrete Stokes
// Equations", Peters, Reichelt, Reusken, Chapter 3.3.
// The preconditioner Apc for A must be "good" (MG-like) to guarantee
// convergence.
//=============================================================================
template <class ApcT, class SpcT, InexactUzawaApcMethodT ApcMeth= APC_OTHER>
  class InexactUzawaCL: public StokesSolverBaseCL
{
  private:
    ApcT& Apc_;
    SpcT& Spc_;
    double innerreduction_;
    int    innermaxiter_;

  public:
    InexactUzawaCL(ApcT& Apc, SpcT& Spc, int outer_iter, double outer_tol,
        double innerreduction= 0.3, int innermaxiter= 500)
        :StokesSolverBaseCL( outer_iter, outer_tol),
         Apc_( Apc), Spc_( Spc), innerreduction_( innerreduction), innermaxiter_( innermaxiter)
    {}
    template<class Mat, class Vec, class ExT>
    inline
    void Solve(const Mat& A, const Mat& B, const Mat& C, Vec& v, Vec& p,
          const Vec& b, const Vec& c, const ExT& vel_ex, const ExT& pr_ex);

#ifdef _PAR
    inline
    void Solve(const MatrixCL& A, const MatrixCL& B, const MatrixCL& C, VectorCL& v, VectorCL& p,
          const VectorCL& b, const VectorCL& c, const ExchangeCL& vel_ex, const ExchangeCL& pr_ex);
    inline
    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, const MLMatrixCL& C, VectorCL& v, VectorCL& p,
          const VectorCL& b, const VectorCL& c, const ExchangeCL& vel_ex, const ExchangeCL& pr_ex);
#endif
    inline
    void Solve(const MatrixCL& A, const MatrixCL& B, const MatrixCL& C, VectorCL& v, VectorCL& p,
          const VectorCL& b, const VectorCL& c, const DummyExchangeCL& vel_ex, const DummyExchangeCL& pr_ex);
    inline
    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, const MLMatrixCL& C, VectorCL& v, VectorCL& p,
          const VectorCL& b, const VectorCL& c, const DummyExchangeCL& vel_ex, const DummyExchangeCL& pr_ex);
};

template <class ApcT, class SpcT, InexactUzawaApcMethodT Apcmeth>
template <class Mat, class Vec, class ExT>
  inline void
  InexactUzawaCL<ApcT, SpcT, Apcmeth>::Solve( const Mat& A, const Mat& B, const Mat& C, Vec& v, Vec& p,
          const Vec& b, const Vec& c, const ExT& vel_ex, const ExT& pr_ex)
{
    res_=  tol_;
    iter_= maxiter_;
    InexactUzawa( A, B, C, v, p, b, c, vel_ex, pr_ex, Apc_, Spc_, iter_, res_, Apcmeth, innerreduction_, innermaxiter_, rel_, output_);
}

template <class ApcT, class SpcT, InexactUzawaApcMethodT Apcmeth>
  inline void
  InexactUzawaCL<ApcT, SpcT, Apcmeth>::Solve( const MatrixCL& A, const MatrixCL& B, const MatrixCL& C,
    VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c, const DummyExchangeCL& vel_ex, const DummyExchangeCL& pr_ex)
{
    Solve<>(A, B, C, v, p, b, c, vel_ex, pr_ex);
}

template <class ApcT, class SpcT, InexactUzawaApcMethodT Apcmeth>
  inline void
  InexactUzawaCL<ApcT, SpcT, Apcmeth>::Solve( const MLMatrixCL& A, const MLMatrixCL& B, const MLMatrixCL& C,
    VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c, const DummyExchangeCL& vel_ex, const DummyExchangeCL& pr_ex)
{
    Solve<>(A, B, C, v, p, b, c, vel_ex, pr_ex);
}

#ifdef _PAR
template <class ApcT, class SpcT, InexactUzawaApcMethodT Apcmeth>
  inline void
  InexactUzawaCL<ApcT, SpcT, Apcmeth>::Solve( const MatrixCL& A, const MatrixCL& B, const MatrixCL& C,
    VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c, const ExchangeCL& vel_ex, const ExchangeCL& pr_ex)
{
    Solve<>(A, B, C, v, p, b, c, vel_ex, pr_ex);
}

template <class ApcT, class SpcT, InexactUzawaApcMethodT Apcmeth>
  inline void
  InexactUzawaCL<ApcT, SpcT, Apcmeth>::Solve( const MLMatrixCL& A, const MLMatrixCL& B, const MLMatrixCL& C,
    VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c, const ExchangeCL& vel_ex, const ExchangeCL& pr_ex)
{
    Solve<>(A, B, C, v, p, b, c, vel_ex, pr_ex);
}
#endif

// Use a Krylow-method (from num/krylovsolver.h) with the standard-interface of
// the Oseensolvers in this file.
template <class SolverT>
class BlockMatrixSolverCL: public StokesSolverBaseCL
{
  private:
    SolverT& solver_;

  public:
    BlockMatrixSolverCL( SolverT& solver)
        : StokesSolverBaseCL(-1, -1.0), solver_( solver) {}

// We overwrite these functions.
    void   SetTol     (double tol) { solver_.SetTol( tol); }
    void   SetMaxIter (int iter)   { solver_.SetMaxIter( iter); }
    void   SetRelError(bool rel)   { solver_.SetRelError( rel); }
    void   SetOutput  (std::ostream* output) {solver_.SetOutput( output); }

    double GetTol     () const { return solver_.GetTol(); }
    int    GetMaxIter () const { return solver_.GetMaxIter(); }
    double GetResid   () const { return solver_.GetResid(); }
    int    GetIter    () const { return solver_.GetIter(); }
    bool   GetRelError() const { return solver_.GetRelError(); }
#ifdef _PAR
    template <typename Mat, typename Vec>
    void
    Solve(const Mat& A, const Mat& B, const Mat& C, Vec& v, Vec& p, const Vec& b, const Vec& c, const ExchangeCL& vel_ex, const ExchangeCL& pr_ex) {
        BlockMatrixBaseCL<Mat> M( &A, MUL, &B, TRANSP_MUL, &B, MUL, &C, MUL);
        VectorCL rhs( M.num_rows());
        rhs[std::slice( 0, M.num_rows( 0), 1)]= b;
        rhs[std::slice( M.num_rows( 0), M.num_rows( 1), 1)]= c;
        VectorCL x( M.num_cols());
        x[std::slice( 0, M.num_cols( 0), 1)]= v;
        x[std::slice( M.num_cols( 0), M.num_cols( 1), 1)]= p;
        ExchangeBlockCL exBlock;
        exBlock.AttachTo(vel_ex);
        exBlock.AttachTo(pr_ex);
        solver_.Solve( M, x, rhs, exBlock);
        v= x[std::slice( 0, M.num_cols( 0), 1)];
        p= x[std::slice( M.num_cols( 0), M.num_cols( 1), 1)];
    }
#endif

    template <typename Mat, typename Vec>
    void
    Solve(const Mat& A, const Mat& B, const Mat& C, Vec& v, Vec& p, const Vec& b, const Vec& c, const DummyExchangeCL& vel_ex, const DummyExchangeCL& pr_ex) {
        BlockMatrixBaseCL<Mat> M( &A, MUL, &B, TRANSP_MUL, &B, MUL, &C, MUL );
        VectorCL rhs( M.num_rows());
        rhs[std::slice( 0, M.num_rows( 0), 1)]= b;
        rhs[std::slice( M.num_rows( 0), M.num_rows( 1), 1)]= c;
        VectorCL x( M.num_cols());
        x[std::slice( 0, M.num_cols( 0), 1)]= v;
        x[std::slice( M.num_cols( 0), M.num_cols( 1), 1)]= p;
        DummyExchangeBlockCL exBlock;
        exBlock.AttachTo(vel_ex);
        exBlock.AttachTo(pr_ex);
        solver_.Solve( M, x, rhs, exBlock);
        v= x[std::slice( 0, M.num_cols( 0), 1)];
        p= x[std::slice( M.num_cols( 0), M.num_cols( 1), 1)];
    }
#ifdef _PAR
    void Solve( const MatrixCL& A, const MatrixCL& B, const MatrixCL& C, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c, const ExchangeCL& vel_ex, const ExchangeCL& pr_ex) {
        Solve<>(A, B, C, v, p, b, c, vel_ex, pr_ex);
    }
    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, const MLMatrixCL& C, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c, const ExchangeCL& vel_ex, const ExchangeCL& pr_ex) {
        Solve<>(A, B, C, v, p, b, c, vel_ex, pr_ex);
    }
#endif
    void Solve( const MatrixCL& A, const MatrixCL& B, const MatrixCL& C, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c, const DummyExchangeCL& vel_ex, const DummyExchangeCL& pr_ex) {
        Solve<>(A, B, C, v, p, b, c, vel_ex, pr_ex);
    }

    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, const MLMatrixCL& C, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c, const DummyExchangeCL& vel_ex, const DummyExchangeCL& pr_ex) {
        Solve<>(A, B, C, v, p, b, c, vel_ex, pr_ex);
    }

};

/*
// StokesSolverBaseCL can be used as a preconditioner for the methods in solver.h
class StokesSolverAsPreCL
{
  private:
    StokesSolverBaseCL& solver_;
    mutable std::ostream* output_;

  public:
    StokesSolverAsPreCL( StokesSolverBaseCL& solver, int max_iter, std::ostream* output= 0)
        : solver_( solver), output_( output) {solver.SetMaxIter(max_iter);}

    void
    Apply(const BlockMatrixCL& M, VectorCL& x, const VectorCL& rhs) const {
        VectorCL v(M.num_cols( 0)), p(M.num_cols( 1)),
                 b(rhs[std::slice( 0, M.num_rows( 0), 1)]),
                 c(rhs[std::slice( M.num_rows( 0), M.num_rows( 1), 1)]);
        solver_.Solve(*M.GetBlock(0), *M.GetBlock(1), v, p, b, c);
        x[std::slice( 0, M.num_cols( 0), 1)]= v;
        x[std::slice( M.num_cols( 0), M.num_cols( 1), 1)]= p;
        if (output_ != 0)
            *output_<< "StokesSolverAsPreCL: iterations: " << solver_.GetIter()
                    << "\trelative residual: " << solver_.GetResid() << std::endl;
    }
    void
    Apply(const MLBlockMatrixCL& M, VectorCL& x, const VectorCL& rhs) const {
        VectorCL v(M.num_cols( 0)), p(M.num_cols( 1)),
                 b(rhs[std::slice( 0, M.num_rows( 0), 1)]),
                 c(rhs[std::slice( M.num_rows( 0), M.num_rows( 1), 1)]);
        solver_.Solve(*M.GetBlock(0), *M.GetBlock(1), v, p, b, c);
        x[std::slice( 0, M.num_cols( 0), 1)]= v;
        x[std::slice( M.num_cols( 0), M.num_cols( 1), 1)]= p;
        if (output_ != 0)
            *output_<< "StokesSolverAsPreCL: iterations: " << solver_.GetIter()
                    << "\trelative residual: " << solver_.GetResid() << std::endl;
    }
};*/

//=============================================================================
//  SchurComplMatrixCL
//=============================================================================
template<typename, typename, typename>
  class SchurComplMatrixCL;

template<typename T1, typename T2, typename T3>
  VectorCL operator*(const SchurComplMatrixCL<T1, T2, T3>&, const VectorCL&);

template<typename PoissonSolverT, typename MatT, typename ExVCL>
class SchurComplMatrixCL
/// This class can multiply B*A^(-1)*B * x.
{
  private:
    PoissonSolverT& solver_;
    const MatT& A_;
    const MatT& B_;
    //const MatT& C_;
    const ExVCL& ex_;

  public:
    SchurComplMatrixCL(PoissonSolverT& solver, const MatT& A, const MatT& B, const ExVCL& ex)
      : solver_(solver), A_( A), B_( B), ex_(ex) {}

    friend VectorCL
    operator*<>(const SchurComplMatrixCL<PoissonSolverT, MatT, ExVCL>&, const VectorCL&);

    const ExVCL& GetVelEx() const { return ex_; }
};

template<class PoissonSolverT, typename MatT, typename ExVCL>
  VectorCL operator*(const SchurComplMatrixCL<PoissonSolverT, MatT, ExVCL>& M, const VectorCL& v)
{
    VectorCL x(  M.A_.num_cols() );
    M.solver_.Solve( M.A_, x, transp_mul( M.B_, v), M.ex_);
//     IF_MASTER
//       std::cout << "> inner iterations: " << M.solver_.GetIter()
//                 << "\tresidual: " << M.solver_.GetResid()
//                 << std::endl;
    return M.B_*x;
}


// ***************************************************************************
/// \brief approximate Schur complement matrix class
// ***************************************************************************
template<typename, typename, typename>
  class ApproximateSchurComplMatrixCL;

template<typename T1, typename T2, typename T3>
  VectorCL operator*(const ApproximateSchurComplMatrixCL<T1, T2, T3>&, const VectorCL&);

template<typename APC, typename MatT, typename ExVCL>
class ApproximateSchurComplMatrixCL
/// Assume M is a approximate of A. Then instead of performing B*A^(-1)*B^T * v this class can
/// perform the multiplication B*M^(-1)*B^T * v
{
  private:
    const MatT& A_;
    APC& Apc_;
    const MatT& B_;
    const MatT& C_;
    const ExVCL& ex_;

  public:
    ApproximateSchurComplMatrixCL(const MatT& A, APC& Apc, const MatT& B, const MatT& C, const ExVCL& ex)
      : A_( A), Apc_( Apc), B_( B), C_( C ), ex_(ex) {}

    friend VectorCL
    operator*<>(const ApproximateSchurComplMatrixCL<APC, MatT, ExVCL>&, const VectorCL&);
    const ExVCL& GetVelEx() const { return ex_; }
    size_t Version() const  { return A_.Version(); } ///< Get modification version number
};

template<class APC, typename MatT, typename ExVCL>
  VectorCL operator*(const ApproximateSchurComplMatrixCL<APC, MatT, ExVCL>& M, const VectorCL& v)
{
    VectorCL x( 0.0, M.B_.num_cols());
    VectorCL r= transp_mul( M.B_, v);
    M.Apc_.Apply( M.A_, x, r, M.ex_);
    if (!M.Apc_.RetAcc())
        M.ex_.Accumulate(x);    

    VectorCL result( M.B_*x );
    result -= M.C_*v;
    return result;
}

/// \todo: this routine is not (yet) suited for non-empty (2,2)-block (e.g. pressure stabilization)
//-----------------------------------------------------------------------------
// UzawaPCG: The basic scheme is identical to PCG, but two extra recursions for
//     InexactUzawa are performed to avoid an evaluation of a preconditioner
//     there.
//     Also, the matrix is fixed: BQ_A^{-1}B^T and we assume x=zbar=zhat=0
//     upon entry to this function.
//     See "Fast iterative solvers for Stokes equation", Peters, Reichelt,
//     Reusken, Ch. 3.3 Remark 3.
//-----------------------------------------------------------------------------
template <typename APC, typename Mat, typename Vec, typename SPC>
bool
UzawaPCG(const APC& Apc, const Mat& A, const Mat& B,
    Vec& x, Vec& zbar, Vec& zhat, const Vec& b, const SPC& M,
    int& max_iter, double& tol)
{
    const size_t n= x.size();
    Vec p(n), z(n), q(n), d(n), e(n), r= b;
    Vec q1( B.num_cols()), q2( B.num_cols());
    double rho, rho_1= 0.0, resid= norm_sq( r);

    tol*= tol;
    if (resid<=tol)
    {
        tol= std::sqrt( resid);
        max_iter= 0;
        return true;
    }

    for (int i=1; i<=max_iter; ++i)
    {
        M.Apply( B/*dummy*/, z, r, DummyExchangeCL(), DummyExchangeCL());
        rho= dot( r, z);
        if (i == 1)
            p= z;
        else
            z_xpay(p, z, (rho/rho_1), p); // p= z + (rho/rho_1)*p;

        // q= A*p;
        q1= transp_mul( B, p);
        q2= 0.0;
        Apc.Apply( A, q2, q1, DummyExchangeCL());
        q= B*q2;

        const double alpha= rho/dot( p, q);
        axpy(alpha, p, x);                // x+= alpha*p;
        axpy(-alpha, q, r);               // r-= alpha*q;
        zbar+= alpha*q1;
        zhat+= alpha*q2;
        resid= norm_sq( r);
        if (resid<=tol)
        {
            tol= std::sqrt( resid);
            max_iter= i;
            return true;
        }
        rho_1= rho;
    }
    tol= std::sqrt(resid);
    return false;
}

//=============================================================================
// Inexact Uzawa-method from "Fast Iterative Solvers for Discrete Stokes
// Equations", Peters, Reichelt, Reusken, Chapter 3.3.
// The preconditioner Apc for A must be "good" (MG-like) to guarantee
// convergence.
// Due to theory (see paper), we should use 0.2 < innerred < 0.7.
//=============================================================================
template <typename Mat, typename Vec, typename PC1, typename PC2, typename ExVCL, typename ExPCL>
  bool InexactUzawa(const Mat& A, const Mat& B, const Mat& C, Vec& xu_acc, Vec& xp_acc, const Vec& f, const Vec& g,
                    const ExVCL& exV, const ExPCL& exP, PC1& Apc, PC2& Spc,
                    int& max_iter, double& tol, InexactUzawaApcMethodT apcmeth,
                    double innerred, int innermaxiter, bool rel, std::ostream* output)
/// \param[in]     A            Coefficient matrix for velocities
/// \param[in]     B            Coefficient matrix for coupling velocity and pressure
/// \param[in,out] xu_acc IN:   initial vector for velocity, OUT: velocity result
/// \param[in,out] xp_acc IN:   initial vector for pressure, OUT: pressure result
/// \param[in]     f            upper rhs
/// \param[in]     g            lower rhs
/// \param[in]     exV          Exchange class for velocity
/// \param[in]     exP          Exchange class for pressure
/// \param[in]     Apc          Preconditioner for A
/// \param[in]     Spc          Preconditioner for Schur-Complement
/// \param[in,out] max_iter     IN: maximal iterations, OUT: used iterations
/// \param[in,out] tol          IN: tolerance for the residual, OUT: residual
/// \param[in,out] usedinnerit  number of accumulated steps for inner solver
/// \param[in]     apcmeth      Characteristics of the preconditioner for the A-block
/// \param[in]     innerred     inner reduction
/// \param[in]     innermaxiter maximal inner solver steps
/// \param[in,out] output       Debug information
{
    if (Apc.NeedDiag())
        Apc.SetDiag(A, exV);
    if (Spc.NeedDiag())
        Spc.SetDiag(B, exP);
    VectorCL ru( f - A*xu_acc - transp_mul( B, xp_acc));
    VectorCL rp( g - B*xu_acc - C*xp_acc);
    VectorCL w( f.size());
    VectorCL z( g.size());
    VectorCL z2( g.size());
    VectorCL zbar( f.size());
    VectorCL zhat( f.size());
    VectorCL du( f.size());
    VectorCL c( g.size());
    ApproximateSchurComplMatrixCL<PC1, Mat, ExVCL>* asc = apcmeth == APC_SYM_LINEAR ? 0 :
            new ApproximateSchurComplMatrixCL<PC1, Mat, ExVCL>( A, Apc, B, C, exV);
    double innertol;
    int inneriter;
    int pr_iter_cumulative= 0;
    double res_u = exV.Norm_sq( ru, false);
    double res_p = exP.Norm_sq( rp, false);
    double resid0= std::sqrt( res_u + res_p);
    double resid= resid0;
    if (output)
        (*output) << "   o InexactUzawa 0: res " << resid
                  << ", res-impuls " << std::sqrt(res_u)
                  << ", res-mass " << std::sqrt(res_p)
                  << '\n';
    if (rel) tol*= resid0;
    if (resid <= tol) { // The fixed point iteration between levelset and Stokes
        tol= resid;     // equation uses this to determine convergence.
        max_iter= 0;
        if (asc!=0) delete asc;
        return true;
    }
    for (int k= 1; k <= max_iter; ++k) {
        w= 0.0;
        Apc.Apply( A, w, ru, exV);
        if (!Apc.RetAcc())
            exV.Accumulate(w);
        c= B*w - rp;            // w is accumulated
        z= 0.0;
        z2= 0.0;
        inneriter= innermaxiter;
        switch (apcmeth) {
             case APC_SYM_LINEAR:         // this case has not been implemented yet!
#ifdef _PAR
        throw DROPSErrCL("InexactUzawa: parallel implementation of \"APC_SYM_LINEAR\" is not available");
#endif
        throw DROPSErrCL("InexactUzawa: implementation of \"APC_SYM_LINEAR\" is not available for ghost penalty stabilization");
                 zbar= 0.0;
                 zhat= 0.0;
                 innertol= innerred*norm( c);
                 UzawaPCG( Apc, A, B, z, zbar, zhat, c, Spc, inneriter, innertol);
                 break;
            case APC_SYM:
                innertol= innerred*exP.Norm(c, false);
                PCG( *asc, z, c, exP, Spc, inneriter, innertol);
                break;
            default:
                std::cout << "WARNING: InexactUzawa: Unknown apcmeth; using GMRes.\n";
            // fall through
            case APC_OTHER:
                innertol= innerred; // GMRES can do relative tolerances.
#ifdef _PAR
                ModGMRES( *asc, z, c, exP, Spc, /*restart*/ inneriter,
                    inneriter, innertol, /*relative errors*/ true, /*use MGS*/false,
                    LeftPreconditioning);
#else
                GMRES( *asc, z, c, exP, Spc, /*restart*/ inneriter, inneriter, innertol,
                                /*relative errors*/ true, /*don't check 2-norm*/ false);
#endif
                break;
        }
        if (apcmeth != APC_SYM_LINEAR) {
            zbar= transp_mul( B, z);        // z is accumulated
            zhat= 0.0;
            Apc.Apply( A, zhat, zbar, exV);
            if (!Apc.RetAcc())
                exV.Accumulate(zhat);
        }
        pr_iter_cumulative+= inneriter;
        du= w - zhat;           // w and zhat are accumulated, so du is accumulated as well
        xu_acc+= du;
        xp_acc+= z;             // z is in every case accumulated, because it is a result of a solver
        ru-= A*du + zbar;
        rp= g - B*xu_acc - C*xp_acc;
        //rp-= B*du + C*z;
        res_u= exV.Norm_sq( ru, false);
        res_p = exP.Norm_sq( rp, false);
        resid= std::sqrt( res_u + res_p);
        if (output)
            (*output) << "   o InexactUzawa "<<k<<": res " << resid
                      << ", res-impuls " << std::sqrt(res_u)
                      << ", res-mass " << std::sqrt(res_p)
                      << ", inner_iter "<<inneriter
                      << '\n';

        if (resid <= tol) { // absolute errors
            tol= resid;
            max_iter= k;
            if (asc!=0) delete asc;
            return true;
        }
    }
    tol= resid;
    if (asc!=0) delete asc;
    std::cout << "===> Warning, InexactUzawa stopped without reaching tolerance!" << std::endl;
    return false;
}

/// \todo: this routine is not (yet) suited for non-empty (2,2)-block (e.g. pressure stabilization)
//-----------------------------------------------------------------------------
// UzawaPCG: The basic scheme is identical to PCG, but two extra recursions for
//     InexactUzawa are performed to avoid an evaluation of a preconditioner
//     there.
//     Also, the matrix is fixed: BQ_A^{-1}B^T and we assume x=zbar=zhat=0
//     upon entry to this function.
//     See "Fast iterative solvers for Stokes equation", Peters, Reichelt,
//     Reusken, Ch. 3.3 Remark 3.
//-----------------------------------------------------------------------------

template <typename Mat, typename Vec, typename PC1, typename PC2>
bool
UzawaCGEff(const Mat& A, const Mat& B, Vec& xu, Vec& xp, const Vec& f, const Vec& g,
    PC1& Apc, PC2& Spc,
    int& max_iter, double& tol)
{
    double err= std::sqrt( norm_sq( f - ( A*xu + transp_mul( B, xp))) + norm_sq( g - B*xu));
    const double err0= err;
    Vec rbaru( f - (A*xu + transp_mul(  B, xp)));
    Vec rbarp( g - B*xu);
    Vec ru( f.size());
    Apc.Apply( A, ru, rbaru);
    Vec rp( B*ru - rbarp);
    Vec a( f.size()), b( f.size()), s( f.size()), pu( f.size()), qu( f.size());
    Vec z( g.size()), pp( g.size()), qp( g.size()), t( g.size());
    double alpha= 0.0, initialbeta=0.0, beta= 0.0, beta0= 0.0, beta1= 0.0;
    for (int i= 0; i < max_iter; ++i) {
        z= 0.0;
        Spc.Apply( B, z, rp);
        a= A*ru;
        beta1= dot(a, ru) - dot(rbaru,ru) + dot(z,rp);
        if (i==0) initialbeta= beta1;
        if (beta1 <= 0.0) {throw DROPSErrCL( "UzawaCGEff: Matrix is not spd.\n");}
        // This is for fair comparisons of different solvers:
        err= std::sqrt( norm_sq( f - (A*xu + transp_mul( B, xp))) + norm_sq( g - B*xu));
        std::cout << "relative residual (2-norm): " << err/err0
                  << "\t(problem norm): " << std::sqrt( beta1/initialbeta) << '\n';
//        if (beta1/initialbeta <= tol*tol) {
//            tol= std::sqrt( beta1/initialbeta);
        if (err/err0 <= tol) {
            tol= err/err0;
            max_iter= i;
            return true;
        }
        if (i != 0) {
            beta= beta1/beta0;
            z_xpay( pu, ru, beta, pu); // pu= ru + beta*pu;
            z_xpay( pp, z, beta, pp);  // pp= z  + beta*pp;
            z_xpay( b, a, beta, b);    // b= a + beta*b;
        }
        else {
            pu= ru;
            pp= z;
            b= a; // A*pu;
        }
        qu= b + transp_mul( B, pp);
        qp= B*pu;
        s= 0.0;
        Apc.Apply( A, s, qu);
        z_xpay( t, B*s, -1.0, qp); // t= B*s - qp;
        alpha= beta1/(dot( s, b) - dot(qu, pu) + dot(t, pp));
        axpy( alpha, pu, xu); // xu+= alpha*pu;
        axpy( alpha, pp, xp); // xp+= alpha*pp;
        axpy( -alpha, s, ru); // ru-= alpha*s;
        axpy( -alpha, t, rp); // rp-= alpha*t;
        axpy( -alpha, qu, rbaru);// rbaru-= alpha*qu;
        // rbarp-= alpha*qp; rbarp is not needed in this algo.
        beta0= beta1;
    }
    tol= std::sqrt( beta1);
    return false;
}

/// \todo: this routine is not (yet) suited for non-empty (2,2)-block (e.g. pressure stabilization)
//-----------------------------------------------------------------------------
// CG: The return value indicates convergence within max_iter (input)
// iterations (true), or no convergence within max_iter iterations (false).
// max_iter - number of iterations performed before tolerance was reached
//      tol - residual after the final iteration
//-----------------------------------------------------------------------------

template <typename Mat, typename Vec, typename PC1, typename PC2>
bool
UzawaCG(const Mat& A, const Mat& B, Vec& u, Vec& p, const Vec& b, const Vec& c,
        PC1& Apc, PC2& Spc, int& max_iter, double& tol)
{
    double err= std::sqrt( norm_sq( b - (A*u + transp_mul( B, p))) + norm_sq( c - B*u));
    const double err0= err;
    Vec ru( b - ( A*u + transp_mul(  B, p)));
    Vec rp( c - B*u);
    Vec s1( b.size()); // This is r2u...
    Apc.Apply( A, s1, ru);
    Vec s2( B*s1 - rp);
    Vec r2p( c.size());
    Spc.Apply( B, r2p, s2);
    double rho0= dot( s1, VectorCL( A*s1 - ru)) + dot( r2p, s2);
    const double initialrho= rho0;
//    std::cout << "UzawaCG: rho: " << rho0 << '\n';
    if (rho0<=0.0) throw DROPSErrCL("UzawaCG: Matrix is not spd.\n");
//    tol*= tol*rho0*rho0; // For now, measure the relative error.
//    if (rho0<=tol) {
//        tol= std::sqrt( rho0);
//        max_iter= 0;
//        return true;
//    }
    Vec pu= s1; // s1 is r2u.
    Vec pp= r2p;
    Vec qu( A*pu + transp_mul( B, pp));
    Vec qp= B*pu;
    double rho1= 0.0;
    Vec t1( b.size());
    Vec t2( c.size());
    for (int i= 1; i<=max_iter; ++i) {
        Apc.Apply( A, t1, qu);
        z_xpay( t2, B*t1, -1.0, qp); // t2= B*t1 - qp;
        const double alpha= rho0/( dot(pu, VectorCL( A*t1 - qu)) + dot( pp, t2));
        axpy(alpha, pu, u);  // u+= alpha*pu;
        axpy(alpha, pp, p);  // p+= alpha*pp;
        axpy( -alpha, qu, ru);
        axpy( -alpha, qp, rp);
        s1= 0.0;
        Apc.Apply( A, s1, ru);
        z_xpay( s2, B*s1, -1.0, rp); // s2= B*s1 - rp;
        //axpy( -alpha, t1, s1); // kann die beiden oberen Zeilen ersetzen,
        //axpy( -alpha, t2, s2); // Algorithmus wird schneller, aber Matrix bleibt nicht spd
        r2p= 0.0;
        Spc.Apply( B, r2p, s2);
        rho1= dot( s1, VectorCL( A*s1 - ru)) + dot( r2p, s2);
//        std::cout << "UzawaCG: rho: " << rho1 << '\n';
        if (rho1<=0.0) throw DROPSErrCL("UzawaCG: Matrix is not spd.\n");
        // This is for fair comparisons of different solvers:
        err= std::sqrt( norm_sq( b - (A*u + transp_mul( B, p))) + norm_sq( c - B*u));
        std::cout << "relative residual (2-norm): " << err/err0
                << "\t(problem norm): " << std::sqrt( rho1/initialrho)<< '\n';
//        if (rho1 <= tol) {
//            tol= std::sqrt( rho1);
        if (err <= tol*err0) {
            tol= err/err0;
            max_iter= i;
            return true;
        }
        const double beta= rho1/rho0;
        z_xpay( pu, s1, beta, pu); // pu= s1 + beta*pu; // s1 is r2u.
        z_xpay( pp, r2p, beta, pp); // pp= r2p + beta*pp;
        qu= (A*pu) + transp_mul( B, pp);
        qp= (B*pu);
        rho0= rho1;
        t1= 0.0;
    }
    tol= std::sqrt( rho1);
    return false;
}

/// \todo: this routine is not (yet) suited for non-empty (2,2)-block (e.g. pressure stabilization)
template <typename PC1, typename PC2>
class UzawaCGSolverEffCL : public StokesSolverBaseCL
{
  private:
    PC1& Apc_;
    PC2& Spc_;

  public:
    UzawaCGSolverEffCL (PC1& Apc, PC2& Spc, int maxiter, double tol)
        : StokesSolverBaseCL(maxiter, tol), Apc_( Apc), Spc_( Spc) {}

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c) {
        res_=  tol_;
        iter_= maxiter_;
        UzawaCGEff( A, B, v, p, b, c, Apc_, Spc_, iter_, res_);
    }
    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c) {
        res_=  tol_;
        iter_= maxiter_;
        UzawaCGEff( A, B, v, p, b, c, Apc_, Spc_, iter_, res_);
    }
};

/// \todo: this routine is not (yet) suited for non-empty (2,2)-block (e.g. pressure stabilization)
template <typename PC1, typename PC2>
class UzawaCGSolverCL : public StokesSolverBaseCL
{
  private:
    PC1& Apc_;
    PC2& Spc_;

  public:
    UzawaCGSolverCL (PC1& Apc, PC2& Spc, int maxiter, double tol)
        : StokesSolverBaseCL(maxiter, tol), Apc_( Apc), Spc_( Spc) {}

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c) {
        res_=  tol_;
        iter_= maxiter_;
        UzawaCG( A, B, v, p, b, c, Apc_, Spc_, iter_, res_);
    }
    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c) {
        res_=  tol_;
        iter_= maxiter_;
        UzawaCG( A, B, v, p, b, c, Apc_, Spc_, iter_, res_);
    }
};

//important for ScaledMGPreCL
template <class ProlongationT>
double
EigenValueMaxMG(const MLMatrixCL& A, const ProlongationT& P, VectorCL& x, int iter, double  tol= 1e-3)
{
    MLMatrixCL::const_iterator finest= A.GetFinestIter();
    typename ProlongationT::const_iterator finestP= P.GetFinestIter();
    Uint   sm   =  1; // how many smoothing steps?
    int    lvl  = -1; // how many levels? (-1=all)
    double omega= 1.; // relaxation parameter for smoother
    double l= -1.0;
    double l_old= -1.0;
    VectorCL tmp( x.size());
    VectorCL z( x.size());
//    JORsmoothCL smoother( omega); // Jacobi
//    GSsmoothCL smoother( omega); // Gauss-Seidel
//    SGSsmoothCL smoother( omega); // symmetric Gauss-Seidel
//    SORsmoothCL smoother( omega); // Gauss-Seidel with over-relaxation
    SSORsmoothCL smoother( omega); // symmetric Gauss-Seidel with over-relaxation
//    CGSolverCL  solver( 200, tol); //CG-Verfahren
    SSORPcCL directpc; PCGSolverCL<SSORPcCL> solver( directpc, 200, 1e-15);
    x/= norm( x);
    std::cout << "EigenValueMaxMG:\n";
    for (int i= 0; i<iter; ++i) {
        tmp= 0.0;
        MGM( A.begin(), finest, finestP, tmp, A*x, smoother, sm, solver, lvl, -1);
        z= x - tmp;
        l= dot( x, z);
        std::cout << "iteration: " << i  << "\tlambda: " << l << "\trelative_change= : " << (i==0 ? -1 : std::fabs( (l-l_old)/l_old)) << '\n';
        if (i > 0 && std::fabs( (l-l_old)/l_old) < tol) break;
        l_old= l;
        x= z/norm( z);
    }
    std::cout << "maximal value for lambda: " << l << '\n';
    return l;
}

//scaled MG as preconditioner
template <class ProlongationT= MLMatrixCL>
class ScaledMGPreCL
{
  private:
    const ProlongationT& P_;
    Uint iter_;
    double s_;
    mutable SSORPcCL directpc_;
    mutable PCGSolverCL<SSORPcCL> solver_;
    Uint sm_; // how many smoothing steps?
    int lvl_; // how many levels? (-1=all)
    double omega_; // relaxation parameter for smoother
    mutable SSORsmoothCL smoother_;  // Symmetric-Gauss-Seidel with over-relaxation

  public:
    ScaledMGPreCL( const ProlongationT& P, Uint iter, double s= 1.0)
        : P_(P), iter_( iter), s_( s), solver_( directpc_, 200, 1e-12), sm_( 1),
         lvl_( -1), omega_( 1.0), smoother_( omega_)
    {}

    inline void
    Apply( const MLMatrixCL& A, VectorCL& x, const VectorCL& r) const;
    inline void
    Apply( const MatrixCL&, VectorCL&, const VectorCL&) const
    {
        throw DROPSErrCL( "ScaledMGPreCL::Apply: need multilevel data structure\n");
    }
};

template <class ProlongationT>
inline void
ScaledMGPreCL<ProlongationT>::Apply( const MLMatrixCL& A, VectorCL& x, const VectorCL& r) const
{
    x= 0.0;
    MLMatrixCL::const_iterator finest = A.GetFinestIter();
    typename ProlongationT::const_iterator finestP= P_.GetFinestIter();
    const double oldres= norm(r - A*x);
    for (Uint i= 0; i < iter_; ++i) {
        VectorCL r2( r - A*x);
        VectorCL dx( x.size());
        MGM( A.begin(), finest, finestP, dx, r2, smoother_, sm_, solver_, lvl_, -1);
        x+= dx*s_;
    }
    const double res= norm(r - A*x);
    std::cout << "ScaledMGPreCL: it: " << iter_ << "\treduction: " << (oldres==0.0 ? res : res/oldres) << '\n';
}

template<class SmootherT, class PVelT= MLMatrixCL, class PPrT= MLMatrixCL>
class StokesMGSolverCL: public StokesSolverBaseCL
{
  private:
    PVelT PVel_;
    PPrT  PPr_;
    const MLMatrixCL    &prM_;
    const SmootherT&    smoother_;
    StokesSolverBaseCL& directSolver_;
    Uint  smoothSteps_;
    int   usedLevels_;
    MLMatrixCL BT_;
    size_t BVersion_;

    void UpdateBT( const MLMatrixCL& B)
    {
        BT_.resize( B.size());
        MLMatrixCL::iterator BT = BT_.begin();
        for ( MLMatrixCL::const_iterator it = B.begin(); it != B.end(); ++it)
        {
            transpose (*it, *BT);
            ++BT;
        }
        BVersion_ = B.Version();
    }

  public:
    StokesMGSolverCL( const MLMatrixCL& prM, const SmootherT& smoother, StokesSolverBaseCL& ds,
                      Uint iter_vel, double tol, bool rel= false, Uint sm = 2, int lvl = -1)
      : StokesSolverBaseCL(iter_vel, tol, rel), prM_(prM), smoother_(smoother), directSolver_(ds),
        smoothSteps_(sm), usedLevels_(lvl), BVersion_( 0) {}
    ~StokesMGSolverCL() {}

    PVelT* GetPVel() { return &PVel_;}
    PPrT*  GetPPr()  { return &PPr_; }
#ifdef _PAR
    void
    Solve(const MatrixCL& /*A*/, const MatrixCL& /*B*/, const MatrixCL&, VectorCL&, VectorCL&, const VectorCL&, const VectorCL&, const ExchangeCL&, const ExchangeCL&)
    {
        throw DROPSErrCL( "StokesMGSolverCL::Solve: need multilevel data structure\n");
    }
    void
    Solve(const MLMatrixCL& /*A*/, const MLMatrixCL& /*B*/, const MLMatrixCL&, VectorCL&, VectorCL&, const VectorCL&, const VectorCL&, const ExchangeCL&, const ExchangeCL&)
    {
        throw DROPSErrCL( "StokesMGSolverCL::Solve: not implemented in parallel\n");
    }
#endif
    void
    Solve(const MatrixCL& /*A*/, const MatrixCL& /*B*/, const MatrixCL&, VectorCL&, VectorCL&, const VectorCL&, const VectorCL&, const DummyExchangeCL&, const DummyExchangeCL&)
    {
        throw DROPSErrCL( "StokesMGSolverCL::Solve: need multilevel data structure\n");
    }

    void
    Solve(const MLMatrixCL& A, const MLMatrixCL& B, const MLMatrixCL& C, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c, const DummyExchangeCL&, const DummyExchangeCL&) {
// define MG parameters for the first diagonal blockS
        if (B.Version() != BVersion_) UpdateBT( B);
        int nit=maxiter_;
        double actualtol = 1;
        Uint   wc   = 1;   // how many W-cycle steps? (1=V-cycle)

// define initial approximation
        MLMatrixCL::const_iterator A_end   ( A.GetFinestIter());
        MLMatrixCL::const_iterator B_end   ( B.GetFinestIter());
        MLMatrixCL::const_iterator BT_end  ( BT_.GetFinestIter());
        MLMatrixCL::const_iterator C_end   ( C.GetFinestIter());
        MLMatrixCL::const_iterator prM_end ( prM_.GetFinestIter());
        typename PVelT::const_iterator PVel( PVel_.GetFinestIter());
        typename PPrT::const_iterator  PPr ( PPr_.GetFinestIter());

        const double runorm0= norm_sq( A * v + transp_mul( B, p) - b);
        const double rpnorm0= norm_sq( B * v - c);
        const double resid0= std::sqrt(runorm0+rpnorm0);
        double resid;
        for (int j=0; j<maxiter_; ++j)
        {
            StokesMGM( A.begin(), A_end, B_end, BT_end, C_end, prM_end, PVel, PPr, v, p, b, c, smoother_, smoothSteps_, wc, directSolver_, usedLevels_, -1);
            const double runorm= norm_sq( A * v + transp_mul(B, p ) - b);
            const double rpnorm= norm_sq( B * v - c);
            resid= std::sqrt(runorm+rpnorm);
            if (rel_)
                actualtol= resid/resid0;
            else
                actualtol= resid;
            std::cout << "P2P1:StokesMGSolverCL: residual = " << actualtol << std::endl;
            if (actualtol<=tol_)
            {
                nit= j+1;
                break;
            }
        }
        std::cout << "StokesMGM: actual residual = " << actualtol
                  << "  after " << nit << " iterations " << std::endl;
        res_ = actualtol;
        iter_= nit;
   }
};

///\brief Preconditioner for the (Navier-) Stokes operator
///
/// Vertex-based multiplicative Schwartz-method. The implementation is based on the Vanka-smoother class for
/// Stokes-multigrid in Drops. The effect of Apply is one smoother-iteration with 0 as starting value.
class VankaPreCL
{
  private:
    PVankaSmootherCL smoother_;
    mutable MatrixCL BT_;
    mutable size_t BVersion_;

    void UpdateBT (const MatrixCL& B) const {
        transpose (B, BT_);
        BVersion_ = B.Version();
    }

  public:
    VankaPreCL (const MLIdxDescCL* idx= 0) : smoother_( 0, 1., idx), BVersion_( 0) {}

    void Setidx (const MLIdxDescCL* idx) { smoother_.Setidx( idx); }

    template<typename ExT>
    void
    Apply (const BlockMatrixCL& M, VectorCL& x, const VectorCL& rhs, const ExT& ex) const {
        VectorCL v(M.num_cols( 0)), p(M.num_cols( 1)),
                 b(rhs[std::slice( 0, M.num_rows( 0), 1)]),
                 c(rhs[std::slice( M.num_rows( 0), M.num_rows( 1), 1)]);
        if (M.GetBlock( 1)->Version() != BVersion_) UpdateBT( *M.GetBlock( 1));
        smoother_.Apply( *M.GetBlock( 0), *M.GetBlock( 1), BT_, /*dummy*/ *M.GetBlock( 0), v, p, b, c, ex.GetEx(0), ex.GetEx(1));
        x[std::slice( 0, M.num_cols( 0), 1)]= v;
        x[std::slice( M.num_cols( 0), M.num_cols( 1), 1)]= p;
    }
    template<typename ExT>
    void
    Apply (const MLBlockMatrixCL& M, VectorCL& x, const VectorCL& rhs, const ExT& ex) const {
        BlockMatrixCL MM( &M.GetBlock( 0)->GetFinest(), MUL, &M.GetBlock( 1)->GetFinest(), TRANSP_MUL,
                          &M.GetBlock( 1)->GetFinest(), MUL);
        this->Apply( MM, x, rhs, ex);
    }
    /// \name Parallel preconditioner setup ...
    //@{
    bool NeedDiag() const { return false; }
    bool RetAcc()   const { return true; }
    template<typename Mat, typename ExT>
    void SetDiag( const Mat&, const ExT&) {} // just for consistency
    void SetDiag( const VectorCL&)        {} // just for consistency
    //@}
};

///\brief Preconditioner for the Schur-complement of the (Navier-) Stokes operator
///
/// Vertex-based multiplicative Schwartz-method. The implementation is based on the Vanka-smoother class for
/// Stokes-multigrid in Drops. The effect of Apply is one smoother-iteration with 0 as starting value. The
/// right-hand side for the momentum equations is set to 0. Thus, the action of the Schur-complement is approximated.
class VankaSchurPreCL: public SchurPreBaseCL
{
  private:
    PVankaSmootherCL smoother_;
    mutable MatrixCL BT_;
    mutable size_t BVersion_;
    mutable BlockMatrixCL M;

    void UpdateBT (const MatrixCL& B) const {
        transpose (B, BT_);
        BVersion_ = B.Version();
    }

  public:
    VankaSchurPreCL (const MLIdxDescCL* idx= 0)
        : SchurPreBaseCL( 0, 0), smoother_( 0, 1., idx), BVersion_( 0), M( 0, MUL, 0, TRANSP_MUL, 0, MUL) {}

    void Setidx (const MLIdxDescCL* idx)   { smoother_.Setidx( idx); }
    void SetAB   (const MatrixCL* A, const MatrixCL* B) { M.SetBlock( 0, A); M.SetBlock( 1, B); M.SetBlock( 2, B); }

    template<typename Mat,  typename Vec, typename ExT>
    void
    Apply (const Mat& /*S*/, Vec& x, const Vec& rhs, const ExT& vel_ex, const ExT& p_ex) const {
        Vec v(M.num_cols( 0)), p(M.num_cols( 1)), b( M.num_rows( 0));
        if (M.GetBlock( 1)->Version() != BVersion_) UpdateBT( *M.GetBlock( 1));
        smoother_.Apply( *M.GetBlock( 0), *M.GetBlock( 1), BT_, /*dummy*/ *M.GetBlock( 0), v, p, b, rhs, vel_ex, p_ex);
        x= p;
    }
#ifdef _PAR
    void Apply(const MatrixCL& A, VectorCL& x, const VectorCL& b, const ExchangeCL& vel_ex, const ExchangeCL& p_ex) const
    { Apply<>(A, x, b, vel_ex, p_ex);}
    void Apply(const MLMatrixCL& A, VectorCL& x, const VectorCL& b, const ExchangeCL& vel_ex, const ExchangeCL& p_ex) const
    { Apply<>(A, x, b, vel_ex, p_ex);}
#endif
    void Apply(const MatrixCL& A,   VectorCL& x, const VectorCL& b, const DummyExchangeCL& vel_ex, const DummyExchangeCL& p_ex) const
    { Apply<>(A, x, b, vel_ex, p_ex);}
    void Apply(const MLMatrixCL& A, VectorCL& x, const VectorCL& b, const DummyExchangeCL& vel_ex, const DummyExchangeCL& p_ex) const
    { Apply<>(A, x, b, vel_ex, p_ex);}

    /// \name Parallel preconditioner setup ...
    //@{
    bool NeedDiag() const { return false; }
    bool RetAcc()   const { return true; }
    template<typename Mat, typename ExT>
    void SetDiag( const Mat&, const ExT&) {} // just for consistency
    void SetDiag( const VectorCL&)        {} // just for consistency
    //@}
};

} // end of namespace DROPS

#endif
