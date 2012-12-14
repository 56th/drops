/// \file oseenprecond.h
/// \brief preconditioners for the oseen problem and the schur complement
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Helmut Jarausch, Volker Reichelt; SC RWTH Aachen:

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

#include "num/precond.h"
#include "num/MGsolver.h"
#include "num/spblockmat.h"


#ifndef OSEENPRECOND_H_
#define OSEENPRECOND_H_

namespace DROPS
{

/// fwd decl from num/stokessolver.h
template<typename, typename>
class ApproximateSchurComplMatrixCL;

#ifdef _PAR
/// fwd decl from num/parstokessolver.h
template<typename, typename, typename>
class ParApproximateSchurComplMatrixCL;
#endif

/// base class for Schur complement preconditioners
class SchurPreBaseCL: public PreBaseCL
{
  protected:
    double kA_,   ///< scaling factor for pressure stiffness matrix or equivalent
           kM_;   ///< scaling factor for pressure mass matrix

  public:
    SchurPreBaseCL( double kA, double kM, std::ostream* output=0) : PreBaseCL( output), kA_( kA), kM_( kM) {}
    virtual ~SchurPreBaseCL() {}
    void SetWeights( double kA, double kM) { kA_ = kA; kM_ = kM; }

    virtual void Apply( const MatrixCL& A,   VectorCL& x, const VectorCL& b) const = 0;
    virtual void Apply( const MLMatrixCL& A, VectorCL& x, const VectorCL& b) const = 0;
    // dummy Apply routine for inexact Uzawa. The matrix parameter in the Apply is not used, anyway.
    template<typename PcT, typename Mat>
    void Apply( const ApproximateSchurComplMatrixCL<PcT, Mat>&, VectorCL& x, const VectorCL& b) const { const MatrixCL dummy; Apply( dummy, x, b); }
};

//**************************************************************************
// Preconditioner for the instationary Stokes-equations.
// cf. "Iterative Techniques For Time Dependent Stokes Problems",
//     James H. Bramble, Joseph E. Pasciak, January 1994
//
// A Poisson-problem with natural boundary-conditions for the pressure is
// solved via 1 SSOR-step, a problem with the mass-matrix aswell.
// The constants kA_, kM_ have to be chosen according to h and dt, see Theorem 4.1
// of the above paper.
// kA_ = theta/Re and kM_ = 1/dt will do a good job,
// where Re is proportional to the ratio density/viscosity.
//
// A_ is the pressure-Poisson-Matrix for natural boundary-conditions, M_ the
// pressure-mass-matrix.
//**************************************************************************
class ISPreCL : public SchurPreBaseCL
{
  private:
    MatrixCL& A_;
    MatrixCL& M_;
    SSORPcCL  ssor_;

  public:
    ISPreCL( MatrixCL& A_pr, MatrixCL& M_pr,
        double kA= 0., double kM= 1., double om= 1.)
        : SchurPreBaseCL( kA, kM), A_( A_pr), M_( M_pr), ssor_( om)  {}
    ISPreCL( MLMatrixCL& A_pr, MLMatrixCL& M_pr,
             double kA= 0., double kM= 1., double om= 1.)
    : SchurPreBaseCL( kA, kM), A_( A_pr.GetFinest()), M_( M_pr.GetFinest()), ssor_( om)  {}

    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec& p, const Vec& c) const;
    void Apply(const MatrixCL& A,   VectorCL& x, const VectorCL& b) const { Apply<>( A, x, b); }
    void Apply(const MLMatrixCL& A, VectorCL& x, const VectorCL& b) const { Apply<>( A, x, b); }
};


//**************************************************************************
// Preconditioner for the instationary Stokes-equations.
// Confer ISPreCL for details. This preconditioner uses a few CG-iterations
// to solve the linear systems.
// It is well suited for InexactUzawa-Solvers.
//**************************************************************************
#ifndef _PAR
template <class SolverT>
class ISNonlinearPreCL : public SchurPreBaseCL
{
  private:
    MatrixCL&  A_;
    MatrixCL&  M_;
    SolverT&   solver_;

  public:
    ISNonlinearPreCL(SolverT& solver, MatrixCL& A_pr, MatrixCL& M_pr,
        double kA= 0., double kM= 1.)
        : SchurPreBaseCL( kA, kM), A_( A_pr), M_( M_pr),
          solver_( solver)  {}

    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec& p, const Vec& c) const;
    void Apply(const MatrixCL& A,   VectorCL& x, const VectorCL& b) const { Apply<>( A, x, b); }
    void Apply(const MLMatrixCL& A, VectorCL& x, const VectorCL& b) const { Apply<>( A, x, b); }
};
#else
template <typename ASolverT, typename MSolverT>
class ISNonlinearPreCL : public SchurPreBaseCL
{
  private:
    MatrixCL&  A_;
    MatrixCL&  M_;
    ASolverT& Asolver_;
    MSolverT& Msolver_;
    mutable typename ASolverT::PrecondT PcA_;
    mutable typename MSolverT::PrecondT PcM_;

  public:
    ISNonlinearPreCL(ASolverT& Asolver, MSolverT& Msolver, MatrixCL& A_pr, MatrixCL& M_pr,
        double kA= 0., double kM= 1.)
        : SchurPreBaseCL( kA, kM), A_( A_pr), M_( M_pr),
          Asolver_(Asolver), Msolver_(Msolver), PcA_(Asolver_.GetPC()), PcM_(Msolver_.GetPC())  {}

    /// \brief Apply preconditioner
    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec& p, const Vec& c) const;
    void Apply(const MatrixCL& A,   VectorCL& x, const VectorCL& b) const { Apply<>( A, x, b); }
    void Apply(const MLMatrixCL& A, VectorCL& x, const VectorCL& b) const { Apply<>( A, x, b); }

    /// \brief preconditionied vector is accumulated after "Apply"
    inline bool RetAcc() const {
        return true;
    }

    /// \brief Check if preconditioners needs diagonal of the matrices
    inline bool NeedDiag() const {
        return Asolver_.GetPC().NeedDiag() && Msolver_.GetPC().NeedDiag();
    }

    /// \brief Set diagonal of the preconditioner of the solvers (the matrices are known by this class)
    template <typename Mat>
    void SetDiag(const Mat&)
    {
        PcA_.SetDiag(A_);
        PcM_.SetDiag(M_);
    }


};
#endif

//**************************************************************************
// Preconditioner for the instationary Stokes-equations.
// cf. "Iterative Techniques For Time Dependent Stokes Problems",
//     James H. Bramble, Joseph E. Pasciak, January 1994
//
// Confer ISPreCL for details regarding preconditioning of S. This
// preconditioner uses multigrid-solvers.
//
//**************************************************************************
template< class ProlongationT>
class ISMGPreCL : public SchurPreBaseCL
{
  private:
    const Uint sm; // how many smoothing steps?
    const int lvl; // how many levels? (-1=all)
    const double omega; // relaxation parameter for smoother
    SSORsmoothCL smoother;  // Symmetric-Gauss-Seidel with over-relaxation
    SSORPcCL directpc;
    mutable PCGSolverCL<SSORPcCL> solver;

    DROPS::MLMatrixCL& Apr_;
    DROPS::MLMatrixCL& Mpr_;
    ProlongationT      P_;
    DROPS::Uint iter_prA_;
    DROPS::Uint iter_prM_;
    mutable std::vector<DROPS::VectorCL> ones_;

    void MaybeInitOnes() const;

  public:
    ISMGPreCL(DROPS::MLMatrixCL& A_pr, DROPS::MLMatrixCL& M_pr,
                    double kA, double kM, DROPS::Uint iter_prA=1,
                    DROPS::Uint iter_prM = 1)
        : SchurPreBaseCL( kA, kM), sm( 1), lvl( -1), omega( 1.0), smoother( omega), solver( directpc, 200, 1e-12),
          Apr_( A_pr), Mpr_( M_pr), iter_prA_( iter_prA), iter_prM_( iter_prM), ones_(0)
    {}

    template <typename Mat, typename Vec>
    void
    Apply(const Mat& /*A*/, Vec& p, const Vec& c) const;
    void Apply(const MatrixCL& A,   VectorCL& x, const VectorCL& b) const { Apply<>( A, x, b); }
    void Apply(const MLMatrixCL& A, VectorCL& x, const VectorCL& b) const { Apply<>( A, x, b); }

    ProlongationT* GetProlongation() { return &P_; }
};


//**************************************************************************
// Preconditioner for the instationary Stokes-equations.
// It uses BB^T instead of a Laplacian on the pressure space and is
// therefore suited for P1X-elements;
// cf. "Uniform Preconditioners for a Parameter Dependent Saddle Point
// Problem with Application to Generalized Stokes Interface Equations",
// Olshanskii, Peters, Reusken, 2005
//**************************************************************************
class ISBBTPreCL : public SchurPreBaseCL
{
  private:
    const MatrixCL*  B_;
    mutable MatrixCL*  Bs_;                                     ///< scaled Matrix B
    mutable size_t Bversion_;
    const MatrixCL*  M_, *Mvel_;

    double     tolA_, tolM_;                                    ///< tolerances of the solvers
    mutable VectorCL Dprsqrtinv_;                               ///< diag(M)^{-1/2}
#ifndef _PAR
    typedef NEGSPcCL SPcT_;
    SPcT_            spc_;
    JACPcCL          jacpc_;
    mutable PCGNESolverCL<SPcT_> solver_;
    mutable PCGSolverCL<JACPcCL> solver2_;
#else
    mutable CompositeMatrixCL BBT_;
    //typedef ParJacNEG0CL    PCSolver1T;                         ///< type of the preconditioner for solver 1
    typedef ParNEGSPcCL    PCSolver1T;                         ///< type of the preconditioner for solver 1
    typedef ParJac0CL       PCSolver2T;                         ///< type of the preconditioner for solver 2
    PCSolver1T PCsolver1_;
    PCSolver2T PCsolver2_;
    mutable ParPCGNESolverCL<PCSolver1T> solver_;                 ///< solver for BB^T
    mutable ParPCGSolverCL<PCSolver2T> solver2_;                ///< solver for M

    const IdxDescCL* vel_idx_;                                  ///< Accessing ExchangeCL for velocity
#endif
    const IdxDescCL* pr_idx_;                                   ///< Accessing ExchangeCL for pressure; also used to determine, how to represent the kernel of BB^T in case of pure Dirichlet-BCs.
    double regularize_;                                         ///< If regularize_==0. no regularization is performed. Otherwise, a column is attached to Bs.
    void Update () const;                                       ///< Updating the diagonal matrices D and Dprsqrtinv

  public:
#ifndef _PAR
    ISBBTPreCL (const MatrixCL* B, const MatrixCL* M_pr, const MatrixCL* Mvel,
        const IdxDescCL& pr_idx,
        double kA= 0., double kM= 1., double tolA= 1e-2, double tolM= 1e-2, double regularize= 0.)
        : SchurPreBaseCL( kA, kM), B_( B), Bs_( 0), Bversion_( 0),
          M_( M_pr), Mvel_( Mvel), tolA_(tolA), tolM_(tolM),
          solver_( spc_, 500, tolA_, /*relative*/ true),
          solver2_( jacpc_, 500, tolM_, /*relative*/ true),
          pr_idx_( &pr_idx), regularize_( regularize) {}

    ISBBTPreCL (const ISBBTPreCL& pc)
        : SchurPreBaseCL( pc.kA_, pc.kM_), B_( pc.B_), Bs_( pc.Bs_ == 0 ? 0 : new MatrixCL( *pc.Bs_)),
          Bversion_( pc.Bversion_),
          M_( pc.M_), Mvel_( pc.Mvel_),
          tolA_(pc.tolA_), tolM_(pc.tolM_),
          Dprsqrtinv_( pc.Dprsqrtinv_),
          spc_( pc.spc_),
          solver_( spc_, 500, tolA_, /*relative*/ true),
          solver2_( jacpc_, 500, tolM_, /*relative*/ true),
          pr_idx_( pc.pr_idx_), regularize_( pc.regularize_) {}
#else
    ISBBTPreCL (const MatrixCL* B, const MatrixCL* M_pr, const MatrixCL* Mvel,
        const IdxDescCL& pr_idx, const IdxDescCL& vel_idx,
        double kA= 0., double kM= 1., double tolA= 1e-2, double tolM= 1e-2, double regularize= 0.)
        : SchurPreBaseCL( kA, kM), B_( B), Bs_( 0), Bversion_( 0),
          M_( M_pr), Mvel_( Mvel), tolA_(tolA), tolM_(tolM),
          BBT_( 0, TRANSP_MUL, 0, MUL, vel_idx, pr_idx),
          PCsolver1_( pr_idx, vel_idx), PCsolver2_(pr_idx),
          solver_( 800, tolA_, vel_idx, pr_idx, PCsolver1_, /*relative*/ true),
          solver2_( 500, tolM_, pr_idx, PCsolver2_, /*relative*/ true),
          vel_idx_( &vel_idx), pr_idx_( &pr_idx), regularize_( regularize) {}
    ISBBTPreCL (const ISBBTPreCL& pc)
        : SchurPreBaseCL( pc.kA_, pc.kM_), B_( pc.B_), Bs_( pc.Bs_ == 0 ? 0 : new MatrixCL( *pc.Bs_)),
          Bversion_( pc.Bversion_),
          M_( pc.M_), Mvel_( pc.Mvel_),
          tolA_(pc.tolA_), tolM_(pc.tolM_),
          Dprsqrtinv_( pc.Dprsqrtinv_),
          BBT_( Bs_, TRANSP_MUL, Bs_, MUL, *pc.vel_idx_, *pc.pr_idx_),
          PCsolver1_( *pc.pr_idx_, *pc.vel_idx_), PCsolver2_( *pc.pr_idx_),
          solver_( 800, tolA_, *pc.vel_idx_, *pc.pr_idx_, PCsolver1_, /*relative*/ true),
          solver2_( 500, tolM_, *pc.pr_idx_, PCsolver2_, /*relative*/ true),
          vel_idx_( pc.vel_idx_), pr_idx_( pc.pr_idx_), regularize_( pc.regularize_){}

    /// \name Parallel preconditioner setup ...
    //@{
    bool NeedDiag() const { return false; }
    void SetDiag(const VectorCL&) {}        // just for consistency
    template<typename Mat>
    void SetDiag(const Mat&) {}             // just for consistency
    bool RetAcc()   const { return true; }
    //@}
#endif

    ISBBTPreCL& operator= (const ISBBTPreCL&) {
        throw DROPSErrCL( "ISBBTPreCL::operator= is not permitted.\n");
    }

    ~ISBBTPreCL () { delete Bs_; }

    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec& p, const Vec& c) const;
    void Apply(const MatrixCL& A,   VectorCL& x, const VectorCL& b) const { Apply<>( A, x, b); }
    void Apply(const MLMatrixCL& A, VectorCL& x, const VectorCL& b) const { Apply<>( A, x, b); }

    void SetMatrices (const MatrixCL* B, const MatrixCL* Mvel, const MatrixCL* M, const IdxDescCL* pr_idx) {
        B_= B;
        Mvel_= Mvel;
        M_= M;
        pr_idx_= pr_idx;
        Bversion_ = 0;
    }
};

template <typename Mat, typename Vec>
void ISBBTPreCL::Apply(const Mat&, Vec& p, const Vec& c) const
{
    if (B_->Version() != Bversion_)
        Update();

    p= 0.0;
    if (kA_ != 0.0) {
//#ifndef _PAR
        solver_.Solve( *Bs_, p, VectorCL( Dprsqrtinv_*c));
//#else
//        solver_.Solve( BBT_, p, VectorCL( Dprsqrtinv_*c));
//#endif
//            std::cout << "ISBBTPreCL p: iterations: " << solver_.GetIter()
//                       << "\tresidual: " <<  solver_.GetResid();
        if (solver_.GetIter() == solver_.GetMaxIter()){
            std::cout << "ISBBTPreCL::Apply: BBT-solve: " << solver_.GetIter()
                    << '\t' << solver_.GetResid() << '\n';
        }
        p= kA_*(Dprsqrtinv_*p);
    }
    if (kM_ != 0.0) {
        Vec p2_( c.size());
        solver2_.Solve( *M_, p2_, c);
//            std::cout << "\tISBBTPreCL p2: iterations: " << solver2_.GetIter()
//                       << "\tresidual: " <<  solver2_.GetResid()
//                       << '\n';
        if (solver2_.GetIter() == solver2_.GetMaxIter()){
            std::cout << "ISBBTPreCL::Apply: M-solve: " << solver2_.GetIter()
                    << '\t' << solver2_.GetResid() << '\n';
        }

        p+= kM_*p2_;
    }
}

//**************************************************************************
// Preconditioner for the instationary (Navier-) Stokes-equations.
// It is a scaled version of the Min-Commutator-PC of Elman and can be used
// with P1X-elements.
//**************************************************************************

class MinCommPreCL : public SchurPreBaseCL
{
  private:
    const MatrixCL* A_, *B_, *Mvel_, *M_;
    mutable MatrixCL* Bs_;
    mutable size_t Aversion_, Bversion_, Mvelversion_, Mversion_;
    mutable VectorCL Dprsqrtinv_, Dvelsqrtinv_;
    double  tol_;

#ifndef _PAR
    typedef NEGSPcCL SPcT_;
    SPcT_ spc_;
    mutable PCGNESolverCL<SPcT_> solver_;
#else
    mutable CompositeMatrixCL BBT_;
    typedef ParJacNEG0CL    PCSolver1T;                         ///< type of the preconditioner for solver 1
    PCSolver1T PCsolver1_;
    mutable ParPCGSolverCL<PCSolver1T> solver_;                 ///< solver for BB^T
    const IdxDescCL* vel_idx_;                                  ///< Accessing ExchangeCL for velocity
#endif
    const IdxDescCL* pr_idx_;                                   ///< Used to determine, how to represent the kernel of BB^T in case of pure Dirichlet-BCs.
    double regularize_;

    void Update () const;

  public:
#ifndef _PAR
    MinCommPreCL (const MatrixCL* A, MatrixCL* B, MatrixCL* Mvel, MatrixCL* M_pr, const IdxDescCL& pr_idx,
                  double tol=1e-2, double regularize= 0.0)
        : SchurPreBaseCL( 0, 0), A_( A), B_( B), Mvel_( Mvel), M_( M_pr), Bs_( 0),
          Aversion_( 0), Bversion_( 0), Mvelversion_( 0), Mversion_( 0),
          tol_(tol),
          spc_( /*symmetric GS*/ true), solver_( spc_, 200, tol_, /*relative*/ true),
          pr_idx_( &pr_idx), regularize_( regularize) {}

    MinCommPreCL (const MinCommPreCL & pc)
        : SchurPreBaseCL( pc.kA_, pc.kM_), A_( pc.A_), B_( pc.B_), Mvel_( pc.Mvel_), M_( pc.M_),
          Bs_( pc.Bs_ == 0 ? 0 : new MatrixCL( *pc.Bs_)),
          Aversion_( pc.Aversion_), Bversion_( pc.Bversion_), Mvelversion_( pc.Mvelversion_),
          Mversion_( pc.Mversion_),
          Dprsqrtinv_( pc.Dprsqrtinv_), Dvelsqrtinv_( pc.Dvelsqrtinv_), tol_(pc.tol_),
          spc_( pc.spc_), solver_( spc_, 200, tol_, /*relative*/ true),
          pr_idx_( pc.pr_idx_), regularize_( pc.regularize_) {}
#else
    MinCommPreCL (const MatrixCL* A, MatrixCL* B, MatrixCL* Mvel, MatrixCL* M_pr, const IdxDescCL& pr_idx, const IdxDescCL& vel_idx,
                  double tol=1e-2, double regularize= 0.0)
        : SchurPreBaseCL( 0, 0), A_( A), B_( B), Mvel_( Mvel), M_( M_pr), Bs_( 0),
          Aversion_( 0), Bversion_( 0), Mvelversion_( 0), Mversion_( 0),
          tol_(tol), BBT_( 0, TRANSP_MUL, 0, MUL, vel_idx, pr_idx),
          PCsolver1_( pr_idx), solver_( 800, tol_, pr_idx, PCsolver1_, /*relative*/ true),
          vel_idx_( &vel_idx), pr_idx_( &pr_idx), regularize_( regularize) {}

    MinCommPreCL (const MinCommPreCL & pc)
        : SchurPreBaseCL( pc.kA_, pc.kM_), A_( pc.A_), B_( pc.B_), Mvel_( pc.Mvel_), M_( pc.M_),
          Bs_( pc.Bs_ == 0 ? 0 : new MatrixCL( *pc.Bs_)),
          Aversion_( pc.Aversion_), Bversion_( pc.Bversion_), Mvelversion_( pc.Mvelversion_),
          Mversion_( pc.Mversion_),
          Dprsqrtinv_( pc.Dprsqrtinv_), Dvelsqrtinv_( pc.Dvelsqrtinv_), tol_(pc.tol_), BBT_( 0, TRANSP_MUL, 0, MUL, *pc.vel_idx_, *pc.pr_idx_),
          PCsolver1_( *pc.pr_idx_), solver_( 800, tol_, *pc.pr_idx_, PCsolver1_, /*relative*/ true),
          vel_idx_( pc.vel_idx_), pr_idx_( pc.pr_idx_), regularize_( pc.regularize_) {}
#endif

    MinCommPreCL& operator= (const MinCommPreCL&) {
        throw DROPSErrCL( "MinCommPreCL::operator= is not permitted.\n");
    }

    ~MinCommPreCL () { delete Bs_; }

    template <typename Mat, typename Vec>
    void Apply (const Mat&, Vec& x, const Vec& b) const;
    void Apply(const MatrixCL& A,   VectorCL& x, const VectorCL& b) const { Apply<>( A, x, b); }
    void Apply(const MLMatrixCL& A, VectorCL& x, const VectorCL& b) const { Apply<>( A, x, b); }

    void SetMatrixA  (const MatrixCL* A) { A_= A; }
    void SetMatrices (const MatrixCL* A, const MatrixCL* B, const MatrixCL* Mvel, const MatrixCL* M,
                      const IdxDescCL* pr_idx) {
        A_= A;
        B_= B;
        Mvel_= Mvel;
        M_= M;
        pr_idx_= pr_idx;
        Aversion_ = Bversion_ = Mvelversion_ = Mversion_ = 0;
    }
    /// \name Parallel preconditioner setup ...
    //@{
    bool NeedDiag() const { return false; }
    void SetDiag(const VectorCL&) {}        // just for consistency
    template<typename Mat>
    void SetDiag(const Mat&) {}             // just for consistency
    bool RetAcc()   const { return true; }
    //@}
};

template <typename Mat, typename Vec>
  void
  MinCommPreCL::Apply (const Mat&, Vec& x, const Vec& b) const
{
    if ((A_->Version() != Aversion_) || (Mvel_->Version() != Mvelversion_) || (B_->Version() != Bversion_))
        Update();

    VectorCL y( b.size());
#ifndef _PAR
    solver_.Solve( *Bs_, y, VectorCL( Dprsqrtinv_*b));
#else
    solver_.Solve( BBT_, y, VectorCL( Dprsqrtinv_*b));
#endif
    if (solver_.GetIter() == solver_.GetMaxIter())
        std::cout << "MinCommPreCL::Apply: 1st BBT-solve: " << solver_.GetIter()
                  << '\t' << solver_.GetResid() << '\n';
    y*= Dprsqrtinv_;
    VectorCL z( Dprsqrtinv_*((*B_)*VectorCL( Dvelsqrtinv_*Dvelsqrtinv_*
        ( (*A_)*VectorCL( Dvelsqrtinv_*Dvelsqrtinv_*transp_mul( *B_, y)) ))));
    VectorCL t( b.size());
#ifndef _PAR
    solver_.Solve( *Bs_, t, z);
#else
    solver_.Solve( BBT_, t, z);
#endif
    if (solver_.GetIter() == solver_.GetMaxIter())
        std::cout << "MinCommPreCL::Apply: 2nd BBT-solve: " << solver_.GetIter()
                  << '\t' << solver_.GetResid() << '\n';
    x= Dprsqrtinv_*t;
}


//**************************************************************************
// Preconditioner for the instationary (Navier-) Stokes-equations.
// It uses the approximate Schur complement B diag(L)^{-1} B^T
// with L the upper left block of the saddle point system, and can be used
// with P1X-elements.
//**************************************************************************
class BDinvBTPreCL: public SchurPreBaseCL
{
  private:
    const MatrixCL* L_, *B_, *Mvel_, *M_;
    mutable MatrixCL* Bs_;
    mutable size_t Lversion_, Bversion_, Mvelversion_, Mversion_;
    mutable VectorCL Dprsqrtinv_, Dvelinv_, DSchurinv_;
    double  tol_;
#ifndef _PAR
    mutable DiagPcCL diagVelPc_, diagSchurPc_;
    typedef ApproximateSchurComplMatrixCL<DiagPcCL,MatrixCL> AppSchurComplMatrixT;
#else
    mutable ParDiagPcCL diagVelPc_, diagSchurPc_;
    typedef ParApproximateSchurComplMatrixCL<ParDiagPcCL,MatrixCL,ExchangeCL> AppSchurComplMatrixT;
#endif
    mutable AppSchurComplMatrixT *BDinvBT_;
#ifndef _PAR
    mutable PCGSolverCL<DiagPcCL> solver_;
#else
    mutable ParPCGSolverCL<ParDiagPcCL> solver_;
    const IdxDescCL* vel_idx_;
#endif
    const IdxDescCL* pr_idx_;                                   ///< Used to determine, how to represent the kernel of BB^T in case of pure Dirichlet-BCs.
    double regularize_;
    bool lumped_;

    void Update () const;

  public:
#ifndef _PAR
    BDinvBTPreCL (const MatrixCL* L, MatrixCL* B, MatrixCL* M_vel, MatrixCL* M_pr, const IdxDescCL& pr_idx,
                  double tol=1e-2, double regularize= 0.0)
        : SchurPreBaseCL( 0, 0), L_( L), B_( B), Mvel_( M_vel), M_( M_pr), Bs_( 0),
          Lversion_( 0), Bversion_( 0), Mvelversion_( 0), Mversion_( 0), tol_(tol),
          diagVelPc_(Dvelinv_), diagSchurPc_(DSchurinv_), BDinvBT_(0),
          solver_( diagSchurPc_, 200,  tol_, /*relative*/ true), pr_idx_( &pr_idx),
          regularize_( regularize), lumped_(false) {}

    BDinvBTPreCL (const BDinvBTPreCL & pc)
        : SchurPreBaseCL( pc.kA_, pc.kM_), L_( pc.L_), B_( pc.B_), Mvel_( pc.Mvel_), M_( pc.M_),
          Bs_( pc.Bs_ == 0 ? 0 : new MatrixCL( *pc.Bs_)),
          Lversion_( pc.Lversion_), Bversion_( pc.Bversion_), Mversion_( pc.Mversion_),
          Dprsqrtinv_( pc.Dprsqrtinv_), Dvelinv_( pc.Dvelinv_), DSchurinv_( pc.DSchurinv_), tol_(pc.tol_),
          diagVelPc_( Dvelinv_), diagSchurPc_( DSchurinv_), BDinvBT_(0),
          solver_( diagSchurPc_, 200, tol_, /*relative*/ true), pr_idx_( pc.pr_idx_),
          regularize_( pc.regularize_), lumped_( pc.lumped_) {}
#else
    BDinvBTPreCL (const MatrixCL* L, MatrixCL* B, MatrixCL* M_vel, MatrixCL* M_pr, const IdxDescCL& vel_idx, const IdxDescCL& pr_idx,
                  double tol=1e-2, double regularize= 0.0)
        : SchurPreBaseCL( 0, 0), L_( L), B_( B), Mvel_( M_vel), M_( M_pr), Bs_( 0),
          Lversion_( 0), Bversion_( 0), Mvelversion_( 0), Mversion_( 0), tol_(tol),
          diagVelPc_(vel_idx, Dvelinv_), diagSchurPc_(pr_idx, DSchurinv_), BDinvBT_(0),
          solver_( 200, tol_, pr_idx, diagSchurPc_, /*relative*/ true), vel_idx_( &vel_idx), pr_idx_( &pr_idx),
          regularize_( regularize), lumped_(false) {}

    BDinvBTPreCL (const BDinvBTPreCL & pc)
        : SchurPreBaseCL( pc.kA_, pc.kM_), L_( pc.L_), B_( pc.B_), Mvel_( pc.Mvel_), M_( pc.M_),
          Bs_( pc.Bs_ == 0 ? 0 : new MatrixCL( *pc.Bs_)),
          Lversion_( pc.Lversion_), Bversion_( pc.Bversion_), Mversion_( pc.Mversion_),
          Dprsqrtinv_( pc.Dprsqrtinv_), Dvelinv_( pc.Dvelinv_), DSchurinv_( pc.DSchurinv_), tol_(pc.tol_),
          diagVelPc_( *pc.vel_idx_, Dvelinv_), diagSchurPc_( *pc.pr_idx_, DSchurinv_), BDinvBT_(0),
          solver_( 200, tol_, *pc.pr_idx_, diagSchurPc_, /*relative*/ true), vel_idx_( pc.vel_idx_), pr_idx_( pc.pr_idx_),
          regularize_( pc.regularize_), lumped_( pc.lumped_) {}

    /// \name Parallel preconditioner setup ...
    //@{
    bool NeedDiag() const { return false; }
    void SetDiag(const VectorCL&) {}        // just for consistency
    template<typename Mat>
    void SetDiag(const Mat&) {}             // just for consistency
    bool RetAcc()   const { return true; }
    const ExchangeCL& GetEx() const
    {
        return pr_idx_->GetEx();
    }
    //@}
#endif

    BDinvBTPreCL& operator= (const BDinvBTPreCL&) {
        throw DROPSErrCL( "BDinvBTPreCL::operator= is not permitted.\n");
    }

    ~BDinvBTPreCL ();

    template <typename Mat, typename Vec>
    void Apply (const Mat&, Vec& x, const Vec& b) const;
    void Apply(const MatrixCL& A,   VectorCL& x, const VectorCL& b) const { Apply<>( A, x, b); }
    void Apply(const MLMatrixCL& A, VectorCL& x, const VectorCL& b) const { Apply<>( A, x, b); }

    void SetMatrixA  (const MatrixCL* L) { L_= L; }
    void SetMatrices (const MatrixCL* L, const MatrixCL* B, const MatrixCL* M_vel, const MatrixCL* M_pr,
                      const IdxDescCL* pr_idx) {
        L_= L;
        B_= B;
        Mvel_= M_vel;
        M_= M_pr;
        pr_idx_= pr_idx;
        Lversion_ = Bversion_ = Mvelversion_= Mversion_ = 0;
    }
    /// Returns true, if the lumped diag of the velocity mass matrix is used, and false, if the diagonal of the velocity convection-diffusion-reaction matrix is considered.
    bool UsesMassLumping() const { return lumped_; }
    /// For lump==true, the lumped diag of the velocity mass matrix is used. Otherwise the diagonal of the velocity convection-diffusion-reaction matrix is considered.
    void SetMassLumping( bool lump) { lumped_= lump; }
    /// If lumping is switched on, the lumped diag of the velocity mass matrix is returned. Otherwise the diagonal of the velocity convection-diffusion-reaction matrix is returned.
    VectorCL GetVelDiag() const {
        if ((L_->Version() != Lversion_) || (Mvel_->Version() != Mvelversion_))
            Update();
        return VectorCL(1.0/Dvelinv_);
    }
};

template <typename Mat, typename Vec>
  void
  BDinvBTPreCL::Apply (const Mat&, Vec& x, const Vec& b) const
{
    if ((L_->Version() != Lversion_) || (Mvel_->Version() != Mvelversion_) || (M_->Version() != Mversion_) || (B_->Version() != Bversion_))
        Update();

    Vec y( b.size());
    solver_.Solve( *BDinvBT_, y, Vec( Dprsqrtinv_*b));
    if (solver_.GetIter() == solver_.GetMaxIter())
        std::cout << "BDinvBTPreCL::Apply: BLBT-solve: " << solver_.GetIter()
                  << '\t' << solver_.GetResid() << '\n';
    x= Dprsqrtinv_*y;
}

/// Upper block-triangular preconditioning strategy in BlockPreCL
struct UpperBlockPreCL
{
    template <class PC1T, class PC2T, class Mat, class Vec>
    static void
    Apply (const PC1T& pc1, const PC2T& pc2, const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c) {
        pc2.Apply( /*dummy*/ B, p, c);
        p*= -1.;
        Vec b2( b);
        b2-= transp_mul( B, p);
        pc1.Apply( A, v, b2);
    }
};

/// Block-diagonal preconditioning strategy in BlockPreCL
struct DiagBlockPreCL
{
    template <class PC1T, class PC2T, class Mat, class Vec>
    static void
    Apply (const PC1T& pc1, const PC2T& pc2, const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c) {
        pc1.Apply( A, v, b);
        pc2.Apply( /*dummy*/ B, p, c);
        p*= -1.;
   }
};

/// Block-diagonal preconditioning strategy in BlockPreCL if an spd preconditioner is required
struct DiagSpdBlockPreCL
{
    template <class PC1T, class PC2T, class Mat, class Vec>
    static void
    Apply (const PC1T& pc1, const PC2T& pc2, const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c) {
        pc1.Apply( A, v, b);
        pc2.Apply( /*dummy*/ B, p, c);
   }
};

/// Lower block-triangular preconditioning strategy in BlockPreCL
struct LowerBlockPreCL
{
    template <class PC1T, class PC2T, class Mat, class Vec>
    static void
    Apply (const PC1T& pc1, const PC2T& pc2, const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c) {
        pc1.Apply( A, v, b);
#ifdef _PAR
        Assert(pc1.RetAcc(), DROPSErrCL("LowerBlockPreCL::Apply: Accumulation is missing"), DebugParallelNumC);
#endif
        Vec c2( B*v - c);
        pc2.Apply( /*dummy*/ B, p, c2);
    }
};

/// SIMPLER type preconditioning strategy in BlockPreCL
struct SIMPLERBlockPreCL
{
    template <class PC1T, class PC2T, class Mat, class Vec>
    static void
    Apply (PC1T& pc1, PC2T& pc2, const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c) {
        Vec D( pc2.GetVelDiag()), dp(p);

        // find some initial pressure for the SIMPLER preconditioner (initial p=0 would result in SIMPLE)
#ifdef _PAR
        pc2.Apply( B, p, Vec(B*Vec(pc1.GetEx().GetAccumulate(b)/D) - c));
        if (!pc2.RetAcc())
            pc2.GetEx().Accumulate(p);
#else
        pc2.Apply( B, p, Vec(B*Vec(b/D) - c));
#endif
        // from here on it's the SIMPLE preconditioner
        pc1.Apply( A, v, Vec(b - transp_mul( B, p)));
#ifdef _PAR
        if (!pc1.RetAcc()) pc1.GetEx().Accumulate(v);
#endif

        pc2.Apply( /*dummy*/ B, dp, Vec(B*v - c));

#ifdef _PAR
        if (!pc2.RetAcc()) pc2.GetEx().Accumulate(dp);
        v-= pc1.GetEx().GetAccumulate(Vec(transp_mul( B, dp)/D));
#else
        v-= transp_mul( B, dp)/D;
#endif
        p+= dp;
   }
};

// With BlockShapeT= DiagBlockPreCL, this is the diagonal PC ( pc1^(-1) 0 \\ 0 pc2^(-1) ),
// else it is block-triangular:
// Upper... ( pc1^(-1) B^T \\ 0 pc2^(-1) ), resp. lower: ( pc1^(-1) 0 \\ B pc2^(-1) )
template <class PC1T, class PC2T, class BlockShapeT= DiagBlockPreCL>
class BlockPreCL
{
  private:
    PC1T& pc1_; // Preconditioner for A.
    PC2T& pc2_; // Preconditioner for S.

  public:
    BlockPreCL (PC1T& pc1, PC2T& pc2)
        : pc1_( pc1), pc2_( pc2) {}

    template <typename Mat, typename Vec>
    void
    Apply(const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c) const {
        BlockShapeT::Apply( pc1_, pc2_, A, B, v, p, b, c);
    }

    template <typename Mat, typename Vec>
    void
    Apply(const BlockMatrixBaseCL<Mat>& A, Vec& x, const Vec& b) const {
        VectorCL b0( b[std::slice( 0, A.num_rows( 0), 1)]);
        VectorCL b1( b[std::slice( A.num_rows( 0), A.num_rows( 1), 1)]);
        VectorCL x0( A.num_cols( 0));
        VectorCL x1( A.num_cols( 1));
        BlockShapeT::Apply( pc1_, pc2_, *A.GetBlock( 0), *A.GetBlock( 2), x0, x1, b0, b1);
        x[std::slice( 0, A.num_cols( 0), 1)]= x0;
        x[std::slice( A.num_cols( 0), A.num_cols( 1), 1)]= x1;
    }

#ifdef _PAR
    /// \brief Check if the preconditioned vector is accumulated
    bool RetAcc() const {
        Assert( pc1_.RetAcc()==pc2_.RetAcc(), DROPSErrCL("BlockPreCL::RetAcc: Preconditioners do not match"),
                DebugNumericC);
        return pc1_.RetAcc();
    }

    /// \brief Check if the diagonal of the matrix needs to be computed
    bool NeedDiag() const { return pc1_.NeedDiag() || pc2_.NeedDiag(); }

    /// \brief Set accumulated diagonal of a matrix, that is needed by most of the preconditioners
    template<typename Mat>
    void SetDiag(const Mat& A)
    {
        pc1_.SetDiag(*A.GetBlock( 0));
        pc2_.SetDiag(/*dummy*/ *(A.GetBlock( 3)!=0 ? A.GetBlock( 3) : A.GetBlock( 1)));
    }
    const PC1T& GetPC1() const { return pc1_; }
          PC1T& GetPC1()       { return pc1_; }
    const PC2T& GetPC2() const { return pc2_; }
          PC2T& GetPC2()       { return pc2_; }

#endif
};


template <typename Mat, typename Vec>
void ISPreCL::Apply(const Mat&, Vec& p, const Vec& c) const
{
//    double new_res;
//    double old_res= norm( c);
    ssor_.Apply( A_, p, c);
//    std::cout << " residual: " <<  (new_res= norm( A_*p - c)) << '\t';
//    std::cout << " reduction: " << new_res/old_res << '\t';
    p*= kA_;
//    double mnew_res;
//    double mold_res= norm( c);
    Vec p2_( c.size());
    ssor_.Apply( M_, p2_, c);
//    std::cout << " residual: " <<  (mnew_res= norm( M_*p2_ - c)) << '\t';
//    std::cout << " reduction: " << mnew_res/mold_res << '\n';
    p+= kM_*p2_;
}

#ifndef _PAR
template <class SolverT>
template <typename Mat, typename Vec>
void ISNonlinearPreCL<SolverT>::Apply(const Mat&, Vec& p, const Vec& c) const
{
    p= 0.0;
    if (kA_ != 0.0) {
        solver_.Solve( A_, p, c);
        if (solver_.GetIter() == solver_.GetMaxIter())
        std::cout << "ISNonlinearPreCL::Apply (1st solve: iterations: " << solver_.GetIter()
                  << "\tresidual: " <<  solver_.GetResid() << '\n';
        p*= kA_;
    }
    if ( kM_ != 0.0) {
        Vec p2_( c.size());
        solver_.Solve( M_, p2_, c);
        if (solver_.GetIter() == solver_.GetMaxIter())
        std::cout << "ISNonlinearPreCL::Apply (2nd solve: iterations: " << solver_.GetIter()
                  << "\tresidual: " <<  solver_.GetResid() << '\n';
        p+= kM_*p2_;
    }
}
#else
template <typename ASolverT, typename MSolverT>
template <typename Mat, typename Vec>
void ISNonlinearPreCL<ASolverT, MSolverT>::Apply(const Mat&, Vec& p, const Vec& c) const
{
    p= 0.0;
    if (kA_ != 0.0) {
        Asolver_.Solve( A_, p, c);
//       std::cout << "ISNonlinearPreCL p: iterations: " << solver_.GetIter()
//                 << "\tresidual: " <<  solver_.GetResid()<<std::endl;
        p*= kA_;
    }
    if ( kM_ != 0.0) {
        Vec p2_( c.size());
        Msolver_.Solve( M_, p2_, c);
//           std::cout << "\t p2: iterations: " << solver_.GetIter()
//                     << "\tresidual: " <<  solver_.GetResid()
//                     << '\n';
        p+= kM_*p2_;
    }
}
#endif

/* with orthogonalization
template <typename Mat, typename Vec>
void ISNonlinearPreCL::Apply(const Mat&, Vec& p, const Vec& c) const
{
    VectorCL e( 1., p.size());
    const double ee= p.size(); // = dot( e, e);
    VectorCL c2( (dot(c, e)/ee)*e);
    VectorCL c3( c - c2);
std::cout << "norm( c): " << norm( c) << "\tnorm( e): " << norm( e)
          << "\tdot(c, e)/norm( c)/norm( e): " << dot(c, e)/norm( c)/ee << '\n';
    p= 0.0;
    solver_.Solve( A_, p, c3);
    p+= c2;
    std::cout << "ISNonlinearPreCL p: iterations: " << solver_.GetIter()
              << "\tresidual: " <<  solver_.GetResid();
    p*= kA_;
    if ( kM_ != 0.0) {
        Vec p2_( c.size());
        solver_.Solve( M_, p2_, c);
std::cout << "norm( p2): " << norm( p2_);
        std::cout << "\t p2: iterations: " << solver_.GetIter()
                  << "\tresidual: " <<  solver_.GetResid()
                  << '\n';
        p+= kM_*p2_;
    }
}
*/

template<class SmootherCL, class DirectSolverCL, class ProlongationIterT>
void
MGMPr(const std::vector<VectorCL>::const_iterator& ones,
      const MLMatrixCL::const_iterator& begin, const MLMatrixCL::const_iterator& fine,
      ProlongationIterT P, VectorCL& x, const VectorCL& b,
      const SmootherCL& Smoother, const Uint smoothSteps,
      DirectSolverCL& Solver, const int numLevel, const int numUnknDirect)
// Multigrid method, V-cycle. If numLevel==0 or #Unknowns <= numUnknDirect,
// the direct solver Solver is used.
// If one of the parameters is -1, it will be neglected.
// If MLMatrixCL.begin() has been reached, the direct solver is used too.
// Concerning the stabilization see Hackbusch Multigrid-Methods and Applications;
// Basically we project on the orthogonal complement of the kernel of A before
// the coarse-grid correction.
{
    MLMatrixCL::const_iterator coarse = fine;
    ProlongationIterT coarseP= P;
    if(  ( numLevel==-1      ? false : numLevel==0 )
       ||( numUnknDirect==-1 ? false : x.size() <= static_cast<Uint>(numUnknDirect) )
       || fine==begin)
    { // use direct solver
        Solver.Solve( *fine, x, b);
        x-= dot( *ones, x);
        return;
    }
    --coarse;
    --coarseP;
    // presmoothing
    for (Uint i=0; i<smoothSteps; ++i) Smoother.Apply( *fine, x, b);
    // restriction of defect
    VectorCL d( transp_mul( *P, VectorCL( b - (*fine)*x)));
    d-= dot( *(ones-1), d);
    VectorCL e( d.size());
    // calculate coarse grid correction
    MGMPr( ones-1, begin, coarse, coarseP, e, d, Smoother, smoothSteps, Solver, (numLevel==-1 ? -1 : numLevel-1), numUnknDirect);
    // add coarse grid correction
    x+= (*P) * e;
    // postsmoothing
    for (Uint i=0; i<smoothSteps; ++i) Smoother.Apply( *fine, x, b);
    // This projection could probably be avoided, but it is cheap and with it,
    // we are on the safe side.
    x-= dot( *ones, x);
}

template<class ProlongationT>
void ISMGPreCL<ProlongationT>::MaybeInitOnes() const
{
    if (Mpr_.size() == ones_.size()) return;
    // Compute projection on constant pressure function only once.
    Uint i= 0;
    ones_.resize(0); // clear all
    ones_.resize(Mpr_.size());
    for (MLMatrixCL::const_iterator it= Mpr_.begin(); it != Mpr_.end(); ++it, ++i) {
        ones_[i].resize( it->num_cols(), 1.0/it->num_cols());
    }
}

template<class ProlongationT>
template <typename Mat, typename Vec>
void
ISMGPreCL<ProlongationT>::Apply(const Mat& /*A*/, Vec& p, const Vec& c) const
{
    MaybeInitOnes();
    p= 0.0;
    const Vec c2_( c - dot( ones_.back(), c));
    typename ProlongationT::const_iterator finestP = --P_.end();
//    double new_res= (Apr_.back().A.Data*p - c).norm();
//    double old_res;
//    std::cout << "Pressure: iterations: " << iter_prA_ <<'\t';
    for (DROPS::Uint i=0; i<iter_prA_; ++i) {
        DROPS::MGMPr( ones_.end()-1, Apr_.begin(), --Apr_.end(), finestP, p, c2_, smoother, sm, solver, lvl, -1);
//        old_res= new_res;
//        std::cout << " residual: " <<  (new_res= (Apr_.back().A.Data*p - c).norm()) << '\t';
//        std::cout << " reduction: " << new_res/old_res << '\n';
    }
    p*= kA_;
//    std::cout << " residual: " <<  (Apr_.back().A.Data*p - c).norm() << '\t';

    Vec p2( p.size());
    for (DROPS::Uint i=0; i<iter_prM_; ++i)
        DROPS::MGM( Mpr_.begin(), --Mpr_.end(), finestP, p2, c, smoother, sm, solver, lvl, -1);
//    std::cout << "Mass: iterations: " << iter_prM_ << '\t'
//              << " residual: " <<  (Mpr_.back().A.Data*p2 - c).norm() << '\n';

    p+= kM_*p2;
}


} // end of namespace DROPS

#endif /* OSEENPRECOND_H_ */
