/// \file precond.h
/// \brief preconditioner for Ax = b
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

#ifndef PRECOND_H_
#define PRECOND_H_

#include "misc/container.h"
#include "num/spmat.h"
#include "num/spblockmat.h"
#ifdef _PAR
#include "parallel/exchange.h"
#endif

namespace DROPS{

//*****************************************************************************
//
//  Gauss-Seidel type methods for preconditioning
//
//*****************************************************************************

//=============================================================================
//  Template magic for the selection of the preconditioner
//=============================================================================

// Available methods
enum PreMethGS
{
    P_JAC,     //  0 Jacobi
    P_GS,      //  1 Gauss-Seidel
    P_SGS,     //  2 symmetric Gauss-Seidel
    P_SGS0,    //  3 symmetric Gauss-Seidel with initial vector 0
    P_JOR,     //  4 P_JAC with over-relaxation
    P_SOR,     //  5 P_GS with over-relaxation
    P_SSOR,    //  6 P_SGS with over-relaxation
    P_SSOR0,   //  7 P_SGS0 with over-relaxation
    P_SGS0_D,  //  8 P_SGS0 using SparseMatDiagCL
    P_SSOR0_D, //  9 P_SSOR0 using SparseMatDiagCL
    P_DUMMY,   // 10 identity
    P_GS0,     // 11 Gauss-Seidel with initial vector 0
    P_JAC0     // 12 Jacobi with initial vector 0
};

// Base methods
enum PreBaseGS { PB_JAC, PB_GS, PB_SGS, PB_SGS0, PB_DUMMY, PB_GS0, PB_JAC0 };

// Properties of the methods
template <PreMethGS PM> struct PreTraitsCL
{
    static const PreBaseGS BaseMeth= PreBaseGS(PM<8 ? PM%4
                                     : (PM==10 ? PB_DUMMY : (PM==11 ? PB_GS0 : (PM==12 ? PB_JAC0 : PB_SGS0))) );
    static const bool      HasOmega= (PM>=4 && PM<8) || PM==9 || PM==11 || PM==12;
    static const bool      HasDiag=  PM==8 || PM==9;
};

// Used to make a distinct type from each method
template <PreBaseGS> class PreDummyCL {};


//=============================================================================
//  Implementation of the methods
//=============================================================================

// One step of the Jacobi method with start vector x
template <bool HasOmega, typename Vec>
void
SolveGSstep(const PreDummyCL<PB_JAC>&, const MatrixCL& A, Vec& x, const Vec& b, double omega)
{
    const size_t n= A.num_rows();
    Vec          y(x.size());

    for (size_t i=0, nz=0; i<n; ++i)
    {
        double aii=0, sum= b[i];
        for (const size_t end= A.row_beg(i+1); nz<end; ++nz)
            if (A.col_ind(nz) != i)
                sum-= A.val(nz)*x[A.col_ind(nz)];
            else
                aii= A.val(nz);
        if (HasOmega)
            y[i]= (1.-omega)*x[i]+omega*sum/aii;
        else
            y[i]= sum/aii;
    }

    std::swap(x,y);
}

// One step of the Jacobi method with start vector 0
template <bool HasOmega, typename Vec>
void
SolveGSstep(const PreDummyCL<PB_JAC0>&, const MatrixCL& A, Vec& x, const Vec& b, double omega)
{
    const size_t n= A.num_rows();
    size_t nz;

    for (size_t i= 0; i < n; ++i) {
        nz= A.row_beg( i);
        for (const size_t end= A.row_beg( i+1); A.col_ind( nz) != i && nz < end; ++nz) ; // empty loop
        if (HasOmega)
            x[i]= /*(1.-omega)*x[i]=0  + */ omega*b[i]/A.val( nz);
        else
            x[i]= b[i]/A.val( nz);
    }
}

// One step of the Gauss-Seidel/SOR method with start vector x
template <bool HasOmega, typename Vec>
void
SolveGSstep(const PreDummyCL<PB_GS>&, const MatrixCL& A, Vec& x, const Vec& b, double omega)
{
    const size_t n= A.num_rows();
    double aii, sum;

    for (size_t i=0, nz=0; i<n; ++i) {
        sum= b[i];
        const size_t end= A.row_beg( i+1);
        for (; A.col_ind( nz) != i; ++nz) // This is safe: Without diagonal entry, Gauss-Seidel would explode anyway.
            sum-= A.val( nz)*x[A.col_ind( nz)];
        aii= A.val( nz++);
        for (; nz<end; ++nz)
            sum-= A.val( nz)*x[A.col_ind( nz)];
        if (HasOmega)
            x[i]= (1.-omega)*x[i]+omega*sum/aii;
        else
            x[i]= sum/aii;
    }
}

// One step of the Gauss-Seidel/SOR method with start vector x
template <bool HasOmega, typename Vec>
void
SolveGSstep(const PreDummyCL<PB_GS0>&, const MatrixCL& A, Vec& x, const Vec& b, double omega)
{
    const size_t n= A.num_rows();
    double aii, sum;

    for (size_t i=0, nz=0; i<n; ++i) {
        sum= b[i];
        const size_t end= A.row_beg( i+1);
        for (; A.col_ind( nz) != i; ++nz) // This is safe: Without diagonal entry, Gauss-Seidel would explode anyway.
            sum-= A.val( nz)*x[A.col_ind( nz)];
        aii= A.val( nz);
        nz= end;
        if (HasOmega)
            x[i]= (1.-omega)*x[i]+omega*sum/aii;
        else
            x[i]= sum/aii;
    }
}


// One step of the Symmetric-Gauss-Seidel/SSOR method with start vector x
template <bool HasOmega, typename Vec>
void
SolveGSstep(const PreDummyCL<PB_SGS>&, const MatrixCL& A, Vec& x, const Vec& b, double omega)
{
    const size_t n= A.num_rows();
    double aii, sum;

    for (size_t i=0, nz=0; i<n; ++i) {
        sum= b[i];
        const size_t end= A.row_beg( i+1);
        for (; A.col_ind( nz) != i; ++nz) // This is safe: Without diagonal entry, Gauss-Seidel would explode anyway.
            sum-= A.val( nz)*x[A.col_ind( nz)];
        aii= A.val( nz++);
        for (; nz<end; ++nz)
            sum-= A.val( nz)*x[A.col_ind( nz)];
        if (HasOmega)
            x[i]= (1.-omega)*x[i]+omega*sum/aii;
        else
            x[i]= sum/aii;
    }
    for (size_t i= n, nz= A.row_beg( n); i>0; ) { // This is safe: Without diagonal entry, Gauss-Seidel would explode anyway.
        --i;
        double aii, sum= b[i];
        const size_t beg= A.row_beg( i);
        for (; A.col_ind( --nz) != i; ) {
            sum-= A.val( nz)*x[A.col_ind( nz)];
        }
        aii= A.val( nz);
        for (; nz>beg; ) {
            --nz;
            sum-= A.val( nz)*x[A.col_ind( nz)];
        }
        if (HasOmega)
            x[i]= (1.-omega)*x[i]+omega*sum/aii;
        else
            x[i]= sum/aii;
    }
}


// One step of the Symmetric-Gauss-Seidel/SSOR method with start vector 0
template <bool HasOmega, typename Vec>
void
SolveGSstep(const PreDummyCL<PB_SGS0>&, const MatrixCL& A, Vec& x, const Vec& b, double omega)
{
    const size_t n= A.num_rows();

    for (size_t i=0; i<n; ++i)
    {
        double sum= b[i];
        size_t j= A.row_beg(i);
        for ( ; A.col_ind(j) < i; ++j)
            sum-= A.val(j)*x[A.col_ind(j)];
        if (HasOmega)
            x[i]= omega*sum/A.val(j);
        else
            x[i]= sum/A.val(j);
    }

    for (size_t i=n; i>0; )
    {
        --i;
        double sum= 0;
        size_t j= A.row_beg(i+1)-1;
        for ( ; A.col_ind(j) > i; --j)
            sum-= A.val(j)*x[A.col_ind(j)];
        if (HasOmega)
            x[i]= (2.-omega)*x[i]+omega*sum/A.val(j);
        else
            x[i]+= sum/A.val(j);
    }
}


// One step of the Symmetric-Gauss-Seidel/SSOR method with start vector 0,
// uses SparseMatDiagCL for the location of the diagonal
template <bool HasOmega, typename Vec>
void
SolveGSstep(const PreDummyCL<PB_SGS0>&, const MatrixCL& A, Vec& x, const Vec& b, const SparseMatDiagCL& diag, double omega)
{
    const size_t n= A.num_rows();

    for (size_t i=0; i<n; ++i)
    {
        double sum= b[i];
        for (size_t j= A.row_beg(i); j < diag[i]; ++j)
            sum-= A.val(j)*x[A.col_ind(j)];
        if (HasOmega)
            x[i]= omega*sum/A.val(diag[i]);
        else
            x[i]= sum/A.val(diag[i]);
    }

    for (size_t i=n; i>0; )
    {
        --i;
        double sum= 0;
        for (size_t j= A.row_beg(i+1)-1; j > diag[i]; --j)
            sum-= A.val(j)*x[A.col_ind(j)];
        if (HasOmega)
            x[i]= (2.-omega)*x[i]+omega*sum/A.val(diag[i]);
        else
            x[i]+= sum/A.val(diag[i]);
    }
}

template <bool HasOmega, typename  Vec, PreBaseGS PBT>
void
SolveGSstep(const PreDummyCL<PBT>& pd, const MLMatrixCL& M, Vec& x, const Vec& b, double omega)
{
    SolveGSstep<HasOmega, Vec>( pd, M.GetFinest(), x, b, omega);
}

template <bool HasOmega, typename  Vec, PreBaseGS PBT>
void
SolveGSstep(const PreDummyCL<PBT>& pd, const MLMatrixCL& M, Vec& x, const Vec& b)
{
    SolveGSstep<HasOmega, Vec>( pd, M.GetFinest(), x, b);
}

template <bool HasOmega, typename  Vec, PreBaseGS PBT>
void
SolveGSstep(const PreDummyCL<PBT>& pd, const MLMatrixCL& A, Vec& x, const Vec& b, const SparseMatDiagCL& diag, double omega)
{
    SolveGSstep<HasOmega, Vec>( pd, A.GetFinest(), x, b, diag, omega);
}

template <bool HasOmega, typename  Vec, PreBaseGS PBT>
void
SolveGSstep(const PreDummyCL<PBT>& pd, const MLMatrixCL& A, Vec& x, const Vec& b, const SparseMatDiagCL& diag)
{
    SolveGSstep<HasOmega, Vec>( pd, A.GetFinest(), x, b, diag);
}
//=============================================================================
//  Preconditioner classes
//=============================================================================

// TODO: Init ueberdenken.

// Preconditioners without own matrix
template <PreMethGS PM, bool HasDiag= PreTraitsCL<PM>::HasDiag> class PreGSCL;

// Simple preconditioners
template <PreMethGS PM>
class PreGSCL<PM,false>
{
  private:
    double _omega;

  public:
    PreGSCL (double om= 1.0) : _omega(om) {}

    template <typename Mat, typename Vec>
    void Apply(const Mat& A, Vec& x, const Vec& b) const
    {
        SolveGSstep<PreTraitsCL<PM>::HasOmega,Vec>(PreDummyCL<PreTraitsCL<PM>::BaseMeth>(), A, x, b, _omega);
    }
};


// Preconditioner with SparseMatDiagCL
template <PreMethGS PM>
class PreGSCL<PM,true>
{
  private:
    const SparseMatDiagCL* _diag;
    double                 _omega;

  public:
    PreGSCL (double om= 1.0) : _diag(0), _omega(om) {}
    PreGSCL (const PreGSCL& p) : _diag(p._diag ? new SparseMatDiagCL(*(p._diag)) : 0), _omega(p._omega) {}
    ~PreGSCL() { delete _diag; }

    void Init(const MatrixCL& A)
    {
        delete _diag; _diag=new SparseMatDiagCL(A);
    }

    template <typename Vec>
    void Apply(const MatrixCL& A, Vec& x, const Vec& b) const
    {
        SolveGSstep<PreTraitsCL<PM>::HasOmega,Vec>(PreDummyCL<PreTraitsCL<PM>::BaseMeth>(), A, x, b, *_diag, _omega);
    }
    template <typename Vec>
    void Apply(const MLMatrixCL& A, Vec& x, const Vec& b) const
    {
        SolveGSstep<PreTraitsCL<PM>::HasOmega,Vec>(PreDummyCL<PreTraitsCL<PM>::BaseMeth>(), A.GetFinest(), x, b, *_diag, _omega);
    }
};


// Preconditioners with own matrix
template <PreMethGS PM, bool HasDiag= PreTraitsCL<PM>::HasDiag>
class PreGSOwnMatCL;

// Simple preconditioners
template <PreMethGS PM>
class PreGSOwnMatCL<PM,false>
{
  private:
    const MatrixCL& _M;
    double          _omega;

  public:
    PreGSOwnMatCL (const MatrixCL& M, double om= 1.0) : _M(M), _omega(om) {}

    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec& x, const Vec& b) const
    {
        SolveGSstep<PreTraitsCL<PM>::HasOmega,Vec>(PreDummyCL<PreTraitsCL<PM>::BaseMeth>(), _M, x, b, _omega);
    }
};


// Preconditioner with SparseMatDiagCL
template <PreMethGS PM>
class PreGSOwnMatCL<PM,true>
{
  private:
    const MatrixCL&        _M;
    const SparseMatDiagCL* _diag;
    double                 _omega;

  public:
    PreGSOwnMatCL (const MatrixCL& M, double om= 1.0) : _M(M), _diag(0), _omega(om) {}
    PreGSOwnMatCL (const PreGSOwnMatCL&); // not defined
    ~PreGSOwnMatCL() { delete _diag; }

    void Init(const MatrixCL& A)
    {
        delete _diag; _diag=new SparseMatDiagCL(A);
    }

    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec& x, const Vec& b) const
    {
        SolveGSstep<PreTraitsCL<PM>::HasOmega,Vec>(PreDummyCL<PreTraitsCL<PM>::BaseMeth>(), _M, x, b, *_diag, _omega);
    }
};

class MultiSSORPcCL
// do multiple SSOR-steps
{
  private:
    double _omega;
    int   _num;

  public:
    MultiSSORPcCL(double om= 1.0, int num= 1) : _omega(om), _num(num) {}

    template <typename Mat, typename Vec>
    void Apply(const Mat& A, Vec& x, const Vec& b) const
    {
        // one SSOR0-step
        SolveGSstep<PreTraitsCL<P_SSOR0>::HasOmega,Vec>(PreDummyCL<PreTraitsCL<P_SSOR0>::BaseMeth>(), A, x, b, _omega);
        // _num-1 SSOR-steps
        for (int i=1; i<_num; ++i)
            SolveGSstep<PreTraitsCL<P_SSOR>::HasOmega,Vec>(PreDummyCL<PreTraitsCL<P_SSOR>::BaseMeth>(), A, x, b, _omega);
    }
};

class DummyPcCL
{
  public:
    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec& x, const Vec& b) const
    {
        x = b;
    }
    template <typename Mat, typename Vec>
    void ApplyTranspose(const Mat&, Vec& x, const Vec& b) const
    {
        x = b;
    }
    template <typename Mat, typename Vec>
    Vec transp_mul(const Mat&, const Vec& b) const
    {
        return b;
    }
};

/// \brief Apply a diagonal-matrix given as a vector.
class DiagPcCL
{
  private:
    const VectorCL& D_;

  public:
    DiagPcCL (const VectorCL& D) : D_( D) {}

    template <typename Mat, typename Vec>
    void Apply (const Mat&, Vec& x, const Vec& b) const
    {
        Assert( D_.size()==b.size(), DROPSErrCL("DiagPcCL: incompatible dimensions"), DebugNumericC);
        x= D_*b;
    }
};

/// \brief (Symmetric) Gauss-Seidel preconditioner (i.e. start-vector 0) for A*A^T (Normal Equations).
class NEGSPcCL
{
  private:
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
    NEGSPcCL (bool symmetric= true) : symmetric_( symmetric), Aaddr_( 0), Aversion_( 0) {}

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
void NEGSPcCL::Update (const Mat& A) const
{
    if (&A == Aaddr_ && Aversion_ == A.Version()) return;
    Aaddr_= &A;
    Aversion_= A.Version();

    D_.resize( A.num_rows());
    D_= BBTDiag( A);
    y_.resize( A.num_cols());
}

template <typename Mat, typename Vec>
void NEGSPcCL::ForwardGS(const Mat& A, Vec& x, const Vec& b) const
{
    // x= 0.; // implied, but superfluous, because we can assign the x-values below, not update.
    y_= 0.;
    double t;
    for (size_t i= 0; i < A.num_rows(); ++i) {
        t= (b[i] - mul_row( A, y_, i))/D_[i];
        x[i]= t;
        add_row_to_vec( A, t, y_, i); // y+= t* (i-th row of A)
    }
}

template <typename Mat, typename Vec>
void NEGSPcCL::BackwardGS(const Mat& A, Vec& x, const Vec& b) const
{
    // x= 0.; // implied, but superfluous, because we can assign the x-values below, not update.
    y_= 0.;
    for (size_t i= b.size() - 1; i < b.size(); --i) {
        x[i]= (b[i] - mul_row( A, y_, i))/D_[i];
        add_row_to_vec( A, x[i], y_, i); // y+= t* (i-th row of A)
    }
}

template <typename Mat, typename Vec>
void NEGSPcCL::Apply(const Mat& A, Vec& x, const Vec& b) const
{
    Update( A);

    ForwardGS( A, x, b);
    if (!symmetric_) return;
    BackwardGS( A, x, VectorCL( D_*x));
}

template <typename Mat, typename Vec>
void NEGSPcCL::ApplyTranspose(const Mat& A, Vec& x, const Vec& b) const
{
    Update( A);

    BackwardGS( A, x, b);
    if (!symmetric_) return;
    ForwardGS( A, x, VectorCL( D_*x));
}

template <typename Mat, typename Vec>
void NEGSPcCL::ForwardMulGS(const Mat& A, Vec& x, const Vec& b) const
{
    // x= 0.; // implied, but superfluous, because we can assign the x-values below, not update.
    y_= 0.;
    for (size_t i= 0 ; i < b.size(); ++i) {
        add_row_to_vec( A, b[i], y_, i); // y+= b[i]* (i-th row of A)
        x[i]= mul_row( A, y_, i);
    }
}

template <typename Mat, typename Vec>
void NEGSPcCL::BackwardMulGS(const Mat& A, Vec& x, const Vec& b) const
{
    // x= 0.; // implied, but superfluous, because we can assign the x-values below, not update.
    y_= 0.;
    for (size_t i= b.size() - 1 ; i < b.size(); --i) {
        add_row_to_vec( A, b[i], y_, i); // y+= b[i]* (i-th row of A)
        x[i]= mul_row( A, y_, i);
    }
}

template <typename Mat, typename Vec>
Vec NEGSPcCL::mul (const Mat& A, const Vec& b) const
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
Vec NEGSPcCL::transp_mul(const Mat& A, const Vec& b) const
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

#ifdef _PAR
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
    
    bool RetAcc() const {return false;}
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
#endif

//=============================================================================
//  Typedefs
//=============================================================================

typedef PreGSCL<P_SSOR>    SSORsmoothCL;
typedef PreGSCL<P_SOR>     SORsmoothCL;
typedef PreGSCL<P_SGS>     SGSsmoothCL;
typedef PreGSCL<P_JOR>     JORsmoothCL;
typedef PreGSCL<P_GS>      GSsmoothCL;
typedef PreGSCL<P_JAC0>    JACPcCL;
typedef PreGSCL<P_SGS0>    SGSPcCL;
typedef PreGSCL<P_SSOR0>   SSORPcCL;
typedef PreGSCL<P_SSOR0_D> SSORDiagPcCL;
typedef PreGSCL<P_GS0>     GSPcCL;


/// \brief base class for "expensive" preconditioners.
///
/// Preconditioner base class for Schur complement preconditioners (see SchurPreBaseCL)
/// and SolverAsPreCL, as those are computational expensive compared to the additional virtual function call.
/// Anyway, Jacobi and Gauss-Seidel type preconditioners are not derived from this class to avoid the virtual function overhead.
class PreBaseCL
{
  protected:
    mutable std::ostream* output_;
    mutable Uint iter_;

    PreBaseCL( std::ostream* output= 0)
      : output_( output), iter_(0) {}
    virtual ~PreBaseCL() {}

    void AddIter( int iter) const { iter_+= iter; }

  public:
    virtual void Apply( const MatrixCL& A,   VectorCL& x, const VectorCL& b) const = 0;
    virtual void Apply( const MLMatrixCL& A, VectorCL& x, const VectorCL& b) const = 0;

    /// reset iteration counter
    void ResetIter()    const { iter_= 0; }
    /// return total number of iterations
    Uint GetTotalIter() const { return iter_; }
};


//=============================================================================
// Krylov-methods as preconditioner.
// Needed for, e.g., InexactUzawa.
//=============================================================================
/// Wrapper to use a solver as preconditioner. Needed for, e.g., inexact Uzawa.
template <class SolverT>
class SolverAsPreCL: public PreBaseCL
{
  private:
    SolverT& solver_;

  public:
    SolverAsPreCL( SolverT& solver, std::ostream* output= 0)
        : PreBaseCL( output), solver_( solver) {}
    /// return solver object
    SolverT& GetSolver()             { return solver_; }
    /// return solver object
    const SolverT& GetSolver() const { return solver_; }
#ifdef _PAR
    /// \brief Check if return preconditioned vectors are accumulated after calling Apply
    bool RetAcc() const   { return true; }
    /// \brief Check if the diagonal of the matrix is needed
    bool NeedDiag() const { return solver_.GetPC().NeedDiag(); }
    /// \brief Set diagonal to the preconditioner of the solver
    template<typename Mat>
    void SetDiag(const Mat& A)
    {
        solver_.GetPC().SetDiag(A);
    }
    const ExchangeCL& GetEx() const
    {
        return solver_.GetEx();
    }
#endif

    template <typename Mat, typename Vec>
    void
    Apply(const Mat& A, Vec& x, const Vec& b) const {
        x= 0.0;
        solver_.Solve( A, x, b);
//         if (solver_.GetIter()==solver_.GetMaxIter())
//           IF_MASTER
//             std::cout << "===> Warning: Cannot solve inner system!\n";
        if (output_ != 0)
            *output_ << "SolverAsPreCL: iterations: " << solver_.GetIter()
                     << "\trelative residual: " << solver_.GetResid() << std::endl;
        AddIter( solver_.GetIter());
    }
    void Apply(const MatrixCL& A,   VectorCL& x, const VectorCL& b) const { Apply<>( A, x, b); }
    void Apply(const MLMatrixCL& A, VectorCL& x, const VectorCL& b) const { Apply<>( A, x, b); }
};

#ifdef _PAR
//***************************************************************************
// implemented parallel methods for preconditioning
//***************************************************************************
template <typename Mat, typename Vec>
void Jacobi(const Mat& A, const Vec& Diag, Vec& x, const Vec& b, const double omega, const bool HasOmega);

template <typename Mat, typename Vec>
void Jacobi0(const Mat&, const Vec& Diag, Vec& x, const Vec& b, const double omega, const bool HasOmega);

template <typename Mat, typename Vec>
void SSOR0(const Mat& A, const Vec& Diag, Vec& x, const Vec& b, const double omega, const bool HasOmega);


//***************************************************************************
// preconditioning-classes
//***************************************************************************

// ***************************************************************************
/// \brief Base class for parallel preconditioning classes
// ***************************************************************************
class ParPreconditioningBaseCL
/** This class can be used as a base class for parallel preconditioners.
    Therefore this class stores a reference to an IdxDescCL to get a reference
    to the ExchangeCL. Since, the index within the problem classes are lying all
    the time at the same position in memory, we store a reference to the
    IdxDescCL instead of storing a reference to the ExchanegCL.
    This class can create an accumulated diagonal out of a matrix by its own.
    This diagonal is only created, if the corresponding matrix has been
    changed.*/
{
  protected:
    const IdxDescCL& idx_;                ///< index for accessing the ExchangeCL
    VectorCL   diag_;                     ///< accumulated diagonal of corresponding matrix
    size_t     mat_version_;              ///< version of the corresponding matrix
    bool       check_mat_version_;        ///< check matrix version before apply preconditioning (only if DebugNumericC is set)

  public:
    ParPreconditioningBaseCL(const IdxDescCL& idx)
      : idx_(idx), diag_(0), mat_version_(0), check_mat_version_(true) {}
    // default c-tor and de-tor

    /// \brief Set accumulated diagonal of a matrix, that is needed by most of the preconditioners
    template<typename Mat>
    void SetDiag(const Mat& A)
    {
        Comment( (A.Version() == mat_version_ ? "SetDiag, Reusing OLD diagonal\n" : "SetDiag: Creating NEW diagonal\n"), DebugNumericC);
        // Just set diagonal, if the diagonal has changed
        if (A.Version() == mat_version_)
            return;
        if (diag_.size() != A.num_rows())
            diag_.resize(A.num_rows());
        // accumulate diagonal of the matrix
        diag_= A.GetDiag();
        GetEx().Accumulate(diag_);
        // remeber version of the matrix
        mat_version_= A.Version();
    }

    /// \brief Set accumulated diagonal of a matry by the given accumulated diagonal
    void SetDiag(const VectorCL& diag)
    /** \pre diag must be accumulated*/
        { if (diag.size()!=diag_.size()) diag_.resize(diag.size()); diag_=diag; }
    /// \brief Get constant reference on accumulated diagonal of corresponding matrix
    const VectorCL& GetDiag() const { return diag_; }
    /// \brief Get reference on accumulated diagonal of corresponding matrix
    VectorCL&       GetDiag()       { return diag_; }
    /// \brief Get constant reference on exchange class
    const ExchangeCL& GetEx() const { return idx_.GetEx(); }
    /// \brief Check matrix version before apply preconditioner in Debug-Mode (DebugNumericC) default
    void CheckMatVersion() { check_mat_version_=true; }
    /// \brief Don't check matrix version before apply preconditioner
    void DoNotCheckMatVersion() { check_mat_version_=false; }
};

// ***************************************************************************
/// \brief Class for performing a preconditioning step with the identity matrix
// ***************************************************************************
class ParDummyPcCL : public ParPreconditioningBaseCL
{
  private:
    typedef ParPreconditioningBaseCL base_;

  public:
    ParDummyPcCL(const IdxDescCL& idx) : base_(idx) {}

    /// \brief Check if return preconditioned vectors are accumulated after calling Apply
    bool RetAcc()   const { return false; }
    /// \brief Check if the diagonal of the matrix is needed
    bool NeedDiag() const { return false; }
    /// \name Set diagonal of the matrix for consistency
    //@{
    void SetDiag(const VectorCL&) {}        // just for consistency
    template<typename Mat>
    void SetDiag(const Mat&) {}             // just for consistency
    //@}

    /// \brief Apply preconditioner: x <- b
    template <typename Mat, typename Vec>
    void Apply(const Mat& /*A*/, Vec &x, const Vec& b) const
    {
        x=b;
    }
    /// \brief Apply preconditioner of A^T: x <- b
    template <typename Mat, typename Vec>
    void transp_Apply(const Mat&, Vec &x, const Vec& b) const
    {
        x=b;
    }

};

// ***************************************************************************
/// \brief Class for performing one Step of the Jacobi-Iteration
// ***************************************************************************
class ParJacCL : public ParPreconditioningBaseCL
{
  protected:
    typedef ParPreconditioningBaseCL base_;
    using base_::diag_;
    using base_::mat_version_;
    using base_::check_mat_version_;
    using base_::idx_;

  protected:
    double  omega_;                                 // overrelaxion-parameter

  public:
    ParJacCL (const IdxDescCL& idx, double omega=1) : base_(idx), omega_(omega) {}

    /// \brief Check if return preconditioned vectors are accumulated after calling Apply
    bool RetAcc() const   { return true; }
    /// \brief Check if the diagonal of the matrix is needed
    bool NeedDiag() const { return true; }
    /// \brief Get overrelaxation parameter
    double GetOmega() const {return omega_;}

    /// \brief Apply preconditioner: one step of the Jacobi-iteration
    template <typename Mat, typename Vec>
    void Apply(const Mat& A, Vec &x, const Vec& b) const
    {
        Assert(mat_version_==A.Version() || !check_mat_version_,
               DROPSErrCL("ParJacCL::Apply: Diagonal of actual matrix has not been set"),
               DebugNumericC);
        Jacobi(A, diag_, x, b, omega_, std::fabs(omega_-1.)>DoubleEpsC, idx_.GetEx());
    }
};

// ********************************************************************************
/// \brief Class for performing one Step of the Jacobi-Iteration with startvector 0
// ********************************************************************************
class ParJac0CL : public ParJacCL
{
  protected:
    typedef ParJacCL base_;
    using base_::mat_version_;
    using base_::diag_;
    using base_::omega_;
    using base_::check_mat_version_;

  public:

    ParJac0CL (const IdxDescCL& idx, double omega=1) : base_(idx, omega) {}

    /// \brief Check if return preconditioned vectors are accumulated after calling Apply
    bool RetAcc() const   { return false; }

    /// \brief Apply preconditioner: one step of the Jacobi-iteration with start vector 0
    template <typename Mat, typename Vec>
    void Apply(const Mat& A, Vec &x, const Vec& b) const
    {
        Assert(mat_version_==A.Version() || !check_mat_version_,
               DROPSErrCL("ParJac0CL::Apply: Diagonal of actual matrix has not been set"),
               DebugNumericC);
        Jacobi0(A, diag_, x, b, omega_, std::fabs(omega_-1.)>DoubleEpsC);
    }

    /// \brief Apply preconditioner of A^T: one step of the Jacobi-iteration with start vector 0
    template <typename Mat, typename Vec>
    void transp_Apply(const Mat& A, Vec &x, const Vec& b) const
    {
        Assert(mat_version_==A.Version()  || !check_mat_version_,
               DROPSErrCL("ParJac0CL::transp_Apply: Diagonal of actual matrix has not been set"),
               DebugNumericC);
        Jacobi0(A, diag_, x, b, omega_, std::fabs(omega_-1.)>DoubleEpsC);
    }
};

// ********************************************************************************
/// \brief Class for performing one step of the Jacobi-Iteration on a matrix A*A^T
///   with startvector 0
// ********************************************************************************
class ParJacNEG0CL : public ParJac0CL
{
  protected:
    typedef ParJac0CL base_;        ///< base class
    ExchangeMatrixCL  exMat_;       ///< handling of accumulating matrix-entries

    /// \brief Determine the diagonal of A*A^T
    inline void MySetDiag(const MatrixCL& A, const ExchangeCL& RowEx, const ExchangeCL& ColEx);

  public:
    /// \brief Constructor
    ParJacNEG0CL (const IdxDescCL& idx, double omega=1) : base_(idx, omega) {}

    /// \name Compute diagonal of matrix AA^T
    //@{
    /// \brief Set known diagonal of A*A^T as diagonal
    void SetDiag(const VectorCL& d) { base_::diag_.resize(d.size()); base_::diag_=d; }

    /// \brief Determine diagonal of a matrix A is not implemented for all MatTs
    template<typename MatT>
    void SetDiag(const MatT&){
        throw DROPSErrCL("ParJacNEG0CL::SetDiag: Not defined for that matrix type.");
    }
    //@}
};

inline void ParJacNEG0CL::MySetDiag(const MatrixCL& A, const ExchangeCL& RowEx, const ExchangeCL& ColEx)
/// This function determines the diagonal of the matrix A*A^T and stores this diagonal
/// in the base class. Therefore, the matrix A has to be accumulated. This is done by
/// the class MatrixExchangeCL.
///
/// Conservatively, we assume that the pattern of the matrix has changed, so the MatrixExchangeCL
/// is created each time, this function is called.
///
/// \todo(par) Determine when to create the communication-pattern more precisely
/// \param A     Compute diagonal of the matrix A*A^T
/// \param RowEx ExchangeCL according to the RowIdx of matrix A
/// \param ColEx ExchangeCL according to the ColIdx of matrix A
{
    // Accumulate matrix
    exMat_.BuildCommPattern(A, RowEx, ColEx);
    MatrixCL Aacc(exMat_.Accumulate(A));
    // Determine diagonal of AA^T
    base_::diag_.resize(Aacc.num_rows());
    for (size_t i = 0; i < Aacc.num_rows(); ++i)
        for (size_t nz = Aacc.row_beg(i); nz < Aacc.row_beg(i + 1); ++nz)
            base_::diag_[i] += Aacc.val(nz) * A.val(nz);
    RowEx.Accumulate(base_::diag_);
    mat_version_= A.Version();
}

/// \brief (Specialization) SetDiag for CompositeMatrixBaseCL<MatrixCL, MatrixCL>
template <>
inline void ParJacNEG0CL::SetDiag(const CompositeMatrixBaseCL<MatrixCL, MatrixCL>& BBT)
{
    MySetDiag(*BBT.GetBlock1(), BBT.GetEx1(), BBT.GetEx0());
}


// ********************************************************************************
/// \brief Class for performing one Step of the Jacobi-Iteration with startvector 0 and own matrix
// ********************************************************************************
class ParJac0OwnMatCL : public ParJac0CL
{
  private:
    typedef ParJac0CL base_;
    using base_::mat_version_;
    using base_::diag_;
    using base_::omega_;
    using base_::check_mat_version_;

  public:
    ParJac0OwnMatCL (const IdxDescCL& idx, double omega=1) : base_(idx, omega) {}

    /// \brief Apply preconditioner: one step of the Jacobi-iteration with start vector 0
    template <typename Mat, typename Vec>
    void Apply(const Mat& A, Vec &x, const Vec& b) const
    {
        Jacobi0(A, diag_, x, b, omega_, std::fabs(omega_-1.)>DoubleEpsC);
    }

    /// \brief Apply preconditioner of A^T: one step of the Jacobi-iteration with start vector 0
    template <typename Mat, typename Vec>
    void transp_Apply(const Mat& A, Vec &x, const Vec& b) const
    {
        Jacobi0(A, diag_, x, b, omega_, std::fabs(omega_-1.)>DoubleEpsC);
    }

    /// \brief Check if the diagonal of the matrix is needed
    bool NeedDiag() const { return false; }

    /// \brief Set Diag makes nothing, because own matrix is used
    template<typename Mat>
    void SetDiag(const Mat&) {}

    /// \brief Set new diagonal
    template<typename Mat>
    void SetNewDiag(const Mat& M) { base_::SetDiag(M); }

    /// \brief Set new diagonal (given as an accumulated vector)
    void SetNewDiag(const VectorCL& v) { base_::SetDiag(v); }
};

// ********************************************************************************
/// \brief Class for performing an accumulation as preconditioning
// ********************************************************************************
class ParAccPcCL : public ParDummyPcCL
{
  private:
    typedef ParDummyPcCL base_;

  public:
    ParAccPcCL(const IdxDescCL& idx) : base_(idx) {}

    /// \brief Check if the diagonal of the matrix is needed
    bool NeedDiag() const { return false; }
    /// \brief Check if return preconditioned vectors are accumulated after calling Apply
    bool RetAcc() const { return true; }

    /// \brief Apply preconditioner: x <- Accumulate(b)
    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec &x, const Vec& b) const
    {
        x= base_::GetEx().GetAccumulate(b);
    }
};

// ********************************************************************************
/// \brief Class for performing one Step of the blocked inexact SSOR0 step
// ********************************************************************************
class ParSSOR0CL : public ParJacCL
{
  private:
    typedef ParJacCL base_;
    using base_::diag_;
    using base_::mat_version_;
    using base_::omega_;
    using base_::check_mat_version_;

  public:
    ParSSOR0CL (const IdxDescCL& idx, double omega=1) : base_(idx, omega) {}

    /// \brief Check if return preconditioned vectors are accumulated after calling Apply
    bool RetAcc() const { return false; }
    /// \brief Check if the diagonal of the matrix is needed
    bool NeedDiag() const { return true; }

    /// \brief Apply preconditioner: one step of the incomplete SSOR-iteration
    template <typename Mat, typename Vec>
    void Apply(const Mat& A, Vec &x, const Vec& b) const
    {
        Assert(mat_version_==A.Version() || !check_mat_version_,
               DROPSErrCL("ParJac0CL::Apply: Diagonal of actual matrix has not been set"),
               DebugNumericC);
        SSOR0(A, diag_, x, b, omega_, std::fabs(omega_-1.)>DoubleEpsC);
    }
};

// ********************************************************************************
/// \brief serial CG as preconditioner
// ********************************************************************************
class ParCGPreCL
{
  private:
    int iters_;
    double tol_;

  public:
    ParCGPreCL(int Steps, double tol) : iters_(Steps), tol_(tol) {}

    inline bool RetAcc() const {return true;}
    inline bool NeedDiag() const {return false;}
    void SetDiag(VectorCL */*diag*/) {}

    template <typename Mat, typename Vec>
    void Apply(const Mat& A, Vec &x, const Vec& b) const
    {
        int iter= iters_;
        double tol=tol_;
        x=0.;
        PCG(A, x, b, SSORPcCL(), iter, tol);
        if (ProcCL::MyRank()==0)
            std::cout << "  -- CG as Precond used " << iter <<" steps\n";
    }
};

// ********************************************************************************
/// \brief serial GMRES as preconditioner
// ********************************************************************************
class ParGMResPreCL
{
  private:
    int iters_;
    int restart_;
    double tol_;

 public:
    ParGMResPreCL(int Restart, int Steps, double tol) : iters_(Steps), restart_(Restart), tol_(tol) {}

    inline bool RetAcc() const {return true;}
    inline bool NeedDiag() const {return false;}
    void SetDiag(VectorCL */*diag*/) {}

    template <typename Mat, typename Vec>
    void Apply(const Mat& A, Vec &x, const Vec& b) const
    {
        int iter= iters_;
        double tol=tol_;
        GMRES(A, x, b, MultiSSORPcCL(1,1), restart_, iter, tol);
        if (ProcCL::MyRank()==0)
            std::cout << "  -- GMRES as Precond used " << iter <<" steps\n";
    }
};


// ********************************************************************************
///\brief Apply a diagonal-matrix given as a vector.
// ********************************************************************************
class ParDiagPcCL : public ParPreconditioningBaseCL
{
  private:
    typedef ParPreconditioningBaseCL base_;
    const VectorCL& D_;

  public:
    ParDiagPcCL (const IdxDescCL& idx, const VectorCL& D) : base_(idx), D_( D) {}

    template <typename Mat, typename Vec>
    void Apply (const Mat&, Vec& x, const Vec& b) const
    {
        Assert( D_.size()==b.size(), DROPSErrCL("DiagPcCL: incompatible dimensions"), DebugNumericC);
        x= D_*b;
    }

    inline bool RetAcc() const {return false;}
    inline bool NeedDiag() const {return false;}
    /// \brief Set Diag makes nothing, because own matrix is used
    template<typename Mat>
    void SetDiag(const Mat&) {}
};

// ********************************************************************************
///\brief Multilevel Smoother
// ********************************************************************************

template <class SmootherT>
class MLSmootherCL : public MLDataCL<SmootherT> {
  public:
    MLSmootherCL( const MLIdxDescCL& idx, const double omega = 1.0) {
        for ( MLIdxDescCL::const_iterator it = idx.begin(); it != idx.end(); ++it){
            this->push_back( SmootherT(*it, omega));
        }
    }
    void SetDiag( const MLMatrixCL& A) {
        Assert( A.size() == this->size(), DROPSErrCL ("MLSmootherCL::SetDiag: dimensions do not fit\n"), DebugNumericC);
        MLMatrixCL::const_iterator Ait = A.begin();
        for ( typename MLSmootherCL::iterator Sit = this->begin(); Sit != this->end(); ++Sit, ++Ait)
            Sit->SetDiag(*Ait);
    }
};

//***************************************************************************
// Implementations of the methods
//***************************************************************************

/// \brief One step of a Jacobi-iteration
template <typename Mat, typename Vec, typename ExT>
void Jacobi(const Mat& A, const Vec& Diag, Vec& x, const Vec& b, const double omega, const bool HasOmega, const ExT& ex)
        /// \param[in]     A        local distributed coefficients-matrix of the linear equation system
        /// \param[in]     Diag     accumulated form of the diagonalelements of A
        /// \param[in,out] x        startvector (accumulated) and result (distributed)
        /// \param[in]     b        rhs (distributed form)
        /// \param[in]     omega    overrelaxion-parameter
        /// \param[in]     HasOmega flag, if this routine takes care of an overrelaxion-parameter
{
    const size_t n= A.num_rows();
    Vec          y(x.size());

    for (size_t i=0, nz=0; i<n; ++i)
    {
        double sum= b[i];
        for (const size_t end= A.row_beg(i+1); nz<end; ++nz)
            if (A.col_ind(nz) != i)
                sum-= A.val(nz)*x[A.col_ind(nz)];

        y[i] = sum/Diag[i];
    }
    ex.Accumulate(y);

    if (HasOmega)
        y = (1.-omega)*x+omega*y;

    std::swap(x,y);
}

/// \brief One Step of a Jacobi-iteration with startvector 0
template <typename Mat, typename Vec>
void Jacobi0(const Mat&, const Vec& Diag, Vec& x, const Vec& b, const double omega, const bool HasOmega)
        /// \param[in]     Diag     accumulated form of the diagonalelements of A
        /// \param[out]    x        result (distributed form)
        /// \param[in]     b        rhs (distributed form)
        /// \param[in]     omega    overrelaxion-parameter
        /// \param[in]     HasOmega flag, if this routine takes care of an overrelaxion-parameter

{
    const size_t n = x.size();
    for (size_t i=0; i<n; ++i)
    {
        if (HasOmega)
            x[i] = omega * b[i] / Diag[i];
        else
            x[i] = b[i] / Diag[i];
    }
}

// One step of the Symmetric-Gauss-Seidel/SSOR method with start vector 0
template <typename Mat, typename Vec>
void SSOR0(const Mat& A, const Vec& Diag, Vec& x, const Vec& b, const double omega, const bool HasOmega)
{
    const size_t n= A.num_rows();

    for (size_t i=0; i<n; ++i)
    {
        double sum= b[i];
        size_t j= A.row_beg(i);
        for ( ; A.col_ind(j) < i; ++j)
            sum-= A.val(j)*x[A.col_ind(j)];
        if (HasOmega)
            x[i]= omega*sum/Diag[i];
        else
            x[i]= sum/Diag[i];
    }

    for (size_t i=n; i>0; )
    {
        --i;
        double sum= 0;
        size_t j= A.row_beg(i+1)-1;
        for ( ; A.col_ind(j) > i; --j)
            sum-= A.val(j)*x[A.col_ind(j)];
        if (HasOmega)
            x[i]= (2.-omega)*x[i]+omega*sum/Diag[i];
        else
            x[i]+= sum/Diag[i];
    }
}
#endif

} // end of namespace DROPS

#endif /* PRECOND_H_ */
