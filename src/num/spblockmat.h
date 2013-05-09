/// \file spblockmat.h
/// \brief block and composed matrices
/// \author LNM RWTH Aachen: Sven Gross, Joerg Peters; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_SPBLOCKMAT_H
#define DROPS_SPBLOCKMAT_H

#include "num/spmat.h"
#include "misc/container.h"
#include "misc/problem.h"

namespace DROPS{

//*****************************************************************************
//
/// \brief 2x2-Blockmatrix: ( A & B \\ C & D)
//
//*****************************************************************************
enum BlockMatrixOperationT { MUL, TRANSP_MUL };

template <class MatT>
class BlockMatrixBaseCL
{
  public:
    typedef BlockMatrixOperationT OperationT;

  private:
    const MatT* block_[4];
    OperationT operation_[4];

    bool block_num_rows(size_t b, size_t& nr) const;
    bool block_num_cols(size_t b, size_t& nc) const;

  public:
    BlockMatrixBaseCL( const MatT* A, OperationT Aop, const MatT* B, OperationT Bop,
        const MatT* C, OperationT Cop, const MatT* D= 0, OperationT Dop= MUL);

    /// \brief number of rows of a block
    size_t num_rows(size_t) const;
    /// \brief number of columns of a block
    size_t num_cols(size_t) const;
    /// \brief number of overall rows
    size_t num_rows() const { return this->num_rows( 0) + this->num_rows( 1); }
    /// \brief number of overall columns
    size_t num_cols() const { return this->num_cols( 0) + this->num_cols( 1); }

    /// \name Getters and Setters
    //@{
    const MatT* GetBlock( size_t b) const { return block_[b]; }
    void        SetBlock( size_t b, const MatT* m) { block_[b]= m; }

    OperationT GetOperation( size_t b) const { return operation_[b]; }
    OperationT GetTransposeOperation( size_t b) const {
        return operation_[b] == MUL ? TRANSP_MUL : MUL;
    }
    BlockMatrixBaseCL GetTranspose() const;
    VectorBaseCL<typename MatT::value_type> GetDiag() const;
    //@}
};

template <class MatT>
BlockMatrixBaseCL<MatT>::BlockMatrixBaseCL( const MatT* A, OperationT Aop,
    const MatT* B, OperationT Bop, const MatT* C, OperationT Cop,
    const MatT* D, OperationT Dop)
{
    block_[0]= A; operation_[0]= Aop;
    block_[1]= B; operation_[1]= Bop;
    block_[2]= C; operation_[2]= Cop;
    block_[3]= D; operation_[3]= Dop;
}

template <class MatT>
  bool BlockMatrixBaseCL<MatT>::block_num_rows(size_t b, size_t& nr) const
{
    if (block_[b] == 0) return false;
    switch (operation_[b]) {
      case MUL:        nr= block_[b]->num_rows(); return true;
      case TRANSP_MUL: nr= block_[b]->num_cols(); return true;
      default:
        Comment("BlockMatrixBaseCL::block_num_rows: No such operation.\n", DebugNumericC);
        return false;
    }
}

template <class MatT>
  bool BlockMatrixBaseCL<MatT>::block_num_cols(size_t b, size_t& nc) const
{
    if (block_[b] == 0) return false;
    switch (operation_[b]) {
      case MUL:        nc= block_[b]->num_cols(); return true;
      case TRANSP_MUL: nc= block_[b]->num_rows(); return true;
      default:
        Comment("BlockMatrixBaseCL::block_num_cols: No such operation.\n", DebugNumericC);
        return false;
    }
}

template <class MatT>
  size_t BlockMatrixBaseCL<MatT>::num_rows(size_t block_row) const
{
    size_t ret= 0;
    __UNUSED__ bool block_found;
    switch (block_row) {
      case 0:
        block_found= block_num_rows( 0, ret) || block_num_rows( 1, ret);
        break;
      case 1:
        block_found= block_num_rows( 2, ret) || block_num_rows( 3, ret);
        break;
      default:
        Comment("BlockMatrixBaseCL::num_rows: No such block_row.\n", DebugNumericC);
        return 0;
    }
    Assert( block_found, "BlockMatrixBaseCL::num_rows: All pointers are 0.\n", DebugNumericC);
    return ret;
}

template <class MatT>
  size_t BlockMatrixBaseCL<MatT>::num_cols(size_t block_col) const
{
    size_t ret= 0;
    __UNUSED__ bool block_found;
    switch (block_col) {
      case 0:
        block_found= block_num_cols( 0, ret) || block_num_cols( 2, ret);
        break;
      case 1:
        block_found= block_num_cols( 1, ret) || block_num_cols( 3, ret);
        break;
      default:
        Comment("BlockMatrixBaseCL::num_cols: No such block_col.\n", DebugNumericC);
        return 0;
    }
    Assert( block_found, "BlockMatrixBaseCL::num_cols: All pointers are 0.\n", DebugNumericC);
    return ret;
}

template <class MatT>
  BlockMatrixBaseCL<MatT> BlockMatrixBaseCL<MatT>::GetTranspose() const
{
    return BlockMatrixBaseCL<MatT>( block_[0], GetTransposeOperation( 0),
        block_[2], GetTransposeOperation( 2),
        block_[1], GetTransposeOperation( 1),
        block_[3], GetTransposeOperation( 3));
}

template <typename MatT>
  VectorBaseCL<typename MatT::value_type> BlockMatrixBaseCL<MatT>::GetDiag() const
{
    const size_t n=num_rows();
    Assert(n==num_cols(), "SparseMatBaseCL::GetDiag: no square Matrix", DebugParallelC);
    VectorBaseCL<typename MatT::value_type> diag(n);

    if (block_[0])
        diag[std::slice(0,num_rows(0),1)] = block_[0]->GetDiag();
    else
        diag[std::slice(0,num_rows(0),1)] = 1.;

    if (block_[3])
        diag[std::slice(num_rows(0),num_rows(1),1)] = block_[3]->GetDiag();
    else
        diag[std::slice(num_rows(0),num_rows(1),1)] = 1.;

    return diag;
}

template <typename Mat, typename Vec>
  Vec operator* (const BlockMatrixBaseCL<Mat>& A, const Vec& x)
{
    Vec x0( x[std::slice( 0, A.num_cols( 0), 1)]),
        x1( x[std::slice( A.num_cols( 0), A.num_cols( 1), 1)]);
    Vec r0( A.num_rows( 0)),
        r1( A.num_rows( 1));
    const Mat* mat;

    if ( (mat= A.GetBlock( 0)) != 0) {
        switch( A.GetOperation( 0)) {
          case MUL:
            r0= (*mat)*x0; break;
          case TRANSP_MUL:
            r0= transp_mul( *mat, x0); break;
        }
    }
    if ( (mat= A.GetBlock( 1)) != 0) {
        switch( A.GetOperation( 1)) {
          case MUL:
            r0+= (*mat)*x1; break;
          case TRANSP_MUL:
            r0+= transp_mul( *mat, x1); break;
        }
    }
    if ( (mat= A.GetBlock( 2)) != 0) {
        switch( A.GetOperation( 2)) {
          case MUL:
            r1= (*mat)*x0; break;
          case TRANSP_MUL:
            r1= transp_mul( *mat, x0); break;
        }
    }
    if ( (mat= A.GetBlock( 3)) != 0) {
        switch( A.GetOperation( 3)) {
          case MUL:
            r1+= (*mat)*x1; break;
          case TRANSP_MUL:
            r1+= transp_mul( *mat, x1); break;
        }
    }
    Vec ret( A.num_rows());
    ret[std::slice( 0, A.num_rows( 0), 1)]= r0;
    ret[std::slice( A.num_rows( 0), A.num_rows( 1), 1)]= r1;
    return ret;
}

template <typename Mat, typename Vec>
  Vec transp_mul(const BlockMatrixBaseCL<Mat>& A, const Vec& x)
{
    return A.GetTranspose()*x;
}


//*****************************************************************************
//
/// \brief adapter: 2x2-Blockmatrix: ( A & B \\ C & D) to Matrix
//
//*****************************************************************************
template <typename Mat>
Mat BuildMatrix( BlockMatrixBaseCL<Mat>& A)
{
    size_t rows = A.num_rows(), cols = A.num_cols();
    Mat ret;
    SparseMatBuilderCL<typename Mat::value_type> B(&ret, rows, cols);

    for (int block = 0; block <4; ++block) {
        if (!A.GetBlock(block)) continue;
        size_t row_shift = 0, col_shift = 0;
        switch (block) {
            case 0 : { row_shift = 0;             col_shift = 0;             } break;
            case 1 : { row_shift = 0;             col_shift = A.num_cols(0); } break;
            case 2 : { row_shift = A.num_rows(0); col_shift = 0;             } break;
            case 3 : { row_shift = A.num_rows(0); col_shift = A.num_cols(0); } break;
            default: throw DROPSErrCL("wrong block number");
        }
        for (size_t row = 0; row < A.GetBlock(block)->num_rows(); ++row)
            for (size_t nz = A.GetBlock(block)->row_beg(row); nz < A.GetBlock(block)->row_beg(row+1); ++ nz)
                if (A.GetOperation(block) == MUL)
                    B(row+row_shift, A.GetBlock(block)->col_ind(nz)+col_shift) = A.GetBlock(block)->val(nz);
                else
                    B(A.GetBlock(block)->col_ind(nz)+row_shift, row+col_shift) = A.GetBlock(block)->val(nz);
    }

    B.Build();
    return ret;
}

//*****************************************************************************
//
///  \brief Composition of 2 matrices
/** The matrices are applied as block1_*block0_*v */
//
//*****************************************************************************
template <class MatT0, class MatT1, class ExT>
class CompositeMatrixBaseCL
{
  public:
    typedef BlockMatrixOperationT OperationT;

  private:
    const ExT& ex0_;      // ExchangeCL for handling resulting vectors of block[0]*v
    const ExT& ex1_;      // ExchangeCL for handling resulting vectors of block_[1]*block[0]*v

    const MatT0* block0_;
    const MatT1* block1_;
    OperationT operation_[2];

  public:
    CompositeMatrixBaseCL( const MatT0* A, OperationT Aop, const ExT& ex0, const MatT1* B, OperationT Bop, const ExT& ex1)
      : ex0_(ex0), ex1_(ex1)
    {
        block0_= A; operation_[0]= Aop;
        block1_= B; operation_[1]= Bop;
    }

    /// \brief Get number of rows
    size_t num_rows() const;
    /// \brief Get number of columns
    size_t num_cols() const;
    /// \brief Get number of rows of the second matrix
    size_t intermediate_dim() const;
    /// \name Get a reference on the ExchangeCLs
    //@{
    const ExT& GetEx0() const { return ex0_; }
    const ExT& GetEx1() const { return ex1_; }
    VectorCL GetDiag() const { throw DROPSErrCL("CompositeMatrixBaseCL:GetDiag(): Not implemented"); return VectorCL(); };
    //@}

    /// \name Getters and Setters
    //@{
    const MatT0* GetBlock0 () const { return block0_; }
    const MatT1* GetBlock1 () const { return block1_; }
    void SetBlock0 (const MatT0* p) { block0_= p; }
    void SetBlock1 (const MatT1* p) { block1_= p; }

    OperationT GetOperation( size_t b) const { return operation_[b]; }
    OperationT GetTransposeOperation( size_t b) const {
        return operation_[b] == MUL ? TRANSP_MUL : MUL;
    }
    CompositeMatrixBaseCL<MatT1, MatT0, ExT> GetTranspose() const;
    size_t Version() const { return block0_->Version(); }
    //@}
};

template <class MatT0, class MatT1, class ExT>
  size_t CompositeMatrixBaseCL<MatT0, MatT1, ExT>::num_rows() const
{
    switch (operation_[1]) {
      case MUL:        return block1_->num_rows();
      case TRANSP_MUL: return block1_->num_cols();
      default:
        Comment("CompositeMatrixBaseCL::num_rows: No such operation.\n", DebugNumericC);
        return (size_t)-1;
    }
}

template <class MatT0, class MatT1, class ExT>
  size_t CompositeMatrixBaseCL<MatT0, MatT1, ExT>::num_cols() const
{
    switch (operation_[0]) {
      case MUL:        return block0_->num_cols();
      case TRANSP_MUL: return block0_->num_rows();
      default:
        Comment("CompositeMatrixBaseCL::num_cols: No such operation.\n", DebugNumericC);
        return (size_t)-1;
    }
}

template <class MatT0, class MatT1, class ExT>
  size_t CompositeMatrixBaseCL<MatT0, MatT1, ExT>::intermediate_dim() const
{
    switch (operation_[0]) {
      case MUL:        return block0_->num_rows();
      case TRANSP_MUL: return block0_->num_cols();
      default:
        Comment("CompositeMatrixBaseCL::intermediate_dim: No such operation.\n", DebugNumericC);
        return (size_t)-1;
    }
}


template <class MatT0, class MatT1, class ExT>
  CompositeMatrixBaseCL<MatT1, MatT0, ExT> CompositeMatrixBaseCL<MatT0, MatT1, ExT>::GetTranspose() const
{
    return CompositeMatrixBaseCL<MatT1,MatT0, ExT>( block1_, GetTransposeOperation( 1), ex1_,
                                               block0_, GetTransposeOperation( 0), ex0_);
}

template <typename _MatT0, typename _MatT1, class ExT, typename _VecEntry>
  VectorBaseCL<_VecEntry> operator*(const CompositeMatrixBaseCL<_MatT0, _MatT1, ExT>& A,
                                    const VectorBaseCL<_VecEntry>& x)
{
    VectorBaseCL<_VecEntry> tmp( A.intermediate_dim());
    switch( A.GetOperation( 0)) {
      case MUL:
        tmp= (*A.GetBlock0())*x; break;
      case TRANSP_MUL:
        tmp= transp_mul( *A.GetBlock0(), x); break;
    }
    A.GetEx0().Accumulate(tmp);
    VectorBaseCL<_VecEntry> ret( A.num_cols());
    switch( A.GetOperation( 1)) {
      case MUL:
        ret= (*A.GetBlock1())*tmp; break;
      case TRANSP_MUL:
        ret= transp_mul( *A.GetBlock1(), tmp); break;
    }
    return ret;
}

template <typename _MatT0, typename _MatT1, class ExT, typename _VecEntry>
VectorBaseCL<_VecEntry>
transp_mul(const CompositeMatrixBaseCL<_MatT0, _MatT1, ExT>& A,
    const VectorBaseCL<_VecEntry>& x)
{
    return A.GetTranspose()*x;
}

//=============================================================================
//  Typedefs
//=============================================================================

typedef BlockMatrixBaseCL<MatrixCL>               BlockMatrixCL;
typedef BlockMatrixBaseCL<MLMatrixCL>             MLBlockMatrixCL;
#ifndef _PAR
typedef CompositeMatrixBaseCL<MatrixCL, MatrixCL, DummyExchangeCL> CompositeMatrixCL;
#else
typedef CompositeMatrixBaseCL<MatrixCL, MatrixCL, ExchangeCL> CompositeMatrixCL;
#endif

} // end of namespace DROPS

#endif
