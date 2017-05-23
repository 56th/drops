/// \file oseenprecond.cpp
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

#include "num/oseenprecond.h"
#include "num/oseensolver.h"

namespace DROPS
{

#ifdef _PAR
void BDinvBTPreCL::Update(const ExchangeCL& vel_ex, const ExchangeCL& pr_ex) const
{
    std::cout << "BDinvBTPreCL::Update: old/new versions: " << Lversion_  << '/' << L_->Version()
        << '\t' << Bversion_ << '/' << B_->Version() << '\t' << Mversion_ << '/' << M_->Version()
        << '\t' << Mvelversion_ << '/' << Mvel_->Version() << '\n';
    delete Bs_;
    Bs_= new MatrixCL( *B_);
    Lversion_= L_->Version();
    Bversion_= B_->Version();
    Mversion_= M_->Version();
    Mvelversion_= Mvel_->Version();
    Cversion_ = C_->Version();

    Dvelinv_.resize( Mvel_->num_rows());
    if (lumped_)
        Dvelinv_ = 1.0/vel_ex.GetAccumulate(VectorCL(LumpInRows(*Mvel_)));
    else
        Dvelinv_ = 1.0/vel_ex.GetAccumulate(L_->GetDiag());

    // The lumped P2-mass-matrix has negative entries. Hence, CGNE cannot be used for solving. However, B Dvelinv_ B^T has only positive diagonal entries
    VectorCL Dprsqrt( std::sqrt( pr_ex.GetAccumulate(M_->GetDiag())));
    Dprsqrtinv_.resize( M_->num_rows());
    Dprsqrtinv_= 1.0/Dprsqrt;

    ScaleRows( *Bs_, Dprsqrtinv_);
    DSchurinv_.resize( Dprsqrt.size());

    VectorCL D(DSchurinv_.size());
    ExchangeMatrixCL exMat_;
    exMat_.BuildCommPattern(*Bs_, pr_ex, vel_ex);
    MatrixCL BSacc(exMat_.Accumulate(*Bs_));

    ScaleCols( BSacc, Dvelinv_);
    for (size_t i = 0; i < BSacc.num_rows(); ++i)
        for (size_t nz = BSacc.row_beg(i); nz < BSacc.row_beg(i + 1); ++nz)
            D[i] += BSacc.val(nz) * Bs_->val(nz);
    pr_idx_->GetEx().Accumulate(D);
    DSchurinv_ = 1.0/D;

    delete BDinvBT_;

    BDinvBT_= new AppSchurComplMatrixT( *L_, diagVelPc_, *Bs_, *C_, vel_ex);
//    SchurPc_.SetDiag(*BDinvBT_, D, pr_ex);
}
#endif

void BDinvBTPreCL::Update(const DummyExchangeCL& vel_ex, const DummyExchangeCL& pr_ex) const
{
    std::cout << "BDinvBTPreCL::Update: old/new versions: " << Lversion_  << '/' << L_->Version()
        << '\t' << Bversion_ << '/' << B_->Version() << '\t' << Mversion_ << '/' << M_->Version()
        << '\t' << Mvelversion_ << '/' << Mvel_->Version() << '\n';
    delete Bs_;
    delete Cs_;
    Bs_= new MatrixCL( *B_);
    Cs_= new MatrixCL( *C_);
    Lversion_= L_->Version();
    Bversion_= B_->Version();
    Mversion_= M_->Version();
    Mvelversion_= Mvel_->Version();
    Cversion_ = C_->Version();

    Dvelinv_.resize( Mvel_->num_rows());
    if (lumped_)
        Dvelinv_ = 1.0/vel_ex.GetAccumulate(VectorCL(LumpInRows(*Mvel_)));
    else
        Dvelinv_ = 1.0/vel_ex.GetAccumulate(L_->GetDiag());

    // The lumped P2-mass-matrix has negative entries. Hence, CGNE cannot be used for solving. However, B Dvelinv_ B^T has only positive diagonal entries
    VectorCL Dprsqrt( std::sqrt( pr_ex.GetAccumulate(M_->GetDiag())));
    Dprsqrtinv_.resize( M_->num_rows());
    Dprsqrtinv_= 1.0/Dprsqrt;

    ScaleRows( *Bs_, Dprsqrtinv_);
    ScaleRows( *Cs_, Dprsqrtinv_);
    ScaleCols( *Cs_, Dprsqrtinv_);
    DSchurinv_.resize( Dprsqrt.size());

    DSchurinv_= 1.0/( Bs_->GetSchurDiag(Dvelinv_) - Cs_->GetDiag() );
    delete SerBDinvBT_;

    SerBDinvBT_= new SerAppSchurComplMatrixT( *L_, diagVelPc_, *Bs_, *Cs_, vel_ex);
#ifndef _PAR
    if (regularize_ != 0.) {
        NEGSPcCL spc;
        Regularize( *Bs_, *pr_idx_, Dprsqrt, spc, regularize_, vel_ex, pr_ex);
        // Add a row to Dvelinv_ corresponding to the new column in Bs.
        VectorCL tmp= Dvelinv_;
        Dvelinv_.resize(Dvelinv_.size() + 1);
        Dvelinv_[std::slice( 0, tmp.size(), 1)]= tmp;
        Dvelinv_[tmp.size()]= 1.;
    }
#endif
//    SchurPc_.SetDiag(*SerBDinvBT_, Bs_->GetSchurDiag(Dvelinv_), pr_ex);
}

BDinvBTPreCL::~BDinvBTPreCL()
{
    delete Bs_;
    delete Cs_;
#ifdef _PAR
    delete BDinvBT_;
#endif
    delete SerBDinvBT_;
}


} // end of namespace DROPS
