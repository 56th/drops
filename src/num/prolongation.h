/// \file prolongation.h
/// \brief prolongation and restriction routines for (parallel) multigrid solver
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef PROLONGATION_H_
#define PROLONGATION_H_

#include "geom/simplex.h"
#include "geom/multigrid.h"
#include "misc/problem.h"
#ifdef _PAR
#include "parallel/exchange.h"
#endif


namespace DROPS{

template <class ValueT>
class ProlongationCL
{
  public:
    typedef ValueT value_type;

  private:
    IdxDescCL* coarse_;
    IdxDescCL* fine_;

    const MultiGridCL& mg_;

    VectorCL scale_;                   ///< weights for prolongation/restriction

#ifdef _PAR
    typedef LocalP2CL<value_type> ChildDataT;
    typedef DROPS_STD_UNORDERED_MAP<const TetraCL*, ChildDataT > DataT;

    mutable DataT data_;

    class HandlerParentDataCL;
    void Restrict(const TetraCL* sit, const VecDescCL& fine, const bool doinjection=false) const;
    void Prolong( const TetraCL* sit, const LocalP2CL<ValueT>& coarse, VecDescCL& fine) const;
#endif

    std::vector<IdxT> CollectChildUnknownsP2(const TetraCL& t, const Uint f_idx) const;

    void BuildP2ProlongationMatrix(const IdxDescCL& coarse, const IdxDescCL& fine);
    void BuildP1ProlongationMatrix(const IdxDescCL& coarse, const IdxDescCL& fine);

    MatrixCL prolongation_;
    MatrixCL injection_;

  public:
    ProlongationCL(const MultiGridCL& mg) : coarse_(0), fine_(0), mg_(mg){}

    void Create(IdxDescCL* coarse, IdxDescCL* fine);

    VectorCL operator* (const VectorCL& vec) const;

    VectorCL mytransp_mul(const VectorCL& vec, const bool doinjection = false) const;

    VectorCL restrict_vec(const VectorCL& vec) const { return mytransp_mul(vec, true);}
};


template <class ValueT>
VectorCL transp_mul(const ProlongationCL<ValueT>& p, const VectorCL& v){
    return p.mytransp_mul(v);
}

template<class ValueT>
void SetupProlongationMatrix(const MultiGridCL& mg, MLDataCL<ProlongationCL<ValueT> >& P, MLIdxDescCL* ColIdx, MLIdxDescCL* RowIdx)
{
    typedef MLDataCL<ProlongationCL<ValueT> > ProlongT;
    P.clear();
    for (size_t i=0; i< ColIdx->size(); ++i)
        P.push_back(ProlongationCL<ValueT>(mg));
    typename ProlongT::iterator itProlong = ++P.begin();
    MLIdxDescCL::iterator itcIdx = ColIdx->begin();
    MLIdxDescCL::iterator itfIdx = ++RowIdx->begin();
    for (size_t lvl=1; lvl < P.size(); ++lvl, ++itProlong, ++itcIdx, ++itfIdx)
        itProlong->Create(&*itcIdx, &*itfIdx);
}

/// \brief Observes the MultiGridCL-changes by AdapTriangCL to repair the prolongation for velocity.
template<class ValueT>
class UpdateProlongationCL : public MGObserverCL
{
  private:
    const MultiGridCL& MG_;
    MLDataCL<ProlongationCL<ValueT> >  *P_;
    MLIdxDescCL *ColIdx_, *RowIdx_;

  public:
    UpdateProlongationCL( const MultiGridCL& MG, MLDataCL<ProlongationCL<ValueT> >* P, MLIdxDescCL* ColIdx, MLIdxDescCL* RowIdx)
        : MG_( MG), P_( P), ColIdx_( ColIdx), RowIdx_( RowIdx) { post_refine_sequence (); }

    void pre_refine  () {}
    void post_refine () {}

    void pre_refine_sequence  () {}
    void post_refine_sequence () {
        if (P_ != 0) {
            P_->clear();
            SetupProlongationMatrix( MG_, *P_, ColIdx_, RowIdx_);
        }
    }
    const IdxDescCL* GetIdxDesc() const { return (const IdxDescCL*)0; }
#ifdef _PAR
    const VectorCL*  GetVector()  const { return 0; }
    void swap( IdxDescCL&, VectorCL&) {}
#endif
};

} //end of namespace DROPS

#include "num/prolongation.tpp"

#endif /* PROLONGATION_H_ */
