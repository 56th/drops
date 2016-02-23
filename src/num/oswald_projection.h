/// \file oswald_projection.h
/// \brief Oswald-projection (arithmetic mean in each dof).
/// \author Joerg Grande

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
 * Copyright 2016 Joerg Grande, Aachen, Germany
*/

#ifndef DROPS_OSWALD_PROJECTION_H
#define DROPS_OSWALD_PROJECTION_H

#include "misc/problem.h"
#include "num/accumulator.h"
#include "num/fe.h"

namespace DROPS
{

///\brief Accumulator to compute the Oswald-projection (arithmetic mean in each dof) of the data in LocalP2T.
/// LocalP2T is not a short-hand for LocalP2CL<...>. It must provide:
///    * a typedef value_type,
///    * the number of components of value_type (static int num_components),
///    * an operator[](size_t i) returning the local value in the P2-dof i,
///    * the function set_tetra(const TetraCL&) to setup the local data on a tetra,
///    * the predicate invalid_p(size_t i).
template <typename LocalP2T>
class OswaldProjectionP2AccuCL : public TetraAccumulatorCL
{
  private:
    LocalP2T loc_;
    std::valarray<double>* n_;
    VecDescCL& avg_;

    LocalNumbP2CL numg;

    void set_n (std::valarray<double>* n) { n_= n; } // The clones must refer to the n_ of thread 0.

  public:
    OswaldProjectionP2AccuCL (LocalP2T loc, VecDescCL& avg)
        : loc_( loc), n_( 0), avg_( avg) {}

    virtual void begin_accumulation   () {
        n_= new std::valarray<double>( avg_.Data.size()/loc_.num_components);
    }
    virtual void finalize_accumulation() {
        delete n_;
    }

    virtual void visit (const TetraCL& t) {
        loc_.set_tetra( &t);
        numg.assign_indices_only( t, *avg_.RowIdx);
        for (Uint i= 0; i < 10; ++i) {
            if (!numg.WithUnknowns( i) || loc_.invalid_p (i))
                continue;
            const IdxT dof= numg.num[i];
            double& n= n_[0][dof/loc_.num_components]; // This assumes that the local gradient is in the components dof/3..dof/3 + 2.
            n+= 1.;
            typedef typename LocalP2T::value_type value_type;
            const value_type& oldavg= DoFHelperCL<value_type, VectorCL>::get( avg_.Data, dof);
            DoFHelperCL<value_type, VectorCL>::set( avg_.Data, dof, ((n - 1.)/n)*oldavg + (1./n)*loc_[i]);
        }
    }

    virtual TetraAccumulatorCL* clone (int /*clone_id*/) {
        OswaldProjectionP2AccuCL* p= new OswaldProjectionP2AccuCL( loc_, avg_);
        p->set_n( n_);
        return p;
    }
};

} // end of namespace DROPS

#endif