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

#include "geom/subtriangulation.h"
#include "misc/problem.h"
#include "num/accumulator.h"
#include "num/fe.h"
#include "num/lattice-eval.h"

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
    std::valarray<double>* n_,
                         * n_invalid_;
    bool check_averaging_;
    VecDescCL& avg_;

    LocalNumbP2CL numg;

    const VecDescCL*   ls;      // a P2-level-set function
    const BndDataCL<>* lsetbnd; // boundary data for the level set function
    std::valarray<double> ls_loc;
    const PrincipalLatticeCL* lat;

    OswaldProjectionP2AccuCL& set_n (std::valarray<double>* n) { // The clones must refer to the n_ of thread 0.
        n_= n;
        return *this;
    }
    OswaldProjectionP2AccuCL& set_n_invalid (std::valarray<double>* n) { // The clones must refer to the n_invalid of thread 0.
        n_invalid_= n;
        return *this;
    }

  public:
    OswaldProjectionP2AccuCL (LocalP2T loc, VecDescCL& avg)
        : loc_( loc), n_( 0), n_invalid_( 0), check_averaging_( false), avg_( avg), ls( 0), lsetbnd( 0), lat( 0) {}

    OswaldProjectionP2AccuCL& set_check_averaging (bool b= true) {
        check_averaging_= b;
        return *this;
    }

    OswaldProjectionP2AccuCL& set_level_set_function (const VecDescCL* lsarg, const BndDataCL<>* lsetbndarg, const PrincipalLatticeCL* latarg) {
        ls= lsarg;
        lsetbnd= lsetbndarg;
        lat= latarg;
        ls_loc.resize ( lat ? lat->vertex_size () : 0);
        return *this;
    }

    virtual void begin_accumulation   () {
        n_= new std::valarray<double>( avg_.Data.size()/loc_.num_components);
        if (check_averaging_)
            n_invalid_= new std::valarray<double>( avg_.Data.size()/loc_.num_components);
    }
    virtual void finalize_accumulation() {
        loc_.finalize_accumulation ();
        if (check_averaging_)
            for (size_t i= 0; i < n_->size (); ++i)
                if (n_[0][i] == 0 && n_invalid_[0][i] > 0)
                    std::cerr << "OswaldProjectionP2AccuCL::finalize_accumulation: No local value for " << i << "; invalid_p: " << n_invalid_[0][i] << ".\n";
        delete n_;
        delete n_invalid_;
    }

    virtual void visit (const TetraCL& t) {
        if (ls != 0) {
            LocalP2CL<> locp2_ls( t, *ls, *lsetbnd);
            evaluate_on_vertexes( locp2_ls, *lat, Addr( ls_loc));
            if (equal_signs( ls_loc))
                return;
        }
        loc_.set_tetra( &t);
        numg.assign_indices_only( t, *avg_.RowIdx);
        for (Uint i= 0; i < 10; ++i) {
            if (!numg.WithUnknowns( i))
                continue;
            const IdxT dof= numg.num[i];
        if (loc_.invalid_p (i)) {
            if (check_averaging_)
                ++n_invalid_[0][dof/loc_.num_components];
            continue;
        }
            double& n= n_[0][dof/loc_.num_components]; // This assumes that the local gradient is in the components dof/3..dof/3 + 2.
            n+= 1.;
            typedef typename LocalP2T::value_type value_type;
            const value_type& oldavg= DoFHelperCL<value_type, VectorCL>::get( avg_.Data, dof);
            DoFHelperCL<value_type, VectorCL>::set( avg_.Data, dof, ((n - 1.)/n)*oldavg + (1./n)*loc_[i]);
        }
    }

    virtual TetraAccumulatorCL* clone (int /*clone_id*/) {
        OswaldProjectionP2AccuCL* p= new OswaldProjectionP2AccuCL( loc_, avg_);
        p->set_n( n_)
          .set_n_invalid (n_invalid_)
          .set_check_averaging (check_averaging_)
          .set_level_set_function (ls, lsetbnd, lat);
       return p;
    }
};

} // end of namespace DROPS

#endif
