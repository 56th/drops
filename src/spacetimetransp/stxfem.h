/// \file stxfem.h
/// \brief space time xfem classes: IdxDescCL
/// \author LNM RWTH Aachen: Christoph Lehrenfeld

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

#ifndef DROPS_STXFEM_H
#define DROPS_STXFEM_H

#include "misc/problem.h"
#include "num/spacetime_geom.h"
#include "num/spacetime_quad.h"


#ifdef _PAR
#  include "parallel/parallel.h"   // for parallel reductions
/* #  include "parallel/interface.h"  // for accumulation of vectors */
#endif

namespace DROPS
{
namespace STXFEM 
{ // namespace for space time xfem helper functions 

void UpdateSTXNumbering( IdxDescCL* Idx, const MultiGridCL& mg, const VecDescCL& lset_old, 
                         const VecDescCL& lset_new, 
                         instat_scalar_fun_ptr lset_fpt, const TimeInterval ti,
                         const BndDataCL<>& lsetbnd, 
                         bool NumberingChanged, double vmax);


enum TIMEDIRECTION
{
    PAST = 0,
    FUTURE =1,
    NOTIME =2
};

template <TIMEDIRECTION timedir>
void GetTrace( const VecDescCL& stsol, VecDescCL& sol, const MultiGridCL& mg);
void GetFutureTrace( const VecDescCL& stsol, VecDescCL& sol, const MultiGridCL& mg);
void GetPastTrace( const VecDescCL& stsol, VecDescCL& sol, const MultiGridCL& mg); 


}//end of namespace DROPS::STXFEM

}//end of namespace DROPS

#endif
