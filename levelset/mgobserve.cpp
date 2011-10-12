/// \file mgobserve.cpp
/// \brief Observer-base-class for MultiGridCL-changes through AdaptriangCL
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
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#include "levelset/mgobserve.h"
#ifdef _PAR
#  include "parallel/DiST.h"
#endif

namespace DROPS{

/** Merge the known and received DoF values into the new vector \see VecDescCL.

    \todo Do we need the priority PrioHasUnk?
*/
void MGObserverCL::CopyVecElements( VecDescCL& new_vec_desc, const std::vector<Usint>& dims)
{
    const Uint       new_idx = new_vec_desc.RowIdx->GetIdx();   // new index
    const Uint       old_idx = GetIdxDesc()->GetIdx();          // old index
          VectorCL&  new_vec = new_vec_desc.Data;             // new vector
    const VectorCL&  old_vec = *GetVector();                    // old vector
    const std::vector<double>& recv_buf= GetRecvBuffer();       // receive buffer

    // Create a WholeRemoteDataIteratorCL
    DiST::LevelListCL lvls;
    DiST::PrioListCL prios; prios.push_back( PrioMaster); 
    prios.push_back( PrioHasUnk);   // do we need this?
    DiST::Helper::WholeRemoteDataIteratorCL it( dims, lvls, prios, false);

    // iterate over all simplices and copy the DoF values
    for ( ; it!=it.GetEnd(); ++it){
        DiST::TransferableCL& simplex= it->second.GetLocalObject();
        UnknownHandleCL& Unknowns= simplex.Unknowns;
        if ( Unknowns.Exist( new_idx)){
            Assert( Unknowns.Exist( old_idx), 
                DROPSErrCL("LevelsetRepairCL::post_migrate: Unknowns do not exist before migration"),
                DebugParallelNumC);

            // copy data from known DoF values or from the receive buffer
            double const * in= !Unknowns.UnkRecieved( old_idx) 
                                ? Addr(old_vec)  + Unknowns(old_idx)
                                : Addr(recv_buf) + Unknowns(old_idx);
            // ... into the new vector
            double * out = Addr(new_vec)+Unknowns(new_idx);
            std::copy( in, in + new_vec_desc.RowIdx->NumUnknownsSimplex( simplex), out);
        }
    }
}

}   // end of namespace DROPS
