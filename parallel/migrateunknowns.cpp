/// \file migrateunknowns.cpp
/// \brief Observer-base-class migrating unknowns
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross; SC RWTH Aachen: Oliver Fortmeier

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

#include "parallel/migrateunknowns.h"
#include "levelset/mgobserve.h"

namespace DROPS{

/** Get the index description provided by the observer. 
    \todo If this function is called very often, e.g. by the UnknownHandlerCL::Pack,
          it would be nice to have it inlined. (Store a local copy of the index description?)
*/
const IdxDescCL* MigrateFECL::GetIdxDesc() const
{
    return obs_->GetIdxDesc();
}

/** Get the vector of the old data by the corresponding observer.
    \todo If this function is called very often, e.g. by the UnknownHandlerCL::Pack,
          it would be nice to have it inlined. (Store a local copy of the vector?)
*/
const VectorCL* MigrateFECL::GetVector() const
{
    return obs_->GetVector();
}

/** Merge the known and received DoF values into the new vector new_vec_desc.
    \todo Do we need the priority PrioHasUnk?
*/
void MigrateFECL::CopyVecElements( VecDescCL& new_vec_desc, const std::vector<Usint>& dims)
{
    const Uint       new_idx = new_vec_desc.RowIdx->GetIdx();   // new index
    const Uint       old_idx = GetIdxDesc()->GetIdx();          // old index
          VectorCL&  new_vec = new_vec_desc.Data;               // new vector
    const VectorCL&  old_vec = *GetVector();                    // old vector

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
                                : Addr(recvBuf_) + Unknowns(old_idx);
            // ... into the new vector
            double * out = Addr(new_vec)+Unknowns(new_idx);
            std::copy( in, in + new_vec_desc.RowIdx->NumUnknownsSimplex( simplex), out);
            // .. this unknown has been handled, so forget about the received-status
            Unknowns.ResetUnkRecieved( old_idx);
        }
    }
}

/** Clear all entries which are given in the receive buffer and reserve 
    some memory.
*/
void MigrateFECL::pre_migrate()
{
    Assert( obs_->GetVector(), 
        DROPSErrCL("MigrateFECL::pre_migrate: Vector data not set"), 
        DebugParallelNumC);
    recvBuf_.clear();
    recvBuf_.reserve( GetIdxDesc()->NumUnknowns());
}

/** After the migration has take place, build the local vector of 
    accumulated vector indices. Afterwards, fill this vector by all the
    DoF which are either already known or just received.
    \post a new index number is created for the FE function.
*/
void MigrateFECL::post_migrate()
{
    // Abbrevations
    const IdxDescCL& old_idx= *GetIdxDesc();

    // Create a new numbering
    VecDescCL loc_vec;
    IdxDescCL loc_idx( old_idx.GetFE());
    const Uint lvl = old_idx.TriangLevel();

    loc_idx.CreateNumbering( lvl, *mg_, old_idx);

    // Assign the new index (and allocate memory for the new vector)
    loc_vec.SetIdx( &loc_idx);

    // Copy the DoF values into the new vector
    DiST::Helper::WholeRemoteDataIteratorCL::DimListT dims; 
    if ( old_idx.NumUnknownsVertex())  dims.push_back( 0);
    if ( old_idx.NumUnknownsEdge())    dims.push_back( 1);
    if ( old_idx.NumUnknownsFace())    dims.push_back( 2);
    if ( old_idx.NumUnknownsTetra())   dims.push_back( 3);
    
    CopyVecElements( loc_vec, dims);

    obs_->swap( loc_idx, loc_vec.Data);
    recvBuf_.clear();
}


/** Handle all FE types where the Vector is provided by the ObservedVectorsCL.
    Furthermore, store a reference to the \see MultiGridCL to create numberings.
*/
void ObservedMigrateFECL::Init( ObservedVectorsCL& obsvec, MultiGridCL& mg)
{
    for ( ObservedVectorsCL::iterator obs= obsvec.begin(); obs!=obsvec.end(); ++obs){
        (*obs)->pre_refine_sequence();
        if ( (*obs)->GetVector()){
            push_back( MigrateFECL( *obs, mg));
        }
    }
}

}   // end of namespace DROPS
