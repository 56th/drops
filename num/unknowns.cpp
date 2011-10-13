/// \file unknowns.cpp
/// \brief Implementation of the mapping from simplices to indices to
///    linear-algebra data-structures.)
/// \author LNM RWTH Aachen: Sven Gross, Joerg Peters, Volker Reichelt; SC RWTH Aachen:

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

#include "num/unknowns.h"
#ifdef _PAR
#  include "parallel/DiST.h"
#  include "parallel/mpistream.h"
#  include "levelset/mgobserve.h"
#endif
#include <algorithm>


namespace DROPS
{


UnknownIdxCL::UnknownIdxCL( const UnknownIdxCL& orig)
    : _Idx( orig._Idx), received_(orig.received_) {}


UnknownIdxCL& UnknownIdxCL::operator=( const UnknownIdxCL& rhs)
{
    if(&rhs == this) return *this;
    _Idx= rhs._Idx;
    return *this;
}

#ifdef _PAR

/** This function put information---if existent---onto the sendstream. Therefore,
    it used the singleton ObservedVectorsCL. This function assumes the following
    two conditions:
    \cond it is assumed, that all processes store the "same" observers, i.e.,
          same ordering of the vectors and indices
    \param sendstream where to put the data
    \param t          the transferable object which is used to determine the number of unknowns
*/
void UnknownHandleCL::Pack( DiST::Helper::MPIostreamCL& sendstream, const DiST::TransferableCL& t) const
{
    const double NoDataC= std::numeric_limits<double>::quiet_NaN();
    if ( Exist()){
        sendstream << true;
        const ObservedVectorsCL& observers= ObservedVectorsCL::Instance();
        // write value of dofs onto the stream for all observed vectors
        for ( ObservedVectorsCL::const_iterator it( observers.begin()); it!=observers.end(); ++it){
            IdxDescCL const *  idx_desc= (*it)->GetIdxDesc();
            VectorCL  const *  vec     = (*it)->GetVector();
            // The vector and the index description is available (assume that this is true
            // on receiving process, too).
            if ( vec && idx_desc){
                const Uint idx= idx_desc->GetIdx();
                // check if there are unknowns
                if ( Exist( idx)){
                    const IdxT dof= (*this)(idx);
                    for ( Uint i=0; i<idx_desc->NumUnknownsSimplex( t); ++i){  // edges must have the same number of unknowns than vertices!
                        sendstream << (*vec)[dof+i];
                    }
                }
                else{
                    sendstream << NoDataC;    // flag, that the proc does not have a dof for that idx
                }
            }   // end of (vec && idx_desc)
        }   // end of observer loop
    }
    else{           // There are no DOF values available
        sendstream << false;
    }
}

/** If any DOF information are received, handle these DOFs. That is, put a copy of the
    DoF value into the corresponding receive buffer in the \see MGObserverCL. Furthermore,
    mark that DoFs have been received.

    \cond no DoFs of respectivly each type are known before the migration takes place.

    \param recvstream the stream where all information can be found
    \param t          the transferable object which is used to determine the number of unknowns
*/
void UnknownHandleCL::UnPack( DiST::Helper::MPIistreamCL& recvstream, const DiST::TransferableCL& t)
{
    bool unk_recv= false;
    const double NoDataC= std::numeric_limits<double>::quiet_NaN(); // flag for sending/receiving no data
    double val_recv;                                                // buffer for receiving a DoF value
    recvstream >> unk_recv;

    if ( unk_recv){
        ObservedVectorsCL& observers= ObservedVectorsCL::Instance();
        for ( ObservedVectorsCL::const_iterator it( observers.begin()); it!=observers.end(); ++it){
            IdxDescCL const *    idx_desc= (*it)->GetIdxDesc();
            std::vector<double>& recvbuf = (*it)->GetRecvBuffer();
            // The vector and the index description is available
            if ( (*it)->GetVector() && idx_desc){
                const Uint idx= idx_desc->GetIdx();
                recvstream >> val_recv;
                if ( val_recv!=NoDataC){    // OK, we are receiving 'good' data
                    // check for errors
                    Assert( !Exist(idx),
                        DROPSErrCL("UnknownHandleCL::UnPack: Merging of received dof is not possible"),
                        DebugParallelNumC);

                    // Now, prepare for receiving data ...
                    Prepare( idx);
                    // ... remember that this is a received unknown
                    SetUnkRecieved( idx);
                    // ... remember where the dofs are put
                    (*this)(idx)= recvbuf.size();
                    // ... and the store the received dofs in the buffer
                    recvbuf.push_back( val_recv);   // first value has already been received
                    for ( int i=0; i<(int)idx_desc->NumUnknownsSimplex( t)-1; ++i){
                        recvstream >> val_recv;
                        recvbuf.push_back( val_recv);
                    }
                }   // end of receiving the index idx
            }   // end if (vec && idx_desc)
        }   // end of the observer loop
    }
}
#endif

} // end of namespace DROPS
