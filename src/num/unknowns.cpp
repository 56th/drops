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
#  include "DiST/DiST.h"
#  include "DiST/mpistream.h"
#  include "parallel/migrateunknowns.h"
#endif
#include <algorithm>


namespace DROPS
{


UnknownIdxCL::UnknownIdxCL( const UnknownIdxCL& orig)
    : _Idx( orig._Idx)
#ifdef _PAR
      , received_(orig.received_)
#endif
{}


UnknownIdxCL& UnknownIdxCL::operator=( const UnknownIdxCL& rhs)
{
    if(&rhs == this) return *this;
    _Idx= rhs._Idx;
#ifdef _PAR
    received_= rhs.received_;
#endif
    return *this;
}

void UnknownHandleCL::DebugInfo( std::ostream& os) const
{
#ifdef _PAR
    os << "in triang ";
    for (size_t i= 0; i<mayHaveUnk_.size(); ++i)
        if (mayHaveUnk_[i]) os << i << " ";
    os<< ", ";
#endif
    if (!_unk)
        os << "no DoFs stored.\n";
    else {
        os << "(sysnum, dof) = ";
        for (Uint i=0; i< _unk->GetNumSystems(); ++i)
            if (Exist(i))
                os << "( " << i << ", " << _unk->GetIdx(i) << ")   ";
            else
                os << "( " << i << ", invalid)   ";
        os << '\n';
    }
}

#ifdef _PAR

/** This function put information---if existent---onto the sendstream. Therefore,
    it used the singleton ObservedVectorsCL. This function assumes the following
    condition:
    \cond it is assumed, that all processes store the "same" observers, i.e.,
          same ordering of the vectors and indices
    \todo
    \param sendstream where to put the data
    \param t          the transferable object which is used to determine the number of unknowns
*/
void UnknownHandleCL::Pack( DiST::MPIostreamCL& sendstream, const DiST::TransferableCL& t) const
{
//const DiST::GeomIdCL gid( 1, MakePoint3D(0.5, 0.5, 0.5), 0);
//bool report= gid == t.GetGID();
//if (report) DebugInfo( cdebug << ">> Pack " << gid);
    const Uint NoIdxFlag= static_cast<Uint>(-1);                            // flag for sending/receiving no data
    if ( Exist()){
        sendstream << true;
        const ObservedMigrateFECL& observers= ObservedMigrateFECL::Instance();
        // write value of dofs onto the stream for all observed vectors
        for ( ObservedMigrateFECL::const_iterator it( observers.begin()); it!=observers.end(); ++it){
            IdxDescCL const *  idx_desc= it->GetIdxDesc();
            VectorCL  const *  vec     = it->GetVector();
            const Uint idx             = idx_desc->GetIdx();
            // if there are unknowns, transfer them, otherwise, send NoDataC
            if ( Exist( idx)){
                const IdxT dof= (*this)(idx);
                sendstream << idx;
                for ( Uint i=0; i<idx_desc->NumUnknownsSimplex( t); ++i){  // edges must have the same number of unknowns than vertices!
                    sendstream << (*vec)[dof+i];
                }
            }
            else{
                sendstream << NoIdxFlag;        // flag, that the proc does not have a dof for that idx
            }
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
void UnknownHandleCL::UnPack( DiST::MPIistreamCL& recvstream, const DiST::TransferableCL& t)
{
//const DiST::GeomIdCL gid( 1, MakePoint3D(0.5, 0.5, 0.5), 0);
//bool report= gid == t.GetGID();
//if (report) DebugInfo( cdebug << "<< Unpack " << gid);
    bool unk_recv= false;
    const Uint NoIdxFlag= static_cast<Uint>(-1);                    // flag for sending/receiving no data
    Uint idx_recv= NoIdxFlag;
    double val_recv;                                                // buffer for receiving a DoF value
    recvstream >> unk_recv;

    if ( unk_recv){
        ObservedMigrateFECL& observers= ObservedMigrateFECL::Instance();
        for ( ObservedMigrateFECL::iterator it( observers.begin()); it!=observers.end(); ++it){
            IdxDescCL const *    idx_desc= it->GetIdxDesc();
            std::vector<double>& recvbuf = it->GetRecvBuffer();
            const Uint idx               = idx_desc->GetIdx();
            // Read the data ...
            recvstream >> idx_recv;
            if ( idx_recv!=NoIdxFlag){    // OK, we are receiving 'good' data
                // check for errors
                Assert( idx==idx_recv,
                    DROPSErrCL("UnknownHandleCL::UnPack: Mismatch in received data!"),
                    DebugParallelNumC);

                if (Exist(idx)) { // DOF already exists, so we just skip the received data
                    for ( int i=0; i<(int)idx_desc->NumUnknownsSimplex( t); ++i)
                        recvstream >> val_recv;
                    continue;
                }
                // Now, prepare for receiving data ...
                Prepare( idx);
                // ... remember that this is a received unknown
                SetUnkReceived( idx);
                // ... remember where the dofs are put
                (*this)(idx)= recvbuf.size();
                // ... and store the received dofs in the buffer
                for ( int i=0; i<(int)idx_desc->NumUnknownsSimplex( t); ++i){
                    recvstream >> val_recv;
                    recvbuf.push_back( val_recv);
                }
            }   // end of receiving the index idx
        }   // end of the observer loop
    }
}
#endif

} // end of namespace DROPS
