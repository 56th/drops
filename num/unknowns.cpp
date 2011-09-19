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
#endif
#include <algorithm>


namespace DROPS
{


UnknownIdxCL::UnknownIdxCL( const UnknownIdxCL& orig)
    : _Idx( orig._Idx) {}


UnknownIdxCL& UnknownIdxCL::operator=( const UnknownIdxCL& rhs)
{
    if(&rhs == this) return *this;
    _Idx= rhs._Idx;
    return *this;
}

#ifdef _PAR

/** For all system numbers available in _unk, put the system number and the 
    value of the DOF onto the stream sendstream. If there are no DOFs of a 
    system number, then std::numeric_limits<Uint>::max() is written 
    onto the stream.
    \param sendstream where to put all information
*/
void UnknownIdxCL::Pack( DiST::Helper::SendStreamCL& sendstream) const
{
    for ( Uint sysnum=0; sysnum<GetNumSystems(); ++sysnum){
        if ( _Idx[sysnum]!=NoIdx)
            sendstream << sysnum; /// \TODO: Write DOF values on the stream
        else
            sendstream << std::numeric_limits<Uint>::max();
    }
}

/** Get all DOF from the recvstream and store them in the receive buffer.
    Furthermore, all DOF values are marked as received.
    \param recvstream the stream where all information can be found
*/
void UnknownIdxCL::UnPack( DiST::Helper::RecvStreamCL& recvstream)
{
    Uint recv_sysnum;       // received system number
    double recv_DOFvalue;   // received DOF value
    // Allocate memory for receive flags
    received_.resize( GetNumSystems(), false);

    // Get all information from the stream
    for ( Uint sysnum=0; sysnum<GetNumSystems(); ++sysnum){
        recvstream >> recv_sysnum;
        // there are some information about DOF, so unpack them
        if ( recv_sysnum!=std::numeric_limits<Uint>::max()){
            Assert( recv_sysnum==sysnum /* && sysnum= XXX.GetIdx() */, 
                DROPSErrCL("UnknownHandleCL::UnPack: Missing information about a system number"), 
                DebugUnknownsC | DebugParallelC);
            /// \TODO: Receive DOF values and store them!
            recvstream >> recv_DOFvalue;
            SetUnkRecv( sysnum);
        }
    }
}

/// \brief Write UnknownIdxCL onto a send stream
inline DiST::Helper::SendStreamCL& operator<< ( DiST::Helper::SendStreamCL& sendstream, const UnknownIdxCL& unk)
{
    unk.Pack( sendstream);
    return sendstream;
}

/// \brief Read UnknownIdxCL from a receive stream
inline DiST::Helper::RecvStreamCL& operator>> ( DiST::Helper::RecvStreamCL& recvstream, UnknownIdxCL& unk)
{
    unk.UnPack( recvstream);
    return recvstream;
}

/** If there are any DOF information available, then put this onto the
    the send stream.
    \param sendstream where to put the data
*/
void UnknownHandleCL::Pack( DiST::Helper::SendStreamCL& sendstream) const
{
    if ( Exist()){  // There are DOF values available
        sendstream << _unk->GetNumSystems() << *_unk;
    }
    else{           // There are no DOF values available
        sendstream << Uint(0);
    }
}

/** If any DOF information are received, handle these DOFs.
    \param recvstream the stream where all information can be found    
*/
void UnknownHandleCL::UnPack( DiST::Helper::RecvStreamCL& recvstream)
{
    // get the number of system numbers written onto stream by Pack
    Uint num_sysnums;
    recvstream >> num_sysnums;

    // if there are any DOF information, receive them.
    if ( num_sysnums){
        Assert( _unk==0, DROPSErrCL("UnknownHandleCL::UnPack: UnknownIdxCL already constructed"), DebugUnknownsC | DebugParallelC);
        // Allocate memory for storing all 
        _unk= new UnknownIdxCL( num_sysnums);
        recvstream >> *_unk;
    }
}
#endif

} // end of namespace DROPS
