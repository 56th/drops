/// \file DiST.tpp
/// \brief The DiST (short for distributed simplex type) module is responsible for the distributed geometric data structure on a parallel machine.
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Yuanjun Zhang; SC RWTH Aachen: Oliver Fortmeier, Daniel Medina Cardona.

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

namespace DROPS {

namespace DiST{


// I N T E R F A C E  C L
// ----------------------

inline MPIostreamCL& operator<< (MPIostreamCL& os,
    const InterfaceCL::MessagesCL& m)
{
    os << m.numData;
    os.write( Addr( m.messages), m.messages.size());
#   if DROPSDebugC & DebugDiSTC
        // append delimiting char to find inconsistent gather/scatter routines
        os << '|';
#   endif
    return os;
}

template <typename HandlerT>
void InterfaceCL::GatherData( HandlerT& handler, const iterator& begin,
    const iterator& end, CommPhase phase)
/** This function also allocates the memory for sending data to the owner processes. */
{
    for ( iterator it( begin); it != end; ++it) {
        const int owner= it->second.GetOwnerProc();
        // Generates the sendbuf_ for owner, if not already there.
        // Needed for correctness, because the owner expects a (n empty) message.
        SendStreamCL& buf= sendbuf_[owner];
        if (phase==fromowner && !it->second.AmIOwner()) // Check, if this process has to gather data.
            continue;
        // Check if gather wants to put something on the send stream and write
        // the GeomIdCL, and the associated sub-sendstream on the stream.
        SendStreamCL tmp_buf( binary_);
        if (handler.Gather( it->second.GetLocalObject(), tmp_buf))
            buf << it->first << tmp_buf;
    }
}

template <typename HandlerT>
bool InterfaceCL::ScatterData (HandlerT& handler)
/// \return local accumulated AND of all scatter calls
{
    bool result= true;
    for (RecvListT::iterator it( recvbuf_.begin()); it != recvbuf_.end(); ++it)
        result= result && ScatterData( handler, it->second);
    return result;
}

template <typename HandlerT, typename IStreamT>
bool InterfaceCL::ScatterData( HandlerT& handler, IStreamT& recv)
{
    bool result= true;
    GeomIdCL gid;
    size_t numData;
    while ((recv >> gid).good()) {
        recv >> numData;
#       if DROPSDebugC & DebugDiSTC
        if (!recv) {
            cdebug << "error while reading object " << gid << std::endl;
            throw DROPSErrCL("InterfaceCL::ScatterData: Receive stream is broken!");
        }
#       endif
        RemoteDataCL& rd= InfoCL::Instance().GetRemoteData( gid);
        const bool scatter_result= handler.Scatter( rd.GetLocalObject(), numData, recv);
        result= result && scatter_result;
#       if DROPSDebugC & DebugDiSTC
            // check for delimiter
            char delim= '|';
            recv >> delim;
            Assert( delim=='|', ErrorCL("InterfaceCL::ScatterData: "
                "incomplete receive while reading object ", gid), DebugDiSTC);
#       endif
     }
     return result;
}


template <typename HandlerT>
bool InterfaceCL::Perform( HandlerT& handler, CommPhase phase)
{
    SetupCommunicationStructure(); ///\todo Call of SetupCommunicationStructure in constructor?

    // Gather the data
    GatherData( handler, begin_from_, end_, phase);

    // communicate data
    if (phase == toowner) {
        loc_recv_toowner_.setBinary( binary_);
        loc_send_toowner_.setBinary( binary_);
    }
    ExchangeData( phase);

    // scatter the data
    bool result= ScatterData( handler);
    if (phase == toowner) {
        result= result && ScatterData( handler, loc_recv_toowner_);
        loc_recv_toowner_.clear();
        loc_recv_toowner_.setbuf( 0, 0);
        loc_send_toowner_.clear();
        loc_send_toowner_.clearbuffer();
    }
    // clear receive buffer
    recvbuf_.clear();

    return result;
}


template <typename HandlerT>
bool InterfaceCL::Communicate( HandlerT& handler)
/** Basically, call the Gather function of the handler for all entities covered by the interface
    and having the \a from priority. Then transfer the data to the respective owners which bundle
    the data and send the data bundles to the copies having the \a to priority. After receiving
    the data, call the Scatter function on the copies.
    \return Local accumulated AND of all Scatter calls (without reduction).
*/
{
    return Perform( handler, bothPhases);
}

template <typename HandlerT>
bool InterfaceCL::InformOwners( HandlerT& handler)
{
    return Perform( handler, toowner);
}

template <typename HandlerT>
bool InterfaceCL::InformCopies( HandlerT& handler)
{
    return Perform( handler, fromowner);
}


template <typename HandlerT>
bool InterfaceCL::ExecuteLocal( HandlerT& handler)
{
    return ExecuteLocal( handler, begin_from_, end_);
}

template <typename ExecuteHandlerT, typename IteratorT>
bool InterfaceCL::ExecuteLocal( ExecuteHandlerT& handler, const IteratorT& begin, const IteratorT& end)
{
    bool result=true;
    for (IteratorT it(begin); it!=end; ++it) {
        const bool execute_result= handler( it->second.GetLocalObject());
        result= result && execute_result;
    }
    return result;
}


// T R A N S F E R  C L A S S
//---------------------------

// template specialization for tetras
template <>
void DiST::TransferCL::ReceiveSimplices<TetraCL>( DiST::RecvStreamCL& recvstream, size_t num);

/** Receive \a num simplices of type SimplexT from the \a recvstream and update the
    corresponding remote data list
*/
template <typename SimplexT>
void DiST::TransferCL::ReceiveSimplices( DiST::RecvStreamCL& recvstream, size_t num)
{
    for ( size_t i=0; i<num; ++i) {
        SimplexT stmp;                                  // temporary to receive a single simplex
        DiST::RemoteDataCL::ProcListT procList; // temporary to receive process list

        // receive simplex
        stmp.UnPack( recvstream);
        Assert( stmp.GetDim()==GetDim<SimplexT>(), DROPSErrCL("Mismatch in dimension of a received simplex!"), DebugDiSTC);
        // receive proc/prio list
        recvstream >> procList;
        CreateSimplex<SimplexT>( stmp, procList); // creates simplex and remote data list entry
    }
}

// I N F O  C L A S S
// ------------------

const RemoteDataCL& InfoCL::GetRemoteData( const GeomIdCL& h) const
{
    RemoteDataListCL::const_iterator it= remoteData_[h.dim].find( h);
    Assert( it!=remoteData_[h.dim].end(), ErrorCL( "InfoCL::GetRemoteData (const): Simplex not registered: ", h, NoGID), DebugDiSTC);
    return it->second;
}

RemoteDataCL& InfoCL::GetRemoteData( const GeomIdCL& h)
{
    RemoteDataListCL::iterator it= remoteData_[h.dim].find( h);
    Assert( it!=remoteData_[h.dim].end(), ErrorCL( "InfoCL::GetRemoteData: Simplex not registered: ", h, NoGID), DebugDiSTC);
    return it->second;
}

InfoCL& InfoCL::Instance( MultiGridCL* mg)
{
    Assert( instance_!=0 || mg!=0, DROPSErrCL("InfoCL::GetInstance: Cannot initialize InfoCL without a pointer to a multigrid"), DebugDiSTC);
    return instance_==0 ? *( instance_= new InfoCL(mg)) : *instance_;
}

InfoCL* InfoCL::InstancePtr()
{
    Assert( instance_!=0, DROPSErrCL("InfoCL::GetInstancePtr: InfoCL is not initialized"), DebugDiSTC);
    return instance_;
}


}   // end of namespace DiST
}   // end of namespace DROPS

