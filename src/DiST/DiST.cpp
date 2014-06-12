/// \file DiST.cpp
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

#include "DiST/DiST.h"
#include "geom/simplex.h"
#include "geom/multigrid.h"

#include <set>

namespace DROPS{
namespace DiST{



std::ostream& operator<< ( std::ostream& os, const SimplexTransferInfoCL::ProcSetT& pl)
// for debugging
{
    SimplexTransferInfoCL::ProcSetT::const_iterator it( pl.begin());
    if (it!=pl.end()) {
        os << '(' << it->first << ',' << PriorityToString( it->second) << ')';
        ++it;
    }
    for (; it!=pl.end(); ++it) {
        os << ", (" << it->first << ',' << PriorityToString( it->second) << ')';
    }
    return os;
}

// S I M P L E X  T R A N S F E R  I N F O  C L
//---------------------------------------------

void SimplexTransferInfoCL::AddProc( int p, Priority prio, UpdatePolicyE changeLocalPrio)
{
    ProcSetT::iterator it= postProcs_.find(p);
    if (it == postProcs_.end()) // not in proc list, yet
        postProcs_[p]= prio;
    else if (p != ProcCL::MyRank()) // merge prios
        it->second= merge_prio( it->second, prio);
    else { // p == me
        switch (changeLocalPrio) {
            case overwrite: // set prio
                it->second= prio; break;
            case merge: // merge with present prio
                it->second= merge_prio( it->second, prio); break;
            default:; // do nothing
        }
    }
}

void SimplexTransferInfoCL::AddProcSet( const ProcSetT& procs)
{
    for (ProcSetT::const_iterator it= procs.begin(), end= procs.end(); it!=end; ++it) {
        const int proc= it->first;
        if (postProcs_.find(proc) == postProcs_.end()) { // not in proc list, yet
            Priority prio= rd_.GetPrio(proc);
            if (prio==NoPrio) // not found in remote data
                prio= PrioMaster;
            postProcs_[proc]= prio;
        }
    }
}

int SimplexTransferInfoCL::GetPostProc( Priority prio) const
{
    for (ProcSetT::const_iterator it= postProcs_.begin(), end= postProcs_.end(); it!=end; ++it)
        if (it->second==prio)
          return it->first;
    return -1;
}

void SimplexTransferInfoCL::ComputeSendToProcs( bool tetra)
/// For tetrahedra (\a tetra == true) we have to treat the special case that the object should be sent also to processes which already have a remote copy.
{
    const int me= ProcCL::MyRank();
    if (tetra) {
        if (postProcs_.size()==2 && rd_.GetNumProcs()==2) { // Ma/Gh before and after transfer: only send to remote post procs with same prio
            procsToSend_.clear();
            Priority myprio= rd_.GetLocalPrio();
            const int p= GetPostProc(myprio);
            if (p!=-1 && p!=me)
                procsToSend_[p]= myprio;
        } else { // send to all remote post procs
            procsToSend_= postProcs_;
            procsToSend_.erase( me);
        }
    } else { // for all sub-simplices (non-tetra): procsToSend_ = postProcs_ - remote data proc list
        procsToSend_= postProcs_;
        for (RemoteDataCL::ProcListT::const_iterator it= rd_.GetProcListBegin(), end= rd_.GetProcListEnd(); it!=end; ++it)
            procsToSend_.erase( it->proc);
    }
}


// I N T E R F A C E  C L
// ----------------------

void Print( const std::set<int>& set, std::ostream& os)
{
    for( std::set<int>::const_iterator it=set.begin(), end=set.end(); it!=end; ++it)
        os << *it << " ";
    os << std::endl;
}

void InterfaceCL::SetupCommunicationStructure()
/** For a more detailed description, see InterfaceCL::Communicate, where the use of
    ownerRecvFrom_, ownerSendTo_ and IRecvFromOwners_ is given. */
{
    ownerRecvFrom_.clear();
    ownerSendTo_.clear();
    IRecvFromOwners_.clear();
    // note: ISendToOwner is not needed, created on the fly in GatherData(...)

    typedef RemoteDataCL::ProcList_const_iterator ProcList_const_iterator;

    // fill ownerRecvFrom_ and ownerSendTo_:
    // each owner has to determine the process ranks to receive data from and send data to
    for ( iterator it = begin_; it != end_; ++it) {
        const int owner= it->second.GetOwnerProc();
        if ( owner==ProcCL::MyRank()){
            ProcList_const_iterator pit=it->second.GetProcListBegin();
            for ( ; pit!=it->second.GetProcListEnd(); ++pit){
                if (from_.contains( pit->prio))
                    ownerRecvFrom_.insert( pit->proc);
                if (to_.contains( pit->prio))
                    ownerSendTo_.insert( pit->proc);
            }
        }
    }
    // fill IRecvFromOwners_:
    // during owner->to communication, each process p has to determine the owner processes to receive data from
    for ( iterator it = begin_to_; it != end_; ++it) {
        if (to_.contains( it->second.GetLocalPrio()))
            IRecvFromOwners_.insert( it->second.GetOwnerProc());
    }
//    std::cout << "[" << ProcCL::MyRank() << "] ownerRecvFrom  = "; Print(ownerRecvFrom_, std::cout);
//    std::cout << "[" << ProcCL::MyRank() << "] ownerSendTo    = "; Print(ownerSendTo_, std::cout);
//    std::cout << "[" << ProcCL::MyRank() << "] IRecvFromOwner = "; Print(IRecvFromOwners_, std::cout);
}

void InterfaceCL::SendData (SendListT& sendbuf, std::vector<ProcCL::RequestT>& req, int tag)
{
    req.reserve( sendbuf.size());
    for (SendListT::iterator it= sendbuf.begin(); it != sendbuf.end(); ++it) {
        if ( it->first != ProcCL::MyRank()) {
            // initiate the MPI call for sending the data (in all other cases, the data are just copied)
            req.push_back( it->second.Isend( it->first, tag));
        }
    }
}

/// \brief Reads a stream (gid1, numdata1, data1, gid2, numdata2, data2, ...) until the stream fail()s.
/// The data is put in collect[gid].
template <class IStreamT>
void collect_streams (IStreamT& recv, InterfaceCL::CollectDataT& collect)
{
    GeomIdCL gid;

    while ((recv >> gid).good()) {
        RefMPIistreamCL gid_data( 0, 0, recv.isBinary());
        recv >> gid_data;
        collect[gid].append( gid_data.begin(), gid_data.end());
    }
}

void InterfaceCL::to_owner (std::vector<ProcCL::RequestT>& reqFirstSend, CollectDataT& collect)
{
    const int myrank= ProcCL::MyRank();
    const int firstSendTag= 5; // tag for sending in phase (1)

    // Phase (1): Send data to owning processes
    //-----------------------------------------
    SendData( sendbuf_, reqFirstSend, firstSendTag);

    // Phase (2): Owning process receives data
    //----------------------------------------
    for (InterfaceCL::ProcSetT::const_iterator sender= ownerRecvFrom_.begin(); sender != ownerRecvFrom_.end(); ++sender) {
        if (*sender != myrank) {
            RecvStreamCL locrecvbuf( binary_);
            collect_streams( locrecvbuf.Recv( *sender, firstSendTag), collect);
        }
        else {
            SendStreamCL& locsendbuf= sendbuf_[myrank];
            RefMPIistreamCL locrecvbuf( locsendbuf.begin(),
                locsendbuf.cur() - locsendbuf.begin(), binary_);
            collect_streams( locrecvbuf, collect);
        }
    }
}

void InterfaceCL::from_owner (std::vector<ProcCL::RequestT>& reqSecondSend, SendListT& sendstreams)
{
    const int myrank= ProcCL::MyRank();
    const int secondSendTag= 6; // tag for sending in phase (4)

    // Phase (4): send buffers to non-owner copies
    // --------------------------------------------
    SendData( sendstreams, reqSecondSend, secondSendTag);

    // Phase (5): Receive data from owning processes
    // ----------------------------------------------
    for (InterfaceCL::ProcSetT::const_iterator sender= IRecvFromOwners_.begin();
        sender != IRecvFromOwners_.end(); ++sender)
        if (*sender != myrank)
            recvbuf_[*sender].Recv( *sender, secondSendTag);
        else
            recvbuf_[myrank].clearbuffer( sendstreams[myrank].str());
}

void InterfaceCL::ExchangeData (CommPhase phase)
/** The interface communication has four phases.
    (1) The data collected by the gather routine is sent to processes which are owner of
    at least one entity. Let p denote a process who is owner of an entity.
    (2) Then, p receives data from the processes given in the member variable
    ownerRecvFrom_ and orders the data according the GID of the corresponding entity.
    (3) After receiving the data from the processes owning a non-owner copy, p generates
    a send stream containing the following data:
    GID_1 numData_1 data_1, GID_2 numData_2 data_2, ..., GID_n numData_n data_n NoGID.
    Here, numData_i describes the number of copies of the entity with GID_i which have
    sent some data to the owner.
    (4) The data which are collected in phase (3) is sent to the processes who own
    at least one entity covered by the "to"-interface.
    (5) For receiving the data, each process receives data from the previously computed
    list IRecvFromOwners_.
    \param phase How to communicate: bothPhases: copies -> owner -> copies; toowner: copies -> owner; fromowner: owner -> copies

    \todo Split up the function to control the communication better and
          to make it easier to read
    \todo Test for one-way communication is missing since no other priorities as master
          are given so far.
    \todo Instead of numData_i data_i, one should simply write data_i as substream with: sendbuf << data_i.
*/
{
    const int myrank= ProcCL::MyRank();

    // Phase (1): Send data to owning processes
    // Phase (2): Owning process receives data
    //----------------------------------------
    std::vector<ProcCL::RequestT> reqFirstSend;
    CollectDataT collect; // Collect the data for one GID from all senders
    if (phase==bothPhases || phase==toowner)
        to_owner( reqFirstSend, collect);
    else { // Local operation
        if (ownerRecvFrom_.count( myrank) > 0) {
            SendStreamCL& locsendbuf= sendbuf_[myrank];
            RefMPIistreamCL locrecvbuf( locsendbuf.begin(), locsendbuf.cur() - locsendbuf.begin(), binary_) ;
            collect_streams( locrecvbuf, collect);
        }
    }

    // Phase (3): Generate and fill the send buffers
    //-------------------------------------
    // Generate a SendStream for each receiver. Some receivers,
    // which are waiting for a message, might receive an empty message.
    // MPI permits this. (receivers == ownerSendTo_)
    // Note, that the message must be sent nonetheless, as the receiver posts a Recv.
    SendListT sendstreams;
    if (phase == fromowner || phase == bothPhases)
        for (ProcSetT::const_iterator receiver= ownerSendTo_.begin(); receiver != ownerSendTo_.end(); ++receiver)
            sendstreams[*receiver].setBinary( binary_);

    // For each collected GID, put the GID, number of copies, where gather was
    // called, and the gathered data into a stream buffer.
    typedef RemoteDataCL::ProcList_const_iterator PL_IterT;
    for (CollectDataT::iterator it= collect.begin(); it != collect.end(); ++it) {
        RemoteDataCL& rd= InfoCL::Instance().GetRemoteData( it->first);
        for (PL_IterT pit= rd.GetProcListBegin(); pit != rd.GetProcListEnd(); ++pit) {
            if (to_.contains( pit->prio)) {
                const int receiver= pit->proc;
                if (phase==toowner) { // Data will be referenced by loc_recv_toowner_ in Phase (4) and (5) and ScatterData.
                    if (receiver == myrank)
                        loc_send_toowner_ << it->first << it->second;
                }
                else {
                    Assert( sendstreams.count( receiver) > 0,
                        DROPSErrCL("InterfaceCL::Communicate: Missing sendbuffer"),
                        DebugDiSTC);
                    sendstreams[receiver] << it->first << it->second;
                }
            }
        }
    }
    collect.clear();

    // Phase (4): Send buffers to non-owner copies
    // Phase (5): Receive data from owning processes
    // ----------------------------------------------
    std::vector<ProcCL::RequestT> reqSecondSend;
    if (phase == fromowner || phase == bothPhases)
        from_owner( reqSecondSend, sendstreams);
    else
        if (IRecvFromOwners_.count( ProcCL::MyRank()) > 0)
            loc_recv_toowner_.setbuf( loc_send_toowner_.begin(), loc_send_toowner_.cur() - loc_send_toowner_.begin());

    // Wait until all messages have left me before deleting the buffers
    if (!reqFirstSend.empty())
        ProcCL::WaitAll( reqFirstSend);
    if (!reqSecondSend.empty())
        ProcCL::WaitAll( reqSecondSend);
    // clear all send streams
    sendbuf_.clear();
    sendstreams.clear();
}

// M O D I F Y  C L A S S
//-----------------------

/// handler for the DiST interface to merge proc lists
class ModifyCL::MergeProcListHandlerCL
{
  private:
    ModifyCL& mod_;

  public:
    MergeProcListHandlerCL( ModifyCL& mod) : mod_(mod) {}

    bool Gather( const TransferableCL& t, SendStreamCL& s)
    {
        ModifyCL::UpdateListT& ul= mod_.entsToUpdt_[t.GetDim()];
        ModifyCL::UpdateIterator it= ul.find( &t);
        if (it==ul.end())
            return false;
        const RemoteDataCL::ProcListT proclist( it->second.GetPostProcs().begin(), it->second.GetPostProcs().end());
        s << proclist;
        return true;
    }

    bool Scatter( TransferableCL& t, const size_t numData, MPIistreamCL& r){
        ModifyCL::UpdateListT& ul= mod_.entsToUpdt_[t.GetDim()];
        ModifyCL::UpdateIterator it= ul.find( &t);
        if (it==ul.end())
            return false;
        for (size_t i=0; i<numData; ++i) {
            RemoteDataCL::ProcListT proclist;
            r >> proclist;
            for (RemoteDataCL::ProcListT::const_iterator pit= proclist.begin(), pend= proclist.end(); pit!=pend; ++pit)
                it->second.AddProc( pit->proc, pit->prio, SimplexTransferInfoCL::merge);
        }
        return true;
    }

    void Call()
    // do interface comm
    {
        std::vector<GeomIdCL> updateVec;

        // collect all objects to communicate
        for (int dim=0; dim<4; ++dim)
            for (ModifyCL::UpdateIterator it= mod_.UpdateBegin(dim), end= mod_.UpdateEnd(dim); it!=end; ++it)
                updateVec.push_back( it->first->GetGID());

        InterfaceCL comm( updateVec.begin(), updateVec.end(), mod_.binary_);
        comm.Communicate( *this);
    }
};

/// handler for the DiST interface to communicate update lists
class ModifyCL::CommToUpdateHandlerCL
{
  private:
    ModifyCL& mod_;

  public:
    CommToUpdateHandlerCL( ModifyCL& mod) : mod_(mod) {}

    bool Gather( const TransferableCL& t, SendStreamCL& )
    {
        const Usint dim= t.GetDim();
        ModifyCL::UpdateListT&   ul = mod_.entsToUpdt_[dim];
        ModifyCL::UpdateIterator it = ul.find( &t);

        return it != ul.end();
    }

    bool Scatter( TransferableCL& t, const size_t, MPIistreamCL& )
    {
        const Usint dim = t.GetDim();
        if (dim==GetDim<TetraCL>()) { // tetra
            ModifyCL::UpdateListT& ul= mod_.entsToUpdt_[dim];
            ModifyCL::UpdateIterator it= ul.find( &t);
            if (it==ul.end()) { // not already in update list
                it= mod_.AddSimplexToUpdate( dim, &t, false );
                // add (local proc,local prio) to update list, otherwise local Ma/Gh copy will be lost
                it->second.AddProc( ProcCL::MyRank(), t.GetPrio());
            }
        } else // non-tetra
            mod_.AddSimplexToUpdate( dim, &t, false);
        return true;
    }

    void Call()
    // do interface comm
    {
        InterfaceCL::DimListT dimlist;
        for (int dim=0; dim<4; ++dim)
            dimlist.push_back( dim);

        const PrioListT   allPrios;
        const LevelListCL allLvls;
        // communicate over all objects
        InterfaceCL comm( allLvls, allPrios, allPrios, dimlist, /*dist*/ true, mod_.binary_);
        comm.Communicate( *this);
    }
};

void ModifyCL::Init()
{
    Assert( !modifiable_, DROPSErrCL("ModifyCL::Init: Class is already in the modifiable mode"), DebugDiSTC);
    modifiable_= true;
    entsToUpdt_= new UpdateListT[4];
}

void ModifyCL::Finalize()
{
    Assert( modifiable_, DROPSErrCL("TransferCL::Finalize: Class is not in the modifiable mode, call Init() first!"), DebugDiSTC);
    // algorithm based on migration algorithm of FMDB
    CreateUpdateList();
    const bool foundDel= AssignPostProcs();
    UpdateRemoteData();

    if (foundDel) { // unregister/remove unused simplices (depending on del_)
        if (del_) InfoCL::Instance().PreMultiGridMod();
        DeleteUnusedSimplices( del_);
        if (del_) InfoCL::Instance().PostMultiGridMod();
    }
    // Now we are ready to remove the update lists
    delete[] entsToUpdt_;
    entsToUpdt_= 0;
    modifiable_= false;
}

void ModifyCL::ChangePrio( const TransferableCL& t, Priority prio)
{
    Assert( modifiable_, DROPSErrCL("ModifyCL::ChangePrio: Class is not in the modifiable mode, call Init() first!"), DebugDiSTC);

    UpdateIterator it= AddSimplexToUpdate( t.GetDim(), &t, false);
    it->second.AddProc( ProcCL::MyRank(), prio, /*changeLocalPrio*/SimplexTransferInfoCL::overwrite);
}

void ModifyCL::Delete( const TransferableCL& t)
{
    Assert( modifiable_, DROPSErrCL("ModifyCL::Delete: Class is not in the modifiable mode, call Init() first!"), DebugDiSTC);

    AddSimplexToUpdate( t.GetDim(), &t, false);
}

void ModifyCL::Keep( const TransferableCL& t)
{
    Assert( modifiable_, DROPSErrCL("ModifyCL::Keep: Class is not in the modifiable mode, call Init() first!"), DebugDiSTC);

    UpdateIterator it= AddSimplexToUpdate( t.GetDim(), &t, false);
    RemoteDataCL& rd= it->second.GetRemoteData();
    for (RemoteDataCL::ProcList_const_iterator pit= rd.GetProcListBegin(), pend= rd.GetProcListEnd(); pit!= pend; ++pit)
        it->second.AddProc( pit->proc, pit->prio, /*changeLocalPrio*/SimplexTransferInfoCL::keep);
}

ModifyCL::UpdateIterator ModifyCL::AddSimplexToUpdate( int dim, const TransferableCL* t, bool updateSubs)
{
    UpdateListT& simplices= entsToUpdt_[dim];
    // map::insert only inserts, if key is not found in the map, otherwise iterator of the found element is returned
    return simplices.insert( std::make_pair( t, SimplexTransferInfoCL( InfoCL::Instance().GetRemoteData( *t), updateSubs))).first;
}

void ModifyCL::CreateUpdateList()
{
    UpdateListT& tetras= entsToUpdt_[3];
    // insert all subsimplices of tetras to update in respective update lists
    for (UpdateIterator it= tetras.begin(), end= tetras.end(); it!=end; ++it)
        if (it->second.UpdateSubs()) {
            TetraCL* t;
            simplex_cast<TetraCL>( *it->first, t);

            for ( Uint i=0; i<NumVertsC; ++i)
                AddSimplexToUpdate( 0, t->GetVertex(i), false);
            for ( Uint i=0; i<NumEdgesC; ++i)
                AddSimplexToUpdate( 1, t->GetEdge(i), false);
            for ( Uint i=0; i<NumFacesC; ++i)
                AddSimplexToUpdate( 2, t->GetFace(i), false);
        }
    // communicate simplices to update
    CommToUpdateHandlerCL commToUpdate(*this);
    commToUpdate.Call();
}

bool ModifyCL::AssignPostProcs()
/// returns whether there are entities to be removed.
{
    const ProcSetT* postProcs;
    UpdateListT &verts= entsToUpdt_[0],
        &edges= entsToUpdt_[1],
        &faces= entsToUpdt_[2],
        &tetras= entsToUpdt_[3];
    int me= ProcCL::MyRank();
    RemoteDataListCL& remoteTetras= InfoCL::Instance().GetRemoteList<TetraCL>();

    for (RemoteDataListCL::iterator it= remoteTetras.begin(), end= remoteTetras.end(); it!=end; ++it) {
        TetraCL* t;
        simplex_cast<TetraCL>( it->second.GetLocalObject(), t);
        const UpdateIterator tit= tetras.find( it->second.GetLocalObjectPtr()),
                tend= tetras.end();

        // get tetra's newProcs set
        ProcSetT tmpProcSet;
        if (tit!=tend) { // tetra to be updated
            postProcs= &(tit->second.GetPostProcs());
        } else { // tetra will not be moved  ->  newProcs = (me,tetraPrio)
            tmpProcSet.insert( std::make_pair( me, t->GetPrio()));
            postProcs= &tmpProcSet;
        }

        // update proc sets of tetra's sub-simplices
        for ( Uint i=0; i<NumVertsC; ++i) {
            const UpdateIterator sit= verts.find( t->GetVertex(i)),
                    send= verts.end();
            if (sit!=send) // simplex to be updated
                sit->second.AddProcSet( *postProcs);
            }
        for ( Uint i=0; i<NumEdgesC; ++i) {
            const UpdateIterator sit= edges.find( t->GetEdge(i)),
                    send= edges.end();
            if (sit!=send) // simplex to be updated
                sit->second.AddProcSet( *postProcs);
            }
        for ( Uint i=0; i<NumFacesC; ++i) {
            const UpdateIterator sit= faces.find( t->GetFace(i)),
                    send= faces.end();
            if (sit!=send) // simplex to be updated
                sit->second.AddProcSet( *postProcs);
        }
    }
    // merge proc/prio lists. Note: local prio wins as in SimplexTransferInfoCL::AddProc (otherwise ChangePrio does not work properly)
    MergeProcListHandlerCL mergeProcList( *this);
    mergeProcList.Call();

    // now simplices to remove and proc/prio lists for sending can be determined
    bool foundDel= false; // is there any entity to be deleted?
    for (int dim=0; dim<4; ++dim)
        for (UpdateIterator it= UpdateBegin(dim), end= UpdateEnd(dim); it!=end; ++it) {
            SimplexTransferInfoCL& info= it->second;
            if (!info.WillBeOnProc(me)) {
                foundDel= true;
                info.SetRemoveMark();
            }
            info.ComputeSendToProcs( dim==3); // tetras have to be treated in a special way
        }
    return foundDel;
}

void ModifyCL::UpdateRemoteData()
{
    const RemoteDataCL::LoadVecT& loadOfProc= InfoCL::Instance().GetLoadVector();
    for (int dim=0; dim<4; ++dim)
        for (UpdateListT::iterator it= entsToUpdt_[dim].begin(), end= entsToUpdt_[dim].end(); it!=end; ++it) {
            SimplexTransferInfoCL& sti= it->second;
            if (!sti.WillBeRemoved()) {
                // update remote data
                RemoteDataCL::ProcListT proclist( sti.GetPostProcs().begin(), sti.GetPostProcs().end());
                sti.GetRemoteData().SetProcList( proclist);
                sti.GetRemoteData().UpdateOwner(loadOfProc);
            }
        }
}

void ModifyCL::DeleteUnusedSimplices( bool del)
{
    // remove entries from remote data list
    for (int dim=3; dim>=0; --dim) {
        RemoteDataListCL& rdl= InfoCL::Instance().GetRemoteList( dim);
        for (UpdateListT::iterator it= entsToUpdt_[dim].begin(), end= entsToUpdt_[dim].end(); it!=end; ++it)
            if (it->second.WillBeRemoved()) {
                TransferableCL& t= it->second.GetRemoteData().GetLocalObject();
                t.SetRemoveMark();
                rdl.Unregister( t);
            }
    }
    if (!del) return;

    // now remove simplices from multigrid
    SimplexFactoryCL& sf= mg_.GetSimplexFactory();
    for (int lvl= mg_.GetLastLevel(); lvl>=0; --lvl) {
        sf.DestroyMarkedTetras(lvl);
        sf.DestroyMarkedVEFs(lvl);
    }
}



// T R A N S F E R  C L A S S
//---------------------------

void TransferCL::Init()
{
    base::Init();
}

void TransferCL::Finalize()
/// \todo Delete Sendbuffer as last operation in this function?
{
    Assert( modifiable_, DROPSErrCL("TransferCL::Finalize: Class is not in the modifiable mode, call Init() first!"), DebugDiSTC);
    // algorithm based on migration algorithm of FMDB
    CreateUpdateList();
    AssignPostProcs();
    FillSendBuffer();

    std::vector<ProcCL::RequestT> req; req.reserve(sendBuffer_.size());
    // non-blocking sends
    for (SendBufT::const_iterator it(sendBuffer_.begin()), end(sendBuffer_.end()); it!=end; ++it)
        req.push_back( it->second->Isend( it->first));

    InfoCL::Instance().PreMultiGridMod();
    Receive();

    // Wait until each message has left me, before deleting the send buffer
    ProcCL::WaitAll(req);

    // All messages are delivered, so free all allocated memory
    for (SendBufT::iterator it= sendBuffer_.begin(), end= sendBuffer_.end(); it!=end; ++it)
        delete it->second;
    sendBuffer_.clear();
    req.clear();

    DeleteUnusedSimplices( del_);
    InfoCL::Instance().PostMultiGridMod();
    UpdateOwners();
    // Now we are ready to remove the update lists
    delete[] entsToUpdt_;
    entsToUpdt_= 0;
    modifiable_= false;
}

void TransferCL::Transfer( const TetraCL& t, int toProc, Priority prio, bool del)
{
    Assert( modifiable_, DROPSErrCL("TransferCL::MarkForTransfer: Class is not in the modifiable mode, call Init() first!"), DebugDiSTC);

    UpdateIterator it= AddSimplexToUpdate( 3, &t, /*updateSubs*/true);
    SimplexTransferInfoCL& info= it->second;
    /// add/merge proc/prio
    info.AddProc( toProc, prio);
    if (!del) { // keep local object, if not already there
        prio= info.GetRemoteData().GetLocalPrio();
        if (info.GetPostProc(prio)==-1) // prio not found
            info.AddProc( ProcCL::MyRank(), prio);
    } else
        info.RemoveProc( ProcCL::MyRank());

    // Code below was copied from ParMultiGridCL::Transfer. Please use ParMultiGridCL::Transfer directly!
//    if (t.IsRegularlyRef() && t.IsMaster() && prio==PrioMaster)
//        t.UnCommitRegRefMark();

}

struct TransferCL::CompByDescLevelCL {
    bool operator() (const UpdateEntryT& u, const UpdateEntryT& v)
    {
        return u.first->GetLevel() > v.first->GetLevel();
    }
};

TransferCL::SortedListT* TransferCL::SortUpdateTetras()
{
    const int dim= GetDim<TetraCL>();
    SortedListT* sortedTetras= new SortedListT( entsToUpdt_[dim].begin(), entsToUpdt_[dim].end());
    sortedTetras->sort( CompByDescLevelCL());
    return sortedTetras;
}

void TransferCL::UpdateSendRemoteData( const TransferableCL& t, SimplexTransferInfoCL& sti)
{
    RemoteDataCL& rd= sti.GetRemoteData();
    // update remote data
    RemoteDataCL::ProcListT proclist( sti.GetPostProcs().begin(), sti.GetPostProcs().end());
    const bool isTetra= t.GetDim()==GetDim<TetraCL>();
    if (sti.GetBroadcaster()==ProcCL::MyRank() || isTetra) {
        // write simplex and remote data to all procs in SendToProcs
        SendStreamCL buf( binary_);
        buf << t << proclist;
        for (ProcSetT::const_iterator sit= sti.GetSendToProcs().begin(), send= sti.GetSendToProcs().end(); sit!=send; ++sit) {
            (*sendBuffer_[sit->first]).write( buf.begin(), buf.cur() - buf.begin());
            if (isTetra)
                (*sendBuffer_[sit->first]) << int(rd.GetLocalPrio());
        }
    }
    if (!sti.WillBeRemoved())
        rd.SetProcList( proclist);
}

void TransferCL::SendAddedData( const TransferableCL& t, SimplexTransferInfoCL& sti)
{
    // send to every proc except me
    ProcSetT sendTo= sti.GetPostProcs();
    sendTo.erase( ProcCL::MyRank());
    // write geom id and DOF data to all procs in sendTo
    SendStreamCL buf( binary_);
    t.Unknowns.Pack( buf << t.GetGID(), t);
    for (ProcSetT::const_iterator sit= sendTo.begin(), send= sendTo.end(); sit!=send; ++sit) {
        (*sendBuffer_[sit->first]).write( buf.begin(), buf.cur() - buf.begin());
    }
}

void TransferCL::UpdateOwners()
{
    const RemoteDataCL::LoadVecT& loadOfProc= InfoCL::Instance().GetLoadVector();
    InfoCL& info= InfoCL::Instance();
    for (int dim=0; dim<4; ++dim)
        for (RemoteDataListCL::iterator it= info.GetRemoteList(dim).begin(), end= info.GetRemoteList(dim).end(); it!=end; ++it)
            it->second.UpdateOwner(loadOfProc);
}

/** Generate for each process that gets some data a buffer containing all
    information to be sent.
*/
void TransferCL::FillSendBuffer()
{
    int me= ProcCL::MyRank();
    typedef std::map<int,SArrayCL<size_t,5> > msgCountT;
    msgCountT numMsg; // number of simplices of dimension 0,..,3 and number of added data to be sent to certain proc
    // first count number of messages to be sent
    for (int dim=0; dim<4; ++dim)
        for (UpdateListT::const_iterator it= entsToUpdt_[dim].begin(), end= entsToUpdt_[dim].end(); it!=end; ++it) {
            const SimplexTransferInfoCL& sti= it->second;
            if (sti.GetBroadcaster()==me || dim==GetDim<TetraCL>()) {
                for (ProcSetT::const_iterator pit= sti.GetSendToProcs().begin(), pend= sti.GetSendToProcs().end(); pit!=pend; ++pit)
                    numMsg[pit->first][dim]++;
            }
            for (ProcSetT::const_iterator pit= sti.GetPostProcs().begin(), pend= sti.GetPostProcs().end(); pit!=pend; ++pit)
                if (pit->first != me)
                    numMsg[pit->first][4]++; // count number of added data
        }
    // allocate send buffers
    for (msgCountT::iterator it= numMsg.begin(), end= numMsg.end(); it!=end; ++it) {
        SendStreamCL* buf= new SendStreamCL( binary_);
        sendBuffer_[it->first]= buf;
        // write number of verts/edges/faces/tetras/addedData to be sent
        for (int dim=0; dim<5; ++dim)
            (*buf) << it->second[dim];
    }
    // sort tetras to update by descending level (so that children are received first to enable parents to set the Children_ array)
    SortedListT* sortedTetras= SortUpdateTetras();
    // now update remote data and write simplices and corresponding remote data to send buffers. First for verts, edges, faces...
    for (int dim=0; dim<3; ++dim)
        for (UpdateListT::iterator it= entsToUpdt_[dim].begin(), end= entsToUpdt_[dim].end(); it!=end; ++it)
            UpdateSendRemoteData( *it->first, it->second);
    // ...now for tetra in sorted order
    for (SortedListT::iterator it= sortedTetras->begin(), end= sortedTetras->end(); it!=end; ++it)
        UpdateSendRemoteData( *it->first, it->second);
    delete sortedTetras;

    // now send added data (e.g., numerical data)
    for (int dim=0; dim<4; ++dim)
        for (UpdateListT::iterator it= entsToUpdt_[dim].begin(), end= entsToUpdt_[dim].end(); it!=end; ++it)
            SendAddedData( *it->first, it->second);
}

void TransferCL::Receive()
{
    // compute send pattern (procs I send to get a "1")
    std::valarray<int> send(0, ProcCL::Size());
    for ( SendBufT::const_iterator it(sendBuffer_.begin()), end(sendBuffer_.end()); it!=end; ++it){
        send[ it->first]= 1;
    }
    // get the receive pattern (procs I receive from have a nonzero entry)
    std::valarray<int> IreceiveFrom= ProcCL::Alltoall( send);
    // get the number of procs to receive from
    int numRecvFrom= 0;
    for ( int p=0; p < ProcCL::Size(); ++p)
        if ( IreceiveFrom[p]>0)
            ++numRecvFrom;
    // receive from other processes
    std::vector<RecvStreamCL*> recv( numRecvFrom);
    std::vector<SArrayCL<size_t,5> >   numMsg( numRecvFrom);
    for ( int p=0, i=0; p < ProcCL::Size(); ++p)
        if ( IreceiveFrom[p]>0) {
            recv[i]= new RecvStreamCL( binary_);
            recv[i]->Recv( p);
            // receive number of verts/edges/faces/tetras
            for (int dim=0; dim<5; ++dim)
                (*recv[i]) >> numMsg[i][dim];
            ++i;
        }
    // Receive vertices
    for (int i=0; i<numRecvFrom; ++i)
        ReceiveSimplices<VertexCL>( *recv[i], numMsg[i][0]);
    // Receive edges
    for (int i=0; i<numRecvFrom; ++i)
        ReceiveSimplices<EdgeCL>  ( *recv[i], numMsg[i][1]);
    // Receive faces
    for (int i=0; i<numRecvFrom; ++i)
        ReceiveSimplices<FaceCL>  ( *recv[i], numMsg[i][2]);
    // receive tetras
    for (int i=0; i<numRecvFrom; ++i)
        ReceiveSimplices<TetraCL> ( *recv[i], numMsg[i][3]);
    // receive added data
    for (int i=0; i<numRecvFrom; ++i) {
        ReceiveAddedData ( *recv[i], numMsg[i][4]);
        // current receive stream not needed anymore
        delete recv[i];
    }
}

void TransferCL::ReceiveAddedData( RecvStreamCL& recvstream, size_t num)
{
    GeomIdCL gid;
    for (size_t i=0; i<num; ++i) {
        // read gid and added data
        recvstream >> gid;
        TransferableCL& t= InfoCL::Instance().GetRemoteData(gid).GetLocalObject();
        t.Unknowns.UnPack( recvstream, t);
    }
}

template <>
VertexCL& TransferCL::CreateSimplex<VertexCL>( const VertexCL& s, const RemoteDataCL::ProcListT& pl)
{
    return mg_.GetSimplexFactory().MakeCopy( s, pl);
}

template <>
EdgeCL& TransferCL::CreateSimplex<EdgeCL>    ( const EdgeCL& s,   const RemoteDataCL::ProcListT& pl)
{
    return mg_.GetSimplexFactory().MakeCopy( s, pl);
}

template <>
FaceCL& TransferCL::CreateSimplex<FaceCL>    ( const FaceCL& s,   const RemoteDataCL::ProcListT& pl)
{
    return mg_.GetSimplexFactory().MakeCopy( s, pl);
}

template <>
TetraCL& TransferCL::CreateSimplex<TetraCL>  ( const TetraCL& s,  const RemoteDataCL::ProcListT& pl)
{
    TetraCL& t= mg_.GetSimplexFactory().MakeCopy( s, pl);
    // now link to faces
    for ( Uint i=0; i<NumFacesC; ++i)
        const_cast<FaceCL*>(t.GetFace(i))->LinkTetra(&t);
    return t;
}

template <>
void DiST::TransferCL::ReceiveSimplices<TetraCL>( DiST::RecvStreamCL& recvstream, size_t num)
{
    for ( size_t i=0; i<num; ++i) {
        TetraCL stmp;                                   // temporary to receive a single simplex
        DiST::RemoteDataCL::ProcListT procList; // temporary to receive process list

        // receive simplex
        stmp.UnPack( recvstream);
        Assert( stmp.GetDim()==GetDim<TetraCL>(), DROPSErrCL("Mismatch in dimension of a received simplex!"), DebugDiSTC);
        // receive proc/prio list
        recvstream >> procList;
        int formerPrio;
        recvstream >> formerPrio;
        bool modMFR= true;
        TetraCL* tp= 0;
        if (!InfoCL::Instance().Exists(stmp.GetGID())) {
            tp= &CreateSimplex<TetraCL>( stmp, procList); // creates simplex and remote data list entry
            if (formerPrio==PrioGhost && procList.size()==1) // there will be a merge afterwards by master received from another proc who will take care of the edge MFRs
            modMFR= false;
        } else { // merge with existing tetra
            tp= InfoCL::Instance().GetTetra( stmp.GetGID());
            tp->Merge(stmp);
            modMFR= formerPrio==PrioMaster;
        }

        // TODO: hier experimentell, sollte per Handler ausgelagert werden (aus ParMultiGridCL::HandlerTObjMkCons).
        // Note: tetra can be received twice (Ma/Gh copy merged on a third proc), avoid double increment of edge MFRs.
        if (tp->IsRegularlyRef() && tp->IsMaster() && modMFR)
        {
            tp->CommitRegRefMark();
            // nun wird auf allen Edges _AccMFR:=_MFR gesetzt, um Unkonsistenzen bei vorher verteilt
            // und nun nur noch lokal gespeicherten Edges zu vermeiden. Ein abschliessenden
            // AccumulateMFR in XferEnd() setzt auf den verteilt gespeicherten Edges dann die
            // richtigen _AccMFR-Werte.
            for (TetraCL::const_EdgePIterator it= tp->GetEdgesBegin(), end= tp->GetEdgesEnd(); it!=end; ++it)
                (*it)->SetAccMFR( (*it)->GetMFR());
        }

    }
}

// I N F O  C L A S S
// ------------------
InfoCL* InfoCL::instance_= 0;

void InfoCL::PreMultiGridMod()
{
    mg_->PrepareModify();
}

void InfoCL::PostMultiGridMod()
{
    mg_->FinalizeModify();
    mg_->ClearTriangCache();
}

VertexCL* InfoCL::GetVertex( const GeomIdCL& gid) const
{
    VertexCL* v=0;
    simplex_cast<VertexCL>( GetRemoteData(gid).GetLocalObject(), v);
    return v;
}

EdgeCL* InfoCL::GetEdge(   const GeomIdCL& gid) const
{
    EdgeCL* v=0;
    simplex_cast<EdgeCL>( GetRemoteData(gid).GetLocalObject(), v);
    return v;
}

FaceCL* InfoCL::GetFace(   const GeomIdCL& gid) const
{
    FaceCL* v=0;
    simplex_cast<FaceCL>( GetRemoteData(gid).GetLocalObject(), v);
    return v;
}

TetraCL* InfoCL::GetTetra(  const GeomIdCL& gid) const
{
    TetraCL* v=0;
    simplex_cast<TetraCL>( GetRemoteData(gid).GetLocalObject(), v);
    return v;
}

/** Also checks the multigrid for sanity*/
bool InfoCL::IsSane( std::ostream& os) const
{
    bool sane=true;
    // Check remote data lists
    for ( size_t i=0; i<4; ++i) {
        sane= sane && remoteData_[i].IsSane( os);
    }
//    sane= sane && mg_->IsSane(os);
    return sane;
}

void InfoCL::DebugInfo( std::ostream& os) const
{
    for ( size_t i=0; i<4; ++i) {
        os << "dim = " << i << ":\t";
        remoteData_[i].DebugInfo( os);
    }
}

void InfoCL::SizeInfo( std::ostream& os) const
{
    os << "Remote data on proc " << ProcCL::MyRank() << ":";
    for ( size_t i=0; i<4; ++i){
        os << '\t' << remoteData_[i].size() << " of dim = " << i;
    }
    os << std::endl;
}

/** Show debug information of simplex given by its geometric id.
    \param gid geometric id of simplex
    \param os output stream
*/
void InfoCL::ShowSimplex( const GeomIdCL& gid, std::ostream& os) const
{
    if (Exists(gid)) {
        GetRemoteData(gid).GetLocalObject().DebugInfo( os);
        return;
    }
}

}   // end of namespace DiST
}   // end of namespace DROPS
