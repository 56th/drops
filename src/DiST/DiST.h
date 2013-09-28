/// \file DiST.h
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

#ifndef DROPS_DIST_H
#define DROPS_DIST_H

#include "misc/container.h"
#include "misc/utils.h"
#include "num/unknowns.h"
#include "DiST/mpistream.h"
#include "DiST/geomid.h"
#include "DiST/remotedata.h"

#include <list>
#include <set>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>

#ifdef _PAR
#  include "parallel/parallel.h"
#endif

namespace DROPS
{
// fwd decl
// Classes for simplices in the triangulation
class VertexCL;
class EdgeCL;
class FaceCL;
class TetraCL;
class SimplexFactoryCL;
class MultiGridCL;


/// \brief Namespace for all functions constituting a parallel multigrid
/** Everything in this namespace acts as an interface to DROPS. */
namespace DiST{

class SimplexTransferInfoCL
{
  public:
    typedef std::map<int,Priority> ProcSetT;
    /// enum to control how local priority entries are changed by AddProc
    enum UpdatePolicyE {
        keep, overwrite, merge
    };
  private:
    RemoteDataCL& rd_;       ///< remote data of simplex on local proc
    ProcSetT postProcs_;             ///< set of procs/prios where object will be after transfer
    ProcSetT procsToSend_;           ///< set of procs/prios where object has to be newly created
    bool RemoveMark_;                ///< whether simplex should be deleted locally after transfer
    bool UpdateSubs_;                ///< whether sumbsimplices should be updated (e.g. true for MarkForTransfer of a tetra, but false for ChangePrio of a tetra)
    int  proc_bc_;                   ///< proc rank of broadcaster

    /// \brief Merge two priorities
    inline Priority merge_prio( Priority a, Priority b) const { return std::max(a,b); }

  public:
    SimplexTransferInfoCL( RemoteDataCL& rd, bool updateSubs)
        : rd_(rd), RemoveMark_(false), UpdateSubs_(updateSubs), proc_bc_(rd.GetOwnerProc()) {}
    SimplexTransferInfoCL( const SimplexTransferInfoCL& sti)
        : rd_(sti.rd_), postProcs_(sti.postProcs_), procsToSend_(sti.procsToSend_), RemoveMark_(sti.RemoveMark_), UpdateSubs_(sti.UpdateSubs_), proc_bc_(sti.proc_bc_) {}

    /// \brief add/merge new \a proc/\a prio entry into set of post procs.
    /// If \a proc is the current process, the priority of an existing entry will not be changed unless \a changeLocalPrio==true.
    void AddProc( int proc, Priority prio, UpdatePolicyE changeLocalPrio=keep);
    /// \brief Adds new proc/prio entries without changing existing entries.
    void AddProcSet( const ProcSetT& procs);
    /// \brief Remove \a proc entry from set of post procs
    void RemoveProc( int proc) { postProcs_.erase(proc); }
    /// return whether simplex will be on proc \a p after transfer
    bool WillBeOnProc( int p) const { return postProcs_.find(p) != postProcs_.end(); }
    /// return process rank where simplex will have a certain \a prio after transfer (return -1 if such a rank doesn't exist)
    int GetPostProc( Priority prio) const;
    /// return set of procs/prios after transfer
    const ProcSetT& GetPostProcs() const { return postProcs_; }
    /// return set of procs/prios where simplex has to be sent to
    const ProcSetT& GetSendToProcs() const { return procsToSend_; }
    /// compute set of procs/prios where simplex has to be sent to (call after post proc list is complete)
    void ComputeSendToProcs( bool tetra);
    /// mark simplex for deletion
    void SetRemoveMark()       { RemoveMark_= true; }
    /// return whether simplex will be removed
    bool WillBeRemoved() const { return RemoveMark_; }
    /// return whether simplex will update its sub simplices (only makes sense for tetrahedron)
    bool UpdateSubs()    const { return UpdateSubs_; }
    /// return rank of broadcaster
    int  GetBroadcaster() const { return proc_bc_; }
    /// return simplex' remote data
    const RemoteDataCL& GetRemoteData() const { return rd_; }
    /// return simplex' remote data
    RemoteDataCL&       GetRemoteData()       { return rd_; }
};

/// \brief Enables communication across certain process boundaries.
class InterfaceCL
/// A process boundary includes all geometric entities which are shared by at least
/// two processes. Communication can be one-way or two-way. Basically, for each
/// entity specified by a level, dimension and a priority (in to_) a function gather
/// is called. The gathered data are collected by the owner of the entity. Afterwards,
/// the owner informs all copies of the entity to call a scatter function.
///
/// DiST::InterfaceCL is a replacement for DDD interfaces.
/// \todo
/// - Implementation of iterators, cf TriangCL
/// - Interface for grid level (or all levels) and triangulation level (?)
/// - When to call the possible expensive function SetupCommunicationStructure?
///   Maybe, an observer pattern would be nice here.
{
  public:
    typedef WholeRemoteDataIteratorCL iterator;    ///< type for iterator over elements in the interface
    typedef iterator::DimListT DimListT;                   ///< type for storing the involved simplices
    typedef iterator::GIDIteratorT GIDIteratorT;           ///< type for iterator over GID's
    typedef std::map<int, SendStreamCL> SendListT; ///< type for storing data to be sent (proc -> data)
    typedef std::map<int, RecvStreamCL> RecvListT; ///< type for receiving data (proc -> data)
    typedef std::set<int> ProcSetT;                        ///< type for a set of proccessor numbers

    /// \brief Helper types for ExchangeData and the function collect_streams in Dist.cpp
    ///@{
    class MessagesCL {
      private:
        size_t            numData;  ///< number of messages
        std::vector<char> messages; ///< concatenation of the messages as sequence of char

      public:
        MessagesCL () : numData( 0) {}
        void append (const char* begin, const char* end)
            { messages.insert( messages.end(), begin, end); ++numData; }
        friend MPIostreamCL& operator<< (MPIostreamCL& os, const InterfaceCL::MessagesCL& msg);
    };
    typedef DROPS_STD_UNORDERED_MAP<GeomIdCL, MessagesCL, Hashing > CollectDataT;
    ///@}

  private:
    enum CommPhase {            ///< Communication phases
        bothPhases,             ///< owner and copies gather and scatter data, i.e., copies -> owner -> copies
        toowner,                ///< copies send to owners, i.e., copies -> owner
        fromowner               ///< owner send to copies, i.e., owner -> copies
    };

    RecvListT recvbuf_;         ///< received data (after call of function Communicate())
    SendListT sendbuf_;         ///< data to be sent (filled in GatherData)
    SendStreamCL    loc_send_toowner_; ///< used in Perform, phase==toowner, to avoid a copy of the local message
    RefMPIistreamCL loc_recv_toowner_; ///< used in Perform, phase==toowner, to avoid a copy of the local message

    PrioListT from_;            ///< list of priority on sender side, the interface operates on
    PrioListT to_;              ///< list of priority on receiver side, the interface operates on

    bool      binary_;          ///< transfer data binary or in ASCII
    ProcSetT  ownerRecvFrom_;   ///< list of processes, the owner receives data from
    ProcSetT  ownerSendTo_;     ///< list of processes, the owner sends data to
    ProcSetT  IRecvFromOwners_; ///< list of owner processes, I have to receive data from

    iterator  begin_from_,      ///< iterators defining sequence of interface elements
              begin_to_,
              begin_,
              end_;

    /// \brief Call the gather handler for each entity covered by the iterators [begin, end).
    template <typename HandlerT>
    void GatherData( HandlerT&, const iterator& begin, const iterator& end, CommPhase phase);
    /// \brief MPI Isend of the streams in sendbuf.
    void SendData (SendListT& sendbuf, std::vector<ProcCL::RequestT>& req, int tag);
    /// \brief Phase (1) and (2) of ExchangeData: The owner acquires the information.
    void to_owner (std::vector<ProcCL::RequestT>&, CollectDataT&);
    /// \brief Phase (4) and (5) of ExchangeData: The owner distributes the accumulated information.
    void from_owner (std::vector<ProcCL::RequestT>&, SendListT&);
    /// \brief For a stream [GID, tail), call handler(GID, numData, tail).
    template <typename HandlerT, typename IStreamT>
    bool ScatterData( HandlerT& handler, IStreamT& recv);
    /// \brief Call the scatter handler for each entity whose GID is given in the istream_.
    template <typename HandlerT>
    bool ScatterData( HandlerT&);
    /// \brief Handles the communication between the processes
    void ExchangeData( CommPhase phase);
    /// \brief set up the ownerRecvFrom_ and IRecvFromOwners_
    void SetupCommunicationStructure();
    /// \brief do the specified communication
    template <typename HandlerT>
    bool Perform( HandlerT&, CommPhase);

  public:
    /// \brief Construct an interface for communication across process boundaries
    /** Therefore, the simplices can be grouped according to the parameters:
        \param lvl    the level of the simplices
        \param from   priorities on "gather side"
        \param to     priorities on "scatter side"
        \param dims   dimensions of the simplices
        \param binary transfer of data in binary or ASCII mode
    */
    InterfaceCL( const LevelListCL& lvl, const PrioListT& from, const PrioListT& to,
                 const DimListT& dims, const bool dist=true, bool binary= use_binaryMPIstreams)
        : from_(from), to_(to), binary_(binary),
          begin_from_( dims, lvl, from, dist),
          begin_to_( dims, lvl, to, dist), begin_( dims, lvl, /*all prios*/PrioListT(), dist), end_( begin_from_.GetEnd()) {}

    InterfaceCL( GIDIteratorT begin, GIDIteratorT end, bool binary= use_binaryMPIstreams)
        : binary_(binary), begin_from_( begin, end), begin_to_( begin_from_), begin_(begin_from_), end_( begin_from_.GetEnd()) {}

    /// \name ExecuteLocal, functions for calling an execute handler for all
    ///       entities covered by the interface
    /// \tparam ExecuteHandlerT entities h of this type must provide a function header of <br>
    ///         bool h( TransferableCL& t)
    /// \return Local accumulated AND for all calls of the ExecuteHandlers (without global reduction)
    //@{
    /// \brief Call \a handler for each entity covered by this interface
    template <typename ExecuteHandlerT>
    bool ExecuteLocal( ExecuteHandlerT& handler);

    /// \brief Call \a handler for each entity covered by the iterators [begin, end)
    /// \tparam IteratorT each iterator it of this type must provide:
    ///         it->first :  GeomIdCL
    ///         it->second : RemoteDataCL
    template <typename ExecuteHandlerT, typename IteratorT>
    bool ExecuteLocal( ExecuteHandlerT& handler, const IteratorT& begin, const IteratorT& end);
    //@}

    /// \name Do The interface communication
    /// \param HandlerT A handler h is used to gather and scatter data on
    ///     TransferableCL's. Therefore, the handler must provide the members
    ///    Gather( DiST::TransferableCL&, DiST::SendStreamCL&)
    ///    Scatter( DiST::TransferableCL&, const size_t&, DiST::RecvStreamCL&);
    //@{
    /// \brief Do the interface communication on the interface specified by the
    ///   the constructor
    template <typename HandlerT>
    bool Communicate( HandlerT&);

    /// \brief Do the interface communication from non-owners to owners
    template <typename HandlerT>
    bool InformOwners( HandlerT&);

    /// \brief Do the interface communication from owners to non-owners
    template <typename HandlerT>
    bool InformCopies( HandlerT&);
    //@}
};

/// \brief Enables modifications of the distributed data such as priority changes and
///        unregistering/deleting objects.
class ModifyCL
/// All modifications should be embraced by calls to Init() and Finalize().
/// ModifyCL is a replacement for the DDD Prio module.
{
  protected:
    // handlers for interface comm
    class MergeProcListHandlerCL;
    class CommToUpdateHandlerCL;
    /// \brief type for remembering which simplices need to be updated, plus for each some temporary information needed for the transfer
    typedef std::map<const TransferableCL*, SimplexTransferInfoCL> UpdateListT;
    typedef std::pair<const TransferableCL*,SimplexTransferInfoCL> UpdateEntryT;
    typedef UpdateListT::iterator                                          UpdateIterator;
    typedef SimplexTransferInfoCL::ProcSetT                        ProcSetT;

    bool          modifiable_;        ///< for checking if Init and Finalize are called
    bool          del_;               ///< physically remove simplices
    bool          binary_;            ///< transfer data in binary format
    MultiGridCL&  mg_;                ///< Multigrid
    UpdateListT*  entsToUpdt_;        ///< simplices to be updated

    /// \name iterator over simplices of dimension \a dim to be updated
    //@{
    UpdateIterator UpdateBegin(int dim) { return entsToUpdt_[dim].begin(); }
    UpdateIterator UpdateEnd  (int dim) { return entsToUpdt_[dim].end(); }
    //@}
    /// adds simplex of dimension \a dim to the respective update list
    UpdateIterator AddSimplexToUpdate( int dim, const TransferableCL*, bool updateSubs);

    /// \brief Create list of simplices to be updated (needs interface comm.)
    void CreateUpdateList();
    /// \brief Assign new procs (needs interface comm.) and determine which simplices will be deleted after transfer
    bool AssignPostProcs();
    /// \brief Update the remote data (incl. owner) to account for changed prios
    void UpdateRemoteData();
    /// \brief Delete all simplices that are not needed any more.
    void DeleteUnusedSimplices( bool del);

  public:
    /// \brief Constructor with a given multigrid (\a mg) and decision if communication should be done \a binary.
    /// For \a del = true all unused simplices are removed from the multigrid, otherwise the RemoveMark is set and the simplex is unregistered from the DiST module.
    ModifyCL( MultiGridCL& mg, bool del= true, bool binary= use_binaryMPIstreams)
        : modifiable_(false), del_(del), binary_( binary), mg_(mg), entsToUpdt_(0) {}

    /// \brief Call Init() before any modifications
    void Init();
    /// \brief Call Finalize() after any modifications (initiates the communication)
    void Finalize();
    /// \brief Change the priority of a locally stored object
    void ChangePrio( const TransferableCL&, Priority);
    /// \brief Consider a local simplex for deletion.
    void Delete( const TransferableCL&);
    /// \brief Keep a local simplex (even though all adjacent tetras have called "Delete()").
    void Keep( const TransferableCL&);
};

/// \brief Enables transfer of objects.
class TransferCL : public ModifyCL
/// All transfer operations should be embraced by calls to Init() and Finalize().
/// TransferCL is a replacement for the DDD Xfer module.
{
  public:
    typedef ModifyCL                                base;
    typedef TransferCL                              self;

  private:
    /// \brief type for sending information to other processes
    typedef std::map<int, SendStreamCL*> SendBufT;
    /// \brief type for sorted tetras to update
    typedef std::list<UpdateEntryT>              SortedListT;
    // comparison functor to sort by descending level
    struct CompByDescLevelCL;

  private:
    SendBufT      sendBuffer_;        ///< all information to be sent

    /// create list of tetras to update, sorted by descending order
    SortedListT* SortUpdateTetras();
    /// Update and send remote data of a given simplex.
    void UpdateSendRemoteData( const TransferableCL&, SimplexTransferInfoCL&);
    /// Send added data of a given simplex.
    void SendAddedData( const TransferableCL&, SimplexTransferInfoCL&);
    /// \brief Update ownership (in remote data) for all registered objects
    void UpdateOwners();
    /// \brief allocate and fill send buffers, update remote data
    void FillSendBuffer();
    /// \brief Receive and create simplices and remote data
    void Receive();
    /// \brief Receive vertices, edges, faces, tetras (helper for Receive())
    template <typename SimplexT>
    void ReceiveSimplices( RecvStreamCL&, size_t);
    /// \brief Receive added data
    void ReceiveAddedData( RecvStreamCL&, size_t);
    /// \brief Create simplex using the SimplexFactoryCL
    template <typename SimplexT>
    SimplexT& CreateSimplex( const SimplexT&, const RemoteDataCL::ProcListT&);

  public:
    /// \brief Constructor with a given multigrid (\a mg) and decision if the transfer should be done \a binary.
    /// For \a del = true all unused simplices are removed from the multigrid, otherwise the RemoveMark is set and the simplex is unregistered from the DiST module.
    TransferCL( MultiGridCL& mg, bool del= true, bool binary= use_binaryMPIstreams)
        : base( mg, del, binary) {}
    /// \brief To be called before marking tetrahedra for transfer
    void Init();
    /// \brief To be called after marking tetrahedra for transfer (initiates the communication)
    void Finalize();
    /// \brief Mark a tetrahedron for transfer
    void Transfer( const TetraCL&, int toProc, Priority prio, bool del= true);
};


// fwd declaration
class InfoInitCL;

/// \brief The entry point for the DiST module storing a RemoteDataListCL for each simplex type.
class InfoCL
/// Implemented using the singleton pattern.
{
private:
	// all DiST modules are friends to allow access to RemoteDataListCL objects
    friend class InfoInitCL;
    friend class ModifyCL;
    friend class TransferCL;
    friend class OldTransferCL;
    friend class InterfaceCL;
    friend class MultiGridCL;
    friend class SimplexFactoryCL;

  private:
    /// \brief RemoteDataListCL
    /** For each simplex type, a single list is stored where list i contains all simplices of dimension i. */
    RemoteDataListCL remoteData_[4];
    RemoteDataCL::LoadVecT loadOfProc_;

    MultiGridCL* mg_; // do we need it?

    /// \name singleton pattern
    //@{
    static InfoCL* instance_;
    InfoCL( MultiGridCL* mg) : loadOfProc_( ProcCL::Size()), mg_(mg) {}
    //@}

    InfoCL( const InfoCL&);           // copy ctr not implemented
    InfoCL& operator=(const InfoCL&); // assignment op not implemented

    void PreMultiGridMod();  ///< call before modification of multigrid
    void PostMultiGridMod(); ///< call after modification of multigrid

  public:
    /// \name Implementing the singleton pattern. So access to this class is only provided via the get instance functions.
    //@{
    static inline InfoCL& Instance( MultiGridCL* mg=0);
    static inline InfoCL* InstancePtr();
    static inline void Destroy() { if (instance_) delete instance_; }
    //@}

    /// \brief Access the RemoteDataCL of a given simplex
    const RemoteDataCL& GetRemoteData( const TransferableCL& t) const { return GetRemoteData( t.GetGID()); }
    /// \brief Access the RemoteDataCL of a given simplex (private?)
    RemoteDataCL& GetRemoteData( const TransferableCL& t) { return GetRemoteData( t.GetGID()); }
    inline const RemoteDataCL& GetRemoteData( const GeomIdCL&) const; ///< Get remote data by GID
    inline       RemoteDataCL& GetRemoteData( const GeomIdCL&);       ///< Get remote data by GID

    /// \brief Access the RemoteDataListCL of a given dimension
    const RemoteDataListCL& GetRemoteList( int dim) const { return remoteData_[dim]; }
    /// \brief Access the RemoteDataListCL of a given dimension
    RemoteDataListCL& GetRemoteList( int dim) { return remoteData_[dim]; }
    /// \brief Access the RemoteDataListCL of a given simplex type
    template <typename SimplexT>
    const RemoteDataListCL& GetRemoteList() const { return remoteData_[GetDim<SimplexT>()]; }
    /// \brief Access the RemoteDataListCL of a given simplex type (private?)
    template <typename SimplexT>
    RemoteDataListCL& GetRemoteList() { return remoteData_[GetDim<SimplexT>()]; }

    /// \brief Access loads of processes used to determine the ownership in the remote data.
    const RemoteDataCL::LoadVecT& GetLoadVector() const { return loadOfProc_; }

    /// \brief Is this GID registered?
    inline bool Exists( const GeomIdCL& gid) const { return gid.dim < 4 ? remoteData_[gid.dim].Exists(gid) : false; }

    /// \name Direct access to a simplex by its GID
    /// \todo Make these functions inline
    //@{
    VertexCL* GetVertex( const GeomIdCL&) const;
    EdgeCL*   GetEdge  ( const GeomIdCL&) const;
    FaceCL*   GetFace  ( const GeomIdCL&) const;
    TetraCL*  GetTetra ( const GeomIdCL&) const;
    //@}

    /// \name Debugging
    //@{
    bool IsSane( std::ostream&) const;
    void DebugInfo( std::ostream&) const;
    void SizeInfo( std::ostream&) const;
    void ShowSimplex( const GeomIdCL&, std::ostream&) const;     ///< Show the debug information of a simplex given by its geometric id
    //@}
};

}   // end of namespace DiST
}   // end of namespace DROPS

#include "DiST/DiST.tpp"
#endif
