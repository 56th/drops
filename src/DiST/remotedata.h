/// \file remotedata.h
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

#ifndef REMOTEDATA_H_
#define REMOTEDATA_H_

#include "parallel/parallel.h"
#include "DiST/geomid.h"
#include "DiST/mpistream.h"
#include "num/unknowns.h"

namespace DROPS
{

/// \enum Priority Priorities for distributed objects
/// \todo Put these priorities in the DiST namespace
enum Priority
{
    NoPrio          = -1,    ///< invalid, used for cases where no priority has been assigned
    PrioNeutral     = 0,     ///< no priority
    PrioKilledGhost = 1,     ///< prio of ghost tetras marked for removal
    PrioVGhost      = 2,     ///< for subs of overlapping master tetra (don't constitute a proc bnd)
    PrioGhost       = 3,     ///< for ghost tetras and subs only owned by ghost tetras (skipped by the public iterators)
    PrioMaster      = 4      ///< master copy of tetras and their subs
};

inline std::string PriorityToString( const Priority);

namespace DiST
{

/// \brief Priority list. By definition, an empty list contains all priorities.
class PrioListCL: public std::vector<Priority>
{
  public:
    PrioListCL() : std::vector<Priority>() {}
    bool contains( const Priority prio) const { return empty() || is_in( begin(), end(), prio); }
};

/// \brief type for priority list
typedef PrioListCL PrioListT;

/// \brief Level list. By definition, an empty list contains all levels.
class LevelListCL: public std::vector<Uint>
{
  public:
    explicit inline LevelListCL( Uint max_level); ///< creates level list with levels 0, ..., max_level
                    LevelListCL() {}              ///< creates empty list, containing all levels by definition

    bool contains( Uint level) const { return empty() || is_in( begin(), end(), level); }
};

// fwd declaration
class TransferableCL;
class ModifyCL;
class TransferCL;
class RemoteDataListCL;
class RemoteDataListIteratorCL;

/// error class for exception handling
class ErrorCL: public DROPSErrCL
{
  protected:
    GeomIdCL gid_;   ///< geometric id of associated simplex

  public:
    ErrorCL()
      : DROPSErrCL(), gid_(NoGID) {}
    ErrorCL( const std::string& mesg, const GeomIdCL& gid= NoGID)
      : DROPSErrCL(mesg), gid_(gid) {}
    /// print one piece of data
    template<class T>
    ErrorCL(const std::string& mesg, const T& data, const GeomIdCL& gid= NoGID);

    std::ostream& what  (std::ostream&) const;
    void handle() const;
};

/// \brief For each distributed entity, a list of process ranks (and corresponding priority) is stored.
class RemoteDataCL
/// A list of process ranks and corresponding priorities, which is to be stored for each distributed entity. This
/// list contains the calling process as first entry.
/// \todo Make it easier to iterate over ProcListT, maybe a macro as DROPS_FOR_TRIANG_CONST_TETRA
{
  public:
    struct ProcListEntryCL {
        int      proc;                                              ///< process rank
        Priority prio;                                              ///< priority on this process

        ProcListEntryCL()                                  : proc(-1), prio(PrioNeutral)   {}
        ProcListEntryCL( int rank, Priority priority)      : proc(rank), prio(priority)    {}
        ProcListEntryCL( const std::pair<int,Priority>& p) : proc(p.first), prio(p.second) {}

        bool operator== (const ProcListEntryCL& p) const { return p.proc==proc && p.prio==prio; }
        bool operator== (int rank) const { return rank==proc; }
    };

    typedef std::vector< ProcListEntryCL> ProcListT;                ///< list containing couplings of process ranks and priorities
    typedef ProcListT::iterator           ProcList_iterator;        ///< iterator for the list ProcListT
    typedef ProcListT::const_iterator     ProcList_const_iterator;  ///< constant iterator for the list ProcListT
    typedef std::valarray<Uint>           LoadVecT;                 ///< load vector for owner update

    friend class DiST::ModifyCL;
    friend class DiST::TransferCL;

  private:
    TransferableCL* localObj_;      ///< Pointer to local element
    ProcListT       procList_;      ///< processes where the object is stored + corresponding priorities
    int             owner_;         ///< process which is responsible for the communication

    // sort proc list such that the local entry is the first
    void MakeConsistent();
    // for Modify/TransferCL
    void SetProcList( const ProcListT& pl) { procList_= pl; MakeConsistent(); }

  public:
    RemoteDataCL() : localObj_(0), procList_(), owner_(-1) {}
    /// \brief only for debugging purpose!
    RemoteDataCL( const ProcListT& v) : localObj_(0), procList_(v), owner_(-1) { MakeConsistent(); }
    /// \brief Construct a regular RemoteDataCL. Note: owner has to be updated afterwards.
    RemoteDataCL( TransferableCL* tp, const ProcListT& v) : localObj_(tp), procList_(v), owner_(-1) { MakeConsistent(); }
    /// \brief Assumes, that only the calling process owns a copy
    RemoteDataCL( TransferableCL* tp, Priority prio) : localObj_(tp), procList_(1), owner_(ProcCL::MyRank()) { procList_[0]=ProcListEntryCL(ProcCL::MyRank(), prio); }

    /// \name Get iterator to begin and end
    //@{
    ProcList_iterator       GetProcListBegin()       { return procList_.begin(); }
    ProcList_iterator       GetProcListEnd  ()       { return procList_.end(); }
    ProcList_const_iterator GetProcListBegin() const { return procList_.begin(); }
    ProcList_const_iterator GetProcListEnd  () const { return procList_.end(); }
    //@}

          Uint            GetNumProcs( Priority prio=NoPrio) const;             ///< Get number of processes owning the distributed entity
    const TransferableCL& GetLocalObject()    const { return *localObj_; }      ///< Get a constant reference to the local stored entity
          TransferableCL& GetLocalObject()          { return *localObj_; }      ///< Get a reference to the local stored entity
    const TransferableCL* GetLocalObjectPtr() const { return localObj_; }       ///< Get a constant pointer to the local stored entity
          TransferableCL* GetLocalObjectPtr()       { return localObj_; }       ///< Get a pointer to the local stored entity

    void Identify( const TransferableCL& parent, const PrioListCL&);            ///< Make local entity a distributed object with same processes as the parent entity (note: no communication required)

    /// \brief Ask for the owner of the entity
    /** The owner definitively needs a copy with priority>=master!
        2 possible ways to decide: (a) compute consistently by GID or (b) use process with smallest work load (cf. FMDB)
    */
           int      GetOwnerProc() const { return owner_; }
           void     UpdateOwner( const LoadVecT& load);                         ///< Update owner: process with smallest load
           bool     AmIOwner()     const { return owner_==ProcCL::MyRank(); }   ///< Ask if the calling process is the owner
    inline Priority GetLocalPrio() const;                                       ///< Get the priority of the local stored entity
           Priority GetPrio(int rank) const;                                    ///< Get the priority of process with given \a rank. Returns NoPrio, if simplex is not stored on \a rank.
    inline bool     PrioExists(Priority prio) const;                            ///< Check that there exists a remote object with priority >= \a prio
           bool     IsLocal()      const { return GetNumProcs()==1; }           ///< Check if the simplex is only stored on this process
    inline bool     IsDistributed( Priority prio=NoPrio) const;                 ///< Check if a simplex is distributed and at least two remote objects have a priority>=prio
    inline bool     IsOnProcBnd()  const;                                       ///< Check if the simplex is located on a process boundary

    /// \brief Debugging
    //@{
    bool IsSane( std::ostream&) const;
    void DebugInfo( std::ostream&) const;
    //@}
};

/// \brief MPI-IO of ProcListT
///@{
MPIostreamCL& operator<< (MPIostreamCL& os, const RemoteDataCL::ProcListT& pl);
MPIistreamCL& operator>> (MPIistreamCL& is,       RemoteDataCL::ProcListT& pl);
///@}

/// \brief Container to store a RemoteDataCL object for each distributed geometric entity of type SimplexT.
/** This constitutes the main data structure for the DiST module. However, this is not designed to act as an
    interface to drops. Each access should be done via the functions in the DiST namespace.
    \todo
    - Iterators should be separately implemented. They should provide the functionality to iterate
      through specialized data, e.g., iterate over all simplices of level l with priority Master or Ghost.
*/
class RemoteDataListCL : public DROPS_STD_UNORDERED_MAP<GeomIdCL, RemoteDataCL, Hashing>
{
  public:
    typedef DROPS_STD_UNORDERED_MAP<GeomIdCL, RemoteDataCL, Hashing> base;
    typedef base::iterator                                           iterator;
    typedef base::const_iterator                                     const_iterator;

  private:
    friend class RemoteDataListIteratorCL;

    // handler for interface comm
    class DebugHandlerCL;

  public:
    using base::end;
    using base::begin;

    /** \brief Create an iterator that goes through all elements of the RemoteDataList, which
               have certain properties.
        \param[in] Levels Iterate through elements, whose level is contained in level list 'Levels'.
        \param[in] a      Iterate through elements, whose priority matches any of the priorities listed in 'a'.
        \param[in] dist   Iterate only through distributed elements.
    */
    inline RemoteDataListIteratorCL beginLvlPrioIt( const LevelListCL& Levels=LevelListCL(), const PrioListT& a=PrioListT(), bool dist=false);
    /// \brief end-Iterator
    inline RemoteDataListIteratorCL endLvlPrioIt();

    /// \brief Put a simplex with a given priority into the list
    inline void Register( TransferableCL&, Priority prio=PrioMaster);
    /// \brief Put a simplex with a given proc/prio list into the list
    inline void Register( TransferableCL&, const RemoteDataCL::ProcListT&);
    /// \brief Remove simplex from list
    inline void Unregister( TransferableCL&);
    /// \brief Is this GID registered?
    inline bool Exists( const GeomIdCL& gid) const { return find(gid)!=end(); }

    /// \name Checking and debugging functions
    //@{
    bool IsSane( std::ostream&) const;
    void DebugInfo( std::ostream&) const;
    void SizeInfo( std::ostream&) const;
    //@}
};

/// \brief Iterator for the class RemoteDataListCL which iterates through selected objects
class RemoteDataListIteratorCL
{
  public:
    typedef RemoteDataListCL  ListT;                        ///< type of the remote data list
    typedef RemoteDataListCL::iterator   iterator;          ///< iterator of the list
    typedef RemoteDataListCL::reference  reference;         ///< reference to the RemoteDataCL
#if defined(__INTEL_COMPILER) || (__GNUC__>=4)
    // std::tr1::unordered_map does not define a type 'pointer'
    typedef std::pair<const GeomIdCL, RemoteDataCL>* pointer; ///< pointer to the RemoteDataCL
#else
    typedef RemoteDataListCL::pointer    pointer;             ///< pointer to the RemoteDataCL
#endif

  protected:
    ListT*      list_;            ///< Pointer to an object of type RemoteDataListCL
    iterator    pos_;             ///< Shows the position of the iterator at the RemoteDataListCL object
    LevelListCL levels_;          ///< Iterate over the elements of RemoteDataListCL, which level is contained in 'levels_'
    PrioListT   prios_;           ///< that belong to any of the priorities contained in the vector prios_
    bool        distributed_;     ///< and is at least located at two processors.

    RemoteDataListIteratorCL( RemoteDataListCL* list) // only for use with derived classes
        : list_(list), pos_(list->begin()), prios_(PrioListT()), distributed_(false) {}

    /// \brief Check if the element "elem" is covered by this iterator class
    inline bool contains();

  public:

    /// \brief Constructor that creates a RemoteDataListIteratorCL object from an
    ///        iterator of RemoteDataListCL. (Used for end().)
    RemoteDataListIteratorCL( const iterator& it)
        : pos_(it), levels_(LevelListCL()), prios_(PrioListT()), distributed_(false) {}

    /// \brief Standard constructor, it receives a pointer to the container and the level and priority list.
    inline RemoteDataListIteratorCL( RemoteDataListCL*, const LevelListCL&, const PrioListT&, bool dist);

    /// \brief Copy constructor
    inline RemoteDataListIteratorCL(const RemoteDataListIteratorCL& Lit)
        : list_(Lit.list_), pos_(Lit.pos_), levels_(Lit.levels_), prios_(Lit.prios_), distributed_(Lit.distributed_) {}

    /// \brief return a reference to element of container, this iterator points to
    reference operator* () const { return  *pos_; }
    /// \brief return a pointer to element of container, this iterator points to
    pointer operator-> () const { return &*pos_; }

    /// \brief Go to the next element (prefix)
    inline RemoteDataListIteratorCL& operator ++ ();
    /// \brief Go to the next element (suffix)
    inline RemoteDataListIteratorCL  operator ++ (int);

    /// \brief assignment operator
    inline RemoteDataListIteratorCL& operator= ( const RemoteDataListIteratorCL&);

    /// \name Comparing two iterators
    //@{
    /// \brief Test, if two iterators are the same by comparing the iterators onto the container
    friend bool operator== ( const RemoteDataListIteratorCL& It0, const RemoteDataListIteratorCL& It1);
    /// \brief Test, if two iterators are not the same by comparing the iterators onto the container
    friend bool operator!= ( const RemoteDataListIteratorCL& It0, const RemoteDataListIteratorCL& It1);
    //@}
};

/** Just compare the position*/
inline bool operator== ( const RemoteDataListIteratorCL& It0, const RemoteDataListIteratorCL& It1)
{
    return It0.pos_==It1.pos_;
}

/** Just compare the position*/
inline bool operator!= ( const RemoteDataListIteratorCL& It0, const RemoteDataListIteratorCL& It1)
{
    return !(It0==It1);
}

class WholeRemoteDataIteratorCL: public RemoteDataListIteratorCL
/// \brief Iterate through selected objects of different dimensions stored in remote data lists.
///
/// Two different possibilities for definition:
/// 1) for a fixed level list, priority list, dimension list,
/// 2) for a given set of geometric id's (GID).
{
  public:
    typedef std::vector<Usint>      DimListT;     ///< type for storing the dimensions of the involved simplices
    typedef std::vector<GeomIdCL>   GIDVecT;      ///< type for storing geometric id's (GID)
    typedef GIDVecT::const_iterator GIDIteratorT; ///< type for GID iterators

  private:
    typedef RemoteDataListIteratorCL base;

    bool useGIDVec_;                              ///< use GID vector or levels/dims/prios for definition of iterator sequence
    GIDIteratorT itGID_;                          ///< iterator pointing to current GID
    GIDIteratorT endGID_;                         ///< end iterator for GID's
    DimListT dims_;                               ///< list of simplex dimensions to be considered
    size_t dimIdx_;                               ///< index of current dimension

    /// \brief Set iterator corresponding to the current GID. Returns dimension (0..4) of current GID.
    inline Usint FindGID();

  public:
    /// \brief Standard constructor defining sequence by dimensions, levels and priorities.
    WholeRemoteDataIteratorCL( const DimListT&, const LevelListCL&, const PrioListT&, bool dist);
    /// \brief Standard constructor defining sequence by a vector of GID's
    WholeRemoteDataIteratorCL( GIDIteratorT begin, GIDIteratorT end);
    /// \brief Constructor that creates a WholeRemoteDataIteratorCL object from an
    ///        iterator of RemoteDataListCL. (Used for end().)
    WholeRemoteDataIteratorCL( const ListT::iterator& it) : base(it) {}

    /// \brief Go to the next element (prefix)
    inline WholeRemoteDataIteratorCL& operator ++ ();
    /// \brief Go to the next element (suffix)
    inline WholeRemoteDataIteratorCL  operator ++ (int);

    /// \brief assignment operator
    inline WholeRemoteDataIteratorCL& operator= ( const WholeRemoteDataIteratorCL&);
    /// \brief get end iterator
    ListT::iterator GetEnd() const { return list_[3].end(); }
};

/// \brief base class for all distributed data managed by the DiST module, i.e., all simplex classes
class TransferableCL
{
  public:
    typedef RemoteDataCL::ProcList_const_iterator ProcList_const_iterator;

  protected:
    GeomIdCL gid_;           ///< Each entity stores a geometric id

    const RemoteDataCL& GetRemoteData() const;                               ///< Access its RemoteDataCL
          RemoteDataCL& GetRemoteData();                                     ///< Access its RemoteDataCL
    virtual void UpdateGID() = 0;

  public:
    /// \brief Uninitialized object
    TransferableCL() : gid_( NoGID) {}
    /// \brief Copy a transferable object
    TransferableCL( const TransferableCL& t) : gid_(t.gid_), Unknowns(t.Unknowns) {}
    TransferableCL( const GeomIdCL& h) : gid_(h) {}
    TransferableCL( Uint lvl, const Point3DCL& p, const Usint dim) : gid_( lvl, p, dim) {}
    virtual ~TransferableCL() {}

    /// \brief Put the simplex on an outgoing stream
    virtual void Pack (MPIostreamCL&) const = 0;
    /// \brief Init the simplex by an incoming stream
    virtual void UnPack (MPIistreamCL&) = 0;

    //// \name Information each simplex provides concerning distributed computing
    //@{
    const GeomIdCL& GetGID()      const { return gid_; }                             ///< Ask for the geometric id
    const Point3DCL&        GetBary()     const { return gid_.bary; }                        ///< Ask for the barycenter of the simplex
    Usint                   GetDim()      const { return gid_.dim; }                         ///< Ask for the dimension of the simplex
    Uint                    GetLevel()    const { return gid_.level; }                       ///< Ask for the level
    size_t                  GetHashId()   const { Hashing h; return h(GetGID()); }   ///< Ask for the hash of the GID
    Priority                GetPrio()     const { return GetRemoteData().GetLocalPrio(); }   ///< Ask for the priority
    bool                    IsMaster()    const { return GetPrio()>=PrioMaster; }            ///< Check if simplex is a master copy
    bool                    IsLocal()     const { return GetNumDist()==1; }                  ///< Check if the simplex is local
    Uint                    GetNumDist ( Priority prio=NoPrio) const                         ///< Get number of procs on which the simplex with priority >= prio is stored
        { return GetRemoteData().GetNumProcs(prio); }
    bool                    IsDistributed( Priority prio=NoPrio) const                       ///< Check a copy with priority >= prio exists on another process
        { return GetRemoteData().IsDistributed(prio); }
    int                     GetOwner()    const { return GetRemoteData().GetOwnerProc(); }   ///< Get the process rank of the owning process
    bool                    AmIOwner()    const { return GetRemoteData().AmIOwner(); }       ///< Check if I am the owner
    virtual bool            IsOnProcBnd() const { return GetRemoteData().IsOnProcBnd(); }    ///< Check if the simplex is located at a process boundary (overloaded by tetras)
    ProcList_const_iterator GetProcListBegin() const { return GetRemoteData().GetProcListBegin(); }
    ProcList_const_iterator GetProcListEnd()   const { return GetRemoteData().GetProcListEnd(); }
    //@}
    void Identify( const TransferableCL& parent, const PrioListCL& prios)                    ///< Make simplex a distributed object with same processes as the parent simplex (notes: no communication required, assigning PrioNeutral on all processes)
        { GetRemoteData().Identify( parent, prios); }

    UnknownHandleCL Unknowns;

    virtual void DebugInfo( std::ostream&) const;                                            ///< Write debug information to a stream
    virtual bool IsMarkedForRemovement() const = 0;                                          ///< check if marked for removal
    virtual void SetRemoveMark        () const = 0;                                          ///< set mark for removal
};

/// \brief Use operator << to put data on a stream
inline MPIostreamCL& operator<< (MPIostreamCL& os, const TransferableCL& t)
{
    t.Pack( os);
    return os;
}

/// \brief Use operator >> to get data out of stream
inline MPIistreamCL& operator>> (MPIistreamCL& is, TransferableCL& t)
{
    t.UnPack( is);
    return is;
}


} // end of namespace DiST
} // end of namespace DROPS

#include "DiST/remotedata.tpp"

#endif /* REMOTEDATA_H_ */
