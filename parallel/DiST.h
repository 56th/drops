/// \file DiST.h
/// \brief The DiST (short for distributed simplex type) module is responsible for the distributed geometric data structure on a parallel machine.
/// It should replace the DDD library in the near future.
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
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#ifndef DROPS_DIST_H
#define DROPS_DIST_H

#include "misc/container.h"
#include "misc/utils.h"
#include "num/unknowns.h"
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



/// \enum Priority Priorities for distributed objects
/// \todo Put these priorities in the DiST namespace
enum Priority
{
    NoPrio          = -1,    ///< invalid, used for cases where no priority has been assigned
    PrioNeutral     = 0,     ///< no priority
    PrioKilledGhost = 1,     ///< prio of ghost tetras marked for removal
    PrioVGhost      = 2,     ///< for subs of overlapping master tetra (don't constitute a proc bnd)
    PrioGhost       = 3,     ///< for ghost tetras and subs only owned by ghost tetras (skipped by the public iterators)
    PrioMaster      = 4,     ///< master copy of tetras and their subs
    PrioHasUnk      = 5      ///< simplices having potentially unknowns
};

inline std::string PriorityToString( const Priority);

/// \brief Namespace for all functions constituting a parallel multigrid
/** Everything in this namespace acts as an interface to DROPS. */
namespace DiST{

/// \brief Priority list. By definition, an empty list contains all priorities.
class PrioListCL: public std::vector<Priority>
{
  public:
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

/// \name Ask for the dimension of a simplex
//@{
template <typename SimplexT> inline Usint GetDim();
template <> inline Usint GetDim<VertexCL>() { return 0; }
template <> inline Usint GetDim<EdgeCL>()   { return 1; }
template <> inline Usint GetDim<FaceCL>()   { return 2; }
template <> inline Usint GetDim<TetraCL>()  { return 3; }
//@}

/// \brief Helper functions for constituting the parallel multigrid
/** Everything in this namespace should not act as an interface to DROPS. */
namespace Helper{

// fwd declaration
class RemoteDataListCL;
class RemoteDataListIteratorCL;

/// \brief Assign each simplex an unique geometric id
struct GeomIdCL
{
    Uint        level;      ///< level, the simplex occurs first
    Point3DCL   bary;       ///< barycenter of the simplex
    Usint       dim;        ///< Dimension of the simplex, i.e., 0 - vertex, 1 - edge, 2 - face, 3 - tetrahedron, 4 - uninitialized

    GeomIdCL() : level((Uint)(-1)), bary(), dim(4) {}
    GeomIdCL(Uint lvl, const Point3DCL& p, Usint dimension) : level(lvl), bary(p), dim(dimension) {}
    GeomIdCL( const GeomIdCL& h) : level(h.level), bary(h.bary), dim(h.dim) {}
    template <typename SimplexT>
    GeomIdCL(Uint lvl, const SimplexT& s) : level(lvl), bary( GetBaryCenter(s)), dim(GetDim<SimplexT>()) {}
    bool operator== (const GeomIdCL& h) const { return h.level == level && h.bary == bary;}
    bool operator!= (const GeomIdCL& h) const { return !(h==*this); }
};

inline std::ostream& operator << ( std::ostream& os, const GeomIdCL& h)
{
    static char scode[]= "VEFT"; // simplex code for each dimension
    os << scode[h.dim] << h.level << " (" << h.bary << ')';
    return os;
}
const GeomIdCL NoGID= GeomIdCL( (Uint)(-1), Point3DCL(), 4);  // Dummy id, if simplex does not exist

#if __GNUC__ >= 4 || DROPS_WIN
struct Hashing : std::unary_function<GeomIdCL, size_t>
{
    size_t operator()(const GeomIdCL& h) const
    {
        size_t seed = 0;
        DROPS_STD_HASH<int> inthasher;
        DROPS_STD_HASH<double> doublehasher;
        // see boost::hash_combine
        seed ^= inthasher(h.level) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= doublehasher(h.bary[0]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= doublehasher(h.bary[1]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= doublehasher(h.bary[2]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};
#endif

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

/// \brief Helper union for writing and reading numbers in streams
/** \todo Maybe we need a memcopy and this is not save?*/
template <typename T>
union ToBinary
{
    T    value;
    char binary[sizeof(T)];
};

/// \brief Outgoing stream.
class SendStreamCL : public std::ostringstream
/** This stream is employed as a mean of transport for sending data (integers,
    vertices, elements, etc.) out of the process towards other processes. <p>
    It is derived from the ostringstream class of the Standard IOstream Library
    and it is characterized by saving all kind of data as a string. <p>
    It's counterpart is the class RecvStreamCL.
    \todo Is there a way to get rid of the copy sendbuf_?
*/
{
public:
    typedef std::ostringstream base;

  private:
    bool binary_;           ///< flag for binary sending/receiving
    std::string sendbuf_;   ///< buffering the sending data

  public:
    SendStreamCL( const bool binary=true) : std::ostringstream(), binary_(binary) {}
    inline bool isBinary() const {return binary_;}
    /// \brief Non-blocking send to process 'dest'.
    inline ProcCL::RequestT Isend(int dest, int tag=5);
    /// \brief Clear buffer.
    /** We do not really clear the buffer, we just move the inserting pointer (put pointer) to the beginning. */
    void clearbuffer() { seekp(0); sendbuf_.clear(); }
};

/// \brief Incoming stream.
class RecvStreamCL : public std::istringstream
/** This stream is employed as a mean of transport for receiving incoming data
    (integers, vertices, elements, etc.) sent from other process. <p>
    It is derived from the istringstream class of the Standard IOstream Library
    and it is characterized by saving all kind of data as a string. <p>
    It's counterpart is the class SendStreamCL.
*/
{
public:
    typedef std::istringstream base;

  private:
    bool binary_;   ///< flag for binary sending/receiving

  public:
    /// @param[in] binary is true if the stream store the data in binary; in ASCII otherwise.
    RecvStreamCL( const bool binary=true)
        : std::istringstream( binary ? std::ios_base::binary : std::ios_base::in), binary_(binary) {}
    RecvStreamCL( const SendStreamCL& s)
        : std::istringstream( s.str()), binary_(s.isBinary()) {}
    /// \brief Gives back how many bytes of data contains our stream at the time.
    inline int getbufsize() const { return ((int) const_cast<RecvStreamCL*>(this)->tellg())+1; }
    inline bool isBinary() const {return binary_;}
    /// \brief Blocking receive from process 'source'.
    void Recv(int source, int tag=5);
    /// \brief Go back to the beginning of the stream.
    inline void resetbuffer() { this->seekg(0); }
};

/// \brief operator << for SendStreamCL
template<typename T>
SendStreamCL& operator<<( SendStreamCL& os, const T& t)
{
    if (os.isBinary()) {
        ToBinary<T> bin;
        bin.value= t;
        os.write( bin.binary, sizeof(T));
    } else {
        SendStreamCL::base& oss= dynamic_cast<SendStreamCL::base&>(os);
        oss << t << ' ';
    }
    return os;
}

/// \brief operator >> for RecvStreamCL
template<typename T>
RecvStreamCL& operator>>( RecvStreamCL& is, T& t)
{
    if (is.isBinary()) {
        ToBinary<T> bin;
        is.read( bin.binary, sizeof(T));
        t= bin.value;
    } else {
        RecvStreamCL::base& iss= dynamic_cast<RecvStreamCL::base&>(is);
        iss >> t;
    }
    return is;
}

/// \brief Use operator << to put a GeomIdCL object on a stream
SendStreamCL& operator<< ( SendStreamCL&, const GeomIdCL&);
/// \brief Use operator << to put a Point3DCL object on a stream
SendStreamCL& operator<< ( SendStreamCL&, const Point3DCL&);

/// \brief Use operator >> to get a GeomIdCL object out of stream
RecvStreamCL& operator>> ( RecvStreamCL&, GeomIdCL&);
/// \brief Use operator >> to get a Point3DCL object out of stream
RecvStreamCL& operator>> ( RecvStreamCL&, Point3DCL&);


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

          Uint            GetNumProcs()       const { return procList_.size(); }///< Get number of processes owning the distributed entity
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
           Uint     GetNumPrio(Priority prio) const;                            ///< Get the number of processes with priority >= \a prio
           bool     IsLocal()      const { return GetNumProcs()==1; }           ///< Check if the simplex is only stored on this process
    inline bool     IsDistributed( Priority prio=NoPrio) const;                 ///< Check if a simplex is distributed and at least one remote object has a priority>=prio
    inline bool     IsOnProcBnd()  const;                                       ///< Check if the simplex is located on a process boundary

    /// \brief Debugging
    //@{
    bool IsSane( std::ostream&) const;
    void DebugInfo( std::ostream&) const;
    //@}
};

/// \brief Use operator << to put a RemoteDataCL object on a stream
SendStreamCL& operator<< ( SendStreamCL&, const GeomIdCL&);

/// \brief Use operator >> to get a RemoteDataCL object out of stream
RecvStreamCL& operator>> ( RecvStreamCL&, GeomIdCL&);

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
    inline RemoteDataListIteratorCL( RemoteDataListCL*, const LevelListCL&, const PrioListT&, bool);

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
    inline WholeRemoteDataIteratorCL( const DimListT&, const LevelListCL&, const PrioListT&, bool dist);
    /// \brief Standard constructor defining sequence by a vector of GID's
    inline WholeRemoteDataIteratorCL( GIDIteratorT begin, GIDIteratorT end);
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

class SimplexTransferInfoCL
{
  public:
    typedef std::map<int,Priority> ProcSetT;

  private:
    Helper::RemoteDataCL& rd_;       ///< remote data of simplex on local proc
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
    void AddProc( int proc, Priority prio, bool changeLocalPrio=false);
    /// \brief Adds new proc/prio entries without changing existing entries.
    void AddProcSet( const ProcSetT& procs);
    /// return whether simplex will be on proc \a p after transfer
    bool WillBeOnProc( int p) const { return postProcs_.find(p) != postProcs_.end(); }
    /// return set of procs/prios after transfer
    const ProcSetT& GetPostProcs() const { return postProcs_; }
    /// return set of procs/prios where simplex has to be sent to
    const ProcSetT& GetSendToProcs() const { return procsToSend_; }
    /// compute set of procs/prios where simplex has to be sent to (call after post proc list is complete)
    void ComputeSendToProcs();
    /// mark simplex for deletion
    void SetRemoveMark()       { RemoveMark_= true; }
    /// return whether simplex will be removed
    bool WillBeRemoved() const { return RemoveMark_; }
    /// return whether simplex will update its sub simplices (only makes sense for tetrahedron)
    bool UpdateSubs()    const { return UpdateSubs_; }
    /// return rank of broadcaster
    int  GetBroadcaster() const { return proc_bc_; }
    /// return simplex' remote data
    const Helper::RemoteDataCL& GetRemoteData() const { return rd_; }
    /// return simplex' remote data
    Helper::RemoteDataCL&       GetRemoteData()       { return rd_; }
};

}   // end of namespace Helper


/// \brief base class for all distributed data managed by the DiST module, i.e., all simplex classes
class TransferableCL
{
  public:
    typedef Helper::RemoteDataCL::ProcList_const_iterator ProcList_const_iterator;

  protected:
    Helper::GeomIdCL gid_;           ///< Each entity stores a geometric id

    inline const Helper::RemoteDataCL& GetRemoteData() const;                               ///< Access its RemoteDataCL
    inline       Helper::RemoteDataCL& GetRemoteData();                                     ///< Access its RemoteDataCL

  public:
    /// \brief Uninitialized object
    TransferableCL() : gid_( Helper::NoGID) {}
    /// \brief Copy a transferable object
    TransferableCL( const TransferableCL& t) : gid_(t.gid_), Unknowns(t.Unknowns) {}
    TransferableCL( const Helper::GeomIdCL& h) : gid_(h) {}
    TransferableCL( Uint lvl, const Point3DCL& p, const Usint dim) : gid_( lvl, p, dim) {}
    virtual ~TransferableCL() {}

    /// \brief Put the simplex on an outgoing stream
    virtual void Pack( Helper::SendStreamCL&) const = 0;
    /// \brief Init the simplex by an incoming stream
    virtual void UnPack( Helper::RecvStreamCL&) = 0;

    //// \name Information each simplex provides concerning distributed computing
    //@{
    const Helper::GeomIdCL& GetGID()      const { return gid_; }                             ///< Ask for the geometric id
    const Point3DCL&        GetBary()     const { return gid_.bary; }                        ///< Ask for the barycenter of the simplex
    Usint                   GetDim()      const { return gid_.dim; }                         ///< Ask for the dimension of the simplex
    Uint                    GetLevel()    const { return gid_.level; }                       ///< Ask for the level
    size_t                  GetHashId()   const { Helper::Hashing h; return h(GetGID()); }   ///< Ask for the hash of the GID
    Priority                GetPrio()     const { return GetRemoteData().GetLocalPrio(); }   ///< Ask for the priority
    bool                    IsMaster()    const { return GetPrio()>=PrioMaster; }            ///< Check if simplex is a master copy
    bool                    MayStoreUnk() const { return GetPrio()==PrioMaster; }            ///< Check for ability of storing unknowns due to priority \todo Is PrioMaster correct?
    bool                    IsLocal()     const { return GetNumDist()==1; }                  ///< Check if the simplex is local
    Uint                    GetNumDist () const { return GetRemoteData().GetNumProcs(); }    ///< Get number of procs on which the simplex is stored
    bool                    IsDistributed( Priority prio=NoPrio) const                       ///< Check a copy with priority >= prio exists on another process
        { return GetRemoteData().IsDistributed(prio); }
    int                     GetOwner()    const { return GetRemoteData().GetOwnerProc(); }   ///< Get the process rank of the owning process
    bool                    AmIOwner()    const { return GetRemoteData().AmIOwner(); }       ///< Check if I am the owner
    virtual bool            IsOnProcBnd() const { return GetRemoteData().IsOnProcBnd(); }    ///< Check if the simplex is located at a process boundary (overloaded by tetras)
    virtual bool            IsExclusive( const Priority) const = 0;                          ///< Check if simplex is exclusive
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
inline Helper::SendStreamCL& operator<< ( Helper::SendStreamCL& os, const TransferableCL& t)
{
    t.Pack( os);
    return os;
}

/// \brief Use operator >> to get data out of stream
inline Helper::RecvStreamCL& operator>> ( Helper::RecvStreamCL& is, TransferableCL& t)
{
    t.UnPack( is);
    return is;
}

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
    typedef Helper::WholeRemoteDataIteratorCL iterator;         ///< type for iterator over elements in the interface
    typedef iterator::DimListT DimListT;                        ///< type for storing the involved simplices
    typedef iterator::GIDIteratorT GIDIteratorT;                ///< type for iterator over GID's
    typedef std::map<int, Helper::SendStreamCL*> SendListT;     ///< type for storing data to be sent (proc -> data)
    typedef std::map<int, Helper::RecvStreamCL*> RecvListT;     ///< type for receiving data (proc -> data)

  private:
    enum CommPhase {                ///< Communication phases
        bothPhases,                 ///< owner and copies gather and scatter data, i.e., copies -> owner -> copies
        toowner,                    ///< copies send to owners, i.e., copies -> owner
        fromowner                   ///< owner send to copies, i.e., owner -> copies
    };
    RecvListT      recvbuf_;        ///< received data (after call of function Communicate())
    SendListT      sendbuf_;        ///< data to be sent
    PrioListT      from_;           ///< list of priority on sender side, the interface operates on
    PrioListT      to_;             ///< list of priority on receiver side, the interface operates on
    bool           binary_;         ///< transfer data binary or in ASCII
    std::set<int>  ownerRecvFrom_;  ///< list of processes, the owner receives data from
    std::set<int>  ownerSendTo_;    ///< list of processes, the owner sends data to
    std::set<int>  IRecvFromOwners_;///< list of owner processes, I have to receive data from
    iterator       begin_from_,     ///< iterators defining sequence of interface elements
                   begin_to_,
                   begin_,
                   end_;

    /// \brief Call the gather handler for each entity covered by the iterators [begin, end).
    template <typename HandlerT>
    void GatherData( HandlerT&, const iterator& begin, const iterator& end, CommPhase phase);
    /// \brief Call the scatter handler for each entity whose GID is given in the istream_.
    template <typename HandlerT>
    bool ScatterData( HandlerT&);
    /// \brief Handles the communication between the processes
    void Communicate( CommPhase phase);
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
                 const DimListT& dims, const bool dist=true, const bool binary=true)
        : from_(from), to_(to), binary_(binary),
          begin_from_( dims, lvl, from, dist),
          begin_to_( dims, lvl, to, dist), begin_( dims, lvl, /*all prios*/PrioListT(), dist), end_( begin_from_.GetEnd()) {}

    InterfaceCL( GIDIteratorT begin, GIDIteratorT end, const bool binary=true)
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
    ///    Gather( DiST::TransferableCL&, DiST::Helper::SendStreamCL&)
    ///    Scatter( DiST::TransferableCL&, const size_t&, DiST::Helper::RecvStreamCL&);
    //@{
    /// \brief Do the interface communication on the interface specified by the
    ///   the constructor
    template <typename HandlerT>
    bool PerformInterfaceComm( HandlerT&);

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
    typedef std::map<const TransferableCL*, Helper::SimplexTransferInfoCL> UpdateListT;
    typedef std::pair<const TransferableCL*,Helper::SimplexTransferInfoCL> UpdateEntryT;
    typedef UpdateListT::iterator                                          UpdateIterator;
    typedef Helper::SimplexTransferInfoCL::ProcSetT                        ProcSetT;

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
    /// \brief Update the remote data to account for changed prios
    void UpdateRemoteData();
    /// \brief Update ownership (in remote data)
    void UpdateOwners();
    /// \brief Delete all simplices that are not needed any more.
    void DeleteUnusedSimplices( bool del);

  public:
    /// \brief Constructor with a given multigrid (\a mg) and decision if communication should be done \a binary.
    /// For \a del = true all unused simplices are removed from the multigrid, otherwise the RemoveMark is set and the simplex is unregistered from the DiST module.
    ModifyCL( MultiGridCL& mg, bool del= true, bool binary=true)
        : modifiable_(false), del_(del), binary_( binary), mg_(mg), entsToUpdt_(0) {}

    /// \brief Call Init() before any modifications
    void Init();
    /// \brief Call Finalize() after any modifications (initiates the communication)
    void Finalize();
    /// \brief Change the priority of a locally stored object
    void ChangePrio( const TransferableCL&, Priority);
    /// \brief Consider a local simplex for deletion.
    void Delete( const TransferableCL&);
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
    typedef std::map<int, Helper::SendStreamCL*> SendBufT;
    /// \brief type for sorted tetras to update
    typedef std::list<UpdateEntryT>              SortedListT;
    // comparison functor to sort by descending level
    struct CompByDescLevelCL;

  private:
    SendBufT      sendBuffer_;        ///< all information to be sent

    /// create list of tetras to update, sorted by descending order
    SortedListT* SortUpdateTetras();
    /// Update and send remote data of a given simplex.
    void UpdateSendRemoteData( const TransferableCL&, Helper::SimplexTransferInfoCL&);
    /// \brief allocate and fill send buffers, update remote data
    void FillSendBuffer();
    /// \brief Receive and create simplices and remote data
    void Receive();
    /// \brief Receive vertices, edges, faces, tetras (helper for Receive()),
    template <typename SimplexT>
    void ReceiveSimplices( Helper::RecvStreamCL&, size_t);
    /// \brief Create simplex using the SimplexFactoryCL
    template <typename SimplexT>
    SimplexT& CreateSimplex( const SimplexT&, const Helper::RemoteDataCL::ProcListT&);

  public:
    /// \brief Constructor with a given multigrid (\a mg) and decision if the transfer should be done \a binary.
    /// For \a del = true all unused simplices are removed from the multigrid, otherwise the RemoveMark is set and the simplex is unregistered from the DiST module.
    TransferCL( MultiGridCL& mg, bool del= true, bool binary=true)
        : base( mg, del, binary) {}
    /// \brief To be called before marking tetrahedra for transfer
    void Init();
    /// \brief To be called after marking tetrahedra for transfer (initiates the communication)
    void Finalize();
    /// \brief Mark a tetrahedron for transfer
    void MarkForTransfer( const TetraCL&, int toProc, Priority prio, bool del= true);
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
	Helper::RemoteDataListCL remoteData_[4];
    Helper::RemoteDataCL::LoadVecT loadOfProc_;

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
    const Helper::RemoteDataCL& GetRemoteData( const TransferableCL& t) const { return GetRemoteData( t.GetGID()); }
    /// \brief Access the RemoteDataCL of a given simplex (private?)
    Helper::RemoteDataCL& GetRemoteData( const TransferableCL& t) { return GetRemoteData( t.GetGID()); }
    inline const Helper::RemoteDataCL& GetRemoteData( const Helper::GeomIdCL&) const; ///< Get remote data by GID
    inline       Helper::RemoteDataCL& GetRemoteData( const Helper::GeomIdCL&);       ///< Get remote data by GID

    /// \brief Access the RemoteDataListCL of a given dimension
    const Helper::RemoteDataListCL& GetRemoteList( int dim) const { return remoteData_[dim]; }
    /// \brief Access the RemoteDataListCL of a given dimension
    Helper::RemoteDataListCL& GetRemoteList( int dim) { return remoteData_[dim]; }
    /// \brief Access the RemoteDataListCL of a given simplex type
    template <typename SimplexT>
    const Helper::RemoteDataListCL& GetRemoteList() const { return remoteData_[GetDim<SimplexT>()]; }
    /// \brief Access the RemoteDataListCL of a given simplex type (private?)
    template <typename SimplexT>
    Helper::RemoteDataListCL& GetRemoteList() { return remoteData_[GetDim<SimplexT>()]; }

    /// \brief Access loads of processes used to determine the ownership in the remote data.
    const Helper::RemoteDataCL::LoadVecT& GetLoadVector() const { return loadOfProc_; }

    /// \brief Is this GID registered?
    inline bool Exists( const Helper::GeomIdCL& gid) const { return gid.dim < 4 ? remoteData_[gid.dim].Exists(gid) : false; }

    /// \name Direct access to a simplex by its GID
    /// \todo Make these functions inline
    //@{
    VertexCL* GetVertex( const Helper::GeomIdCL&) const;
    EdgeCL*   GetEdge  ( const Helper::GeomIdCL&) const;
    FaceCL*   GetFace  ( const Helper::GeomIdCL&) const;
    TetraCL*  GetTetra ( const Helper::GeomIdCL&) const;
    //@}

    /// \name Debugging
    //@{
    bool IsSane( std::ostream&) const;
    void DebugInfo( std::ostream&) const;
    void SizeInfo( std::ostream&) const;
    void ShowSimplex( const size_t hash_val, std::ostream&) const;     ///< Show the debug information of a simplex by a given hash value
    //@}
};

}   // end of namespace DiST
}   // end of namespace DROPS

#include "parallel/DiST.tpp"
#endif
