/// \file remotedata.cpp
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

#include "DiST/remotedata.h"
#include "DiST/DiST.h"

namespace DROPS
{

namespace DiST
{

// E R R O R   C L A S S
//-----------------------

std::ostream& ErrorCL::what(std::ostream& out) const
{
    out << "["<<ProcCL::MyRank()<<"] ";
    out << _ErrMesg << std::endl;
    if (gid_!=NoGID) {
        if (!InfoCL::Instance().Exists(gid_))
            out << "simplex " << gid_ << " not known to DiST!" << std::endl;
        else {
            const RemoteDataCL& rd= InfoCL::Instance().GetRemoteData(gid_);
            rd.DebugInfo(out);
            rd.GetLocalObject().DebugInfo(out);
        }
    }
    return out;
}

void ErrorCL::handle() const
{
    ProcCL::RecoverStdOstreams();
    what(std::cout);
    std::cout.flush();
    std::abort();
}

// R E M O T E  D A T A  C L
//--------------------------

std::ostream& operator<< ( std::ostream& os, const RemoteDataCL::ProcListT& pl)
// for debugging
{
    RemoteDataCL::ProcList_const_iterator it( pl.begin());
    if (it!=pl.end()) {
        os << '(' << it->proc << ',' << PriorityToString( it->prio) << ')';
        ++it;
    }
    for (; it!=pl.end(); ++it) {
        os << ", (" << it->proc << ',' << PriorityToString( it->prio) << ')';
    }
    return os;
}

void RemoteDataCL::MakeConsistent()
{
    const int me= ProcCL::MyRank();
    // find local entry
    ProcListT::iterator it= GetProcListBegin(),
            end= GetProcListEnd();
    for (; it->proc!=me && it!=end; ++it) {}
    if (it==end) {
        std::cerr << "RemoteDataCL::MakeConsistent: found no entry for local object in proc list, my rank = " << me << std::endl;
        this->DebugInfo(std::cerr);
        throw ErrorCL( "RemoteDataCL::MakeConsistent: found no entry for local object in proc list", GetLocalObject().GetGID());
    }
    // local entry at first position
    std::swap( *GetProcListBegin(), *it);
}

void RemoteDataCL::UpdateOwner( const LoadVecT& load)
/// Owner is proc with smallest load. If more than one proc with same minimal load exist, the one with minimal rank is taken.
{
    double minLoad= std::numeric_limits<double>::max();
    owner_= std::numeric_limits<int>::max();
    const bool masterExists= PrioExists(PrioMaster);
    for (ProcListT::const_iterator it= GetProcListBegin(), end= GetProcListEnd(); it!=end; ++it) {
        if (masterExists && it->prio < PrioMaster)
            continue;
        const int proc= it->proc;
        if (load[proc] < minLoad) {
            minLoad= load[proc];
            owner_= proc;
        } else if ((load[proc] == minLoad) && (proc < owner_))
            owner_= proc;
    }
}

Priority RemoteDataCL::GetPrio(int rank) const
{
    for (ProcListT::const_iterator it= GetProcListBegin(), end= GetProcListEnd(); it!=end; ++it)
        if (it->proc==rank)
            return it->prio;

    return NoPrio;
}

void RemoteDataCL::Identify( const TransferableCL& parent, const PrioListCL& prios)
/// \pre Should only be called for local objects (IsLocal()==true).
/// \pre Should be called on all procs storing the local object.
/// \post Processor list is a copy of parent's one except for the priorities. These are set to PrioMaster on all processes.
/// \param prios List of priorities that should be considered for \a parent
{
    Assert( IsLocal(), DROPSErrCL("RemoteDataCL::Identify: object is already distributed"), DebugDiSTC);
    procList_.resize(0); // remove local entry
    for (ProcListT::const_iterator it= parent.GetProcListBegin(), end= parent.GetProcListEnd(); it!=end; ++it)
        if ( prios.contains(it->prio))
            procList_.push_back( ProcListEntryCL( it->proc, PrioMaster));
    // proc list is consistent due to consistency of parent's proc list.
    // update owner
    UpdateOwner( InfoCL::Instance().GetLoadVector());
}

void RemoteDataCL::DebugInfo( std::ostream& os) const
{
    os << "Simplex " << localObj_->GetGID() << ", owner = " << owner_ << '\n'
       << " o stored on " << procList_ << std::endl;
}

bool RemoteDataCL::IsSane( std::ostream& os) const
{
    bool sane= true;
    // Check if first entry belongs to me
    if ( GetProcListBegin()->proc!=ProcCL::MyRank()){
        os << "First entry in the process list is not me!" << std::endl;
        sane= false;
    }
    // Check that owner was set properly (does not check for consistency on other processes)
    if ( owner_<0 || owner_>=ProcCL::Size()) {
        os << "Found invalid owner with rank " << owner_ << std::endl;
        sane= false;
    }
    return sane;
}

MPIostreamCL& operator<< (MPIostreamCL& os, const RemoteDataCL::ProcListT& pl)
{
    if (os.isBinary())
        write_char_array( os, reinterpret_cast<const char*>( &pl[0]),
                          pl.size()*sizeof( RemoteDataCL::ProcListEntryCL));
    else {
        os << pl.size();
        for (RemoteDataCL::ProcListT::const_iterator it= pl.begin(), end= pl.end(); it != end; ++it)
            os << static_cast<int>( it->proc) << static_cast<int>( it->prio);
    }
    return os;
}

MPIistreamCL& operator>> (MPIistreamCL& is, RemoteDataCL::ProcListT& pl)
{
    pl.clear();
    if (is.isBinary()) {
        std::streamsize num_char;
        is >> num_char;
        pl.resize( num_char/sizeof( RemoteDataCL::ProcListEntryCL));
        is.read( reinterpret_cast<char*>( &pl[0]), num_char);
    }
    else {
        size_t num;
        is >> num;
        pl.reserve( num);
        int proc, prio;
        for (size_t i= 0; i < num; ++i) {
            is >> proc >> prio;
            pl.push_back( RemoteDataCL::ProcListEntryCL( proc, Priority( prio)));
        }
    }
    return is;
}

Uint RemoteDataCL::GetNumProcs( Priority prio) const
{
    if (prio == NoPrio)
        return procList_.size();
    Uint num= 0;
    for (ProcList_const_iterator it=GetProcListBegin(); it!=GetProcListEnd(); ++it){
        if ( it->prio>=prio)
        ++num;
    }
    return num;
}


// R E M O T E  D A T A  L I S T  C L
//-----------------------------------

/// handler for the DiST interface to check, whether proc/prio lists are the same on different processes
class RemoteDataListCL::DebugHandlerCL
{
  private:
    std::ostream& os_;

  public:
    DebugHandlerCL( std::ostream& os) : os_(os) {}

    bool Gather( const TransferableCL& t, SendStreamCL& s)
    { // write proc/prio list
        RemoteDataCL::ProcListT proclist( t.GetProcListBegin(), t.GetProcListEnd());
        s << proclist << t.GetOwner();
        return true;
    }

    bool Scatter( TransferableCL& t, const size_t numData, MPIistreamCL& r)
    { // read proc/prio lists and compare with own
        RemoteDataCL::ProcListT myplist( t.GetProcListBegin(), t.GetProcListEnd());
        for (size_t i=0; i<numData; ++i) {
            RemoteDataCL::ProcListT plist;
            int owner;
            r >> plist >> owner;
            if (owner != t.GetOwner()) {
                os_ << "Inconsistent owner, mine = " << t.GetOwner() << ", remote = " << owner << std::endl;
                t.DebugInfo( os_);
                return false;
            }
            if (myplist.size() != plist.size()) {
                os_ << "Differing size of remote data for simplex " << t.GetGID() << std::endl;
                t.DebugInfo( os_);
                return false;
            }
            for (RemoteDataCL::ProcListT::const_iterator pit= plist.begin(), pend= plist.end(); pit!=pend; ++pit)
                if (!is_in( myplist.begin(), myplist.end(), *pit)) {
                    os_ << "Missing remote data entry ( " << pit->proc << ", " << pit->prio << " ) for simplex " << t.GetGID() << std::endl;
                    t.DebugInfo( os_);
                    return false;
                }
        }
        return true;
    }

    bool Check( Usint dim)
    // do interface comm
    {
        InterfaceCL::DimListT dimlist;
        dimlist.push_back( dim);

        const PrioListT   allPrios;
        const LevelListCL allLvls;
        // communicate over all objects
        InterfaceCL comm( allLvls, allPrios, allPrios, dimlist);
        return comm.Communicate( *this);
    }
};

void RemoteDataListCL::DebugInfo( std::ostream& os) const
{
    os << "#Elements in RemoteDataList: " << this->size() << '\n';
    for ( const_iterator it(begin()); it!=end(); ++it){
        it->second.DebugInfo(os);
    }
}

bool RemoteDataListCL::IsSane( std::ostream& os) const
{
    bool sane= true;
    if ( !this->empty()){
        // Check if each entry in this list has the same dimension
        const Usint dim= this->begin()->first.dim;
        for ( const_iterator it(this->begin()); it!=this->end(); ++it){
            if ( it->first.dim!=dim){
                os << "Dimension in RemoteDataListCL is not consistent!" << std::endl;
                it->second.DebugInfo( os);
                sane= false;
            }
            sane =  sane && it->second.IsSane( os);
        }
        // Check that remote data is the same on different processes
        DebugHandlerCL debug( os);
        sane= sane && debug.Check(dim);
    }
    return sane;
}

WholeRemoteDataIteratorCL::WholeRemoteDataIteratorCL( const DimListT& dims, const LevelListCL& lvl, const PrioListT& prios, bool dist)
    : base( &InfoCL::InstancePtr()->GetRemoteList(0)), useGIDVec_(false), dims_(dims)
{
    levels_= lvl;
    prios_= prios;
    distributed_= dist;
    // Iterate through all elements of the unordered_map until an element is covered by the iterator
    for (dimIdx_= 0; dimIdx_< dims_.size(); ++dimIdx_) {
        const Usint dim= dims_[dimIdx_];
        for (pos_=list_[dim].begin(); pos_!=list_[dim].end(); ++pos_) {
            if ( contains())
                return;
        }
    }
    // for empty sequence, we have pos_ == GetEnd()
    pos_= GetEnd();
}

WholeRemoteDataIteratorCL::WholeRemoteDataIteratorCL( GIDIteratorT begin, GIDIteratorT end)
    : base( &InfoCL::InstancePtr()->GetRemoteList(0)), useGIDVec_(true), itGID_(begin), endGID_(end)
{
    // search for first element of GID vector
    FindGID(); // sets pos_
    // for empty sequence, we have pos_ == GetEnd()
}

// T R A N S F E R A B L E  C L
// ----------------------------

const RemoteDataCL& TransferableCL::GetRemoteData() const
{
    return InfoCL::Instance().GetRemoteData( *this);
}

RemoteDataCL& TransferableCL::GetRemoteData()
{
    return InfoCL::Instance().GetRemoteData( *this);
}

void TransferableCL::DebugInfo( std::ostream& os) const
{
    os << " parallel information:\n"
       << " o priority: " << PriorityToString( GetPrio()) << "\n";
    if ( IsLocal()){
        os << " o only stored locally\n";
    }
    else{
        os << " o " << (IsOnProcBnd() ? "is" : "is not") << " located at a process boundary\n"
           << " o owner: " << GetRemoteData().GetOwnerProc() << '\n'
           << " o processes: ";
        for ( DiST::RemoteDataCL::ProcList_const_iterator p(GetRemoteData().GetProcListBegin());
              p!= GetRemoteData().GetProcListEnd(); ++p){
                  os << p->proc << " (" << PriorityToString( p->prio) << ") ";
        }
        os << '\n';
    }
}

} // end of namespace DiST
} // end of namespace DROPS
