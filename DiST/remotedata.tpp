/// \file remotedata.tpp
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

std::string PriorityToString( const Priority prio)
{
    switch (prio){
    case NoPrio         : return "invalid";
    case PrioNeutral    : return "neutral";
    case PrioKilledGhost: return "killed ghost";
    case PrioVGhost     : return "vertical ghost";
    case PrioGhost      : return "ghost";
    case PrioMaster     : return "master";
    }
    std::stringstream ret;
    ret << "undefined" << static_cast<int>(prio);
    return ret.str();
}

namespace DiST{
class InfoCL;
// L E V E L  L I S T  C L A S S
//------------------------------

LevelListCL::LevelListCL( Uint max_level)
{
    reserve( max_level+1);
    for (Uint lvl=0; lvl<=max_level; ++lvl)
        push_back( lvl);
}


template<class T>
ErrorCL::ErrorCL(const std::string& mesg, const T& data, const GeomIdCL& gid)
  : DROPSErrCL(mesg), gid_(gid)
{
    std::ostringstream oss;
    oss << data;
    _ErrMesg+= oss.str();
}


// R E M O T E  D A T A  C L A S S
//--------------------------------

/** Returns the priority of the first element in the procList_.
    \pre the first element in the list of the RemoteDataCL must be the process itself. */
Priority RemoteDataCL::GetLocalPrio() const
{
    Assert( !procList_.empty() && procList_.front().proc==ProcCL::MyRank(), DROPSErrCL("RemoteDataCL: Inconsistency"), DebugDiSTC);
    return procList_.front().prio;
}

bool RemoteDataCL::PrioExists(Priority prio) const
{
    for (ProcListT::const_iterator it= GetProcListBegin(), end= GetProcListEnd(); it!=end; ++it)
        if (it->prio >= prio)
            return true;
    return false;
}

/** Checks, if the simplex is not stored only local and at least two non VGhost copies exist. */
bool RemoteDataCL::IsOnProcBnd() const
{
    if ( IsLocal())
        return false;
    size_t counter=0;
    for ( ProcList_const_iterator it(GetProcListBegin()); it!=GetProcListEnd(); ++it){
        if ( it->prio > PrioVGhost)
            if (++counter >= 2)
                return true;
    }
    return false;
}

bool RemoteDataCL::IsDistributed( Priority prio) const
{
    Uint num= 0;
    for (ProcList_const_iterator it=GetProcListBegin(); it!=GetProcListEnd(); ++it){
        if ( it->prio>=prio)
            if (++num >= 2)
                return true;
    }
    return false;
}

// R E M O T E  D A T A  L I S T  C L A S S
// ----------------------------------------

/** The following specification of entities can be formulated:
    \param[in] Levels Iterate through elements of level list 'Levels'.
    \param[in] a      Iterate through elements, whose priority matches any of the priorities listed in 'a'.
    \param[in] dist   Iterate only through distributed elements.
*/
RemoteDataListIteratorCL RemoteDataListCL::beginLvlPrioIt ( const LevelListCL& Levels, const PrioListT& a, bool dist)
{
    return RemoteDataListIteratorCL(this, Levels, a, dist);
}

/// \brief end-Iterator
RemoteDataListIteratorCL RemoteDataListCL::endLvlPrioIt ()
{
    return RemoteDataListIteratorCL(this->end());
}


/** Note that this function does not update some inter-process information and does not
    invoke any communication.*/
void RemoteDataListCL::Register( TransferableCL& t, Priority prio)
{
    Assert( this->find(t.GetGID())==this->end(), ErrorCL("RemoteDataListCL::Register: Simplex already known: ", t.GetGID()), DebugDiSTC);
    Assert( empty() || begin()->first.dim==t.GetDim(), ErrorCL("RemoteDataListCL::Register: wrong dimension of ", t.GetGID(), NoGID), DebugDiSTC);
    (*this)[t.GetGID()]= RemoteDataCL( &t, prio);
}

/** Note that this function does not update some inter-process information and does not
    invoke any communication.*/
void RemoteDataListCL::Register( TransferableCL& t, const RemoteDataCL::ProcListT& pl)
{
    Assert( this->find(t.GetGID())==this->end(), ErrorCL("RemoteDataListCL::Register: Simplex already known", t.GetGID()), DebugDiSTC);
    Assert( empty() || begin()->first.dim==t.GetDim(), ErrorCL("RemoteDataListCL::Register: wrong dimension of ", t.GetGID(), NoGID), DebugDiSTC);
    (*this)[t.GetGID()]= RemoteDataCL( &t, pl);
}

/** Note that this function does not update some inter-process information and does not
    invoke any communication.*/
void RemoteDataListCL::Unregister( TransferableCL& t)
{
    Assert( this->find(t.GetGID())!=this->end(), ErrorCL("RemoteDataListCL::Unregister: Simplex unknown", t.GetGID()), DebugDiSTC);
    this->erase(t.GetGID());
}

// R E M O T E  D A T A  L I S T  I T E R A T O R  C L A S S
//----------------------------------------------------------

bool RemoteDataListIteratorCL::contains()
{
    ///Check that the number or processors is greater 2 if distributed, else continue with the next element of the unordered_map
    if (!distributed_ || !pos_->second.IsLocal())
    {
        ///Check that the level of the current element is contained in 'levels_', else continue with the next element of the unordered_map
        if (levels_.contains( pos_->first.level))
        {
            ///Check if the priority of the Simplex is any of the priorities in the Priority List
            return prios_.contains( pos_->second.GetLocalPrio());
        }
    }
    return false;
}

RemoteDataListIteratorCL::RemoteDataListIteratorCL( RemoteDataListCL* list, const LevelListCL& lvl, const PrioListT& prios, bool dist)
    : list_(list), levels_(lvl), prios_(prios), distributed_(dist)
{
    // Iterate through all elements of the unordered_map until an element is covered by the iterator
    for ( pos_=list_->begin(); pos_!=list_->end(); pos_++){
        if ( contains())
            return;
    }
    // for empty sequence, we have pos == list_->end
}


/** Go to the next element (prefix)
*/
RemoteDataListIteratorCL& RemoteDataListIteratorCL::operator++ ()
{
    Assert( pos_!=list_->end(), DROPSErrCL("RemoteDataListIteratorCL::operator ++: iterator already points to the end"), DebugDiSTC);
    for (++pos_; pos_!=list_->end(); ++pos_)
    {
        if ( contains())
            return *this;
    }
    return *this;
}

RemoteDataListIteratorCL RemoteDataListIteratorCL::operator++ ( int)
{
    RemoteDataListIteratorCL tmp(*this);
    ++(*this);
    return tmp;
}

inline RemoteDataListIteratorCL& RemoteDataListIteratorCL::operator= ( const RemoteDataListIteratorCL& it)
{
    list_       = it.list_;
    pos_        = it.pos_;
    levels_     = it.levels_;
    prios_      = it.prios_;
    distributed_= it.distributed_;
    return *this;
}


// W H O L E  R E M O T E  D A T A  L I S T  I T E R A T O R  C L A S S
//---------------------------------------------------------------------

Usint WholeRemoteDataIteratorCL::FindGID()
{
    const Usint dim= itGID_!=endGID_ ? itGID_->dim : 3;
    pos_= itGID_!=endGID_ ? list_[dim].find( *itGID_) : GetEnd();
    return dim;
}

WholeRemoteDataIteratorCL& WholeRemoteDataIteratorCL::operator++ ()
{
    Assert( pos_!=GetEnd(), DROPSErrCL("RemoteDataListIteratorCL::operator ++: iterator already points to the end"), DebugDiSTC);
    if (useGIDVec_) {
        ++itGID_;
        FindGID();
    } else {
        ++pos_;
        for (;;) { // ever
            for (; pos_!=list_[dims_[dimIdx_]].end(); ++pos_)
            {
                if ( contains())
                    return *this;
            }
            ++dimIdx_;
            if (dimIdx_==dims_.size()) {
                pos_= GetEnd();
                return *this;
            }
            // otherwise pos_ should point to beginning of next dimension
            pos_= list_[dims_[dimIdx_]].begin();
        }
    }
    return *this;
}

WholeRemoteDataIteratorCL WholeRemoteDataIteratorCL::operator++ ( int)
{
    WholeRemoteDataIteratorCL tmp(*this);
    ++(*this);
    return tmp;
}

inline WholeRemoteDataIteratorCL& WholeRemoteDataIteratorCL::operator= ( const WholeRemoteDataIteratorCL& orig)
{
    dynamic_cast<base&>(*this)= orig;
    useGIDVec_ = orig.useGIDVec_;
    itGID_     = orig.itGID_;
    endGID_    = orig.endGID_;
    dims_      = orig.dims_;
    dimIdx_    = orig.dimIdx_;
    return *this;
}

} // end of namespace DiST
} // end of namespace DiST
