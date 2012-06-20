/// \file unknowns.h
/// \brief Implementation of the mapping from simplices to indices to
///    linear-algebra data-structures.
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

#ifndef DROPS_UNKNOWNS_H
#define DROPS_UNKNOWNS_H

#include "misc/utils.h"
#include <limits>
#include <vector>

namespace DROPS
{

#ifdef _PAR
// fwd declaration
namespace DiST{
    class TransferableCL;
    class MPIostreamCL;
    class MPIistreamCL;
}
#endif

/// The type of indices used to access the components of vectors,
/// matrices, etc.
typedef Ulint IdxT;

/// \brief Value of all unset indices.
///
/// \note If we ever have as many unknowns, there will of course be a
/// problem -- but probably there would be still one more index and the
/// wrap around would kill us anyways.
const IdxT NoIdx= std::numeric_limits<IdxT>::max();


/// Implementation-detail of UnknownHandleCL.
class UnknownIdxCL
{
  private:
    std::vector<IdxT> _Idx;
#ifdef _PAR
    /// If a vertex or an edge is transfered during the migration algorithm, a temporary
    /// copy of the DOFs is stored in a receive buffer. To decide where the DOF on the
    /// simplex is stored (receive buffer or in the regular vector), we make use of the
    /// flags in the vector received_.
    mutable std::vector<bool> received_;
#endif

  public:
    UnknownIdxCL( Uint numsys) : _Idx( numsys, NoIdx) {}
    UnknownIdxCL() {}
    UnknownIdxCL(const UnknownIdxCL&);
    ~UnknownIdxCL() {}
    UnknownIdxCL& operator=( const UnknownIdxCL&);

    /// get the position in the vector of the DOF of the system number sysnum (writing access)
    IdxT& GetIdx( Uint sysnum)
    {
        Assert( sysnum<GetNumSystems(), DROPSErrCL("UnknownIdxCL: Sysnum out of range"), DebugUnknownsC);
        return _Idx[sysnum];
    }

    /// get the position in the vector of the DOF of the system number sysnum (reading access)
    IdxT  GetIdx( Uint sysnum) const
    {
        Assert( sysnum<GetNumSystems(), DROPSErrCL("UnknownIdxCL: Sysnum out of range"), DebugUnknownsC);
        return _Idx[sysnum];
    }

    /// get the maximal number of system numbers
    Uint GetNumSystems()            const { return _Idx.size(); }

    void resize( Uint size, const IdxT defaultIdx= NoIdx)
    {
        _Idx.resize(size, defaultIdx);
    }

    /// DOF to a new system number
    void push_back(IdxT idx= NoIdx) { _Idx.push_back( idx); }

#ifdef _PAR
    /// mark a DOF of the system number sysnum as received
    void SetUnkRecv( Uint sysnum) const
    {
        if ( received_.size()<=sysnum)
            received_.resize( sysnum+1,false);
        received_[sysnum]=true;
    }

    /// check if a DOF of the system number sysnum has been received
    bool GetUnkRecv( Uint sysnum) const
    {
        if ( received_.size()<=sysnum)
            return false;
        else
            return received_[ sysnum];
    }

    /// mark all DOF as not received, i.e., forget about receive information
    void ResetUnkRecv() { received_.resize(0); }

    /// mark a DOF of the system number sysnum as not received
    void ResetUnkRecv( Uint sysnum)
    {
        if ( sysnum < received_.size())
            received_[ sysnum]= false;
    }

    /// check if any DOF has been received
    bool HasUnkRecv() const
    {
        for ( Uint sysnum=0; sysnum<received_.size(); ++sysnum)
            if ( received_[sysnum])
                return true;
        return false;
    }
#endif
};


/// \brief Maps a simplex and a "sysnum" (system number) on an index for
///     accessing numerical data.
///
/// Every simplex has a public member Unknowns of type UnknownHandleCL,
/// which behaves as a container of indices (for numerical data) that
/// can be accessed via a system number.  This class is only a handle
/// for an UnknownIdxCL.
class UnknownHandleCL
{
  private:
    UnknownIdxCL* _unk;              ///< stores the indices of DoF for different FE spaces
#ifdef _PAR
    std::vector<bool>  mayHaveUnk_;  ///< stores for each triangulation level, whether this DoF may belong to a FE space of this level (important in parallel case)
#endif

  public:
    UnknownHandleCL() : _unk(0) {}
    UnknownHandleCL( const UnknownHandleCL& orig)
#ifdef _PAR
      : mayHaveUnk_(orig.mayHaveUnk_)
#endif
    {
        _unk= orig._unk ? new UnknownIdxCL( *orig._unk)
                        : 0;
    }

    UnknownHandleCL& operator=( const UnknownHandleCL& rhs)
    {
        if (this==&rhs) return *this;
        delete _unk;
        _unk= rhs._unk ? new UnknownIdxCL( *rhs._unk)
                       : 0;
#ifdef _PAR
        mayHaveUnk_= rhs.mayHaveUnk_;
#endif
        return *this;
    }

    ~UnknownHandleCL() { delete _unk; }

    void Init(Uint numsys= 0)
    {
        Assert( _unk==0, DROPSErrCL("UnknownHandleCL: Init was called twice"), DebugUnknownsC);
        _unk= new UnknownIdxCL(numsys);
    }

    void Destroy() { delete _unk; _unk= 0; }

    /// True, iff this instance has already acquired an UnknownIdxCL-object.
    bool Exist()             const { return _unk; }
    /// True, iff the system sysnum exists and has a valid index-entry.
    bool Exist( Uint sysnum) const { return _unk && sysnum<_unk->GetNumSystems() && _unk->GetIdx(sysnum) != NoIdx; }

    /// Effectively deletes the index belonging to system sysnum.
    void Invalidate( Uint sysnum) { _unk->GetIdx(sysnum)= NoIdx; }

    UnknownIdxCL* Get() const { return _unk; }
    /// Retrieves the index for sysnum for writing.
    IdxT&        operator() ( Uint i)       { return _unk->GetIdx(i); }
    /// Retrieves the index for sysnum for reading.
    IdxT         operator() ( Uint i) const { return _unk->GetIdx(i); }

    /// Allocates memory for a system with number sysnum.  Afterwards, an index
    /// can be stored for sysnum.
    /// The initial index is set to NoIdx. Thus, .Exist( sysnum)==false.
    void Prepare( Uint sysnum)
    {
        if (!_unk) _unk= new UnknownIdxCL( sysnum+1);
        else if ( !(sysnum < _unk->GetNumSystems()) )
            _unk->resize( sysnum+1, NoIdx);
    }

    /// \brief Write out some debug information
    void DebugInfo( std::ostream& os) const;

#ifdef _PAR
    /// \brief Returns whether a FE space on triangulation level \a lvl may contain current DoF.
    ///
    /// In the parallel case, only DoF of a master tetra may have unknowns on a certain triangulation level the master tetra belongs to.
    /// Thus, e.g., DoFs on a vertex may contribute to FE spaces on level 0, 2, but not on level 1.
    bool InTriangLevel( int lvl) const
    {
        if (lvl < 0) lvl += mayHaveUnk_.size();
    	return mayHaveUnk_[lvl];
    }
    /// Commit that a FE space on triangulation level \a lvl may contain current DoF
    void EnableUnknowns( Uint lvl)
    {
    	mayHaveUnk_[lvl] = true;
    }
    /// Allocates memory for \a maxlevel levels. Initially, InTriangLevel() is false for all levels.
    void DisableAllUnknowns( int maxlevel)
    {
    	mayHaveUnk_.clear();
    	mayHaveUnk_.resize(maxlevel, false);
    }

    /// Remember that the DOF of the system number sysnum has been received
    void SetUnkReceived( Uint sysnum) const
    {
        Assert(_unk!=0, DROPSErrCL("UnknownHandleCL::SetUnkReceived: missing Prepare()"), DebugUnknownsC | DebugParallelC);
        _unk->SetUnkRecv( sysnum);
    }

    /// Check if the DOF of the system number sysnum has been received
    bool UnkReceived( Uint sysnum) const { return _unk && _unk->GetUnkRecv( sysnum); }

    /// Mark DOF of all system numbers as not received
    void ResetUnkReceived() const
    {
        if (_unk!=0)
            _unk->ResetUnkRecv();
    }

    /// Mark DOF of the system number sysnum as not received
    void ResetUnkReceived( Uint sysnum) const
    {
        if (_unk!=0)
            _unk->ResetUnkRecv( sysnum);
    }

    /// For Debugging Purpose: Check if there is an UnkRecv-Flag
    bool HasUnkReceived() const
    {
        if (_unk==0)
            return false;
        else
            return _unk->HasUnkRecv();
    }

    /// Put all DOFs handled by this handler onto the stream
    void Pack( DiST::MPIostreamCL&, const DiST::TransferableCL&) const;

    /// Get all information about DOF from a stream
    void UnPack( DiST::MPIistreamCL&, const DiST::TransferableCL&);
#else
    bool InTriangLevel( int) const {return true;}
#endif
};

} // end of namespace DROPS

#endif
