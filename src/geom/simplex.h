/// \file simplex.h
/// \brief classes that constitute the simplex
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Eva Loch, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_SIMPLEX_H
#define DROPS_SIMPLEX_H

#include <list>
#include <vector>
#include <algorithm>
#include <iostream>
#include "misc/utils.h"
#include "geom/boundary.h"
#include "geom/topo.h"
#include "num/unknowns.h"
#include "DiST/geomid.h"

#ifdef _PAR
#  include "DiST/DiST.h"
#  include "parallel/parallel.h"
#  include <map>
#endif

namespace DROPS
{

// fwd decl
class BoundaryCL;
class PeriodicEdgesCL;
class SimplexFactoryCL;

// Classes for simplices in the triangulation
class VertexCL;
class EdgeCL;
class FaceCL;
class TetraCL;

// Containers for storage of the simplices
typedef GlobalListCL<VertexCL> MG_VertexContT;
typedef GlobalListCL<EdgeCL>   MG_EdgeContT;
typedef GlobalListCL<FaceCL>   MG_FaceContT;
typedef GlobalListCL<TetraCL>  MG_TetraContT;


//**************************************************************************
// Classes: FaceWrapperCL, RecycleBinCL                                    *
// Purpose: When the triangulation is modified, some simplices will have   *
//          to be removed in one step of the algorithm, but will be needed *
//          again some steps later. In order not to lose the information   *
//          in the 'Unknowns', we store those simplices in a 'RecycleBin'. *
//          To identify a simplex in the 'RecycleBin' we use its vertices. *
//          'Edges' and 'Tetras' know their vertices, but 'faces' do not.  *
//          Therefore, we provide a 'FaceWrapper' that adds its 'Vertices' *
//          to a 'Face'.                                                   *
//**************************************************************************

struct FaceWrapperCL
{
  public:
    const FaceCL*   face;
    const VertexCL* vert1;
    const VertexCL* vert2;

    FaceWrapperCL(const FaceCL* f, const VertexCL* v1, const VertexCL* v2)
        : face(f), vert1(v1), vert2(v2) {}
};


class RecycleBinCL
{
  public:
    typedef std::list<const EdgeCL*>  EdgeContT;
    typedef std::list<FaceWrapperCL>  FaceWrapperContT;
    typedef std::list<const TetraCL*> TetraContT;

  private:
    EdgeContT        Edges_;
    FaceWrapperContT Faces_;
    TetraContT       Tetras_;

  public:
    inline const EdgeCL*  FindEdge  (const VertexCL* v) const;
    inline const FaceCL*  FindFace  (const VertexCL* v1, const VertexCL* v2) const;
    inline const TetraCL* FindTetra (const VertexCL* v1, const VertexCL* v2, const VertexCL* v3) const;

    void Recycle (const EdgeCL* ep)  { Edges_.push_back(ep); } // ??? nicht alles doppelt speichern ???
    void Recycle (const FaceCL* fp, const VertexCL* v1, const VertexCL* v2)
        { Faces_.push_back( FaceWrapperCL(fp,v1,v2) ); }
    void Recycle (const TetraCL* tp) { Tetras_.push_back(tp); }

    void DebugInfo(std::ostream&) const;
};

/*******************************************************************
*   V E R T E X  C L                                               *
*******************************************************************/
/// \brief Represents vertices in the multigrid
/** Contains the geometric part ('Coord_', 'BndVerts_') of a
    point in the multigrid as well as some topological ('Level_')
    information.
    It also stores some algorithmic information like 'RemoveMark_'
    and the RecycleBins.                                          */
/*******************************************************************
*   V E R T E X  C L                                               *
*******************************************************************/
class VertexCL
#ifdef _PAR
    : public DiST::TransferableCL
#endif
{
  public: // friend declarations
    friend class MultiGridCL;
    friend class SimplexFactoryCL;
#ifdef _PAR
    friend class ParMultiGridCL;
    typedef DiST::TransferableCL base;
#else
    UnknownHandleCL Unknowns;       ///< access to the unknowns on the vertex
#endif

  public: // typedefs and members
    typedef std::vector<BndPointCL>::iterator       BndVertIt;
    typedef std::vector<BndPointCL>::const_iterator const_BndVertIt;



  private: // members
#ifndef _PAR
    Point3DCL                Coord_;          ///< global coordinates of the vertex. In parallel mode, the level is stored by the base class
    Uint                     Level_ : 8;      ///< level, where the vertex occurs first. In parallel mode, the level is stored by the base class
#endif
    IdCL<VertexCL>           Id_;             ///< id of the vertex on this proc
    std::vector<BndPointCL>* BndVerts_;       ///< parametrization for each boundary segment the vertex is part of
    RecycleBinCL*            Bin_;            ///< recycle-bin
    mutable bool             RemoveMark_;     ///< flag, if this vertex should be removed

  private: // functions
    /// \name RecycleBin
    //@{
    bool  HasRecycleBin() const { return Bin_; }
    const RecycleBinCL* GetRecycleBin() const { return Bin_; }
    RecycleBinCL* GetCreateRecycleBin() { return HasRecycleBin() ? Bin_ : (Bin_= new RecycleBinCL); }
    //@}

    /// \brief Regular constructor called by the SimplexFactoryCL
    inline VertexCL (const Point3DCL& Point3D, Uint FirstLevel, IdCL<VertexCL> id= IdCL<VertexCL>());

  public:
    /// \brief Copy-constructor.
    inline VertexCL (const VertexCL&);
    /// \brief Constructs an uninitialized vertex.
    inline VertexCL();
    /// \brief Delete a vertex, also deletes the recycle-bin and boundary information
    inline ~VertexCL();

#ifndef _PAR
    /// \brief get level of vertex (=first appearance in the multigrid)
    /** In the parallel version this function is defined in the base class*/
    Uint GetLevel() const { return Level_; }
    DiST::GeomIdCL GetGID() const { return DiST::GeomIdCL(this->GetLevel(), *this);}
#endif

    /// \brief get coordinate of this vertex
    inline const Point3DCL& GetCoord() const;

    /// \brief change the coordinate of the vertex, e.g. ALE method in poisson problem
    void ChangeCoord     (Point3DCL& p);

    /// \brief Check if the vertex can be found in a triangulation level
    bool IsInTriang( Uint TriLevel) const { return  GetLevel() <= TriLevel; }

    /// \name Boundary
    //@{
    inline void AddBnd ( const BndPointCL&);        ///< add boundary-information
    inline void BndSort();                          ///< sort boundary-segments
    bool IsOnBoundary() const { return BndVerts_; } ///< check if this vertex lies on domain boundary
    const_BndVertIt GetBndVertBegin() const { return BndVerts_->begin(); }
    const_BndVertIt GetBndVertEnd  () const { return BndVerts_->end(); }
    /// \brief check if Vertex belongs on a special boundary index given by a boundary vertex
    inline bool HasBnd( const BndPointCL& BndVert) const;
    //@}

    /// \name RemovementMarks
    //@{
    bool IsMarkedForRemovement () const { return RemoveMark_; }     ///< check if vertex is marked for removement
    void SetRemoveMark         () const { RemoveMark_ = true; }     ///< set mark for removement
    void ClearRemoveMark       ()       { RemoveMark_ = false; }    ///< clear mark for removement
    //@}

    /// \name RecycleBin
    //@{
    /// \brief empty recycle-bin
    void DestroyRecycleBin () { delete Bin_; Bin_= 0; }
    /// \brief put a pointer to an edge into the recycle-bin of this vertex
    void Recycle (const EdgeCL* ep) { GetCreateRecycleBin()->Recycle(ep); }
    /// \brief put a pointer to a face into the recycle-bin of this vertex
    void Recycle (const FaceCL* fp, const VertexCL* vp1, const VertexCL* vp2)
        { GetCreateRecycleBin()->Recycle(fp,vp1,vp2); }
    /// \brief put a pointer to a tetra into the recycle-bin of this vertex
    void Recycle (const TetraCL* tp)
        { GetCreateRecycleBin()->Recycle(tp); }
    /// \brief Find an edge in the recycle-bin by the opposite vertex. Returns 0 if no edge is found.
    EdgeCL*  FindEdge  (const VertexCL* v) const
        { return HasRecycleBin() ? const_cast<EdgeCL*>(GetRecycleBin()->FindEdge(v)) : 0; }
    /// \brief Find a face in the recycle-bin by the other vertices. Returns 0 if no face is found.
    FaceCL*  FindFace  (const VertexCL* v1, const VertexCL* v2) const
        { return HasRecycleBin() ? const_cast<FaceCL*>(GetRecycleBin()->FindFace(v1,v2)) : 0; }
    /// \brief Find a tetra in the recycle-bin by the other vertices. Returns 0 if no tetra is found.
    TetraCL* FindTetra (const VertexCL* v1, const VertexCL* v2, const VertexCL* v3) const
        { return HasRecycleBin() ? const_cast<TetraCL*>(GetRecycleBin()->FindTetra(v1,v2,v3)) : 0; }
    //@}

#ifdef _PAR
    /// \name Functions for parallel computing demanded by the base class
    //@{
    /// \brief Put vertex on a stream
    void Pack( DiST::MPIostreamCL&) const;
    /// \brief Generate vertex from a stream
    void UnPack( DiST::MPIistreamCL&);
    /// \brief Generate GID
    void UpdateGID();
    //@}
#endif

    /// \name Debugging
    //@{
    bool IsSane    (std::ostream&, const BoundaryCL& ) const;       ///< check for sanity of this vertex
    void DebugInfo (std::ostream&) const;                           ///< get debug-information
    const IdCL<VertexCL>& GetId() const { return Id_; }             ///< Get id of this vertex (numbered locally on this proc)
   //@}
};


/*******************************************************************
*   E D G E  C L                                                   *
*******************************************************************/
/// \brief Represents an edge in the multigrid
/** The refinement algorithm works by manipulating the marks 'MFR_'
    of the edges. This is possible, because the refinement pattern
    of a face/tetrahedron is determined by the pattern of its edges.
    Marking the edges ensures consistency between the neighbors.  */
/*******************************************************************
*   E D G E  C L                                                   *
*******************************************************************/
class EdgeCL
#ifdef _PAR
    : public DiST::TransferableCL
#endif
{
  public:
    friend class PeriodicEdgesCL;
    friend class SimplexFactoryCL;
#ifdef _PAR
    friend class ParMultiGridCL;
    typedef DiST::TransferableCL base;
#else
    UnknownHandleCL Unknowns;                       ///< access to unknowns on this edge
#endif

  public:
    typedef MG_VertexContT::LevelCont VertContT;    ///< container for vertices
    typedef MG_EdgeContT::LevelCont   EdgeContT;    ///< container for subedges



  private:
#ifndef _PAR
    Uint                   Level_ : 8;    ///< level of the edge (according to owning tetras)
#else
    mutable short int      AccMFR_;       ///< accumulated MFR over all procs
#endif
    SArrayCL<VertexCL*, 2> Vertices_;     ///< "left" and "right" vertex of the edge
    VertexCL*              MidVertex_;    ///< midvertex, if the edge is refined
    SArrayCL<BndIdxT, 2>   Bnd_;          ///< an edge can be found on (max) two boundary-segments
    mutable short int      MFR_;          ///< mark, if the edge should be/is refined (set by refinement-algo)
    short int              localMFR_;     ///< MFR!=localMFR iff edge is on periodic boundary
    mutable bool           RemoveMark_;   ///< mark for removement

    /// \brief constructor called by the SimplexFactoryCL
    inline EdgeCL (VertexCL* vp0, VertexCL* vp1, Uint Level, BndIdxT bnd0= NoBndC, BndIdxT bnd1= NoBndC, short int MFR=0);

  public:
    /// \brief Copy-constructor
    inline EdgeCL (const EdgeCL&);
    /// \brief Constructs an uninitialized edge
    inline EdgeCL();
    // default dtor

    /// \brief Add boundary-information
    void AddBndIdx(BndIdxT idx) { Bnd_[(Bnd_[0]==NoBndC ? 0 : 1)]= idx; }

    /// \name Midvertex
    //@{
    VertexCL*   GetMidVertex   ()             { return MidVertex_; }    ///< get pointer to midvertex
    void        SetMidVertex   (VertexCL* vp) { MidVertex_= vp; }       ///< set pointer to midvertex
    void        RemoveMidVertex()             { MidVertex_= 0; }        ///< remove pointer to midvertex without deleting it
    void        BuildMidVertex (SimplexFactoryCL&, const BoundaryCL&);  ///< create midvertex
    /// \brief Build both subedges
    inline void BuildSubEdges  (SimplexFactoryCL&, const BoundaryCL&, SArrayCL<EdgeCL*,2>&);
    //}

    /// \name Marks
    //@{
#ifndef _PAR
    bool IsMarkedForRef       () const { return MFR_; }                         ///< check if this edge is marked for refinement
    void IncMarkForRef        ()       { ++MFR_; ++localMFR_;}                  ///< increase mark for refinement count
    void DecMarkForRef        ()       { --MFR_; --localMFR_;}                  ///< decrease mark for refinement count
    void ResetMarkForRef      ()       { MFR_= 0; localMFR_= 0; }               ///< remove mark for refinement
#else
    bool IsMarkedForRef       () const { return AccMFR_; }                      ///< check if this edge is marked for refinement
    void IncMarkForRef        () const { ++MFR_; ++AccMFR_; }                   ///< increase the local and accumulated mark for refinement count
    void DecMarkForRef        () const { --MFR_; --AccMFR_; }                   ///< decrease the local and accumulated mark for refinement count
    void ResetMarkForRef      ()       { MFR_= 0; AccMFR_= 0; }                 ///< remove the local and accumulated mark for refinement
    short int GetAccMFR       () const {return AccMFR_;}                        ///< get accumulated mark for refinement
    void SetAccMFR ( short int MFR)    { AccMFR_= MFR; }                        ///< Set accumulated MFR
#endif
    bool IsMarkedForRemovement() const { return RemoveMark_; }                  ///< check if edge is marked for removement
    void SetRemoveMark        () const { RemoveMark_= true; }                   ///< set mark for removement
    void ClearRemoveMark      ()       { RemoveMark_= false; }                  ///< clear mark for removement
    short int GetMFR          () const { return MFR_; }                         ///< get mark for refinement of this proc
    //@}

    void RecycleMe   () const { Vertices_[0]->Recycle(this); }                  ///< put a pointer to this edge into the recycle-bin of the "left" vertex
    void SortVertices()                                                         ///< sort vertices by id
        { if (Vertices_[1]->GetId() < Vertices_[0]->GetId()) std::swap(Vertices_[0],Vertices_[1]); }

#ifndef _PAR
    /// \brief get level of vertex (=first appearance in the multigrid)
    /** In the parallel version this function is defined in the base class*/
    Uint GetLevel() const { return Level_; }
    DiST::GeomIdCL GetGID() const { return DiST::GeomIdCL(this->GetLevel(), *this);}
#else
    /// \name Functions for parallel computing demanded by the base class
    //@{
    /// \brief Put edge on a stream
    void Pack( DiST::MPIostreamCL&) const;
    /// \brief Generate edge from a stream
    void UnPack( DiST::MPIistreamCL&);
    /// \brief Generate GID
    void UpdateGID();
    //@}
#endif

    const VertexCL* GetVertex     (Uint i)             const { return Vertices_[i]; }       ///< get pointer to the "left" or the "right" vertex
    const VertexCL* GetMidVertex  ()                   const { return MidVertex_; }         ///< get midvertex
    /// \brief Get opposite vertex
    const VertexCL* GetNeighbor   (const VertexCL* vp) const { return vp==Vertices_[0] ? Vertices_[1] : Vertices_[0]; }
    /// \brief Check if the edge has a vertex
    bool            HasVertex     (const VertexCL* vp) const { return vp==Vertices_[0] || vp==Vertices_[1]; }
    bool            IsRefined     ()                   const { return MidVertex_; }         ///< check if edge is refined
    bool            IsOnBoundary  ()                   const { return Bnd_[0] != NoBndC; }  ///< check if edge lies on the domain boundary
    const BndIdxT*  GetBndIdxBegin()                   const { return Bnd_.begin(); }
    const BndIdxT*  GetBndIdxEnd  ()                   const
        { return IsOnBoundary() ? (Bnd_[1] == NoBndC ? Bnd_.begin()+1 : Bnd_.end() ) : Bnd_.begin(); }
    bool            IsInTriang    (Uint TriLevel)      const                                ///< check if edge can be found in a triangulation level
        { return GetLevel() == TriLevel || ( GetLevel() < TriLevel && !IsRefined() ); }

    // Debugging
    bool IsSane    (std::ostream&) const;                                                   ///< check for sanity
    void DebugInfo (std::ostream&) const;                                                   ///< get debug-information
};


/*******************************************************************
*   F A C E  C L                                                   *
*******************************************************************/
/// \brief Represents a face in the multigrid
/** It can have neighbors on two levels: two regular ones (one per
    side) on level 'Level_' and two green ones on the next.       */
/*******************************************************************
*   F A C E  C L                                                   *
*******************************************************************/
class FaceCL
#ifdef _PAR
    : public DiST::TransferableCL
#endif
{
  public:
    friend class SimplexFactoryCL;
#ifdef _PAR
    friend class ParMultiGridCL;
    typedef DiST::TransferableCL base;
#else
    UnknownHandleCL Unknowns;       ///< access to unknowns on a face (not yet used)
#endif

  private:
#ifndef _PAR
    Uint                       Level_ : 8;      ///< level of the face (=level according to tetras) (in parallel stored in _dddH.attr)
#endif
    SArrayCL<const TetraCL*,4> Neighbors_;      ///< neighbor tetras of the face (two on the same level, two on finer level)
    BndIdxT                    Bnd_;            ///< boundary-index of this face
    mutable bool               RemoveMark_;     ///< mark for removement

    /// \brief Create a face (serial)
    inline FaceCL (Uint Level, BndIdxT bnd= NoBndC);
    /// \brief Create a face (parallel)
    inline FaceCL (Uint Level, const Point3DCL& v0, const Point3DCL& v1, const Point3DCL& v2, BndIdxT bnd= NoBndC);

  public:
    /// \brief Copy a face
    inline FaceCL (const FaceCL&);
    /// \brief Create an uninitialized face
    inline FaceCL();
    // default dtor

    /// \name RemovementMarks
    //@{
    bool IsMarkedForRemovement() const { return RemoveMark_; }  ///< check if marked for removement
    void SetRemoveMark        () const { RemoveMark_= true; }   ///< set mark for removement
    void ClearRemoveMark      ()       { RemoveMark_= false; }  ///< clear mark for removement
    //@}

    /// \name Neighbors
    //@{
    void LinkTetra  (const TetraCL*);                           ///< link a tetra to face
    void UnlinkTetra(const TetraCL*);                           ///< unlink a tetra of face
    void SetNeighbor(Uint i, TetraCL* tp) { Neighbors_[i]=tp; } ///< Set tetra as neighbor
    inline bool IsLinkedTo (const TetraCL*) const;              ///< check if a tetra is linked to this face
    //@}

    /// \brief  Recycling
    /** put a pointer to this face into the recycle-bin of the first vertex */
    void RecycleMe(VertexCL* vp0, const VertexCL* vp1, const VertexCL* vp2) const
        { vp0->Recycle(this,vp1,vp2); }

#ifndef _PAR
    Uint GetLevel() const { return Level_; }                    ///< get level of the face (stored within the class)
    DiST::GeomIdCL GetGID() const { return DiST::GeomIdCL(this->GetLevel(), *this);}
#else
    /// \name Functions for parallel computing demanded by the base class
    //@{
    /// \brief Put face on a stream
    void Pack( DiST::MPIostreamCL&) const;
    /// \brief Generate face from a stream
    void UnPack( DiST::MPIistreamCL&);
    /// \brief Generate GID
    void UpdateGID();

    /// \brief Check if face is ghost
    bool IsGhost() const { return GetPrio()<PrioMaster; }
    /// \brief Get neighbor proc over the face
    /** Since a (non-VGhost) face is maximal stored by two processes, get the other process*/
    int GetNeighborProc() const
        { return GetPrio()!=PrioVGhost ? (GetRemoteData().GetProcListBegin()+1)->proc : -1; }
    //@}
#endif

    bool        IsOnNextLevel() const { return Neighbors_[2] || Neighbors_[3]; }        ///< check if face can be found in the next level
    bool        IsOnBoundary () const { return Bnd_ != NoBndC; }                        ///< check if face lies on the domain-boundary
    BndIdxT     GetBndIdx    () const { return Bnd_; }                                  ///< get index of the boundary-segment
    inline bool IsRefined    () const;                                                  ///< check if face is refined
    bool        IsInTriang   (Uint TriLevel) const                                      ///< check if face can be found in a triangulation level
      { return GetLevel() == TriLevel || ( GetLevel() < TriLevel && !IsRefined() ); }

    /// \brief Get simplex
    //@{
    inline const VertexCL* GetVertex(Uint)       const;                                 ///< get i'th vertex of the face
    inline const EdgeCL*   GetEdge  (Uint)       const;                                 ///< get i'th edge of the face
    inline const TetraCL*  GetTetra (Uint, Uint) const;                                 ///< get tetra of level and number
    //@}

    /// \name Neighboring tetras
    //@{
           const TetraCL* GetSomeTetra     () const { return Neighbors_[0]; }           ///< return pointer to first neighbor
    inline bool           HasNeighborTetra (const TetraCL*)       const;                ///< check if a tetra is neighbor
    inline const TetraCL* GetNeighborTetra (const TetraCL*)       const;                ///< get neighbor tetra of another tetra
           const TetraCL* GetNeighInTriang (const TetraCL*, Uint) const;                ///< get neighbor tetra in triangulation level
    inline Uint           GetFaceNumInTetra(const TetraCL*)       const;                ///< get number of face within a tetra
           const TetraCL* GetNeighbor      (Uint i) const { return Neighbors_[i];}      ///< get raw tetra-pointer from the array
    //@}

    /// \name Debugging
    //@{
    bool IsSane   (std::ostream&) const;                                                ///< check for sanity
    void DebugInfo(std::ostream&) const;                                                ///< get debug-information
    //@}
};


/*******************************************************************
*   T E T R A  C L                                                 *
*******************************************************************/
/// \brief Represents a tetrahedron in multigrid
/** This is probably the most important data structure of DROPS.
    All major routines that work on the grid (i.e. the refinement
    algorithm, the routines to set up a discretized system, the
    error estimator) "do it tetra by tetra".
    \todo (merge) RefRule_ and RefMark_ as single Usint OK?       */
/*******************************************************************
*   T E T R A  C L                                                 *
*******************************************************************/
class TetraCL
#ifdef _PAR
    : public DiST::TransferableCL
#endif
{
  public:
    friend class MultiGridCL;
    friend class SimplexFactoryCL;
#ifdef _PAR
    friend class ParMultiGridCL;
    typedef DiST::TransferableCL base;
#else
    UnknownHandleCL Unknowns;   ///< access to unknowns on tetras
#endif

    typedef SArrayCL<VertexCL*,NumVertsC>::iterator         VertexPIterator;            ///< iterator of pointers to vertices of this tetra
    typedef SArrayCL<VertexCL*,NumVertsC>::const_iterator   const_VertexPIterator;      ///< const version
    typedef SArrayCL<EdgeCL*,NumEdgesC>::iterator           EdgePIterator;              ///< iterator of pointers to edges of this tetra
    typedef SArrayCL<EdgeCL*,NumEdgesC>::const_iterator     const_EdgePIterator;        ///< const version
    typedef SArrayCL<FaceCL*,NumFacesC>::iterator           FacePIterator;              ///< iterator of pointers to faces of this tetra
    typedef SArrayCL<FaceCL*,NumFacesC>::const_iterator     const_FacePIterator;        ///< const version
    typedef SArrayCL<TetraCL*,MaxChildrenC>::iterator       ChildPIterator;             ///< iterator of pointers to children of this tetra
    typedef SArrayCL<TetraCL*,MaxChildrenC>::const_iterator const_ChildPIterator;       ///< const version
    typedef MG_VertexContT::LevelCont                       VertContT;                  ///< container for verts for linking purpose
    typedef MG_EdgeContT::LevelCont                         EdgeContT;                  ///< container for edges for linking purpose
    typedef MG_FaceContT::LevelCont                         FaceContT;                  ///< container for faces for linking purpose
    typedef MG_TetraContT::LevelCont                        TetraContT;                 ///< container for children for linking purpose

  private:
    // static arrays for linking edges and faces during the refinement algorithm
    static SArrayCL<EdgeCL*, NumAllEdgesC> ePtrs_;                      ///< EdgePointers for linking edges within refinement
    static SArrayCL<FaceCL*, NumAllFacesC> fPtrs_;                      ///< FacePointers for linking faces within refinement

#ifndef _PAR
    Uint                             Level_ : 8;
#endif
    IdCL<TetraCL>                    Id_;                               ///< id-number (locally numbered on one proc)
    Usint                            RefRule_;                          ///< actual refinement of the tetrahedron
    mutable Usint                    RefMark_;                          ///< refinement-mark (e.g. set by the error estimator)

    // subsimplices, parent, children
    SArrayCL<VertexCL*,NumVertsC>    Vertices_;                         ///< container for verts of tetra
    SArrayCL<EdgeCL*,NumEdgesC>      Edges_;                            ///< container for edges of tetra
    SArrayCL<FaceCL*,NumFacesC>      Faces_;                            ///< container for faces of tetra
    TetraCL*                         Parent_;                           ///< container for parent of tetra
    SArrayCL<TetraCL*,MaxChildrenC>* Children_;                         ///< container for children of tetra, for leaves: null-pointer
// ??? TODO: Kann man das ohne Speicherluecken und ohne Fragmentieren hinkriegen ???

    /// \brief constructor of verts and parent; FileBuilderCL has to construct the Id_, too, thus it can optionally be set.
    inline TetraCL (VertexCL*, VertexCL*, VertexCL*, VertexCL*, TetraCL*, IdCL<TetraCL> id= IdCL<TetraCL>());
#ifdef _PAR
    /// \brief constructor of verts and level, if no parent is available; FileBuilderCL has to construct the Id_, too, thus it can optionally be set.
    inline TetraCL (VertexCL*, VertexCL*, VertexCL*, VertexCL*, TetraCL*, Uint, IdCL<TetraCL> id= IdCL<TetraCL>());
#endif

  public:
    /// \brief Copy a tetrahedron
    inline TetraCL (const TetraCL&);
    /// \brief Constructor of an uninitialized tetrahedron
    inline TetraCL();
    /// \brief Delete a tetrahedron
    ~TetraCL() { if (Children_) delete Children_; }

    /// \name Interface for the refinement algorithm
    //@{
    /// \name Access to children, vertices
    //@{
    ChildPIterator GetChildBegin  ()        ///< "Pointer-Iterator" to first child
        { return Children_ ? Children_->begin() : 0; }
    ChildPIterator GetChildEnd    ()        ///< "Pointer-Iterator" to end of children
        { return Children_ ? Children_->begin() + GetRefData().ChildNum : 0; }
    VertexCL*      GetVertMidVert (Uint i)  ///< return pointer to midvertex of edge or vertex of tetra
        { return IsMidVert(i) ? Edges_[EdgeOfMidVert(i)]->GetMidVertex() : Vertices_[i]; }
    //@}

    /// \name Rules and marks
    //@{
    void        SetRefRule          (Uint RefRule) { RefRule_ = RefRule; }      ///< set refinement rule for this tetra
    void        RestrictMark        ();                                         ///< set regular refinement if marked and manipulate MFR on edges
    inline void CommitRegRefMark    () const;                                   ///< increase MFR on edges
    inline void UnCommitRegRefMark  () const;                                   ///< decrease MFR on edges
    inline void Close               ();                                         ///< calculate green closure
    void        ClearAllRemoveMarks ();                                         ///< save all subsimplices of this tetra and all child tetras from removement
    //@}

    /// \name Recycling
    //@{
    void RecycleMe        () const { Vertices_[0]->Recycle(this); }             ///< put pointer to the tetra into recycle-bin of vertex(0)
    void RecycleReusables ();                                                   ///< put all subsimplices of the tetra and children that are needed in next refinement step into recycle-bins
    //@}

    /// \name Building children
    //@{
    void        CollectEdges           (const RefRuleCL&, SimplexFactoryCL&, const BoundaryCL&);        ///< build or unrecycle edges that are needed for refinement
    void        CollectFaces           (const RefRuleCL&, SimplexFactoryCL&);                           ///< build or unrecycle faces that are needed for refinement
    inline void LinkEdges              (const ChildDataCL&);                                            ///< link edges from "ePtrs_" to the tetra according to the ChildDataCL
    inline void LinkFaces              (const ChildDataCL&);                                            ///< link faces from "fPtrs_" to the tetra according to the ChildDataCL
    void CollectAndLinkChildren (const RefRuleCL&, SimplexFactoryCL&);                                  ///< build, unrecycle and link children
    void        UnlinkFromFaces        ()                                                               ///< remove link from faces to the tetra
        { for (Uint face=0; face<NumFacesC; ++face) Faces_[face]->UnlinkTetra(this); }
    void DeleteChildren() { delete Children_; Children_=0; }
    //@}

    /// \name Functions for a builder
    //@{
    void BuildEdges        (SimplexFactoryCL&);                                  ///< build edges
    void BuildAndLinkFaces (SimplexFactoryCL&);                                  ///< build and link faces
    void SetFace           (Uint f, FaceCL* fp) { Faces_[f]= fp; }               ///< set face
    void SetEdge           (Uint e, EdgeCL* ep) { Edges_[e]= ep; }               ///< set edge
    void SetRefMark        (Uint refmark)       { RefMark_= refmark; }           ///< set RefMark
    void SetChild          (Uint, TetraCL*);                                     ///< set a child-pointer; if neccessary the Children_-Array is allocated first
    //@}
    //@}


#ifndef _PAR
    Uint GetLevel() const { return Level_; }                                     ///< return level of tetra
    bool IsGhost() const {return false;}
    DiST::GeomIdCL GetGID() const { return DiST::GeomIdCL(this->GetLevel(), *this);}
#else
    /// \name Functions for parallel computing demanded by the base class
    //@{
    /// \brief Put tetra on a stream
    void Pack( DiST::MPIostreamCL&) const;
    /// \brief Generate tetra from a stream
    void UnPack( DiST::MPIistreamCL&);
    /// \brief Generate GID
    void UpdateGID();
    /// \brief Merge tetra with given tetra. Needed by DiST::TransferCL.
    void Merge( const TetraCL&);

    bool IsProcBnd (Uint face) const { return Faces_[face]->IsOnProcBnd(); }     ///< check if face of tetra is on processor-boundary
    bool IsGhost    () const { return GetPrio()<PrioMaster; }                    ///< check if tetra is ghost
    bool HasGhost   () const;                                                    ///< check if tetra has a ghost-copy somewhere
    //@}
#endif

    const IdCL<TetraCL>& GetId () const { return Id_; }                          ///< get local id
    Uint GetRefMark            () const { return RefMark_; }                     ///< get refinement mark
    Uint GetRefRule            () const { return RefRule_; }                     ///< get refinement rule
    inline const RefRuleCL& GetRefData () const                                  ///< get information about refinement data
        { return DROPS::GetRefRule(this->GetRefRule() & 63); }

// RefMark_ is mutable
    void SetRegRefMark () const { RefMark_= RegRefMarkC; }                       ///< mark tetra for regular refinement
    void SetRemoveMark () const { RefMark_= RemoveMarkC; }                       ///< mark tetra for removement
    void SetNoRefMark  () const { RefMark_= NoRefMarkC; }                        ///< mark tetra for no refinement

    bool IsMarkEqRule   () const { return RefMark_ == RefRule_; }                ///< check if tetra is refined as the mark says
    bool IsUnrefined    () const { return RefRule_ == UnRefRuleC; }              ///< check if tetra is unrefined
    bool IsRegularlyRef () const { return RefRule_ == RegRefRuleC; }             ///< check if tetra is regular refined
    bool IsRegular      () const                                                 ///< check if the tetra is regular
#ifndef _PAR
      { return GetLevel()!=0 ? Parent_->GetRefRule() == RegRefRuleC : true; }
#else
      { return GetLevel()!=0 ? (IsGhost() ? true : Parent_->GetRefRule() == RegRefRuleC ) : true; }
#endif

    bool IsMarkedForRef        () const                                          ///< check if tetra is marked for refinement
      { return RefMark_ != NoRefMarkC && RefMark_ != RemoveMarkC; }
    bool IsMarkedForRegRef     () const { return RefMark_ == RegRefMarkC; }      ///< check if tetra is marked for regular refinement
    bool IsMarkedForRemovement () const { return RefMark_ == RemoveMarkC; }      ///< check if tetra is marked for removement
    bool IsMarkedForNoRef      () const { return RefMark_ == NoRefMarkC; }       ///< check if tetra is marked for no refinement

    /// \name access to subsimplices
    //@{
    const_VertexPIterator GetVertBegin ()   const { return Vertices_.begin(); }
    const_VertexPIterator GetVertEnd   ()   const { return Vertices_.end(); }
    const_EdgePIterator   GetEdgesBegin()   const { return Edges_.begin(); }
    const_EdgePIterator   GetEdgesEnd  ()   const { return Edges_.end(); }
    const_FacePIterator   GetFacesBegin()   const { return Faces_.begin(); }
    const_FacePIterator   GetFacesEnd  ()   const { return Faces_.end(); }
    const VertexCL*       GetVertex(Uint i) const { return Vertices_[i]; }
    const EdgeCL*         GetEdge  (Uint i) const { return Edges_[i]; }
    const FaceCL*         GetFace  (Uint i) const { return Faces_[i]; }
    //@}

    bool           IsBndSeg        (Uint face) const { return Faces_[face]->IsOnBoundary(); }           ///< check if face lies on domain-boundary
    bool           IsNeighbor      (Uint face) const { return Faces_[face]->HasNeighborTetra(this); }   ///< check if this tetra is neigbor to a face
    BndIdxT        GetBndIdx       (Uint face) const { return Faces_[face]->GetBndIdx(); }              ///< get boundary index of a face
    const TetraCL* GetNeighbor     (Uint face) const { return Faces_[face]->GetNeighborTetra(this); }   ///< get pointer to neighbor tetra over face
    const TetraCL* GetNeighInTriang(Uint face, Uint trilevel) const
      { return Faces_[face]->GetNeighInTriang( this, trilevel); }

    /// \name access to parent and children
    //@{
    const TetraCL*       GetParent     ()       const { return Parent_; }
    const_ChildPIterator GetChildBegin ()       const { return Children_ ? Children_->begin() : 0; }
    const_ChildPIterator GetChildEnd   ()       const { return Children_ ? Children_->begin() + GetRefData().ChildNum : 0; }
    const TetraCL*       GetChild      (Uint i) const { return (*Children_)[i]; }
    //@}

    double               GetVolume     () const;                                           ///< get volume of tetra
    double               GetNormal     (Uint face, Point3DCL& normal, double& dir) const;  ///< get normal onto face with direction
    double               GetOuterNormal(Uint face, Point3DCL& normal)              const;  ///< get outer normal onto face
    bool                 IsInTriang    (Uint TriLevel) const                               ///< check if tetra is in triangulation level
      { return GetLevel() == TriLevel || ( GetLevel() < TriLevel && IsUnrefined() ); }

    bool IsSane    (std::ostream&) const;                                                  ///< check for sanity
    void DebugInfo (std::ostream&) const;                                                  ///< get debug-information
};


///\brief A spatial tetra and a time interval t0 < t1.
struct TetraPrismCL
{
    const TetraCL& t;
    double t0,
           t1;

    TetraPrismCL (const TetraCL& targ, double t0arg, double t1arg)
        : t( targ), t0( t0arg), t1( t1arg) {}
};


#ifdef _PAR
/// \brief Cast a transferable object to a given Simplex type.
template <typename SimplexT>
void simplex_cast( const DiST::TransferableCL& t, SimplexT*& s)
{
    Assert( DiST::GetDim<SimplexT>() == t.GetDim(), DiST::ErrorCL("simplex_cast: dimension does not match for ", t.GetGID(), DiST::NoGID), ~0);
    s= dynamic_cast<SimplexT*>(const_cast<DiST::TransferableCL*>(&t));
}
#endif

/// \brief Factory pattern for generating simplices
/** This class is responsible for generating the simplices, store them in the
    right list, and register them (in the parallel version).
*/
class SimplexFactoryCL
{
private:
    /// \name Where to store the simplices
    //@{
    MG_VertexContT& vertices_;
    MG_EdgeContT&   edges_;
    MG_FaceContT&   faces_;
    MG_TetraContT&  tetras_;
    //@}

public:
#ifdef _PAR
    typedef DiST::RemoteDataCL::ProcListT ProcListT;
#endif
    /// \brief Constructor
    SimplexFactoryCL( MG_VertexContT& verts, MG_EdgeContT& edges, MG_FaceContT& faces, MG_TetraContT& tetras)
        : vertices_(verts), edges_(edges), faces_(faces), tetras_(tetras) {}
    ///\name Factory methods to generate all types of simplices. Basically, each function calls the constructor of the simplex, put the simplex
    ///      into the MG_Container, and calls the Register function of the RemoteDataListCL
    //{@
    /// \brief Create a vertex by coordinate and level; FileBuilderCL has to construct the Id_, too, thus it can optionally be set.
    VertexCL& MakeVertex( const Point3DCL& Point3D, Uint Level, IdCL<VertexCL> id= IdCL<VertexCL>(), const bool donotRegister=false);
    /// \brief Create an edge
    EdgeCL& MakeEdge( VertexCL* vp0, VertexCL* vp1, Uint Level, BndIdxT bnd0= NoBndC, BndIdxT bnd1= NoBndC, short int MFR=0, const bool donotRegister=false);
    /// \brief Create a face (normally, the serial version)
    FaceCL& MakeFace ( Uint Level, BndIdxT bnd= NoBndC, const bool donotRegister=false);
    /// \brief Create a face (normally, the parallel version)
    FaceCL& MakeFace ( Uint Level, const Point3DCL& v0, const Point3DCL& v1, const Point3DCL& v2, BndIdxT bnd= NoBndC, const bool donotRegister=false);
    /// \brief Create a tetrahedron of verts and parent; FileBuilderCL has to construct the Id_, too, thus it can optionally be set.
    TetraCL& MakeTetra( VertexCL*, VertexCL*, VertexCL*, VertexCL*, TetraCL*, IdCL<TetraCL> id= IdCL<TetraCL>(), const bool donotRegister=false);
#ifdef _PAR
    /// \brief Create a tetrahedron of verts and level, if no parent is available; FileBuilderCL has to construct the Id_, too, thus it can optionally be set.
    TetraCL& MakeTetra( VertexCL*, VertexCL*, VertexCL*, VertexCL*, TetraCL*, Uint, IdCL<TetraCL> id= IdCL<TetraCL>(), const bool donotRegister=false);

    /// \brief Make a copy of a simplex
    template <typename SimplexT>
    SimplexT& MakeCopy( const SimplexT&, const ProcListT&);
    //@}
#endif

    ///\name Factory methods to destroy simplices which are marked for removement. Such a simplex is erased from the MG_Container.
    //{@
    ///\brief Destroy tetras of a given level after unlinking them from their faces.
    void DestroyMarkedTetras( Uint Level);
    void DestroyMarkedVEFs  ( Uint Level);
    //@}
};


/// \name barycentric center of simplices
//@{
inline
Point3DCL GetBaryCenter(const VertexCL& v) { return v.GetCoord(); }
Point3DCL GetBaryCenter(const EdgeCL&);
Point3DCL GetBaryCenter(const FaceCL&);
Point3DCL GetBaryCenter(const TetraCL&);
Point3DCL GetBaryCenter(const TetraCL& t, Uint face);
//@}

/// \name "safe" calculation of barycentric center of simplices. The compiler attributes/pragmas ensure the correct order of summation
//@{
#if _MSC_VER > 1400
# pragma float_control(push)
# pragma float_control(precise, on)
#endif

Point3DCL
#if GCC_VERSION > 40305 && !__INTEL_COMPILER
__attribute__((__optimize__("no-associative-math")))
#endif
    ComputeBaryCenter(const Point3DCL& v0, const Point3DCL& v1);

Point3DCL
#if GCC_VERSION > 40305 && !__INTEL_COMPILER
__attribute__((__optimize__("no-associative-math")))
#endif
    ComputeBaryCenter(const Point3DCL& v0, const Point3DCL& v1, const Point3DCL& v2);

Point3DCL
#if GCC_VERSION > 40305 && !__INTEL_COMPILER
__attribute__((__optimize__("no-associative-math")))
#endif
    ComputeBaryCenter(const Point3DCL& v0, const Point3DCL& v1, const Point3DCL& v2, const Point3DCL& v3);

#if _MSC_VER > 1400
# pragma float_control(pop)
#endif
//@}

Point3DCL ComputeBaryCenter(const VertexCL& s);
Point3DCL ComputeBaryCenter(const EdgeCL& s);
Point3DCL ComputeBaryCenter(const FaceCL& s);
Point3DCL ComputeBaryCenter(const TetraCL& s);

Point3DCL GetWorldCoord(const TetraCL&, const SVectorCL<3>&);
Point3DCL GetWorldCoord(const TetraCL&, Uint face, const SVectorCL<2>&);
// barycentric coordinates:
Point3DCL GetWorldCoord(const TetraCL&, const SVectorCL<4>&);
Point3DCL GetRefCoord(const TetraCL&, const SVectorCL<4>&);

SVectorCL<3> FaceToTetraCoord(const TetraCL& t, Uint f, SVectorCL<2> c);

/// \brief Maps world-coordinates p to barycentric coordinates of the tetra t.
class World2BaryCoordCL
{
  private:
    QRDecompCL<4> qr_;

  public:
    World2BaryCoordCL (const TetraCL& t) { assign( t); }
    World2BaryCoordCL (const Point3DCL* coordVerts);
    World2BaryCoordCL () {};

    void assign (const TetraCL& t);

    /// \brief Maps a point p to M\inv*(p, 1)^T.
    BaryCoordCL operator() (const Point3DCL& p) const;

    /// \brief Maps a direction v to M\inv*(v, 0)^T. Given v= q - p, map_direction (v) = map(q) - map(p). 
    BaryCoordCL map_direction (const Point3DCL& p) const;
};

///\brief Compute world-coordinates from a tetra and barycentric coordinates.
/// This class is similar to GetWorldCoord( tetra, barycoord). It is better suited for many consecutive evaluations. It can be used in std-algorithms.
class Bary2WorldCoordCL
{
  private:
    SMatrixCL<3,4> mat_;

  public:
    typedef Point3DCL value_type;

    Bary2WorldCoordCL (const TetraCL& tet) : mat_( Uninitialized) { assign( tet); }
    Bary2WorldCoordCL () : mat_( Uninitialized) {}

    void assign (const TetraCL& tet);

    Point3DCL operator() (const BaryCoordCL& b) const { return mat_*b; }
};

///\brief Compute (space-time) world-coordinates from a tetra-prism and STCoord-coordinates.
/// This class is suited for many consecutive evaluations. It can be used in std-algorithms.
class STCoord2WorldCoordCL
{
  private:
    const Bary2WorldCoordCL space_mapper_;
    const double t0_,
                 dt_;

  public:
    typedef Point4DCL value_type;

    STCoord2WorldCoordCL (const TetraPrismCL& prism)
        : space_mapper_( prism.t), t0_( prism.t0), dt_( prism.t1 - prism.t0) {}

    Point3DCL space (const STCoordCL& b) const { return space_mapper_( b.x_bary); }
    double    time  (const STCoordCL& b) const { return t0_ + dt_*b.t_ref; }
    Point4DCL operator() (const STCoordCL& b) const {
        const Point3DCL& x( space( b));
        return MakePoint4D( x[0], x[1], x[2], time( b));
    }
};


inline SMatrixCL<4,4> parent_to_child_bary (Ubyte ch)
{
    SMatrixCL<4,4> ret( Uninitialized);
    const size_t pos= ch*16;
    for (Uint i= 0; i < 16; ++i)
        ret[i]= parent_to_child_bary_ar[pos + i];
    return ret;
}

inline SMatrixCL<4,4> child_to_parent_bary (Ubyte ch)
{
    SMatrixCL<4,4> ret( Uninitialized);
    const size_t pos= ch*16;
    for (Uint i= 0; i < 16; ++i)
        ret[i]= 0.5*child_to_parent_bary_ar[pos + i];
    return ret;
}


/// \brief Collect the faces of the the children of a refined tetra that refine a given face.
///
/// If the tetra is unrefined, the given face is returned.
/// \todo The mapping (refrule, face) -> (list of (child, face of child)) could be
///     precomputed and put in a table in geom/topo.{h,cpp}
/// \param p The parent tetra
/// \param f Face number in the parent (0..3)
/// \param childfaces Collects the pointers to the child faces that refine face f of p
void ComputeChildFacesOfFace (const TetraCL& p, Uint f, std::vector<const FaceCL*>& childfaces);

} // end of namespace DROPS

#include "../geom/simplex.tpp"
#endif
