/// \file multigrid.h
/// \brief classes that constitute the multigrid
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

/// TODO: Use information hiding, access control and const-qualification more
///       extensively to avoid accidental changes of the multigrid structure.

#ifndef DROPS_MULTIGRID_H
#define DROPS_MULTIGRID_H

#include "geom/simplex.h"
#include "num/bndData.h"

namespace DROPS
{

//**************************************************************************
// Classes that constitute a multigrid and helpers                         *
//**************************************************************************

class MultiGridCL;
class MGBuilderCL;

template <class SimplexT>
struct TriangFillCL;

template <class SimplexT>
class TriangCL
{
  public:
    typedef std::vector<SimplexT*> LevelCont;

    typedef SimplexT**                     ptr_iterator;
    typedef const SimplexT**         const_ptr_iterator;

    typedef ptr_iter<SimplexT>             iterator;
    typedef ptr_iter<const SimplexT> const_iterator;

  private:
    mutable std::vector<LevelCont> triang_;
    MultiGridCL&                   mg_;

    inline void MaybeCreate (int lvl) const;

  public:
    TriangCL (MultiGridCL& mg);

    void clear () { triang_.clear(); }
    size_t size  (int lvl= -1) const
        { MaybeCreate( lvl); return triang_[StdIndex( lvl)].size() - 1; }

    ptr_iterator begin (int lvl= -1)
        { MaybeCreate( lvl); return &*triang_[StdIndex( lvl)].begin(); }
    ptr_iterator end   (int lvl= -1)
        {
            MaybeCreate( lvl);
            return &*(triang_[StdIndex( lvl)].end() - 1);
        }

    const_ptr_iterator begin (int lvl= -1) const
        { MaybeCreate( lvl); return const_cast<const_ptr_iterator>( &*triang_[StdIndex( lvl)].begin()); }
    const_ptr_iterator end   (int lvl= -1) const
        {
            MaybeCreate( lvl);
            return const_cast<const_ptr_iterator>( &*( triang_[StdIndex( lvl)].end() - 1));
        }

    ///@{ Cave: The returned level-container contains a zero-pointer as last element, which serves as end-iterator for the end()-functions in this class.
    LevelCont&       operator[] (int lvl)
        { MaybeCreate( lvl); return triang_[StdIndex( lvl)]; }
    const LevelCont& operator[] (int lvl) const
        { MaybeCreate( lvl); return triang_[StdIndex( lvl)]; }
    ///@}

    inline int  StdIndex    (int lvl) const;

};

typedef  TriangCL<VertexCL> TriangVertexCL;
typedef  TriangCL<EdgeCL>   TriangEdgeCL;
typedef  TriangCL<FaceCL>   TriangFaceCL;
typedef  TriangCL<TetraCL>  TriangTetraCL;

/// \brief Type of functions used to identify points on periodic boundaries,
///     that share the same dof.
typedef bool (*match_fun)(const Point3DCL&, const Point3DCL&);


class BoundaryCL
/// \brief stores boundary segments and information on periodic boundaries (if some exist)
{
  friend class MGBuilderCL;

  public:
    enum BndType {
        Per1Bnd= 1,    ///< periodic boundary 1
        Per2Bnd= 2,    ///< periodic boundary 2
        OtherBnd= 0    ///< non-periodic boundary
    };

    typedef std::vector<BndSegCL*> SegPtrCont;
    typedef std::vector<BndType>   BndTypeCont;

  private:
    SegPtrCont          Bnd_;
    mutable BndTypeCont BndType_;
    mutable match_fun   match_;

  public:
    BoundaryCL() : match_(0) {}
    /// deletes the objects pointed to in Bnd_.
    ~BoundaryCL();

    const BndSegCL* GetBndSeg(BndIdxT idx)  const { return Bnd_[idx]; }
    BndIdxT         GetNumBndSeg()          const { return Bnd_.size(); }
    BndType         GetBndType(BndIdxT idx) const { return !BndType_.empty() ? BndType_[idx] : OtherBnd; }

    void      SetPeriodicBnd( const BndTypeCont& type, match_fun) const;
    match_fun GetMatchFun() const { return match_; }
    bool      Matching ( const Point3DCL& p, const Point3DCL& q) const { return match_(p,q); }
    bool      HasPeriodicBnd() const { return match_; }
};


#ifdef _PAR
// fwd declaration
class LbIteratorCL;
namespace DiST{
class InfoCL;
}
#endif

class ColorClassesCL; ///< forward declaration of the partitioning of the tetras in a triangulation into color classes

class MultiGridCL
{

  friend class MGBuilderCL;
#ifdef _PAR
  friend class ParMultiGridCL;
  friend class LbIteratorCL;
  friend class DiST::InfoCL;
#endif

  public:
    typedef MG_VertexContT VertexCont;
    typedef MG_EdgeContT   EdgeCont;
    typedef MG_FaceContT   FaceCont;
    typedef MG_TetraContT  TetraCont;

    typedef VertexCont::LevelCont VertexLevelCont;
    typedef EdgeCont::LevelCont   EdgeLevelCont;
    typedef FaceCont::LevelCont   FaceLevelCont;
    typedef TetraCont::LevelCont  TetraLevelCont;

    typedef VertexCont::LevelIterator             VertexIterator;
    typedef EdgeCont::LevelIterator               EdgeIterator;
    typedef FaceCont::LevelIterator               FaceIterator;
    typedef TetraCont::LevelIterator              TetraIterator;
    typedef VertexCont::const_LevelIterator const_VertexIterator;
    typedef EdgeCont::const_LevelIterator   const_EdgeIterator;
    typedef FaceCont::const_LevelIterator   const_FaceIterator;
    typedef TetraCont::const_LevelIterator  const_TetraIterator;

    typedef TriangVertexCL::iterator             TriangVertexIteratorCL;
    typedef TriangEdgeCL::iterator               TriangEdgeIteratorCL;
    typedef TriangFaceCL::iterator               TriangFaceIteratorCL;
    typedef TriangTetraCL::iterator              TriangTetraIteratorCL;
    typedef TriangVertexCL::const_iterator const_TriangVertexIteratorCL;
    typedef TriangEdgeCL::const_iterator   const_TriangEdgeIteratorCL;
    typedef TriangFaceCL::const_iterator   const_TriangFaceIteratorCL;
    typedef TriangTetraCL::const_iterator  const_TriangTetraIteratorCL;

  private:
    BoundaryCL Bnd_;
    VertexCont Vertices_;
    EdgeCont   Edges_;
    FaceCont   Faces_;
    TetraCont  Tetras_;

    TriangVertexCL TriangVertex_;
    TriangEdgeCL   TriangEdge_;
    TriangFaceCL   TriangFace_;
    TriangTetraCL  TriangTetra_;

    size_t           version_;                      ///< each modification of the multigrid increments this number
    SimplexFactoryCL factory_;                      ///< factory for generating simplices

    mutable std::map<int, ColorClassesCL*> colors_; // map: level -> Color-classes of the tetra for that level

#ifdef _PAR
    bool killedGhostTetra_;                         // are there ghost tetras, that are marked for removement
    bool IsLevelEmpty(Uint lvl)
        { return Vertices_[lvl].empty() && Edges_[lvl].empty() && Faces_[lvl].empty() && Tetras_[lvl].empty(); }
#endif

    void PrepareModify   () { Vertices_.PrepareModify(); Edges_.PrepareModify(); Faces_.PrepareModify(); Tetras_.PrepareModify(); }
    void FinalizeModify  () { Vertices_.FinalizeModify(); Edges_.FinalizeModify(); Faces_.FinalizeModify(); Tetras_.FinalizeModify(); IncrementVersion(); }
    void AppendLevel     () { Vertices_.AppendLevel(); Edges_.AppendLevel(); Faces_.AppendLevel(); Tetras_.AppendLevel(); }
    void RemoveLastLevel () { Vertices_.RemoveLastLevel(); Edges_.RemoveLastLevel(); Faces_.RemoveLastLevel(); Tetras_.RemoveLastLevel(); }

    void IncrementVersion() { ++version_; }                    ///< Increment version of the multigrid
    void ClearTriangCache ();

    void RestrictMarks (Uint Level) { std::for_each( Tetras_[Level].begin(), Tetras_[Level].end(), std::mem_fun_ref(&TetraCL::RestrictMark)); }
    void CloseGrid     (Uint);
    void UnrefineGrid  (Uint);
    void RefineGrid    (Uint);

  public:
    MultiGridCL (const MGBuilderCL& Builder);
    MultiGridCL (const MultiGridCL&); // Dummy
    // default ctor


    ~MultiGridCL()
    { 
#ifdef _PAR            
        DiST::InfoCL::Instance().Destroy();
#endif 
    }

    const BoundaryCL& GetBnd     () const { return Bnd_; }
    const VertexCont& GetVertices() const { return Vertices_; }
    const EdgeCont&   GetEdges   () const { return Edges_; }
    const FaceCont&   GetFaces   () const { return Faces_; }
    const TetraCont&  GetTetras  () const { return Tetras_; }

    const TriangVertexCL& GetTriangVertex () const { return TriangVertex_; }
    const TriangEdgeCL&   GetTriangEdge   () const { return TriangEdge_; }
    const TriangFaceCL&   GetTriangFace   () const { return TriangFace_; }
    const TriangTetraCL&  GetTriangTetra  () const { return TriangTetra_; }

    VertexIterator GetVerticesBegin (int Level=-1) { return Vertices_.level_begin( Level); }
    VertexIterator GetVerticesEnd   (int Level=-1) { return Vertices_.level_end( Level); }
    EdgeIterator   GetEdgesBegin    (int Level=-1)  { return Edges_.level_begin( Level); }
    EdgeIterator   GetEdgesEnd      (int Level=-1)  { return Edges_.level_end( Level); }
    FaceIterator   GetFacesBegin    (int Level=-1) { return Faces_.level_begin( Level); }
    FaceIterator   GetFacesEnd      (int Level=-1) { return Faces_.level_end( Level); }
    TetraIterator  GetTetrasBegin   (int Level=-1) { return Tetras_.level_begin( Level); }
    TetraIterator  GetTetrasEnd     (int Level=-1) { return Tetras_.level_end( Level); }
    const_VertexIterator GetVerticesBegin (int Level=-1) const { return Vertices_.level_begin( Level); }
    const_VertexIterator GetVerticesEnd   (int Level=-1) const { return Vertices_.level_end( Level); }
    const_EdgeIterator   GetEdgesBegin    (int Level=-1) const { return Edges_.level_begin( Level); }
    const_EdgeIterator   GetEdgesEnd      (int Level=-1) const { return Edges_.level_end( Level); }
    const_FaceIterator   GetFacesBegin    (int Level=-1) const { return Faces_.level_begin( Level); }
    const_FaceIterator   GetFacesEnd      (int Level=-1) const { return Faces_.level_end( Level); }
    const_TetraIterator  GetTetrasBegin   (int Level=-1) const { return Tetras_.level_begin( Level); }
    const_TetraIterator  GetTetrasEnd     (int Level=-1) const { return Tetras_.level_end( Level); }

    VertexIterator GetAllVertexBegin (int= -1     ) { return Vertices_.begin(); }
    VertexIterator GetAllVertexEnd   (int Level=-1) { return Vertices_.level_end( Level); }
    EdgeIterator   GetAllEdgeBegin   (int= -1     ) { return Edges_.begin(); }
    EdgeIterator   GetAllEdgeEnd     (int Level=-1) { return Edges_.level_end( Level); }
    FaceIterator   GetAllFaceBegin   (int= -1     ) { return Faces_.begin(); }
    FaceIterator   GetAllFaceEnd     (int Level=-1) { return Faces_.level_end( Level); }
    TetraIterator  GetAllTetraBegin  (int= -1     ) { return Tetras_.begin(); }
    TetraIterator  GetAllTetraEnd    (int Level=-1) { return Tetras_.level_end( Level); }
    const_VertexIterator GetAllVertexBegin (int= -1     ) const { return Vertices_.begin(); }
    const_VertexIterator GetAllVertexEnd   (int Level=-1) const { return Vertices_.level_end( Level); }
    const_EdgeIterator   GetAllEdgeBegin   (int= -1     ) const  { return Edges_.begin(); }
    const_EdgeIterator   GetAllEdgeEnd     (int Level=-1) const  { return Edges_.level_end( Level); }
    const_FaceIterator   GetAllFaceBegin   (int= -1     ) const { return Faces_.begin(); }
    const_FaceIterator   GetAllFaceEnd     (int Level=-1) const { return Faces_.level_end( Level); }
    const_TetraIterator  GetAllTetraBegin  (int= -1     ) const { return Tetras_.begin(); }
    const_TetraIterator  GetAllTetraEnd    (int Level=-1) const { return Tetras_.level_end( Level); }

    TriangVertexIteratorCL GetTriangVertexBegin (int Level=-1) { return TriangVertex_.begin( Level); }
    TriangVertexIteratorCL GetTriangVertexEnd   (int Level=-1) { return TriangVertex_.end( Level); }
    TriangEdgeIteratorCL   GetTriangEdgeBegin   (int Level=-1) { return TriangEdge_.begin( Level); }
    TriangEdgeIteratorCL   GetTriangEdgeEnd     (int Level=-1) { return TriangEdge_.end( Level); }
    TriangFaceIteratorCL   GetTriangFaceBegin   (int Level=-1) { return TriangFace_.begin( Level); }
    TriangFaceIteratorCL   GetTriangFaceEnd     (int Level=-1) { return TriangFace_.end( Level); }
    TriangTetraIteratorCL  GetTriangTetraBegin  (int Level=-1) { return TriangTetra_.begin( Level); }
    TriangTetraIteratorCL  GetTriangTetraEnd    (int Level=-1) { return TriangTetra_.end( Level); }
    const_TriangVertexIteratorCL GetTriangVertexBegin (int Level=-1) const { return TriangVertex_.begin( Level); }
    const_TriangVertexIteratorCL GetTriangVertexEnd   (int Level=-1) const { return TriangVertex_.end( Level); }
    const_TriangEdgeIteratorCL   GetTriangEdgeBegin   (int Level=-1) const { return TriangEdge_.begin( Level); }
    const_TriangEdgeIteratorCL   GetTriangEdgeEnd     (int Level=-1) const { return TriangEdge_.end( Level); }
    const_TriangFaceIteratorCL   GetTriangFaceBegin   (int Level=-1) const { return TriangFace_.begin( Level); }
    const_TriangFaceIteratorCL   GetTriangFaceEnd     (int Level=-1) const { return TriangFace_.end( Level); }
    const_TriangTetraIteratorCL  GetTriangTetraBegin  (int Level=-1) const { return TriangTetra_.begin( Level); }
    const_TriangTetraIteratorCL  GetTriangTetraEnd    (int Level=-1) const { return TriangTetra_.end( Level); }

    Uint GetLastLevel() const { return Tetras_.GetNumLevel()-1; }
    Uint GetNumLevel () const { return Tetras_.GetNumLevel(); }

    size_t GetVersion() const { return version_; }              ///< Get version of the multigrid

    void Refine();                                              // in parallel mode, this function uses a parallel version for refinement!

    void Scale( double);
    void Transform( Point3DCL (*mapping)(const Point3DCL&));
    void MakeConsistentNumbering();
    void SplitMultiBoundaryTetras();                            ///< Tetras adjacent to more than one boundary-segment are subdivided into four tetras using the barycenter. This method must be called prior to Refine or MakeConsistentNumbering.
    void SizeInfo(std::ostream&);                               // all procs have to call this function in parallel mode!
    void ElemInfo(std::ostream&, int Level= -1) const;          // all procs have to call this function in parallel mode
    void DebugInfo(std::ostream&, int Level=-1) const;          ///< Put all vertices, edges, faces, and tetras on the stream
#ifdef _PAR
    Uint GetNumDistributedObjects() const;                      // get number of distributed objects
    Uint GetNumTriangTetra(int Level=-1);                       // get number of tetras of a given level
    Uint GetNumTriangFace(int Level=-1);                        // get number of faces of a given level
    Uint GetNumDistributedFaces(int Level=-1);                  // get number of faces on processor boundary
    SimplexFactoryCL& GetSimplexFactory() { return factory_; }
    void MakeConsistentHashes();
#endif

    const ColorClassesCL& GetColorClasses (int Level, const BndCondCL& Bnd) const;

    bool IsSane (std::ostream&, int Level=-1) const;
};


class PeriodicEdgesCL
/// \brief handles edges on periodic boundaries.
///
/// This class is used by the refinement algorithm in MultiGridCL to accumulate the MFR counters on linked periodic edges.
/// This assures that periodic boundaries are matching after refinement.
{
  public:
    typedef std::pair<EdgeCL*,EdgeCL*>  IdentifiedEdgesT;
    typedef std::list<IdentifiedEdgesT> PerEdgeContT;
    typedef PerEdgeContT::iterator      iterator;
    typedef MultiGridCL::EdgeIterator   EdgeIterator;

  private:
    PerEdgeContT      list_;
    MultiGridCL&      mg_;

    /// recompute data structure
    void Recompute( EdgeIterator begin, EdgeIterator end);
    /// accumulate local MFR counters of periodic edges and store the sum in the MFR counter
    void Accumulate();
    /// delete all data
    void Shrink();

  public:
    PeriodicEdgesCL( MultiGridCL& mg) : mg_(mg) {}
    // standard dtor

    BoundaryCL::BndType GetBndType( const EdgeCL& e) const;
    void AccumulateMFR( int lvl);
    /// print out list of identified edges for debugging
    void DebugInfo(std::ostream&);
};

/// \brief Storage of independent set of tetrahedra for assembling
class ColorClassesCL
{
  public:
    typedef std::vector<const TetraCL*> ColorClassT;
    typedef std::vector<ColorClassT>::const_iterator const_iterator;

  private:
    std::vector<ColorClassT> colors_;

  public:
    ColorClassesCL () {}
    ColorClassesCL (const MultiGridCL& mg, Uint lvl, const BndCondCL& Bnd) {
        compute_color_classes( const_cast<MultiGridCL&>( mg), lvl, Bnd);
    }

    void compute_color_classes (MultiGridCL& mg, Uint lvl, const BndCondCL& Bnd);
    /// \brief Put all tetras in the same class --> perfect parallelization.
    /// To avoid update races in code using this, you must only update tetra-specific data.
    void make_single_color_class(MultiGridCL::const_TriangTetraIteratorCL begin,
                                 MultiGridCL::const_TriangTetraIteratorCL end);

    size_t num_colors () const { return colors_.size(); }
    const_iterator begin () const { return colors_.begin(); }
    const_iterator end   () const { return colors_.end(); }
};


template <class SimplexT>
struct TriangFillCL
{
  static void // not defined
  fill (MultiGridCL& mg, typename TriangCL<SimplexT>::LevelCont& c, int lvl);
};

template <>
struct TriangFillCL<VertexCL>
{
    static void fill (MultiGridCL& mg, TriangCL<VertexCL>::LevelCont& c, int lvl);
};

template <>
struct TriangFillCL<EdgeCL>
{
    static void fill (MultiGridCL& mg, TriangCL<EdgeCL>::LevelCont& c, int lvl);
};

template <>
struct TriangFillCL<FaceCL>
{
    static void fill (MultiGridCL& mg, TriangCL<FaceCL>::LevelCont& c, int lvl);
};

template <>
struct TriangFillCL<TetraCL>
{
    static void fill (MultiGridCL& mg, TriangCL<TetraCL>::LevelCont& c, int lvl);
};


#define DROPS_FOR_TRIANG_VERTEX( mg, lvl, it) \
for (DROPS::TriangVertexCL::iterator it( mg.GetTriangVertexBegin( lvl)), end__( mg.GetTriangVertexEnd( lvl)); it != end__; ++it)

#define DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) \
for (DROPS::TriangVertexCL::const_iterator it( mg.GetTriangVertexBegin( lvl)), end__( mg.GetTriangVertexEnd( lvl)); it != end__; ++it)


#define DROPS_FOR_TRIANG_EDGE( mg, lvl, it) \
for (DROPS::TriangEdgeCL::iterator it( mg.GetTriangEdgeBegin( lvl)), end__( mg.GetTriangEdgeEnd( lvl)); it != end__; ++it)

#define DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it) \
for (DROPS::TriangEdgeCL::const_iterator it( mg.GetTriangEdgeBegin( lvl)), end__( mg.GetTriangEdgeEnd( lvl)); it != end__; ++it)


#define DROPS_FOR_TRIANG_FACE( mg, lvl, it) \
for (DROPS::TriangFaceCL::iterator it( mg.GetTriangFaceBegin( lvl)), end__( mg.GetTriangFaceEnd( lvl)); it != end__; ++it)

#define DROPS_FOR_TRIANG_CONST_FACE( mg, lvl, it) \
for (DROPS::TriangFaceCL::const_iterator it( mg.GetTriangFaceBegin( lvl)), end__( mg.GetTriangFaceEnd( lvl)); it != end__; ++it)


#define DROPS_FOR_TRIANG_TETRA( mg, lvl, it) \
for (DROPS::TriangTetraCL::iterator it( mg.GetTriangTetraBegin( lvl)), end__( mg.GetTriangTetraEnd( lvl)); it != end__; ++it)

#define DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) \
for (DROPS::TriangTetraCL::const_iterator it( mg.GetTriangTetraBegin( lvl)), end__( mg.GetTriangTetraEnd( lvl)); it != end__; ++it)


class MGBuilderCL
{
  protected:
    Uint parnumLevel_; /// \todo There is only one user of this and of set_par_numlevel... remove it?

    MultiGridCL::VertexCont& GetVertices(MultiGridCL* MG_) const { return MG_->Vertices_; }
    MultiGridCL::EdgeCont&   GetEdges   (MultiGridCL* MG_) const { return MG_->Edges_; }
    MultiGridCL::FaceCont&   GetFaces   (MultiGridCL* MG_) const { return MG_->Faces_; }
    MultiGridCL::TetraCont&  GetTetras  (MultiGridCL* MG_) const { return MG_->Tetras_; }
    BoundaryCL::SegPtrCont&  GetBnd     (MultiGridCL* MG_) const { return MG_->Bnd_.Bnd_; }
    void PrepareModify  (MultiGridCL* MG_) const { MG_->PrepareModify(); }
    void FinalizeModify (MultiGridCL* MG_) const { MG_->FinalizeModify(); }
    void AppendLevel    (MultiGridCL* MG_) const { MG_->AppendLevel(); }
    void RemoveLastLevel(MultiGridCL* MG_) const { MG_->RemoveLastLevel(); }

  public:
    MGBuilderCL (Uint parnumLevel= 1);
    virtual ~MGBuilderCL()  {}

    void set_par_numlevel (Uint parnumLevel) { parnumLevel_=  parnumLevel; }

    virtual void buildBoundary(MultiGridCL* MG_) const = 0;

    /** In order to create a multigrid with MPI, the following strategy is
    used. Only the master process creates the multigrid and all other
    processes have to create an "empty" multigrid, i.e. they only create the
    level views and the boundary. */
    ///\brief Used on the master-process (and serially).
    virtual void build_ser_impl(MultiGridCL* mgp) const = 0;
    ///\brief Used on the non-masters in parallel. A default implementation is provided.
    virtual void build_par_impl(MultiGridCL* mgp) const;
    void build(MultiGridCL* mgp) const {
        if (MASTER)
            build_ser_impl( mgp);
        else
            build_par_impl( mgp);
    }
};

class LocatorCL;

/// \brief Class for combining a tetrahedra with its bary center
class LocationCL
{
  private:
    const TetraCL* Tetra_;
    SVectorCL<4>   Coord_;

  public:
    LocationCL()
        : Tetra_(0), Coord_() {}
    LocationCL(const TetraCL* t, const SVectorCL<4>& p)
        : Tetra_(t), Coord_(p) {}
    LocationCL(const LocationCL& loc)
        : Tetra_(loc.Tetra_), Coord_(loc.Coord_) {}

#ifndef _PAR
    bool IsValid() const                        ///< Check if the tetrahedra is set
        { return Tetra_; }
#else
    bool IsValid(Uint lvl) const                ///< Check if the tetrahedra is set
        { return Tetra_ && Tetra_->IsInTriang(lvl)/* && Tetra_->MayStoreUnk()*/; }
#endif
    const TetraCL& GetTetra() const             ///< Get a reference to the tetrahedra
        { return *Tetra_; }
    const SVectorCL<4>& GetBaryCoord() const    ///< Get bary center coordinates of the tetrahedra
        { return Coord_; }

    friend class LocatorCL;
};

/// \brief Find a tetrahedra that surrounds a given Point
class LocatorCL
{
  private:
    static bool InTetra(const SVectorCL<4>& b, double tol= 0)
        { return b[0] >= -tol && b[0] <= 1.+tol
              && b[1] >= -tol && b[1] <= 1.+tol
              && b[2] >= -tol && b[2] <= 1.+tol
              && b[3] >= -tol && b[3] <= 1.+tol; }

    static void MakeMatrix(const TetraCL& t, SMatrixCL<4,4>& M)
    {
        for (Uint j=0; j<4; ++j)
        {
            for (Uint i=0; i<3; ++i)
                M(i, j)= t.GetVertex(j)->GetCoord()[i];
            M(3, j)= 1.;
        }
    }

  public:
    // default ctor, copy-ctor, dtor, assignment-op

    /// \brief Locate a point with a given tetrahedra
    static void
    LocateInTetra(LocationCL&, Uint, const Point3DCL&, double tol= 0);
    /// \brief Find the tetrahedra that surounds a point
    static void
    Locate(LocationCL&, const MultiGridCL& MG, int, const Point3DCL&, double tol= 1e-14);
};
// inline functions

template <class SimplexT>
  TriangCL<SimplexT>::TriangCL (MultiGridCL& mg)
      : triang_( mg.GetNumLevel()), mg_( mg)
{}

template <class SimplexT>
  inline int
  TriangCL<SimplexT>::StdIndex(int lvl) const
{
    return lvl >= 0 ? lvl : lvl + mg_.GetNumLevel();
}

template <class SimplexT>
  inline void
  TriangCL<SimplexT>::MaybeCreate(int lvl) const
{
    const int level= StdIndex( lvl);
    Assert ( level >= 0 && level < static_cast<int>( mg_.GetNumLevel()),
        DROPSErrCL( "TriangCL::MaybeCreate: Wrong level."), DebugContainerC);
#   pragma omp critical(TriangCL_MaybeCreate)
  {
    if (triang_.size() != mg_.GetNumLevel()) {
        triang_.clear();
        triang_.resize( mg_.GetNumLevel());
    }
    if (triang_[level].empty()) {
        TriangFillCL<SimplexT>::fill( mg_, triang_[level], level);
        // Append a zero-pointer as explicit end-iterator of the sequence.
        triang_[level].push_back( 0);
    }
  }
}


void circumcircle(const TetraCL& t, Point3DCL& c, double& r);
void circumcircle(const TetraCL& t, Uint face, Point3DCL& c, double& r);



inline void GetTrafo( SMatrixCL<3,3>& T, const TetraCL & t)
{
    const Point3DCL& pt0= t.GetVertex(0)->GetCoord();
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            T(j,i)= t.GetVertex(i+1)->GetCoord()[j] - pt0[j];
}

/// calculates the transpose of the transformation  Tetra -> RefTetra
inline void GetTrafoTr( SMatrixCL<3,3>& T, double& det, const Point3DCL pt[4])
{
    double M[3][3];
    const Point3DCL& pt0= pt[0];
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            M[j][i]= pt[i+1][j] - pt0[j];
    det=   M[0][0] * (M[1][1]*M[2][2] - M[1][2]*M[2][1])
         - M[0][1] * (M[1][0]*M[2][2] - M[1][2]*M[2][0])
         + M[0][2] * (M[1][0]*M[2][1] - M[1][1]*M[2][0]);

    T(0,0)= (M[1][1]*M[2][2] - M[1][2]*M[2][1])/det;
    T(0,1)= (M[2][0]*M[1][2] - M[1][0]*M[2][2])/det;
    T(0,2)= (M[1][0]*M[2][1] - M[2][0]*M[1][1])/det;
    T(1,0)= (M[2][1]*M[0][2] - M[0][1]*M[2][2])/det;
    T(1,1)= (M[0][0]*M[2][2] - M[2][0]*M[0][2])/det;
    T(1,2)= (M[2][0]*M[0][1] - M[0][0]*M[2][1])/det;
    T(2,0)= (M[0][1]*M[1][2] - M[1][1]*M[0][2])/det;
    T(2,1)= (M[1][0]*M[0][2] - M[0][0]*M[1][2])/det;
    T(2,2)= (M[0][0]*M[1][1] - M[1][0]*M[0][1])/det;
}

/// calculates the transpose of the transformation  Tetra -> RefTetra
inline void GetTrafoTr( SMatrixCL<3,3>& T, double& det, const TetraCL& t)
{
    double M[3][3];
    const Point3DCL& pt0= t.GetVertex(0)->GetCoord();
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            M[j][i]= t.GetVertex(i+1)->GetCoord()[j] - pt0[j];
    det=   M[0][0] * (M[1][1]*M[2][2] - M[1][2]*M[2][1])
         - M[0][1] * (M[1][0]*M[2][2] - M[1][2]*M[2][0])
         + M[0][2] * (M[1][0]*M[2][1] - M[1][1]*M[2][0]);

    T(0,0)= (M[1][1]*M[2][2] - M[1][2]*M[2][1])/det;
    T(0,1)= (M[2][0]*M[1][2] - M[1][0]*M[2][2])/det;
    T(0,2)= (M[1][0]*M[2][1] - M[2][0]*M[1][1])/det;
    T(1,0)= (M[2][1]*M[0][2] - M[0][1]*M[2][2])/det;
    T(1,1)= (M[0][0]*M[2][2] - M[2][0]*M[0][2])/det;
    T(1,2)= (M[2][0]*M[0][1] - M[0][0]*M[2][1])/det;
    T(2,0)= (M[0][1]*M[1][2] - M[1][1]*M[0][2])/det;
    T(2,1)= (M[1][0]*M[0][2] - M[0][0]*M[1][2])/det;
    T(2,2)= (M[0][0]*M[1][1] - M[1][0]*M[0][1])/det;
}


void MarkAll (MultiGridCL&);
void UnMarkAll (MultiGridCL&);


class ParamCL; // forward declaration for read_PeriodicBoundaries.

/// \brief Read PeriodicMatching and the periodic boundary-segments from P and insert them into mg.Bnd_.
/// The key PeriodicMatching is optional; ommitting it or setting it to the empty string disables periodic matching.
/// The default for all boundary-segments is OtherBnd.
/// All other keys are interpreted as boundary-segment indices.
/// The values have the form "Per1BC" or "Per2BC".
void read_PeriodicBoundaries (MultiGridCL& mg, const ParamCL& P);

} // end of namespace DROPS

#endif
