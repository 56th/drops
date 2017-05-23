/// \file subtriangulation.h
/// \brief Triangulation of a principal-lattice of a tetra adapted to a piecewise linear level-set function.
/// \author LNM RWTH Aachen: Joerg Grande, Liang Zhang; SC RWTH Aachen:

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
 * Copyright 2011 LNM/SC RWTH Aachen, Germany
*/

#ifndef DROPS_SUBTRIANGULATION_H
#define DROPS_SUBTRIANGULATION_H

#include "misc/container.h"
#include "geom/principallattice.h"
#include "geom/reftetracut.h"

#include <vector>
#include <valarray>
//#include <tr1/unordered_map>


namespace DROPS {

/// \brief Common types used by TetraPartitionCL, SPatchCL and their helpers
namespace LatticePartitionTypesNS {

typedef SArrayCL<Uint, 5>          PentaT; ///< representation of a penta of the partition via its vertices: index in the vertex_-array
typedef std::vector<PentaT>        PentaContT;
typedef PentaContT::const_iterator const_penta_iterator;

typedef SArrayCL<Uint, 4>          TetraT; ///< representation of a tetra of the partition via its vertices: index in the vertex_-array
typedef std::vector<TetraT>        TetraContT;
typedef TetraContT::const_iterator const_tetra_iterator;

typedef SArrayCL<Uint, 3>             TriangleT; ///< representation of a triangle of the interface via its vertices: index in the vertex_-array
typedef std::vector<TriangleT>        TriangleContT;
typedef TriangleContT::const_iterator const_triangle_iterator;

typedef std::vector<BaryCoordCL>    VertexContT;
typedef VertexContT::const_iterator const_vertex_iterator;
typedef VertexContT::      iterator       vertex_iterator;

typedef std::vector<STCoordCL>        STVertexContT;
typedef STVertexContT::const_iterator const_stvertex_iterator;
typedef STVertexContT::      iterator       stvertex_iterator;

} // end of namespace DROPS::LatticePartitionTypesNS

///\brief sums the signs and returns fabs(sum) == \# dof. If true is returned, the interface does not intersect the domain of f up to subgrid-resolution.
template <class GridFunT>
  inline bool
  equal_signs (const GridFunT& f);


class TetraPartitionCL; ///< forward declaration for output-routines

///\brief declaration of debug output (neccessary due to friend declaration in TetraPartitionCL)
std::ostream&
operator<< (std::ostream&, const TetraPartitionCL&);

///\brief declaration of .vtu output (neccessary due to friend declaration in TetraPartitionCL)
void
write_paraview_vtu (std::ostream&, const TetraPartitionCL&, TetraSignEnum= AllTetraC);


///\brief Partition the principal lattice of a tetra t (n intervals on each edge) according to the sign of a levelset function ls.
///
/// The sub-tetras, which are cut by the interface are triangulated to match the interface. The values of the level set function on the vertices of the principal lattice must be prescribed. The sequence of all tetrahedra contains the negative tetrahedra as initial subsequence.
///
/// For the vertexes, there are two properties, which can be selected by the policy template-parameters of make_partition():
/// VertexPartitionPolicyT: Determines, how the vertices are ordered with respect to the sign of the levelset function. This is important for the generation of composite quadrature rules.
///     Unordered:   First the vertexes from the principal lattice, then all proper cut-vertexes.
///     Sorted:      First the negative vertexes, then the zero vertexes, then the positive vertexes. The vertexes of the negative and of the positive tetras potentially overlap in the zero vertexes.
///     Partitioned: Like sorted, and all zero vertexes are duplicated, such that the vertexes of the negative tetras and of the positive tetras are disjoint.
/// Unordered is fastest and appropriate for quadrature rules that use no vertexes of the triangulation. Quadrature rules that use vertexes need Sorted, if a continuous integrand is prescribed (and integrated on the negative or positive or both domains), and Partitioned, if the integrand is discontinuous.
///
/// VertexCutMergingPolicyT: Determines, how often cut vertexes are stored.
///     Duplicate: A cut-vertex is added for each tetra, on which it is discovered; fast, but leads to more vertices, which in turn leads to more dof for quadrature rules that use the vertexes of the partition.
///     Merge: The edge cuts are memoized for each edge and added only once -- leads to the minimal amount of cut vertexes.
class TetraPartitionCL
{
  public:
    typedef LatticePartitionTypesNS::TetraT               TetraT;
    typedef LatticePartitionTypesNS::TetraContT           TetraContT;
    typedef LatticePartitionTypesNS::const_tetra_iterator const_tetra_iterator;

    typedef LatticePartitionTypesNS::VertexContT           VertexContT;
    typedef LatticePartitionTypesNS::const_vertex_iterator const_vertex_iterator;

  private:
    TetraContT tetras_;          ///< All tetras of the partition.
    Uint       pos_tetra_begin_; ///< begin of the subsequence of positive tetras

    VertexContT vertexes_;         ///< All vertices of the partition. 
    Uint        pos_vertex_begin_; ///< begin of the subsequence of vertexes of positive tetras
    Uint        neg_vertex_end_;   ///< end of the subsequence of of vertexes of negative tetras

    template <template <Uint Dim> class VertexCutMergingPolicyT>
      const TetraT ///< Create a single sub-tetra and its vertexes
      make_sub_tetra (const RefTetraPartitionCL::TetraT& ref_tet, const PrincipalLatticeCL::TetraT& lattice_tet,
        const double lset[4], Uint lattice_num_vertexes, VertexCutMergingPolicyT<3>& edgecut);

  public:
    TetraPartitionCL () : pos_tetra_begin_( 0), pos_vertex_begin_( 0), neg_vertex_end_( 0) {} ///< Empty default-cut

    ///\brief Computes the partition of the principal lattice with num_intervals on each edge of the reference-tetra given the level set values in ls.
    template <class VertexPartitionPolicyT,
              template <Uint Dim> class VertexCutMergingPolicyT>
    void make_partition (const PrincipalLatticeCL& lat, const std::valarray<double>& ls);

    Uint tetra_size  (TetraSignEnum s= AllTetraC) const ///< number of tetras with given sign
         { return tetra_end( s) - tetra_begin( s); }
    Uint vertex_size (TetraSignEnum s= AllTetraC) const ///< number of vertexes used by the tetras of the corresponding sign; depending on VertexCutMergingPolicyT, interface vertexes occur multiple times.
         { return vertex_end( s) - vertex_begin( s); }

    int sign (const_tetra_iterator t) const { return t < tetra_end( NegTetraC) ? -1 : 1; } ///< Sign of the tetra, to which t points

    /// Random-access to the tetras: all tetras, or negative and positive tetras separately, see TetraSignEnum.
    ///@{
    const_tetra_iterator tetra_begin (TetraSignEnum s= AllTetraC) const
        { return tetras_.begin() + (s == PosTetraC ? pos_tetra_begin_ : 0); }
    const_tetra_iterator tetra_end   (TetraSignEnum s= AllTetraC) const
        { return s == NegTetraC ? tetras_.begin() + pos_tetra_begin_ : tetras_.end(); }
    ///@}
    /// Random-access to the vertices.
    ///@{
    const_vertex_iterator vertex_begin (TetraSignEnum s= AllTetraC) const
        { return vertexes_.begin() + (s == PosTetraC ? pos_vertex_begin_ : 0); }
    const_vertex_iterator vertex_end   (TetraSignEnum s= AllTetraC) const
        { return s == NegTetraC ? vertexes_.begin() + neg_vertex_end_ : vertexes_.end(); }
    ///@}

    friend std::ostream& operator<< (std::ostream&, const TetraPartitionCL&); ///< Debug-output to a stream (dumps all members)
    friend void write_paraview_vtu (std::ostream&, const TetraPartitionCL&, TetraSignEnum);  ///< Debug-output to a stream: VTU-format for paraview.
};


/// @{ forward declarations from simplex.h, which is only included in subtriangulation.cpp.
class TetraCL;
struct TetraPrismCL;
class Bary2WorldCoordCL;
class STCoord2WorldCoordCL;
/// @}


template <Uint Dim>
struct DimensionTraitsCL
{
};

template <>
struct DimensionTraitsCL<3>
{
    typedef LatticePartitionTypesNS::TriangleT               FacetT;
    typedef LatticePartitionTypesNS::TriangleContT           FacetContT;
    typedef LatticePartitionTypesNS::const_triangle_iterator const_facet_iterator;

    typedef BaryCoordCL                                    VertexT;
    typedef LatticePartitionTypesNS::VertexContT           VertexContT;
    typedef LatticePartitionTypesNS::const_vertex_iterator const_vertex_iterator;
    typedef LatticePartitionTypesNS::      vertex_iterator       vertex_iterator;

    typedef PrincipalLatticeCL         LatticeT;
    typedef PrincipalLatticeCL::TetraT LatticeBodyT;

    typedef RefPatchCL<3>         RefPatchT;
    typedef RefPatchCL<3>::FacetT RefPatchFacetT;

    typedef TetraCL WorldBodyT;

    typedef Point3DCL                        WorldVertexT;
    typedef std::vector<Point3DCL>           WorldVertexContT;
    typedef WorldVertexContT::const_iterator const_world_vertex_iterator;

    typedef Bary2WorldCoordCL VertexToWorldVertexMapperT;

    typedef WorldVertexContT NormalContT;
    typedef NormalContT::const_iterator const_normal_iterator;
};

template <>
struct DimensionTraitsCL<4>
{
    typedef LatticePartitionTypesNS::TetraT               FacetT;
    typedef LatticePartitionTypesNS::TetraContT           FacetContT;
    typedef LatticePartitionTypesNS::const_tetra_iterator const_facet_iterator;

    typedef STCoordCL                                        VertexT;
    typedef LatticePartitionTypesNS::STVertexContT           VertexContT;
    typedef LatticePartitionTypesNS::const_stvertex_iterator const_vertex_iterator;
    typedef LatticePartitionTypesNS::      stvertex_iterator       vertex_iterator;

    typedef TetraPrismLatticeCL         LatticeT;
    typedef TetraPrismLatticeCL::PentaT LatticeBodyT;

    typedef RefPatchCL<4>         RefPatchT;
    typedef RefPatchCL<4>::FacetT RefPatchFacetT;

    typedef TetraPrismCL WorldBodyT;

    typedef Point4DCL                        WorldVertexT;
    typedef std::vector<Point4DCL>           WorldVertexContT;
    typedef WorldVertexContT::const_iterator const_world_vertex_iterator;

    typedef STCoord2WorldCoordCL VertexToWorldVertexMapperT;

    typedef WorldVertexContT NormalContT;
    typedef NormalContT::const_iterator const_normal_iterator;
};

template <Uint Dim>
class SPatchCL; ///< forward declaration

///\brief Debug-output to a stream: VTU-format for paraview (fwd-declaration neccessary due to friend declaration in TetraPartitionCL)
void write_paraview_vtu (std::ostream&, const SPatchCL<3>&);

///\brief Partition the bodies of LatticeT according to the sign of a levelset function ls. This class computes the simplicial facets of the resulting piecewise linear interface.
/// VertexCutMergingPolicyT: Determines, how often cut vertexes are stored.
///     Duplicate: A cut-vertex is added for each tetra, on which it is discovered; fast, but leads to more vertices, which in turn leads to more dof for quadrature rules that use the vertexes of the partition.
///     Merge: The edge cuts are memoized for each edge and added only once -- leads to the minimal amount of cut vertexes.
template <Uint Dim>
class SPatchCL
{
  public:
    typedef DimensionTraitsCL<Dim> DimTraitsT;

    typedef typename DimTraitsT::FacetT               FacetT;
    typedef typename DimTraitsT::FacetContT           FacetContT;
    typedef typename DimTraitsT::const_facet_iterator const_facet_iterator;

    typedef typename DimTraitsT::VertexT               VertexT;
    typedef typename DimTraitsT::VertexContT           VertexContT;
    typedef typename DimTraitsT::const_vertex_iterator const_vertex_iterator;

    typedef typename DimTraitsT::LatticeT        LatticeT;
    typedef typename DimTraitsT::LatticeBodyT LatticeBodyT;

    typedef typename DimTraitsT::RefPatchT      RefPatchT;
    typedef typename DimTraitsT::RefPatchFacetT RefPatchFacetT;

    typedef typename DimTraitsT::WorldBodyT WorldBodyT;

    typedef typename DimTraitsT::WorldVertexT                WorldVertexT;
    typedef typename DimTraitsT::WorldVertexContT            WorldVertexContT;
    typedef typename DimTraitsT::const_world_vertex_iterator const_world_vertex_iterator;

    typedef typename DimTraitsT::VertexToWorldVertexMapperT VertexToWorldVertexMapperT;

    typedef typename DimTraitsT::NormalContT NormalContT;
    typedef typename DimTraitsT::const_normal_iterator const_normal_iterator;

    typedef std::vector<double>::const_iterator const_absdet_iterator;

  private:
    typedef std::pair<Uint, Uint> RenumberVertexPairT; ///< Helper type to handle zero-vertexes

    FacetContT        facets_;            ///< All facets of the interface.
    std::vector<bool> is_boundary_facet_; ///< True, iff the facet is a facet of one of the simplexes of the principal lattice.

    VertexContT      vertexes_;
    mutable WorldVertexContT world_vertexes_;

    mutable NormalContT         normals_;
    mutable std::vector<double> absdets_;

    template <template <Uint> class VertexCutMergingPolicyT>
      const FacetT ///< Create a single sub-facet and its vertexes
      make_sub_facet (const RefPatchFacetT& ref_tri, const LatticeBodyT& lattice_tet,
        const LatticeT& lattice, const double lset[Dim],
        std::vector<Uint>& copied_vertexes, std::vector<RenumberVertexPairT>& renumber_zero_verts,
        VertexCutMergingPolicyT<Dim>& edgecut);

  public:
    /// Empty default-interface

    ///\brief Computes the piecewise triangular interface for the principal lattice with num_intervals on each edge of the reference-tetra given the level set values in ls.
    template <template <Uint> class VertexCutMergingPolicyT>
    void make_patch (const LatticeT& lat, const std::valarray<double>& ls);

    /// True, iff the facet is a facet of one of the bodies of the principal lattice.
    ///@{
    bool is_boundary_facet (Uint i) const { return is_boundary_facet_[i]; }
    bool is_boundary_facet (const_facet_iterator it) const { return is_boundary_facet_[it - facets_.begin()]; }
    ///@}

    Uint facet_size  () const ///< number of triangles
         { return facets_.size(); }
    Uint vertex_size () const ///< number of vertexes
         { return vertexes_.size(); }
    bool empty () const { return facets_.empty(); } ///< True, iff there is no surface patch
    void clear (); ///< reset to empty default-state

    /// Random-access to the tetras and vertices.
    ///@{
    const_facet_iterator facet_begin () const { return facets_.begin(); }
    const_facet_iterator facet_end   () const { return facets_.end(); }
    const_vertex_iterator vertex_begin () const { return vertexes_.begin(); }
    const_vertex_iterator vertex_end   () const { return vertexes_.end(); }
    ///@}

    bool world_vertex_empty () const { return world_vertexes_.empty(); }
    void compute_world_vertexes (const WorldBodyT& wb) const;
    const_world_vertex_iterator world_vertex_begin () const { return world_vertexes_.begin(); }
    const_world_vertex_iterator world_vertex_end   () const { return world_vertexes_.end(); }

    bool absdets_empty () const { return absdets_.empty(); }
    void compute_absdets (const WorldBodyT& wb) const;
    const_absdet_iterator absdet_begin () const { return absdets_.begin(); }
    const_absdet_iterator absdet_end   () const { return absdets_.end(); }

    bool normal_empty () const { return normals_.empty(); }
    // Also computes absdets.
    void compute_normals (const WorldBodyT& wb) const;
    const_normal_iterator normal_begin () const { return normals_.begin(); }
    const_normal_iterator normal_end   () const { return normals_.end(); }

#   pragma GCC diagnostic ignored "-Wnon-template-friend"
    friend void write_paraview_vtu (std::ostream&, const SPatchCL<Dim>&);
//#   pragma GCC diagnostic pop
};

typedef SPatchCL<3> SurfacePatchCL;

/// \brief In cases that a triangle on the boundary is cut by the interface of two fluids, we want to partition this cut triangle into subtriangles for integration.
/// This class partitions the cut triangle by the interface on the boundary to sub-triangles
/// Todo: TetraSignEnum is confusing for understanding the code, rename to TriangleSignEnum maybe
/// If you add/change functionalities or refactor code of this class, please check the unit test case: tests/bndTrianglePartition.cpp
class BndTriangPartitionCL
{
  public:
    typedef LatticePartitionTypesNS::TriangleT               TriangleT;
    typedef LatticePartitionTypesNS::TriangleContT           TriangleContT;
    typedef LatticePartitionTypesNS::const_triangle_iterator const_triangle_iterator;

    typedef LatticePartitionTypesNS::VertexContT           VertexContT;
    typedef LatticePartitionTypesNS::const_vertex_iterator const_vertex_iterator;

  private:
    TriangleContT     triangles_;               ///< All triangles of the cut boundary face in the order: negative triangles, positive triangles
    Uint       pos_triangles_begin_;            ///< begin of the subsequence of positive triangles

    VertexContT vertexes_;                      /// Store all the vertices in the tetra and vertices generated by the edge cut on the BndTriangle

    template <class VertexCutMergingPolicyT>
      const TriangleT ///< Create a single sub-triangle and its vertexes
      make_sub_triangle (const RefTrianglePartitionCL::TriangleT& ref_tri, const PrincipalLatticeCL::TetraT& lattice_tet,
        const double lset[4], Uint lattice_num_vertexes,  VertexCutMergingPolicyT& edgecut);

  public:
    /// Empty default-interface

    ///\brief partition a cut face on the special boundary 
    template <class VertexPartitionPolicyT, template <Uint Dim> class VertexCutMergingPolicyT>
    void make_partition2D (const PrincipalLatticeCL& lat, Uint face, const std::valarray<double>& ls);

    Uint triangle_size  (TetraSignEnum s= AllTetraC) const ///< number of triangles
         {  return triangle_end(s)-triangle_begin(s); }
    Uint vertex_size () const                              ///< number of vertexes
         { return   vertex_end()-vertex_begin(); }

    /// Random-access to the triangles and vertices, order of vertices and triangles: first the negative, then the positve
    ///@{
    const_triangle_iterator triangle_begin (TetraSignEnum s= AllTetraC) const 
        { return triangles_.begin() + (s == PosTetraC ? pos_triangles_begin_ : 0); }
    const_triangle_iterator triangle_end   (TetraSignEnum s= AllTetraC) const 
        { return s == NegTetraC ? triangles_.begin() + pos_triangles_begin_ : triangles_.end(); }

    // Vertices are not sorted like TetraPartitionCL
    const_vertex_iterator vertex_begin () const
        { return vertexes_.begin(); }
    const_vertex_iterator vertex_end   () const
        { return vertexes_.end();   }
    ///@}
};


/// \brief Vertices are not ordered with respect to the sign of the levelset function: First the vertexes from the principal lattice, then all proper cut-vertexes.
class UnorderedVertexPolicyCL
{
  private:
    typedef LatticePartitionTypesNS::TetraContT  TetraContT;
    typedef LatticePartitionTypesNS::VertexContT VertexContT;

  public:
    UnorderedVertexPolicyCL (const PrincipalLatticeCL&, const std::valarray<double>&,
        TetraContT::iterator, TetraContT::iterator, Uint) {}

    /// \brief Append the cut_vertexes to vertexes.
    /// pos_vertex_begin_ and neg_vertex_end_ are set to zero to indicate, that the vertexes are not ordered with respect to the sign of the level set function
    void sort_vertexes (VertexContT&, VertexContT&, Uint&, Uint&);
};

/// \brief Vertices are ordered with respect to the sign of the levelset function as follows: First the negative vertexes, then the zero vertexes, then the positive vertexes. The vertexes of the negative and of the positive tetras potentially overlap in the zero vertexes.
class SortedVertexPolicyCL
{
  private:
    typedef LatticePartitionTypesNS::TetraContT  TetraContT;
    typedef LatticePartitionTypesNS::VertexContT VertexContT;

    const PrincipalLatticeCL&    lat_;
    const std::valarray<double>& ls_;
    const TetraContT::iterator   tetra_begin_,
                                 tetra_end_;

  public:
    SortedVertexPolicyCL (const PrincipalLatticeCL& lat, const std::valarray<double>& ls,
        TetraContT::iterator tetra_begin, TetraContT::iterator tetra_end, Uint)
        : lat_( lat), ls_( ls), tetra_begin_( tetra_begin), tetra_end_(tetra_end) {}

    /// \brief Sort the vertexes and update the vertex numbers in the tetras.
    void sort_vertexes (VertexContT& vertexes, VertexContT& cut_vertexes,
        Uint& pos_vertex_begin, Uint& neg_vertex_end);
};

/// \brief Vertices are ordered with respect to the sign of the levelset function as follows: First the negative vertexes, then the zero vertexes, then the positive vertexes. Opposite to SortedVertexPolicyC, all zero vertexes are duplicated, such that the vertexes of the negative tetras and of the positive tetras are disjoint.
class PartitionedVertexPolicyCL
{
  private:
    typedef LatticePartitionTypesNS::TetraContT  TetraContT;
    typedef LatticePartitionTypesNS::VertexContT VertexContT;

    SortedVertexPolicyCL       pol_;
    const TetraContT::iterator tetra_begin_,
                               tetra_end_;
    Uint                       pos_tetra_begin_;

  public:
    PartitionedVertexPolicyCL (const PrincipalLatticeCL& lat, const std::valarray<double>& ls,
        TetraContT::iterator tetra_begin, TetraContT::iterator tetra_end, Uint pos_tetra_begin)
        : pol_( lat, ls, tetra_begin, tetra_end, pos_tetra_begin),
          tetra_begin_( tetra_begin), tetra_end_( tetra_end), pos_tetra_begin_( pos_tetra_begin) {}

    /// \brief Sort the vertexes and update the vertex numbers in the tetras: Special care must be taken for the duplicated vertexes
    void sort_vertexes (VertexContT& vertexes, VertexContT& cut_vertexes,
        Uint& pos_vertex_begin, Uint& neg_vertex_end);
};

///\brief A cut-vertex is added to the list of vertexes for each tetra, on which it is discovered; fast, but leads to more vertices, which in turn leads to more dof for quadrature rules that use the vertexes of the partition.
template <Uint Dim>
class DuplicateCutPolicyCL
{
  private:
    typedef DimensionTraitsCL<Dim> DimTraitsT;

    typedef typename DimTraitsT::VertexT               VertexT;
    typedef typename DimTraitsT::VertexContT           VertexContT;
    typedef typename DimTraitsT::const_vertex_iterator const_vertex_iterator;

    typedef typename DimTraitsT::LatticeT        LatticeT;

    const typename LatticeT::const_vertex_iterator lattice_vertexes_;
    VertexContT vertexes_;

  public:
    DuplicateCutPolicyCL (typename LatticeT::const_vertex_iterator lattice_vertexes)
        : lattice_vertexes_( lattice_vertexes) {}

    ///\brief Add the cut vertex and return its number.
    Uint operator() (Uint v0, Uint v1, double ls0, double ls1) {
        const double edge_bary1_cut= ls0/(ls0 - ls1); // the root of the level set function on the edge
        vertexes_.push_back( ConvexComb( edge_bary1_cut, lattice_vertexes_[v0], lattice_vertexes_[v1]));
        return vertexes_.size() - 1;
    }

    VertexContT& cut_vertex_container () { return vertexes_; }
};

///\brief A cut-vertex is added to the list of vertexes only by the first tetra, on which it is discovered: cuts are memoized for each edge.
template <Uint Dim>
class MergeCutPolicyCL
{
  private:
    typedef DimensionTraitsCL<Dim> DimTraitsT;

    typedef typename DimTraitsT::VertexT               VertexT;
    typedef typename DimTraitsT::VertexContT           VertexContT;
    typedef typename DimTraitsT::const_vertex_iterator const_vertex_iterator;

    typedef typename DimTraitsT::LatticeT        LatticeT;

    struct UintPairHasherCL
    {
        size_t operator() (const std::pair<Uint, Uint>& p) const
            { return p.first << (4*sizeof(Uint)) ^ p.second; } // for less than 2^(sizeof(Uint)/2) vertices this is a bijection
    };

    typedef std::pair<Uint, Uint> EdgeT;
    typedef DROPS_STD_UNORDERED_MAP<EdgeT, Uint, UintPairHasherCL> EdgeToCutMapT;

    const typename LatticeT::const_vertex_iterator lattice_vertexes_;
    VertexContT vertexes_;
    EdgeToCutMapT edge_to_cut_;

  public:
    MergeCutPolicyCL (typename LatticeT::const_vertex_iterator lattice_vertexes)
        : lattice_vertexes_( lattice_vertexes) {}

    ///\brief Return the number of the cut vertex, if it is already memoized, otherwise add it and return its number.
    Uint operator() (Uint v0, Uint v1, double ls0, double ls1) {
        const EdgeT e= v0 < v1 ? std::make_pair( v0, v1) : std::make_pair( v1, v0);
        typename EdgeToCutMapT::const_iterator e_it= edge_to_cut_.find( e);
        if (e_it != edge_to_cut_.end())
            return e_it->second;
        else {
            const double edge_bary1_cut= ls0/(ls0 - ls1); // the root of the level set function on the edge
            vertexes_.push_back( ConvexComb( edge_bary1_cut, lattice_vertexes_[v0], lattice_vertexes_[v1]));
            return edge_to_cut_[e]= vertexes_.size() - 1;
        }
    }

    VertexContT& cut_vertex_container () { return vertexes_; }
};

} // end of namespace DROPS

#include "geom/subtriangulation.tpp"

#endif
