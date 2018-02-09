/// \file reftetracut.h
/// \brief Triangulation of the reference tetraeder adapted to a linear level-set function.
/// \author LNM RWTH Aachen: Joerg Grande;

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
 * Copyright 2011,2012 LNM, RWTH Aachen, Germany
*/

#ifndef DROPS_REFTETRACUT_H
#define DROPS_REFTETRACUT_H

#include "misc/container.h"
#include "misc/staticinit.h"
#include "geom/signtraits.h"

namespace DROPS {

template <Uint Dim>
class RefPatchCL; ///< forward declaration for RefPatchBuilderCL.

template <Uint Dim>
class RefPatchBuilderCL
{
  public:
    typedef SArrayCL<Ubyte, Dim> FacetT; ///< the vertices of a facet of the cut: the body's vertices are denoted by 0..Dim, the edge-cuts by edge-num + Dim + 1.

    /// \brief Called by generic StaticInit in RefPatchCL.
    static void StaticInit (RefPatchCL<Dim> instances[SignTraitsCL<Dim>::num_pattern]); ///< undefined

    ///\brief maximal number of facets in any RefPatchCL.
    enum { max_num_facets= 0 };

    ///\brief Assign a sign pattern on the vertices; returns the value of RefPatchCL<Dim>::empty(). The generic version is undefined.
    static bool assign (const SignTraitsCL<Dim>&, RefPatchCL<Dim>&);
};

template <>
class RefPatchBuilderCL<3>
{
  public:
    typedef SArrayCL<Ubyte, 3> FacetT;

    enum { max_num_facets= 2 }; // Number of facets: 0, 1, or 2. The maximum comes from a quadrilateral.

  private:
    static FacetT MakeTriangle (Ubyte v0, Ubyte v1, Ubyte v2) { return MakeSArray( v0, v1, v2); }

  public:
    static void StaticInit ( RefPatchCL<3> instances[SignTraitsCL<3>::num_pattern]);
    static bool assign (const SignTraitsCL<3>&, RefPatchCL<3>&);
};

template <>
class RefPatchBuilderCL<4>
{
  public:
    typedef SArrayCL<Ubyte, 4> FacetT;

    enum { max_num_facets= 3 }; // Number of facets: 0, 1, 2 or 3. The maximum comes from a triangular prism.

  private:
    ///\brief The extended vertex-number (see SignTraitsCL) on the reference-penta, given the extended vertex-number n (in 0..9) on the facet tetra.
    static Ubyte PentaCutByTetraCut (Ubyte tetra, Ubyte n) {
        return n < 4 ? RefPenta::VertByTetraVert( tetra, n)
                     : RefPenta::EdgeByTetraEdge( tetra, n - 4) + 5;
    }
    static FacetT MakeTetra (Ubyte v0, Ubyte v1, Ubyte v2, Ubyte v3) { return MakeSArray( v0, v1, v2, v3); }
    static bool cone_construction (const SignTraitsCL<4>& cut, RefPatchCL<4>& p, Ubyte f); ///< produce tetras as the convex hulls of v=cut( 0) and the cut of the facet-tetra f. Returns, whether the cut with f is on the boundary of the f.

  public:
    static void StaticInit ( RefPatchCL<4> instances[SignTraitsCL<4>::num_pattern]);
    static bool assign (const SignTraitsCL<4>&, RefPatchCL<4>&);
};

///\brief The facets of the intersection of the reference-tetra with a linear levelset-function.
///
/// The class memoizes used sign-patterns if the triangulations are accessed via the instance( ls)-function. Individual instances may still be constructed (useful for debugging).
template <Uint Dim>
class RefPatchCL
{
  public:
    typedef RefPatchBuilderCL<Dim> BuilderT;
    typedef typename BuilderT::FacetT FacetT;
    typedef const FacetT*             const_facet_iterator;
    typedef       FacetT*             facet_iterator;

    enum { max_num_facets= BuilderT::max_num_facets };

    friend class RefPatchBuilderCL<Dim>;

    /// \brief Initializes RefPatchCL::instance_array_ (see below)  by calling RefPatchCL::assign() for all non-zero sign patterns.
    static void StaticInit () { BuilderT::StaticInit( RefPatchCL<Dim>::instance_array_); }
    static void StaticDestruct () {}

  private:
    FacetT facet_[max_num_facets]; ///< the facets
    Ubyte size_;                   ///< number of facets
    Ubyte is_boundary_facet_;      ///< true if the facet is one of the body's facets.

    static RefPatchCL instance_array_[SignTraitsCL<Dim>::num_pattern];

  public:
    RefPatchCL () : size_( static_cast<Ubyte>( -1)), is_boundary_facet_( 0) {} ///< Uninitialized default state
    RefPatchCL (const SignTraitsCL<Dim>& cut) { assign( cut); } ///< Initialize with sign pattern on the vertices
    bool assign (const SignTraitsCL<Dim>& cut) { return BuilderT::assign( cut, *this); } ///< Assign a sign pattern on the vertices; returns the value of empty(). The assignment is delegated to the builder, which ensures the right actions depending on Dim.

    bool  is_initialized () const { return size_ <= max_num_facets; } ///< True after assign(...)

    ///@{ Recommended access to the facets for a given sign-pattern; memoizes the result. The functions throw an error for the 0-sign-pattern.
    static inline const RefPatchCL& instance (const byte   ls[Dim + 1]);
    static inline const RefPatchCL& instance (const double ls[Dim + 1]);
    ///@}

    bool is_boundary_facet () const { return is_boundary_facet_ == 1; } ///< true, iff the facet is one of the tetra's faces.

    bool  empty () const { return size_ == 0; } ///< true, iff the area of the intersection is 0.
    size_t size () const { return size_; }      ///< Number of facets in 0..max_num_facets

    ///@{ Random-access to the facets
    const_facet_iterator facet_begin () const { return facet_; }
    const_facet_iterator facet_end   () const { return facet_ + size_; }
    ///@}
};

template <Uint Dim>
RefPatchCL<Dim> RefPatchCL<Dim>::instance_array_[SignTraitsCL<Dim>::num_pattern];

extern template RefPatchCL<3> RefPatchCL<3>::instance_array_[SignTraitsCL<3>::num_pattern];
extern template RefPatchCL<4> RefPatchCL<4>::instance_array_[SignTraitsCL<4>::num_pattern];

template <Uint Dim>
inline const RefPatchCL<Dim>&
RefPatchCL<Dim>::instance (const byte ls[Dim + 1])
{
    const byte idx= SignTraitsCL<Dim>::pattern_idx ( ls);
    if (idx == 0)
        throw DROPSErrCL( "RefPatchCL::instance: found a full body in the zero level set, grid is too coarse!");
    return RefPatchCL<Dim>::instance_array_[idx + SignTraitsCL<Dim>::zero_pattern_offset];
}

template <Uint Dim>
inline const RefPatchCL<Dim>&
RefPatchCL<Dim>::instance (const double ls[Dim + 1])
{
    byte ls_byte[Dim + 1];
    std::transform( ls + 0, ls + Dim + 1, ls_byte + 0, DROPS::sign);
    return RefPatchCL<Dim>::instance( ls_byte);
}

typedef RefPatchCL<3> RefTetraPatchCL;
typedef RefPatchCL<4> RefPentaPatchCL;


typedef SignTraitsCL<3> SignPatternTraitCL;


///\brief The tetras partition the positive and negative part of the reference-tetra with respect to a linear levelset-function ls.
///
/// The class memoizes used sign-patterns if the triangulations are accessed via the instance( ls)-function. Individual instances may still be constructed (useful for debugging).
///
/// Layout of the tetra-sequence: tetras_ has at most 6 members, of which at most 3 are positive and 3 are negative.
/// [tetras_ ... <= ... begin_ ... <= ... tetras_ + 3 ... <= ... end_ ... <= ...tetras_ + 6]
/// [begin_..tetras_+3) are the negative tetras, [tetras_ +3..end_) are the positive tetras.
class RefTetraPartitionCL
{
  public:
    typedef SArrayCL<Ubyte, 4> TetraT; ///< representation of a tetra of the partition via its vertices: (0..3): original vertices; (4..9): original edges + 4
    typedef const TetraT* const_tetra_iterator;
    typedef       TetraT*       tetra_iterator;

    /// \brief Initializes RefTetraPatchCL::instance_array_ (see below)  by calling RefTetraPartitionCL::assign() for all non-zero sign patterns.
    static void StaticInit ();
    static void StaticDestruct () {}

  private:
    TetraT tetras_[6]; ///< at most six tetras
    tetra_iterator begin_;
    tetra_iterator end_;

    TetraT MakeTetra (Ubyte v0, Ubyte v1, Ubyte v2, Ubyte v3) const { return  MakeSArray( v0, v1, v2, v3); }

    void AddTetra (Ubyte v0, Ubyte v1, Ubyte v2, Ubyte v3, int sign) ///< The sequences grow away from tetras_+3
        { (sign == -1 ? *--begin_ : *end_++)= MakeTetra( v0, v1, v2, v3); }
     ///\brief e,f,g are assumed to be the equally-oriented edges of the quadrilateral faces.
    void AddPrism (Ubyte e0, Ubyte e1, Ubyte f0, Ubyte f1, Ubyte g0, Ubyte g1, int sign) {
        AddTetra( e0, e1, f1, g1, sign);
        AddTetra( e0, f0, f1, g1, sign);
        AddTetra( e0, f0, g0, g1, sign);
    }

    ///\brief If the intersection is quadrilateral, this returns the first uncut edge.
    Ubyte first_uncut_edge (const SignPatternTraitCL& cut) const { return cut[0] == 1 ? 0 : (cut[1] == 2 ? 1 : 2); }
    Ubyte some_non_zero_vertex (const SignPatternTraitCL& cut) const;

    static RefTetraPartitionCL instance_array_[81]; // 81 = 3^4 = all possible sign-patterns on the vertices

  public:
    RefTetraPartitionCL () : begin_( tetras_ + 3), end_( tetras_) {} ///< Uninitialized default state
    RefTetraPartitionCL (const SignPatternTraitCL& cut) { assign( cut); } ///< Initialize with sign pattern on the vertices
    bool assign (const SignPatternTraitCL& cut); ///< Assign a sign pattern on the vertices; returns the value of is_uncut()
    bool is_initialized () const { return begin_ < end_; } ///< True after assign(...)

    ///@{ Recommended access to the triangles for a given sign-pattern; memoizes the result. The functions throw an error for the 0-sign-pattern.
    static inline const RefTetraPartitionCL& instance (const byte   ls[4]);
    static inline const RefTetraPartitionCL& instance (const double ls[4]);
    ///@}

    bool is_uncut () const { return end_ == begin_ + 1; } ///< True, iff the partition has exactly one tetra
    int sign (const_tetra_iterator t) const { return t < tetras_ + 3 ? -1 : 1; } ///< Sign of the tetra, to which t points

    Ubyte tetra_size (TetraSignEnum s= AllTetraC) const ///< number of tetras with given sign
         { return tetra_end( s) - tetra_begin( s); }

    ///@{ Random-access to the tetras: all tetras, or negative and positive tetras separately, see TetraSignEnum
    const_tetra_iterator tetra_begin (TetraSignEnum s= AllTetraC) const
        { return s == PosTetraC ? tetras_ + 3 : begin_; }
    const_tetra_iterator tetra_end (TetraSignEnum s= AllTetraC) const
        { return s == NegTetraC ? tetras_ + 3 : end_; }
    ///@}

    friend std::ostream& operator<< (std::ostream&, const RefTetraPartitionCL&); ///< Debug-output to a stream (dumps all members)
};

inline const RefTetraPartitionCL&
RefTetraPartitionCL::instance (const byte ls[4])
{
    const byte idx= SignTraitsCL<3>::pattern_idx ( ls);
    if (idx == 0)
        throw DROPSErrCL( "RefTetraPartitionCL::instance: found 3-dim. zero level set, grid is too coarse!");
    return instance_array_[idx + 40];
}

inline const RefTetraPartitionCL&
RefTetraPartitionCL::instance (const double ls[4])
{
    byte ls_byte[4];
    std::transform( ls + 0, ls + 4, ls_byte + 0, DROPS::sign);
    return instance( ls_byte);
}

/// \brief In cases that a triangle on the boundary is cut by the interface of two fluids, we want to partition this cut triangle into subtriangles for integration.
/// This class partitions the positive and negative part of the reference-triangle (a face of the tetra) with respect to a linear levelset-function.
/// If you add/change functionalities or refactor code of this class, please check the unit test case: tests/bndTrianglePartition.cpp
class RefTrianglePartitionCL
{
  public:
    typedef SArrayCL<Ubyte, 3> TriangleT;     ///< the vertices of a triangle of the cut: the tetra's vertices are denoted by 0..3, the edge-cuts by edge-num + 4, which is in 4..9.
    typedef const TriangleT* const_triangle_iterator;
    typedef       TriangleT*       triangle_iterator;
    
  private:
    TriangleT triangle_[3];      ///< at most three triangles
    Ubyte size_;                 ///< number of triangles
    int sign_[3];                ///< sign of the triangles

    TriangleT MakeTriangle (Ubyte v0, Ubyte v1, Ubyte v2) const { return MakeSArray( v0, v1, v2); }

  public:
    /// Setting the level-set value of the opposite vertex to the triangle to 0.  
    /// With RefTetraPartitionCL we can first make a partition of the refTetra in which the refTriangle lies, 
    /// then we just map this refTetra partition to the face we are interested in.
    RefTrianglePartitionCL(byte ls[4], Ubyte VertexNum); 
    size_t size () const { return size_; }      ///< Number of triangles, 0, 1, or 2
    int sign (const_triangle_iterator t) const { return sign_[ t- triangle_begin() ]; } ///< Sign of the triangle, to which t points

    ///@{ Random-access to the triangles
    const_triangle_iterator triangle_begin () const { return triangle_; }
    const_triangle_iterator triangle_end   () const { return triangle_ + size_; }
    ///@}
};

} // end of namespace DROPS

#endif
