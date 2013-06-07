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

namespace DROPS {

///\brief Represents the reference tetra, which is cut by a linear level set function ls. The values of the latter are prescribed on the vertices.
class SignPatternTraitCL
{
  private:
    Ubyte num_root_vert_;  ///< number of vertices, where the level set function is zero.
    Ubyte num_root_;       ///< number of roots of the level set function; invariant: num_root_vert <= num_root
    byte sign_[4];         ///< Sign of the level set function of the vertices; \f$\in\{-1,0,1\}\f$
    Ubyte cut_simplex_[4]; ///< local number with respect to the reference tetra of the object on the cut: [0,num_root_vert): vertex numbers in (0..3); [num_root_vert, num_root): edge numbers in (0..5). Both parts are sorted in increasing order.
    Ubyte cut_simplex_rep_[4]; ///< local number of the object on the cut: (0..9)

    void compute_cuts ();

  public:
    SignPatternTraitCL () : num_root_vert_( 4), num_root_( 4) {} ///< Uninitialized default state
    SignPatternTraitCL (const byte   ls[4]) { assign( ls); } ///< Assign the sign pattern on the vertices.
    SignPatternTraitCL (const double ls[4]) { assign( ls); } ///< Assign a sign pattern on the vertices.
    void assign (const byte   ls[4]); ///< Assign a sign pattern on the vertices.
    void assign (const double ls[4]); ///< Assign a sign pattern on the vertices.

    bool empty () const { return num_root_ == 0; } ///< True, iff there is no intersection.
    bool is_2d () const { return num_root_ > 2; }  ///< True, iff the intersection has positive area.
    bool is_3d () const { return num_root_vert_ == 4; }
    bool no_zero_vertex () const { return num_root_vert_ == 0; } ///< True, iff there is no vertex, in which ls vanishes.

    Ubyte num_cut_simplexes () const { return num_root_; } ///< Number of edges and vertices with a root of ls.
    Ubyte num_zero_vertexes () const { return num_root_vert_; } ///< Number of vertices of the tetra that are roots of ls.

    byte sign (int i) const { return sign_[i]; } ///< -1,0,1; sign of vertex i.

    /// Return local number of edges/verts with a root of ls. For edges, [] returns a edge number in 0..5, and () returns an extended vertex number in 4..9.
    ///@{
    Ubyte operator[] (int i) const { return cut_simplex_[i]; }
    Ubyte operator() (int i) const { return cut_simplex_rep_[i]; }
    ///@}

    friend std::ostream& operator<< (std::ostream&, const SignPatternTraitCL&); ///< Debug-output to a stream (dumps all members)
};

inline Ubyte
num_triangles (const SignPatternTraitCL& cut)
{
    return cut.is_2d() ? cut.num_cut_simplexes() - 2 : 0;
}

///\brief The triangles of the intersection of the reference-tetra with a linear levelset-function.
///
/// The class memoizes used sign-patterns if the triangulations are accessed via the instance( ls)-function. Individual instances may still be constructed (useful for debugging).
class RefTetraPatchCL
{
  public:
    typedef SArrayCL<Ubyte, 3> TriangleT; ///< the vertices of a triangle of the cut: the tetra's vertices are denoted by 0..3, the edge-cuts by edge-num + 4, which is in 4..9.
    typedef const TriangleT* const_triangle_iterator;
    typedef       TriangleT*       triangle_iterator;

    /// \brief Initializes RefTetraPatchCL::instance_array_ (see below)  by calling RefTetraPatchCL::assign() for all non-zero sign patterns.
    static void StaticInit ();
    static void StaticDestruct () {}

  private:
    TriangleT triangle_[2];      ///< at most two triangles
    Ubyte size_;                 ///< number of triangles
    Ubyte is_boundary_triangle_; ///< true if the triangle is one of the tetra's faces.

    TriangleT MakeTriangle (Ubyte v0, Ubyte v1, Ubyte v2) const { return MakeSArray( v0, v1, v2); }

    static RefTetraPatchCL instance_array_[81]; // 81 = 3^4 = all possible sign-patterns on the vertices

  public:
    RefTetraPatchCL () : size_( static_cast<Ubyte>( -1)), is_boundary_triangle_( 0) {} ///< Uninitialized default state
    RefTetraPatchCL (const SignPatternTraitCL& cut) { assign( cut); } ///< Initialize with sign pattern on the vertices
    bool assign (const SignPatternTraitCL& cut); ///< Assign a sign pattern on the vertices; returns the value of empty()

    bool  is_initialized () const { return size_ <= 2; } ///< True after assign(...)

    ///@{ Recommended access to the triangles for a given sign-pattern; memoizes the result. The functions throw an error for the 0-sign-pattern.
    static inline const RefTetraPatchCL& instance (const byte   ls[4]);
    static inline const RefTetraPatchCL& instance (const double ls[4]);
    ///@}

    bool is_boundary_triangle () const { return is_boundary_triangle_ == 1; } ///< true, iff the triangle is one of the tetra's faces.

    bool  empty () const { return size_ == 0; } ///< true, iff the area of the intersection is 0.
    size_t size () const { return size_; }      ///< Number of triangles, 0, 1, or 2

    ///@{ Random-access to the triangles
    const_triangle_iterator triangle_begin () const { return triangle_; }
    const_triangle_iterator triangle_end   () const { return triangle_ + size_; }
    ///@}
};

///\brief Return a signed array-index for the possible 3^4 sign-patterns on the vertices of a tetra.
/// The index ranges from [-40..40].
inline byte instance_idx (const byte ls[4])
{
    return  27*ls[0] + 9*ls[1] + 3*ls[2] + ls[3];
}

///\brief Return a signed array-index for the possible 3^4 sign-patterns on the vertices of a tetra.
/// The index ranges from [-40..40].
inline byte instance_idx (const double ls[4])
{
    return  27*sign( ls[0]) + 9*sign( ls[1]) + 3*sign( ls[2]) + sign( ls[3]);
}

inline const RefTetraPatchCL&
RefTetraPatchCL::instance (const byte ls[4])
{
    const byte idx= instance_idx ( ls);
    if (idx == 0)
        throw DROPSErrCL( "RefTetraPatchCL::instance: found 3-dim. zero level set, grid is too coarse!");
    return instance_array_[idx + 40];
}

inline const RefTetraPatchCL&
RefTetraPatchCL::instance (const double ls[4])
{
    byte ls_byte[4];
    std::transform( ls + 0, ls + 4, ls_byte + 0, DROPS::sign);
    return instance( ls_byte);
}

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
    const byte idx= instance_idx ( ls);
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

/// \brief creat a partition of a boundary triangle face cut by interface
//Setting level-set value of the opposite vertex of the face to 0,  with the rule of RefTetraPartitionCL we can first make a partition of the tetra which is the only tetra has this face;
//Then we just map this tetra partition to the face we are interested in.   
class RefTrianglePartitionCL
{
  public:
    typedef SArrayCL<Ubyte, 3> TriangleT;            ///< the vertices of a triangle of the cut: the tetra's vertices are denoted by 0..3, the edge-cuts by edge-num + 4, which is in 4..9.
    typedef const TriangleT* const_triangle_iterator;
    typedef       TriangleT*       triangle_iterator;
	
  private:
    TriangleT triangle_[3];      ///< at most three triangles
    Ubyte size_;                 ///< number of triangles
	int sign_[3];                ///< sign  the triangles

    TriangleT MakeTriangle (Ubyte v0, Ubyte v1, Ubyte v2) const { return MakeSArray( v0, v1, v2); }
	
	public:
    //parameters: level set sign of the tetra which has the face with ls[VertexNum] = 0;
    RefTrianglePartitionCL(const byte ls[4], Ubyte VertexNum); 
    size_t size () const { return size_; }      ///< Number of triangles, 0, 1, or 2
    int sign (const_triangle_iterator t) const { return sign_[ t- triangle_begin() ]; } ///< Sign of the triangle, to which t points

    ///@{ Random-access to the triangles
    const_triangle_iterator triangle_begin () const { return triangle_; }
    const_triangle_iterator triangle_end   () const { return triangle_ + size_; }
    ///@}
};

} // end of namespace DROPS

#endif
