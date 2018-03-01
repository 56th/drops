/// \file signtraits.h
/// \brief Properties of the reference-d-simplex, which is cut by a linear level set
/// \author LNM RWTH Aachen: Joerg Grande

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

#ifndef DROPS_SIGNTRAITS_H
#define DROPS_SIGNTRAITS_H

#include "geom/topo.h"

namespace DROPS
{

template <Uint Dim= 3>
class SignTraitsCL;

template <Uint Dim>
std::ostream&
operator<< (std::ostream& out, const SignTraitsCL<Dim>& c);

///\brief Represents the reference simplex of dimension Dim, which is cut by a linear level set function ls. The values of the latter are prescribed on the vertices.
template <Uint Dim>
class SignTraitsCL
{
  public:
    enum { NumVerts= Dim + 1,
           NumEdges= (Dim*(Dim + 1))/2,
// dimension vs. max. number of roots:
// 1 2, the 0-pattern
// 2 3, the 0-pattern
// 3 4, 2 vertexes on each side (and 0-pattern)
// 4 6, 2 verts on one side, three on the other: 10 - 1 - 3
// in general:
// d < 3: see above,
// else: distribute vertexes as evenly as possible on both sides:
//       d odd:  (d+1 over 2) - 2*((d+1)/2 over 2)              = ((d+1)/2)^2
//       d even: (d+1 over 2) - (d/2 over 2) - (d/2 + 1 over 2) = d*(d+2)/4
           max_num_root= Dim < 4 ? Dim + 1 : (Dim & 1u ? NumVerts*NumVerts : Dim*(Dim + 2))/4,

           num_pattern= Dim== 1 ? 9 : Dim == 2 ? 27 : Dim == 3 ? 81 : Dim == 4 ? 243 : 0, // 3^(Dim+1) sign patterns.
           zero_pattern_offset= (num_pattern - 1)/2
    };

    /// \brief enumerate all sign patterns: returns an index in [-(num_pattern-1)/2..num_pattern-1)/2]. The zero pattern is at index 0.
    /// Dimensions > 4 are currently not supported (byte --> short)
    static byte pattern_idx (const byte ls[Dim + 1]) {
        return  Dim == 1 ? ls[0] + 3*ls[1] :
                Dim == 2 ? ls[0] + 3*ls[1] + 9*ls[2] :
                Dim == 3 ? ls[0] + 3*ls[1] + 9*ls[2] + 27*ls[3] :
                Dim == 4 ? ls[0] + 3*ls[1] + 9*ls[2] + 27*ls[3] + 81*ls[4]:
                0;
    }

  private:
    Ubyte num_root_vert_; ///< number of vertices, where the level set function is zero.
    Ubyte num_root_;      ///< number of roots of the level set function; invariant: num_root_vert <= num_root <= max_num_root.
    byte sign_[NumVerts]; ///< Sign of the level set function of the vertices; \f$\in\{-1,0,1\}\f$

    ///\brief local number with respect to the reference simplex of the object on the cut: [0,num_root_vert): vertex numbers in (0..NumVerts-1); [num_root_vert, num_root): edge numbers in (0..NumEdges-1). Both parts are sorted in increasing order.
    Ubyte cut_simplex_    [max_num_root];
    Ubyte cut_simplex_rep_[max_num_root]; ///< local number of the object on the cut: (0..NumVerts+NumEdges-1)

    void compute_cuts ();

  public:
    SignTraitsCL () : num_root_vert_( NumVerts), num_root_( NumVerts) {} ///< Uninitialized default state
    SignTraitsCL (const byte   ls[NumVerts]) { assign( ls); }            ///< Assign the sign pattern on the vertices.
    SignTraitsCL (const double ls[NumVerts]) { assign( ls); }            ///< Assign the sign pattern on the vertices.
    void assign (const byte   ls[NumVerts]); ///< Assign a sign pattern on the vertices.
    void assign (const double ls[NumVerts]); ///< Assign a sign pattern on the vertices.

    bool empty () const { return num_root_ == 0; } ///< True, iff there is no intersection.

    bool has_codim_le_1 () const { return num_root_ >= Dim; }        ///< True, iff the intersection has non-zero (Dim-1)-dimensional measure.
    bool has_codim_0 () const { return num_root_vert_ == NumVerts; } ///< True, exactly for the zero-pattern.
    bool no_zero_vertex () const { return num_root_vert_ == 0; }     ///< True, iff there is no vertex, in which ls vanishes.

    Ubyte num_cut_simplexes () const { return num_root_; }      ///< Number of edges and vertices with a root of ls.
    Ubyte num_zero_vertexes () const { return num_root_vert_; } ///< Number of vertices of the simplex that are roots of ls.

    byte sign (int i) const { return sign_[i]; } ///< -1,0,1; sign of vertex i.

    /// Return local number of edges/verts with a root of ls. For edges, [] returns a edge number in 0..NumEdges-1, and () returns an extended vertex number in NumVerts..NumVerts+NumEdges-1.
    ///@{
    Ubyte operator[] (int i) const { return cut_simplex_[i]; }
    Ubyte operator() (int i) const { return cut_simplex_rep_[i]; }
    ///@}
    friend std::ostream& operator<< <Dim>(std::ostream&, const SignTraitsCL<Dim>&); ///< Debug-output to a stream (dumps all members)
};

template <Uint Dim>
void
SignTraitsCL<Dim>::compute_cuts ()
{
    for (Ubyte i= 0; i < NumVerts; ++i)
        if (sign( i) == 0)
            cut_simplex_[num_root_vert_++]= i;
    num_root_= num_root_vert_;
    for (Ubyte i= 0; i < NumEdges; ++i)
        if (sign( VertOfEdge<Dim>( i, 0))*sign( VertOfEdge<Dim>( i, 1)) == -1)
            cut_simplex_[num_root_++]= i;
    std::memcpy( cut_simplex_rep_, cut_simplex_, num_root_*sizeof(byte));
    for (int i= num_root_vert_; i < num_root_; ++i)
        cut_simplex_rep_[i]+= NumVerts;
}

template <Uint Dim>
void
SignTraitsCL<Dim>::assign (const byte ls[NumVerts])
{
    num_root_vert_= num_root_= 0;

    byte sum= 0;
    for (Ubyte i= 0; i < NumVerts; ++i)
        sum+= (sign_[i]= ls[i]);
    if (std::abs( sum) == static_cast<int>( NumVerts)) // optimize the case of uncut tetras
        return;

    compute_cuts ();
}

template <Uint Dim>
inline void
SignTraitsCL<Dim>::assign (const double ls[NumVerts])
{
    byte ls_byte[NumVerts];
    std::transform( ls + 0, ls + NumVerts, ls_byte + 0, DROPS::sign);
    this->assign( ls_byte);

}

template <Uint Dim>
std::ostream&
operator<< (std::ostream& out, const SignTraitsCL<Dim>& c)
{
    const int NumVerts= Dim + 1;
    out << static_cast<int>( c.num_root_vert_) << ' ' << static_cast<int>( c.num_root_) << '\n';
    for (int i= 0; i < NumVerts; ++i)
        out << static_cast<int>( c.sign_[i]) << ' ';
    out << '\n';
    for (int i= 0; i < NumVerts; ++i)
        out << static_cast<int>( c.cut_simplex_[i]) << ' ';
    out << '\n';
    for (int i= 0; i < NumVerts; ++i)
        out << static_cast<int>( c.cut_simplex_rep_[i]) << ' ';
    return out << '\n';
}

} // end of namespace DROPS
#endif