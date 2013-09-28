/// \file principallattice.cpp
/// \brief The principal lattice on the reference tetra obtained by slicing with planes parallel to the faces.
/// \author LNM RWTH Aachen: Joerg Grande; SC RWTH Aachen:

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
 * Copyright 2011,2012 LNM/SC RWTH Aachen, Germany
*/

#include "geom/principallattice.h"

namespace DROPS {

PrincipalLatticeCL::PrincipalLatticeCacheCL PrincipalLatticeCL::cache_;

PrincipalLatticeCL::TetraContT PrincipalLatticeCL::tetra_;

PrincipalLatticeCL::PrincipalLatticeCL (Uint n)
    : n_( n)
{
    vertex_.reserve( vertex_size());
    // Compute vertices of the first tetra in the Kuhn-triangulation of the reference cube in rlex-order of their cartesian coordinates (x,y,z).
    for (Uint x= 0; x <= num_intervals(); ++x)
        for (Uint y= 0; y <= x; ++y)
            for (Uint z= 0; z <= y; ++z) {
                vertex_.push_back( MakeBaryCoord( num_intervals() - x, x - y, y - z, z)); // (x,y,z) in barycentric coordinates
                vertex_.back()/= num_intervals();
            }
}

void PrincipalLatticeCL::create_tetras (Uint xbegin, Uint xend)
{
    // All 6 permutations of (0,1,2)
    const Uint perm[][3]= { {0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0} };
    // compute tetras
    TetraT t( Uninitialized);
    SArrayCL< SArrayCL<Uint,3>, 4> tet;
    // For all candidate starting points:
    tetra_.reserve( xend*xend*xend);
    for (Uint i= xbegin; i < xend; ++i)
        for (Uint j= 0; j <= i; ++j)
            for (Uint k= 0; k <= j; ++k) {
                tet[0][0]= i; tet[0][1]= j; tet[0][2]= k;
                // For all 6 permutations
                for (Uint l= 0; l < 6; ++l) {
                    const Uint* const p= perm[l];
                    // construct tetra and check whether it is inside the reference tetra
                    tet[1]= tet[0];
                    ++tet[1][p[0]];
                    if (!in_ref_tetra( tet[1], xend)) continue;
                    tet[2]= tet[1];
                    ++tet[2][p[1]];
                    if (!in_ref_tetra( tet[2], xend)) continue;
                    tet[3]= tet[2];
                    ++tet[3][p[2]];
                    if (!in_ref_tetra( tet[3], xend)) continue;
                    t[0]= vertex_index( tet[0][0], tet[0][1], tet[0][2]);
                    t[1]= vertex_index( tet[1][0], tet[1][1], tet[1][2]);
                    t[2]= vertex_index( tet[2][0], tet[2][1], tet[2][2]);
                    t[3]= vertex_index( tet[3][0], tet[3][1], tet[3][2]);
                    tetra_.push_back( t);
                }
            }
}

const PrincipalLatticeCL* PrincipalLatticeCL::memoize (Uint n)
{
    if (n-1 >= cache_.size()) {
        const Uint oldsize= cache_.size();
        cache_.resize( n);
        create_tetras( oldsize, n);
    }
    return cache_[n-1]= new PrincipalLatticeCL( n);
}

const PrincipalLatticeCL& PrincipalLatticeCL::instance (Uint n)
{
    const PrincipalLatticeCL* tmp;
#   pragma omp critical
    tmp= is_memoized( n) ? read_cache( n) : memoize( n);
    return *tmp;
}

const size_t p1_dof_on_lattice_2[4]= { 0, 4, 7, 9 };
const size_t p2_dof_on_lattice_2[10]= {
    0, 4, 7, 9,      // vertexes
    1, 2, 5, 3, 6, 8 // edges
};


TetraPrismLatticeCL::TetraPrismLatticeCacheCL TetraPrismLatticeCL::cache_;

TetraPrismLatticeCL::TetraPrismLatticeCL (Uint n, Uint m)
    : m_( m), pl_( PrincipalLatticeCL::instance( n))
{
    vertex_.reserve( vertex_size());
    // Compute vertices
    for (Uint t= 0; t <= num_time_intervals(); ++t)
        for (PrincipalLatticeCL::const_vertex_iterator x= pl_.vertex_begin(); x != pl_.vertex_end(); ++x)
            vertex_.push_back( MakeSTCoord( *x, static_cast<double>( t)/num_time_intervals()));

    penta_.reserve( penta_size());
    for (Uint t= 0; t < num_time_intervals(); ++t)
        for (PrincipalLatticeCL::const_tetra_iterator tet= pl_.tetra_begin(); tet != pl_.tetra_end(); ++tet)
            AddPrism( *tet, t);
}

void TetraPrismLatticeCL::AddPrism (const TetraT& tet, Uint t)
{
/// This is the iterated one construction.
/// First, v=(tet0, t2) is chosen as tip. There are two facets of the tetrahedral prism not containing v: f0 is the tetra (tet0..tet3) x t; f1 is the triangular prism tp=(tet1, tet2, tet3)x(t, t2).
/// tp is tetrahedralized with the cone construction with tip v2=(tet1, t2). The facets of tp not containing v2 are the triangle ff=(tet1..tet3) x t and the quadrilateral q=(tet2, tet3) x (t, t2), which is triangulated by connecting (tet2, t2) with (tet3,t).
/// v           v2        (tet2,t2) (tet3,t2)
/// |           |         |         |
/// (tet0,  t)  (tet1,t)  (tet2,t)  (tet3,t)

    const Uint t2= t + 1;
    penta_.push_back( // v, f0
        MakeSArray( vertex_index( tet[0], t),
                    vertex_index( tet[1], t),
                    vertex_index( tet[2], t),
                    vertex_index( tet[3], t),
                    vertex_index( tet[0], t2)));
    penta_.push_back( // v, v2, ff
        MakeSArray( vertex_index( tet[1], t),
                    vertex_index( tet[2], t),
                    vertex_index( tet[3], t),
                    vertex_index( tet[0], t2),
                    vertex_index( tet[1], t2)));
    penta_.push_back( // v, v2, (tet2,t), (tet3,t), (tet2,t2)
        MakeSArray( vertex_index( tet[2], t),
                    vertex_index( tet[3], t),
                    vertex_index( tet[0], t2),
                    vertex_index( tet[1], t2),
                    vertex_index( tet[2], t2)));
    penta_.push_back( // v, v2, (tet3,t), (tet2,t2), (tet3,t2)
        MakeSArray( vertex_index( tet[3], t),
                    vertex_index( tet[0], t2),
                    vertex_index( tet[1], t2),
                    vertex_index( tet[2], t2),
                    vertex_index( tet[3], t2)));
///\todo The above sequence of indices allows a more compact storage scheme, as the penta-indices overlap perfectly: instead of storing 4*5 indices we get away with the 8-tuple (tet0..tet3) x t, (tet0..tet3) x t2.
}

const TetraPrismLatticeCL* TetraPrismLatticeCL::memoize (Uint n, Uint m)
{
    if (n-1 >= cache_.size())
        cache_.resize( n);
    std::vector<const TetraPrismLatticeCL*>& v( cache_[n-1]);
    if (m-1 >= v.size())
        v.resize( m);
    return v[m-1]= new TetraPrismLatticeCL( n, m);
}

const TetraPrismLatticeCL& TetraPrismLatticeCL::instance (Uint n, Uint m)
{
    const TetraPrismLatticeCL* tmp;
//#   pragma omp critical
    tmp= is_memoized( n, m) ? read_cache( n, m) : memoize( n, m);
    return *tmp;
}

} // end of namespace DROPS
