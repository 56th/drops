/// \file principallattice.h
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
 * Copyright 2011 LNM/SC RWTH Aachen, Germany
*/

#ifndef DROPS_PRINCIPALLATTICE_H
#define DROPS_PRINCIPALLATTICE_H

#include "misc/container.h"
#include <vector>

namespace DROPS {

/// \brief The principal lattice on the reference tetra.
///
/// The barycentric coordinates of the vertices are computed. The tetras are computed as 4-tuples to the container of vertices.
/// Access to the concrete objects via instance( n), where n is the number of intervals used to subdivide the edges of the reference tetra.
class PrincipalLatticeCL
{
  public:
    typedef SArrayCL<Uint, 4> TetraT; ///< Represents a tetra as indices in the vertex-container.

    typedef std::vector<BaryCoordCL> VertexContT;
    typedef VertexContT::const_iterator const_vertex_iterator;
    typedef std::vector<TetraT> TetraContT;
    typedef TetraContT::const_iterator const_tetra_iterator;

  private:
    /// \brief cache for computed lattices
    ///@{
    class PrincipalLatticeCacheCL : public std::vector<const PrincipalLatticeCL*>
    {
      public:
    	~PrincipalLatticeCacheCL() { for ( size_t i = 0; i< this->size(); ++i) if ((*this)[i]) delete (*this)[i];; }
    };

    static PrincipalLatticeCacheCL cache_;
    static bool is_memoized (Uint n)
        { return n-1 < cache_.size() && cache_[n-1]; }
    static const PrincipalLatticeCL* read_cache (Uint n) { return cache_[n-1]; }
    static const PrincipalLatticeCL* memoize (Uint n);
    ///@}

    ///\brief computation of the tetras
    ///@{
    static TetraContT  tetra_;  ///< All tetras of the triangulation as tuple of indices to vertex_; the list of tetras of a finer subdivision has the list of tetras of a coarser subdivision as common initial subsequence. Therefore tetra_ is static.

    static inline bool in_ref_tetra (const SArrayCL<Uint, 3>& v, Uint num_interval)
        { return num_interval >= v[0] && v[0] >= v[1] && v[1] >= v[2]; }
    static inline Uint vertex_index (Uint i, Uint j, Uint k) { ///< Enumerates the integer-lattice points in the cone generated by the first tetra of the Kuhn-triangulation
        return (i > 0 ? ((i+2)*(i+1)*i)/6 : 0) + (j > 0 ? ((j+1)*j)/2 : 0) + k; // (i>0 ? \binom(i-1 + 3, 3) : 0) + (j>0 ? \binom(j-1 +2, 2) : 0) + k
    }
    ///\brief Creates the tetras in the slice [xbegin..xend] of the first tetra of the Kuhn-triangulation
    static void create_tetras (Uint xbegin, Uint xend);
    ///@}

    Uint n_; ///< number of intervals for the edges
    VertexContT vertex_; ///< All vertices of the lattice as barycentric coordinates

    PrincipalLatticeCL (Uint n);

  public:
    ///\brief number of intervals on each edge
    Uint num_intervals () const { return n_; }
    ///\brief number of vertices in the lattice
    Uint vertex_size () const {
        const Uint nv= num_intervals() + 1;
        return (nv*(nv + 1)*(nv + 2))/6; // \binom(num_intervals() + 3, 3)
    }
    ///\brief number of tetras in the triangulation
    Uint tetra_size  () const { return n_*n_*n_; }

    ///\brief Get baryCoord of a vertex with index number
	BaryCoordCL GetBaryCoord(Uint vertex_num) const { return vertex_[vertex_num];}
    ///\brief Access to vertexes and tetras as sequences (random access iterators)
    ///@{
    const_vertex_iterator vertex_begin ()  const { return vertex_.begin(); }
    const_vertex_iterator vertex_end   ()  const { return vertex_.end(); }
    const_tetra_iterator tetra_begin ()  const { return tetra_.begin(); }
    const_tetra_iterator tetra_end   ()  const { return tetra_.begin() + tetra_size(); }
    ///@}

    ///\brief Access the principal lattice with n intervals on each edge (singleton pattern)
    static const PrincipalLatticeCL& instance (Uint n);
};

extern const size_t p1_dof_on_lattice_2[4];  ///< For vertex i (in 0..3) as counted in topo.h, p1_dof_on_lattice_2[i] is the number of the vertex in the principal lattice of order 2.
extern const size_t p2_dof_on_lattice_2[10]; ///< For a P2-dof i (numbered from 0..9: vertexes, then edges) as counted in topo.h, p2_dof_on_lattice_2[i] is the number of the vertex in the principal lattice of order 2.

} // end of namespace DROPS

#endif