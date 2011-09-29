/// \file renumber.h
/// \brief bandwith reduction and downwind numbering
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

#ifndef DROPS_RENUMBER_H
#define DROPS_RENUMBER_H

#include "num/spmat.h"

namespace DROPS
{

//=============================================================================
//  Reverse Cuthill-McKee ordering
//=============================================================================

/// \brief The first component is the degree of the vertex, the second is its number.
///
/// This is used by reverse_cuthill_mckee and its helpers.
typedef std::pair<size_t, size_t> DegVertT;
const size_t NoVert= std::numeric_limits<size_t>::max();

/// \brief Traverse the graph starting at v in breadth first preorder.
///
/// This is a helper of reverse_cuthill_mckee.  Vertices in each level
/// are numbered by increasing degree.
/// \param v The first vertex.
/// \param M Interpreted as adjacency matrix of the graph; see reverse_cuthill_mckee. Only verticed with p[v] != NoIdx are considered.
/// \param p Numbering of the vertices.
/// \param idx Next number to be used.
/// \param degree The degree of each vertex
template <typename T>
  void
  breadth_first (size_t v, const SparseMatBaseCL<T>& M, PermutationT& p, size_t& idx,
    std::vector<size_t>& degree);

/// \brief Find a vertex with high eccentricity.
///
/// This is a helper of reverse_cuthill_mckee. It computes (heuristically)
/// a start-vertex, the distance of which to at least one other vertex
/// is close to the diameter of of the graph. In a sense, such nodes
/// are on the boundary of the graph.
///
/// After a breadth first traversal a min-degree vertex from the last
/// level is chosen and the procedure is repeated if its eccentricity
/// is larger than that of the previous candidate.
///
/// \param V The vertices of the (sub-) graph together with their degree.
/// \param M Interpreted as adjacency matrix of the graph; see reverse_cuthill_mckee.
/// \param max_iter Maximal number of iteration. Though this heuristic will terminate
///     after at least |V| steps, this should be a small number like 5.
template <typename T>
  size_t
  pseudo_peripheral_node (std::vector<DegVertT>& V, const SparseMatBaseCL<T>& M,
    size_t& max_iter);

/// \brief Compute a permutation such that M has (in most cases) near
/// minimal bandwith.
///
/// M_in is interpreted as adjacency matrix of a graph with vertices i
/// from [0...M.num_rows()).  Edge (i, j) exists, iff M_ij != 0.
///
/// From the linear algebra standpoint that means that x_i depends on
/// x_j, so these edges are actually the incoming edges of x_i. We use
/// this definition because we can traverse the incoming edges (rows)
/// much faster than the outgoing edges (columns). Thus, this routine
/// computes the rcm-permutation of M^T. Note, that the bandwith of A and
/// A^T is identical, so this probably does not matter. If desired the
/// function can compute the transpose of M and find the corresponding
/// rcm-permutation. This needs as much storage as an assembly of M.
///
/// \param M_in Interpreted as adjacency matrix of the graph.
/// \param p Contains for each unknown i its new number p[i].
/// \param use_indegree true: Compute the rcm-permutation of M^T
///     false: Compute the rcm-permutation of M; this is memory intensive.
template <typename T>
void
reverse_cuthill_mckee (const SparseMatBaseCL<T>& M_in, PermutationT& p,
    bool use_indegree= true);


//=============================================================================
//  Downwind numbering
//=============================================================================

/// \brief Produce a downwind numbering for the graph of a convection matrix.
/// The square matrix M is interpreted as graph with vertices 0..M.nom_cols()-1.
/// There is a directed edge (i,j), iff M_ij > 0.
/// The strongly connected components of the graph are computed with Tarjan's algorithm.
/// If there are no cycles in the graph, all components contain exactly one vertex.
/// The algo discovers the components in downwind direction (actually in upwind direction
/// of the graph of M^T).
/// The numbering within a component has no furhter meaning.
///
/// Permuting M with the permutation returned makes the positive part of M block-lower-triangular,
/// where the blocks correspond to the connected components. For an acyclic graph, the matrix is
/// lower triangular.
/// This makes one Gauss-Seidel sweep an excellent preconditioner.
class TarjanDownwindCL
{
  private:
    PermutationT        new_index;       ///< new_index[i] is the new number of vertex i.
    std::vector<size_t> t_index;         ///< index in Tarjan's algorithm
    std::vector<size_t> low_link;        ///< the smallest index encountered in the current component
    std::vector<size_t> component_;      ///< component[i] is the number of the component of i.
    std::vector<size_t> component_size_; ///< number of elements in component i.

    size_t idx;                ///< Counter for next available vertex index
    size_t tidx;               ///< Counter for the next available index in Tarjan's algorithm
    std::vector<size_t> stack; ///< Stack, which holds unfinished, preceding connected components.
    std::vector<bool> is_on_stack;

    /// \brief Perform a depth-first search to discover the connected component of v.
    template < class T>
      void number_connected_component_of (size_t v, const SparseMatBaseCL<T>& M);

  public:
    /// \brief Compute a downwind numbering for the convection-matrix M.
    template <class T>
      const PermutationT& number_connected_components (const SparseMatBaseCL<T>& M);

    /// \brief Get the computed permutation
    const PermutationT&        permutation    () const { return new_index; }
    /// \brief Number of strongly connected components.
    size_t                     num_components () const { return component_size_.size(); }
    /// \brief component_map()[i] is the number of the component, to which vertex i belongs.
    const std::vector<size_t>& component_map  () const { return component_; }
    /// \brief component(c) contains all vertices in component c. The vertices are in ascending order.
    std::vector<size_t> component  (size_t c) const;
    /// \brief the number of vertices in component i.
    const std::vector<size_t>& component_size () const { return component_size_; }

    /// \brief High-level overview of the numbering proccess.
    inline void stats (std::ostream& os) const;
};

} // end of namspace DROPS

#include "num/renumber.tpp"

#endif
