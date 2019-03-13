/// \file subtriangulation.tpp
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
 * Copyright 2016 LNM/SC RWTH Aachen, Germany
*/

#include "geom/topo.h"

#include <iterator>

namespace DROPS
{

template <class GridFunT>
  inline bool
  equal_signs (const GridFunT& f)
{
    int sum= 0;
    for (Uint i= 0; i < f.size(); ++i)
        sum+= sign( f[i]);

    return static_cast<size_t>( std::abs( sum)) == f.size(); // std::abs return a signed type --> silence g++ warning
}

    template <class GridFunT>
    double  distance (const GridFunT& f)
    {
        double min_val= 10000;
        if(equal_signs( f))
        {
            for (Uint i= 0; i < f.size(); ++i)
                min_val= std::min( std::fabs(f[i]),min_val);
        }
        else
            min_val=0;
        return min_val;
    }
///\brief Write the sign of the levelset function src to the sequence dst.
inline void
copy_levelset_sign (const std::valarray<double>& src, std::valarray<byte>& dst)
{
    dst.resize( src.size());
    std::transform( Addr( src), Addr( src) + src.size(), Addr( dst), sign);
}

///\brief Copies the level set values and the corresponding signs from the global vectors to an array specific to the given lattice-tetra
template <Uint Dim>
inline void
copy_local_level_set_values (const std::valarray<double>& ls,
                             const std::valarray<byte>&   ls_sign,
                             const typename DimensionTraitsCL<Dim>::LatticeBodyT& lattice_tet,
                             double loc_ls[Dim + 1],
                             byte   loc_ls_sign[Dim + 1])
{
    for (Uint i= 0; i < Dim + 1; ++i) {
        loc_ls_sign[i]= ls_sign[lattice_tet[i]];
        loc_ls     [i]= ls     [lattice_tet[i]];
    }
}

///\brief Check if a sub-tetra in the lattice has the face with FaceNum, and copy the local level set values and signs;
// This function is created to check if a sub-tetra has a face which is a cut face w.r.t. the interface
inline bool
on_cutBnd(const std::valarray<double>& ls, const std::valarray<byte>& ls_sign, const PrincipalLatticeCL& lattice, const PrincipalLatticeCL::TetraT& lattice_tet, 
double loc_ls[4], byte loc_ls_sign[4], Uint FaceNum)
{
    bool OnCutBnd = false;
    int counter = 0;
    for (Uint i= 0; i < 4; ++i)
    {
        loc_ls_sign[i]= ls_sign[lattice_tet[i]];
        loc_ls     [i]= ls     [lattice_tet[i]];
        //check if the vertex i is on the face 	
        if( lattice.GetBaryCoord(lattice_tet[i])[FaceNum] < 10e-9)
        {
            counter++;
        }
    }
    if (counter == 3 )
        OnCutBnd = true;
    return OnCutBnd; 
}

template <template <Uint Dim> class VertexCutMergingPolicyT>
  const TetraPartitionCL::TetraT
  TetraPartitionCL::make_sub_tetra (
    const RefTetraPartitionCL::TetraT& ref_tet, const PrincipalLatticeCL::TetraT& lattice_tet,
    const double lset[4], Uint lattice_num_vertexes,
    VertexCutMergingPolicyT<3>& edgecut)
{
    Uint loc_vert_num;
    TetraT tet;

    for (Uint j= 0; j < 4; ++j) {
        loc_vert_num= ref_tet[j];
        if (loc_vert_num < 4) // zero vertex
            tet[j]= lattice_tet[loc_vert_num];
        else { // genuine edge cut
            const Ubyte v0= VertOfEdge( loc_vert_num - 4, 0),
                        v1= VertOfEdge( loc_vert_num - 4, 1);
            tet[j]= lattice_num_vertexes + edgecut( lattice_tet[v0], lattice_tet[v1], lset[v0], lset[v1]);
        }
    }
    return tet;
}

template <class VertexPartitionPolicyT,
          template <Uint Dim> class VertexCutMergingPolicyT>
  void
  TetraPartitionCL::make_partition (const PrincipalLatticeCL& lat, const std::valarray<double>& ls)
{
    const Uint lattice_num_vertexes= lat.vertex_size();

    tetras_.clear();
    pos_tetra_begin_= 0;
    vertexes_.clear();

    std::valarray<byte> ls_sign;
    copy_levelset_sign( ls, ls_sign);

    VertexCutMergingPolicyT<3> edgecut( lat.vertex_begin()); // Stores the genuine cuts.

    TetraContT loc_tetras; // temporary container for the positive tetras.
    double loc_ls[4];
    byte   loc_ls_sign[4];
    for (PrincipalLatticeCL::const_tetra_iterator lattice_tet= lat.tetra_begin(), lattice_end= lat.tetra_end(); lattice_tet != lattice_end; ++lattice_tet) {
        copy_local_level_set_values<3>( ls, ls_sign, *lattice_tet, loc_ls, loc_ls_sign);
        const RefTetraPartitionCL& cut= RefTetraPartitionCL::instance( loc_ls_sign);
        for (RefTetraPartitionCL::const_tetra_iterator it= cut.tetra_begin(), end= cut.tetra_end(); it != end; ++it)
            (cut.sign( it) == -1 ? tetras_ : loc_tetras).push_back( make_sub_tetra(
                *it, *lattice_tet, loc_ls, lattice_num_vertexes, edgecut));
    }
    pos_tetra_begin_= tetras_.size();
    std::copy( loc_tetras. begin(), loc_tetras.end(), std::back_inserter( tetras_));

    VertexPartitionPolicyT vertex_order_policy( lat, ls, tetras_.begin(), tetras_.end(), pos_tetra_begin_);
    vertex_order_policy.sort_vertexes( vertexes_, edgecut.cut_vertex_container(), pos_vertex_begin_, neg_vertex_end_);
}

template <Uint Dim>
template <template <Uint> class VertexCutMergingPolicyT>
  const typename SPatchCL<Dim>::FacetT
  SPatchCL<Dim>::make_sub_facet (const RefPatchFacetT& ref_tri, const LatticeBodyT& lattice_tet,
    const LatticeT& lattice, const double lset[Dim+1],
    std::vector<Uint>& copied_vertexes, std::vector<RenumberVertexPairT>& zero_vertex_uses,
    VertexCutMergingPolicyT<Dim>& edgecut)
{
    const Uint NotCopiedC= static_cast<Uint>( -1);
    FacetT tri;
    Uint loc_vert_num;

    for (Uint j= 0; j < Dim; ++j) {
        loc_vert_num= ref_tri[j];
        if (loc_vert_num < Dim + 1) { // zero vertex
            if (copied_vertexes.empty()) // lazy initialization of the lookup-table for copied zero-vertexes
                copied_vertexes.resize( lattice.vertex_size(), NotCopiedC);
            const Uint lattice_vert_num= lattice_tet[loc_vert_num];
            if (copied_vertexes[lattice_vert_num] == NotCopiedC) {
                copied_vertexes[lattice_vert_num]= vertexes_.size();
                vertexes_.push_back( lattice.vertex_begin()[lattice_vert_num]);
            }
            tri[j]= copied_vertexes[lattice_vert_num];
            zero_vertex_uses.push_back( std::make_pair( facets_.size(), j));
        }
        else { // genuine edge vertex
            const Ubyte v0= VertOfEdge<Dim>( loc_vert_num - (Dim + 1), 0),
                        v1= VertOfEdge<Dim>( loc_vert_num - (Dim + 1), 1);
            tri[j]= edgecut( lattice_tet[v0], lattice_tet[v1], lset[v0], lset[v1]);
        }
    }
    return tri;
}

template <Uint Dim>
  void
  SPatchCL<Dim>::clear ()
{
    facets_.clear();
    is_boundary_facet_.clear();
    vertexes_.clear();
    world_vertexes_.clear();
    normals_.clear();
    absdets_.clear();

}

template <Uint Dim>
template <template <Uint> class VertexCutMergingPolicyT>
  void
  SPatchCL<Dim>::make_patch (const LatticeT& lat, const std::valarray<double>& ls)
{
    this->clear();

    std::valarray<byte> ls_sign;
    copy_levelset_sign( ls, ls_sign);

    VertexCutMergingPolicyT<Dim> edgecut( lat.vertex_begin());

    std::vector<Uint> copied_vertexes; // lookup-table for copied zero-vertexes from the lattice
    std::vector<RenumberVertexPairT> zero_vertex_uses; // used to renumber the zero_vertexes

    double loc_ls[Dim + 1];
    byte   loc_ls_sign[Dim + 1];
    for (typename LatticeT::const_body_iterator lattice_tet= lat.body_begin(), lattice_end= lat.body_end();
        lattice_tet != lattice_end; ++lattice_tet) {
        copy_local_level_set_values<Dim>( ls, ls_sign, *lattice_tet, loc_ls, loc_ls_sign);
        const RefPatchT& cut= RefPatchT::instance( loc_ls_sign);
        if (cut.empty()) continue;
        for (typename RefPatchT::const_facet_iterator it= cut.facet_begin(), end= cut.facet_end(); it != end; ++it) {
            facets_.push_back( make_sub_facet(
                *it, *lattice_tet, lat, loc_ls, copied_vertexes, zero_vertex_uses, edgecut));
            is_boundary_facet_.push_back( cut.is_boundary_facet());
        }
    }

    // Renumber the zero vertexes in the triangles, as they will be at offset #edge-cuts in vertexes_.
    const Uint num_genuine_cuts= edgecut.cut_vertex_container().size();
    for (std::vector<RenumberVertexPairT>::const_iterator it= zero_vertex_uses.begin(); it != zero_vertex_uses.end(); ++it)
        facets_[it->first][it->second]+= num_genuine_cuts;
    // Prepend the edge-cuts to vertexes_
    std::copy( vertexes_.begin(), vertexes_.end(), std::back_inserter( edgecut.cut_vertex_container()));
    vertexes_.swap( edgecut.cut_vertex_container());
}


//***************************************************************************************************
//                                        BndTriangPartitionCL
//***************************************************************************************************
template <class VertexCutMergingPolicyT>
  const BndTriangPartitionCL::TriangleT 
  BndTriangPartitionCL::make_sub_triangle (const RefTrianglePartitionCL::TriangleT& ref_tri, const PrincipalLatticeCL::TetraT& lattice_tet,
        const double lset[4], Uint lattice_num_vertexes, VertexCutMergingPolicyT& edgecut)
{
    TriangleT tri;
    Uint loc_vert_num;
    for (Uint j= 0; j < 3; ++j) {
        loc_vert_num= ref_tri[j];
        if (loc_vert_num < 4)   // on the origianl vertices
            tri[j]= lattice_tet[loc_vert_num];
        else {                  // genuine edge cut
            const Ubyte v0= VertOfEdge( loc_vert_num - 4, 0),
                        v1= VertOfEdge( loc_vert_num - 4, 1);
            tri[j]= lattice_num_vertexes + edgecut( lattice_tet[v0], lattice_tet[v1], lset[v0], lset[v1]); // add the cut vertex if not find and return its number 
        }
    }
    return tri;
}

template <class VertexPartitionPolicyT, template <Uint Dim> class VertexCutMergingPolicyT>
  void 
  BndTriangPartitionCL::make_partition2D (const PrincipalLatticeCL& lat, Uint face, const std::valarray<double>& ls)
{
    const Uint lattice_num_vertexes= lat.vertex_size();

    triangles_.clear();
    pos_triangles_begin_ = 0;
    vertexes_.clear();

    std::valarray<byte> ls_sign;
    copy_levelset_sign( ls, ls_sign);

    VertexCutMergingPolicyT<3> edgecut( lat.vertex_begin()); // Stores the genuine cuts.

    TriangleContT  loc_triangles; // temporary container for the positive triangle.
    double loc_ls[4];
    byte   loc_ls_sign[4];
    for (PrincipalLatticeCL::const_tetra_iterator lattice_tet= lat.tetra_begin(), lattice_end= lat.tetra_end(); lattice_tet != lattice_end; ++lattice_tet) {
        
        // check if a sub_tetra in the lattice (with 2 intervals) is a sub_tetra which has the face
        bool OnCutBnd = on_cutBnd( ls, ls_sign, lat, *lattice_tet, loc_ls, loc_ls_sign, face);
        if (OnCutBnd)
        {
            const RefTrianglePartitionCL cut(loc_ls_sign, Ubyte(face));                          //casting from unsigned int to unsigned char?
            for (RefTrianglePartitionCL::const_triangle_iterator it= cut.triangle_begin(), end= cut.triangle_end(); it != end; ++it)
                (cut.sign( it) == -1 ? triangles_ : loc_triangles).push_back( make_sub_triangle(
                    *it, *lattice_tet, loc_ls, lattice_num_vertexes, edgecut));
        }
    }
    pos_triangles_begin_= triangles_.size(); 

    if (loc_triangles.end()-loc_triangles.begin()!=0)
        std::copy( loc_triangles. begin(), loc_triangles.end(), std::back_inserter( triangles_)); //insert the positive triangles

    //unsorted vertexes, other sorted policy could be added if needed in the future
    std::copy( lat.vertex_begin(), lat.vertex_end(), std::back_inserter( vertexes_));             //insert all vertices in the principal lattice;
    VertexContT cut_vertexes=edgecut.cut_vertex_container();
    std::copy( cut_vertexes.begin(), cut_vertexes.end(), std::back_inserter( vertexes_));         //insert all vertices generated by edge cut;
}

} // end of namespace DROPS
