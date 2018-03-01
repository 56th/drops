/// \file reftetracut.cpp
/// \brief Triangulation of the reference tetraeder adapted to a linear level-set function.
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

#include "geom/reftetracut.h"
#include "geom/topo.h"

#include <iostream>
#include <cstring>

namespace DROPS {

inline Ubyte
num_triangles (const SignTraitsCL<3>& cut)
{
    return cut.has_codim_le_1() ? cut.num_cut_simplexes() - 2 : 0;
}


void RefPatchBuilderCL<3>::StaticInit (RefPatchCL<3> instance_array[SignTraitsCL<3>::num_pattern])
{
    byte ls[4];
    for (ls[0]= -1; ls[0] < 2; ++ls[0])
        for (ls[1]= -1; ls[1] < 2; ++ls[1])
            for (ls[2]= -1; ls[2] < 2; ++ls[2])
                for (ls[3]= -1; ls[3] < 2; ++ls[3]) {
                    if ( ls[0] == 0 && ls[1] == 0 && ls[2] == 0 && ls[3] == 0)
                        continue;
                    instance_array[SignTraitsCL<3>::pattern_idx( ls) + SignTraitsCL<3>::zero_pattern_offset].assign( ls);
                }
}

bool
RefPatchBuilderCL<3>::assign (const SignTraitsCL<3>& cut, RefPatchCL<3>& p)
{
    for (p.size_= 0; p.size_ < num_triangles( cut); ++p.size_)
        p.facet_[p.size_]= RefPatchBuilderCL<3>::MakeTriangle( cut(p.size_), cut(p.size_ + 1), cut(p.size_ + 2));
    p.is_boundary_facet_= cut.num_zero_vertexes() == 3 ? 1 : 0;
    return p.empty();
}

template class RefPatchCL<3>;

namespace {
StaticInitializerCL<RefPatchCL<3> > RefPatch3_initializer_;

} // end of anonymous namespace


void RefPatchBuilderCL<4>::StaticInit (RefPatchCL<4> instance_array[SignTraitsCL<4>::num_pattern])
{
    byte ls[5];
    for (ls[0]= -1; ls[0] < 2; ++ls[0])
        for (ls[1]= -1; ls[1] < 2; ++ls[1])
            for (ls[2]= -1; ls[2] < 2; ++ls[2])
                for (ls[3]= -1; ls[3] < 2; ++ls[3])
                    for (ls[4]= -1; ls[4] < 2; ++ls[4]) {
                        if ( ls[0] == 0 && ls[1] == 0 && ls[2] == 0 && ls[3] == 0 && ls[4] == 0)
                            continue;
                        instance_array[SignTraitsCL<4>::pattern_idx( ls) + SignTraitsCL<4>::zero_pattern_offset].assign( ls);
                    }
}

inline bool
RefPatchBuilderCL<4>::cone_construction (const SignTraitsCL<4>& cut, RefPatchCL<4>& p, Ubyte f)
{
    const Ubyte v= cut( 0); // The tip of the cones.
    byte f_ls[4]; // level set signs on the facet f of the penta.
    for (Uint i= 0; i < 4; ++i)
        f_ls[i]= cut.sign( RefPenta::VertOfTetra( f, i));
//     const RefTetraPatchCL& pp= RefTetraPatchCL::instance( f_ls);
//     for (RefTetraPatchCL::const_facet_iterator it= pp.facet_begin(), end= pp.facet_end(); it != end; ++it) {
//         tetra_[size_++]= MakeTetra( v,
//                                     PentaCutByTetraCut( f, it[0][0]),
//                                     PentaCutByTetraCut( f, it[0][1])),
//                                     PentaCutByTetraCut( f, it[0][2]));
//     return p.is_boundary_facet();
    SignPatternTraitCL f_cut( f_ls);
    for (Uint i= 0; i < num_triangles( f_cut); ++i)
        p.facet_[p.size_++]= MakeTetra( v,
                                        PentaCutByTetraCut( f, f_cut( i)),
                                        PentaCutByTetraCut( f, f_cut( i + 1)),
                                        PentaCutByTetraCut( f, f_cut( i + 2)));
    return f_cut.num_zero_vertexes() == 3 ? 1 : 0;
}

bool
RefPatchBuilderCL<4>::assign (const SignTraitsCL<4>& cut, RefPatchCL<4>& p)
{
    p.size_= 0;
    p.is_boundary_facet_= 0;

    if (cut.empty() || !cut.has_codim_le_1())
        return p.empty();
    // From here on, there are at least 4 roots of the level set function.

    if (cut.no_zero_vertex()) { //cut(0) is an edge, there are two facets f0, f1 of the penta not containing it.
        // is_on_boundary_ is false (otherwise, at least four zero-vertexes are needed).
        for (Uint v= 0; v < 2; ++v)
            cone_construction( cut, p, RefPenta::OppTetra( RefPenta::VertOfEdge( cut[0], v)));
    }
    else //cut(0) is a vertex, there is one facet f of the penta not containing it.
        p.is_boundary_facet_= cone_construction( cut, p, RefPenta::OppTetra( cut[0]));

    return p.empty();
}

template class RefPatchCL<4>;

namespace {
StaticInitializerCL<RefPatchCL<4> > RefPatch4_initializer_;

} // end of anonymous namespace


RefTetraPartitionCL RefTetraPartitionCL::instance_array_[SignTraitsCL<3>::num_pattern];

void RefTetraPartitionCL::StaticInit ()
{
    byte ls[4];
    for (ls[0]= -1; ls[0] < 2; ++ls[0])
        for (ls[1]= -1; ls[1] < 2; ++ls[1])
            for (ls[2]= -1; ls[2] < 2; ++ls[2])
                for (ls[3]= -1; ls[3] < 2; ++ls[3]) {
                    if ( ls[0] == 0 && ls[1] == 0 && ls[2] == 0 && ls[3] == 0)
                        continue;
                    RefTetraPartitionCL::instance_array_[SignTraitsCL<3>::pattern_idx( ls) + SignTraitsCL<3>::zero_pattern_offset].assign( ls);
                }
}

namespace {
StaticInitializerCL<RefTetraPartitionCL> RefTetraPartition_initializer_;

} // end of anonymous namespace

Ubyte
RefTetraPartitionCL::some_non_zero_vertex (const SignPatternTraitCL& cut) const
{
    Ubyte v;
    for (v= 0; v < cut.num_zero_vertexes() && cut[v] == v; ++v)
        /*empty body*/;
    return v;
}

bool
RefTetraPartitionCL::assign (const SignPatternTraitCL& cut)
{
    end_= begin_= tetras_ + 3;

    if (cut.empty()) { // Most common case: no cut.
        AddTetra( 0, 1, 2, 3, cut.sign( 0));
    }
    else if (cut.no_zero_vertex()) { // next common case: cuts without vertices on the zero level
        if (cut.num_cut_simplexes() == 3) { // triangular cut: a tetra and a remaining prism
            const Ubyte v= VertByEdge( cut[0], cut[1]);
            AddTetra( v, cut(0), cut(1), cut(2), cut.sign( v));
            AddPrism( OppVertOfEdge( cut[0], v), cut(0),
                      OppVertOfEdge( cut[1], v), cut(1),
                      OppVertOfEdge( cut[2], v), cut(2),
                      -cut.sign( v));
        }
        else if (cut.num_cut_simplexes() == 4) { // quadrilateral cut: two prisms
            const Ubyte e= first_uncut_edge( cut);
            const Ubyte f= OppEdge( e);
            AddPrism( VertOfEdge(e, 0), VertOfEdge( e, 1),
                      cut(0), cut(1),
                      cut(2), cut(3),
                      cut.sign(VertOfEdge( e, 0)));
            AddPrism( VertOfEdge(f, 0), VertOfEdge( f, 1),
                      cut(0), cut(2),
                      cut(1), cut(3),
                      cut.sign(VertOfEdge( f, 0)));
        }
    }
    else if (cut.num_cut_simplexes() > cut.num_zero_vertexes()) { // next common case: there are cut edges, and also 1 or 2 vertices of the tetra with value 0 (the latter as we are in the else-part of cut.no_zero_vertex())
        if (cut.num_zero_vertexes() == 1) { // triangular cut through a vertex: a tetra and a remaining pyramid with quadrilateral base
            const Ubyte e= cut[1], f= cut[2];
            const Ubyte v= VertByEdge( e, f);
            AddTetra( v, cut(0), cut(1), cut(2), cut.sign( v));
            const Ubyte opp_v_in_e= v == VertOfEdge( e, 0) ? VertOfEdge( e, 1) : VertOfEdge( e, 0);
            const Ubyte opp_v_in_f= v == VertOfEdge( f, 0) ? VertOfEdge( f, 1) : VertOfEdge( f, 0);
            // the pyramid
            AddTetra( cut(0), cut(1), opp_v_in_f, opp_v_in_e, -cut.sign( v));
            AddTetra( cut(0), cut(1), opp_v_in_f, cut(2), -cut.sign( v));
        }
        else if (cut.num_zero_vertexes() == 2) { // triangular cut through 2 vertexes: two tetras
            const Ubyte e= OppEdge( EdgeByVert( cut[0], cut[1]));
            const Ubyte v0= VertOfEdge( e, 0), v1= VertOfEdge( e, 1);
            AddTetra( cut(0), cut(1), v0, cut(2), cut.sign( v0));
            AddTetra( cut(0), cut(1), v1, cut(2), cut.sign( v1));
        }
    }
    else // remaining cases: 1, 2 or 3 cuts, which are vertices of the tetra
        AddTetra( 0, 1, 2, 3, cut.sign( some_non_zero_vertex( cut)));
    return is_uncut();
}

std::ostream&
operator<< (std::ostream& out, const RefTetraPartitionCL& c)
{
    out << c.end_ - c.begin_ << ' ' << c.tetras_ - c.begin_ << '\n';
    for (Uint i= 0; i < c.end_ - c.begin_; ++i) {
        for (Uint j= 0; j < 4; ++j)
            out << static_cast<Uint>( c.begin_[i][j]) << ' ';
        out << '\n';
    }
    return out;
}

RefTrianglePartitionCL::RefTrianglePartitionCL(byte ls[4], Ubyte VertexNum)
{
    ls[VertexNum] = 0;
    Ubyte vert[3]; //Used to store indices of vertices of the triangle
    const RefTetraPartitionCL& RefTetraPart= RefTetraPartitionCL::instance(ls);
    size_=RefTetraPart.tetra_size(AllTetraC);
    if(size_> 3)
        throw DROPSErrCL("Size of the triangles of RefTrianglePartitionCL is bigger than 3");
    Uint index=0;
    for (RefTetraPartitionCL::const_tetra_iterator it= RefTetraPart.tetra_begin(), end=RefTetraPart.tetra_end(); it != end; ++it) {
        Uint j=0;
        for(Uint i=0; i<4; ++i) {
            if((*it)[i]!=VertexNum) {
                vert[j]=(*it)[i];
                ++j;
                if (j > 3)
                    throw DROPSErrCL("Index of vertex of the triangle is bigger than 3");
            }	 
        }
        //make a triangle and store its sign;
        triangle_[index] = MakeTriangle(vert[0], vert[1], vert[2]);
        sign_[index] = RefTetraPart.sign(it);
        ++index;
    }
}

} // end of namespace DROPS
