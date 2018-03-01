/// \file topo.h
/// \brief Topology of the reference tetrahedron and the refinement rules
/// \author LNM RWTH Aachen: Sven Gross, Joerg Peters, Volker Reichelt; SC RWTH Aachen:

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

#ifndef DROPS_TOPO_H
#define DROPS_TOPO_H

#include "misc/utils.h"

namespace DROPS
{

/// Connectivity of the reference tetrahedron and its refinements.

// Constants
const Uint NumVertsC      = 4;
const Uint NumFacesC      = 4;
const Uint NumEdgesC      = 6;
const Uint MaxChildrenC   = 8;
const Uint NumAllVertsC   = 10;
const Uint NumAllEdgesC   = 45;
//const Uint MaxInnerEdgeC  = 13;
const Uint MaxEdgesC      = 25;
const Uint NumObviousEdgesC= 18;
const Uint NumAllFacesC   = 102;
//const Uint MaxInnerFacesC = 8;
//const Uint MaxOuterFacesC = 16;
const Uint MaxFacesC      = 24;
const Uint NumAllChildrenC= 94;
const Uint NumAllSubFacesPerFaceC= 13;

const Uint NoRefMarkC      = 0;
const Uint RegRefMarkC     = 63;
const Uint GreenRegRefMarkC= 127; ///< RegRefMarkC plus an extra bit
const Uint RemoveMarkC     = 64;

const Uint UnRefRuleC = 0;
const Uint RegRefRuleC= 63;


/// \brief Vertices of a given face, ordered as follows: \link #VertOfFaceAr (more...)\endlink
///
/// -# 4 faces of the tetrahedron itself; Face i does not contain vertex i
/// -# For faces 0 till 3, their subfaces as follows (all on one line):
///    -# four faces of the regular refinement ordered lexicographically
///    -# the six subfaces obtained by dividing one edge, sorted according
///       to the divided vertex and then according to the nondivided vertex.
///    -# the three subfaces from dividing two of the faces' edges
///    .
/// -# All inner subfaces, ordered lexicographically
extern const Ubyte VertOfFaceAr[NumAllFacesC][3];

/// \brief Vertices of a given edge (45 equals 10 over 2). \link #VertOfEdgeAr (more...)\endlink
///
/// - 6 edges of the father
/// - 12 subedges of these edges
/// - 12 connections on the faces: vertex with opposite midvertex
/// - 12 connections on the faces: two adjacent midvertices
/// - 3 diagonals (in this order so that the regular pattern has edge 30-42)
extern const Ubyte VertOfEdgeAr[NumAllEdgesC][2];

/// Edge with given vertices
extern const byte EdgeByVertAr[NumAllVertsC][NumAllVertsC];

/// Vertex as intersection of two edges; -1, if the edges are opposites; only the upper triangular part is meaningful.
extern const byte VertByEdgeAr[NumEdgesC][NumEdgesC];

/// Edges of a given Face
const Ubyte EdgeOfFaceAr[NumFacesC][3] = { {2,4,5}, {1,3,5}, {0,3,4}, {0,1,2} };

/// Faces intersecting in a given edge
const Ubyte FaceOfEdgeAr[NumEdgesC][2] = { {2,3}, {1,3}, {0,3}, {1,2}, {0,2}, {0,1} };

/// Orientation of cross-product when iterating over VertOfFace for a given face; 1 stands for outer normal
const byte OrientOfFaceAr[NumFacesC] = { 1, -1, 1, -1 };


/// Vertices of a given edge
inline Ubyte VertOfEdge    (Ubyte edge, Ubyte num) { return VertOfEdgeAr[edge][num]; }
/// Edge with given vertices
inline byte  EdgeByVert    (Ubyte v0, Ubyte v1) { return EdgeByVertAr[v0][v1]; }
/// Edge opposing a given edge
inline Ubyte OppEdge       (Ubyte edge) { return 5-edge; }
/// Opposite vertex of edge
inline Ubyte OppVertOfEdge (Ubyte edge, Ubyte vert)
    { return vert == VertOfEdge( edge, 0) ? VertOfEdge( edge, 1) : VertOfEdge( edge, 0); }
/// The intersection of two edges (-1 if the edges are opposites)
inline byte VertByEdge (Ubyte e0, Ubyte e1)
    { if (e0 > e1) std::swap( e0, e1); return VertByEdgeAr[e0][e1]; }
/// Vertices of a given face
inline Ubyte VertOfFace    (Ubyte face, Ubyte num) { return VertOfFaceAr[face][num]; }
/// Face with given vertices
inline Ubyte FaceByVert    (Ubyte v0, Ubyte v1, Ubyte v2) { return 6-v0-v1-v2; }
/// Faces intersecting in a given vertex
inline Ubyte FaceOfVert    (Ubyte vert, Ubyte num) { return num < vert ? num : num+1; }
/// Edges of a given Face
inline Ubyte EdgeOfFace    (Ubyte face, Ubyte num) { return EdgeOfFaceAr[face][num]; }
/// Returns, on which face of the parent chedge lies and also determines,
/// whether the barycenter of chedge lies at v0 (0), v1 (1) or v2 (2) of the face.
/// As to why this works, take a look at topo.h, specifically the sorting
/// of VertOfEdgeAr.
inline void WhichEdgeInFace(Uint chedge, Uint& parface, Uint& pos)
{ parface= chedge<30 ? chedge/3 -6 : chedge/3 -10; pos= chedge%3; }
/// Vertex opposing a given face
inline Ubyte OppVert       (Ubyte face) { return face; }
/// Face opposing a given vertex
inline Ubyte OppFace       (Ubyte vert) { return vert; }
/// Faces intersecting in a given edge
inline Ubyte FaceOfEdge    (Ubyte edge, Ubyte num) { return FaceOfEdgeAr[edge][num]; }
/// Orientation of cross-product when iterating over VertOfFace for a given face
inline byte OrientOfFace (Ubyte face) { return OrientOfFaceAr[face]; }
/// Number of the midvertex of a given edge
inline Ubyte MidVertOfEdge (Ubyte edge) { return edge+NumVertsC; }
/// returns the edge to which this midvertex belongs
inline Ubyte EdgeOfMidVert (Ubyte vert) { return vert-NumVertsC; }
/// true, if the given vertex is a midvertex of an edge
inline bool IsMidVert     (Ubyte vert) { return vert >= NumVertsC; }
/// true, if the edge is an edge of the parent-tetra
inline bool IsParentEdge(Ubyte edge) { return edge < NumEdgesC; }
/// returns, if edge is a sub-edge of an edge in the parent-tetrahedron
inline bool IsSubEdge     (Ubyte edge) { return edge < NumEdgesC+12 && edge > NumEdgesC-1; }
/// First or second subedge of a given edge
inline Ubyte SubEdge       (Ubyte edge, Ubyte num) { return NumEdgesC+2*edge+num; }
/// returns the number of the parent-edge
inline Ubyte ParentEdge    (Ubyte subedge) { return (subedge-NumEdgesC)/2; }
/// returns, whether subedge is subedge zero or one in the parent-edge
inline Ubyte NumOfSubEdge  (Ubyte subedge) { return subedge%2; }
/// true, if subedge is not an edge of the parent, not a subedge of a parent-edge and not a diagonal
inline bool IsSubInParFace(Ubyte subedge) { return subedge > 17 && subedge < 42; }
/// returns the number of the parent face for subedges, if IsSubInParFace is true
inline Ubyte ParFaceOfEdge(Ubyte subedge) { return (subedge - (subedge<30 ? 18 : 30))/3; }
/// true, if subedge is a diagonal of a tetrahedron
inline bool IsDiagonal    (Ubyte subedge) { return subedge > 41; }
/// true, if the face is a also face of the parent-tetra
inline bool IsParentFace(Ubyte face) { return face < NumFacesC; }
/// true, if the face lies in the interior of the parent-tetra
inline bool IsInnerFace(Ubyte face) { return face > 55; }
/// true, if the face is a subface of a parent-face, but not the parent-face itself
inline bool IsSubFace(Ubyte face) { return face >= NumFacesC && face <= 55; }
/// returns the number of the parent-face for faces, for which IsSubFace() is true
inline Ubyte ParentFace(Ubyte subface) { return (subface-NumFacesC)/NumAllSubFacesPerFaceC; }
/// returns the number of subface num of face 'face'; '13' is a kludge to refer to the
/// unrefined face 'face' itself
inline Ubyte SubFaceOfFace   (Ubyte face, Ubyte num) { return num == 13 ? face : 4 + face*NumAllSubFacesPerFaceC + num; }
/// returns the number of a subface relative to the numbering of subfaces in a face; thus a return value
/// of 13 refers to the unrefined parent-face itself
inline Ubyte NumOfSubFace(Ubyte subface) { return (subface-NumFacesC)%NumAllSubFacesPerFaceC; }
/// returns the number of subface num of face 'face'; '13' is a kludge to refer to the
/// unrefined face 'face' itself
inline Ubyte ChildOfFace   (Ubyte face, Ubyte num) { return num == 13 ? face : 4 + face*NumAllSubFacesPerFaceC + num; }
/// True, iff refinement-rule refpat subdivides the edge 'edge'
inline bool RefinesEdge(Ubyte refpat, Ubyte edge) { return (refpat>>edge) & 1; }
/// True, iff refinement-rule refpat subdivides the face 'face'
inline bool RefinesFace(Ubyte refpat, Ubyte face)
{
    return    (refpat>>EdgeOfFace(face, 0) & 1)
           || (refpat>>EdgeOfFace(face, 1) & 1)
           || (refpat>>EdgeOfFace(face, 2) & 1);
}


/// Refinement-rules

/// Describes a child in a refinement rule
struct ChildDataCL
{
    /// Number of each vertex
    Ubyte Vertices[NumVertsC];
    /// Number of each edge
    Ubyte Edges[NumEdgesC];
    /// Number of each face
    Ubyte Faces[NumFacesC];
};

/// Describes neighbourhood relations of children and faces in a refinement rule
struct RefRuleCL
{
    /// How many edges are used for the children (parent edges inclusive)
    Ubyte EdgeNum;
    /// VertIndices of the edges
    byte  Edges[MaxEdgesC];

    /// How many faces are used for the children (parent faces inclusive)
    Ubyte FaceNum;
    /// VertIndices of the faces
    byte  Faces[MaxFacesC];

    /// How many children does this rule generate
    Ubyte ChildNum;
    /// ChildData index of the children
    byte  Children[MaxChildrenC];
};

/// All refinement-rules
extern const byte RefRuleAr[];
extern const byte ChildDataAr[];

/// Obtain a refinement-rule for the edge-pattern 'refpattern'
inline const RefRuleCL& GetRefRule(Ubyte refpattern)
    { return *reinterpret_cast<const RefRuleCL*>( RefRuleAr + sizeof( RefRuleCL)*refpattern); }

inline const ChildDataCL& GetChildData(Ubyte VertIndex)
    { return *reinterpret_cast<const ChildDataCL*>( ChildDataAr + sizeof( ChildDataCL)*VertIndex); }


/// Transformation of barycentric coordinates (children <-> parent)

/// Data for the transformation matrices of barycentric coordinates from parent to child.
/// The accessor function SMatrixCL<4,4> parent_to_child_bary (Ubyte ch) is provided in geom/simplex.h.
extern const byte parent_to_child_bary_ar[];

/// Data for the transformation matrices of barycentric coordinates from child to parent.
/// Note, that the values must be scaled by 0.5.
/// The accessor function SMatrixCL<4,4> child_to_parent_bary (Ubyte ch) is provided in geom/simplex.h.
extern const byte child_to_parent_bary_ar[];


/// connectivity of the reference-edge

namespace RefEdge {

const Uint NumVertsC= 2;

/// Vertices of a given edge
inline Ubyte VertOfEdge (Ubyte, Ubyte num) { return num; }

/// Vertex opposing a given vertex
/// Note, that this is different from the n-dimensional cases, n > 2, as vertexes and facets are the same in 1d. Therefore, they should not be numbered differently. This is a natural ordering, which also could have been used for tetras. However, changing this now would be painful.
inline Ubyte OppVert (Ubyte vert) { return 1 - vert; }

} // end of namespace DROPS::RefEdge


/// connectivity of the reference-triangle

namespace RefTri {

const Uint NumVertsC= 3;
const Uint NumEdgesC= 3;

/// Vertices of a given edge
const Ubyte vert_of_edge_ar[NumEdgesC][2]= {
    {0,1}, {0,2}, {1,2}
};
inline Ubyte VertOfEdge (Ubyte edge, Ubyte num) { return vert_of_edge_ar[edge][num]; }

/// Edge opposing a given vertex
/// Note, that this is different from the higher-dimensional cases as edges and facets are the same in 2d. Therefore, they should not be numbered differently. This is a natural ordering, which also could have been used for tetras. However, changing this now would be painful.
inline Ubyte OppEdge (Ubyte vert) { return 2 - vert; }

} // end of namespace DROPS::RefTri


/// connectivity of the reference-pentatope

namespace RefPenta {

const Uint NumVertsC= 5;
const Uint NumEdgesC= 10;
const Uint NumTrisC=  10;
const Uint NumTetrasC= 5;

/// Vertices of a given edge
const Ubyte vert_of_edge_ar[NumEdgesC][2]= {
    {0,1}, {0,2}, {1,2}, {0,3}, {1,3}, {2,3}, // The first six edges are from the reference-tetra
    {0,4}, {1,4}, {2,4}, {3,4}
};
inline Ubyte VertOfEdge (Ubyte edge, Ubyte num) { return vert_of_edge_ar[edge][num]; }

/// Edge given by two vertices
const byte edge_by_vert_ar[NumVertsC][NumVertsC]= {
    {-1,  0,  1,  3,  6},
    { 0, -1,  2,  4,  7},
    { 1,  2, -1,  5,  8},
    { 3,  4,  5, -1,  9},
    { 6,  7,  8,  9, -1}
};
inline byte EdgeByVert (Ubyte v0, Ubyte v1) { return edge_by_vert_ar[v0][v1]; }

/// Vertices of a given tetra
const Ubyte vert_of_tetra_ar[NumVertsC][4]= {
    {1,2,3,4}, {0,2,3,4}, {0,1,3,4}, {0,1,2,4}, {0,1,2,3}
};
inline Ubyte VertOfTetra (Ubyte tetra, Ubyte num) { return vert_of_tetra_ar[tetra][num]; }

/// Tetra opposing a given vertex
inline Ubyte OppTetra (Ubyte vert) { return vert; }

/// Vertex of the penta defined by vertex vert of the facet tetra.
inline Ubyte VertByTetraVert (Ubyte tetra, Ubyte vert) { return vert + (tetra > vert ? 0 : 1); }
/// Edge of the penta defined by edge of the facet tetra.
inline Ubyte EdgeByTetraEdge (Ubyte tetra, Ubyte edge)
{
    return EdgeByVert( VertByTetraVert( tetra, DROPS::VertOfEdge( edge, 0)),
                       VertByTetraVert( tetra, DROPS::VertOfEdge( edge, 1)));
}

} // end of namespace DROPS::RefPenta


/// dimension-independent connectivity functions

/// Vertices of a given edge
template <Uint Dim>
inline Ubyte VertOfEdge (Ubyte edge, Ubyte num); ///< not defined, but full specializations for Dim=1..4 are.

template <>
inline Ubyte VertOfEdge<1> (Ubyte, Ubyte num)
{ return RefEdge::VertOfEdge( 0, num); }

template <>
inline Ubyte VertOfEdge<2> (Ubyte edge, Ubyte num)
{ return RefTri::VertOfEdge( edge, num); }

template <>
inline Ubyte VertOfEdge<3> (Ubyte edge, Ubyte num)
{ return VertOfEdge( edge, num); }

template <>
inline Ubyte VertOfEdge<4> (Ubyte edge, Ubyte num)
{ return RefPenta::VertOfEdge( edge, num); }

} // end  of namespace DROPS
#endif
