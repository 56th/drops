/// \file simplex.cpp
/// \brief classes that constitute the simplex
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Eva Loch, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

/// Remarks: We should use the const-qualifier to make it difficult to
///          accidentally change the multigrid structure from anywhere
///          outside of the multigrid algorithms.
///          Thus the pointer to user data structures should probably be
///          a pointer to mutable.

#ifdef _PAR
#include "parallel/parmultigrid.h"
#include "parallel/parallel.h"
#endif

#include "geom/multigrid.h"
#include "num/gauss.h"
#include <iterator>

//for curved
#include "geom/deformation.h"
#include "geom/isoparamP2.h" 

namespace DROPS
{

//
// static members of TetraCL
//
SArrayCL<EdgeCL*, NumAllEdgesC> TetraCL::ePtrs_(static_cast<EdgeCL*>(0));
SArrayCL<FaceCL*, NumAllFacesC> TetraCL::fPtrs_(static_cast<FaceCL*>(0));

// V E R T E X  C L
// ----------------

#ifdef _PAR

/** Puts the GID, mark for removement, and boundary information onto the stream.
*/
void VertexCL::Pack( DiST::MPIostreamCL& ostrstream) const
{
    ostrstream << GetGID() << RemoveMark_;
    if ( IsOnBoundary()) {
        ostrstream << BndVerts_->size();
        for ( const_BndVertIt it( GetBndVertBegin()), end(GetBndVertEnd()); it!=end; ++it) {
            ostrstream << it->GetBndIdx() << it->GetCoord2D()[0] << it->GetCoord2D()[1];
        }
    }
    else {
        ostrstream << size_t(0);
    }
}

/** Reads the GID, mark for removement, and boundary information from the stream.
*/
void VertexCL::UnPack( DiST::MPIistreamCL& istrstream)
{
    size_t numBnd=0;
    istrstream >> gid_ >> RemoveMark_ >> numBnd;
    BndIdxT bidx;
    Point2DCL p2d;
    // delete old bnd verts if they exist
    if (BndVerts_) {
        delete BndVerts_;
        BndVerts_= 0;
    }
    for ( size_t i=0; i<numBnd; ++i) {
        istrstream >> bidx >> p2d[0] >> p2d[1];
        AddBnd( BndPointCL(bidx, p2d));
    }
}

void VertexCL::UpdateGID()
{
    if (DiST::InfoCL::Instance().Exists(this->GetGID()))
        throw DROPSErrCL("VertexCL::UpdateGID(): object is already registered");
    DiST::GeomIdCL hash(this->GetLevel(), *this);
    this->gid_ = hash;
}

#endif

void VertexCL::ChangeCoord (__UNUSED__ Point3DCL& p)
{
#ifdef _PAR
    throw DROPSErrCL("VertexCL::ChangeCoord in parallel not implemented, yet");
#else
    Coord_ = p;
#endif
}

// E D G E  C L
// ------------

// Assumes that the next Level with respect to ep exists.
void EdgeCL::BuildMidVertex(SimplexFactoryCL& factory, const BoundaryCL& Bnd)
// TODO: For nonlinear boundaries, we must treat the case, in which an edge lies in two
//       boundary-segments in a different manner: Project to the common "edge" of the
//       boundary-segments!
// TODO: Due to MeshReader-Boundaries, which have no 2D-parameterization, we calculate
//       the barycenter of the edge directly. This, of course, breaks nonlinear boundary
//       segments.
{
    const VertexCL* const vp0 ( GetVertex( 0));
    const VertexCL* const vp1 ( GetVertex( 1));

    if (IsOnBoundary()) {
        const BndIdxT bndidx= *GetBndIdxBegin();
        const BndPointCL& bndvert0= *std::find_if( vp0->GetBndVertBegin(), vp0->GetBndVertEnd(),
                                                   BndPointSegEqCL( bndidx));
        const BndPointCL& bndvert1= *std::find_if( vp1->GetBndVertBegin(), vp1->GetBndVertEnd(),
                                                   BndPointSegEqCL( bndidx));
        BndPairCL bndpair= Bnd.GetBndSeg( bndidx)->MidProject( bndvert0, bndvert1);
        // XXX: Revise this for nonlinear boundary-segments.
        bndpair.second= GetBaryCenter( *this);
        VertexCL& newMidVert= factory.MakeVertex( bndpair.second, GetLevel()+1);
        SetMidVertex( &newMidVert);
        newMidVert.AddBnd( BndPointCL( bndidx, bndpair.first));
        if ( std::distance( GetBndIdxBegin(), GetBndIdxEnd()) == 2 ) {
            const BndIdxT bndidx= *(GetBndIdxBegin() + 1);
            const BndPointCL& bndvert0= *std::find_if( vp0->GetBndVertBegin(), vp0->GetBndVertEnd(),
                                                       BndPointSegEqCL( bndidx));
            const BndPointCL& bndvert1= *std::find_if( vp1->GetBndVertBegin(), vp1->GetBndVertEnd(),
                                                       BndPointSegEqCL( bndidx));
            BndPairCL bndpair1= Bnd.GetBndSeg( bndidx)->MidProject( bndvert0, bndvert1);
            newMidVert.AddBnd( BndPointCL( bndidx, bndpair1.first));
            newMidVert.BndSort();
//            Assert( bndpair.second == bndpair1.second, DROPSErrCL("BuildMidVertex: Projection leads to different 3D-coords."), ~0 );
        }
    }
    else {
        VertexCL& newMidVert= factory.MakeVertex( BaryCenter( vp0->GetCoord(), vp1->GetCoord()), GetLevel() + 1);
        SetMidVertex( &newMidVert);
    }
#ifdef _PAR
//    throw DROPSErrCL("EdgeCL::BuildMidVertex: Identify is missing!");
    // new created Vertex must be identified with DDD
    if ( IsOnProcBnd())
        ParMultiGridCL::Instance().IdentifyVertex( this);
#endif
}

Point3DCL GetBaryCenter(const EdgeCL& e)
{
    return (e.GetVertex(0)->GetCoord() + e.GetVertex(1)->GetCoord() )*0.5;
}

#ifdef _PAR

/** Put the GID, both vertices, boundary information, accumulated MFR, and mark for removement
    onto the stream
*/
void EdgeCL::Pack( DiST::MPIostreamCL& ostrstream) const
{
    ostrstream << GetGID()
               << GetVertex(0)->GetGID()
               << GetVertex(1)->GetGID()
               << (GetMidVertex()==0 ? DiST::NoGID : GetMidVertex()->GetGID())
               << Bnd_[0] << Bnd_[1] << AccMFR_ << RemoveMark_;
}

/** Reads the GID, both vertices, boundary information, accumulated MFR, and mark for removement
    from the stream
*/
void EdgeCL::UnPack( DiST::MPIistreamCL& istrstream)
{
    DiST::GeomIdCL vertex0, vertex1, midVertex;
    istrstream >> gid_ >> vertex0 >> vertex1 >> midVertex;
    istrstream >>  Bnd_[0] >>  Bnd_[1] >> AccMFR_ >> RemoveMark_;
    Vertices_[0]= DiST::InfoCL::Instance().GetVertex(vertex0);
    Vertices_[1]= DiST::InfoCL::Instance().GetVertex(vertex1);
    if ( DiST::InfoCL::Instance().Exists(midVertex))
        MidVertex_= DiST::InfoCL::Instance().GetVertex(midVertex);
}

void EdgeCL::UpdateGID()
{
    if (DiST::InfoCL::Instance().Exists(this->GetGID()))
        throw DROPSErrCL("EdgeCL::UpdateGID(): object is already registered");
    DiST::GeomIdCL hash(this->GetLevel(), *this);
    this->gid_ = hash;
}


#endif

// F A C E  C L
// ------------

// Assumptions:
// - a green child and its parent with a common face are both stored on the same side,
//   i.e. 0 and 2 lies on one side and 1 and 3 on the other one

void FaceCL::LinkTetra(const TetraCL* tp)
{
    int offset=0;

    if (tp->GetLevel()==GetLevel()) // tetra on same level
    {
#ifndef _PAR
        if (Neighbors_[0] && Neighbors_[0]!=tp) // in sequential version:  if _Neighbors[0]!=0 then  _Neighbors[0]!=tp is always true
            offset=1;
#else
        if (Neighbors_[2]==0 && Neighbors_[3]==0) {
            if (Neighbors_[0] && Neighbors_[0]!=tp) // in sequential version:  if Neighbors_[0]!=0 then  Neighbors_[0]!=tp is always true
                offset=1;
        }
        else { // during transfer, the children are unpacked before their parents. Hence, we have to find the position where one of the children is stored on the next level.
        	if (tp->GetChildBegin()==tp->GetChildEnd()) { // no children known to tetra -> use position where no children are linked
        		offset= Neighbors_[2] ? 1 : 0;
                Assert( Neighbors_[offset+2]==0, DiST::ErrorCL("FaceCL::LinkTetra: Occupied child position ", offset+2, GetGID()), DebugRefineEasyC);
        	} else if (!is_in( tp->GetChildBegin(), tp->GetChildEnd(), Neighbors_[2]))
            	if (Neighbors_[3]==0 || is_in( tp->GetChildBegin(), tp->GetChildEnd(), Neighbors_[3]))
            		offset=1;
            Assert( Neighbors_[offset+2]==0 || is_in( tp->GetChildBegin(), tp->GetChildEnd(), Neighbors_[offset+2]), DiST::ErrorCL("FaceCL::LinkTetra: Wrong child at position ", offset+2, GetGID()), DebugRefineEasyC);
        }
#endif
    }
#ifndef _PAR
    else                            // green child of parent
    {
        Assert(tp->GetLevel() == GetLevel()+1, DROPSErrCL("FaceCL::LinkTetra: Illegal level of green tetra"), DebugRefineEasyC);
        // tetra is stored on the same side as the parent
        offset= Neighbors_[0]==tp->GetParent() ? 2 : 3;
    }
#else
    else {                          // green child of parent
        Assert(tp->GetLevel() == GetLevel()+1, DROPSErrCL("FaceCL::LinkTetra: Illegal level of green tetra"), DebugRefineEasyC);
        if (tp->GetParent()) {
        // tetra is stored on the same side as the parent
            offset= Neighbors_[0]==tp->GetParent() ? 2 : 3;
            Assert( tp->GetParent()==Neighbors_[offset-2], DiST::ErrorCL("FaceCL::LinkTetra: Wrong parent at position ", offset-2, GetGID()), DebugRefineEasyC);
        } else { // during transfer, the children are unpacked before their parents. Hence, we need an empty position without parent.
            offset= (Neighbors_[0]==0 && (Neighbors_[2]==0 || Neighbors_[2]==tp)) ? 2 : 3;
            Assert( !Neighbors_[offset-2], DiST::ErrorCL("FaceCL::LinkTetra: Occupied parent position ", offset-2, GetGID()), DebugRefineEasyC);
        }
    }
#endif
    Assert(!Neighbors_[offset] || Neighbors_[offset]==tp, DROPSErrCL("FaceCL::LinkTetra: Link occupied by another tetra!"), DebugRefineEasyC);
    Neighbors_[offset]= tp;
}


void FaceCL::UnlinkTetra(const TetraCL* tp)
/// \todo (of) Kann es passieren, dass ein Tetraeder nicht angelinkt ist und
///   dennoch entfernt werden soll? Assert ersteinmal stehen gelassen ...
/// If tetra is not linked do nothing
{
    int i;

    // find tetra in neighbor list
    for ( i=0; i<4; ++i )
        if (Neighbors_[i] == tp) break;

    // if tetra is not linked, do nothing.
    if (i==4) return;

    Assert(i<4, DROPSErrCL("FaceCL::UnlinkTetra: No such tetra."), DebugRefineEasyC);

    if (i == 0)
    {
		i= 1;
#ifndef _PAR
        Neighbors_[0]= Neighbors_[1];
        Neighbors_[2]= Neighbors_[3];
        Neighbors_[3]= 0;
#else
        Neighbors_[0]= Neighbors_[1];
        std::swap( Neighbors_[2], Neighbors_[3]);
#endif
    }
    Neighbors_[i]= 0;
}

const TetraCL* FaceCL::GetNeighInTriang(const TetraCL* tp, Uint TriLevel) const
/** Get a pointer to the neighbor tetrahedra of a tetraeder \a tp in the triangulation level
    specified by parameter \a TriLevel.
    \pre tetrahedra and face must exist in triangulation level and face must not
     lie on boundary
*/
{
    Assert( tp->IsInTriang(TriLevel) && IsInTriang(TriLevel),
            DROPSErrCL("FaceCL::GetNeighInTriang: Face or Tetra not in triangulation!"), DebugRefineEasyC);
    Assert( !IsOnBoundary(),
            DROPSErrCL("FaceCL::GetNeighInTriang: Face of tetra lies on boundary!"), DebugRefineEasyC);

    const Uint oppSide= Neighbors_[ tp->GetLevel()==GetLevel() ? 0 : 2 ]==tp;

    return Neighbors_[oppSide]->IsInTriang(TriLevel) ?
           Neighbors_[oppSide] : Neighbors_[oppSide+2];
}


#ifdef _PAR

// parallel functions for faces
// ----------------------------

/** Puts the GID, all neighbors, boundary information, and mark for removement onto the stream
*/
void FaceCL::Pack( DiST::MPIostreamCL& ostrstream) const
{
    ostrstream << GetGID();
    for ( Uint i=0; i<4; ++i) {
        ostrstream << ( GetNeighbor(i)==0 ? DiST::NoGID : GetNeighbor(i)->GetGID());
    }
    ostrstream << Bnd_ << RemoveMark_;
}

/** Reads the GID, all neighbors, boundary information, and mark for removement from the stream
*/
void FaceCL::UnPack( DiST::MPIistreamCL& istrstream)
{
    DiST::GeomIdCL tmpNeigh;
    istrstream >> gid_;
    for ( size_t i=0; i<4; ++i)
        istrstream >> tmpNeigh;
    // neighboring tetras are linked later by TransferCL::CreateSimplex<TetraCL>(...)
    istrstream >> Bnd_ >> RemoveMark_;
}

void FaceCL::UpdateGID()
{
    if (DiST::InfoCL::Instance().Exists(this->GetGID()))
        throw DROPSErrCL("FaceCL::UpdateGID(): object is already registered");
    DiST::GeomIdCL hash(this->GetLevel(), *this);
    this->gid_ = hash;
}


#endif

Point3DCL GetBaryCenter(const FaceCL& f)
{
    const TetraCL* const tp= f.GetSomeTetra();
    const Uint face= f.GetFaceNumInTetra(tp);

    return ( tp->GetVertex( VertOfFace(face, 0))->GetCoord()
           + tp->GetVertex( VertOfFace(face, 1))->GetCoord()
           + tp->GetVertex( VertOfFace(face, 2))->GetCoord() )/3.0;
}

// T E T R A  C L
// --------------

Point3DCL GetBaryCenter(const TetraCL& t)
{
    return 0.25*( t.GetVertex(0)->GetCoord() + t.GetVertex(1)->GetCoord()
                + t.GetVertex(2)->GetCoord() + t.GetVertex(3)->GetCoord() );
}

Point3DCL GetBaryCenter(const TetraCL& t, Uint face)
{
    return ( t.GetVertex(VertOfFace(face, 0))->GetCoord()
            +t.GetVertex(VertOfFace(face, 1))->GetCoord()
            +t.GetVertex(VertOfFace(face, 2))->GetCoord() )/3.;
}

Point3DCL GetWorldCoord(const TetraCL& t, const SVectorCL<3>& c)
{
    static MeshDeformationCL & m = MeshDeformationCL::getInstance();
    if (m.IsUsed())
    {
        if (!m.IsTetraCurved(t))
            return GetWorldCoord(m.GetLocalP1Deformation(t),c);
        else    
            return GetWorldCoord(m.GetLocalP2Deformation(t),c);
    }
    else
        return (1. -c[0] -c[1] -c[2])*t.GetVertex(0)->GetCoord()
            +c[0]*t.GetVertex(1)->GetCoord()
            +c[1]*t.GetVertex(2)->GetCoord()
            +c[2]*t.GetVertex(3)->GetCoord();
}

Point3DCL GetWorldCoord(const TetraCL& t, const SVectorCL<4>& c)
{
    static MeshDeformationCL & m = MeshDeformationCL::getInstance();
    if (m.IsUsed())
    {
        if (!m.IsTetraCurved(t))
            return GetWorldCoord(m.GetLocalP1Deformation(t),c);
        else    
            return GetWorldCoord(m.GetLocalP2Deformation(t),c);
    }
    else
        return c[0]*t.GetVertex(0)->GetCoord()
            +c[1]*t.GetVertex(1)->GetCoord()
            +c[2]*t.GetVertex(2)->GetCoord()
            +c[3]*t.GetVertex(3)->GetCoord();
}

Point3DCL GetRefCoord(const TetraCL& t, const SVectorCL<4>& c)
{
        return c[0]*t.GetVertex(0)->GetCoord()
            +c[1]*t.GetVertex(1)->GetCoord()
            +c[2]*t.GetVertex(2)->GetCoord()
            +c[3]*t.GetVertex(3)->GetCoord();
}

Point3DCL GetWorldCoord(const TetraCL& t, Uint face, const SVectorCL<2>& c)
{
    static MeshDeformationCL & m = MeshDeformationCL::getInstance();
    if (m.IsUsed())
    {
        if (!m.IsTetraCurved(t))
            return GetWorldCoord(m.GetLocalP1Deformation(t),FaceToTetraCoord(t,face,c));
        else    
            return GetWorldCoord(m.GetLocalP2Deformation(t),FaceToTetraCoord(t,face,c));
    }
    else
        return (1. -c[0] -c[1])*t.GetVertex(VertOfFace(face, 0))->GetCoord()
            +c[0]*t.GetVertex(VertOfFace(face, 1))->GetCoord()
            +c[1]*t.GetVertex(VertOfFace(face, 2))->GetCoord();
}

SVectorCL<3> FaceToTetraCoord(__UNUSED__ const TetraCL& t, Uint f, SVectorCL<2> c)
{
    SVectorCL<3> ret(0.);
    switch(f)
    {
      case 0:
        ret[0]= 1 -c[0] -c[1]; ret[1]= c[0]; ret[2]= c[1]; break;
      case 1:
        ret[1]= c[0]; ret[2]= c[1]; break;
      case 2:
        ret[0]= c[0]; ret[2]= c[1]; break;
      case 3:
        ret[0]= c[0]; ret[1]= c[1]; break;
      default: throw DROPSErrCL("FaceToTetraCoord: illegal face-number.");
    }
    Assert( (GetWorldCoord(t,f,c)-GetWorldCoord(t, ret)).norm_sq() < 1.e-15, DROPSErrCL("FaceToTetraCoord: inconsistent mapping!"), DebugNumericC);
    return ret;
}


World2BaryCoordCL::World2BaryCoordCL (const TetraCL& t)
{
    SMatrixCL<4,4>& m= qr_.GetMatrix();
    for (Uint v= 0; v < 4; ++v) {
        const Point3DCL& p= t.GetVertex( v)->GetCoord();
        for (Uint i= 0; i < 3; ++i) m( i, v)= p[i];
    }
    for (Uint j= 0; j < 4; ++j) m( 3, j)= 1.;
    qr_.prepare_solve();
}

World2BaryCoordCL::World2BaryCoordCL(const Point3DCL* coordVerts)
{
    SMatrixCL<4,4>& m= qr_.GetMatrix();
    for (Uint v= 0; v < 4; ++v) {
        const Point3DCL& p= coordVerts[v];
        for (Uint i= 0; i < 3; ++i) m( i, v)= p[i];
    }
    for (Uint j= 0; j < 4; ++j) m( 3, j)= 1.;
    qr_.prepare_solve();
}

BaryCoordCL World2BaryCoordCL::operator() (const Point3DCL& p) const
{
    BaryCoordCL r( MakeBaryCoord( p[0], p[1], p[2], 1.));
    qr_.Solve( r);
    return r;
}

void ComputeChildFacesOfFace (const TetraCL& p, Uint f, std::vector<const FaceCL*>& childfaces)
{
    if (p.IsUnrefined()) {
        childfaces.push_back( p.GetFace( f));
        return;
    }

    bool done= false; // Used to shortcut the search for child faces, if the child face coincides with the parent face.
    const RefRuleCL& ref= p.GetRefData();
    for (Uint c= 0; c < ref.ChildNum && !done; ++c) { // Loop over all children
        const ChildDataCL& ch= GetChildData( ref.Children[c]);
        for (Uint i= 0; i < NumFacesC; ++i) {
            const Uint childFace= ch.Faces[i];
            if (IsParentFace( childFace)) {
                done= true;
                childfaces.push_back( p.GetChild( c)->GetFace( i)); // The child-face is identical to the face f of p.
            }
            else if (IsSubFace( childFace) && f == ParentFace( childFace))
                childfaces.push_back( p.GetChild( c)->GetFace( i)); // This child-face refines the face f of p.
        }
    }
}

double TetraCL::GetVolume () const
{
    Point3DCL v[3];
    const Point3DCL p0( Vertices_[0]->GetCoord() );
    for (Uint i=1; i<NumVertsC; ++i)
        v[i-1] = Vertices_[i]->GetCoord() - p0;
    return std::fabs(  v[0][0] * (v[1][1]*v[2][2] - v[1][2]*v[2][1])
                     - v[0][1] * (v[1][0]*v[2][2] - v[1][2]*v[2][0])
                     + v[0][2] * (v[1][0]*v[2][1] - v[1][1]*v[2][0]) ) / 6.0;
}

double TetraCL::GetNormal(Uint face, Point3DCL& normal, double& dir) const
/**
dir is set to 1.0 if the normal points out of the tetra,
else it is set to -1.0.
Normal has unit length, but the length of the cross-product is returned,
which is useful for integration over that face.
If the triangulation is consistently numbered, both tetras on a face will
return the same normal in "normal" (, of course "dir" will be different).
*/
{
    const VertexCL* v[3];
    for (Uint i=0; i<3; ++i)
        v[i]= GetVertex( VertOfFace(face, i) );
    const VertexCL* const opvert= GetVertex( OppVert(face) );

    cross_product(normal, v[1]->GetCoord()-v[0]->GetCoord(),
                          v[2]->GetCoord()-v[0]->GetCoord());
    dir= inner_prod(opvert->GetCoord() - v[0]->GetCoord(), normal) < 0.0 ? 1. : -1.;
    const double absdet2D= normal.norm();
    normal/= absdet2D;
    return absdet2D;
}

double TetraCL::GetOuterNormal(Uint face, Point3DCL& normal) const
/**
Returns the length of the cross-product;
"normal" is set to the unit outward normal of face "face"
*/
{
    double dir;
    const double absdet2D= GetNormal(face, normal, dir);
    normal*= dir;
    return absdet2D;
}

void TetraCL::BuildEdges( SimplexFactoryCL& factory)
{
    for (Uint edge=0; edge<NumEdgesC; ++edge)
    {
        VertexCL* const v0= Vertices_[VertOfEdge(edge,0)];
        VertexCL* const v1= Vertices_[VertOfEdge(edge,1)];
        if ( !(Edges_[edge]= v0->FindEdge(v1) ) )
        {
            std::vector<BndPointCL> commonBndVerts;
            if ( v0->IsOnBoundary() && v1->IsOnBoundary() )
            {
                commonBndVerts.reserve(2);
                std::set_intersection( v0->GetBndVertBegin(), v0->GetBndVertEnd(),
                                       v1->GetBndVertBegin(), v1->GetBndVertEnd(),
                                       std::back_inserter(commonBndVerts), BndPointSegLessCL() );
            }
            EdgeCL* newEdge=0;
            switch( commonBndVerts.size() )
            {
              case 0:
                  newEdge= &factory.MakeEdge( v0, v1, GetLevel());
                  break;
              case 1:
                  newEdge= &factory.MakeEdge( v0, v1, GetLevel(), commonBndVerts[0].GetBndIdx());
                  break;
              case 2:
                  newEdge= &factory.MakeEdge( v0, v1, GetLevel(), commonBndVerts[0].GetBndIdx(), commonBndVerts[1].GetBndIdx());
                  break;
              default:
                v0->DebugInfo( std::cout); std::cout << std::endl;
                v1->DebugInfo( std::cout); std::cout << std::endl;
                throw DROPSErrCL("TetraCL::BuildEdges: Found edge on more than two BndSegs!");
            }
            Edges_[edge]= newEdge;
            Edges_[edge]->RecycleMe();
        }
    }
}

BndIdxT GetCommonBndSeg(const VertexCL* v0, const VertexCL* v1, const VertexCL* v2)
{
    if (!(v0->IsOnBoundary() && v1->IsOnBoundary() && v2->IsOnBoundary() ))
        return NoBndC;
    std::list<BndPointCL> intersec01, intersec012;
    std::set_intersection( v0->GetBndVertBegin(), v0->GetBndVertEnd(),
                           v1->GetBndVertBegin(), v1->GetBndVertEnd(),
                           std::back_inserter(intersec01), BndPointSegLessCL() );
    std::set_intersection( v2->GetBndVertBegin(), v2->GetBndVertEnd(),
                           intersec01.begin(), intersec01.end(),
                           std::back_inserter(intersec012), BndPointSegLessCL() );
    if (intersec012.empty() )
        return NoBndC;
    Assert( intersec012.size()==1,
        DROPSErrCL("GetCommonBndSeg: found more than one BndSeg connected to three different boundary vertices"), ~0);
    return intersec012.begin()->GetBndIdx();
}

void TetraCL::BuildAndLinkFaces ( SimplexFactoryCL& factory)  // used by XXXBuilderCL
{
    for (Uint face=0; face<NumFacesC; ++face)
    {
        VertexCL* const vp0= GetVertMidVert(VertOfFace(face, 0));
        VertexCL* const vp1= GetVertMidVert(VertOfFace(face, 1));
        VertexCL* const vp2= GetVertMidVert(VertOfFace(face, 2));

        if ( !(Faces_[face]= vp0->FindFace(vp1,vp2) ) )
        {
#ifndef _PAR
            FaceCL& newFace=factory.MakeFace( GetLevel(), GetCommonBndSeg(vp0, vp1, vp2) );
#else
            FaceCL& newFace=factory.MakeFace( GetLevel(), vp0->GetCoord(), vp1->GetCoord(), vp2->GetCoord(), GetCommonBndSeg(vp0, vp1, vp2));
#endif
            Faces_[face]= &newFace;
            Faces_[face]->RecycleMe(vp0, vp1, vp2);
        }
        Faces_[face]->LinkTetra(this);
    }
}

void TetraCL::SetChild(Uint c, TetraCL* cp)
{
    if (!Children_) Children_= new SArrayCL<TetraCL*, MaxChildrenC>;
    (*Children_)[c]= cp;
}


// member functions for r e f i n e m e n t

void TetraCL::RecycleReusables()
/**
is called, if the refinement rule has changed.
It recycles and rescues simplices, that will be reused:
- Edges of children, that are not edges of the parent, if they are used by the new rule
- Faces of children, that are not faces of the parent,  -- " --
- Children                                           ,  -- " --
*/
{
    const RefRuleCL& myRule= GetRefData();
    const RefRuleCL& newRule= DROPS::GetRefRule(GetRefMark() & 63);

    // prepare a list of the faces that must be recycled
    std::list<byte> commonEdges;
    std::list<byte> commonFaces;
    std::list<byte> commonChildren;

    typedef std::list<byte>::iterator SetIterT;

    std::set_intersection( myRule.Edges,  myRule.Edges +  myRule.EdgeNum,
                          newRule.Edges, newRule.Edges + newRule.EdgeNum,
                          std::back_inserter(commonEdges) );
    std::set_intersection( myRule.Faces,  myRule.Faces +  myRule.FaceNum,
                          newRule.Faces, newRule.Faces + newRule.FaceNum,
                          std::back_inserter(commonFaces) );
    std::set_intersection( myRule.Children,  myRule.Children +  myRule.ChildNum,
                          newRule.Children, newRule.Children + newRule.ChildNum,
                          std::back_inserter(commonChildren) );

#ifdef _PAR
    Assert(Children_, DiST::ErrorCL( "TetraCL::RecycleReusables: no children array", this->GetGID()), DebugDiSTC);
    ParMultiGridCL& pmg= ParMultiGridCL::Instance();
#endif
    for (Uint ch=0; ch<myRule.ChildNum; ++ch)
    {
        const ChildDataCL childdat= GetChildData(myRule.Children[ch]);
        TetraCL* const child= (*Children_)[ch];
        // recycle and rescue common edges
        for (Uint edge= 0; edge<NumEdgesC; ++edge)
        {
            if ( IsParentEdge(childdat.Edges[edge]) ) continue;
            SetIterT it= std::lower_bound( commonEdges.begin(), commonEdges.end(), childdat.Edges[edge]);
            if (it != commonEdges.end() && *it == childdat.Edges[edge])
            {
#ifdef _PAR
                // take care that sub simplices are not unregistered by DiST::ModifyCL
                pmg.Keep( child->Vertices_[VertOfEdge(edge,0)]);
                pmg.Keep( child->Vertices_[VertOfEdge(edge,1)]);
                pmg.Keep( child->Edges_[edge]);
#endif
                child->Vertices_[VertOfEdge(edge,0)]->ClearRemoveMark();
                child->Vertices_[VertOfEdge(edge,1)]->ClearRemoveMark();
                child->Edges_[edge]->ClearRemoveMark();
                child->Edges_[edge]->RecycleMe();
                commonEdges.erase(it);  // because edge is now already recycled
            }
        }
        // recycle and rescue common faces
        for (Uint face=0; face<NumFacesC; ++face)
        {
            if ( IsParentFace(childdat.Faces[face]) ) continue;
            SetIterT it= std::lower_bound(commonFaces.begin(), commonFaces.end(), childdat.Faces[face]);
            if (it != commonFaces.end() && *it == childdat.Faces[face])
            {
                VertexCL* const vp0= child->Vertices_[VertOfFace(face, 0)];
                VertexCL* const vp1= child->Vertices_[VertOfFace(face, 1)];
                VertexCL* const vp2= child->Vertices_[VertOfFace(face, 2)];

#ifdef _PAR
                // take care that sub simplices are not unregistered by DiST::ModifyCL
                pmg.Keep(child->Faces_[face]);
#endif
                child->Faces_[face]->ClearRemoveMark();
                child->Faces_[face]->RecycleMe(vp0,vp1,vp2);
                commonFaces.erase(it);  // because face is now already recycled
            }
        }
    }
    // recycle, rescue and unlink children
    for (Uint ch=0; ch<myRule.ChildNum; ++ch)
    {
        if ( std::binary_search(commonChildren.begin(), commonChildren.end(), myRule.Children[ch]) )
        {
            TetraCL* const child= (*Children_)[ch];

            child->SetNoRefMark();
            child->UnlinkFromFaces();
            child->RecycleMe();
        }
    }
}


void TetraCL::ClearAllRemoveMarks()
/**
is called if MarkEqRule(). It safes all sub simplices from removement.
*/
{
    const RefRuleCL& myRule= GetRefData();

    for (Uint ch=0; ch<myRule.ChildNum; ++ch)
    {
        TetraCL* const child= (*Children_)[ch];

        for (VertexPIterator vertPIt(child->Vertices_.begin()); vertPIt!=child->Vertices_.end(); ++vertPIt)
            (*vertPIt)->ClearRemoveMark();
        for (Uint edge= 0; edge<NumEdgesC; ++edge)
            child->Edges_[edge]->ClearRemoveMark();
        for (Uint face= 0; face<NumFacesC; ++face)
            child->Faces_[face]->ClearRemoveMark();
    }
}


void TetraCL::RestrictMark()
{
    if ( IsUnrefined() )
    {
        if ( IsMarkedForRegRef() && IsRegular() )
            CommitRegRefMark();
        if ( GetLevel() == 0 && IsMarkedForRemovement() )
            SetNoRefMark();
    }
    else
    {
#ifdef _PAR
        // if tetra has a ghost tetra, the marks are computed on the ghost tetra
        if ( HasGhost() ) return;
#endif
        if ( IsRegularlyRef() )
        {
            bool keepAnyChild= false;
            for (ChildPIterator childPIt=GetChildBegin(); childPIt!=GetChildEnd(); ++childPIt)
                if ( (*childPIt)->IsMarkedForRemovement() )
                    (*childPIt)->SetNoRefMark();
                else keepAnyChild= true;
            if (!keepAnyChild)
            {
                SetNoRefMark();
#ifdef _PAR
                if ( IsGhost() ) return;
#endif
                UnCommitRegRefMark();
            }
        }
        else
        { // tetra is irregularly refined
            bool setregrefmark= false;
            for (ChildPIterator chp= GetChildBegin(), end= GetChildEnd(); chp!=end; ++chp)
            {
                if (!setregrefmark) {
                    if ( (*chp)->IsMarkedForRef() )
                        setregrefmark= true;
                    else
                    {
                        for (const_EdgePIterator ep= (*chp)->GetEdgesBegin(), edgeend= (*chp)->GetEdgesEnd(); !setregrefmark && ep!=edgeend; ++ep)
                            if ( (*ep)->IsMarkedForRef() && (*ep)->GetLevel()!=GetLevel() )
                            // parent edges are ignored
                            {
                                setregrefmark= true;
                                break;
                            }
                    }
                }
                (*chp)->SetNoRefMark();
            }
            if (setregrefmark)
            {
                SetRegRefMark();
#ifdef _PAR
                if ( IsGhost() ) return;
#endif
                CommitRegRefMark();
            }
            else
                SetNoRefMark();
        }
    }
}

void TetraCL::CollectEdges (const RefRuleCL& refrule,
                            SimplexFactoryCL& factory,
                            const BoundaryCL& Bnd)
/**
The edges for new refinement are stored in the static TetraCL::ePtrs array.
First look for them in the recycle bin (maybe they were created before),
if the edge cannot be found, create it.
*/
/// \todo (of): Ist auf verschiedenen Prozessen tatsaechlich die Reihenfolge der Subedges eindeutig???
{
    const Uint nextLevel= GetLevel()+1;

    // Collect obvious edges
    for (Uint edge=0; edge<NumEdgesC; ++edge)
    {
        VertexCL* const vp0= Vertices_[VertOfEdge(edge, 0)];
        VertexCL* const vp1= Vertices_[VertOfEdge(edge, 1)];
        EdgeCL*   const ep = Edges_[edge];

        if ( ep->IsMarkedForRef() )
        {
            if ( ep->IsRefined() )
            {
                ePtrs_[SubEdge(edge, 0)]= vp0->FindEdge(ep->GetMidVertex());
                ePtrs_[SubEdge(edge, 1)]= ep->GetMidVertex()->FindEdge(vp1);
                Assert(ePtrs_[SubEdge(edge, 0)], DROPSErrCL("CollectEdges: SubEdge 0 not found."), DebugRefineEasyC);
                Assert(ePtrs_[SubEdge(edge, 1)], DROPSErrCL("CollectEdges: SubEdge 1 not found."), DebugRefineEasyC);
            }
            else
            {
                SArrayCL<EdgeCL*,2> subedge;
                ep->BuildSubEdges(factory, Bnd, subedge);
                for (size_t i=0; i<2; ++i) {
                    ePtrs_[SubEdge(edge, i)]= subedge[i];
                    subedge[i]->RecycleMe();
                }
#ifdef _PAR
                // if new edges are created on the proc-boundary, identify them
                if ( ep->IsOnProcBnd())
                {
                    ParMultiGridCL::Instance().IdentifyEdge( ePtrs_[SubEdge(edge, 0)], ep);
                    ParMultiGridCL::Instance().IdentifyEdge( ePtrs_[SubEdge(edge, 1)], ep);
                }
#endif
            }
        }
        else
        {
            ePtrs_[edge]= ep;
        }
    }

    // Collect inner edges
    byte const* tmp= std::lower_bound(refrule.Edges, refrule.Edges+refrule.EdgeNum, static_cast<byte>(NumObviousEdgesC) );
        // pointer to first inner edge entry
    for (Uint i= tmp - refrule.Edges; i<refrule.EdgeNum; ++i)
    {
        const Uint      edge = refrule.Edges[i];
        VertexCL* const vp0  = GetVertMidVert(VertOfEdge(edge, 0));
        VertexCL* const vp1  = GetVertMidVert(VertOfEdge(edge, 1));

        if ( !(ePtrs_[edge]=vp0->FindEdge(vp1)) )
        {
            EdgeCL* newEdge=0;
            if ( IsDiagonal(edge) )
                newEdge= &factory.MakeEdge( vp0, vp1, nextLevel);
            else // lies on parent face
            {
                newEdge= &factory.MakeEdge( vp0, vp1, nextLevel, Faces_[ParFaceOfEdge(edge)]->GetBndIdx());
#ifdef _PAR
                // if new faces are created on the proc-boundary, identify them
                if ( GetFace(ParFaceOfEdge(edge))->IsOnProcBnd() ) //Identify
                    ParMultiGridCL::Instance().IdentifyEdge( newEdge, GetFace(ParFaceOfEdge(edge)));
#endif
            }
            ePtrs_[edge] = newEdge;
            ePtrs_[edge]->RecycleMe();
        }
    }
}

void TetraCL::CollectFaces (const RefRuleCL& refrule, SimplexFactoryCL& factory)
/**
The faces for new refinement are stored in the static TetraCL::fPtrs array.
First look for them in the recycle bin (maybe they were created before),
if the face cannot be found, create it and link boundary, if necessary.
*/
{
    const Uint nextLevel= GetLevel()+1;

    for (Uint i=0; i<refrule.FaceNum; ++i)
    {
        const Uint face= refrule.Faces[i];

        if (IsParentFace(face))
            fPtrs_[face]= Faces_[face];
        else
        {
                  VertexCL* const vp0= GetVertMidVert(VertOfFace(face, 0));
            const VertexCL* const vp1= GetVertMidVert(VertOfFace(face, 1));
            const VertexCL* const vp2= GetVertMidVert(VertOfFace(face, 2));
            if (!(fPtrs_[face]= vp0->FindFace(vp1, vp2) ) )
            {
                FaceCL* newFace;
                if ( IsSubFace(face) )
                {
#ifndef _PAR
                    newFace= &factory.MakeFace( nextLevel, GetFace(ParentFace(face))->GetBndIdx());
#else
                    newFace= &factory.MakeFace( nextLevel, vp0->GetCoord(), vp1->GetCoord(), vp2->GetCoord(), GetFace(ParentFace(face))->GetBndIdx());
                    // if new faces are created on the proc-boundary, identify them
                    if ( GetFace(ParentFace(face))->IsOnProcBnd())  //Identify
                        ParMultiGridCL::Instance().IdentifyFace( newFace, GetFace(ParentFace(face)));
#endif
                }
                else{
#ifndef _PAR
                    newFace= &factory.MakeFace( nextLevel);
#else
                    newFace= &factory.MakeFace( nextLevel, vp0->GetCoord(), vp1->GetCoord(), vp2->GetCoord());
#endif
                }

                fPtrs_[face] = newFace;
                fPtrs_[face]->RecycleMe(vp0, vp1, vp2);
            }
        }
    }
}

void TetraCL::CollectAndLinkChildren (const RefRuleCL& refrule, SimplexFactoryCL& factory)
/**
The child tetras for new refinement are stored in the Children_ array.
First look for them in the recycle bin (maybe they are still left from the old rule),
if the child cannot be found, create it.
*/
{
    if ( !Children_ ) Children_= new SArrayCL<TetraCL*, MaxChildrenC>;
    Uint ChildNum= refrule.ChildNum;
    Uint ch;
    for (ch=0; ch < ChildNum; ++ch)
    {
        const ChildDataCL childdat= GetChildData(refrule.Children[ch]);
        VertexCL* const vp0= GetVertMidVert(childdat.Vertices[0]);
        VertexCL* const vp1= GetVertMidVert(childdat.Vertices[1]);
        VertexCL* const vp2= GetVertMidVert(childdat.Vertices[2]);
        VertexCL* const vp3= GetVertMidVert(childdat.Vertices[3]);

        if (!( (*Children_)[ch]= vp0->FindTetra(vp1, vp2, vp3) ))
        {
            (*Children_)[ch] = &factory.MakeTetra( vp0, vp1, vp2, vp3, this);
        }
        (*Children_)[ch]->LinkEdges(childdat);
        (*Children_)[ch]->LinkFaces(childdat);
    }

    for (; ch < MaxChildrenC; ++ch)
        (*Children_)[ch]= 0;
}

#ifdef _PAR

/**
*/
void TetraCL::Pack( DiST::MPIostreamCL& ostrstream) const
{
    ostrstream << GetGID()
               << RefRule_ << RefMark_;
    for ( Uint i=0; i<NumVertsC; ++i) {
        ostrstream << GetVertex(i)->GetGID();
    }
    for ( Uint i=0; i<NumEdgesC; ++i) {
        ostrstream << GetEdge(i)->GetGID();
    }
    for ( Uint i=0; i<NumFacesC; ++i) {
        ostrstream << GetFace(i)->GetGID();
    }
    ostrstream << ( GetParent()==0 ? DiST::NoGID : GetParent()->GetGID());

    Uint numChildren= std::distance(GetChildBegin(), GetChildEnd());
    ostrstream << numChildren;
    for ( const_ChildPIterator chp(GetChildBegin()); chp!=GetChildEnd(); ++chp) {
        ostrstream << (*chp)->GetGID();
    }
}

/**
*/
void TetraCL::UnPack( DiST::MPIistreamCL& istrstream)
{
    DiST::GeomIdCL tmp;

    istrstream >> gid_;
    istrstream >> RefRule_ >> RefMark_;
    for ( Uint i=0; i<NumVertsC; ++i) {
        istrstream >> tmp;
        Vertices_[i]= DiST::InfoCL::Instance().GetVertex(tmp);
    }
    for ( Uint i=0; i<NumEdgesC; ++i) {
        istrstream >> tmp;
        Edges_[i]= DiST::InfoCL::Instance().GetEdge(tmp);
    }
    for ( Uint i=0; i<NumFacesC; ++i) {
        istrstream >> tmp;
        Faces_[i]= DiST::InfoCL::Instance().GetFace(tmp);
    }
    istrstream >> tmp;
    if (DiST::InfoCL::Instance().Exists(tmp))
        Parent_= DiST::InfoCL::Instance().GetTetra(tmp);
    else
        Parent_= 0; // will be set by parent itself in DiST::TransferCL::CreateSimplex<TetraCL>(...)
    Uint numChildren, ch= 0;
    istrstream >> numChildren;
    if (numChildren>0) {
        if ( !Children_ ) Children_= new SArrayCL<TetraCL*, MaxChildrenC>;
        bool delChildren= false;
        for (; ch<numChildren; ++ch) {
            istrstream >> tmp;
            if (DiST::InfoCL::Instance().Exists(tmp))
                (*Children_)[ch]= DiST::InfoCL::Instance().GetTetra(tmp);
            else // children not on processor: delete Children_ array
                delChildren= true;
        }
        if (delChildren)
        	DeleteChildren();
        else {
			for (; ch < MaxChildrenC; ++ch)
				(*Children_)[ch]= 0;
        }
    }
}

void TetraCL::UpdateGID()
{
    if (DiST::InfoCL::Instance().Exists(this->GetGID()))
        throw DROPSErrCL("TetraCL::UpdateGID(): object is already registered");
    DiST::GeomIdCL hash(this->GetLevel(), *this);
    this->gid_ = hash;
}


/**
 * - For former HasGhost, set children array and midvertex pointers properly.
 * - For former Ghost, set parent pointer properly.
 */
void TetraCL::Merge( const TetraCL& t)
{
	if (!Children_) { // former HasGhost: set children array properly. Midvertex pointers will be set by ParMultiGridCL::AdaptMidVertex().
		Children_= new SArrayCL<TetraCL*, MaxChildrenC>;
		std::copy( t.Children_->begin(), t.Children_->end(), Children_->begin());
	} else if (!Parent_) // former ghost: set parent pointer properly
		Parent_= t.Parent_;
}

#endif


/**********************************************************************************************************
*
*    D e b u g   f u n c t i o n s   to verify the validity of a multigrid
*
**********************************************************************************************************/
void RecycleBinCL::DebugInfo(std::ostream& os) const
{
    os << "RecycleBinCL:\nrecycled Edges: ";
    for (EdgeContT::const_iterator it= Edges_.begin(), theend= Edges_.end(); it!=theend;++it)
        (*it)->DebugInfo(os);
    os << "recycled Faces: ";
    for (FaceWrapperContT::const_iterator it= Faces_.begin(), theend= Faces_.end(); it!=theend;++it)
        it->face->DebugInfo(os);
    os << "recycled Tetras: ";
    for (TetraContT::const_iterator it= Tetras_.begin(), theend= Tetras_.end(); it!=theend;++it)
        (*it)->DebugInfo(os);
    os << std::endl;
}

/// \todo: Checke die BoundarySegments!!

bool VertexCL::IsSane(std::ostream& os, const BoundaryCL& Bnd) const
/**
Check for:
<ol>
 <li> boundary descriptions map to the same coordinates</li>
 <li> recycle bin is empty</li>
 <li> removemark not set</li>
 <li> simplex known to DiST </li>
</ol>
*/
{
    bool sane= true;
    // Check, if all boundary descriptions map to the same coordinates
    if (BndVerts_) {
        for (std::vector<BndPointCL>::const_iterator bIt(BndVerts_->begin()); bIt != BndVerts_->end(); ++bIt) {
            if (dynamic_cast<const MeshBoundaryCL*>( Bnd.GetBndSeg(bIt->GetBndIdx())))
                continue; // We ignore MeshBoundaryCL as it does not have 2D-coordinates...
            double diff=(Bnd.GetBndSeg(bIt->GetBndIdx())->Map(bIt->GetCoord2D()) - GetCoord()).norm();
            if ( diff > DoubleEpsC)
            {
                sane= false;
                os << "BndSegCL description " << bIt->GetBndIdx()
                   << " does not match the coordinates. "<< "Difference is "<< diff<< " " << std::endl;
                os << "  Mapping gives: " << Bnd.GetBndSeg(bIt->GetBndIdx())->Map(bIt->GetCoord2D()) << ' ';
                os << "stored: " << GetCoord() << ' ';
            }
        }
    }
    // Check, that the refinement algorithm did not miss any RecycleBins
    if ( HasRecycleBin() ) {
        sane= false;
        os << "Clear your RecycleBin!";
    }
    // Check if the removemark is not set
    if (IsMarkedForRemovement()) {
        sane=false;
        os << "Vertex is marked for removement. ";
    }

#ifdef _PAR
    // Check if the vertex is known to a remote data list
    if ( DiST::InfoCL::Instance().GetRemoteList<VertexCL>().find( GetGID())==DiST::InfoCL::Instance().GetRemoteList<VertexCL>().end()) {
        sane= false;
        os << "Vertex is not known to a remote data list. ";
    }
#endif

    if (!sane) os << std::endl;
    return sane;
}


void VertexCL::DebugInfo(std::ostream& os) const
{
    os << "VertexCL: "<< GetGID() << '\n'
       << " level " << GetLevel() << ", coord " << GetCoord() << '\n'
       << ' ' << ( IsMarkedForRemovement() ? "is" : "is not") << " marked for removement\n";
    if ( IsOnBoundary()) {
        os << " boundary vertices: ";
        for (const_BndVertIt biter(BndVerts_->begin()); biter!=BndVerts_->end(); ++biter)
            os << biter->GetBndIdx()<< "  "
               << biter->GetCoord2D() << "    ";
    }
    else{
        os << " no boundary vertices found";
    }
    if (HasRecycleBin())
    {
        os << "\n recylce bin exists and is:";
        Bin_->DebugInfo(os);
    }
    else{
        os << "\n no recycle bin found\n";
    }
    os << "Has Unknowns " << Unknowns.Exist()
#ifdef _PAR
       << ", received any Unknowns " << Unknowns.HasUnkReceived()
#endif
       << std::endl;
#ifdef _PAR
    base::DebugInfo( os);
#endif
    os << std::endl;
}

bool EdgeCL::IsSane(std::ostream& os) const
/**
Check for:
<ol>
 <li> both vertices are on domain boundary, if the edge lies on it</li>
 <li> removemark not set
 <li> if the mark for refinement is greater than zero, a midvertex must exist
 <li> in parallel, if AccMFR is equal to MFR, the edge shouldn't be distributed</li>
 <li> correct GID </li>
 <li> simplex known to DiST </li>
</ol>
/// \todo Needs an update for the non-DDD interface
*/
{
    bool sane= true;
    // Check if the boundary information matches with the vertices
    if ( IsOnBoundary() && !(GetVertex(0)->IsOnBoundary() && GetVertex(1)->IsOnBoundary())) {
        sane= false;
        os << "One of the vertices is on no boundary even though the edge is. ";
    }
    if (sane) {
        for (const BndIdxT *it= GetBndIdxBegin(), *end= GetBndIdxEnd(); it!=end; ++it) {
            if (!is_in_if( GetVertex(0)->GetBndVertBegin(), GetVertex(0)->GetBndVertEnd(), BndPointSegEqCL(*it) ) ) {
               sane= false;
               os << "BndIdx " << *it << " is not among the boundaries of vertex 0. ";
            }
            if (!is_in_if( GetVertex(1)->GetBndVertBegin(), GetVertex(1)->GetBndVertEnd(), BndPointSegEqCL(*it) ) ) {
               sane= false;
               os << "BndIdx " << *it << " is not among the boundaries of vertex 1. ";
            }
        }
    }

    // Check if the remove mark is not set
    if (IsMarkedForRemovement()) {
        sane=false;
        os << "Edge is marked for removement\n";
    }
    // check if an edge, that is marked for refinement stores a pointer to the
    // midvertex
#ifndef _PAR
    if ( GetMFR()>0 && !GetMidVertex())
#else
    if (GetPrio()>=PrioGhost && GetAccMFR()>0 && !GetMidVertex())
#endif
    {
        sane= false;
        os << "Refined edge does not store a pointer to a midvertex. ";
    }
#ifdef _PAR
    // check, if the MFR and AccMFR is right for a not distributed edge
    if ( IsLocal() && AccMFR_!=MFR_)
    {
        sane= false;
        os << "Inconsistent MFR for undistributed edge. ";
    }

    // Check if the GID is computed correctly
    DiST::GeomIdCL gid( GetLevel(), *this);
    if ( gid!=GetGID()) {
        sane= false;
        os << "GID does not match. ";
    }

    // Check if the edge is known to a remote data list
    if ( DiST::InfoCL::Instance().GetRemoteList<EdgeCL>().find( GetGID())==DiST::InfoCL::Instance().GetRemoteList<EdgeCL>().end()) {
        sane= false;
        os << "Edge is not known to a remote data list. ";
    }
#endif
    if (!sane) os << std::endl;
    return sane;
}

void EdgeCL::DebugInfo (std::ostream& os) const
{
    os << "EdgeCL: " << GetGID() << '\n'
#ifdef _PAR
       << " level: " << GetLevel() << ", barycenter " << GetBary() << '\n'
#endif
       << " vertices: " << GetVertex(0)->GetGID() << ", " << GetVertex(1)->GetGID() << '\n'
       << ' ' << ( IsMarkedForRemovement() ? "is" : "is not") << " marked for removement\n"
       << " mark for refinement: "
#ifdef _PAR
       << "acc.= " << GetAccMFR() << ", local= "
#endif
       << GetMFR();
    if (IsRefined())
       os << ", midvertex " << GetMidVertex()->GetGID();
    os << '\n';
    if ( IsOnBoundary() ) {
        os << " boundary indices: ";
        for (const BndIdxT *it= GetBndIdxBegin(), *end= GetBndIdxEnd(); it!=end; ++it)
            os << *it << ' ';
        os << '\n';
    }
    else{
        os << " not on a boundary\n";
    }
#ifdef _PAR
    base::DebugInfo( os);
#endif
    os << std::endl;
}

bool FaceCL::IsSane(std::ostream& os) const
/**
Check for:
<ol>
 <li> both sides of the face point to a tetra, boundary segment or lies on proc boundary</li>
 <li> removemark not set</li>
 <li> check if neighbor tetras have me as face</li>
 <li> check if levels of neighbor tetras are correct
 <li> check if all three vertices lies on the same boundary as me, if the face lies on a domain boundary</li>
 <li> correct GID </li>
 <li> simplex known to DiST </li>
</ol>
*/
{
    bool sane= true;
    // Check if neighbors exist
    if ( GetNeighbor(0)==0) {
        sane= false;
        os << "No tetra is linked" << std::endl;
    }
    // check, that both sides point to a tetra or a BndSeg in my level
#ifndef _PAR
    if ((Neighbors_[1]!=0) == IsOnBoundary() )
    {
        sane= false;
        os << "A tetra/boundary is missing/superfluous. ";
    }
#else
    if ( (GetNeighbor(1)!=0) == IsOnBoundary() && !IsOnProcBnd() && GetPrio()!=PrioVGhost) {
        sane= false;
        if (IsOnBoundary() )
            os << "Two tetras share this face eventhough it's on bnd. ";
        else
            os << "Second tetra missing eventhough face is not on bnd. ";
    }
#endif
    for (int i=0; i<2; ++i)
    	if (GetNeighbor(i) && GetNeighbor(i+2) && GetNeighbor(i) != GetNeighbor(i+2)->GetParent()/*!is_in( GetNeighbor(i)->GetChildBegin(), GetNeighbor(i)->GetChildEnd(), GetNeighbor(i+2))*/) {
    		sane= false;
    		os << "Unrelated tetra and child at position " << i << ". ";
    	}
    if (IsMarkedForRemovement()) {
        sane=false;
        os << "Face is marked for removement\n";
    }

    // check, that linked neighbor tetras have me as face
    for(Uint i=0; i<NumFacesC; ++i) {
        if ( GetNeighbor(i)!=0 ) {
            if ( !is_in( GetNeighbor(i)->GetFacesBegin(), GetNeighbor(i)->GetFacesEnd(), this)) {
                sane= false;
                os << "Found linked tetra, that lacks me as face:\n";
                GetNeighbor(i)->DebugInfo(os);
            }
            if (GetNeighbor(i)->GetLevel() != GetLevel() + i/2) {
                sane= false;
                os << "Found linked tetra on wrong level:\n";
                GetNeighbor(i)->DebugInfo(os);
            }
        }
    }
    // Check boundary information
    if ( IsOnBoundary()) {
        if (!(GetVertex(0)->IsOnBoundary() && GetVertex(1)->IsOnBoundary() && GetVertex(2)->IsOnBoundary())) {
            sane= false;
            os << "One of the vertices is on no boundary even though the face is. ";
        }
        else if (GetBndIdx() != GetCommonBndSeg(GetVertex(0), GetVertex(1), GetVertex(2)) ) {
            sane= false;
            os << "BndIdx " << GetBndIdx() << " is not among the common boundaries of the face. ";
        }
    }

#ifdef _PAR
    // Check if the GID is computed correctly
    DiST::GeomIdCL gid( GetLevel(), *this);
    if ( gid!=GetGID()) {
        sane= false;
        os << "GID does not match. ";
    }

    // Check if the face is known to a remote data list
    if ( DiST::InfoCL::Instance().GetRemoteList<FaceCL>().find( GetGID())==DiST::InfoCL::Instance().GetRemoteList<FaceCL>().end()) {
        sane= false;
        os << "Face is not known to a remote data list. ";
    }
#endif

    if (!sane) os << std::endl;
    return sane;
}

void FaceCL::DebugInfo(std::ostream& os) const
{
    os << "FaceCL: " << GetGID() << '\n'
#ifdef _PAR
       << " level: " << GetLevel() << ", barycenter " << GetBary() << '\n'
#endif
       << ' ' << ( IsMarkedForRemovement() ? "is" : "is not") << " marked for removement\n";
    if ( IsOnBoundary()) {
        os << " is on boundary of index " << GetBndIdx();
    }
    else{
        os << " is not on a boundary";
    }
    os << "\n vertices: ";
    for ( Uint i=0; i<3; ++i) {
        os << GetVertex(i)->GetGID() << ' ';
    }
    os << "\n edges: ";
    for ( Uint i=0; i<3; ++i) {
        os << GetEdge(i)->GetGID() << ' ';
    }
    os << "\n neighbor tetras:";
    for ( int neigh=0; neigh<(IsOnNextLevel() ? 4 : 2); ++neigh) {
        if ( neigh==0) {
            os << "\n  o on level " << GetLevel() << ": ";
        }
        if ( neigh==2) {
            os << "\n  o on level " << GetLevel()+1 << ": ";
        }
        if ( GetNeighbor(neigh)) {
            os << GetNeighbor(neigh)->GetGID() << ' ';
        }
        else { // Neighbor tetra is not accessible ...
            if ( (neigh==1 || neigh==3) && IsOnBoundary())
                os << "boundary " << GetBndIdx() << ' ';
#ifdef _PAR
            else if ( GetPrio()!=PrioVGhost && !IsLocal())
                os << "on process " << GetNeighborProc() << ' ';
            else if ( GetPrio()==PrioVGhost && !IsLocal())
                os << "not known due to VGhost priority ";
#endif
            else
                os << "is missing ";
        }
    }
    os << '\n';

#ifdef _PAR
    base::DebugInfo( os);
#endif
    os << std::endl;
}

bool TetraCL::IsSane(std::ostream& os) const
/**
Check for:
<ol>
 <li> volume of children sums up to my volume</li>
 <li> Master exists, if I am Ghost</li>
 <li> all children are masters</li>
 <li> If a ghost exists on another process, I have no children stored <li>
 <li> ghost tetras must have children</li>
 <li> At most two copies exist </li>
 <li> Only Ma/Gh pairs allowed </li>
 <li> Pointer to parent exists (if master copy and level!=0)</li>
 <li> neighbor-connections are right </li>
 <li> whether the ordering of the vertices in each edge is induced by the ordering of the vertices in the tetra</li>
 <li> whether the vertices of opposing edges contain all four vertices of the tetra </li>
 <li> vertices, edges and faces have the right priority </li>
 <li> refinement rule matchs to refinement of edges </li>
 <li> correct GID </li>
 <li> simplex known to DiST </li>
 <li> remote data points to the right simplex </li>
</ol>
*/
{
    bool sane= true;

    // Check, if the volume of my children adds up to my volume
    if ( Children_) {
        double vol = 0.0;
        for ( const_ChildPIterator ch(GetChildBegin()); ch!=GetChildEnd(); ++ch)
            vol += (*ch)->GetVolume();
        if ( std::fabs(GetVolume() - vol)>DoubleEpsC ) {
            sane= false;
            os << "Volume of children does not sum up to parent's volume. ";
        }
#ifdef _PAR
        // if a ghost copy exists, no children must be stored on this process
        if ( HasGhost()) {
            sane= false;
            os << "Ghost exists, but children are stored locally. ";
        }
        // Check if all children are not ghost
        for ( Uint ch= 0, numCh= GetRefData().ChildNum; ch<numCh; ++ch) {
            if ( GetChild(ch)->IsGhost() ) {
                sane= false;
                os << "Child " << ch << " should be a Master. ";
            }
        }
#endif
    }

#ifdef _PAR
    // check if a master copy exists, if this is a ghost
    if ( IsGhost() && IsLocal()) {
        sane= false;
        os << "Found ghost tetra without master. ";
    }
    // Check that a ghost has children
    if ( IsGhost() && !Children_) {
        sane= false;
        os << "Ghost with no child!\n";
    }
    // Check if at most two copies are generated
    if ( GetNumDist()>2) {
        sane= false;
        os << "Too many copies exist. ";
    }
    // Check that only Ma/Gh pairs are present
    if ( GetNumDist()==2 && GetPrio()==(++GetRemoteData().GetProcListBegin())->prio) {
        sane= false;
        os << "Copies should not have same priorities. ";
    }
    // Check if parent exists for a master copy
    if ( IsMaster() && GetLevel()!=0 && !GetParent() ) {
        sane= false;
        os << "Parent is missing!";
    }
#endif

    // Check, if the neighbor-connections are right:
    // If it is another tetra, check, if the numbering
    // across the face is consistent
    // If it is a boundary-segment, check,
    // if the three vertices belong to this boundary-segment
    for ( Uint face=0; face<NumFacesC; ++face) {
        if ( IsNeighbor(face)
#ifdef _PAR
                    && !IsProcBnd(face)
#endif
           )
        {
            const TetraCL* const np = GetNeighbor(face);
            Uint nface= Faces_[face]->GetFaceNumInTetra(np);
            std::vector<Uint> vertnum0;
            std::vector<Uint> vertnum1;
            vertnum0.reserve(3);
            vertnum1.reserve(3);
            for (Uint i=0; i<3; ++i)
            {
              vertnum0.push_back(Vertices_[VertOfFace(face, i)]->GetId().GetIdent());
              vertnum1.push_back(np->Vertices_[VertOfFace(nface,i)]->GetId().GetIdent());
            }
            // Since we assume that the initial triangulation is numbered consistently
            // there is no need to sort the vertices in the face, as they should be ordered
            // in exactly the same way in both tetras;
//                    std::sort(vertnum0.begin(), vertnum0.end());
//                    std::sort(vertnum1.begin(), vertnum1.end());
            if (vertnum0 != vertnum1)
            {
                sane= false;
                os << "Neighborhood across face " << face << " is screwed or the ordering induced\n"
                   << "by the tetras is not the same, which is a BadThing (TM). ";
            }
        }
    }
    // Check, whether the ordering of the vertices in each edge is induced by the ordering of
    // the vertices in the tetra
    for (Uint edge=0; edge <NumEdgesC; ++edge)
    {
        if (   Edges_[edge]->GetVertex(0) != Vertices_[VertOfEdge(edge, 0)]
            || Edges_[edge]->GetVertex(1) != Vertices_[VertOfEdge(edge, 1)])
        {
            sane= false;
            os << "Edge " << edge << " and this tetra are not numbered consistently. ";
        }
    }
    // Check, whether the vertices of opposing edges contain all four vertices of the tetra
    std::vector<const VertexCL*> tverts;
    tverts.reserve(4);
    tverts.push_back(Vertices_[0]); tverts.push_back(Vertices_[1]);
    tverts.push_back(Vertices_[2]); tverts.push_back(Vertices_[3]);
    std::sort(tverts.begin(), tverts.end());
    for (Uint edge=0, mask= 1; edge<NumEdgesC; ++edge, mask*=2)
    {
        std::vector<const VertexCL*> everts;
        everts.reserve(4);
        everts.push_back(Edges_[edge]->GetVertex(0)); everts.push_back(Edges_[edge]->GetVertex(1));
        everts.push_back(Edges_[OppEdge(edge)]->GetVertex(0));
        everts.push_back(Edges_[OppEdge(edge)]->GetVertex(1));
        std::sort(everts.begin(), everts.end());
        if (tverts!=everts)
        {
            sane= false;
            os << "Edge " << edge << " and its opposite " << Uint(OppEdge(edge))
               << " have wrong vertices with respect to tetra " << GetId().GetIdent() << ". ";
        }
        if ( static_cast<bool>(RefRule_ & mask) != GetEdge(edge)->IsMarkedForRef() )
        {
            sane= false;
            os << "Refinement rule does not match refinement of " << (GetEdge(edge)->IsMarkedForRef() ? "" : "un") <<"refined edge " << edge << ". ";
        }
    }
#ifdef _PAR
    const Priority prio= HasGhost() ? PrioVGhost : GetPrio();
    for (Uint i= 0; i<NumVertsC; ++i) {
        if (GetVertex(i)->GetPrio() < prio) {
            sane=  false;
            os << "Vertex has wrong priority " << PriorityToString(GetVertex(i)->GetPrio()) << ". ";
        }
    }
    for (Uint i= 0; i<NumEdgesC; ++i) {
        if (GetEdge(i)->GetPrio() < prio) {
            sane=  false;
            os << "Edge has wrong priority " << PriorityToString(GetEdge(i)->GetPrio()) << ". ";
        }
    }
    for (Uint i= 0; i<NumFacesC; ++i) {
        if (GetFace(i)->GetPrio() < prio) {
            sane=  false;
            os << "Face has wrong priority " << PriorityToString(GetFace(i)->GetPrio()) << ". ";
        }
    }

    if (IsMarkedForRemovement() || (IsGhost() && IsMarkedForNoRef())) {
        sane=false;
        os << "Tetra is marked for removement\n";
    }
#endif

#ifdef _PAR
    // Check if the hash is computed correctly
    DiST::GeomIdCL gid( GetLevel(), *this);
    if ( gid!=GetGID()) {
        sane= false;
        os << "GID does not match. ";
    }

    // Check if the tetra is known to a remote data list
    if ( !DiST::InfoCL::Instance().Exists( GetGID())) {
        sane= false;
        os << "Tetra is not known to a remote data list. ";
    }
    // Check if the local object pointer in the remote data links to me
    else if (DiST::InfoCL::Instance().GetTetra( GetGID()) != this)
    {
        sane= false;
        os << "Inconsistent local object in remote data. ";
    }
#endif

    if (!sane) os << std::endl;
    return sane;
}


void TetraCL::DebugInfo (std::ostream& os) const
{
    os << "TetraCL: " << GetGID() << '\n'
#ifdef _PAR
       << " level: " << GetLevel() << ", barycenter " << GetBary() << '\n'
#endif
       << " RefRule " << GetRefRule() << ", RefMark " << GetRefMark() << '\n'
       << ' ' << ( IsMarkedForRemovement() ? "is" : "is not") << " marked for removement\n";
    os << " vertices: ";
    for ( Uint i=0; i<NumVertsC; ++i) {
        os << GetVertex(i)->GetGID() << ' ';
    }
    os << "\n edges: ";
    for ( Uint i=0; i<NumEdgesC; ++i) {
        os << GetEdge(i)->GetGID() << ' ';
    }
    os << "\n faces: ";
    for ( Uint i=0; i<NumFacesC; ++i) {
        os << GetFace(i)->GetGID() << ' ';
    }
    os << "\n parent: ";
    if ( GetParent()) {
        os << GetParent()->GetGID() << '\n';
    }
    else { // Parent is gone ...
        if ( GetLevel()==0)
            os << "no parent due to level\n";
#ifdef _PAR
        else if ( IsGhost())
            os << "parent is stored on process " << (++GetRemoteData().GetProcListBegin())->proc << '\n';
#endif
        else
            os << " the parent is vanished\n";
    }
    os << " children: ";
    if ( IsUnrefined()) {
        os << "tetra is childless";
    }
    else {
        for ( const_ChildPIterator ch(GetChildBegin()); ch!=GetChildEnd(); ++ch) {
            os << (*ch)->GetGID() << ' ';
        }
    }

    // Neighbors
    os << "\n neighbors:\n";
    if ( IsGhost())
        os << "  o not known, because tetra is ghost\n";
    else {
        for (Uint i=0; i<NumFacesC; ++i) {
            os << "  o over face " << GetFace(i)->GetGID() << ": ";
            if ( GetFace(i)->IsOnBoundary())
                os << "boundary " << GetFace(i)->GetBndIdx() << '\n';
#ifdef _PAR
            else if ( GetFace(i)->IsOnProcBnd() )
                os << "process "<<GetFace(i)->GetNeighborProc() << '\n';
#endif
            else if( IsNeighbor(i))
                os << GetNeighbor(i)->GetGID() << '\n';
            else
                os << " not there.    ";
        }
    }

#ifdef _PAR
    base::DebugInfo( os);
#endif

    os << std::endl;
}

void circumcircle(const TetraCL& t, Point3DCL& c, double& r)
{
    SMatrixCL<3,3> M;
    const Point3DCL p0= t.GetVertex(0)->GetCoord();
    const double p0_sq= p0.norm_sq();

    for (int i=0; i<3; ++i)
    {
        const Point3DCL p= t.GetVertex(i+1)->GetCoord();
        const Point3DCL p2= p - p0;
        for (Uint j=0; j<3; ++j)
            M(i, j)= p2[j];
        c[i]= 0.5*(p.norm_sq()-p0_sq);
    }
    gauss_pivot(M, c);

    if (DebugRefineEasyC
        && (std::fabs( (c-p0).norm() - (c-t.GetVertex(1)->GetCoord()).norm()) > DoubleEpsC
            || std::fabs( (c-p0).norm() - (c-t.GetVertex(2)->GetCoord()).norm()) > DoubleEpsC
            || std::fabs( (c-p0).norm() - (c-t.GetVertex(3)->GetCoord()).norm()) > DoubleEpsC) )
    {
        std::cout << "circumcircle: Didn't find circumcenter. c: " << c << '\n';
        t.DebugInfo(std::cout);
        throw DROPSErrCL();
    }
    r= (c-p0).norm();
}

void circumcircle(const TetraCL& t, Uint face, Point3DCL& c, double& r)
{
    const Point3DCL& v0= t.GetVertex( VertOfFace(face, 0) )->GetCoord();
    const Point3DCL& p1= t.GetVertex( VertOfFace(face, 1) )->GetCoord() - v0;
    const Point3DCL& p2= t.GetVertex( VertOfFace(face, 2) )->GetCoord() - v0;
    double m[2];
    const double p1_sq= p1.norm_sq();
    const double p2_sq= p2.norm_sq();
    const double p1p2= inner_prod(p1,p2);
    const double det= p1_sq*p2_sq - p1p2*p1p2;

    m[0]= 0.5*p2_sq*(p1_sq - p1p2)/det;
    m[1]= 0.5*p1_sq*(p2_sq - p1p2)/det;
    c= v0 + m[0]*p1 + m[1]*p2;

    if (DebugRefineEasyC
        && (   std::fabs( (c-v0).norm() - (c-(p1+v0)).norm()) > DoubleEpsC
            || std::fabs( (c-v0).norm() - (c-(p2+v0)).norm()) > DoubleEpsC) )
    {
        std::cout << "circumcircle: Didn't find circumcenter. c: " << c << " face: " << face << '\n';
        t.DebugInfo(std::cout);
        throw DROPSErrCL();
    }
    r= (c-v0).norm();
}


// S I M P L E X  F A C T O R Y  C L A S S
//----------------------------------------

VertexCL& SimplexFactoryCL::MakeVertex( const Point3DCL& Point3D, Uint Level, IdCL<VertexCL> id, __UNUSED__ const bool donotRegister)
{
    vertices_[Level].push_back( VertexCL(Point3D, Level, id));
    VertexCL& ret= vertices_[Level].back();
#ifdef _PAR
    if (!donotRegister)
        DiST::InfoCL::Instance().GetRemoteList<VertexCL>().Register( ret);
#endif

    return ret;
}

#ifdef _PAR
template <>
VertexCL& SimplexFactoryCL::MakeCopy<VertexCL>( const VertexCL& v, const ProcListT& pl)
{
    vertices_[v.GetLevel()].push_back( v);
    VertexCL& ret= vertices_[v.GetLevel()].back();
#ifdef _PAR
    DiST::InfoCL::Instance().GetRemoteList<VertexCL>().Register( ret, pl);
#endif
    return ret;
}
#endif

EdgeCL& SimplexFactoryCL::MakeEdge( VertexCL* vp0, VertexCL* vp1, Uint Level, BndIdxT bnd0, BndIdxT bnd1, short int MFR, __UNUSED__ const bool donotRegister)
{
    edges_[Level].push_back( EdgeCL(vp0, vp1, Level, bnd0, bnd1, MFR));
    EdgeCL& ret= edges_[Level].back();
#ifdef _PAR
    if (!donotRegister)
        DiST::InfoCL::Instance().GetRemoteList<EdgeCL>().Register( ret);
#endif
    return ret;
}

#ifdef _PAR
template <>
EdgeCL& SimplexFactoryCL::MakeCopy<EdgeCL>( const EdgeCL& e, const ProcListT& pl)
{
    edges_[e.GetLevel()].push_back( e);
    EdgeCL& ret= edges_[e.GetLevel()].back();
#ifdef _PAR
    DiST::InfoCL::Instance().GetRemoteList<EdgeCL>().Register( ret, pl);
#endif
    return ret;
}
#endif

/** Normally used by the serial version*/
FaceCL& SimplexFactoryCL::MakeFace ( Uint Level, BndIdxT bnd, __UNUSED__ const bool donotRegister)
{
    faces_[Level].push_back( FaceCL(Level, bnd));
    FaceCL& ret= faces_[Level].back();
#ifdef _PAR
    throw DROPSErrCL("SimplexFactoryCL::MakeFace: Constructor only suitable for the serial version");
#endif
    return ret;
}

/** Normally used by the parallel version*/
FaceCL& SimplexFactoryCL::MakeFace ( Uint Level, const Point3DCL& v0, const Point3DCL& v1, const Point3DCL& v2, BndIdxT bnd, __UNUSED__ const bool donotRegister)
{
    faces_[Level].push_back( FaceCL( Level, v0, v1, v2, bnd));
    FaceCL& ret= faces_[Level].back();
#ifdef _PAR
    if (!donotRegister)
        DiST::InfoCL::Instance().GetRemoteList<FaceCL>().Register( ret);
#endif
    return ret;
}

#ifdef _PAR
template <>
FaceCL& SimplexFactoryCL::MakeCopy<FaceCL>( const FaceCL& f, const ProcListT& pl)
{
    faces_[f.GetLevel()].push_back( f);
    FaceCL& ret= faces_[f.GetLevel()].back();
#ifdef _PAR
    DiST::InfoCL::Instance().GetRemoteList<FaceCL>().Register( ret, pl);
#endif
    return ret;
}
#endif

/** Normally used by the serial version*/
TetraCL& SimplexFactoryCL::MakeTetra( VertexCL* vp0, VertexCL* vp1, VertexCL* vp2, VertexCL* vp3, TetraCL* Parent, IdCL<TetraCL> id, __UNUSED__ const bool donotRegister)
{
    const Uint Level= (Parent==0 ? 0 : Parent->GetLevel()+1);
    tetras_[Level].push_back( TetraCL(vp0, vp1, vp2, vp3, Parent, id));
    TetraCL& ret= tetras_[Level].back();
#ifdef _PAR
    if (!donotRegister)
        DiST::InfoCL::Instance().GetRemoteList<TetraCL>().Register( ret);
#endif
    return ret;
}

#ifdef _PAR
/** Normally used by the parallel version*/
TetraCL& SimplexFactoryCL::MakeTetra( VertexCL* vp0, VertexCL* vp1, VertexCL* vp2, VertexCL* vp3, TetraCL* Parent, Uint Level, IdCL<TetraCL> id, __UNUSED__ const bool donotRegister)
{
    tetras_[Level].push_back( TetraCL(vp0, vp1, vp2, vp3, Parent, Level, id));
    TetraCL& ret= tetras_[Level].back();
#ifdef _PAR
    if (!donotRegister)
        DiST::InfoCL::Instance().GetRemoteList<TetraCL>().Register( ret);
#endif
    return ret;
}
#endif

#ifdef _PAR
template <>
TetraCL& SimplexFactoryCL::MakeCopy<TetraCL>( const TetraCL& t, const ProcListT& pl)
{
    tetras_[t.GetLevel()].push_back( t);
    TetraCL& ret= tetras_[t.GetLevel()].back();
#ifdef _PAR
    DiST::InfoCL::Instance().GetRemoteList<TetraCL>().Register( ret, pl);
#endif
    for ( TetraCL::const_ChildPIterator chp(ret.GetChildBegin()); chp!=ret.GetChildEnd(); ++chp)
        (*chp)->Parent_= &ret;
    return ret;
}
#endif

void SimplexFactoryCL::DestroyMarkedTetras( Uint Level)
{
    MG_TetraContT::LevelCont& tets= tetras_[Level];
    for (MultiGridCL::TetraIterator it= tets.begin(), end= tets.end(); it!=end; )
    {
        if ( it->IsMarkedForRemovement() )
        {
            it->UnlinkFromFaces();
            tets.erase(it++);
        }
        else
            ++it;
    }
}

void SimplexFactoryCL::DestroyMarkedVEFs( Uint Level)
{
    vertices_[Level].remove_if( std::mem_fun_ref( &VertexCL::IsMarkedForRemovement) );
    edges_   [Level].remove_if( std::mem_fun_ref(   &EdgeCL::IsMarkedForRemovement) );
    faces_   [Level].remove_if( std::mem_fun_ref(   &FaceCL::IsMarkedForRemovement) );
}

Point3DCL ComputeBaryCenter(const Point3DCL& v0, const Point3DCL& v1)
{
    Point3DCL ret;
    for (int i=0; i<3; ++i)
        ret[i] = 0.5*(v0[i] + v1[i]);
    return ret;
}

Point3DCL ComputeBaryCenter(const Point3DCL& v0, const Point3DCL& v1, const Point3DCL& v2)
{
    Point3DCL ret;
    for (int i=0; i<3; ++i)
        ret[i] = 1.0/3.0*(v0[i] + v1[i] + v2[i]);
    return ret;
}

Point3DCL ComputeBaryCenter(const Point3DCL& v0, const Point3DCL& v1, const Point3DCL& v2, const Point3DCL& v3)
{
    Point3DCL ret;
    for (int i=0; i<3; ++i)
        ret[i] = 0.25*(v0[i] + v1[i] + v2[i] + v3[i]);
    return ret;
}

Point3DCL ComputeBaryCenter(const VertexCL& s)
{
    return s.GetCoord();
}

Point3DCL ComputeBaryCenter(const EdgeCL& s)
{
    return ComputeBaryCenter(s.GetVertex(0)->GetCoord(), s.GetVertex(1)->GetCoord());
}

Point3DCL ComputeBaryCenter(const FaceCL& f)
{
    const TetraCL* const tp= f.GetSomeTetra();
    const Uint face= f.GetFaceNumInTetra(tp);

    return ComputeBaryCenter( tp->GetVertex( VertOfFace(face, 0))->GetCoord(),
           tp->GetVertex( VertOfFace(face, 1))->GetCoord(),
           tp->GetVertex( VertOfFace(face, 2))->GetCoord());
}

Point3DCL ComputeBaryCenter(const TetraCL& s)
{
    return ComputeBaryCenter(s.GetVertex(0)->GetCoord(), s.GetVertex(1)->GetCoord(), s.GetVertex(2)->GetCoord(), s.GetVertex(3)->GetCoord());
}


} // end of namespace DROPS
