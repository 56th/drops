/// \file simplex.tpp
/// \brief Template and inline definitions of simplices
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Yuanjun Zhang; SC RWTH Aachen: Oliver Fortmeier, Daniel Medina Cardona.

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

namespace DROPS {

// R E C Y C L E  B I N  C L
//--------------------------

inline const EdgeCL* RecycleBinCL::FindEdge (const VertexCL* v) const
{
    for (EdgeContT::const_iterator it= Edges_.begin(), end= Edges_.end(); it!=end; ++it)
        if ( (*it)->GetVertex(1) == v ) return *it;
    return 0;
}


inline const FaceCL* RecycleBinCL::FindFace (const VertexCL* v1, const VertexCL* v2) const
{
    for(FaceWrapperContT::const_iterator it= Faces_.begin(), end= Faces_.end(); it!=end; ++it)
        if ( it->vert1 == v1 && it->vert2 == v2 ) return it->face;
    return 0;
}


inline const TetraCL* RecycleBinCL::FindTetra (const VertexCL* v1, const VertexCL* v2, const VertexCL* v3) const
{
    for(TetraContT::const_iterator it= Tetras_.begin(), end= Tetras_.end(); it!=end; ++it)
        if ( (*it)->GetVertex(1) == v1 && (*it)->GetVertex(2) == v2 && (*it)->GetVertex(3) == v3 )
            return *it;
    return 0;
}


// V E R T E X  C L
// ----------------

/** Most common constructor for vertices. Construct a vertex located at Coord which occurs
    first on level FirstLevel.
*/
VertexCL::VertexCL (const Point3DCL& Coord, Uint FirstLevel, IdCL<VertexCL> id) :
#ifndef _PAR
    Coord_(Coord), Level_(FirstLevel),
#else
    base( FirstLevel, Coord, /*dim*/ 0),
#endif
    Id_( id), BndVerts_(0), Bin_(0), RemoveMark_( false)
{}

/** Copy a vertex.
    <p><b>Danger!!! Copying simplices might corrupt the multigrid structure!!!</b></p>
*/
VertexCL::VertexCL (const VertexCL& v) :
#ifndef _PAR
     Unknowns( v.Unknowns), Coord_( v.Coord_), Level_ ( v.Level_),
#else
    base( v),
#endif
    Id_(v.Id_), BndVerts_( v.BndVerts_ ? new std::vector<BndPointCL>(*v.BndVerts_) : 0),
    Bin_(v.Bin_ ? new RecycleBinCL(*v.Bin_) : 0),
    RemoveMark_(v.RemoveMark_)
{ }

/** Normally, this constructor is only called while receiving vertices.
*/
VertexCL::VertexCL() :
#ifndef _PAR
    Coord_(0.0), Level_(0),
#else
    base(),
#endif
    Id_(), BndVerts_(0), Bin_(0), RemoveMark_(false)
{}

VertexCL::~VertexCL()
{
    delete BndVerts_;
    DestroyRecycleBin();
}

const Point3DCL& VertexCL::GetCoord() const
{
#ifndef _PAR
    return Coord_;
#else
    return base::GetBary();
#endif
}

void VertexCL::AddBnd (const BndPointCL& BndVert)
{
    if (BndVerts_)
        BndVerts_->push_back(BndVert);
    else
        BndVerts_= new std::vector<BndPointCL>(1,BndVert);
}

void VertexCL::BndSort ()
{
    std::sort( BndVerts_->begin(), BndVerts_->end(), BndPointSegLessCL());
}

bool VertexCL::HasBnd(const BndPointCL& BndVert) const
{
    if (!IsOnBoundary()) return false;
        for (const_BndVertIt it(GetBndVertBegin()), end(GetBndVertEnd()); it!=end; ++it)
            if (it->GetBndIdx()==BndVert.GetBndIdx()) return true;
    return false;
}

// E D G E  C L
// ------------

EdgeCL::EdgeCL (VertexCL* vp0, VertexCL* vp1, Uint Level, BndIdxT bnd0, BndIdxT bnd1, short int MFR) :
#ifndef _PAR
    Level_( Level),
#else
    base( Level, ComputeBaryCenter(vp0->GetCoord(), vp1->GetCoord()), /*dim*/ 1), AccMFR_(MFR),
#endif
    MidVertex_(0), MFR_(MFR), localMFR_(MFR), RemoveMark_(false)
{
    Vertices_[0]= vp0; Vertices_[1]= vp1;
    Bnd_[0]= bnd0; Bnd_[1]= bnd1;
}

/** Danger!!! Copying simplices might corrupt the multigrid structure!!! */
EdgeCL::EdgeCL (const EdgeCL& e) :
#ifndef _PAR
    Unknowns(e.Unknowns), Level_( e.Level_),
#else
    base( e), AccMFR_( e.AccMFR_),
#endif
    Vertices_( e.Vertices_), MidVertex_( e.MidVertex_), Bnd_( e.Bnd_),
    MFR_( e.MFR_), localMFR_( e.localMFR_), RemoveMark_( e.RemoveMark_)
{ }

/** Normally, used to receive an edge. */
EdgeCL::EdgeCL() :
#ifndef _PAR
    Level_(0),
#else
    base(), AccMFR_(-1),
#endif
    Vertices_(static_cast<VertexCL*>(0)), MidVertex_(0), Bnd_(NoBndC),
    MFR_(0), RemoveMark_(false)
{ }

/** The new subedges are stored in e1 and e2*/
void EdgeCL::BuildSubEdges( SimplexFactoryCL& factory, const BoundaryCL& Bnd, SArrayCL<EdgeCL*,2>& e)
{
    BuildMidVertex( factory, Bnd);
    e[0]= &factory.MakeEdge( Vertices_[0], GetMidVertex(), GetLevel()+1, Bnd_[0], Bnd_[1]);
    e[1]= &factory.MakeEdge( GetMidVertex(), Vertices_[1], GetLevel()+1, Bnd_[0], Bnd_[1]);
}


// F A C E  C L
// ------------

/** This constructor is normally used by the serial version of DROPS */
FaceCL::FaceCL ( __UNUSED__ Uint Level, __UNUSED__ BndIdxT bnd)
#ifndef _PAR
    : Level_(Level), Bnd_(bnd), RemoveMark_(false) {}
#else
    { throw DROPSErrCL("FaceCL::FaceCL: This constructor cannot be used by the parallel version");}
#endif

/** This constructor is normally used by the parallel version of DROPS */
FaceCL::FaceCL ( Uint Level, __UNUSED__ const Point3DCL& v0, __UNUSED__ const Point3DCL& v1, __UNUSED__ const Point3DCL& v2, BndIdxT bnd) :
#ifndef _PAR
    Level_(Level),
#else
    base( Level, ComputeBaryCenter(v0, v1, v2), /*dim*/ 2),
#endif
    Bnd_(bnd), RemoveMark_(false)
{ }

/** Danger!!! Copying simplices might corrupt the multigrid structure!!! */
FaceCL::FaceCL (const FaceCL& f) :
#ifndef _PAR
    Unknowns(f.Unknowns), Level_(f.Level_),
#else
    base( f.GetGID()),
#endif
    Neighbors_(f.Neighbors_), Bnd_(f.Bnd_),RemoveMark_(f.RemoveMark_)
{ }

/** Normally, this constructor is used for receiving a face. */
FaceCL::FaceCL() :
#ifndef _PAR
    Level_(0),
#else
    base(),
#endif
    Neighbors_(static_cast<const TetraCL*>(0)), Bnd_(NoBndC), RemoveMark_(false)
{ }

bool FaceCL::IsLinkedTo( const TetraCL* tp) const
{
    return GetLevel()==tp->GetLevel() ? Neighbors_[0]==tp || Neighbors_[1]==tp
                                      : Neighbors_[2]==tp || Neighbors_[3]==tp;
}

bool FaceCL::IsRefined () const
{
    const TetraCL* const tp= GetSomeTetra();
    return RefinesFace( tp->GetRefRule(), GetFaceNumInTetra(tp) );
}

bool FaceCL::HasNeighborTetra (const TetraCL* tp) const
{
    if ( IsOnBoundary() ) return false;
    if ( tp->GetLevel() == GetLevel() ) // sequential version: always true, parallel: neighbor may be stored on a different proc
        return Neighbors_[0]==tp ? Neighbors_[1] : Neighbors_[0];
    Assert( IsOnNextLevel() && tp->GetLevel() == GetLevel()+1,
        DROPSErrCL("FaceCL::HasNeighborTetra: Illegal Level."), DebugRefineEasyC);
    return Neighbors_[2]==tp ? Neighbors_[3] : Neighbors_[2];
}

const TetraCL* FaceCL::GetNeighborTetra (const TetraCL* tp) const
{
    if ( tp->GetLevel()==GetLevel() )
        return Neighbors_[0]==tp ? Neighbors_[1] : Neighbors_[0];
    Assert( IsOnNextLevel() && tp->GetLevel() == GetLevel()+1,
        DROPSErrCL("FaceCL::GetNeighborTetra: Illegal Level."), DebugRefineEasyC);
    return Neighbors_[2]==tp ? Neighbors_[3] : Neighbors_[2];
}

Uint FaceCL::GetFaceNumInTetra (const TetraCL* tp) const
{
    for (Uint face=0; face<NumFacesC; ++face)
        if (tp->GetFace(face) == this) return face;
    throw DROPSErrCL("FaceCL::GetFaceNumInTetra: I'm not face of given tetra!");
}

const VertexCL* FaceCL::GetVertex (Uint v) const
{
    const TetraCL* const tp= GetSomeTetra();
    return tp->GetVertex( VertOfFace(GetFaceNumInTetra(tp), v) );
}

const EdgeCL* FaceCL::GetEdge (Uint e) const
{
    const TetraCL* const tp= GetSomeTetra();
    return tp->GetEdge( EdgeOfFace(GetFaceNumInTetra(tp), e) );
}

const TetraCL* FaceCL::GetTetra (Uint Level, Uint side) const
{
    if ( Level == GetLevel() ) return Neighbors_[side];
    else
    {
        Assert(IsOnNextLevel(), DROPSErrCL("FaceCL::GetTetra: No such level."), DebugRefineEasyC);
        return Neighbors_[side+2];
    }
}


// T E T R A  C L
// --------------

/** Constructs a tetraedra by given vertices and a parent-tetrahedra. If no parent is given, this constructor
    assumes that the level is 0. This must not be the case in the parallel version.
*/
TetraCL::TetraCL (VertexCL* vp0, VertexCL* vp1, VertexCL* vp2, VertexCL* vp3, TetraCL* Parent, IdCL<TetraCL> id) :
#ifndef _PAR
    Level_(Parent==0 ? 0 : Parent->GetLevel()+1),
#else
    base(Parent==0 ? 0 : Parent->GetLevel()+1, 0.25*( vp0->GetCoord() + vp1->GetCoord()+ vp2->GetCoord() + vp3->GetCoord()), /*dim*/ 3),
#endif
    Id_(id), RefRule_(UnRefRuleC), RefMark_(NoRefMarkC),
    Parent_(Parent), Children_(0)
{
    Vertices_[0] = vp0; Vertices_[1] = vp1;
    Vertices_[2] = vp2; Vertices_[3] = vp3;
}

/** \todo: do we need lvl? */
#ifdef _PAR
TetraCL::TetraCL (VertexCL* vp0, VertexCL* vp1, VertexCL* vp2, VertexCL* vp3, TetraCL* Parent, __UNUSED__ Uint lvl, IdCL<TetraCL> id)
    : base( Parent==0 ? 0 : Parent->GetLevel()+1, ComputeBaryCenter( vp0->GetCoord(), vp1->GetCoord(), vp2->GetCoord(), vp3->GetCoord()), /*dim*/ 3),
      Id_(id), RefRule_(UnRefRuleC), RefMark_(NoRefMarkC),
      Parent_( Parent), Children_(0)
{
    Assert(!Parent && Parent->GetLevel()!=lvl-1, DROPSErrCL("TetraCL::TetraCL: Parent and given level does not match"), DebugRefineEasyC);
    Vertices_[0] = vp0; Vertices_[1] = vp1;
    Vertices_[2] = vp2; Vertices_[3] = vp3;
}
#endif

/** Danger!!! Copying simplices might corrupt the multigrid structure!!!*/
TetraCL::TetraCL (const TetraCL& T) :
#ifndef _PAR
    Unknowns(T.Unknowns), Level_(T.Level_),
#else
    base( T),
#endif
    Id_(T.Id_), RefRule_(T.RefRule_), RefMark_(T.RefMark_),
    Vertices_(T.Vertices_), Edges_(T.Edges_),
    Faces_(T.Faces_), Parent_(T.Parent_),
    Children_(T.Children_ ? new SArrayCL<TetraCL*,MaxChildrenC> (*T.Children_) : 0)
{ }

/** Normally used to receive tetras*/
TetraCL::TetraCL() :
#ifndef _PAR
    Level_(0),
#else
    base(),
#endif
    Id_(), RefRule_(UnRefRuleC),RefMark_(NoRefMarkC),
    Vertices_(static_cast<VertexCL*>(0)),
    Edges_(static_cast<EdgeCL*>(0)),Faces_(static_cast<FaceCL*>(0)),
    Parent_(0), Children_(0)
{ }


inline void TetraCL::CommitRegRefMark() const
{
    for (const_EdgePIterator epiter(Edges_.begin()); epiter!=Edges_.end(); ++epiter)
        (*epiter)->IncMarkForRef();
}


inline void TetraCL::UnCommitRegRefMark() const
{
    for (const_EdgePIterator epiter(Edges_.begin()); epiter!=Edges_.end(); ++epiter)
        (*epiter)->DecMarkForRef();
}

inline void TetraCL::Close()
{
    Uint mask= 1;
    Uint newrefmark= 0;

    for (const_EdgePIterator edgep( Edges_.begin() ); edgep!=Edges_.end(); ++edgep, mask<<=1)
        if ( (*edgep)->IsMarkedForRef() ) newrefmark|= mask;
    switch (newrefmark)
    {
      case RegRefMarkC:
        RefMark_= GreenRegRefMarkC;
        return;

      case NoRefMarkC:
      // if this tetra was marked for removement, it can be removed,
      // i. e. we leave the mark for removement, so that RestrictMarks
      // will catch it in the next coarser level
        if (RefMark_!=RemoveMarkC) RefMark_= newrefmark;
        return;

      default:
        RefMark_= newrefmark;
        return;
    }
}

inline void TetraCL::LinkEdges (const ChildDataCL& childdat)
{
    for (Uint edge=0; edge<NumEdgesC; ++edge)
    {
#ifdef _PAR
        Assert(!ePtrs_[childdat.Edges[edge]]->IsMarkedForRemovement(), DiST::ErrorCL("TetraCL::LinkEdges: link edge that is marked for removement", ePtrs_[childdat.Edges[edge]]->GetGID(), GetGID()), ~0);
#endif
        Edges_[edge]= ePtrs_[childdat.Edges[edge]];
    }
}

inline void TetraCL::LinkFaces (const ChildDataCL& childdat)
{
    for (Uint face=0; face<NumFacesC; ++face)
    {
        Faces_[face]= fPtrs_[childdat.Faces[face]];
        Faces_[face]->LinkTetra(this);
    }
}

#ifdef _PAR
inline bool TetraCL::HasGhost () const
/// Checks, if there is a ghost copy of this tetrahedron
/// \return false if there this tetra is a ghost, is local or this is master and the ghost copy on
///         another proc is marked for removement
{
    if ( IsGhost() || IsLocal())
        return false;
    else if ( (++GetRemoteData().GetProcListBegin())->prio==PrioKilledGhost ){
        return false;
    }
    return true;
}

#endif
}
