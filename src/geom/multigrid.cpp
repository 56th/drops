/// \file multigrid.cpp
/// \brief classes that constitute the multigrid
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
#include "misc/params.h"
#include "misc/singletonmap.h"
#include "misc/problem.h"
#include "num/gauss.h"
#include <iterator>
#include <set>

namespace DROPS
{

BoundaryCL::~BoundaryCL()
{
    for (SegPtrCont::iterator It=Bnd_.begin(); It!=Bnd_.end(); ++It)
        delete *It;
}

void BoundaryCL::SetPeriodicBnd( const BndTypeCont& type, match_fun match) const
{
    if (type.size()!=GetNumBndSeg())
        throw DROPSErrCL("BoundaryCL::SetPeriodicBnd: inconsistent vector size!");
#ifdef _PAR
    for (size_t i=0; i<type.size(); ++i){
        if (type[i]!=OtherBnd){
            throw DROPSErrCL("No periodic boundary conditions implemented in the parallel version, yet");
        }
    }
#endif
    BndType_= type;
    match_= match;
}

BoundaryCL::BndType PeriodicEdgesCL::GetBndType( const EdgeCL& e) const
{
    BoundaryCL::BndType type= BoundaryCL::OtherBnd;
    for (const BndIdxT *bndIt= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); bndIt!=end; ++bndIt)
        type= std::max( type, mg_.GetBnd().GetBndType(*bndIt));
    return type;
}

void PeriodicEdgesCL::Accumulate()
{
    // initialize MFR counters on all Per1 edges
    for (iterator It( list_.begin()), End(list_.end()); It!=End; ++It)
        It->first->MFR_= It->first->localMFR_;
    // compute sum in Per1 MFR counters
    for (iterator It( list_.begin()), End(list_.end()); It!=End; ++It)
        It->first->MFR_+= It->second->localMFR_;
    // copy Per1 MFR counter to Per2 MFR counter
    for (iterator It( list_.begin()), End(list_.end()); It!=End; ++It)
        It->second->MFR_= It->first->MFR_;
}

void PeriodicEdgesCL::Recompute( EdgeIterator begin, EdgeIterator end)
{
    typedef std::list<EdgeCL*> psetT;
    psetT s1, s2;
    // collect all objects on Per1/Per2 bnds in s1, s2 resp.
    for (EdgeIterator it= begin; it!=end; ++it)
        if (it->IsOnBoundary())
        {
            BoundaryCL::BndType type= GetBndType( *it);
            if (type==BoundaryCL::Per1Bnd)
                s1.push_back( &*it);
            else if (type==BoundaryCL::Per2Bnd)
                s2.push_back( &*it);
        }
    // now we have s1.size() <= s2.size()
    // match objects in s1 and s2
    const BoundaryCL& bnd= mg_.GetBnd();
    for (psetT::iterator it1= s1.begin(), end1= s1.end(); it1!=end1; ++it1)
    {
        // search corresponding object in s2
        for (psetT::iterator it2= s2.begin(), end2= s2.end(); it2!=end2; )
            if (bnd.Matching( GetBaryCenter( **it1), GetBaryCenter( **it2)) )
            {
                // store pair in list_
                list_.push_back( IdentifiedEdgesT( *it1, *it2));
                // remove it2 from s2
                s2.erase( it2++);
            }
            else it2++;
    }
    if (!s2.empty())
        throw DROPSErrCL( "PeriodicEdgesCL::Recompute: Periodic boundaries do not match!");
}

void PeriodicEdgesCL::DebugInfo( std::ostream& os)
{
    int num= 0;
    for (PerEdgeContT::iterator it= list_.begin(), end=  list_.end(); it!=end; ++it, ++num)
    {
        it->first->DebugInfo( os);
        os << "\t\t<-- " << num << " -->\n";
        it->second->DebugInfo( os);
        os << "===================================================================\n";
    }
    os << num << " identified edges found.\n\n";
}


void PeriodicEdgesCL::Shrink()
{
    list_.clear();
}

void PeriodicEdgesCL::AccumulateMFR( int lvl)
{
    if (!mg_.GetBnd().HasPeriodicBnd()) return;
    Shrink();
    for (int i=0; i<=lvl; ++i)
        Recompute( mg_.GetEdgesBegin(i), mg_.GetEdgesEnd(i));
//std::cout << " \n>>> After Recompute:\n"; DebugInfo( std::cout);
    Accumulate();
//std::cout << " \n>>> After Accumulate:\n"; DebugInfo( std::cout);
    Shrink();
}


MultiGridCL::MultiGridCL (const MGBuilderCL& Builder)
    : TriangVertex_( *this), TriangEdge_( *this), TriangFace_( *this), TriangTetra_( *this), version_(0),
    factory_( Vertices_, Edges_, Faces_, Tetras_)
{
#ifdef _PAR
    DiST::InfoCL::Instance( this);  // tell InfoCL about the multigrid before(!) building the grid
    ParMultiGridCL::Instance().AttachTo( *this);
#endif
    Builder.build(this);
    FinalizeModify();
#ifdef _PAR
    ParMultiGridCL::Instance().MarkSimplicesForUnknowns();
#endif
}

void MultiGridCL::ClearTriangCache ()
{
    TriangVertex_.clear();
    TriangEdge_.clear();
    TriangFace_.clear();
    TriangTetra_.clear();
#ifdef _PAR
    ParMultiGridCL::Instance().MarkSimplicesForUnknowns();
#endif

    for (std::map<int, ColorClassesCL*>::iterator it= colors_.begin(), end= colors_.end(); it != end; ++it)
        delete it->second;
    colors_.clear();
}

void MultiGridCL::CloseGrid(Uint Level)
{
    Comment("Closing grid " << Level << "." << std::endl, DebugRefineEasyC);

    for (TetraIterator tIt(Tetras_[Level].begin()), tEnd(Tetras_[Level].end()); tIt!=tEnd; ++tIt)
    {
#ifdef _PAR
        //AllComment("Now closing tetra " << tIt->GetGID() << std::endl, DebugRefineHardC);
#else
        Comment("Now closing tetra " << tIt->GetId().GetIdent() << std::endl, DebugRefineHardC);
#endif
        if ( tIt->IsRegular() && !tIt->IsMarkedForRegRef() )
            tIt->Close();
    }
    Comment("Closing grid " << Level << " done." << std::endl, DebugRefineEasyC);
}

#ifdef _PAR
/// \brief Adaptor to rescue simplices on level 0. Needed by MultiGridCL::UnrefineGrid().
template<class SimplexT>
class KeepLevel0_fun {
  private:
    ParMultiGridCL& pmg_;
  public:
    KeepLevel0_fun() : pmg_(ParMultiGridCL::Instance()) {}
    void operator() (SimplexT* sp)
    {
        if (sp->GetLevel()==0) {
            sp->ClearRemoveMark();
            pmg_.Keep( sp);
        }
    }
};
#endif

void MultiGridCL::UnrefineGrid (Uint Level)
{
    Comment("Unrefining grid " << Level << "." << std::endl, DebugRefineEasyC);

    const Uint nextLevel(Level+1);
    std::for_each(Vertices_[nextLevel].begin(), Vertices_[nextLevel].end(), std::mem_fun_ref(&VertexCL::SetRemoveMark));
    std::for_each(Edges_[nextLevel].begin(),    Edges_[nextLevel].end(),    std::mem_fun_ref(&EdgeCL::SetRemoveMark));
    std::for_each(Faces_[nextLevel].begin(),    Faces_[nextLevel].end(),    std::mem_fun_ref(&FaceCL::SetRemoveMark));
#ifdef _PAR
    // mark all subsimplices on level 0 for removement. All subsimplices that
    // are still needed will be rescued within the refinement algorithm.
    if (Level==0) {
        std::for_each(Vertices_[0].begin(), Vertices_[0].end(), std::mem_fun_ref(&VertexCL::SetRemoveMark));
        std::for_each(Edges_[0].begin(),    Edges_[0].end(),    std::mem_fun_ref(&EdgeCL::SetRemoveMark));
        std::for_each(Faces_[0].begin(),    Faces_[0].end(),    std::mem_fun_ref(&FaceCL::SetRemoveMark));
    }
    ParMultiGridCL& pmg= ParMultiGridCL::Instance();
    typedef std::vector<TetraIterator> GhostContT;
    GhostContT killedGhosts;
    pmg.ModifyBegin();
#endif

    for (TetraIterator tIt(Tetras_[Level].begin()), tEnd(Tetras_[Level].end()); tIt!=tEnd; ++tIt)
    {
#ifndef _PAR
        Comment("inspecting children of tetra " << tIt->GetId().GetIdent() << "." << std::endl, DebugRefineHardC);
#else
        if ( tIt->IsGhost() ? !tIt->IsMarkedForNoRef() : !tIt->IsMarkedForRemovement() ) {
            // Maybe some sub simplices on level 0 have to be rescued, so get rid of their RemoveMarks and keep them during Modify
            std::for_each( tIt->GetVertBegin(), tIt->GetVertEnd(),   KeepLevel0_fun<VertexCL>());
            std::for_each( tIt->GetEdgesBegin(), tIt->GetEdgesEnd(), KeepLevel0_fun<EdgeCL>());
            std::for_each( tIt->GetFacesBegin(), tIt->GetFacesEnd(), KeepLevel0_fun<FaceCL>());
        }
        if (tIt->HasGhost())
            continue;
#endif
        if ( !tIt->IsUnrefined() ){
            if ( tIt->IsMarkEqRule() )
                tIt->ClearAllRemoveMarks();
            else
            {
                std::for_each(tIt->GetChildBegin(), tIt->GetChildEnd(), std::mem_fun(&TetraCL::SetRemoveMark));
                if ( !tIt->IsMarkedForNoRef() ) tIt->RecycleReusables();
#ifdef _PAR
                // if tetra is ghost and will have no children on this proc after unref, we can delete this tetra
                if ( tIt->IsGhost() && tIt->IsMarkedForNoRef())
                {
                    killedGhosts.push_back(tIt);
                    killedGhostTetra_= true;
                }
#endif
            }
        }
    }

#ifndef _PAR
    Comment("Now physically unlinking and removing superfluous tetras." << std::endl, DebugRefineEasyC);
    factory_.DestroyMarkedTetras(nextLevel);
    Comment("Now adapting midvertex pointers on level " << Level << ". " << std::endl, DebugRefineEasyC);
    for (EdgeIterator eIt(Edges_[Level].begin()), eEnd(Edges_[Level].end()); eIt!=eEnd; ++eIt)
        if ( (eIt->IsRefined() && eIt->GetMidVertex()->IsMarkedForRemovement()) )
            eIt->RemoveMidVertex();
    Comment("Now removing superfluous faces, edges and vertices." << std::endl, DebugRefineEasyC);
    factory_.DestroyMarkedVEFs(nextLevel);
#else
    // rescue subs that are owned by ghost tetras on next level
    pmg.TreatGhosts( nextLevel);
    // rescue verts in a special way, because they can be found in different levels
    pmg.RescueGhostVerts( 0);
    // if ghost tetras were killed during the refinement algorithm,
    // remove simplices on Level 0 once in the last call of UnrefineGrid()
    if (killedGhostTetra_ && Level==GetLastLevel()-1) {
        for_each_if( Faces_[0].begin(), Faces_[0].end(),
                Delete_fun<FaceCL>(), std::mem_fun_ref(&FaceCL::IsMarkedForRemovement) );
        for_each_if( Edges_[0].begin(), Edges_[0].end(),
                Delete_fun<EdgeCL>(), std::mem_fun_ref(&EdgeCL::IsMarkedForRemovement) );
        for_each_if( Vertices_[0].begin(), Vertices_[0].end(),
                Delete_fun<VertexCL>(), std::mem_fun_ref(&VertexCL::IsMarkedForRemovement) );
    }
    // delete marked simplices on next level
    for_each_if( Tetras_[nextLevel].begin(), Tetras_[nextLevel].end(),
            Delete_fun<TetraCL>(), std::mem_fun_ref(&TetraCL::IsMarkedForRemovement) );
    for_each_if( Faces_[nextLevel].begin(), Faces_[nextLevel].end(),
            Delete_fun<FaceCL>(), std::mem_fun_ref(&FaceCL::IsMarkedForRemovement) );
    for_each_if( Edges_[nextLevel].begin(), Edges_[nextLevel].end(),
            Delete_fun<EdgeCL>(), std::mem_fun_ref(&EdgeCL::IsMarkedForRemovement) );
    for_each_if( Vertices_[nextLevel].begin(), Vertices_[nextLevel].end(),
            Delete_fun<VertexCL>(), std::mem_fun_ref(&VertexCL::IsMarkedForRemovement) );
    // delete killed ghosts
    Delete_fun<TetraCL> del;
    for (GhostContT::iterator it=killedGhosts.begin(); it!=killedGhosts.end(); ++it)
        del(**it);
    // now DiST can internally delete references and information!
    pmg.ModifyEnd();
    GetSimplexFactory().DestroyMarkedTetras(nextLevel); // this is needed, for some reason, even though ModifyCL::Finalize() should do it
    // remove killed ghosts
    for (GhostContT::iterator it=killedGhosts.begin(); it!=killedGhosts.end(); ++it) {
        (*it)->UnlinkFromFaces();
        Tetras_[Level].erase(*it);
    }
    killedGhosts.clear();

    Comment("Now adapting midvertex pointers on level " << Level << ". " << std::endl, DebugRefineEasyC);
    for (EdgeIterator eIt(Edges_[Level].begin()), eEnd(Edges_[Level].end()); eIt!=eEnd; ++eIt)
        if ( (eIt->IsRefined() && !eIt->IsMarkedForRef()) )
            eIt->RemoveMidVertex();
#endif
    Comment("Unrefining grid " << Level << " done." << std::endl, DebugRefineEasyC);
}


void MultiGridCL::RefineGrid (Uint Level)
{
    Comment("Refining grid " << Level << std::endl, DebugRefineEasyC);

#ifdef _PAR
    ParMultiGridCL::Instance().IdentifyBegin();
#endif

    const Uint nextLevel(Level+1);
    if ( Level==GetLastLevel() ) AppendLevel();

    for (TetraIterator tIt(Tetras_[Level].begin()), tEnd(Tetras_[Level].end()); tIt!=tEnd; ++tIt)
    {
        if ( tIt->IsMarkEqRule() ) continue;

        tIt->SetRefRule( tIt->GetRefMark() );
        if ( tIt->IsMarkedForNoRef() )
        {
#ifndef _PAR
            Comment("refining " << tIt->GetId().GetIdent() << " with rule 0." << std::endl, DebugRefineHardC);
#else
            //AllComment("refining " << tIt->GetGID() << " with rule 0." << std::endl, DebugRefineHardC);
#endif
            if ( tIt->Children_ )
                { delete tIt->Children_; tIt->Children_=0; }
        }
        else
#ifdef _PAR
            if ( !tIt->HasGhost() ) // refinement will be done on ghost tetra!
#endif
            {
                const RefRuleCL& refrule( tIt->GetRefData() );
#ifdef _PAR
                //AllComment("refining " << tIt->GetGID() << " with rule " << tIt->GetRefRule() << "." << std::endl, DebugRefineHardC);
#else
                Comment("refining " << tIt->GetId().GetIdent() << " with rule " << tIt->GetRefRule() << "." << std::endl, DebugRefineHardC);
#endif
                tIt->CollectEdges           (refrule, factory_, Bnd_);
                tIt->CollectFaces           (refrule, factory_);
                tIt->CollectAndLinkChildren (refrule, factory_);
            }
    }
    for (Uint lvl= 0; lvl <= nextLevel; ++lvl)
        std::for_each( Vertices_[lvl].begin(), Vertices_[lvl].end(),
            std::mem_fun_ref( &VertexCL::DestroyRecycleBin));

#ifdef _PAR
    ParMultiGridCL::Instance().IdentifyEnd();
#endif
//    DiST::InfoCL::Instance().IsSane( std::cerr);
    Comment("Refinement of grid " << Level << " done." << std::endl, DebugRefineEasyC);
//    if (DROPSDebugC & DebugRefineHardC) if ( !IsSane(cdebug, Level) ) cdebug << std::endl;
}


void MultiGridCL::Refine()
{
#ifndef _PAR
    PeriodicEdgesCL perEdges( *this);
#endif
    PrepareModify();
#ifdef _PAR
    ParMultiGridCL& pmg= ParMultiGridCL::Instance();
    killedGhostTetra_= false;
    pmg.AdjustLevel();          // all procs must have the same number of levels
#endif

    const int tmpLastLevel( GetLastLevel() );


    for (int Level=tmpLastLevel; Level>=0; --Level)
    {
        RestrictMarks(Level);
#ifndef _PAR
        perEdges.AccumulateMFR( Level);
#else
        // calc marks over proc boundaries
        pmg.CommunicateRefMarks( Level );
        pmg.AccumulateMFR( Level );
#endif
        CloseGrid(Level);
    }

    for (int Level=0; Level<=tmpLastLevel; ++Level)
    {
#ifndef _PAR
        if ( Tetras_[Level].empty() ) continue;
#endif
        if (Level)
            CloseGrid(Level);
        if ( Level != tmpLastLevel )
            UnrefineGrid(Level);
        RefineGrid(Level);
    }

#ifdef _PAR
    pmg.AdaptPrioOnSubs();

    for (Uint l=0; l<GetLastLevel(); ++l)
    {
        Vertices_[l].remove_if( std::mem_fun_ref(&VertexCL::IsMarkedForRemovement) );
        Edges_[l].remove_if( std::mem_fun_ref(&EdgeCL::IsMarkedForRemovement) );
        Faces_[l].remove_if( std::mem_fun_ref(&FaceCL::IsMarkedForRemovement) );
    }
    killedGhostTetra_= false;

    while ( GetLastLevel()>0 && IsLevelEmpty(GetLastLevel()))
        RemoveLastLevel();
    pmg.AdjustLevel();
#else
    while ( Tetras_[GetLastLevel()].empty() ) RemoveLastLevel();
#endif

    FinalizeModify();
    ClearTriangCache();

    std::for_each( GetAllVertexBegin(), GetAllVertexEnd(),
        std::mem_fun_ref( &VertexCL::DestroyRecycleBin));
}


void MultiGridCL::Scale( __UNUSED__ double s)
{
#ifndef _PAR
    for (VertexIterator it= GetAllVertexBegin(), end= GetAllVertexEnd(); it!=end; ++it)
        it->Coord_*= s;
#else
    throw DROPSErrCL("MultiGridCL::Transform: Not implemented, yet. Sorry");
#endif
}

void MultiGridCL::Transform( __UNUSED__ Point3DCL (*mapping)(const Point3DCL&))
{
#ifndef _PAR
    for (VertexIterator it= GetAllVertexBegin(), end= GetAllVertexEnd(); it!=end; ++it)
        it->Coord_= mapping(it->Coord_);
#else
    throw DROPSErrCL("MultiGridCL::Transform: Not implemented, yet. Sorry");
#endif
}

class VertPtrLessCL : public std::binary_function<const VertexCL*, const VertexCL* , bool>
{
  public:
    bool operator() (const VertexCL* v0, const VertexCL* v1)
        { return v0->GetId() < v1->GetId(); }
};

#ifdef _PAR
void MultiGridCL::MakeConsistentHashes()
{
    for (FaceIterator sit = GetFacesBegin(0); sit != GetFacesEnd(0); ++sit){
        sit->UpdateGID();
        DiST::InfoCL::Instance().GetRemoteList<FaceCL>().Register( *sit);
    }
    for (TetraIterator sit = GetTetrasBegin(0); sit != GetTetrasEnd(0); ++sit){
        sit->UpdateGID();
        DiST::InfoCL::Instance().GetRemoteList<TetraCL>().Register( *sit);
    }
}
#endif

void MultiGridCL::MakeConsistentNumbering()
// Applicable only before the first call to Refine()
// Rearranges the Vertexorder in Tetras and Edges, so that it is the one induced by
// the global vertex-numbering in level 0
{
    // correct vertex-order in the edges
    std::for_each (GetEdgesBegin(0), GetEdgesEnd(0), std::mem_fun_ref(&EdgeCL::SortVertices));

    for (TetraIterator sit= GetTetrasBegin(0), theend= GetTetrasEnd(0); sit!=theend; ++sit)
    {
        VertexCL* vp[NumVertsC];
        std::copy(sit->Vertices_.begin(), sit->Vertices_.end(), vp+0);
        // correct vertex-order in tetras
        std::sort( sit->Vertices_.begin(), sit->Vertices_.end(), VertPtrLessCL() );

        // sort edge-pointers according to new vertex-order
        EdgeCL* ep[NumEdgesC];
        std::copy(sit->Edges_.begin(), sit->Edges_.end(), ep+0);
        for (Uint edge=0; edge<NumEdgesC; ++edge)
        {
            const Uint v0= std::distance( sit->GetVertBegin(),
                               std::find(sit->GetVertBegin(), sit->GetVertEnd(), ep[edge]->GetVertex(0)) );
            const Uint v1= std::distance( sit->GetVertBegin(),
                               std::find(sit->GetVertBegin(), sit->GetVertEnd(), ep[edge]->GetVertex(1)) );
            sit->Edges_[EdgeByVert(v0, v1)]= ep[edge];
        }

        // sort face-pointers according to new vertex-order
        FaceCL* fp[NumFacesC];
        std::copy(sit->Faces_.begin(), sit->Faces_.end(), fp);
        for (Uint face=0; face<NumFacesC; ++face)
        {
            const Uint v0= std::distance( sit->GetVertBegin(),
                               std::find(sit->GetVertBegin(), sit->GetVertEnd(), vp[VertOfFace(face, 0)]) );
            const Uint v1= std::distance( sit->GetVertBegin(),
                               std::find(sit->GetVertBegin(), sit->GetVertEnd(), vp[VertOfFace(face, 1)]) );
            const Uint v2= std::distance( sit->GetVertBegin(),
                               std::find(sit->GetVertBegin(), sit->GetVertEnd(), vp[VertOfFace(face, 2)]) );
            sit->Faces_[FaceByVert(v0, v1, v2)]= fp[face];
        }
    }
#ifdef _PAR
    MakeConsistentHashes();
#endif
}


void SetAllEdges (TetraCL* tp, EdgeCL* e0, EdgeCL* e1, EdgeCL* e2, EdgeCL* e3, EdgeCL* e4, EdgeCL* e5)
{
    tp->SetEdge( 0, e0);
    tp->SetEdge( 1, e1);
    tp->SetEdge( 2, e2);
    tp->SetEdge( 3, e3);
    tp->SetEdge( 4, e4);
    tp->SetEdge( 5, e5);
}

void SetAllFaces (TetraCL* tp, FaceCL* f0, FaceCL* f1, FaceCL* f2, FaceCL* f3)
{
    tp->SetFace( 0, f0);
    tp->SetFace( 1, f1);
    tp->SetFace( 2, f2);
    tp->SetFace( 3, f3);
}

bool HasMultipleBndSegs (const TetraCL& t)
{
    Uint numbnd= 0;
    for (Uint i= 0; numbnd < 2 && i < NumFacesC; ++i)
        if (t.IsBndSeg( i)) ++numbnd;

    return numbnd > 1;
}

void MultiGridCL::SplitMultiBoundaryTetras()
{
    PrepareModify();
    // Uint count= 0;

    // Note that new tetras can be appended to Tetras_[0] in this loop; it is assumed that
    // the pointers and iterators to the present tetras are not invalidated by appending
    // to Tetras_[0]. (This assumption must of course hold for the refinement algorithm to
    // work at all.) The new tetras never require further splitting.
    for (TetraLevelCont::iterator t= Tetras_[0].begin(), theend= Tetras_[0].end(); t != theend; ) {
        if (!HasMultipleBndSegs( *t)) {
            ++t;
            continue;
        }
        // ++count;

        // t: (v0 v1 v2 v3)

        // One new vertex: The barycenter b is in level 0 in the interior of \Omega; store it and its address.
        Point3DCL b( GetBaryCenter( *t));
        VertexCL* bp= &factory_.MakeVertex( b, /*first level*/ 0);

        // Four new edges: (v0 b) (v1 b) (v2 b) (b v3); they have no midvertices and no boundary descriptions
        EdgeCL* ep[4];
        ep[0]= &factory_.MakeEdge( t->Vertices_[0], bp, /*level*/ 0);
        ep[1]= &factory_.MakeEdge( t->Vertices_[1], bp, /*level*/ 0);
        ep[2]= &factory_.MakeEdge( t->Vertices_[2], bp, /*level*/ 0);
        ep[3]= &factory_.MakeEdge( bp, t->Vertices_[3], /*level*/ 0);

        // Six new faces: (v0 v1 b) (v0 v2 b) (v0 b v3) (v1 v2 b) (v1 b v3) (v2 b v3); they have no boundary descriptions
        FaceCL* fp[6];
#ifndef _PAR
        for (int i= 0; i < 6; ++i) {
            fp[i]= &factory_.MakeFace(  /*level*/ 0);
        }
#else
        VertexCL *v0= t->Vertices_[0], *v1=t->Vertices_[1], *v2=t->Vertices_[2], *v3=t->Vertices_[3];
        fp[0]= &factory_.MakeFace( 0, v0->GetCoord(), v1->GetCoord(), bp->GetCoord());
        fp[1]= &factory_.MakeFace( 0, v0->GetCoord(), v2->GetCoord(), bp->GetCoord());
        fp[2]= &factory_.MakeFace( 0, v0->GetCoord(), bp->GetCoord(), v3->GetCoord());
        fp[3]= &factory_.MakeFace( 0, v1->GetCoord(), v2->GetCoord(), bp->GetCoord());
        fp[4]= &factory_.MakeFace( 0, v1->GetCoord(), bp->GetCoord(), v3->GetCoord());
        fp[5]= &factory_.MakeFace( 0, v2->GetCoord(), bp->GetCoord(), v3->GetCoord());
#endif

        // Four new tetras: (v0 v1 v2 b) (v0 v1 b v3) (v0 v2 b v3) (v1 v2 b v3)
        TetraCL* tp[4];
        tp[0]= &factory_.MakeTetra( t->Vertices_[0], t->Vertices_[1], t->Vertices_[2], bp, /*parent*/ 0);
        tp[1]= &factory_.MakeTetra( t->Vertices_[0], t->Vertices_[1], bp, t->Vertices_[3], /*parent*/ 0);
        tp[2]= &factory_.MakeTetra( t->Vertices_[0], t->Vertices_[2], bp, t->Vertices_[3], /*parent*/ 0);
        tp[3]= &factory_.MakeTetra( t->Vertices_[1], t->Vertices_[2], bp, t->Vertices_[3], /*parent*/ 0);
        // Set the edge-pointers
        SetAllEdges( tp[0], t->Edges_[EdgeByVert(0, 1)], t->Edges_[EdgeByVert(0, 2)], t->Edges_[EdgeByVert(1, 2)], ep[0], ep[1], ep[2]);
        SetAllEdges( tp[1], t->Edges_[EdgeByVert(0, 1)], ep[0], ep[1], t->Edges_[EdgeByVert(0, 3)], t->Edges_[EdgeByVert(1, 3)], ep[3]);
        SetAllEdges( tp[2], t->Edges_[EdgeByVert(0, 2)], ep[0], ep[2], t->Edges_[EdgeByVert(0, 3)], t->Edges_[EdgeByVert(2, 3)], ep[3]);
        SetAllEdges( tp[3], t->Edges_[EdgeByVert(1, 2)], ep[1], ep[2], t->Edges_[EdgeByVert(1, 3)], t->Edges_[EdgeByVert(2, 3)], ep[3]);
        // Set the face-pointers
        SetAllFaces( tp[0], fp[3], fp[1], fp[0], t->Faces_[FaceByVert( 0, 1, 2)]);
        SetAllFaces( tp[1], fp[4], fp[2], t->Faces_[FaceByVert( 0, 1, 3)], fp[0]);
        SetAllFaces( tp[2], fp[5], fp[2], t->Faces_[FaceByVert( 0, 2, 3)], fp[1]);
        SetAllFaces( tp[3], fp[5], fp[4], t->Faces_[FaceByVert( 1, 2, 3)], fp[3]);

        // Set tetra-pointers of the new faces
        fp[0]->SetNeighbor( 0, tp[0]); fp[0]->SetNeighbor( 1, tp[1]);
        fp[1]->SetNeighbor( 0, tp[0]); fp[1]->SetNeighbor( 1, tp[2]);
        fp[2]->SetNeighbor( 0, tp[1]); fp[2]->SetNeighbor( 1, tp[2]);
        fp[3]->SetNeighbor( 0, tp[0]); fp[3]->SetNeighbor( 1, tp[3]);
        fp[4]->SetNeighbor( 0, tp[1]); fp[4]->SetNeighbor( 1, tp[3]);
        fp[5]->SetNeighbor( 0, tp[2]); fp[5]->SetNeighbor( 1, tp[3]);

        // Set tetra-pointers of the faces of t to the corresponding new tetra
        if (t->GetFace( 0)->GetNeighbor( 0) == &*t)
            t->Faces_[0]->SetNeighbor( 0, tp[3]);
        else
            t->Faces_[0]->SetNeighbor( 1, tp[3]);
        if (t->GetFace( 1)->GetNeighbor( 0) == &*t)
            t->Faces_[1]->SetNeighbor( 0, tp[2]);
        else
            t->Faces_[1]->SetNeighbor( 1, tp[2]);
        if (t->GetFace( 2)->GetNeighbor( 0) == &*t)
            t->Faces_[2]->SetNeighbor( 0, tp[1]);
        else
            t->Faces_[2]->SetNeighbor( 1, tp[1]);
        if (t->GetFace( 3)->GetNeighbor( 0) == &*t)
            t->Faces_[3]->SetNeighbor( 0, tp[0]);
        else
            t->Faces_[3]->SetNeighbor( 1, tp[0]);

        // Remove *t (now unused), increment t *before* erasing
        TetraLevelCont::iterator tmp= t;
        ++t;
#ifdef _PAR
        DiST::InfoCL::Instance().GetRemoteList<TetraCL>().Unregister(*tmp);
#endif
        Tetras_[0].erase( tmp);
    }

    FinalizeModify();
    ClearTriangCache();

    // std::cerr << "Split " << count << " tetras.\n";
}
#ifndef _PAR
class EdgeByVertLessCL : public std::binary_function<const EdgeCL*, const EdgeCL* , bool>
{
  public:
    bool operator() (const EdgeCL* e0, const EdgeCL* e1)
        { return    e0->GetVertex(1)->GetId() == e1->GetVertex(1)->GetId()
                 ?  e0->GetVertex(0)->GetId() <  e1->GetVertex(0)->GetId()
                 :  e0->GetVertex(1)->GetId() <  e1->GetVertex(1)->GetId(); }
};

class EdgeEqualCL : public std::binary_function<const EdgeCL*, const EdgeCL*, bool>
{
  public:
    bool operator() (const EdgeCL* e0, const EdgeCL* e1)
        { return    e0->GetVertex(1)->GetId() == e1->GetVertex(1)->GetId()
                 && e0->GetVertex(0)->GetId() == e1->GetVertex(0)->GetId(); }
};
#else

class EdgeByVertLessCL : public std::binary_function<const EdgeCL*, const EdgeCL* , bool>
{
  public:
    bool operator() (const EdgeCL* e0, const EdgeCL* e1)
        { return    e0->GetVertex(1)->GetGID() == e1->GetVertex(1)->GetGID()
                 ?  e0->GetVertex(0)->GetGID() <  e1->GetVertex(0)->GetGID()
                 :  e0->GetVertex(1)->GetGID() <  e1->GetVertex(1)->GetGID(); }
};


class EdgeEqualCL : public std::binary_function<const EdgeCL*, const EdgeCL*, bool>
{
  public:
    bool operator() (const EdgeCL* e0, const EdgeCL* e1)
        { return    e0->GetVertex(1)->GetGID() == e1->GetVertex(1)->GetGID()
                 && e0->GetVertex(0)->GetGID() == e1->GetVertex(0)->GetGID(); }
};
#endif

bool MultiGridCL::IsSane (std::ostream& os, int Level) const
{
    bool sane=true;

    if (Level==-1)
    {
        // Check all levels
        for (int lvl=0; lvl<=static_cast<int>(GetLastLevel()); ++lvl)
            if ( !IsSane(os, lvl) )
                sane = false;

        // Check if all stored vertices, edges and faces are needed by at least one tetra
        std::set<const VertexCL*> neededVertices;
        std::set<const EdgeCL*>   neededEdges;
        std::set<const FaceCL*>   neededFaces;
        for ( const_TetraIterator sit( GetAllTetraBegin()); sit!=GetAllTetraEnd(); ++sit){
            for ( Uint i=0; i<NumVertsC; ++i)
                neededVertices.insert(sit->GetVertex(i));
            for ( Uint i=0; i<NumEdgesC; ++i)
                neededEdges.insert(sit->GetEdge(i));
            for ( Uint i=0; i<NumFacesC; ++i)
                neededFaces.insert(sit->GetFace(i));
        }
        for (const_VertexIterator sit( GetAllVertexBegin()); sit != GetAllVertexEnd( Level); ++sit){
            if ( neededVertices.find(&*sit)==neededVertices.end()){
                sane=false;
                os << "Not needed vertex:\n";
                sit->DebugInfo(os);
            }
        }
        for (const_EdgeIterator sit( GetAllEdgeBegin()); sit != GetAllEdgeEnd( Level); ++sit)
            if ( neededEdges.find(&*sit)==neededEdges.end()){
                sane=false;
                os << "Not needed edge:\n";
                sit->DebugInfo(os);
            }
        for (const_FaceIterator sit( GetAllFaceBegin()); sit != GetAllFaceEnd( Level); ++sit)
            if ( neededFaces.find(&*sit)==neededFaces.end()){
                sane=false;
                os << "Not needed face:\n";
                sit->DebugInfo(os);
            }
    }
    else
    {
        // Check Vertices
        for (const_VertexIterator vIt( GetVerticesBegin( Level));
             vIt != GetVerticesEnd( Level); ++vIt)
        {
            if ( int(vIt->GetLevel())!=Level )
            {
                sane=false;
                os <<"Wrong Level (should be "<<Level<<") for\n";
                vIt->DebugInfo(os);
            }
            if ( !vIt->IsSane(os, Bnd_) )
            {
                sane=false;
                vIt->DebugInfo(os);
            }
        }
        // Check Edges
        for (const_EdgeIterator eIt( GetEdgesBegin( Level));
             eIt!=GetEdgesEnd( Level); ++eIt)
        {
            if ( int(eIt->GetLevel())!=Level )
            {
                sane=false;
                os <<"Wrong Level (should be "<<Level<<") for\n";
                eIt->DebugInfo(os);
            }

            if ( !eIt->IsSane(os) )
            {
                sane=false;
                eIt->DebugInfo(os);
            }
        }
        // An edge connecting two vertices should be unique in its level
        // This is memory-expensive!
        std::list<const EdgeCL*> elist;
        ref_to_ptr<const EdgeCL> conv;
        std::transform( GetEdgesBegin( Level), GetEdgesEnd( Level),
            std::back_inserter( elist), conv);
        elist.sort( EdgeByVertLessCL());
        if (std::adjacent_find( elist.begin(), elist.end(), EdgeEqualCL()) != elist.end() )
        {
            sane = false;
            os << "Found an edge more than once in level " << Level << ".\n";
            (*std::adjacent_find( elist.begin(), elist.end(), EdgeEqualCL()))->DebugInfo( os);
        }
        // Check Faces
        for (const_FaceIterator It( GetFacesBegin( Level));
             It!=GetFacesEnd( Level); ++It)
        {
            if ( It->GetLevel() != static_cast<Uint>(Level) )
            {
                sane=false;
                os <<"Wrong Level (should be "<<Level<<") for\n";
                It->DebugInfo(os);
            }
            if ( !It->IsSane(os) )
            {
                sane = false;
                It->DebugInfo(os);
            }
        }
        // Check Tetras
        for (const_TetraIterator tIt( GetTetrasBegin( Level));
             tIt!=GetTetrasEnd( Level); ++tIt)
        {
            if ( tIt->GetLevel() != static_cast<Uint>(Level) )
            {
                sane=false;
                os <<"Wrong Level (should be "<<Level<<") for\n";
                tIt->DebugInfo(os);
            }
            if ( !tIt->IsSane(os) )
            {
                sane = false;
                tIt->DebugInfo(os);
            }
        }
    }
    return sane;
}


void MultiGridCL::SizeInfo(std::ostream& os)
{
#ifndef _PAR
    size_t numVerts= GetVertices().size(),
           numEdges= GetEdges().size(),
           numFaces= GetFaces().size(),
           numTetras= GetTetras().size(),
           numTetrasRef= numTetras - std::distance( GetTriangTetraBegin(), GetTriangTetraEnd());
    os << numVerts  << " Verts, "
       << numEdges  << " Edges, "
       << numFaces  << " Faces, "
       << numTetras << " Tetras"
       << std::endl;
#else
    int  elems[5],
        *recvbuf=0;

    if (ProcCL::IamMaster())
        recvbuf = new int[5*ProcCL::Size()];

    elems[0] = GetVertices().size(); elems[1]=GetEdges().size();
    elems[2] = GetFaces().size();    elems[3]=GetTetras().size();
    elems[4] = elems[3] - std::distance( GetTriangTetraBegin(), GetTriangTetraEnd());

    ProcCL::Gather(elems, recvbuf, 5, ProcCL::Master());

    Uint numVerts=0, numEdges=0, numFaces=0, numTetras=0, numTetrasRef=0;
    if (ProcCL::IamMaster()){
        for (int i=0; i<ProcCL::Size(); ++i){
            numVerts  += recvbuf[i*5+0];
            numEdges  += recvbuf[i*5+1];
            numFaces  += recvbuf[i*5+2];
            numTetras += recvbuf[i*5+3];
            numTetrasRef += recvbuf[i*5+4];
        }

        if ( ProcCL::Size()<8){
            for (int i=0; i<ProcCL::Size(); ++i){
                os << "     On Proc "<<i<<" are: "
                   << recvbuf[i*5+0] << " Verts, "
                   << recvbuf[i*5+1] << " Edges, "
                   << recvbuf[i*5+2] << " Faces, "
                   << recvbuf[i*5+3] << " Tetras"
                   << '\n';
            }
        }
        os << "  Accumulated: "
           << numVerts << " Verts, "
           << numEdges << " Edges, "
           << numFaces << " Faces, "
           << numTetras << " Tetras"
           << std::endl;
    }
    delete[] recvbuf;
#endif
    IF_MASTER
    {
        // print out memory usage.
        // before manipulating stream, remember previous precision and format flags
        const int prec= os.precision();
        const std::ios_base::fmtflags ff= os.flags();
        // print only one digit after decimal point
        os.precision(1);
        os.setf( std::ios_base::fixed);

        size_t vMem= numVerts*sizeof(VertexCL),
               eMem= numEdges*sizeof(EdgeCL),
               fMem= numFaces*sizeof(FaceCL),
               tMem= numTetras*sizeof(TetraCL) + numTetrasRef*8*sizeof(TetraCL*),
                   // also account for Children_ arrays which are allocated for all refined tetras
               Mem= vMem + eMem + fMem + tMem;
        double MemMB= double(Mem)/1024/1024;
        os << "Memory used for geometry: " << MemMB << " MB ("
           << (double(vMem)/Mem*100) << "% verts, "
           << (double(eMem)/Mem*100) << "% edges, "
           << (double(fMem)/Mem*100) << "% faces, "
           << (double(tMem)/Mem*100) << "% tetras)\n";
        // restore precision and format flags
        os.precision(prec);
        os.flags( ff);
    }
}

void MultiGridCL::ElemInfo(std::ostream& os, int Level) const
{
    double hmax= -1, hmin= 1e99,
           rmax= -1, rmin= 1e99;
    DROPS_FOR_TRIANG_CONST_TETRA( (*this), Level, It) {
        double loc_max= -1, loc_min= 1e99;
        for (Uint i=0; i<3; ++i)
        {
            Point3DCL pi= It->GetVertex(i)->GetCoord();
            for (Uint j=i+1; j<4; ++j)
            {
                const double h= (It->GetVertex(j)->GetCoord() - pi).norm();
                if (h < loc_min) loc_min= h;
                if (h > loc_max) loc_max= h;
            }
        }
        if (loc_min < hmin) hmin= loc_min;
        if (loc_max > hmax) hmax= loc_max;
        const double ratio= loc_max/loc_min;
        if (ratio < rmin) rmin= ratio;
        if (ratio > rmax) rmax= ratio;
    }
#ifdef _PAR
    hmin = ProcCL::GlobalMin(hmin, ProcCL::Master());
    rmin = ProcCL::GlobalMin(rmin, ProcCL::Master());
    hmax = ProcCL::GlobalMax(hmax, ProcCL::Master());
    rmax = ProcCL::GlobalMax(rmax, ProcCL::Master());
#endif
    IF_MASTER
      os << hmin << " <= h <= " << hmax << '\t'
         << rmin << " <= h_max/h_min <= " << rmax << std::endl;
}

void MultiGridCL::DebugInfo(std::ostream& os, int Level) const
{
    for (const_VertexIterator sit(GetVerticesBegin(Level)), end(GetVerticesEnd(Level)); sit!=end; ++sit)
        sit->DebugInfo( os);
    for (const_EdgeIterator sit(GetEdgesBegin(Level)), end(GetEdgesEnd(Level)); sit!=end; ++sit)
        sit->DebugInfo( os);
    for (const_FaceIterator sit(GetFacesBegin(Level)), end(GetFacesEnd(Level)); sit!=end; ++sit)
        sit->DebugInfo( os);
    for (const_TetraIterator sit(GetTetrasBegin(Level)), end(GetTetrasEnd(Level)); sit!=end; ++sit)
        sit->DebugInfo( os);
}

#ifdef _PAR
/// \brief Get number of distributed objects on local processor
Uint MultiGridCL::GetNumDistributedObjects() const
/** Count vertices, edges, faces and tetrahedra, that are stored on at least two
    processors. */
{
    Uint numdistVert=0, numdistEdge=0, numdistFace=0, numdistTetra=0;
    for (const_VertexIterator sit(GetVerticesBegin()), end(GetVerticesEnd()); sit!=end; ++sit)
        if (!sit->IsLocal()) ++numdistVert;
    for (const_EdgeIterator sit(GetEdgesBegin()), end(GetEdgesEnd()); sit!=end; ++sit)
        if (!sit->IsLocal()) ++numdistEdge;
    for (const_FaceIterator sit(GetFacesBegin()), end(GetFacesEnd()); sit!=end; ++sit)
        if (!sit->IsLocal()) ++numdistFace;
    for (const_TetraIterator sit(GetTetrasBegin()), end(GetTetrasEnd()); sit!=end; ++sit)
        if (!sit->IsLocal()) ++numdistTetra;

    return numdistVert+numdistEdge+numdistFace+numdistTetra;
}

/// \brief Get number of tetrahedra of a given level
Uint MultiGridCL::GetNumTriangTetra(int Level)
{
    Uint numTetra=0;
    DROPS_FOR_TRIANG_TETRA( (*this), Level, It) {
        ++numTetra;
    }
    return numTetra;
}

/// \brief Get number of faces of a given level
Uint MultiGridCL::GetNumTriangFace(int Level)
{
    Uint numFace=0;
    DROPS_FOR_TRIANG_FACE( (*this), Level, It) {
        ++numFace;
    }
    return numFace;
}

/// \brief Get number of faces on processor boundary
Uint MultiGridCL::GetNumDistributedFaces(int Level)
{
    Uint numdistFace=0;
    DROPS_FOR_TRIANG_FACE( (*this), Level, It)
        if( It->IsOnProcBnd() )
            ++numdistFace;
    return numdistFace;
}
#endif

void
TriangFillCL<VertexCL>::fill (MultiGridCL& mg, TriangCL<VertexCL>::LevelCont& c, int lvl)
{
    for (MultiGridCL::VertexIterator it= mg.GetAllVertexBegin( lvl),
         theend= mg.GetAllVertexEnd( lvl); it != theend; ++it)
        if (it->IsInTriang( lvl)
#ifdef _PAR
            && it->Unknowns.InTriangLevel(lvl)
#endif
           )
            c.push_back( &*it);
    TriangCL<VertexCL>::LevelCont tmp= c;
    c.swap( tmp);
}

void
TriangFillCL<EdgeCL>::fill (MultiGridCL& mg, TriangCL<EdgeCL>::LevelCont& c, int lvl)
{
    for (MultiGridCL::EdgeIterator it= mg.GetAllEdgeBegin( lvl),
         theend= mg.GetAllEdgeEnd( lvl); it != theend; ++it)
        if (it->IsInTriang( lvl)
  #ifdef _PAR
            && it->Unknowns.InTriangLevel(lvl)
  #endif
           )
            c.push_back( &*it);
    TriangCL<EdgeCL>::LevelCont tmp= c;
    c.swap( tmp);
}

void
TriangFillCL<FaceCL>::fill (MultiGridCL& mg, TriangCL<FaceCL>::LevelCont& c, int lvl)
{
    for (MultiGridCL::FaceIterator it= mg.GetAllFaceBegin( lvl),
         theend= mg.GetAllFaceEnd( lvl); it != theend; ++it)
           if (it->IsInTriang( lvl)
  #ifdef _PAR
               && it->Unknowns.InTriangLevel(lvl)
  #endif
           )
            c.push_back( &*it);
    TriangCL<FaceCL>::LevelCont tmp= c;
    c.swap( tmp);
}

void
TriangFillCL<TetraCL>::fill (MultiGridCL& mg, TriangCL<TetraCL>::LevelCont& c, int lvl)
{
    for (MultiGridCL::TetraIterator it= mg.GetAllTetraBegin( lvl),
         theend= mg.GetAllTetraEnd( lvl); it != theend; ++it)
        if (it->IsInTriang( lvl)
 #ifdef _PAR
            && it->IsMaster()
 #endif
           ) c.push_back( &*it);
    TriangCL<TetraCL>::LevelCont tmp= c;
    c.swap( tmp);
}


MGBuilderCL::MGBuilderCL (Uint parnumLevel)
    : parnumLevel_( parnumLevel)
{}

void MGBuilderCL::build_par_impl(MultiGridCL* mgp) const
{
    for (Uint i= 0; i < parnumLevel_; ++i)
        AppendLevel( mgp);

    // Create boundary
    buildBoundary( mgp);
}


void
LocatorCL::LocateInTetra(LocationCL& loc, Uint trilevel, const Point3DCL&p, double tol)
// Assumes, that p lies in loc.Tetra_ and that loc.Tetra_ contains the barycentric
// coordinates of p therein. If these prerequisites are not met, this function might
// loop forever or lie to you. You have been warned!
// Searches p in the children of loc.Tetra_ up to triangulation-level trilevel
{
    const TetraCL*& t= loc.Tetra_;
    SVectorCL<4>& b= loc.Coord_;
    SMatrixCL<4,4> M;

    for (Uint lvl=t->GetLevel(); lvl<trilevel && !t->IsUnrefined(); ++lvl)
        {
            // Adjust relative tolerances on finer grid, so that the absolute
            // tolerances stay the same.
            tol *= 2;
            for (TetraCL::const_ChildPIterator it=t->GetChildBegin(), theend=t->GetChildEnd(); it!=theend; ++it)
            {
                MakeMatrix(**it, M);
                std::copy(p.begin(), p.end(), b.begin());
                b[3]= 1.;
                gauss_pivot(M, b);
                if ( InTetra(b, tol) )
                {
                    t= *it;
                    break;
                }
            }
        }
}

void
LocatorCL::Locate(LocationCL& loc, const MultiGridCL& MG, int trilevel, const Point3DCL& p, double tol)
/// \todo this only works for triangulations of polygonal domains, which resolve the geometry of the domain exactly (on level 0).
/// \todo this only works for FE-functions living on the finest level
{
    SVectorCL<4>& b= loc.Coord_;
    SMatrixCL<4,4> M;
#ifndef _PAR
    const Uint search_level=0;
#else
    const Uint search_level=MG.GetLastLevel()-1;
#endif

    for (MultiGridCL::const_TetraIterator it= MG.GetTetrasBegin(search_level), theend= MG.GetTetrasEnd(search_level); it!=theend; ++it)
    {
        MakeMatrix(*it, M);
        std::copy(p.begin(), p.end(), b.begin());
        b[3]= 1.;
        gauss_pivot(M, b);
        if ( InTetra(b, tol) )
        {
            loc.Tetra_= &*it;
            LocateInTetra(loc, MG.GetTriangTetra().StdIndex( trilevel), p, tol);
            return;
        }
    }
    loc.Tetra_= 0; std::fill(b.begin(), b.end(), 0.);
}

void MarkAll (DROPS::MultiGridCL& mg)
{
    DROPS_FOR_TRIANG_TETRA( mg, /*default-level*/-1, It)
        It->SetRegRefMark();
}


void UnMarkAll (DROPS::MultiGridCL& mg)
{
     DROPS_FOR_TRIANG_TETRA( mg, /*default-level*/-1, It)
     {
#ifdef _PAR
         if (!It->IsMaster()) std::cerr <<"Marking non-master tetra for removement!!!\n";
#endif
         It->SetRemoveMark();
     }
}


///\brief Helper to ColorClassesCL.
/// b Is a (dynamic) bitset of used colors, global_b is one such bitset for each p1-dof.
class used_colors_CL
{
  private:
    ///\brief count trailing zeros in the byte c.
    static Uint ctz (unsigned char c) {
        if (c == 0) // Short cut for 0 which could also be treated naturally by the following binary search.
            return 8;
        Uint n;
#ifdef __GNUC__
         n= __builtin_ctz( c);
#else
        n= 0;
        if ((c & 15) == 0) { // 2^4 - 1
            n+= 4;
            c>>= 4;
        }
        if ((c & 3) == 0) {  // 2^2 - 1
            n+= 2;
            c>>= 2;
        }
        if ((c & 1) == 0) {  // 2^1 - 1
            n+= 1;
            c>>= 1;
        }
#endif
        return n;
    }

    Uint size;
    unsigned char* b;

    const unsigned char* global_b;

  public:
    used_colors_CL (Uint nb, const unsigned char* global_barg)
        : size( nb), b( new unsigned char[size]), global_b( global_barg) {}
    ~used_colors_CL () {
        delete[] b;
    }

    void assign( const TetraCL& t, Uint sys) {
        size_t idx= t.GetVertex( 0)->Unknowns( sys);
        std::copy( global_b + idx*size, global_b + (idx + 1)*size, b);
        for (Uint i= 1; i < 4; ++i) {
            idx= t.GetVertex( i)->Unknowns( sys);
            for (Uint j= 0; j < size; ++j)
                b[j]|= global_b[idx*size + j];
        }
    }
    Uint get_free_color () const {
        Uint byte;
        for (byte= 0; byte < size && b[byte] == 0xff /*a byte with all bits set*/; ++byte)
            ;
        return ctz( ~b[byte] & (b[byte] + 1)) + (byte << 3u); // The argument of ctz has exactly one bit set in the position of the least significant zero bit of b[byte] (provided b != 0xff, which cannot happen).
    }
};

void ColorClassesCL::compute_color_classes (MultiGridCL& mg, Uint lvl, const BndCondCL& Bnd)
{
#   ifdef _PAR
        ParTimerCL timer;
#   else
        TimerCL timer;
#   endif
        timer.Start();

    colors_.resize( 0);

    // strip all Dirichlet bc from Bnd
    BndCondCL BndNoDirichlet= Bnd;
    BndNoDirichlet.StripDirichletBC();

    IdxDescCL p2idx( P2_FE, BndNoDirichlet);
    p2idx.CreateNumbering( lvl, mg);
    const Uint sys= p2idx.GetIdx();
    const size_t n= p2idx.NumUnknowns();
    VecDescBaseCL<VectorBaseCL<unsigned short> > numtetra( &p2idx);
    DROPS_FOR_TRIANG_TETRA( mg, lvl, it) {
        for (Uint i= 0; i < 4; ++i)
            ++numtetra.Data[it->GetVertex( i)->Unknowns( sys)];
        for (Uint i= 0; i < 6; ++i)
            ++numtetra.Data[it->GetEdge( i)->Unknowns( sys)];
    }
    Uint maxcolors= 0,
         tmp;
    // Use the inclusion-exclusion principle to count all tetras which intersect each tetra
    DROPS_FOR_TRIANG_TETRA( mg, lvl, it) {
        tmp= 4*2 - 1; // Tetras at the faces (upper bound), the tetra itself
        for (Uint i= 0; i < 4; ++i)
            tmp+= numtetra.Data[it->GetVertex( i)->Unknowns( sys)];
        for (Uint i= 0; i < 6; ++i)
            tmp-=numtetra.Data[it->GetEdge( i)->Unknowns( sys)];
        maxcolors= std::max( maxcolors, tmp);
    }
    ++maxcolors;
    const Uint nb= ((maxcolors + 7u) & ~7u) >> 3u; // number of bytes with at least maxcolors bits.
//     std::cerr << "\nmaxcolors: " << maxcolors << " nb: " << nb << " old maxcolors: " << 1 + 4*(*std::max_element( Addr( numtetra.Data), Addr( numtetra.Data) + n) - 1) << ".\n";
    numtetra.Reset();

    std::vector<unsigned char> color_bitsets( n*nb);
    colors_.resize( maxcolors);

    size_t idx;
    Uint c,    // The color in [0, maxcolors).
         byte; // The byte in which bit number c is.
    unsigned char bit_in_byte; // The bitmask with the bit corresponding to c set.
    used_colors_CL used_colors( nb, Addr( color_bitsets));
    DROPS_FOR_TRIANG_TETRA( mg, lvl, it) {
        used_colors.assign( *it, sys);
        c= used_colors.get_free_color();
        colors_[c].push_back( &*it);

        byte=         (c & ~7u) >> 3u;
        bit_in_byte=  1u << (c & 7u);
        for (Uint i= 0; i < 4; ++i) {
            idx= it->GetVertex( i)->Unknowns( sys);
            color_bitsets[idx*nb + byte] |= bit_in_byte;
        }
    }

    p2idx.DeleteNumbering( mg);

    Uint i;
//     size_t tetracheck= 0;
    // tetra sorting for better memory access pattern
    for (i= 0; i < colors_.size() && !colors_[i].empty(); ++i) {
        std::sort( colors_[i].begin(), colors_[i].end());
//         tetracheck+= colors_[i].size();
    }
    colors_.resize( i);
//     if (tetracheck != std::distance (mg.GetTriangTetraBegin( lvl), mg.GetTriangTetraEnd( lvl))) {
//         std::cerr << "tetracheck: " << tetracheck << " numtetra: " << std::distance (mg.GetTriangTetraBegin( lvl), mg.GetTriangTetraEnd( lvl)) << ".\n";
//         throw DROPSErrCL( "ColorClassesCL::my_compute_color_classes: messed up tetras.\n");
//     }

    timer.Stop();
    const double duration= timer.GetTime();
    std::cout << "ColorClassesCL::compute_color_classes: " << duration << " seconds, " << num_colors() << " colors.\n";
}

void ColorClassesCL::make_single_color_class (MultiGridCL::const_TriangTetraIteratorCL begin,
                                              MultiGridCL::const_TriangTetraIteratorCL end)
{
#   ifdef _PAR
        ParTimerCL timer;
#   else
        TimerCL timer;
#   endif
        timer.Start();

    colors_.resize( 1);
    colors_[0].clear();
    colors_[0].reserve( std::distance( begin, end));
    for (; begin != end; ++begin)
        colors_[0].push_back( &*begin);

    // tetra sorting for better memory access pattern
    std::sort( colors_[0].begin(), colors_[0].end());

    timer.Stop();
    std::cout << "ColorClassesCL::make_single_color_class: Creation of the tetra-coloring took " << timer.GetTime() << " seconds.\n";
}

const ColorClassesCL& MultiGridCL::GetColorClasses (int Level, const BndCondCL& Bnd) const
{
    if (Level < 0)
        Level+= GetNumLevel();

    if (colors_.find( Level) == colors_.end())
        colors_[Level]= new ColorClassesCL( *this, Level, Bnd);

    return *colors_[Level];
}

void read_PeriodicBoundaries (MultiGridCL& mg, const ParamCL& P)
{
    const BoundaryCL& bnd= mg.GetBnd();
    const BndIdxT num_bnd= bnd.GetNumBndSeg();

    match_fun mfun= 0;
    BoundaryCL::BndTypeCont bnd_type( num_bnd, BoundaryCL::OtherBnd);

    // Try to read and set PeriodicMatching.
    const ParamCL::ptree_type* child= 0;
    try {
        child= &P.get_child( "PeriodicMatching");
    } catch (DROPSParamErrCL e) {}
    if (child != 0)
        try {
            const std::string s= child->get_value<std::string>();
            if (s != "")
                mfun= SingletonMapCL<match_fun>::getInstance()[s];
        } catch (DROPSErrCL e) {
            std:: cerr << "read_PeriodicBoundaries: While processing 'PeriodicMatching'...\n";
            throw e;
        }

    // Read data for the boundary segments.
    for (ParamCL::ptree_const_iterator_type it= P.begin(), end= P.end(); it != end; ++it) {
        const std::string key= it->first;
        if (key == std::string( "PeriodicMatching"))
            continue;

        BndIdxT i;
        std::istringstream iss( key);
        iss >> i;
        if (!iss) // As BndIdxT is unsigned, the 2nd test is redundant.
            throw DROPSErrCL( "read_PeriodicBoundaries: Invalid boundary segment '" + key + "'.\n");

        BoundaryCL::BndType type= BoundaryCL::OtherBnd;
        try {
            const std::string s= it->second.get_value<std::string>();
            if (s == "Per1BC")
                type= BoundaryCL::Per1Bnd;
            else if (s == "Per2BC")
                type= BoundaryCL::Per2Bnd;
            else {
                if (s.empty())
                    throw DROPSParamErrCL("Specify PerBnd type without brackets '[...]'");
                else
                    throw DROPSParamErrCL("Unknown PerBnd type: "+s);
            }
        } catch (DROPSParamErrCL e) {
            std:: cerr << "read_PeriodicBoundaries: While processing key '" << key << "'...\n";
            throw e;
        }
        if (type != BoundaryCL::Per1Bnd && type != BoundaryCL::Per2Bnd && type != BoundaryCL::OtherBnd)
            throw DROPSErrCL( "read_PeriodicBoundaries: Key '" + key + "' specifies an invalid BndType.\n");
        bnd_type[i]= type;
    }

    // Enter the data to bnd
    bnd.SetPeriodicBnd( bnd_type, mfun);
}


} // end of namespace DROPS
