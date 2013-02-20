/// \file parmultigrid.cpp
/// \brief handling of a parallel multigrid
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

#include "parallel/parallel.h"
#include "parallel/parmultigrid.h"
#include <iostream>
#include <limits>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

namespace DROPS
{
/****************************************************************************
* C O N S T R U C T O R S                                                   *
*****************************************************************************
* Constructor and Destructor of the parallel multigrid                       *
* AttachTo: assign the given Multigrid to the parallel mutigrid             *
*           or tell the Multigrid about a VectorDescriber-Class             *
****************************************************************************/
/// Init the parallel stuff
ParMultiGridCL::ParMultiGridCL()
    : mg_(0), modify_(0), transfer_(0)
{
    // for Identify, we need all priorities except PrioVGhost
    priosId_.push_back(PrioNeutral); priosId_.push_back(PrioKilledGhost);
    priosId_.push_back(PrioGhost);   priosId_.push_back(PrioMaster);
}

/// \brief Assign the MultiGridCL to this ParMultiGridCL
void ParMultiGridCL::AttachTo( MultiGridCL& mg)
{
    // Store reference to the multigrid
    mg_= &mg;
}

/// \brief Marks all tetras on last level!
void ParMultiGridCL::MarkAll()
{
    for (MultiGridCL::TriangTetraIteratorCL it=mg_->GetTriangTetraBegin(); it!=mg_->GetTriangTetraEnd(); ++it)
        it->SetRegRefMark();
}

/// \brief Refine the Multigrid by calling the procedure ParRefine from the MulitGridCL
void ParMultiGridCL::Refine()
{
    mg_->Refine();
}

/// \brief This function assures, that every proc has the same number of level.
/// This is important in the refinement algorithm.
/** Make sure that the containers are modifiable*/
void ParMultiGridCL::AdjustLevel()
{
    int myLastLevel  =mg_->GetLastLevel();
    int allLastLevel = ProcCL::GlobalMax(myLastLevel);       // this procedure is from parallel.h

    // all procs should have the same number of levels!
    for (; myLastLevel<allLastLevel; ++myLastLevel)
    {
        mg_->AppendLevel();
    }
}

/****************************************************************************
* C H E C K I N G   F O R   S A N I T Y                                          *
*****************************************************************************
*   Functions for checking the sanity of the parallel structures                 *
****************************************************************************/
class ParMultiGridCL::SanityCheckCL
{
  private:
    int           level_;    ///< level to check
    std::ostream& os_;      ///< where to report errors

  public:
    SanityCheckCL( int level, std::ostream& os) : level_(level), os_(os) {}

    ///\name Handler for DiST::InterfaceCL
    //@{
    bool Gather( DiST::TransferableCL& t, DiST::SendStreamCL& send)
    {
        if ( t.GetDim()==1){        // edge
            EdgeCL* ep= 0;
            simplex_cast( t, ep);
            for ( Uint i=0; i<2; ++i)
                send << ep->GetVertex(i)->GetGID();
        }
        else if ( t.GetDim()==2){   // face
            FaceCL* fp= 0;
            simplex_cast( t, fp);
            for ( Uint i=0; i<3; ++i)
                send << fp->GetVertex( i)->GetGID();
            for ( Uint i=0; i<3; ++i)
                send << fp->GetEdge( i)->GetGID();
        }
        return true;
    }

    bool Scatter( DiST::TransferableCL& t, const size_t& numData, DiST::MPIistreamCL& recv)
    {
        DiST::GeomIdCL gid;
        bool sane= true;
        if ( t.GetDim()==1){        // edge
            EdgeCL* ep= 0;
            simplex_cast( t, ep);
            for ( size_t i=0; i<numData; ++i){
                for ( Uint i=0; i<2; ++i){
                    recv >> gid;
                    if ( gid!=ep->GetVertex( i)->GetGID()){
                        os_ << "Vertex " << i << " of edge does not match for local edge\n";
                        ep->DebugInfo( os_);
                        sane= false;
                    }
                }
            }
        }
        else if ( t.GetDim()==2){   // face
            FaceCL* fp= 0;
            simplex_cast( t, fp);
            for ( size_t i=0; i<numData; ++i){
                for ( Uint i=0; i<3; ++i){
                    recv >> gid;
                    if ( gid!=fp->GetVertex( i)->GetGID()){
                        os_ << "Vertex " << i << " of face does not match for local face\n";
                        fp->DebugInfo( os_);
                        sane= false;
                    }
                }
                for ( Uint i=0; i<3; ++i){
                    recv >> gid;
                    if ( gid!=fp->GetEdge( i)->GetGID()){
                        os_ << "Edge " << i << " of face does not match for local face\n";
                        fp->DebugInfo( os_);
                        sane= false;
                    }
                }
            }
        }
        return sane;
    }

    bool Call()
    // do interface comm
    {
        DiST::InterfaceCL::DimListT dimlist;
        dimlist.push_back( 1); dimlist.push_back( 2);

        const DiST::PrioListT   allPrios;
        DiST::LevelListCL Lvls; Lvls.push_back( level_);
        // communicate over all objects
        DiST::InterfaceCL comm( Lvls, allPrios, allPrios, dimlist);
        return comm.Communicate( *this);
    }
    //@}
};

/// \brief Check for sanity on a given level or on all levels (\a Level = -1)
/** For a given level \a Level,
 *  - check if the edges and faces have the same subsimplices (check with GID).
    - check multigrid for sanity.
    - check local MFR counter on edges.
    For \a Level = -1, check above items on each level, and additionally
    - check that all procs have the same number of levels.
    - check that the number of simplicies in the multigrid and those registered by DiST are the same.
    */

bool ParMultiGridCL::IsSane(std::ostream& os, int Level) const
{
    bool sane= true;

    // Checking all levels
    if (Level==-1)
    {
        Comment("Checking number of levels and simplices in multigrid and DiST" << std::endl,DebugParallelC);
        Uint maxLevel= ProcCL::GlobalMax( mg_->GetLastLevel());
        if (maxLevel != mg_->GetLastLevel() )
        {
            sane= false;
            os << "Local MultiGrid has too few levels: " << mg_->GetLastLevel()
                    << " instead of " << maxLevel <<std::endl;
        }
        // check that number of simplices in the multigrid and those registered by DiST are the same.
        DiST::InfoCL& info= DiST::InfoCL::Instance();
        std::vector<size_t> numMG(4);
        std::string simplex[4]= {"vert", "edge", "face", "tetra"};
        numMG[0]= mg_->GetVertices().size();
        numMG[1]= mg_->GetEdges().size();
        numMG[2]= mg_->GetFaces().size();
        numMG[3]= mg_->GetTetras().size();
        for (int dim=0; dim<4; ++dim)
            if (numMG[dim] != info.GetRemoteList(dim).size())
            {
                sane= false;
                os << "Inconsistent number of " << simplex[dim] << "s:\t" << numMG[dim] << " in multigrid != "
                   << info.GetRemoteList(dim).size() << " registered by DiST" << std::endl;
            }
        // make checks on each level
        for (Uint lvl=0; lvl<=maxLevel; ++lvl)
        {
            sane= ParMultiGridCL::IsSane( os, lvl) && sane;
        }
    }
    // checking on level
    else
    {
        Comment("Checking level " << Level << std::endl,DebugParallelC);
        sane= mg_->IsSane( os, Level);
        Comment("Checking inter-processor dependencies on level " << Level << std::endl, DebugParallelC);
        SanityCheckCL sanitycheck( Level, os);
        sane= sanitycheck.Call() && sane;
        Comment("Checking edge MFR counters on level " << Level << std::endl,DebugParallelC);
        sane= CheckMFR( Level, os) && sane;
    }
    return sane;
}

/// \brief Check local MFR counters of edges on a given level or on all levels (\a Level = -1)
bool ParMultiGridCL::CheckMFR( int Level, std::ostream& os) const
{
    bool sane= true;
    if (Level==-1) { // all levels
        for (Uint lvl=0; lvl<=mg_->GetLastLevel();++lvl)
            sane= ParMultiGridCL::CheckMFR( lvl, os) && sane;
    }
    else {
        DROPS_STD_UNORDERED_MAP<DiST::GeomIdCL,short int,DiST::Hashing> edgeMFR;

        for (MultiGridCL::const_TetraIterator it= mg_->GetTetrasBegin(Level), end= mg_->GetTetrasEnd(); it!=end; ++it)
            if (it->IsRegularlyRef() && it->IsMaster()) {
                // has marked its edges' MFR, so do the same here
                for (int e=0; e<6; ++e)
                    edgeMFR[ it->GetEdge(e)->GetGID()]++;
            }
        // now compare with edges' local MFRs
        for (MultiGridCL::const_EdgeIterator it= mg_->GetEdgesBegin(Level), end= mg_->GetEdgesEnd(); it!=end; ++it) {
            const DiST::GeomIdCL gid= it->GetGID();
            if (it->GetMFR() != edgeMFR[gid]) {
                sane= false;
                os << "Wrong MFR on edge " << gid << ": local MFR = " << it->GetMFR() << ", but should be " << edgeMFR[gid] << std::endl;
                it->DebugInfo(os);
            }

        }
    }
    return sane;
}

/****************************************************************************
* F U N C T I O N S   F O R   I N T E R F A C E S                           *
*****************************************************************************
* The following functions are handlers to call the DDD-Interfaces correctly *
****************************************************************************/

/// \brief accumulate the marks for refinement on Edges on Level (-1==all levels)
/** Provides Gather/Scatter for DiST::InterfaceCL
*/
class ParMultiGridCL::HandlerAccMFRCL
{
  private:
    int level_;

  public:
    HandlerAccMFRCL( int level) : level_(level) {}
    /// \brief Set AccMFR to zero and put my MFR into the message
    bool Gather( DiST::TransferableCL& t, DiST::SendStreamCL& send)
    {
        EdgeCL* ep; simplex_cast( t, ep);
        if ( !ep->IsMarkedForRemovement()){
            ep->SetAccMFR( 0);
            send << ep->GetMFR();
        }
        else
            return false;
        return true;

    }
    /// \brief Add received MFR to the accumulated MFR
    bool Scatter( DiST::TransferableCL& t, const size_t& numData, DiST::MPIistreamCL& recv)
    {
        EdgeCL* ep; simplex_cast( t, ep);
        ep->AccMFR_= 0;
        short int recvMFR=-1;
        for ( size_t i=0; i<numData; ++i){
            recv >> recvMFR;
            ep->AccMFR_+= recvMFR;
        }
        return true;
    }
    /// \brief Actually perform the communication
    void Call()
    {
        const DiST::PrioListCL allPrios;
        DiST::LevelListCL Lvls;
        if ( level_!=-1)
            Lvls.push_back( level_);
        DiST::InterfaceCL::DimListT dimlist; dimlist.push_back( DiST::GetDim<EdgeCL>());
        DiST::InterfaceCL interf( Lvls, allPrios, allPrios, dimlist);
        if ( !interf.Communicate( *this))
            throw DROPSErrCL("HandlerAccMFRCL::Call: Interface communication is broken");
    }
};

/// \brief accumulate the marks for refinement on Edges on Level (-1==all levels)
void ParMultiGridCL::AccumulateMFR( int Level)
{
    HandlerAccMFRCL handlerMFR( Level);
    handlerMFR.Call();
}

class ParMultiGridCL::HandlerRefMarkCL
{
private:
    int level_;
public:
    HandlerRefMarkCL( int level) : level_(level) {}
    /// \brief For tetra on proc bnd, only the ghost copy called the function TetraCL::RestrictMark. Therefore,
    ///     send the mark to the master copy.
    bool Gather( DiST::TransferableCL& t, DiST::SendStreamCL& send)
    {
        TetraCL* tp; simplex_cast( t, tp);

        Assert( t.GetPrio()==PrioGhost, DROPSErrCL("HandlerRefMarkCL::Gather: called for non-ghost"), DebugParallelC);
        send << tp->IsMarkedForRegRef();
        return true;
    }
    /// \brief This is called by the master copy. The corresponding ghost copy definitely put in the
    ///     message, whether the tetrahedron is marked for regular refinement
    bool Scatter( DiST::TransferableCL& t, __UNUSED__ const size_t& numData, DiST::MPIistreamCL& recv)
    {
        TetraCL* tp; simplex_cast( t, tp);

        Assert( t.GetPrio()==PrioMaster, DROPSErrCL("HandlerRefMarkCL::Gather: called for non-master"), DebugParallelC);
        Assert( numData==1, DROPSErrCL("HandlerRefMarkCL::Scatter: Received data from more than one ghost"), DebugParallelC);
        Assert( !tp->IsUnrefined(), DROPSErrCL("HandlerRefMarkCL::Scatter: Master has Ghost, eventhough, tetrahedron is unrefined!"), DebugParallelC);

        bool MarkForRegRef=false;
        recv >> MarkForRegRef;
        // Do the work that cannot be done in TetraCL::RestrictMark
        if ( tp->IsRegularlyRef()){
            if ( !MarkForRegRef){
                tp->SetNoRefMark();
                tp->UnCommitRegRefMark();
            }
        }
        else{ // tetra is irregularly refined
            if (MarkForRegRef){
                tp->SetRegRefMark();
                tp->CommitRegRefMark();
            }
            else{
                tp->SetNoRefMark();
            }
        }
        return true;
    }
    /// \brief Actually perform the communication
    void Call()
    {
        DiST::PrioListCL prioGhost; prioGhost.push_back( PrioGhost);
        DiST::PrioListCL prioMaster; prioMaster.push_back( PrioMaster);
        DiST::InterfaceCL::DimListT dimlist; dimlist.push_back( DiST::GetDim<TetraCL>());
        DiST::LevelListCL Lvls;
        if ( level_!=-1)
            Lvls.push_back( level_);
        DiST::InterfaceCL interf( Lvls, prioGhost, prioMaster, dimlist);
        if ( !interf.Communicate( *this))
            throw DROPSErrCL("HandlerRefMarkCL::Call: Interface communication is broken");
    }
};

/// \brief Communicate the Refinement Marks from all procs
/** Interface-function (TetraIF) calls Gather- and ScatterTetraRestrictMarks */
void ParMultiGridCL::CommunicateRefMarks( Uint Level)
{
	HandlerRefMarkCL handler( (int)Level);
    handler.Call();
}

class ParMultiGridCL::RescueGhostVertsCL
{
private:
    Uint level_;
    Uint lastLevel_;

public:
    RescueGhostVertsCL( Uint level, Uint lastLevel) : level_(level), lastLevel_(lastLevel) {}
    bool operator() ( DiST::TransferableCL& t)
    {
        TetraCL* tp; simplex_cast( t, tp);
        if ( !tp->IsGhost())
            return true;
        if ( !tp->IsMarkedForNoRef())
            std::for_each( tp->GetVertBegin(), tp->GetVertEnd(), std::mem_fun( &VertexCL::ClearRemoveMark) );
        return true;
    }
    void Call()
    {
        DiST::PrioListCL prioGhost; prioGhost.push_back( PrioGhost);
        DiST::InterfaceCL::DimListT dimlist; dimlist.push_back( DiST::GetDim<TetraCL>());
        // Vertices on  <level> can only be owned by ghost tetras on  <Level> to <LastLevel-1>.
        DiST::LevelListCL Lvls;
        for ( Uint lvl= level_; lvl<lastLevel_; ++lvl)
            Lvls.push_back( lvl);
        DiST::InterfaceCL interf( Lvls, prioGhost, prioGhost, dimlist);
        interf.ExecuteLocal( *this);
    }
};

/// \brief Treat Ghost-Vertices, so they are not deleted
/** Vertices have to be treated in a special way for removement, because they can be found in other
    levels than the ghost tetra (which is deleted).
*/
void ParMultiGridCL::RescueGhostVerts( Uint Level)
{
    RescueGhostVertsCL rescueGhostVerts( Level, mg_->GetLastLevel());
    rescueGhostVerts.Call();
}


class ParMultiGridCL::TreatHasGhostsCL
{
protected:
    int level_;
    DiST::ModifyCL& modify_;

    void ChangePrioForAllSubs( TetraCL& t, Priority prio)
    {
//        std::cout << "[" << ProcCL::MyRank() << "]: Change Prio of Tetra " << t.GetGID() << " from "
//            << PriorityToString(t.GetPrio()) << " to " << PriorityToString( prio) << std::endl;
        for (TetraCL::const_VertexPIterator it= t.GetVertBegin(), end= t.GetVertEnd(); it!=end; ++it){
//            std::cout << "[" << ProcCL::MyRank() << "]:   - for vertex " << (*it)->GetGID()
//                      << " from " << PriorityToString((*it)->GetPrio()) << " to "
//                      << PriorityToString( prio) << std::endl;
            modify_.ChangePrio( **it, prio);
        }
        for (TetraCL::const_EdgePIterator it= t.GetEdgesBegin(), end= t.GetEdgesEnd(); it!=end; ++it){
//            std::cout << "[" << ProcCL::MyRank() << "]:   - for edge " << (*it)->GetGID()
//                      << " from " << PriorityToString((*it)->GetPrio()) << " to "
//                      << PriorityToString( prio) << std::endl;
            modify_.ChangePrio( **it, prio);
        }
        for (TetraCL::const_FacePIterator it= t.GetFacesBegin(), end= t.GetFacesEnd(); it!=end; ++it){
//            std::cout << "[" << ProcCL::MyRank() << "]:   - for face " << (*it)->GetGID()
//                      << " from " << PriorityToString((*it)->GetPrio()) << " to "
//                      << PriorityToString( prio) << std::endl;
            modify_.ChangePrio( **it, prio);
        }
    }

    void ClearRemoveMarkForAllSubs( TetraCL& t)
    {
        std::for_each( t.GetVertBegin(),  t.GetVertEnd(), std::mem_fun( &VertexCL::ClearRemoveMark) );
        std::for_each( t.GetEdgesBegin(), t.GetEdgesEnd(), std::mem_fun( &EdgeCL::ClearRemoveMark) );
        std::for_each( t.GetFacesBegin(), t.GetFacesEnd(), std::mem_fun( &FaceCL::ClearRemoveMark) );
    }

public:
    TreatHasGhostsCL( int level, DiST::ModifyCL& modify) : level_(level), modify_(modify) {}
    bool operator() ( DiST::TransferableCL& t)
    {
        TetraCL* tp; simplex_cast( t, tp);
        Assert( !tp->IsGhost(), DROPSErrCL("TreatHasGhostsCL::operator(): Called for a ghost tetrahedron"), DebugParallelC);
        // Check if this tetrahedron needs to be handled
        if ( !tp->HasGhost()) return true;                          // there is no ghost copy, so everything is already done
        if ( tp->IsGhost() && tp->IsMarkedForNoRef()) return true;  // this tetrahedron will be deleted anyway, so ignore it

        // this tetrahedron won't be deleted, so rescue the sub-simplices ...
        ClearRemoveMarkForAllSubs( *tp);
        // and set their priority to PrioVGhost
        ChangePrioForAllSubs( *tp, PrioVGhost);
        // Delete all children, these are stored by the ghost copy
        tp->DeleteChildren();
        return true;
    }
    void Call()
    {
        DiST::PrioListCL prioMaster; prioMaster.push_back( PrioMaster);
        DiST::InterfaceCL::DimListT dimlist; dimlist.push_back( DiST::GetDim<TetraCL>());
        DiST::LevelListCL Lvls;
        if ( level_!=-1)
            Lvls.push_back( level_);
        DiST::InterfaceCL interf( Lvls, prioMaster, prioMaster, dimlist, false);
        interf.ExecuteLocal( *this);
    }
};

/// \brief Treat subsimplices that of HasGhosts
/** Mark subsimplices of HasGhosts as VGhost and rescue them <p>
    Interface-function (TetraIF) calls ExecHasGhostC*/
void ParMultiGridCL::TreatHasGhosts( int Level)
{
    TreatHasGhostsCL treatHasGhost( Level, *modify_);
    treatHasGhost.Call();
}

class ParMultiGridCL::RescueGhostsCL : public ParMultiGridCL::TreatHasGhostsCL
{
public:
    typedef ParMultiGridCL::TreatHasGhostsCL base;
    RescueGhostsCL( int level, DiST::ModifyCL& modify) : base( level, modify) {}
    /// \brief This is called for ghost tetrahedra. If the tetrahedron won't be deleted,
    ///     rescue all the sub-simplices and set their priority to PrioGhost
    bool operator() ( DiST::TransferableCL& t)
    {
        TetraCL* tp; simplex_cast( t, tp);
        Assert( tp->IsGhost(), DROPSErrCL("RescueGhostsCL::operator(): Called for a non-ghost tetrahedron"), DebugParallelC);
        if (!tp->IsMarkedForNoRef()){
            ChangePrioForAllSubs( *tp, PrioGhost);
            ClearRemoveMarkForAllSubs( *tp);
        }
        // remove link to parent
        tp->Parent_=0;
        return true;
    }
    void Call()
    {
        DiST::PrioListCL prioGhost; prioGhost.push_back( PrioGhost);
        DiST::InterfaceCL::DimListT dimlist; dimlist.push_back( DiST::GetDim<TetraCL>());
        DiST::LevelListCL Lvls;
        if ( level_!=-1)
            Lvls.push_back( level_);
        DiST::InterfaceCL interf( Lvls, prioGhost, prioGhost, dimlist);
        interf.ExecuteLocal( *this);
    }
};

/// \brief Treat Ghosts, so they are not deleted
/** Removes the del-mark on all ghost on level k+1 without NoRefMark, so they are not deleted.
    This proecude markes them also as Ghost <p>
    Interface-function (TetraIF) calls ExecGhostRescue */
void ParMultiGridCL::TreatGhosts( int Level)
{
    RescueGhostsCL rescueGhostSubs( Level, *modify_);
    rescueGhostSubs.Call();
}


class ParMultiGridCL::RescueMasterCL : public ParMultiGridCL::TreatHasGhostsCL
{
public:
    typedef ParMultiGridCL::TreatHasGhostsCL base;
    RescueMasterCL( int level, DiST::ModifyCL& modify)
        : base( level, modify) {}
    /// \brief This is called for master tetrahedra. If the tetrahedron won't be deleted or is handled by a TreatHasGhostsCL,
    ///     rescue all the sub-simplices and set their priority to PrioGhost
    bool operator() ( DiST::TransferableCL& t)
    {
        TetraCL* tp; simplex_cast( t, tp);
        Assert( !tp->IsGhost(), DROPSErrCL("RescueMasterCL::operator(): Called for a ghost tetrahedron"), DebugParallelC);
        if ( !(tp->IsMarkedForRemovement() || tp->HasGhost())){
            ChangePrioForAllSubs( *tp, PrioMaster);
            ClearRemoveMarkForAllSubs( *tp);
        }
        return true;
    }
    void Call()
    {
        DiST::PrioListCL prioMaster; prioMaster.push_back( PrioMaster);
        DiST::InterfaceCL::DimListT dimlist; dimlist.push_back( DiST::GetDim<TetraCL>());
        const DiST::LevelListCL AllLvls;
        DiST::InterfaceCL interf( AllLvls, prioMaster, prioMaster, dimlist, false);
        interf.ExecuteLocal( *this);
    }
};

/// \brief Set prios of all subsimplices right
/** This uses TreatHasGhosts(), TreatGhosts() and RescueSubs() */
void ParMultiGridCL::AdaptPrioOnSubs()
{
    Comment( "Adapting prios on subsimplices...\n", DebugParallelC);
    ModifyBegin();
    TreatHasGhostsCL treatHasGhost(-1, *modify_);
    RescueGhostsCL rescueGhostSubs(-1, *modify_);
    RescueMasterCL rescueMasterSubs(-1, *modify_);
    treatHasGhost.Call();
    rescueGhostSubs.Call();
    rescueMasterSubs.Call(); // for all master tetras
    // additionally, we have to treat all ghost tetras which are not marked for no refinement
//    for (MultiGridCL::TetraIterator it= mg_->GetAllTetraBegin(), end= mg_->GetAllTetraEnd(); it!=end; ++it)
//        if (!it->IsMarkedForRemovement() && it->IsGhost() && !it->IsMarkedForNoRef())
//            rescueMasterSubs(*it);
    ModifyEnd();
    Comment( "AdaptPrioOnSubs done.\n", DebugParallelC);
}

class ParMultiGridCL::AdaptMidVertexCL
{
  public:
    AdaptMidVertexCL() {}
    /** \brief Delete mid-vertex of PrioVGhost edges and set mid-vertex of the other edges, if necessary. Set MFR of local edges, if necessary. */
    void operator() ( EdgeCL& e)
    {
        if (e.GetPrio()==PrioVGhost)
            e.RemoveMidVertex();
        else if ( !e.GetMidVertex() && e.IsMarkedForRef()) {
            const DiST::GeomIdCL midVertGID( e.GetLevel()+1, e.GetGID().bary, DiST::GetDim<VertexCL>());
            e.SetMidVertex( DiST::InfoCL::Instance().GetVertex(midVertGID));
        }
        if (e.IsLocal() && e.GetMFR() != e.GetAccMFR()) // set local MFR = AccMFR
            e.MFR_= e.AccMFR_;
    }
};


/// \brief Adapt mid-vertex pointers on edges
void ParMultiGridCL::AdaptMidVertex()
{
    AdaptMidVertexCL adaptMidVerts;
    std::for_each( mg_->GetAllEdgeBegin(), mg_->GetAllEdgeEnd(), adaptMidVerts);
}

/// \brief Change the priority of a tetra
template<>
  void ParMultiGridCL::PrioChange<TetraCL>(TetraCL* const Tp, Priority Prio)
{
    Assert( modify_, DROPSErrCL("ParMultiGridCL::PrioChange: There must be an active Transfer or Modify module to run this procedure"), DebugParallelC);
    Assert(!(Prio==PrioGhost && Tp->IsMaster() && !Tp->IsLocal() ),
             DROPSErrCL("ParMultiGridCL::PrioChange: illegal prio for T"),
             DebugParallelC
          );
    modify_->ChangePrio( *Tp, Prio);
}

/****************************************************************************
* B E G I N   A N D   E N D   T R A N S F E R                               *
*****************************************************************************
* Call these functions to start and end transfer-commands                   *
* They are preparing and finishing everything for the transfer.             *
****************************************************************************/
/// \brief prepare everything for the transfers
///
/// Call this everytime before using a transfer command!
void ParMultiGridCL::TransferBegin(int Level)
{
    Assert( !transfer_, DROPSErrCL("ParMultiGridCL::TransferBegin: already called TransferBegin"), DebugParallelC);
    Assert( !modify_, DROPSErrCL("ParMultiGridCL::TransferBegin: cannot use Transfer module inside Modify module, call ModifyEnd() before"), DebugParallelC);
    modify_= transfer_= new DiST::TransferCL( *mg_, true, true);

    // all procs should have the same number of levels!
    mg_->PrepareModify();
    AdjustLevel();
    mg_->FinalizeModify();

    level_ = (Level==-1) ? mg_->GetLastLevel() : Level;

    Comment("- Starting Transfer"<<std::endl, DebugParallelC);
    transfer_->Init();
}

/// \brief End the transfer phase by sending simplices, delete unused simplices and eventually rescue unknowns
/** */
void ParMultiGridCL::TransferEnd()
{
    Assert( transfer_ , DROPSErrCL("ParMultiGridCL::TransferEnd: not in Transfer mode"), DebugParallelC);
    transfer_->Finalize();

    delete transfer_;    modify_= transfer_= 0;

    ModifyBegin();

    Comment("  * Rescue HasGhosts"<<std::endl,DebugParallelC);
    TreatHasGhosts();       // mark subs of HasGhosts as VGhost, rescue subs
    Comment("  * Rescue Ghosts"<<std::endl,DebugParallelC);
    TreatGhosts();          // rescue subs of Ghosts, mark as Ghost
    Comment("  * Rescue Masters"<<std::endl,DebugParallelC);
    RescueMasterCL rescueMasterSubs(-1, *modify_);
    rescueMasterSubs.Call();
    ModifyEnd();

    // Accumulate Ref-counter on edges
    Comment("  * Accumulate MFR"<<std::endl,DebugParallelC);
    AccumulateMFR();

    // Adapt midvertex pointers on edges
    Comment("  * Adapting Midvertex on Edges"<<std::endl,DebugParallelC);
    AdaptMidVertex();

    // important: clear triang cache, otherwise migration of unknowns not working properly
    mg_->ClearTriangCache();

    Comment("- Transfer finished"<<std::endl,DebugParallelC);
}

void ParMultiGridCL::MarkSimplicesForUnknowns()
{
    for (MultiGridCL::VertexIterator it= mg_->GetAllVertexBegin(); it != mg_->GetAllVertexEnd(); ++it)
        it->Unknowns.DisableAllUnknowns( mg_->GetNumLevel());

    for (MultiGridCL::EdgeIterator it= mg_->GetAllEdgeBegin(); it != mg_->GetAllEdgeEnd(); ++it)
        it->Unknowns.DisableAllUnknowns( mg_->GetNumLevel());

    for (MultiGridCL::FaceIterator it= mg_->GetAllFaceBegin(); it != mg_->GetAllFaceEnd(); ++it)
        it->Unknowns.DisableAllUnknowns( mg_->GetNumLevel());

    for (MultiGridCL::TetraIterator it= mg_->GetAllTetraBegin(); it != mg_->GetAllTetraEnd(); ++it)
        it->Unknowns.DisableAllUnknowns( mg_->GetNumLevel());

    for (Uint lvl= 0; lvl <= mg_->GetLastLevel(); ++lvl)
    {
        for (MultiGridCL::TriangTetraIteratorCL tit(mg_->GetTriangTetraBegin(lvl));
             tit!=mg_->GetTriangTetraEnd(lvl); ++tit)
        {
        // master tetras in current triang level are able to store unknowns on their sub-simplices
            for (TetraCL::const_VertexPIterator it(tit->GetVertBegin()); it!=tit->GetVertEnd(); ++it)
                (**it).Unknowns.EnableUnknowns(lvl);

            for (TetraCL::const_EdgePIterator it(tit->GetEdgesBegin()); it!=tit->GetEdgesEnd(); ++it)
              (**it).Unknowns.EnableUnknowns(lvl);

            for (TetraCL::const_FacePIterator it(tit->GetFacesBegin()); it!=tit->GetFacesEnd(); ++it)
              (**it).Unknowns.EnableUnknowns(lvl);

            tit->Unknowns.EnableUnknowns(lvl);
        }
    }
}


/****************************************************************************
* I D E N T I F Y  -  F U N C T I O N S                                     *
*****************************************************************************
*   The following functions identify Vertices, Edges or Faces. If procs     *
*   create the same subsimplices independently from each other, this has    *
*   to be done to update the remote data lists accordingly.                 *
****************************************************************************/
/// \brief Identify a vertex by parent edge
void ParMultiGridCL::IdentifyVertex( const EdgeCL* Parent)
{
    VertexCL* vp= const_cast<VertexCL*>( Parent->GetMidVertex());
    vp->Identify( *Parent, priosId_);
}

/// \brief Identify an edge by parent edge
void ParMultiGridCL::IdentifyEdge( EdgeCL* Me, const EdgeCL* Parent)
{
    Me->Identify( *Parent, priosId_);
}

/// \brief Identify an edge by parent face and two vertices
void ParMultiGridCL::IdentifyEdge( EdgeCL* Me, const FaceCL* Parent)
{
    Me->Identify( *Parent, priosId_);
}

/// \brief Identify a face with parent face and a number
void ParMultiGridCL::IdentifyFace( FaceCL* Me, const FaceCL* Parent)
{
    Me->Identify( *Parent, priosId_);
}


/// \brief Send a Tetra
void ParMultiGridCL::Transfer(TetraCL &t, int dest, Priority prio, bool del)
{
    transfer_->Transfer( t, dest, prio, del);
    if (t.IsRegularlyRef() && t.IsMaster() && prio==PrioMaster)
        t.UnCommitRegRefMark();
}


/****************************************************************************
* G E T  M G                                                                *
****************************************************************************/
/// \brief receive a reference to the stored MultiGrid
MultiGridCL& ParMultiGridCL::GetMG()
{
    Assert(mg_!=0, DROPSErrCL("ParMultiGridCL: GetMG: No MultiGrid is assigned"), DebugParallelC);
    return *mg_;
}


/****************************************************************************
* D E B U G I N F O                                                         *
****************************************************************************/
/// \brief Show simplex by GID
/** Search on proc (or all procs) for a simplex of given GID and show DebugInfo*/
void ParMultiGridCL::Show( const DiST::GeomIdCL& gid, char *mesg, int proc) const
{
    for (Uint l= 0; l<= mg_->GetLastLevel(); ++l)
    {
        ShowSimplex( mg_->GetVerticesBegin(l), mg_->GetVerticesEnd(l), gid, mesg, proc);
        ShowSimplex( mg_->GetEdgesBegin(l), mg_->GetEdgesEnd(l), gid, mesg, proc);
        ShowSimplex( mg_->GetFacesBegin(l), mg_->GetFacesEnd(l), gid, mesg, proc);
        ShowSimplex( mg_->GetTetrasBegin(l), mg_->GetTetrasEnd(l), gid, mesg, proc);
    }

    if (proc != -1 && proc == ProcCL::MyRank() )
    {
        const DiST::RemoteDataListCL& rdl=
            DiST::InfoCL::Instance().GetRemoteList( gid.dim);
        DiST::RemoteDataListCL::const_iterator it= rdl.find( gid);
        if ( it==rdl.end()){
            std::cerr << "...stored only locally." << std::endl;
        }
        else{
            DiST::RemoteDataCL::ProcList_const_iterator pit=
                it->second.GetProcListBegin();
            for ( ; pit!=it->second.GetProcListEnd(); ++pit){
                std::cerr << "...stored on proc " << pit->proc
                          << " with prio " << PriorityToString(pit->prio)
                          << std::endl;
            }
        }
    }
}


/// \brief Writes all vertices, edges, faces and tetrahedra onto the stream
void ParMultiGridCL::DebugInfo(std::ostream & os) const
{
    os << "I have:\n";
    os << mg_->GetVertices().size() << " vertices, " << mg_->GetEdges().size() << " edges, "
       << mg_->GetFaces().size() << " faces," << mg_->GetTetras().size() << " tetras" << std::endl;

    os << "\nThe Vertices are:\n";
    MultiGridCL::const_VertexIterator vit=mg_->GetAllVertexBegin();
    for (; vit!=mg_->GetAllVertexEnd(); ++vit){
        vit->DebugInfo(os);
//         if (vit->Unknowns.Exist())
//         {
//             for (size_t i=0; i<_VecDesc.size(); ++i)
//             {
//                 const Uint idx=_VecDesc[i]->RowIdx->GetIdx();
//                 if (vit->Unknowns.Exist(idx))
//                 {
//                     IdxT sysnum=vit->Unknowns(idx);
//                     if (!vit->Unknowns.Unkreceived(idx))
//                         os << " Unknowns of idx "<<idx<<" at pos "<<sysnum<<" out of Vector: "<< _VecDesc[i]->Data[sysnum]<<std::endl;
//                 }
//                 else
//                     os << " No Unknowns of idx "<<idx<<std::endl;
//             }
//         }
    }

    os << "\nThe Edges are:\n";
    MultiGridCL::const_EdgeIterator eit=mg_->GetAllEdgeBegin();
    for (; eit!=mg_->GetAllEdgeEnd(); ++eit)
    {
        eit->DebugInfo(os);
//         if (eit->Unknowns.Exist())
//         {
//             for (size_t i=0; i<_VecDesc.size(); ++i)
//             {
//                 const Uint idx=_VecDesc[i]->RowIdx->GetIdx();
//                 if (eit->Unknowns.Exist(idx))
//                 {
//                     IdxT sysnum=eit->Unknowns(idx);
//                     if (!eit->Unknowns.Unkreceived(idx))
//                         os << " Unknowns of idx "<<idx<<" out of Vector: "<< _VecDesc[i]->Data[sysnum]<<std::endl;
//                 }
//                 else
//                     os << " No Unknowns of idx "<<idx<<std::endl;
//             }
//         }
    }

    os << "\nThe Faces are:\n";
    MultiGridCL::const_FaceIterator fit=mg_->GetAllFaceBegin();
    for (; fit!=mg_->GetAllFaceEnd(); ++fit)
        fit->DebugInfo(os);

    os << "\nThe Tetras are:\n";
    MultiGridCL::const_TetraIterator sit(mg_->GetAllTetraBegin());
    for (; sit!=mg_->GetAllTetraEnd(); ++sit)
        sit->DebugInfo(os);

    os << "\n";
}

/// \brief Calculate balance of tetras over procs in the last triangulation level
/** This function compares the largest number of tetras in the last triangulation level
    with the number of the smallest number of tetras and returns the ratio.
    (1 is perfect, 0 is totally unbalanced) */
double ParMultiGridCL::GetBalance() const
{
    int myTetras = mg_->TriangTetra_.size();    // all procs count their tetras
    int maxTetras = ProcCL::GlobalMax(myTetras);
    int minTetras = ProcCL::GlobalMin(myTetras);

    return std::fabs(1- (double)(maxTetras)/(double)(minTetras));
}

bool ParMultiGridCL::ConsCheck( std::ostream& os) const
{
    const bool cons= DROPS::ProcCL::Check( DROPS::DiST::InfoCL::Instance().IsSane( os));
    if ( cons )
        std::cout << " DiST-module seems to be alright!" << std::endl;
    else
        std::cout << " DiST-module seems to be broken!" << std::endl;
    return cons;
}

void PrintMG(const DROPS::ParMultiGridCL& pmg, int type)
/** Write out all information about vertices, edges, faces and tetrahedra that
    can be get via the DebugInfo member function of these classes
    \param pmg  the parallel multigrid
    \param type after refinement or migration
*/
{
    const int me=DROPS::ProcCL::MyRank();
    static int REFnum=0;
    static int MIGnum=0;
    char filename[30];

    if (type==REF)
        std::sprintf(filename, "output/%img__REF_%i.mg",me,REFnum++);
    else if (type == MIG)
        std::sprintf(filename, "output/%img__MIG_%i.mg",me,MIGnum++);
    if (me==0)
        std::cout << " - Writing multigrid into: " << filename<< " (for proc 0)"<<std::endl;

    std::ofstream file(filename);
    pmg.DebugInfo(file);
    file.close();
}

bool CheckParMultiGrid()
/**  Check parallel and sequential multigrid on each processor
*/
{
    char dat[30];
    std::ostream output (std::cout.rdbuf());
    std::sprintf(dat,"output/sane%i.chk",DROPS::ProcCL::MyRank());
    std::ofstream checkfile(dat);
    if (!checkfile){
        IF_MASTER
          std::cout << "Cannot open file "<<dat<<" to write sanity check output. Using std::cout"<<std::endl;
    }
    else{
        output.rdbuf(checkfile.rdbuf());
    }
    const ParMultiGridCL& pmg= ParMultiGridCL::Instance();
    bool pmg_sane = pmg.IsSane(output),
         mg_sane  = pmg.GetMG().IsSane(output),
         dist_sane= pmg.ConsCheck(output);
    bool sane     = ProcCL::Check(pmg_sane && mg_sane && dist_sane);
    if (!sane)
        throw DROPSErrCL("CheckParMultiGrid: Multigrid is not sane!");
    return sane;
}

} // end of namespace DROPS
