/// \file parmultigrid.cpp
/// \brief handling of a parallel multigrid
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

#include "parallel/parallel.h"
#include "parallel/parmultigrid.h"
#include "num/solver.h"
#include <iostream>
#include <limits>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

namespace DROPS
{
/****************************************************************************
* I N I T I A L   S T A T I C   A S S I G N M E N T                         *
*****************************************************************************
* initial assignment of the static members of ParMultiGridCL                *
****************************************************************************/
ParMultiGridCL* ParMultiGridCL::instance_ = 0;

/****************************************************************************
* W R A P P E R                                                             *
*****************************************************************************
* Converts C++ functions to C functions, so DDD can use them correctly      *
****************************************************************************/
/*
extern "C" void DeleteObjC(void* buffer, size_t size, int ddd_type) {ParMultiGridCL::DeleteObj(buffer,size,ddd_type);}

extern "C" void HandlerVDeleteC(OBJT obj) {ParMultiGridCL::HandlerDelete<VertexCL>(obj);}
extern "C" void HandlerEDeleteC(OBJT obj) {ParMultiGridCL::HandlerDelete<EdgeCL>(obj);}
extern "C" void HandlerFDeleteC(OBJT obj) {ParMultiGridCL::HandlerDelete<FaceCL>(obj);}
extern "C" void HandlerTDeleteC(OBJT obj) {ParMultiGridCL::HandlerDelete<TetraCL>(obj);}

extern "C" OBJT HandlerVConstructorC(size_t s, PrioT p, ATTRT l) {return ParMultiGridCL::HandlerVConstructor(s,p,l);}
extern "C" OBJT HandlerEConstructorC(size_t s, PrioT p, ATTRT l) {return ParMultiGridCL::HandlerEConstructor(s,p,l);}
extern "C" OBJT HandlerFConstructorC(size_t s, PrioT p, ATTRT l) {return ParMultiGridCL::HandlerFConstructor(s,p,l);}
extern "C" OBJT HandlerTConstructorC(size_t s, PrioT p, ATTRT l) {return ParMultiGridCL::HandlerTConstructor(s,p,l);}

extern "C" void HandlerVXferC( OBJT o, PROCT to, PrioT p) {ParMultiGridCL::HandlerVXfer(o,to,p);}
extern "C" void HandlerEXferC( OBJT o, PROCT to, PrioT p) {ParMultiGridCL::HandlerEXfer(o,to,p);}
extern "C" void HandlerFXferC( OBJT o, PROCT to, PrioT p) {ParMultiGridCL::HandlerFXfer(o,to,p);}
extern "C" void HandlerTXferC( OBJT o, PROCT to, PrioT p) {ParMultiGridCL::HandlerTXfer(o,to,p);}

extern "C" void HandlerVGatherC(OBJT o, int i, TypeT t, void* d) {ParMultiGridCL::HandlerVGather(o,i,t,d);}
extern "C" void HandlerEGatherC(OBJT o, int i, TypeT t, void* d) {ParMultiGridCL::HandlerEGather( o,i,t,d);}
extern "C" void HandlerTGatherC(OBJT o, int i, TypeT t, void* d) {ParMultiGridCL::HandlerTGather( o,i,t,d);}

extern "C" void HandlerVScatterC(OBJT o, int i, TypeT t, void* d, int n) {ParMultiGridCL::HandlerVScatter(o,i,t,d,n);}
extern "C" void HandlerEScatterC(OBJT o, int i, TypeT t, void* d, int n) {ParMultiGridCL::HandlerEScatter(o,i,t,d,n);}
extern "C" void HandlerTScatterC(OBJT o, int i, TypeT t, void* d, int n) {ParMultiGridCL::HandlerTScatter(o,i,t,d,n);}

extern "C" void HandlerTUpdateC(OBJT o)              {ParMultiGridCL::HandlerTUpdate(o);}
extern "C" void HandlerTObjMkConsC( OBJT o, int i)   {ParMultiGridCL::HandlerTObjMkCons(o,i);}
extern "C" void HandlerTSetPrioC(OBJT o, PrioT p) {ParMultiGridCL::HandlerTSetPrio(o,p);}

extern "C" int GatherEdgeMFRC(OBJT o, void* d) {return ParMultiGridCL::GatherEdgeMFR(o,d);}
extern "C" int ScatterEdgeMFRC(OBJT o, void* d) {return ParMultiGridCL::ScatterEdgeMFR(o,d);}

extern "C" int GatherTetraRestrictMarksC( OBJT o, void* b) {return ParMultiGridCL::GatherTetraRestrictMarks(o,b);}
extern "C" int ScatterTetraRestrictMarksC(OBJT o, void* b) {return ParMultiGridCL::ScatterTetraRestrictMarks(o,b);}

extern "C" int ExecGhostRescueC(         OBJT o) {return ParMultiGridCL::ExecGhostRescue(o);}
extern "C" int ExecGhVertRescueC(        OBJT o) {return ParMultiGridCL::ExecGhVertRescue(o);}
extern "C" int ExecHasGhostC(            OBJT o) {return ParMultiGridCL::ExecHasGhost(o);}
extern "C" int ExecAdaptVGhostMidVertexC(OBJT o) {return ParMultiGridCL::ExecAdaptVGhostMidVertex(o);}

extern "C" int GatherEdgeSaneC( OBJT o, void* d, PROCT p, ATTRT a) {return ParMultiGridCL::GatherEdgeSane(o,d,p,a);}
extern "C" int ScatterEdgeSaneC(OBJT o, void* d, PROCT p, ATTRT a) {return ParMultiGridCL::ScatterEdgeSane(o,d,p,a);}
extern "C" int GatherFaceSaneC( OBJT o, void* d, PROCT p, ATTRT a) {return ParMultiGridCL::GatherFaceSane(o,d,p,a);}
extern "C" int ScatterFaceSaneC(OBJT o, void* d, PROCT p, ATTRT a) {return ParMultiGridCL::ScatterFaceSane(o,d,p,a);}

extern "C" int GatherUnknownsRefC (OBJT o, void* b) { return ParMultiGridCL::GatherUnknownsRef(o,b); }
extern "C" int ScatterUnknownsRefC(OBJT o, void* b) { return ParMultiGridCL::ScatterUnknownsRef(o,b); }

extern "C" int GatherUnknownsMigVC (OBJT o, void* b) { return ParMultiGridCL::GatherUnknownsMigV(o,b); }
extern "C" int ScatterUnknownsMigVC(OBJT o, void* b) { return ParMultiGridCL::ScatterUnknownsMigV(o,b); }
extern "C" int GatherUnknownsMigEC (OBJT o, void* b) { return ParMultiGridCL::GatherUnknownsMigE(o,b); }
extern "C" int ScatterUnknownsMigEC(OBJT o, void* b) { return ParMultiGridCL::ScatterUnknownsMigE(o,b); }

extern "C" int GatherInterpolValuesVC (OBJT o, void* b) { return ParMultiGridCL::GatherInterpolValues<VertexCL>(o,b); }
extern "C" int ScatterInterpolValuesVC(OBJT o, void* b) { return ParMultiGridCL::ScatterInterpolValues<VertexCL>(o,b); }
extern "C" int GatherInterpolValuesEC (OBJT o, void* b) { return ParMultiGridCL::GatherInterpolValues<EdgeCL>(o,b); }
extern "C" int ScatterInterpolValuesEC(OBJT o, void* b) { return ParMultiGridCL::ScatterInterpolValues<EdgeCL>(o,b); }
*/

/****************************************************************************
* C O N S T R U C T O R S                                                   *
*****************************************************************************
* Constructor an Destructor of the parallel multigrid                       *
* AttachTo: assign the given Multigrid to the parallel mutigrid             *
*           or tell the Multigrid about a VectorDescriber-Class             *
****************************************************************************/
/// \brief Constructor with number of Vector-Describer-Classes
///
/// Init the parallel stuff
ParMultiGridCL::ParMultiGridCL()
    : mg_(0), modify_(0), transfer_(0)
{
    Assert( instance_==0, DROPSErrCL("ParMultiGridCL: Constructor is called twice"), DebugParallelC);

    // Init the containers for the unknowns
    _RecvBuf.resize(0);

    // There are no active transfers neither Prio-Change-enviroment while creating this class
    _UnkOnSimplex[0]= _UnkOnSimplex[1]= _UnkOnSimplex[2]= false;

    // for Identify, we need all priorities except PrioVGhost
    priosId_.push_back(PrioNeutral); priosId_.push_back(PrioKilledGhost);
    priosId_.push_back(PrioGhost);   priosId_.push_back(PrioMaster);
}

ParMultiGridCL::~ParMultiGridCL()
{
//  if (_UnkOnSimplex)
//      delete[] _UnkOnSimplex;
}

/// \brief Assign the MultiGridCL to this ParMultiGridCL
void ParMultiGridCL::AttachTo( MultiGridCL& mg)
{
    // Store referenz to the Multigrid and to the simplex containers
    mg_= &mg;
//    _VertCont=  &mg_->_Vertices;
//    _EdgeCont=  &mg_->_Edges;
//    _FaceCont=  &mg_->_Faces;
//    _TetraCont= &mg_->_Tetras;
}

/// \brief Store a pointer to a scalar boundary condition, that belongs to an Index
template<>
  void ParMultiGridCL::AttachTo<BndDataCL<double> >(const IdxDescCL* idxDesc, const BndDataCL<double>* bndp)
/** Scalar boundary conditions are given. Hence, a pointer to the boundary
    consitions are stored _ScalBnd in at the same position where the VecDescCL
    can be found.
    \pre The VecDesc, corresponding to \a idxDesc,  must have been set by
         AttachTo(VecDescCL*) before calling this routine
    \param idxDesc IdxDescCL matching to \a bndp
    \param bndp    pointer to the correspondint BndDataCL  */
{
    size_t vecPos= GetStorePos( idxDesc);
    Assert( vecPos!=_VecDesc.size() && vecPos<_ScalBnd.size(),
            DROPSErrCL("ParMultiGridCL::AttachTo<BndDataCL<double>>: VecDesc is not known so far, set it with AttachTo before calling this routine"),
            DebugParallelNumC);
    _ScalBnd[vecPos]= bndp;
}

/// \brief Store a pointer to a vectorial boundary condition, that belongs to an Index
template<>
  void ParMultiGridCL::AttachTo<BndDataCL<Point3DCL> >(const IdxDescCL* idxDesc, const BndDataCL<Point3DCL>* bndp)
/** Vectorial boundary conditions are given. Hence, a pointer to the boundary
    consitions are stored _VecBnd in at the same position where the VecDescCL
    can be found.
    \pre The VecDesc, corresponding to \a idxDesc,  must have been set by
         AttachTo(VecDescCL*) before calling this routine
    \param idxDesc IdxDescCL matching to \a bndp
    \param bndp    pointer to the correspondint BndDataCL  */
{
    size_t vecPos= GetStorePos( idxDesc);
    Assert( vecPos!=_VecDesc.size() && vecPos<_VecBnd.size(),
            DROPSErrCL("ParMultiGridCL::AttachTo<BndDataCL<double>>: VecDesc is not known so far, set it with AttachTo before calling this routine"),
            DebugParallelNumC);
    _VecBnd[vecPos]= bndp;
}

/// \brief Delete all information about the VecDescCL
void ParMultiGridCL::DeleteVecDesc()
{
    _VecDesc.resize( 0);
    _VecBnd.resize( 0);
    _ScalBnd.resize( 0);
    _UnkOnSimplex[0]=false;
    _UnkOnSimplex[1]= false;
    _UnkOnSimplex[2]= false;
}

/// \brief Get scalar boundary condition to a known VecDesCL
template<>
  const BndDataCL<double>* ParMultiGridCL::GetBndCond<BndDataCL<double> >( const IdxDescCL* idxDesc)
/// First, find the position, where the VecDescCL is stored. Second, the pointer
/// to the boundary conditions is returned
{
    size_t vecPos= GetStorePos( idxDesc);
    Assert(vecPos<_VecDesc.size() && vecPos<_ScalBnd.size() && _ScalBnd[vecPos]!=0,
           DROPSErrCL("ParMultiGridCL::GetBndCond<BndDataCL<double>>: BC not set"),
           DebugParallelNumC);
    return _ScalBnd[vecPos];
}

/// \brief Get vectorial boundary condition to a known VecDesCL
template<>
  const BndDataCL<Point3DCL>* ParMultiGridCL::GetBndCond<BndDataCL<Point3DCL> >( const IdxDescCL* idxDesc)
/// First, find the position, where the VecDescCL is stored. Second, the pointer
/// to the boundary conditions is returned
{
    size_t vecPos= GetStorePos( idxDesc);
    Assert(vecPos<_VecDesc.size() && vecPos<_VecBnd.size() && _VecBnd[vecPos]!=0,
           DROPSErrCL("ParMultiGridCL::GetBndCond<BndDataCL<Point3DCL>>: BC not set"),
           DebugParallelNumC);
    return _VecBnd[vecPos];
}

/// \brief Handle Unknowns after a refine operation
void ParMultiGridCL::HandleUnknownsAfterRefine(/*const std::vector<VecDescCL*> &newVecDesc*/)
/**
    This procedure performs the following operations:
    <ol>
     <li>Transfer unknowns from killed ghost tetras to the master tetras</li>
     <li>Delete all subsimplices, that are marked for removement by the refinement algorithm</li>
     <li>Unsubscribe the killed ghost tetras from the DDD-System</li>
     <li>Delete ghost tetras that are not needed any more</li>
     <li>Copy received values into the receive buffer</li>
    </ol>
    \pre  The VecDescCL's, that describes the unknowns before the refinement algorithm has been performed,
          has to be attached to the ParMultiGridCL by the procedure 'AttachTo'.
    \post All subsimplices, that stores unknowns has an accumulated value of its unknowns
*/
{
    Assert(_RecvBuf.size()==0, DROPSErrCL("ParMultiGridCL::HandleUnknownsAfterRefine: receive Buffer is not empty!"), DebugParallelNumC);
    if (!mg_->UnknownsForRefine())
        return;
    Assert(VecDescRecv(), DROPSErrCL("ParMultiGridCL::HandleUnknownsAfterRefine: missing vector describers"), DebugParallelNumC);

    // Allocate mem for receive buffer
    IdxT numKnownUnknowns=0;
    for (size_t i=0; i<_VecDesc.size(); ++i)
        numKnownUnknowns+=_VecDesc[i]->Data.size();
    _RecvBuf.resize(numKnownUnknowns);

    // Transfer the unknowns to the right processor
    /// \todo DiST: Implement Transfer of unknowns
//    DynamicDataInterfaceCL::IFOneway( _GhToMaTetra_IF,                               // transfer unknowns from killed ghost-tetra to corresponding master tetra
//                  IF_FORWARD,                                    // transfer in one direction
//                  NumberOfUnknownsOnTetra()*sizeof(TransferUnkT),// number of (double+bool) values
//                  &GatherUnknownsRefC, &ScatterUnknownsRefC);    // handler for gather and scatter

    // Right now no more unused simplices are needed any more!
    mg_->PrepareModify();

    // kill all subsimplices, that has not beed done by refinement algorithm
    for (Uint l=0; l<=mg_->GetLastLevel(); ++l)
    {
        mg_->Vertices_[l].remove_if( std::mem_fun_ref(&VertexCL::IsMarkedForRemovement) );
        mg_->Edges_[l].remove_if( std::mem_fun_ref(&EdgeCL::IsMarkedForRemovement) );
        mg_->Faces_[l].remove_if( std::mem_fun_ref(&FaceCL::IsMarkedForRemovement) );
    }

    // Kill Ghosts tetras that are not needed any more
    if ( mg_->KilledGhosts())
    {
        std::list<MultiGridCL::TetraIterator>::iterator end(mg_->toDelGhosts_.end());
        // unsubscribe them from the DiST-System
        DiST::ModifyCL modify( *mg_, true, true);
        modify.Init();
        for (std::list<MultiGridCL::TetraIterator>::iterator it(mg_->toDelGhosts_.begin()); it!=end; ++it)
            modify.Delete( **it);
        modify.Finalize();
        // information about ghost, that should be deleted are not needed any more
        mg_->toDelGhosts_.resize(0);
    }

    // there are no more ghost tetras for deleting
    mg_->killedGhostTetra_= false;

    // Adjust level
    while ( mg_->GetLastLevel()>0 && mg_->Tetras_[mg_->GetLastLevel()].empty() )
        mg_->RemoveLastLevel();
    AdjustLevel();

    // no more modifications has to be done
    mg_->FinalizeModify();
    mg_->ClearTriangCache();
}

/// \brief Fill a new vector describer class with datas (Has to be called after the refinement an migration algorithm!)
void ParMultiGridCL::HandleNewIdx(IdxDescCL* oldIdxDesc, VecDescCL* newVecDesc)
/** Within the refinement or migration algorithm a single processor
    can get new unknowns and give other unknowns to other processors. So each
    processor put the received unknowns in a buffer. In order to get a fully filled
    data vector, this routine builds this vector by taking values out of the
    receive buffer or the "old" vector.*/
/*  This routine iterates over all vertices, edges and tetras and calls the routine PutData.
    PutData collects the data and put it into the new vector at the right position.*/
{
    Assert(VecDescRecv(), DROPSErrCL("ParMultiGridCL::HandleNewIdx: No Indices received before transfer"), DebugParallelNumC);

    // old and new index
    const Uint old_idx= oldIdxDesc->GetIdx(),
               new_idx= newVecDesc->RowIdx->GetIdx();

    // find the right index within _VecDesc
    size_t vecIdx= GetStorePos( oldIdxDesc);
    Assert(vecIdx!=_VecDesc.size(),
           DROPSErrCL("ParMultiGridCL::HandleNewIdx: The appropriate VecDesc has not been set"),
           DebugParallelNumC);

    // create the new data vector
    newVecDesc->Data.resize(newVecDesc->RowIdx->NumUnknowns());

    // abbrevations for accessing the data vectors
    VectorCL *new_data= &(newVecDesc->Data);
    const VectorCL* const old_data= &(_VecDesc[vecIdx]->Data);

    // error checking
    Assert(newVecDesc->RowIdx->NumUnknownsVertex() == oldIdxDesc->NumUnknownsVertex(),
           DROPSErrCL("ParMultiGridCL::HandleNewIdx: number of unknowns on vertices does not match"), DebugParallelNumC);
    Assert(newVecDesc->RowIdx->NumUnknownsEdge()==oldIdxDesc->NumUnknownsEdge(),
           DROPSErrCL("ParMultiGridCL::HandleNewIdx: number of unknowns on edges does not match"), DebugParallelNumC);
    Assert(newVecDesc->RowIdx->NumUnknownsTetra()==oldIdxDesc->NumUnknownsTetra(),
           DROPSErrCL("ParMultiGridCL::HandleNewIdx: number of unknowns on tetras does not match"), DebugParallelNumC);

    // Put unknowns on vertices into the new vector
    for (MultiGridCL::const_VertexIterator sit=mg_->GetAllVertexBegin();
            sit!=mg_->GetAllVertexEnd(); ++sit)
    {
        PutData(sit, old_data, new_data, old_idx, new_idx, oldIdxDesc);
    }

    // Put unknowns on edges into the new vector and interpolate very special midvertices
    if (oldIdxDesc->NumUnknownsEdge()){
        for (MultiGridCL::const_EdgeIterator sit=mg_->GetAllEdgeBegin();
                sit!=mg_->GetAllEdgeEnd(); ++sit)
        {
            if (oldIdxDesc->NumUnknownsVertex()==1)
                PutData(sit, old_data, new_data, old_idx, new_idx, oldIdxDesc, GetBndCond<BndDataCL<double> >(oldIdxDesc));
            else
                PutData(sit, old_data, new_data, old_idx, new_idx, oldIdxDesc, GetBndCond<BndDataCL<Point3DCL> >(oldIdxDesc));
        }
    }
}

/// \brief Send interpolated values from processers that can interpolate a value to processors that could not interpolate a value
void ParMultiGridCL::CompleteRepair(VecDescCL* /*newVec*/)
/** After RepairAfterRefine[P1|P2] call this function to exchange interpolated
    values.
    \todo DiST: Implement me
    \param newVec the interpolated VecDescCL
*/
{
/*
    Assert(newVec->RowIdx->NumUnknownsEdge()==0 || newVec->RowIdx->NumUnknownsEdge()==newVec->RowIdx->NumUnknownsVertex(),
           DROPSErrCL("ParMultiGridCL::CompleteRepair: Unknowns on vertices and edges must be the same"),
           DebugParallelNumC);

    if (newVec->RowIdx->NumUnknownsVertex()){
        _actualVec=newVec;             // gather and scatter functions has to know about the index and values
        DynamicDataInterfaceCL::IFExchange(InterfaceCL<VertexCL>::GetIF(),                              // exchange datas over distributed vertices
                    _actualVec->RowIdx->NumUnknownsVertex() * sizeof(TransferUnkT), // number of datas to be exchanged
                    GatherInterpolValuesVC,                                         // how to gather datas
                    ScatterInterpolValuesVC                                         // how to scatter datas
                    );
        _actualVec=0;                           // fortget about the actual index and values
    }

    if (newVec->RowIdx->NumUnknownsEdge()){
        _actualVec=newVec;
        DynamicDataInterfaceCL::IFExchange(InterfaceCL<EdgeCL>::GetIF(),
                    _actualVec->RowIdx->NumUnknownsEdge() * sizeof(TransferUnkT),
                    GatherInterpolValuesEC,
                    ScatterInterpolValuesEC
                    );
        _actualVec=0;
    }
*/
}

/// \brief Put datas on a vertex from the old vector into a new vector
void ParMultiGridCL::PutData(MultiGridCL::const_VertexIterator& sit,
                             const VectorCL* const old_data, VectorCL* new_data,
                             const Uint old_idx, const Uint new_idx,
                             const IdxDescCL* idxDesc)
/** This routine puts the unknowns on a vertex according to an index into
    a new data vector. Therefore datas are taken from the receive buffer, if the
    data has been stored on another proc before the refinement and migration, or
    out of the "old" vector, if the calling proc has already owned these data
    before the refinement.
    \param sit         Iterator onto a simplex
    \param old_data    Pointer to the old data vector
    \param new_data    Pointer to the new data vector
    \param old_idx     old index
    \param new_idx     new index
    \param idxDesc     describer of index, just used for getting information
                       about number of unknowns on simplices*/
{
    const Uint numUnknowns=idxDesc->NumUnknownsVertex();
    // check if there are unknowns to copy. (They have to exist on the new and old ParMultiGrid)
    if (sit->Unknowns.Exist()
        && sit->Unknowns.Exist(new_idx)
        && sit->Unknowns.Exist(old_idx) )
    {
        const IdxT new_sysnum = sit->Unknowns(new_idx),                 // position, where to find the new unknown(s)
                   old_sysnum = sit->Unknowns(old_idx);                 // position, where to find the old unknown(s)
        if (sit->Unknowns.UnkReceived(old_idx))                         // this is an "new" unknown
        {
            for (Uint i=0; i<numUnknowns; ++i)
                (*new_data)[new_sysnum+i]= _RecvBuf[old_sysnum+i];
            sit->Unknowns.SetUnkReceived(new_idx);                      // set received flag as a flag, that this unknown is an old one
        }
        else
        {
            for (Uint i=0; i<numUnknowns; ++i)
                (*new_data)[new_sysnum+i]=  (*old_data)[old_sysnum+i];
            sit->Unknowns.SetUnkReceived(new_idx);                      // set received flag as a flag, that this unknown is an old one
        }
    }
}

/// \brief Save unknowns on edges, that are deleted within migration
void ParMultiGridCL::RescueUnknownsOnEdges()
/** Within the migration it may happen, that an edge is removed, but the
    midvertex stays on the processor. If the midvertex has unknowns after the
    migration on this processor and do not has unknowns before the refinement,
    the interpolation must be done, before the edge is deleted. So this is done
    here. <br>
    We distinguish between two cases:
    <ol>
     <li>P1-finite elements: The value on the midvertex is interpolated by the
         values on the two vertices of the edge. If a value is given by
         Dirichlet BC, this must be checked.</li>
     <li>P2-finite elements: The value on the edge is set to the midvertex. Here
         no check for Dirichlet BC is needed.</li>
    </ol>
    \todo DiST: Implement me
*/
{
/*
    Assert(UnknownsOnSimplices(),
           DROPSErrCL("ParMultiGridCL::RescueUnknownsOnEdges: Called without unknowns"),
           DebugParallelNumC);

    // check all edges, if the edge is removed, but its midvertex stays on this
    // processor and the midvertex has no value
    for (int l=0; l<=_level; ++l){
        for (MultiGridCL::EdgeIterator sit((*_EdgeCont)[l].begin()), tEnd((*_EdgeCont)[l].end()); sit!=tEnd; ++sit){
            if (sit->IsRefined()
                && sit->IsMarkedForRemovement()
                && !sit->GetMidVertex()->IsMarkedForRemovement()
                && !sit->GetMidVertex()->Unknowns.Exist())
            {

                // rescue unknowns for all known indices
                for (size_t idx_type=0; idx_type<_VecDesc.size(); ++idx_type){
                    const VecDescCL* sol          = _VecDesc[idx_type];                 // Pointer to the VecDescCL
                    const Uint       idx          = sol->RowIdx->GetIdx();              // index of the unknowns
                    const Uint       numUnkOnEdge = sol->RowIdx->NumUnknownsEdge(),     // unknowns on edges
                                     numUnkOnVert = sol->RowIdx->NumUnknownsVertex();

                    if (numUnkOnEdge==0 && numUnkOnVert){                                 // Rescue unknowns for P1 Elements
                        if (numUnkOnVert==1){   // scalar unknowns
                            const BndDataCL<double>* bnd = GetBndCond<BndDataCL<double> >(sol->RowIdx);
                            double new_dof;
                            if(!bnd->IsOnDirBnd(*sit->GetMidVertex()) ){
                                if ( LinearInterpolation(*sit, idx, bnd, sol->Data, new_dof) ){
                                    PutDofIntoRecvBuffer(*sit->GetMidVertex(), idx, new_dof);
                                    sit->GetMidVertex()->Unknowns.SetUnkreceived(idx);
                                }
                                else
                                    throw DROPSErrCL("ParMultiGridCL::RescueUnknownsOnEdges: Cannot interpolate onto midvertex");
                            }
                        }
                        else{                   // vectorial unknowns
                            const BndDataCL<Point3DCL>* bnd = GetBndCond<BndDataCL<Point3DCL> >(sol->RowIdx);;
                            Point3DCL new_dof;
                            if(!bnd->IsOnDirBnd(*sit->GetMidVertex()) ){
                                if ( LinearInterpolation(*sit, idx, bnd, sol->Data, new_dof) ){
                                    PutDofIntoRecvBuffer(*sit->GetMidVertex(), idx, new_dof);
                                    sit->GetMidVertex()->Unknowns.SetUnkreceived(idx);
                                }
                                else
                                    throw DROPSErrCL("ParMultiGridCL::RescueUnknownsOnEdges: Cannot interpolate onto midvertex");
                            }
                        }
                    }
                    else{                                                               // Rescue unknowns for P2 Elements
                        if (sit->Unknowns.Exist()){
                            if (numUnkOnEdge==1)
                                PutDofIntoRecvBuffer(*sit->GetMidVertex(), idx,
                                    GetDofOutOfVector<EdgeCL, VectorCL, double>()(*sit, idx, sol->Data));
                            else
                                PutDofIntoRecvBuffer(*sit->GetMidVertex(), idx,
                                    GetDofOutOfVector<EdgeCL, VectorCL, Point3DCL>()(*sit, idx, sol->Data));
                            sit->GetMidVertex()->Unknowns.SetUnkreceived(idx);
                        }
                    }
                }
            }
        }
    }
*/
}

/// \brief Delete receive Buffer
/** After a refinement and a miragtion the receive buffer is not needed any more */
void ParMultiGridCL::DeleteRecvBuffer()
{
    _RecvBuf.resize(0);
    _RecvBufPos=0;
}

/// \brief Get overall number of unknowns on a tetrahedron
/** Check number of unknowns on tetra and all its subsimplices for number of
    unknowns according to vecdesc */
Uint ParMultiGridCL::NumberOfUnknownsOnTetra()
{
    Uint result=0;
    for (VecDescPCT::const_iterator it(_VecDesc.begin()), end(_VecDesc.end()); it!=end; ++it)
        result +=  NumVertsC * (*it)->RowIdx->NumUnknownsVertex()
                 + NumEdgesC * (*it)->RowIdx->NumUnknownsEdge()
                 + 1         * (*it)->RowIdx->NumUnknownsTetra();
    return result;
}

/// \brief Get number of unknowns on a vertex
Uint ParMultiGridCL::NumberOfUnknownsOnVertex()
{
    Uint result=0;
    for (VecDescPCT::const_iterator it(_VecDesc.begin()), end(_VecDesc.end()); it!=end; ++it)
        result +=(*it)->RowIdx->NumUnknownsVertex();
    return result;
}

/// \brief Get number of unknowns on an edge
Uint ParMultiGridCL::NumberOfUnknownsOnEdge()
{
    Uint result=0;
    for (VecDescPCT::const_iterator it(_VecDesc.begin()), end(_VecDesc.end()); it!=end; ++it)
        result += (*it)->RowIdx->NumUnknownsEdge();
    return result;
}

/// \brief Enlarge the receive buffer
void ParMultiGridCL::EnlargeReceiveBuffer()
{
    if (_RecvBuf.size()>0){
        BufferCT tmpBuf(_RecvBuf);
        _RecvBuf.resize(2*_RecvBuf.size());
        std::copy(tmpBuf.begin(), tmpBuf.end(), _RecvBuf.begin());
        Comment("["<<ProcCL::MyRank()<<"]===> Enlarge receive buffer from "<<tmpBuf.size()<<" to "<<_RecvBuf.size()<<"!"<<std::endl, DebugParallelC);
    }
    else{
        _RecvBuf.resize(1024);
        Comment("["<<ProcCL::MyRank()<<"]===> Create receive buffer of size "<<_RecvBuf.size()<<"!"<<std::endl, DebugParallelC);
    }
}

/// \brief Gather unknowns on ghost tetras for sending these to master the tetra
int ParMultiGridCL::GatherUnknownsRef (OBJT obj, void* buf)
/** Get all values on a ghost tetrahedron, that will be deleted, and put these values into the
    given buffer.
*/
{
    TetraCL* const tp= ddd_cast<TetraCL*>(obj);

    if (tp->GetPrio()!=PrioKilledGhost || !tp->IsMarkedForNoRef())
        return 1;

    TransferUnkT* buffer    = static_cast<TransferUnkT*>(buf);
    Uint          buffer_pos= 0;

    for (VecDescPCT::const_iterator it(_VecDesc.begin()); it!=_VecDesc.end(); ++it)
    {
        const Uint idx=(*it)->RowIdx->GetIdx(); // this should be the index before the refinement algorithm has been started

        if ((*it)->RowIdx->NumUnknownsVertex())
        { // collect unknowns on vertices
            for (TetraCL::const_VertexPIterator sit(tp->GetVertBegin());
                 sit!=tp->GetVertEnd();
                 ++sit,
                 buffer_pos+=(*it)->RowIdx->NumUnknownsVertex())
            {
                if ((*sit)->Unknowns.Exist() && (*sit)->Unknowns.Exist(idx)){    // check for existance
                    for (Uint i=0; i<(*it)->RowIdx->NumUnknownsVertex(); ++i){
                        buffer[buffer_pos+i].val = (*it)->Data[(*sit)->Unknowns(idx)+i];
                        buffer[buffer_pos+i].idx = idx;
                        buffer[buffer_pos+i].mark= true;
                    }
                }
                else{
                    buffer[buffer_pos].mark=false;
                }
            }
        }
        if ((*it)->RowIdx->NumUnknownsEdge())
        { // collect unknowns on edges, if they exists
            for (TetraCL::const_EdgePIterator sit(tp->GetEdgesBegin());
                 sit!=tp->GetEdgesEnd();
                 ++sit,
                buffer_pos+=(*it)->RowIdx->NumUnknownsEdge())
            {
                if ((*sit)->Unknowns.Exist() && (*sit)->Unknowns.Exist(idx))    // check for existance
                    for (Uint i=0; i<(*it)->RowIdx->NumUnknownsVertex(); ++i){
                        buffer[buffer_pos+i].val = (*it)->Data[(*sit)->Unknowns(idx)+i];
                        buffer[buffer_pos+i].idx = idx;
                        buffer[buffer_pos+i].mark= true;
                    }
                else{
                    buffer[buffer_pos].mark=false;
                }
            }
        }
        if ((*it)->RowIdx->NumUnknownsTetra())
        { // collect unknowns on the tetra itselfe
            if (tp->Unknowns.Exist() && tp->Unknowns.Exist(idx))                // check for existance
                for (Uint i=0; i<(*it)->RowIdx->NumUnknownsTetra(); ++i){
                    buffer[buffer_pos+i].val = (*it)->Data[tp->Unknowns(idx)+i];
                    buffer[buffer_pos+i].idx = idx;
                    buffer[buffer_pos+i].mark= true;
                }
            else{
                buffer[buffer_pos].mark=false;
            }
            buffer_pos+=(*it)->RowIdx->NumUnknownsTetra();
        }
    }
    Assert(buffer_pos==NumberOfUnknownsOnTetra(),
           DROPSErrCL("ParMultiGridCL::GatherUnknowns: I haven't check enough places for unknowns!"),
           DebugParallelNumC);
    return 0;
}

/// \brief Scatter unknowns on master tetras, which have been send from killed ghost tetra
int ParMultiGridCL::ScatterUnknownsRef(OBJT obj, void* buf)
/** This procedure puts all unknowns that can live on subsimplices or the tetrahedron itselfe
    into a receive buffer, if the unknowns are not known so far
*/
{
    TetraCL* const tp= ddd_cast<TetraCL*>(obj);

    if (tp->GetPrio()!=PrioMaster)
        return 1;

    TransferUnkT* buffer    = static_cast<TransferUnkT*>(buf);
    Uint          buffer_pos= 0;

    // for all known indices
    for (size_t index_type=0; index_type<_VecDesc.size(); ++index_type)
    {
        // index number and unknowns on simplices
        const Uint idx          = _VecDesc[index_type]->RowIdx->GetIdx(),
                   lvl          = _VecDesc[index_type]->RowIdx->TriangLevel(),
                   numUnkOnVert = _VecDesc[index_type]->RowIdx->NumUnknownsVertex(),
                   numUnkOnEdge = _VecDesc[index_type]->RowIdx->NumUnknownsEdge(),
                   numUnkOnTetra= _VecDesc[index_type]->RowIdx->NumUnknownsTetra();

        // Check if there are unknowns on vertices and receive them
        if (numUnkOnVert)
        {
            // iterate over all vertices of this tetra
            for (TetraCL::const_VertexPIterator sit(tp->GetVertBegin()); sit!=tp->GetVertEnd(); ++sit)
            {
                // put the received unknowns into _RecvBuf, if this vertex should store unknowns (i.e. is master)
                // and if the unknown is not known so far.
                if (!(*sit)->Unknowns.Exist(idx)                        // this vert stores no "old" unknown
                     && !(*sit)->Unknowns.UnkReceived(idx)              // this vert has no unknowns received yet
                     && buffer[buffer_pos].mark                         // there have been unknowns on sender side
                     && buffer[buffer_pos].idx==idx                     // and right index
                   )
                {
                    if (_RecvBufPos+numUnkOnVert>_RecvBuf.size())
                        EnlargeReceiveBuffer();
                    (*sit)->Unknowns.Prepare(idx);                      // create UnknownIdxCL
                    (*sit)->Unknowns(idx)= _RecvBufPos;                 // remeber position, where the unknowns has been put

                    for (Uint i=0; i<numUnkOnVert; ++i)                 // store all unknowns
                        _RecvBuf[_RecvBufPos+i]=buffer[buffer_pos+i].val;
                    _RecvBufPos+=numUnkOnVert;

                    (*sit)->Unknowns.SetUnkReceived(idx);               // this vertex has received unknowns
                }
                buffer_pos+=numUnkOnVert;                               // unknowns on vertices has been handled
            }
        }

        // for documentation look at vertices above
        if (numUnkOnEdge)
        {
            for (TetraCL::const_EdgePIterator sit(tp->GetEdgesBegin()); sit!=tp->GetEdgesEnd(); ++sit)
            {
                if ((*sit)->Unknowns.InTriangLevel(lvl)
                     && !(*sit)->Unknowns.Exist(idx)
                     && !(*sit)->Unknowns.UnkReceived(idx)
                     && buffer[buffer_pos].mark
                     && buffer[buffer_pos].idx==idx
                   )
                {
                    if (_RecvBufPos + numUnkOnEdge>_RecvBuf.size())
                        EnlargeReceiveBuffer();
                    (*sit)->Unknowns.Prepare(idx);
                    (*sit)->Unknowns(idx)= _RecvBufPos;

                    for (Uint i=0; i<numUnkOnEdge; ++i)
                        _RecvBuf[_RecvBufPos+i]=buffer[buffer_pos+i].val;
                    _RecvBufPos+=numUnkOnEdge;

                    (*sit)->Unknowns.SetUnkReceived(idx);
                }
                buffer_pos+=numUnkOnEdge;
            }
        }

        // receive unknowns on tetra itselfe
        if (numUnkOnTetra)
        {
            if (tp->Unknowns.InTriangLevel(lvl)
                    && !tp->Unknowns.Exist(idx)
                    && !tp->Unknowns.UnkReceived(idx)
                    && buffer[buffer_pos].mark
                    && buffer[buffer_pos].idx==idx
               )
            {
                if (_RecvBufPos + numUnkOnTetra>_RecvBuf.size())
                    EnlargeReceiveBuffer();

                tp->Unknowns.Prepare(idx);
                tp->Unknowns(idx)= _RecvBufPos;

                for (Uint i=0; i<numUnkOnTetra; ++i)
                    _RecvBuf[_RecvBufPos+i] = buffer[buffer_pos+i].val;
                _RecvBufPos+=numUnkOnTetra;

                tp->Unknowns.SetUnkReceived(idx);
            }
            buffer_pos+=numUnkOnTetra;
        }
    }

    Assert(buffer_pos==NumberOfUnknownsOnTetra(),
           DROPSErrCL("ParMultiGridCL::ScatterUnknowns: I haven't check enough places for unknowns!"),
           DebugParallelNumC);
    return 0;
}

int ParMultiGridCL::GatherUnknownsMigV (OBJT obj, void* buf)
{
    VertexCL* const sp= ddd_cast<VertexCL*>(obj);
    TransferUnkT* buffer = static_cast<TransferUnkT*>(buf);
    int bufferpos=0;

    for (size_t index_type=0; index_type<_VecDesc.size(); ++index_type)
    {
        const Uint idx          = _VecDesc[index_type]->RowIdx->GetIdx(),
                   numUnkOnVert = _VecDesc[index_type]->RowIdx->NumUnknownsVertex();
        if (sp->Unknowns.Exist(idx))
        {
            buffer[bufferpos].mark=true;
            buffer[bufferpos].idx=idx;
            const Uint sysnum= sp->Unknowns(idx);
            if (!sp->Unknowns.UnkReceived(idx))
                for (Uint i=0; i<numUnkOnVert; ++i)
                    buffer[bufferpos++].val= _VecDesc[index_type]->Data[sysnum+i];
            else //sp->Unknowns.Unkreceived(idx)
                for (Uint i=0; i<numUnkOnVert; ++i)
                    buffer[bufferpos++].val= _RecvBuf[sysnum+i];
        }
        else
        {
            buffer[bufferpos].mark=false;
            bufferpos+=numUnkOnVert;
        }
    }
    return 0;
}

int ParMultiGridCL::ScatterUnknownsMigV(OBJT obj, void* buf)
{
    VertexCL* const sp= ddd_cast<VertexCL*>(obj);
    TransferUnkT* buffer = static_cast<TransferUnkT*>(buf);
    int bufferpos=0;

    for (size_t index_type=0; index_type<_VecDesc.size(); ++index_type)
    {
        const Uint idx          = _VecDesc[index_type]->RowIdx->GetIdx(),
                   numUnkOnVert = _VecDesc[index_type]->RowIdx->NumUnknownsVertex();
        if (sp->Unknowns.Exist(idx))        // do nothing, unknowns allready known
            bufferpos += numUnkOnVert;
        else if (buffer[bufferpos].mark && buffer[bufferpos].idx==idx)    // new unknowns has been sent to right index
        {
            if (numUnkOnVert+_RecvBufPos>_RecvBuf.size())
                EnlargeReceiveBuffer();

            sp->Unknowns.Prepare(idx);
            sp->Unknowns(idx)=_RecvBufPos;

            for (Uint i=0; i<numUnkOnVert; ++i)
                _RecvBuf[_RecvBufPos+i]=buffer[bufferpos++].val;
            _RecvBufPos+=numUnkOnVert;

            sp->Unknowns.SetUnkReceived(idx);
        }
        else                                // not known unknowns and no information received
            bufferpos+=numUnkOnVert;
    }
    return 0;
}

int ParMultiGridCL::GatherUnknownsMigE (OBJT obj, void* buf)
{
    EdgeCL* const sp= ddd_cast<EdgeCL*>(obj);
    TransferUnkT* buffer = static_cast<TransferUnkT*>(buf);
    int bufferpos=0;
    for (size_t index_type=0; index_type<_VecDesc.size(); ++index_type)
    {
        const Uint idx          = _VecDesc[index_type]->RowIdx->GetIdx(),
                   numUnkOnEdge = _VecDesc[index_type]->RowIdx->NumUnknownsEdge();
        if (!numUnkOnEdge)
            continue;
        if (sp->Unknowns.Exist(idx))
        {
            buffer[bufferpos].mark=true;
            buffer[bufferpos].idx=idx;
            const Uint sysnum= sp->Unknowns(idx);
            if (!sp->Unknowns.UnkReceived(idx))
                for (Uint i=0; i<numUnkOnEdge; ++i)
                    buffer[bufferpos++].val= _VecDesc[index_type]->Data[sysnum+i];
            else //sp->Unknowns.Unkreceived(idx)
                for (Uint i=0; i<numUnkOnEdge; ++i)
                    buffer[bufferpos++].val= _RecvBuf[sysnum+i];
        }
        else
        {
            buffer[bufferpos].mark=false;
            buffer[bufferpos].idx=idx;
            bufferpos+=numUnkOnEdge;
        }
    }
    return 0;
}

int ParMultiGridCL::ScatterUnknownsMigE(OBJT obj, void* buf)
{
    EdgeCL* const sp= ddd_cast<EdgeCL*>(obj);
    TransferUnkT* buffer = static_cast<TransferUnkT*>(buf);
    int bufferpos=0;

    for (size_t index_type=0; index_type<_VecDesc.size(); ++index_type)
    {
        const Uint idx          = _VecDesc[index_type]->RowIdx->GetIdx(),
                   numUnkOnEdge = _VecDesc[index_type]->RowIdx->NumUnknownsEdge();
        if (!numUnkOnEdge)
            continue;
        if (sp->Unknowns.Exist(idx))        // do nothing, unknowns allready known
            bufferpos += numUnkOnEdge;
        else if (buffer[bufferpos].mark && buffer[bufferpos].idx==idx)    // new unknowns has been sent to right index
        {
            if (numUnkOnEdge+_RecvBufPos>_RecvBuf.size())
                EnlargeReceiveBuffer();
            sp->Unknowns.Prepare(idx);
            sp->Unknowns(idx)=_RecvBufPos;
            sp->Unknowns.SetUnkReceived(idx);

            for (Uint i=0; i<numUnkOnEdge; ++i)
                _RecvBuf[_RecvBufPos+i]=buffer[bufferpos++].val;
            _RecvBufPos+=numUnkOnEdge;
        }
        else                                // not known unknowns and no information received
            bufferpos+=numUnkOnEdge;
    }
    return 0;
}


/// \brief Clear all XferNew-Marks
void ParMultiGridCL::DelAllUnkRecv()
{
    for (MultiGridCL::const_VertexIterator sit=mg_->GetAllVertexBegin(); sit!=mg_->GetAllVertexEnd(); ++sit){
        sit->Unknowns.ResetUnkReceived();
    }
    for (MultiGridCL::const_EdgeIterator sit=mg_->GetAllEdgeBegin(); sit!=mg_->GetAllEdgeEnd(); ++sit){
        sit->Unknowns.ResetUnkReceived();
    }
    for (MultiGridCL::const_TetraIterator sit=mg_->GetAllTetraBegin(); sit!=mg_->GetAllTetraEnd(); ++sit){
        sit->Unknowns.ResetUnkReceived();
    }
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


/// \brief This function assures, that every proc has the same number of level. This is important in the refinement algorithm
/** Make sure, the containers are modifiable*/
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
    bool Gather( DiST::TransferableCL& t, DiST::Helper::SendStreamCL& send)
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

    bool Scatter( DiST::TransferableCL& t, const size_t& numData, DiST::Helper::MPIistreamCL& recv)
    {
        DiST::Helper::GeomIdCL gid;
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
        DROPS_STD_UNORDERED_MAP<DiST::Helper::GeomIdCL,short int,DiST::Helper::Hashing> edgeMFR;

        for (MultiGridCL::const_TetraIterator it= mg_->GetTetrasBegin(Level), end= mg_->GetTetrasEnd(); it!=end; ++it)
            if (it->IsRegularlyRef() && it->IsMaster()) {
                // has marked its edges' MFR, so do the same here
                for (int e=0; e<6; ++e)
                    edgeMFR[ it->GetEdge(e)->GetGID()]++;
            }
        // now compare with edges' local MFRs
        for (MultiGridCL::const_EdgeIterator it= mg_->GetEdgesBegin(Level), end= mg_->GetEdgesEnd(); it!=end; ++it) {
            const DiST::Helper::GeomIdCL gid= it->GetGID();
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
    bool Gather( DiST::TransferableCL& t, DiST::Helper::SendStreamCL& send)
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
    bool Scatter( DiST::TransferableCL& t, const size_t& numData, DiST::Helper::MPIistreamCL& recv)
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
    bool Gather( DiST::TransferableCL& t, DiST::Helper::SendStreamCL& send)
    {
        TetraCL* tp; simplex_cast( t, tp);

        Assert( t.GetPrio()==PrioGhost, DROPSErrCL("HandlerRefMarkCL::Gather: called for non-ghost"), DebugParallelC);
        send << tp->IsMarkedForRegRef();
        return true;
    }
    /// \brief This is called by the master copy. The corresponding ghost copy definitely put in the
    ///     message, whether the tetrahedron is marked for regular refinement
    bool Scatter( DiST::TransferableCL& t, __UNUSED__ const size_t& numData, DiST::Helper::MPIistreamCL& recv)
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
            const DiST::Helper::GeomIdCL midVertGID( e.GetLevel()+1, e.GetGID().bary, DiST::GetDim<VertexCL>());
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

    // All tetras, that are marked for removement, should be removed now!
    // All subsimplices are marked for removement
    for (int l=0; l<=level_; ++l)
    {
//        mg_->Tetras_[l].remove_if( std::mem_fun_ref(&TetraCL::IsMarkedForRemovement) );
//        for_each( mg_->GetVerticesBegin(l), mg_->GetVerticesEnd(l), std::mem_fun_ref( &VertexCL::SetRemoveMark ) );
//        for_each( mg_->GetEdgesBegin(l), mg_->GetEdgesEnd(l), std::mem_fun_ref( &EdgeCL::SetRemoveMark ) );
//        for_each( mg_->GetFacesBegin(l), mg_->GetFacesEnd(l), std::mem_fun_ref( &FaceCL::SetRemoveMark ) );
    }

    ModifyBegin();

    Comment("  * Rescue HasGhosts"<<std::endl,DebugParallelC);
    TreatHasGhosts();       // mark subs of HasGhosts as VGhost, rescue subs
    Comment("  * Rescue Ghosts"<<std::endl,DebugParallelC);
    TreatGhosts();          // rescue subs of Ghosts, mark as Ghost
    Comment("  * Rescue Masters"<<std::endl,DebugParallelC);
    RescueMasterCL rescueMasterSubs(-1, *modify_);
    rescueMasterSubs.Call();
    // Remove link to parent, if tetra is ghost (note: already done in TreatGhosts)
//    for (int l=0; l<=level_; ++l){
//        for (MultiGridCL::TetraIterator sit(mg_->Tetras_[l].begin()); sit!=mg_->Tetras_[l].end(); ++sit){
//            if (sit->IsGhost() && sit->GetParent() ){
//                sit->Parent_=0;
//            }
//        }
//    }

//    Comment("  * Tell DiST to delete Objects!"<<std::endl,DebugParallelC);
//    for (MultiGridCL::VertexIterator sit( mg_->GetAllVertexBegin()); sit!=mg_->GetAllVertexEnd(level_); ++sit)
//        if ( sit->IsMarkedForRemovement())
//            modify_->Delete( *sit);
//    for (MultiGridCL::EdgeIterator sit( mg_->GetAllEdgeBegin()); sit!=mg_->GetAllEdgeEnd(level_); ++sit)
//        if ( sit->IsMarkedForRemovement())
//            modify_->Delete( *sit);
//    for (MultiGridCL::FaceIterator sit( mg_->GetAllFaceBegin()); sit!=mg_->GetAllFaceEnd(level_); ++sit)
//        if ( sit->IsMarkedForRemovement())
//            modify_->Delete( *sit);
    ModifyEnd();

    // Accumulate Ref-counter on edges
    Comment("  * Accumulate MFR"<<std::endl,DebugParallelC);
    AccumulateMFR();

    // Adapt midvertex pointers on edges
    Comment("  * Adapting Midvertex on Edges"<<std::endl,DebugParallelC);
    AdaptMidVertex();

    // Rescue unknowns on edges, that are deleted and midvertex stays on processor
    /// \todo DiST: Implement rescue of dof!
/*    if (VecDescRecv()){
        Comment("  * Send unknowns "<<std::endl,DebugParallelC);
        DynamicDataInterfaceCL::IFExchange(AllSimplexIFCL<VertexCL>::GetIF(),               // exchange datas over distributed vertices
                       NumberOfUnknownsOnVertex()* sizeof(TransferUnkT),// number of datas to be exchanged
                       GatherUnknownsMigVC,                             // how to gather datas
                       ScatterUnknownsMigVC                             // how to scatter datas
                      );
        DynamicDataInterfaceCL::IFExchange(AllSimplexIFCL<EdgeCL>::GetIF(),                 // exchange datas over distributed vertices
                       NumberOfUnknownsOnEdge()* sizeof(TransferUnkT),  // number of datas to be exchanged
                       GatherUnknownsMigEC,                             // how to gather datas
                       ScatterUnknownsMigEC                             // how to scatter datas
                      );
        RescueUnknownsOnEdges();
    }
*/

    // now physically delete the simplices with RemoveMark from memory
//    Comment("- Delete all unused Simplices!"<<std::endl,DebugParallelC);
//    mg_->PrepareModify();
//    for (int l=0; l<=level_; ++l)
//    {
//        mg_->Vertices_[l].remove_if( std::mem_fun_ref(&VertexCL::IsMarkedForRemovement) );
//        mg_->Edges_[l].remove_if( std::mem_fun_ref(&EdgeCL::IsMarkedForRemovement) );
//        mg_->Faces_[l].remove_if( std::mem_fun_ref(&FaceCL::IsMarkedForRemovement) );
//    }
//
//    level_=-1; // invalidate level_
//    mg_->FinalizeModify();
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

/// \brief Destroy unknowns on non-master vertices, edges and tetras if there are information about them
void ParMultiGridCL::DeleteUnksOnGhosts(int Level)
/** This procedure iterates over all simplices, that stores unknowns and delete them.
    \param Level level of the simplices (default all levels)
*/
{
    if (_UnkOnSimplex[0])
        for (MultiGridCL::VertexIterator sit(mg_->GetAllVertexBegin(Level)), end(mg_->GetAllVertexEnd(Level)); sit!=end; ++sit)
            if (!sit->IsMaster() && sit->Unknowns.Exist())
                sit->Unknowns.Destroy();
    if (_UnkOnSimplex[1])
        for (MultiGridCL::EdgeIterator sit(mg_->GetAllEdgeBegin(Level)), end(mg_->GetAllEdgeEnd(Level)); sit!=end; ++sit)
            if (!sit->IsMaster() && sit->Unknowns.Exist())
                sit->Unknowns.Destroy();
    if (_UnkOnSimplex[2])
        for (MultiGridCL::TetraIterator sit(mg_->GetAllTetraBegin(Level)), end(mg_->GetAllTetraEnd(Level)); sit!=end; ++sit)
            if (!sit->IsMaster() && sit->Unknowns.Exist())
                sit->Unknowns.Destroy();


/*    if (_VecDescRecv){
        if (Level==-1)
        {
            DynamicDataInterfaceCL::IFExecLocal(NotMasterIF_, &ExecDestroyUnksVC);
            DynamicDataInterfaceCL::IFExecLocal(NotMasterIF_, &ExecDestroyUnksEC);
        }
        else
        {
            DynamicDataInterfaceCL::IFAExecLocal(NotMasterIF_, Level, &ExecDestroyUnksVC);
            DynamicDataInterfaceCL::IFAExecLocal(NotMasterIF_, Level, &ExecDestroyUnksEC);
        }
    }*/
}

template<class SimplexT>
  int ParMultiGridCL::DestroyUnksOnSimplex(OBJT obj)
{
    SimplexT *sp= ddd_cast<SimplexT*>(obj);

    if (sp->Unknowns.Exist())
    {
        sp->Unknowns.Destroy();
        return 0;
    }
    else
        return 1;
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

/****************************************************************************
* D E C L A R E A L L                                                            *
*****************************************************************************
*   This procedure declares the TypeTs: Vertex, Edge, Face, Tetra,        *
*   AddedScal, AddedVec, BndPtr, ChildPtr                                        *
****************************************************************************/
/// \brief Declare all DDD-Types
/// \todo DiST: Remove this function
/*
void ParMultiGridCL::DeclareAll()
{

    VertexCL::Declare();
    EdgeCL::Declare();
    FaceCL::Declare();
    TetraCL::Declare();
    AddedScalCL::Declare();
    AddedVecCL::Declare();
    DeclareBndPtT();
    DeclareChildPtrT();
    Comment("- All Types are declared" << std::endl, DebugParallelC);
}
*/

/****************************************************************************
* D E F I N E A L L                                                             *
*****************************************************************************
*   This procedure defines the TypeTs: Vertex, Edge, Face, Tetra,         *
*   AddedScal, AddedVec, BndPtr, ChildPtr                                        *
****************************************************************************/
/// \brief Define all DDD-Types
/// \todo DiST: Remove this function
/*
void ParMultiGridCL::DefineAll()
{
    VertexCL::Define();
    EdgeCL::Define();
    FaceCL::Define();
    TetraCL::Define();
    AddedScalCL::Define();
    AddedVecCL::Define();
    DefineBndPtT();
    DefineChildPtrT();
    Comment("- All Types are defined" << std::endl, DebugParallelC);
}
*/

/****************************************************************************
* I N I T  I F                                                                      *
*****************************************************************************
* This procedure definies the Interfaces for: Edges, Faces, Tetras          *
* The interfaces for numerical accumulations are also defined here. So the  *
* user does not have to worry about that!                                       *
****************************************************************************/
/// \brief Initialize all Interfaces
/// \todo DiST: Remove this function
/*
void ParMultiGridCL::InitIF()
{
	TypeT  O[4];
    O[0]= VertexCL::GetType();
    O[1]= EdgeCL::GetType();
    O[2]= TetraCL::GetType();
    O[3]= FaceCL::GetType();

    PrioT  A[5], B[4];   // arrays of priorities

    A[0]= B[0] = PrioHasUnk;
    A[1]= B[1] = PrioMaster;
    A[2]= B[2] = PrioVGhost;
    A[3]= B[3] = PrioGhost;
    A[4]=        PrioKilledGhost;

    // interface of edges
    Assert(!_EdgeIF, DROPSErrCL("ParMultiGridCL: InitIF: EdgeIF allready declared"), DebugParallelC);
    _EdgeIF = DynamicDataInterfaceCL::IFDefine(1, O+1, 4, A, 4, B);               // Master, VGhost, Ghost -> Master, VGhost, Ghost
    DynamicDataInterfaceCL::IFSetName( _EdgeIF, (char*)"Edge-IF for AccumulateMFR");

    // interface of faces
    Assert(!_FaceIF, DROPSErrCL("ParMultiGridCL: InitIF: FaceIF allready declared"), DebugParallelC);
    _FaceIF = DynamicDataInterfaceCL::IFDefine(1, O+3, 2, A, 2, B);               // Master Face -> Master Face
    DynamicDataInterfaceCL::IFSetName( _FaceIF, (char*)"Face-IF for Sanity Check");

    // interface of tetras
    Assert(!_TetraIF, DROPSErrCL("ParMultiGridCL: InitIF: TetraIF allready declared"), DebugParallelC);
    _TetraIF = DynamicDataInterfaceCL::IFDefine(1, O+2, 1, A+3, 2, B);    // Ghost Tetra -> Master Tetra
    DynamicDataInterfaceCL::IFSetName( _TetraIF, (char*)"Tetra-IF for RestrictMarks");

    // interface of not master vertices, edges and tetras
    Assert(!NotMasterIF_, DROPSErrCL("ParMultiGridCL: InitIF: NotMasterV_IF allready declared"), DebugParallelC);
    NotMasterIF_= DynamicDataInterfaceCL::IFDefine(3, O, 2, A+2, 2, B+2);   // PrioVGhost, PrioGhost -> PrioVGhost, PrioGhost
    DynamicDataInterfaceCL::IFSetName( NotMasterIF_, (char*)"non master Vertex, Edge, Tetrahedron-IF");

    // Ghost tetras to master-tetras
    _GhToMaTetra_IF= DynamicDataInterfaceCL::IFDefine(1, O+2, 1, A+4, 2, B);    // PrioKilledGhost -> PrioMaster
    DynamicDataInterfaceCL::IFSetName( _GhToMaTetra_IF, (char*)"Killed Ghost to Master Interface");

    // Init the Interfaces for numerical accumulations!
    InterfaceCL<VertexCL>::InitIF();
    InterfaceCL<EdgeCL>::InitIF();
    InterfaceCL<TetraCL>::InitIF();
    AllSimplexIFCL<VertexCL>::InitIF();
    AllSimplexIFCL<EdgeCL>::InitIF();
}
*/

/****************************************************************************
* < S i m p l e x > X f e r                                                 *
*****************************************************************************
* These procedures transfer a simplex from the calling proc to another      *
* proc. XferBegin() and XferEnd() have to call before and after respectively*
****************************************************************************/
/*
/// \brief Send a Vertex
void ParMultiGridCL::VXfer(VertexCL &v, PROCT dest, PrioT prio, bool del)
{
    // Set the removement-mark before calling the transfer command, because this mark is checked, whether additional data is send too!
    if (del) v.SetRemoveMark();
    const HDRT hdr= const_cast<HDRT>(&v._dddH);
    DynamicDataInterfaceCL::XferCopyObj( hdr, dest, prio);
}

/// \brief Send a Edge
void ParMultiGridCL::EXfer(EdgeCL &e, PROCT dest, PrioT prio, bool del)
{
    if (del)
        e.SetRemoveMark();

    const HDRT hdr= const_cast<HDRT>(&e._dddH);
    DynamicDataInterfaceCL::XferCopyObj( hdr, dest, prio);

}

/// \brief Send a Face
void ParMultiGridCL::FXfer(FaceCL &f, PROCT dest, PrioT prio, bool del)
{
    if (del)
        f.SetRemoveMark();

    const HDRT hdr= const_cast<HDRT>(&f._dddH);
    DynamicDataInterfaceCL::XferCopyObj(hdr, dest, prio);
    // del richtig beruecksichtigen!
}
*/

/// \brief Send a Tetra
void ParMultiGridCL::Transfer(TetraCL &t, int dest, Priority prio, bool del)
{
    transfer_->Transfer( t, dest, prio, del);
    if (t.IsRegularlyRef() && t.IsMaster() && prio==PrioMaster)
        t.UnCommitRegRefMark();
}

/****************************************************************************
* H A N D L E R S   F O R   D D D                                           *
*****************************************************************************
* The following handlers tell the DDD-System how to tread the DROPS classes *
* if they are touched by Xfer-Commands or by Interface-Commands             *
****************************************************************************/

//  Deleting an object by DDD
//----------------------------
/// \brief Handle delete of an object
///
/// Set the remove mark and desctuct the DDD-Hdr. The real delete of the object is done within
/// the code.
//template<class SimplexT> void ParMultiGridCL::HandlerDelete( OBJT obj)
//{
//
//    SimplexT* const sp= ddd_cast<SimplexT*>(obj);
//    sp->SetRemoveMark();
//    //  Ich denke, dass man hier den Hdr_Destructor nicht braucht, da er automatisch von DDD aufgerufen wird
//    // ---> Doch, den braucht man ganz ganz unbedingt! (siehe Hdr_Destructor bei DDD!)
//    DynamicDataInterfaceCL::HdrDestructor( &sp->_dddH);
//    AllComment("  * Simplex with GID "<<sp->GetGID()<<" deleted on proc "<<ProcCL::MyRank()<<std::endl,DebugParallelHardC);
//}

//  Construct an object by DDD
//-----------------------------
/// \todo DiST: Do we need the handlers for construction and sending; of doesn't think so
/*
/// \brief Construct a vertex and return it as a DDD-Object
OBJT ParMultiGridCL::HandlerVConstructor( size_t, PrioT, ATTRT level){
    (*_VertCont)[level].push_back( VertexCL());
    return ddd_cast(&(*_VertCont)[level].back());
}

/// \brief Construct a edge and return it as a DDD-Object
OBJT ParMultiGridCL::HandlerEConstructor( size_t, PrioT, ATTRT level){
    (*_EdgeCont)[level].push_back( EdgeCL());
    return ddd_cast(&(*_EdgeCont)[level].back());
}

/// \brief Construct a face and return it as a DDD-Object
OBJT ParMultiGridCL::HandlerFConstructor( size_t, PrioT, ATTRT level){
    (*_FaceCont)[level].push_back( FaceCL());
    return ddd_cast(&(*_FaceCont)[level].back());
}

/// \brief Construct a tera and return it as a DDD-Object
OBJT ParMultiGridCL::HandlerTConstructor( size_t, PrioT, ATTRT level){
    (*_TetraCont)[level].push_back( TetraCL());
    return ddd_cast(&(*_TetraCont)[level].back());
}
*/

//  Transfer an object by DDD
//----------------------------
/// \brief transfer a vertex
///
/// Check if numerical data have to be transfered too. If this case happens, tell DDD that additional
/// data will be transfered. Than DDD calls the Gather and Scatter functions <p>
/// The Recylce-Bin is also destroyed and boundary information are send too
/*
void ParMultiGridCL::HandlerVXfer( OBJT obj, __UNUSED__ PROCT proc, __UNUSED__ PrioT prio)
{
    VertexCL* const vp= ddd_cast<VertexCL*>(obj);
    int numSendScalUnk= 0,
        numSendVecUnk = 0;

    vp->DestroyRecycleBin();

    // if there are bondary-information, then send them too!
    if (vp->IsOnBoundary() )
    	DynamicDataInterfaceCL::XferAddData( vp->_BndVerts->size(), _BndPtT);

    // if there this ParMultiGridCL knowns about unknowns, ther are unknwons on this vertex and this vertex is marked for removement, count the unknowns and give this Information to DDD!
    if (VecDescRecv() && vp->Unknowns.Exist())
    {
        for (size_t i=0; i<_VecDesc.size(); ++i)
        {
            if ((_VecDesc[i]) && (_VecDesc[i])->RowIdx &&  vp->Unknowns.Exist((_VecDesc[i])->RowIdx->GetIdx()) )
            {
                Assert((_VecDesc[i])->RowIdx->NumUnknownsVertex() ==1 || (_VecDesc[i])->RowIdx->NumUnknownsVertex()==3,
                        DROPSErrCL("ParMultiGridCL::HandlerVXfer: Can only send scalar or 3d-vectors"), DebugParallelC);
                if ( (_VecDesc[i])->RowIdx->NumUnknownsVertex() ==1 )
                    ++numSendScalUnk;
                if( (_VecDesc[i])->RowIdx->NumUnknownsVertex() ==3 )
                    ++numSendVecUnk;
            }
        }
        DynamicDataInterfaceCL::XferAddData( numSendScalUnk, AddedScalCL::_dddT );
        DynamicDataInterfaceCL::XferAddData( numSendVecUnk,  AddedVecCL::_dddT  );
    }

    AllComment("  * Vertex with GID " << vp->GetGID() << " from " << ProcCL::MyRank() << " to " << proc << " "
            << "with " << ( !(vp->IsOnBoundary()) ? 0 : vp->_BndVerts->size()) << " boundary points and "
            << "with " << numSendScalUnk << " scalar unknowns and "
            << "with " << numSendVecUnk  << " vectoriel unknowns!" << std::endl, DebugParallelHardC);
}
*/

/// \brief transfer an edge
/** Check if numerical data have to be transfered too. If this case happens, tell DDD that additional
    data will be transfered. Than DDD calls the Gather and Scatter functions <p>*/
/*
void ParMultiGridCL::HandlerEXfer(OBJT obj, __UNUSED__ PROCT proc, PrioT)
{
    EdgeCL* const ep= ddd_cast<EdgeCL*>(obj);

    int numSendScalUnk= 0,
        numSendVecUnk = 0;

    // if ParMultiGridCL knowns about unknowns and there are unknowns on this
    // edge count the unknowns and give this information to DDD
    if (VecDescRecv() && ep->Unknowns.Exist())
    {
        for (size_t i=0; i<_VecDesc.size(); ++i)
        {
            if ((_VecDesc[i]) && (_VecDesc[i])->RowIdx &&  ep->Unknowns.Exist((_VecDesc[i])->RowIdx->GetIdx()) )
            {
                Assert((_VecDesc[i])->RowIdx->NumUnknownsEdge()==1 || (_VecDesc[i])->RowIdx->NumUnknownsEdge()==3,
                        DROPSErrCL("ParMultiGridCL::HandlerVXfer: Can only send scalar or 3d-vectors"), DebugParallelC);

                if ( (_VecDesc[i])->RowIdx->NumUnknownsEdge() ==1 )
                    ++numSendScalUnk;
                if ( (_VecDesc[i])->RowIdx->NumUnknownsEdge() ==3 )
                            ++numSendVecUnk;
            }
        }
        DynamicDataInterfaceCL::XferAddData( numSendScalUnk, AddedScalCL::_dddT );
        DynamicDataInterfaceCL::XferAddData( numSendVecUnk,  AddedVecCL::_dddT  );
    }

    AllComment("  * Edge with GID " << ep->GetGID() << " from " << ProcCL::MyRank() << " to " << proc << " "
            << "with " << numSendScalUnk << " scalar unknowns and "
            << "with " << numSendVecUnk  << " vectoriel unknowns!" << std::endl, DebugParallelHardC);
}
*/

/// \brief transfer a face
/** Nothing is done. <p>*/
/*
void ParMultiGridCL::HandlerFXfer(__UNUSED__ OBJT obj, __UNUSED__ PROCT proc, PrioT)
{
    AllComment("  * Face with GID " << ddd_cast<FaceCL*>(obj)->GetGID() << " from " << ProcCL::MyRank() << " to " << proc << std::endl, DebugParallelHardC);
}
*/
/// \brief transfer a tetra
/** Transfer a tetraeder with all subsimplices. Also pointer to childern are transfered and
    if neccessary numerical data.*/
/*
void ParMultiGridCL::HandlerTXfer(OBJT obj, PROCT proc, PrioT prio)
{
    TetraCL* const tp= ddd_cast<TetraCL*>(obj);
    int numChilds=0;
    if (tp->_Children && !tp->HasGhost() )
    {
    	DynamicDataInterfaceCL::XferAddData( tp->GetRefData().ChildNum, _ChildPtrT);
        numChilds = tp->GetRefData().ChildNum;
    }

    int numSendScalUnk= 0,
        numSendVecUnk = 0;

    // if there this ParMultiGridCL knowns about unknowns, ther are unknwons on this tetra, count the unknowns and give this Information to DDD!
    if (VecDescRecv() && tp->Unknowns.Exist())
    {
        for (size_t i=0; i<_VecDesc.size(); ++i)
        {
            if ((_VecDesc[i]) && (_VecDesc[i])->RowIdx &&  tp->Unknowns.Exist((_VecDesc[i])->RowIdx->GetIdx()) )
            {
                Assert((_VecDesc[i])->RowIdx->NumUnknownsTetra() ==1 || (_VecDesc[i])->RowIdx->NumUnknownsTetra()==3,
                        DROPSErrCL("ParMultiGridCL::HandlerTXfer: Can only send scalar or 3d-vectors"), DebugParallelC);

                if ( (_VecDesc[i])->RowIdx->NumUnknownsTetra()==1 )
                    ++numSendScalUnk;
                if ( (_VecDesc[i])->RowIdx->NumUnknownsTetra()==3 )
                    ++numSendVecUnk;
            }
        }
        DynamicDataInterfaceCL::XferAddData( numSendScalUnk, AddedScalCL::_dddT );
        DynamicDataInterfaceCL::XferAddData( numSendVecUnk,  AddedVecCL::_dddT  );
    }

    // all subsimplices and of the tetra are transferred too.
    for (Uint i=0; i<NumVertsC; ++i)
    	DynamicDataInterfaceCL::XferCopyObj( tp->_Vertices[i]->GetHdr(), proc, prio);
    for (Uint i=0; i<NumEdgesC; ++i)
    	DynamicDataInterfaceCL::XferCopyObj( tp->_Edges[i]->GetHdr(), proc, prio);
    for (Uint i=0; i<NumFacesC; ++i)
    	DynamicDataInterfaceCL::XferCopyObj( tp->_Faces[i]->GetHdr(), proc, prio);

    AllComment("  * Tetra with GID " << tp->GetGID() << " from " << ProcCL::MyRank() << " to " << proc
            << " with " << numSendScalUnk << " scalar unknowns,"
            << " with " << numSendVecUnk << " vectoriel unknowns"
            << " and " << numChilds << " Children" << std::endl, DebugParallelHardC);
}
*/

//  Gather and Scatter-Handlers
//-----------------------------

/// \brief Send additional data with the simplices
/** These procedures put additional data to a message or receive this data.
    This data might be geometrical or numerical, like boundary-information or children-information,
    or the numerical values onto the simplex. */
/*
void ParMultiGridCL::HandlerVGather( OBJT obj, int cnt, TypeT type, void* buf)
{
    VertexCL* const vp= ddd_cast<VertexCL*>(obj);
    Assert(type==AddedScalCL::GetType() || type==AddedVecCL::GetType() || type==_BndPtT,
           DROPSErrCL("ParMultiGridCL: HandlerVGather: Cannot handle this type!"), DebugParallelC);
    // if this vert is on a boundary, add the boundary-points to the message
    if (type == _BndPtT)
    {
        // the buffer is a storage of boundary-pointers
        BndPointCL* buffer= static_cast<BndPointCL*>(buf);
        // put all boundary vertices into the buffer
        for( VertexCL::const_BndVertIt it= vp->GetBndVertBegin(), end= vp->GetBndVertEnd(); it!=end; ++it, ++buffer)
        {
            //std::cout << it->GetCoord2D() << std::endl;
            *buffer= *it;
        }
    }
    // if there are Unknowns on this vert, add these to the message
    else if ( VecDescRecv() && (type == AddedScalCL::GetType() || type == AddedVecCL::GetType()) )
        SendUnknowns(vp,type,buf,cnt);
}
*/
/*
void ParMultiGridCL::HandlerVScatter( OBJT obj, int cnt, TypeT type, void* buf, int newness)
{
    VertexCL* const vp= ddd_cast<VertexCL*>(obj);
    Assert(type==AddedScalCL::GetType() || type==AddedVecCL::GetType() || type==_BndPtT ,
           DROPSErrCL("ParMultiGridCL: HandlerVScatter: Cannot handle this type!"), DebugParallelC);
    // if boundary information are received
    if ( type == _BndPtT )
    {
        const BndPointCL* const buffer= static_cast<BndPointCL*>(buf);
//         vp->_BndVerts= new std::vector<BndPointCL>;         // allocate memory for the boundary-points
        for( int i=0; i<cnt; ++i)
            if (!vp->HasBnd(buffer[i]))
                vp->AddBnd(buffer[i]);
//             vp->_BndVerts->push_back( buffer[i]);           // store the received boundary-points
    }

    // if numerical data are received
    else if ( VecDescRecv() && (type == AddedScalCL::GetType() || type==AddedVecCL::GetType()) )
        RecvUnknowns(vp,type,buf,cnt);
}
*/

/*
void ParMultiGridCL::HandlerEGather( OBJT obj, int cnt, TypeT type, void* buf)
{
    EdgeCL* const ep= ddd_cast<EdgeCL*>(obj);

    Assert(type==AddedScalCL::GetType() || type==AddedVecCL::GetType(),
           DROPSErrCL("ParMultiGridCL: HandlerEGather: Cannot handle this type!"), DebugParallelC);
    if ( VecDescRecv() && (type == AddedScalCL::GetType() || type == AddedVecCL::GetType()) )
        SendUnknowns(ep,type,buf,cnt);
}
*/
/*
void ParMultiGridCL::HandlerEScatter( OBJT obj, int cnt, TypeT type, void* buf, int newness)
{
    EdgeCL* const ep= ddd_cast<EdgeCL*>(obj);
    Assert(type==AddedScalCL::GetType() || type==AddedVecCL::GetType(),
           DROPSErrCL("ParMultiGridCL: HandlerEScatter: Cannot handle this type!"), DebugParallelC);

    if ( VecDescRecv() && (type == AddedScalCL::GetType() || type==AddedVecCL::GetType()) )
        RecvUnknowns(ep,type,buf,cnt);
}
*/
/// \brief Collect additional data for a tetra transfer
/** The additional data may be pointer to children or numerical data within a migration.*/
/*
void ParMultiGridCL::HandlerTGather( OBJT obj, int cnt, TypeT type, void* buf)
{
    TetraCL* const tp= ddd_cast<TetraCL*>(obj);
    Assert(type==_ChildPtrT || type==AddedScalCL::GetType() || type==AddedVecCL::GetType(),
           DROPSErrCL("ParMultiGridCL: HandlerTGather: Cannot handle this type!"), DebugParallelC);
    if (type==_ChildPtrT)           // if this Tetra should send his children
    {
        TetraCL** buffer= static_cast<TetraCL**>(buf);
        for (TetraCL::ChildPIterator it(tp->GetChildBegin()), end(tp->GetChildEnd()); it!=end; ++it, ++buffer)
        {
            *buffer= *it;
        }
    }
    else if ( VecDescRecv() && (type == AddedScalCL::GetType() || type == AddedVecCL::GetType()) )
        SendUnknowns(tp,type,buf,cnt);
}
*/

/// \brief receive additional data for a tetra transfer
/*
void ParMultiGridCL::HandlerTScatter( OBJT obj, int cnt, TypeT type, void* buf, int newness)
{
    TetraCL* const tp= ddd_cast<TetraCL*>(obj);
    Assert(type==_ChildPtrT || type==AddedScalCL::GetType() || type==AddedVecCL::GetType(),
           DROPSErrCL("ParMultiGridCL: HandlerTScatter: Cannot handle this type!"), DebugParallelC);
    if (newness!=XFER_NEW && tp->IsGhost() )
    {
        // This case shouldn't happen, because this case is catched by the migrate function in LoadBalCL
        std::cerr << ">>>>["<<ProcCL::MyRank()<<"] HandlerScatterTetra: ("<<tp->GetGID()<<") Master replaced by Ghost, continuing anyway, workaround enabled! Seems to be an illegal xfer!!!!"<<std::endl;
        tp->GetHdr()->prio= PrioMaster;
    }

    if (type==_ChildPtrT)                           // hiphiphurra, I am mother!
    {
        TetraCL** const buffer= static_cast<TetraCL**>(buf);

        // For Problem Problems/prob2_sun_2_procs.param! (Test-Case for pointer-length in DDD!)
//      if (ProcCL::MyRank()==0 && tp->GetGID()==1728){
//          std::cout << "["<<ProcCL::MyRank()<<"]  sizeof(*buffer)="<<sizeof(*buffer)<<", sizeof(buffer)="<<sizeof(buffer)<<": Tetra "<<tp->GetGID()<<" Soll folgende Kinder empfangen:\n   ";
//          for (int i=0; i<cnt; ++i)
//              std::cout << (tmp[i])->GetGID() << ",  ";
//          std::cout  << std::endl;
//      }

        // Create new children-container if necessary
        if (!tp->_Children)
            tp->_Children= new SArrayCL<TetraCL*, MaxChildrenC>;
        // put received children into children-container!
        for( int i=0; i<cnt; ++i)
            (*(tp->_Children))[i]= buffer[i];
    }
    else if ( VecDescRecv() && (type == AddedScalCL::GetType() || type==AddedVecCL::GetType()) )
        RecvUnknowns(tp,type,buf,cnt);
}
*/


/// \brief Make Edge MFR consistent
/** If a tetra is regular refined and the priority is master and the tetraeder is made on this proc within the actual
    transfer mode, then this tetra increas MFR and set AccMFR to MFR.
    XferEnd() will accumulate the MFRs.
*/
/*
void ParMultiGridCL::HandlerTObjMkCons( OBJT obj, int newness)
{
    TetraCL* const tp= ddd_cast<TetraCL*>(obj);

    if (tp->IsRegularlyRef() && tp->IsMaster() && newness==XFER_NEW)
    {
        tp->CommitRegRefMark();
        // nun wird auf allen Edges _AccMFR:=_MFR gesetzt, um Unkonsistenzen bei vorher verteilt
        // und nun nur noch lokal gespeicherten Edges zu vermeiden. Ein abschliessenden
        // AccumulateMFR in XferEnd() setzt auf den verteilt gespeicherten Edges dann die
        // richtigen _AccMFR-Werte.
        for (TetraCL::const_EdgePIterator it= tp->GetEdgesBegin(), end= tp->GetEdgesEnd(); it!=end; ++it)
            (*it)->_AccMFR= (*it)->_MFR;
    }
}
*/


/// \brief Set Prio of an tetra
/** if priority is set to ghost, the pointer to parent is set to null
    If Prio changes from Ghost to Master, the MFR-Counter on Edges must be increased.
    The accumulated MFR will be set correct within AccumulateMFR().
*/
/*
void ParMultiGridCL::HandlerTSetPrio( OBJT obj, PrioT prio)
{
    TetraCL* const tp= ddd_cast<TetraCL*>(obj);
    Assert(!(prio==PrioGhost && tp->IsMaster() && !DynamicDataInterfaceExtraCL::InfoIsLocal( tp->GetHdr() )),
             DROPSErrCL("HandlerSetPrio: illegal prio for T"),
             DebugParallelC
          );
    if (prio==PrioGhost)
        tp->_Parent= 0;
    if (prio==PrioMaster && tp->GetPrio()==PrioGhost && tp->IsRegularlyRef())
    {   // It may happen, that this routine increases the MFR on a ghost edge, that will be after the transfer
        // not distributed any more. So remember this tetra and repair the MFR on local edges of this tetra
        // within ParMultiGridCL::AccumulateMFR()
        tp->CommitRegRefMark();
        ToHandleTetra_.push_back(tp);
    }
}
*/


//  Update-Handlers
//-----------------------------
/// \brief Update tetra after recieving
/** link tetra tho faces and set prio to ghost, if no parent is found*/
/*
void ParMultiGridCL::HandlerTUpdate( OBJT obj)
{
    TetraCL* const tp= ddd_cast<TetraCL*>(obj);
    // link tetra to his faces
    // NOTE: parent has to be linked first, if this Tetra has PrioMaster, has a parent and this parent hasn't linked before
    if (tp->IsMaster() && tp->GetParent() && !tp->GetParent()->GetFace(0)->IsLinkedTo( tp->GetParent() ) )
         // link parent
        for (Uint i=0; i<NumFacesC; ++i)
            const_cast<FaceCL*>(tp->GetParent()->GetFace(i))->LinkTetra( tp->GetParent() );
    if (tp->IsMaster() && !tp->GetParent() && tp->GetLevel()!=0)
    {
        PrioChange(tp, PrioGhost);
        // This case happens within the migration phase. Imagine the following situation:
        //     Proc A has a master copy of a tetra that should be moved to another proc (tetra is
        //         marked for removement on this Proc A)
        //     Proc B moves a ghost copy of the tetra to proc A.
        // Now proc A does not delete the tetrahedron, because it should be received by proc B. The both priorities
        // (Ghost and Master) are merged. And Master is wrongly the winner. This is corrected here!
//         AllComment("["<<ProcCL::MyRank()<<"] ====> Set Prio of Tetra " << tp->GetGID() <<" to Ghost, because this tetra has no parent on this proc! (of: This should be OK!)"<<std::endl, ~0);
    }

    // now link this Tetra to his faces
    AllComment("  * Tetra with "<<tp->GetGID()<<" is updated (linked to faces) on proc "<<ProcCL::MyRank()<<std::endl,DebugParallelHardC);
    for (Uint i=0; i<NumFacesC; ++i)
    {
        const_cast<FaceCL*>(tp->GetFace(i))->LinkTetra( tp);
        AllComment("    + to face with GID: " << const_cast<FaceCL*>(tp->GetFace(i))->GetGID()<<std::endl, DebugParallelHardC);
    }
}
*/


/****************************************************************************
* D E L E T E  O B J                                                        *
*****************************************************************************
* If DDD tries to delete Objects there occurs an error, because the STL     *
* does not allow that other code delete objects from a list, vector or so   *
* so we uses function to delete the objects!                                *
****************************************************************************/
/// \brief DDD may use this function to delete a simplex, but this must not be happend!
//void ParMultiGridCL::DeleteObj(void * /* buffer*/, size_t /*size*/, int ddd_typ)
//{
//    std::cout << "Deleting Object of type " << ddd_typ << " is still missing!" << std::endl;
//}



/****************************************************************************
* S E T H A N D L E R A L L                                                 *
*****************************************************************************
* This procedure tells DDD how to treat DROPS-stuff                         *
****************************************************************************/
/*
void ParMultiGridCL::SetAllHandler()
{
    // How to construct Simplices
	DynamicDataInterfaceCL::SetHandlerCONSTRUCTOR (VertexCL::GetType(), &HandlerVConstructorC);
	DynamicDataInterfaceCL::SetHandlerCONSTRUCTOR (EdgeCL::GetType(),   &HandlerEConstructorC);
	DynamicDataInterfaceCL::SetHandlerCONSTRUCTOR (FaceCL::GetType(),   &HandlerFConstructorC);
	DynamicDataInterfaceCL::SetHandlerCONSTRUCTOR (TetraCL::GetType(),  &HandlerTConstructorC);

    // How to delete Simplices
	DynamicDataInterfaceCL::SetHandlerDELETE(VertexCL::GetType(), &HandlerVDeleteC);
	DynamicDataInterfaceCL::SetHandlerDELETE(EdgeCL::GetType(),   &HandlerEDeleteC);
	DynamicDataInterfaceCL::SetHandlerDELETE(FaceCL::GetType(),   &HandlerFDeleteC);
	DynamicDataInterfaceCL::SetHandlerDELETE(TetraCL::GetType(),  &HandlerTDeleteC);

//  DynamicDataInterfaceCL::SetHandlerXFERDELETE(VertexCL::GetType(), &HandlerVDeleteC);
//  DynamicDataInterfaceCL::SetHandlerXFERDELETE(EdgeCL::GetType(),   &HandlerEDeleteC);
//  DynamicDataInterfaceCL::SetHandlerXFERDELETE(FaceCL::GetType(),   &HandlerFDeleteC);
//  DynamicDataInterfaceCL::SetHandlerXFERDELETE(TetraCL::GetType(),  &HandlerTDeleteC);



    // How to transfer simplices
	DynamicDataInterfaceCL::SetHandlerXFERCOPY (VertexCL::GetType(), &HandlerVXferC);
	DynamicDataInterfaceCL::SetHandlerXFERCOPY (EdgeCL::GetType(),   &HandlerEXferC);
	DynamicDataInterfaceCL::SetHandlerXFERCOPY (FaceCL::GetType(),   &HandlerFXferC);
	DynamicDataInterfaceCL::SetHandlerXFERCOPY (TetraCL::GetType(),  &HandlerTXferC);

    // How to pack data to a transfer of a simplex
	DynamicDataInterfaceCL::SetHandlerXFERGATHER(VertexCL::GetType(), &HandlerVGatherC);
	DynamicDataInterfaceCL::SetHandlerXFERGATHER(EdgeCL::GetType(),   &HandlerEGatherC);
	DynamicDataInterfaceCL::SetHandlerXFERGATHER(TetraCL::GetType(),  &HandlerTGatherC);

    // How to unpack data from a transfer of a simplex
	DynamicDataInterfaceCL::SetHandlerXFERSCATTER(VertexCL::GetType(), &HandlerVScatterC);
	DynamicDataInterfaceCL::SetHandlerXFERSCATTER(EdgeCL::GetType(),   &HandlerEScatterC);
	DynamicDataInterfaceCL::SetHandlerXFERSCATTER(TetraCL::GetType(),  &HandlerTScatterC);

    // How to handle Tetra right after the transfer!
	DynamicDataInterfaceCL::SetHandlerUPDATE     (TetraCL::GetType(), &HandlerTUpdateC);
	DynamicDataInterfaceCL::SetHandlerOBJMKCONS  (TetraCL::GetType(), &HandlerTObjMkConsC);
	DynamicDataInterfaceCL::SetHandlerSETPRIORITY(TetraCL::GetType(), &HandlerTSetPrioC);

    Comment("- All Handlers are set" << std::endl, DebugParallelC);
}
*/

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
/// \brief Show an simplex by GID
/** Search on proc (or all procs) for a simplex of given GID and show DebugInfo*/
void ParMultiGridCL::Show( const DiST::Helper::GeomIdCL& gid, char *mesg, int proc)
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
        const DiST::Helper::RemoteDataListCL& rdl=
            DiST::InfoCL::Instance().GetRemoteList( gid.dim);
        DiST::Helper::RemoteDataListCL::const_iterator it= rdl.find( gid);
        if ( it==rdl.end()){
            std::cerr << "...stored only locally." << std::endl;
        }
        else{
            DiST::Helper::RemoteDataCL::ProcList_const_iterator pit=
                it->second.GetProcListBegin();
            for ( ; pit!=it->second.GetProcListEnd(); ++pit){
                std::cerr << "...stored on proc " << pit->proc
                          << " with prio " << PriorityToString(pit->prio)
                          << std::endl;
            }
        }
    }
}


/// \brief Writes all vertices, edges, faces and tetraeders onto the ostream
/// \todo DiST: Implement me!
void ParMultiGridCL::DebugInfo(std::ostream &) const
{
/*
    os << "I have:\n";
    os << _VertCont->size() << " vertices, " << _EdgeCont->size() << " edges, " << _FaceCont->size() << " faces,"
            << _TetraCont->size() << " tetras" << std::endl;

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
*/
}

/// \brief Calc Balance over procs over tetras in the last triangulation level
/** This proc compares the highest number of tetras in the last triangulation level
    with the number of the lowest number of tetras and return the ratio
    (1 is perfect, 0 is unbalanced) */
double ParMultiGridCL::GetBalance()
{
    int myTetras = mg_->TriangTetra_.size();                                    // all procs count their tetras
    int maxTetras = ProcCL::GlobalMax(myTetras);
    int minTetras = ProcCL::GlobalMin(myTetras);

    return std::fabs(1- (double)(maxTetras)/(double)(minTetras));
}

/// \brief Display all types that are defined and declared
/*
void ParMultiGridCL::ShowTypes() const
{
	DynamicDataInterfaceCL::TypeDisplay(VertexCL::GetType());
	DynamicDataInterfaceCL::TypeDisplay(EdgeCL::GetType());
	DynamicDataInterfaceCL::TypeDisplay(FaceCL::GetType());
	DynamicDataInterfaceCL::TypeDisplay(TetraCL::GetType());
	DynamicDataInterfaceCL::TypeDisplay(AddedScalCL::GetType());
	DynamicDataInterfaceCL::TypeDisplay(AddedVecCL::GetType());
	DynamicDataInterfaceCL::TypeDisplay(_BndPtT);
	DynamicDataInterfaceCL::TypeDisplay(_ChildPtrT);
}
*/

/// \brief Display all interfaces used by ParMultiGridCL
///
/// Show all interfaces, that the DDD-System knows. This procedure uses the DDD-function DynamicDataInterfaceCL::IFDisplayAll
/*
void ParMultiGridCL::ShowInterfaces() const
{
	DynamicDataInterfaceCL::IFDisplayAll();
}
*/

/// \brief DDD-Consistency-Check
/// \todo DiST: Implement a consistency check!
/*
void ParMultiGridCL::ConsCheck()
{
	DynamicDataInterfaceCL::ConsCheck();
}
*/

/// \brief Get the size of the receive buffer
size_t ParMultiGridCL::GetRecvBufferSize()
{
    return _RecvBuf.size();
}

/****************************************************************************
* P A R A L L E L   F U N C T I O N S   F R O M   M U L T I G R I D         *
*****************************************************************************
* Declare and Define the simplices and of the Types declared in             *
* parmultigrid.h                                                            *
****************************************************************************/
/*
void VertexCL::Declare(){
    Assert(!_dddT, DROPSErrCL("VertexCL: Declare: Type allready declared"), DebugParallelC);
    _dddT = DynamicDataInterfaceCL::TypeDeclare((char*)"Vertex");
}

void EdgeCL::Declare(){
    Assert(!_dddT, DROPSErrCL("EdgeCL: Declare: Type allready declared"), DebugParallelC);
    _dddT = DynamicDataInterfaceCL::TypeDeclare((char*)"Edge");
}

void FaceCL::Declare(){
    Assert(!_dddT, DROPSErrCL("FaceCL: Declare: Type allready declared"), DebugParallelC);
    _dddT = DynamicDataInterfaceCL::TypeDeclare((char*)"Face");
}

void TetraCL::Declare(){
    Assert(!_dddT, DROPSErrCL("TetaCL: Declare: Type allready declared"), DebugParallelC);
    _dddT = DynamicDataInterfaceCL::TypeDeclare((char*)"Tetraeder");
}

void ParMultiGridCL::DeclareBndPtT(){
    Assert(!_BndPtT, DROPSErrCL("ParMultiGridCL: Declare: BndPtT-Type allready declared"), DebugParallelC);
    _BndPtT = DynamicDataInterfaceCL::TypeDeclare((char*)"Boundary-Points");
}

void ParMultiGridCL::DeclareChildPtrT(){
    Assert(!_ChildPtrT, DROPSErrCL("ParMultiGridCL: Declare: ChildPtrT-Type allready declared"), DebugParallelC);
    _ChildPtrT = DynamicDataInterfaceCL::TypeDeclare((char*)"Tetraeder-Pointer");
}

void AddedScalCL::Declare(){
    Assert(!_dddT, DROPSErrCL("AddedDataCL: Declare: Type allready declared"), DebugParallelC);
    _dddT = DynamicDataInterfaceCL::TypeDeclare((char*)"Scalar-Unknown-Transfer-Type");
}

void AddedVecCL::Declare(){
    Assert(!_dddT, DROPSErrCL("AddedDataCL: Declare: Type allready declared"), DebugParallelC);
    _dddT = DynamicDataInterfaceCL::TypeDeclare((char*)"Vector-Unknown-Transfer-Type");
}


void VertexCL::Define()
{
    VertexCL* v = 0;                                                // example VertexCL for telling DDD some information on this class
    DynamicDataInterfaceCL::TypeDefineVertex( _dddT, v,                                       // type and an example
                    EL_LDATA,  &v->_Id, sizeof(v->_Id),
                    EL_GDATA,  &v->_Coord, sizeof(v->_Coord),
                    EL_LDATA,  &v->_BndVerts, sizeof(v->_BndVerts),
                    EL_LDATA,  &v->_RemoveMark, sizeof(v->_RemoveMark),
                    EL_DDDHDR, &v->_dddH,                               // offset header and other members of VertexCL
                    EL_LDATA,  &v->Unknowns, sizeof(v->Unknowns),
                    EL_END,    v+1);
//        EL_DATAPTR, &v->_BndVerts, sizeof(v->_BndVerts),
// WARNUNG: obige Zeile fuehrt dazu, dass der Zeiger _BndVerts ueberschrieben wird; das fuehrt zu Fehlern in HandlerVertexScatter, falls das Zielobjekt schon
//          auf dem Prozess existiert!
}


void EdgeCL::Define()
{
    EdgeCL* e =0;
    DynamicDataInterfaceCL::TypeDefineEdge( _dddT, e,
                    EL_DDDHDR, &e->_dddH,
                    EL_OBJPTR, &e->_Vertices, sizeof(e->_Vertices), VertexCL::GetType(),
                    EL_OBJPTR, &e->_MidVertex, sizeof(e->_MidVertex), VertexCL::GetType(),
                    EL_GDATA,  &e->_Bnd, sizeof(e->_Bnd),
                    EL_LDATA,  &e->_MFR, sizeof(e->_MFR),               // Wird mit Interface ... akkumuliert!
                    EL_GDATA,  &e->_AccMFR, sizeof(e->_AccMFR),
                    EL_LDATA,  &e->_RemoveMark, sizeof(e->_RemoveMark),
                    EL_END,    e+1);
    // neglect other members for now:
    // Unknowns
    // eventuell noch:        EL_GDATA,  &e->_AccMFR, sizeof(e->_AccMFR),
}


void FaceCL::Define()
{
    FaceCL* f=0;
    DynamicDataInterfaceCL::TypeDefineFace( _dddT, f,
                    EL_DDDHDR, &f->_dddH,
                    EL_GDATA,  &f->_Bnd, sizeof(f->_Bnd),
                    EL_LDATA,  &f->_RemoveMark, sizeof(f->_RemoveMark),
                    EL_END,    f+1);
    // eventuell noch:        EL_OBJPTR, &f->_Neighbors, sizeof(f->_Neighbors), TetraCL::GetType(),
    // neglect other members for now:
    // Unknowns, _Neighbors
    // _Neighbors sind i.A. auf jedem Proc in anderer Reihenfolge gespeichert!
    //    -> Tetras tragen sich nach Xfer selbst ein.
}


void TetraCL::Define()
{
    TetraCL* t=0;
    DynamicDataInterfaceCL::TypeDefineTetra( _dddT, t,
                    EL_DDDHDR, &t->_dddH,
                    EL_GDATA,  &t->_RefRule, sizeof(t->_RefRule),
                    EL_GDATA,  &t->_RefMark, sizeof(t->_RefMark),
                    EL_OBJPTR, &t->_Vertices, sizeof(t->_Vertices), VertexCL::GetType(),
                    EL_OBJPTR, &t->_Edges, sizeof(t->_Edges), EdgeCL::GetType(),
                    EL_OBJPTR, &t->_Faces, sizeof(t->_Faces), FaceCL::GetType(),
                    EL_OBJPTR, &t->_Parent, sizeof(t->_Parent), TetraCL::GetType(),
                    EL_LDATA,  &t->Unknowns, sizeof(t->Unknowns),
                    EL_END,    t+1);
    // neglect other members for now:
    // _Children, Unknowns
}

void ParMultiGridCL::DefineBndPtT()
{
    BndPointCL* b=0;
    DynamicDataInterfaceCL::TypeDefineBndPtT( _BndPtT, b,
                    EL_GDATA,  &b->_BndIdx,  sizeof(b->_BndIdx),        // The Index and
                    EL_GDATA,  &b->_Coord2D, sizeof(b->_Coord2D),       // Coordinate of the Boundary-Point shopuld be submittet
                    EL_END,    b+1);
}

void ParMultiGridCL::DefineChildPtrT()
{
    TetraCL** chp=0;
    DynamicDataInterfaceCL::TypeDefineChildPtrT( _ChildPtrT, chp,
                    EL_OBJPTR, &*chp, sizeof(*chp), TetraCL::GetType(),
                    EL_END,    chp+1);
    //delete chp;
}

void AddedScalCL::Define()
{
    AddedScalCL *add=0;
    DynamicDataInterfaceCL::TypeDefineAddedScal( _dddT, add,
                    EL_GDATA, &add->idxVecDesc_, sizeof(add->idxVecDesc_),
                    EL_GDATA, &add->data_, sizeof(add->data_),
                    EL_END,   add+1);
}

void AddedVecCL::Define()
{
    AddedVecCL *add=0;
    DynamicDataInterfaceCL::TypeDefineAddedVec( _dddT, add,
                    EL_GDATA, &add->idxVecDesc_, sizeof(add->idxVecDesc_),
                    EL_GDATA, &add->data_, sizeof(add->data_),
                    EL_END,   add+1);
}
*/

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

bool CheckParMultiGrid(const ParMultiGridCL& pmg)
/**  Check parallel and sequential multigrid on each processor
    \param pmg The parallel multigrid
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
    bool pmg_sane = pmg.IsSane(output),
         mg_sane  = pmg.GetMG().IsSane(output);
    bool sane     = ProcCL::Check(pmg_sane && mg_sane);
    if (!sane)
        throw DROPSErrCL("CheckParMultiGrid: Multigrid is not sane!");
    return sane;
}

} // end of namespace DROPS
