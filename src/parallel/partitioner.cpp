/// \file partitionerclass.cpp
/// \brief Implementation of the methods of the classes inside partitioner.h file
/// \author LNM RWTH Aachen: ; SC RWTH Aachen: Oliver Fortmeier, Alin Bastea
/// Begin: 18.01.2010

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

#include "parallel/partitioner.h"

namespace DROPS
{

//Implementation of the methods of the GraphST structure
//========================================================
/** Clear memory reserved for all arrays*/
GraphST::~GraphST()
{
    Clear();
}

/// \brief Method used to allocate memory for the vtxdist vector in the struct
void GraphST::ResizeVtxDist()
{
    if (vtxdist!=0) delete[] vtxdist;
    vtxdist = new int[ProcCL::Size()+1];
}

/**This method will allocate memory for the arrays:
   xadj, adjncy, vwgt, adjwgt, and part
   Additionally, the values for number of adjacenies and vertices as well as the parameters geom and ncon are set.
   \param numadj  number of neighbors
   \param myVerts number of verticies
   \param geom    flag for geometry
   \param ncon    number of weight conditions per vertex
 */
void GraphST::Resize(int numadj, int numverts, bool hasgeom, int numncon)
{
    myVerts= numverts;
    myAdjs= numadj;
    geom= hasgeom;
    ncon= numncon;
    if (xadj!=0) delete[] xadj;
    xadj = new idxtype[myVerts+1];
    if (adjncy!=0) delete[] adjncy;
    adjncy = new idxtype[myAdjs];
    if (vwgt!=0) delete[] vwgt;
    vwgt = new idxtype[myVerts*ncon];
    if (adjwgt!=0) delete[] adjwgt;
    adjwgt = new idxtype[myAdjs];
    if (geom)
    {
        if (xyz!=0) delete[] xyz;
        xyz = new float[3*myVerts];
    }
    if (part!=0) delete[] part;
    part = new idxtype[myVerts];
}

/// \brief This method will liberate the memory used for storing the arrays in the structure
void GraphST::Clear()
{
    /** Free all allocated arrays */
    if (xadj!=0)    delete[] xadj;
    if (adjncy!=0)  delete[] adjncy;
    if (vtxdist!=0) delete[] vtxdist;
    if (part!=0)    delete[] part;
    if (vwgt!=0)    delete[] vwgt;
    if (adjwgt!=0)  delete[] adjwgt;
    if (xyz!=0)     delete[] xyz;
    xadj=0; adjncy=0; vtxdist=0; part=0;
    vwgt=0; adjwgt=0; xyz=0;
}
//End of the implementations of the methods in the structure GraphST
//============================================================================================

//Implementation of the methods of the parent class partitioner
//===========================================================================================

/// \brief Constructor of the base partitioner class (allocates memory for the structure attribute)
PartitionerCL::PartitionerCL( PartMethod meth)
    : time_(0), meth_(meth)
{
    allArrays_ = new GraphST();
}

/// \brief Destructor of the base partitioner class (empties memory)
PartitionerCL::~PartitionerCL()
{
    delete allArrays_;
}

/// \brief Getter/setter for the structure member inside the base partitioner class
GraphST& PartitionerCL::GetGraph()
{
    return *allArrays_;
}

/**Factory method that based on the partOption parameter allocates memory for the chosen type of object
  (it implements the Factory design pattern)
 \param partitioner type of partitioner used( 0 - Metis, 1 - Zoltan, 2 - Scotch)
        quality     quality of the partitioning
        method      partition method invoked by the partitioner
 */
PartitionerCL* PartitionerCL::newPartitioner(Partitioner partitioner, float quality, PartMethod method)
{
    PartitionerCL* pa = NULL;
     switch (partitioner)
     {
         case metis:
         {
             pa = new ParMetisCL( method, quality);
             break;
         }
         case zoltan:
         {
#ifdef _ZOLTAN
             pa = new ZoltanCL( method);
#else
             throw DROPSErrCL("PartitionerCL::newPartitioner: Zoltan library is not included");
#endif
             break;
         }
         case scotch:
         {
#ifdef _SCOTCH
             pa = new ScotchCL( method);
#else
             throw DROPSErrCL("PartitionerCL::newPartitioner: Scotch library is not included");
#endif
             break;
         }
     }
     return pa;
}

/** Use own processor rank as destination*/
void PartitionerCL::ApplyIdentity()
{
    std::fill(GetGraph().part, GetGraph().part+GetGraph().myVerts, ProcCL::MyRank());
}


//Implementation of the methods in the derived class ParMetisCL
//=================================================================================================

/** Constructor of the derived class ParMetisCL (parent PartitionerCL)
 \param ubvec       type of partitioner used
        quality     quality of the partitioning
        meth        quality parameter
 */
ParMetisCL::ParMetisCL(PartMethod meth, float ubvec, float quality) : PartitionerCL(meth)
{
    GetGraph().ubvec = ubvec;
    quality_ = quality;
}

/** Serial partitioning
 \param master      rank of the master thread
 */
void ParMetisCL::PartGraphSer(int master)
{
    if ( ProcCL::MyRank()!=master)
        return;

    Comment("- Start calculate LoadBalance with METIS-"<<(meth_ == KWay ? "Kway" : "Recursive")<<std::endl, DebugLoadBalC);

    if (GetGraph().myVerts != GetGraph().vtxdist[ProcCL::Size()])
        Comment("LoadBalCL: PartGraphSer: This procedure is called by a proc, that does not have all nodes!"<<std::endl, DebugLoadBalC);

    int    wgtflag    = 3,                                      ///< Weights on vertices and adjacencies are given
           numflag    = 0,                                      ///< Numbering of verts starts by 0 (C-Style)
           nparts     = ProcCL::Size(),                         ///< Number of subdomains (per proc one)
           n          = GetGraph().myVerts,                     ///< Number of vertices
           options[5] = {0,0,0,0,0};                            ///< Default options
    TimerCL timer; timer.Reset();
    // Depending on the method chosen KWay or Recursive the appropriate metis function is called
    // if number of balance conditions is not equal to 1 use the multi constrained versions
    if ( GetGraph().ncon==1){
        if (meth_ == KWay || meth_== Identity)
            METIS_PartGraphKway(      &n, GetGraph().xadj, GetGraph().adjncy, GetGraph().vwgt, GetGraph().adjwgt, &wgtflag, &numflag,  &nparts, options,&GetGraph().edgecut, GetGraph().part);
        else if (meth_==Recursive)
            METIS_PartGraphRecursive( &n, GetGraph().xadj, GetGraph().adjncy, GetGraph().vwgt, GetGraph().adjwgt, &wgtflag, &numflag,  &nparts, options,&GetGraph().edgecut, GetGraph().part);
    }
    else{
        if (GetGraph().ncon>15)
            throw DROPSErrCL("ParMetisCL::PartGraphSer: Too many constrains are given");
        wgtflag = 1;    // weights on vertices have to be provided, so this is the flag for specifying additional weights on edges
        if (meth_ == KWay || meth_== Identity){
            std::valarray<float> ubvec( GetGraph().ubvec, GetGraph().ncon);
            METIS_mCPartGraphKway(      &n, &GetGraph().ncon, GetGraph().xadj, GetGraph().adjncy, GetGraph().vwgt, GetGraph().adjwgt, &wgtflag, &numflag,  &nparts, Addr(ubvec), options, &GetGraph().edgecut, GetGraph().part);
        }
        else if (meth_==Recursive){
            METIS_mCPartGraphRecursive( &n, &GetGraph().ncon, GetGraph().xadj, GetGraph().adjncy, GetGraph().vwgt, GetGraph().adjwgt, &wgtflag, &numflag,  &nparts, options, &GetGraph().edgecut,GetGraph().part);
        }
    }
    timer.Stop(); time_= timer.GetTime();
    Comment("  * Number of Edgecut: "<<GetGraph().edgecut<<std::endl, DebugLoadBalC);
}

/// \brief Parallel partitioning
void ParMetisCL::PartGraphPar()
/** This procedure uses ParMetis in order to compute a partitioning of the dual
    reduced graph, that has been set up with the member function
    CreateDualRedGraph. The previous distribution of the vertices among the
    processors is used.
    \pre Setup of the graph
*/
{
    Comment("- Start calculate LoadBalance with ParMETIS-"<<(meth_ == Adaptive ? "AdaptiveRepart" : (meth_ == KWay ? "KWay" : "Identity"))<<std::endl,
            DebugLoadBalC);

    Assert( GetGraph().xadj && GetGraph().adjncy && GetGraph().vwgt && GetGraph().adjwgt,
            DROPSErrCL("ParMetisCL::PartGraphPar: Graph has not been set up. Maybe use CreateDualRedGraph before calling this routine"),
            DebugLoadBalC);

    int    wgtflag = 3,                 // Weights on vertices and adjacencies are given
           numflag = 0,                 // numbering of verts starts by 0 (C-Style)
           nparts  = ProcCL::Size(),    // number of subdomains (per proc one)
           options[5] = {0,0,0,0,0};    // default options and no debug information
    float  itr     = quality_;          // how much an exchange costs
    std::valarray<float> ubvec( GetGraph().ubvec, GetGraph().ncon); // allowed inbalance
    std::valarray<float> tpwgts( 1.f/(float)nparts, GetGraph().ncon*ProcCL::Size());

    if (GetGraph().part == 0)
        GetGraph().part = new idxtype[GetGraph().myVerts];

    MPI_Comm comm = MPI_COMM_WORLD;
    //Depending on the method choosen Adaptive, KWay or Identity the apropiate parmetis function is called
    ParTimerCL timer; timer.Reset();
    switch (meth_)
    {
        case Adaptive:
            ParMETIS_V3_AdaptiveRepart(
                    GetGraph().vtxdist, GetGraph().xadj, GetGraph().adjncy, GetGraph().vwgt, GetGraph().vwgt,
                    GetGraph().adjwgt, &wgtflag, &numflag, &GetGraph().ncon, &nparts, Addr(tpwgts),
                    Addr(ubvec), &itr, options, &GetGraph().edgecut, GetGraph().part, &comm);
            break;
        case KWay:
            ParMETIS_V3_PartKway(
                    GetGraph().vtxdist, GetGraph().xadj, GetGraph().adjncy, GetGraph().vwgt, GetGraph().adjwgt,
                    &wgtflag, &numflag, &GetGraph().ncon, &nparts, Addr(tpwgts), Addr(ubvec), options,
                    &GetGraph().edgecut, GetGraph().part, &comm);
            break;
        case Identity:
            ApplyIdentity();
            break;
        case NoMig:
            break;
        case Recursive: // Fall through
        default:
            throw DROPSErrCL("ParMetisCL::PartGraphPar: Unknown partitioning method");
    }
    timer.Stop(); time_= timer.GetTime();
}


#ifdef _ZOLTAN
/// \brief Constructor of the derived class ZoltanCL (parent class PartitionerCL)
ZoltanCL::ZoltanCL( PartMethod meth)
    : PartitionerCL( meth)
{
    zz = Zoltan_Create( ProcCL::GetComm());
}

/// \brief DEstructor of the derived class ZoltanCL (parent class PartitionerCL)
ZoltanCL::~ZoltanCL()
{
    Zoltan_Destroy(&zz);
}

/** Querry function required by Zoltan
 \param data      pointer to the GraphST attribute
        ierr      zoltan error type
 */
int ZoltanCL::get_number_of_vertices(void *data, int *ierr)
{
    GraphST *graph = (GraphST *)data;
    *ierr = ZOLTAN_OK;
     return graph->myVerts;
}

/** Querry function required by Zoltan
 \param data        pointer to the GraphST attribute
        globalID    global id of the actual vertex
        localID     local id of the actual vertex
        obj_wgts    weights of the vertices
        ierr        zoltan error type
 */
void ZoltanCL::get_vertex_list(void *data, int, int, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int, float *obj_wgts, int *ierr)
{
    int i;
    GraphST *graph = (GraphST *)data;
    *ierr = ZOLTAN_OK;
    // Return the IDs of our verticies, and the weights.
    for (i=0; i<graph->myVerts; i++)
    {
      globalID[i] = i+ graph->myfirstVert;
      localID[i] = i;
      obj_wgts[i] = graph->vwgt[i];
    }
}

/** Querry function required by Zoltan
 \param data        pointer to the GraphST attribute
        sizeGID     dimension size of the global id of the actual vertex
        sizeLID     dimension size of the local id of the actual vertex
        num_obj     number of vertices
        localID     local id of the actual vertex
        num_edges   number of edges in the graph
        nbor_GID    global id of the neighboring vertices
        nbor_Proc   processor ot which the neighbor is assigned
        wgt_dim     dimension size of the weights of the edges
        ewgts       weights of the edges
        ierr        Zoltan error type
 */
void ZoltanCL::get_edge_list(void *data, int sizeGID, int sizeLID, int num_obj, ZOLTAN_ID_PTR , ZOLTAN_ID_PTR localID, int *num_edges, ZOLTAN_ID_PTR nborGID, int *nborProc, int wgt_dim, float *ewgts, int *ierr)//returns the list of edges of each processor
{
    int i, j, from, to;
    int *nextNbor, *nextProc;
    float *nextWgt;

    GraphST *graph = (GraphST *)data;
    *ierr = ZOLTAN_OK;

    if ( (sizeGID != 1) || (sizeLID != 1) || (num_obj != graph->myVerts)|| (wgt_dim != 1))
    {
        *ierr = ZOLTAN_FATAL;
        return;
    }

      nextNbor = (int *)nborGID;
      nextProc = nborProc;
      nextWgt = ewgts;

      for (i=0; i < num_obj; i++)
      {
        to = graph->xadj[localID[i]+1];
        from = graph->xadj[localID[i]];
        if ((to - from) != num_edges[i])
        {
            *ierr = ZOLTAN_FATAL;
            return;
        }
        for (j=from; j < to; j++)
        {
          *nextNbor++ = graph->adjncy[j];//?????
          *nextProc++ = graph->GetProc(graph->adjncy[j]);
          *nextWgt++ = graph->adjwgt[j];
        }
      }
}

/** Querry function required by Zoltan
 \param data        pointer to the GraphST attribute
        sizeGID     dimension size of the global id of the actual vertex
        sizeLID     dimension size of the local id of the actual vertex
        num_obj     number of vertices
        num_Edges   number of edges
        ierr        Zoltan error type
 */
void ZoltanCL::get_num_edges_list(void *data, int sizeGID, int sizeLID, int num_obj, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int *numEdges, int *ierr)// returns the total number of edges
{
    int i;
      GraphST *graph = (GraphST *)data;
      if ( (sizeGID != 1) || (sizeLID != 1) || (num_obj != graph->myVerts))
      {
        *ierr = ZOLTAN_FATAL;
        return;
      }
      for (i=0;  i < num_obj ; i++)
        numEdges[i] = graph->xadj[i+1] - graph->xadj[i];
      *ierr = ZOLTAN_OK;
}

/// \brief Method that sets some parameters inside Zoltan so that the communication between drops and the partitioner is correct
void ZoltanCL::CreateGraph()
{
    if ( GetGraph().ncon>1){
        throw DROPSErrCL("ZoltanCL::CreateGraph: Cannot handle multiple vertex weights by a Zoltan partitioner");
    }

    Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM","1");
    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
    Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1");
    Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");

    Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, &GetGraph());
    Zoltan_Set_Obj_List_Fn(zz, get_vertex_list,  &GetGraph());
    Zoltan_Set_Num_Edges_Multi_Fn(zz, get_num_edges_list,  &GetGraph());
    Zoltan_Set_Edge_List_Multi_Fn(zz, get_edge_list,  &GetGraph());
}

/// \brief Serial partitioning
void ZoltanCL::PartGraphSer( int)
{
    Comment("- Start calculate LoadBalanace with Zoltan-Partition"<<std::endl, DebugLoadBalC);

    if (GetGraph().myVerts != GetGraph().vtxdist[DynamicDataInterfaceCL::InfoProcs()])
        Comment("LoadBalCL: PartGraphSer: This procedure is called by a proc, that does not have all nodes!"<<std::endl, DebugLoadBalC);
    int i;
    int changes ; // Set to 1 or .TRUE. if the decomposition was changed by the load-balancing method; 0 or .FALSE. otherwise.
    int numGidEntries; //the number of array entries used to describe a single global ID
    int numLidEntries; //the number of array entries used to describe a single local ID
    int numImport;
    ZOLTAN_ID_PTR importGlobalGids;
    ZOLTAN_ID_PTR importLocalGids;
    int *importProcs;
    int *importToPart;
    int numExport;
    ZOLTAN_ID_PTR exportGlobalGids;
    ZOLTAN_ID_PTR exportLocalGids;
    int *exportProcs;
    int *exportToPart;
    int rc;

    Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");//Partition "from scratch," not taking into account the current data distribution;
                                                //this option is recommended for static load balancing
    ParTimerCL timer; timer.Reset();
    rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
          &changes,     /* 1 if partitioning was changed, 0 otherwise */
          &numGidEntries,  /* Number of integers used for a global ID */
          &numLidEntries,  /* Number of integers used for a local ID */
          &numImport,      /* Number of vertices to be sent to me */
          &importGlobalGids,  /* Global IDs of vertices to be sent to me */
          &importLocalGids,   /* Local IDs of vertices to be sent to me */
          &importProcs,    /* Process rank for source of each incoming vertex */
          &importToPart,   /* New partition for each incoming vertex */
          &numExport,      /* Number of vertices I must send to other processes*/
          &exportGlobalGids,  /* Global IDs of the vertices I must send */
          &exportLocalGids,   /* Local IDs of the vertices I must send */
          &exportProcs,    /* Process to which I send each of the vertices */
          &exportToPart);
    timer.Stop(); time_= timer.GetTime();

    if (rc != ZOLTAN_OK)
        throw DROPSErrCL("ZoltanCL::PartGraphSer: Zoltan returned an error!");

    for (i=0; i < GetGraph().myVerts; i++)
        GetGraph().part[i] = ProcCL::MyRank();

    for (i=0; i < numExport; i++)
        GetGraph().part[exportLocalGids[i]] = exportToPart[i];
}

/// \brief Parallel partitioning
void ZoltanCL::PartGraphPar()
{
    Comment("- Start calculate LoadBalanace with Zoltan-Repartition"<<std::endl,DebugLoadBalC);

    Assert( GetGraph().xadj && GetGraph().adjncy && GetGraph().vwgt && GetGraph().adjwgt,
            DROPSErrCL("LoadBalCL::PartGraphPar: Graph has not been set up. Maybe use CreateDualRedGraph before calling this routine"),
            DebugLoadBalC);
    int i;
    int changes; // Set to 1 or .TRUE. if the decomposition was changed by the load-balancing method; 0 or .FALSE. otherwise.
    int numGidEntries; //the number of array entries used to describe a single global ID
    int numLidEntries; //the number of array entries used to describe a single local ID
    int numImport;
    ZOLTAN_ID_PTR importGlobalGids;
    ZOLTAN_ID_PTR importLocalGids;
    int *importProcs;
    int *importToPart;
    int numExport;
    ZOLTAN_ID_PTR exportGlobalGids;
    ZOLTAN_ID_PTR exportLocalGids;
    int *exportProcs;
    int *exportToPart;
    int rc;

    if ( meth_==Adaptive){
        Zoltan_Set_Param(zz, "LB_APPROACH", "REPARTITION"); //Partition but take into account current data distribution to keep data migration low;
                                                            //this option is recommended for dynamic load balancing
    }
    else {
        Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
    }

    ParTimerCL timer; timer.Reset();
    if ( meth_==Identity){
        ApplyIdentity();
    }
    else if ( meth_!=NoMig){
        rc = Zoltan_LB_Partition(zz,/* input (all remaining fields are output) */
                &changes,           /* 1 if partitioning was changed, 0 otherwise */
                &numGidEntries,     /* Number of integers used for a global ID */
                &numLidEntries,     /* Number of integers used for a local ID */
                &numImport,         /* Number of vertices to be sent to me */
                &importGlobalGids,  /* Global IDs of vertices to be sent to me */
                &importLocalGids,   /* Local IDs of vertices to be sent to me */
                &importProcs,       /* Process rank for source of each incoming vertex */
                &importToPart,      /* New partition for each incoming vertex */
                &numExport,         /* Number of vertices I must send to other processes*/
                &exportGlobalGids,  /* Global IDs of the vertices I must send */
                &exportLocalGids,   /* Local IDs of the vertices I must send */
                &exportProcs,       /* Process to which I send each of the vertices */
                &exportToPart);
    }
    timer.Stop(); time_= timer.GetTime();

    if (rc != ZOLTAN_OK)
        throw DROPSErrCL("ZoltanCL::PartGraphPar: Zoltan returned an error!");

    for (i=0; i < GetGraph().myVerts; i++)
        GetGraph().part[i] = ProcCL::MyRank();

    for (i=0; i < numExport; i++)
        GetGraph().part[exportLocalGids[i]] = exportToPart[i];
}
# endif

#ifdef _SCOTCH
ScotchCL::ScotchCL(PartMethod meth)
    : PartitionerCL(meth)
{}

void ScotchCL::CreateGraph()
{
    if ( GetGraph().ncon>1){
        throw DROPSErrCL("ZoltanCL::CreateGraph: Cannot handle multiple vertex weights by a Scotch partitioner");
    }
}

void ScotchCL::PartGraphPar()
{
    if ( meth_==Identity){
        ParTimerCL timer; timer.Reset();
        ApplyIdentity();
        timer.Stop(); time_= timer.GetTime();
        return;
    }

    SCOTCH_Dgraph   graph;
    SCOTCH_Arch     archdat;
    SCOTCH_Strat    strategy;
    SCOTCH_Dmapping mappdat;
    std::valarray<SCOTCH_Num> velotab( 1, ProcCL::Size());

    SCOTCH_dgraphInit( &graph, ProcCL::GetComm());
    if ( SCOTCH_dgraphBuild( &graph, 0, GetGraph().myVerts, GetGraph().myVerts, GetGraph().xadj, GetGraph().xadj+1, GetGraph().vwgt, 0, GetGraph().myAdjs, GetGraph().myAdjs, GetGraph().adjncy, 0, GetGraph().adjwgt)!=0)
        throw DROPSErrCL("ScotchCL::PartGraphPar: Cannot build SCOTCH graph");

    Assert( SCOTCH_dgraphCheck( &graph)==0, DROPSErrCL("ScotchCL::PartGraphPar: Graph is not consistent"), DebugLoadBalC);

    SCOTCH_stratInit( &strategy);

    if ( SCOTCH_archCmpltw (&archdat, ProcCL::Size(), Addr(velotab))!=0)
        throw DROPSErrCL("ScotchCL::PartGraphPar: Cannot compute ltw");
    if ( SCOTCH_dgraphMapInit (&graph, &mappdat, &archdat, GetGraph().part)!=0)
        throw DROPSErrCL("ScotchCL::PartGraphPar: Cannot assign architecture");

    ParTimerCL timer; timer.Reset();
    SCOTCH_dgraphMapCompute ( &graph, &mappdat, &strategy);
    timer.Stop(); time_= timer.GetTime();

    SCOTCH_dgraphMapExit ( &graph, &mappdat);
    SCOTCH_archExit (&archdat);
    SCOTCH_stratExit( &strategy);
    SCOTCH_dgraphExit(  &graph);
}

void ScotchCL::PartGraphSer( int master)
{
    if ( ProcCL::MyRank()!=master)
        return;
    SCOTCH_Graph graph;
    SCOTCH_graphInit( &graph);
    SCOTCH_graphBuild( &graph, 0, GetGraph().myVerts, GetGraph().xadj, GetGraph().xadj+1, GetGraph().vwgt, 0, GetGraph().myAdjs, GetGraph().adjncy, GetGraph().adjwgt);
    Assert( SCOTCH_graphCheck( &graph)==0, DROPSErrCL("ScotchCL::PartGraphSer: Graph is not consistent"), DebugLoadBalC);
    SCOTCH_Strat strategy;
    SCOTCH_stratInit( &strategy);
    TimerCL timer; timer.Reset();
    SCOTCH_graphPart( &graph, ProcCL::Size(), &strategy, GetGraph().part);
    timer.Stop(); time_= timer.GetTime();
    SCOTCH_stratExit( &strategy);
    SCOTCH_graphExit( &graph);
}
#endif

}//end of namespace
