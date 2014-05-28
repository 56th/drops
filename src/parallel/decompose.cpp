/// \file decompose.cpp
/// \brief determining a partitioning of a distributed hierarchy of triangulations.
/// \author LNM RWTH Aachen: Patrick Esser, Sven Gross; SC RWTH Aachen: Oliver Fortmeier

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
*/

#include "parallel/decompose.h"
#include "parallel/parallel.h"
#include "num/interfacePatch.h"

namespace DROPS {

// L B  I T E R A T O R  C L
//--------------------------

/** Take tetra out of the coarsest level of the multigrid an check if it is in
    the load balancing set
    \param mg the MultiGrid
    \param TriLevel if greater than zero the level to be partitioned, else
     the finest level
*/
LbIteratorCL LbIteratorCL::makeBegin( MultiGridCL& mg, int TriLevel)
{
    TriLevel = TriLevel<0 ? mg.GetLastLevel() : TriLevel;
    Uint non_empty_level=0;
    while (mg.GetTetras().IsLevelEmpty(non_empty_level) && non_empty_level<mg.GetTetras().GetNumLevel()-1)
        ++non_empty_level;
    LbIteratorCL ret( &mg, mg.GetTetrasBegin(non_empty_level), non_empty_level, TriLevel);
    return ( !mg.GetTetras().IsLevelEmpty(non_empty_level) && ret.IsInLbSet()) ? ret : ++ret;
}

/** iterator to the last tetra in LoadBalSet
    \param mg the MultiGrid
    \param TriLevel if greater than zero the level to be partitioned, else
     the finest level
*/
LbIteratorCL LbIteratorCL::makeEnd( MultiGridCL& mg, int TriLevel)
{
    TriLevel = TriLevel<0 ? mg.GetLastLevel() : TriLevel;
    return LbIteratorCL( &mg, mg.GetTetrasEnd(TriLevel), TriLevel, TriLevel);
}


// G R A P H  C L
//---------------

/** Only allocate memory for vtxdist.*/
GraphCL::GraphCL()
    : xadj_(0), adjncy_(0), vtxdist_( ProcCL::Size()+1), vwgt_(0), adjwgt_(0)
{}

/** Determine the first vertex each process is responsible for. */
void GraphCL::buildVtxDist( graph_index_type num_vert)
{
    num_procs_w_verts_= 1;
    if (num_vert==0) {
        num_procs_w_verts_= 0;
        // ParMetis does not allow procs without any graph vertex.
        // As a work-around, in this case an artificial graph vertex with zero weight is added.
        num_vert= 1;
    }
    num_procs_w_verts_= ProcCL::GlobalSum( num_procs_w_verts_);

    IndexArray verts= ProcCL::Gather(num_vert, -1);
    vtxdist_[0]=0;
    for ( int p=0; p<ProcCL::Size(); ++p)
        vtxdist_[p+1]= vtxdist_[p]+verts[p];
}


/** Allocate memory to represent the graph given by the number of
    adjacencies and vertices.
    \param num_adj number of adjnacencies
    \param num_Vert number of vertices
*/
void GraphCL::resize( int num_adj, int num_vert)
{
    if (num_adj==0) {
        // ParMetis does not allow procs without any graph edge.
        // As a work-around, in this case an artificial graph edge with zero weight is added, which connects graph vertex 0 with ithself.
        num_adj= 1;
    }
    if (num_vert==0)
        num_vert= 1;
    xadj_.resize( num_vert+1);
    adjncy_.resize( num_adj);
    vwgt_.resize( num_vert);
    adjwgt_.resize( num_adj);
    part_.resize( num_vert);
}


// V E R T E X  W E I G H T E R S
//-------------------------------

/// \brief Class for weighting vertices. For more details about the weighting function, see Diss. Fortmeier
/** This is only the base class, that assigns each vertex the weight 1. */
class VertexWeighterCL
{
public:
    /// \brief constructor
    VertexWeighterCL() {}
    /// \brief destructor
    virtual ~VertexWeighterCL() {}
    /// \brief get the weight of the vertex (which is given by the tetrahedron refered by an LBIteratorCL)
    virtual int getWeight( const TetraCL&) const { return 1; }
};

/// \brief Weight each vertex by the number of represented tetrahedra located on a
///        specified triangulation level
class VertexWeighterTriangCL : public VertexWeighterCL
{
protected:
    Uint triang_level_;         ///< the triangulation level

public:
    /// \brief constructor
    VertexWeighterTriangCL( const Uint triang_level)
        : VertexWeighterCL(), triang_level_( triang_level) {}
    /// \brief destructor
    virtual ~VertexWeighterTriangCL() {}
    /// \brief get the weight of the vertex
    virtual int getWeight( const TetraCL& t) const
    {
        graph_index_type vert_weight= 1;
        if ( !t.IsUnrefined()){
            TetraCL::const_ChildPIterator ch= t.GetChildBegin();
            for ( ; ch!=t.GetChildEnd(); ++ch){
                if ( (*ch)->IsInTriang( triang_level_)){
                    ++vert_weight;
                }
            }
        }
        return vert_weight;
    }
};

/// \brief Weight tetrahedra intersected by the interface by a factor rho_I more than
///        not intersected tetrahedra
/** \todo Check me, if DOF handling works with DiST! */
class VertexWeighterTwoPhaseCL : public VertexWeighterTriangCL
{
protected:
    int                factor_;     ///< factor: how much is an intersected tetrahedron expensive
    const VecDescCL&   lset_;       ///< level set function
    const BndDataCL<>& lsetbnd_;    ///< level set boundary conditions
    /// \brief check if on each vertex and edge of a given tetrahedron a level set unknown is given
    inline bool CheckForLsetUnk( const TetraCL& t) const;

public:
    /// \brief constructor that needs access to the level set fucntion
    VertexWeighterTwoPhaseCL( const Uint triang_level, const int factor, const VecDescCL& lset, const BndDataCL<>& lsetbnd)
        : VertexWeighterTriangCL( triang_level), factor_(factor),
          lset_(lset), lsetbnd_(lsetbnd) {}
    /// \brief destructor
    ~VertexWeighterTwoPhaseCL() {}
    /// \brief get the weight
    inline int getWeight( const TetraCL& t) const;
};

/** Check if level set unknowns exists on all vertices and edges of a tetrahedron */
bool VertexWeighterTwoPhaseCL::CheckForLsetUnk( const TetraCL& t) const
{
    const Uint idx= lset_.RowIdx->GetIdx();
    // check vertices
    for ( TetraCL::const_VertexPIterator sit=t.GetVertBegin(); sit!=t.GetVertEnd(); ++sit){
        if ( !(*sit)->Unknowns.Exist(idx))
            return false;
    }
    // check edges
    for ( TetraCL::const_EdgePIterator sit= t.GetEdgesBegin(); sit!=t.GetEdgesEnd(); ++sit){
        if ( !(*sit)->Unknowns.Exist(idx))
            return false;
    }
    return true;
}

/** Determine the weight of a vertex by rho_I*num_intersected + num_nut_intersected children.
    \todo CHECK ME! (I cannot check this model, since no unknowns are working right now in DiST!
*/
int VertexWeighterTwoPhaseCL::getWeight( const TetraCL& t) const
{
    InterfacePatchCL patch;                 // used to check for an intersection
    int intersected=0, notintersected=0;    // number of intersected and not-intersected tetras

    if ( t.IsUnrefined()){                      // unrefined tetra
        if ( CheckForLsetUnk(t)){               // only apply patch if level set unknowns are available
            patch.Init( t, lset_, lsetbnd_);
            if ( patch.Intersects())
                intersected= 1;                 // increase number of intersected tetras
            else
                notintersected= 1;              // increase number of not-intersected tetras
        }
        else{
            notintersected= 1;                  // tetras without level set unknowns count for not-intersected tetras
        }
    }
    else{                                       // refined tetra (docu: see above ;-) )
        TetraCL::const_ChildPIterator ch=t.GetChildBegin();
        for ( ; ch!=t.GetChildEnd(); ++ch){
            if ( CheckForLsetUnk(**ch)){
                patch.Init( **ch, lset_, lsetbnd_);
                if ( patch.Intersects())
                    ++intersected;
                else
                    ++notintersected;
            }
        }
    }
    return intersected*factor_+notintersected;  // compute the weight
}

// E D G E   W E I G H T E R S
//----------------------------

/// \brief Base class for weighting edges. This class weights each edge by 1
class EdgeWeighterCL
{
public:
    /// \brief Map each neighbor to a weight
    typedef std::map<graph_index_type, graph_index_type> NeighWeightMapT;
    /// \brief Iterator of the map that maps each neighbor to a weight
    typedef NeighWeightMapT::const_iterator const_neigh_iterator;

protected:
    NeighWeightMapT  adj_weight_;   ///< map neighbor index to the weight 1
    graph_index_type myid_;         ///< id of the calling vertex

public:
    /// \brief constructor
    EdgeWeighterCL() {}
    /// \brief destructor
    virtual ~EdgeWeighterCL() {}
    /// \brief Init the datastructures for a single vertex with a given ID
    virtual void init( const graph_index_type myGlobalId)
        { adj_weight_.clear(); myid_= myGlobalId; }
    /// \brief Finalize the weighting process of a single vertex
    virtual void finalize() {}
    /// \brief Get the weight that corresponds to a given neighbor
    int getWeight( graph_index_type globalNeighId) { return adj_weight_[ globalNeighId]; }
    /// \brief Update the neighbor list and the weight
    virtual void update( graph_index_type globalNeighId, const FaceCL* const)
        { if ( globalNeighId!=myid_) adj_weight_[globalNeighId]= 1; }
    /// \name iterate over all neighbors
    //@{
    const_neigh_iterator begin() const { return adj_weight_.begin();  }
    const_neigh_iterator end() const { return adj_weight_.end();  }
    //@}
};

/// \brief Weight the edges by the number of faces among the two tetrahedron families
///        represented by the two adjacent vertices
class EdgeWeighterTriangCL : public EdgeWeighterCL
{
public:
    /// \brief constructor
    EdgeWeighterTriangCL() : EdgeWeighterCL() {}
    /// \brief destructor
    ~EdgeWeighterTriangCL() {}
    /// \brief remember the neighbor
    /// \pre face must not be refined!
    void update( graph_index_type globalNeighId, __UNUSED__ const FaceCL* const fp)
    {
        Assert( !fp->IsRefined(), DROPSErrCL("EdgeWeighterTriangCL::update: Face must not be refined"), DebugLoadBalC);
        if ( globalNeighId!=myid_) ++adj_weight_[globalNeighId];
    }
};

/// \brief Weight edges by the number of DoF located at the faces among the two tetrahedron families
///        represented by the two adjacent vertices
/** \todo Check me, if DOF handling works with DiST! */
class EdgeWeighterDOFCL : public EdgeWeighterCL
{
protected:
    typedef EdgeWeighterCL base;
    typedef std::vector< std::set<IdxT> > DOFContainer;     ///< type for number of common dof among two tetrahedron families
    std::vector<const IdxDescCL*> idx_;                     ///< index describer of DoF that should be considered
    std::map<graph_index_type, DOFContainer> dofs_;         ///< number of common dof among two tetrahedron families

public:
    /// \brief constructor with a given set of dof that should be considered
    EdgeWeighterDOFCL(const ObservedMigrateFECL& obs)
        : base()
    {
        ObservedMigrateFECL::const_iterator it= obs.begin();
        for ( ; it!=obs.end(); ++it)
            idx_.push_back( it->GetIdxDesc());
    }
    /// \brief destructor
    ~EdgeWeighterDOFCL() {}
    /// \brief init for a given vertex
    void init( const graph_index_type myGlobalId)
        { base::init( myGlobalId); dofs_.clear(); }

    /// \brief Compute the weights
    void finalize();

    /// \brief update the weight for a given face
    void update( graph_index_type globalNeighId, const FaceCL* const fp);
};

/** Call this function before(!) call of getWeight. This function
    assumes the same number on unknowns---if exists--- on edges
    than on vertices.
*/
void EdgeWeighterDOFCL::finalize()
{
    std::map<graph_index_type, graph_index_type>::iterator it= adj_weight_.begin();
    for ( ; it!=adj_weight_.end(); ++it){
        int numUnk=0;
        for ( size_t i=0; i<dofs_[it->first].size(); ++i)
            numUnk +=  dofs_[it->first][i].size() * idx_[i]->NumUnknownsVertex();
        it->second= numUnk;
    }
}

/** Consider the dof information located on the face for weighting.
    \pre face must not be refined
*/
void EdgeWeighterDOFCL::update( graph_index_type globalNeighId, const FaceCL* const fp)
{
    Assert( !fp->IsRefined(), DROPSErrCL("EdgeWeighterDOFCL::update: Face must not be refined"), DebugLoadBalC);
    if ( dofs_.find( globalNeighId)==dofs_.end()){
        dofs_[ globalNeighId]= DOFContainer( idx_.size());
    }
    for ( Uint v=0; v<3; ++v){
        for ( Uint i=0; i<idx_.size(); ++i){
            if ( fp->GetVertex(v)->Unknowns.Exist( idx_[i]->GetIdx())){
                const IdxT DOF= fp->GetVertex(v)->Unknowns(idx_[i]->GetIdx());
                dofs_[ globalNeighId][i].insert( DOF);
                if ( idx_[i]->IsExtended( DOF))
                    dofs_[ globalNeighId][i].insert( idx_[i]->GetXidx()[DOF]);
            }
        }
    }
    for ( Uint e=0; e<3; ++e){
        for ( Uint i=0; i<idx_.size(); ++i){
            if ( fp->GetEdge(e)->Unknowns.Exist( idx_[i]->GetIdx())){
                const IdxT DOF= fp->GetEdge(e)->Unknowns(idx_[i]->GetIdx());
                dofs_[ globalNeighId][i].insert( DOF);
                if ( idx_[i]->IsExtended( DOF))
                    dofs_[ globalNeighId][i].insert( idx_[i]->GetXidx()[DOF]);
            }
        }
    }
}


// G R A P H  B U I L D E R S
//---------------------------

/// \brief Abstract base class for all graph builders
class BaseGraphBuilderCL
{
protected:
    MultiGridCL&     mg_;               ///< the multigrid
    Uint             triang_level_;     ///< finest triangulation level that should be considered
    const IdxDescCL& vertexIdx_;        ///< index describer for numbering the tetrahedra
    graph_index_type num_vert_;         ///< local number of graph vertices
    graph_index_type my_first_vert_;    ///< my first vertex number

    /// \brief get the triangulation level
    Uint getTriangLevel() const { return triang_level_; }
    /// \brief get the index of the graph numbering
    Uint getVertexIdx() const { return vertexIdx_.GetIdx(); }
    /// \brief Locally number the tetrahedra to obtain a numbering of the vertices for the graph
    void buildVertexNumbering();

public:
    /// \brief constructor
    BaseGraphBuilderCL( MultiGridCL& mg, Uint triang_level, const IdxDescCL& vIdx)
        : mg_( mg), triang_level_( triang_level), vertexIdx_( vIdx),
          num_vert_(0),  my_first_vert_(0) {}
    /// \brief destrcutor
    virtual ~BaseGraphBuilderCL() {}
    /// \brief build the graph
    virtual void build( BaseGraphCL*&, VertexWeighterCL&, EdgeWeighterCL&) = 0;
    /// \brief get the global number of a vertex that corresponds to a tetrahedron
    graph_index_type getGlobalVertexId( const TetraCL& t) const
        { return t.Unknowns( getVertexIdx()) + my_first_vert_; }
    /// \brief Get a pointer to the graph
    virtual BaseGraphCL* getGraph() { return 0; }
    /// \brief Get a constant pointer to the graph
    virtual BaseGraphCL const * getGraph() const { return 0; }
};

/** Make use of the index description vertexIdx_ to assign a number to the
    tetrahedra, i.e., vertices of the graph
*/
void BaseGraphBuilderCL::buildVertexNumbering()
{
    num_vert_=0;
    const Uint idx= vertexIdx_.GetIdx();
    LbIteratorCL it = LbIteratorCL::makeBegin( mg_, triang_level_),
                 end= LbIteratorCL::makeEnd( mg_, triang_level_);
    TetraCL::ChildPIterator ch;
    for ( ; it!=end; ++it, ++num_vert_){
        it->Unknowns.Prepare( idx);
        it->Unknowns( idx)= num_vert_;
        for ( ch=it->GetChildBegin(); ch!=it->GetChildEnd(); ++ch){
            if ( (*ch)->IsInTriang( triang_level_)){
                (*ch)->Unknowns.Prepare( idx);
                (*ch)->Unknowns( idx)= num_vert_;
            }
        }
    }
}

/// \brief Class for building a standard graph
class GraphBuilderCL : public BaseGraphBuilderCL
{
private:
    class CommunicateAdjacencyCL;       ///< Interface communication class for communicating vertex numbers across process boundaries

    GraphCL*         graph_;            ///< the graph to be build
    IdxDescCL        neighVertexIdx_;   ///< index describer, to store inter process neighbors on faces
    graph_index_type num_adj_;          ///< local number of adjacencies

    /// \brief Communicate the local vertices among the processes
    void buildVtxDist( const graph_index_type myVerts)
        { getGraph()->buildVtxDist( myVerts); my_first_vert_= getGraph()->vtxdist()[ProcCL::MyRank()]; }
    /// \brief Get the index of storing neighbors on faces
    Uint getNeighVertexIdx() const { return neighVertexIdx_.GetIdx(); }
    /// \brief Determine the number of local adjacencies
    void determineNumberAdjacencies();
    /// \brief For a given tetrahedron, put the neighbors and the weights in the graph data structure
    void buildNeighbors( const TetraCL&, graph_index_type&, EdgeWeighterCL&);
    /// \brief Fill the data into the graph data structure
    void buildVerticesAndEdges( VertexWeighterCL&, EdgeWeighterCL&);

public:
    /// \brief constructor
    GraphBuilderCL( MultiGridCL& mg, Uint triang_level, const IdxDescCL& vIdx)
        : BaseGraphBuilderCL( mg, triang_level, vIdx),
          graph_( 0), neighVertexIdx_( P1D_FE), num_adj_(0) {}

    /// \brief build the graph
    void build( BaseGraphCL*&, VertexWeighterCL&, EdgeWeighterCL&);
    /// \brief get a pointer to the graph
    GraphCL* getGraph() { return graph_; }
    /// \brief get a constant pointer to the graph
    GraphCL const * getGraph() const { return graph_; }
};


/** Transfer of vertex number among process boundaries (via faces) and, also, delete them.*/
class GraphBuilderCL::CommunicateAdjacencyCL
{
private:
    GraphBuilderCL&  gb_;               ///< the corresponding graph builder to get vertex numbers
    const IdxDescCL& neighVertexIdx_;   ///< where to store the neighbor on faces

public:
    /// \brief constructor
    CommunicateAdjacencyCL( GraphBuilderCL& gb, const IdxDescCL& neighVertIdx)
        : gb_(gb), neighVertexIdx_( neighVertIdx) {}
    /// \brief destructor
    ~CommunicateAdjacencyCL() {}

    /// \brief Collect the adjacent vertex number on the sender side
    bool Gather( const DiST::TransferableCL& t, DiST::SendStreamCL& s)
    {
        FaceCL const * fp; simplex_cast( t, fp);   // transform t to a face
        if ( !fp->IsInTriang( gb_.getTriangLevel())){
            return false;
        }
        TetraCL const * tp= fp->GetSomeTetra();
        if ( tp->HasGhost())
            tp= fp->GetNeighborTetra( tp);
        Assert( tp->Unknowns.Exist( gb_.getVertexIdx()), DROPSErrCL("CommunicateAdjacencyCL::Gather: Missing LB number"), DebugLoadBalC);
        s << gb_.getGlobalVertexId( *tp);
        return true;
    }

    /// \brief Assign the adjacent vertex number from the sender to the face
    bool Scatter( DiST::TransferableCL& t, const size_t numData, DiST::MPIistreamCL& r)
    {
        FaceCL * fp; simplex_cast( t, fp);   // transform t to a face
        graph_index_type remote_neigh;
        TetraCL const * tp= fp->GetSomeTetra();
        if ( tp->HasGhost())
            tp= fp->GetNeighborTetra( tp);
        const graph_index_type local_lb= gb_.getGlobalVertexId(*tp);

        for ( size_t i=0; i<numData; ++i){
            r >> remote_neigh;
            if ( local_lb!=remote_neigh){
                fp->Unknowns.Prepare( neighVertexIdx_.GetIdx());
                fp->Unknowns( neighVertexIdx_.GetIdx())= remote_neigh;
            }
        }
        return true;
    }

    /// \brief Perform the interface communication
    void Call()
    {
        DiST::InterfaceCL::DimListT dimlist; dimlist.push_back( 2);
        DiST::PrioListT Prios; Prios.push_back(PrioMaster);
        const DiST::LevelListCL allLevels;
        DiST::InterfaceCL comm( allLevels, Prios, Prios, dimlist);
        comm.Communicate( *this);
    }

    /// \name Delete the assigned LB numbers
    //@{
    void RemoveLbNumbersOnFaces()
    {
        DiST::InterfaceCL::DimListT dimlist; dimlist.push_back( 2);
        DiST::PrioListT Prios; Prios.push_back(PrioMaster);
        const DiST::LevelListCL allLevels;
        DiST::InterfaceCL comm( allLevels, Prios, Prios, dimlist);
        comm.ExecuteLocal( *this);
    }
    bool operator() ( DiST::TransferableCL& t)
    {
        FaceCL* fp; simplex_cast( t, fp);   // transform t to a face
        if ( fp->Unknowns.Exist( neighVertexIdx_.GetIdx()))
            fp->Unknowns.Invalidate( neighVertexIdx_.GetIdx());
        return true;
    }
    //@}
};


void GraphBuilderCL::determineNumberAdjacencies()
{
    num_adj_= 0;                        // number of adjacencies
    std::set<graph_index_type> adj;     // neighbors
    TetraCL::ChildPIterator ch;         // iterator over children

    // iterate over all tetrahedra referenced by the LbIteratorCL
    LbIteratorCL it= LbIteratorCL::makeBegin( mg_, triang_level_);
    for ( ; it!=LbIteratorCL::makeEnd( mg_, triang_level_); ++it){
        if ( it->IsUnrefined()){                            // unrefined tetrahedron
            for ( Uint face=0; face<NumFacesC; ++face){     // for each face, there might be a neighbor
                if ( !it->IsBndSeg(face)){
                    ++num_adj_;
                }
            }
        }
        else{                                               // refined tetrahedron, so check children
            adj.clear();
            for ( ch= it->GetChildBegin(); ch!=it->GetChildEnd(); ++ch){
                if ( (*ch)->IsInTriang( triang_level_)){
                    for ( Uint face=0; face<NumFacesC; ++face){
                        FaceCL const * fp= (*ch)->GetFace(face);
                        if ( fp->IsOnBoundary()){       // no neighbor
                            continue;
                        }
                        else if ( fp->IsOnProcBnd()){   // neighbor on other process, so it is stored on the face
                            Assert( fp->Unknowns.Exist( getNeighVertexIdx()), DROPSErrCL("GraphBuilderCL::determineNumberAdjacencies: Adjacency across process boundary not known"), DebugLoadBalC);
                            adj.insert( fp->Unknowns( getNeighVertexIdx()));
                        }
                        else{                           // local neighbor, so ask for the vertex number
                            TetraCL const * neigh= (*ch)->GetNeighInTriang( face, triang_level_);
                            Assert( neigh, DROPSErrCL("GraphBuilderCL::buildAdjNum: No local neighbor found"), DebugLoadBalC);
                            Assert( neigh->Unknowns.Exist( getVertexIdx()), DROPSErrCL("GraphBuilderCL::determineNumberAdjacencies: Tetrahedron without load balancing number"), DebugLoadBalC);
                            adj.insert( getGlobalVertexId( *neigh));
                        }
                    }
                }
            }
            adj.erase( getGlobalVertexId( *it));    // do not count the edge to myself
            num_adj_+= adj.size();                  // increase the number of neighbors
        }
    }
}

/** The EdgeWeighterCL is used to store all neighbors and the corresponding weight. The
    adjk_counter is always (in and out) the position where to store the next adjanency.
*/
void GraphBuilderCL::buildNeighbors( const TetraCL& t, graph_index_type& adj_count, EdgeWeighterCL& ew)
{
    ew.init( getGlobalVertexId( t));
    // collect neighbors and weights for an unrefined tetrahedron
    if ( t.IsUnrefined()){
        for ( Uint face= 0; face<NumFacesC; ++face){        // neighbors have a common face
            if ( !t.IsBndSeg(face)){
                FaceCL const * fp = t.GetFace(face);
                // get global id of the neighbor vertex
                graph_index_type neigh_id= -1;
                if ( fp->IsOnProcBnd()){                    // neighbor on other process, so the number is stored on the face
                    Assert( fp->Unknowns.Exist( getNeighVertexIdx()), DROPSErrCL("GraphBuilderCL::buildNeighbors: (unref tetra) Adjacency across process boundary not known"), DebugLoadBalC);
                    neigh_id= fp->Unknowns( getNeighVertexIdx());
                }
                else{                                       // ask the local stored neighbor for the vertex number
                    const TetraCL* const neigh= t.GetNeighbor( face);
                    Assert( neigh, DROPSErrCL("GraphBuilderCL::buildNeighbors: (unref tetra) Missing neighbor"), DebugLoadBalC);
                    Assert( neigh->Unknowns.Exist( getVertexIdx()), DROPSErrCL("GraphBuilderCL::buildNeighbors: Tetrahedron without load balancing number"), DebugLoadBalC);
                    neigh_id= getGlobalVertexId( *neigh);
                }
                Assert( neigh_id>=0, DROPSErrCL("GraphBuilderCL::buildNeighbors: (unref tetra) No neighbor found!"), DebugLoadBalC);
                ew.update( neigh_id, fp);                   // remember the neighbor (and update the weight)
            }
        }
    }
    else {  // collect neighbors and weights for an refined tetrahedron
        TetraCL::const_ChildPIterator ch( t.GetChildBegin());
        for ( ; ch!=t.GetChildEnd(); ++ch){                 // check all children for neighbors
            if ( (*ch)->IsInTriang( triang_level_)){        // for each children, see documentation above!
                for ( Uint face=0; face<NumFacesC; ++face){
                    if ( !(*ch)->IsBndSeg( face)){
                        FaceCL const * fp = (*ch)->GetFace(face);
                        graph_index_type neigh_id= -1;
                        if ( fp->IsOnProcBnd()){
                            Assert( fp->Unknowns.Exist( getNeighVertexIdx()), DROPSErrCL("GraphBuilderCL::buildNeighbors: (ref tetra) Adjacency across process boundary not known"), DebugLoadBalC);
                            neigh_id= fp->Unknowns( getNeighVertexIdx());
                        }
                        else {
                            const TetraCL* const neigh= (*ch)->GetNeighInTriang( face, triang_level_);
                            Assert( neigh, DROPSErrCL("GraphBuilderCL::buildNeighbors: Missing neighbor"), DebugLoadBalC);
                            Assert( neigh->Unknowns.Exist( getVertexIdx()), DROPSErrCL("GraphBuilderCL::buildNeighbors: Tetrahedron without load balancing number"), DebugLoadBalC);
                            neigh_id= getGlobalVertexId( *neigh);
                        }
                        Assert( neigh_id>=0, DROPSErrCL("GraphBuilderCL::buildNeighbors: (ref tetra) No neighbor found!"), DebugLoadBalC);
                        ew.update( neigh_id, fp);
                    }
                }
            }
        }
    }
    ew.finalize();  // determine the weights

    // Store the neighbors and the weights in the data structure of the graph
    EdgeWeighterCL::const_neigh_iterator it= ew.begin();
    for ( ; it!=ew.end(); ++it, ++adj_count){
        getGraph()->adjncy()[ adj_count]= it->first;
        getGraph()->adjwgt()[ adj_count]= it->second;
    }
}

/** Iterate over all vertices and put their weight and the neighbors in the
    corresponding fields of the graph data structure.
*/
void GraphBuilderCL::buildVerticesAndEdges( VertexWeighterCL& vw, EdgeWeighterCL& ew)
{
    graph_index_type adj_count=0, vert_count= 0;    // position where to store the vertex weight and the adjacency
    LbIteratorCL it=LbIteratorCL::makeBegin( mg_, triang_level_);
    for ( ; it!= LbIteratorCL::makeEnd( mg_, triang_level_); ++it, ++vert_count){
        getGraph()->xadj()[vert_count]= adj_count;
        buildNeighbors( *it, adj_count, ew);
        getGraph()->vwght()[vert_count]= vw.getWeight( *it);
    }
    getGraph()->xadj()[vert_count]= adj_count;
}

/** Call the member functions to construct a graph. */
void GraphBuilderCL::build( BaseGraphCL*& graph, VertexWeighterCL& vw, EdgeWeighterCL& ew)
{
    graph_= new GraphCL();
    graph= graph_;
    buildVertexNumbering();
    buildVtxDist( num_vert_);
    CommunicateAdjacencyCL comm( *this, neighVertexIdx_);
    comm.Call();
    determineNumberAdjacencies();
    getGraph()->resize( num_adj_, num_vert_);
    buildVerticesAndEdges( vw, ew);
    comm.RemoveLbNumbersOnFaces();
    Comment( "Built a graph with " << getGraph()->vtxdist()[ProcCL::Size()] << " vertices!"<< std::endl, DebugLoadBalC);
}

// M E T I S  P A R T I T I O N E R  C L
//--------------------------------------

/** Graph is given on all processes. */
void MetisPartitionerCL::doParallelPartition()
{
    graph_index_type wgtflag    = 3,                 // Weights on vertices and adjacencies are given
                     numflag    = 0,                 // numbering of verts starts by 0 (C-Style)
                     ncon       = 1,                 // number of conditions
                     nparts     = ProcCL::Size(),    // number of sub-domains (per proc one)
                     options[5] = {0,0,0,0,0};       // default options and no debug information
    graph_real_type  itr        = 1000.0,           // how much an exchange costs
                     ubvec      = 1.05;             // allowed imbalance
    std::valarray<graph_real_type> tpwgts( 1.f/(graph_real_type)nparts, ProcCL::Size());
    MPI_Comm comm = MPI_COMM_WORLD;

    int result= 0;

    switch (method_) {
        case  3: std::cout << "bisection method of graph partitioning is not implemented, yet, using default instead" << std::endl;
        case -1: case 1: result = ParMETIS_V3_AdaptiveRepart(
                Addr(graph().vtxdist()), Addr(graph().xadj()), Addr(graph().adjncy()),
                Addr(graph().vwght()), Addr(graph().vwght()), Addr(graph().adjwgt()),
                &wgtflag, &numflag, &ncon, &nparts, Addr(tpwgts),
                &ubvec, &itr, options, &objective_, graph().part(), &comm); break;
        case  2: result = ParMETIS_V3_PartKway(
                Addr(graph().vtxdist()), Addr(graph().xadj()), Addr(graph().adjncy()),
                Addr(graph().vwght()), Addr(graph().adjwgt()),
                &wgtflag, &numflag, &ncon, &nparts, Addr(tpwgts),
                &ubvec, options, &objective_, graph().part(), &comm); break;
    }
    if ( result != METIS_OK){
        throw DROPSErrCL("Error in ParMETIS-function ParMETIS_V3_AdaptiveRepart");
    }
}

/** Graph is given on a single process. */
void MetisPartitionerCL::doSerialPartition()
{
    if ( !ProcCL::IamMaster()) return;

    graph_index_type nparts     = ProcCL::Size(),             // number of sub-domains (per proc one)
                     n          = graph().get_num_verts(),    // number of vertices
                     ncon       = 1,                          // number of conditions per vertex
                    *vsize      = 0;                          // default
    graph_index_type options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    graph_real_type    ubvec= 1.01;
    std::valarray<graph_real_type> tpwgts( 1.f/(graph_real_type)nparts, ProcCL::Size()*ncon);

    int result = 0;

    switch (method_) {
      case  1: std::cout << "adaptive re-computation of graph partitioning is not implemented, yet, using default instead" << std::endl;
      case -1: case 2: result= METIS_PartGraphKway(
              &n, &ncon, Addr(graph().xadj()), Addr(graph().adjncy()),
              Addr(graph().vwght()), vsize, Addr(graph().adjwgt()), &nparts, Addr(tpwgts), &ubvec, options,
              &objective_, graph().part()); break;
      case  3: result = METIS_PartGraphRecursive(
              &n, &ncon, Addr(graph().xadj()), Addr(graph().adjncy()),
              Addr(graph().vwght()), vsize, Addr(graph().adjwgt()), &nparts, Addr(tpwgts), &ubvec, options,
              &objective_, graph().part()); break;
    }
    if ( result != METIS_OK){
        throw DROPSErrCL("Error in METIS-function METIS_PartGraphKway");
    }
}

/** Perform the graph partitioning. */
void MetisPartitionerCL::doPartition()
{
    if ( graph().isSerial())
        doSerialPartition();
    else
        doParallelPartition();
}



/** Construct a class to determine a partitioning of a distributed triangulation
    hierarchy which is given by \a mg.
    \param mg the parallel tetrahedral hierarchy
*/
PartitioningCL::PartitioningCL( MultiGridCL& mg, int triang_level)
    : graph_( 0), partitioner_(0), mg_(mg),
      triang_level_( triang_level<0 ? mg.GetLastLevel() : triang_level),
      vertexIdx_( P0_FE)
{}

/** Determine a decomposition. \anchor DetermineDecompositionCL_make
    \param method this parameter specifies the strategy which is used to determine
           a partitioning. Therefore, four digits are used: P GT EW VW with
           <ul>
             <li> P the partitioner, i.e., identity(0), Metis(1), Zoltan(2), Scotch(3), Mondriaan(4 - not implemented so far)... </li>
             <li> GT type of the graph, i.e., graph(0) or hypergraph(1) </li>
             <li> EW method to weight the edges, i.e., unity(0), number of faces(1) or number of DOF(2) </li>
             <li> VW method to weight the vertices, i.e., unity(0), number of children on finest triangulation(1),
                  number of DOF(2), intersected tetras cause more work (3)</li>
           </ul>
*/
void PartitioningCL::make( const int method, int rho_I, const VecDescCL* lset, const BndDataCL<>* lsetbnd, const ObservedMigrateFECL* obs)
{
    EdgeWeighterCL* ew= 0;
    VertexWeighterCL *vw=0;
    BaseGraphBuilderCL* gb=0;

    // make the graph builder
    switch ( (method%1000)/100){
    case 0: gb= new GraphBuilderCL( mg_, triang_level_, vertexIdx_); break;
    default: throw DROPSErrCL("PartitioningCL::make: Unknown graph model");
    }

    // make the edge weighter
    switch ( (method%100)/10){
    case 0: ew= new EdgeWeighterCL(); break;
    case 1: ew= new EdgeWeighterTriangCL(); break;
    case 2:
        if ( obs==0){
            throw DROPSErrCL("PartitioningCL: To weight edges by DOF information, you have to provide DOF information");
        }
        ew= new EdgeWeighterDOFCL( *obs);
        break;
    default: throw DROPSErrCL("PartitioningCL::make: Unknown edge weighting method");
    }

    // make the vertex weighter
    switch ( method%10){
    case 0: vw= new VertexWeighterCL(); break;
    case 1: vw= new VertexWeighterTriangCL( triang_level_); break;
    case 3:
        if (lset==0 || lsetbnd==0){
            throw DROPSErrCL("PartitioningCL: To weight vertices by intersection metric, you have to provide the level set function!");
        }
        vw= new VertexWeighterTwoPhaseCL( triang_level_, rho_I, *lset, *lsetbnd);
        break;
    default: throw DROPSErrCL("PartitioningCL::make: Unknown vertex weighting method");
    }

    // build the graph
    gb->build( graph_, *vw, *ew);

    // make the partitioner
    switch (method/1000){
    case 0: partitioner_= new IdentityPartitionerCL( *graph_); break;
    case 1: partitioner_= new MetisPartitionerCL( *dynamic_cast<GraphCL*>(graph_)); break;
    default: throw DROPSErrCL("PartitioningCL::make: Unknown partitioner");
    }

    // determine a partitioning
    partitioner_->doPartition();

    // free memory
    delete ew;
    delete vw;
    delete gb;
}

/** Free the memory of the graph, partitioner and delete the graph vertex numbering.*/
void PartitioningCL::clear()
{
    delete graph_;       graph_=0;
    delete partitioner_; partitioner_=0;
    vertexIdx_.DeleteNumbering( mg_);
}

}   // end of namespace DROPS
