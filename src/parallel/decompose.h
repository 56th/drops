/// \file decompose.h
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

#include "geom/multigrid.h"
#include "misc/utils.h"
#include "misc/problem.h"
#include "parallel/migrateunknowns.h"
#include <vector>

#ifndef DROPS_DECOMPOSE_H
#define DROPS_DECOMPOSE_H

namespace DROPS{

/*****************************************************************************
This file contains the following classes:
    - LbIteratorCL, iterate over tetrahedra that represent a vertex in the
      graph, i.e.,
      * unrefined on a specified level
      * refined and children in the specified level
    - Graph classes for storing standard (and maybe in future hyper-) graphs
      * BaseGraphCL describes an interface for all graphs
      * GraphCL a standard graph used by, e.g., METIS
    - Partitioner classes to determine graph partitionings
      * BasePartitionerCL describes and interface for all partitioners
      * MetisPartitionerCL, the interface to the METIS family
==> - PartitioningCL class using all the methods above to determine
      a partitioning of a MultiGridCL.
*****************************************************************************/


/// \brief Class for iterating through the multi-nodes of this proc to set up a
/// reduced, dual graph
/** For numbering the tetrahedra for ParMETIS, we have to go through special
    tetrahedra. This class helps to touch only the right ones.
*/
class LbIteratorCL
{
  public:
    typedef MultiGridCL::TetraIterator IteratorT;     ///< Iteration through the tetrahedra of the multigrid

  private:
    MultiGridCL* mg_;                                 ///< pointer to the multigrid, that is used for iterating
    IteratorT    pos_;                                ///< The actual position of this iterator
    Uint         Level_,                              ///< actual level
                 TriLevel_;                           ///< the iterator should only go through this triangulation-level
                                                      ///< (normaly set to the last level)
  public:
    /// \brief Constructor
    LbIteratorCL (MultiGridCL* mg, IteratorT pos, Uint Level, Uint TriLevel) :
        mg_(mg), pos_(pos), Level_(Level), TriLevel_(TriLevel) {}

    /// \brief Copyconstructor
    LbIteratorCL (const LbIteratorCL &LbI) :
        mg_(LbI.mg_), pos_(LbI.pos_), Level_(LbI.Level_), TriLevel_(LbI.TriLevel_) {}

    // default destructor

    // Test if a tetra is in a loadbalancing set
    inline bool IsInLbSet() const                     ///< Test, if the tetrahedron, on which this iterator points, is in the loadbalancing set
        { return IsInLbSet(*pos_); }
    inline bool IsInLbSet( const TetraCL&) const;     ///< Test if a tetrahedron is in the loadbalancing set, see above

    // access of the tetraeder
    TetraCL& operator *  () const { return  *pos_; }  ///< return a reference onto the tetra, this iterator points to
    TetraCL* operator -> () const { return &*pos_; }  ///< return a pointer to a tetra

    // iterate through the elements of the LoadBalance-Set
    inline LbIteratorCL& operator ++ ();              ///< Go to the next multi-node tetra (prefix)
    inline LbIteratorCL  operator ++ (int);           ///< Go to the next multi-node tetra (suffix)

    // assignment
    LbIteratorCL& operator = (const LbIteratorCL& it) ///< assignment operator
    {
        mg_ = it.mg_; pos_=it.pos_;
        Level_=it.Level_; TriLevel_ = it.TriLevel_;
        return *this;
    }

    /// \name Factory method
    //@{
    static LbIteratorCL makeBegin( MultiGridCL&, int TriLevel= -1);
    static LbIteratorCL makeEnd(   MultiGridCL&, int TriLevel= -1);
    //@}

    /// \name comparing two iterators
    //@{
    /// \brief Test, if two iterators are the same by comparing the iterators onto the tetras
    friend bool operator == (const LbIteratorCL& It0, const LbIteratorCL& It1)
        { return It0.pos_ == It1.pos_; }
    /// \brief Test, if two iterators are not the same by comparing the iterators onto the tetras
    friend bool operator != (const LbIteratorCL& It0, const LbIteratorCL& It1)
        { return !(It0 == It1); }
    //@}
};


/// \brief Base class representing a (hyper)graph
/** This class is only the base class and should not be instanciated. However,
    this class also stores the partitioning. */
class BaseGraphCL
{
protected:
    std::vector<graph_index_type> part_;             ///< obtained partitioning
    int num_procs_w_verts_;

public:
    BaseGraphCL() : part_(0), num_procs_w_verts_(0) {}
    virtual ~BaseGraphCL() {}
    /// \brief
    bool isSerial() const { return num_procs_w_verts_ <= 1; }
    /// \brief Get partition number of a vertex
    graph_index_type getPartition( const graph_index_type v) const { return part_[v]; }
    /// \brief Get number of vertices stored on this process
    inline size_t get_num_verts() const { return part_.size(); }
    /// \brief Getter for partitioners
    graph_index_type* part() { return Addr(part_); }
};

/// \brief Class that represents a concrete graph for METIS and SCOTCH
/** Note: for procs without any graph vertex, an artificial vertex with weight 0
 *  and an artificial edge with weight zero is created to circumvent METIS' restrictions.
 */
class GraphCL : public BaseGraphCL
{
private:
    typedef std::vector<graph_index_type> IndexArray;     ///< type of storing an index array
    typedef std::vector<graph_index_type> WeightArray;    ///< type of storing an array of weights

    IndexArray  xadj_,                              ///< Starting index in the array adjncy, where node[i] writes its neighbors
                adjncy_,                            ///< adjacencies of the nodes
                vtxdist_;                           ///< number of nodes, that is stored by all procs
    WeightArray vwgt_,                              ///< weight of the Nodes
                adjwgt_;                            ///< weight of the edges

public:
    /// \brief Constructor
    GraphCL();
    /// \brief build the array vtxdist
    void buildVtxDist( graph_index_type num_vert);
    /// \brief resize the graph
    void resize( int num_adj, int num_vert);
    /// \brief free memory
    void clear();
    /// \brief Get number of adjacencies stored on this process
    inline size_t get_num_adj() const { return adjwgt_.size(); }
    /// \brief Get the first vertex stored by this proc
    inline graph_index_type get_first_vert() const { return vtxdist_[ ProcCL::MyRank()]; }
    /// \brief get the process storing a given vertex (by the global id)
    /// \todo check me!
    inline int getProc( graph_index_type globalidx)
        { return std::distance( vtxdist_.begin(), std::lower_bound(vtxdist_.begin(), vtxdist_.end(), globalidx)); }

    /// \name Getters and setters
    //@{
    IndexArray&  xadj()    { return xadj_; }
    IndexArray&  adjncy()  { return adjncy_; }
    IndexArray&  vtxdist() { return vtxdist_; }
    WeightArray& vwght()   { return vwgt_; }
    WeightArray& adjwgt()  { return adjwgt_; }
    //@}
};

/// \brief Base class representing the partitioners
/** The derived classes are responsible for partitioning a
    (hyper)graph. This class is an abstract class.*/
class BasePartitionerCL
{
protected:
    graph_index_type objective_; ///< resulting objective of the partitioner
    double balance_;             ///< resulting balancing constraint of the partitioner

public:
    BasePartitionerCL() : objective_(0), balance_(0.0) {}
    virtual ~BasePartitionerCL() {}
    /// \brief do the partitioning ( the result is stored in the graph)
    virtual void doPartition() = 0;
};

/// \brief Apply identity as partitioning. i.e., all tetrahedra stay on their process
class IdentityPartitionerCL : public BasePartitionerCL
{
private:
    BaseGraphCL& graph_;
    BaseGraphCL& graph() { return graph_; }

public:
    IdentityPartitionerCL( BaseGraphCL& graph)
        : BasePartitionerCL(), graph_( graph) {}
    void doPartition()
        { std::fill( graph().part(), graph().part()+graph().get_num_verts(), ProcCL::MyRank()); }
};

/// \brief Using (Par)Metis for partitioning a graph
/**
 * possible parameters for method:
 *   1: adaptive re-computation of graph partitioning (adaptive)
 *   2: multilevel method (KWay)
 *   3: bisection (recursive)
 *  -1: default (parallel: method 1, serial: method 2
 */
class MetisPartitionerCL : public BasePartitionerCL
{
private:
    GraphCL& graph_;       // the corresponding graph

    int method_;           // tells which method should be used to compute graph partition problem

    /// \brief Get the graph
    GraphCL& graph() { return graph_; }

    /// \brief Compute a partition for a distributed graph, i.e., use of ParMETIS_V3_AdaptiveRepart
    void doParallelPartition();
    /// \brief Compute a partition for a sequential graph, i.e., use of METIS_PartGraphKway
    void doSerialPartition();

public:
    MetisPartitionerCL( GraphCL& graph, int method = -1)
        : BasePartitionerCL(), graph_(graph), method_(method) {}
    /// \brief Use (Par)Metis functions to partition the given \a graph
    void doPartition();
};


/// \brief Using Zoltan for partitioning a graph
/** \todo Implement me. See parallel/partitioner.cpp for a reference implementation */
class ZoltanPartitionerCL : public BasePartitionerCL
{
    // THIS CLASS IS NOT IMPLEMENTED, YET.
};

/// \brief Using pt-Scotch for partitioning a graph
/** \todo Implement me. See parallel/partitioner.cpp for a reference implementation */
class ScotchPartitionerCL : public BasePartitionerCL
{
    // THIS CLASS IS NOT IMPLEMENTED, YET.
};


/// \brief Class for determining a partitioning of a distributed hierarchy of triangulations.
///        \anchor DetermineDecompositionCL
/** This class provides an interface to determine a partitioning of a distributed
    hierarchy of triangulations. Therefore, this class internally uses a graph model which
    can be specified by the method parameter given by the factory method "make".
*/
class PartitioningCL
{
private:
    BaseGraphCL*       graph_;          ///< storing a (hyper)graph
    BasePartitionerCL* partitioner_;    ///< storing a partitioner
    MultiGridCL&       mg_;             ///< store a reference to the tetrahedral hierarchy
    Uint               triang_level_;   ///< triangulation level that should be decomposed
    IdxDescCL          vertexIdx_;      ///< Used to number tetrahedra corresponding to graph vertices

public:
    /// \brief Constructor
    PartitioningCL( MultiGridCL& mg, int triang_level=-1);
    /// \brief Destructor cleans everything up (and deletes the partitioning)
    ~PartitioningCL() { clear(); }
    /// \brief Determine a partitioning
    void make( const int method, int rho_I=11,
               const VecDescCL* lset=0, const BndDataCL<>* lsetbnd=0,
               const ObservedMigrateFECL* obs= 0);
    /// \brief get the destination process of a tetrahedron
    int getDestination( const TetraCL& t) const
        { return graph_->getPartition( t.Unknowns( vertexIdx_.GetIdx())); }
    /// \brief Clean up and free the index used to number the tetrahedra
    void clear();
};

}

#include "parallel/decompose.tpp"

#endif

