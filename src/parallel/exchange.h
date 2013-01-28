/// \file exchange.h
/// \brief handling of a parallel distributed vectors and distributed matrices
/// \author LNM RWTH Aachen: ; SC RWTH Aachen: Oliver Fortmeier

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

/// These classes do not use the DDD-Interfaces. After the lists
/// are created no geometric data are needed to do the
/// accumulation in opposite to the DDD-Interface. And this class
/// split the send and the receive, so other work can be done
/// between these commands.

#ifndef DROPS_EXCHANGE_H
#define DROPS_EXCHANGE_H

#ifdef _PAR
#include "parallel/parallel.h"
#endif
#include <list>
#include <vector>
#include "num/spmat.h"
#include "geom/multigrid.h"
#include "misc/problem.h"

namespace DROPS{

class DummyExchangeCL
{
  public:
    /// \brief Accumulate a vector
    void Accumulate( VectorCL&, const Ulint= 0) const {}
    /// \brief Accumulate a vector of vectors
    void Accumulate( std::vector<VectorCL>&) const {}
    /// \brief Get an accumulated copy of a vector
    VectorCL GetAccumulate (const VectorCL& u) const {return u;}
    /// \brief Get accumulated version of a vector of vectors
    std::vector<VectorCL> GetAccumulate( const std::vector<VectorCL>& u) const { return u; }

    /// \brief Parallel inner product without final reduction over all processes
    double LocalDot( const VectorCL& x, bool, const VectorCL& y, bool, VectorCL* x_acc=0, VectorCL* y_acc=0) const
    {
        if (x_acc != 0)
            *x_acc = x;
        if (y_acc != 0)
            *y_acc = y;

        return KahanInnerProd( x, y, double());
    }
    /// \brief Parallel inner product of vectors with final reduction over all processes
    ///        x_acc must not point to x and y_acc must not point to y
    double ParDot( const VectorCL& x, bool x_is_acc, const VectorCL& y, bool y_is_acc, VectorCL* x_acc=0, VectorCL* y_acc=0) const
    {
        return LocalDot(x, x_is_acc, y, y_is_acc, x_acc, y_acc);
    }
    /// \brief Parallel squared Euclidian norm without final reduction over all processes
    double LocalNorm_sq( const VectorCL& x, bool is_acc, VectorCL* x_acc=0) const{
        return LocalDot(x, is_acc, x, is_acc, x_acc);
    }
    /// \brief Parallel squared Euclidian norm with final reduction over all processes
    double Norm_sq( const VectorCL& x, bool is_acc, VectorCL* x_acc=0) const{
        return LocalNorm_sq(x, is_acc, x_acc);
    }
    /// \brief Parallel Euclidian norm with final reduction over all processes
    double Norm( const VectorCL& x, bool is_acc, VectorCL* x_acc=0) const{
        return std::sqrt(LocalNorm_sq(x, is_acc, x_acc));
    }
};

class DummyExchangeBlockCL : public DummyExchangeCL
{
  public:
    typedef std::vector<const DummyExchangeCL*>    ExchangeCLCT;    ///< Container for ExchangeCL

  private:
    ExchangeCLCT         exchange_;     ///< store all index describers to access ExchangeCLs

    /// \brief Update of datastructure, i.e. blockoffset_
    void Update() {}

  public:
    DummyExchangeBlockCL() {}

    /// \brief Attach an index describer
    void AttachTo(const DummyExchangeCL& ex) {
        exchange_.push_back(&ex);
    }
    /// \brief Ask for number of handled blocks
    size_t GetNumBlocks() const { return exchange_.size(); }
    /// \brief Ask for an ExchangeCL
    const DummyExchangeCL& GetEx( size_t i) const { return *exchange_[i]; }
};

#ifdef _PAR

/// fwd declaration
class ExchangeBuilderCL;
class ExchangeMatrixCL;

/******************************************************************************
* S E N D  N U M  D A T A  C L                                                *
******************************************************************************/
/// \brief Class for sending distributed vector entries to neighbor processes
/// \tparam T type of the entries in the vector to be sent.
template <typename T>
class SendNumDataCL
{
  public:
    typedef T value_type;

  private:
    friend class ExchangeBuilderCL;     ///< for generating the data structures
    friend class ExchangeMatrixCL;      ///< for generating the data structures
    int toproc_;                        ///< to whom the data are send
    int minlengthvec_;                  ///< minimal length of the vector (for debugging)
    ProcCL::DatatypeT mpidatatype_;     ///< derived from MPI datatype indexed

    /// \brief Create the MPI datatype
    void CreateDataType(const int, const int bl, const int ad[]);
    /// \brief Free the MPI datatype
    void freeType();

  public:
    /// \brief Construct a class for sending numerical data
    SendNumDataCL( int toproc) : toproc_(toproc), mpidatatype_( ProcCL::NullDataType) {}
    /// \brief Delete this class, i.e., free the MPI datatype
    ~SendNumDataCL() { freeType(); }
    /// \brief Ask for the receiver-rank
    int GetReceiver() const { return toproc_; }
    /// \brief Send the data (nonblocking, asynchronous)
    template <typename VectorT>
    inline ProcCL::RequestT Isend(const VectorT&, int tag, Ulint offset) const;
     // Send data to "toProc_" (nonblocking, asynchronous)
    inline ProcCL::RequestT Isend(const double*, int tag, Ulint offset) const;
};

/******************************************************************************
* R E C V  N U M  D A T A  C L                                                *
******************************************************************************/
/// \brief Class for receiving and accumulating data
/// \tparam T type of the entries in the vector to be received.
template <typename T>
class RecvNumDataCL
{
  public:
    typedef T value_type;

  private:
    friend class ExchangeBuilderCL; ///< for generating the data structures
    int fromproc_;                  ///< who is the sender
    std::vector<IdxT> sysnums_;     ///< where to add the received data

  public:
    /// \brief Construct a class for receiving and accumulating data
    RecvNumDataCL( int fromproc) : fromproc_(fromproc) {}
    // default copy- and destructor

    /// \brief Ask for the sender
    int GetSender() const { return fromproc_; }
    /// \brief Receive data
    inline ProcCL::RequestT Irecv( int tag, VectorBaseCL<T>& recvBuf, Ulint offset) const;
    /// \brief Accumulate the received data (adding values from \a recvBuf to \a x)
    inline void Accumulate( VectorBaseCL<T>& x, Ulint offsetV, const VectorBaseCL<T>& recvBuf, Ulint offsetRecv) const;
    /// \brief Assign the received data ( by taking values from \a recvBuf)
    inline void Assign( VectorBaseCL<T>& x, Ulint offsetV, const VectorBaseCL<T>& recvBuf, Ulint offsetRecv) const;
};


// fwd declaration
class ExchangeBlockCL;

/******************************************************************************
* E X C H A N G E  C L                                                        *
******************************************************************************/
/// \brief Accumulate vectors, compute norm, and determine inner product of
///        parallel distributed vectors
/** This class is one of the core classes to handle distributed vectors. These
    vectors can be given in two different forms:

    (a) Accumulated form: Each process, storing a dof, knows the complete value
    of the dof.
    (b) Distributed form: Each process, storing a dof, knows only a part of the
    value of the dof. The global value is then given by the sum over all
    partial values on the processes.

    This class provides three main functionalities:

    - (1) Transforming a vector of type (b) to type (a).
    - (2) Performing an inner product of two vectors
    - (3) Getting information about the number of a dof on another process.

    (1) Transforming a vector of form (b) to (a) is called "accumulation of a
    vector." This task is the main functionality of this class and is
    implemented by the function <b>void Accumulate( VectorCL&) const</b>.

    The accumulation process is two-fold. In the first communication phase,
    each process holding a partial value sends this value to the DoF owner of the
    corresponding simplex (vertex or edge). Afterwards, the DoF owner computes the
    sum. In the second communication phase, the DoF owner informs the copies about
    the global value. Note that the DoF owner may be different from the geometric
    owner used in DiST.

    (2) The second core functionality of this class is determining the inner
    product of two vectors, given by the function <b>
    double ParDot( const VectorCL& x, bool isXacc,
            const VectorCL& y, bool isYacc,
            VectorCL* x_acc, VectorCL* y_acc) const</b>.
    We refer to the description of this function for detailed
    information about the function parameters. Note, that the sums are
    determined by the Kahan sum formula.

    (3) The third functionality of this class is to get information about a
    dof on another process. Therefore, we implemented the function <b>
    DOFInfoList_const_iterator GetProcListBegin( IdxT localdof) const</b>
    and<b>
    DOFInfoList_const_iterator GetProcListEnd( IdxT localdof) const</b>.
    With these functions, one can access a list over all all processes
    storing the local dof as well and get information about their local dof.

    \todo Right now, if different positive number of unknowns exists for vertices
    edges, faces or tetrahedra, this class does not work correctly. And, in
    particular, if unknowns do not exist on vertices but on another type of
    simplices, this class breaks down as well.
    \todo Implement communication pattern with a single communication phase
    \todo Test with extended dofs
*/
class ExchangeCL
{
  public:
    typedef std::vector<IdxT>                  IdxVecT;
    // neighbor information
    typedef DROPS_STD_UNORDERED_SET<int>       NeighListT;      ///< \todo Correct data-type?
    typedef NeighListT::const_iterator         NeighListT_const_iterator;
    // distributed DOF information
    typedef std::map<int,IdxT>                 DOFInfoT;
    typedef std::vector< DOFInfoT >            DOFProcListT;
    typedef DOFInfoT::const_iterator           DOFInfoList_const_iterator;

  private:
    // internal data containers
    typedef std::list< SendNumDataCL<double> > SendListT;       ///< type for storing send information
    typedef std::list< RecvNumDataCL<double> > RecvListT;       ///< type for storing receive information
    typedef std::vector<VectorCL>              BufferListT;     ///< list of buffers for MPI Recv
    typedef std::vector<ProcCL::RequestT>      RequestListT;    ///< list of MPI Requests

  private:
    friend class ExchangeBuilderCL;
    friend class ExchangeBlockCL;
    friend class IdxDescCL;

    bool viaowner_;                                         ///< flag, if the communication is twofold

    SendListT sendListPhase1_;                              ///< send information for the first communication phase
    SendListT sendListPhase2_;                              ///< send information for the second communication phase
    RecvListT recvListPhase1_;                              ///< receive information for the first communication phase
    RecvListT recvListPhase2_;                              ///< receive information for the second communication phase

    /// \brief Receive buffers two vectors x and y
    mutable BufferListT xBuf_, yBuf_;

    NeighListT   neighs_;                                   ///< neighbor processes
    DOFProcListT dofProcList_;                              ///< storing information about distributed dof

    /// \brief Call MPI_Isend and MPI_Irecv
    void InitComm(int Phase, const VectorCL&, ProcCL::RequestT* sendreq, ProcCL::RequestT* recvreq, BufferListT& buf, int tag=1001, Ulint offset=0) const;
    /// \brief Let the DoF owner do all accumulations
    void DoAllAccumulations( VectorCL&, const BufferListT& buf, const Ulint offset=0) const;
    /// \brief Let the copies assign the accumulated values
    void DoAllAssigning( VectorCL& v, const BufferListT& buf, const Ulint offset=0) const;

    /// \brief Do the computation of an inner product with different number of accumulations
    //@{
    double LocalDotNoAccumulation(const VectorCL&, const VectorCL&) const;
    double LocalDotOneAccumulation(const VectorCL& x, const VectorCL& y, VectorCL* y_acc) const;
    double LocalNormSQAccumulation( const VectorCL& x, VectorCL* x_acc) const;
    double LocalDotTwoAccumulations(const VectorCL& x, const VectorCL& y, VectorCL* x_acc, VectorCL* y_acc) const;
    //@}

  public:
    /// \todo Replace LocalIndex, DistrIndex, OwnerDistrIndex by Owner vector, storing the DoF owner for distributed dof, or -1 for non-distributed dof.
    IdxVecT LocalIndex;                 ///< list of local dof only locally existing
    IdxVecT DistrIndex;                 ///< list of local dof existing on at least two processes
    IdxVecT OwnerDistrIndex;            ///< list of local dof, this process owns the corresponding simplex

    ExchangeCL( bool viaowner=true) : viaowner_(viaowner) {}

    /// \brief Factory method. This function calls an ExchangeBuilderCL to build this class
    void CreateList( const MultiGridCL& mg, IdxDescCL*, bool, bool);
    /// \brief Clear the memory
    void clear();
    /// \brief Get size of the vector this class is responsible for
    size_t GetNum() const { return LocalIndex.size()+DistrIndex.size(); }
    /// \brief Get the number of local and owned DoF's
    size_t GetNumOwned() const { return LocalIndex.size()+OwnerDistrIndex.size(); }
    /// \brief Check if the communication is performed via the DoF owner
    bool CommViaOwner() const { return viaowner_; }

    /// \brief Accumulate a vector
    void Accumulate( VectorCL&, const Ulint offset=0) const;
    /// \brief Accumulate a vector of vectors
    void Accumulate( std::vector<VectorCL>&) const;
    /// \brief Get an accumulated copy of a vector
    VectorCL GetAccumulate (const VectorCL&) const;
    /// \brief Get accumulated version of a vector of vectors
    std::vector<VectorCL> GetAccumulate( const std::vector<VectorCL>&) const;

    /// \brief Parallel inner product without final reduction over all processes
    double LocalDot( const VectorCL&, bool, const VectorCL&, bool, VectorCL* x_acc=0, VectorCL* y_acc=0) const;
    /// \brief Parallel inner product of vectors with final reduction over all processes
    ///        x_acc must not point to x and y_acc must not point to y
    double ParDot( const VectorCL&, bool, const VectorCL&, bool, VectorCL* x_acc=0, VectorCL* y_acc=0) const;
    /// \brief Parallel squared Euclidian norm without final reduction over all processes
    double LocalNorm_sq( const VectorCL&, bool, VectorCL* x_acc=0) const;
    /// \brief Parallel squared Euclidian norm with final reduction over all processes
    double Norm_sq( const VectorCL&, bool, VectorCL* x_acc=0) const;
    /// \brief Parallel Euclidian norm with final reduction over all processes
    double Norm( const VectorCL&, bool, VectorCL* x_acc=0) const;

    /// \name Get information about neighbor processes
    //@{
    Uint GetNumNeighs() const                                           ///< Get number of neighbor processes
        { return neighs_.size(); }
    NeighListT_const_iterator GetNeighBegin() const                     ///< Get iterator to the first neighbor process
        { return neighs_.begin(); }
    NeighListT_const_iterator GetNeighEnd() const                       ///< Get iterator behind the last neighbor process
        { return neighs_.end(); }
    //@}

    /// \name Get information about a given dof on other processes
    //@{
    DOFInfoList_const_iterator GetProcListBegin( IdxT localdof) const   ///< Get begin iterator of the process list of \a localdof
        { return dofProcList_[localdof].begin(); }
    DOFInfoList_const_iterator GetProcListEnd( IdxT localdof) const     ///< Get end iterator of the process list of \a localdof
        { return dofProcList_[localdof].end(); }
    bool IsDist(IdxT localdof) const                                    ///< Check if \a localdof is distributed
        { return !dofProcList_[localdof].empty(); }
    inline bool IsOnProc( IdxT localdof, int proc) const;               ///< Check if \a proc owns a copy of \a localdof
    size_t GetNumProcs( IdxT localdof) const                            ///< Get the number of processes owning a copy of \a localdof
        { return dofProcList_[localdof].size(); }
    inline IdxT GetExternalIdxFromProc( const IdxT localdof, int proc) const;   ///< Get dof number of proc
    //@}

    /// \todo Make more efficient?
    bool AmIOwner( IdxT dof) const {
    	if ( GetNumProcs( dof) == 0)
    		return true;
    	for ( size_t i = 0; i< OwnerDistrIndex.size(); ++i)
    		if (OwnerDistrIndex[i] == dof)
    			return true;

    	return false;
    }
};


/****************************************************************************
* E X C H A N G E  B L O C K  C L A S S                                     *
****************************************************************************/
/// \brief Handle exchange of all numerical data for a blocked vector, i.e.
///    vectors, that have multiple IdxDescCL
/** This class handles the exchange of a blocked vector containing multiple
    describers. This is used to perform a blocked version of an iterative
    solver. For instance, GCR can be used to solve the Oseen problem.
 */
/****************************************************************************
* E X C H A N G E  B L O C K  C L A S S                                     *
****************************************************************************/
class ExchangeBlockCL
{
  public:
    typedef std::vector<const ExchangeCL*>         ExchangeCLCT;    ///< Container for ExchangeCL
    typedef std::vector<IdxT>                      BlockOffsetCT;   ///< Container of starting index of block elements
    typedef std::vector<ExchangeCL::RequestListT>  RequestListT;    ///< List of list of MPI requests
    typedef std::vector<ExchangeCL::BufferListT*>  BufferListT;     ///< List of list of buffers for MPI Recv

  private:
    ExchangeCLCT         exchange_;     ///< store all ExchangeCLs
    BlockOffsetCT        blockOffset_;  ///< store the length of vectors

    mutable BufferListT xBuf_, yBuf_;

    /// \brief Call MPI_Isend and MPI_Irecv
    void InitComm(int Phase, const VectorCL&, RequestListT& sendreq, RequestListT& recvreq, BufferListT& buf, int tag=1001) const;
    /// \brief Let the DoF owner do all accumulations
    void DoAllAccumulations( VectorCL&, const BufferListT&) const;
    /// \brief Let the copies assign the accumulated values
    void DoAllAssigning( VectorCL& v, const BufferListT&) const;

    /// \name Do the computation of an inner product with different number of accumulations
    //@{
    double LocalDotNoAccumulation(const VectorCL&, const VectorCL&) const;
    double LocalDotOneAccumulation(const VectorCL& x, const VectorCL& y, VectorCL* y_acc) const;
    double LocalNormSQAccumulation( const VectorCL& x, VectorCL* x_acc) const;
    double LocalDotTwoAccumulations(const VectorCL& x, const VectorCL& y, VectorCL* x_acc, VectorCL* y_acc) const;
    //@}

    /// \brief Update of data structure, i.e. blockoffset_
    void Update();

  public:
    ExchangeBlockCL() {}

    /// \brief Attach an ExchangeCL
    void AttachTo(const ExchangeCL&);
    /// \brief Ask for number of handled blocks
    size_t GetNumBlocks() const { return exchange_.size(); }
    /// \brief Ask for length of vectors that can be accumulated
    IdxT GetNum() const { return blockOffset_.back(); }
    /// \brief Ask for an ExchangeCL
    const ExchangeCL& GetEx( size_t i) const { return *exchange_[i]; }
    /// \brief Get the offset for block \a i
    size_t GetOffset( const size_t i) const { return blockOffset_[i]; }

    /// \brief Accumulate a vector
    void Accumulate( VectorCL&) const;
    /// \brief Accumulate a vector of vectors
    void Accumulate( std::vector<VectorCL>&) const;
    /// \brief Get an accumulated copy of a vector
    VectorCL GetAccumulate (const VectorCL&) const;
    /// \brief Get accumulated version of a vector of vectors
    std::vector<VectorCL> GetAccumulate( const std::vector<VectorCL>&) const;

    /// \brief Parallel inner product without final reduction over all processes
    double LocalDot( const VectorCL&, bool, const VectorCL&, bool, VectorCL* x_acc=0, VectorCL* y_acc=0) const;
    /// \brief Parallel inner product of vectors with final reduction over all processes.
    ///        x_acc must not point to x and y_acc must not point to y
    double ParDot( const VectorCL&, bool, const VectorCL&, bool, VectorCL* x_acc=0, VectorCL* y_acc=0) const;
    /// \brief Parallel squared Euclidian norm without final reduction over all processes
    double LocalNorm_sq( const VectorCL&, bool, VectorCL* x_acc=0) const;
    /// \brief Parallel squared Euclidian norm with final reduction over all processes
    double Norm_sq( const VectorCL&, bool, VectorCL* x_acc=0) const;
    /// \brief Parallel Euclidian norm with final reduction over all processes
    double Norm( const VectorCL&, bool, VectorCL* x_acc=0) const;
};


/****************************************************************************
* E X C H A N G E  M A T R I X  C L A S S                                   *
****************************************************************************/
/// \brief Handle the accumulation of a sparse matrix (MatrixCL)
/** This class is capable of determining the communication pattern for
    accumulating a sparse matrix, and performing the accumulation.
    \todo(par) Develop an "accure" version of accumulation
    \todo TestMe
 */
/****************************************************************************
* E X C H A N G E  M A T R I X  C L A S S                                   *
****************************************************************************/
class ExchangeMatrixCL
{
  public:
    typedef std::vector<int>       ProcNumCT;       ///< Container for storing neighbor processes
    typedef ProcNumCT::iterator    ProcNum_iter;    ///< iterator of ProcNumCT
    typedef std::vector<size_t>    CouplingCT;      ///< Container of distributed matrix elements

  private:

    /// each element of ExList_ handles the send-operation with a single neighbor processor
    std::vector< SendNumDataCL<double> > ExList_;
    /// Buffer for receiving elements
    std::vector<VectorCL>              RecvBuf_;
    /// Where to add/store received non-zeroes
    std::vector<CouplingCT>            Coupl_;
    /// flag, if non-zero is not stored on local processor
    static size_t NoIdx_;

    /// Determine the intersection of two processor lists
    inline ProcNumCT Intersect( const ExchangeCL& RowEx, const ExchangeCL& ColEx, const size_t i, const size_t j) const;

    /// Determine the position, where a nonzero is stored
    inline size_t GetPosInVal(const size_t row, const size_t col, const MatrixCL& mat)
    /// if the non-zero (row,col) is not stored by the local processor, this function
    /// returns NoIdx_
    {
        Assert( row<mat.num_rows() && col<mat.num_cols(), DROPSErrCL("ExchangeMatrixCL::GetPosInVal: Row or col out of bounds"), DebugParallelNumC);
        const size_t *pos= std::lower_bound( mat.GetFirstCol(row), mat.GetFirstCol(row+1), col);
        return (pos != mat.GetFirstCol(row+1) && *pos==col) ? pos-mat.GetFirstCol(0) : NoIdx_;
    }

  public:
    // default constructors and destructors

    /// \brief Reset
    void Clear() { ExList_.clear(); RecvBuf_.clear(); Coupl_.clear(); }

    /// \brief Determine the communication pattern for accumulating a matrix
    void BuildCommPattern(const MatrixCL& mat, const IdxDescCL& RowIdx, const IdxDescCL& ColIdx){
        BuildCommPattern(mat, RowIdx.GetEx(), ColIdx.GetEx());
    }

    /// \brief Determine the communication pattern for accumulating a matrix
    void BuildCommPattern(const MatrixCL&, const ExchangeCL& RowEx, const ExchangeCL& ColEx);

    /// \brief Accumulate a matrix
    MatrixCL Accumulate(const MatrixCL&);
};

/******************************************************************************
* E X C H A N G E  B U I L D E R  C L                                         *
******************************************************************************/
/// \brief Build an ExchangeCL for a given MultiGridCL and IdxDescCL
class ExchangeBuilderCL
{
  private:
    typedef SArrayCL<int,3> int3T;
    typedef std::vector<int3T> tmpRecvT;    ///< for each process, store (proc, remote dof, remote extended dof) temporarily

    ExchangeCL&               ex_;      ///< ExchangeCL to  be build
    const MultiGridCL&        mg_;      ///< The underlying multigrid
    IdxDescCL&                rowidx_;  ///< corresponding index describer
    DiST::InterfaceCL*        interf_;  ///< used for internal interface communication
    static int                NoInt_;   ///< flag for undefined send position \todo negative number?

    // Copy constructor and assignment is not allowed and not implemented
    ExchangeBuilderCL( const ExchangeBuilderCL&);
    ExchangeBuilderCL& operator=( const ExchangeBuilderCL&);

    /// \brief Build the data structures needed to communicate via the DoF owner
    ///        process
    void buildViaOwner();
    /// \brief Build the data structures where all processes communicate
    ///        directly with their neighbors
    void buildNtoN();
    /// \name Helper functions for the build
    //{@
    /// \brief Determine DOF owner among all procs, which hold the local dof
    static int GetDOFOwner( const tmpRecvT& rcvTmp);
    /// \brief Clear all information in the ExchangeCL
    void clearEx();
    /// \brief Build index lists
    void BuildIndexLists();
    /// \brief Allocate memory for the receive buffers
    void BuildRecvBuffer();
    //@}

    // handlers
    class HandlerDOFExchangeCL;
    class HandlerDOFSendCL;
    class HandlerDOFRecvCL;
    class HandlerDOFNtoNSendCL;
    class HandlerDOFNtoNRecvCL;
    class HandlerDOFIndexCL;

  public:
    ExchangeBuilderCL( ExchangeCL& ex, const MultiGridCL& mg, IdxDescCL *RowIdx);
    ~ExchangeBuilderCL() { if (interf_) delete interf_; interf_=0; }

    /// \brief Build the ExchangeCL specified in the constructor
    void build() { ex_.viaowner_ ? buildViaOwner() : buildNtoN(); }
};

#endif // parallel

} // end of namespace DROPS

// File, where the inline an template-functions are declared
#include "parallel/exchange.tpp"

#endif
