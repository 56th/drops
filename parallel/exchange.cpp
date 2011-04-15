/// \file exchange.cpp
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

#ifdef _PAR
#  include "parallel/exchange.h"
#  include "parallel/parmultigrid.h"
#endif
#include <iomanip>
#include <map>
#include <limits>

namespace DROPS{

void ExchangeCL::InitComm(
    int Phase, const VectorCL& v,
    ProcCL::RequestT* sendreq, ProcCL::RequestT* recvreq,
    BufferListT& buf, int tag, Ulint offset) const
/** Call of the non-blocking MPI send and receive operations.
    \param Phase Flag, if the send and receive calls should be performed
        for the 1 or 2 phase.
    \param v This vector is used as a send buffer. So make sure, that its
        entries are not modified before calling the Waitall function for
        the send requests.
    \param sendreq list of requests, where to store the requests return by
        the MPI Isend function
    \param recvreq list of requests, where to store the requests return by
        the MPI Irecv function
    \param buf list of buffers to store the received values. Make sure, not
        to modify this buffer before calling the Waitall function for the
        receive requests
    \param tag used tag for the send and receive operation
    \param offset This offset is used to access elements in the vector \a v.
*/
{
    // set iterators for sending
    SendListT::const_iterator sendit=
        (Phase==1) ? sendListPhase1_.begin() : sendListPhase2_.begin();
    SendListT::const_iterator sendend=
        (Phase==1) ? sendListPhase1_.end() : sendListPhase2_.end();
    size_t i=0;
    // init send
    for ( ; sendit!=sendend; ++sendit, ++i){
        sendreq[i]= sendit->Isend( v, tag, offset);
    }

    // set iterators for receiving
    RecvListT::const_iterator recvit=
        (Phase==1) ? recvListPhase1_.begin() : recvListPhase2_.begin();
    RecvListT::const_iterator recvend=
        (Phase==1) ? recvListPhase1_.end() : recvListPhase2_.end();
    BufferListT::iterator bufit= buf.begin();
    // init recv
    i=0;
    for ( ; recvit!=recvend; ++recvit, ++bufit, ++i){
        recvreq[i]= recvit->Irecv( tag, *bufit, 0);
    }
}

void ExchangeCL::DoAllAccumulations(
    VectorCL& v, const BufferListT& buf, const Ulint offset) const
/** Owner does all the accumulation with information about all neighbors.
    \pre Receive for the phase I has to be completed and if v is used as a send
         buffer, the sending has to be completed as well.
    \param v Store the accumulated values in this vector
    \param buf Receive buffer where the values of the copies can be found
    \param offset This offset is used to access elements in the vector \a v.
*/
{
    RecvListT::const_iterator recvit= recvListPhase1_.begin();
    BufferListT::const_iterator bufit= buf.begin();
    for ( ; recvit!=recvListPhase1_.end(); ++recvit, ++bufit){
        recvit->Accumulate( v, offset, *bufit, 0);
    }
}

void ExchangeCL::DoAllAssigning(
    VectorCL& v, const BufferListT& buf, const Ulint offset) const
/** Copies does the assigning of the values determined by the owner
    \pre Receive for the phase II has to be completed and if v is used as a send
         buffer, the sending has to be completed as well.
    \param v Store the accumulated values in this vector
    \param buf Receive buffer where the values of the owner can be found
    \param offset This offset is used to access elements in the vector \a v.
*/
{
    RecvListT::const_iterator recvit=  recvListPhase2_.begin();
    BufferListT::const_iterator bufit= buf.begin();
    for ( ; recvit!=recvListPhase2_.end(); ++recvit, ++bufit){
        recvit->Assign( v, offset, *bufit, 0);
    }
}

double ExchangeCL::LocalDotNoAccumulation(const VectorCL& x, const VectorCL& y) const
/** Both vectors \a x and \a y are given in accumulated form.
    \return inner product of x and y without reducing this value over
        all processes.
*/
{
    double sum=0;
    sum=KahanInnerProd( x, y, LocalIndex.begin(), LocalIndex.end(), double());
    sum=KahanInnerProd( x, y, OwnerDistrIndex.begin(), OwnerDistrIndex.end(), sum);
    return sum;
}

double ExchangeCL::LocalDotOneAccumulation(
    const VectorCL& x, const VectorCL& y, VectorCL* y_acc) const
/** The vector \a x is given in accumulated form and \a y is given in distributed form.
    \return inner product of x and y without reducing this value over
        all processes.
*/
{
    double result= std::numeric_limits<double>::max();
    bool y_acc_created=false;

    // Check if memory for the accumulated form of y is provided
    if ( y_acc){
        *y_acc= y;
    }
    else{
        y_acc_created= true;
        y_acc= new VectorCL(y);
    }

    // Determine, if a second communication step (and assigning) is neccessary
    const bool doSecondComm= !y_acc_created && viaowner_;

    // number of send and receive for phase 1 and 2
    const size_t num_sr_1=
        sendListPhase1_.size() + recvListPhase1_.size();
    const size_t num_sr_2=
         sendListPhase2_.size() + recvListPhase2_.size();
    RequestListT req(num_sr_1+num_sr_2);

    // initiate the communication for phase I
    InitComm( 1, y, Addr(req), Addr(req)+sendListPhase1_.size(), yBuf_);

    // While communicating, do product on local elements
    const double sum1=KahanInnerProd( x, y, LocalIndex.begin(), LocalIndex.end(), double());

    // Do accumulation on owners, therefore, first, wait until all
    // messages are received by the owner
    ProcCL::WaitAll( recvListPhase1_.size(), Addr(req)+sendListPhase1_.size());
    DoAllAccumulations( *y_acc, yBuf_);

    // Init Second transfer phase only if accumulated version of y is requested.
    // Here, the buffer for y can be re-used.
    if ( doSecondComm)
        InitComm( 2, *y_acc, Addr(req)+num_sr_1, Addr(req)+num_sr_1+sendListPhase2_.size(), yBuf_);

    // While communication accumulated values, do product on distributed elements
    result= KahanInnerProd( x, *y_acc, OwnerDistrIndex.begin(), OwnerDistrIndex.end(), sum1);

    // Before touching the memory of y_acc and returning y, wait
    // until send operation and, eventually, receive operation are completed.
    ProcCL::WaitAll( sendListPhase1_.size(), Addr(req));
    if ( doSecondComm){
        ProcCL::WaitAll( num_sr_2, Addr(req)+num_sr_1);
        DoAllAssigning( *y_acc, yBuf_);
    }

    // Free memory, if allocated here
    if ( y_acc_created){
        delete y_acc; y_acc=0;
    }

    return result;
}

double ExchangeCL::LocalNormSQAccumulation( const VectorCL& x, VectorCL* x_acc) const
/** The vector \a x is given in distributed form.
    \return inner product of x and x without reducing this value over
        all processes.
*/
{
    double result= std::numeric_limits<double>::max();
    bool x_acc_created=false;

    // Check if memory for the accumulated form of y is provided
    if ( x_acc){
        *x_acc= x;
    }
    else{
        x_acc_created= true;
        x_acc= new VectorCL(x);
    }
    const bool doSecondComm= !x_acc_created && viaowner_;

    // number of send and receive for phase 1 and 2
    const size_t num_sr_1=
        sendListPhase1_.size() + recvListPhase1_.size();
    const size_t num_sr_2=
         sendListPhase2_.size() + recvListPhase2_.size();
    RequestListT req(num_sr_1+num_sr_2);

    // initiate the communication for phase I
    InitComm( 1, x, Addr(req), Addr(req)+sendListPhase1_.size(), xBuf_, 1001);

    // While communicating, do product on local elements
    const double sum1=KahanInnerProd( x, x, LocalIndex.begin(), LocalIndex.end(), double());

    // Do accumulation on owners
    ProcCL::WaitAll( recvListPhase1_.size(), Addr(req)+sendListPhase1_.size());
    DoAllAccumulations( *x_acc, xBuf_);

    // Init Second transfer phase only if accumulated version of x is requested
    if ( doSecondComm)
        InitComm( 2, *x_acc, Addr(req)+num_sr_1, Addr(req)+num_sr_1+sendListPhase2_.size(), xBuf_, 1002);

    // While communication accumulated values, do product on distributed elements
    result= KahanInnerProd( *x_acc, *x_acc, OwnerDistrIndex.begin(), OwnerDistrIndex.end(), sum1);

    // Before touching the memory of y_acc and returning y, wait
    // until send and received are done.
    ProcCL::WaitAll( sendListPhase1_.size(), Addr(req));
    if ( doSecondComm){
        ProcCL::WaitAll( num_sr_2, Addr(req)+num_sr_1);
        DoAllAssigning( *x_acc, xBuf_);
    }

    // Free memory, if allocated here
    if ( x_acc_created){
        delete x_acc; x_acc=0;
    }

    return result;
}

double ExchangeCL::LocalDotTwoAccumulations(
    const VectorCL& x, const VectorCL& y,
    VectorCL* x_acc, VectorCL* y_acc) const
/** Both vectors \a x and \a y are given in distributed form.
    \return inner product of x and y without reducing this value over
        all processes.
*/
{
    double result= std::numeric_limits<double>::max();
    bool x_acc_created=false, y_acc_created=false;

    // Check if memory for the accumulated form of x and/or y is provided
    if ( x_acc){
        *x_acc= x;
    }
    else{
        x_acc_created= true;
        x_acc= new VectorCL(x);
    }

    if ( y_acc){
        *y_acc= y;
    }
    else{
        y_acc_created= true;
        y_acc= new VectorCL(y);
    }

    // number of send and receive for phase 1 and 2
    const size_t num_sr_1=
        sendListPhase1_.size() + recvListPhase1_.size();
    const size_t num_sr_2=
         sendListPhase2_.size() + recvListPhase2_.size();
    // Request list for first and second phase
    RequestListT reqX(num_sr_1+num_sr_2);
    RequestListT reqY(num_sr_1+num_sr_2);

    // Initiate the first communication phase for x and y
    InitComm( 1, x, Addr(reqX), Addr(reqX)+sendListPhase1_.size(), xBuf_, 1001);
    InitComm( 1, y, Addr(reqY), Addr(reqY)+sendListPhase1_.size(), yBuf_, 1002);

    // While communicating, do product on local elements
    const double sum1=KahanInnerProd( x, y, LocalIndex.begin(), LocalIndex.end(), double());

    // Do accumulation on owners, check before, if the receive has been performed.
    // Since *x_acc and *y_acc is not used as sendbuffer, we do not have to wait
    // for completing the sending
    ProcCL::WaitAll( recvListPhase1_.size(), Addr(reqX)+sendListPhase1_.size());
    DoAllAccumulations( *x_acc, xBuf_);
    ProcCL::WaitAll( recvListPhase1_.size(), Addr(reqY)+sendListPhase1_.size());
    DoAllAccumulations( *y_acc, yBuf_);

    // Init Second transfer phase
    if ( !x_acc_created && viaowner_)
        InitComm( 2, *x_acc, Addr(reqX)+num_sr_1, Addr(reqX)+num_sr_1+sendListPhase2_.size(), xBuf_, 1003);
    if ( !y_acc_created && viaowner_)
        InitComm( 2, *y_acc, Addr(reqY)+num_sr_1, Addr(reqY)+num_sr_1+sendListPhase2_.size(), xBuf_, 1004);

    // While communication accumulated values, do product on distributed elements
    result= KahanInnerProd( *x_acc, *y_acc, OwnerDistrIndex.begin(), OwnerDistrIndex.end(), sum1);

    // Before touching the memory of x_acc (y_acc) and returning x (y), wait
    // until send and received are done.
    ProcCL::WaitAll( sendListPhase1_.size(), Addr(reqX));
    ProcCL::WaitAll( sendListPhase1_.size(), Addr(reqY));
    if ( !x_acc_created && viaowner_){
        ProcCL::WaitAll( num_sr_2, Addr(reqX)+num_sr_1);
        DoAllAssigning( *x_acc, xBuf_);
    }
    if ( !y_acc_created && viaowner_){
        ProcCL::WaitAll( num_sr_2, Addr(reqY)+num_sr_1);
        DoAllAssigning( *y_acc, yBuf_);
    }
    // Free memory of x_acc and y_acc if allocated here
    if ( x_acc_created){
        delete x_acc; x_acc=0;
    }
    if ( y_acc_created){
        delete y_acc; y_acc=0;
    }

    return result;
}

void ExchangeCL::Accumulate( VectorCL& v, const Ulint offset) const
/** Transform the vector \a v given in distributed form into the accumulated form.
    \param v Vector given in distributed form
    \param offset This offset is used to access elements of v. In most cases,
        this parameter should be set to 0.
*/
{
    // number of send and receive for phase 1 and 2
    const size_t num_sr_1=
        sendListPhase1_.size() + recvListPhase1_.size();
    const size_t num_sr_2=
         sendListPhase2_.size() + recvListPhase2_.size();
    RequestListT req(num_sr_1+num_sr_2);

    // first communication phase
    InitComm( 1, v, Addr(req), Addr(req)+sendListPhase1_.size(), xBuf_, 1001, offset);

    // Wait until all data are received
    ProcCL::WaitAll( num_sr_1, Addr(req));

    // Accumulate values by owners
    DoAllAccumulations(v, xBuf_, offset);

    // second communication phase
    if ( viaowner_){
        InitComm( 2, v, Addr(req)+num_sr_1, Addr(req)+num_sr_1+sendListPhase2_.size(), xBuf_, 1002, offset);
        ProcCL::WaitAll( num_sr_2, Addr(req)+num_sr_1);
        DoAllAssigning( v, xBuf_, offset);
    }
}

void ExchangeCL::Accumulate( std::vector<VectorCL>& v) const
/// \todo test me
{
    // number of send and receive for phase 1 and 2
    const size_t num_sr_1=
        sendListPhase1_.size() + recvListPhase1_.size();
    const size_t num_sr_2=
         sendListPhase2_.size() + recvListPhase2_.size();

    std::vector<BufferListT> recvbuf(v.size(), xBuf_);
    std::vector<RequestListT> req( v.size(), RequestListT( num_sr_1+num_sr_2));

    // init send and receive for the first phase
    for ( size_t i=0; i<v.size(); ++i){
        InitComm( 1, v[i], Addr(req[i]), Addr(req[i])+sendListPhase1_.size(), recvbuf[i], 1001+i);
    }

    for ( size_t i=0; i<v.size(); ++i){
        ProcCL::WaitAll( num_sr_1, Addr(req[i]));
    }

    for ( size_t i=0; i<v.size(); ++i){
        DoAllAccumulations( v[i], recvbuf[i]);
    }

    if ( viaowner_){
        for ( size_t i=0; i<v.size(); ++i){
            InitComm( 2, v[i], Addr(req[i])+num_sr_1, Addr(req[i])+num_sr_1+sendListPhase2_.size(), recvbuf[i], 2001+i);
        }

        for ( size_t i=0; i<v.size(); ++i){
            ProcCL::WaitAll( num_sr_2, Addr(req[i])+num_sr_1);
        }

        for ( size_t i=0; i<v.size(); ++i){
            DoAllAssigning( v[i], recvbuf[i]);
        }
    }
}

VectorCL ExchangeCL::GetAccumulate (const VectorCL& v) const
/** Get an accumulated copy of the distributed vector \a v. */
{
    VectorCL v_acc(v);
    Accumulate(v_acc);
    return v_acc;
}

std::vector<VectorCL> ExchangeCL::GetAccumulate( const std::vector<VectorCL>& v) const
/** Get an accumulated copy of the distributed vectors given in \a v. */
{
    std::vector<VectorCL> v_acc( v);
    Accumulate( v_acc);
    return v_acc;
}

double ExchangeCL::LocalDot(
    const VectorCL& x, bool isXacc,
    const VectorCL& y, bool isYacc,
    VectorCL* x_acc, VectorCL* y_acc) const
/** Do an inner product without the global reduction over all processes. That is,
    each process only stores a part of the inner product. For a description of the
    parameters, we refer to this member function double ExchangeCL::ParDot(...).
*/
{
    double result= std::numeric_limits<double>::max();
    if ( isXacc && isYacc){
        result= LocalDotNoAccumulation( x, y);
        if ( x_acc) *x_acc= x;
        if ( y_acc) *y_acc= y;
    }
    if ( isXacc && !isYacc){
        result= LocalDotOneAccumulation( x, y, y_acc);
        if ( x_acc) *x_acc= x;
    }
    if ( isYacc && !isXacc){
        result= LocalDotOneAccumulation( y, x, x_acc);
        if ( y_acc) *y_acc= y;
    }
    if ( !isYacc && !isXacc){
        result = LocalDotTwoAccumulations( x, y, x_acc, y_acc);
    }
    return result;
}

double ExchangeCL::ParDot(
    const VectorCL& x, bool isXacc,
    const VectorCL& y, bool isYacc,
    VectorCL* x_acc, VectorCL* y_acc) const
/** Compute the inner product of two vectors \a x and \a y and returns this value.
    The parameter \a isXacc indicates whether \a x is accumulated or distributed,
    i.e., of form (b) or (a). The parameter \a isYacc provides the same information
    for the vector \a y. <br>
    If pointers \a x_acc or \a y_acc are given, on output, the referenced vectors
    containing a copy of the accumulated form of \a x or \a y.
*/
{
    return ProcCL::GlobalSum( LocalDot( x, isXacc, y, isYacc, x_acc, y_acc));
}

double ExchangeCL::LocalNorm_sq( const VectorCL& x, bool isXacc, VectorCL* x_acc) const
/** This function computes the squared Euklidian norm of a vector \a x without
    performing the reduction via all processes. For more detailed information, we refer
    to double ExchangeCL::Norm(...). */
{
    double result= std::numeric_limits<double>::max();
    if ( isXacc){
        result= LocalDotNoAccumulation( x, x);
        if ( x_acc) *x_acc=x;
    }
    else{
        result= LocalNormSQAccumulation( x, x_acc);
    }
    return result;
}

double ExchangeCL::Norm_sq( const VectorCL& x, bool isXacc, VectorCL* x_acc) const
/** This function computes the squared Euklidian norm of a vector \a x. For
    more detailed information, we refer to double ExchangeCL::Norm(...). */
{
    return ProcCL::GlobalSum( LocalNorm_sq(x, isXacc, x_acc));
}

double ExchangeCL::Norm( const VectorCL& x, bool isXacc, VectorCL* x_acc) const
/** Determine the Euclidian norm of a vector \a x. The flag isXacc indicates, if the
    vector x is provided in accumulated (b) or distributed (a) form. If a pointer
    x_acc is given, on exits, the referenced vector contains the accumulated version
    of \a x. This function returns the value of the Euclidian norm of \a x.
*/
{
    return std::sqrt(Norm_sq(x, isXacc, x_acc));
}

void ExchangeCL::CreateList( const MultiGridCL& mg, IdxDescCL* rowidx, bool, bool)
/** Build the internal data structures to be able to provide the
    functionality of this class.
    \param mg The underlying multigrid
    \param rowidx the corresponding description of the finite element function
*/
{
    ExchangeBuilderCL builder( *this, mg, rowidx);
    builder.build();
}

void ExchangeCL::clear()
/** Free memory used by this class. */
{
    sendListPhase1_.clear();
    sendListPhase2_.clear();
    recvListPhase1_.clear();
    recvListPhase2_.clear();
    xBuf_.clear();
    yBuf_.clear();
    neighs_.clear();
    dofProcList_.clear();
    LocalIndex.clear();
    DistrIndex.clear();
    OwnerDistrIndex.clear();
}

/****************************************************************************
* E X C H A N G E  B L O C K  C L A S S                                     *
****************************************************************************/

void ExchangeBlockCL::AttachTo(const IdxDescCL& idx)
/** Beside attaching the index describer to the known indices, this functions fill
    the offsets array and link the receive buffer to allocated memory in the
    corresponding ExchangeCL s.
    \param[in] idx new index description
*/
{
    idxDesc_.push_back(&idx);
    Update();
}

void ExchangeBlockCL::Update()
/** Update the blockOffset_, where the sizes of the block vectors can be found. And,
    init the buffers for receiving. This function has to be called each time, the
    corresponding index describers has been changed
*/
{
    blockOffset_.resize( idxDesc_.size()+1);

    // fill block offsets
    blockOffset_[0]=0;
    for (size_t i=1; i<blockOffset_.size(); ++i){
        blockOffset_[i]= blockOffset_[i-1] + idxDesc_[i-1]->NumUnknowns();
    }

    // Set buffers
    xBuf_.resize( GetNumBlocks());
    yBuf_.resize( GetNumBlocks());
    for ( size_t i=0; i<GetNumBlocks(); ++i){
        xBuf_[i]= &(idxDesc_[i]->GetEx().xBuf_);
        yBuf_[i]= &(idxDesc_[i]->GetEx().yBuf_);
    }

    // Check if all ExchangeCLs have the same communication pattern
    bool comm= idxDesc_[0]->GetEx().CommViaOwner();
    for ( size_t i=1; i<GetNumBlocks(); ++i){
        if ( comm!=idxDesc_[i]->GetEx().CommViaOwner())
            throw DROPSErrCL("ExchangeBlockCL::Update: All ExchangeCLs must have \
                             the same communication pattern");
    }
}

void ExchangeBlockCL::InitComm(
    int Phase, const VectorCL& v,
    RequestListT& sendreq, RequestListT& recvreq,
    BufferListT& buf, int tag) const
/** For a detailed description, see ExchangeCL::InitComm. */
{
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        const ExchangeCL& ex= idxDesc_[j]->GetEx();
        // set iterators for sending
        ExchangeCL::SendListT::const_iterator sendit=
            (Phase==1) ? ex.sendListPhase1_.begin() : ex.sendListPhase2_.begin();
        ExchangeCL::SendListT::const_iterator sendend=
            (Phase==1) ? ex.sendListPhase1_.end() : ex.sendListPhase2_.end();
        size_t i=0;
        // init send
        for ( ; sendit!=sendend; ++sendit, ++i){
            sendreq[j][i]= sendit->Isend( v, tag+j, blockOffset_[j]);
        }

        // set iterators for receiving
        ExchangeCL::RecvListT::const_iterator recvit=
            (Phase==1) ? ex.recvListPhase1_.begin() : ex.recvListPhase2_.begin();
        ExchangeCL::RecvListT::const_iterator recvend=
            (Phase==1) ? ex.recvListPhase1_.end() : ex.recvListPhase2_.end();
        ExchangeCL::BufferListT::iterator bufit= buf[j]->begin();
        // init recv
        i=0;
        for ( ; recvit!=recvend; ++recvit, ++bufit, ++i){
            recvreq[j][i]= recvit->Irecv( tag+j, *bufit, 0);
        }

    }
}

void ExchangeBlockCL::DoAllAccumulations( VectorCL& v, const BufferListT& buf) const
/** For a detailed description, see ExchangeCL::DoAllAccumulations. */
{
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        const ExchangeCL& ex= idxDesc_[j]->GetEx();
        ex.DoAllAccumulations( v, *buf[j], blockOffset_[j]);
    }
}

void ExchangeBlockCL::DoAllAssigning( VectorCL& v, const BufferListT& buf) const
/** For a detailed description, see ExchangeCL::DoAllAssigning. */
{
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        const ExchangeCL& ex= idxDesc_[j]->GetEx();
        ex.DoAllAssigning( v, *buf[j], blockOffset_[j]);
    }
}

double ExchangeBlockCL::LocalDotNoAccumulation(const VectorCL& x, const VectorCL& y) const
{
    double sum=0.0;
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        const ExchangeCL& ex= idxDesc_[j]->GetEx();
        sum= KahanInnerProd( x, y, ex.LocalIndex.begin(), ex.LocalIndex.end(),
            sum, blockOffset_[j]);
        sum= KahanInnerProd( x, y, ex.OwnerDistrIndex.begin(), ex.OwnerDistrIndex.end(),
            sum, blockOffset_[j]);
    }
    return sum;
}

double ExchangeBlockCL::LocalDotOneAccumulation(
    const VectorCL& x, const VectorCL& y, VectorCL* y_acc) const
/** The vector x is given in accumulated form and y is given in distributed form.
    For a detailed description, see ExchangeCL::LocalDotOneAccumulation.
*/
{
    double result= 0;
    bool y_acc_created=false;

    // Check if memory for the accumulated form of y is provided
    if ( y_acc){
        *y_acc= y;
    }
    else{
        y_acc_created= true;
        y_acc= new VectorCL(y);
    }
    const bool doSecondComm= !y_acc_created && idxDesc_[0]->GetEx().CommViaOwner();

    RequestListT sendreq1( GetNumBlocks()), sendreq2( GetNumBlocks()), recvreq( GetNumBlocks());
    // Send distributed entries to the owners
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        const ExchangeCL& ex= idxDesc_[j]->GetEx();
        sendreq1[j].resize( ex.sendListPhase1_.size());
        recvreq[j].resize( ex.recvListPhase1_.size());
        ex.InitComm( 1, y, Addr(sendreq1[j]), Addr(recvreq[j]), *yBuf_[j],
            1001+j, blockOffset_[j]);
    }

    for ( size_t j=0; j<GetNumBlocks(); ++j){
        const ExchangeCL& ex= idxDesc_[j]->GetEx();
        result= KahanInnerProd( x, y, ex.LocalIndex.begin(),
            ex.LocalIndex.end(), result, blockOffset_[j]);
    }

    // Do accumulation on owners, therefore, first, wait until all
    // messages are received by the owner
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        ProcCL::WaitAll( recvreq[j]);
    }
    DoAllAccumulations( *y_acc, yBuf_);

    // Init Second transfer phase only if accumulated version of y is requested.
    // Here, the buffer for y can be re-used.
    if ( doSecondComm){
        for ( size_t j=0; j<GetNumBlocks(); ++j){
            const ExchangeCL& ex= idxDesc_[j]->GetEx();
            sendreq2[j].resize( ex.sendListPhase2_.size());
            recvreq[j].resize( ex.recvListPhase2_.size());
            ex.InitComm( 2, *y_acc, Addr(sendreq2[j]), Addr(recvreq[j]), *yBuf_[j],
                2001+j, blockOffset_[j]);
        }
    }

    // While communication accumulated values, do product on distributed elements
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        const ExchangeCL& ex= idxDesc_[j]->GetEx();
        result= KahanInnerProd( x, *y_acc, ex.OwnerDistrIndex.begin(),
            ex.OwnerDistrIndex.end(), result, blockOffset_[j]);
    }

    // Before touching the memory of y_acc and returning y, wait
    // until send operation and, eventually, receive operation are completed.
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        ProcCL::WaitAll( sendreq1[j]);
    }
    if ( doSecondComm){
        for ( size_t j=0; j<GetNumBlocks(); ++j){
            ProcCL::WaitAll( sendreq2[j]);
            ProcCL::WaitAll( recvreq[j]);
        }
        DoAllAssigning( *y_acc, yBuf_);
    }

    // Free memory, if allocated here
    if ( y_acc_created){
        delete y_acc; y_acc=0;
    }

    return result;
}

double ExchangeBlockCL::LocalNormSQAccumulation( const VectorCL& x, VectorCL* x_acc) const
/** For a detailed description, see ExchangeCL::LocalNormSQAccumulation. */
{
    double result= 0;
    bool x_acc_created=false;

    // Check if memory for the accumulated form of y is provided
    if ( x_acc){
        *x_acc= x;
    }
    else{
        x_acc_created= true;
        x_acc= new VectorCL(x);
    }
    const bool doSecondComm= !x_acc_created && idxDesc_[0]->GetEx().CommViaOwner();

    RequestListT sendreq1( GetNumBlocks()), sendreq2( GetNumBlocks()), recvreq( GetNumBlocks());
    // Send distributed entries to the owners
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        const ExchangeCL& ex= idxDesc_[j]->GetEx();
        sendreq1[j].resize( ex.sendListPhase1_.size());
        recvreq[j].resize( ex.recvListPhase1_.size());
        ex.InitComm( 1, x, Addr(sendreq1[j]), Addr(recvreq[j]), *xBuf_[j],
            1001+j, blockOffset_[j]);
    }

    // While communicating, do product on local elements
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        const ExchangeCL& ex= idxDesc_[j]->GetEx();
        result= KahanInnerProd( x, x, ex.LocalIndex.begin(),
            ex.LocalIndex.end(), result, blockOffset_[j]);
    }

    // Do accumulation on owners
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        ProcCL::WaitAll( recvreq[j]);
    }
    DoAllAccumulations( *x_acc, xBuf_);

    // Init Second transfer phase only if accumulated version of x is requested
    if ( doSecondComm){
        for ( size_t j=0; j<GetNumBlocks(); ++j){
            const ExchangeCL& ex= idxDesc_[j]->GetEx();
            sendreq2[j].resize( ex.sendListPhase2_.size());
            recvreq[j].resize( ex.recvListPhase2_.size());
            ex.InitComm( 2, *x_acc, Addr(sendreq2[j]), Addr(recvreq[j]), *xBuf_[j],
                2001+j, blockOffset_[j]);
        }
    }

    // While communication accumulated values, do product on distributed elements
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        const ExchangeCL& ex= idxDesc_[j]->GetEx();
        result= KahanInnerProd( *x_acc, *x_acc, ex.OwnerDistrIndex.begin(),
            ex.OwnerDistrIndex.end(), result, blockOffset_[j]);
    }

    // Before touching the memory of x_acc, wait
    // until send and received are done.
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        ProcCL::WaitAll( sendreq1[j]);
    }

    if ( doSecondComm){
        for ( size_t j=0; j<GetNumBlocks(); ++j){
            ProcCL::WaitAll( sendreq2[j]);
            ProcCL::WaitAll( recvreq[j]);
        }
        DoAllAssigning( *x_acc, xBuf_);
    }

    // Free memory, if allocated here
    if ( x_acc_created){
        delete x_acc; x_acc=0;
    }

    return result;
}

double ExchangeBlockCL::LocalDotTwoAccumulations(
    const VectorCL& x, const VectorCL& y, VectorCL* x_acc, VectorCL* y_acc) const
/** For a detailed description, see ExchangeCL::LocalDotTwoAccumulations. */
{
    double result= 0;
    bool x_acc_created=false, y_acc_created=false;

    // Check if memory for the accumulated form of x and/or y is provided
    if ( x_acc){
        *x_acc= x;
    }
    else{
        x_acc_created= true;
        x_acc= new VectorCL(x);
    }

    if ( y_acc){
        *y_acc= y;
    }
    else{
        y_acc_created= true;
        y_acc= new VectorCL(y);
    }

    RequestListT sendreqX1( GetNumBlocks()), sendreqX2( GetNumBlocks()), recvreqX( GetNumBlocks());
    RequestListT sendreqY1( GetNumBlocks()), sendreqY2( GetNumBlocks()), recvreqY( GetNumBlocks());

    // Initiate the first communication phase for x and y
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        const ExchangeCL& ex= idxDesc_[j]->GetEx();
        sendreqX1[j].resize( ex.sendListPhase1_.size());
        recvreqX[j].resize( ex.recvListPhase1_.size());
        ex.InitComm( 1, x, Addr(sendreqX1[j]), Addr(recvreqX[j]), *xBuf_[j],
            1001+j, blockOffset_[j]);
        sendreqY1[j].resize( ex.sendListPhase1_.size());
        recvreqY[j].resize( ex.recvListPhase1_.size());
        ex.InitComm( 1, y, Addr(sendreqY1[j]), Addr(recvreqY[j]), *yBuf_[j],
            1001+j+GetNumBlocks(), blockOffset_[j]);
    }

    // While communicating, do product on local elements
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        const ExchangeCL& ex= idxDesc_[j]->GetEx();
        result=KahanInnerProd( x, y, ex.LocalIndex.begin(),
            ex.LocalIndex.end(), result, blockOffset_[j]);
    }

    // Do accumulation on owners, check before, if the receive has been performed.
    // Since *x_acc and *y_acc is not used as sendbuffer, we do not have to wait
    // for completing the sending
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        ProcCL::WaitAll( recvreqX[j]);
    }
    DoAllAccumulations( *x_acc, xBuf_);
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        ProcCL::WaitAll( recvreqY[j]);
    }
    DoAllAccumulations( *y_acc, yBuf_);

    // Init Second transfer phase
    if ( !x_acc_created && idxDesc_[0]->GetEx().CommViaOwner()){
        for ( size_t j=0; j<GetNumBlocks(); ++j){
            const ExchangeCL& ex= idxDesc_[j]->GetEx();
            sendreqX2[j].resize( ex.sendListPhase2_.size());
            recvreqX[j].resize( ex.recvListPhase2_.size());
            ex.InitComm( 2, *x_acc, Addr(sendreqX2[j]), Addr(recvreqX[j]), *xBuf_[j],
                2001+j, blockOffset_[j]);
        }
    }
    if ( !y_acc_created && idxDesc_[0]->GetEx().CommViaOwner()){
        for ( size_t j=0; j<GetNumBlocks(); ++j){
            const ExchangeCL& ex= idxDesc_[j]->GetEx();
            sendreqY2[j].resize( ex.sendListPhase2_.size());
            recvreqY[j].resize( ex.recvListPhase2_.size());
            ex.InitComm( 2, *y_acc, Addr(sendreqY2[j]), Addr(recvreqY[j]), *yBuf_[j],
                2001+j+GetNumBlocks(), blockOffset_[j]);
        }
    }

    // While communication accumulated values, do product on distributed elements
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        const ExchangeCL& ex= idxDesc_[j]->GetEx();
        result= KahanInnerProd( *x_acc, *y_acc, ex.OwnerDistrIndex.begin(),
            ex.OwnerDistrIndex.end(), result, blockOffset_[j]);
    }

    // Before touching the memory of x_acc (y_acc) and returning x (y), wait
    // until send and received are done.
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        ProcCL::WaitAll( sendreqX1[j]);
        ProcCL::WaitAll( sendreqY1[j]);
    }
    if ( !x_acc_created && idxDesc_[0]->GetEx().CommViaOwner()){
        for ( size_t j=0; j<GetNumBlocks(); ++j){
            ProcCL::WaitAll( sendreqX2[j]);
            ProcCL::WaitAll( recvreqX[j]);
        }
        DoAllAssigning( *x_acc, xBuf_);
    }
    if ( !y_acc_created && idxDesc_[0]->GetEx().CommViaOwner()){
        for ( size_t j=0; j<GetNumBlocks(); ++j){
            ProcCL::WaitAll( sendreqY2[j]);
            ProcCL::WaitAll( recvreqY[j]);
        }
        DoAllAssigning( *y_acc, yBuf_);
    }
    // Free memory of x_acc and y_acc if allocated here
    if ( x_acc_created){
        delete x_acc; x_acc=0;
    }
    if ( y_acc_created){
        delete y_acc; y_acc=0;
    }

    return result;
}

void ExchangeBlockCL::Accumulate( VectorCL& v) const
/** For a detailed description, see ExchangeCL::Accumulate. */
{
    RequestListT sendreq( GetNumBlocks()), recvreq( GetNumBlocks());
    // Send distributed entries to the owners
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        const ExchangeCL& ex= idxDesc_[j]->GetEx();
        sendreq[j].resize( ex.sendListPhase1_.size());
        recvreq[j].resize( ex.recvListPhase1_.size());
        ex.InitComm( 1, v, Addr(sendreq[j]), Addr(recvreq[j]), *xBuf_[j], 1001+j, blockOffset_[j]);
    }
    // Wait for all receive and send operations are completed
    for ( size_t j=0; j<GetNumBlocks(); ++j){
        ProcCL::WaitAll( sendreq[j]);
        ProcCL::WaitAll( recvreq[j]);
    }
    // Do all accumulations
    DoAllAccumulations( v, xBuf_);

    if ( idxDesc_[0]->GetEx().CommViaOwner()){
        // Send back to copies
        for ( size_t j=0; j<GetNumBlocks(); ++j){
            const ExchangeCL& ex= idxDesc_[j]->GetEx();
            sendreq[j].resize( ex.sendListPhase2_.size());
            recvreq[j].resize( ex.recvListPhase2_.size());
            ex.InitComm( 2, v, Addr(sendreq[j]), Addr(recvreq[j]), *xBuf_[j], 2001+j, blockOffset_[j]);
        }

        // Wait for all sends receives to be completed
        // Wait for all receive and send operations are completed
        for ( size_t j=0; j<GetNumBlocks(); ++j){
            ProcCL::WaitAll( sendreq[j]);
            ProcCL::WaitAll( recvreq[j]);
        }

        // Do all assignings
        DoAllAssigning( v, xBuf_);
    }
}

void ExchangeBlockCL::Accumulate( std::vector<VectorCL>& v) const
{
    for (size_t i = 0; i< v.size(); ++i)
        Accumulate(v[i]);
}

VectorCL ExchangeBlockCL::GetAccumulate (const VectorCL& v) const
/** For a detailed description, see ExchangeCL::GetAccumulate. */
{
    VectorCL v_acc(v);
    Accumulate( v_acc);
    return v_acc;
}

std::vector<VectorCL> ExchangeBlockCL::GetAccumulate( const std::vector<VectorCL>& v) const
/** Get an accumulated copy of the distributed vectors given in \a v. */
{
    std::vector<VectorCL> v_acc( v);
    Accumulate( v_acc);
    return v_acc;
}

double ExchangeBlockCL::LocalDot(
    const VectorCL& x, bool isXacc,
    const VectorCL& y, bool isYacc,
    VectorCL* x_acc, VectorCL* y_acc) const
/** Do an inner product without the global reduction over all processes. That is,
    each process only stores a part of the inner product. For a description of the
    parameters, we refer to \a double ExchangeCL::ParDot(...).
*/
{
    double result= std::numeric_limits<double>::max();
    if ( isXacc && isYacc){
        result= LocalDotNoAccumulation( x, y);
        if ( x_acc) *x_acc= x;
        if ( y_acc) *y_acc= y;
    }
    if ( isXacc && !isYacc){
        result= LocalDotOneAccumulation( x, y, y_acc);
        if ( x_acc) *x_acc= x;
    }
    if ( isYacc && !isXacc){
        result= LocalDotOneAccumulation( y, x, x_acc);
        if ( y_acc) *y_acc= y;
    }
    if ( !isYacc && !isXacc){
        result = LocalDotTwoAccumulations( x, y, x_acc, y_acc);
    }
    return result;
}

double ExchangeBlockCL::ParDot(
    const VectorCL& x, bool isXacc,
    const VectorCL& y, bool isYacc,
    VectorCL* x_acc, VectorCL* y_acc) const
/** Compute the inner product of two vectors \a x and \a y and returns this value.
    The parameter \a isXacc indicates whether \a x is accumulated or distributed,
    i.e., of form (b) or (a). The parameter \a isYacc provides the same information
    for the vector \a y.
    If pointers \a x_acc or \a y_acc are given, on output, the referenced vectors
    containing a copy of the accumulated form of \a x or \a y.
*/
{
    return ProcCL::GlobalSum( LocalDot( x, isXacc, y, isYacc, x_acc, y_acc));
}

double ExchangeBlockCL::LocalNorm_sq( const VectorCL& x, bool isXacc, VectorCL* x_acc) const
/** This function computes the squared Euklidian norm of a vector \a x without
    performing the reduction via all processes. For more detailed information, we refer
    to double ExchangeCL::Norm(...). */
{
    double result= std::numeric_limits<double>::max();
    if ( isXacc){
        result= LocalDotNoAccumulation( x, x);
        if ( x_acc) *x_acc=x;
    }
    else{
        result= LocalNormSQAccumulation( x, x_acc);
    }
    return result;
}

double ExchangeBlockCL::Norm_sq( const VectorCL& x, bool isXacc, VectorCL* x_acc) const
/** This function computes the squared Euklidian norm of a vector \a x. For
    more detailed information, we refer to double ExchangeCL::Norm(...). */
{
    return ProcCL::GlobalSum( LocalNorm_sq(x, isXacc, x_acc));
}

double ExchangeBlockCL::Norm( const VectorCL& x, bool isXacc, VectorCL* x_acc) const
/** Determine the Euclidian norm of a vector \a x. The flag isXacc indicates, if the
    vector x is provided in accumulated (b) or distributed (a) form. If a pointer
    x_acc is given, on exits, the referenced vector contains the accumulated version
    of \a x. This function returns the value of the Euclidian norm of \a x.
*/
{
    return std::sqrt(Norm_sq(x, isXacc, x_acc));
}


/******************************************************************************
* E X C H A N G E  B U I L D E R  C L                                         *
******************************************************************************/

// STATIC MEMBER INITIALIZATION OF EXCHANGE BUILDER CL
//----------------------------------------------------
int ExchangeBuilderCL::NoInt_= std::numeric_limits<int>::max();

ExchangeBuilderCL::ExchangeBuilderCL(
    ExchangeCL& ex, const MultiGridCL& mg, IdxDescCL *RowIdx)
    : ex_(ex), mg_(mg), rowidx_(*RowIdx), interf_(0)
{
    // init the interface
    DiST::PrioListT priolist;
    priolist.push_back( PrioMaster);
    DiST::InterfaceCL::DimListT dimlist;
    if ( rowidx_.NumUnknownsVertex())
        dimlist.push_back( DiST::GetDim<VertexCL>());
    if ( rowidx_.NumUnknownsEdge())
        dimlist.push_back( DiST::GetDim<EdgeCL>());
    if ( rowidx_.NumUnknownsFace())
        dimlist.push_back( DiST::GetDim<FaceCL>());
    if ( rowidx_.NumUnknownsTetra())
        dimlist.push_back( DiST::GetDim<TetraCL>());
    DiST::LevelListCL lvls( rowidx_.TriangLevel()); // levels 0,..,TriangLevel

    interf_ = new DiST::InterfaceCL( lvls, priolist, priolist, dimlist, true);
}

void ExchangeBuilderCL::clearEx()
{
    ex_.clear();
}

void ExchangeBuilderCL::BuildIndexLists()
/** This function fills the lists LocalIndex, DistrIndex, OwnerDistrIndex,
    and dofProcList_ of the ExchangeCL.
*/
{
    ex_.dofProcList_.resize( rowidx_.NumUnknowns());

    // Communicate dof positions via process boundaries, after the
    // call, OwnerDistrIndex and dofProcList_ is set up.
    HandlerDOIndexCL handlerIndex( ex_, rowidx_);
    interf_->PerformInterfaceComm( handlerIndex);

    // sort the list of OwnerDistrIndex for better memory access pattern
    std::sort( ex_.OwnerDistrIndex.begin(), ex_.OwnerDistrIndex.end());

    // determine the indices of local dof
    ExchangeCL::IdxVecT& LocalIndex= ex_.LocalIndex;
    ExchangeCL::IdxVecT& DistrIndex= ex_.DistrIndex;
    size_t numDistIdx=0;
    ExchangeCL::DOFProcListT::const_iterator it;
    for ( it=ex_.dofProcList_.begin(); it!= ex_.dofProcList_.end(); ++it){
        if ( it->size()>=1)
            ++numDistIdx;
    }
    // reserve memory.
    LocalIndex.reserve( numDistIdx);
    DistrIndex.reserve( rowidx_.NumUnknowns()- numDistIdx);

    // fill LocalIndex and DistrIndex
    for ( IdxT i=0; i<ex_.dofProcList_.size(); ++i){
        if ( ex_.dofProcList_[i].empty())
            LocalIndex.push_back(i);
        else
            DistrIndex.push_back(i);
    }

    // determine neighbors
    for ( size_t i=0; i<DistrIndex.size(); ++i){
        for ( ExchangeCL::DOFInfoList_const_iterator it=ex_.GetProcListBegin( DistrIndex[i]); it!=ex_.GetProcListEnd( DistrIndex[i]); ++it)
            ex_.neighs_.insert( it->first);
    }
}

void ExchangeBuilderCL::buildViaOwner()
{
    // clear all information previously determined
    clearEx();

    // Fill the lists sendListPhase1_ and sendListPhase2_, i.e., the dofs
    // that need to be sent in the first and second communication phase
    HandlerDOFtoOwnerCL handlerToOwner( rowidx_, mg_);
    interf_->InformOwners( handlerToOwner);
    handlerToOwner.buildSendStructures( ex_.sendListPhase1_);
    handlerToOwner.buildRecvStructures( ex_.recvListPhase1_);

    // let the copies know about the send position of all distributed dofs.
    // That is, determine information for the second communication phase
    HandlerDOFFromOwnerCL handlerFromOwner( rowidx_, mg_);
    interf_->InformCopies( handlerFromOwner);
    handlerFromOwner.buildSendStructures( ex_.sendListPhase2_);
    handlerFromOwner.buildRecvStructures( ex_.recvListPhase2_);

    // determine information about distributed dof
    BuildIndexLists();

    // allocate memory for receive buffers
    BuildRecvBuffer();
}

void ExchangeBuilderCL::BuildRecvBuffer()
{
    ex_.xBuf_.resize( std::max( ex_.recvListPhase2_.size(), ex_.recvListPhase1_.size()));
    ex_.yBuf_.resize( ex_.xBuf_.size());

    size_t i=0;
    ExchangeCL::RecvListT::const_iterator it= ex_.recvListPhase2_.begin();
    for ( ; it!=ex_.recvListPhase2_.end(); ++it, ++i){
        ex_.xBuf_[i].resize( it->sysnums_.size());
        ex_.yBuf_[i].resize( it->sysnums_.size());
    }
    i=0;
    it= ex_.recvListPhase1_.begin();
    for ( ; it!=ex_.recvListPhase1_.end(); ++it, ++i){
        ex_.xBuf_[i].resize( std::max( ex_.xBuf_[i].size(), it->sysnums_.size()));
        ex_.yBuf_[i].resize( std::max( ex_.yBuf_[i].size(), it->sysnums_.size()));
    }
}


void ExchangeBuilderCL::buildDirectComm()
{
    clearEx();

    // determine communication structure
    HandlerDOFDirectCommCL handerDOFDirect( rowidx_, mg_);
    interf_->PerformInterfaceComm( handerDOFDirect);
    handerDOFDirect.buildSendStructures( ex_.sendListPhase1_);
    handerDOFDirect.buildRecvStructures( ex_.recvListPhase1_);

    // determine information about distributed dof
    BuildIndexLists();

    // allocate memory for receive buffers
    BuildRecvBuffer();
}


/******************************************************************************
* H A N D L E R  D O F  E X C H A N G E   C L                                 *
******************************************************************************/

void ExchangeBuilderCL::HandlerDOFExchangeCL::collectDOF()
/** Collect all dof which must be send to the owning process and sort them afterwards.
    \todo Collect DOFs on Faces
*/
{
    if ( rowidx_.NumUnknownsVertex())
        DROPS_FOR_TRIANG_CONST_VERTEX( mg_, rowidx_.TriangLevel(), it)
            collectDOFonSimplex( *it);
    if ( rowidx_.NumUnknownsEdge())
        DROPS_FOR_TRIANG_CONST_EDGE( mg_, rowidx_.TriangLevel(), it)
            collectDOFonSimplex( *it);
//    if ( rowidx_.NumUnknownsFace())
//        DROPS_FOR_TRIANG_CONST_FACE( mg_, rowidx_.TriangLevel(), it)
//            this->collectDOF( *it);
    if ( rowidx_.NumUnknownsTetra())
        DROPS_FOR_TRIANG_CONST_TETRA( mg_, rowidx_.TriangLevel(), it)
            collectDOFonSimplex( *it);

    SendDOFListT::iterator it;
    for ( it=sendList_.begin(); it!=sendList_.end(); ++it)
        std::sort( it->second.begin(), it->second.end());
}

void ExchangeBuilderCL::HandlerDOFExchangeCL::buildSendStructures(
    ExchangeCL::SendListT& ex_sendlist)
/** Create the MPI datatype for all neighbors for sending. */
{
    // Create the data types
    SendDOFListT::const_iterator it;
    for ( it=sendList_.begin(); it!=sendList_.end(); ++it){
        const int toproc= it->first;
        if ( toproc!=ProcCL::MyRank()){
            ex_sendlist.push_back( SendNumDataCL<double>(toproc));
            const int count= static_cast<int>(it->second.size());
            const int bl   = static_cast<int>( rowidx_.NumUnknownsVertex());
            const int* ad  = Addr( it->second);
            ex_sendlist.back().CreateDataType( count, bl, ad);
        }
    }
}

void ExchangeBuilderCL::HandlerDOFExchangeCL::buildRecvStructures(
    ExchangeCL::RecvListT& ex_recvlist)
/** Store how to receive data from neighbors. */
{
    // Create the receive ordering
    RecvDofT::const_iterator mit;
    const Uint numUnk= rowidx_.NumUnknownsVertex();
    for ( mit=recvList_.begin(); mit!=recvList_.end(); ++mit){
        const int fromproc= mit->first;
        if ( fromproc!=ProcCL::MyRank()){
            std::cout << "In ... buildRecvStructures, first 10 positions are\n";
            int j=0;
            for (RecvDofT::mapped_type::const_iterator it=mit->second.begin(); j!=10; ++it, ++j){
                std::cout << it->first << ' ';
            }
            std::cout << std::endl;

            ex_recvlist.push_back( RecvNumDataCL<double>(fromproc));
            std::vector<IdxT>& sysnums= ex_recvlist.back().sysnums_;
            sysnums.reserve( numUnk*mit->second.size());
            RecvDofT::mapped_type::const_iterator recvposit= mit->second.begin();
            int i=0;    // for debugging
            for ( ; recvposit!=mit->second.end(); ++recvposit, ++i){
                const IdxT local_dof= recvposit->second;
                Assert( (int)numUnk*i==recvposit->first,
                    DROPSErrCL("ExchangeBuilderCL::buildSendRecv1CommPhase: Missing send position"),
                    DebugParallelNumC);
                for ( Uint j=0; j<numUnk; ++j)
                    sysnums.push_back( local_dof+j);
            }
        }
    }
}


/******************************************************************************
* H A N D L E R  D O F  T O  O W N E R   C L                                  *
******************************************************************************/

void ExchangeBuilderCL::HandlerDOFtoOwnerCL::collectDOFonSimplex( const DiST::TransferableCL& s)
/** Collect the dof that need to be send to the owner in the first communication phase,
    i.e., fill the list sendList_.
*/
{
    const Uint idx= rowidx_.GetIdx();
    // Only collect data on simplices, where the dof are distributed among processes
    if ( !s.IsDistributed( PrioMaster))
        return;
    // Only collect data on simplices, where a dof is given
    if ( s.Unknowns.Exist() && s.Unknowns.Exist(idx)){

        // local dof information
        const IdxT dof= s.Unknowns(idx);
        const bool isExtended= (rowidx_.IsExtended() && rowidx_.GetXidx()[dof]!=NoIdx);
        const IdxT extdof= (isExtended) ? rowidx_.GetXidx()[dof] : NoIdx;

        // dof that need to be sent to the owner
        sendList_[ s.GetOwner()].push_back(static_cast<int>(dof));
        if ( isExtended)
            sendList_[ s.GetOwner()].push_back(static_cast<int>(extdof));
    }
}

bool ExchangeBuilderCL::HandlerDOFtoOwnerCL::Gather( DiST::TransferableCL& t,
    DiST::Helper::SendStreamCL& send)
/** If unknowns on the simplex \a t exists, then put
      (1) my rank
      (2) and the send position of the dof
      (3) the send position of the extended dof
      into the send buffer. Return true.
    Else,
      return false;
*/
{
    const Uint idx= rowidx_.GetIdx();
    if ( t.Unknowns.Exist() && t.Unknowns.Exist( idx) && t.IsDistributed(PrioMaster)){
        send << ProcCL::MyRank();

        const IdxT dof= t.Unknowns(idx);
        const Uint numUnk= rowidx_.NumUnknownsSimplex( t);
        const int firstPos= getSendPos( static_cast<int>(dof), t.GetOwner())*numUnk;
        send << firstPos;
        const bool isExtended= rowidx_.IsExtended() && rowidx_.GetXidx()[dof]!=NoIdx;
        if ( isExtended){
            const IdxT exdof= rowidx_.GetXidx()[dof];
            const int exfirstPos= getSendPos( static_cast<int>(exdof), t.GetOwner())*numUnk;
            send << exfirstPos;
        }
        else{
            send << NoInt_;
        }
        return true;
    }
    else {
        return false;
    }
}

bool ExchangeBuilderCL::HandlerDOFtoOwnerCL::Scatter( DiST::TransferableCL& t,
    const size_t& numData, DiST::Helper::RecvStreamCL& recv)
/** Owners have to create a mapping from the sendposition to the local dofs.*/
{
    const Uint idx= rowidx_.GetIdx();

    Assert (t.Unknowns.Exist() && t.Unknowns.Exist(idx),
        DROPSErrCL("ExchangeBuilderCL::ScatterDOFtoOwner: Unknowns does not exist"),
        DebugParallelNumC);
    Assert( t.AmIOwner(),
        DROPSErrCL("ExchangeBuilderCL::ScatterDOFtoOwner: Only intended to be called by owners"),
        DebugParallelNumC);

    // local dofs
    const IdxT dof= t.Unknowns(idx);
    const bool isExtended= (rowidx_.IsExtended() && rowidx_.GetXidx()[dof]!=NoIdx);
    const IdxT exdof=  isExtended ? rowidx_.GetXidx()[dof] : NoIdx;

    // temporaries for reading the receive stream
    int sender=-1;
    int sendpos=NoInt_;

    for ( size_t data=0; data<numData; ++data){
        recv >> sender;
        // receive non-extended dof
        recv >> sendpos;

        recvList_[sender][sendpos]= dof;
        // receive extended dof
        recv >> sendpos;
        if ( sendpos!=NoInt_){
            Assert( isExtended,
                DROPSErrCL("ExchangeBuilderCL::ScatterDOFtoOwner: DOF on owner is not extended"),
                DebugParallelNumC);
            recvList_[sender][sendpos]= exdof;
        }

    }
    return true;
}


/******************************************************************************
* H A N D L E R  D O F  F R O M  O W N E R   C L                              *
******************************************************************************/

void ExchangeBuilderCL::HandlerDOFFromOwnerCL::collectDOFonSimplex( const DiST::TransferableCL& s)
/** Collect the dof that need to be send from the owner in the second communication phase,
    i.e., fill the list sendList_.
*/
{
    const Uint idx= rowidx_.GetIdx();
    // Only collect data on simplices where DOF are distributed and I am owner of
    if ( !s.IsDistributed(PrioMaster) || !s.AmIOwner())
        return;
    // Only collect data on simplices, where a dof is given
    if ( s.Unknowns.Exist() && s.Unknowns.Exist(idx)){
        // local dof information
        const IdxT dof= s.Unknowns(idx);
        const bool isExtended= (rowidx_.IsExtended() && rowidx_.GetXidx()[dof]!=NoIdx);
        const IdxT extdof= (isExtended) ? rowidx_.GetXidx()[dof] : NoIdx;

        DiST::TransferableCL::ProcList_const_iterator it;
        for ( it=s.GetProcListBegin(); it!=s.GetProcListEnd(); ++it){
            const int toproc= it->proc;
            if ( it->prio==PrioMaster){   // Unknowns exist on toproc
                sendList_[toproc].push_back(static_cast<int>(dof));
                if ( isExtended){
                    sendList_[toproc].push_back(static_cast<int>(extdof));
                }
            }
        }
    }
}

bool ExchangeBuilderCL::HandlerDOFFromOwnerCL::Gather( DiST::TransferableCL& t,
    DiST::Helper::SendStreamCL& send)
/** The owner informs all copies about the sendposition of the dof located at the simplex.
    Therefore, put the rank of the receiving process, sendposition of the dof and extended
    dof into the send buffer. Finalize the stream by NoInt_
*/
{
    Assert( t.AmIOwner(),
        DROPSErrCL("ExchangeBuilderCL::GatherDOFfromOwner: Only intended to be called by owners"),
        DebugParallelNumC);

    const Uint idx= rowidx_.GetIdx();
    if ( t.Unknowns.Exist() && t.Unknowns.Exist( idx) && t.IsDistributed( PrioMaster)){
        // local information
        const IdxT dof= t.Unknowns(idx);
        const bool isExtended= (rowidx_.IsExtended() && rowidx_.GetXidx()[dof]!=NoIdx);
        const IdxT extdof= (isExtended) ? rowidx_.GetXidx()[dof] : NoIdx;
        const Uint numUnk= rowidx_.NumUnknownsSimplex( t);

        // Inform all (master) neighbors
        DiST::TransferableCL::ProcList_const_iterator it;
        for ( it=t.GetProcListBegin(); it!=t.GetProcListEnd(); ++it){
            const int toproc= it->proc;
            send << toproc;
            if ( it->prio==PrioMaster){
                const int firstPos=getSendPos( static_cast<int>(dof), toproc)*numUnk;
                send << firstPos;
                if ( isExtended){
                    const int extFirstPos=
                        getSendPos( static_cast<int>(extdof), toproc)*numUnk;
                    send << extFirstPos;
                }
                else{
                    send << NoInt_;         // dof is not extended
                }
            }
            else{
                send << NoInt_ << NoInt_;   // the copy does not own any dofs
            }
        }
        send << NoInt_;                     // finalize stream
        return true;
    }
    else{
        return false;
    }
}

bool ExchangeBuilderCL::HandlerDOFFromOwnerCL::Scatter( DiST::TransferableCL& t,
    __UNUSED__ const size_t& numData, DiST::Helper::RecvStreamCL& recv)
{
    // local dof
    const Uint idx= rowidx_.GetIdx();
    const IdxT dof= t.Unknowns(idx);
    const bool isExtended= (rowidx_.IsExtended() && rowidx_.GetXidx()[dof]!=NoIdx);
    const IdxT extdof= (isExtended) ? rowidx_.GetXidx()[dof] : NoIdx;

    // Check for wrong calls
    Assert( t.Unknowns.Exist() && t.Unknowns.Exist( idx),
        DROPSErrCL("ExchangeBuilderCL::ScatterDOFfromOwner: No dof available."),
        DebugParallelNumC);
    Assert( numData==1,
        DROPSErrCL("ExchangeBuilderCL::ScatterDOFfromOwner: More than one owner gathered data"),
        DebugParallelNumC);

    // temporaries for receiving
    int receiver= -1,    dummyreceiver=-1,
        sendpos_dof= -1, dummy1= -1,
        sendpos_ext= -1, dummy2= -1;

    // definitively, we have to read the complete receive stream
    recv >> dummyreceiver;
    while (dummyreceiver!=NoInt_){
        recv >> dummy1 >> dummy2;
        if ( dummyreceiver==ProcCL::MyRank()){
            Assert(receiver==-1,
                DROPSErrCL("ExchangeBuilderCL::ScatterDOFfromOwner: Received multiple information. I am confused."),
                DebugParallelNumC);
            receiver= dummyreceiver;
            sendpos_dof= dummy1;
            sendpos_ext= dummy2;
        }
        recv >> dummyreceiver;
    }
    // Check, if we have received valid data
    Assert( sendpos_dof>=0 && sendpos_ext>=0,
        DROPSErrCL("ExchangeBuilderCL::ScatterDOFfromOwner: No data received"),
        DebugParallelNumC);

    // Remember the position of this dof in a send operation from the owner (phase II).
    recvList_[t.GetOwner()][sendpos_dof]= dof;
    if ( isExtended){
        Assert( sendpos_ext!=NoInt_,
            DROPSErrCL("ExchangeBuilderCL::ScatterDOFfromOwner: Owner's dof is not extended"), DebugParallelNumC);
        recvList_[t.GetOwner()][sendpos_ext]= extdof;
    }
    return true;
}


/******************************************************************************
* H A N D L E R  D O F  D I R E C T  C O M M  C L                             *
******************************************************************************/

void ExchangeBuilderCL::HandlerDOFDirectCommCL::collectDOFonSimplex( const DiST::TransferableCL& s)
/** Collect the dof that need to be send to processes which store a copy of \a s,
    i.e., fill the list sendList_.
*/
{
    const Uint idx= rowidx_.GetIdx();
    // Only collect data on distributed simplices
    if ( !s.IsDistributed( PrioMaster))
        return;
    // Only collect data on simplices, where a dof is given
    if ( s.Unknowns.Exist() && s.Unknowns.Exist(idx)){

        // local dof information
        const IdxT dof= s.Unknowns(idx);
        const bool isExtended= (rowidx_.IsExtended() && rowidx_.GetXidx()[dof]!=NoIdx);
        const IdxT extdof= (isExtended) ? rowidx_.GetXidx()[dof] : NoIdx;

        // dof that need to be sent to copies (if there is a master copy)
        DiST::TransferableCL::ProcList_const_iterator it= s.GetProcListBegin();
        for ( ; it!=s.GetProcListEnd(); ++it){
            if ( it->prio== PrioMaster){
                sendList_[ it->proc].push_back(static_cast<int>(dof));
                if ( isExtended)
                    sendList_[ it->proc].push_back(static_cast<int>(extdof));
            }
        }
    }
}

bool ExchangeBuilderCL::HandlerDOFDirectCommCL::Gather( DiST::TransferableCL& t,
    DiST::Helper::SendStreamCL& send)
/** The process informs all copies about the sendposition of the dof located at the simplex.
    Therefore, put the rank of the receiving process, sendposition of the dof and extended
    dof into the send buffer. Finalize the stream by NoInt_
*/
{
    const Uint idx= rowidx_.GetIdx();
    if ( t.Unknowns.Exist() && t.Unknowns.Exist( idx)){
        // local information
        const IdxT dof= t.Unknowns(idx);
        const bool isExtended= (rowidx_.IsExtended() && rowidx_.GetXidx()[dof]!=NoIdx);
        const IdxT extdof= (isExtended) ? rowidx_.GetXidx()[dof] : NoIdx;
        const Uint numUnk= rowidx_.NumUnknownsSimplex( t);

        // First put my_rank on the stream
        send << ProcCL::MyRank();

        // Inform all (master) neighbors
        DiST::TransferableCL::ProcList_const_iterator it;
        for ( it=t.GetProcListBegin(); it!=t.GetProcListEnd(); ++it){
            const int toproc= it->proc;
            send << toproc;
            if ( it->prio==PrioMaster){
                const int firstPos=getSendPos( static_cast<int>(dof), toproc)*numUnk;
                send << firstPos;
                if ( isExtended){
                    const int extFirstPos=
                        getSendPos( static_cast<int>(extdof), toproc)*numUnk;
                    send << extFirstPos;
                }
                else{
                    send << NoInt_;         // dof is not extended
                }
            }
            else{
                send << NoInt_ << NoInt_;   // the copy does not own any dofs
            }
        }
        send << NoInt_;                     // finalize stream
        return true;
    }
    else{
        return false;
    }
}

bool ExchangeBuilderCL::HandlerDOFDirectCommCL::Scatter( DiST::TransferableCL& t,
    const size_t& numData, DiST::Helper::RecvStreamCL& recv)
{
    // local dof
    const Uint idx= rowidx_.GetIdx();
    const IdxT dof= t.Unknowns(idx);
    const bool isExtended= (rowidx_.IsExtended() && rowidx_.GetXidx()[dof]!=NoIdx);
    const IdxT extdof= (isExtended) ? rowidx_.GetXidx()[dof] : NoIdx;

    // Check for wrong calls
    Assert( t.Unknowns.Exist() && t.Unknowns.Exist( idx),
        DROPSErrCL("ExchangeBuilderCL::ScatterDOFfromOwner: No dof available."),
        DebugParallelNumC);

    // temporaries for receiving
    int receiver= -1,    dummyreceiver=-1,
        sendpos_dof= -1, dummy1= -1,
        sendpos_ext= -1, dummy2= -1,
        sender= -1;

    for ( size_t i=0; i<numData; ++i){
        // First get the sender of the dof information
        recv >> sender;

        // definitively, we have to read the complete receive stream
        recv >> dummyreceiver;
        while (dummyreceiver!=NoInt_){
            recv >> dummy1 >> dummy2;
            if ( dummyreceiver==ProcCL::MyRank()){
                Assert(receiver==-1,
                    DROPSErrCL("ExchangeBuilderCL::ScatterDOFfromOwner: Received multiple information. I am confused."),
                    DebugParallelNumC);
                receiver= dummyreceiver;
                sendpos_dof= dummy1;
                sendpos_ext= dummy2;
            }
            recv >> dummyreceiver;
        }

        // Check, if we have received valid data
        Assert( sendpos_dof>=0 && sendpos_ext>=0,
            DROPSErrCL("ExchangeBuilderCL::ScatterDOFfromOwner: No data received"),
            DebugParallelNumC);

        // Remember the position of this dof sent by another process.
        recvList_[sender][sendpos_dof]= dof;
        if ( isExtended){
            Assert( sendpos_ext!=NoInt_,
                DROPSErrCL("ExchangeBuilderCL::ScatterDOFfromOwner: Owner's dof is not extended"), DebugParallelNumC);
            recvList_[sender][sendpos_ext]= extdof;
        }

        // enable error checking
        receiver= -1; sendpos_dof=-1; sendpos_ext=-1;
    }

    return true;
}


/******************************************************************************
* H A N D L E R  D O F  I N D E X   C L                                       *
******************************************************************************/

bool ExchangeBuilderCL::HandlerDOIndexCL::Gather(
    DiST::TransferableCL& t, DiST::Helper::SendStreamCL& send)
/** For generating information about distributed dof, put my rank and the local (extended) dof
    into the buffer. Additionally, if this process is the owner of \a t, remember the dof as
    "owner dof."
*/
{
    const Uint idx= rowidx_.GetIdx();
    Assert( t.GetNumDist()>1,
        DROPSErrCL("ExchangeBuilderCL::GatherLocalDOF: Handler is called for a non-distributed simplex"),
        DebugParallelNumC);

    if ( t.Unknowns.Exist() && t.Unknowns.Exist(idx)){
        // local dof information
        const IdxT dof= t.Unknowns(idx);
        const bool isExtended= (rowidx_.IsExtended() && rowidx_.GetXidx()[dof]!=NoIdx);
        const IdxT extdof= (isExtended) ? rowidx_.GetXidx()[dof] : NoIdx;
        const Uint numUnk= rowidx_.NumUnknownsSimplex(t);

        send << ProcCL::MyRank() << dof << extdof;
        // Additionally, remember this dof as an "owner dof."
        if ( t.AmIOwner()){
            for ( Uint j=0; j<numUnk; ++j)
                ownerDistrIndex_.push_back( dof+j);
            if ( isExtended)
                for ( Uint j=0; j<numUnk; ++j)
                    ownerDistrIndex_.push_back( extdof+j);
        }
        return true;
    }
    return false;
}

bool ExchangeBuilderCL::HandlerDOIndexCL::Scatter(
    DiST::TransferableCL& t, const size_t& numData,
    DiST::Helper::RecvStreamCL& recv)
/** Store dof information on other processes in the dofProcList_. */
{
    const Uint idx= rowidx_.GetIdx();

    Assert( t.Unknowns.Exist() && t.Unknowns.Exist(idx),
        DROPSErrCL("ExchangeBuilderCL::ScatterLocalDOF: Received dof, but I do not store any dof"),
        DebugParallelNumC);

    // local information
    const IdxT dof= t.Unknowns(idx);
    const bool isExtended= (rowidx_.IsExtended() && rowidx_.GetXidx()[dof]!=NoIdx);
    const IdxT extdof= (isExtended) ? rowidx_.GetXidx()[dof] : NoIdx;
    const Uint numUnk= rowidx_.NumUnknownsSimplex(t);

    // temporaries for reading the stream
    int fromproc=-1;
    IdxT remote_dof, remote_extdof;

    // read information
    for ( size_t i=0; i<numData; ++i){
        recv >> fromproc >> remote_dof >> remote_extdof;
        if ( fromproc!=ProcCL::MyRank()){
            for ( Uint j=0; j<numUnk; ++j)
                dofProcList_[dof+j].insert( std::make_pair(fromproc, remote_dof+j));
            if ( isExtended || remote_extdof!=NoIdx){
                Assert(remote_extdof!=NoIdx && isExtended,
                    DROPSErrCL("ExchangeBuilderCL::ScatterLocalDOF: Remote process has not extended its local dof"),
                    DebugParallelNumC);
                for ( Uint j=0; j<numUnk; ++j)
                   dofProcList_[extdof+j].insert( std::make_pair(fromproc, remote_extdof+j));
            }
        }
    }

    return true;
}

// -----------------------------------------
// E X C H A N G E   M A T R I X   C L A S S
// -----------------------------------------

size_t ExchangeMatrixCL::NoIdx_= std::numeric_limits<size_t>::max();

/// \brief Determine the communication pattern for accumulating a matrix
void ExchangeMatrixCL::BuildCommPattern(const MatrixCL& mat,
        const ExchangeCL& RowEx, const ExchangeCL& ColEx)
/** To accumulate a matrix, the non-zeros which are stored by multiple processors, have to
    communicated among the processors. Therefore, each processor determines the non-zeroes
    it has to send to neighbors. And second, each processor determines how handle the received
    non-zeroes from a neighbor processor.
    \todo Are RowEx.GetProcs( i) and ColEx.GetProcs( j) already sorted?
    \param mat distributed matrix
    \param RowEx ExchangeCL that corresponds to row
    \param ColEx ExchangeCL that corresponds to column
*/
{
    // reset
    Clear();

    // Collect information about distributed non-zeroes in the matrix
    typedef std::map<int, std::vector<int>  > AODMap;
    typedef std::map<int, std::vector<IdxT> > RemoteDOFMap;

    AODMap aod;                              // mapping: proc -> array of displacements
    RemoteDOFMap remoteRowDOF, remoteColDOF; // mapping: proc -> remote (row|col)DOFs

    ProcNumCT NZonProcs( ProcCL::Size());   // stores the intersection
    for ( size_t i=0; i<mat.num_rows(); ++i) {
        if ( RowEx.IsDist(i)){
            for ( size_t nz=mat.row_beg(i); nz<mat.row_beg(i+1); ++nz){
                const size_t j= mat.col_ind(nz);
                if ( ColEx.IsDist( j)){     // here, i and j are both distributed
                    // determine all neighbor processors, that owns i *and* j as well
                    ProcNumCT NZonProcs= Intersect( RowEx, ColEx, i, j);

                    for (ProcNum_iter proc= NZonProcs.begin(); proc!=NZonProcs.end(); ++proc){
                        // mark the non-zero, that this non-zero should be sent to neighbor *proc
                        aod[*proc].push_back( (int)nz);
                        // determine the dof number on remote processor *proc
                        remoteRowDOF[*proc].push_back( RowEx.GetExternalIdxFromProc( i, *proc));
                        remoteColDOF[*proc].push_back( ColEx.GetExternalIdxFromProc( j, *proc));
                    }
                }
            }
        }
    }

    // Create MPI-Type for sending
    ExList_.reserve( aod.size());
    for (AODMap::const_iterator it= aod.begin(); it!=aod.end(); ++it){
        ExList_.push_back( SendNumDataCL<double>(it->first));
        ExList_.back().CreateDataType( it->second.size(), 1, Addr(it->second));
    }

    // Send "send-order"
    std::vector<ProcCL::RequestT> req;
    req.reserve( 2*remoteRowDOF.size());
    for (RemoteDOFMap::const_iterator it= remoteRowDOF.begin(); it!=remoteRowDOF.end(); ++it)
        req.push_back( ProcCL::Isend( it->second, it->first, 2211));
    for (RemoteDOFMap::const_iterator it= remoteColDOF.begin(); it!=remoteColDOF.end(); ++it)
        req.push_back( ProcCL::Isend( it->second, it->first, 2212));


    // Create receive sequence
    std::vector<IdxT> recvBufRowDOF, recvBufColDOF;
    RecvBuf_.resize( ExList_.size());
    Coupl_.resize( ExList_.size());
    for (size_t ex=0; ex<ExList_.size(); ++ex){
        // receive sequence
        int messagelength= ProcCL::GetMessageLength<IdxT>( ExList_[ex].GetReceiver(), 2211);
        recvBufRowDOF.resize( messagelength);
        recvBufColDOF.resize( messagelength);
        ProcCL::Recv( Addr(recvBufRowDOF), messagelength, ExList_[ex].GetReceiver(), 2211);
        ProcCL::Recv( Addr(recvBufColDOF), messagelength, ExList_[ex].GetReceiver(), 2212);

        // create sequence
        Coupl_[ex].resize( messagelength);
        for (size_t k=0; k<recvBufRowDOF.size(); ++k)
            Coupl_[ex][k]= GetPosInVal(recvBufRowDOF[k], recvBufColDOF[k], mat);

        // reserve memory for receiving
        RecvBuf_[ex].resize( messagelength);
    }
    // Before cleaning up remoteRowDOF and remoteColDOF, check if all messages has been received
    ProcCL::WaitAll( req);
}

MatrixCL ExchangeMatrixCL::Accumulate(const MatrixCL& mat)
/** According to the communication pattern accumulate the distributed non-zero elements of
    the matrix mat.
    \pre BuildCommPattern must have been called for the pattern of the matrix \a mat
    \param  mat input, distributed matrix
    \return matrix, with accumulated non-zeros
 */
{

    // Make a copy of distributed values
    MatrixCL result( mat);

    // Initialize communication
    std::vector<ProcCL::RequestT> send_req( ExList_.size());
    std::vector<ProcCL::RequestT> recv_req( ExList_.size());
    for (size_t ex=0; ex<ExList_.size(); ++ex){
        send_req[ex]= ExList_[ex].Isend( mat.val(), 2213, 0);
        recv_req[ex]= ProcCL::Irecv( RecvBuf_[ex], ExList_[ex].GetReceiver(), 2213);
    }

    // do accumulation
    for ( size_t ex=0; ex<ExList_.size(); ++ex){
        // wait until non-zeros have been received
        ProcCL::Wait( recv_req[ex]);
        // add received non-zeros
        for ( size_t nz=0; nz<RecvBuf_[ex].size(); ++nz){
            if (Coupl_[ex][nz]!=NoIdx_)
                result.val()[Coupl_[ex][nz]]+= RecvBuf_[ex][nz];
        }
    }

    // wait until send are finished before leaving this routine
    ProcCL::WaitAll(send_req);
    return result;
}

} // end of namespace DROPS
