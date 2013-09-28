/// \file exchange.tpp
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

namespace DROPS{

#ifdef _PAR

template <typename T>
void SendNumDataCL<T>::freeType()
/** If the datatype mpidatatype_ has been commited to MPI unregister
    this type.
*/
{
    if (mpidatatype_!=ProcCL::NullDataType){
        ProcCL::Free(mpidatatype_);
        mpidatatype_= ProcCL::NullDataType;
    }
}

/// \brief Create Datatype for sending
template <typename T>
void SendNumDataCL<T>::CreateDataType(const int count,
        const int blocklength,
        const int array_of_displacements[])
/** This function creates a MPI datatype that is responsible for gathering
    the data.
    \param count number of elements to be send
    \param blocklength length of each block
    \param array_of_displacements position of the elements
*/
{
    if (mpidatatype_!=ProcCL::NullDataType)
        freeType();
    mpidatatype_ = ProcCL::CreateBlockIndexed<T>(count, blocklength, array_of_displacements);
    ProcCL::Commit(mpidatatype_);
    minlengthvec_= array_of_displacements[count-1]+blocklength;
}

template <typename T>
template <typename VectorT>
inline ProcCL::RequestT SendNumDataCL<T>::Isend(
     const VectorT& v, int tag, Ulint offset) const
/** Call the MPI Isend routine via the ProcCL. A request handler is return
    for asking if the vector \a v can be overwritten. This function uses the
    internal MPI datatype mpidatattype_ for sending.
    \param v data to be send
    \param tag a message tag
    \param offset start sending form the offset element in \a v. Used for
                  blocked vectors.
*/
{
    Assert( (int)v.size()-(int)offset>=minlengthvec_,
        DROPSErrCL("SendNumDataCL::Isend: Given vector does not hold enough entries"),
        DebugParallelNumC);
    return ProcCL::Isend( Addr(v)+offset, 1, mpidatatype_, toproc_, tag);
}

template <typename T>
inline ProcCL::RequestT SendNumDataCL<T>::Isend(
     const double* v, int tag, Ulint offset) const
/** Call the MPI Isend routine via the ProcCL. A request handler is return
    for asking if the vector \a v can be overwritten. This function uses the
    internal MPI datatype mpidatattype_ for sending.
    \param v data to be send
    \param tag a message tag
    \param offset start sending form the offset element in \a v. Used for
                  blocked vectors.
*/
{
/*    if ((int)v.size()-(int)offset<minlengthvec_){
        printf("Proc %d, v should contain at least %d elements to be sent to %d\n",
            ProcCL::MyRank(), minlengthvec_, toproc_);
    }
    Assert( (int)v.size()-(int)offset>=minlengthvec_,
        DROPSErrCL("SendNumDataCL::Isend: Given vector does not hold enough entries"),
        DebugParallelNumC);*/
    return ProcCL::Isend( v+offset, 1, mpidatatype_, toproc_, tag);
}

template <typename T>
ProcCL::RequestT RecvNumDataCL<T>::Irecv(int tag, VectorBaseCL<T>& recvBuf,
    Ulint offset) const
/** This procedure receives the datas with an non-blocking non-synchronous MPI
    Receive.
    \param tag     used tag for communication
    \param recvBuf buffer for storing received unknowns
    \param offset  first position, where to store unknowns
    \pre recvBuf has to be big enough to store all data
    \pre no other procedures is allowed to work on the memory (in particular no other Irecv!)
    \return MPI request handle for this receive
*/
{
    Assert( recvBuf.size()+offset>=sysnums_.size(),
        DROPSErrCL("RecvNumDataCL::Irecv: Receive buffer is too small"),
        DebugParallelNumC);
    return ProcCL::Irecv( Addr(recvBuf)+offset, sysnums_.size(), fromproc_, tag);
}

template <typename T>
void RecvNumDataCL<T>::Accumulate(VectorBaseCL<T>& v, Ulint offsetV,
    const VectorBaseCL<T>& recvBuf, Ulint offsetRecv) const
/** This procedure accumulates received data. It assumes, that the data has been
    received.
    \param v          original value, that contains all local unknowns
    \param offsetV    start element of the vector, that contains local unknowns
    \param recvBuf    vector of all received elements
    \param offsetRecv first element in receive buffer
    \pre Communication has to be done before entering this procedure
*/
{
    for ( size_t i=0; i<sysnums_.size(); ++i){
        v[sysnums_[i]+offsetV]+= recvBuf[i+offsetRecv];
    }
}

template <typename T>
void  RecvNumDataCL<T>::Assign(VectorBaseCL<T>& v, Ulint offsetV,
    const VectorBaseCL<T>& recvBuf, Ulint offsetRecv) const
/** This procedure assigns the received data. It assumes, that the data has been
    received.
    \param v          original value, that contains all local unknowns
    \param offsetV    start element of the vector, that contains local unknowns
    \param recvBuf    vector of accumulated elements (received by DoF owner)
    \param offsetRecv first element in receive buffer
    \pre Communication has to be done before entering this procedure
*/
{
    for ( size_t i=0; i<sysnums_.size(); ++i){
        v[sysnums_[i]+offsetV] = recvBuf[i+offsetRecv];
    }
}

bool ExchangeCL::IsOnProc( IdxT localdof, int proc) const
/** Check in dofProcList_, if \a proc owns a copy. */
{
    const DOFInfoT& remoteDOF= dofProcList_[localdof];
    return remoteDOF.find(proc) != remoteDOF.end();
}

inline IdxT ExchangeCL::GetExternalIdxFromProc( const IdxT localdof, int proc) const
/** Get dof number on process \a proc. If the dof is not located on proc, return NoIdx. */
{
    const DOFInfoT& remoteDOF= dofProcList_[localdof];
    DOFInfoList_const_iterator it= remoteDOF.find( proc);
    if ( it != remoteDOF.end())
        return it->second;
    else
        return NoIdx;
}

inline ExchangeMatrixCL::ProcNumCT ExchangeMatrixCL::Intersect(
    const ExchangeCL& RowEx, const ExchangeCL& ColEx, const size_t i, const size_t j) const
{
    ProcNumCT result;
    ExchangeCL::DOFInfoList_const_iterator it_row( RowEx.GetProcListBegin(i)),
        it_col( ColEx.GetProcListBegin(j));
    const ExchangeCL::DOFInfoList_const_iterator end_row( RowEx.GetProcListEnd(i)),
        end_col( ColEx.GetProcListEnd(j));
    while (it_row!=end_row && it_col!=end_col){
        if ( it_row->first == it_col->first){
            result.push_back( it_row->first);
            ++it_row; ++it_col;
        }
        else if ( it_row->first < it_col->first)
            ++it_row;
        else
            ++it_col;
    }
    return result;
}

#endif // parallel

} // end of namespace DROPS
