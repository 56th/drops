/// \file mpistream.h
/// \brief SendStreamCL and RecvStreamCL for transfering objects (simplexes + data) as byte-messages.
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

#include "misc/container.h"
#include "misc/utils.h"
#include "parallel/parallel.h"

#ifndef DROPS_MPISTREAM_H
#define DROPS_MPISTREAM_H

namespace DROPS {

namespace DiST {
namespace Helper {

class GeomIdCL; ///< forward declaration for operator>>/<<

/// \brief Streambuf for MPI-messages based on std::stringbuf.
/// The MPI operations are forwarded to ProcCL.
/// The messages are buffered in a std::stringbuf.
class MPIstringbufCL : public std::stringbuf
{
  public:
    typedef std::stringbuf base_type;

    explicit MPIstringbufCL (std::ios_base::openmode which= std::ios_base::in | std::ios_base::out)
        : base_type( which) {}
    explicit MPIstringbufCL (const std::string& s, std::ios_base::openmode which= std::ios_base::in | std::ios_base::out)
        : base_type( s, which) {}

    /// \brief Non-blocking send to process 'dest'. The buffer may not be modified, until the send is completed.
    ProcCL::RequestT Isend (int dest, int tag);
    /// \brief Blocking receive from process 'source'. The buffer is reset before new data is filled in.
    void Recv (int source, int tag);
    /// \brief Reset the buffer to the empty default-state. This releases the memory of the buffer.
    void clearbuffer () { str( std::string()); }
};

/// \brief Output-stream.
/// This stream is employed for sending data (integers, vertices, elements, etc.)
/// out of the process towards other processes.<p>
/// It is derived from the ostream class of the Standard IOstream Library.
/// The data is interpreted as char-array in memory (isBinary()==true) or via its stream-defined ascii-representation. In ascii, items are terminated by the value of SendRecvStreamAsciiTerminatorC.<p>
/// In ascii, FP-values are written with 17 significant decimal digits. This suffices for 8-byte doubles.
/// It's counterpart is the class RecvStreamCL.<p>
class SendStreamCL : public std::ostream
{
  public:
    typedef std::ostream base_type;

  private:
    bool binary_;        ///< flag for binary sending/receiving
    MPIstringbufCL buf_; ///< string based output-buffer

  public:
    SendStreamCL (const bool binary= true)
        : base_type( &buf_), binary_( binary), buf_( std::ios_base::out)
    { precision( 17); }

    inline bool isBinary() const { return binary_; }
    /// \brief Non-blocking send to process 'dest'.
    inline ProcCL::RequestT Isend(int dest, int tag= 5) { return buf_.Isend( dest, tag); }
    /// \brief Return a copy of the string
    inline std::string str () const { return buf_.str(); }
    /// \brief Reset the buffer to the empty default-state. This releases the memory of the buffer.
    void clearbuffer () { buf_.clearbuffer(); }
};

/// \brief Input stream.
/// This stream is employed for receiving data (integers, vertices, elements, etc.)
/// from other processes.<p>
/// It is derived from the istream class of the Standard IOstream Library.
/// The data is interpreted as char-array in memory (isBinary()==true) or via its stream-defined ascii-representation. In ascii, items are terminated by the value of SendRecvStreamAsciiTerminatorC.<p>
/// It's counterpart is the class SendStreamCL. <p>
class RecvStreamCL : public std::istream
{
  public:
    typedef std::istream base_type;

  private:
    bool binary_;        ///< flag for binary sending/receiving
    MPIstringbufCL buf_; ///< string based output-buffer

  public:
    /// @param[in] binary is true if the stream store the data in binary; in ASCII otherwise.
    RecvStreamCL (const bool binary= true)
        : base_type( &buf_), binary_( binary), buf_( std::ios_base::in)
    { unsetf( std::ios_base::skipws); }

    RecvStreamCL( const SendStreamCL& s)
        : base_type( &buf_), binary_( s.isBinary()), buf_( s.str(), std::ios_base::in) {}

    inline bool isBinary() const { return binary_; }
    /// \brief Blocking receive from process 'source'. Removes any prior content of the buffer.
    void Recv(int source, int tag= 5) { buf_.Recv( source, tag); }
    /// \brief Reset the buffer to the empty default-state. This releases the memory of the buffer.
    void clearbuffer () { buf_.clearbuffer(); }
};

/// \brief The data-item-terminator used in ascii-mode.
const char SendRecvStreamAsciiTerminatorC= ' ';

/// \brief Helper union for writing and reading numbers in streams
/** Technically, this is undefined behavior (reading another member of the union as was previously written). However, all major compilers recognize this as idiom for type-punning. Alternatively, casting the address of *any* object to char* is save by the standard. Maybe go for that.*/
template <typename T>
union ToBinary
{
    T    value;
    char binary[sizeof(T)];
};


/// \brief operator << for SendStreamCL
template<typename T>
SendStreamCL& operator<< (SendStreamCL& os, const T& t)
{
    if (os.isBinary()) {
        ToBinary<T> bin;
        bin.value= t;
        os.write( bin.binary, sizeof(T));
    } else {
        SendStreamCL::base_type& ostr= static_cast<SendStreamCL::base_type&>( os);
        ostr << t << SendRecvStreamAsciiTerminatorC;
    }
    return os;
}

/// \brief operator >> for RecvStreamCL
template<typename T>
RecvStreamCL& operator>> (RecvStreamCL& is, T& t)
{
    if (is.isBinary()) {
        ToBinary<T> bin;
        is.read( bin.binary, sizeof(T));
        t= bin.value;
    } else {
        RecvStreamCL::base_type& istr= static_cast<RecvStreamCL::base_type&>( is);
        istr >> t;
        RecvStreamCL::char_type c;
        istr.get( c);
        if (c != SendRecvStreamAsciiTerminatorC)
            throw DROPSErrCL( "RecvStreamCL& operator>>( RecvStreamCL& is, T& t): Ascii item-terminator not found.\n");
    }
    return is;
}


/// \brief input/output of GeomIdCL
/// @{
SendStreamCL& operator<< (SendStreamCL&, const GeomIdCL&);
RecvStreamCL& operator>> (RecvStreamCL&,       GeomIdCL&);
///@}

/// \brief input/output of Point3DCL
/// @{
SendStreamCL& operator<< (SendStreamCL&, const Point3DCL&);
RecvStreamCL& operator>> (RecvStreamCL&,       Point3DCL&);
///@}


} // end of namespace Helper
} // end of namespace DiSt
} // end of namespace DROPS

#endif