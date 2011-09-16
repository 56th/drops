/// \file mpistream.h
/// \brief SendStreamCL and ReicStreamCL for transfering objects (simplexes + data) as byte-messages.
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
#include "parallel/parallel.h"

#ifndef DROPS_MPISTREAM_H
#define DROPS_MPISTREAM_H

namespace DROPS {

namespace DiST {
namespace Helper {

class GeomIdCL; ///< forward declaration for operator>>/<<

/// \brief Output-streambuf based on std::stringbuf.
class MPIostringbufCL : public std::stringbuf
{
  public:
    typedef std::stringbuf base_type;

    /// \brief Non-blocking send to process 'dest'.
    ProcCL::RequestT Isend(int dest, int tag) {
        const std::streamsize size= pptr() - eback();
        return ProcCL::Isend( eback(), size, dest, tag);
    }
    /// \brief Reset the buffer to the empty default-state. This releases the memory of the buffer.
    void clearbuffer () { str( std::string()); }
};

/// \brief Outgoing stream.
class SendStreamCL : public std::ostream
/** This stream is employed as a means of transport for sending data (integers,
    vertices, elements, etc.) out of the process towards other processes. <p>
    It is derived from the ostream class of the Standard IOstream Library
    and it is characterized by saving all kind of data as a string. <p>
    It's counterpart is the class RecvStreamCL. <p>
*/
{
  public:
    typedef std::ostream base_type;

  private:
    bool binary_;         ///< flag for binary sending/receiving
    MPIostringbufCL buf_; ///< string based output-buffer

  public:
    SendStreamCL( const bool binary= true) : base_type( &buf_), binary_( binary) {}

    inline bool isBinary() const {return binary_;}
    /// \brief Non-blocking send to process 'dest'.
    inline ProcCL::RequestT Isend(int dest, int tag= 5) { return buf_.Isend( dest, tag); }
    /// \brief Return a copy of the string
    inline std::string str () const { return buf_.str(); }
    /// \brief Reset the buffer to the empty default-state. This releases the memory of the buffer.
    void clearbuffer () { buf_.clearbuffer(); }
};

/// \brief Incoming stream.
class RecvStreamCL : public std::istringstream
/** This stream is employed as a mean of transport for receiving incoming data
    (integers, vertices, elements, etc.) sent from other process. <p>
    It is derived from the istringstream class of the Standard IOstream Library
    and it is characterized by saving all kind of data as a string. <p>
    It's counterpart is the class SendStreamCL.
*/
{
public:
    typedef std::istringstream base;

  private:
    bool binary_;   ///< flag for binary sending/receiving

  public:
    /// @param[in] binary is true if the stream store the data in binary; in ASCII otherwise.
    RecvStreamCL( const bool binary=true)
        : std::istringstream( binary ? std::ios_base::binary : std::ios_base::in), binary_(binary) {}
    RecvStreamCL( const SendStreamCL& s)
        : std::istringstream( s.str()), binary_(s.isBinary()) {}
    /// \brief Gives back how many bytes of data contains our stream at the time.
    inline int getbufsize() const { return ((int) const_cast<RecvStreamCL*>(this)->tellg())+1; }
    inline bool isBinary() const {return binary_;}
    /// \brief Blocking receive from process 'source'.
    void Recv(int source, int tag=5);
    /// \brief Go back to the beginning of the stream.
    inline void resetbuffer() { this->seekg(0); }
};


/// \brief Helper union for writing and reading numbers in streams
/** Technically, this is undefined behavior (reading another member of the union as was previously written). All major compilers recognize this as idiom for type-punning. Alternatively, casting the address of *any* object to char* is save by the standard. Maybe go for that.*/
template <typename T>
union ToBinary
{
    T    value;
    char binary[sizeof(T)];
};


/// \brief operator << for SendStreamCL
template<typename T>
SendStreamCL& operator<<( SendStreamCL& os, const T& t)
{
    if (os.isBinary()) {
        ToBinary<T> bin;
        bin.value= t;
        os.write( bin.binary, sizeof(T));
    } else {
        SendStreamCL::base_type& oss= dynamic_cast<SendStreamCL::base_type&>(os);
        oss << t << ' ';
    }
    return os;
}

/// \brief operator >> for RecvStreamCL
template<typename T>
RecvStreamCL& operator>>( RecvStreamCL& is, T& t)
{
    if (is.isBinary()) {
        ToBinary<T> bin;
        is.read( bin.binary, sizeof(T));
        t= bin.value;
    } else {
        RecvStreamCL::base& iss= dynamic_cast<RecvStreamCL::base&>(is);
        iss >> t;
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