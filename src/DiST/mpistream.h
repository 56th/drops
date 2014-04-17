/// \file mpistream.h
/// \brief SendStreamCL and RecvStreamCL for transferring objects (simplexes + data) as byte-messages.
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
 * Copyright 2012 LNM/SC RWTH Aachen, Germany
*/

#include "misc/container.h"
#include "misc/utils.h"
#include "parallel/parallel.h"

#ifndef DROPS_MPISTREAM_H
#define DROPS_MPISTREAM_H

namespace DROPS {

namespace DiST {

/// \brief Default for the MPI-message format used by DiST.
/// True on program startup. Can be reassigned.
extern bool use_binaryMPIstreams;

struct GeomIdCL; ///< forward declaration for operator>>/<<

/// \brief Streambuf for MPI-messages based on std::stringbuf.
/// The MPI operations are forwarded to ProcCL.
/// The messages are buffered in a std::stringbuf.
/// The input/output area can be accessed as char_type-sequences with begin(), cur() and end().
class MPIstringbufCL : public std::stringbuf
{
  public:
    typedef std::stringbuf base_type;
    typedef base_type::char_type   char_type;
    typedef base_type::traits_type traits_type;
    typedef base_type::int_type    int_type;
    typedef base_type::pos_type    pos_type;
    typedef base_type::off_type    off_type;

    explicit MPIstringbufCL (std::ios_base::openmode which= std::ios_base::in | std::ios_base::out)
        : base_type( which) {}
    explicit MPIstringbufCL (const std::string& s, std::ios_base::openmode which= std::ios_base::in | std::ios_base::out)
        : base_type( s, which) {}

    /// \brief Non-blocking send to process 'dest'. The buffer may not be modified, until the send is completed.
    ProcCL::RequestT Isend (int dest, int tag);
    /// \brief Blocking receive from process 'source'. The buffer is reset before new data is filled in.
    void Recv (int source, int tag);
    /// \brief Reset the buffer to contain a copy of s. The prior memory is released.
    void clearbuffer (const std::string& s= std::string()) { str( s); }

    ///\brief access the input/output area (an array of char_type); cur is the current i/o-position.
    ///@{
    char_type* begin(std::ios_base::openmode m) { return m == std::ios_base::in ? eback() : pbase(); }
    char_type* cur  (std::ios_base::openmode m) { return m == std::ios_base::in ? gptr()  : pptr(); }
    char_type* end  (std::ios_base::openmode m) { return m == std::ios_base::in ? egptr() : epptr(); }

    const char_type* begin(std::ios_base::openmode m) const { return m == std::ios_base::in ? eback() : pbase(); }
    const char_type* cur  (std::ios_base::openmode m) const { return m == std::ios_base::in ? gptr()  : pptr(); }
    const char_type* end  (std::ios_base::openmode m) const { return m == std::ios_base::in ? egptr() : epptr(); }
    ///@}
};


/// \brief Base for the MPI-output-streams.
/// This stream defines the stream-format for sending data (integers, vertices, elements, etc.)
/// out of the process towards other processes.
/// Formats are defined by writing output-operators.<p>
/// The data is interpreted as char-array in memory (isBinary()==true) or via its MPIostream-defined ascii-representation.
/// In ascii, items are terminated by the value of SendRecvStreamAsciiTerminatorC.<p>
/// In ascii, FP-values are written with 17 significant decimal digits. This suffices for 8-byte doubles.
class MPIostreamCL : public std::ostream
{
  private:
    bool binary_; ///< flag for binary sending/receiving

    template <class T>
      inline MPIostreamCL& write_fundamental_type (const T& t);

  public:
    typedef std::ostream base_type;

    MPIostreamCL (std::streambuf* buf, bool binary)
        : base_type( buf), binary_( binary) { if (!isBinary()) precision( 17); }

    inline bool isBinary ()             const { return binary_; }
    inline void setBinary (bool binary)       { binary_= binary; }

    /// \brief Output of C++'s arithmetic types
    ///@{
    MPIostreamCL& operator<< (bool val)          { return write_fundamental_type( val); }
    MPIostreamCL& operator<< (short val)         { return write_fundamental_type( val); }
    MPIostreamCL& operator<< (unsigned short val){ return write_fundamental_type( val); }
    MPIostreamCL& operator<< (int val)           { return write_fundamental_type( val); }
    MPIostreamCL& operator<< (unsigned int val)  { return write_fundamental_type( val); }
    MPIostreamCL& operator<< (long val)          { return write_fundamental_type( val); }
    MPIostreamCL& operator<< (unsigned long val) { return write_fundamental_type( val); }
    MPIostreamCL& operator<< (float val)         { return write_fundamental_type( val); }
    MPIostreamCL& operator<< (double val)        { return write_fundamental_type( val); }
    MPIostreamCL& operator<< (long double val)   { return write_fundamental_type( val); }
    MPIostreamCL& operator<< (const void* val)   { return write_fundamental_type( val); }
    ///@}
    /// \brief Output of single char.
    /// These are global functions for std::ostream. There seems to be no good reason for this.
    ///@{
    MPIostreamCL& operator<< (char val)          { return write_fundamental_type( val); }
    MPIostreamCL& operator<< (signed char val)   { return write_fundamental_type( val); }
    MPIostreamCL& operator<< (unsigned char val) { return write_fundamental_type( val); }
    ///@}
    // MPIostreamCL& operator<< (streambuf* sb); // Implement, if needed.

    MPIostreamCL& write (const char* s , std::streamsize n) { base_type::write( s, n); return *this; }
};

/// \brief Input stream.
/// This stream is employed for receiving data (integers, vertices, elements, etc.)
/// from other processes.<p>
/// It is derived from the std::istream.
/// The data is interpreted as char-array in memory (isBinary()==true) or via its stream-defined ascii-representation.
/// In ascii, terminating SendRecvStreamAsciiTerminatorC are consumed.<p>
class MPIistreamCL : public std::istream
{
  private:
    bool binary_;        ///< flag for binary sending/receiving

    template <class T>
      inline MPIistreamCL& read_fundamental_type (T& t);

  public:
    typedef std::istream base_type;

    MPIistreamCL (std::streambuf* buf, bool binary)
        : base_type( buf), binary_( binary) { if (!isBinary()) unsetf( std::ios_base::skipws); }

    inline bool isBinary ()             const { return binary_; }
    inline void setBinary (bool binary)       { binary_= binary; }

    /// \brief Input of C++'s arithmetic types
    ///@{
    MPIistreamCL& operator>> (bool& val)          { return read_fundamental_type( val); }
    MPIistreamCL& operator>> (short& val)         { return read_fundamental_type( val); }
    MPIistreamCL& operator>> (unsigned short& val){ return read_fundamental_type( val); }
    MPIistreamCL& operator>> (int& val)           { return read_fundamental_type( val); }
    MPIistreamCL& operator>> (unsigned int& val)  { return read_fundamental_type( val); }
    MPIistreamCL& operator>> (long& val)          { return read_fundamental_type( val); }
    MPIistreamCL& operator>> (unsigned long& val) { return read_fundamental_type( val); }
    MPIistreamCL& operator>> (float& val)         { return read_fundamental_type( val); }
    MPIistreamCL& operator>> (double& val)        { return read_fundamental_type( val); }
    MPIistreamCL& operator>> (long double& val)   { return read_fundamental_type( val); }
    MPIistreamCL& operator>> (void*& val)         { return read_fundamental_type( val); }
    ///@}
    /// \brief Input of single char.
    /// These are global functions for std::istream. There seems to be no good reason for this.
    ///@{
    MPIistreamCL& operator>> (char& val)          { return read_fundamental_type( val); }
    MPIistreamCL& operator>> (signed char& val)   { return read_fundamental_type( val); }
    MPIistreamCL& operator>> (unsigned char& val) { return read_fundamental_type( val); }
    ///@}
    // MPIistreamCL& operator>> (streambuf* sb); // Implement, if needed.

    MPIistreamCL& read (char* s, std::streamsize n ) { base_type::read( s, n); return *this; }

};


/// \brief MPI-output-stream with MPIstringbufCL-buffer
/// This stream is employed for formatting and buffering outgoing MPI-streams.
/// Copy and assignment have value-semantics, i.e., the buffer contents are copied/assigned.
class SendStreamCL : public MPIostreamCL
{
  public:
    typedef MPIostreamCL base_type;

  private:
    MPIstringbufCL buf_; ///< string based output-buffer

  public:
    explicit SendStreamCL (bool binary= use_binaryMPIstreams)
        : base_type( &buf_, binary), buf_( std::ios_base::out) {}
    SendStreamCL (const SendStreamCL& s);
    SendStreamCL& operator= (const SendStreamCL& s);

    /// \brief Non-blocking send to process 'dest'.
    inline ProcCL::RequestT Isend(int dest, int tag= 5) { return buf_.Isend( dest, tag); }
    /// \brief Return a copy of the string
    inline std::string str () const { return buf_.str(); }
    /// \brief Reset the buffer to the empty default-state. This releases the memory of the buffer.
    void clearbuffer (const std::string& s= std::string()) { buf_.clearbuffer( s); }

    ///\brief access the input/output area (an array of char_type); cur is the current i/o-position.
    ///@{
    char_type* begin() { return buf_.begin( std::ios_base::out); }
    char_type* cur  () { return buf_.cur(   std::ios_base::out); }
    char_type* end  () { return buf_.end(   std::ios_base::out); }

    const char_type* begin() const { return buf_.begin( std::ios_base::out); }
    const char_type* cur  () const { return buf_.cur(   std::ios_base::out); }
    const char_type* end  () const { return buf_.end(   std::ios_base::out); }
    ///@}
};

/// \brief Input stream.
/// This stream is employed for formatting and buffering outgoing MPI-streams.
/// Copy and assignment have value-semantics, i.e., the buffer contents are copied/assigned.
class RecvStreamCL : public MPIistreamCL
{
  public:
    typedef MPIistreamCL base_type;

  private:
    MPIstringbufCL buf_; ///< string based output-buffer

  public:
    /// @param[in] binary is true if the stream store the data in binary; in ASCII otherwise.
    explicit RecvStreamCL (bool binary= use_binaryMPIstreams)
        : base_type( &buf_, binary), buf_( std::ios_base::in) {}
    RecvStreamCL (const RecvStreamCL& s);
    RecvStreamCL& operator= (const RecvStreamCL& s);

    RecvStreamCL( const SendStreamCL& s)
        : base_type( &buf_, s.isBinary()), buf_( s.str(), std::ios_base::in) {}

    /// \brief Blocking receive from process 'source'. Removes any prior content of the buffer.
    RecvStreamCL& Recv(int source, int tag= 5)
        { buf_.Recv( source, tag); return *this; }
    /// \brief Reset the buffer to the empty default-state. This releases the memory of the buffer.
    void clearbuffer (const std::string& s= std::string()) { buf_.clearbuffer( s); }

    ///\brief access the input/output area (an array of char_type); cur is the current i/o-position.
    ///@{
    char_type* begin() { return buf_.begin( std::ios_base::in); }
    char_type* cur  () { return buf_.cur(   std::ios_base::in); }
    char_type* end  () { return buf_.end(   std::ios_base::in); }

    const char_type* begin() const { return buf_.begin( std::ios_base::in); }
    const char_type* cur  () const { return buf_.cur(   std::ios_base::in); }
    const char_type* end  () const { return buf_.end(   std::ios_base::in); }
    ///@}

    friend MPIistreamCL& operator>> (MPIistreamCL& is, RecvStreamCL& sub);
};

/// \brief Streambuf for MPI-messages with reference semantics.
/// The buffer uses an external char_type field or another buffers memory for its own controlled sequences.
/// The sequences have a fixed length.
/// This is a standard-conforming streambuf.
class MPIrefbufCL : public std::streambuf
{
  public:
    typedef std::streambuf base_type;
    typedef base_type::char_type   char_type;
    typedef base_type::traits_type traits_type;
    typedef base_type::int_type    int_type;
    typedef base_type::pos_type    pos_type;
    typedef base_type::off_type    off_type;

  private:
    std::ios_base::openmode mode_;

  protected:
    MPIrefbufCL* setbuf (char_type* b, std::streamsize s);
    pos_type seekoff (off_type off, std::ios_base::seekdir base,
        std::ios_base::openmode m= std::ios_base::in | std::ios_base::out);
    pos_type seekpos (pos_type p, std::ios_base::openmode m= std::ios_base::in | std::ios_base::out)
        { return seekoff( p, std::ios_base::beg, m); }
    std::streamsize showmanyc () { return 0; }
    std::streamsize xsgetn (char_type* p, std::streamsize n);
    std::streamsize xsputnb(const char_type* p, std::streamsize n);

  public:
    explicit MPIrefbufCL (std::streambuf::char_type* b, std::streamsize s,
        std::ios_base::openmode m= std::ios_base::in | std::ios_base::out)
        : mode_( m) { this->setbuf( b, s); }

    MPIrefbufCL* pubsetbuf (char_type* b, std::streamsize s) { return setbuf( b, s); }

    ///\brief access the input/output area (an array of char_type); cur is the current i/o-position.
    ///@{
    char_type* begin(std::ios_base::openmode m) { return m == std::ios_base::in ? eback() : pbase(); }
    char_type* cur  (std::ios_base::openmode m) { return m == std::ios_base::in ? gptr()  : pptr(); }
    char_type* end  (std::ios_base::openmode m) { return m == std::ios_base::in ? egptr() : epptr(); }

    const char_type* begin(std::ios_base::openmode m) const { return m == std::ios_base::in ? eback() : pbase(); }
    const char_type* cur  (std::ios_base::openmode m) const { return m == std::ios_base::in ? gptr()  : pptr(); }
    const char_type* end  (std::ios_base::openmode m) const { return m == std::ios_base::in ? egptr() : epptr(); }
    ///@}
};

/// \brief Output stream with reference semantics.
/// This stream is employed for formatting and buffering outgoing MPI-streams.
class RefMPIostreamCL : public MPIostreamCL
{
  public:
    typedef MPIostreamCL base_type;

  private:
    MPIrefbufCL buf_; ///< reference output-buffer

  public:
    explicit RefMPIostreamCL (char_type* p= 0, std::streamsize n= 0,
                              bool binary= use_binaryMPIstreams)
        : base_type( &buf_, binary), buf_( p, n, std::ios_base::out) {}

    ///\brief access the input/output area (an array of char_type); cur is the current i/o-position.
    ///@{
    char_type* begin() { return buf_.begin( std::ios_base::out); }
    char_type* cur  () { return buf_.cur(   std::ios_base::out); }
    char_type* end  () { return buf_.end(   std::ios_base::out); }

    const char_type* begin() const { return buf_.begin( std::ios_base::out); }
    const char_type* cur  () const { return buf_.cur(   std::ios_base::out); }
    const char_type* end  () const { return buf_.end(   std::ios_base::out); }
    ///@}
};

/// \brief Input stream with reference semantics.
/// This stream is employed for formatting and buffering incoming MPI-streams.
class RefMPIistreamCL : public MPIistreamCL
{
  public:
    typedef MPIistreamCL base_type;

  private:
    MPIrefbufCL buf_;    ///< string based output-buffer

  public:
    /// @param[in] binary is true if the stream store the data in binary; in ASCII otherwise.
    explicit RefMPIistreamCL (char_type* p= 0, std::streamsize n= 0,
                              bool binary= use_binaryMPIstreams)
        : base_type( &buf_, binary), buf_( p, n, std::ios_base::in) {}

    MPIrefbufCL* setbuf( char_type* p, std::streamsize n) { return buf_.pubsetbuf( p, n); }

    ///\brief access the input/output area (an array of char_type); cur is the current i/o-position.
    ///@{
    char_type* begin() { return buf_.begin( std::ios_base::in); }
    char_type* cur  () { return buf_.cur(   std::ios_base::in); }
    char_type* end  () { return buf_.end(   std::ios_base::in); }

    const char_type* begin() const { return buf_.begin( std::ios_base::in); }
    const char_type* cur  () const { return buf_.cur(   std::ios_base::in); }
    const char_type* end  () const { return buf_.end(   std::ios_base::in); }
    ///@}
};

/// \brief file-output-stream with std::filebuf-buffer
/// This stream is employed for formatting and buffering outgoing file-streams.
class SerialOFStreamCL : public MPIostreamCL
{
  public:
    typedef MPIostreamCL base_type;

  private:
    std::filebuf buf_; ///< string based output-buffer

  public:
    explicit SerialOFStreamCL (const char * s, bool binary= use_binaryMPIstreams)
        : base_type( &buf_, binary) { buf_.open(s, std::ios_base::out); }
    explicit SerialOFStreamCL (const std::string s, bool binary= use_binaryMPIstreams)
        : base_type( &buf_, binary) { buf_.open( s.c_str(), std::ios_base::out); }
};

/// \brief file-input-stream with std::filebuf-buffer
/// This stream is employed for formatting and buffering incoming file-streams.
class SerialIFStreamCL : public MPIistreamCL
{
  public:
    typedef MPIistreamCL base_type;

  private:
    std::filebuf buf_; ///< string based output-buffer

  public:
    explicit SerialIFStreamCL (const char * s, bool binary= use_binaryMPIstreams)
        : base_type( &buf_, binary) { buf_.open(s, std::ios_base::in); }
    explicit SerialIFStreamCL (const std::string s, bool binary= use_binaryMPIstreams)
        : base_type( &buf_, binary) { buf_.open( s.c_str(), std::ios_base::in); }
};

/// \brief The data-item-terminator used in ascii-mode.
const char SendRecvStreamAsciiTerminatorC= ' ';


/// \brief input/output of GeomIdCL
/// @{
MPIostreamCL& operator<< (MPIostreamCL&, const GeomIdCL&);
MPIistreamCL& operator>> (MPIistreamCL&,       GeomIdCL&);
///@}

/// \brief input/output of SVectorCL (Point2DCL, Point3DCL, BaryCoordCL)
/// @{
template <Uint rows>
  inline MPIostreamCL&
  operator<< (MPIostreamCL& os, const SVectorCL<rows>& p);

template <Uint rows>
  inline MPIistreamCL&
  operator>> (MPIistreamCL& is, SVectorCL<rows>& p);
///@}

///\brief input/output of MPIstreamCL as a sub-object of an MPIi/ostreamCL.
/// This allows for substructured MPI-messages.
/// The range [begin, cur) of the substream sub is considered for output with operator<<, not [begin, end);
/// In substreams created with operator>>, there holds begin() == cur().
///@{
MPIostreamCL& operator<< (MPIostreamCL& os, const SendStreamCL& sub);
MPIostreamCL& operator<< (MPIostreamCL& os, const RecvStreamCL& sub);
MPIistreamCL& operator>> (MPIistreamCL& is,       RecvStreamCL& sub);

MPIostreamCL&    operator<< (MPIostreamCL&    os, const RefMPIostreamCL& sub);
RecvStreamCL&    operator>> (RecvStreamCL&    is,       RefMPIistreamCL& sub);
RefMPIistreamCL& operator>> (RefMPIistreamCL& is,       RefMPIistreamCL& sub);
///@}

/// \brief input/output of C-strings and std::string.
/// @{
MPIostreamCL& write_char_array (MPIostreamCL& os,
    const MPIostreamCL::char_type* p, std::streamsize n);

inline MPIostreamCL& operator<< (MPIostreamCL& os, const          char* s);
inline MPIostreamCL& operator<< (MPIostreamCL& os, const signed   char* s);
inline MPIostreamCL& operator<< (MPIostreamCL& os, const unsigned char* s);
inline MPIostreamCL& operator<< (MPIostreamCL& os, const std::string& s);

inline MPIistreamCL& operator>> (MPIistreamCL& is, std::string& s);
///@}

/// \brief Writes the size of the array to follow.
/// In ascii-mode, the size is formatted to a fixed length with leading zeros
/// and ends with SendRecvStreamAsciiTerminatorC.
MPIostreamCL& write_array_header (MPIostreamCL& os, std::streamsize n);

/// \brief Append a sub-stream to an SendStreamCL inplace.
/// Equivalent to constructing the MPIostreamCL sub and insert it to an MPIostreamCL os via
/// os << sub. This class spares the copy of the temporary buffer sub by dircectly
/// writing to os.
class SubostreamBuilderCL
{
  private:
    SendStreamCL& os_; ///< Substreams are inserted into this stream.
    SendStreamCL::pos_type        begin_; ///< begin of the substream relative to its parents begin.
    SendStreamCL::pos_type        begin_data_; ///< begin of the data, i. e. after the size-header.

  public:
    SubostreamBuilderCL (SendStreamCL& os) : os_( os) {}
    void init () {
        begin_= os_.tellp();
        write_array_header( os_, 0);
        begin_data_= os_.tellp();
    }
    void append (SendStreamCL& s) {
        os_.write( s.begin(), s.cur() - s.begin());
    }
    void append (RefMPIostreamCL& s) {
        os_.write( s.begin(), s.cur() - s.begin());
    }
    void finalize () {
        const SendStreamCL::pos_type end= os_.tellp();
        os_.seekp( begin_);
        write_array_header( os_, end - begin_data_);
        os_.seekp( end);
    }
    void finalize_skip_empty () {
        if (os_.tellp() == begin_data_) {
            os_.seekp( begin_);
        }
        else
            finalize();
    }
};

} // end of namespace DiST
} // end of namespace DROPS

#include "DiST/mpistream.tpp"

#endif
