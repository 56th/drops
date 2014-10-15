/// \file mpistream.cpp
/// \brief SendStreamCL and RecStreamCL for transferring objects (simplexes + data) as byte-messages.
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

#include "DiST/mpistream.h"
#include "DiST/DiST.h"
#include "geom/simplex.h"
#include "geom/multigrid.h"

namespace DROPS{
namespace DiST{

bool use_binaryMPIstreams= true;


ProcCL::RequestT MPIstringbufCL::Isend (int dest, int tag)
{
    const std::streamsize size= pptr() - pbase();
    return ProcCL::Isend( pbase(), size, dest, tag);
}

void MPIstringbufCL::Recv (int source, int tag)
{
    if (!str().empty())
        throw ErrorCL( "MPIstringbufCL::Recv: Not cleared before reuse.\n");

    const int bufsize= ProcCL::GetMessageLength<char_type>( source, tag);
    str( std::string( bufsize, SendRecvStreamAsciiTerminatorC));
    ProcCL::Recv( eback(), bufsize, source, tag);
}


MPIrefbufCL* MPIrefbufCL::setbuf (char_type* b, std::streamsize s)
{
    if (mode_ & std::ios_base::in)
        setg( b, b, b + s);
    if (mode_ & std::ios_base::out)
        setp( b, b + s);
    return this;
}

MPIrefbufCL::pos_type MPIrefbufCL::seekoff (off_type off,
    std::ios_base::seekdir base, std::ios_base::openmode m)
{
    pos_type ret=  pos_type( off_type( -1));
    const bool testin=  (std::ios_base::in  & mode_ & m) != 0;
    const bool testout= (std::ios_base::out & mode_ & m) != 0;
    pos_type new_g_pos, new_p_pos;
    bool ok= true;
    if (testin) {
        if (base == std::ios_base::beg)
            new_g_pos= off;
        else if (base == std::ios_base::cur)
            new_g_pos=  (gptr() - eback()) + off;
        else
            new_g_pos= (egptr() - eback()) + off;
        char_type* const new_gptr= eback() + new_g_pos;
        if (new_gptr < eback() || new_gptr > egptr())
            ok= false;
    }
    if (testout) {
        if (base == std::ios_base::beg)
            new_p_pos= off;
        else if (base == std::ios_base::cur)
            new_p_pos= ( pptr() - pbase()) + off;
        else
            new_p_pos= (epptr() - pbase()) + off;
        char_type* const new_pptr= pbase() + new_p_pos;
        if (new_pptr < pbase() || new_pptr > epptr())
            ok= false;
    }
    if (ok) {
        if (testin) {
            ret= new_g_pos;
            setg( eback(), eback() + new_g_pos, egptr());
        }
        if (testout) {
            ret= new_p_pos;
            setp( pbase(), epptr());
            pbump( new_p_pos);
        }
    }
    return ret;
}

std::streamsize MPIrefbufCL::xsgetn (char_type* p, std::streamsize n)
{
    const std::streamsize m= std::min( n, egptr() - gptr());
    std::memcpy( p, gptr(), m*sizeof(char_type));
    gbump( m);
    // if ( gptr() == egptr())
    //     uflow();
    return m;
}

std::streamsize MPIrefbufCL::xsputnb(const char_type* p, std::streamsize n)
{
    const std::streamsize m= std::min( n, epptr() - pptr());
    std::memcpy( pptr(), p, m*sizeof(char_type));
    pbump( m);
    // if (pptr() == epptr())
    //     overflow();
    return m;
}

MPIostreamCL& operator<< (MPIostreamCL& os, const GeomIdCL& h)
{
    return os << h.level << h.bary << h.dim;
}

MPIistreamCL& operator>> (MPIistreamCL& is, GeomIdCL& h)
{
    return is >> h.level >> h.bary >> h.dim;
}


std::string terminated_array_header( std::streamsize n)
{
    const std::streamsize terminated_header_size= std::numeric_limits<std::streamsize>::digits10 + 2;
    std::ostringstream s;
    s.width( terminated_header_size - 1);
    s.fill( '0');
    s << n << SendRecvStreamAsciiTerminatorC;
    return s.str();
}

MPIostreamCL& write_array_header (MPIostreamCL& os, std::streamsize n)
{
    if (os.isBinary())
        os << n;
    else {
        const std::string& s= terminated_array_header( n);
        os.write( &s[0], s.size());
    }
    return os;
}

MPIostreamCL& write_char_array (MPIostreamCL& os,
    const MPIostreamCL::char_type* p, std::streamsize n)
{
    write_array_header( os, n);
    os.write( p, n);
    return os;
}

MPIostreamCL& operator<< (MPIostreamCL& os, const SendStreamCL& sub)
{
    return write_char_array( os, sub.begin(), sub.cur() - sub.begin());
}

MPIostreamCL& operator<< (MPIostreamCL& os, const RecvStreamCL& sub)
{
    return write_char_array( os, sub.begin(), sub.cur() - sub.begin());
}

MPIistreamCL& operator>> (MPIistreamCL& is, RecvStreamCL& sub)
{
    std::streamsize n;
    is >> n;
    sub.buf_.str( std::string( n, SendRecvStreamAsciiTerminatorC));
    is.read( sub.begin(), n);
    return is;
}

MPIostreamCL& operator<< (MPIostreamCL& os, const RefMPIostreamCL& sub)
{
    return write_char_array( os, sub.begin(), sub.cur() - sub.begin());
}


RecvStreamCL& operator>> (RecvStreamCL& is, RefMPIistreamCL& sub)
{
    std::streamsize n;
    is >> n;
    sub.setbuf( is.cur(), n);
    is.seekg( n, std::ios_base::cur);
    return is;
}

RefMPIistreamCL& operator>> (RefMPIistreamCL& is, RefMPIistreamCL& sub)
{
    std::streamsize n;
    is >> n;
    sub.setbuf( is.cur(), n);
    is.seekg( n, std::ios_base::cur);
    return is;
}

SendStreamCL::SendStreamCL (const SendStreamCL& s)
    : std::ios(), base_type( &buf_, binary), buf_( s.buf_.str(), std::ios_base::out)
{
    seekp( const_cast<SendStreamCL&>( s).tellp());
    copyfmt( s);
    clear( s.rdstate());
}

SendStreamCL& SendStreamCL::operator= (const SendStreamCL& s)
{
    if (&s == this)
        return *this;
    setBinary( s.isBinary());
    clear();
    clearbuffer( s.buf_.str());
    seekp( const_cast<SendStreamCL&>( s).tellp());
    copyfmt( s);
    clear( s.rdstate());
    return *this;
}

RecvStreamCL::RecvStreamCL (const RecvStreamCL& s)
    : std::ios(), base_type( &buf_, binary), buf_( s.buf_.str(), std::ios_base::in)
{
    seekg( const_cast<RecvStreamCL&>( s).tellg());
    copyfmt( s);
    clear( s.rdstate());

}

RecvStreamCL& RecvStreamCL::operator= (const RecvStreamCL& s)
{
    if (&s == this)
        return *this;
    setBinary( s.isBinary());
    clear();
    clearbuffer( s.buf_.str());
    seekg( const_cast<RecvStreamCL&>( s).tellg());
    copyfmt( s);
    clear( s.rdstate());
    return *this;
}


} // end of namespace DiST
} // end of namespace DROPS
