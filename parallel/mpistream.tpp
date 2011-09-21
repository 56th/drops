/// \file mpistream.tpp
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

namespace DROPS{
namespace DiST{
namespace Helper{

/// \brief operator << for SendStreamCL
template<typename T>
MPIostreamCL& operator<< (MPIostreamCL& os, const T& t)
{
    if (os.isBinary())
        os.write( reinterpret_cast<const char*>( &t), sizeof( T));
    else
        static_cast<MPIostreamCL::base_type&>( os) << t << SendRecvStreamAsciiTerminatorC;

    return os;
}

/// \brief operator >> for RecvStreamCL
template<typename T>
MPIistreamCL& operator>> (MPIistreamCL& is, T& t)
{
    if (is.isBinary())
        is.read( reinterpret_cast<char*>( &t), sizeof(T));
    else {
        MPIistreamCL::base_type& istr= static_cast<MPIistreamCL::base_type&>( is);
        istr >> t;
        MPIistreamCL::char_type c;
        istr.get( c);
        if (c != SendRecvStreamAsciiTerminatorC)
            throw DROPSErrCL( "MPIistreamCL& operator>>( MPIistreamCL& is, T& t):"
                "Ascii item-terminator not found.\n");
    }
    return is;
}


inline MPIostreamCL& operator<< (MPIostreamCL& os, const Point3DCL& h)
{
    os << h[0] << h[1] << h[2];
    return os;
}

inline MPIistreamCL& operator>> (MPIistreamCL& is, Point3DCL& h)
{
    is >> h[0] >> h[1] >> h[2];
    return is;
}


MPIostreamCL& write_char_array (MPIostreamCL& os,
    const MPIostreamCL::char_type* p, std::streamsize n);

inline MPIostreamCL& operator<< (MPIostreamCL& os, const char* s)
{
    return write_char_array( os, s, std::strlen( s));
}

inline MPIostreamCL& operator<< (MPIostreamCL& os, const signed char* s)
{
    const char* const p= reinterpret_cast<const char*>( s);
    return write_char_array( os, p, std::strlen( p));
}

inline MPIostreamCL& operator<< (MPIostreamCL& os, const unsigned char* s)
{
    const char* const p= reinterpret_cast<const char*>( s);
    return write_char_array( os, p, std::strlen( p));
}

inline MPIostreamCL& operator<< (MPIostreamCL& os, const std::string& s)
{
    return write_char_array( os, &s[0], s.size());
}

inline MPIistreamCL& operator>> (MPIistreamCL& is, std::string& s)
{
    std::streamsize n;
    is >> n;
    s.resize( 0);
    s.resize( n);
    is.read( &s[0], n);
    return is;
}


} // end of namespace Helper
} // end of namespace DiSt
} // end of namespace DROPS
