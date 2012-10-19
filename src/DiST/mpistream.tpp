/// \file mpistream.tpp
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

namespace DROPS{
namespace DiST{

template<typename T>
inline MPIostreamCL& MPIostreamCL::write_fundamental_type (const T& t)
{
    if (isBinary())
        this->write( reinterpret_cast<const char*>( &t), sizeof( T));
    else
        static_cast<MPIostreamCL::base_type&>( *this) << t << SendRecvStreamAsciiTerminatorC;

    return *this;
}

/// \brief operator >> for RecvStreamCL
template<typename T>
inline MPIistreamCL& MPIistreamCL::read_fundamental_type (T& t)
{
    if (isBinary())
        this->read( reinterpret_cast<char*>( &t), sizeof(T));
    else {
        MPIistreamCL::base_type& istr= static_cast<MPIistreamCL::base_type&>( *this);
        istr >> t;
        MPIistreamCL::char_type c= SendRecvStreamAsciiTerminatorC; // Initialisation: If the stream is not good(), nothing is read. Do not fail due to not reading at all. For extra credit: Check, if !good() is due to the previous read.
        istr.get( c);
        if (c != SendRecvStreamAsciiTerminatorC)
            throw DROPSErrCL( "MPIistreamCL& operator>>( MPIistreamCL& is, T& t):"
                "Ascii item-terminator not found.\n");
    }
    return *this;
}

template <Uint rows>
  inline MPIostreamCL&
  operator<< (MPIostreamCL& os, const SVectorCL<rows>& p)
{
    if (os.isBinary())
        os.write( reinterpret_cast<const char*>( &p), sizeof( SVectorCL<rows>));
    else
        for (Uint i= 0; i < rows; ++i)
            os << p[i];
    return os;
}

template <Uint rows>
  inline MPIistreamCL&
  operator>> (MPIistreamCL& is, SVectorCL<rows>& p)
{
    if (is.isBinary())
        is.read( reinterpret_cast<char*>( &p), sizeof( SVectorCL<rows>));
    else
        for (Uint i= 0; i < rows; ++i)
            is >> p[i];
    return is;
}


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


} // end of namespace DiST
} // end of namespace DROPS
