/// \file mpistream.cpp
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

#include "parallel/DiST.h"
#include "geom/simplex.h"
#include "geom/multigrid.h"

namespace DROPS{
namespace DiST{
namespace Helper{



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


SendStreamCL& operator<< (SendStreamCL& os, const GeomIdCL& h)
{
    os << h.level << h.bary[0] << h.bary[1] << h.bary[2] << h.dim;
    return os;
}

RecvStreamCL& operator>> (RecvStreamCL& is, GeomIdCL& h)
{
    is >> h.level >> h.bary[0] >> h.bary[1] >> h.bary[2] >> h.dim;
    return is;
}

SendStreamCL& operator<< (SendStreamCL& os, const Point3DCL& h)
{
    os << h[0] << h[1] << h[2];
    return os;
}

RecvStreamCL& operator>> ( RecvStreamCL& is, Point3DCL& h)
{
    is >> h[0] >> h[1] >> h[2];
    return is;
}

} // end of namespace Helper
} // end of namespace DiST
} // end of namespace DROPS
