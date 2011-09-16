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

/** This member function is just a mask of the already in DROPS implemented
    blocking receive, the only difference is that we always receive objects
    of datatype MPI_CHAR.
    \param source from whom should I receive
    \param tag    the tag which has been used to send the message
*/
void RecvStreamCL::Recv(int source, int tag)
{
    int bufsize = ProcCL::GetMessageLength<char>( source, tag);
    std::string temp(bufsize,' ');
    ProcCL::Recv(&temp[0], bufsize, source, tag);
    this->str( str() + temp);
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
