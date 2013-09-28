/// \file fe.cpp
/// \brief Description of various finite-element functions
/// \author LNM RWTH Aachen: Joer Grande, Sven Gross, Volker Reichelt ; SC RWTH Aachen:

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

#include "num/fe.h"

namespace DROPS
{

const double FE_P1CL::_gradient[4][3]=
    { {-1., -1., -1.}, {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.} };

const double FE_P1DCL::_gradient[4][3]=
    { {3., 3., 3.}, {-3., 0., 0.}, {0., -3., 0.}, {0., 0., -3.} };

const double FE_P1DCL::_vertexvalue[4][4]=
    { {-2., 1., 1., 1.}, {1., -2., 1., 1.}, {1., 1., -2., 1.}, {1., 1., 1., -2.} };

const double FE_P2CL::_D2H[10][3][3]= {
    { {4., 4., 4.}, {4., 4., 4.}, {4., 4., 4.} },
    { {4., 0., 0.}, {0., 0., 0.}, {0., 0., 0.} },
    { {0., 0., 0.}, {0., 4., 0.}, {0., 0., 0.} },
    { {0., 0., 0.}, {0., 0., 0.}, {0., 0., 4.} },
    { {-8., -4., -4.}, {-4., 0., 0.}, {-4., 0., 0.} },
    { {0., -4., 0.}, {-4., -8., -4.}, {0., -4., 0.} },
    { {0., 4., 0.}, {4., 0., 0.}, {0., 0., 0.} },
    { {0., 0., -4.}, {0., 0., -4.}, {-4., -4., -8.} },
    { {0., 0., 4.}, {0., 0., 0.}, {4., 0., 0.} },
    { {0., 0., 0.}, {0., 0., 4.}, {0., 4., 0.} } };

void
FE_P2CL::ApplyAll(Uint numpt, const BaryCoordCL* const pt, std::valarray<double>* v)
{
    for (Uint i= 0; i < 10; ++i)
        v[i].resize( numpt);
    for (Uint i= 0; i < numpt; ++i) {
        v[0][i]= H0( pt[i]);
        v[1][i]= H1( pt[i]);
        v[2][i]= H2( pt[i]);
        v[3][i]= H3( pt[i]);
        v[4][i]= H4( pt[i]);
        v[5][i]= H5( pt[i]);
        v[6][i]= H6( pt[i]);
        v[7][i]= H7( pt[i]);
        v[8][i]= H8( pt[i]);
        v[9][i]= H9( pt[i]);
    }
}

} // end of namespace DROPS
