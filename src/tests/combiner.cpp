/// \file combiner.cpp
/// \brief tests implementation for the StrategyCombinerCL
/// \author LNM RWTH Aachen: Matthias Kirchhart; SC RWTH Aachen: 

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
 * Copyright 2014 LNM/SC RWTH Aachen, Germany
*/

#include <iostream>

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "levelset/adaptriang.h"
#include "out/output.h"
#include "out/vtkOut.h"

using namespace DROPS;

double pos1 = -0.5;
double pos2 =  0.5;

double lset1( const Point3DCL &x, double );
double lset2( const Point3DCL &x, double );

int main()
{
    try
    {
    const double L = 1; 
    Point3DCL orig(-L), e1, e2, e3;
    e1[0] = e2[1] = e3[2] = 2*L;

    const int n = 1;
    BrickBuilderCL builder( orig, e1, e2, e3, n, n, n );
    MultiGridCL MG( builder );

    StrategyCombinerCL    marker;
    DistMarkingStrategyCL marker1( lset1, 0.1, 0, 5 ); 
    DistMarkingStrategyCL marker2( lset2, 0.1, 0, 5 );

    marker.push_back( marker1 );
    marker.push_back( marker2 );


    AdapTriangCL adap( MG, &marker );
    adap.MakeInitialTriang();

    VTKOutCL vtk( MG, "combiner", 21, ".", "combiner", "combiner", true );
    vtk.Write( 0 );

    for ( int i = 1; i <= 20; ++i )
    {
        pos1 += 0.05;
        pos2 -= 0.05;
        adap.UpdateTriang();
        vtk.Write( i );
    }
    }
    catch ( DROPSErrCL err )
    {
        err.handle();
    }
}

inline
double lset1( const Point3DCL &x, double )
{
    return x[0] + pos1;
}

inline
double lset2( const Point3DCL &x, double )
{
    return x[0] + pos2;
}

