/// \file triang.cpp
/// \brief tests iterators on triangulation
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Volker Reichelt; SC RWTH Aachen:

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

#include "misc/utils.h"
#include "geom/multigrid.h"
#include "geom/builder.h"

using namespace DROPS;

enum  OutputModeT { SILENT, NOISY };


void MarkDrop(DROPS::MultiGridCL& mg, int maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;
    DROPS_FOR_TRIANG_TETRA( mg, maxLevel, It) {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.2,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRegRefMark();
    }
}

void UnMarkDrop(DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte( 0.5);
    DROPS_FOR_TRIANG_TETRA( mg, maxLevel, It) {
        if ( (GetBaryCenter( *It)-Mitte).norm() <= std::max( 0.2, 1.5*std::pow( It->GetVolume(), 1.0/3.0)) )
            It->SetRemoveMark();
    }
}


Uint
Old( DROPS::MultiGridCL& mg)
{
    std::cout << "\n-----------------------------------------------------------------"
                 "\niterator:\n";

    Uint tmp= 0;
    DROPS_FOR_TRIANG_VERTEX( mg, /*default level*/ -1, it) {
        tmp+= it->GetId().GetIdent();
    }
    DROPS_FOR_TRIANG_EDGE( mg, /*default level*/ -1, it) {
        ++tmp;
    }
    DROPS_FOR_TRIANG_FACE( mg, /*default level*/-1, it) {
        ++tmp;
    }
    DROPS_FOR_TRIANG_TETRA( mg, /*default level*/-1, it) {
        tmp+= it->GetId().GetIdent();
    }
    return tmp;
}

Uint
New(DROPS::MultiGridCL&, const TriangVertexCL& vt, const TriangEdgeCL& et,
    const TriangFaceCL& ft, const TriangTetraCL& tt)
{
    std::cout << "\n-----------------------------------------------------------------"
                 "\nptr_iterator and tmp-variable:\n";

    Uint tmp= 0;
    for (TriangVertexCL::const_ptr_iterator sit__( vt.begin()),
         end( vt.end()); sit__ != end; ++sit__) {
        const VertexCL* const sit( *sit__);
        tmp+= sit->GetId().GetIdent();
    }
    for (TriangEdgeCL::const_ptr_iterator sit( et.begin()),
         end( et.end()); sit != end; ++sit) {
        ++tmp;
    }
    for (TriangFaceCL::const_ptr_iterator sit( ft.begin()),
         end( ft.end()); sit != end; ++sit) {
        ++tmp;
    }
    for (TriangTetraCL::const_ptr_iterator sit__( tt.begin()),
         end( tt.end()); sit__ != end; ++sit__) {
        const TetraCL* const sit( *sit__);
        tmp+= sit->GetId().GetIdent();
    }
    for (TriangTetraCL::const_ptr_iterator sit__( tt.begin()),
         end( tt.end()); sit__ != end; ++sit__) {
        const TetraCL* const sit( *sit__);
        tmp+= sit->GetId().GetIdent();
    }
    return tmp;
}


int main ()
{
  try {
    DROPS::BrickBuilderCL brick(DROPS::std_basis<3>(0),
                                DROPS::std_basis<3>(1),
                                DROPS::std_basis<3>(2),
                                DROPS::std_basis<3>(3),
                                30, 30, 30);
    DROPS::MultiGridCL mg( brick);
    mg.SizeInfo( std::cout);
    MarkDrop( mg, -1);
    mg.Refine();
    MarkDrop( mg, -1);
    mg.Refine();
    MarkDrop( mg, -1);
    mg.Refine();
    mg.SizeInfo( std::cout);

    TriangVertexCL vt( mg);
    TriangEdgeCL et( mg);
    TriangFaceCL ft( mg);
    TriangTetraCL tt( mg);

    TimerCL time;
    time.Start();
    Uint q0= Old( mg);
    time.Stop();
    std::cout << "value1: " << q0 << "\ttime: " << time.GetTime() << " seconds"
        << std::endl;
    time.Reset();
    time.Start();
    Uint q3= Old( mg);
    time.Stop();
    std::cout << "value2: " << q3 << "\ttime: " << time.GetTime() << " seconds"
        << std::endl;
    time.Reset();
    time.Start();
    Uint q1= New( mg, vt, et, ft, tt);
    time.Stop();
    std::cout << "value3: " << q1 << "\ttime: " << time.GetTime() << " seconds"
        << std::endl;
    time.Reset();
    time.Start();
    Uint q2= New( mg, vt, et, ft, tt);
    time.Stop();
    std::cout << "value4: " << q2 << "\ttime: " << time.GetTime() << " seconds"
        << std::endl;
    time.Reset();
    std::cout << "tt.size: " << tt.size() << std::endl;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
