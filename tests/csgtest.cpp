/// \file csgtest.cpp
/// \brief test CSG for level sets.
/// \author LNM RWTH Aachen: Joerg Grande

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

#include <fstream>

#include "geom/csg.h"
#include "misc/utils.h"
#include "geom/multigrid.h"
#include "levelset/surfacetension.h"
#include "levelset/levelset.h"
#include "geom/builder.h"
#include "num/fe.h"
#include "misc/problem.h"
#include "out/vtkOut.h"


using namespace DROPS;

//dummy
double sigmaf (const Point3DCL&, double) { return 0.; }

const CSG::BodyCL* thebody;

inline double csg_fun (const Point3DCL& x, double t)
{
    return (*thebody)( x, t);
}

int Test (MultiGridCL& mg)
{
    std::ifstream jsonfile( "csgbody.json");
    ParamCL p;
    jsonfile >> p;
    thebody= CSG::body_builder( p);

    SurfaceTensionCL sf( sigmaf);   // dummy class
    LsetBndDataCL lsbnd( 6);
    LevelsetP2CL lset( mg, lsbnd, sf);
    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
    lset.Init( csg_fun);

    // writer for vtk-format
    VTKOutCL vtkwriter( mg, "DROPS data", 1,
                        ".", "nnn", 
                        "nnn", /* <- time file name */
                        true, 0, -1, 0);
    vtkwriter.Register( make_VTKScalar( lset.GetSolution(), "level-set") );
    vtkwriter.Write( 0.);

    delete thebody;
    return 0;
}

int main ()
{
  try {
    BrickBuilderCL mgb( Point3DCL( -1.), 2*std_basis<3>( 1), 2*std_basis<3>( 2), 2*std_basis<3>( 3), 8, 8, 8);
    MultiGridCL mg( mgb);
    return  Test( mg);
  }
  catch (DROPSErrCL err) { err.handle(); }
}
