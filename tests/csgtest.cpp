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
#include "misc/funcmap.h"
#include "misc/utils.h"
#include "geom/multigrid.h"
#include "levelset/surfacetension.h"
#include "levelset/levelset.h"
#include "geom/builder.h"
#include "num/fe.h"
#include "misc/problem.h"
#include "out/vtkOut.h"


using namespace DROPS;

/// The level set is a heart with an inward and an outward cusp on the p[2]-axis; contained in \f$(-1.2,1.2)\times(-0.7,0.7)\times(-1.1,1.3)\f$.
double suess (const Point3DCL& p, double)
{
    return std::pow( std::pow( p[0], 2) + 2.25*std::pow( p[1], 2) + std::pow( p[2], 2) - 1., 3) - std::pow( p[0], 2)*std::pow( p[2], 3) - 9./80.*std::pow( p[1], 2)*std::pow( p[2], 3);
}
RegisterScalarFunction reg_suess( "suess", &suess);

//dummy
double sigmaf (const Point3DCL&, double) { return 0.; }

const CSG::BodyCL* thebody;

inline double csg_fun (const Point3DCL& x, double t)
{
    return (*thebody)( x, t);
}

int TestExamples (MultiGridCL& mg)
{
    SurfaceTensionCL sf( sigmaf);   // dummy class
    LsetBndDataCL lsbnd( 6);
    LevelsetP2CL lset( mg, lsbnd, sf);
    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);

    std::ifstream jsonfile( "../geom/csg-examples.json");
    ParamCL p;
    jsonfile >> p;
    const size_t num= std::distance( p.begin(), p.end());

   // writer for vtk-format
    VTKOutCL vtkwriter( mg, "DROPS data", num,
                        ".", "csg-examples", 
                        "csg-examples", /* <- time file name */
                        true, 0, -1, 0);
    vtkwriter.Register( make_VTKScalar( lset.GetSolution(), "level-set") );

    size_t i= 0;
    for (ParamCL::ptree_const_iterator_type it= p.begin(); it != p.end(); ++it, ++i) {
        std::cout << "Processing example \"" << it->first <<"\"." << std::endl;
        thebody= CSG::body_builder( it->second);
        lset.Init( csg_fun);
        vtkwriter.Write( i);
        delete thebody;
    }
    return 0;
}

int main ()
{
  try {
    BrickBuilderCL mgb( Point3DCL( -1.), 2*std_basis<3>( 1), 2*std_basis<3>( 2), 2*std_basis<3>( 3), 16, 16, 16);
    MultiGridCL mg( mgb);
    return  TestExamples( mg);
  }
  catch (DROPSErrCL err) { err.handle(); }
}
