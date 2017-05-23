/// \file show_domain
/// \brief Read a domain and a VTK-writer from the command line and output the domain and its boundary segments.
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
 * Copyright 2013 LNM, Germany
*/

#include "geom/multigrid.h"
#include "out/output.h"
#include "out/vtkOut.h"
#include "geom/builder.h"
#include "misc/params.h"
#include "misc/dynamicload.h"

// #include "misc/funcmap.h"
//
// bool mymatch (const DROPS::Point3DCL&, const DROPS::Point3DCL&)
// { return true; }
//
// static DROPS::RegisterMatchingFunction regmymatch("mymatch",  mymatch);


int main (int argc, char** argv)
{
  try
  {

    DROPS::ParamCL P;

    DROPS::read_parameter_file_from_cmdline( P, argc, argv, "../../param/tests/show_domain/show_domain.json");
    P.put_if_unset<std::string>( "VTK.TimeFileName", P.get<std::string>("VTK.VTKName"));
    std::cout << P << std::endl;

    DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"), P.get<std::vector<std::string> >("General.DynamicLibs") );

    std::unique_ptr<DROPS::MGBuilderCL> builder( DROPS::make_MGBuilder( P) );
    DROPS::MultiGridCL mg( *builder);
    const DROPS::ParamCL::ptree_type* ch= 0;
    try {
        ch= &P.get_child( "Mesh.Periodicity");
    }
    catch (DROPS::DROPSParamErrCL) {}
    if (ch)
        read_PeriodicBoundaries( mg, *ch);
    std::cout << DROPS::SanityMGOutCL(mg) << std::endl;

    DROPS::IdxDescCL p2idx( DROPS::P2_FE);
    const DROPS::Uint lvl= mg.GetLastLevel();
    p2idx.CreateNumbering( lvl, mg);
    DROPS::VecDescCL bndvec( &p2idx);
    DROPS::BndDataCL<> nobnd( 0);

    DROPS::VTKOutCL vtk(
        /*mg*/          mg,
        /*dataname*/    "DROPS data",
        /*numsteps*/    1,
        /*dirname*/     P.get<std::string>("VTK.VTKDir"),
        /*filename*/    P.get<std::string>("VTK.VTKName"),
        /*pvdfilename*/ P.get<std::string>("VTK.TimeFileName"),
        /*binary*/      P.get<int>("VTK.Binary"),
        /*onlyP1*/      P.get<int>("VTK.UseOnlyP1"),
        /*P2DG*/        false,
        /*level*/       mg.GetLastLevel(),
        /*reusepvd*/    P.get<int>("VTK.ReUseTimeFile"),
        /*usedeformed*/ P.get<int>("VTK.UseDeformation"));
    vtk.Register( make_VTKScalar( make_P2Eval( mg, nobnd, bndvec), "boundary_segments"));

    bndvec.Data= -1.;
    const DROPS::Uint idx= p2idx.GetIdx();
    // Write the smallest boundary-segment number into bndvec.
    DROPS_FOR_TRIANG_VERTEX( mg, lvl, it)
        if (it->IsOnBoundary())
            bndvec.Data[it->Unknowns( idx)]= it->GetBndVertBegin()->GetBndIdx();
    DROPS_FOR_TRIANG_EDGE( mg, lvl, it)
        if (it->IsOnBoundary())
            bndvec.Data[it->Unknowns( idx)]= *it->GetBndIdxBegin();

    vtk.Write( 0);
    return 0;
  }
  catch (DROPS::DROPSErrCL& err) { err.handle(); }
}

