/// \file dist_ref.cpp
/// \brief test of refinement
/// \author LNM RWTH Aachen: Patrick Esser, Sven Gross; SC RWTH Aachen: Oliver Fortmeier

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

#include "parallel/parallel.h"
#include "misc/container.h"
#include "DiST/DiST.h"
#include "parallel/partime.h"
#include "parallel/loadbal.h"
#include "geom/simplex.h"
#include "geom/builder.h"
#include "geom/multigrid.h"
#include "misc/problem.h"
#include "levelset/levelset.h"
#include "out/vtkOut.h"
#include <iostream>
#include <sstream>

using namespace std;

namespace DROPS{

/// \brief Build a brick and tell parallel info class about the multigrid
void BuildBrick( MultiGridCL*& mg)
{
    Point3DCL origin, e1, e2, e3;
    e1[0]= e2[1]= e3[2]= 1.;
    Uint ref[3]= { 2, 2, 2};
    BrickBuilderCL builder( origin, e1, e2, e3, ref[0], ref[1], ref[2]);
    mg = new MultiGridCL( builder);
}

double DistToSphere( const Point3DCL& p, double)
{
    Point3DCL origin( 0.5);
    return (origin-p).norm() - 0.25;
}

void RefineBrick( MultiGridCL& mg)
{
//    MarkAll( mg);
    MarkInterface( DistToSphere, 0.1, mg);
    mg.Refine();
}

void CheckDiST( const MultiGridCL& mg, std::ostream& os)
{
    if ( DROPS::ProcCL::Check( DROPS::DiST::InfoCL::Instance().IsSane( os)))
        std::cout << " DiST-module seems to be alright!" << std::endl;
    else
        std::cout << " DiST-module seems to be broken!" << std::endl;
    if ( DROPS::ProcCL::Check( mg.IsSane( os)))
        std::cout << " multigrid seems to be alright!" << std::endl;
    else
        std::cout << " multigrid seems to be broken!" << std::endl;
}

}

int main( int argc, char **argv)
{
#ifdef _PAR
    DROPS::ProcCL::Instance(&argc, &argv);
#endif
    try {

/*
        DROPS::Point3DCL origin, e1, e2, e3;
        e1[0]= e2[1]= e3[2]= 1.;
        DROPS::BrickBuilderCL emptybrick( origin, e1, e2, e3);
        DROPS::MultiGridCL mg( brick);
        mg.PrepareModify();
        DROPS::VertexCL& vert=mg.GetSimplexFactory().MakeVertex( origin, 0);
        mg.FinalizeModify();
        DROPS::DiST::ModifyCL modify( mg);
        modify.Init();
        modify.Identify(vert, vert.GetPrio());
        modify.Finalize();

        std::cout << "Vertex is local: " << vert.IsLocal() << std::endl;
        return 0;
*/

        DROPS::MultiGridCL* mg= 0;
        DROPS::BuildBrick( mg);
        std::cout << "=====================================\ninitial migration\n";
        DROPS::LoadBalCL lb( *mg);
        lb.DoMigration();
        //DROPS::LoadBalHandlerCL lb( *mg, DROPS::metis);     // loadbalancing
        //lb.DoInitDistribution( DROPS::ProcCL::Master());    // distribute initial grid
        const int num_ref= 5;
        // writer for vtk-format
        DROPS::VTKOutCL vtkwriter( *mg, "dist_ref", num_ref+1, "vtk", "dist_ref", "dist_ref", true);
        vtkwriter.Write(0);
        // refinement
        std::string filename("sane.chk");
        DROPS::ProcCL::AppendProcNum(filename);
        std::ofstream sanityfile( filename.c_str());
        for (int i=0; i<num_ref; ++i) {
            if (i!=0) {
                std::cout << "=====================================\nmigration " << i+1 << "\n";
                sanityfile<< "=====================================\nmigration " << i+1 << "\n";
                lb.DoMigration();
            }
            mg->SizeInfo( std::cout);
            DROPS::CheckDiST( *mg, sanityfile);
            std::cout << "=====================================\nrefinement " << i+1 << std::endl;
            sanityfile<< "=====================================\nrefinement " << i+1 << std::endl;
            DROPS::RefineBrick( *mg);
            mg->SizeInfo( std::cout);
            DROPS::CheckDiST( *mg, sanityfile);
            vtkwriter.Write(i+1);
        }

        sanityfile<< "=====================================\nDiST debug info" << std::endl;
        const DROPS::DiST::InfoCL& info= DROPS::DiST::InfoCL::Instance();
        info.SizeInfo(sanityfile);
        info.Instance().DebugInfo(sanityfile);
        sanityfile<< "=====================================\nmultigrid debug info" << std::endl;
        mg->SizeInfo(sanityfile);
        mg->DebugInfo(sanityfile);
        delete mg; mg=0;
    }
    catch (DROPS::DROPSErrCL err) {err.handle();}
    return 0;
}
