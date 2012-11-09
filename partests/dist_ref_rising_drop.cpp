/// \file dist_ref_rising_drop.cpp
/// \brief test of refinement when simulating a rising drop. No fluid dynamics are calculated,
///        only the refinements are performed
/// \author LNM RWTH Aachen:; SC RWTH Aachen: Oliver Fortmeier

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
#include "misc/params.h"
#include "levelset/levelset.h"
#include "out/vtkOut.h"
#include <iostream>
#include <sstream>

DROPS::ParamCL P;

namespace DROPS{

/// \brief Build a brick specified by Brick.Mesh of the parameter file
void BuildBrick( MultiGridCL*& mg)
{
    int nx, ny, nz;
    double dx, dy, dz;
    std::string mesh( P.get<string>("Brick.Mesh","1x10x1@4x4x4")), delim("x@");
    size_t idx;
    while ((idx= mesh.find_first_of( delim)) != std::string::npos )
        mesh[idx]= ' ';
    std::istringstream brick_info( mesh);
    brick_info >> dx >> dy >> dz >> nx >> ny >> nz;
    Point3DCL orig, px, py, pz;
    px[0]= dx; py[1]= dy; pz[2]= dz;
    MGBuilderCL* builder=0;
    if ( ProcCL::IamMaster()){
        builder = new BrickBuilderCL( orig, px, py, pz, nx, ny, nz);
    }
    else{
        builder = new DROPS::EmptyBrickBuilderCL( orig, px, py, pz);
    }
    mg = new MultiGridCL( *builder);
    delete builder;
}


/// \name Global variables to determine the distance to the interface
//@{
Point3DCL origin;
Point3DCL drop_vel;
double radius;
//@}

double DistToSphere( const Point3DCL& p, double)
{
    return (origin-p).norm() - 0.25;
}

void RefineBrick( MultiGridCL& mg)
{
    const int f_level= P.get<int>("Exp.FinestLevel", 2),        // finest level
              c_level= P.get<int>("Exp.CoarsestLevel", 0);      // coarsest level
    const int min_ref_num= f_level - c_level;                   // minimal number of refinements
    bool marked= true;
    for (int i= 0; i < 2*min_ref_num && marked; ++i) {
        marked= MarkInterface( DistToSphere, P.get<double>("Exp.RefWidth",0.1), mg, (Uint)f_level, (Uint)c_level);
        std::cout << "Mark interface for " << i << '/' << (2*min_ref_num) << " time, marked: " << marked << std::endl;
        mg.Refine();
    }
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
    DROPS::ProcInitCL procinit(&argc, &argv);
#endif
    try {
        std::ifstream param;
        if (argc!=2)
        {
            std::cout << "Using default parameter file: risingdroplet.json\n";
            param.open( "risingdroplet.json");
        }
        else{
            std::cout << "Opening file " << argv[1] << "\n";
            param.open( argv[1]);
        }
        if (!param)
        {
            std::cerr << "error while opening parameter file\n";
            return 1;
        }
        param >> P;
        param.close();
        std::cout << P << std::endl;


        DROPS::MultiGridCL* mg= 0;
        DROPS::BuildBrick( mg);
        std::cout << "=====================================\ninitial migration\n";
        DROPS::LoadBalCL lb( *mg);  // loadbalancing
        lb.DoMigration( );          // distribute initial grid

        // writer for vtk-format
        DROPS::VTKOutCL *vtkwriter=0;
        if (P.get("VTK.VTKOut", 0)!=0){
            vtkwriter = new DROPS::VTKOutCL(*mg, "DROPS data", P.get<int>("Time.NumSteps")/P.get("VTK.VTKOut", 0)+1,
                                P.get<std::string>("VTK.VTKDir"), P.get<std::string>("VTK.VTKName"), P.get<std::string>("VTK.VTKName"),
                                P.get<int>("VTK.Binary"));
            vtkwriter->Write(0);
        }

        // file for sanity
        std::ofstream* sanityfile=0;
        if ( P.get<int>("Exp.CheckSanity", 1)!=0){
            std::string filename("sane.chk");
            DROPS::ProcCL::AppendProcNum(filename);
            sanityfile= new std::ofstream( filename.c_str());
        }
        // file for timings
        std::ostringstream oss;
        oss << "timings" << DROPS::ProcCL::Size() << ".txt";
        std::ofstream timings( oss.str().c_str());

        // Set initial drop
        DROPS::origin  = P.get<DROPS::Point3DCL>("Exp.Origin");
        DROPS::drop_vel= P.get<DROPS::Point3DCL>("Exp.DropVel");
        DROPS::radius  = P.get<double>(          "Exp.DropRadius", 0.1);;

        DROPS::ParTimerCL timer;

        for (int i=0; i<P.get<int>("Exp.NumSteps", 10); ++i) {
            timer.Reset();
            if (i!=0) {
                std::cout << "=====================================\n- migration " << i+1 << "\n";
                if ( P.get<int>("Exp.CheckSanity", 1)!=0)
                    (*sanityfile)<< "=====================================\nmigration " << i+1 << "\n";
                lb.DoMigration();
            }
            mg->SizeInfo( std::cout);
            if ( P.get<int>("Exp.CheckSanity", 1)!=0){
                DROPS::CheckDiST( *mg, *sanityfile);
                (*sanityfile)<< "=====================================\nrefinement " << i+1 << std::endl;
            }
            std::cout << "- refinement " << i+1 << std::endl;
            DROPS::RefineBrick( *mg);
            DROPS::origin += DROPS::drop_vel;
            mg->SizeInfo( std::cout);
            std::cout << "MG has " << mg->GetNumLevel() << " levels\n";
            if ( P.get<int>("Exp.CheckSanity", 1)!=0)
                DROPS::CheckDiST( *mg, *sanityfile);
            if (P.get("VTK.VTKOut", 0)!=0){
                vtkwriter->Write(i+1);
            }
            timer.Stop();
            std::cout << "====> Step " << i << " took " << timer.GetTime() << " second.\n";
            timings << i << ' ' << timer.GetTime() << std::endl;;
        }

        if ( P.get<int>("Exp.CheckSanity", 1)!=0){
            (*sanityfile)<< "=====================================\nDiST debug info" << std::endl;
            const DROPS::DiST::InfoCL& info= DROPS::DiST::InfoCL::Instance();
            info.SizeInfo( *sanityfile);
            info.Instance().DebugInfo( *sanityfile);
            (*sanityfile) << "=====================================\nmultigrid debug info" << std::endl;
            mg->SizeInfo(*sanityfile);
            mg->DebugInfo(*sanityfile);
        }
        mg->SizeInfo( std::cout);
        delete mg; mg=0;
        if (sanityfile) delete sanityfile;
        if (vtkwriter) delete vtkwriter;
        timings << std::endl;
    }
    catch (DROPS::DROPSErrCL& err) {err.handle();}
    return 0;
}
