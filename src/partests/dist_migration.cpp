/// \file dist_migration.cpp
/// \brief test of migration with unknowns
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
#include "stokes/instatstokes2phase.h"
#include "geom/geomselect.h"
#include "out/vtkOut.h"
#include "misc/dynamicload.h"
#include <iostream>
#include <sstream>

using namespace std;

namespace DROPS{

/// \brief Build a brick and tell parallel info class about the multigrid
void BuildBrick( MultiGridCL*& mg)
{
    Point3DCL origin, e1, e2, e3;
    e1[0]= e2[1]= e3[2]= 1.;
    Uint ref[3]= { 1, 1, 1};
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

void SetDOFLset( const MultiGridCL& mg, LevelsetP2CL& ls)
{
    const Uint idx = ls.Phi.RowIdx->GetIdx();
    VectorCL&  Data= ls.Phi.Data;

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, mg.GetLastLevel(), it){
        if ( it->Unknowns.Exist( idx))
            Data[ it->Unknowns(idx)]= (double)it->GetHashId()/1e20;
    }
    DROPS_FOR_TRIANG_CONST_EDGE( mg, mg.GetLastLevel(), it){
        if ( it->Unknowns.Exist( idx))
            Data[ it->Unknowns(idx)]= (double)it->GetHashId()/1e20;
    }
}


void SetDOFStokes( const MultiGridCL& mg, InstatStokes2PhaseP2P1CL& stokes)
{
    const Uint vel_idx = stokes.v.RowIdx->GetIdx();
    const Uint pr_idx  = stokes.p.RowIdx->GetIdx();
    VectorCL&  vel_data= stokes.v.Data;
    VectorCL&  pr_data = stokes.p.Data;

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, mg.GetLastLevel(), it){
        const double hash= it->GetHashId();
        if ( it->Unknowns.Exist( vel_idx)){
            const IdxT DOF= it->Unknowns( vel_idx);
            vel_data[ DOF]  = hash/2e20;
            vel_data[ DOF+1]= hash/3e20;
            vel_data[ DOF+2]= hash/4e20;
        }
        if ( it->Unknowns.Exist( pr_idx)){
            const IdxT DOF= it->Unknowns( pr_idx);
            pr_data[ DOF]  = hash/1e19;
        }
    }
    DROPS_FOR_TRIANG_CONST_EDGE( mg, mg.GetLastLevel(), it){
        const double hash= it->GetHashId();
        if ( it->Unknowns.Exist( vel_idx)){
            const IdxT DOF= it->Unknowns( vel_idx);
            vel_data[ DOF]  = hash/2e20;
            vel_data[ DOF+1]= hash/3e20;
            vel_data[ DOF+2]= hash/4e20;
        }
    }
}

void Migrate( LoadBalCL& lb, ObservedVectorsCL& obs)
{
    // Do the migration process (including the unknowns)
    obs.notify_pre_refmig_sequence( lb.GetMG());
    obs.notify_pre_migrate();
    lb.DoMigration();
    obs.notify_post_migrate();
    obs.notify_post_refmig_sequence();
}


void CheckDOFLset( const MultiGridCL& mg, LevelsetP2CL& ls)
{
    const Uint idx = ls.Phi.RowIdx->GetIdx();
    const VectorCL&  Data= ls.Phi.Data;

    std::vector<DiST::TransferableCL const *> fail_elems;
    // Check all vertices and edges for correct values
    DROPS_FOR_TRIANG_CONST_VERTEX( mg, mg.GetLastLevel(), it){
        if ( it->Unknowns.Exist( idx)){
            const IdxT DOF= it->Unknowns(idx);
            if ( std::abs(Data[DOF]-(double)it->GetHashId()/1e20)>DoubleEpsC){
                fail_elems.push_back( &*it);
            }
        }
    }
    DROPS_FOR_TRIANG_CONST_EDGE( mg, mg.GetLastLevel(), it){
        if ( it->Unknowns.Exist( idx)){
            const IdxT DOF= it->Unknowns(idx);
            if ( std::abs(Data[DOF]-(double)it->GetHashId()/1e20)>DoubleEpsC){
                fail_elems.push_back( &*it);
            }
        }
    }

    bool correct= ProcCL::Check( fail_elems.empty(), ProcCL::Master());
    if ( correct){
        std::cout << "Migration with unknowns is correct for level set function" << std::endl;
    }
    else{
        std::cout << "Migration with unknowns is NOT correct for level set function" << std::endl;
    }
}

void CheckDOFStokes( const MultiGridCL& mg, InstatStokes2PhaseP2P1CL& stokes)
{
    const Uint vel_idx = stokes.v.RowIdx->GetIdx();
    const Uint pr_idx  = stokes.p.RowIdx->GetIdx();
    VectorCL&  vel_data= stokes.v.Data;
    VectorCL&  pr_data = stokes.p.Data;
    typedef std::vector<DiST::TransferableCL const *> failContT;
    failContT fail_elems_vel, fail_elems_pr;
    // Check all vertices and edges for correct values
    DROPS_FOR_TRIANG_CONST_VERTEX( mg, mg.GetLastLevel(), it){
        const double hash= it->GetHashId();
        if ( it->Unknowns.Exist( vel_idx)){
            const IdxT DOF= it->Unknowns( vel_idx);
            if (   std::abs(vel_data[DOF]   - hash/2e20)>DoubleEpsC
                || std::abs(vel_data[DOF+1] - hash/3e20)>DoubleEpsC
                || std::abs(vel_data[DOF+2] - hash/4e20)>DoubleEpsC
                ){
                fail_elems_vel.push_back( &*it);
            }
        }
        if ( it->Unknowns.Exist( pr_idx)){
            const IdxT DOF= it->Unknowns( pr_idx);
            if ( std::abs(pr_data[DOF]-(double)it->GetHashId()/1e19)>DoubleEpsC){
                fail_elems_pr.push_back( &*it);
            }
        }
    }
    DROPS_FOR_TRIANG_CONST_EDGE( mg, mg.GetLastLevel(), it){
        if ( it->Unknowns.Exist( vel_idx)){
            const double hash= it->GetHashId();
            const IdxT DOF= it->Unknowns( vel_idx);
            if (   std::abs(vel_data[DOF]   - hash/2e20)>DoubleEpsC
                || std::abs(vel_data[DOF+1] - hash/3e20)>DoubleEpsC
                || std::abs(vel_data[DOF+2] - hash/4e20)>DoubleEpsC
                ){
                fail_elems_vel.push_back( &*it);
            }
        }
    }

    bool pr_correct= ProcCL::Check( fail_elems_pr.empty(), ProcCL::Master());
    bool vel_correct= ProcCL::Check( fail_elems_vel.empty(), ProcCL::Master());
    if ( pr_correct)
        std::cout << "Migration with unknowns is correct for pressure" << std::endl;
    else
        std::cout << "Migration with unknowns is NOT correct for pressure, failed on 0 in " << fail_elems_pr.size() << " cases\n";
    if ( vel_correct)
        std::cout << "Migration with unknowns is correct for velocity" << std::endl;
    else
        std::cout << "Migration with unknowns is NOT correct for velocity, failed on [0] in " << fail_elems_vel.size() << " cases\n";
}


double sigmaf (const Point3DCL&, double) { return 0; }

/** Assign a value to the DOF and migrate tetrahedra in a 'wild' order.
    Afterwards, check if the DOF are still right.
*/
void CheckMigration( LoadBalCL& lb)
{
    std::cout << "Checking migration with unknowns ..." << std::endl;

    std::vector<std::string> libs;
    libs.push_back("misc/libmisc-scalarFunctions");
    libs.push_back("misc/libmisc-vectorFunctions");
//    libs.push_back("levelset/liblevelset-twophaseCoeff");
    dynamicLoad( "../", libs);

    MultiGridCL& mg= lb.GetMG();

    // Create boundary conditions
    std::string perbndtypestr;
    std::string zerobndfun, zerovelbndfun;
    for( size_t i= 1; i<=mg.GetBnd().GetNumBndSeg(); ++i) {
        zerobndfun    += "Zero";
        zerovelbndfun += "ZeroVel";
        if (i!=mg.GetBnd().GetNumBndSeg()){
            zerobndfun   += "!";
            zerovelbndfun+= "!";
        }
    }
    match_fun periodic_match = 0;

    // Create all the stuff which is necessary to make one level set function
    LsetBndDataCL* lsetbnddata= 0;
    BuildBoundaryData( &mg, lsetbnddata, perbndtypestr, zerobndfun, periodic_match);
    SurfaceTensionCL sft( DROPS::sigmaf);
    LevelsetP2CL & lset( * LevelsetP2CL::Create( mg, *lsetbnddata, sft) );
    lset.CreateNumbering( mg.GetLastLevel());
    MLIdxDescCL* lidx= &lset.idx;


    // Create all the stuff which is necessary to make a Stokes class
    typedef BndDataCL<DROPS::Point3DCL> VelBndDataCL;
    typedef BndDataCL<double>    PrBndDataCL;
    VelBndDataCL *velbnddata = 0;
    PrBndDataCL *prbnddata = 0;
    BuildBoundaryData( &mg, prbnddata, perbndtypestr, zerobndfun, periodic_match);
    BuildBoundaryData( &mg, velbnddata, perbndtypestr, zerovelbndfun, periodic_match);
    StokesBndDataCL bnddata(*velbnddata,*prbnddata);
    TwoPhaseFlowCoeffCL tpf( 0, 0, 0, 0, 0, Point3DCL());
    InstatStokes2PhaseP2P1CL Stokes( mg, tpf, bnddata, DROPS::P1_FE, -0.1);
    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;
    Stokes.CreateNumberingVel( mg.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  mg.GetLastLevel(), pidx, &lset);
    Stokes.SetIdx();
    Stokes.v.SetIdx  ( vidx);
    Stokes.p.SetIdx  ( pidx);


    // make the repair classes and let them be observed
    LevelsetRepairCL lset_repair( lset);
    VelocityRepairCL velrepair( Stokes);
    PressureRepairCL prrepair( Stokes, lset);

    // let the level set, velocity and pressure FE function be observed
    ObservedVectorsCL obs;
    obs.push_back( &lset_repair);
    obs.push_back( &velrepair);
    obs.push_back( &prrepair);

    std::cout << "Set values for the pressure, the velocity and the level set function" << std::endl;
    SetDOFLset( mg, lset);
    SetDOFStokes( mg, Stokes);

    // do the migration
    std::cout << "Migrate the tetrahedra with the following DOF information\n"
              << " level set " << lidx->GetIdx() << '\n'
              << " velocity " << Stokes.v.RowIdx->GetIdx() << '\n'
              << " pressure " << Stokes.p.RowIdx->GetIdx() << '\n'
              << std::endl;
//    const DiST::GeomIdCL gid( 3, MakePoint3D(0.75, 0.25, 0.375), 0);
//    VertexCL *v= DiST::InfoCL::Instance().Exists(gid) ? DiST::InfoCL::Instance().GetVertex(gid) : 0;
//    IF_MASTER if (v) { v->Unknowns.DebugInfo(cdebug << gid << ": "); v->DebugInfo(cdebug); }
    Migrate( lb, obs);
    std::cout << "After migration:\n"
              << " level set " << lidx->GetIdx() << '\n'
              << " velocity " << Stokes.v.RowIdx->GetIdx() << '\n'
              << " pressure " << Stokes.p.RowIdx->GetIdx() << '\n'
              << std::endl;
//    v= DiST::InfoCL::Instance().Exists(gid) ? DiST::InfoCL::Instance().GetVertex(gid) : 0;
//    IF_MASTER if (v) { v->Unknowns.DebugInfo(cdebug << gid << ": "); v->DebugInfo(cdebug); }

    // Check for correctness
    CheckDOFLset( mg, lset);
    CheckDOFStokes( mg, Stokes);
    delete lsetbnddata; lsetbnddata=0;
    delete velbnddata; velbnddata=0;
    delete prbnddata; prbnddata=0;
    delete &lset;
    std::cout << "Check migration with unknowns performed" << std::endl;
}

}

int main( int argc, char **argv)
{
#ifdef _PAR
    DROPS::ProcCL::Instance(&argc, &argv);
#endif
    try {
        DROPS::MultiGridCL* mg= 0;
        DROPS::BuildBrick( mg);
        std::cout << "=====================================\ninitial migration\n";
        DROPS::LoadBalCL lb( *mg);  // loadbalancing
        lb.DoMigration( );          // distribute initial grid
        const int num_ref= 5;
        // writer for vtk-format
        DROPS::VTKOutCL vtkwriter( *mg, "dist_ref", num_ref+1, "vtk", "dist_mig", "dist_mig", true);
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

        DROPS::CheckMigration( lb);

        delete mg; mg=0;
    }
    catch (DROPS::DROPSErrCL err) {err.handle();}
    return 0;
}
