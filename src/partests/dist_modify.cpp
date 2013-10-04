/// \file dist_modify.cpp
/// \brief test remote data for DiST
/// \author LNM RWTH Aachen: ; SC RWTH Aachen: Oliver Fortmeier, Daniel Medina Cardona

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
#include "geom/simplex.h"
#include "geom/builder.h"
#include "geom/multigrid.h"
#include <iostream>
#include <sstream>

using namespace std;

namespace DROPS{

/// \brief Build a brick and tell parallel info class about the multigrid
void BuildBrick( MultiGridCL*& mg)
{
    Point3DCL origin, e1, e2, e3;
    e1[0]= e2[1]= e3[2]= 1.;
    Uint ref[3]= { 4, 4, 4};
    BrickBuilderCL builder( origin, e1, e2, e3, ref[0], ref[1], ref[2]);
    mg = new MultiGridCL( builder);
}


/// \brief Process sender sends its i-th tetrahedron to process sender+i mod P
void SendTetras( MultiGridCL& mg, int sender, const bool binary=true)
{
    const int P= ProcCL::Size();
    sender= sender%P;
    std::cout << "Proc " << sender << " is sending its tetras around..." << std::endl;
    DiST::TransferCL transfer( mg, true, binary);
    transfer.Init();
    if ( ProcCL::MyRank()==sender){
        int            toProc= (sender+1)%P;
        const Priority prio  = PrioMaster;
        const bool     todel = true;
        for ( MultiGridCL::const_TetraIterator sit( mg.GetTetrasBegin()); sit!=mg.GetTetrasEnd(); ++sit) {
            transfer.Transfer( *sit, toProc, prio, todel);
            toProc= (toProc+1)%P;
        }
    }
    transfer.Finalize();
}

/// \brief Change prios for all vertices
void ChangePrios( MultiGridCL& mg, int sender, const bool binary=true)
{
    const int me= ProcCL::MyRank();
    const Priority prio= static_cast<Priority>( PrioMaster + (me+sender)%ProcCL::Size());
    std::cerr << "Proc " << me << " is changing its verts' priority to " << PriorityToString(prio) << std::endl;
    DiST::ModifyCL mod( mg, true, binary);
    mod.Init();
    for ( MultiGridCL::const_VertexIterator sit( mg.GetVerticesBegin()); sit!=mg.GetVerticesEnd(); ++sit)
        mod.ChangePrio( *sit, prio);
    mod.Finalize();
}

} // end of namespace DROPS

int main( int argc, char **argv)
{
#ifdef _PAR
    DROPS::ProcCL::Instance(&argc, &argv);
#endif
    try {
        bool binary=true;
        if ( argc==2){
            const int arg_binary=atoi(argv[1]);
            if ( arg_binary==0)
                binary=false;
        }
        else{
            std::cout << "usage: " << argv[0] << " 0 : for transferring information ASCII based\n"
                      << "usage: " << argv[0] << " 1 : for transferring information binary based (default)"
                      << std::endl;
        }
        DROPS::MultiGridCL* mg= 0;
        DROPS::BuildBrick( mg);
        const DROPS::DiST::InfoCL& info= DROPS::DiST::InfoCL::Instance();
        std::string filename("sane.chk");
        DROPS::ProcCL::AppendProcNum(filename);
        std::ofstream sanityfile( filename.c_str());

        for (int sender=0; sender < DROPS::ProcCL::Size(); ++sender) {
            DROPS::SendTetras( *mg, sender, binary);
            DROPS::ChangePrios( *mg, sender, binary);
            mg->SizeInfo( std::cout);
            if ( DROPS::ProcCL::Check( info.IsSane( sanityfile))){
                std::cout << " DiST-module seems to be alright!" << std::endl;
            }
            else{
                std::cout << " DiST-module seems to be broken!" << std::endl;
            }
            sanityfile << "\n\n===========================\n\n";
        }
        info.SizeInfo(sanityfile);
        info.Instance().DebugInfo(sanityfile);
    }
    catch (DROPS::DROPSErrCL err) {err.handle();}
    return 0;
}
