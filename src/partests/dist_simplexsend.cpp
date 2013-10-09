/// \file dist_simplexsend.cpp
/// \brief test streaming for DiST
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
#include <iostream>
#include <sstream>

using namespace std;
using namespace DROPS;

void CheckDiST( std::ostream& os)
{
    if ( DROPS::ProcCL::Check( DROPS::DiST::InfoCL::Instance().IsSane( os)))
        std::cout << " DiST-module seems to be alright!" << std::endl;
    else
        std::cout << " DiST-module seems to be broken!" << std::endl;
}

int main( int argc, char **argv)
{
#ifdef _PAR
    DROPS::ProcCL::Instance(&argc, &argv);
#endif
    try {
        if (DROPS::ProcCL::Size()!=3){
            std::cout <<  "usage: " << argv[0] << " with 3 processes" << std::endl;
            throw DROPS::DROPSErrCL("Start program with 3 processes");
        }

        TetraBuilderCL *builder=0;
        MultiGridCL *mg= 0;

        // Process 0 constructs some simplices and transfers these to processes 1 and 2
        if ( DROPS::ProcCL::IamMaster()){
            builder= new TetraBuilderCL(0);
            mg = new MultiGridCL( *builder);
            DiST::TransferCL transfer( *mg, true, true);

            transfer.Init();
            transfer.Transfer( *mg->GetTetrasBegin(), 1, PrioMaster, false);
            transfer.Transfer( *mg->GetTetrasBegin(), 2, PrioMaster, false);
            transfer.Finalize();
        }
        else{
            builder= new TetraBuilderCL(0);
            mg = new MultiGridCL( *builder);
            DiST::TransferCL transfer( *mg, true, true);

            transfer.Init();
            transfer.Finalize();

            for ( MultiGridCL::const_VertexIterator sit(mg->GetVerticesBegin()); sit!=mg->GetVerticesEnd(); ++sit){
                std::cerr << ProcCL::MyRank() << ": Received Vertex " << sit->GetGID() << '\n';
            }
            for ( MultiGridCL::const_EdgeIterator sit(mg->GetEdgesBegin()); sit!=mg->GetEdgesEnd(); ++sit){
                std::cerr << ProcCL::MyRank() << ": Received Edge " << sit->GetGID() << '\n'
                          << " with Vertices " << sit->GetVertex(0)->GetGID()
                          << " and " <<sit->GetVertex(1)->GetGID()
                          << '\n';

            }
            for ( MultiGridCL::const_FaceIterator sit(mg->GetFacesBegin()); sit!=mg->GetFacesEnd(); ++sit){
                std::cerr << ProcCL::MyRank() << ": Received Face " << sit->GetGID() << '\n';
            }
            for ( MultiGridCL::const_TetraIterator sit(mg->GetTetrasBegin()); sit!=mg->GetTetrasEnd(); ++sit){
                std::cerr << ProcCL::MyRank() << ": Received Tetra " << sit->GetGID() << '\n';
                for (Uint i=0; i<NumVertsC; ++i){
                    std::cerr << " vertex(" << i << "): " << sit->GetVertex(i)->GetGID() << '\n';
                }
                for (Uint i=0; i<NumEdgesC; ++i){
                    std::cerr << " edge(" << i << "): " << sit->GetEdge(i)->GetGID() << '\n';
                }
                for (Uint i=0; i<NumFacesC; ++i){
                    std::cerr << " faces(" << i << "): " << sit->GetFace(i)->GetGID() << '\n';
                    std::cerr << "  Linked to: tetras:";
                    for ( size_t i=0; i<4; ++i){
                        if ( sit->GetFace(i)->GetNeighbor(i)==0){
                            std::cerr << "NoTetra ";
                        }
                        else{
                            std::cerr << sit->GetFace(i)->GetNeighbor(i)->GetGID()<< ' ';
                        }
                    }
                    std::cerr << std::endl;

                }
            }
        }
        CheckDiST( std::cerr);

    }
    catch (DROPS::DROPSErrCL err) {err.handle();}
}
