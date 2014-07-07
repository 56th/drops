/// \file dist_interface.cpp
/// \brief test of interfaces for DiST
/// \author LNM RWTH Aachen: Patrick Esser; SC RWTH Aachen: Oliver Fortmeier

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
    Uint ref[3]= { 2, 2, 2};
    BrickBuilderCL builder( origin, e1, e2, e3, ref[0], ref[1], ref[2]);
    mg = new MultiGridCL( builder);
}

/// \brief Process 0 sends its i-th tetrahedron to process i
void SendTetras( MultiGridCL& mg, const bool binary=true)
{
    DiST::TransferCL transfer( mg, true, binary);
    transfer.Init();
    if ( ProcCL::IamMaster()){
        int            toProc= 0;
        const Priority prio  = PrioMaster;
        const bool     todel = true;
        for ( MultiGridCL::const_TetraIterator sit( mg.GetTetrasBegin()); sit!=mg.GetTetrasEnd(); ++sit){
            transfer.Transfer( *sit, toProc, prio, todel);
            toProc= (toProc+1)%ProcCL::Size();
        }
    }
    transfer.Finalize();
}

/// \brief Verts get PrioMaster/PrioGhost
void MixPrios( MultiGridCL& mg, const bool binary=true)
{
    DiST::ModifyCL mod( mg, true, binary);
    mod.Init();
    if ( ProcCL::IamMaster()){
        bool gh= true;
        for (MultiGridCL::VertexIterator it= mg.GetAllVertexBegin(), end= mg.GetAllVertexEnd(); it!=end; ++it, gh=!gh)
            mod.ChangePrio( *it, gh ? PrioGhost : PrioMaster);
    }
    mod.Finalize();
}

DiST::PrioListT prioFrom, prioTo;
DiST::InterfaceCL::DimListT dims;
DiST::LevelListCL lvls;

bool Executer( DiST::TransferableCL& t)
{
    const bool rightPrio= prioFrom.contains( t.GetPrio());
    const bool rightLvl=  lvls.contains( t.GetLevel());
    const bool rightDim = std::find(dims.begin(), dims.end(), t.GetDim())!=dims.end();
    if (  !rightPrio || !rightDim || !rightLvl) {
        printf( "Executer of Transferable dim %d, lvl %d, bary (%f,%f,%f) is called, but must not\n",
            t.GetDim(), t.GetLevel(), t.GetBary()[0], t.GetBary()[1], t.GetBary()[2]);
        return false;
    }
    return true;
}

struct BaryHandlerCL
{
    bool Gather( DiST::TransferableCL& t, DiST::SendStreamCL& send)
    {
        send << t.GetBary()[0] << t.GetBary()[1] << t.GetBary()[2];
        return true;
    }

    bool Scatter( DiST::TransferableCL& t, size_t numData, DiST::MPIistreamCL& recv)
    {
        Point3DCL bary;
        const DiST::RemoteDataCL& rd= DiST::InfoCL::Instance().GetRemoteData(t);
        size_t numFrom= 0;
        for (DiST::RemoteDataCL::ProcList_const_iterator it= rd.GetProcListBegin(), end= rd.GetProcListEnd(); it!=end; ++it)
            if (prioFrom.contains(it->prio))
                numFrom++;
        bool correct= numData==numFrom;

        for ( size_t i=0; i<numData; ++i){
            recv >> bary[0] >> bary[1] >> bary[2];
            if ( !(t.GetBary()==bary)){
                correct= false;
            }
        }
        return correct;
    }
};

struct ProcRankToOwnerHandlerCL
{
    bool Gather( DiST::TransferableCL&, DiST::SendStreamCL& send)
    {
        send << ProcCL::MyRank();
        return true;
    }

    bool Scatter( DiST::TransferableCL& t, size_t numData, DiST::MPIistreamCL& recv)
    {
        std::vector<int> ranks( numData, -1);
        for ( size_t i=0; i<numData; ++i){
            recv >> ranks[i];
        }
        if ( !t.AmIOwner()){
            return false;
        }
        return true;
    }
};

struct ProcRankFromOwnerHandlerCL
{
    bool Gather( DiST::TransferableCL& t, DiST::SendStreamCL& send)
    {
        if ( !t.AmIOwner()){
            printf("Proc %d, called GatherProcRankFromOwner for a non-owning simplex\n",
                ProcCL::MyRank());
        }
        send << ProcCL::MyRank();
        return true;
    }

    bool Scatter( DiST::TransferableCL& t, size_t numData, DiST::MPIistreamCL& recv)
    {
        bool correct=true;
        std::vector<int> ranks( numData, -1);
        for ( size_t i=0; i<numData; ++i){
            recv >> ranks[i];
        }

        if ( ranks.size()!=1){
            correct=false;
        }
        if( ranks[0]!=t.GetOwner()){
            correct= false;
        }
        return correct;
    }
};


void CheckInterface( const MultiGridCL&, const bool binary=true)
{
    DiST::InterfaceCL interf( lvls, prioFrom, prioTo, dims, binary);

    if ( !ProcCL::Check(interf.ExecuteLocal( Executer))){
        std::cout << "LocalExecute returns false" << std::endl;
    }
    else{
        std::cout << "LocalExecute returns true" << std::endl;
    }
    BaryHandlerCL BaryHandler;
    if ( !ProcCL::Check(interf.Communicate( BaryHandler))){
        std::cout << "In interface communication: Scatter returns false" << std::endl;
    }
    else{
        std::cout << "Interface communication seems to be OK" << std::endl;
    }
    ProcRankToOwnerHandlerCL ProcRankToOwnerHandler;
    if ( !ProcCL::Check( interf.InformOwners( ProcRankToOwnerHandler))){
        std::cout << "Inform owners seems to be broken" << std::endl;
    }
    else{
        std::cout << "Inform owners seems to be OK" << std::endl;
    }
    ProcRankFromOwnerHandlerCL ProcRankFromOwnerHandler;
    if ( !ProcCL::Check( interf.InformCopies( ProcRankFromOwnerHandler))){
        std::cout << "Inform copies seems to be broken" << std::endl;
    }
    else{
        std::cout << "Inform copies seems to be OK" << std::endl;
    }
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
        DROPS::SendTetras( *mg, binary);
        mg->SizeInfo( std::cout);
        DROPS::MixPrios( *mg, binary);

        DROPS::dims.push_back(0);
        DROPS::dims.push_back(1);
        DROPS::dims.push_back(2);
//        dims.push_back(3);

        std::cout << "\nFirst checking interface for prioFrom == prioTo ...\n";
        DROPS::prioFrom.push_back( DROPS::PrioMaster);
        DROPS::prioTo= DROPS::prioFrom; // Ma -> Ma
        CheckInterface( *mg, binary);

        std::cout << "\nNow checking interface for prioFrom != prioTo ...\n";
        DROPS::prioFrom.clear();
        DROPS::prioFrom.push_back( DROPS::PrioGhost); // Gh -> Ma
        CheckInterface( *mg, binary);

        std::string filename("sane.chk");
        DROPS::ProcCL::AppendProcNum(filename);
        std::ofstream sanityfile( filename.c_str());
        const DROPS::DiST::InfoCL& info= DROPS::DiST::InfoCL::Instance();
        info.SizeInfo(sanityfile);
        info.Instance().DebugInfo(sanityfile);
        sanityfile << "\n\n===========================\n\n";

        if ( DROPS::ProcCL::Check( DROPS::DiST::InfoCL::Instance().IsSane( sanityfile))){
            std::cout << " DiST-module seems to be alright!" << std::endl;
        }
        else{
            std::cout << " DiST-module seems to be broken!" << std::endl;
        }
        delete mg; mg=0;

    }
    catch (DROPS::DROPSErrCL err) {err.handle();}
    return 0;
}
