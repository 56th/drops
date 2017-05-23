/// \file dist_remotedatatest.cpp
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
#include <iostream>
#include <sstream>

using namespace std;

int main( int argc, char **argv)
{
#ifdef _PAR
    DROPS::ProcCL::Instance(&argc, &argv);
    // The following code should only be executed by the master
    if (!DROPS::ProcCL::IamMaster()) return 0;
#endif
    try {
        using DROPS::DiST::GeomIdCL;
        using DROPS::VertexCL;
        using DROPS::DiST::RemoteDataCL;
        using DROPS::DiST::RemoteDataListCL;
        using DROPS::DiST::RemoteDataListIteratorCL;
        typedef RemoteDataCL::ProcListEntryCL ProcEntryT;

        ProcEntryT p0(0,DROPS::PrioNeutral);
        ProcEntryT p1(0,DROPS::PrioKilledGhost);
        ProcEntryT p2(0,DROPS::PrioVGhost);
        ProcEntryT p3(0,DROPS::PrioGhost);
        ProcEntryT p4(0,DROPS::PrioMaster);

        typedef vector<ProcEntryT> ProcVecT;
        ProcVecT v0, v1, v2, v3, v4, v5;
        v0.push_back(p0);
        v1.push_back(p1);
        v1.push_back(p1);
        v2.push_back(p2);
        v3.push_back(p3);
        v4.push_back(p4);
        v4.push_back(p1);
        v5.push_back(p1);

        RemoteDataCL data0(v0);
        RemoteDataCL data1(v1);
        RemoteDataCL data2(v2);
        RemoteDataCL data3(v3);
        RemoteDataCL data4(v4);
        RemoteDataCL data5(v5);

        RemoteDataListCL list;
        typedef RemoteDataListCL::value_type ValuePair;

        list.insert(ValuePair(GeomIdCL(0, DROPS::MakePoint3D( 0., 0., 0.), 0), data0));
        list.insert(ValuePair(GeomIdCL(1, DROPS::MakePoint3D( 0., 0., 0.), 0), data1));
        list.insert(ValuePair(GeomIdCL(2, DROPS::MakePoint3D( 1., 1., 1.), 0), data2));
        list.insert(ValuePair(GeomIdCL(2, DROPS::MakePoint3D( 0., 0., 0.), 0), data3));
        list.insert(ValuePair(GeomIdCL(1, DROPS::MakePoint3D( 0., 0., 1.), 0), data4));
        list.insert(ValuePair(GeomIdCL(1, DROPS::MakePoint3D( 1., 0., 1.), 0), data5));

        // this should throw an error, because the hash already exists in the list
//        list.insert( ValuePair(hashlist[4], data));
//        RemoteDataCL<VertexCL> searchedData= list[hashlist[1]];
        // check if searchedData matches!

        std::cerr << "The container has the following vertices:" << std::endl;
        for ( RemoteDataListCL::iterator itr = list.begin(); itr!=list.end(); ++itr)
        {
            std::cerr << "Coordinates: " << (*itr).first.bary[0]         <<
                                    "  " << (*itr).first.bary[1]         <<
                                    "  " << (*itr).first.bary[2]         <<
                              " Level: " << (*itr).first.level           <<
                          " LocalPrio: " << (*itr).second.GetLocalPrio() << std::endl;
        }

        DROPS::DiST::PrioListT prioList;
        prioList.push_back( DROPS::PrioMaster);
        prioList.push_back( DROPS::PrioKilledGhost);
        DROPS::Uint prioLevel = 1;
        DROPS::DiST::LevelListCL lvlList;
        lvlList.push_back( prioLevel);

        std::cerr<< "The requested properties are:" << std::endl;
        std::cerr<< "+I want vertices of level: " << prioLevel << std::endl;
        std::cerr<< "+With any of the following priorities: ";
        for( DROPS::DiST::PrioListT::iterator itp = prioList.begin();itp!=prioList.end();++itp)
        {
            std::cerr<< "   " << *itp;
        }

        std::cerr << "\nThe vertices with the requested properties are:" << std::endl;
        for (RemoteDataListIteratorCL it = list.beginLvlPrioIt(lvlList,prioList,true); it!=list.endLvlPrioIt(); ++it)
        {
            std::cerr << "Coordinates: " << (*it).first.bary[0]         <<
                                    "  " << (*it).first.bary[1]         <<
                                    "  " << (*it).first.bary[2]         <<
                              " Level: " << (*it).first.level           <<
                          " LocalPrio: " << (*it).second.GetLocalPrio() << std::endl;
            if (prioLevel != (*it).first.level)
                std::cerr << "Level does not match!\n";
            if (!prioList.contains( (*it).second.GetLocalPrio()))
                std::cerr << "Prio does not match!\n";
        }
    }
    catch (DROPS::DROPSErrCL err) {err.handle();}
    return 0;
}
