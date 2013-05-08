/// \file reftest.cpp
/// \brief tests refinement algorithm
/// \author LNM RWTH Aachen: Joerg Grande, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "out/output.h"
#include <fstream>
#include <cmath>
#include <cstdio>

void MarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    long ct= 0;

    for (DROPS::MultiGridCL::TetraIterator It(mg.GetAllTetraBegin(maxLevel)),
             ItEnd(mg.GetAllTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( It->IsUnrefined() && (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
        {
            It->SetRegRefMark(); ++ct;
        }
    }
    std::cout << ct <<" Tetras wurden markiert\n";
}


void UnMarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TetraIterator It(mg.GetAllTetraBegin(maxLevel)),
             ItEnd(mg.GetAllTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( It->IsUnrefined() && (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRemoveMark();
    }
}


int main()
{
  try
  {
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.0;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 2, 2, 2);
//    DROPS::BBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 2, 2, 2);
//    DROPS::LBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 2, 2);
    DROPS::MultiGridCL mg(brick);

//        std::cout << DROPS::DumpMGCL(mg) << std::endl << "hallo" << std::endl;
    mg.MakeConsistentNumbering();

//        std::cout << DROPS::DumpMGCL(mg);
    std::cout << DROPS::SanityMGOutCL(mg) << std::endl;
    mg.SizeInfo(std::cout);
    for (DROPS::Uint i=0; i<2; ++i)
    {
        DROPS::MarkAll(mg);
        std::cout << i+1 << ". Totalverfeinerung ------------------"<<std::endl;
        mg.Refine();
        std::cout << DROPS::SanityMGOutCL(mg) << std::endl;
        std::cout<<i+1<<".Totalverf.: "; mg.SizeInfo(std::cout);
    }
//    std::cout << DROPS::SanityMGOutCL(mg);
    for (DROPS::Uint i=0; i<4; ++i)
    {
        MarkDrop(mg, std::min(mg.GetLastLevel(),10u));
//        DROPS::MarkAll(mg);
//        MarkSome(mg);
//        std::cout << DROPS::SanityMGOutCL(mg);
//        std::cout << DROPS::DumpMGCL(mg);
        std::cout << i+1 << ". Tropfenverfeinerung ------------------"<<std::endl;
        mg.Refine();
        std::cout<<i+1<<".Tropfenverf.: "; mg.SizeInfo(std::cout);
        std::cout << DROPS::SanityMGOutCL(mg) << std::endl;
//        DebugIt(&mg);
        char str[20];
        std::ofstream ofs;
        std::sprintf(str, "drop%i.off", i+1);
        ofs.open(str);
        ofs << DROPS::GeomMGOutCL(mg, -1, false, 0.5) << std::endl;
        ofs.close();
        std::sprintf(str, "bary%i.mg", i+1);
        ofs.open(str);
        for (DROPS::MultiGridCL::TriangTetraIteratorCL it= mg.GetTriangTetraBegin(), end= mg.GetTriangTetraEnd(); it!=end; ++it)
            ofs << GetBaryCenter(*it) << '\n';
        ofs.close();
        std::sprintf(str, "verts%i.mg", i+1);
        ofs.open(str);
        for (DROPS::MultiGridCL::TriangVertexIteratorCL it= mg.GetTriangVertexBegin(), end= mg.GetTriangVertexEnd(); it!=end; ++it)
            ofs << it->GetCoord() << '\n';
        ofs.close();
    }

//    std::cout << DROPS::SanityMGOutCL(mg) << std::endl;
//    mg.SizeInfo(std::cout);

//    int wait;
//    std::cout << "Press a key: " << std::flush;
//    std::cin>>wait;

    for (DROPS::Uint i=0; i<6; ++i)
    {
        std::cout << i+1 << ". Tropfenentfeinerung ------------------"<< std::endl;
        UnMarkDrop(mg, std::min(mg.GetLastLevel(),10u));
//        DROPS::UnMarkAll(mg);
        mg.Refine();
        std::cout << DROPS::SanityMGOutCL(mg)  << std::endl;
        std::cout<<i+1<<".Tropfenentf.: "; mg.SizeInfo(std::cout);
    }
//    DROPS::UnMarkAll(mg);
//    for (DROPS::MultiGridCL::TetraIterator it= mg.GetAllTetraBegin(); it!=mg.GetAllTetraEnd(); ++it)
//        if (it->IsMarkedForRemovement()) std::cout << "MFR ";
//    std::cout << DROPS::DumpMGCL(mg);

    for (DROPS::Uint i=0; i<3; ++i)
    {
        DROPS::UnMarkAll(mg);
        std::cout << i+1 << ". Totalentfeinerung ------------------"<<std::endl;
        mg.Refine();
        std::cout << DROPS::SanityMGOutCL(mg) << std::endl;
        std::cout<<i+1<<".Totalentf.: "; mg.SizeInfo(std::cout);
    }

    {std::ofstream os("ttt2.off");
    os << DROPS::GeomMGOutCL(mg, -1, false) << std::endl;
    os.close();}
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
