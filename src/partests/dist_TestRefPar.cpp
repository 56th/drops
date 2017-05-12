/// \file dist_TestRefPar.cpp
/// \brief DiST variant of Olli's DDD refinement/coarsening test case partests/TestRefPar
/// \author LNM RWTH Aachen: Sven Gross; SC RWTH Aachen: Oliver Fortmeier

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

// "parallel" Header-Files
#include "parallel/parmultigrid.h"
#include "DiST/DiST.h"
#include "parallel/partime.h"
#include "parallel/loadbal.h"

// geometry Header-Files
#include "geom/builder.h"
#include "geom/multigrid.h"

// Ausgabe in geomview-format
#include "out/output.h"

// Parameterdatei
#include "misc/params.h"

// Standard-Header-Files fuer Ausgaben
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

using namespace std;

/****************************************************************************
* G L O B A L E  V A R I A B L E N                                          *
****************************************************************************/
// Zeiten, die gemessen werden sollen
enum TimePart{
    T_Ref,
    T_SetupGraph,
    T_CalcDist,
    T_Migration,
    T_Check,
    T_print
};

// Tabelle, in der die Zeiten stehen
DROPS::TimeStoreCL Times(6);

// Parameter-Klasse
DROPS::ParamCL C;

/****************************************************************************
    * S E T   D E S C R I B E R   F O R   T I M E S T O R E  C L                *
****************************************************************************/
void SetDescriber()
{
    Times.SetDescriber(T_Ref, "Refinement");
    Times.SetDescriber(T_SetupGraph, "Setup loadbalancing graph");
    Times.SetDescriber(T_CalcDist, "Calculate Distribution");
    Times.SetDescriber(T_Migration, "Migration");
    Times.SetDescriber(T_Check, "Checking MG");
    Times.SetDescriber(T_print, "Printing");
    Times.SetCounterDescriber("Moved MultiNodes");
}

/****************************************************************************
* F I L E  H A N D L I N G                                                  *
****************************************************************************/
void PrintGEO(const DROPS::ParMultiGridCL& pmg)
{
    const int me=DROPS::ProcCL::MyRank();
    static int num=0;
    char filename[30];
    sprintf (filename, "output/geo_%i_GEOM_%i.geo",me,num++);
    ofstream file(filename);
    file << DROPS::GeomMGOutCL(pmg.GetMG(),-1,false,0.08,0.15) << std::endl;
    file.close();
}

/****************************************************************************
* C H E C K  P A R  M U L T I  G R I D                                      *
*****************************************************************************
*   Checkt, ob die parallele Verteilung und die MultiGrid-Struktur gesund   *
*   zu sein scheint. Zudem wird der Check von DDD aufgerufen                *
****************************************************************************/
void CheckParMultiGrid(DROPS::ParMultiGridCL& pmg, int type)
{
    if (type==DROPS::MIG && !C.get<int>("Misc.CheckAfterMig")) return;
    if (type==DROPS::REF && !C.get<int>("Misc.CheckAfterRef")) return;

    DROPS::ParTimerCL time;
    double duration;

    std::cout << "  - Check of parallel MultiGrid... ";

    char dat[30];
    sprintf(dat,"output/sane%i.chk", DROPS::ProcCL::MyRank());
    ofstream check(dat);
    time.Reset();

    check << "\n======== DiST check ========\n";
    if ( DROPS::ProcCL::Check( DROPS::DiST::InfoCL::Instance().IsSane( check)))
        std::cout << " DiST-module seems to be alright!" << std::endl;
    else
        std::cout << " DiST-module seems to be broken!" << std::endl;
    check << "\n======== ParMultiGrid check ========\n";
    bool pmg_sane = pmg.IsSane(check);
    check << "\n======== MultiGrid check ========\n";
    bool mg_sane  = pmg.GetMG().IsSane(check);

    std::cout << "  - Check of parallel MultiGrid... ";
    if (DROPS::ProcCL::Check(pmg_sane && mg_sane)){
         std::cout << "OK\n";
    }
    else{
        // Always exit on error
        std::cout << "not OK!!!\n";
        std::cout << "EXIT: Error found in multigrid\n";
        exit(-1);
    }
    time.Stop();
    duration = time.GetMaxTime();
    Times.AddTime(T_Check,duration);
    if (C.get<int>("Misc.PrintTime")) std::cout << "       --> "<<duration<<" sec\n";
    check.close();

    DROPS::DiST::GeomIdCL gid(1,DROPS::MakePoint3D(0.40625, 0.46875, 0.4375),3); // T1 (0.40625 0.46875 0.4375 )
    DROPS::DiST::InfoCL::Instance().ShowSimplex( gid, cdebug);
}

/****************************************************************************
* M A R K I N G   S T R A T E G I E S                                       *
*****************************************************************************
*   Setze Markierungen auf den Tetraedern. Zum Verfeinern und Vergroebern.  *
****************************************************************************/
void MarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
         ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRegRefMark();
    }
}

void UnMarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
         ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRemoveMark();
    }
}

void MarkCorner (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Corner; Corner[0]=0.; Corner[1]=0.; Corner[2]=0.;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
         ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Corner).norm()<=0.3)
            It->SetRegRefMark();
    }
}

bool MarkAround(DROPS::MultiGridCL& mg, const DROPS::Point3DCL& p, double rad, int maxLevel= -1)
{
    bool mod=false;
    if (maxLevel==-1)
        maxLevel = mg.GetLastLevel()+1;
    for (DROPS::MultiGridCL::TriangTetraIteratorCL it= mg.GetTriangTetraBegin(),
         end= mg.GetTriangTetraEnd(); it!=end; ++it)
    {
        if ((int)it->GetLevel()<maxLevel && (GetBaryCenter(*it)-p).norm()<=rad )
        {
            it->SetRegRefMark();
            mod=true;
        }
    }
    return mod;
}

bool UnMarkAround(DROPS::MultiGridCL& mg, const DROPS::Point3DCL& p, double rad)
{
    bool mod=false;
    for (DROPS::MultiGridCL::TriangTetraIteratorCL it= mg.GetTriangTetraBegin(),
         end= mg.GetTriangTetraEnd(); it!=end; ++it)
    {
        if ((GetBaryCenter(*it)-p).norm()<=rad )
        {
            it->SetRemoveMark();
            mod=true;
        }
    }
    return mod;
}

bool UnMarkForGhostKill (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
// search for a ghost tetra and unmark all children
{
    int done=1;
    if (DROPS::ProcCL::MyRank()!=0)
        DROPS::ProcCL::Recv(&done, 1, DROPS::ProcCL::MyRank()-1, 563738);
    else
        done=0;
    if (!done)
    {
        for (DROPS::MultiGridCL::const_TetraIterator  It(mg.GetTetrasBegin(maxLevel-1)),
            ItEnd(mg.GetTetrasEnd(maxLevel-1)); It!=ItEnd && !done; ++It)
        {
            if (It->IsGhost() && It->IsRegularlyRef()){
                for (DROPS::TetraCL::const_ChildPIterator ch(It->GetChildBegin()),
                    chEnd(It->GetChildEnd()); ch!=chEnd; ++ch)
                    (*ch)->SetRemoveMark();
                std::cout << "Tetra "<<It->GetGID()<<" marked for ghost-kill by proc "<<DROPS::ProcCL::MyRank()<<std::endl;
                done=1;
            }
        }
    }
    if (DROPS::ProcCL::MyRank()<DROPS::ProcCL::Size()-1)
        DROPS::ProcCL::Send(&done, 1, DROPS::ProcCL::MyRank()+1, 563738);


    return DROPS::ProcCL::GlobalOr(done);
}

DROPS::MultiGridCL* CreateInitGrid(int master= 0)
{
    using namespace DROPS;
    MultiGridCL *mg;
    const int size = ProcCL::Size();
    DROPS::ParTimerCL time;
    double duration;

    Point3DCL e1(0.0), e2(0.0), e3(0.0), orig(0.0);

    if(C.get<int>("Refine.InitCond")==0)
    {
        const Point3DCL brk_dim = C.get<Point3DCL>("Brick.dim"),
                        brk_orig= C.get<Point3DCL>("Brick.orig");
        e1[0]=brk_dim[0]; e2[1]=brk_dim[1]; e3[2]= brk_dim[2];
        BrickBuilderCL brick(brk_orig, e1, e2, e3, C.get<int>("Brick.BasicRefX"), C.get<int>("Brick.BasicRefY"), C.get<int>("Brick.BasicRefZ"));
        mg = new DROPS::MultiGridCL(brick);
    }
    else if (C.get<int>("Refine.InitCond")==1)
    {
        e1[0]=e2[1]=e3[2]= 1.;
        if (ProcCL::MyRank()==master)
        {
            BrickBuilderCL builder(orig, e1, e2, e3, 4, 4, 4);
            FileBuilderCL fileBuilder(C.get<std::string>("Misc.InitPrefix"), &builder);
            mg = new DROPS::MultiGridCL(fileBuilder);

            std::ofstream serSanity("sanity.txt");

            std::cout << "\n \n MultiGrid mit "<<mg->GetNumLevel()<<" Leveln aus Datei gelesen\n \n";
            serSanity << SanityMGOutCL(*mg) << std::endl;
        }
        else
        {
            BrickBuilderCL builder(orig, e1, e2, e3, 4, 4, 4);
            builder.set_par_numlevel( C.get<int>("Refine.Refined")+1);
            mg = new DROPS::MultiGridCL(builder);
        }
    }
    else
    {
        throw DROPSErrCL("Unknown init condition");
    }

    // now distribute the grid on master to other procs
//    pmg.AttachTo(*mg);

    LoadBalCL lb(*mg);
    if (size>1)
    {
        time.Reset();
        lb.DoMigration();
        time.Stop();
        if (C.get<int>("Misc.PrintTime")){
            duration = time.GetMaxTime();
            std::cout << "       --> "<<duration<<" sec\n";
        }
        Times.IncCounter(lb.GetNumMovedMultiNodes());
    }

    return mg;
}

void DoMigration( DROPS::LoadBalCL &LoadBal)
{
    DROPS::ParTimerCL timer;
    timer.Reset();
    LoadBal.DoMigration();
    timer.Stop();
    if (C.get<int>("Misc.PrintTime")){
        std::cout << " Migration took "<<timer.GetMaxTime()<<" sec\n";
    }
    Times.IncCounter(LoadBal.GetNumMovedMultiNodes());
}

using namespace DROPS;
/****************************************************************************
* M A I N                                                                   *
****************************************************************************/
int main(int argc, char* argv[])
{
    DROPS::ProcCL::Instance(&argc, &argv);
    try
    {
        const char line[] = "----------------------------------------------------------------------------------";
        const char dline[]= "==================================================================================";
        SetDescriber();

        DROPS::ParTimerCL alltime, time;
        double duration;

        const int me= DROPS::ProcCL::MyRank();

        // Parameter file einlesen ...
        DROPS::read_parameter_file_from_cmdline( C, argc, argv, "../../param/partests/dist_TestRefPar/Ref.json");
        std::cout << C << std::endl;

        const bool printTime= C.get<int>("Misc.PrintTime"),
                   printSize= C.get<int>("Misc.PrintSize"),
                   printPMG=  C.get<int>("Misc.PrintPMG"),
                   printGEO=  C.get<int>("Misc.PrintGEO");
        if (printTime)
            DROPS::ParTimerCL::TestBandwidth(std::cout);

        cout << dline << endl << " + Erstelle initiales Gitter (Wuerfel der Laenge 1) auf Prozessor 0 ...\n";

        DROPS::MultiGridCL &mg = *CreateInitGrid();
        DROPS::ParMultiGridCL& pmg= DROPS::ParMultiGridCL::Instance();

        if (printSize){
            cout << "  - Verteilung der Elemente:\n";
            mg.SizeInfo(cout);
        }
        if (printPMG){
            cout << " + Schreibe Debug-Informationen in ein File ... ";
            PrintMG(pmg);
            cout << " OK\n";
        }
        if (printGEO){
            cout << " + Schreibe das Multigrid im Geomview-Format in ein File ... ";
            PrintGEO(pmg);
            cout << " OK\n";
        }

        CheckParMultiGrid(pmg,REF);

		const int markall   =  C.get<int>("Refine.All"),
		          markdrop  =  C.get<int>("Refine.Drop"),
		          markcorner=  C.get<int>("Refine.Corner"),
		          markingproc= C.get<int>("Refine.MarkingProc");
        cout << dline << endl << " Verfeinere das Gitter nun " << markall << " mal global, " << markdrop
				<< " mal in der Mitte um den Tropfen\n und " << markcorner << " mal um der Ecke (0,0,0)\n"
				<< " Es wird die Strategie ";
		switch (C.get<int>("LoadBalancing.RefineStrategy")){
			case 0 : cout << "No Loadbalancing ";break;
			case 1 : cout << "AdaptiveRefine "; break;
			case 2 : cout << "PartKWay "; break;
			default: cout << "Unbekannte Strategy ...\n EXIT"; exit(0);
		}
		cout << "verwendet. Es markiert der Prozessor " << markingproc << "\n" << dline << endl;

		int movedRefNodes=0, movedCoarseNodes=0;

//        switch (C.get<int>("Refine.Strategy"))
//        {
//            case 0:  numrefs= markall+markdrop+markcorner; break;
//            case 1:  numrefs=5; break;
//            case 2:  numrefs=markall; break;
//            case 3:  numrefs=markall; break;
//            default: throw DROPSErrCL("Specify the refinement strategy!");
//        }

        DROPS::LoadBalCL LoadBal(mg);
        for (int ref=0; ref<markall+markdrop+markcorner; ++ref)
        {
            DROPS::Point3DCL e, e1;
            bool killedghost=false;

            switch (C.get<int>("Refine.Strategy"))
            {
            case 0:
                cout << " + Refine " << (ref) << " : ";
                if (ref < markall){
                    cout << "all ...\n";
                    if (markingproc==-1 || markingproc==me)
                        DROPS::MarkAll(mg);
                }
                else if (ref < markdrop+markall){
                    cout << "drop ...\n";
                    if (markingproc==-1 || markingproc==me)
                        MarkDrop(mg, mg.GetLastLevel());
                }
                else{
                    cout << "corner ...\n";
                    if (markingproc==-1 || markingproc==me)
                        MarkCorner(mg, mg.GetLastLevel());
                }
            break;
            case 1:
                e[0]=1.; e[1]=2.; e[2]=0.5;
                e1[0]=0.;  e1[1]=0.; e1[2]=0.5;
                switch (ref)
                {
                    case 0:
                        std::cout << "Mark all "<<std::endl;
                        MarkAll(mg);
                        break;
                    case 1:
                        std::cout << "Mark all"<<std::endl;
                        MarkAll(mg);
    //                     MarkAround(mg, e1, 0.5);
                        break;
                    case 2:
                        std::cout << "Mark around "<<e<<std::endl;
                        MarkAround(mg, e, 0.5);
                        break;
                    case 3:
                        std::cout << "UnMark around "<<e<<std::endl;
                        UnMarkAround(mg, e, 0.6);
                        break;
                    case 4:
                        std::cout << "UnMark for ghost tetra kill"<<std::endl;
                        killedghost=UnMarkForGhostKill(mg, mg.GetLastLevel());
                        killedghost= ProcCL::GlobalOr(killedghost);
                        if (ProcCL::IamMaster() && killedghost)
                            std::cout << "A ghost tetra will be killed"<<std::endl;
                        break;
                    default:
                        std::cout << "I do not know this case!\n";
                }
            break;
            case 3:
                if (ref%2==0)
                    MarkAll(mg);
                else
                    UnMarkAll(mg);
            break;
            }   // end of switch C.Strategy


            time.Reset(); pmg.Refine(); time.Stop();
            duration = time.GetMaxTime();
            Times.AddTime(T_Ref,duration);
            if (printTime) std::cout << "       --> "<<duration<<" sec\n";

            if (printPMG){
                cout << "  - Schreibe Debug-Informationen in ein File ... ";
                PrintMG(pmg,REF);
                cout << " OK\n";
            }

            CheckParMultiGrid(pmg,REF);

            DoMigration( LoadBal);
            movedRefNodes += LoadBal.GetNumMovedMultiNodes();

            if (printPMG){
                cout << "  - Schreibe Debug-Informationen in ein File ... ";
                PrintMG(pmg,MIG);
                cout << " OK\n";
            }
            if (printGEO){
                cout << "  - Schreibe das Multigrid im Geomview-Format in ein File ... ";
                PrintGEO(pmg);
                cout << " OK\n";
            }
            if (printSize){
                cout << "  - Verteilung der Elemente:\n";
                mg.SizeInfo(cout);
            }

            CheckParMultiGrid(pmg,MIG);

            if (ref!=markall+markdrop+markcorner-1) cout << line << endl;
        }

        if (C.get<int>("LoadBalancing.MiddleMig")){
            cout <<dline<<endl<< " + Last-Verteilung zwischen dem Verfeinern und Vergroebern ...\n";
            DoMigration( LoadBal);
            movedRefNodes += LoadBal.GetNumMovedMultiNodes();
            CheckParMultiGrid(pmg,MIG);
        }

        const int coarseall    = C.get<int>("Coarsen.All"),
                  coarsedrop   = C.get<int>("Coarsen.Drop"),
                  unmarkingproc= C.get<int>("Coarsen.MarkingProc");
		cout <<dline<<endl << " Vergroebere nun das Gitter zunaechst " << coarsedrop
				<< " mal um den Tropfen herum und dann " << coarseall << " ueberall\n Es wird die Strategie ";
		switch (C.get<int>("LoadBalancing.CoarseStrategy")){
			case 0 : cout << "No Loadbalancing ";break;
			case 1 : cout << "AdaptiveRefine "; break;
			case 2 : cout << "PartKWay "; break;
			default: cout << "Unbekannte Strategy ...\n EXIT"; exit(0);
		}
		cout << "verwendet. Es markiert der Prozessor " << unmarkingproc << "\n" << dline << endl;

        for (int ref =0; ref<coarsedrop+coarseall; ++ref)
        {
            if (ref < coarsedrop){
                cout << " + Coarse drop (" << ref << ") ... \n";
                if (unmarkingproc==-1 || unmarkingproc==me)
                    UnMarkDrop(mg, mg.GetLastLevel());
            }
            else {
                cout << " + Coarse all (" << ref << ") ... \n";
                if (unmarkingproc==-1 || unmarkingproc==me){
                    DROPS::UnMarkAll(mg);
                }
            }

            time.Reset(); pmg.Refine(); time.Stop();
            duration = time.GetMaxTime();
            Times.AddTime(T_Ref,duration);
            if (printTime) std::cout << "       --> "<<duration<<" sec\n";
            if (printPMG)
            {
                cout << "  - Schreibe Debug-Informationen in ein File ... ";
                PrintMG(pmg, REF);
                cout << " OK\n";
            }
            if (printGEO){
                cout << "  - Schreibe das Multigrid im Geomview-Format in ein File ... ";
                PrintGEO(pmg);
                cout << " OK\n";
            }

            CheckParMultiGrid(pmg,REF);

            DoMigration( LoadBal);

            movedCoarseNodes += LoadBal.GetNumMovedMultiNodes();

            if (printSize){
                cout << "  - Verteilung der Elemente:\n";
                mg.SizeInfo(cout);
            }
            if (printPMG)
            {
                cout << "  - Schreibe Debug-Informationen in ein File ... ";
                PrintMG(pmg, MIG);
                cout << " OK\n";
            }

            CheckParMultiGrid(pmg,MIG);

            if (ref!=coarsedrop+coarseall-1) cout << line << endl;
        }

		cout << dline<< endl;

        if (printTime)
            Times.Print(cout);
		cout << "Moved Multinodes for refinement: " << movedRefNodes << endl
	  		 << "Moved Multinodes for coarsening: " << movedCoarseNodes << endl
		   	 << dline << endl << "Shuting down ...\n";

    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
    return 0;
}

