/// \file loadbal.cpp
/// \brief Loadbalancing of tetrahedal multigrids
/// \author LNM RWTH Aachen: Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier, Timo Henrich

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

#include "parallel/loadbal.h"
#include "parallel/parallel.h"
#include "parallel/migrateunknowns.h"
#include "num/interfacePatch.h"
#include "misc/scopetimer.h"
#include <iomanip>

namespace DROPS{

/****************************************************************************
* L O A D  B A L  C L A S S                                                 *
****************************************************************************/

/** Iterate over the tetrahedra that are in the loadbalancing set and tell
    DiST to which process these tetras belong. For a detailed description see
    Diploma thesis of Sven Gross.
*/
void LoadBalCL::Migrate( const PartitioningCL& detdecomp)
{
#if DROPSDebugC
/*  Maybe, you want to use this for tracking a specific tetrahedron (also it blows up the code)
	Point3DCL p= MakePoint3D( 0.6875, 0.625, 0.0625),
	    q= MakePoint3D( 0.5, 0.75, 0.25),
	    r= MakePoint3D( 0.8125, 0.9375, 0.375);
	DiST::GeomIdCL observe1(2,p,3), observe2(0,q,3), observe3(2,r,3);
*/
#endif
    Comment("- Start Migrating"<<std::endl, DebugLoadBalC);

    ParMultiGridCL& pmg= ParMultiGridCL::Instance();
    pmg.TransferBegin();
    int me = ProcCL::MyRank();
    movedNodes_= 0;
    int dest;

    LbIteratorCL it= LbIteratorCL::makeBegin( *mg_, TriLevel_);
    for ( ; it!=LbIteratorCL::makeEnd( *mg_, TriLevel_); ++it){

        // determine the destination
        dest= detdecomp.getDestination( *it);
        if (dest==me) continue;
        ++movedNodes_;

#if DROPSDebugC
/*
        if ( it->GetGID()==observe1 || it->GetGID()==observe2 || it->GetGID()==observe3)
            cdebug << " ===> Transfer des Tetras mit GID "<<it->GetGID() << " nach " << dest << " als ";
*/
#endif

        if (it->IsUnrefined() )
        { // E1-Xfer: parent transfer without children
            pmg.Transfer( *it, dest, PrioMaster, true);
#if DROPSDebugC
/*
            if ( it->GetGID()==observe1 || it->GetGID()==observe2 || it->GetGID()==observe3)
                std::cerr << "E1-Xfer mit delete =1 und PrioMaster" << std::endl;
*/
#endif
        }
        else
        { // E2-Xfer: parent transfer with children (M1/M2-Xfer)
#if DROPSDebugC
/*
            if ( it->GetGID()==observe1 || it->GetGID()==observe2 || it->GetGID()==observe3)
                std::cerr << "E2-Xfer mit delete ="<< (it->GetPrio()==PrioGhost)
                        << " und PrioGhost" << std::endl;
*/
#endif
            pmg.Transfer( *it, dest, PrioGhost, it->GetPrio()==PrioGhost);
        }

        if (!it->IsUnrefined() )
        {
            for (TetraCL::ChildPIterator ch(it->GetChildBegin()), chend(it->GetChildEnd()); ch!=chend; ++ch)
            {
                if ((*ch)->IsUnrefined() || (*ch)->HasGhost() )
                { // M1-Xfer: master transfer of unrefined child
                    pmg.Transfer( **ch, dest, PrioMaster, true);
#if DROPSDebugC
/*
                    if ( it->GetGID()==observe1 || (*ch)->GetGID()==observe1 || (*ch)->GetGID()==observe2 || (*ch)->GetGID()==observe3)
                        cdebug <<" ===> Transfer des Tetras mit GID "<< (*ch)->GetGID()
                                << " als Kind von " << it->GetGID() << " nach " << dest
                                << " als M1-Xfer mit delete =1 und PrioMaster" << std::endl;
*/
#endif
                }
                else
                { // M2-Xfer: master transfer of refined child. Maybe additional E2-Xfer for this child (as a parent of its children).
                    const bool E2Xfer= it.IsInLbSet( **ch) && detdecomp.getDestination( **ch)!=me;

                    pmg.Transfer( **ch, dest, PrioMaster, E2Xfer);
                    if (!E2Xfer)
                    {
                        pmg.PrioChange(*ch,PrioGhost);
                    }
#if DROPSDebugC
/*
                    if ( (*ch)->GetGID()==observe1 || (*ch)->GetGID()==observe2 || (*ch)->GetGID()==observe3)
                        cdebug << " ===> Transfer des Tetras mit GID "<< (*ch)->GetGID()
                                << " als Kind von " << it->GetGID() << " nach " << dest
                                << " als M2-Xfer mit delete =" << E2Xfer
                                << " und PrioMaster und ChangePrio to Prio"<< (E2Xfer?"Master":"Ghost") << std::endl;
*/
#endif
                }
            }
        }
    }
    pmg.TransferEnd();
    movedNodes_= ProcCL::GlobalSum( movedNodes_);
}

/** This function encapsulates all necessary steps to perform a load balancing
    step. So it creates a graph representing the hierarchy of tetrahedral grids,
    determine a partitioning of the graph and migrate the tetrahedra, accordingly.
*/
void LoadBalCL::DoMigration()
{
    // if no migration should be performed or only 1 process is used, don't do it ;-)
    if ( method_/1000==0 || ProcCL::Size()<2){
        std::cout << "Skip migration" << std::endl;
        return;
    }

    Comment( "Perform load balancing step:\n - Determine a decomposition\n", DebugLoadBalC);
    PartitioningCL detdecomp( *mg_, TriLevel_>0 ? TriLevel_ : mg_->GetLastLevel());
    detdecomp.make( method_, rho_I_, lset_, lsetbnd_, &ObservedMigrateFECL::Instance());

    Comment( " - Migrate tetrahedra\n", DebugLoadBalC);
    Migrate( detdecomp);

    Comment( " - Cleaning up\n", DebugLoadBalC);
    detdecomp.clear();
}


}   // end of namespace DROPS
