/// \file adaptriang.cpp
/// \brief adaptive triangulation based on position of the interface provided by the levelset function
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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
#include "adaptriang.h"
namespace DROPS
{

void AdapTriangCL::MakeInitialTriang()
{
    if ( ! marker_ )
    {
        throw DROPSErrCL( "Attempt to initialise a mesh without a strategy.\n"
                          "Here: file __FILE__, line __LINE__.\n" );
    }

#ifndef _PAR
    TimerCL time;
#else
    ParTimerCL time;
#endif

    time.Reset();
    time.Start();
    const Uint min_ref_num = marker_->GetFineLevel();
    Uint i;
    bool modified = true;
    for (i=0; i<2*min_ref_num && modified; ++i)
        modified=ModifyGridStep(true);

    time.Stop();
    const double duration=time.GetTime();

    std::cout << "MakeInitialTriang: " << i
                << " refinements in " << duration << " seconds\n"
                << "last level: " << mg_.GetLastLevel() << '\n';
    mg_.SizeInfo( std::cout );
}

bool AdapTriangCL::ModifyGridStep( bool lb )
/** One step of grid change
    \param lb Do a load-balancing?
    \return true if modifications were necessary,
    false, if nothing changed. */
{
    if ( ! marker_ )
    {
        throw DROPSErrCL( "Attempt to mark a mesh without a strategy.\n"
                          "Here: file __FILE__, line __LINE__.\n" );
    }

    marker_->SetUnmodified();

    TetraAccumulatorTupleCL accutuple;
    accutuple.push_back( marker_ );
    accutuple( mg_.GetTriangTetraBegin(), mg_.GetTriangTetraEnd() );

    bool modified = marker_->modified();
#ifdef _PAR
    modified = ProcCL::GlobalOr(modified);
#endif

    if (modified || lb) {
        observer_.notify_pre_refine();
        mg_.Refine();
        observer_.notify_post_refine();
#ifdef _PAR

        if (lb) {
            // Do the migration process (including the unknowns)
            observer_.notify_pre_migrate();
            lb_.DoMigration();
            Assert( CheckParMultiGrid(), DROPSErrCL("AdapTriangCL::ModifyGridStep: Failure in DiST ConsCheck"),
                   DebugParallelC|DebugParallelNumC|DebugLoadBalC);
            observer_.notify_post_migrate();
        }
#endif
    }
    return modified;
}

bool AdapTriangCL::UpdateTriang()
{
    if ( ! marker_ )
    {
        throw DROPSErrCL( "Attempt to update a mesh without a strategy.\n"
                          "Here: file __FILE__, line __LINE__.\n" );
    }

#ifndef _PAR
    TimerCL time;
#else
    ParTimerCL time;
#endif

    double duration;

    marker_->SetUnmodified();
    const int min_ref_num = marker_->GetFineLevel() - marker_->GetCoarseLevel();
    int i;
    observer_.notify_pre_refmig_sequence( GetMG() );
    for ( i = 0; i < 2*min_ref_num; ++i )
    {
        if (!ModifyGridStep(true))
        {
            break;
        }
    }
    observer_.notify_post_refmig_sequence();

    time.Stop();
    duration= time.GetTime();
    std::cout << "UpdateTriang: " << i
              << " refinements/interpolations in " << duration << " seconds\n"
              << "last level: " << mg_.GetLastLevel() << '\n';
    mg_.SizeInfo( std::cout );

    return WasModified();
}

} // end of namespace DROPS

