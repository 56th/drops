/// \file adaptriang.h
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

#ifndef DROPS_ADAPTRIANG_H
#define DROPS_ADAPTRIANG_H

#include "levelset/mgobserve.h"
#include "levelset/levelset.h"
#include "levelset/marking_strategy.h"
#include <vector>

#ifdef _PAR
# include "parallel/parallel.h"
# include "parallel/parmultigrid.h"
# include "parallel/loadbal.h"
# include "parallel/logger.h"
#endif

namespace DROPS
{

/*!
 * \brief Adaptive triangulation based on a given marking strategy.
 *
 * This class implements adaptive triangulations based on arbitrary marking
 * strategy. A marking strategy is an object of type MarkingStrategyCL. It
 * marks the leafs of a multigrid triangulation for removal, i.e., coarsening,
 * or refinement.
 *
 * Given such a marking strategy, this class does everything that's needed to
 * adapt a given multigrid hierarchy. The AdapTriangCL only stores a pointer to
 * a marking strategy, it does *not* call delete on that pointer.
 */
class AdapTriangCL
{
  private:
    MultiGridCL& mg_;
    MarkingStrategyCL *marker_;

#ifdef _PAR
    ParMultiGridCL*  pmg_;                                  ///< Reference to the parallel multigrid
    LoadBalCL lb_;                                          ///< The load balancing class used to determine a partitioning
#endif

    typedef ObservedVectorsCL ObserverContT;                ///< type for observing the multigrid dependent FE functions
    ObserverContT  observer_;                               ///< stores handlers to manipulate FE-functions due to grid changes (refinement, migration)

    /// \brief One step of the grid change
    bool ModifyGridStep( bool lb = true );

  public:
    AdapTriangCL(  MultiGridCL& mg, MarkingStrategyCL *marker = 0, __UNUSED__ int lbStrategy = 1011)
        : mg_( mg ), marker_( marker )
#ifdef _PAR
        , pmg_( ParMultiGridCL::InstancePtr()), lb_( mg_)
#endif
    {

#ifdef _PAR
        if (lbStrategy>=0){
            lb_.DoMigration();
            lb_.SetMethod( lbStrategy);
        }
        else // negative lbStrategy avoids initial load balancing
            lb_.SetMethod( -lbStrategy);
#endif
    }

#ifdef _PAR
    /// \brief Get a reference onto the parallel MultiGrid
    ParMultiGridCL& GetPMG() { return *pmg_; }

    /// \brief Get a constant reference onto the parallel MultiGrid
    const ParMultiGridCL& GetPMG() const { return *pmg_; }

    /// \brief Get a reference onto the LoadBalHandlerCL
    LoadBalCL& GetLb() { return lb_; }

    /// \brief Get a constant reference onto the LoadBalHandlerCL
    const LoadBalCL& GetLb() const { return lb_; }
#endif

    /// \brief Make initial triangulation.
    void MakeInitialTriang();

    /// \brief Coupling of all necessary steps update the triangulation according to the marking strategy.
    /** This function updates the triangulation according to the marking strategy.
        Therefore this function uses the marking strategy to mark the tetras.
        Afterwards it refinesthem and balances the number of tetras over the
        processors. Also the numerical data is interpolated to the new
        triangulation.
        \return WasModified() */
    bool UpdateTriang();


    /// \name Get a reference onto the MultiGrid
    //@{
    const MultiGridCL& GetMG() const { return mg_; }
    MultiGridCL&       GetMG()       { return mg_; }
    //@}

    /// \brief Push back a handler for FE-functions to apply changes due to grid modifications
    void push_back (MGObserverCL* o)
    {
        observer_.push_back( o);
    }

    /// \brief Specify the marking strategy to use.
    void set_marking_strategy( MarkingStrategyCL *s )
    {
        marker_ = s;
    }

    /// \brief Get the marking strategy currently in use.
    MarkingStrategyCL* get_marking_strategy() const
    {
        return marker_;
    }

    bool WasModified()
    {
        if ( ! marker_ ) return false;
#ifndef _PAR
        return marker_->modified();
#else
        return ProcCL::GlobalOr(marker_->modified() || lb_.GetNumMovedMultiNodes() > 0);
#endif
    }
};

} // end of namespace DROPS

#endif

