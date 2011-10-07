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
#include <vector>

#ifdef _PAR
# include "parallel/parallel.h"
# include "parallel/parmultigrid.h"
# include "parallel/loadbal.h"
# include "parallel/partitioner.h"
# include "parallel/logger.h"
#endif

namespace DROPS
{

/// \brief Adaptive triangulation based on position of the interface provided by the levelset function
class AdapTriangCL
{
  private:
    MultiGridCL& mg_;
#ifdef _PAR
    ParMultiGridCL*  pmg_;                                  ///< Reference to the parallel multigrid
    LoadBalHandlerCL lb_;                                   ///< Reference to the load balancing class
#endif

    double width_;                                          ///< width of the refined grid
    int c_level_, f_level_;                                 ///< coarsest and finest level of the grid
    bool modified_;                                         ///< flag if the grid has been modified

    typedef ObservedVectorsCL ObserverContT;                ///< type for observing the multigrid dependent FE functions
    ObserverContT& observer_;                               ///< stores handlers to manipulate FE-functions due to grid changes (refinement, migration)

    /// \name Evaluate a function on a simplex
    //@{
    template <class DistFctT>
    double GetValue( const DistFctT& dist, const VertexCL& v) { return dist.val( v); }
    template <class DistFctT>
    double GetValue( const DistFctT& dist, const EdgeCL& e)   { return dist.val( e); }
    template <class DistFctT>
    double GetValue( const DistFctT& dist, const TetraCL& t)  { return dist.val( t, 0.25, 0.25, 0.25); }
    double GetValue( scalar_fun_ptr dist, const VertexCL& v)  { return dist( v.GetCoord() ); }
    double GetValue( scalar_fun_ptr dist, const EdgeCL& e)    { return dist( GetBaryCenter( e) ); }
    double GetValue( scalar_fun_ptr dist, const TetraCL& t)   { return dist( GetBaryCenter( t) ); }
    //@}

    /// \brief On step of the grid change
    template <class DistFctT>
    bool ModifyGridStep( DistFctT&, bool lb=true);

  public:
    AdapTriangCL(  MultiGridCL& mg, double width, int c_level, int f_level, __UNUSED__ int lbStrategy = 1, __UNUSED__ int partitioner=1)
      : mg_( mg),
#ifdef _PAR
      pmg_( ParMultiGridCL::InstancePtr()), lb_( mg_, Partitioner(partitioner)),
#endif
      width_(width), c_level_(c_level), f_level_(f_level), modified_(false), observer_(ObserverContT::Instance())
      {
        Assert( 0<=c_level && c_level<=f_level, "AdapTriangCL: Levels are cheesy.\n", ~0);
#ifdef _PAR
        pmg_->AttachTo( mg_);
        if (lbStrategy>=0)
            lb_.DoInitDistribution( ProcCL::Master());
        else
            lbStrategy*=-1;
        switch ( lbStrategy) {
            case 0 : lb_.SetStrategy( NoMig);     break;
            case 1 : lb_.SetStrategy( Adaptive);  break;
            case 2 : lb_.SetStrategy( Recursive); break;
        }
#endif
      }

#ifdef _PAR
    /// \brief Get a reference onto the parallel MultiGrid
    ParMultiGridCL& GetPMG() { return *pmg_; }

    /// \brief Get a constant reference onto the parallel MultiGrid
    const ParMultiGridCL& GetPMG() const { return *pmg_; }

    /// \brief Get a reference onto the LoadBalHandlerCL
    LoadBalHandlerCL& GetLb() { return lb_; }

    /// \brief Get a constant reference onto the LoadBalHandlerCL
    const LoadBalHandlerCL& GetLb() const { return lb_; }
#endif

    void SetWidth       (double width) { width_  = width; }
    void SetCoarseLevel (int c_level)  { c_level_= c_level; }
    void SetFineLevel   (int f_level)  { f_level_= f_level; }

    /// \brief Make initial triangulation according to a distance function
    template <class DistFctT>
      void MakeInitialTriang (DistFctT&);

    /// \brief Coupling of all necessary steps update the triangulation according to a levelset function
    void UpdateTriang (const LevelsetP2CL& ls);

    /// \brief Check if the triangulation has been modified within last update
    bool WasModified () const { return modified_; }

    /// \name Get a reference onto the MultiGrid
    //@{
    const MultiGridCL& GetMG() const { return mg_; }
    MultiGridCL&       GetMG()       { return mg_; }
    //@}

    /// \brief Push back a handler for FE-functions to apply changes due to grid modifications
    void push_back (MGObserverCL* o)
    {
        observer_.push_back( o);
#ifdef _PAR
        observer_.back()->SetPMG( *pmg_);
#endif
    }
};

} // end of namespace DROPS

#include "levelset/adaptriang.tpp"

#endif
