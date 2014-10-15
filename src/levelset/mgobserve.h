/// \file mgobserve.h
/// \brief Observer-base-class for MultiGridCL-changes through AdaptriangCL
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_MGOBSERVER_H
#define DROPS_MGOBSERVER_H

#include "misc/utils.h"
#ifdef _PAR
#  include "parallel/migrateunknowns.h"
#endif

namespace DROPS
{

class IdxDescCL;
class MultiGridCL;

/// \brief Observer-base-class for the observer-pattern
///
/// AdapTriangCL calls these methods around multigrid-changes. These can be used
/// to save information prior to grid changes or to repair functions on a multigrid
/// after such changes. Specific behavior must be added through subclassing.
/// \todo (merge) Put a reference on MultiGridCL (or ParMultiGridCL as member?) to this class?
class MGObserverCL
{
  public:
    virtual ~MGObserverCL () {};

    /// Called immediately before MultiGridCL::Refine() in AdapTriangCL::ModifyGridStep.
    virtual void pre_refine  ()= 0;
    /// Called immediately after MultiGridCL::Refine() in AdapTriangCL::ModifyGridStep.
    virtual void post_refine ()= 0;

    /// Called at the beginning of AdapTriangCL::UpdateTriang().
    virtual void pre_refine_sequence  ()= 0;
    /// Called at the end of AdapTriangCL::UpdateTriang().
    virtual void post_refine_sequence ()= 0;
    /// Get a pointer to the FE space.
    virtual const IdxDescCL* GetIdxDesc() const= 0;

#ifdef _PAR
    /// \name Used for migration handling. Observer is skipped if GetVector() returns null_ptr.
    // @{
    /// Get a pointer to an (accumulated) vector. Return value null_ptr indicates that this observer should be skipped for migration handling.
    virtual const VectorCL* GetVector() const= 0;
    /// Swap FE space and vector.
    virtual void swap( IdxDescCL& idx, VectorCL& v)= 0;
    // @}
#endif
};

///\brief Container for MGObserverCL objects, e.g., solution vectors which should be adapted after refinement or (in the parallel case) migration
class ObservedVectorsCL: public std::vector<MGObserverCL*>
{
  private:
    typedef std::vector<MGObserverCL*> base;

  public:
    /// \name Call handlers (MGObserverCL) to manipulate FE-functions
    // @{
    /// \brief Tell Observer, that MG will be refined (and migration will be performed)
    void notify_pre_refine () {
        for (iterator obs= begin(); obs != end(); ++obs)
            (*obs)->pre_refine();
    }

#ifdef _PAR
    /// \brief Tell Observer that a migration will be performed
    void notify_pre_migrate() {
        ObservedMigrateFECL::Instance().notify_pre_migrate();
    }

    /// \brief Tell Observer that a migration has been performed
    void notify_post_migrate() {
        ObservedMigrateFECL::Instance().notify_post_migrate();
    }
#endif
    /// \brief Tell Observer, that MG has been refined (and migration has been performed)
    void notify_post_refine () {
        for (iterator obs= begin(); obs != end(); ++obs)
            (*obs)->post_refine();
    }

    /// \brief Tell Observer, that a sequence of refinements (and migrations) will take place
    void notify_pre_refmig_sequence( __UNUSED__ MultiGridCL& mg) {
        for (iterator obs= begin(); obs != end(); ++obs){
            (*obs)->pre_refine_sequence();
        }
#ifdef _PAR
        ObservedMigrateFECL::Instance().Init( *this, mg);
#endif
    }

    /// \brief Tell Observer, that a sequence of refinements (and migrations) has taken place
    void notify_post_refmig_sequence() {
        for (iterator obs= begin(); obs != end(); ++obs)
            (*obs)->post_refine_sequence();
#ifdef _PAR
        ObservedMigrateFECL::Instance().clear();
#endif
    }
    //@}
};

} // end of namespace DROPS

#endif
