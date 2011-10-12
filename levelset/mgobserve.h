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
#  include "parallel/parmultigrid.h"
#endif

namespace DROPS
{

/// \brief Observer-base-class for the observer-pattern
///
/// AdapTriangCL calls these methods around multigrid-changes. These can be used
/// to save information prior to grid changes or to repair functions on a multigrid
/// after such changes. Specific behavior must be added through subclassing.
/// \todo (merge) Put a reference on MultiGridCL (or ParMultiGridCL as member?) to this class?
/// \todo Number of index (velocity, pressure, levelset) for the ParMultiGridCL should be more comfortable.
class MGObserverCL
{
#ifdef _PAR
  protected:
    std::vector<double> recvBuf_;   ///< during migration, buffer the received values here
    
    /// \brief Copy vector elements from the old describer and receive buffer into the new vector
    void CopyVecElements( VecDescCL& new_vec_desc, const std::vector<Usint>& dims);

      /// \todo OF: Do we need still the ParMultiGridCL?
    ParMultiGridCL* pmg_;       ///< All parallel Observers need the ParMultiGridCL
#endif

  public:
#ifdef _PAR
    MGObserverCL(ParMultiGridCL& pmg) : pmg_(&pmg) {}
    MGObserverCL() : pmg_(0) {}

    /// Get parallel multigrid.
    ParMultiGridCL& GetPMG() { return *pmg_; }
    /// Set parallel multigrid.
    void SetPMG( ParMultiGridCL& pmg) { pmg_= &pmg; }
#endif
    virtual ~MGObserverCL () {};

    /// Called immediately before MultiGridCL::Refine() in AdapTriangCL::ModifyGridStep.
    virtual void pre_refine  ()= 0;
    /// Called immediately after MultiGridCL::Refine() in AdapTriangCL::ModifyGridStep.
    virtual void post_refine ()= 0;

    /// Called at the beginning of AdapTriangCL::UpdateTriang().
    virtual void pre_refine_sequence  ()= 0;
    /// Called at the end of AdapTriangCL::UpdateTriang().
    virtual void post_refine_sequence ()= 0;
#ifdef _PAR
    /// \name Used for migration handling. Observer is skipped if GetVector() returns null_ptr.
    // @{
    /// \name gatter and setter for the receive buffer and exteded indices
    //@{
    const std::vector<double>& GetRecvBuffer() const { return recvBuf_; }
          std::vector<double>& GetRecvBuffer()       { return recvBuf_; }
    //@}
    /// Get a pointer to the FE space.
    virtual const IdxDescCL* GetIdxDesc() const= 0;
    /// Get a pointer to an (accumulated) vector. Return value null_ptr indicates that this observer should be skipped for migration handling.
    virtual const VectorCL* GetVector() const= 0;
    /// Swap FE space and vector.
    virtual void swap( IdxDescCL& idx, VectorCL& v)= 0;

    virtual void pre_migrate () = 0;
    virtual void post_migrate() = 0;
    // @}
#endif
};

///\brief Container for MGObserverCL objects, e.g., solution vectors which should be adapted after refinement or (in the parallel case) migration
class ObservedVectorsCL: public std::vector<MGObserverCL*>
{
  private:
    typedef std::vector<MGObserverCL*> base;
    // Singleton: make std ctr's and assignment operator private
    ObservedVectorsCL() : base() {}
    ObservedVectorsCL( const ObservedVectorsCL&);            // copy ctr not defined
    ObservedVectorsCL& operator=( const ObservedVectorsCL);  // assignment operator not defined

  public:
    /// Get instance of Singleton
    static ObservedVectorsCL& Instance()
    {
        static ObservedVectorsCL instance;
        return instance;
    }

    /// \name Call handlers (MGObserverCL) to manipulate FE-functions
    // @{
    /// \brief Tell Observer, that MG will be refined (and migration will be performed)
    void notify_pre_refine () {
        for (iterator obs= begin(); obs != end(); ++obs)
            (*obs)->pre_refine();
    }

    /// \brief Tell Observer, that MG has been refined (and migration has been performed)
    void notify_post_refine () {
        for (iterator obs= begin(); obs != end(); ++obs)
            (*obs)->post_refine();
    }

    /// \brief Tell Observer, that MG will be refined (and migration will be performed)
    void notify_pre_migrate () {
        for (iterator obs= begin(); obs != end(); ++obs){
            (*obs)->pre_migrate();
        }
    }

    /// \brief Tell Observer, that MG has been refined (and migration has been performed)
    void notify_post_migrate () {
        for (iterator obs= begin(); obs != end(); ++obs){
            (*obs)->post_migrate();
        }
    }

    /// \brief Tell Observer, that a sequence of refinements (and migrations) will take place
    void notify_pre_refmig_sequence() {
        for (iterator obs= begin(); obs != end(); ++obs)
            (*obs)->pre_refine_sequence();
    }

    /// \brief Tell Observer, that a sequence of refinements (and migrations) has taken place
    void notify_post_refmig_sequence() {
        for (iterator obs= begin(); obs != end(); ++obs)
            (*obs)->post_refine_sequence();
#ifdef _PAR
        ParMultiGridCL::Instance().DelAllUnkRecv();
#endif
    }
    //@}
};

} // end of namespace DROPS

#endif
