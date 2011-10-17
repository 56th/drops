/// \file migrateunknowns.h
/// \brief Observer-base-class migrating unknowns
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_MIGRATEUNKNOWNS_H
#define DROPS_MIGRATEUNKNOWNS_H

#include "misc/problem.h"
#include "geom/multigrid.h"
#include "misc/utils.h"

namespace DROPS{

// fwd declaration
class MGObserverCL;
class ObservedVectorsCL;

/** \brief Class for handling the migration of one FE type such as level set, velocity, etc.

    Since this class also creates a new numbering, a reference to the MultiGridCL is used. And,
    a pointer to the MGObserverCL is also stored to perform the swap operation
*/
class MigrateFECL
{
  private:
    std::vector<double> recvBuf_;   ///< during migration, buffer the received values here
    MGObserverCL*       obs_;       ///< the corresponding MGObserverCL (for getting the index, the data and the swap function)
    const IdxDescCL*    old_idx_;   ///< index from MGObserverCL (to avoid frequent virtual function calls)
    const VectorCL*     old_vec_;   ///< vector from MGObserverCL (to avoid frequent virtual function calls)
    MultiGridCL*        mg_;        ///< MultiGridCL for creating a numbering


    /// \brief Copy vector elements from the old describer and receive buffer into the new vector
    void CopyVecElements( VecDescCL& new_vec_desc, const std::vector<Usint>& dims);

public:
    MigrateFECL( MGObserverCL* obs, MultiGridCL& mg);

    /// \brief reserve memory for the receive buffer
    void pre_migrate();

    /// \brief handle the received data
    void post_migrate();

    /// \name getters
    // @{
    const IdxDescCL* GetIdxDesc() const  { return old_idx_; }
    const VectorCL*  GetVector() const   { return old_vec_; }
    std::vector<double>& GetRecvBuffer() { return recvBuf_; }
    // @}

};

/** \brief This is a list of all FE types that are handled during the migration.

    This class is initialized by the ObservedVectorsCL
*/
class ObservedMigrateFECL : public std::vector<MigrateFECL>
{
private:
    typedef std::vector<MigrateFECL> base;       ///< base class

    /// \brief Singleton: make std ctr's and assignment operator private
    explicit ObservedMigrateFECL() : base() {}
    ObservedMigrateFECL( const ObservedMigrateFECL&);            // copy ctr not defined
    ObservedMigrateFECL& operator=( const ObservedMigrateFECL);  // assignment operator not defined

public:
    /// Get instance of Singleton
    static ObservedMigrateFECL& Instance()
    {
        static ObservedMigrateFECL instance;
        return instance;
    }

    /// \brief Init the list of observed FE types by the list of the ObservedVectorsCL
    void Init( ObservedVectorsCL&, MultiGridCL& mg);

    /// \brief notify all observed FE classes that migration will take place
    void notify_pre_migrate(){
        for ( iterator it(begin()); it!=end(); ++it){
            it->pre_migrate();
        }
    }

    /// \brief notify all observed FE classes that migration has been performed
    void notify_post_migrate(){
        for ( iterator it(begin()); it!=end(); ++it){
            it->post_migrate();
        }
    }
};

}   // end of namespace DROPS

#endif
