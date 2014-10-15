/// \file loadbal.h
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

#ifndef DROPS_LOADBAL_H
#define DROPS_LOADBAL_H
#include "parallel/parallel.h"
#include "parallel/partime.h"
#include "parallel/parmultigrid.h"
#include "geom/multigrid.h"
#include "misc/problem.h"
#include "DiST/DiST.h"
#include "parallel/decompose.h"
#include <map>
#include <set>
#include <iostream>
#include <fstream>

namespace DROPS{

/// \brief Class for performing a load balancing step that consists of determining
///        a decomposition of tetrahedra and migrate the tetrahedra.
/** This class makes intensive use of the class \ref PartitioningCL to
    determine a partitioning. Right now, the Triangulation- and Two-phase graph,
    see Diss. Fortmeier
    As Partitioners, ParMetis, Zoltan and Scotch are available. All choices to
    determine a partitioning are explained in detail at \ref DetermineDecompositionCL_make.
*/
class LoadBalCL
{
  private:
    MultiGridCL*             mg_;           ///< Pointer to the multigrid
    int                      TriLevel_;     ///< Triangulation level, on which the LoadBalance should be made, normally set to LastTriangLevel (static for HandlerGather)
    int                      method_;       ///< method used to determine a decomposition, see PartitioningCL::make for a detailed description
    const VecDescCL*         lset_;         ///< Eventually use information about interface for loadbalancing
    const BndDataCL<>*       lsetbnd_;      ///< Eventually use information about interface for loadbalancing
    int                      rho_I_;        ///< weight factor of intersected tetrahedra
    size_t                   movedNodes_;   ///< number of multi nodes which are transferred during the migration

    /// \brief Migrate the tetrahedra according to the decomposition
    void Migrate( const PartitioningCL&);

  public:
    /// \brief Constructor
    LoadBalCL(MultiGridCL& mg, int TriLevel=-1, int method=1011)
        : mg_( &mg), TriLevel_( TriLevel),
          method_(method), lset_(0), lsetbnd_(0), rho_I_(11)
    {}
    /// \brief Destructor
    ~LoadBalCL() {}

    /// \brief Perform the load balancing
    void DoMigration();

    /// \name Set the strategies
    //@{
    void SetMethod( const int method)
        { method_= method; }
    /// \brief Set level set, so the two-phase graph model can be used
    void SetLset( const VecDescCL& lset, const BndDataCL<>& lsetbnd, const int rho_I=11)
        { rho_I_= rho_I; lset_=&lset; lsetbnd_=&lsetbnd; }
    //@}

    /// \brief Get the number of moved multi nodes
    size_t GetNumMovedMultiNodes() const { return movedNodes_; }

    /// \name Get the MultiGrid
    //@{
    const MultiGridCL& GetMG() const { return *mg_; }
          MultiGridCL& GetMG()       { return *mg_; }
    //@}
};

} // end of namespace DROPS

#endif // _LOADBAL_H_
