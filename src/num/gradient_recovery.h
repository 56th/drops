/// \file gradient_recovery.h
/// \brief Compute a Lipschitz approximation of the gradient of a finite element function.
/// \author Joerg Grande

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
 * Copyright 2014, 2015, 2016 Joerg Grande, Aachen, Germany
*/


#include "misc/problem.h"

#ifndef DROPS_GRADIENT_RECOVERY_H
#define DROPS_GRADIENT_RECOVERY_H

namespace DROPS
{

class MultiGridCL;
template <class>
  class BndDataCL;

// void ppr_gradient_recovery (const MultiGridCL& mg, const VecDescCL& f, const BndDataCL<>& fbnd, VecDescCL& grad);

/// \brief Compute a continuous approximation of the gradient of the P2-FE f by averaging in the tetras neighboring a dof.
/// grad must be a Point3D-valued P2-FE without Dirichlet-boundary conditions.
void averaging_P2_gradient_recovery (const MultiGridCL& mg, const VecDescCL& f, const BndDataCL<double>& fbnd, VecDescCL& grad);

} // end of namespace DROPS

#endif