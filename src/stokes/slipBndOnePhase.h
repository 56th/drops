/// \file slipBnd.h
/// \brief classes that add Nitsche terms for problems with slip or symmetry boundary segments.
/// \author LNM RWTH Aachen: Liang Zhang; SC RWTH Aachen:

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
 * Copyright 2016 LNM/SC RWTH Aachen, Germany
*/

#ifndef DROPS_SLIPBND_H
#define DROPS_SLIPBND_H

#include <vector>
#include "misc/problem.h"
#include "num/fe.h"
#include "num/discretize.h"
#include "geom/principallattice.h"


#ifdef _PAR
#  include "parallel/exchange.h"
#endif

namespace DROPS
{
//forward declaration
class StokesBndDataCL;
class LocalSystem1DataCL;

/// \brief Update the local system 1 (nocut) with respect to slip Bnd and symmetric Bnd;
class SlipBndSystem1OnePhaseP2CL
{
  private:
    const StokesBndDataCL& BndData_;
    double mu_;
    double beta_;            //Slip coefficient, beta_=0 for symmetric Bnd;
    const double alpha_;     //Coefficient for Nitsche method
    LocalP1CL<Point3DCL> GradRef[10];
    public:
    SlipBndSystem1OnePhaseP2CL(const StokesBndDataCL& BndData, const double alpha=0)
     : BndData_(BndData), alpha_(alpha)
    {
        P2DiscCL::GetGradientsOnRef( GradRef);
    }
    void setMu(double mu){mu_=mu;}
    void setBeta(double beta){beta_=beta;}
    /// update local system 1
    void setup(const TetraCL& tet, const SMatrixCL<3,3>& T, SMatrixCL<3,3> Ak[10][10]);
    /// for Slip boundary condition (the slip wall is moving)
    void setupRhs(const TetraCL& tet, Point3DCL loc_b[10], double t);  
};

/// \brief Due to the weak imposition of bu * n = 0 with Nitsche's method, setup the integral of (bv * bn) * q on the slip boundary for uncut element
class SlipBndSystem2OnePhaseCL
{
 private:
    const StokesBndDataCL& BndData_;
  public:
    SlipBndSystem2OnePhaseCL(const StokesBndDataCL& BndData)
     : BndData_(BndData) {}
    void setupB(const TetraCL& tet, SMatrixCL<1, 3> loc_b[10][4]);
};

}

#endif
