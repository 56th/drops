/// \file stokesCoeff.h
/// \brief  coefficient class for stokes problems
/// \author LNM RWTH Aachen: Eva Loch, Yuanjun Zhang, Thorolf Schulte

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

#ifndef STOKESCOEFF_H_
#define STOKESCOEFF_H_

#include "misc/funcmap.h"
#include "misc/container.h"
#include <sstream>
#include "misc/params.h"

namespace DROPS{

class StokesFlowCoeffCL
{
  public:
  //reaction
    static InstatScalarFunction q;
  //source term
    static InstatVectorFunction f;
  //solution of velocity
    static InstatVectorFunction LsgVel;
  //solution of Jacobi-matrix of exact solution for velocity
    static InstatMatrixFunction DLsgVel;
    //solution of pressure
    static InstatScalarFunction LsgPr;

    const double rho, nu;
    const DROPS::Point3DCL g;

    StokesFlowCoeffCL( const DROPS::ParamCL& P)
      : rho( P.get<double>("NavStokes.Coeff.Dens")),
        nu( P.get<double>("NavStokes.Coeff.Visc")),
        g( P.get<DROPS::Point3DCL>("NavStokes.Coeff.Gravity")){
        DROPS::InScaMap & scamap = DROPS::InScaMap::getInstance();
        q = scamap[P.get<std::string>("NavStokes.Coeff.Reaction")];
        DROPS::InVecMap & vecmap = DROPS::InVecMap::getInstance();
        f = vecmap[P.get<std::string>("NavStokes.Coeff.Source")];
        if( P.get<std::string>("NavStokes.Coeff.Solution_Vel").compare("None")!=0)
            LsgVel = vecmap[P.get<std::string>("NavStokes.Coeff.Solution_Vel")];
        else
            LsgVel = NULL;
        DROPS::InMatMap & matmap = DROPS::InMatMap::getInstance();
        if( P.get<std::string>("NavStokes.Coeff.Solution_DVel").compare("None")!=0)
            DLsgVel = matmap[P.get<std::string>("NavStokes.Coeff.Solution_DVel")];
        else
            DLsgVel = NULL;
        DROPS::InScaMap & inscamap = DROPS::InScaMap::getInstance();
        if( P.get<std::string>("NavStokes.Coeff.Solution_Pr").compare("None")!=0)
            LsgPr = inscamap[P.get<std::string>("NavStokes.Coeff.Solution_Pr")];
        else
            LsgPr = NULL;
    }



};

InstatScalarFunction StokesFlowCoeffCL::q;
InstatVectorFunction StokesFlowCoeffCL::f;
InstatVectorFunction StokesFlowCoeffCL::LsgVel;
InstatMatrixFunction StokesFlowCoeffCL::DLsgVel;
InstatScalarFunction StokesFlowCoeffCL::LsgPr;
}//end of namespace

#endif
