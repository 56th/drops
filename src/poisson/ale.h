/// \file  ale.h
/// \brief classes that move the grids for ale method
/// \author LNM RWTH Aachen: Liang Zhang; SC RWTH Aachen;

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

#ifndef DROPS_ALE_H
#define DROPS_ALE_H

#include "misc/container.h"
#include "misc/params.h"
#include "misc/funcmap.h"
#include "poisson/poissonParam.h"
#include <sstream>

namespace DROPS{

class ALECL{
///A class to handle ALE method for free surface one phase scalar problem;
    private:
    const bool IfALE_;
    const ParamCL Para_;
    MultiGridCL& mg_;
    const double dt_;        //step size of time integration
    const double Ly_;        //Height of the reference interface

  public:
    //free surface function
    instat_scalar_fun_ptr interface_;

    ALECL(ParamCL P, MultiGridCL& mg):
    IfALE_(P.get<int>("ALE.wavy")),
    Para_(P), mg_(mg), dt_(P.get<int>("Time.NumSteps")!=0 ? P.get<double>("Time.FinalTime")/P.get<int>("Time.NumSteps") : 0),
    Ly_(norm(P.get<DROPS::Point3DCL>("Mesh.E2")))
    {
        DROPS::InScaMap & scamap = DROPS::InScaMap::getInstance();
        interface_ = scamap[P.get<std::string>("ALE.Interface")];
    }
    bool GetALE() const {return IfALE_;}
    //Initialize the grids
    void InitGrid() const;
    //Scale the grids according to the free surface functions
    void MovGrid(double t) const;
};

}
#endif
