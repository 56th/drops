/// \file poissonCoeff.cpp
/// \brief boundary and source functions for the poisson-type problems
/// \author LNM RWTH Aachen: Christoph Lehrenfeld, Joerg Grande, Thorolf Schulte, Liang Zhang

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

#include "poisson/poissonParam.h"

//======================================================================================================================
//                                  static members declarations of poissonCoeff class
//======================================================================================================================

namespace DROPS{
ParamCL PoissonCoeffCL::C_;
double PoissonCoeffCL::dx_;
double PoissonCoeffCL::dy_;
int PoissonCoeffCL::nx_;
int PoissonCoeffCL::ny_;
double PoissonCoeffCL::dt_;
scalar_tetra_function PoissonCoeffCL::q;
double PoissonCoeffCL::alpha;
scalar_tetra_function PoissonCoeffCL::f;
scalar_tetra_function PoissonCoeffCL::Solution;
instat_scalar_fun_ptr PoissonCoeffCL::InitialCondition;
vector_tetra_function PoissonCoeffCL::Vel;
instat_scalar_fun_ptr PoissonCoeffCL::interface;
}//end of namespace
