/// \file scalarFunctions.cpp
/// \brief collections of general scalar functions (like zero, one, etc..). No problem-specific functions!
/// \author LNM RWTH Aachen: Martin Horsky; SC RWTH Aachen:
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
#include "misc/container.h"
#include "misc/funcmap.h"

#ifndef SCALARFUNCTIONS_H_
#define SCALARFUNCTIONS_H_

namespace DROPS
{

//========================================================================
//                         General Functions
//========================================================================
/// returning zero
double Zero( const Point3DCL&, double) { return 0.; }
/// returning zero (as scalar_tetra_function)
double ZeroTet( const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) { return 0.; }
/// returning one
double One( const Point3DCL&, double) { return 1.; }
/// returning one (as scalar_tetra_function)
double OneTet( const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) { return 1.; }


//========================================================================
//                   Registrierung der Funktionen
//========================================================================
static RegisterScalarFunction regscazero("Zero", Zero);
static RegisterScalarFunction regscaone("One", One);
static RegisterScalarTetraFunction regscazerotet("Zero", ZeroTet);
static RegisterScalarTetraFunction regscaonetet("One", OneTet);

}//end namespace DROPS
#endif /* SCALARFUNCTIONS_H_ */
