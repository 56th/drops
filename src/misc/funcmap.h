/// \file funcmap.h
/// \brief global map for boundary and coefficient functions
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

#ifndef FUNCMAP_H
#define FUNCMAP_H

#include <map>
#include <string>
#include "misc/singletonmap.h"
#include "num/discretize.h"


namespace DROPS
{


typedef SingletonMapCL< InstatVectorFunction> InVecMap;
typedef SingletonMapCL< InstatScalarFunction> InScaMap;
typedef SingletonMapCL< vector_tetra_function> VecTetMap;
typedef SingletonMapCL< scalar_tetra_function> ScaTetMap;
typedef SingletonMapCL< MatchFunction> MatchMap;
typedef SingletonMapCL< InstatMatrixFunction> InMatMap;

typedef MapRegisterCL< InstatVectorFunction> RegisterVectorFunction;
typedef MapRegisterCL< InstatScalarFunction> RegisterScalarFunction;
typedef MapRegisterCL< vector_tetra_function> RegisterVectorTetraFunction;
typedef MapRegisterCL< scalar_tetra_function> RegisterScalarTetraFunction;
typedef MapRegisterCL< MatchFunction> RegisterMatchingFunction;
typedef MapRegisterCL< InstatMatrixFunction> RegisterMatrixFunction;


} //end of namespace DROPS

#endif
