/// \file bndmap.h
/// \brief global map for boundary functions
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

#ifndef BNDMAP_H
#define BNDMAP_H

#include <map>
#include <string>
#include "misc/singletonmap.h"
#include "num/discretize.h"


namespace DROPS
{


typedef SingletonMapCL<DROPS::instat_vector_fun_ptr> InVecMap;
typedef SingletonMapCL<DROPS::instat_scalar_fun_ptr> InScaMap;
typedef SingletonMapCL<DROPS::vector_tetra_function> VecTetMap;
typedef SingletonMapCL<DROPS::scalar_tetra_function> ScaTetMap;
typedef SingletonMapCL<DROPS::match_fun> MatchMap;
typedef SingletonMapCL<DROPS::instat_matrix_fun_ptr> InMatMap;


Point3DCL TestFunction(const Point3DCL& , double);

typedef MapRegisterCL< instat_vector_fun_ptr> RegisterVectorFunction;
typedef MapRegisterCL< vector_tetra_function> RegisterVectorTetraFunction;

typedef MapRegisterCL< instat_scalar_fun_ptr> RegisterScalarFunction;
typedef MapRegisterCL< scalar_tetra_function> RegisterScalarTetraFunction;


typedef MapRegisterCL< match_fun> RegisterMatchingFunction;

typedef MapRegisterCL< instat_matrix_fun_ptr> RegisterMatrixFunction;


} //end of namespace DROPS

#endif
