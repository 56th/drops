/// \file funcmap.cpp
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

#include<map>
#include "misc/funcmap.h"
#include "misc/singletonmap.h"

namespace DROPS
{
    template class SingletonMapCL<DROPS::instat_scalar_fun_ptr>;
    template class SingletonMapCL<DROPS::instat_vector_fun_ptr>;
    template class SingletonMapCL<DROPS::scalar_tetra_function>;
    template class SingletonMapCL<DROPS::vector_tetra_function>;
    template class SingletonMapCL<DROPS::instat_matrix_fun_ptr>;
    template class SingletonMapCL<DROPS::match_fun>;
} //end of namespace DROPS
