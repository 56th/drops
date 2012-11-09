/// \file velFunctions.cpp
/// \brief collections of general vector functions (like zero, one, etc..). No problem-specific functions!
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

#ifndef VECTORFUNCTIONS_H_
#define VECTORFUNCTIONS_H_

//========================================================================
//                         General Functions
//========================================================================
/// returns vector of zero velocities
DROPS::Point3DCL ZeroVel( const DROPS::Point3DCL&, double) { return DROPS::Point3DCL(0.); }
/// returns vector of zero velocities (as vector_tetra_function)
DROPS::Point3DCL ZeroVelTet( const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) { return DROPS::Point3DCL(0.); }

/// returns unity vector in direction \a D
template <int D>
DROPS::Point3DCL UnitVel( const DROPS::Point3DCL&, double) { DROPS::Point3DCL p(0.); p[D]=1.0; return p; }
/// returns unity vector in direction \a D (as vector_tetra_function)
template <int D>
DROPS::Point3DCL UnitVelTet( const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) { DROPS::Point3DCL p(0.); p[D]=1.0; return p; }

//========================================================================
//            Registration of functions in the func-container
//========================================================================
static DROPS::RegisterVectorFunction regvelzerovel("ZeroVel", ZeroVel);
static DROPS::RegisterVectorFunction regvelunitvelx("UnitVelx", UnitVel<0>);
static DROPS::RegisterVectorFunction regvelunitvely("UnitVely", UnitVel<1>);
static DROPS::RegisterVectorFunction regvelunitvelz("UnitVelz", UnitVel<2>);
static DROPS::RegisterVectorTetraFunction regvelzeroveltet("ZeroVel", ZeroVelTet);
static DROPS::RegisterVectorTetraFunction regvelunitvelxtet("UnitVelx", UnitVelTet<0>);
static DROPS::RegisterVectorTetraFunction regvelunitvelytet("UnitVely", UnitVelTet<1>);
static DROPS::RegisterVectorTetraFunction regvelunitvelztet("UnitVelz", UnitVelTet<2>);

#endif /* VECTORFUNCTIONS_H_ */
