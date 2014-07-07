/// \file spacetime_quad.cpp
/// \brief helper functions to construct integraiton rules on (cutted) 4d-geometries (consisting of pentatopes)
/// \author LNM RWTH Aachen: Christoph Lehrenfeld

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
 * Copyright 2012 LNM/SC RWTH Aachen, Germany
*/

#include "num/spacetime_geom.h"
#include "num/spacetime_quad.tpp"
#include "num/discretize.h"

namespace DROPS
{

// explicit template instantiations for Quad3_4DDataCL, Quad3DataCL and the corresponding tuple:

template 
void Gather4DIntegrationPoints<QuadRule::Volume>(const std::vector<PentatopeCL> &, GridFunctionCL<Point4DCL> & );
template 
void Gather4DIntegrationWeights<QuadRule::Volume>(const std::vector<PentatopeCL> &, GridFunctionCL<double> &);
template 
void Gather4DIntegrationPoints<QuadRule::Surface>(const std::vector<Tetra4DCL> &, GridFunctionCL<Point4DCL> & );
template 
void Gather4DIntegrationWeights<QuadRule::Surface>(const std::vector<Tetra4DCL> &, GridFunctionCL<double> &);
template 
void Gather4DNormals<QuadRule::Surface>(const std::vector<Tetra4DCL> &, GridFunctionCL<Point4DCL> &);
template 
void Gather4DNu<QuadRule::Surface>(const std::vector<Tetra4DCL> &, GridFunctionCL<double> &);
// -
template 
SArrayCL<Point4DCL,QuadRule::Volume::NumNodesC> Transform4DIntegrationPoints<QuadRule::Volume>(const PentatopeCL & penta);
template 
SArrayCL<Point4DCL,QuadRule::Surface::NumNodesC> Transform4DIntegrationPoints<QuadRule::Surface>(const Tetra4DCL & penta);
template 
SArrayCL<double,QuadRule::Volume::NumNodesC> Transform4DIntegrationWeights<QuadRule::Volume>(const PentatopeCL & penta);
template 
SArrayCL<double,QuadRule::Surface::NumNodesC> Transform4DIntegrationWeights<QuadRule::Surface>(const Tetra4DCL & tet);
// -


template class CompositeSTQuadCL<QuadRule >;
/* template instationations for T = double and T = Point3DCL */ 
template double CompositeSTQuadCL<QuadRule >::QuadOnPart ( const GridFunctionCL<double> & f, 
                                                                         bool posPart) const;
template Point3DCL CompositeSTQuadCL<QuadRule>::QuadOnPart ( const GridFunctionCL<Point3DCL> & f, 
                                                                           bool posPart) const;
template double CompositeSTQuadCL<QuadRule >::QuadOnInterface ( const GridFunctionCL<double> & f) const;
template Point3DCL CompositeSTQuadCL<QuadRule>::QuadOnInterface ( const GridFunctionCL<Point3DCL> & f) const;
template Point4DCL CompositeSTQuadCL<QuadRule>::QuadOnInterface ( const GridFunctionCL<Point4DCL> & f) const;

template GridFunctionCL<double> CompositeSTQuadCL<QuadRule>::EvalLinearOnPart ( const LocalP2CL<double>& fold, const LocalP2CL<double>& fnew, bool posPart) const;

template GridFunctionCL<double> CompositeSTQuadCL<QuadRule>::EvalLinearOnInterface ( const LocalP2CL<double>& fold, const LocalP2CL<double>& fnew) const;

template GridFunctionCL<double> CompositeSTQuadCL<QuadRule>::EvalLinear ( const LocalP2CL<double> & fold, const LocalP2CL<double> & fnew, const GridFunctionCL<Point4DCL> & points) const;

template GridFunctionCL<Point3DCL> CompositeSTQuadCL<QuadRule>::EvalLinearOnPart ( const LocalP2CL<Point3DCL>& fold, const LocalP2CL<Point3DCL>& fnew, bool posPart) const;

template GridFunctionCL<Point3DCL> CompositeSTQuadCL<QuadRule>::EvalLinearOnInterface ( const LocalP2CL<Point3DCL>& fold, const LocalP2CL<Point3DCL>& fnew) const;

template GridFunctionCL<Point3DCL> CompositeSTQuadCL<QuadRule>::EvalLinear ( const LocalP2CL<Point3DCL> & fold, const LocalP2CL<Point3DCL> & fnew, const GridFunctionCL<Point4DCL> & points) const;


} // end of namespace DROPS

