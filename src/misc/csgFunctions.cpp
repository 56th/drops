/// \file csgFunctions.cpp
/// \brief make csg functions accessible through function container
/// \author LNM RWTH Aachen: Christoph Lehrenfeld; SC RWTH Aachen:
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
 * Copyright 2013 LNM/SC RWTH Aachen, Germany
*/
#include "misc/container.h"
#include "misc/funcmap.h"
#include "misc/params.h"
#include "geom/csg.h"

#ifndef CSGFUNCTIONS_H_
#define CSGFUNCTIONS_H_

extern DROPS::ParamCL P;

namespace DROPS
{

// -------------------------CSG EXAMPLE(S)------------------------------ // 

double csg_fun (const Point3DCL& x, double t)
{
    static bool first = true;
    static const CSG::BodyCL* thebody;
#pragma omp critical(csg_fun)
    if (first)
    {
        std::ifstream jsonfile( (P.get<std::string>("CSG.Library")).c_str());
        ParamCL p;
        jsonfile >> p;
        thebody= CSG::body_builder( p.get_child( P.get<std::string>("CSG.Geometry") ) );
        first = false;
    }
    
    return (*thebody)( x, t);
}

static RegisterScalarFunction regsca_pdistz("csg_fun", csg_fun);

}//end namespace DROPS
#endif /* CSGFUNCTIONS_H_ */
