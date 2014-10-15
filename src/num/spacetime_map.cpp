/// \file spacetime_map.cpp
/// \brief classes that constitute mappings from R4 to R4
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

#include "num/spacetime_map.h"

namespace DROPS
{

Point4DCL SpaceTimeIdentity::Map(const Point4DCL &p ) const { return p;}

Point4DCL TPSpaceTimeMapping::Map(const Point4DCL & p) const
{
    double tres = tatb.first + p[3] * (tatb.second - tatb.first);
    Point3DCL sp(a);
    for (int i = 0; i < 3; ++i)
        sp += p[i] * diff[i];
    return Point4DCL(sp,tres);
}

} // end of namespace DROPS

