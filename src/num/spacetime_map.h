/// \file spacetime_map.h
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

#ifndef DROPS_SPACETIME_MAP_H
#define DROPS_SPACETIME_MAP_H

// #include <list>
#include <vector>
// #include <algorithm>
/* #include <iostream> */
// #include <iterator>
#include "misc/utils.h"
#include "geom/boundary.h"
#include "geom/topo.h"
#include "geom/simplex.h"
#include "geom/principallattice.h"
#include "num/spacetime_geom.h"
//#include "num/unknowns.h"
// 


namespace DROPS
{

class SpaceTimeMapping
{
public:
    virtual Point4DCL Map(const Point4DCL &) const = 0;
};

class SpaceTimeIdentity : public SpaceTimeMapping
{
private:
    SpaceTimeIdentity(){;}
public:
    virtual Point4DCL Map(const Point4DCL &p ) const;
    static SpaceTimeIdentity& getInstance()
    { 
        static SpaceTimeIdentity instance;
        return instance;
    }
};



class TPSpaceTimeMapping : public SpaceTimeMapping
{
private:
    /* const TetraCL & tet; */
    const TimeInterval tatb;
    Point3DCL a;
    Point3DCL diff[3];
public:
    TPSpaceTimeMapping(const TetraCL & tet, const TimeInterval & tatb_in)
        :/*tet(tet_in),*/tatb(tatb_in) 
    {
        a = tet.GetVertex(0)->GetCoord();
        for (int i = 1; i < 4; ++i)
            diff[i-1] = tet.GetVertex(i)->GetCoord() - a;
    };

    virtual Point4DCL Map(const Point4DCL & p) const;
};


} // end of namespace DROPS

#endif
