/// \file locator.h
/// \brief Find a tetra containing a given point by walking the coarse mesh and descending the refinement hierarchy.
/// \author LNM RWTH Aachen: Joerg Grande

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
 * Copyright 2014 LNM RWTH Aachen, Germany
*/

#ifndef DROPS_LOCATOR_H
#define DROPS_LOCATOR_H

#include "geom/simplex.h"

namespace DROPS {

class MultiGridCL;

/// \brief Given a point x in R^3, find a tetra tet on level lvl that contains x.
class MyLocatorCL
{
  private:
    double eps;
    bool greedy;

    const MultiGridCL& mg;
    Uint               lvl;

    World2BaryCoordCL w2b;

    Point3DCL      x;
    const TetraCL* tet;
    BaryCoordCL    xb;

    double min_bary_coeff (const BaryCoordCL& b) const;
    void locate_in_tetra ();

  public:
    MyLocatorCL(const MultiGridCL& MG, Uint level, bool greedyarg= false)
        : eps( 1e-10), greedy( greedyarg), mg( MG), lvl( level), tet( 0) {}
    // default copy-ctor, dtor, assignment-op

    void set_epsilon (double neweps) { eps= neweps; }
    void locate (const Point3DCL& x);

    bool               is_valid  () const { return tet != 0; }
    const TetraCL&     get_tetra () const { return *tet; }
    const BaryCoordCL& get_bary  () const { return xb; }
    const Point3DCL&   get_point () const { return x; }
};

} // end of namespace DROPS

#endif

