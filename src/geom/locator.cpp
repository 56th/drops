/// \file locator.cpp
/// \brief 
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

#include "geom/locator.h"
#include "geom/multigrid.h"

namespace DROPS {

inline double
MyLocatorCL::min_bary_coeff (const BaryCoordCL& b) const
{
    return *std::min_element( b.begin(), b.end());
}

void
MyLocatorCL::locate_in_tetra ()
// Assumes, that x lies in tet and xb contains the barycentric coordinates of x in tet. If this prerequisite is not met, this function might
// loop forever or lie to you. You have been warned!
// Searches x in the children of tet up to triangulation-level lvl.
{
    double mincoeff;
    const TetraCL* mtet; // Tentative tetra for boundary cases.
    double mc; // Tentative mincoeff for boundary cases.
    BaryCoordCL mb; // Tentative xb for boundary cases.

    for (Uint l= tet->GetLevel(); l < lvl; ++l) {
        if (tet->IsUnrefined())
            break;

        mincoeff= -1.;
        mtet= 0;
        mc= -eps;
        for (TetraCL::const_ChildPIterator it= tet->GetChildBegin(), theend= tet->GetChildEnd(); it != theend; ++it) {
            w2b.assign( **it);
            xb= w2b( x);
            mincoeff= min_bary_coeff( xb);
            if (mincoeff >= 0. || (greedy == true && mincoeff >= -eps)) { // containment.
                tet= *it;
                break;
            }
            else if (mincoeff > mc) { // close to containent. Continue searching a better match.
                mc= mincoeff;
                mtet= *it;
                mb= xb;
            }
        }
        if (mincoeff < 0. && mtet != 0) { // Take the closest match for a containing tetra and try it.
            mincoeff= mc;
            tet= mtet;
            xb= mb;
        }
        if (mincoeff <= -eps)
            throw DROPSErrCL( "MyLocatorCL::LocateInTetra: Lost inclusion.\n");
    }
}

// This only works for triangulations of polygonal domains, which resolve the geometry of the domain exactly (on level 0).
void
MyLocatorCL::locate (const Point3DCL& xarg)
{
    x= xarg;
    tet= 0;

    const TetraCL* mtet= 0;
    double mc= -eps;
    BaryCoordCL mb;

    double mincoeff;
    DROPS_FOR_TRIANG_CONST_TETRA( mg, 0, it) {
        w2b.assign( *it);
        xb= w2b( x);
        mincoeff= min_bary_coeff( xb);
        if (mincoeff >= 0. || (greedy == true && mincoeff >= -eps)) { // containment.
            tet= &*it;
            break;
        }
        else if (mincoeff > mc) {// close to containent. Continue searching a better match.
            mc= mincoeff;
            mtet= &*it;
            mb= xb;
        }
    }
    if (tet == 0 && mtet != 0) {// Take the closest match for a containing tetra and try it.
        tet= mtet;
        xb= mb;
    }
    if (tet == 0)
        throw DROPSErrCL( "MyLocatorCL::Locate: Point not found on level 0.\n");

    locate_in_tetra();
}

} // end of namespace DROPS
