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

#include <sstream>
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
        if (mincoeff <= -eps) {
            std::ostringstream os;
            os << "MyLocatorCL::LocateInTetra: Lost inclusion; mincoeff: " << mincoeff << ".\n";
            throw DROPSErrCL( os.str());
        }
    }
}


void
StructuredCoarseTriangLocatorCL::compute_cubes (Point3DCL x)
{
        x-= orig;
        qr.Solve( x);
        for (Uint i= 0; i < 3; ++i) {
            const double q= x[i]*n[i],
                         m= std::floor( q),
                         M= std::ceil(  q);
            if (std::abs( q - m) < 1e-13) {
                idx[i][0]= m < n[i] ? m     : std::numeric_limits<Uint>::max();
                idx[i][1]= m > 0    ? m - 1 : std::numeric_limits<Uint>::max();
            }
            else if (std::abs( M - q) < 1e-13) {
                idx[i][0]= M > 0    ? M - 1 : std::numeric_limits<Uint>::max();
                idx[i][1]= M < n[i] ? M     : std::numeric_limits<Uint>::max();
            }
            else {
                idx[i][0]= m >= 0 && m < n[i] ? m : std::numeric_limits<Uint>::max();
                idx[i][1]= std::numeric_limits<Uint>::max();
            }
            if (idx[i][0] == std::numeric_limits<Uint>::max())
                std::swap( idx[i][0], idx[i][1]);
        }
}

bool
StructuredCoarseTriangLocatorCL::read_json (const ParamCL& P, const MultiGridCL& mg)
{
        ta.clear();

        if (P.get<std::string>("Type") != std::string("BrickBuilder")) {
            return false;
        }
        n[0]= P.get<Uint>( "N1");
        n[1]= P.get<Uint>( "N2");
        n[2]= P.get<Uint>( "N3");
        ta.resize( 6*n[0]*n[1]*n[2]);

        orig= P.get<Point3DCL>( "Origin");

        SMatrixCL<3, 3>& M= qr.GetMatrix();
        M.col( 0, P.get<Point3DCL>( "E1"));
        M.col( 1, P.get<Point3DCL>( "E2"));
        M.col( 2, P.get<Point3DCL>( "E3"));
        qr.prepare_solve();

        DROPS_FOR_TRIANG_CONST_TETRA( mg, 0, it) {
            const Point3DCL b= GetBaryCenter( *it);
            compute_cubes( b);
            const Uint pos= std::find( Addr( ta) + t_idx(idx[0][0], idx[1][0], idx[2][0], 0),
                                       Addr( ta) + t_idx(idx[0][0], idx[1][0], idx[2][0], 6),
                                       static_cast<const TetraCL*>( 0))
                            - Addr( ta);
             ta[pos]= &*it; // For the barycenters, the should be no ambiguity.
        }
        return true;
}

void
StructuredCoarseTriangLocatorCL::locate_tetras (const Point3DCL& x)
{
        tets.clear();
        compute_cubes( x);
        for (Uint i= 0; i < 2; ++i) {
            for (Uint j= 0; j < 2; ++j) {
                for (Uint k= 0; k < 2; ++k) {
                    if (idx[0][i] != std::numeric_limits<Uint>::max() &&
                        idx[1][j] != std::numeric_limits<Uint>::max() &&
                        idx[2][k] != std::numeric_limits<Uint>::max()    )
                        tets.insert( tets.end(),
                                     Addr( ta) + t_idx(idx[0][i], idx[1][j], idx[2][k], 0),
                                     Addr( ta) + t_idx(idx[0][i], idx[1][j], idx[2][k], 6));
                }
            }
        }
}


// This only works for triangulations of polygonal domains, which resolve the geometry of the domain exactly (on level 0).
void
MyLocatorCL::locate (const Point3DCL& xarg)
{
    x= xarg;
    tet= 0;

    if (structured_coarse_locator.is_initialized()) {
        locate_structured();
        locate_in_tetra();
        return;
    }

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
    if (tet == 0) {
        std::ostringstream os;
        os << "MyLocatorCL::locate: Point not found on level 0.; x: " << x << ".\n";
        throw DROPSErrCL( os.str());
    }

    locate_in_tetra();
}

void
MyLocatorCL::locate_structured ()
{
    const TetraCL* mtet= 0;
    double mc= -eps;
    BaryCoordCL mb;

    structured_coarse_locator.locate_tetras ( x);
    const std::vector<const TetraCL*>& tetras= structured_coarse_locator.get_tetras();
    double mincoeff;
    for (std::vector<const TetraCL*>::const_iterator it= tetras.begin(); it != tetras.end(); ++it) {
        w2b.assign( **it);
        xb= w2b( x);
        mincoeff= min_bary_coeff( xb);
        if (mincoeff >= 0. || (greedy == true && mincoeff >= -eps)) { // containment.
            tet= *it;
            break;
        }
        else if (mincoeff > mc) {// close to containent. Continue searching a better match.
            mc= mincoeff;
            mtet= *it;
            mb= xb;
        }
    }
    if (tet == 0 && mtet != 0) {// Take the closest match for a containing tetra and try it.
        tet= mtet;
        xb= mb;
    }
    if (tet == 0) {
        std::ostringstream os;
        os << "MyLocatorCL::locate_structured: Point not found on level 0.; x: " << x << ".\n";
        throw DROPSErrCL( os.str());
    }
if ((GetWorldCoord( *tet, xb) - x).norm() > 1e-12 || min_bary_coeff( xb) <= -eps) {
    std::cerr << "x: " << x << " xb: " << xb << " diff: " << (GetWorldCoord( *tet, xb) - x).norm() << std::endl;
}

}

} // end of namespace DROPS
