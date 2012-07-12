/// \file isoparamP2.h
/// \brief classes and helper function for isoparametrically(P2) curved tets
/// \author LNM RWTH Aachen: Christoph Lehrenfeld, Liang Zhang

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

#ifndef DROPS_ISOPARAMP2_H
#define DROPS_ISOPARAMP2_H

//#include "geom/simplex.h"
#include "misc/container.h"

namespace DROPS
{

class TetraCL;

// Description of a isoparametrically curved tetrahedra. 
// Stores reference to uncurved tetrahedra and (6) additional
// points.
class CurvedTetraCL
{
  public:
    const TetraCL& uncurved_;   // < original tetrahedra
    Point3DCL* extrapoint_;
    ~CurvedTetraCL();
    ///puts tetra (uncurved) and additional points together and (important) 
    ///makes the ordering fitting to the standard local numbering (see static int between)
    CurvedTetraCL(const TetraCL& tetra, Point3DCL* new_extrapoints);    
    const Point3DCL& GetPoint(int i) const;    
};

/** calculates the transpose of the transformation  Tetra -> RefTetra at
 *  each point (important for isoparametric elements).
 *  if \f$ \Phi \f$ denotes the trafo from \f$ T_ref \f$ to \f$ T \f$
 *  then the result is \f$ T = (\nabla \Phi)^{-T} \f$ and 
 *  det \f$ = det(\nabla \Phi) \f$
 *  Inputs are p the point on the reference triangle and pt, the TEN points 
 *  of the second order curved tetraeder as a CurvedTetraCL. Note that the
 *  points are ordered such that
 *        4 lies in between 0 and 1
 *        5                 0     2
 *        6                 1     2
 *        7                 0     3
 *        8                 1     3
 *        9                 2     3
 */
void GetTrafoTrAtPoint( SMatrixCL<3,3>& T, double& det, const Point3DCL& p, const CurvedTetraCL & ct);

/** calculates the physical coord of a reference point of a curved tetra
 *  It is important the the points in pt have the right numbering, i.e.
 *  they have to be planar-consistent, which means that for planar elements the
 *  following must hold:
 *  point 4 lies in between 0 and 1
 *        5                 0     2
 *        6                 1     2
 *        7                 0     3
 *        8                 1     3
 *        9                 2     3
 */
Point3DCL GetWorldCoord( const CurvedTetraCL & ct, const Point3DCL& p);
Point3DCL GetWorldCoord( const CurvedTetraCL & ct, const BaryCoordCL& p);

} // end of namespace DROPS

#endif
