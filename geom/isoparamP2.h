/// \file isoparamP2.h
/// \brief just (inlined) helper function for isoparametrically(P2) curved tets
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
#include "num/discretize.h"

namespace DROPS
{

/** calculates the transpose of the transformation  Tetra -> RefTetra at
 *  each point (important for isoparametric elements).
 *  if \f$ \Phi \f$ denotes the trafo from \f$ T_ref \f$ to \f$ T \f$
 *  then the result is \f$ T = (\nabla \Phi)^{-T} \f$ and 
 *  det \f$ = det(\nabla \Phi) \f$
 *  Inputs are p the point on the reference triangle and pt, the TEN points 
 *  of the second order curved tetraeder as a CurvedTetraCL.  */
inline void GetTrafoTrAtPoint( SMatrixCL<3,3>& T, double& det, const Point3DCL& p, const LocalP2CL<Point3DCL> & ct)
{
    double M[3][3];  
    for (int j=0; j<3; j++)  
        for (int k=0; k<3; k++)  
            M[j][k] = 0.;
    Point3DCL point(0.);
    Point3DCL gradphi(0.);
    for (int i=0; i<10; i++){
        gradphi = FE_P2CL::DHRef(i,p[0],p[1],p[2]);
        point = ct[i];
        for (int j=0; j<3; j++)  
            for (int k=0; k<3; k++)  
                M[j][k] += point[j] * gradphi[k];
    }
    
    det=   M[0][0] * (M[1][1]*M[2][2] - M[1][2]*M[2][1])
        - M[0][1] * (M[1][0]*M[2][2] - M[1][2]*M[2][0])
        + M[0][2] * (M[1][0]*M[2][1] - M[1][1]*M[2][0]);

    T(0,0)= (M[1][1]*M[2][2] - M[1][2]*M[2][1])/det;
    T(0,1)= (M[2][0]*M[1][2] - M[1][0]*M[2][2])/det;
    T(0,2)= (M[1][0]*M[2][1] - M[2][0]*M[1][1])/det;
    T(1,0)= (M[2][1]*M[0][2] - M[0][1]*M[2][2])/det;
    T(1,1)= (M[0][0]*M[2][2] - M[2][0]*M[0][2])/det;
    T(1,2)= (M[2][0]*M[0][1] - M[0][0]*M[2][1])/det;
    T(2,0)= (M[0][1]*M[1][2] - M[1][1]*M[0][2])/det;
    T(2,1)= (M[1][0]*M[0][2] - M[0][0]*M[1][2])/det;
    T(2,2)= (M[0][0]*M[1][1] - M[1][0]*M[0][1])/det;      
}

/** calculates the physical coord of a reference point of a curved tetra */
inline Point3DCL GetWorldCoord( const LocalP2CL<Point3DCL> & ct, const Point3DCL& p)
{
    Point3DCL point(0.);
    double val;
    for (int i=0; i<10; i++){
        val = FE_P2CL::H(i,p[0],p[1],p[2]);
        point += val * ct[i];
    }
    return point;
}

inline Point3DCL GetWorldCoord( const LocalP2CL<Point3DCL> & ct, const BaryCoordCL& p)
{
    std::cout << " here " << std::endl;
    Point3DCL point(0.);
    double val;
    for (int i=0; i<10; i++){
        val = FE_P2CL::H(i,p[1],p[2],p[3]);
        point += val * ct[i];
    }
    return point;
}


} // end of namespace DROPS

#endif
