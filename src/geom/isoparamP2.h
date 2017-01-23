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
 * Copyright 2016 LNM/SC RWTH Aachen, Germany
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
 *  if \f$ \Phi \f$ denotes the transformation from \f$ T_ref \f$ to \f$ T \f$
 *  then the result is \f$ T = (\nabla \Phi)^{-T} \f$ and 
 *  det \f$ = det(\nabla \Phi) \f$
 *  Inputs are p,  the coordinates on the reference tetrahedron 
 *  and ct, the ten points on the curved tetrahedron        */
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

/** calculates the transpose of the transformation  Tetra -> RefTetra at
 *  each point (important for isoparametric elements).
 *  if \f$ \Phi \f$ denotes the transformation from \f$ T_ref \f$ to \f$ T \f$
 *  then the result is \f$ T = (\nabla \Phi)^{-T} \f$ and 
 *  det \f$ = det(\nabla \Phi) \f$
 *  Inputs are p,  the barycentric coordinates on the reference tetrahedron
 *  and ct, the ten points on the curved tetrahedron*/
inline void GetTrafoTrAtPoint( SMatrixCL<3,3>& T, double& det, const BaryCoordCL& p, const LocalP2CL<Point3DCL> & ct)
{
    double M[3][3];  
    for (int j=0; j<3; j++)  
        for (int k=0; k<3; k++)  
            M[j][k] = 0.;
    Point3DCL point(0.);
    Point3DCL gradphi(0.);
    for (int i=0; i<10; i++){
        gradphi = FE_P2CL::DHRef(i,p[1],p[2],p[3]);
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
template<class LocalFE >
inline Point3DCL GetWorldCoord( const LocalFE & ct, const Point3DCL& p)
{
    Point3DCL point(0.);
    double val;
    for (Uint i=0; i<LocalFE::FETYPE::NumDoFC; i++){
        val = LocalFE::FETYPE::H(i,p[0],p[1],p[2]);
        point += val * ct[i];
    }
    return point;
}

template<class LocalFE >
inline Point3DCL GetWorldCoord( const LocalFE & ct, const BaryCoordCL& p)
{
    Point3DCL point(0.);
    double val;
    for (Uint i=0; i<LocalFE::FETYPE::NumDoFC; i++){
        val = LocalFE::FETYPE::H(i,p[1],p[2],p[3]);
        point += val * ct[i];
    }
    return point;
}

template<class LocalFE, class QuadCL_double, class QuadCL_mat>
inline void GetTrafoAsQuad( const LocalFE & ct, QuadCL_double & adet, QuadCL_mat & T)
{
    for (Uint i=0; i<QuadCL_double::DataClass::NumNodesC; i++){
        const BaryCoordCL & b(QuadCL_double::DataClass::Node[i]);
        GetTrafoTrAtPoint( T[i], adet[i], b, ct);
        adet[i] = std::abs(adet[i]);
    }
}
//*************************************************************************************************
//*                                  Two dimensional transformation
//*************************************************************************************************

/** calculates the transpose of the transformation  Triangle -> RefTriangle at
 *  each point (important for isoparametric elements). 
 *  Inputs are p, the coordinates on the reference tetrahedron in which the triangle lies,
 *  and ct, the ten points on the curved tetrahedron in which the triangle lies,
 *  and the barycentric coordinates of the vertices of the triangle.  */
inline void Get2DTrafoTrAtPoint( SMatrixCL<3,3>& T, double& adet, Point3DCL& Nout, const Point3DCL& p, const LocalP2CL<Point3DCL> & ct, BaryCoordCL Bary[3])
{
    double M[3][3];  
    for (int j=0; j<3; j++)  
        for (int k=0; k<3; k++)  
            M[j][k] = 0.;
    Point3DCL point(0.);
    Point3DCL gradphi(0.);
    // Two orthonormal vectors on the reference face;
    Point3DCL s, t, v1, v2;
    for (int i=0; i<3; i++){
        s[i] = Bary[1][i+1] - Bary[0][i+1];
        t[i] = Bary[2][i+1] - Bary[0][i+1];
    }
    s /= s.norm();
    double temp= t[0]*s[0] + t[1]*s[1] + t[2]*s[2];
    t = t - temp * s;
    t /=t.norm();
    BaryCoordCL pt(p[0] * Bary [0] +  p[1] * Bary [1] +  p[2] * Bary [2]);
    for (int i=0; i<10; i++){
        gradphi = FE_P2CL::DHRef(i,pt[1],pt[2],pt[3]);
        point = ct[i];
        for (int j=0; j<3; j++)  
            for (int k=0; k<3; k++)  
                M[j][k] += point[j] * gradphi[k];
    }
    for (int i=0; i<3; i++){
        v1[i] = M[i][0]*s[0] + M[i][1]*s[1] + M[i][2]*s[2];
        v2[i] = M[i][0]*t[0] + M[i][1]*t[1] + M[i][2]*t[2];
    } 
        
    double det=   M[0][0] * (M[1][1]*M[2][2] - M[1][2]*M[2][1])
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
    cross_product(Nout, v1, v2);
    adet = Nout.norm();
    Nout /= adet; 
}

/** For isoparametric FEM, the transformation matrix and the determinant for the transformation from the reference triangle to a triangle is not a constant.
 *  This function computes the transformation matrix T, the determinant adet and also the outer normal vector on the triangle. */
template<class LocalFE, class Quad2DCL_double, class Quad2DCL_Point3DCL, class QuadCL_mat>
inline void Get2DTrafoAsQuad( const LocalFE & ct, BaryCoordCL Bary[3], Quad2DCL_Point3DCL& Nout, Quad2DCL_double & adet, QuadCL_mat & T)
{
    for (Uint i=0; i<Quad2DCL_double::DataClass::NumNodesC; i++){
        const Point3DCL & b(Quad2DCL_double::DataClass::Node[i]);
        Get2DTrafoTrAtPoint( T[i], adet[i], Nout[i], b, ct, Bary);
    }
}


} // end of namespace DROPS

#endif
