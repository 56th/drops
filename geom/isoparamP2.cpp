/// \file isoparamP2.cpp
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

#include "geom/simplex.h"
#include "geom/isoparamP2.h"
#include "num/fe.h"

namespace DROPS
{

// describes ordering of edges w.r.t. vertices
static int between[6][2]=
    {{0,1},
     {0,2},
     {1,2},
     {0,3},
     {1,3},
     {2,3}};

CurvedTetraCL::~CurvedTetraCL(){
    if (extrapoint_) delete [] extrapoint_;
}

CurvedTetraCL::CurvedTetraCL(const TetraCL& tetra, Point3DCL* new_extrapoints):uncurved_(tetra){
    extrapoint_ = new Point3DCL[6];
    Point3DCL tmp(0.);
    Point3DCL diff(0.);
    for (int i = 0; i < 6; i++){
        tmp = 0.5 * uncurved_.GetVertex(between[i][0])->GetCoord();
        tmp += 0.5 * uncurved_.GetVertex(between[i][1])->GetCoord();
        //find closest of new_extrapoints
        double mindist=1e99;
        int argmin=-1;
        double normdiff = 1e99;
        for (int j = 0; j < 6; j++){
            diff = tmp - new_extrapoints[j];
            normdiff = diff.norm();
            if (mindist>normdiff){
                argmin = j;
                mindist = normdiff;
            }
        }
        if (argmin == -1) 
            throw DROPSErrCL( "CurvedTetraCL::Constructor - no point has a minimal distance?!");          
        extrapoint_[i] = new_extrapoints[argmin];
    }
}
    
const Point3DCL& CurvedTetraCL::GetPoint(int i) const{
    if (i<4){
        return uncurved_.GetVertex(i)->GetCoord(); 
    }
      
    if (i<10){
        return extrapoint_[i-4];
    }else {
        throw DROPSErrCL( "CurvedTetraCL::GetPoint(i) , i should be smaller than 10!");          
    }
}

Point3DCL GetWorldCoord( const CurvedTetraCL & ct, const Point3DCL& p)
{
    Point3DCL point(0.);
    double val;
    for (int i=0; i<10; i++){
        val = FE_P2CL::H(i,p[0],p[1],p[2]);
        point += val * ct.GetPoint(i);
    }
    return point;
}

Point3DCL GetWorldCoord( const CurvedTetraCL & ct, const BaryCoordCL& p)
{
    std::cout << " here " << std::endl;
    Point3DCL point(0.);
    double val;
    for (int i=0; i<10; i++){
        val = FE_P2CL::H(i,p[1],p[2],p[3]);
        point += val * ct.GetPoint(i);
    }
    return point;
}

void GetTrafoTrAtPoint( SMatrixCL<3,3>& T, double& det, const Point3DCL& p, const CurvedTetraCL & ct)
{
    double M[3][3];  
    for (int j=0; j<3; j++)  
        for (int k=0; k<3; k++)  
            M[j][k] = 0.;
    Point3DCL point(0.);
    Point3DCL gradphi(0.);
    for (int i=0; i<10; i++){
        gradphi = FE_P2CL::DHRef(i,p[0],p[1],p[2]);
        point = ct.GetPoint(i);
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


} // end of namespace DROPS
