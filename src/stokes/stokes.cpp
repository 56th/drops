/// \file stokes.cpp
/// \brief classes that constitute the stokes-problem
/// \author LNM RWTH Aachen: Liang Zhang

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

#include "stokes.h"

namespace DROPS
{
/// Setup the integral of (bv * bn) * q on the slip bounary for uncut element 
void SlipBndSystem2OnePhaseCL::setupB(const TetraCL& tet, SMatrixCL<1, 3> loc_b[10][4])
{

    for (Uint k =0; k< 4; ++k) //Go throught all faces of a tet
    {
        LocalP2CL<double> phiVelP2[6];   //local basis for velocity
        LocalP1CL<double> phiPrP1[3];    //local basis for pressure
        Quad5_2DCL<double> pr2Dj;
        Quad5_2DCL<double> vel2Di;
        BaryCoordCL bary[3];
        Point3DCL normal;
        Uint unknownIdx[6];
        
        if( BndData_.Vel.GetBC(*tet.GetFace(k))==Slip0BC || BndData_.Vel.GetBC(*tet.GetFace(k))==SlipBC || BndData_.Vel.GetBC(*tet.GetFace(k))==SymmBC){
            const FaceCL& face = *tet.GetFace(k);            //Get a face on a slip boundary 
            double absdet = FuncDet2D(face.GetVertex(1)->GetCoord()-face.GetVertex(0)->GetCoord(),
                                      face.GetVertex(2)->GetCoord()-face.GetVertex(0)->GetCoord());
            tet.GetOuterNormal(k, normal);
            for (Uint i= 0; i<3; ++i)  
            {
                unknownIdx[i]   = VertOfFace(k, i);          // index for Vertex
                unknownIdx[i+3] = EdgeOfFace(k, i) + 4;      // ndex for Edge
                bary[i][unknownIdx[i]]=1;
                phiPrP1[i][unknownIdx[i]]=1;
            }

            for(Uint i=0; i<6; ++i)
                phiVelP2[i][unknownIdx[i]] = 1;
                
            for(Uint i=0; i<6; ++i){
                vel2Di.assign(phiVelP2[i], bary);  
                for(Uint j=0; j<3; ++j){
                    pr2Dj.assign(phiPrP1[j], bary);  
                    Quad5_2DCL<double> quad2D(pr2Dj * vel2Di);
                    //const double qq= quad2D.quad(absdet);
                    //loc_b[unknownIdx[i]][unknownIdx[j]](0, 0)-= qq*normal[0];
                    //loc_b[unknownIdx[i]][unknownIdx[j]](0, 1)-= qq*normal[1];
                    //loc_b[unknownIdx[i]][unknownIdx[j]](0, 2)-= qq*normal[2];  
                    loc_b[unknownIdx[i]][unknownIdx[j]] -= quad2D.quad(absdet)*SMatrixCL<1,3>(normal); 
                }
            }
        }
    }

}


void SlipBndSystem2OnePhaseCL::setupB(const TetraCL& tet, SMatrixCL<1, 3> loc_b[10][4], instat_vector_fun_ptr bnd_normal) 
///< curved case
{
    MeshDeformationCL& md = MeshDeformationCL::getInstance();
    Quad5_2DCL<double> adet;
    Quad5_2DCL< SMatrixCL<3, 3> > Tq;
    for (Uint k =0; k< 4; ++k) //Go throught all faces of a tet
    {
        LocalP2CL<double> phiVelP2[6];   //local basis for velocity
        LocalP1CL<double> phiPrP1[3];    //local basis for pressure
        Quad5_2DCL<double> pr2Dj;
        Quad5_2DCL<double> vel2Di;
        Quad5_2DCL<Point3DCL> Nout;
        Point3DCL normal;
        Uint unknownIdx[6];

        BaryCoordCL bary[3];
        if( BndData_.Vel.GetBC(*tet.GetFace(k))==Slip0BC || BndData_.Vel.GetBC(*tet.GetFace(k))==SlipBC || BndData_.Vel.GetBC(*tet.GetFace(k))==SymmBC){
            const FaceCL& face = *tet.GetFace(k);            //Get a face on a special boundary 
            double absdet = FuncDet2D(face.GetVertex(1)->GetCoord()-face.GetVertex(0)->GetCoord(),
                                      face.GetVertex(2)->GetCoord()-face.GetVertex(0)->GetCoord());
            tet.GetOuterNormal(k, normal);
            for (Uint i= 0; i<3; ++i){
                unknownIdx[i]   = VertOfFace(k, i);          // i is index for Vertex
                unknownIdx[i+3] = EdgeOfFace(k, i) + 4;      // i is index for Edge
                bary[i][unknownIdx[i]]=1;
                phiPrP1[i][unknownIdx[i]]=1;
            }
            if(md.IsTetraCurved(tet)){
                Get2DTrafoAsQuad(md.GetLocalP2Deformation(tet), bary, Nout, adet, Tq);
            }
            else{
                adet = Quad5_2DCL<double> (absdet);
                Nout = Quad5_2DCL<Point3DCL>(normal);	
            }
            //Nout.assign(tet, bary, bnd_normal);
            for(Uint i=0; i<6; ++i)
                phiVelP2[i][unknownIdx[i]] = 1;
            for(Uint i=0; i<6; ++i){
                vel2Di.assign(phiVelP2[i], bary);  
                for(Uint j=0; j<3; ++j){
                    pr2Dj.assign(phiPrP1[j], bary);  
                    Quad5_2DCL<double> temp(pr2Dj * vel2Di); 
                    Quad5_2DCL<Point3DCL> quad2D( temp * Nout * adet);
                    loc_b[unknownIdx[i]][unknownIdx[j]]-= SMatrixCL<1,3>( quad2D.quad(1.));
                }
            }
        }
    }
}

}