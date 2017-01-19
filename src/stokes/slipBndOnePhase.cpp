/// \file slipBndOnephase.cpp
/// \brief classes that add Nitsche terms for problems with slip or symmetry boundary segments in one phase.
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

#include "slipBndOnePhase.h"
#include "stokes/stokes.h"

namespace DROPS{

//***********************************************************************
//                  SlipBndSystem1OnePhaseP2CL
//***********************************************************************
void SlipBndSystem1OnePhaseP2CL::setup(const TetraCL& tet, const SMatrixCL<3,3>& T, SMatrixCL<3,3> Ak[10][10])
{ 
    LocalP1CL<Point3DCL> Grad[10];
    P2DiscCL::GetGradients( Grad, GradRef, T);
    for (Uint k =0; k< 4; ++k) //Go through all faces of a tet
    {
        Point3DCL normal;
        SMatrixCL<3, 3> dm[10][10];
        LocalP2CL<double> phi[10]; 
        Quad5_2DCL<double> mass2Dj;
        Quad5_2DCL<double> mass2Di;
        Quad5_2DCL<double> Grad2Dj;   // \nabla phi_j* n
        Quad5_2DCL<double> Grad2Di;   // \nabla phi_i* n
        BaryCoordCL bary[3];
        bool symmBC=false;
        if( BndData_.Vel.IsOnSlipBnd(*tet.GetFace(k)) || BndData_.Vel.IsOnSymmBnd(*tet.GetFace(k))){
            if(BndData_.Vel.IsOnSymmBnd(*tet.GetFace(k)))
                symmBC=true;
            const FaceCL& face = *tet.GetFace(k);
            double absdet = FuncDet2D(face.GetVertex(1)->GetCoord()-face.GetVertex(0)->GetCoord(),
                                      face.GetVertex(2)->GetCoord()-face.GetVertex(0)->GetCoord()); 
            double h= std::sqrt(absdet);
            tet.GetOuterNormal(k, normal);
            for (Uint i= 0; i<3; ++i)    
            {
                bary[i][VertOfFace(k, i)]=1;
            }
            for(Uint i=0; i<10; ++i)
                phi[i][i] = 1;
                
            double temp = symmBC? 0.: beta_;
            for(Uint i=0; i<10; ++i){
                LocalP1CL<double> Gradin(dot( normal, Grad[i]));
                mass2Di.assign(phi[i], bary);
                Grad2Di.assign(Gradin, bary);
                for(Uint j=0; j<=i; ++j){
                    LocalP1CL<double> Gradjn(dot( normal, Grad[j]));
                    mass2Dj.assign(phi[j], bary); 
                    Grad2Dj.assign(Gradjn, bary);
                    Quad5_2DCL<double> mass2D(mass2Dj * mass2Di); 
                    Quad5_2DCL<double> Grad2D(Grad2Di * mass2Dj + Grad2Dj * mass2Di); 
                    // three additional terms
                    dm[j][i](0, 0)=dm[j][i](1, 1) = dm[j][i](2, 2) = temp * mass2D.quad(absdet);
                    dm[j][i]     += (alpha_/h * mu_ - temp) * mass2D.quad(absdet) * SMatrixCL<3,3> (outer_product(normal, normal));
                    // if(BndData_.Vel.GetBC(*tet.GetFace(k))!= SymmBC) not necessary
                    dm[j][i]     -=  2. * mu_ * Grad2D.quad(absdet) * SMatrixCL<3,3> (outer_product(normal, normal));  
                    Ak[j][i] += dm[j][i];
                    if (i != j){
                        assign_transpose( dm[i][j], dm[j][i]);
                        Ak[i][j] += dm[i][j];
                    }	
                }
            }
        }
    }
}

void SlipBndSystem1OnePhaseP2CL::setupRhs(const TetraCL& tet, Point3DCL loc_b[10], double t)
{
    for (Uint k =0; k< 4; ++k){ //Go throught all faces of a tet
        Point3DCL normal;
        Uint unknownIdx[6];
        LocalP2CL<double> phi[6]; 
        Quad5_2DCL<double> locP2[6];
        Quad5_2DCL<Point3DCL> WallVel;
        BaryCoordCL bary[3];
        if( BndData_.Vel.IsOnMovSlipBnd(*tet.GetFace(k))){
            const FaceCL& face = *tet.GetFace(k);
            double absdet = FuncDet2D(face.GetVertex(1)->GetCoord()-face.GetVertex(0)->GetCoord(),
                                      face.GetVertex(2)->GetCoord()-face.GetVertex(0)->GetCoord()); 
            tet.GetOuterNormal(k, normal);
            for (Uint i= 0; i<3; ++i)    
            {
                unknownIdx[i]   = VertOfFace(k, i);      // i is index for Vertex
                unknownIdx[i+3] = EdgeOfFace(k, i) + 4;  // i is index for Edge
                bary[i][unknownIdx[i]]=1;
            }
            typedef StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
            bnd_val_fun bf= BndData_.Vel.GetBndSeg(face.GetBndIdx()).GetBndFun();
            WallVel.assign(tet, bary, bf, t);
            for(Uint i=0; i<6; ++i)
                phi[i][unknownIdx[i]] = 1;
                
            for(Uint i=0; i<6; ++i){
                    locP2[i].assign(phi[i], bary);
            }

            for(Uint i=0; i<6; ++i){//setup right hand side	
                    Quad5_2DCL<Point3DCL> WallVelRhs( locP2[i]* WallVel);
                    loc_b[unknownIdx[i]] += (beta_ * WallVelRhs.quad(absdet) - beta_ * SMatrixCL<3,3> (outer_product(normal, normal)) * WallVelRhs.quad(absdet));
            }
        }
    }
}

//***********************************************************************
//                  SlipBndSystem2OnePhaseCL
//***********************************************************************
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
        
        if( BndData_.Vel.IsOnSlipBnd(*tet.GetFace(k)) || BndData_.Vel.IsOnSymmBnd(*tet.GetFace(k))){
            const FaceCL& face = *tet.GetFace(k);          
            double absdet = FuncDet2D(face.GetVertex(1)->GetCoord()-face.GetVertex(0)->GetCoord(),
                                      face.GetVertex(2)->GetCoord()-face.GetVertex(0)->GetCoord());
            tet.GetOuterNormal(k, normal);
            for (Uint i= 0; i<3; ++i)  
            {
                unknownIdx[i]   = VertOfFace(k, i);
                unknownIdx[i+3] = EdgeOfFace(k, i) + 4;
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
                    loc_b[unknownIdx[i]][unknownIdx[j]] -= quad2D.quad(absdet)*SMatrixCL<1,3>(normal); 
                }
            }
        }
    }

}

}
