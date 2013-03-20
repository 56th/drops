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
 * Copyright 2013 LNM/SC RWTH Aachen, Germany
*/

#include "stokes.h"

namespace DROPS
{

//P2_P1 multiplication
void SpecialBndHandleSystem2OnePhaseCL::setupB(const TetraCL& tet, SMatrixCL<1, 3> loc_b[10][4])
{

	for (Uint k =0; k< 4; ++k) //Go throught all faces of a tet
	{
		LocalP2CL<double> phiP2[6];   //local basis for velocity
		LocalP1CL<double> phiP1[3];   //local basis for pressure
		Quad5_2DCL<double> mass2Dj;
		Quad5_2DCL<double> mass2Di;

		BaryCoordCL bary[3];
		if( BndData_.Vel.GetBC(*tet.GetFace(k))==SlipBC || BndData_.Vel.GetBC(*tet.GetFace(k))==SymmBC){
			const FaceCL& face = *tet.GetFace(k);
            double absdet = FuncDet2D(	face.GetVertex(1)->GetCoord()-face.GetVertex(0)->GetCoord(),
                                           	face.GetVertex(2)->GetCoord()-face.GetVertex(0)->GetCoord());  
			tet.GetOuterNormal(k, normal);
			for (Uint i= 0; i<3; ++i) //m is index for Vertex or Edge
			{
				unknownIdx[i]   = VertOfFace(k, i);
				unknownIdx[i+3] = EdgeOfFace(k, i) + 4;
				bary[i][unknownIdx[i]]=1;
				phiP1[i][unknownIdx[i]]=1;
			}

			for(Uint i=0; i<6; ++i)
				phiP2[i][unknownIdx[i]] = 1;
				
			for(Uint i=0; i<6; ++i){
				for(Uint j=0; j<3; ++j){					
					mass2Dj.assign(phiP2[i], bary);  //
					mass2Di.assign(phiP1[j], bary);
					Quad5_2DCL<double> mass2D(mass2Dj * mass2Di); //
					loc_b[unknownIdx[i]][unknownIdx[j]](0, 0)+= mass2D.quad(absdet)*normal[0];
					loc_b[unknownIdx[i]][unknownIdx[j]](0, 1)+= mass2D.quad(absdet)*normal[1];
					loc_b[unknownIdx[i]][unknownIdx[j]](0, 2)+= mass2D.quad(absdet)*normal[2];
				}
			}
		}
	}	
	
}

}