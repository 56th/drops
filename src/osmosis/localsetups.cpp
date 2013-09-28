
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
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

namespace DROPS
{


bool IsRegBaryCoord(const SArrayCL<BaryCoordCL,4>& T)
{
    for(int i= 0; i < 4; ++i)  
        for(int j= 0; j < 4; ++j) {
            if (/*std::isnan(T[i][j])|| std::isinf(T[i][j])|| */ T[i][j] >1.|| T[i][j] <0.) {
                std::cout << "Irregular coordinate!\n";
                return false;
        }
    }
    return true;
}
//====================================================
//
// Setup routines in each tetra
//
//====================================================

/// compute element vector for integrator: \n
/// \f$ f(v) = \int_T f v dx \f$ (for unstabilized fem)\n 
void ComputeRhsElementVector( double locf[4], LocalConvDiffReacCoefficients& local_cdcoef, TransformedP1FiniteElement& transformedfel)
{
    for(int i= 0; i < 4; ++i)
        locf[i]= Quad3CL<>(local_cdcoef.GetSourceAsQuad3()*transformedfel.GetGridfunctions().GetShapeAsQuad3(i)).quad(transformedfel.GetAbsDeterminant());
}

void SetupLocalTwoPhaseRhs(TransformedP1FiniteElement& transformedfel, InterfaceTetraCL& cut, 
    Elvec4& f, LocalConvDiffReacCoefficients& local_coefs) 
{
    cut.ComputeSubTets();
    
    for (Uint k=0; k< cut.GetNumTetra(); ++k){
        const SArrayCL<BaryCoordCL,4>& T =cut.GetTetra(k);
        if (!IsRegBaryCoord(T))  continue;
        double Vol = transformedfel.GetAbsDeterminant()*VolFrac(T);
        /*if (std::isnan(Vol) || std::isinf(Vol)){
            std::cout<<" SetupLocalTwoPhaseRhs: M Support of XFEM is too small.\t";
            continue;
        }*/
        transformedfel.SetSubTetra(T);
        Quad3CL<> qrhs(transformedfel.GetTetra(), local_coefs.GetSourceAsFuntionPointer(), local_coefs.GetTime(), transformedfel.GetNodes());  
        
        for(Uint i= 0; i < 4; ++i) {
            double tmp;
            tmp = Quad3CL<>(qrhs*transformedfel.GetBaseShapeAsQuad3CL(i)).quad(Vol);
            /*if (std::isnan(tmp) || std::isinf(tmp)) {
                for( Uint j=0; j<4; ++j)
                    f[j]= 0.;
                break;
            }*/
            if (k< cut.GetNumNegTetra())
                f[i]+= tmp;
        }
    }     
}


/// compute the following element matrices for an element which is completely in one phase only: \n
/// \f$ M(u,v) = \int u v \f$ with h inverse of the henry coefficient of the domain \n
/// \f$ A(u,v) = \int d \nabla u \nabla v \f$ with d the diffusion coefficient \n
/// \f$ C(u,v) = 0 \n
/// New variable: pointsInTime: -1 - only old time, 0 - both, 1 - only new time
void SetupLocalOnePhaseSystem(TransformedP1FiniteElement& transformedfel, ConvDiffElementMatrices & elmats, LocalConvDiffReacCoefficients& local_coefs,	int pointsInTime)
{
    //const SMatrixCL<3,4> G = transformedfel.GetDShape();
    const SMatrixCL<4,4> GTG = transformedfel.GetGramShape();
    const double absdet = transformedfel.GetAbsDeterminant();
    //Only one phase, so diffusion coeff always of negative part
    const double d = local_coefs.GetDiffusionCoef();
	
    //const Quad3CL<Point3DCL> & q3_u = local_coefs.GetVelocityAsQuad3();
    for(int i= 0; i < 4; ++i) {
        //Quad3CL<> & phi_i = transformedfel.GetGridfunctions().GetShapeAsQuad3(i);
        for(int j= 0; j < i; ++j) {
            // Quad3CL<> & phi_j = transformedfel.GetGridfunctions().GetShapeAsQuad3(j);
        	if (pointsInTime >= 0) {
        		elmats.M[j][i]= P1DiscCL::GetMass( i, j)*absdet;
        		elmats.M[i][j]= elmats.M[j][i];
        		elmats.A[j][i]= d*GTG( i, j)/6.0*absdet;
        		elmats.A[i][j]= elmats.A[j][i];
        	}
        	if (pointsInTime <= 0) {
        		elmats.C[i][j]= d*GTG( i, j)/6.0*absdet;
            	elmats.C[j][i]= elmats.C[i][j];
        	}
        }
        if (pointsInTime >= 0) {
            elmats.M[i][i]= P1DiscCL::GetMass( i, i)*absdet;
            elmats.A[i][i]= d*GTG( i, i)/6.0*absdet;
        }
        if (pointsInTime <= 0) {
            elmats.C[i][i]= d*GTG( i, i)/6.0*absdet;
        }
    }

}

// Couplings between the basis functions (at the same time step or different time steps) when the tetra is cut by only one interface.
// pointsInTime: -1 - only old, 1 - only new
// imho no scenario for BOTH
void SetupLocalOneInterfaceSystem( TransformedP1FiniteElement& transformedfel, InterfaceTetraCL& cut, 
      ConvDiffElementMatrices & elmats, LocalConvDiffReacCoefficients& local_coefs, int pointsInTime)
{
    //const SMatrixCL<3,4> G = transformedfel.GetDShape();
    const SMatrixCL<4,4> GTG = transformedfel.GetGramShape();
    const double absdet = transformedfel.GetAbsDeterminant();

    double d = local_coefs.GetDiffusionCoef();

    cut.ComputeSubTets();
    Uint NumTets=cut.GetNumNegTetra(); //# of neg subtetras
    
    for (Uint k=0; k< NumTets; ++k){
           
        const SArrayCL<BaryCoordCL,4>& T =cut.GetTetra(k);
        
        if (!IsRegBaryCoord(T))  continue;
        
        transformedfel.SetSubTetra(T);
        
        double Vol = absdet*VolFrac(T);
        /*if (std::isnan(Vol)|| std::isinf(Vol)){
            std::cout<<"Vol " <<VolFrac(T)<<"\n";
            std::cout<<"SetupLocalOneInterfaceSystem: Support is too small.\t";
            continue;
        }*/
        Quad3CL<Point3DCL> q3_u(local_coefs.GetVelocityAsLocalP2(), transformedfel.GetNodes());
        bool irreg = false; 

        if (pointsInTime == 0)
        	std::cout << "SetupLocalOneInterfaceSystem: pointsInTime = 0: CHECK" << std::endl;

        for(int i= 0; i < 4; ++i) {
            //Quad3CL<> qp1(transformedfel.GetGridfunctions().GetShapeAsLocalP1(i), transformedfel.GetNodes());
            for(int j= 0; j < 4; ++j) {
              
                Quad3CL<> qM(transformedfel.GetGridfunctions().GetProductShapeAsLocalP2(i,j),transformedfel.GetNodes());
                double iM, iA, iC;
                if (pointsInTime > 0) {
                	iM = qM.quad(Vol);
                	iA = Vol* GTG( j, i)/6. * d;
                	iC = 0.0;
                }
                else {
                	iM = 0.0;
                	iA = 0.0;
                	iC = Vol* GTG( j, i)/6. * d;
                }
              
                /*if (std::isnan(iM)|| std::isinf(iM)||std::isnan(iA)|| std::isinf(iA)||std::isnan(iC)|| std::isinf(iC)) {
                    elmats.ResetAll();
                    irreg = true;
                    break;
                }*/

                elmats.M[i][j]+= iM;
                elmats.A[i][j]+= iA;
                elmats.C[i][j]+= iC;

            }
            if (irreg) break;
        }
    }

}


void SetupLocalOnePhaseMassMatrix( double locM[4][4], TransformedP1FiniteElement & transfp1fel) //, bool pPart)
{
    //const double h = pPart ? 1. : 1./H;
    
    for(int i= 0; i < 4; ++i) {
        for(int j= 0; j < i; ++j) {
            locM[j][i]= P1DiscCL::GetMass( i, j)*transfp1fel.GetAbsDeterminant();
            locM[i][j]= locM[j][i];
        }
        locM[i][i]= P1DiscCL::GetMass( i, i)*transfp1fel.GetAbsDeterminant();
    }
}

/// compute the mass matrix for the case that the tetrahedron is cut ONLY for the old OR new time. \n
/// computes the off-diagonal block of \f$ M(u,v) = \int h u v \f$ with h inverse of the henry coefficient of the domain \n
/// \f$ M_{1,2}[i,j] = \int h \phi_i^{new} \phi_j^{old} \f$ \n
/// or
/// \f$ M_{2,1}[i,j] = \int h \phi_i^{old} \phi_j^{new} \f$ \n
/// OSMOSIS: Only M_neg needed!
/*void SetupLocalOneInterfaceMassMatrix( InterfaceTetraCL& cut, double M_n[4][4], double M_p[4][4],
    TransformedP1FiniteElement & transfp1fel, const double H, bool sign[4], bool jumps, bool pPart)*/
void SetupLocalOneInterfaceMassMatrix( InterfaceTetraCL& cut, double M[4][4], TransformedP1FiniteElement & transfp1fel)
{
    cut.ComputeSubTets();
    Uint NumNegTets=cut.GetNumNegTetra(); /// # of neg subtetras
    // D, H are constant if T \cap Gamma_old is empty
    //double h = pPart ? 1. : 1./H;
    std::memset( M,0, 4*4*sizeof(double));
    
 //   bool sign[4];

 //   for(int i= 0; i < 4; ++i) {
 //       sign[i]= (cut.GetSign(i) == 1);
 //   }

    for (Uint k=0; k< NumNegTets; ++k){
        const SArrayCL<BaryCoordCL,4>& T =cut.GetTetra(k);
        transfp1fel.SetSubTetra(T);
        if (!IsRegBaryCoord(T)) continue;
        
        double Vol = transfp1fel.GetAbsDeterminant()*VolFrac(T);
        /*if (std::isnan(Vol)|| std::isinf(Vol)){
            std::cout<<"Vol " <<VolFrac(T)<<"\n";
            std::cout<<" Support of XFEM is too small.\t";
            continue;
        }*/
        bool irreg=false;
        for(int i= 0; i < 4; ++i) {
            for(int j= 0; j < 4; ++j) {
                double iM = 0;
                //Quad3CL<> qM(transfp1fel.GetBaseShapeAsQuad3CL(j)*transfp1fel.GetBaseShapeAsQuad3CL(i));
                Quad3CL<> qM(transfp1fel.GetGridfunctions().GetProductShapeAsLocalP2(i,j),transfp1fel.GetNodes());
                iM = qM.quad(Vol);
                /*if (std::isnan(iM)|| std::isinf(iM)) {
                    std::memset( M,0, 4*4*sizeof(double));
                    irreg=true;
                    break; //leave the j-loop 
                }*/

                M[i][j]+= iM;
            }
            if (irreg) break;  //leave the i-loop
        }
    }
}

/// compute the mass matrix for the case that the tetrahedron is cut for the old and new time. \n
/// computes \f$ M(u,v) = \int u v \f$ of the domain (at new time) \n
/// note that \f$ \phi_k \f$  are basis function for the transformed unknowns H*u.
void SetupLocalTwoInterfacesMassMatrix( InterfaceTetraCL& cut, __UNUSED__ InterfaceTetraCL& oldcut, 
    double M22[4][4], TransformedP1FiniteElement & transfp1fel, LocalP2CL<>& lp2_oldlset)
{
    std::memset( M22, 0, 4*4*sizeof(double));
    
    double absdet = transfp1fel.GetAbsDeterminant();

    cut.ComputeSubTets();
    Uint NumTets=cut.GetNumNegTetra(); /// # of neg. subtetras

    for (Uint k=0; k< NumTets; ++k){
        const SArrayCL<BaryCoordCL,4>& T =  cut.GetTetra(k);
        if (!IsRegBaryCoord(T)) continue;
        
        transfp1fel.SetSubTetra(T);
        
        InterfaceTetraCL suboldcut;
        suboldcut.Init(T, lp2_oldlset,0.);

        double VolT = absdet*VolFrac(T);

        bool nPartOld = false;
		for (int i=0; i < 9; i++)
			if (suboldcut.GetSign( i) != 1)
				nPartOld = true;

		if (nPartOld)
		{
			// Meinem Verstaendis nach: Sowohl der Subtetra der einen Levelset als auch der anderen sind
			// im negativen Bereich -> Tet gehoert zum Integrationsgebiet!
			for(int i= 0; i < 4; ++i) {
				bool irreg = false;
				for(int j= 0; j < 4; ++j) {
					double iM = 0;
					Quad3CL<> qM(transfp1fel.GetBaseShapeAsQuad3CL(i)*transfp1fel.GetBaseShapeAsQuad3CL(j));
					iM = qM.quad(VolT);

					/*if (std::isnan(iM)|| std::isinf(iM)) {
					///> TODO: If a local value in a child tetrahedron is irregular, ignore this tetra
					///> The contribution of other tetra must be preserved
					///> For example tmpM21[4][4]
						irreg=true;
						break; //leave the j-loop
					}*/

					M22[i][j]+= iM;
				}
				if (irreg) break;
			}
		}
    }
}

Quad3CL<Point3DCL> GetImprovedTetraNormal(TetraCL& t, VecDescCL& lsetphi, LsetBndDataCL& lsbnd)
{
	//TODO: GradRef static?
	LocalP1CL<Point3DCL> GradRef[10], GradP2[10];
	P2DiscCL::GetGradientsOnRef( GradRef);
	SMatrixCL<3,3> T;
	//double dummy;
	//GetTrafoTr( T, dummy, tet);

	//Transoformation matrix, function: GetTrafoTr
	//but without needing the tetra
	double M[3][3];
	const Point3DCL& pt0= t.GetVertex(0)->GetCoord();
	for(int i=0; i<3; ++i)
		for(int j=0; j<3; ++j)
			M[j][i]= t.GetVertex(i+1)->GetCoord()[j] - pt0[j];
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

	Quad3CL<Point3DCL> GradQ[10]; // and their gradients
	P2DiscCL::GetGradients( GradP2, GradRef, T);
	Quad3CL<Point3DCL> normal;
	LocalP2CL<> lset(t, lsetphi, lsbnd);
	for (int v=0; v<10; ++v)
	{
		GradQ[v].assign( GradP2[v]);
		normal += lset[v]*GradQ[v];
	}

	for (int i =0; i<Quad3DataCL::NumNodesC; i++) if (normal[i].norm()>1e-8) normal[i]/= normal[i].norm();

	return normal;
}

}//end of namespace DROPS
