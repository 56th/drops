/// \file osmosisSetup.cpp
/// \brief Classes that constitute a 1-phase-osmosis-problem.
/// \author Trung Hieu Nguyen, Thorolf Schulte, Christoph Lehrenfeld IGPM

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
#include "osmosis/osmosisSetup.h"
#include "osmosis/localsetups.cpp"
#include "osmosis/tests.cpp"
#include <iostream>
#include <fstream>

enum TZONES { ONLYOLD = -1, BOTH = 0, ONLYNEW = 1 };

namespace DROPS
{
//====================================================
//
// OsmosisP1CL
//
//====================================================

/// Initialize P1 function 
void OsmosisP1CL::Init (instat_scalar_fun_ptr cn, double t)
{
    conc.t = t;
    oldconc.t = t;
    const Uint lvl= conc.GetLevel(),
        concidx= conc.RowIdx->GetIdx(),
        oldconcidx= oldconc.RowIdx->GetIdx();


    DROPS_FOR_TRIANG_VERTEX( MG_, lvl, it) {
        if (it->Unknowns.Exist( concidx))
            conc.Data[it->Unknowns( concidx)]= cn( it->GetCoord(), t);
        if (it->Unknowns.Exist( oldconcidx))
            oldconc.Data[it->Unknowns( oldconcidx)]= cn( it->GetCoord(), t);
    }

	DoInterfaceStep(t);
    
}

void OsmosisP1CL::InitWithScaling (instat_scalar_fun_ptr cn, double totalconcentration, double t)
{
    Init(cn,t);

    double absdet;
    InterfaceTetraCL patch;

    double con = 0.0;
    LocalP2CL<double> ones( 1.);
    const Uint lvl= conc.GetLevel();
    DROPS_FOR_TRIANG_TETRA( MG_, lvl, it) {
        LocalP1CL<> lp1_cn( *it, conc, Bnd_);
        LocalP2CL<> lp2_cn( lp1_cn );
        absdet= std::abs( it->GetVolume()*6.);
        patch.Init( *it, lset_.Phi,0.);
        for (int ch=0; ch<8; ++ch)
        {
            // compute volume and concentration
            patch.ComputeCutForChild(ch);
            //Volume += patch.quad( ones, absdet, false);
            //c_avrg += std::abs(patch.quad( lp2_cn, absdet, false));
            con  += std::abs(patch.quad( lp2_cn, absdet, false));
        }
    }
    conc.Data *= totalconcentration/con;
    oldconc.Data *= totalconcentration/con;
	DoInterfaceStep(t);
}

///Calls 
/// -InitStep (Setup of linear osmosis diffusion system)
/// -DoStep (Solution of linear system)
/// -DoInterfaceStep (Solution of interface system)
/// -CommitStep (Updating of the vectors)
void OsmosisP1CL::DoStep (double new_t)
{
    VectorCL rhs( oldconc.Data.size());
    oldt_= t_;
    t_= new_t;

	//Velocity field
	DoInterfaceStep(new_t);
	//Levelset
	DoLevelsetStep(new_t);
    // //Diffusion
	InitStep( rhs);
	DoStep( rhs);


	static int CoupActive = P.get<int>("Coupling.Active", 0);

	if (CoupActive)
	{

		VecDescCL prevVn, prevConc;
		VectorCL prevLset;
		prevVn.SetIdx(&Velidx);
		prevConc.SetIdx(&idx);
		double res(0.0);
		int i;

		static int CoupMaxIter = P.get<int>("Coupling.maxIter", 10);
		static double CoupTol = P.get<double>("Coupling.Tol", 1e-6);
		static double alpha = P.get<double>("Coupling.Damping", 1.0);
		std::cout << "Coupling started...\n";

		for (i = 0; i < CoupMaxIter; ++i)
		{
			//this->CheckDivergence();

			prevVn.Data = Vn_.Data;
			prevConc.Data = conc.Data;
			prevLset.resize(lset_.Phi.Data.size(), 0.0);
			prevLset = lset_.Phi.Data;

			//Velocity field
			DoInterfaceStep(new_t);
			//Levelset
			DoLevelsetStep(new_t);
			//Diffusion
			InitStep( rhs);
			DoStep( rhs);

			//DAMPING
			conc.Data *= alpha;
			conc.Data += (1.0-alpha) *prevConc.Data;

			Vn_.Data *= alpha;
			Vn_.Data += (1.0-alpha) * prevVn.Data;

			lset_.Phi.Data *= alpha;
			lset_.Phi.Data += (1.0-alpha) * prevLset;

			res = norm(prevVn.Data - Vn_.Data)
				+ norm(prevConc.Data - conc.Data)
				+ norm(prevLset - lset_.Phi.Data);

			if (res < CoupTol)
				break;

			std::cout << "\t Coupling(Vn): res = " << norm(prevVn.Data - Vn_.Data) << ", iter = " << i <<"\n";
			std::cout << "\t Coupling(conc): res = " << norm(prevConc.Data - conc.Data) << ", iter = " << i <<"\n";
			std::cout << "\t Coupling(lset): res = " << norm(prevLset - lset_.Phi.Data) << ", iter = " << i <<"\n";
			std::cout << "\t Coupling: res = " << res << ", iter = " << i <<"\n";

		}
		std::cout << "Coupling: res = " << res << ", iter = " << i <<"\n";
	}

    CommitStep();

	//REPARAM
	static int param = P.get("Reparam.Active", 0);
	static int method = P.get("Reparam.method", 13);
	static bool periodic = 0;
	if (param)
		lset_.Reparam(method, periodic);

}

///Sets up the l.h.s. matrix and r.h.s. vectors
///result is the member MLMatrixCL L_ and the vector rhs
void OsmosisP1CL::InitStep (VectorCL& rhs)
{
    SetTwoStepIdx();

    //Setting up: integral over OLD Omega with OLD u (used for rhs ONLY)
    //A better implementation would directly assemble r.h.s. without a matrix assembly before
    SetupInstatMixedMassMatrix( M, b, oldt_);
    rhs.resize(oldconc.Data.size());

    rhs =  M.Data*oldconc.Data + (1-theta_)*dt_*b.Data;

    SetNewIdx();
    //M: int over NEW Omega u*v
    //A: int over NEW Omega mu*grad(u)*grad(v)
    //C: int over OLD Omega mu*grad(u)*grad(v)
    //u = u(t^n+1)
    SetupInstatSystem( A, M, C, b, t_);
    rhs = rhs + theta_*dt_*b.Data;

    MLMatrixCL L1, L2;
    //Full Matrix: M + alpha_tilde*A + (1-alpha_tilde)*C
    L1.LinComb( theta_, A.Data, (1-theta_), C.Data);
    A.Data.clear();
    C.Data.clear();

    L_.LinComb( 1., M.Data, dt_, L1);
    L1.clear();
    L2.clear();
    M.Data.clear();

}

///Solve the linear equations which were set up in OsmosisP1CL::InitStep
void OsmosisP1CL::DoStep (const VectorCL& rhs)
{
    std::cout << "Diffusion: Before solve res = " << norm( L_*conc.Data - rhs) << std::endl;
    gm_.Solve( L_, conc.Data, rhs, conc.RowIdx->GetEx());
    L_.clear();
    std::cout << "Diffusion: res = " << gm_.GetResid() << ", iter = " << gm_.GetIter()<<"\n";
}

void OsmosisP1CL::DoInterfaceStep(double& t)
{
	//Data structure
	MVn.Data.clear();
	MVn.SetIdx(&Velidx, &Velidx);
	bVn.SetIdx(&Velidx);
	Vn_.SetIdx(&Velidx);

	//Velocity field
	SetupInterfaceVelSystem(MVn, bVn, Velidx, t);

	//1 = GM, 2 = CG
	static int intSolver = P.get<int>("Solver.Interface");

	std::cout << "Velocity: Before solve res = " << norm( MVn.Data*Vn_.Data - bVn.Data) << std::endl;

	switch (intSolver)
	{
	case 1: default:
		gm_.Solve( MVn.Data, Vn_.Data, bVn.Data, Vn_.RowIdx->GetEx());
		std::cout << "(GMRes)Velocity: res = " << gm_.GetResid() << ", iter = " << gm_.GetIter()<<"\n";
		break;
	case 2:
		int iter;
		double res;
		cg_.Solve( MVn.Data, Vn_.Data, bVn.Data, Vn_.RowIdx->GetEx(), iter, res);
		std::cout << "(CG)Velocity: res = " << res << ", iter = " << iter<<"\n";
		break;
	}

}

void OsmosisP1CL::DoLevelsetStep(double& )
{
	//Solve Levelset problem
	const_DiscVelP1SolCL Vn_eval(&Vn_, &Bnd_v_, &MG_);
	lset_.SetupSystem(Vn_eval, dt_);

    const double theta_lset_ = 0.5;
	LsetL_.LinComb( 1./dt_, lset_.E, theta_lset_, lset_.H);
	VectorCL lset_rhs_((1./dt_)*(lset_.E*oldlset_.Phi.Data) - (1.-theta_lset_)*(lset_.H*oldlset_.Phi.Data));

	std::cout << "Levelset: Before solve res = " << norm( LsetL_*lset_.Phi.Data - lset_rhs_) << std::endl;
	gm_.Solve( LsetL_, lset_.Phi.Data, lset_rhs_, lset_.Phi.RowIdx->GetEx());
	std::cout << "Levelset: res = " << gm_.GetResid() << ", iter = " << gm_.GetIter()<<"\n";
	LsetL_.clear();

}


void OsmosisP1CL::CheckDivergence()
{
	double maxdiv = 0.0, detTet;
	Point3DCL Grad[4];

	DROPS_FOR_TRIANG_TETRA(MG_, MG_.GetLastLevel(), sit)
	{
		P1DiscCL::GetGradients( Grad, detTet, *sit);
		LocalP1CL<Point3DCL> locVn(*sit, Vn_, Bnd_v_);

		double tmpdiv = 0.0;
		for (int i=0; i < 3; ++i)
		{
			for (int j=0; j < 4; ++j)
				tmpdiv += locVn[j][i]*Grad[j][i];
		}

		if (tmpdiv > maxdiv)
			maxdiv = tmpdiv;
	}

	std::cout << "MAXDIV: " << maxdiv << "\t MAXDIV^-1: " << 1./maxdiv << std::endl;
	if (1./this->dt_ < maxdiv)
	{
		std::cout << "---------------------------------\n"
					 "WARNING: Divergence check failed!\n"
					 "---------------------------------" << std::endl;
	}
}

void OsmosisP1CL::SetupInterfaceVelSystem(MatDescCL& matM, VecDescCL& b, IdxDescCL& idx, double )
{
	VecDescCL cLset = lset_.Phi;
	InterfaceTriangleCL tri;
	double det, detTet, hstab;
	Quad5_2DCL<double> q[4], m, qk, qu;
	Quad5_2DCL<Point3DCL> PGrad, n, e[3];
	Point3DCL eInit;
	LocalNumbP1CL locIdx;
	Point3DCL Grad[4];
	bool partOfInterface;
	static double alpha = P.get<double>("Osmosis.alpha");
	static double beta = P.get<double>("Osmosis.beta");

	LocalP1CL<> p1[4];
	p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; // P1-Basis-Functions

	eInit[0] = 1.;
	e[0] = eInit;
	eInit[0] = 0.; eInit[1] = 1.;
	e[1] = eInit;
	eInit[1] = 0.; eInit[2] = 1.;
	e[2] = eInit;


	static bool stabon = P.get("Stabilisation.Active", 0);
	static double stabexp = P.get("Stabilisation.Exp", 0.0);
	static double stabsca = P.get("Stabilisation.Scaling", 1.0);

	matM.Data.clear();
	const IdxT num_unks = idx.NumUnknowns();

	MatrixBuilderCL M(&matM.Data, num_unks,  num_unks);

	DROPS_FOR_TRIANG_TETRA(MG_, MG_.GetLastLevel(), sit)
	{

		P1DiscCL::GetGradients( Grad, detTet, *sit);
		detTet = std::abs(detTet);
		// Scaling for stabilisation is done on interface term
		// because we want to stabilise with the smallest h
		// which occurs at the interface if adap is used
		if (stabon)
		{
			const double h = cbrt(detTet);
			hstab = 1./(stabsca*pow(h, stabexp));
		}
		else
			hstab = 1.;

		locIdx.assign( *sit, *matM.RowIdx, Bnd_v_);

		//compute triangles that form the approximation of the interface
		tri.Init(*sit, cLset, 0.);
		partOfInterface = tri.Intersects();

		//Insert stabilisation parameters
		if (stabon)
			SetupInterfaceStabilisation(*sit, detTet, M, locIdx, Grad);

		if (partOfInterface)
		{
			//loop over the 8 sub-tetras (childs) that were constructed
			for (int ch= 0; ch < 8; ++ch)
			{
				//continue if there exists no intersection with the interface on the child
				if (!tri.ComputeForChild( ch))
					continue;
				//switched to child now

				det = tri.GetAbsDet();
				//loop over triangles that form the patch (maximum of 2)
				for (int t= 0; t < tri.GetNumTriangles(); ++t)
				{

					//u v = (I - n n^T) Grad v
					n = tri.GetImprovedNormal(t);

					if (t == 1) //Quadrilateral
						det *= tri.GetAreaFrac();
					if (t > 1)
						std::cout << "Quadrilateral divided in > 2 triangles, should not happen\n";

					//first we need quadrature objects for the ansatz functions
					for (int i=0; i < 4; ++i)
						q[i].assign(p1[i], &tri.GetBary(t));
					//then we can integrate
					double intgrl;

					//mass matrix
					for (int i= 0; i < 4; ++i)
					{

						LocalP1CL<> p1conc(*sit, conc, Bnd_);
						Quad5_2DCL<> qconc(p1conc, &tri.GetBary(t));
						LocalP1CL<Point3DCL> p1loc(Grad[i]);
						Quad5_2DCL<Point3DCL> p1grad(p1loc, &tri.GetBary(t));

						for (int k = 0; k < 3; k++)
						{
							// alpha * Grad_Gamma id * Grad_Gamma v
							PGrad = (e[k] - dot(e[k], n)*n);
							qk = dot(PGrad, p1grad);

							//beta * u * v * n
							qu = qconc*q[i]*dot(e[k], n);

							b.Data[locIdx.num[i] + k] += beta*hstab*qu.quad(det) - alpha*hstab*qk.quad(det) ;
						}

						m = q[i]*q[i];
						intgrl = m.quad(det);
						for (int k = 0; k < 3; k++)
							M(locIdx.num[i] + k, locIdx.num[i] + k) += hstab*intgrl;

						for (int j = 0; j < i; ++j)
						{
							m = q[j]*q[i];
							intgrl = m.quad(det);
							for (int k = 0; k < 3; k++)
							{
								M(locIdx.num[i] + k, locIdx.num[j] + k) += hstab*intgrl;
								M(locIdx.num[j] + k, locIdx.num[i] + k) += hstab*intgrl;
							}
						}//elem mass done
					}//mass matrix done
				}//triangles of child done
			}//childs done
		}
	}//tetras done
	M.Build();
}

void OsmosisP1CL::SetupInterfaceStabilisation(TetraCL& sit, double detTet, MatrixBuilderCL& matM, LocalNumbP1CL& locIdx, Point3DCL Grad[4])
{
	double gtgint;
	static double stabaniso = P.get("Stabilisation.Anisotropic", 0);
	LocalP1CL<> p1[4];
	p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; // P1-Basis-Functions
	static P1FEGridfunctions p1feq;
	TransformedP1FiniteElement & transfp1fel = *new TransformedP1FiniteElement (p1feq);
	transfp1fel.SetTetra(sit);
	const SMatrixCL<4,4> GTG = transfp1fel.GetGramShape();

	Quad3CL<Point3DCL> nTet = GetImprovedTetraNormal(sit, lset_.Phi, Bnd_ls_);

    LocalP2CL<> lsetp2(sit, lset_.Phi, Bnd_ls_);
    
    Quad3CL<double> lsetq(lsetp2);
    Quad3CL<double> weight (1.0);

	static bool weighteddiff = P.get<int>("ZeroVelAtBnd.Active",0) 
        && P.get<int>("ZeroVelAtBnd.WeightedDiffusion",1);

    
    if (weighteddiff)
    {
        static double weighteddiffdist = P.get<double>("ZeroVelAtBnd.WeightedDiffDistance",0.25);
        static double weighteddifffact = P.get<double>("ZeroVelAtBnd.WeightedDiffFactor",0.001);

        for (Uint i = 0; i < Quad3DataCL::NumNodesC; ++i)
            if (lsetq[i] > weighteddiffdist)
                weight[i] = weighteddifffact;
    }

	//if (!partOfInterface)
	for(int i= 0; i < 4; ++i)
	{
		Quad3CL<> m3d(p1[i]);


		Quad3CL<Point3DCL> GradStar[4];
		Quad3CL<Point3DCL> Gradj;
		Quad3CL<double> GradiGradj;
		if (stabaniso)
		{
			//calculate Grad^Star phi_i = n*n^T*Grad phi_i
			for (int j=0; j < 4; ++j)
			{
				Gradj = Grad[j];
				GradStar[j] = nTet*dot(Gradj, nTet);
			}
			GradiGradj = dot(GradStar[i], GradStar[i]) * weight;
			gtgint = GradiGradj.quad(detTet);
		}
		else
			gtgint = GTG( i, i)/6.0*detTet;

		for(int k = 0; k < 3; ++k)
		{
            if (locIdx.num[i] != NoIdx)
                matM(locIdx.num[i] + k, locIdx.num[i] + k)	+= gtgint;
		}
		for(int j= 0; j < i; ++j)
		{
			if (stabaniso)
			{
				GradiGradj = dot(GradStar[i], GradStar[j]) * weight;
				gtgint = GradiGradj.quad(detTet);
			}
			else
				gtgint = GTG( i, j)/6.0*detTet;

			for(int k = 0; k < 3; ++k)
			{
                if (locIdx.num[i] != NoIdx && locIdx.num[j] != NoIdx)
                {    
                    matM(locIdx.num[i] + k, locIdx.num[j] + k)	+= gtgint;
                    matM(locIdx.num[j] + k, locIdx.num[i] + k)	+= gtgint;
                }
			}
		}
	}
	delete &transfp1fel;
}

/// Sets matrices to be of size n_new x n_old, s.t. bilinearform-applications A(uold,v)
/// make sense - This is just used for r.h.s. vectors
void OsmosisP1CL::SetTwoStepIdx()
{
    conc.SetIdx( &idx);
    b.SetIdx( &idx);

    M.Data.clear();
    M.SetIdx( &idx, &oldidx);
    A.Data.clear();
    A.SetIdx(  &idx, &oldidx);
    C.Data.clear();
    C.SetIdx(  &idx, &oldidx);
}

/// Sets matrices to be of size n_new x n_new
void OsmosisP1CL::SetNewIdx()
{
    conc.SetIdx( &idx);
    b.SetIdx( &idx);

    M.Data.clear();
    M.SetIdx( &idx, &idx);
    A.Data.clear();
    A.SetIdx( &idx, &idx);
    C.Data.clear();
    C.SetIdx( &idx, &idx);
}

///change numberings and vectors
void OsmosisP1CL::CommitStep ()
{
    oldlset_.Phi.Data= lset_.Phi.Data;
    oldconc.SetIdx(&idx);
    oldconc.Data= conc.Data;
    oldVn_.SetIdx(&oldVelidx);
    oldVn_.Data= Vn_.Data;
}


/// Setup of all volume integral - Bi- and Linearforms
void OsmosisP1CL::SetupInstatSystem(MatrixCL& matA,
    MatrixCL& matM, MatrixCL& matC, VecDescCL *b,
    IdxDescCL& RowIdx, const double time) const
{

    if (b!=0) b->Data= 0.;
    matM.clear();
    matA.clear();
    matC.clear();
    const IdxT num_unks=  RowIdx.NumUnknowns();
    MatrixBuilderCL A(&matA, num_unks,  num_unks),// diffusion wrt NEW OMEGA
                    M(&matM, num_unks,  num_unks),// mass matrix
                    C(&matC, num_unks,  num_unks);// diffusion wrt OLD OMEGA

    const Uint lvl= RowIdx.TriangLevel();
    LocalNumbP1CL n;

    P1FEGridfunctions p1feq;

    TransformedP1FiniteElement & transfp1fel = *new TransformedP1FiniteElement (p1feq);

    GlobalConvDiffReacCoefficients global_cdcoef(D_, GetVelocity() , c_, f_ ,time);
	
    ConvDiffElementMatrices elmats;
    Elvec4 f; 

    DROPS_FOR_TRIANG_TETRA( MG_, lvl, sit) {
        std::memset( f, 0, 4*sizeof(double));
        transfp1fel.SetTetra(*sit);

        n.assign( *sit, RowIdx, Bnd_);
        InterfaceTetraCL cut;
        InterfaceTetraCL oldcut;
        cut.Init( *sit, lset_.Phi,0.);
        oldcut.Init( *sit, oldlset_.Phi, 0.);
        const bool nocut=!cut.Intersects();
        const bool oldnocut=!oldcut.Intersects();
                

        LocalConvDiffReacCoefficients local_cdcoef(global_cdcoef,*sit);
		bool nPartNew = false;
		for (int i=0; i < 10; i++)
			if (cut.GetSign( i) != 1)
			{
				nPartNew = true;
				break;
			}
		bool nPartOld = false;
		for (int i=0; i < 10; i++)
			if (oldcut.GetSign( i) != 1)
			{
				nPartOld = true;
				break;
			}

        //transfp1fel.SetTetra(*sit);
        transfp1fel.SetLocal(*sit,local_cdcoef,!nPartNew);
        elmats.ResetAll();

        if (nocut) // tetra is NOT intersected by the NEW interface
        {
        	if (oldnocut){ //and NOT intersected by the OLD interface
                if (nPartOld)
        			SetupLocalOnePhaseSystem (transfp1fel, elmats, local_cdcoef, ONLYOLD);
        		if (nPartNew)
        		{
        			SetupLocalOnePhaseSystem (transfp1fel, elmats, local_cdcoef, BOTH);
        			ComputeRhsElementVector( f, local_cdcoef, transfp1fel);
        		}
        	}
        	else { //tetra is NOT intersected by the NEW interface but IS intersected by OLD.
        		if (nPartNew) // if tetra IS inside interface wrt NEW Omega, then assemble M and A
        		{
        			SetupLocalOnePhaseSystem (transfp1fel, elmats, local_cdcoef, ONLYNEW);
        			ComputeRhsElementVector( f, local_cdcoef, transfp1fel);
        		}
        		// and for old interface: compute C only
        		SetupLocalOneInterfaceSystem(transfp1fel, oldcut, elmats, local_cdcoef, ONLYOLD);
        	}
        }
        else // tetra IS intersected by the NEW interface
        {
        	//M and A wrt NEW cut, no C computation!
            SetupLocalOneInterfaceSystem(transfp1fel, cut, elmats, local_cdcoef, ONLYNEW);
            SetupLocalTwoPhaseRhs(transfp1fel, cut, f, local_cdcoef);
        	if (oldnocut) { //NOT intersected by old interface
        		if (nPartOld) // tetra was inside interface at OLD time -> compute full C
        			SetupLocalOnePhaseSystem (transfp1fel, elmats, local_cdcoef, ONLYOLD);
        	}
        	else { //both interfaces intersect tetra.
        		//compute C only with OLD interface
        		SetupLocalOneInterfaceSystem(transfp1fel, oldcut, elmats, local_cdcoef, ONLYOLD);
        	}

        }


        // assemble couplings between standard basis functions
        for(int i= 0; i < 4; ++i)
            if (n.WithUnknowns( i)){
                for(int j= 0; j < 4; ++j)
                    if (n.WithUnknowns( j)) {
                        M( n.num[i], n.num[j])+= elmats.M[i][j];
                        A( n.num[i], n.num[j])+= elmats.A[i][j];
                        C( n.num[i], n.num[j])+= elmats.C[i][j];
                    }
                if (b!=0) b->Data[n.num[i]]+= f[i];
            }

    }//end of drops triang loop


    A.Build();
    M.Build();
    C.Build();
    delete &transfp1fel;
}

/// Setup of all volume integral - Bi- and Linearforms (not Nitsche yet, this is in SetupNitscheSystem)
/// - For multilevel
void OsmosisP1CL::SetupInstatSystem (MLMatDescCL& matA,
    MLMatDescCL& matM, MLMatDescCL& matC, VecDescCL& b, const double time) const
{
    MLMatrixCL::iterator itA = --matA.Data.end();
    MLMatrixCL::iterator itM = --matM.Data.end();
    MLMatrixCL::iterator itC = --matC.Data.end();
    MLIdxDescCL::iterator it = --matA.RowIdx->end();
    SetupInstatSystem (*itA, *itM, *itC, &b, *it, time);
    for (size_t num= 1; num < matA.Data.size(); ++num, --itA, --itM, --itC, --it)
        SetupInstatSystem (*itA, *itM, *itC, &b, *it, time);
}


/// Couplings between basis functions wrt old interfaces, s.t. Bilinearform-Applications
/// M(uold,v) make sense also for the new time step (and the functions therein (like v))
// This is only used as a matrix application. So actually there is no need to setting up the matrix!
void OsmosisP1CL::SetupInstatMixedMassMatrix( MatrixCL& matM, VecDescCL* b, IdxDescCL& RowIdx, IdxDescCL& ColIdx,
    const double time) const
{
    if (b!= 0) b->Data = 0.;
    matM.clear();
    //const ExtIdxDescCL& Xidx= RowIdx.GetXidx();
    //const ExtIdxDescCL& oldXidx= ColIdx.GetXidx();
    const IdxT num_unks=  RowIdx.NumUnknowns(),
               num_cols=  ColIdx.NumUnknowns();
    MatrixBuilderCL M(&matM, num_unks,  num_cols);//mass matrix
    const MultiGridCL& mg= this->GetMG();
    BndDataCL<> Bndlset(mg.GetBnd().GetNumBndSeg());                
    const Uint lvl= RowIdx.TriangLevel();
    LocalNumbP1CL n;
    bool nocut; // oldsign[4],

    P1FEGridfunctions p1feq;
    TransformedP1FiniteElement & transfp1fel = *new TransformedP1FiniteElement (p1feq);    

    GlobalConvDiffReacCoefficients global_cdcoef(D_, GetVelocity() , c_, f_, time);

    Elvec4 f; 

    DROPS_FOR_TRIANG_TETRA( MG_, lvl, sit) { 
        std::memset( f, 0, 4*sizeof(double));
        transfp1fel.SetTetra(*sit);

        n.assign( *sit, RowIdx, Bnd_);
        InterfaceTetraCL oldcut;
        oldcut.Init( *sit, oldlset_.Phi,0.);
        nocut=!oldcut.Intersects();

        LocalConvDiffReacCoefficients local_cdcoef(global_cdcoef,*sit);
        // HERE: OLD INTERFACE!!
        bool nPart = false;
        for (int i=0; i < 10; i++)
        	if (oldcut.GetSign( i) != 1)
        		nPart = true;

        transfp1fel.SetLocal(*sit,local_cdcoef, !nPart);

        Elmat4x4 M_P1NEW_P1OLD_n;
        std::memset( M_P1NEW_P1OLD_n, 0, 4*4*sizeof(double));

        if(nocut)
        { //old interface does not cut
        	if (nPart)
        	{
        		SetupLocalOnePhaseMassMatrix ( M_P1NEW_P1OLD_n, transfp1fel);
        		ComputeRhsElementVector( f, local_cdcoef, transfp1fel);
        	}

		}
		else
		{ //old interface does cut
			SetupLocalOneInterfaceMassMatrix( oldcut, M_P1NEW_P1OLD_n, transfp1fel);
			SetupLocalTwoPhaseRhs(transfp1fel, oldcut, f, local_cdcoef);
		}
        
        
        for(int i= 0; i < 4; ++i)
            if (n.WithUnknowns( i)){
                for(int j= 0; j < 4; ++j)
                {
                    if (n.WithUnknowns( j)) {
                        M( n.num[i], n.num[j])+= M_P1NEW_P1OLD_n[i][j];
                    }
                }
                if (b!= 0) b->Data[n.num[i]]+= f[i];
            }
            else
                throw DROPSErrCL("this is not happening!");
    }
    M.Build();
    delete &transfp1fel;
}

void OsmosisP1CL::SetupInstatMixedMassMatrix(MLMatDescCL& matM, VecDescCL& b, const double time) const
{
    MLMatrixCL::iterator itM = --matM.Data.end();
    MLIdxDescCL::iterator it_row = --matM.RowIdx->end();
    MLIdxDescCL::iterator it_col = --matM.ColIdx->end();
    SetupInstatMixedMassMatrix(*itM, &b, *it_row, *it_col, time);
    for (size_t num= 1; num < matM.Data.size(); ++num, --itM, --it_row, --it_col)
        SetupInstatMixedMassMatrix (*itM, &b, *it_row, *it_col, time);
}


//*****************************************************************************
//                               OsmosisRepairCL
//						FUNKTIONIERT AKTUELL NICHT
//					   NEUE RepairAfterRefine NÖTIG
//*****************************************************************************
void
OsmosisRepairCL::post_refine ()
{
	throw DROPSErrCL( "OsmosisRepairCL::post_refine: Sorry, not yet implemented.");
    VecDescCL loc_conc;
    IdxDescCL loc_cidx( P1X_FE, c_.GetBndData(), c_.GetXFEMOmitBound());
    VecDescCL loc_oldconc;
    IdxDescCL loc_oldcidx( P1X_FE, c_.GetBndData(), c_.GetXFEMOmitBound());
    VecDescCL& conc= c_.conc;
    VecDescCL& oldconc= c_.oldconc;

    loc_cidx.CreateNumbering( mylvl, c_.GetMG(), c_.GetBndData(), &(c_.GetLevelset()),&(c_.GetLevelsetBnd()));
    loc_oldcidx.CreateNumbering( mylvl, c_.GetMG(), c_.GetBndData(), &(c_.GetOldLevelset()),&(c_.GetLevelsetBnd()));
    loc_conc.SetIdx( &loc_cidx);
    loc_oldconc.SetIdx( &loc_oldcidx);
    //RepairAfterRefineP1( c_.GetSolution( conc), loc_conc);
    //RepairAfterRefineP1( c_.GetSolution( oldconc), loc_oldconc);
    double t = conc.t;
    conc.Clear(t);
    c_.DeleteNumbering( &c_.idx);
    c_.idx.GetFinest().swap( loc_cidx);
    conc.SetIdx( &c_.idx);
    conc.Data= loc_conc.Data;

    double oldt = oldconc.t;
    oldconc.Clear(oldt);
    c_.DeleteNumbering( &c_.oldidx);
    c_.oldidx.GetFinest().swap( loc_oldcidx);
    oldconc.SetIdx( &c_.oldidx);
    oldconc.Data= loc_oldconc.Data;
}

void
  OsmosisRepairCL::pre_refine_sequence ()
{
    oldp1xrepair_= std::unique_ptr<P1XRepairCL>( new P1XRepairCL( c_.GetMG(), c_.oldconc));
}

void
  OsmosisRepairCL::post_refine_sequence ()
{
     c_.CreateNumbering(mylvl, &c_.idx, &c_.oldidx, &c_.Velidx, &c_.oldVelidx, (c_.GetLevelset()), (c_.GetOldLevelset()));
    (*oldp1xrepair_)();
     oldp1xrepair_.reset();
     c_.conc.SetIdx( &c_.idx);
}

const IdxDescCL* OsmosisRepairCL::GetIdxDesc() const {
  throw DROPSErrCL( "OsmosisRepairCL::GetIdxDesc: Sorry, not yet implemented.");
  return 0;
}


const IdxDescCL* VelOsmosisRepairCL::GetIdxDesc() const {
  throw DROPSErrCL( "VelOsmosisRepairCL::GetIdxDesc: Sorry, not yet implemented.");
  return 0;
}

void
  VelOsmosisRepairCL::post_refine ()
{
	throw DROPSErrCL( "VVelOsmosisRepairCL::post_refine: Sorry, not yet implemented.");
    VelVecDescCL loc_v;
    IdxDescCL    loc_vidx(vecP2_FE);
    VelVecDescCL& v= v_;
    Uint LastLevel= mg_.GetLastLevel();
    loc_vidx.CreateNumbering( mg_.GetLastLevel(), mg_, Bnd_v_);
    if (LastLevel != v.RowIdx->TriangLevel()) {
        std::cout << "LastLevel: " << LastLevel
                  << " old v->TriangLevel: " << v.RowIdx->TriangLevel() << std::endl;
        throw DROPSErrCL( "VelOsmosisRepairCL::post_refine: Sorry, not yet implemented.");
    }
    loc_v.SetIdx( &loc_vidx);

#ifdef _PAR
    GetPMG().HandleNewIdx(&v.vel_idx, &loc_v);
#endif
    //RepairAfterRefineP2( const_DiscVelSolCL( &v, &Bnd_v_, &mg_), loc_v);
#ifdef _PAR
    GetPMG().CompleteRepair( &loc_v);
#endif

    double t = v.t;
    v.Clear(t);
    (*v.RowIdx).DeleteNumbering( mg_);

    vidx_.swap( loc_vidx);
    v.SetIdx( &vidx_);
    v.Data= loc_v.Data;
}

} // end of namespace DROPS
