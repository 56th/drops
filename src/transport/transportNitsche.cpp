/// \file transportNitsche.cpp
/// \brief Classes that constitute a 2-phase-transport-problem with Nitsche-XFEM discretization.
/// \author Trung Hieu Nguyen (small fixes: Martin Horsky, Christoph Lehrenfeld, Joerg Grande) IGPM

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
#include "transport/transportNitsche.h"
#include "transport/localsetups.cpp"
#include <iostream>
#include <fstream>
namespace DROPS
{
//====================================================
//
//TranportP1XCL
//
//====================================================

/// Initialize P1X function
/// (different initial values inside and outside second phase
/// are provided by cn and cp)
void TransportP1XCL::Init (instat_scalar_fun_ptr cn, instat_scalar_fun_ptr cp, double t)
{
    ct.t = t;
    oldct.t = t;
    const Uint lvl= ct.GetLevel(),
               ctidx= ct.RowIdx->GetIdx(),
               oldctidx= oldct.RowIdx->GetIdx();
    const IdxDescCL& idx1 = idx.GetFinest();
    const IdxDescCL& idx2 = oldidx.GetFinest();
    const ExtIdxDescCL& Xidx= idx1.GetXidx();
    const ExtIdxDescCL& oldXidx= idx2.GetXidx();

    DROPS_FOR_TRIANG_VERTEX( MG_, lvl, it) {
        if (it->Unknowns.Exist( ctidx)){
            bool nPart = lset_.Data[it->Unknowns( lset_.RowIdx->GetIdx())] <= 0.;
            if (nPart)
                ct.Data[it->Unknowns( ctidx)]= H_*cn( it->GetCoord(), t);
            else
                ct.Data[it->Unknowns( ctidx)]= cp( it->GetCoord(), t);
            if (Xidx[it->Unknowns(ctidx)]==NoIdx) continue; //no xfem-enrichment function on this vertex

            // xfem coefficients are set s.t. discontinuity is realized across
            // the interface
            ct.Data[Xidx[it->Unknowns( ctidx)]]=cp( it->GetCoord(), t)- H_*cn( it->GetCoord(), t);
        }
    }
    // loop is called a second time, as  oldctidx and ctidx may differ
    DROPS_FOR_TRIANG_VERTEX( MG_, lvl, it) {
        if (it->Unknowns.Exist( oldctidx)){
            bool nPart= oldlset_.Data[it->Unknowns( oldlset_.RowIdx->GetIdx())] <= 0.;
            if (nPart)
                oldct.Data[it->Unknowns( oldctidx)]= H_*cn( it->GetCoord(), t);
            else
                oldct.Data[it->Unknowns( oldctidx)]= cp( it->GetCoord(), t);
            if (oldXidx[it->Unknowns(oldctidx)]==NoIdx) continue; //no xfem-enrichment function on this vertex
            // xfem coefficients are set s.t. discontinuity is realized across
            // the interface
            oldct.Data[oldXidx[it->Unknowns( oldctidx)]]=cp( it->GetCoord(), t)- H_*cn( it->GetCoord(), t);
        }
    }

    //TransformWithScaling(ct, c, 1.0/GetHenry(true), 1.0/GetHenry(false));
    // \todo: as soon as c_in and c_out are members of masstransp they should be initialized as well (?)

}

/// Transform from a concentration (P1X) to another scaled
/// (different scalings inside and outside second phase
/// can be provided by scalingp and scalingn) concentration (also P1X)
void TransportP1XCL::TransformWithScaling (const VecDescCL& concin, VecDescCL& concout, double scalingp, double scalingn)
{
	std::cerr << "TransportP1XCL::TransformWithScaling is not doing the correct thing!" << std::endl;
	getchar();
    concout.t = concin.t;
    const Uint lvl= concin.GetLevel(),
               ctidx= concin.RowIdx->GetIdx(); //concin and concout have same indices
    const IdxDescCL& idx1 = idx.GetFinest();
    const ExtIdxDescCL& Xidx= idx1.GetXidx();
    double fac = 0.;
    double ofac = 0.;
    DROPS_FOR_TRIANG_VERTEX( MG_, lvl, it) {
        if (it->Unknowns.Exist( ctidx)){
           bool nPart = lset_.Data[it->Unknowns( lset_.RowIdx->GetIdx())] <= 0.;
            if (nPart){
              fac = scalingn;
              ofac = scalingp;
            }
            else{
              fac = scalingp;
              ofac = scalingn;
            }

            concout.Data[it->Unknowns( ctidx)]= fac * concin.Data[it->Unknowns( ctidx)];
            if (Xidx[it->Unknowns(ctidx)]==NoIdx)
              continue; //no xfem-enrichment function on this vertex
            else
              concout.Data[Xidx[it->Unknowns( ctidx)]] = ofac * concin.Data[Xidx[it->Unknowns( ctidx)]];
        }
    }
}


/// Compute Mean Drop Concentration, i.e. mean
/// concentration in second phase (negative sign):
/// integral over concentration / volume of second phase
double TransportP1XCL::MeanDropConcentration()
{
    VecDescCL cn (&idx);
    GetSolutionOnPart(cn, false, false);
    double absdet;
    InterfaceTetraCL patch;

    double c_avrg= 0., Volume= 0.;
    LocalP2CL<double> ones( 1.);
    const Uint lvl= ct.GetLevel();
    DROPS_FOR_TRIANG_TETRA( MG_, lvl, it) {
        LocalP1CL<> lp1_cn( *it, cn, Bnd_);
        LocalP2CL<> lp2_cn( lp1_cn );
        absdet= std::abs( it->GetVolume()*6.);
        patch.Init( *it, lset_,0.);
        for (int ch=0; ch<8; ++ch)
        {
            // compute volume and concentration
            patch.ComputeCutForChild(ch);
            Volume+= patch.quad( ones, absdet, false);
            c_avrg+= patch.quad( lp2_cn, absdet, false);
        }
    }
    c_avrg/= Volume;
    return c_avrg;
}


double TransportP1XCL::CheckSolution(instat_scalar_fun_ptr Lsgn, instat_scalar_fun_ptr Lsgp,
        instat_vector_fun_ptr Gradn, instat_vector_fun_ptr Gradp, double time)
{
    VecDescCL cn (&idx);
    GetSolutionOnPart(cn, false, false);
    VecDescCL cp (&idx);
    GetSolutionOnPart(cp, true, false);

    double errl2p = 0, errl2n = 0;
    double errl1p = 0, errl1n = 0;
    double errh1p = 0, errh1n = 0;
    InterfaceTetraCL patch;

    LocalP2CL<double> ones( 1.);
    Point3DCL G[4];  // gradients of P1 basis functions

    std::cout << "Difference to exact solution:" << std::endl;

    const Uint lvl= ct.GetLevel();

    typedef Quad5CL<> QuadT;
    typedef Quad5CL<Point3DCL> QuadVecT;
    DROPS_FOR_TRIANG_TETRA( MG_, lvl, it)
    {
        LocalP1CL<double> lp1_p (*it, cp, Bnd_);
        LocalP1CL<double> lp1_n (*it, cn, Bnd_);
        LocalP1CL<double> lp1_psol (*it, Lsgp);
        LocalP1CL<double> lp1_nsol (*it, Lsgn);
        SMatrixCL<3,3> M;
        double det;
        GetTrafoTr(M,det,*it);
        double absdet= std::fabs(det);
        P1DiscCL::GetGradients( G, M);
        patch.Init( *it, lset_,0.);
        if (patch.Intersects()){
          patch.ComputeSubTets();
          Uint NumTets=patch.GetNumTetra(); //# of subtetras
          for (Uint k=0; k<NumTets; ++k)
          {
            bool pPart = (k>=patch.GetNumNegTetra());
            const SArrayCL<BaryCoordCL,4>& T =patch.GetTetra(k);
            BaryCoordCL* nodes = QuadT::TransformNodes(T);
            QuadT q_sol= QuadT( *it, pPart ? Lsgp : Lsgn, time, nodes);
            QuadT q_dsol= QuadT(pPart? lp1_p : lp1_n, nodes);
            QuadT q_diff(q_dsol - q_sol);
            QuadT q_diff2(q_diff * q_diff);
            QuadT q_diffabs( std::abs(q_diff));
            const LocalP1CL<>& P1( pPart? lp1_p : lp1_n);
            Point3DCL grad;
            for (int i=0; i<4; ++i)
                grad+= P1[i]*G[i];
            QuadVecT q_gradSol( *it, pPart ? Gradp : Gradn, time, nodes),
                    q_grad( grad),
                    q_diffGrad( q_grad - q_gradSol);
            QuadT q_diffGrad2( dot( q_diffGrad, q_diffGrad));

            double Vol = absdet*VolFrac(T);
            if (pPart) {
              errl2p += q_diff2.quad(Vol);
              errl1p += q_diffabs.quad(Vol);
              errh1p += q_diffGrad2.quad(Vol);
            } else {
              errl2n += q_diff2.quad(Vol);
              errl1n += q_diffabs.quad(Vol);
              errh1n += q_diffGrad2.quad(Vol);
            }
            delete[] nodes;
          }

        }
        else
        {
          bool pPart= (patch.GetSign( 0) == 1);
          QuadT q_sol= QuadT( *it, pPart ? Lsgp : Lsgn, time);
          QuadT q_dsol= QuadT(pPart? lp1_p : lp1_n);
          QuadT q_diff(q_dsol - q_sol);
          QuadT q_diff2(q_diff * q_diff);
          QuadT q_diffabs( std::abs(q_diff));
          const LocalP1CL<>& P1( pPart? lp1_p : lp1_n);
          Point3DCL grad;
          for (int i=0; i<4; ++i)
              grad+= P1[i]*G[i];
          QuadVecT q_gradSol( *it, pPart ? Gradp : Gradn, time),
                  q_grad( grad),
                  q_diffGrad( q_grad - q_gradSol);
          QuadT q_diffGrad2( dot( q_diffGrad, q_diffGrad));

          if (pPart) {
            errl2p += q_diff2.quad(absdet);
            errl1p += q_diffabs.quad(absdet);
            errh1p += q_diffGrad2.quad(absdet);
          } else {
            errl2n += q_diff2.quad(absdet);
            errl1n += q_diffabs.quad(absdet);
            errh1n += q_diffGrad2.quad(absdet);
          }
        }
    }
    const double errl2 = std::sqrt(errl2n + errl2p);
    const double errl1 = errl1n + errl1p;
    const double errh1 = std::sqrt(errl2n + errl2p + errh1n + errh1p);
    std::cout << "errL2p = " << std::sqrt(errl2p) << "\t";
    std::cout << "errL2n = " << std::sqrt(errl2n) << "\t";
    std::cout << "errL2  = " << errl2 << std::endl;
    std::cout << "errL1p = " << errl1p << "\t";
    std::cout << "errL1n = " << errl1n << "\t";
    std::cout << "errL1  = " << errl1 << std::endl;
    std::cout << "errH1p = " << std::sqrt(errh1p) << "\t";
    std::cout << "errH1n = " << std::sqrt(errh1n) << "\t";
    std::cout << "errH1  = " << errh1 << std::endl;
    return errl2;
}


///Calls
/// -InitStep (Setup of linear system)
/// -DoStep (Solution of linear system)
/// -CommitStep (Updating of the vectors)
void TransportP1XCL::DoStep (double new_t)
{
    VectorCL rhs( ct.Data.size());
    oldt_= t_;
    t_= new_t;
    InitStep( rhs);
    DoStep( rhs);
    CommitStep();
}

/// Sets matrices to be of size n_new x n_old, s.t. bilinearform-applications A(uold,v)
/// make sense - This is just used for r.h.s. vectors
void TransportP1XCL::SetTwoStepIdx()
{
    UpdateXNumbering( lset_, false, true);
    std::cout<<"# old concentration unknowns: " << oldidx.NumUnknowns() << "# new concentration unknowns: " << idx.NumUnknowns()<< "\n";
    ct.SetIdx( &idx);
    c.SetIdx( &idx);
    cplM.SetIdx( &idx);
    cplA.SetIdx( &idx);
    cplC.SetIdx( &idx);
    b.SetIdx( &idx);

    M.Data.clear();
    M.SetIdx( &idx, &oldidx);
    A.Data.clear();
    A.SetIdx(  &idx, &oldidx);
    C.Data.clear();
    C.SetIdx(  &idx, &oldidx);
    NA.Data.clear();
    NA.SetIdx( &idx, &oldidx);
}

/// Sets matrices to be of size n_new x n_new
void TransportP1XCL::SetNewIdx()
{
    c.SetIdx( &idx);
    cplM.SetIdx( &idx);
    cplA.SetIdx( &idx);
    cplC.SetIdx( &idx);
    b.SetIdx( &idx);

    M.Data.clear();
    M.SetIdx( &idx, &idx);
    A.Data.clear();
    A.SetIdx( &idx, &idx);
    C.Data.clear();
    C.SetIdx( &idx, &idx);
    NA.Data.clear();
    NA.SetIdx( &idx, &idx);
}

///Sets up the l.h.s. matrix and r.h.s. vectors
///result is the member MLMatrixCL L_ and the vector rhs
void TransportP1XCL::InitStep (VectorCL& rhs)
{
    SetTwoStepIdx();
    SetupInstatMixedMassMatrix( M, cplM, oldt_);
    rhs.resize(cplM.Data.size());
    rhs =  M.Data*oldct.Data - cplM.Data;
    SetNewIdx();
    SetupInstatSystem( A, cplA, M, cplM, C, cplC, b, t_);
    SetupNitscheSystem( NA);
    rhs += cplM.Data + dt_*theta_*(b.Data + cplA.Data + cplC.Data);

    MLMatrixCL L1, L2;
    L1.LinComb( theta_, A.Data, theta_, C.Data);
    A.Data.clear();
    C.Data.clear();
    L2.LinComb( 1., L1, theta_, NA.Data);
    NA.Data.clear();
    L1.clear();
    L_.LinComb( 1., M.Data, dt_, L2);
    L2.clear();
    M.Data.clear();
}

///Solve the linear equations which were set up in TransportP1XCL::InitStep
void TransportP1XCL::DoStep (const VectorCL& rhs)
{
    std::cout << "Before solve: res = " << norm( L_*ct.Data - rhs) << std::endl;
    gm_.Solve( L_, ct.Data, rhs, ct.RowIdx->GetEx());
    L_.clear();
    std::cout << "res = " << gm_.GetResid() << ", iter = " << gm_.GetIter()<<"\n";// <<", norm of ct = " << norm(ct.Data)<< std::endl;
}

///change numberings and vectors
void TransportP1XCL::CommitStep ()
{
    oldlset_.Data= lset_.Data;
    UpdateXNumbering( oldlset_, false, false);
    oldct.SetIdx(&oldidx);
    oldct.Data= ct.Data;
}




/// Setup of all volume integral - Bi- and Linearforms (not Nitsche yet, this is in SetupNitscheSystem)
/**
 * - For one  level only \n
 * Profiling of SetupInstatSystem (made on 2011/12/02) on a typical example (CL): \n
 *               LocalOneInterfaceSetup : 58.0 %   \n
 * Local (on each element) Preparations : 13.0 %   \n
 *                   LocalOnePhaseSetup : 13.0 %   \n
 *                       Matrices Build : 11.5 %   \n
 *  Add ElementMatricesToGlobalMatrices :  2.5 %   \n
 *                   Local Rhs(P1 only) :  1.5 %   \n
 *                  Global Preparations :  0.5 %
 **/

void TransportP1XCL::SetupInstatSystem(MatrixCL& matA, VecDescCL *cplA,
    MatrixCL& matM, VecDescCL *cplM, MatrixCL& matC, VecDescCL *cplC, VecDescCL *b,
    IdxDescCL& RowIdx, const double time) const
{
    const ExtIdxDescCL& Xidx= RowIdx.GetXidx();
    if (b!=0) b->Data= 0.;
    if (cplM!=0){
        cplM->Data= 0.;
        cplA->Data= 0.;
        cplC->Data= 0.;
    }
    matM.clear();
    matA.clear();
    matC.clear();
    const IdxT num_unks=  RowIdx.NumUnknowns();
    MatrixBuilderCL A(&matA, num_unks,  num_unks),// diffusion
                    M(&matM, num_unks,  num_unks),// mass matrix
                    C(&matC, num_unks,  num_unks);// convection

    const Uint lvl= RowIdx.TriangLevel();
    LocalNumbP1CL n;
    bool sign[4];

    P1FEGridfunctions p1feq;

    TransformedP1FiniteElement & transfp1fel = sdstab_?
        *new StabilizedTransformedP1FiniteElement (p1feq,sdstab_) :
        *new TransformedP1FiniteElement (p1feq);

    GlobalConvDiffReacCoefficients global_cdcoef(D_,H_, GetVelocity() , c_, f_ ,time);

    ConvDiffElementMatrices elmats;
    ConvDiffElementVectors elvecs;

    DROPS_FOR_TRIANG_TETRA( MG_, lvl, sit) {

        transfp1fel.SetTetra(*sit);

        n.assign( *sit, RowIdx, Bndt_);
        InterfaceTetraCL cut;
        cut.Init( *sit, lset_,0.);
        bool nocut=!cut.Intersects();

        LocalConvDiffReacCoefficients local_cdcoef(global_cdcoef,*sit);
        bool pPart= (cut.GetSign( 0) == 1);
        //transfp1fel.SetTetra(*sit);
        transfp1fel.SetLocal(*sit,local_cdcoef,pPart);
        elmats.ResetAll();
        elvecs.ResetAll();

        ComputeRhsElementVector( elvecs.f, local_cdcoef, transfp1fel);

        if (nocut) // tetra is not intersected by the interface
        {
            // couplings between standard basis functions
            SetupLocalOnePhaseSystem (transfp1fel, elmats, local_cdcoef, pPart);
        }
        else{
            // compute element matrix for standard basis functions and XFEM basis functions
            SetupLocalOneInterfaceSystem(transfp1fel, cut, elmats, local_cdcoef);
            SetupLocalTwoPhaseRhs(transfp1fel, cut, elvecs, local_cdcoef);
            elmats.SetUnsignedAsSumOfSigned();
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
                     else if (cplM !=0) {
                        const double val= Bndt_.GetBndFun( n.bndnum[j])( sit->GetVertex( j)->GetCoord(), time);
                        cplM->Data[n.num[i]]-= elmats.M[i][j]*val;
                        cplA->Data[n.num[i]]-= elmats.A[i][j]*val;
                        cplC->Data[n.num[i]]-= elmats.C[i][j]*val;
                    }
                if (b!=0) b->Data[n.num[i]]+= elvecs.f[i];
            }
        if (nocut) continue; // no XFEM basis functions
        // assemble couplings between standard basis functions and XFEM basis functions
        for(int i= 0; i < 4; ++i)
            sign[i]= (cut.GetSign(i) == 1);
        for(int i= 0; i < 4; ++i)
            if(n.WithUnknowns(i)){
                const IdxT xidx_i= Xidx[n.num[i]];
                for(int j= 0; j < 4; ++j)
                    if(n.WithUnknowns(j)){
                        const IdxT xidx_j= Xidx[n.num[j]];
                        if (xidx_j!=NoIdx){
                            M( n.num[i], xidx_j)+= sign[j]? -elmats.M_n[i][j]: elmats.M_p[i][j];
                            A( n.num[i], xidx_j)+= sign[j]? -elmats.A_n[i][j]: elmats.A_p[i][j];
                            C( n.num[i], xidx_j)+= sign[j]? -elmats.C_n[i][j]: elmats.C_p[i][j];
                        }
                        if (xidx_i!=NoIdx){
                            M( xidx_i, n.num[j])+= sign[i]? -elmats.M_n[i][j]: elmats.M_p[i][j];
                            A( xidx_i, n.num[j])+= sign[i]? -elmats.A_n[i][j]: elmats.A_p[i][j];
                            C( xidx_i, n.num[j])+= sign[i]? -elmats.C_n[i][j]: elmats.C_p[i][j];
                        }
                        if ((xidx_i!=NoIdx) && (xidx_j!=NoIdx) && (sign[i]==sign[j])){
                            M( xidx_i, xidx_j)+= sign[j]? elmats.M_n[i][j]: elmats.M_p[i][j];
                            A( xidx_i, xidx_j)+= sign[j]? elmats.A_n[i][j]: elmats.A_p[i][j];
                            C( xidx_i, xidx_j)+= sign[j]? elmats.C_n[i][j]: elmats.C_p[i][j];
                        }
                    } else if (cplM !=0) {
                        if (xidx_i!=NoIdx){
                            const double val= Bndt_.GetBndFun( n.bndnum[j])( sit->GetVertex( j)->GetCoord(), time);
                            cplM->Data[xidx_i]-= (sign[i]? -elmats.M_n[i][j]: elmats.M_p[i][j])*val;
                            cplA->Data[xidx_i]-= (sign[i]? -elmats.A_n[i][j]: elmats.A_p[i][j])*val;
                            cplC->Data[xidx_i]-= (sign[i]? -elmats.C_n[i][j]: elmats.C_p[i][j])*val;
                        }
                    }
                if((xidx_i!=NoIdx) && (b!=0))
                    b->Data[xidx_i] +=sign[i] ?  - elvecs.f_n[i] :elvecs.f_p[i];
            }

    }

    A.Build();
    M.Build();
    C.Build();
    delete &transfp1fel;
}

/// Setup of all volume integral - Bi- and Linearforms (not Nitsche yet, this is in SetupNitscheSystem)
/// - For multilevel
void TransportP1XCL::SetupInstatSystem (MLMatDescCL& matA, VecDescCL& cplA,
    MLMatDescCL& matM, VecDescCL& cplM, MLMatDescCL& matC, VecDescCL& cplC, VecDescCL& b, const double time) const
{
    MLMatrixCL::iterator itA = --matA.Data.end();
    MLMatrixCL::iterator itM = --matM.Data.end();
    MLMatrixCL::iterator itC = --matC.Data.end();
    MLIdxDescCL::iterator it = --matA.RowIdx->end();
    SetupInstatSystem (*itA, &cplA, *itM, &cplM, *itC, &cplC, &b, *it, time);
    --itA;
    --itM;
    --itC;
    --it;
    for (size_t num= 1; num < matA.Data.size(); ++num, --itA, --itM, --itC, --it)
        SetupInstatSystem (*itA, 0, *itM, 0, *itC, 0, 0, *it, time);
}

void TransportP1XCL::GetSolutionOnPart( VecDescCL& ct_part, bool pPart, bool Is_ct)
{
    const Uint lvl= ct.RowIdx->TriangLevel(),
        idxnum= ct.RowIdx->GetIdx(),
        phiidx= lset_.RowIdx->GetIdx();
    const IdxDescCL& idx1= idx.GetFinest();
    const ExtIdxDescCL& Xidx= idx1.GetXidx();
    const MultiGridCL& mg= this->GetMG();

    ct_part.SetIdx( ct.RowIdx);
    VectorCL& cp= ct_part.Data;
    cp= ct.Data; //
    // add extended part, s.t. all information (seen from one side)
    // is available in terms of a P1 representation
    for( MultiGridCL::const_TriangVertexIteratorCL it= mg.GetTriangVertexBegin(lvl),
        end= mg.GetTriangVertexEnd(lvl); it != end; ++it)
    {
        if (!it->Unknowns.Exist( idxnum)) continue;
        const IdxT nr= it->Unknowns(idxnum);
        if (Xidx[nr]==NoIdx) continue;

        const bool sign= InterfaceTetraCL::Sign( lset_.Data[it->Unknowns(phiidx)]) == 1;
        if (pPart==sign) continue; // extended hat function ==0 on this part
        //different signs due to the definition of the xfem-enrichment functions
        if (pPart)
            cp[nr]+= ct.Data[Xidx[nr]];
        else
            cp[nr]-= ct.Data[Xidx[nr]];
    }
//    if (!Is_ct && !pPart) cp/=H_;
    if (!Is_ct && (GetHenry(pPart)!=1.0)) cp/=GetHenry(pPart);
}

///Assembles the Nitsche Bilinearform. Gathers the weighting functions and calls the Local NitscheSetup for each
///intersected tetrahedron
void TransportP1XCL::SetupNitscheSystem( MatrixCL& matA, IdxDescCL& RowIdx/*, bool new_time */) const
{
    const ExtIdxDescCL& Xidx= RowIdx.GetXidx();
    matA.clear();
    const IdxT num_unks=  RowIdx.NumUnknowns();
    MatrixBuilderCL A(&matA, num_unks,  num_unks);
    const Uint lvl= RowIdx.TriangLevel();
    LocalP1CL<Point3DCL> GradRef[10], Grad[10];
    P2DiscCL::GetGradientsOnRef( GradRef);
    LocalNumbP1CL ln;
    SMatrixCL<3,3> T;
    double det,VolP, VolN, kappa[2], h;
    int sign[4]={0};
    const MultiGridCL& mg= this->GetMG();
    BndDataCL<> Bndlset(mg.GetBnd().GetNumBndSeg());

    DROPS_FOR_TRIANG_TETRA( MG_, /*default level*/lvl, it)
    {
        InterfaceTetraCL patch;
        patch.Init( *it, lset_,0.);
        InterfaceTriangleCL triangle;
        triangle.Init( *it, lset_,0.);
        if (!patch.Intersects()) continue;
        for(int i= 0; i < 4; ++i)
            sign[i]= patch.GetSign(i);
        ln.assign( *it, RowIdx, Bndt_);
        GetTrafoTr( T, det, *it);
        P2DiscCL::GetGradients( Grad, GradRef, T); //TODO CL: kann weg!
        Point3DCL G[4];
        P1DiscCL::GetGradients( G, det, *it);
        const double h3= it->GetVolume()*6;
        h= cbrt( h3);
        kappa[0]=kappa[1]=0.;
        VolP=VolN=0.;
//         LocalP2CL<> lp2_lset(*it, lset_, Bndlset);
//         LocalP1CL<Point3DCL> lp1_grad_lset;
//         lp1_grad_lset *=0.;
//         for (int i=0; i<10; ++i)
//             lp1_grad_lset+= lp2_lset[i]*Grad[i];
//         LocalP2CL<Point3DCL> lp2_grad_lset(lp1_grad_lset);
        patch.ComputeSubTets();
        Uint NumTets=patch.GetNumTetra(); /// # of subtetras

        for (Uint k=0; k< NumTets; ++k){
            bool pPart= (k>=patch.GetNumNegTetra());
            const SArrayCL<BaryCoordCL,4>& TT =  patch.GetTetra(k);
            if (!IsRegBaryCoord(TT)) continue;
            if (pPart) VolP+= VolFrac(TT);
            else  VolN+= VolFrac(TT);
        }
        kappa[0]= VolP;
        kappa[1]= 1.-kappa[0];
        for (int ch= 0; ch < 8; ++ch)
        {
            triangle.ComputeForChild( ch);
            for (int t= 0; t < triangle.GetNumTriangles(); ++t) {
                const BaryCoordCL * p = &triangle.GetBary( t);
                //Quad5_2DCL<Point3DCL> n(lp2_grad_lset, p);
                Quad5_2DCL<Point3DCL> n(triangle.GetNormal(), p);
                //for (int i =0; i<Quad5_2DDataCL::NumNodesC; i++) if (n[i].norm()>1e-8) n[i]/= n[i].norm();
                SetupLocalNitscheSystem( p, Xidx, n, G, ln, A, triangle.GetAbsDet( t), D_, H_, kappa, lambda_, h, sign/*, new_time*/);
                }
        } // Ende der for-Schleife ueber die Kinder
    }
    A.Build();
}

void TransportP1XCL::SetupNitscheSystem (MLMatDescCL& matA) const
{
    MLMatrixCL::iterator itA = --matA.Data.end();
    MLIdxDescCL::iterator it = --matA.RowIdx->end();
    SetupNitscheSystem (*itA, *it);
    for (size_t num= 1; num < matA.Data.size(); ++num, --itA, --it)
        SetupNitscheSystem(*itA, *it);
}

/// Couplings between basis functions wrt old and new interfaces, s.t. Bilinearform-Applications
/// M(uold,v) make sense also for the new time step (and the functions therein (like v))
// This is only used as a matrix application. So actually there is no need to setting up the matrix!
void TransportP1XCL::SetupInstatMixedMassMatrix( MatrixCL& matM, VecDescCL* cplM,
                                                 IdxDescCL& RowIdx, IdxDescCL& ColIdx,
                                                 const double time) const
{
    if (cplM!=0) cplM->Data= 0.;
    matM.clear();
    const ExtIdxDescCL& Xidx= RowIdx.GetXidx();
    const ExtIdxDescCL& oldXidx= ColIdx.GetXidx();
    const IdxT num_unks=  RowIdx.NumUnknowns(),
        num_cols=  ColIdx.NumUnknowns();
    MatrixBuilderCL M(&matM, num_unks,  num_cols);//mass matrix
    const MultiGridCL& mg= this->GetMG();
    BndDataCL<> Bndlset(mg.GetBnd().GetNumBndSeg());
    const Uint lvl= RowIdx.TriangLevel();
    LocalNumbP1CL n;
    bool sign[4], oldsign[4], no_newcut, no_oldcut;
    // The 16 products of the P1-shape-functions

    P1FEGridfunctions p1feq;
    TransformedP1FiniteElement & transfp1fel = sdstab_?
        *new StabilizedTransformedP1FiniteElement (p1feq,sdstab_) :
        *new TransformedP1FiniteElement (p1feq);
    GlobalConvDiffReacCoefficients global_cdcoef(D_,H_, GetVelocity() , f_, c_, time);

    DROPS_FOR_TRIANG_TETRA( MG_, lvl, sit) {
        n.assign( *sit, RowIdx, Bndt_);
        InterfaceTetraCL cut, oldcut;
        cut.Init( *sit, lset_,0.);
        oldcut.Init( *sit, oldlset_,0.);
        no_newcut=!cut.Intersects();
        no_oldcut=!oldcut.Intersects();

        LocalConvDiffReacCoefficients local_cdcoef(global_cdcoef,*sit);
        bool pPart_old= (oldcut.GetSign( 0) == 1);
        bool pPart_new= (cut.GetSign( 0) == 1);
        transfp1fel.SetLocal(*sit,local_cdcoef,pPart_new);

        Elmat4x4 M_P1NEW_P1OLD,                   ///< (test FEM, shape FEM)
            M_P1NEW_XOLD,                    ///< (test FEM, shape old XFEM) [at least old interface]
            M_XNEW_P1OLD,                    ///< (test new XFEM, shape FEM) [at least new interface]
            M_XNEW_XOLD;                     ///< (test new XFEM, shape old XFEM) [two interfaces]

        std::memset( M_P1NEW_P1OLD, 0, 4*4*sizeof(double));
        std::memset( M_XNEW_XOLD  , 0, 4*4*sizeof(double));
        std::memset( M_XNEW_P1OLD , 0, 4*4*sizeof(double));
        std::memset( M_P1NEW_XOLD , 0, 4*4*sizeof(double));

        //  for debug purposes you should use this variant instead of the active one, s.t. only the method
        //  SetupLocalTwoInterfacesMassMatrix is
        //  involved in the setup of cutted elements and not SetupLocalOneInterfaceMassMatrix (...) as well...
        //  if(!no_oldcut||!no_newcut){ //new or old interface does not cut
        //    LocalP2CL<> lp2_oldlset(*sit, oldlset_, Bndlset);
        //    // couplings between XFEM basis functions wrt old and new interfaces
        //    SetupLocalTwoInterfacesMassMatrix( cut, oldcut, M_P1NEW_P1OLD, M_P1NEW_XOLD,
        //                                       M_XNEW_P1OLD, M_XNEW_XOLD, transfp1fel, H_, lp2_oldlset);
        //  }
        //  else
        //    SetupLocalOnePhaseMassMatrix ( M_P1NEW_P1OLD, transfp1fel, H_, pPart_new);

        if(no_oldcut){ //old interface does not cut
            if (no_newcut){ //new and old interface do not cut
                SetupLocalOnePhaseMassMatrix ( M_P1NEW_P1OLD, transfp1fel, H_, pPart_new);
            }
            else
            { //new interface does cut, old does not
                Elmat4x4 M_XNEW_P1OLD_n, M_XNEW_P1OLD_p;
                std::memset( M_XNEW_P1OLD_n,0, 4*4*sizeof(double));
                std::memset( M_XNEW_P1OLD_p,0, 4*4*sizeof(double));
                SetupLocalOneInterfaceMassMatrix( cut, /*cut_is_new_cut*/ true,
                                                  M_XNEW_P1OLD_n, M_XNEW_P1OLD_p,
                                                  transfp1fel, sign, H_, /*pPart_nocut*/ pPart_old);
                for(int i= 0; i < 4; ++i){
                    for(int j= 0; j < 4; ++j){
                        M_XNEW_P1OLD[i][j]= sign[i]? -M_XNEW_P1OLD_n[i][j] : M_XNEW_P1OLD_p[i][j];
                        M_P1NEW_P1OLD[i][j]= M_XNEW_P1OLD_n[i][j] + M_XNEW_P1OLD_p[i][j];
                    }
                }
            }
        }
        // the old interface cuts the tetra
        else
        {
            Elmat4x4 M_P1NEW_XOLD_n, M_P1NEW_XOLD_p;
            std::memset( M_P1NEW_XOLD_n,0, 4*4*sizeof(double));
            std::memset( M_P1NEW_XOLD_p,0, 4*4*sizeof(double));
            // couplings between standard basis functions and XFEM basis functions wrt old interface
            if (no_newcut){
                SetupLocalOneInterfaceMassMatrix( oldcut, /*cut_is_new_cut*/ false,
                                                  M_P1NEW_XOLD_n,  M_P1NEW_XOLD_p,
                                                  transfp1fel,oldsign, H_, /*pPart_nocut*/ pPart_new);
                for(int i= 0; i < 4; ++i){
                    for(int j= 0; j < 4; ++j){
                        M_P1NEW_P1OLD[i][j]= M_P1NEW_XOLD_n[i][j] +  M_P1NEW_XOLD_p[i][j];
                        M_P1NEW_XOLD[i][j]= oldsign[j] ? - M_P1NEW_XOLD_n[i][j] :   M_P1NEW_XOLD_p[i][j];
                    }
                }
            }
            // both interfaces cut the tetra
            else {
                LocalP2CL<> lp2_oldlset(*sit, oldlset_, Bndlset);
                // couplings between XFEM basis functions wrt old and new interfaces
                SetupLocalTwoInterfacesMassMatrix( cut, oldcut, M_P1NEW_P1OLD, M_P1NEW_XOLD,
                                                   M_XNEW_P1OLD, M_XNEW_XOLD, transfp1fel,
                                                   H_, lp2_oldlset);
            }
        }


        for(int i= 0; i < 4; ++i)
            if (n.WithUnknowns( i)){
                for(int j= 0; j < 4; ++j)
                    if (n.WithUnknowns( j)) {
                        M( n.num[i], n.num[j])+= M_P1NEW_P1OLD[i][j];
                    }
                    else if (cplM!=0){
                        const double val= Bndt_.GetBndFun( n.bndnum[j])( sit->GetVertex( j)->GetCoord(), time);
                        cplM->Data[n.num[i]]-= M_P1NEW_P1OLD[i][j]*val;
                    }
            }
        if (no_newcut && no_oldcut) continue;
        for(int i= 0; i < 4; ++i)
            if(n.WithUnknowns(i)){
                const IdxT xidx_i= Xidx[n.num[i]];
                for(int j= 0; j < 4; ++j)
                    if(n.WithUnknowns(j)){
                        const IdxT xidx_j= oldXidx[n.num[j]];
                        if (xidx_j!=NoIdx) // at least old cuts
                            M( n.num[i], xidx_j)+= M_P1NEW_XOLD[i][j];
                        if (xidx_i!=NoIdx) // at least new cuts
                            M( xidx_i, n.num[j])+= M_XNEW_P1OLD[i][j];
                        if (xidx_i!=NoIdx && xidx_j!=NoIdx) //both cut
                            M( xidx_i, xidx_j)+= M_XNEW_XOLD[i][j];
                    }
            }
    }
    M.Build();
    delete &transfp1fel;
}

void TransportP1XCL::SetupInstatMixedMassMatrix(MLMatDescCL& matM,
						VecDescCL& cplM,
						const double time) const
{
    MLMatrixCL::iterator itM = --matM.Data.end();
    MLIdxDescCL::iterator it_row = --matM.RowIdx->end();
    MLIdxDescCL::iterator it_col = --matM.ColIdx->end();
    SetupInstatMixedMassMatrix(*itM, &cplM, *it_row, *it_col, time);
    for (size_t num= 1; num < matM.Data.size(); ++num, --itM, --it_row, --it_col)
        SetupInstatMixedMassMatrix (*itM, 0, *it_row, *it_col, time);
}



/// compute \f$ \left( \int_{\Gamma} [c_T]^2 dx \right)^{\frac12} \f$
double TransportP1XCL::Interface_L2error() const
{
    const IdxDescCL &RowIdx = idx.GetFinest();
    const ExtIdxDescCL& Xidx= RowIdx.GetXidx();
    const Uint lvl= RowIdx.TriangLevel();
    LocalNumbP1CL ln;
    double err_sq=0.;

    DROPS_FOR_TRIANG_TETRA( MG_, /*default level*/lvl, it)
    {
        InterfaceTriangleCL triangle;
        triangle.Init( *it, lset_,0.);
        if (!triangle.Intersects()) continue;
        ln.assign( *it, RowIdx, Bndt_);

        for (int ch= 0; ch < 8; ++ch) {
            triangle.ComputeForChild( ch);
            for (int t= 0; t < triangle.GetNumTriangles(); ++t) {
                static Quad5_2DCL<>  p1[4];
                Quad5_2DCL<>jump_on_Gamma;
                jump_on_Gamma*=0.;
                P1DiscCL::GetP1Basis( p1, &triangle.GetBary( t));
                double det= triangle.GetAbsDet( t);
                for(int i= 0; i < 4; ++i)
                    if(ln.WithUnknowns(i)){
                    const IdxT xidx_i= Xidx[ln.num[i]];
                    if (xidx_i!=NoIdx)
                        jump_on_Gamma+= p1[i]*ct.Data[xidx_i];
                    }
                err_sq+=  Quad5_2DCL<>(jump_on_Gamma*jump_on_Gamma).quad(det);
            }
        }
    }
    return std::sqrt(err_sq);
}

//*****************************************************************************
//                               TransportXRepairCL
//*****************************************************************************
void
TransportXRepairCL::pre_refine ()
{
    p1oldctrepair_= std::unique_ptr<RepairP1CL<double>::type >(
        new RepairP1CL<double>::type( c_.GetMG(), c_.oldct, c_.GetBndData()));
    p1ctrepair_= std::unique_ptr<RepairP1CL<double>::type >(
        new RepairP1CL<double>::type( c_.GetMG(), c_.ct, c_.GetBndData()));
}

void
TransportXRepairCL::post_refine ()
{
    VecDescCL loc_ct;
    IdxDescCL loc_cidx( P1X_FE, c_.GetBndData(), c_.GetXFEMOmitBound());
    VecDescCL loc_oldct;
    IdxDescCL loc_oldcidx( P1X_FE, c_.GetBndData(), c_.GetXFEMOmitBound());
    VecDescCL& ct= c_.ct;
    VecDescCL& oldct= c_.oldct;

    loc_cidx.CreateNumbering( mylvl, c_.GetMG(), c_.GetBndData(), &(c_.GetLevelset()),&(c_.GetLevelsetBnd()));
    loc_oldcidx.CreateNumbering( mylvl, c_.GetMG(), c_.GetBndData(),  &(c_.GetOldLevelset()),&(c_.GetLevelsetBnd()));
    loc_ct.SetIdx( &loc_cidx);
    loc_oldct.SetIdx( &loc_oldcidx);
    p1oldctrepair_->repair( loc_oldct);
    p1oldctrepair_->repair( loc_ct);
    double t = ct.t;
    ct.Clear(t);
    c_.DeleteNumbering( &c_.idx);
    c_.idx.GetFinest().swap( loc_cidx);
    ct.SetIdx( &c_.idx);
    ct.Data= loc_ct.Data;

    double oldt = oldct.t;
    oldct.Clear(oldt);
    c_.DeleteNumbering( &c_.oldidx);
    c_.oldidx.GetFinest().swap( loc_oldcidx);
    oldct.SetIdx( &c_.oldidx);
    oldct.Data= loc_oldct.Data;
}

void
  TransportXRepairCL::pre_refine_sequence ()
{
    oldp1xrepair_= std::unique_ptr<P1XRepairCL>( new P1XRepairCL( c_.GetMG(), c_.oldct));
}

void
  TransportXRepairCL::post_refine_sequence ()
{
     c_.CreateNumbering(mylvl, &c_.idx, &c_.oldidx, (c_.GetLevelset()), (c_.GetOldLevelset()));
    (*oldp1xrepair_)();
     oldp1xrepair_.reset();
     c_.ct.SetIdx( &c_.idx);
}

const IdxDescCL* TransportXRepairCL::GetIdxDesc() const {
  throw DROPSErrCL( "TransportXRepairCL::GetIdxDesc: Sorry, not yet implemented.");
  return 0;
}

const VectorCL* TransportXRepairCL::GetVector() const {
  throw DROPSErrCL( "TransportXRepairCL::GetVector: Sorry, not yet implemented.");
  return 0;
}

void TransportXRepairCL::swap( IdxDescCL&, VectorCL&) {
  throw DROPSErrCL( "TransportXRepairCL::swap: Sorry, not yet implemented.");
}


const IdxDescCL* VelTranspRepairCL::GetIdxDesc() const {
  throw DROPSErrCL( "VelTranspRepairCL::GetIdxDesc: Sorry, not yet implemented.");
  return 0;
}

const VectorCL* VelTranspRepairCL::GetVector() const {
  throw DROPSErrCL( "VelTranspRepairCL::GetVector: Sorry, not yet implemented.");
  return 0;
}

void VelTranspRepairCL::swap( IdxDescCL&, VectorCL&) {
  throw DROPSErrCL( "VelTranspRepairCL::swap: Sorry, not yet implemented.");
}


void
  VelTranspRepairCL::pre_refine ()
{
    p2repair_= std::unique_ptr<RepairP2CL<Point3DCL>::type >(
        new RepairP2CL<Point3DCL>::type( mg_, v_, Bnd_v_));
}

void
  VelTranspRepairCL::post_refine ()
{
    VelVecDescCL loc_v;
    IdxDescCL    loc_vidx(vecP2_FE);
    VelVecDescCL& v= v_;
    Uint LastLevel= mg_.GetLastLevel();
    loc_vidx.CreateNumbering( mg_.GetLastLevel(), mg_, Bnd_v_);
    if (LastLevel != v.RowIdx->TriangLevel()) {
        std::cout << "LastLevel: " << LastLevel
                  << " old v->TriangLevel: " << v.RowIdx->TriangLevel() << std::endl;
        throw DROPSErrCL( "VelTranspRepairCL::post_refine: Sorry, not yet implemented.");
    }
    loc_v.SetIdx( &loc_vidx);
/*
    const Uint old_idx= v_.RowIdx->GetIdx();
    const Uint idx= loc_vidx.GetIdx();
    std::cout << "&v = " << &v << std::endl;
    std::cout << "&Bnd_v_ = " << &Bnd_v_ << std::endl;
    std::cout << "&mg_ = " << &mg_ << std::endl;
    const_DiscVelSolCL p2eval_temp( &v, &Bnd_v_, &mg_);
    Uint nedg = 0;
    Uint tl = v_.GetLevel();
    mg_.IsSane(std::cout,tl);
    std::ofstream fouta("aa.out");
    std::ofstream foutb("bb.out");
    std::cout << "dist= " << std::distance(mg_.GetAllEdgeBegin( tl),mg_.GetAllEdgeEnd( tl)) << std::endl;
        for (MultiGridCL::const_EdgeIterator sit= mg_.GetAllEdgeBegin( tl),
         theend= mg_.GetAllEdgeEnd( tl); sit!=theend; ++sit) {

          std::cout << "\r edge= " << nedg++ << "\t" << std::flush ;
          std::cout << "sit = " << &*sit << "\t" << std::flush ;
          std::cout << "sit->GetMidVertex() = " << &*sit->GetMidVertex() << "\t" << std::flush ;

          std::cout << "sit->GetVertex(0) = " << &*sit->GetVertex(0) << "\t" << std::flush ;
          std::cout << "sit->GetVertex(1) = " << &*sit->GetVertex(1) << "\t" << std::flush ;
          if (&*sit->GetMidVertex()){
            std::cout << "coord = " << sit->GetMidVertex()->GetCoord() << std::endl << std::flush;
            fouta << sit->GetMidVertex()->GetCoord() << std::endl;
          }
          else
          {
            if (sit->IsRefined())
              std::cout << "|||||||||||||||||||||||||||||||||||||||||" << std::endl << std::flush;
            std::cout << "coord = NONONO " << std::endl << std::flush;
            static int asdf = 0;
            if (asdf++ == 0){
              foutb << sit->GetVertex(0)->GetCoord() << std::endl;
              foutb << sit->GetVertex(1)->GetCoord() << std::endl << std::endl;
            }
          }
          if (p2eval_temp.IsDefinedOn(*sit)){

            std::cout <<  "IsDefinedOn " << std::endl;
//            std::cout <<  p2eval_temp.val( *sit) << std::endl;

          }

          if (sit->IsRefined()
              && sit->GetMidVertex()->Unknowns.Exist()
              && !sit->GetMidVertex()->Unknowns.Exist( old_idx)
              && sit->GetMidVertex()->Unknowns.Exist( idx)) {
                std::cout << " -nerv start" << std::endl << std::flush;

  //return c[0] * FE_P2CL::H0( v1) + c[1] * FE_P2CL::H1( v1) + c[2] * FE_P2CL::H2( v1);
              if (p2eval_temp.IsDefinedOn(*sit)){
                std::cout <<  "p2eval_temp " << std::endl << std::flush;
                std::cout <<  p2eval_temp.val( *sit) << std::endl << std::flush;

              }

  //              std::cout <<  *sit->GetMidVertex() << std::endl;
       //         f.SetDoF( *sit->GetMidVertex(), old_f.val( *sit));
                std::cout << " -nerv end" << std::endl << std::flush;
          }
          else if (sit->Unknowns.Exist()
                   && sit->Unknowns.Exist( old_idx)
                   && sit->Unknowns.Exist( idx)) {
         //         f.SetDoF( *sit, old_f.val( *sit));
        //          ++counter3;
          }
         }
  */
    p2repair_->repair( loc_v);

    double t = v.t;
    v.Clear(t);
    (*v.RowIdx).DeleteNumbering( mg_);

    vidx_.swap( loc_vidx);
    v.SetIdx( &vidx_);
    v.Data= loc_v.Data;
}

} // end of namespace DROPS
