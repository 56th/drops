/// \file
/// \brief Discretization for PDEs on an interface.
/// \author LNM RWTH Aachen: Joerg Grande

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

/// Implementation of Setup-Routines and class-methods.

#include "surfactant/ifacetransp.h"
#include "levelset/levelset.h"
#include "num/spmat.h"
#include "misc/scopetimer.h"
#include "ifacetransp.h"
#include <cstring>
#include <cmath>


namespace DROPS {

void Extend (const MultiGridCL& mg, const VecDescCL& x, VecDescCL& xext)
{
    const Uint xidx(x.RowIdx->GetIdx()),
        xextidx(xext.RowIdx->GetIdx()),
        lvl(x.RowIdx->TriangLevel());
    xext.Data= 0.;
    const int stride= x.RowIdx->NumUnknownsVertex();

    DROPS_FOR_TRIANG_CONST_VERTEX(mg, lvl, it) {
        if (it->Unknowns.Exist(xidx) && it->Unknowns.Exist(xextidx))
            for (int k=0; k<stride; ++k)
                xext.Data[it->Unknowns(xextidx)+k]= x.Data[it->Unknowns(xidx)+k];
    }
    if (x.RowIdx->GetFE() == P1IF_FE)
        return;

    // For P2IF_FE, also fixup the edge-dofs.
    DROPS_FOR_TRIANG_CONST_EDGE(mg, lvl, it) {
        if (it->Unknowns.Exist(xidx) && it->Unknowns.Exist(xextidx))
            for (int k=0; k<stride; ++k)
                xext.Data[it->Unknowns(xextidx)+k]= x.Data[it->Unknowns(xidx)+k];
    }
}

    void ExtendP2 (const MultiGridCL& mg, const VecDescCL& x, VecDescCL& xext)
//extend the piecewise linear function to a piecewise P2 function in a trivial way
    {
        const Uint xidx( x.RowIdx->GetIdx()),
                xextidx( xext.RowIdx->GetIdx()),
                lvl( x.RowIdx->TriangLevel());
        xext.Data= 1000;

        DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
            if (it->Unknowns.Exist( xidx) && it->Unknowns.Exist( xextidx))
                xext.Data[it->Unknowns( xextidx)]= x.Data[it->Unknowns( xidx)];
        }
        DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it) {
            xext.Data[it->Unknowns( xextidx)]=(xext.Data[it->GetVertex(0)->Unknowns( xextidx)]
                                               +xext.Data[it->GetVertex(1)->Unknowns( xextidx)])/2;

        }

    }

///\todo extend to vector-valued FE
void Restrict (const MultiGridCL& mg, const VecDescCL& xext, VecDescCL& x)
{
    const Uint xidx(x.RowIdx->GetIdx()),
        xextidx(xext.RowIdx->GetIdx()),
        lvl(x.RowIdx->TriangLevel());

    DROPS_FOR_TRIANG_CONST_VERTEX(mg, lvl, it) {
        if (it->Unknowns.Exist(xidx) && it->Unknowns.Exist(xextidx))
            x.Data[it->Unknowns(xidx)]= xext.Data[it->Unknowns(xextidx)];
    }
    if (x.RowIdx->GetFE() == P1IF_FE)
        return;

    // For P2IF_FE, also fixup the edge-dofs.
    DROPS_FOR_TRIANG_CONST_EDGE(mg, lvl, it) {
        if (it->Unknowns.Exist(xidx) && it->Unknowns.Exist(xextidx))
            x.Data[it->Unknowns(xidx)]= xext.Data[it->Unknowns(xextidx)];
    }
}


    void MassMultiply(const MultiGridCL& mg, const VecDescCL& xext, VecDescCL& x, const VecDescCL& ls, const BndDataCL<>& lsetbnd)
    {
        const Uint xidx( x.RowIdx->GetIdx()),
        //   xextidx( xext.RowIdx->GetIdx()),
                lvl( x.RowIdx->TriangLevel());

        x.Data=0.0;

        //   IdxT num[4];
        LocalP1CL<> p1[4];
        p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; // P1-Basis-Functions
        Quad5_2DCL<double> q[4], m, qf;



        InterfaceTriangleCL triangle;
        DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
            triangle.Init( *it, ls, lsetbnd);
            if (triangle.Intersects()) { // We are at the phase boundary.
                //  GetLocalNumbP1NoBnd( num, *it, *x.RowIdx);
                LocalP2CL<> p2(*it,xext,lsetbnd);
                for (int ch= 0; ch < 8; ++ch)
                    //int ch=8;
                {
                    triangle.ComputeForChild( ch);
                    for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri)
                    {
                        double det=triangle.GetAbsDet( tri);
                        qf.assign( p2, &triangle.GetBary( tri));
                        for (int i= 0; i < 4; ++i)
                            q[i].assign( p1[i], &triangle.GetBary( tri));

                        for (int i= 0; i < 4; ++i) {
                            //  if (num[i] == NoIdx) continue;
                            if (it->GetVertex(i)->Unknowns.Exist( xidx))
                            {
                                m= qf*q[i];
                                //  x.Data[num[i]]+=m.quad(det);
                                x.Data[it->GetVertex(i)->Unknowns( xidx)]+=m.quad( det);
                            }
                            // if (it->GetVertex(0)->Unknowns.Exist( xidx))
                            // 	 x.Data[it->GetVertex(0)->Unknowns( xidx)]= xext.Data[it->Unknowns( xextidx)];
                        }
                    }
                    //SetupInterfaceRhsP1OnTriangle( p1, q, v->Data, num,*it, &triangle.GetBary( tri), triangle.GetAbsDet( tri), f);
                }
            }
        }


    }

    void MassMultiply(const MultiGridCL& mg, const VecDescCL& xext, VecDescCL& x, const VecDescCL& ls, const BndDataCL<>& lsetbnd,VecDescCL& v,const BndDataCL<Point3DCL>& vbnd,double dt)
    {
        const Uint xidx( x.RowIdx->GetIdx()),
        //   xextidx( xext.RowIdx->GetIdx()),
                lvl( x.RowIdx->TriangLevel());

        x.Data=0.0;

        //   IdxT num[4];
        LocalP1CL<> p1[4];
        p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; // P1-Basis-Functions
        Quad5_2DCL<double> q[4], m, qf;
        Quad5_2DCL<Point3DCL> qv;

        //  P2EvalCL func_xext = make_P2Eval(mg,lsetbnd,xext);//!!!

        InterfaceTriangleCL triangle;
        DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
            triangle.Init( *it, ls, lsetbnd);
            if (triangle.Intersects()) { // We are at the phase boundary.
                //  GetLocalNumbP1NoBnd( num, *it, *x.RowIdx);
                LocalP2CL<Point3DCL> pv2(*it,v,vbnd);

                LocalP2CL<> p2(*it,xext,lsetbnd);
                World2BaryCoordCL w2b_gy(*it);//search neighbours!!!!!!!!!
                for (int ch= 0; ch < 8; ++ch)
                    //int ch=8;
                {
                    triangle.ComputeForChild( ch);
                    for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri)
                    {
                        double det=triangle.GetAbsDet( tri);
                        const BaryCoordCL* baryco=&triangle.GetBary( tri);
                        qv.assign( pv2, baryco);
                        Point3DCL pt[3];
                        for (Uint k= 0; k < 3; ++k)
                            pt[k]= GetWorldCoord(*it, baryco[k]);
                        for (Uint i= 0; i < Quad5_2DDataCL::NumNodesC; ++i)
                        {

                            BaryCoordCL bary0 = w2b_gy(Quad5_2DDataCL::Node[i][0]*pt[0]+Quad5_2DDataCL::Node[i][1]*pt[1]+Quad5_2DDataCL::Node[i][2]*pt[2]-qv[i]*dt);
                            //   for (Uint k= 0; k < 4; ++k)
                            //	   bary0[k]=std::min(std::max(bary0[k],0.0),1.0);///?????????
                            qf[i]= p2(bary0) ;
                        }
                        for (int i= 0; i < 4; ++i)
                            q[i].assign( p1[i], &triangle.GetBary( tri));

                        for (int i= 0; i < 4; ++i) {
                            //  if (num[i] == NoIdx) continue;
                            if (it->GetVertex(i)->Unknowns.Exist( xidx))
                            {
                                m= qf*q[i];
                                //  x.Data[num[i]]+=m.quad(det);
                                x.Data[it->GetVertex(i)->Unknowns( xidx)]+=m.quad( det);
                            }
                            // if (it->GetVertex(0)->Unknowns.Exist( xidx))
                            // 	 x.Data[it->GetVertex(0)->Unknowns( xidx)]= xext.Data[it->Unknowns( xextidx)];
                        }
                    }
                    //SetupInterfaceRhsP1OnTriangle( p1, q, v->Data, num,*it, &triangle.GetBary( tri), triangle.GetAbsDet( tri), f);
                }
            }
        }


    }
    void MassMultiply1(const MultiGridCL& mg, const VecDescCL& xext, VecDescCL& x, const VecDescCL& ls, const BndDataCL<>& lsetbnd,VecDescCL& v,const BndDataCL<Point3DCL>& vbnd,double dt)
// instat_vector_fun_ptr normal,double t
    {
        const Uint xidx( x.RowIdx->GetIdx()),
        //   xextidx( xext.RowIdx->GetIdx()),
                lvl( x.RowIdx->TriangLevel());

        x.Data=0.0;

        //   IdxT num[4];
        LocalP1CL<> p1[4];
        p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; // P1-Basis-Functions
        Quad5_2DCL<double> q[4], m, qf;
        Quad5_2DCL<Point3DCL> qv;
        Point3DCL normal;
        //   Quad5_2DCL<Point3DCL> qnormal;

        //  P2EvalCL func_xext = make_P2Eval(mg,lsetbnd,xext);//!!!

        InterfaceTriangleCL triangle;
        DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
            triangle.Init( *it, ls, lsetbnd);
            if (triangle.Intersects()) { // We are at the phase boundary.
                //  GetLocalNumbP1NoBnd( num, *it, *x.RowIdx);
                //   LocalP2CL<Point3DCL> normal_lvl(*it,normal,t);
                LocalP2CL<Point3DCL> pv2(*it,v,vbnd);
                //  std::cout<<pv2[0]<<" : ";
                // pv2=dot(pv2,normal_lvl)*normal_lvl;
                //  std::cout<<normal_lvl[0]<<"    "<<pv2[0]<<std::endl;
                LocalP2CL<> p2(*it,xext,lsetbnd);
                World2BaryCoordCL w2b_gy(*it);//need search neighbours!!!!!!!!!
                for (int ch= 0; ch < 8; ++ch)
                    //int ch=8;
                {
                    triangle.ComputeForChild( ch);
                    for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri)
                    {
                        double det=triangle.GetAbsDet( tri);
                        const BaryCoordCL* baryco=&triangle.GetBary( tri);
                        normal=triangle.GetNormal();
                        qv.assign( pv2, baryco);
                        // 	   qnormal.assign( normal_lvl, baryco);
                        // qv=dot(qv,qnormal)*qnormal;
                        qv=dot(normal,qv)*normal;
                        Point3DCL pt[3];
                        for (Uint k= 0; k < 3; ++k)
                            pt[k]= GetWorldCoord(*it, baryco[k]);
                        for (Uint i= 0; i < Quad5_2DDataCL::NumNodesC; ++i)
                        {

                            BaryCoordCL bary0 = w2b_gy(Quad5_2DDataCL::Node[i][0]*pt[0]+Quad5_2DDataCL::Node[i][1]*pt[1]+Quad5_2DDataCL::Node[i][2]*pt[2]-qv[i]*dt);
                            qf[i]= p2(bary0) ;
                            Uint k=0;
                            while((bary0[0]<0||bary0[0]>1||bary0[1]<0||bary0[1]>1||bary0[2]<0||bary0[2]>1||bary0[3]<0||bary0[3]>1)&&k<5)
                            {
                                if(k<4){
                                    World2BaryCoordCL w2b_gy(*it->GetNeighbor(k));
                                    bary0 = w2b_gy(Quad5_2DDataCL::Node[i][0]*pt[0]+Quad5_2DDataCL::Node[i][1]*pt[1]+Quad5_2DDataCL::Node[i][2]*pt[2]-qv[i]*dt);
                                }

                                k++;
                            }
                            if(k>1&&k<5)
                            {   LocalP2CL<> p2(*it->GetNeighbor(k-1),xext,lsetbnd);
                                qf[i]= p2(bary0) ;
                                //   std::cout<<"o";
                            }
                            // else
                            //   {if(k>=5)	    std::cout<<"*";}
                            /*for (Uint k= 0; k < 4; ++k)
                            {if(bary0[k]<0||bary0[k]>1)
                                std::cout<<"*";
                            bary0[k]=std::min(std::max(bary0[k],0.0),1.0);
                            }///?????????}*/

                        }
                        for (int i= 0; i < 4; ++i)
                            q[i].assign( p1[i], &triangle.GetBary( tri));

                        for (int i= 0; i < 4; ++i) {
                            //  if (num[i] == NoIdx) continue;
                            if (it->GetVertex(i)->Unknowns.Exist( xidx))
                            {
                                m= qf*q[i];
                                //  x.Data[num[i]]+=m.quad(det);
                                x.Data[it->GetVertex(i)->Unknowns( xidx)]+=m.quad( det);
                            }
                            // if (it->GetVertex(0)->Unknowns.Exist( xidx))
                            // 	 x.Data[it->GetVertex(0)->Unknowns( xidx)]= xext.Data[it->Unknowns( xextidx)];
                        }
                    }
                    //SetupInterfaceRhsP1OnTriangle( p1, q, v->Data, num,*it, &triangle.GetBary( tri), triangle.GetAbsDet( tri), f);
                }
            }
        }


    }


void GetLocalNumbInterface(IdxT* Numb, const TetraCL& s, const IdxDescCL& idx)
{
    const Uint sys= idx.GetIdx();
    for (Uint i= 0; i < 4; ++i)
        Numb[i]= s.GetVertex(i)->Unknowns.Exist(sys) ? s.GetVertex(i)->Unknowns(sys) : NoIdx;
    if (idx.GetFE() == P1IF_FE)
        return;

    for(Uint i= 0; i < 6; ++i)
        Numb[i+4]= s.GetEdge(i)->Unknowns.Exist(sys) ? s.GetEdge(i)->Unknowns(sys) : NoIdx;
}


InterfaceCommonDataP2CL::InterfaceCommonDataP2CL (const VecDescCL& ls_arg, const BndDataCL<>& lsetbnd_arg,
    const QuaQuaMapperCL& quaquaarg, const PrincipalLatticeCL& lat_arg)
    : ls(&ls_arg), lsetbnd(&lsetbnd_arg), lat(&lat_arg),  compute_quaddomains_(true), ls_loc(lat->vertex_size()),
      quaqua(quaquaarg)
{
    qdom_projected.compute_absdets(true);
    P2DiscCL::GetGradientsOnRef(gradrefp2);
    for (Uint i= 0; i < 10 ; ++i)
        p2[i][i]= 1.; // P2-Basis-Functions
}

InterfaceCommonDataDeformP2CL::InterfaceCommonDataDeformP2CL (const VecDescCL& ls_arg, const BndDataCL<>& lsetbnd_arg,
        VecDescCL& Psi_vdarg, const PrincipalLatticeCL& lat_arg)
    : ls(&ls_arg), lsetbnd(&lsetbnd_arg), lat(&lat_arg), Psi_vd (&Psi_vdarg), ls_loc(lat->vertex_size()),
      Phi (this)
{
    P2DiscCL::GetGradientsOnRef(gradrefp2);
    for (Uint i= 0; i < 10 ; ++i)
        p2[i][i]= 1.; // P2-Basis-Functions
}


    void update_global_matrix_P1 (MatrixBuilderCL& M, const double coup[4][4], const IdxT numr[4], const IdxT numc[4])
    {
        for (int i= 0; i < 4; ++i)
            if (numr[i] != NoIdx)
                for (int j= 0; j < 4; ++j)
                    if (numc[j] != NoIdx)
                        M( numr[i], numc[j])+= coup[i][j];
    }

void SetupInterfaceMassP1 (const MultiGridCL& mg, MatDescCL* matM, const VecDescCL& ls, const BndDataCL<>& lsetbnd, double alpha)
{
    //ScopeTimerCL timer("SetupInterfaceMassP1");

    TetraAccumulatorTupleCL accus;
    InterfaceCommonDataP1CL cdata(ls, lsetbnd);
    accus.push_back(&cdata);
    InterfaceMatrixAccuCL<LocalInterfaceMassP1CL, InterfaceCommonDataP1CL> accu(matM, LocalInterfaceMassP1CL(alpha), cdata);
    accus.push_back(&accu);
    const IdxDescCL* RowIdx= matM->RowIdx;
    accumulate(accus, mg, RowIdx->TriangLevel(), RowIdx->GetBndInfo());

    // WriteToFile(matM->Data, "mass.txt", "mass");
}

void SetupInterfaceMassP1OnTriangle (const LocalP1CL<> p1[4], Quad5_2DCL<> q[4],
        const BaryCoordCL triangle[3], double det, double coup[4][4])
{
    for (int i= 0; i < 4; ++i)
        q[i].assign(p1[i], triangle);

    Quad5_2DCL<> m;
    for (int i= 0; i < 4; ++i) {
        m= q[i]*q[i];
        coup[i][i]+= m.quad(det);
        for(int j= 0; j < i; ++j) {
            m= (q[j]*q[i]);
            coup[i][j]+= m.quad(det);
            coup[j][i]= coup[i][j];
        }
    }
}

void SetupInterfaceMassP1X (const MultiGridCL& MG, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd, const double alpha[2])
/// \todo: merge with  SetupInterfaceMassP1 by using accumulator pattern (difference only in update_global)
{
    const IdxT num_rows= mat->RowIdx->NumUnknowns();
    const IdxT num_cols= mat->ColIdx->NumUnknowns();
    const bool xfemr= mat->RowIdx->IsExtended(),
            xfemc= mat->ColIdx->IsExtended();
    ExtIdxDescCL &xidxr= mat->RowIdx->GetXidx(),
            &xidxc= mat->ColIdx->GetXidx();
    MatrixBuilderCL M(&mat->Data, num_rows, num_cols);

    // compute coefficients
    double Kneg, Kpos, Kboth, K[4];
    if (alpha) {
        Kpos= alpha[0];
        Kneg= alpha[1];
    } else
        Kneg= Kpos= 1.;
    Kboth= Kneg+Kpos;

    const Uint lvl= mat->GetRowLevel();
    IdxT numr[4], numc[4], xnumr[4], xnumc[4];
    bool is_pos[4];

    double coup[4][4];
    LocalP1CL<> p1[4];
    p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; // P1-Basis-Functions
    Quad5_2DCL<double> q[4], m;

    InterfaceTriangleCL triangle;
    DROPS_FOR_TRIANG_CONST_TETRA(MG, lvl, it) {
        // setup local matrix
        triangle.Init(*it, ls, lsetbnd);
        if (triangle.Intersects()) { // We are at the phase boundary.
            std::memset(coup, 0, 4*4*sizeof(double));

            for (int ch= 0; ch < 8; ++ch) {
                if (!triangle.ComputeForChild(ch)) // no patch for this child
                    continue;

                for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri)
                    SetupInterfaceMassP1OnTriangle(p1, q, &triangle.GetBary(tri), triangle.GetAbsDet(tri), coup);
            }

            // update global matrix

            GetLocalNumbP1NoBnd(numr, *it, *mat->RowIdx);
            GetLocalNumbP1NoBnd(numc, *it, *mat->ColIdx);
            for(int k= 0; k < 4; ++k) {
                xnumr[k]= xfemr ? xidxr[numr[k]] : NoIdx;
                xnumc[k]= xfemc ? xidxc[numc[k]] : NoIdx;
                is_pos[k]= triangle.GetSign(k)==1;
                K[k]= is_pos[k] ? -Kneg : Kpos;
            }

            for(int i= 0; i < 4; ++i) {// assemble row Numb[i]
                if (numr[i] == NoIdx) continue;
                for(int j= 0; j < 4; ++j) {
                    if (numc[j] == NoIdx) continue;
                    M(numr[i],   numc[j])+= Kboth*coup[j][i];
                    if (xnumc[j] == NoIdx) continue;
                    M(numr[i],   xnumc[j])+= K[j]*coup[j][i];
                }
                if (xnumr[i] == NoIdx) continue;
                for(int j= 0; j < 4; ++j) {
                    if (numc[j] == NoIdx) continue;
                    M(xnumr[i],   numc[j])+= K[i]*coup[j][i];
                    if (xnumc[j] == NoIdx || is_pos[i]!= is_pos[j]) continue;
                    M(xnumr[i],   xnumc[j])+= std::abs(K[i])*coup[j][i];
                }
            }
        }
    }

    M.Build();
}

void SetupInterfaceSorptionP1X (const MultiGridCL& MG, const VecDescCL& ls, const BndDataCL<>& lsetbnd,
        MatDescCL* R, MatDescCL* C, MatDescCL* R_i, MatDescCL* C_i, const IdxDescCL* mass_idx, const IdxDescCL* surf_idx, const double k_a[2], const double k_d[2])
{
    R->SetIdx  (mass_idx, mass_idx);
    C->SetIdx  (mass_idx, surf_idx);
    R_i->SetIdx(surf_idx, surf_idx);
    C_i->SetIdx(surf_idx, mass_idx);

    SetupInterfaceMassP1(MG, R_i, ls, lsetbnd, k_d[0]+k_d[1]);
    SetupInterfaceMassP1X(MG, R, ls, lsetbnd, k_a);

    const double mk_a[2]= { -k_a[0], -k_a[1] },
                 mk_d[2]= { -k_d[0], -k_d[1] };
    SetupInterfaceMassP1X(MG, C, ls, lsetbnd, mk_d);
    SetupInterfaceMassP1X(MG, C_i, ls, lsetbnd, mk_a);
}

void SetupMixedMassP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd)
{
    SetupInterfaceMassP1(mg, mat, ls, lsetbnd);
}

void SetupLBP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd, double D)
{
    //ScopeTimerCL timer("SetuLBP1");

    TetraAccumulatorTupleCL accus;
    InterfaceCommonDataP1CL cdata(ls, lsetbnd);
    accus.push_back(&cdata);
    InterfaceMatrixAccuCL<LocalLaplaceBeltramiP1CL, InterfaceCommonDataP1CL> accu(mat, LocalLaplaceBeltramiP1CL(D), cdata);
    accus.push_back(&accu);
    const IdxDescCL* RowIdx= mat->RowIdx;
    accumulate(accus, mg, RowIdx->TriangLevel(), RowIdx->GetBndInfo());

    // WriteToFile(mat->Data, "lb.txt", "lb");
}

 /// P1
void SetupInterfaceRhsP1OnTriangle (const LocalP1CL<> p1[4],
    Quad5_2DCL<> q[4],VectorCL& v, const IdxT Numb[4],
    const TetraCL& t, const BaryCoordCL triangle[3], double det,
    InstatScalarFunction f,double time=0)
{
    for (int i= 0; i < 4; ++i)
        q[i].assign(p1[i], triangle);
    Quad5_2DCL<> qf(t, triangle, f,time), r;

    for (int i= 0; i < 4; ++i) {
        if (Numb[i] == NoIdx) continue;
        r= qf*q[i];
        v[Numb[i]]+= r.quad(det);
    }
}


void SetupInterfaceRhsP1 (const MultiGridCL& mg, VecDescCL* v,
    const VecDescCL& ls, const BndDataCL<>& lsetbnd, InstatScalarFunction f, double t)
{
    const IdxT num_unks= v->RowIdx->NumUnknowns();
    const Uint lvl = v->GetLevel();

    v->Clear(0);
    IdxT num[4];

    std::cout << "entering SetupInterfaceRhsP1: " << num_unks << " dof... ";

    LocalP1CL<> p1[4];
    p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; // P1-Basis-Functions
    Quad5_2DCL<double> q[4], m;

    InterfaceTriangleCL triangle;

    DROPS_FOR_TRIANG_CONST_TETRA(mg, lvl, it) {
        triangle.Init(*it, ls, lsetbnd);
        if (triangle.Intersects()) { // We are at the phase boundary.
            GetLocalNumbP1NoBnd(num, *it, *v->RowIdx);

            for (int ch= 0; ch < 8; ++ch) {
                triangle.ComputeForChild(ch);
                for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri)
                    SetupInterfaceRhsP1OnTriangle(p1, q, v->Data, num,
                        *it, &triangle.GetBary(tri), triangle.GetAbsDet(tri), f,t);
            }
        }
    }
    std::cout << " Rhs set up." << std::endl;
}

// void SetupInterfaceRhsP1 (const MultiGridCL& mg, VecDescCL* v,
//     const VecDescCL& ls, const BndDataCL<>& lsetbnd, InstatScalarFunction f, double t)
// {
//     TetraAccumulatorTupleCL accus;
//     InterfaceCommonDataP1CL cdata(ls, lsetbnd);
//     accus.push_back(&cdata);
//     InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL> loadaccu(v, LocalVectorP1CL(f, v->t), cdata);
//     accus.push_back(&loadaccu);
//     accumulate(accus, mg, v->RowIdx->TriangLevel(), v->RowIdx->GetBndInfo());

//     // WriteToFile(v->Data, "rhs.txt", "Rhs");
// }

/// SetupStokes
/// \brief Setup of the local Stokes system on a tetra intersected by the dividing surface.
class LocalSurfOseen {
private:
    PrincipalLatticeCL const & lat;
    // shape funcs
    LocalP1CL<> P1Hat[4];
    Point3DCL   P1Grad[4];
    LocalP2CL<> P2Hat[10];
    LocalP1CL<Point3DCL> P2GradRef[10], P2Grad[10];
    SMatrixCL<3,3> P2Hess[10];
    // quadratures
    GridFunctionCL<Point3DCL>
            qP2Normal, // normal computed from P2 interpolant
            qPatchNormal, // normals exact to patch of integration
            qPatchOrP2Normal, // either normal computed from P2 interpolant, or normals exact to patch of integration
            qExactNormal, // exact normal to exact surface
            qrotP1[12], qvP1, qvProjP1,
            qP_row[3], // ith element is ith row (= column) of P
            qF, // moment eqn rhs
            qWind; // wind field
    GridFunctionCL<Point3DCL> qLsGrad, qP1Grad[4], qSurfP1Grad[4], qP2Grad[10], qSurfP2Grad[10], qP2HatCrossN[30];
    GridFunctionCL<> qG, qP1Hat[4], qP2Hat[10], qLsGradNorm, qconvP1[4], qdivP1, qP2NormalComp[3], qExactOrP2NormalComp[3], qvProjP1_comp[3], qvP1_comp[3];
    GridFunctionCL<SMatrixCL<3,3>> qP; // projection
    GridFunctionCL<SMatrixCL<3,3>> qP2E[30]; // surface rate-of-strain stress tensor, qP1E[i] = $E_s(\phi_i)$, $\phi_i$ = ith vector (velocity) shape func
    GridFunctionCL<SMatrixCL<3,3>> qHess, qAnalyticHess;
    QuadDomain2DCL  q2Ddomain;
    std::valarray<double> ls_loc;
    SurfacePatchCL spatch;
    QuadDomainCL q3Ddomain;
    GridFunctionCL<Point3DCL> q3Dnormal, q3DP2Grad[10];
    GridFunctionCL<Point3DCL> q3DP1Grad[4];
    LocalP1CL<Point3DCL> getLevelsetGrad(LocalP2CL<> const &);
    SMatrixCL<3, 3>      getLevelsetHess(LocalP2CL<> const &);
public:
    SurfOseenParam* param;
    LocalSurfOseen(SurfOseenParam* param)
        : param(param)
        , lat(PrincipalLatticeCL::instance(param->input.numbOfVirtualSubEdges))
        , ls_loc(lat.vertex_size()) {
        P1DiscCL::GetP1Basis(P1Hat);
        P2DiscCL::GetP2Basis(P2Hat);
        P2DiscCL::GetGradientsOnRef(P2GradRef);
        std::cout << "numb of virtual subedges = " << param->input.numbOfVirtualSubEdges << std::endl;
    }
    Uint num_intervals() const { return lat.num_intervals(); }
    void exportPatchInfo(std::string const & path, const SMatrixCL<3,3>& T, const LocalP2CL<>& ls, const TetraCL& tet);
    void calcIntegrands(const SMatrixCL<3,3>& T, const LocalP2CL<>& ls, const TetraCL& tet); ///< has to be called before any setup method!
    template<typename Interp>
    void calcIntegrands(const SMatrixCL<3,3>& T, const LocalP2CL<>& ls, Interp const & w, const TetraCL& tet); ///< has to be called before any setup method!
    void calc3DIntegrands(const SMatrixCL<3,3>& T, const LocalP2CL<>& ls, const TetraCL& tet); ///< has to be called after calcIntegrands!
    void setupLB_P2(double LB_P2[10][10]); // laplace-beltrami
    void setupLB_P1(double LB_P1[4][4]);
    void setupA_P2_inconsistent(double A_P2[30][30]); // inconsistent method
    void setupA_P2_consistent(double A_P2[30][30]); // explicit contraction of matrices; consistent formulation
    void setupA_P2_stab(double A_P2_stab[10][10], double absdet);
    void setupA_P1(double A_P1[12][12]);
    void setupA_P1_stab(double A_P1_stab[4][4], double absdet);
    void setupN_P2(double N_P2[30][30]); // convection matrix
    void setupAL_P2(double AL_P2[30][30]); // AL / grad-div stab mtx
    void setupRhsAL_P2(double RhsAL_P2[30]); // AL / grad-div stab rhs
    void setupQ_P1P2(double Q_P1P2[4][30]); // rhs matrix for surface curl projection
    void setupB_P1P2(double B_P1P2[4][30]);
    void setupB_P1P1(double B_P1P1[4][12]);
    void setupB_P2P2(double B_P2P2[10][30]);
    void setupB_P2P1(double B_P2P1[10][12]);
    void setupM_P2(double M_P2[10][10]);
    void setupF_P2(double F_P2[30]);
    void setupM_P1 (double M_P1[4][4]);
    void setupG_P1 (double G_P1[4]);
    void setupS_P2 (double S_P2[30][30]);
    void setupS_P1 (double S_P1[12][12]);
    // Navier-Stokes:
    void setupM_P1_stab (double M_P1_stab[4][4], double absdet);
};

LocalP1CL<Point3DCL> LocalSurfOseen::getLevelsetGrad(const LocalP2CL<>& ls) {
    LocalP1CL<Point3DCL> res;
    for(size_t i = 0; i < 10 ; ++i)
        res += ls[i] * P2Grad[i];
    return res;
}

SMatrixCL<3,3> LocalSurfOseen::getLevelsetHess(const LocalP2CL<>& ls) {
    SMatrixCL<3,3> res;
    for(size_t i = 0; i < 10 ; ++i)
        res += ls[i] * P2Hess[i];
    return res;
}

void LocalSurfOseen::calc3DIntegrands(const SMatrixCL<3,3>& T, const LocalP2CL<>& ls, const TetraCL& tet) {
    make_SimpleQuadDomain<Quad5DataCL> (q3Ddomain, AllTetraC);
    auto Normals = getLevelsetGrad(ls);
    resize_and_evaluate_on_vertexes (Normals, q3Ddomain, q3Dnormal);
    // scale normals accordingly to the Euclidean Norm (only consider the ones which make a contribution in the sense of them being big enough... otherwise one has to expect problems with division through small numbers)
    for(Uint i=0; i<q3Dnormal.size(); ++i)
         q3Dnormal[i]= q3Dnormal[i]/q3Dnormal[i].norm();
    for(int j=0; j<10; ++j) {
        resize_and_evaluate_on_vertexes(P2Grad[j], q3Ddomain, q3DP2Grad[j]);
    }
    for(int j=0; j<4; ++j) {
        q3DP1Grad[j].resize(q3Ddomain.vertex_size());
        q3DP1Grad[j]= P1Grad[j];
    }
}

void LocalSurfOseen::exportPatchInfo(std::string const & path, const SMatrixCL<3,3>& T, const LocalP2CL<>& ls, const TetraCL& tet) {
    auto id = tet.GetId().GetIdent();
    std::ofstream
        vertices    (path + "/" + std::to_string(id) + "_vertices.txt"),
        triangles   (path + "/" + std::to_string(id) + "_triangles.txt"),
        normals     (path + "/" + std::to_string(id) + "_normals.txt"),
        tetrahedron(path + "/" + std::to_string(id) + "_tetrahedron.txt");
    evaluate_on_vertexes(ls, lat, Addr(ls_loc));
    spatch.make_patch<MergeCutPolicyCL>(lat, ls_loc);
    // export tet
    for (auto v = tet.GetVertBegin(); v != tet.GetVertEnd(); ++v)
        tetrahedron << (*v)->GetCoord() << '\n';
    spatch.compute_world_vertexes(tet);
    for (auto v = spatch.world_vertex_begin(); v != spatch.world_vertex_end(); ++v)
        vertices << (*v) << '\n';
    // export triangles
    for (auto t = spatch.facet_begin(); t != spatch.facet_end(); ++t)
        triangles << (*t) << '\n';
    // export normals
    spatch.compute_normals(tet);
    for (auto n = spatch.normal_begin(); n != spatch.normal_end(); ++n)
        normals << (*n) << '\n';
}

void LocalSurfOseen::calcIntegrands(const SMatrixCL<3,3>& T, const LocalP2CL<>& ls, const TetraCL& tet) {
    ++param->output.numbOfCutTetras;
    P2DiscCL::GetGradients(P2Grad, P2GradRef, T);
    P2DiscCL::GetHessians(P2Hess, T);
    P1DiscCL::GetGradients(P1Grad, T);
    evaluate_on_vertexes(ls, lat, Addr(ls_loc));
    spatch.make_patch<MergeCutPolicyCL>(lat, ls_loc);
    // the routine takes the information about the tetrahedra and the cutting surface and generates a two-dimensional triangulation of the cut, including the necessary point-positions and weights for the quadrature
    make_CompositeQuad5Domain2D(q2Ddomain, spatch, tet);
    // THIS CODE MAY BE HANDY LATER!!!
    //    std::cout << "numb of quad nodes: " << q2Ddomain.vertex_size() << '\n';
    //    std::cout << "barycentric:\n";
    //    for (auto v = q2Ddomain.vertex_begin(); v != q2Ddomain.vertex_end(); ++v) {
    //        std::cout << (*v) << '\n';
    //    }
    //    SPatchCL<3>::WorldVertexContT world_quad_nodes;
    //    world_quad_nodes.reserve(q2Ddomain.vertex_size());
    //    SPatchCL<3>::VertexToWorldVertexMapperT mapper(tet);
    //    std::transform(q2Ddomain.vertex_begin(), q2Ddomain.vertex_end(), std::back_inserter<SPatchCL<3>::WorldVertexContT>(world_quad_nodes), mapper);
    //    std::cout << "world:\n";
    //    for (auto v = world_quad_nodes.begin(); v != world_quad_nodes.end(); ++v) {
    //        std::cout << (*v) << '\n';
    //    }
    //    std::cout << '\n';
    // DO NOT DELETE
    // compute rhs
    if (param->input.f != nullptr) resize_and_evaluate_on_vertexes(param->input.f, tet, q2Ddomain, param->input.t, qF);
    if (param->input.g != nullptr) resize_and_evaluate_on_vertexes(param->input.g, tet, q2Ddomain, param->input.t, qG);
    // compute normal
    resize_and_evaluate_on_vertexes(getLevelsetGrad(ls), q2Ddomain, qLsGrad);
    qLsGradNorm   = sqrt(dot(qLsGrad, qLsGrad));
    qP2Normal     = qLsGrad / qLsGradNorm;
    size_t numbOfQuadPoints = qP2Normal.size() / spatch.facet_size();
    std::vector<Point3DCL> qExactNormalInit;
    qExactNormalInit.reserve(numbOfQuadPoints * spatch.facet_size());
    spatch.compute_normals(tet);
    auto sgn = [](double val) { return (0. < val) - (val < 0.); };
    size_t patchIndex = 0;
    for (auto n = spatch.normal_begin(); n != spatch.normal_end(); ++n) {
        auto normal = sgn(inner_prod(*n, qP2Normal[patchIndex * numbOfQuadPoints])) * (*n);
        for (size_t i = 0; i < numbOfQuadPoints; ++i) qExactNormalInit.emplace_back(normal);
        ++patchIndex;
    }
    qPatchNormal.resize(qExactNormalInit.size());
    std::copy(qExactNormalInit.begin(), qExactNormalInit.end(), begin(qPatchNormal));
    auto usePatchNormals = param->input.usePatchNormal;
    qPatchOrP2Normal = usePatchNormals ? qPatchNormal : qP2Normal;
    if (param->input.exactNormal != nullptr) {
        resize_and_evaluate_on_vertexes(param->input.exactNormal, tet, q2Ddomain, param->input.t, qExactNormal);
        param->output.normalErrSq.patch += quad_2D(dot(qExactNormal - qPatchNormal, qExactNormal - qPatchNormal), q2Ddomain);
        param->output.normalErrSq.lvset += quad_2D(dot(qExactNormal - qP2Normal, qExactNormal - qP2Normal), q2Ddomain);
    }
    // compute projector and shape matrix
    qP.resize(q2Ddomain.vertex_size());
    qP = eye<3, 3>() - outer_product(qP2Normal, qP2Normal);
    if (param->input.formulation == SurfOseenParam::Formulation::consistent || param->input.exactShape != nullptr) {
        // qHess = qP;
        qHess.resize(q2Ddomain.vertex_size());
        qHess = getLevelsetHess(ls);
        qHess = qP * (qHess / qLsGradNorm) * qP;
        if (param->input.exactShape != nullptr) {
            resize_and_evaluate_on_vertexes(param->input.exactShape, tet, q2Ddomain, 0., qAnalyticHess);
            param->output.shapeErrSq += quad_2D(contract(qAnalyticHess - qHess, qAnalyticHess - qHess), q2Ddomain);
        }
    }
    if (usePatchNormals) qP = eye<3, 3>() - outer_product(qPatchOrP2Normal, qPatchOrP2Normal);
    // provide all components of the normals
    for (int k = 0; k < 3; ++k) {
        qP2NormalComp[k].resize(q2Ddomain.vertex_size());
        ExtractComponent(qP2Normal, qP2NormalComp[k], k);
        qExactOrP2NormalComp[k].resize(q2Ddomain.vertex_size());
        ExtractComponent(qPatchOrP2Normal, qExactOrP2NormalComp[k], k);
    }
    // resize and evaluate P2 basis functions
    for(int j=0; j<10 ;++j) {
        resize_and_evaluate_on_vertexes(P2Hat[j], q2Ddomain, qP2Hat[j]);
        resize_and_evaluate_on_vertexes(P2Grad[j], q2Ddomain, qP2Grad[j]);
        qSurfP2Grad[j] = qP2Grad[j] - dot(qP2Grad[j], qPatchOrP2Normal) * qPatchOrP2Normal;
    }
    // Resize and evaluate of all the 4 P1 Gradient Functions and apply pointwise projections   (P grad \xi_j  for j=1..4)
    for(int j=0; j<4 ;++j) {
        resize_and_evaluate_on_vertexes(P1Hat[j], q2Ddomain, qP1Hat[j]);
        qP1Grad[j].resize(q2Ddomain.vertex_size());
        qP1Grad[j] = P1Grad[j];
        qSurfP1Grad[j] = qP1Grad[j] - dot(qP1Grad[j], qPatchOrP2Normal) * qPatchOrP2Normal;
    }
    // apply pointwise projection to the 3 std basis vectors e_k
    if (param->input.formulation == SurfOseenParam::Formulation::inconsistent)
        for(int k=0; k<3 ;++k) {
            qP_row[k].resize(q2Ddomain.vertex_size());
            Point3DCL e_k = DROPS::std_basis<3>(k+1);
            qP_row[k] = e_k;
            qP_row[k]-= dot(e_k, qPatchOrP2Normal) * qPatchOrP2Normal;
        }
    // compute surface stress tensor
    auto qVectGrad = [&](size_t vecShapeIndex, GridFunctionCL<Point3DCL>* P1OrP2grad) {
        GridFunctionCL<SMatrixCL<3,3>> res(SMatrixCL<3,3>(), q2Ddomain.vertex_size());
        auto scaShapeIndex = vecShapeIndex / 3;
        auto row = vecShapeIndex - 3 * scaShapeIndex;
        for (size_t i = 0; i < res.size(); ++i) {
            SMatrixCL<3, 3> mtx(0.);
            mtx.col(row, P1OrP2grad[scaShapeIndex][i]);
            assign_transpose(res[i], mtx);
        }
        return res;
    };
    // compute velocity shape func cross normal vect
    auto qHatCrossN = [&](size_t vecShapeIndex, GridFunctionCL<>* P1OrP2Hat) {
        GridFunctionCL<Point3DCL> phi(Point3DCL(0., 0., 0.), q2Ddomain.vertex_size());
        auto is = vecShapeIndex / 3; // scalar shape index
        auto in = vecShapeIndex - 3 * is; // nonzero vect component
        for (size_t i = 0; i < phi.size(); ++i)
            phi[i][in] = P1OrP2Hat[is][i];
        return cross_product(qPatchOrP2Normal, phi);
    };
    for (size_t vecShapeIndex = 0; vecShapeIndex < 30; ++vecShapeIndex) {
        qP2E[vecShapeIndex] = qP * sym_part(qVectGrad(vecShapeIndex, qP2Grad)) * qP;
        qP2HatCrossN[vecShapeIndex] = qHatCrossN(vecShapeIndex, qP2Hat);
    }
}

template<typename Interp>
void LocalSurfOseen::calcIntegrands(const SMatrixCL<3,3>& T, const LocalP2CL<>& ls, Interp const & w, const TetraCL& tet) {
    calcIntegrands(T, ls, tet);
    // compute wind
    resize_and_evaluate_on_vertexes(w, q2Ddomain, qWind);
}

void LocalSurfOseen::setupA_P2_consistent(double A_P2[30][30]) {
    for (size_t i = 0; i < 30; ++i) {
        auto is = i / 3; // scalar shape index
        auto in = i - 3 * is; // nonzero vect component
        for (size_t j = i; j < 30; ++j) {
            auto js = j / 3; // scalar shape index
            auto jn = j - 3 * js; // nonzero vect component
            A_P2[i][j] = quad_2D(contract(qP2E[j] - (qP2Hat[js] * qExactOrP2NormalComp[jn]) * qHess, qP2E[i] - (qP2Hat[is] * qExactOrP2NormalComp[in]) * qHess), q2Ddomain);
            A_P2[j][i] = A_P2[i][j];
        }
    }
}

void LocalSurfOseen::setupLB_P2(double LB_P2[10][10]) {
    for (int i = 0; i < 10; ++i)
        for (int j = i; j < 10; ++j) {
            LB_P2[i][j] = quad_2D(dot(qSurfP2Grad[j], qSurfP2Grad[i]), q2Ddomain);
            LB_P2[j][i] = LB_P2[i][j];
        }
}

void LocalSurfOseen::setupLB_P1(double LB_P1[4][4]) {
    for (int i = 0; i < 4; ++i)
        for (int j = i; j < 4; ++j) {
            LB_P1[i][j] = quad_2D(dot(qSurfP1Grad[j], qSurfP1Grad[i]), q2Ddomain);
            LB_P1[j][i] = LB_P1[i][j];
        }
}

void LocalSurfOseen::setupA_P2_inconsistent(double A_P2[30][30]) {
    for (size_t i = 0; i < 30; ++i)
        for (size_t j = i; j < 30; ++j) {
            A_P2[i][j] = quad_2D(contract(qP2E[j], qP2E[i]), q2Ddomain);
            A_P2[j][i] = A_P2[i][j];
        }
}

void LocalSurfOseen::setupAL_P2(double AL_P2[30][30]) {
    for (size_t i = 0; i < 30; ++i)
        for (size_t j = i; j < 30; ++j) {
            AL_P2[i][j] = quad_2D(trace(qP2E[j]) * trace(qP2E[i]), q2Ddomain);
            AL_P2[j][i] = AL_P2[i][j];
        }
}

void LocalSurfOseen::setupRhsAL_P2(double RhsAL_P2[30]) {
    for (size_t i = 0; i < 30; ++i)
        RhsAL_P2[i] = quad_2D(qG * trace(qP2E[i]), q2Ddomain);
}

// TODO: Den Fall fullGrad testen!
void LocalSurfOseen::setupA_P1 (double A_P1[12][12]) {
    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
            for (int k=0; k<3; ++k) {
                Point3DCL e_k= DROPS::std_basis<3>(k+1);
                for (int l=0; l<3; ++l) {
                    Point3DCL e_l= DROPS::std_basis<3>(l+1);
                    A_P1[3 * i + k][3 * j + l] = quad_2D(dot(e_k, qSurfP1Grad[j]) * dot(e_l, qSurfP1Grad[i]) + dot(qSurfP1Grad[i], qSurfP1Grad[j]) * dot(qP_row[k], qP_row[l]), q2Ddomain);
                }
            }
        }
    }
}

void LocalSurfOseen::setupQ_P1P2(double Q_P1P2[4][30]) {
    for (size_t i = 0; i < 4; ++i)
        for (size_t j = 0; j < 30; ++j)
            Q_P1P2[i][j] = quad_2D(dot(qP2HatCrossN[j], qSurfP1Grad[i]), q2Ddomain);
}

// TODO: Den Fall fullGrad testen!
void LocalSurfOseen::setupB_P1P2 (double B_P1P2[4][30])
{
    GridFunctionCL<double> qsurfgrad_k;
    qsurfgrad_k.resize(q2Ddomain.vertex_size());

    // Do all combinations for (i,j) i,j=4 x 30 and corresponding quadrature
    for (int k=0; k<3; ++k) {
        for (int i=0; i < 4; ++i) {
            ExtractComponent(qSurfP1Grad[i], qsurfgrad_k, k);
            for (int j=0; j<10; ++j) {
//                if (!fullGrad) {
                    B_P1P2[i][3*j+k]= quad_2D(qsurfgrad_k*qP2Hat[j], q2Ddomain);
//                }
//                else
//                    B_P1P2[i][3*j+k]= quad_2D(P1Grad[i][k]*qP2Hat[j], q2Ddomain);
            }
        }
    }
}

// TODO: Den Fall fullGrad testen!
void LocalSurfOseen::setupB_P2P2 (double B_P2P2[10][30])
{
    GridFunctionCL<double> qsurfgrad_k;
    qsurfgrad_k.resize(q2Ddomain.vertex_size());

    // Do all combinations for (i,j) i,j=4 x 30 and corresponding quadrature
    for (int k=0; k<3; ++k) {
        for (int i=0; i < 10; ++i) {
            ExtractComponent(qSurfP2Grad[i], qsurfgrad_k, k);
            for (int j=0; j<10; ++j) {
//                if (!fullGrad) {
                    B_P2P2[i][3*j+k]= quad_2D(qsurfgrad_k*qP2Hat[j], q2Ddomain);
//                }
//                else
//                    std::cout << "DOES NOT WORK YET!!!" << std::endl;//B_P2P2[i][3*j+k]= quad_2D(P2Grad[i][k]*qP2Hat[j], q2Ddomain);
            }
        }
    }
}

// TODO: Den Fall fullGrad testen!
void LocalSurfOseen::setupB_P1P1 (double B_P1P1[4][12])
{
    GridFunctionCL<double> qsurfgrad_k;
    qsurfgrad_k.resize(q2Ddomain.vertex_size());

    // Do all combinations for (i,j) i,j=4 x 30 and corresponding quadrature
    for (int k=0; k<3; ++k) {
        for (int i=0; i < 4; ++i) {
            ExtractComponent(qSurfP1Grad[i], qsurfgrad_k, k);
            for (int j=0; j<4; ++j) {
//                if (!fullGrad) {
                    B_P1P1[i][3*j+k]= quad_2D(qsurfgrad_k*qP1Hat[j], q2Ddomain);
//                }
//                else
//                    B_P1P1[i][3*j+k]= quad_2D(P1Grad[i][k]*qP1Hat[j], q2Ddomain);
            }
        }
    }
}

// TODO: Den Fall fullGrad testen!
void LocalSurfOseen::setupB_P2P1 (double B_P2P1[10][12]) {
    GridFunctionCL<double> qsurfgrad_k;
    qsurfgrad_k.resize(q2Ddomain.vertex_size());
    // Do all combinations for (i,j) i,j=4 x 30 and corresponding quadrature
    for (int k=0; k<3; ++k) {
        for (int i=0; i < 10; ++i) {
            ExtractComponent(qSurfP2Grad[i], qsurfgrad_k, k);
            for (int j=0; j<4; ++j)
                B_P2P1[i][3*j+k]= quad_2D(qsurfgrad_k*qP1Hat[j], q2Ddomain);
        }
    }
}

void LocalSurfOseen::setupM_P2 (double M_P2[10][10]) {
     for (int i=0; i < 10; ++i)
         for (int j=0; j<10; ++j)
             M_P2[i][j] = quad_2D(qP2Hat[i]*qP2Hat[j], q2Ddomain);
}

void LocalSurfOseen::setupF_P2(double F_P2[30]) {
    for (size_t i = 0; i < 30; ++i) {
        auto is = i / 3; // scalar shape index
        auto in = i - 3 * is; // nonzero vect component
        auto e_in = DROPS::std_basis<3>(in + 1);
        F_P2[i] = quad_2D(dot(e_in, qF) * qP2Hat[is], q2Ddomain);
    }
}

void LocalSurfOseen::setupN_P2(double N_P2[30][30]) {
    for (size_t i = 0; i < 30; ++i) {
        auto is = i / 3; // scalar shape index
        auto in = i - 3 * is; // nonzero vect component
        for (size_t j = 0; j < 30; ++j) {
            auto js = j / 3; // scalar shape index
            auto jn = j - 3 * js; // nonzero vect component
            N_P2[i][j] = quad_2D(dot(qWind, qP2Grad[is]) * qP2Hat[js] * take(qP, in, jn), q2Ddomain);
        }
    }
}

void LocalSurfOseen::setupM_P1 (double M_P1[4][4]) {
    // Do all combinations for (i,j) i,j=4 x 4 and corresponding quadrature
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
             M_P1[i][j]= quad_2D(qP1Hat[i] * qP1Hat[j], q2Ddomain);
}

void LocalSurfOseen::setupG_P1(double G_P1[4]) {
    for (size_t i = 0; i < 4; ++i)
        G_P1[i] = quad_2D(qG * qP1Hat[i], q2Ddomain);
}

void LocalSurfOseen::setupS_P2 (double S_P2[30][30])
{
    // Do all combinations for (i,j) i,j=30 x 30 and corresponding quadrature
    for (int i=0; i < 10; ++i) {
        for (int j=0; j<10; ++j) {
            for (int k=0; k<3; ++k) {
                for (int l=0; l<3; ++l) {
                     S_P2[3*i+k][3*j+l]= quad_2D(qP2Hat[i]*qP2NormalComp[k]*qP2Hat[j]*qP2NormalComp[l], q2Ddomain);
                }
            }
        }
    }
}

void LocalSurfOseen::setupS_P1 (double S_P1[12][12])
{
    // Do all combinations for (i,j) i,j=30 x 30 and corresponding quadrature
    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
            for (int k=0; k<3; ++k) {
                for (int l=0; l<3; ++l) {
                     S_P1[3*i+k][3*j+l]= quad_2D(qP1Hat[i]*qP2NormalComp[k]*qP1Hat[j]*qP2NormalComp[l], q2Ddomain);
                }
            }
        }
    }
}

void LocalSurfOseen::setupA_P2_stab (double A_P2_stab[10][10], double absdet) {
    for (int i=0; i<10; ++i) {
        for (int j=0; j<10; ++j) {
           A_P2_stab[i][j] = quad(dot(q3Dnormal, q3DP2Grad[i])*dot(q3Dnormal, q3DP2Grad[j]), absdet, q3Ddomain, AllTetraC);
        }
    }
}

void LocalSurfOseen::setupA_P1_stab (double A_P1_stab[4][4], double absdet) {
    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
           A_P1_stab[i][j] = quad(dot(q3Dnormal, q3DP1Grad[i])*dot(q3Dnormal, q3DP1Grad[j]), absdet, q3Ddomain, AllTetraC);
        }
    }
}

void LocalSurfOseen::setupM_P1_stab (double M_P1_stab[4][4], double absdet) {
    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
           M_P1_stab[i][j] = quad(dot(q3DP1Grad[i], q3DP1Grad[j]), absdet, q3Ddomain, AllTetraC);

        }
    }
}

/// \brief Basis Class for Accumulator to set up the matrices for interface Stokes.
class StokesIFAccumulator_P1P1CL : public TetraAccumulatorCL
{
  protected:
    const VecDescCL& lset;
    const LsetBndDataCL& lset_bnd;

    IdxDescCL &P1Idx_, &ScalarP1Idx_;
    //IdxT numP1[4], numScalarP1[4];

    MatrixCL &A_P1_, &A_P1_stab_,  &B_P1P1_, &M_P1_, &S_P1_, &L_P1P1_, &L_P1P1_stab_, &M_ScalarP1_, &A_ScalarP1_stab_;
    MatrixBuilderCL *mA_P1_, *mA_P1_stab_, *mB_P1P1_, *mM_P1_, *mS_P1_, *mL_P1P1_, *mL_P1P1_stab_, *mM_ScalarP1_, *mA_ScalarP1_stab_;
    double locA_P1[12][12], locA_P1_stab[4][4], locB_P1P1[4][12], locM_P1[4][4], locS_P1[12][12], locL_P1P1[4][12], locL_P1P1_stab[4][12];

    LocalSurfOseen localProblem;

    SMatrixCL<3,3> T;
    double det, absdet;
    LocalP2CL<> ls_loc;
    LocalNumbP1CL n, nScalar;

    ///\brief Computes the mapping from local to global data "n", the local matrices in loc and.
    void local_setup (const TetraCL& tet);
    ///\brief Update the global system->
    void update_global_system ();

public:
    StokesIFAccumulator_P1P1CL(const VecDescCL& ls, const LsetBndDataCL& ls_bnd, IdxDescCL& P1FE, IdxDescCL& ScalarP1FE, MatrixCL& A_P1, MatrixCL& A_P1_stab, MatrixCL& B_P1P1, MatrixCL& M_P1, MatrixCL& S_P1, MatrixCL& L_P1P1, MatrixCL& L_P1P1_stab, MatrixCL& M_ScalarP1, MatrixCL& A_ScalarP1_stab, SurfOseenParam* param)
        : lset(ls), lset_bnd(ls_bnd), P1Idx_(P1FE), ScalarP1Idx_(ScalarP1FE), A_P1_(A_P1), A_P1_stab_(A_P1_stab), B_P1P1_(B_P1P1), M_P1_(M_P1), S_P1_(S_P1), L_P1P1_(L_P1P1), L_P1P1_stab_(L_P1P1_stab), M_ScalarP1_(M_ScalarP1), A_ScalarP1_stab_(A_ScalarP1_stab), localProblem(param) {}
    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    void visit (const TetraCL& sit);

    TetraAccumulatorCL* clone (int /*tid*/) { return new StokesIFAccumulator_P1P1CL (*this); }

};

void StokesIFAccumulator_P1P1CL::begin_accumulation ()
{
    std::cout << "entering StokesIF: \n";
    const size_t num_unks_p1 = P1Idx_.NumUnknowns(), num_unks_scalarp1 = ScalarP1Idx_.NumUnknowns();
    mA_P1_= new MatrixBuilderCL(&A_P1_, num_unks_p1, num_unks_p1);
    mA_P1_stab_= new MatrixBuilderCL(&A_P1_stab_, num_unks_p1, num_unks_p1);
    mB_P1P1_= new MatrixBuilderCL(&B_P1P1_, num_unks_scalarp1, num_unks_p1);
    mM_P1_= new MatrixBuilderCL(&M_P1_, num_unks_p1, num_unks_p1);
    mS_P1_= new MatrixBuilderCL(&S_P1_, num_unks_p1, num_unks_p1);
    mL_P1P1_= new MatrixBuilderCL(&L_P1P1_, num_unks_scalarp1, num_unks_p1);
    mL_P1P1_stab_= new MatrixBuilderCL(&L_P1P1_stab_, num_unks_scalarp1, num_unks_p1);
    mM_ScalarP1_= new MatrixBuilderCL(&M_ScalarP1_, num_unks_scalarp1, num_unks_scalarp1);
    mA_ScalarP1_stab_= new MatrixBuilderCL(&A_ScalarP1_stab_, num_unks_scalarp1, num_unks_scalarp1);
}

void StokesIFAccumulator_P1P1CL::finalize_accumulation ()
{
    mB_P1P1_->Build();
    delete mB_P1P1_;
    mL_P1P1_->Build();
    delete mL_P1P1_;
    mL_P1P1_stab_->Build();
    delete mL_P1P1_stab_;
    mA_P1_->Build();
    delete mA_P1_;
    mM_P1_->Build();
    delete mM_P1_;
    mS_P1_->Build();
    delete mS_P1_;
    mA_P1_stab_->Build();
    delete mA_P1_stab_;
    mM_ScalarP1_->Build();
    delete mM_ScalarP1_;
    mA_ScalarP1_stab_->Build();
    delete mA_ScalarP1_stab_;
#ifndef _PAR
    std::cout << "StokesIF_P1P1:\t" << A_P1_.num_nonzeros() << " nonzeros in A, " << A_P1_stab_.num_nonzeros() << " nonzeros in A_stab, "
              << B_P1P1_.num_nonzeros() << " nonzeros in B, " << M_P1_.num_nonzeros() << " nonzeros in M, " << S_P1_.num_nonzeros() << " nonzeros in S, "
              << L_P1P1_.num_nonzeros() << " nonzeros in L, " << L_P1P1_stab_.num_nonzeros() << " nonzeros in L_stab, "
              << M_ScalarP1_.num_nonzeros() << " nonzeros in M_ScalarP1, " << A_ScalarP1_stab_.num_nonzeros() << " nonzeros in A_ScalarP1_stab\n";
#endif
    // std::cout << '\n';
}

void StokesIFAccumulator_P1P1CL::visit (const TetraCL& tet)
{
    ls_loc.assign(tet, lset, lset_bnd);

    if (isInCutMesh(ls_loc))
    {
        local_setup(tet);
        update_global_system();
    }

}

void StokesIFAccumulator_P1P1CL::local_setup (const TetraCL& tet)
{
    GetTrafoTr(T, det, tet);
    absdet= std::fabs(det);

    n.assign(tet, P1Idx_, P1Idx_.GetBndInfo());
    nScalar.assign(tet, ScalarP1Idx_, ScalarP1Idx_.GetBndInfo());

    //GetLocalNumbP1NoBnd(numP1, tet, P1Idx_);
    //GetLocalNumbP1NoBnd(numScalarP1, tet, ScalarP1Idx_);
    localProblem.calcIntegrands(T, ls_loc, tet);
    localProblem.calc3DIntegrands(T, ls_loc, tet);
    localProblem.setupA_P1(locA_P1);
    localProblem.setupA_P1_stab(locA_P1_stab, absdet);
    localProblem.setupB_P1P1(locB_P1P1);
    localProblem.setupM_P1(locM_P1);
    localProblem.setupS_P1(locS_P1);
}

void StokesIFAccumulator_P1P1CL::update_global_system ()
{
    MatrixBuilderCL& mA_P1 = *mA_P1_;
    MatrixBuilderCL& mA_P1_stab= *mA_P1_stab_;
    MatrixBuilderCL& mB_P1P1 = *mB_P1P1_;
    MatrixBuilderCL& mM_P1 = *mM_P1_;
    MatrixBuilderCL& mS_P1 = *mS_P1_;
    MatrixBuilderCL& mL_P1P1 = *mL_P1P1_;
    MatrixBuilderCL& mL_P1P1_stab= *mL_P1P1_stab_;
    MatrixBuilderCL& mM_ScalarP1 = *mM_ScalarP1_;
    MatrixBuilderCL& mA_ScalarP1_stab= *mA_ScalarP1_stab_;

    for(int i= 0; i < 4; ++i) {
        const IdxT ii= n.num[i];
        if (ii==NoIdx || !(n.WithUnknowns(i))) continue;
        for(int j=0; j < 4; ++j) {
            const IdxT jj= n.num[j];
            if (jj==NoIdx || !(n.WithUnknowns(j))) continue;
            for (int k=0; k<3; ++k) {
                for (int l=0; l<3; ++l) {
                    mA_P1(ii+k, jj+l) += locA_P1[3*j+l][3*i+k];
                    mS_P1(ii+k, jj+l) += locS_P1[3*j+l][3*i+k];
                    if(k == l) {
                        mM_P1(ii+k, jj+l) += locM_P1[j][i];
                        mA_P1_stab(ii+k, jj+l) += locA_P1_stab[j][i];
                    } else {
                        mM_P1(ii+k, jj+l) += 0.;
                        mA_P1_stab(ii+k, jj+l) += 0.;
                    }
                }
             }
        }
        for (int j=0; j < 4; ++j) {
            if (nScalar.num[j]==NoIdx || !(nScalar.WithUnknowns(j))) continue;
            for (int k=0; k<3; ++k) {
                mB_P1P1(nScalar.num[j], ii+k) += locB_P1P1[j][3*i+k];
                mL_P1P1(nScalar.num[j], ii+k) += locL_P1P1[j][3*i+k];
                mL_P1P1_stab(nScalar.num[j], ii+k) += locL_P1P1_stab[j][3*i+k];
            }
        }
    }
    for(int i= 0; i < 4; ++i) {
        const IdxT ii= nScalar.num[i];
        if (ii==NoIdx || !(nScalar.WithUnknowns(i))) continue;
        for(int j=0; j <4; ++j) {
            const IdxT jj= nScalar.num[j];
            if (jj==NoIdx || !(nScalar.WithUnknowns(j))) continue;
            mM_ScalarP1(ii, jj) += locM_P1[i][j];
            mA_ScalarP1_stab(ii, jj) += locA_P1_stab[i][j];
        }
    }
}

/// \brief Basis Class for Accumulator to set up the matrices for interface Navier-Stokes.
// TODO: derive from Stokes?
class NavierStokesIFAccumulator_P1P1CL : public TetraAccumulatorCL
{
  protected:
    const VecDescCL& lset;
    const VecDescCL& velocity;
    const LsetBndDataCL& lset_bnd;
    const BndDataCL<Point3DCL>& velocity_bnd;
    IdxDescCL &P1Idx_, &ScalarP1Idx_;
    //IdxT numP1[4], numScalarP1[4];

    MatrixCL &A_P1_, &A_P1_stab_,  &B_P1P1_, &Omega_P1P1_, &N_P1_, &NT_P1_,&M_P1_,&D_P1_, &S_P1_, &L_P1P1_, &L_P1P1_stab_, &M_ScalarP1_, &A_ScalarP1_stab_,  &Schur_normalP1_stab_;
    MatrixBuilderCL *mA_P1_, *mA_P1_stab_, *mB_P1P1_,*mOmega_P1P1_,*mN_P1_,  *mNT_P1_,*mM_P1_,*mD_P1_, *mS_P1_, *mL_P1P1_, *mL_P1P1_stab_, *mM_ScalarP1_, *mA_ScalarP1_stab_, *mSchur_normalP1_stab_;
    double locA_P1[12][12], locA_P1_stab[4][4], locB_P1P1[4][12],locOmega_P1P1[4][12],locN_P1[12][12], locNT_P1[12][12],locM_P1[4][4],locD_P1[4][4], locM_P1_stab[4][4], locS_P1[12][12], locL_P1P1[4][12], locL_P1P1_stab[4][12];

    LocalSurfOseen localProblem;

    SMatrixCL<3,3> T;
    double det, absdet;
    LocalP2CL<> ls_loc;
    LocalP1CL<Point3DCL> v_loc;
    LocalNumbP1CL n, nScalar;

    ///\brief Computes the mapping from local to global data "n", the local matrices in loc and.
    void local_setup (const TetraCL& tet);
    ///\brief Update the global system->
    void update_global_system ();

  public:
    NavierStokesIFAccumulator_P1P1CL(const LevelsetP2CL& ls, const VecDescCL& v, const BndDataCL<Point3DCL>& v_bnd, IdxDescCL& P1FE, IdxDescCL& ScalarP1FE, MatrixCL& A_P1, MatrixCL& A_P1_stab, MatrixCL& B_P1P1, MatrixCL& Omega_P1P1, MatrixCL& N_P1, MatrixCL& NT_P1, MatrixCL& M_P1, MatrixCL& D_P1, MatrixCL& S_P1, MatrixCL& L_P1P1, MatrixCL& L_P1P1_stab, MatrixCL& M_ScalarP1, MatrixCL& A_ScalarP1_stab, MatrixCL& Schur_normalP1_stab, SurfOseenParam* param)
        : lset(ls.Phi), velocity(v), lset_bnd(ls.GetBndData()), velocity_bnd(v_bnd), P1Idx_(P1FE), ScalarP1Idx_(ScalarP1FE), A_P1_(A_P1), A_P1_stab_(A_P1_stab), B_P1P1_(B_P1P1), Omega_P1P1_(Omega_P1P1), N_P1_(N_P1), NT_P1_(NT_P1), M_P1_(M_P1), D_P1_(D_P1), S_P1_(S_P1), L_P1P1_(L_P1P1), L_P1P1_stab_(L_P1P1_stab), M_ScalarP1_(M_ScalarP1), A_ScalarP1_stab_(A_ScalarP1_stab), Schur_normalP1_stab_(Schur_normalP1_stab), localProblem(param) {}
    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    void visit (const TetraCL& sit);

    TetraAccumulatorCL* clone (int /*tid*/) { return new NavierStokesIFAccumulator_P1P1CL (*this); }

};

void NavierStokesIFAccumulator_P1P1CL::begin_accumulation ()
{
    std::cout << "entering NavierStokesIF: \n";
    const size_t num_unks_p1 = P1Idx_.NumUnknowns(), num_unks_scalarp1 = ScalarP1Idx_.NumUnknowns();
    mA_P1_= new MatrixBuilderCL(&A_P1_, num_unks_p1, num_unks_p1);
    mN_P1_= new MatrixBuilderCL(&N_P1_, num_unks_p1, num_unks_p1);
    mNT_P1_= new MatrixBuilderCL(&NT_P1_, num_unks_p1, num_unks_p1);
    mA_P1_stab_= new MatrixBuilderCL(&A_P1_stab_, num_unks_p1, num_unks_p1);
    mB_P1P1_= new MatrixBuilderCL(&B_P1P1_, num_unks_scalarp1, num_unks_p1);
    mOmega_P1P1_= new MatrixBuilderCL(&Omega_P1P1_, num_unks_scalarp1, num_unks_p1);
    mM_P1_= new MatrixBuilderCL(&M_P1_, num_unks_p1, num_unks_p1);
    mD_P1_= new MatrixBuilderCL(&D_P1_, num_unks_p1, num_unks_p1);
    mS_P1_= new MatrixBuilderCL(&S_P1_, num_unks_p1, num_unks_p1);
    mL_P1P1_= new MatrixBuilderCL(&L_P1P1_, num_unks_scalarp1, num_unks_p1);
    mL_P1P1_stab_= new MatrixBuilderCL(&L_P1P1_stab_, num_unks_scalarp1, num_unks_p1);
    mM_ScalarP1_= new MatrixBuilderCL(&M_ScalarP1_, num_unks_scalarp1, num_unks_scalarp1);
    mA_ScalarP1_stab_= new MatrixBuilderCL(&A_ScalarP1_stab_, num_unks_scalarp1, num_unks_scalarp1);
    mSchur_normalP1_stab_= new MatrixBuilderCL(&Schur_normalP1_stab_, num_unks_scalarp1, num_unks_scalarp1);
}

void NavierStokesIFAccumulator_P1P1CL::finalize_accumulation ()
{

    mB_P1P1_->Build();
    delete mB_P1P1_;
    mOmega_P1P1_->Build();
        delete mOmega_P1P1_;
    mL_P1P1_->Build();
    delete mL_P1P1_;
    mL_P1P1_stab_->Build();
    delete mL_P1P1_stab_;
    mA_P1_->Build();
    delete mA_P1_;
    mN_P1_->Build();
    delete mN_P1_;
    mNT_P1_->Build();
    delete mNT_P1_;
    mM_P1_->Build();
    delete mM_P1_;

    mD_P1_->Build();
        delete mD_P1_;

    mS_P1_->Build();
    delete mS_P1_;
    mA_P1_stab_->Build();
    delete mA_P1_stab_;
    mM_ScalarP1_->Build();
    delete mM_ScalarP1_;
    mA_ScalarP1_stab_->Build();
    delete mA_ScalarP1_stab_;
    mSchur_normalP1_stab_->Build();
	delete mSchur_normalP1_stab_;
#ifndef _PAR
    std::cout << "NavierStokesIF_P1P1:\t" << A_P1_.num_nonzeros() << " nonzeros in A, " << A_P1_stab_.num_nonzeros() << " nonzeros in A_stab, " << N_P1_.num_nonzeros() << " nonzeros in N, "
    		<< NT_P1_.num_nonzeros() << " nonzeros in NT, "
    		<< B_P1P1_.num_nonzeros() << " nonzeros in B, " << Omega_P1P1_.num_nonzeros() << " nonzeros in Omega, " << M_P1_.num_nonzeros() << " nonzeros in M, "<< D_P1_.num_nonzeros() << " nonzeros in D, " << S_P1_.num_nonzeros() << " nonzeros in S, "
              << L_P1P1_.num_nonzeros() << " nonzeros in L, " << L_P1P1_stab_.num_nonzeros() << " nonzeros in L_stab, "
              << M_ScalarP1_.num_nonzeros() << " nonzeros in M_ScalarP1, " << A_ScalarP1_stab_.num_nonzeros() << " nonzeros in A_ScalarP1_stab\n";
#endif
    // std::cout << '\n';
}

void NavierStokesIFAccumulator_P1P1CL::visit (const TetraCL& tet)
{
    ls_loc.assign(tet, lset, lset_bnd);
    v_loc.assign(tet, velocity, velocity_bnd);

    if (isInCutMesh(ls_loc))
    {
        local_setup(tet);
        update_global_system();
    }

}

void NavierStokesIFAccumulator_P1P1CL::local_setup (const TetraCL& tet)
{
    GetTrafoTr(T, det, tet);
    absdet= std::fabs(det);

    n.assign(tet, P1Idx_, P1Idx_.GetBndInfo());
    nScalar.assign(tet, ScalarP1Idx_, ScalarP1Idx_.GetBndInfo());

    // GetLocalNumbP1NoBnd(numP1, tet, P1Idx_);
    // GetLocalNumbP1NoBnd(numScalarP1, tet, ScalarP1Idx_);
    // localProblem.exportPatchInfo("../../../MKL-Eigs-for-Sparse-Matrices/output/patch/hOver" + std::to_string(localProblem.num_intervals()), T, ls_loc, tet);

    localProblem.calcIntegrands(T, ls_loc, v_loc, tet);
    localProblem.calc3DIntegrands(T, ls_loc, tet);
    localProblem.setupA_P1(locA_P1);
    localProblem.setupA_P1_stab(locA_P1_stab, absdet);
    localProblem.setupB_P1P1(locB_P1P1);
    localProblem.setupM_P1(locM_P1);
    localProblem.setupM_P1_stab(locM_P1_stab, absdet);
    localProblem.setupS_P1(locS_P1);

}

void NavierStokesIFAccumulator_P1P1CL::update_global_system ()
{
    MatrixBuilderCL& mA_P1 = *mA_P1_;
    MatrixBuilderCL& mA_P1_stab= *mA_P1_stab_;
    MatrixBuilderCL& mB_P1P1 = *mB_P1P1_;
    MatrixBuilderCL& mOmega_P1P1 = *mOmega_P1P1_;
    MatrixBuilderCL& mN_P1 = *mN_P1_;
    MatrixBuilderCL& mNT_P1 = *mNT_P1_;
    MatrixBuilderCL& mM_P1 = *mM_P1_;
    MatrixBuilderCL& mD_P1 = *mD_P1_;
    MatrixBuilderCL& mS_P1 = *mS_P1_;
    MatrixBuilderCL& mL_P1P1 = *mL_P1P1_;
    MatrixBuilderCL& mL_P1P1_stab= *mL_P1P1_stab_;
    MatrixBuilderCL& mM_ScalarP1 = *mM_ScalarP1_;
    MatrixBuilderCL& mA_ScalarP1_stab= *mA_ScalarP1_stab_;
    MatrixBuilderCL& mSchur_normalP1_stab= *mSchur_normalP1_stab_;

int count =0;
    for(int i= 0; i < 4; ++i) {
        const IdxT ii= n.num[i];
        if (ii==NoIdx || !(n.WithUnknowns(i))) continue;
        for(int j=0; j < 4; ++j) {
            const IdxT jj= n.num[j];
            if (jj==NoIdx || !(n.WithUnknowns(j))) continue;
            for (int k=0; k<3; ++k) {
                for (int l=0; l<3; ++l) {
                    mA_P1(ii+k, jj+l) += locA_P1[3*j+l][3*i+k];
                    mN_P1(ii+k, jj+l) += locN_P1[3*j+l][3*i+k];
                    mNT_P1(ii+k, jj+l) += locNT_P1[3*j+l][3*i+k];
                    mS_P1(ii+k, jj+l) += locS_P1[3*j+l][3*i+k];

                    if(k == l) {
                    	mM_P1(ii+k, jj+l) += locM_P1[j][i];

                    	mD_P1(ii+k, jj+l) += locD_P1[j][i];

                    	count++;
                    	mA_P1_stab(ii+k, jj+l) += locA_P1_stab[j][i];
                    } else {
                    	mM_P1(ii+k, jj+l) += 0.;
                    	mD_P1(ii+k, jj+l) += 0.;
                    	mA_P1_stab(ii+k, jj+l) += 0.;
                    }

                }
             }
        }
        for (int j=0; j < 4; ++j) {
            if (nScalar.num[j]==NoIdx || !(nScalar.WithUnknowns(j))) continue;
            for (int k=0; k<3; ++k) {
                mB_P1P1(nScalar.num[j], ii+k) += locB_P1P1[j][3*i+k];
                mOmega_P1P1(nScalar.num[j], ii+k) += locOmega_P1P1[j][3*i+k];
                mL_P1P1(nScalar.num[j], ii+k) += locL_P1P1[j][3*i+k];
                mL_P1P1_stab(nScalar.num[j], ii+k) += locL_P1P1_stab[j][3*i+k];
            }
        }
    }
    for(int i= 0; i < 4; ++i) {
        const IdxT ii= nScalar.num[i];
        if (ii==NoIdx || !(nScalar.WithUnknowns(i))) continue;
        for(int j=0; j <4; ++j) {
            const IdxT jj= nScalar.num[j];
            if (jj==NoIdx || !(nScalar.WithUnknowns(j))) continue;
            mM_ScalarP1(ii, jj) += locM_P1[i][j];
            mA_ScalarP1_stab(ii, jj) += locM_P1_stab[i][j];
            mSchur_normalP1_stab(ii, jj) += locA_P1_stab[i][j];
        }
    }
}

void SetupStokesIF_P1P1(const MultiGridCL& MG_, MatDescCL* A_P1, MatDescCL* A_P1_stab, MatDescCL* B_P1P1, MatDescCL* M_P1, MatDescCL* S_P1, MatDescCL* L_P1P1, MatDescCL* L_P1P1_stab, MatDescCL* M_ScalarP1, MatDescCL* A_ScalarP1_stab, const VecDescCL& lset, const LsetBndDataCL& lset_bnd, SurfOseenParam* param)
{
  ScopeTimerCL scope("SetupStokesIF_P1P1");
  StokesIFAccumulator_P1P1CL accu(lset, lset_bnd, *(A_P1->RowIdx), *(L_P1P1->RowIdx), A_P1->Data, A_P1_stab->Data, B_P1P1->Data, M_P1->Data, S_P1->Data, L_P1P1->Data, L_P1P1_stab->Data, M_ScalarP1->Data, A_ScalarP1_stab->Data, param);
  TetraAccumulatorTupleCL accus;
  //    MaybeAddProgressBar(MG_, "LapBeltr(P2) Setup", accus, RowIdx.TriangLevel());
  accus.push_back(&accu);
  accumulate(accus, MG_, A_P1->GetRowLevel(), A_P1->RowIdx->GetBndInfo());
}

void SetupNavierStokesIF_P1P1(
        const MultiGridCL& MG_,
        MatDescCL* A_P1,
        MatDescCL* A_P1_stab,
        MatDescCL* B_P1P1,
        MatDescCL* Omega_P1P1,
        MatDescCL* N_P1,
        MatDescCL* NT_P1,
        MatDescCL* M_P1,
        MatDescCL* D_P1,
        MatDescCL* S_P1,
        MatDescCL* L_P1P1,
        MatDescCL* L_P1P1_stab,
        MatDescCL* M_ScalarP1,
        MatDescCL* A_ScalarP1_stab,
        MatDescCL* Schur_normalP1_stab, const LevelsetP2CL& lset, const VecDescCL& velocity, const BndDataCL<Point3DCL>& velocity_bnd, SurfOseenParam* param)
{
  ScopeTimerCL scope("SetupNavierStokesIF_P1P1");
  NavierStokesIFAccumulator_P1P1CL accu(lset, velocity, velocity_bnd, *(A_P1->RowIdx), *(L_P1P1->RowIdx), A_P1->Data, A_P1_stab->Data, B_P1P1->Data, Omega_P1P1->Data,  N_P1->Data, NT_P1->Data, M_P1->Data, D_P1->Data, S_P1->Data, L_P1P1->Data, L_P1P1_stab->Data, M_ScalarP1->Data, A_ScalarP1_stab->Data, Schur_normalP1_stab->Data, param);
  TetraAccumulatorTupleCL accus;
  //    MaybeAddProgressBar(MG_, "LapBeltr(P2) Setup", accus, RowIdx.TriangLevel());
  accus.push_back(&accu);
  accumulate(accus, MG_, A_P1->GetRowLevel(), /*A_P1->RowIdx->GetMatchingFunction(),*/ A_P1->RowIdx->GetBndInfo());
}

/// \brief Basis Class for Accumulator to set up the matrices for interface Stokes.
class StokesIFAccumulator_P1P2CL : public TetraAccumulatorCL
{
  protected:
    const VecDescCL& lset;
    const LsetBndDataCL& lset_bnd;

    IdxDescCL &P1Idx_, &ScalarP2Idx_;
    //IdxT numP1[4], numScalarP1[4];

    MatrixCL &A_P1_, &A_P1_stab_,  &B_P2P1_, &M_P1_, &S_P1_, &L_P2P1_, &L_P2P1_stab_, &M_ScalarP2_, &A_ScalarP2_stab_;
    MatrixBuilderCL *mA_P1_, *mA_P1_stab_, *mB_P2P1_, *mM_P1_, *mS_P1_, *mL_P2P1_, *mL_P2P1_stab_, *mM_ScalarP2_, *mA_ScalarP2_stab_;
    double locA_P1[12][12], locA_P1_stab[4][4], locA_P2_stab[10][10], locB_P2P1[10][12], locM_P1[4][4], locM_P2[10][10], locS_P1[12][12], locL_P2P1[10][12], locL_P2P1_stab[10][12];

    LocalSurfOseen localProblem;

    SMatrixCL<3,3> T;
    double det, absdet;
    LocalP2CL<> ls_loc;
    LocalNumbP1CL n;
    LocalNumbP2CL nScalar;

    ///\brief Computes the mapping from local to global data "n", the local matrices in loc and.
    void local_setup (const TetraCL& tet);
    ///\brief Update the global system->
    void update_global_system ();
public:
    StokesIFAccumulator_P1P2CL(const VecDescCL& ls, const LsetBndDataCL& ls_bnd, IdxDescCL& P1FE, IdxDescCL& ScalarP2FE, MatrixCL& A_P1, MatrixCL& A_P1_stab, MatrixCL& B_P2P1, MatrixCL& M_P1, MatrixCL& S_P1, MatrixCL& L_P2P1, MatrixCL& L_P2P1_stab, MatrixCL& M_ScalarP2, MatrixCL& A_ScalarP2_stab, SurfOseenParam* param)
        : lset(ls), lset_bnd(ls_bnd), P1Idx_(P1FE), ScalarP2Idx_(ScalarP2FE), A_P1_(A_P1), A_P1_stab_(A_P1_stab), B_P2P1_(B_P2P1), M_P1_(M_P1), S_P1_(S_P1), L_P2P1_(L_P2P1), L_P2P1_stab_(L_P2P1_stab), M_ScalarP2_(M_ScalarP2), A_ScalarP2_stab_(A_ScalarP2_stab), localProblem(param) {}
    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    void visit (const TetraCL& sit);

    TetraAccumulatorCL* clone (int /*tid*/) { return new StokesIFAccumulator_P1P2CL (*this); }

};

void StokesIFAccumulator_P1P2CL::begin_accumulation ()
{
    std::cout << "entering StokesIF: \n";
    const size_t num_unks_p1 = P1Idx_.NumUnknowns(), num_unks_scalarp2 = ScalarP2Idx_.NumUnknowns();
    mA_P1_= new MatrixBuilderCL(&A_P1_, num_unks_p1, num_unks_p1);
    mA_P1_stab_= new MatrixBuilderCL(&A_P1_stab_, num_unks_p1, num_unks_p1);
    mB_P2P1_= new MatrixBuilderCL(&B_P2P1_, num_unks_scalarp2, num_unks_p1);
    mM_P1_= new MatrixBuilderCL(&M_P1_, num_unks_p1, num_unks_p1);
    mS_P1_= new MatrixBuilderCL(&S_P1_, num_unks_p1, num_unks_p1);
    mL_P2P1_= new MatrixBuilderCL(&L_P2P1_, num_unks_scalarp2, num_unks_p1);
    mL_P2P1_stab_= new MatrixBuilderCL(&L_P2P1_stab_, num_unks_scalarp2, num_unks_p1);
    mM_ScalarP2_= new MatrixBuilderCL(&M_ScalarP2_, num_unks_scalarp2, num_unks_scalarp2);
    mA_ScalarP2_stab_= new MatrixBuilderCL(&A_ScalarP2_stab_, num_unks_scalarp2, num_unks_scalarp2);
}

void StokesIFAccumulator_P1P2CL::finalize_accumulation ()
{
    mB_P2P1_->Build();
    delete mB_P2P1_;
    mL_P2P1_->Build();
    delete mL_P2P1_;
    mL_P2P1_stab_->Build();
    delete mL_P2P1_stab_;
    mA_P1_->Build();
    delete mA_P1_;
    mM_P1_->Build();
    delete mM_P1_;
    mS_P1_->Build();
    delete mS_P1_;
    mA_P1_stab_->Build();
    delete mA_P1_stab_;
    mM_ScalarP2_->Build();
    delete mM_ScalarP2_;
    mA_ScalarP2_stab_->Build();
    delete mA_ScalarP2_stab_;
#ifndef _PAR
    std::cout << "StokesIF_P1P2:\t" << A_P1_.num_nonzeros() << " nonzeros in A, " << A_P1_stab_.num_nonzeros() << " nonzeros in A_stab, "
              << B_P2P1_.num_nonzeros() << " nonzeros in B, " << M_P1_.num_nonzeros() << " nonzeros in M, " << S_P1_.num_nonzeros() << " nonzeros in S, "
              << L_P2P1_.num_nonzeros() << " nonzeros in L, " << L_P2P1_stab_.num_nonzeros() << " nonzeros in L_stab, "
              << M_ScalarP2_.num_nonzeros() << " nonzeros in M_ScalarP1, " << A_ScalarP2_stab_.num_nonzeros() << " nonzeros in A_ScalarP1_stab\n";
#endif
    // std::cout << '\n';
}

void StokesIFAccumulator_P1P2CL::visit (const TetraCL& tet)
{
    ls_loc.assign(tet, lset, lset_bnd);

    if (isInCutMesh(ls_loc))
    {
        local_setup(tet);
        update_global_system();
    }

}

void StokesIFAccumulator_P1P2CL::local_setup (const TetraCL& tet)
{
    GetTrafoTr(T, det, tet);
    absdet= std::fabs(det);

    n.assign(tet, P1Idx_, P1Idx_.GetBndInfo());
    nScalar.assign(tet, ScalarP2Idx_, ScalarP2Idx_.GetBndInfo());

    //GetLocalNumbP1NoBnd(numP1, tet, P1Idx_);
    //GetLocalNumbP1NoBnd(numScalarP1, tet, ScalarP1Idx_);
    localProblem.calcIntegrands(T, ls_loc, tet);
    localProblem.calc3DIntegrands(T, ls_loc, tet);
    localProblem.setupA_P1(locA_P1);
    localProblem.setupA_P1_stab(locA_P1_stab, absdet);
    localProblem.setupA_P2_stab(locA_P2_stab, absdet);
    localProblem.setupB_P2P1(locB_P2P1);
    localProblem.setupM_P1(locM_P1);
    localProblem.setupM_P2(locM_P2);
    localProblem.setupS_P1(locS_P1);
}

void StokesIFAccumulator_P1P2CL::update_global_system ()
{
    MatrixBuilderCL& mA_P1 = *mA_P1_;
    MatrixBuilderCL& mA_P1_stab= *mA_P1_stab_;
    MatrixBuilderCL& mB_P2P1 = *mB_P2P1_;
    MatrixBuilderCL& mM_P1 = *mM_P1_;
    MatrixBuilderCL& mS_P1 = *mS_P1_;
    MatrixBuilderCL& mL_P2P1 = *mL_P2P1_;
    MatrixBuilderCL& mL_P2P1_stab= *mL_P2P1_stab_;
    MatrixBuilderCL& mM_ScalarP2 = *mM_ScalarP2_;
    MatrixBuilderCL& mA_ScalarP2_stab= *mA_ScalarP2_stab_;

    for(int i= 0; i < 4; ++i) {
        const IdxT ii= n.num[i];
        if (ii==NoIdx || !(n.WithUnknowns(i))) continue;
        for(int j=0; j < 4; ++j) {
            const IdxT jj= n.num[j];
            if (jj==NoIdx || !(n.WithUnknowns(j))) continue;
            for (int k=0; k<3; ++k) {
                for (int l=0; l<3; ++l) {
                    mA_P1(ii+k, jj+l) += locA_P1[3*j+l][3*i+k];
                    mS_P1(ii+k, jj+l) += locS_P1[3*j+l][3*i+k];
                    if(k == l) {
                        mM_P1(ii+k, jj+l) += locM_P1[j][i];
                        mA_P1_stab(ii+k, jj+l) += locA_P1_stab[j][i];
                    } else {
                        mM_P1(ii+k, jj+l) += 0.;
                        mA_P1_stab(ii+k, jj+l) += 0.;
                    }
                }
             }
        }
        for (int j=0; j < 10; ++j) {
            if (nScalar.num[j]==NoIdx || !(nScalar.WithUnknowns(j))) continue;
            for (int k=0; k<3; ++k) {
                mB_P2P1(nScalar.num[j], ii+k) += locB_P2P1[j][3*i+k];
                mL_P2P1(nScalar.num[j], ii+k) += locL_P2P1[j][3*i+k];
                mL_P2P1_stab(nScalar.num[j], ii+k) += locL_P2P1_stab[j][3*i+k];
            }
        }
    }
    for(int i= 0; i < 10; ++i) {
        const IdxT ii= nScalar.num[i];
        if (ii==NoIdx || !(nScalar.WithUnknowns(i))) continue;
        for(int j=0; j < 10; ++j) {
            const IdxT jj= nScalar.num[j];
            if (jj==NoIdx || !(nScalar.WithUnknowns(j))) continue;
            mM_ScalarP2(ii, jj) += locM_P2[i][j];
            mA_ScalarP2_stab(ii, jj) += locA_P2_stab[i][j];
        }
    }
}

void SetupStokesIF_P1P2(const MultiGridCL& MG_, MatDescCL* A_P1, MatDescCL* A_P1_stab, MatDescCL* B_P2P1, MatDescCL* M_P1, MatDescCL* S_P1, MatDescCL* L_P2P1, MatDescCL* L_P2P1_stab, MatDescCL* M_ScalarP2, MatDescCL* A_ScalarP2_stab, const VecDescCL& lset, const LsetBndDataCL& lset_bnd, SurfOseenParam* param)
{
  ScopeTimerCL scope("SetupStokesIF_P1P2");
  StokesIFAccumulator_P1P2CL accu(lset, lset_bnd, *(A_P1->RowIdx), *(L_P2P1->RowIdx), A_P1->Data, A_P1_stab->Data, B_P2P1->Data, M_P1->Data, S_P1->Data, L_P2P1->Data, L_P2P1_stab->Data, M_ScalarP2->Data, A_ScalarP2_stab->Data, param);
  TetraAccumulatorTupleCL accus;
  //    MaybeAddProgressBar(MG_, "LapBeltr(P2) Setup", accus, RowIdx.TriangLevel());
  accus.push_back(&accu);
  accumulate(accus, MG_, A_P1->GetRowLevel(), A_P1->RowIdx->GetBndInfo());
}

/// \brief Accumulator to set up the matrices for interface Oseen problem
class SurfOseenAccumulatorP2P1 : public TetraAccumulatorCL {
  private:
    const VecDescCL& lset;
    IdxDescCL& P2Idx_;
    IdxDescCL& ScalarP1Idx_;
    IdxT numP2[10], numScalarP1[4];
    SurfOseenSystem* system;
    using LocalVelocitySetupFunction = void(LocalSurfOseen::*)(double[30][30]);
    LocalVelocitySetupFunction setupA_P2;
    using LocalPressureSetupFunction = void(LocalSurfOseen::*)(double[4][4], double);
    LocalPressureSetupFunction setupC_P1;
    MatrixBuilderCL *mN_P2_, *mAL_P2_, *mA_P2_, *mA_P2_stab_, *mB_P1P2_, *mQ_P1P2_, *mM_P2_, *mS_P2_, *mM_P1_, *mC_P1_, *mA_P1_;
    // *mLB_P2_, *mLB_stab_P2_;
    double locAL_P2[30][30], locN_P2[30][30], locA_P2[30][30], locA_P2_stab[10][10], locB_P1P2[4][30], locQ_P1P2[4][30], locM_P2[10][10], locS_P2[30][30], locM_P1[4][4], locC_P1[4][4], locA_P1[4][4];
    // locLB_P2[10][10];
    double locRhsAL_P2[30], locF_P2[30], locG_P1[4];
    LocalSurfOseen localProblem;
    SMatrixCL<3,3> T;
    double det, absdet;
    LocalP2CL<> ls_loc;
    LocalP2CL<Point3DCL> w_loc;
    ///\brief Computes the mapping from local to global data "n", the local matrices in loc and.
    void local_setup (const TetraCL& tet);
    ///\brief Update the global system->
    void update_global_system ();

public:
    SurfOseenAccumulatorP2P1(const LevelsetP2CL& ls, SurfOseenSystem* system, SurfOseenParam* param)
    : lset(ls.Phi)
    , P2Idx_(*(system->A.RowIdx))
    , ScalarP1Idx_(*(system->B.RowIdx))
    , system(system)
    , localProblem(param) {
        setupA_P2 = param->input.formulation == SurfOseenParam::Formulation::consistent ? &LocalSurfOseen::setupA_P2_consistent : &LocalSurfOseen::setupA_P2_inconsistent;
        setupC_P1 = param->input.stab == SurfOseenParam::PressureVolumestab::normal ? &LocalSurfOseen::setupA_P1_stab : &LocalSurfOseen::setupM_P1_stab;
    }
    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();
    void visit (const TetraCL& sit);
    TetraAccumulatorCL* clone (int /*tid*/) { return new SurfOseenAccumulatorP2P1 (*this); }
};

void SurfOseenAccumulatorP2P1::begin_accumulation() {
    std::cout << "entering StokesIF: \n";
    size_t num_unks_p1_scalar = ScalarP1Idx_.NumUnknowns();
    size_t num_unks_p2 = P2Idx_.NumUnknowns();
    mA_P2_= new MatrixBuilderCL(&system->A.Data, num_unks_p2, num_unks_p2);
    mAL_P2_ = new MatrixBuilderCL(&system->AL.Data, num_unks_p2, num_unks_p2);
    mN_P2_= new MatrixBuilderCL(&system->N.Data, num_unks_p2, num_unks_p2);
    mA_P2_stab_= new MatrixBuilderCL(&system->A_stab.Data, num_unks_p2, num_unks_p2);
    mB_P1P2_= new MatrixBuilderCL(&system->B.Data, num_unks_p1_scalar, num_unks_p2);
    mQ_P1P2_= new MatrixBuilderCL(&system->Q.Data, num_unks_p1_scalar, num_unks_p2);
    mM_P2_= new MatrixBuilderCL(&system->M.Data, num_unks_p2, num_unks_p2);
    mS_P2_= new MatrixBuilderCL(&system->S.Data, num_unks_p2, num_unks_p2);
    mM_P1_= new MatrixBuilderCL(&system->M_p.Data, num_unks_p1_scalar, num_unks_p1_scalar);
    mC_P1_ = new MatrixBuilderCL(&system->C.Data, num_unks_p1_scalar, num_unks_p1_scalar);
    mA_P1_ = new MatrixBuilderCL(&system->A_p.Data, num_unks_p1_scalar, num_unks_p1_scalar);
    // size_t num_unks_p2_scalar = num_unks_p2 / 3;
    // mLB_P2_ = new MatrixBuilderCL(&system->LB.Data, num_unks_p2_scalar, num_unks_p2_scalar);
    // mLB_stab_P2_ = new MatrixBuilderCL(&system->LB_stab.Data, num_unks_p2_scalar, num_unks_p2_scalar);
    system->fRHS.Data.resize(num_unks_p2, 0.);
    system->alRHS.Data.resize(num_unks_p2, 0.);
    system->gRHS.Data.resize(num_unks_p1_scalar, 0.);
}

void SurfOseenAccumulatorP2P1::finalize_accumulation() {
    mA_P2_->Build();
    delete mA_P2_;
    mAL_P2_->Build();
    delete mAL_P2_;
    mN_P2_->Build();
    delete mN_P2_;
    mA_P2_stab_->Build();
    delete mA_P2_stab_;
    mB_P1P2_->Build();
    delete mB_P1P2_;
    mQ_P1P2_->Build();
    delete mQ_P1P2_;
    mM_P2_->Build();
    delete mM_P2_;
    mS_P2_->Build();
    delete mS_P2_;
    mM_P1_->Build();
    delete mM_P1_;
    mC_P1_->Build();
    delete mC_P1_;
    mA_P1_->Build();
    delete mA_P1_;
    // mLB_P2_->Build();
    // delete mLB_P2_;
    // mLB_stab_P2_->Build();
    // delete mLB_stab_P2_;
}

void SurfOseenAccumulatorP2P1::visit(const TetraCL& tet) {
    ls_loc.assign(tet, lset, BndDataCL<>());
    if (isInCutMesh(ls_loc)) {
        w_loc.assign(tet, system->w, BndDataCL<Point3DCL>());
        local_setup(tet);
        update_global_system();
    }
}

void SurfOseenAccumulatorP2P1::local_setup (const TetraCL& tet) {
    GetTrafoTr(T, det, tet);
    absdet = std::fabs(det);
    GetLocalNumbP2NoBnd(numP2, tet, P2Idx_);
    GetLocalNumbP1NoBnd(numScalarP1, tet, ScalarP1Idx_);
    // localProblem.exportPatchInfo("../../../MKL-Eigs-for-Sparse-Matrices/output/patch/hOver" + std::to_string(localProblem.num_intervals()), T, ls_loc, tet);
    localProblem.calcIntegrands(T, ls_loc, w_loc, tet);
    localProblem.calc3DIntegrands(T, ls_loc, tet);
    (localProblem.*setupA_P2)(locA_P2);
    localProblem.setupAL_P2(locAL_P2);
    localProblem.setupN_P2(locN_P2);
    localProblem.setupA_P2_stab(locA_P2_stab, absdet);
    localProblem.setupB_P1P2(locB_P1P2);
    localProblem.setupQ_P1P2(locQ_P1P2);
    localProblem.setupM_P2(locM_P2);
    localProblem.setupS_P2(locS_P2);
    localProblem.setupM_P1(locM_P1);
    (localProblem.*setupC_P1)(locC_P1, absdet);
    localProblem.setupLB_P1(locA_P1);
    // localProblem.setupLB_P2(locLB_P2);
    localProblem.setupF_P2(locF_P2);
    localProblem.setupG_P1(locG_P1);
    localProblem.setupRhsAL_P2(locRhsAL_P2);
}

void SurfOseenAccumulatorP2P1::update_global_system() {
    auto& mA_P2 = *mA_P2_;
    auto& mAL_P2 = *mAL_P2_;
    auto& mN_P2 = *mN_P2_;
    auto& mA_P2_stab = *mA_P2_stab_;
    auto& mB_P1P2 = *mB_P1P2_;
    auto& mQ_P1P2 = *mQ_P1P2_;
    auto& mM_P2 = *mM_P2_;
    auto& mS_P2 = *mS_P2_;
    auto& mM_P1 = *mM_P1_;
    auto& mC_P1  = *mC_P1_;
    auto& mA_P1  = *mA_P1_;
    // auto& mLB_P2 = *mLB_P2_;
    // auto& mLB_P2_stab = *mLB_stab_P2_;
    for(int i = 0; i < 10; ++i) {
        const IdxT ii= numP2[i];
        if (ii==NoIdx) continue;
        // assemble rhs
        for (int k = 0; k < 3; ++k) {
            system->fRHS.Data[ii + k] += locF_P2[3 * i + k];
            system->alRHS.Data[ii + k] += locRhsAL_P2[3 * i + k];
        }
        // assemble mtx
        for(int j = 0; j < 10; ++j) {
            const IdxT jj= numP2[j];
            if (jj==NoIdx) continue;
            // mLB_P2(ii/3, jj/3) += locLB_P2[j][i];
            // mLB_P2_stab(ii/3, jj/3) += locA_P2_stab[j][i];
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l) {
                    mA_P2(ii+k, jj+l) += locA_P2[3*j+l][3*i+k];
                    mAL_P2(ii+k, jj+l) += locAL_P2[3*j+l][3*i+k];
                    mN_P2(ii+k, jj+l) += locN_P2[3*j+l][3*i+k];
                    mS_P2(ii+k, jj+l) += locS_P2[3*j+l][3*i+k];
                    if(k == l) {
                        mM_P2(ii+k, jj+l) += locM_P2[j][i];
                        mA_P2_stab(ii+k, jj+l) += locA_P2_stab[j][i];
                    } else {
                        mM_P2(ii+k, jj+l) += 0.;
                        mA_P2_stab(ii+k, jj+l) += 0.;
                    }
                }
        }
        for (int j=0; j < 4; ++j) {
            if (numScalarP1[j]==NoIdx) continue;
            for (int k=0; k<3; ++k) {
                mB_P1P2(numScalarP1[j], ii+k) += locB_P1P2[j][3*i+k];
                mQ_P1P2(numScalarP1[j], ii+k) += locQ_P1P2[j][3*i+k];
            }
        }
    }
    for(int i = 0; i < 4; ++i) {
        const IdxT ii= numScalarP1[i];
        if (ii==NoIdx) continue;
        // assemble rhs
        system->gRHS.Data[ii] += locG_P1[i];
        // assemble mtx
        for(int j=0; j <4; ++j) {
            const IdxT jj= numScalarP1[j];
            if (jj==NoIdx) continue;
            mM_P1(ii, jj) += locM_P1[i][j];
            mC_P1(ii, jj) += locC_P1[i][j];
            mA_P1(ii, jj) += locA_P1[i][j];
        }
    }
}

void SetupSurfOseen_P2P1(const MultiGridCL& MG_, const LevelsetP2CL& lset, SurfOseenSystem* system, SurfOseenParam* param) {
    ScopeTimerCL scope("SetupSurfOseen_P2P1");
    SurfOseenAccumulatorP2P1 accu(lset, system, param);
    TetraAccumulatorTupleCL accus;
    accus.push_back(&accu);
    accumulate(accus, MG_, system->A.GetRowLevel(), system->A.RowIdx->GetBndInfo());
}

/// \brief Accumulator to set up the matrices for interface Stokes.
class StokesIFAccumulator_P2P2CL : public TetraAccumulatorCL
{
  private:
    const VecDescCL& lset;
    const LsetBndDataCL& lset_bnd;

    IdxDescCL& P2Idx_;
    IdxDescCL& ScalarP2Idx_;
    IdxT numP2[10], numScalarP2[10];

    MatrixCL &A_P2_, &A_P2_stab_, &B_P2P2_, &M_P2_, &S_P2_, &L_P2P2_, &L_P2P2_stab_, &M_ScalarP2_, &A_ScalarP2_stab_;
    MatrixBuilderCL *mA_P2_, *mA_P2_stab_, *mB_P2P2_, *mM_P2_, *mS_P2_, *mL_P2P2_, *mL_P2P2_stab_, *mM_ScalarP2_, *mA_ScalarP2_stab_ ;

    double locA_P2[30][30], locA_P2_stab[10][10], locB_P2P2[10][30], locM_P2[10][10], locS_P2[30][30], locL_P2P2[10][30], locL_P2P2_stab[10][30], locM_ScalarP2[10][10], locA_ScalarP2_stab[10][10];

    LocalSurfOseen localProblem;

    SMatrixCL<3,3> T;
    double det, absdet;
    LocalP2CL<> ls_loc;

    ///\brief Computes the mapping from local to global data "n", the local matrices in loc and.
    void local_setup (const TetraCL& tet);
    ///\brief Update the global system->
    void update_global_system ();

public:
    StokesIFAccumulator_P2P2CL(const VecDescCL& ls, const LsetBndDataCL& ls_bnd, IdxDescCL& P2FE, IdxDescCL& ScalarP2FE, MatrixCL& A_P2, MatrixCL& A_P2_stab, MatrixCL& B_P2P2, MatrixCL& M_P2, MatrixCL& S_P2, MatrixCL& L_P2P2, MatrixCL& L_P2P2_stab, MatrixCL& M_ScalarP2, MatrixCL& A_ScalarP2_stab, SurfOseenParam* param)
        : lset(ls), lset_bnd(ls_bnd), P2Idx_(P2FE), ScalarP2Idx_(ScalarP2FE), A_P2_(A_P2), A_P2_stab_(A_P2_stab), B_P2P2_(B_P2P2), M_P2_(M_P2), S_P2_(S_P2), L_P2P2_(L_P2P2), L_P2P2_stab_(L_P2P2_stab), M_ScalarP2_(M_ScalarP2), A_ScalarP2_stab_(A_ScalarP2_stab), localProblem(param) {}

    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    void visit (const TetraCL& sit);

    TetraAccumulatorCL* clone (int /*tid*/) { return new StokesIFAccumulator_P2P2CL (*this); }
};

void StokesIFAccumulator_P2P2CL::begin_accumulation ()
{
    std::cout << "entering StokesIF: \n";
    const size_t num_unks_scalarp2 = ScalarP2Idx_.NumUnknowns(), num_unks_p2 = P2Idx_.NumUnknowns();
    mA_P2_= new MatrixBuilderCL(&A_P2_, num_unks_p2, num_unks_p2);
    mA_P2_stab_= new MatrixBuilderCL(&A_P2_stab_, num_unks_p2, num_unks_p2);
    mB_P2P2_= new MatrixBuilderCL(&B_P2P2_, num_unks_scalarp2, num_unks_p2);
    mM_P2_= new MatrixBuilderCL(&M_P2_, num_unks_p2, num_unks_p2);
    mS_P2_= new MatrixBuilderCL(&S_P2_, num_unks_p2, num_unks_p2);
    mL_P2P2_= new MatrixBuilderCL(&L_P2P2_, num_unks_scalarp2, num_unks_p2);
    mL_P2P2_stab_= new MatrixBuilderCL(&L_P2P2_stab_, num_unks_scalarp2, num_unks_p2);
    mM_ScalarP2_= new MatrixBuilderCL(&M_ScalarP2_, num_unks_scalarp2, num_unks_scalarp2);
    mA_ScalarP2_stab_= new MatrixBuilderCL(&A_ScalarP2_stab_, num_unks_scalarp2, num_unks_scalarp2);
}

void StokesIFAccumulator_P2P2CL::finalize_accumulation ()
{
    mA_P2_->Build();
    delete mA_P2_;
    mA_P2_stab_->Build();
    delete mA_P2_stab_;
    mB_P2P2_->Build();
    delete mB_P2P2_;
    mM_P2_->Build();
    delete mM_P2_;
    mS_P2_->Build();
    delete mS_P2_;
    mL_P2P2_->Build();
    delete mL_P2P2_;
    mL_P2P2_stab_->Build();
    delete mL_P2P2_stab_;
    mM_ScalarP2_->Build();
    delete mM_ScalarP2_;
    mA_ScalarP2_stab_->Build();
    delete mA_ScalarP2_stab_;

#ifndef _PAR
    std::cout << "StokesIF_P2P2:\t" << A_P2_.num_nonzeros() << " nonzeros in A, " << A_P2_stab_.num_nonzeros() << " nonzeros in A_stab, "
              << B_P2P2_.num_nonzeros() << " nonzeros in B, " << M_P2_.num_nonzeros() << " nonzeros in M, " << S_P2_.num_nonzeros() << " nonzeros in S, "
              << L_P2P2_.num_nonzeros() << " nonzeros in L, " << L_P2P2_stab_.num_nonzeros() << " nonzeros in L_stab, "
              << M_ScalarP2_.num_nonzeros() << " nonzeros in M_ScalarP2, " << A_ScalarP2_stab_.num_nonzeros() << " nonzeros in A_ScalarP2_stab\n";
#endif
    // std::cout << '\n';
}

void StokesIFAccumulator_P2P2CL::visit (const TetraCL& tet)
{
    ls_loc.assign(tet, lset, lset_bnd);

    if (isInCutMesh(ls_loc))
    {
        local_setup(tet);
        update_global_system();
    }

}

void StokesIFAccumulator_P2P2CL::local_setup (const TetraCL& tet)
{
    GetTrafoTr(T, det, tet);
    absdet= std::fabs(det);

    GetLocalNumbP2NoBnd(numP2, tet, P2Idx_);
    GetLocalNumbP2NoBnd(numScalarP2, tet, ScalarP2Idx_);
    localProblem.calcIntegrands(T, ls_loc, tet);
    localProblem.calc3DIntegrands(T, ls_loc, tet);
    localProblem.setupA_P2_inconsistent(locA_P2);
    localProblem.setupA_P2_stab(locA_P2_stab, absdet);
    localProblem.setupB_P2P2(locB_P2P2);
    localProblem.setupM_P2(locM_P2);
    localProblem.setupS_P2(locS_P2);
    localProblem.setupM_P2(locM_ScalarP2);
    localProblem.setupA_P2_stab(locA_ScalarP2_stab, absdet);
}

void StokesIFAccumulator_P2P2CL::update_global_system ()
{
    MatrixBuilderCL& mA_P2 = *mA_P2_;
    MatrixBuilderCL& mA_P2_stab= *mA_P2_stab_;
    MatrixBuilderCL& mB_P2P2 = *mB_P2P2_;
    MatrixBuilderCL& mM_P2 = *mM_P2_;
    MatrixBuilderCL& mS_P2 = *mS_P2_;
    MatrixBuilderCL& mL_P2P2 = *mL_P2P2_;
    MatrixBuilderCL& mL_P2P2_stab= *mL_P2P2_stab_;
    MatrixBuilderCL& mM_ScalarP2 = *mM_ScalarP2_;
    MatrixBuilderCL& mA_ScalarP2_stab= *mA_ScalarP2_stab_;

    for(int i= 0; i < 10; ++i) {
        const IdxT ii= numP2[i];
        if (ii==NoIdx) continue;
        for(int j=0; j <10; ++j) {
            const IdxT jj= numP2[j];
            if (jj==NoIdx) continue;
            for (int k=0; k<3; ++k) {
                for (int l=0; l<3; ++l) {
                    mA_P2(ii+k, jj+l) += locA_P2[3*j+l][3*i+k];
                    mS_P2(ii+k, jj+l) += locS_P2[3*j+l][3*i+k];
                    if(k == l) {
                        mM_P2(ii+k, jj+l) += locM_P2[j][i];
                        mA_P2_stab(ii+k, jj+l) += locA_P2_stab[j][i];
                    } else {
                        mM_P2(ii+k, jj+l) += 0.;
                        mA_P2_stab(ii+k, jj+l) += 0.;
                    }
                }
             }
        }
        for (int j=0; j < 10; ++j) {
            if (numScalarP2[j]==NoIdx) continue;
            for (int k=0; k<3; ++k) {
                mB_P2P2(numScalarP2[j], ii+k) += locB_P2P2[j][3*i+k];
                mL_P2P2(numScalarP2[j], ii+k) += locL_P2P2[j][3*i+k];
                mL_P2P2_stab(numScalarP2[j], ii+k) += locL_P2P2_stab[j][3*i+k];
            }
        }
    }
    for(int i= 0; i < 10; ++i) {
        const IdxT ii= numScalarP2[i];
        if (ii==NoIdx) continue;
        for(int j=0; j < 10; ++j) {
            const IdxT jj= numScalarP2[j];
            if (jj==NoIdx) continue;
            mM_ScalarP2(ii, jj) += locM_ScalarP2[i][j];
            mA_ScalarP2_stab(ii, jj) += locA_ScalarP2_stab[i][j];
        }
    }
}

void SetupStokesIF_P2P2(const MultiGridCL& MG_, MatDescCL* A_P2, MatDescCL* A_P2_stab, MatDescCL* B_P2P2, MatDescCL* M_P2, MatDescCL* S_P2, MatDescCL* L_P2P2, MatDescCL* L_P2P2_stab, MatDescCL* M_ScalarP2, MatDescCL* A_ScalarP2_stab, const VecDescCL& lset, const LsetBndDataCL& lset_bnd, SurfOseenParam* param) {
  ScopeTimerCL scope("SetupStokesIF");
  StokesIFAccumulator_P2P2CL accu(lset, lset_bnd, *(A_P2->RowIdx), *(B_P2P2->RowIdx), A_P2->Data, A_P2_stab->Data, B_P2P2->Data, M_P2->Data, S_P2->Data, L_P2P2->Data, L_P2P2_stab->Data, M_ScalarP2->Data, A_ScalarP2_stab->Data, param);
  TetraAccumulatorTupleCL accus;
  //    MaybeAddProgressBar(MG_, "LapBeltr(P2) Setup", accus, RowIdx.TriangLevel());
  accus.push_back(&accu);
  accumulate(accus, MG_, A_P2->GetRowLevel(), A_P2->RowIdx->GetBndInfo());
}

void SetupInterfaceVectorRhsP1 (const MultiGridCL& mg, VecDescCL* v,
                                const VecDescCL& ls, const BndDataCL<>& lsetbnd, InstatVectorFunction f, double t)
{
    const IdxT num_unks= v->RowIdx->NumUnknowns();
    const Uint lvl = v->GetLevel();

    double totalsurfarea(0.);

    v->Clear(0);

    //IdxT num[10];

    LocalNumbP1CL n;

    std::cout << "entering SetupInterfaceVectorRhsP1: " << num_unks << " dof... ";

    LocalP1CL<> p1[4];
    P1DiscCL::GetP1Basis(p1);

//    LocalP2CL<> p2[10];
//    for(int i=0; i<10; i++)  // P2-Basis-Functions
//      p2[i][i] = 1.;

    Quad5_2DCL<double> q1[4]; //basis function

    // double det, absdet;

    InterfaceTriangleCL triangle;

    Point3DCL e_k;

    DROPS_FOR_TRIANG_CONST_TETRA(mg, lvl, it) {
        triangle.Init(*it, ls, lsetbnd);
        if (triangle.Intersects()) { // We are at the phase boundary.
          n.assign(*it, *v->RowIdx, lsetbnd);
          //GetLocalNumbP2NoBnd(num, *it, *v->RowIdx);

            for (int ch= 0; ch < 8; ++ch) {
                triangle.ComputeForChild(ch);
                for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri) {
                  //P2DiscCL::GetP2Basis(q2, &triangle.GetBary(tri));
                    for (int i= 0; i < 4; ++i)
                      q1[i].assign(p1[i], &triangle.GetBary(tri));

                    Quad5_2DCL<Point3DCL> qf(*it, &triangle.GetBary(tri), f,t);
                    Quad5_2DCL<double> surfarea(1.0), rhs;
                    totalsurfarea += surfarea.quad(triangle.GetAbsDet(tri));

                    for (int i= 0; i < 4; ++i) {
                        if (n.num[i] == NoIdx) continue;
                        //rhs = qf*q2[i];
                        for (int k=0; k<3; ++k) {
                            e_k= DROPS::std_basis<3>(k+1);
                            rhs = dot(e_k, qf)*q1[i];
                            v->Data[n.num[i] + k] += rhs.quad(triangle.GetAbsDet(tri));
                        }
                    }
                }
            }
        }
    }
    std::cout << " VectorRhs set up." << std::endl;
    std::cout << "Total surface area: "<<totalsurfarea<<std::endl;
    
}

void SetupInterfaceVectorRhsP2 (const MultiGridCL& mg, VecDescCL* v,
                                const VecDescCL& ls, const BndDataCL<>& lsetbnd, InstatVectorFunction f)
{
    const IdxT num_unks= v->RowIdx->NumUnknowns();
    const Uint lvl = v->GetLevel();

    double totalsurfarea(0.);

    v->Clear(0);

    //IdxT num[10];

    LocalNumbP2CL n;

    std::cout << "entering SetupInterfaceVectorRhsP2: " << num_unks << " dof... ";

    LocalP2CL<> p2[10];
    P2DiscCL::GetP2Basis(p2);

//    LocalP2CL<> p2[10];
//    for(int i=0; i<10; i++)  // P2-Basis-Functions
//      p2[i][i] = 1.;

    Quad5_2DCL<double> q2[10]; //basis function

    // double det, absdet;

    InterfaceTriangleCL triangle;

    Point3DCL e_k;

    DROPS_FOR_TRIANG_CONST_TETRA(mg, lvl, it) {
        triangle.Init(*it, ls, lsetbnd);
        if (triangle.Intersects()) { // We are at the phase boundary.
          n.assign(*it, *v->RowIdx, lsetbnd);
          //GetLocalNumbP2NoBnd(num, *it, *v->RowIdx);

            for (int ch= 0; ch < 8; ++ch) {
                triangle.ComputeForChild(ch);
                for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri) {
                  //P2DiscCL::GetP2Basis(q2, &triangle.GetBary(tri));
                    for (int i= 0; i < 10; ++i)
                      q2[i].assign(p2[i], &triangle.GetBary(tri));

                    Quad5_2DCL<Point3DCL> qf(*it, &triangle.GetBary(tri), f);
                    Quad5_2DCL<double> surfarea(1.0), rhs;
                    totalsurfarea += surfarea.quad(triangle.GetAbsDet(tri));

                    for (int i= 0; i < 10; ++i) {
                        if (n.num[i] == NoIdx) continue;
                        //rhs = qf*q2[i];
                        for (int k=0; k<3; ++k) {
                            e_k= DROPS::std_basis<3>(k+1);
                            rhs = dot(e_k, qf)*q2[i];
                            v->Data[n.num[i] + k] += rhs.quad(triangle.GetAbsDet(tri));
                        }
                    }
                }
            }
        }
    }
    std::cout << " VectorRhs set up." << std::endl;
    std::cout << "Total surface area: "<<totalsurfarea<<std::endl;
}


void P1Init (InstatScalarFunction icf, VecDescCL& ic, const MultiGridCL& mg, double t)
{
    const Uint lvl= ic.GetLevel(),
               idx= ic.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX(mg, lvl, it) {
        if (it->Unknowns.Exist(idx))
            ic.Data[it->Unknowns(idx)]= icf(it->GetCoord(), t);
    }
    ic.t= t;
}

InterfaceDebugP2CL::InterfaceDebugP2CL (const InterfaceCommonDataP2CL& cdata)
    : cdata_(cdata), to_iface(0), ref_dp(0), ref_abs_det(0), true_area(-1.)
{}

void InterfaceDebugP2CL::begin_accumulation ()
{
    max_dph_err= 0;
    surfacemeasP1 = 0.;
    surfacemeasP2 = 0.;
    max_absdet_err= 0.;
    max_dph2_err= 0;
}

void  InterfaceDebugP2CL::finalize_accumulation()
{
    std::cout << "max_dph_err: " << max_dph_err
        << "\nsurfacemeasP1: " << surfacemeasP1;
    if  (true_area > 0.)
        std::cout << " rel. error: " << std::abs(surfacemeasP1 - true_area)/true_area;
    std::cout << "\nsurfacemeasP2: " << surfacemeasP2;
    if  (true_area > 0.)
        std::cout << " rel. error: " << std::abs(surfacemeasP2 - true_area)/true_area;
    std::cout << "\nmax_absdet_err: " << max_absdet_err
              << "\nmax_dph2_err: " << max_dph2_err << std::endl;
}

void InterfaceDebugP2CL::visit (const TetraCL& t)
{
//     std::cout << "Tetra Id: " << t.GetId().GetIdent() << std::endl;
    const InterfaceCommonDataP2CL& cdata= cdata_.get_clone();
    if (cdata.empty())
        return;

    if (to_iface != 0) {
        const Uint sys= to_iface->RowIdx->GetIdx();
        for (Uint i= 0; i < 4; ++i) {
            cdata.quaqua.set_point(&t, std_basis<4>(i + 1))
                        .base_point();
            const TetraBaryPairT& tb= cdata.quaqua.get_base_point();
            Point3DCL offset= t.GetVertex(i)->GetCoord() - GetWorldCoord(*tb.first, tb.second);
            const size_t dof= t.GetVertex(i)->Unknowns(sys);
            std::copy(Addr(offset), Addr(offset) + 3, &to_iface->Data[dof]);
        }
    }

    for (SurfacePatchCL::const_vertex_iterator it= cdata.surf.vertex_begin(); it != cdata.surf.vertex_end(); ++it) {
        cdata.quaqua.set_point(&t, *it)
             .jacobian();
//         const TetraBaryPairT& tb= cdata.quaqua.get_base_point();
        const Point3DCL& x= GetWorldCoord(t, *it);
//         const Point3DCL& xb= GetWorldCoord(*tb.first, tb.second);
//         std::cout  << "    |x-xb|: " << (x - xb).norm();

        const SMatrixCL<3,3>& dph= cdata.quaqua.get_jacobian();
        SMatrixCL<3,3> diff_dp= ref_dp != 0 ? dph - ref_dp(x, 0.) : SMatrixCL<3,3>();
        const double dph_err= std::sqrt(frobenius_norm_sq(diff_dp));
        max_dph_err= std::max(max_dph_err, dph_err);
//         std::cout  << " |dph -dp|_F: " << dph_err;

        const double absdet= abs_det(t, *it, cdata.quaqua, cdata.surf),
                     absdet_err= ref_abs_det != 0 ? std::abs(absdet - ref_abs_det(t, *it, cdata.surf)) : 0.;
        max_absdet_err= std::max(max_absdet_err, absdet_err);
//         std::cout  << " |\\mu - \\mu^s|: " << absdet_err << std::endl;


        Point3DCL n=x/x.norm();
        SMatrixCL<3,3> diff_dp2 = diff_dp - outer_product(n, transp_mul(diff_dp, n));
        Point3DCL nh;
        Point3DCL v1 = GetWorldCoord(t, cdata.surf.vertex_begin()[1]) - GetWorldCoord(t, cdata.surf.vertex_begin()[0]),
                  v2 = GetWorldCoord(t, cdata.surf.vertex_begin()[2]) - GetWorldCoord(t, cdata.surf.vertex_begin()[0]);
        cross_product(nh, v1, v2);
        nh/= nh.norm();
        diff_dp2 = diff_dp2 - outer_product(diff_dp2*nh, nh);
        const double dph2_err= std::sqrt(frobenius_norm_sq(diff_dp2));
//         std::cout  << " |P(dph -dp)\\hat P|_F: " << dph2_err << std::endl;
        max_dph2_err= std::max(max_dph2_err, dph2_err);
    }
    surfacemeasP1+= quad_2D(std::valarray<double>(1., cdata.qdom.vertex_size()), cdata.qdom);
    surfacemeasP2+= quad_2D(cdata.qdom_projected.absdets(), cdata.qdom);
}

// Works only if the nodes of the quadrature rule are not shared across facets in QuadDomain2DCL.
void ProjectedQuadDomain2DCL::assign (const SurfacePatchCL& p, const QuadDomain2DCL& qdomarg, const QuaQuaMapperCL& quaqua)
{
    qdom= &qdomarg;
    resize(qdom->vertex_size());
    QuadDomain2DCL::const_vertex_iterator qv= qdom->vertex_begin();

    if (!compute_absdets_) {
        for (Uint i= 0; i < qdom->vertex_size(); ++i, ++qv) {
            quaqua.set_point(quaqua.get_tetra(), *qv)
                  .base_point();
            vertexes_[i]= quaqua.get_base_point();
        }
        return;
    }

    QRDecompCL<3,2> qr;
    SMatrixCL<3,2>& M= qr.GetMatrix();
    SMatrixCL<3,2> U;
    Point3DCL tmp;
    SMatrixCL<2,2> Gram;
    const Bary2WorldCoordCL b2w(*quaqua.get_tetra());
    const SurfacePatchCL::const_vertex_iterator pv= p.vertex_begin();

    const Uint nodes_per_facet= qdom->vertex_size()/p.facet_size();
    for (Uint i= 0; i < qdom->vertex_size(); ++i, ++qv) {
        if (i % nodes_per_facet == 0) {
            const SurfacePatchCL::FacetT& facet= p.facet_begin()[i/nodes_per_facet];
            M.col(0, b2w(pv[facet[1]]) - b2w(pv[facet[0]]));
            M.col(1, b2w(pv[facet[2]]) - b2w(pv[facet[0]]));
            qr.prepare_solve();
            for (Uint j= 0; j < 2; ++j) {
                tmp= std_basis<3>(j + 1);
                qr.apply_Q(tmp);
                U.col(j, tmp);
            }
        }
        quaqua.set_point(quaqua.get_tetra(), *qv)
              .jacobian();
        vertexes_[i]= quaqua.get_base_point();
        Gram= GramMatrix(quaqua.get_jacobian()*U);
        absdets_[i]= std::sqrt(Gram(0,0)*Gram(1,1) - Gram(0,1)*Gram(1,0));
    }
}

void gradient_trafo (const TetraCL& tet, const BaryCoordCL& xb, const QuaQuaMapperCL& quaqua, const SurfacePatchCL& p,
                     SMatrixCL<3,3>& W)
{
    // Compute the basepoint b.
    quaqua.set_point(&tet, xb)
          .jacobian();
    const TetraBaryPairT& tb= quaqua.get_base_point();

    // nl(x)
    if (p.normal_empty())
        p.compute_normals(tet);
    Point3DCL nl= p.normal_begin()[0];
    // Evaluate the normal to the interface in b, n(y).
    Point3DCL n= quaqua.local_ls_grad(*tb.first, tb.second);
    n/= n.norm();

    // Dp_h(x)^T
    const SMatrixCL<3,3>& dph= quaqua.get_jacobian();
    SMatrixCL<3,3> dphT;
    assign_transpose(dphT, dph);

    W= dphT + outer_product(nl, n/inner_prod(nl, n) - dph*nl);
}

//////////////
/// =Methods with transport on the domain=

class TransportP2FunctionCL
/// Solve D/Dt u = 0, u(t^0) given, with SDFEM
    {
    public:
        typedef BndDataCL<>    BndDataT;
        typedef P1EvalCL<double, const BndDataT, VecDescCL>       DiscSolCL;
        typedef P1EvalCL<double, const BndDataT, const VecDescCL> const_DiscSolCL;

        const IdxDescCL&    p1idx;

    private:
        MultiGridCL&        MG_;
        double              SD_,    ///< streamline diffusion
                theta_,
                dt_;
        MatrixCL            L_;
        BndDataT            Bnd_;
        GSPcCL              pc_;
        GMResSolverCL<GSPcCL>  gm_;

    public:
        MatrixCL E_old,
                E_,
                H_old,
                H_;

        template<class DiscVelSolT>
        TransportP2FunctionCL (MultiGridCL& mg, const DiscVelSolT& v_old, const IdxDescCL& thep1idx, double dt, double theta= 0.5, double SD= 0., int iter= 1000, double tol= 1e-7)
                : p1idx( thep1idx), MG_( mg), SD_( SD),
                  theta_( theta), dt_( dt), Bnd_( BndDataT( mg.GetBnd().GetNumBndSeg())),
                  gm_( pc_, 500, iter, tol, true) {
            if (theta_ != 1.)
                SetupSystem( v_old, E_old, H_old);
        }

        /// \remarks call SetupSystem \em before calling SetTimeStep!
        template<class DiscVelSolT>
        void SetupSystem (const DiscVelSolT&, MatrixCL& E, MatrixCL& H);
        /// perform one time step
        template <class DiscVelSolT>
        void DoStep (VectorCL& u, const DiscVelSolT& /* new velocity*/);

    };

    template<class DiscVelSolT>
    void TransportP2FunctionCL::SetupSystem (const DiscVelSolT& vel, MatrixCL& E, MatrixCL& H)
// Sets up the stiffness matrices:
// E is of mass matrix type:    E_ij = ( v_j       , v_i + SD * u grad v_i )
// H describes the convection:  H_ij = ( u grad v_j, v_i + SD * u grad v_i )
// where v_i, v_j denote the ansatz functions.
    {
        const IdxT num_unks= p1idx.NumUnknowns();
        const Uint lvl= p1idx.TriangLevel();

        SparseMatBuilderCL<double> bE(&E, num_unks, num_unks),
                bH(&H, num_unks, num_unks);
        IdxT Numb[10];

        std::cout << "entering TransportP1Function::SetupSystem: " << num_unks << "  unknowns. ";
        std::cout << "SD_: " << SD_ << " dt_: " << dt_ << " theta_ : " << theta_ << "\n";

        // fill value part of matrices
        Quad5CL<Point3DCL>  u_loc;
        Quad5CL<Point3DCL> GRef[10];
        Quad5CL<Point3DCL> Grad[10];
        Quad5CL<> u_Grad[10], // fuer u grad v_i
                p2[10];


        double det, absdet, h_T;

        LocalP2CL<> p2dummy;
        for (int i= 0; i < 10; ++i) {
            p2dummy[i]= 1.0;
            p2[i].assign( p2dummy);
            p2dummy[i]= 0.0;
        }

        P2DiscCL::GetGradientsOnRef(GRef);

        DROPS_FOR_TRIANG_CONST_TETRA( const_cast<const MultiGridCL&>( MG_), lvl, sit) {

            GetLocalNumbP2NoBnd( Numb, *sit, p1idx);

            SMatrixCL<3,3> M;
            const Point3DCL& pt0= sit->GetVertex(0)->GetCoord();
            for(Uint i=0; i<3; ++i)
                for(Uint j=0; j<3; ++j)
                    M[j*3+i]= sit->GetVertex(i+1)->GetCoord()[j] - pt0[j];
            // M[j][i]= sit->GetVertex(i+1)->GetCoord()[j] - pt0[j];
            det= M(0,0) * (M(1,1) *M(2,2)  - M(1,2) *M(2,1) )
                 - M(0,1) * (M(1,0) *M(2,2)  - M(1,2) *M(2,0) )
                 + M(0,2) * (M(1,0) *M(2,1)  - M(1,1)*M(2,0) );
            absdet= std::fabs( det);
            h_T= std::pow( absdet, 1./3.);
            P2DiscCL::GetGradients( Grad, GRef, M);

            u_loc.assign( *sit, vel);
            for(int i=0; i<10; ++i)
                u_Grad[i]= dot( Grad[i], u_loc);

            ////////////////////!!!!!!!!!!!!!!!!!!!!!!!!!?????????????????????????????!!!!!
            /// \todo fixed limit for maxV (maxV_limit), any better idea?
            double maxV = 1.; // scaling of SD parameter
            for(int i= 0; i < 10; ++i)    // assemble row Numb[i]
                for(int j= 0; j < 10; ++j) {
                    // E is of mass matrix type:    E_ij = ( v_j       , v_i + SD * u grad v_i )
                    bE( Numb[i], Numb[j])+= P2DiscCL::GetMass(i,j) * absdet
                                            + Quad5CL<>( u_Grad[i]*p2[j]).quad( absdet)*SD_/maxV*h_T;

                    // H describes the convection:  H_ij = ( u grad v_j, v_i + SD * u grad v_i )
                    bH( Numb[i], Numb[j])+= Quad5CL<>( u_Grad[j]*p2[i]).quad( absdet)
                                            + Quad5CL<>(u_Grad[i]*u_Grad[j]).quad( absdet) * SD_/maxV*h_T;
                }
        }
        bE.Build();
        bH.Build();
        std::cout << E.num_nonzeros() << " nonzeros in E, "
                  << H.num_nonzeros() << " nonzeros in H! " << std::endl;
    }

    template <class DiscVelSolT>
    void TransportP2FunctionCL::DoStep (VectorCL& u, const DiscVelSolT& vel)
    {
        SetupSystem( vel, E_, H_);
        L_.clear();
        L_.LinComb( 1./dt_, E_, theta_, H_);
        VectorCL rhs( (1./dt_)*u);

//    WriteToFile( L_, "chartranspL.txt", "L_");

//    WriteToFile( E_old, "chartranspE_old.txt", "E_old");

        if (theta_ != 1.) {
            GMResSolverCL<GSPcCL> gm( gm_);
            VectorCL tmp( rhs.size());
            gm.Solve( E_old, tmp, VectorCL( H_old*u), p1idx);
            std::cout << "TransportP1FunctionCL::DoStep rhs: res = " << gm.GetResid() << ", iter = " << gm.GetIter() << std::endl;
            rhs-= (1. - theta_)*tmp;
        }

        gm_.Solve( L_, u, VectorCL(E_*rhs), p1idx);

        std::cout << "TransportP1FunctionCL::DoStep: res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
    }


    void SurfacePDEP1BaseCL::SetTheta (double theta)
    {
        if (theta >= 0. && theta <= 1.)
            theta_= theta;
    }


void SurfactantP1BaseCL::SetInitialValue (InstatScalarFunction icf, double t)
{
    P1Init (icf, ic, MG_, t);
}

void SurfactantP1BaseCL::SetRhs (InstatScalarFunction rhs)
{
    rhs_fun_= rhs;
}

void SurfactantP1BaseCL::InitTimeStep ()
{
    iface_old.SetIdx( &oldidx_);

    if (oldidx_.NumUnknowns() > 0)
        oldidx_.DeleteNumbering(MG_);
    oldidx_.swap(idx);
    oldic_.resize(ic.Data.size());
    oldic_= ic.Data;
    oldt_= ic.t;

    oldls_.RowIdx= lset_vd_.RowIdx;
    oldls_.Data.resize(lset_vd_.Data.size());
    oldls_.Data= lset_vd_.Data;
    oldls_.t= lset_vd_.t;

    oldv_.SetIdx(v_->RowIdx);
    oldv_.Data= v_->Data;
    oldv_.t= v_->t;
}

//////////////////////////////////////////////
///SurfactantNarrowBandStblP1CL


    void SurfactantNarrowBandStblP1CL:: Update()
    {

        ScopeTimerCL timer( "SurfactantNarrowBandStblP1CL::Update");
        std::cout << "SurfactantNarrowBandStblP1CL::Update:\n";

        IdxDescCL* cidx= ic.RowIdx;
        Mass.Data.clear();
        Mass.SetIdx( cidx, cidx);
        Laplace.Data.clear();
        Laplace.SetIdx( cidx, cidx);

        Volume_stab.Data.clear();
        Volume_stab.SetIdx( cidx, cidx);
        Conv.Data.clear();
        Conv.SetIdx( cidx, cidx);
        Massd.Data.clear();
        Massd.SetIdx( cidx, cidx);
        VecDescCL vd_load( &idx);
        TetraAccumulatorTupleCL accus;
        InterfaceCommonDataP1CL cdata( lset_vd_, lsetbnd_);
        accus.push_back( &cdata);

        InterfaceMatrixAccuCL<LocalInterfaceMassP1CL, InterfaceCommonDataP1CL> mass_accu( &Mass, LocalInterfaceMassP1CL(), cdata, "mass");
        //InterfaceMatrixAccuP1CL<LocalInterfaceMassP1CL> mass_accu( &M, LocalInterfaceMassP1CL(), cdata, "mass");
        accus.push_back( &mass_accu);


        //InterfaceMatrixAccuP1CL<LocalLaplaceBeltramiP1CL> lb_accu( &A, LocalLaplaceBeltramiP1CL( D_), cdata, "Laplace-Beltrami");/// To implement method 1 using surface gradient
        InterfaceMatrixAccuCL<LocalLaplaceBeltramiP1CL, InterfaceCommonDataP1CL> lb_accu( &Laplace, LocalLaplaceBeltramiP1CL( D_), cdata, "Laplace-Beltrami");
        accus.push_back( &lb_accu);

        NarrowBandCommonDataP1CL bdata( lset_vd_, lsetbnd_,width_);
        accus.push_back( &bdata);

        //NarrowBandMatrixAccuP1CL<LocalNormalLaplaceBulkP1CL> sb_accu( &Sb, LocalNormalLaplaceBulkP1CL(D_,dt_,normal_,ic.t), bdata, "NormalLaplaceBulk"); ///4* To implement the stabilization Normal gradient
        NarrowBandMatrixAccuP1CL<LocalNormalLaplaceBulkP1CL> sb_accu( &Volume_stab, LocalNormalLaplaceBulkP1CL(1. ,dt_,normal_,ic.t), bdata, "NormalLaplaceBulk"); ///4* To implement the stabilization Normal gradient
        accus.push_back( &sb_accu);

        //  accus.push_back_acquire( make_wind_dependent_matrixP1_accu<LocalInterfaceConvectionP1CL>( &C,  cdata,  make_P2Eval( MG_, Bnd_v_, *v_), "convection"));
        //  accus.push_back_acquire( make_wind_dependent_matrixP1_accu<LocalInterfaceConvectionSkewP1CL>( &C,  cdata,  make_P2Eval( MG_, Bnd_v_, *v_), "convection"));
        //  accus.push_back_acquire( make_wind_dependent_matrixP1_accu1<LocalInterfaceConvectionSkew1P1CL>( &C,  cdata,  make_P2Eval( MG_, Bnd_v_, *v_),normal_,ic.t, "convection"));
        accus.push_back_acquire( make_wind_dependent_matrixP1_accu<LocalInterfaceConvectionP1CL>( &Conv,  cdata,  make_P2Eval( MG_, Bnd_v_, *v_), "convection"));

        //  accus.push_back_acquire( make_wind_dependent_matrixP1_accu<LocalInterfaceMassDivP1CL>( &Md, cdata,  make_P2Eval( MG_, Bnd_v_, *v_), "massdiv"));
        //  accus.push_back_acquire( make_wind_dependent_matrixP1_accu1<LocalInterfaceMassDivSkewP1CL>( &Md, cdata,  make_P2Eval( MG_, Bnd_v_, *v_),normal_,ic.t, "massdiv"));
        accus.push_back_acquire( make_wind_dependent_matrixP1_accu<LocalInterfaceMassDivP1CL>( &Massd, cdata,  make_P2Eval( MG_, Bnd_v_, *v_), "massdiv"));


        if (rhs_fun_)
            accus.push_back_acquire( new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>( &vd_load,
                                                                                                          LocalVectorP1CL( rhs_fun_, ic.t), cdata, "load"));
//        accus.push_back_acquire( new InterfaceVectorAccuP1CL<LocalVectorP1CL>( &vd_load, LocalVectorP1CL( rhs_fun_, ic.t), cdata, "load"));
        { ScopeTimerCL timer( "SurfactantExtensionP1CL::setup-Matrix");
            accumulate( accus, MG_, cidx->TriangLevel(), cidx->GetBndInfo());
            //accumulate( accus, MG_, cidx->TriangLevel(), cidx->GetMatchingFunction(), cidx->GetBndInfo());
        }
        load.resize( idx.NumUnknowns());
        load= vd_load.Data;

         //    WriteToFile( Mass.Data, "chartranspM.txt", "mass");
         //     WriteToFile( Laplace.Data, "chartranspA.txt", "Laplace-Beltrami");
         //     WriteToFile( Volume_stab.Data, "Stabilization.txt", "Stab-matrix");
        //     WriteToFile( Md.Data,"chartranspMd.txt","mass-div");
        //         WriteToFile( vd_load.Data,"chartranspload.txt","load");
        //    WriteToFile( C.Data,"convection.txt","convection");

        // std::cout << "SurfactantCharTransportP1CL::Update: Finished\n";
    }


// Implement Implicit Euler method
    void SurfactantNarrowBandStblP1CL::InitStep1 (double new_t)
    {
        // ScopeTimerCL timer( "SurfactantNarrowBandStblP1CL::InitStep");
        std::cout << "SurfactantNarrowBandStblP1CL::InitStep:\n";

        ic.t= new_t;
        dt_= ic.t - oldt_;

        //idx.GetXidx().SetBound( width_); //transfer the width_ to CreatNumbering
        idx.CreateNumbering( oldidx_.TriangLevel(), MG_, &lset_vd_, &lsetbnd_,width_); // InitTimeStep deletes oldidx_ and swaps idx and oldidx_.
        //idx_c.CreateNumbering( oldidx_c_.TriangLevel(), MG_, &lset_vd_, &lsetbnd_);

        std::cout << "new NumUnknowns: " << idx.NumUnknowns();

        full_idx.CreateNumbering( idx.TriangLevel(), MG_);
        std::cout << " full NumUnknowns: " << full_idx.NumUnknowns() << std::endl;
        DROPS::VecDescCL rhsext( &full_idx);


        /*
           oldls_.RowIdx= lset_vd_.RowIdx;
           oldls_.Data.resize( lset_vd_.Data.size());
           oldls_.Data= lset_vd_.Data;
           oldls_.t= lset_vd_.t;
       */

        VecDescCL rhs( &oldidx_);
        rhs.Data= (1./dt_)*oldic_;


        DROPS::ExtendP2( MG_, rhs, rhsext);

/*
    temp_ic.Reset();
    temp_ic.Data.resize( rhsext.Data.size());
    temp_ic.t=new_t;
    temp_ic.SetIdx( &full_idx);
    temp_ic.Data=rhsext.Data;
    WriteToFile( temp_ic.Data, "Extended.txt", "mass");*/

        VecDescCL rhs1( &idx);
        //   Restrict( MG_, rhsext, rhs1);
        MassMultiply(MG_,rhsext,rhs1,lset_vd_, lsetbnd_);
        //    MassMultiply(MG_,rhsext,rhs1,lset_vd_, lsetbnd_, *v_,Bnd_v_, dt_);
        //   MassMultiply(MG_,rhsext,rhs1,lset_vd_, lsetbnd_, *v_,Bnd_v_,normal_,ic.t, dt_);//characteristic method
        //  MassMultiply1(MG_,rhsext,rhs1,lset_vd_, lsetbnd_, *v_,Bnd_v_, dt_);
        //   WriteToFile( rhs1.Data, "RHS1.txt", "RHS1");
        rhs1_.resize( rhs1.Data.size());
        rhs1_= rhs1.Data;

        if(new_t/dt_==1)//used in InitStep 3.
        {
            rhsext1.Reset();
            rhsext1.SetIdx(&full_idx);
            rhsext1.Data=(std::sqrt(dt_))*rhsext.Data;// Since we choose a smaller time step dt*0.1
        }

    }

// Implement BDF2 method
//In the first step, we use the backward Euler method
    void SurfactantNarrowBandStblP1CL::InitStep2 (double new_t)
    {
        // ScopeTimerCL timer( "SurfactantcGP1CL::InitStep");
        std::cout << "SurfactantNarrowBandStblP1CL::InitStep:\n";

        ic.t= new_t;
        dt_= ic.t - oldt_;



        std::cout<<"Width of the Narrow Band: "<<2*width_<<std::endl;
        // idx.GetXidx().SetBound( width_); //transfer the width_ to CreatNumbering
        idx.CreateNumbering( oldidx_.TriangLevel(), MG_, &lset_vd_, &lsetbnd_,width_); // InitTimeStep deletes oldidx_ and swaps idx and oldidx_.
        std::cout << "new NumUnknowns: " << idx.NumUnknowns();

        full_idx.CreateNumbering( idx.TriangLevel(), MG_);
        std::cout << " full NumUnknowns: " << full_idx.NumUnknowns() << std::endl;
        DROPS::VecDescCL rhsext( &full_idx);
        /*
           oldls_.RowIdx= lset_vd_.RowIdx;
           oldls_.Data.resize( lset_vd_.Data.size());
           oldls_.Data= lset_vd_.Data;
           oldls_.t= lset_vd_.t;
       */

        iface.SetIdx( &idx);
        for (int i=0; i<idx.NumUnknowns();i++)
        {
            iface.Data[i]=1.;
        }
       /* iface.SetIdx( &oldidx_);
        for (int i=0; i<idx.NumUnknowns();i++)
        {
            iface_old.Data[i]=1.;
        }*/



        VecDescCL rhs( &oldidx_);
        rhs.Data= (1./dt_)*oldic_;
        DROPS::ExtendP2( MG_, rhs, rhsext);
        //  WriteToFile( rhsext.Data, "Extended00.txt", "mass");

/*
    temp_ic.Reset();
    temp_ic.Data.resize( rhsext.Data.size());
    temp_ic.t=new_t;
    temp_ic.SetIdx( &full_idx);
    temp_ic.Data=dt_*rhsext.Data;

    WriteToFile( temp_ic.Data, "Extended.txt", "mass");*/

        //fulltransport_= new TransportP2FunctionCL( MG_, make_P2Eval( MG_, Bnd_v_, oldv_), full_idx, dt_,  /*theta=*/theta_, /*SD=*/ 0.1, /*iter=*/ 2000, /*tol=*/ gm_.GetTol());
        //fulltransport_->DoStep( rhsext.Data, make_P2Eval( MG_, Bnd_v_, *v_));

        if(new_t/dt_==1)
        {
            std::cout<<"test---1"<<std::endl;
            VecDescCL rhs1( &idx);
            //Restrict( MG_, rhsext, rhs1);
            MassMultiply(MG_,rhsext,rhs1,lset_vd_, lsetbnd_);
            //std::cout<<"test---1"<<std::endl;
            rhs1_.resize( rhs1.Data.size());
            rhs1_= rhs1.Data;
        }
        else
        {
            VecDescCL rhs1( &idx);//present time step
            //Restrict( MG_, rhsext, rhs1);
            MassMultiply(MG_,rhsext,rhs1,lset_vd_, lsetbnd_);

            //    fulltransport_->DoStep( rhsext1.Data, make_P2Eval( MG_, Bnd_v_, *v_));
            VecDescCL rhs2( &idx);//last time step
            //Restrict( MG_, rhsext1, rhs2);
            MassMultiply(MG_,rhsext1,rhs2,lset_vd_, lsetbnd_);

            rhs1_.resize( rhs1.Data.size());
            rhs1_= rhs1.Data*2.0-rhs2.Data*0.5;
        }
        rhsext1.Reset();
        rhsext1.SetIdx(&full_idx);
        rhsext1.Data=rhsext.Data;

        //  WriteToFile( rhsext1.Data, "Extended01.txt", "mass");

    }

// Implement BDF2 method
// In the first step, we solve the problem with much smaller time step
    void SurfactantNarrowBandStblP1CL::InitStep3 (double new_t)
    {
        // ScopeTimerCL timer( "SurfactantNarrowBandStblP1CL::InitStep");
        std::cout << "SurfactantNarrowBandStblP1CL::InitStep:\n";

        ic.t= new_t;
        dt_= ic.t - oldt_;

        // idx.GetXidx().SetBound( width_); //transfer the width_ to CreatNumbering
        idx.CreateNumbering( oldidx_.TriangLevel(), MG_, &lset_vd_, &lsetbnd_,width_); // InitTimeStep deletes oldidx_ and swaps idx and oldidx_.
        std::cout << "new NumUnknowns: " << idx.NumUnknowns();

        full_idx.CreateNumbering( idx.TriangLevel(), MG_);
        std::cout << " full NumUnknowns: " << full_idx.NumUnknowns() << std::endl;
        DROPS::VecDescCL rhsext( &full_idx);
        /*
           oldls_.RowIdx= lset_vd_.RowIdx;
           oldls_.Data.resize( lset_vd_.Data.size());
           oldls_.Data= lset_vd_.Data;
           oldls_.t= lset_vd_.t;
       */


        VecDescCL rhs( &oldidx_);
        rhs.Data= (1./dt_)*oldic_;
        DROPS::ExtendP2( MG_, rhs, rhsext);
        //  WriteToFile( rhsext.Data, "Extended00.txt", "mass");

/*
    temp_ic.Reset();
    temp_ic.Data.resize( rhsext.Data.size());
    temp_ic.t=new_t;
    temp_ic.SetIdx( &full_idx);
    temp_ic.Data=dt_*rhsext.Data;

    WriteToFile( temp_ic.Data, "Extended.txt", "mass");*/

        if(new_t/dt_==1)
        {
            std::cout<<"Error occur: No solution provided on the first step, please use Initstep 2 instead"<<std::endl;
        }
        else
        {
            VecDescCL rhs1( &idx);//present time step
            //Restrict( MG_, rhsext, rhs1);
            MassMultiply(MG_,rhsext,rhs1,lset_vd_, lsetbnd_);

            //    fulltransport_->DoStep( rhsext1.Data, make_P2Eval( MG_, Bnd_v_, *v_));
            VecDescCL rhs2( &idx);//last time step
            //Restrict( MG_, rhsext1, rhs2);
            MassMultiply(MG_,rhsext1,rhs2,lset_vd_, lsetbnd_);

            rhs1_.resize( rhs1.Data.size());
            rhs1_= rhs1.Data*2.0-rhs2.Data*0.5;
        }
        rhsext1.Reset();
        rhsext1.SetIdx(&full_idx);
        rhsext1.Data=rhsext.Data;

        //  WriteToFile( rhsext1.Data, "Extended01.txt", "mass");

    }
//Implicit Euler method
    void SurfactantNarrowBandStblP1CL::DoStep1 ()
    {
        ic.SetIdx( &idx);
        std::cout<<"test---2  DoStep1()"<<std::endl;
        Update();
        L_.LinComb( 1./dt_, Mass.Data, 1., Laplace.Data, rho_, Volume_stab.Data, 1., Massd.Data,1.,Conv.Data);
        // if(new_t/(new_t-oldt_)==1)
        WriteToFile( L_, "Matrix.txt", "system");
        //  L_.LinComb( 1./dt_, M.Data, 1., Sb.Data,  1., Md.Data,1.,C.Data);
        // L_.LinComb( 1./dt_, M.Data, 5., A.Data, -4.,Sb.Data, 1., Md.Data,1.,C.Data);
        //  L_.LinComb( 1./dt_, M.Data, 1., A.Data, 1.,Sb.Data, 1., Md.Data);
        //const VectorCL therhs((M.Data*rhs1_) + load);
        const VectorCL therhs(rhs1_ + load);//we multiply M in initstep
        std::cout  <<"  Before solve: res = " << norm(therhs)<<" "<<norm(ic.Data)<<" "<< norm( L_*ic.Data - therhs) << std::endl;
        {
            std::cout<<"test 1"<<std::endl;
            ScopeTimerCL timer( "SurfactantExtensionP1CL::DoStep: Solve");
            gm_.Solve( L_, ic.Data, therhs, ic.RowIdx->GetEx());

            std::cout<<"test 2"<<std::endl;
        }
        std::cout << "SurfactantExtensionP1CL::DoStep: res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
    }

//BDF2 method
    void SurfactantNarrowBandStblP1CL::DoStep2 ()
    {

        ic.SetIdx( &idx);

        Update(); {
            L_.LinComb( 3./(2.*dt_), Mass.Data, 1., Laplace.Data,rho_, Volume_stab.Data, 1., Massd.Data,1.0,Conv.Data);
            // L_.LinComb( 3./(2.*dt_), M.Data, 1., Sb.Data, 1., Md.Data,1.0,C.Data);

            // const VectorCL therhs((M.Data*rhs1_) + load);
            const VectorCL therhs(rhs1_ + load);//we multiply M in initstep
            std::cout  <<"  Before solve: res = " << norm( L_*ic.Data - therhs) << std::endl;

            ScopeTimerCL timer( "SurfactantExtensionP1CL::DoStep: Solve");
            gm_.Solve( L_, ic.Data, therhs, ic.RowIdx->GetEx());
        }
        std::cout << "SurfactantExtensionP1CL::DoStep: res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
    }
    void SurfactantNarrowBandStblP1CL::CommitStep ()
    {
        full_idx.DeleteNumbering( MG_);
    }


//only back forward Euler
    void SurfactantNarrowBandStblP1CL::DoStep0 (double new_t)
    {
        ScopeTimerCL timer( "SurfactantNarrowP1CL::::DoStep");
        InitStep1( new_t);
        DoStep1();/**/
        CommitStep();
        //   WriteToFile( temp_ic.Data, "Extended.txt", "mass");

    }

//BDF2 method
    void SurfactantNarrowBandStblP1CL::DoStep (double new_t)
    {
        ScopeTimerCL timer( "SurfactantExtensionP1CL::::DoStep");
       /* InitStep2( new_t);
        if(new_t/(new_t-oldt_)==1)
        {
            std::cout<<"test---1 DoStep"<<std::endl;//first step~
            DoStep1();
        }
        else
            DoStep2();*//*First step-- Backward Euler*/

          InitStep3(new_t);
          DoStep2();
          CommitStep();/*In the first step, solve the equation by some smaller time stepsize-- */

        //   WriteToFile( temp_ic.Data, "Extended.txt", "mass");

    }


    void CahnHilliardP1BaseCL::SetInitialValue (InstatScalarFunction icmu, InstatScalarFunction icc,InstatScalarFunction ics, double t)
    {
        P1Init ( icc, ic, MG_, t);
        P1Init ( icmu, imu, MG_, t);
        P1Init ( ics, is, MG_, t);

    }

    void CahnHilliardP1BaseCL::SetRhs (InstatScalarFunction rhs3, InstatScalarFunction rhs4)
    {
        rhs_fun3_= rhs3;
        rhs_fun4_= rhs4;

    }

    void CahnHilliardP1BaseCL::InitTimeStep ()
    {
        if (oldoldidx_.NumUnknowns() > 0) {
            oldoldidx_.DeleteNumbering(MG_);
        }


        if (oldidx_.NumUnknowns() > 0) {
            oldoldidx_.swap( oldidx_);
            oldoldic_.resize( oldic_.size());
            oldoldic_= oldic_;
            oldidx_.DeleteNumbering(MG_);
        }

        oldidx_.swap( idx);
        oldic_.resize( ic.Data.size());
        oldic_= ic.Data;
        oldimu_.resize( imu.Data.size());
        oldimu_= imu.Data;
        oldis_.resize( is.Data.size());
        oldis_= is.Data;

        oldt_= ic.t;

        oldls_.RowIdx= lset_vd_.RowIdx;
        oldls_.Data.resize( lset_vd_.Data.size());
        oldls_.Data= lset_vd_.Data;
        oldls_.t= lset_vd_.t;

        oldv_.SetIdx( v_->RowIdx);
        oldv_.Data= v_->Data;
        oldv_.t= v_->t;
    }

    void CahnHilliardcGP1CL::Update()
    {

        IdxDescCL* cidx= ic.RowIdx;
        Mass.Data.clear();
        Mass.SetIdx( cidx, cidx);
        Laplace.Data.clear();
        Laplace.SetIdx( cidx, cidx);
        LaplaceM.Data.clear();
        LaplaceM.SetIdx( cidx, cidx);
        Conv.Data.clear();
        Conv.SetIdx( cidx, cidx);
        Massd.Data.clear();
        Massd.SetIdx( cidx, cidx);
        Volume_stab.Data.clear();
        Volume_stab.SetIdx( cidx, cidx);

        TetraAccumulatorTupleCL accus;
        InterfaceCommonDataP1CL cdata( lset_vd_, lsetbnd_);
        accus.push_back( &cdata);


        InterfaceMatrixAccuCL<LocalInterfaceMassP1CL, InterfaceCommonDataP1CL> mass_accu( &Mass, LocalInterfaceMassP1CL(), cdata, "mass");
        accus.push_back( &mass_accu);
        InterfaceMatrixAccuCL<LocalLaplaceBeltramiP1CL, InterfaceCommonDataP1CL> lb_accu( &Laplace, LocalLaplaceBeltramiP1CL( 1.), cdata, "Laplace-Beltrami");
        accus.push_back( &lb_accu);

      /*  InterfaceMatrixAccuCL<LocalLaplaceMobilityP1CL, InterfaceCommonDataP1CL> lbmob_accu( &LaplaceM, LocalLaplaceMobilityP1CL( make_P1Eval( MG_, Bnd_, ic), normal_, ic.t), cdata, "Laplace-Mobility");
        accus.push_back( &lbmob_accu);*/

        VecDescCL vd_oldic( &oldidx_);
        vd_oldic.Data= oldic_;
        vd_oldic.t= oldt_;
        accus.push_back_acquire( make_concentration_dependent_matrixP1_accu<LocalLaplaceMobilityP1CL>( &LaplaceM,  cdata, normal_, oldt_,  make_P1Eval( MG_, Bnd_, vd_oldic), "Laplace_Mobility"));

        accus.push_back_acquire( make_wind_dependent_matrixP1_accu<LocalInterfaceConvectionP1CL>( &Conv,  cdata,  make_P2Eval( MG_, Bnd_v_, *v_), "convection"));
        accus.push_back_acquire( make_wind_dependent_matrixP1_accu<LocalInterfaceMassDivP1CL>   ( &Massd, cdata,  make_P2Eval( MG_, Bnd_v_, *v_), "massdiv"));

        InterfaceMatrixAccuCL<LocalNormalLaplaceBulkP1CL, InterfaceCommonDataP1CL> normalstab_accu( &Volume_stab, LocalNormalLaplaceBulkP1CL( 1., dt_, normal_, oldt_), cdata, "Normal_stab");
        accus.push_back( &normalstab_accu);

        if (theta_ != 1.0) {
            Mass2.Data.clear();
            Mass2.SetIdx( cidx, cidx);
            InterfaceCommonDataP1CL* oldcdata= new InterfaceCommonDataP1CL( oldls_, lsetbnd_);
            accus.push_back_acquire( oldcdata);
            accus.push_back_acquire( new InterfaceMatrixAccuCL<LocalInterfaceMassP1CL, InterfaceCommonDataP1CL>( &Mass2, LocalInterfaceMassP1CL(), *oldcdata, "old mass"));
        }
        accumulate( accus, MG_, cidx->TriangLevel(), cidx->GetBndInfo());
        MatrixCL Precond3, Precond4;

        std::valarray<double> ones( 1., Mass.Data.num_rows());
        MatrixCL Ident(ones);


        Precond3.LinComb(1.0, Mass.Data, sigma_*dt_, LaplaceM.Data, dt_*rho_, Volume_stab.Data);
        Precond4.LinComb(1.0, Mass.Data, epsilon_*epsilon_, Laplace.Data, epsilon_*epsilon_*rho_, Volume_stab.Data);
        block_pc_.GetPC1().Reset(Precond3);
        block_pc_.GetPC2().Reset(Precond4);

//     WriteToFile( M.Data, "cGcGM.txt", "mass");
//     WriteToFile( A.Data, "cGcGA.txt", "Laplace-Beltrami");
//     WriteToFile( C.Data, "cGcGC.txt", "material derivative");
//     WriteToFile( Md.Data,"cGcGMd.txt","mass-div");

        // std::cout << "SurfactantP1CL::Update: Finished\n";
    }

    void CahnHilliardcGP1CL::InitStep(VectorCL& rhs3, VectorCL& rhs4, double new_t)
    {
        ic.t= new_t;
        imu.t= new_t;
        dt_= new_t - oldt_;
        idx.CreateNumbering( oldidx_.TriangLevel(), MG_, &lset_vd_, &lsetbnd_);
        std::cout << "new NumUnknowns: " << idx.NumUnknowns() << std::endl;
        ic.SetIdx( &idx);
        imu.SetIdx( &idx);
        iface.SetIdx( &idx);
        for (int i=0; i<idx.NumUnknowns();i++)
        {
            iface.Data[i]=1.;
        }


        VecDescCL vd_timeder( &idx),    // right-hand sides from integrals over the old/new interface
                vd_oldtimeder( &idx),
                vd_load3( &idx),
                vd_load4( &idx),
                vd_oldres( &idx),
                vd_oldload3( &idx),
                vd_oldload4( &idx),
                well_potential( &oldidx_),//double-well potential
                vd_well( &idx),
                vd_oldic( &oldidx_);  // the initial data.

        vd_oldic.Data= oldic_;
        vd_oldic.t= oldt_;

        TetraAccumulatorTupleCL accus;
        InterfaceCommonDataP1CL cdata( lset_vd_, lsetbnd_);
        accus.push_back( &cdata);

        InterfaceVectorAccuCL<LocalMatVecP1CL<LocalInterfaceMassP1CL>, InterfaceCommonDataP1CL> mass_accu( &vd_timeder,
           LocalMatVecP1CL<LocalInterfaceMassP1CL>( LocalInterfaceMassP1CL(), &vd_oldic), cdata, "mixed-mass");
        accus.push_back( &mass_accu);

        for (int i=0; i<well_potential.Data.size();i++)
        {
            well_potential.Data[i]=Potential_prime_function(oldic_[i]);//BDF1
        }

        InterfaceVectorAccuCL<LocalMatVecP1CL<LocalInterfaceMassP1CL>, InterfaceCommonDataP1CL> well_accu( &vd_well,
           LocalMatVecP1CL<LocalInterfaceMassP1CL>( LocalInterfaceMassP1CL(), &well_potential), cdata, "chemical_potential");
        accus.push_back( &well_accu);

        if (rhs_fun3_)
            accus.push_back_acquire( new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>( &vd_load3,
                LocalVectorP1CL( rhs_fun3_, new_t), cdata, "load3"));

        if (rhs_fun4_)
            accus.push_back_acquire( new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>( &vd_load4,
                LocalVectorP1CL( rhs_fun4_, new_t), cdata, "load4"));


        if (theta_ == 1.0) {
            accumulate( accus, MG_, idx.TriangLevel(), idx.GetBndInfo());
            rhs3 = VectorCL( theta_*(vd_timeder.Data + dt_*vd_load3.Data));
            rhs4 = VectorCL( (-1.)*vd_well.Data + (S_)*vd_timeder.Data + vd_load4.Data);
            return;
        }

        InterfaceCommonDataP1CL oldcdata( oldls_, lsetbnd_);
        accus.push_back( &oldcdata);

        if (rhs_fun3_)
            accus.push_back_acquire( new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>( &vd_oldload3,
                    LocalVectorP1CL( rhs_fun3_, oldt_), oldcdata, "load3 on old iface"));

        if (rhs_fun4_)
            accus.push_back_acquire( new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>( &vd_oldload4,
                    LocalVectorP1CL( rhs_fun4_, oldt_), oldcdata, "load4 on old iface"));

        InterfaceVectorAccuCL<LocalMatVecP1CL<LocalInterfaceMassP1CL>, InterfaceCommonDataP1CL> old_mass_accu( &vd_oldtimeder,
           LocalMatVecP1CL<LocalInterfaceMassP1CL>( LocalInterfaceMassP1CL(), &vd_oldic), oldcdata, "mixed-mass on old iface");
        accus.push_back( &old_mass_accu);

        InterfaceVectorAccuCL<LocalMatVecP1CL<LocalLaplaceBeltramiP1CL>, InterfaceCommonDataP1CL> old_lb_accu( &vd_oldres,
            LocalMatVecP1CL<LocalLaplaceBeltramiP1CL>( LocalLaplaceBeltramiP1CL( sigma_), &vd_oldic), oldcdata, "Laplace-Beltrami on old iface");
        accus.push_back( &old_lb_accu);

        accus.push_back_acquire( make_wind_dependent_vectorP1_accu<LocalInterfaceConvectionP1CL>( &vd_oldres, &vd_oldic,
           oldcdata,  make_P2Eval( MG_, Bnd_v_, oldv_), "convection on old iface"));
        accus.push_back_acquire( make_wind_dependent_vectorP1_accu<LocalInterfaceMassDivP1CL>   ( &vd_oldres, &vd_oldic,
           oldcdata,  make_P2Eval( MG_, Bnd_v_, oldv_), "mass-div on old iface"));


        accumulate( accus, MG_, idx.TriangLevel(), idx.GetBndInfo());

        rhs3 = VectorCL( theta_*vd_timeder.Data + (1. - theta_)*vd_oldtimeder.Data
                    + dt_*(theta_*vd_load3.Data + (1. - theta_)*(vd_oldload3.Data - vd_oldres.Data)));

        rhs4 = VectorCL((-1.)*vd_well.Data + (S_)*vd_timeder.Data + vd_load4.Data);

    }

    VectorCL CahnHilliardcGP1CL::InitStep3 (double new_t)
    {
        // ScopeTimerCL timer( "SurfactantcGP1CL::InitStep");
        // std::cout << "SurfactantcGP1CL::InitStep:\n";

        ic.t= new_t;
        dt_= new_t - oldt_;
        idx.CreateNumbering( oldidx_.TriangLevel(), MG_, &lset_vd_, &lsetbnd_); // InitTimeStep deletes oldidx_ and swaps idx and oldidx_.
        std::cout << "new NumUnknowns: " << idx.NumUnknowns() << std::endl;
        ic.SetIdx( &idx);

        for (int i=0; i<idx.NumUnknowns();i++)
        {
            iface.Data[i]=1.;
        }



        VecDescCL vd_timeder( &idx),    // right-hand sides from integrals over the old/new interface
                vd_oldtimeder( &idx),
                vd_load( &idx),
                vd_oldres( &idx),
                vd_oldload( &idx),
                vd_oldic( &oldidx_);  // the initial data.
                vd_oldic.Data= oldic_;
                vd_oldic.t= oldt_;

        TetraAccumulatorTupleCL accus;

        InterfaceCommonDataP1CL cdata( lset_vd_, lsetbnd_);
        accus.push_back( &cdata);

        InterfaceVectorAccuCL<LocalMatVecP1CL<LocalInterfaceMassP1CL>, InterfaceCommonDataP1CL> mass_accu( &vd_timeder,
                                                                                                           LocalMatVecP1CL<LocalInterfaceMassP1CL>( LocalInterfaceMassP1CL(), &vd_oldic), cdata, "mixed-mass");
        accus.push_back( &mass_accu);

        if (rhs_fun3_)
            accus.push_back_acquire( new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>( &vd_load, LocalVectorP1CL( rhs_fun3_, new_t), cdata, "load"));

        if (theta_ == 1.0) {
            accumulate( accus, MG_, idx.TriangLevel(), idx.GetBndInfo());
            return VectorCL( theta_*(vd_timeder.Data + dt_*vd_load.Data));
        }

     /*   InterfaceCommonDataP1CL oldcdata( oldls_, lsetbnd_);
        accus.push_back( &oldcdata);

        if (rhs_fun3_)
            accus.push_back_acquire( new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>( &vd_oldload, LocalVectorP1CL( rhs_fun3_, oldt_), oldcdata, "load on old iface"));

        InterfaceVectorAccuCL<LocalMatVecP1CL<LocalInterfaceMassP1CL>, InterfaceCommonDataP1CL> old_mass_accu( &vd_oldtimeder,
                                                                                                               LocalMatVecP1CL<LocalInterfaceMassP1CL>( LocalInterfaceMassP1CL(), &vd_oldic), oldcdata, "mixed-mass on old iface");
        accus.push_back( &old_mass_accu);
        InterfaceVectorAccuCL<LocalMatVecP1CL<LocalLaplaceBeltramiP1CL>, InterfaceCommonDataP1CL> old_lb_accu( &vd_oldres,
                                                                                                               LocalMatVecP1CL<LocalLaplaceBeltramiP1CL>( LocalLaplaceBeltramiP1CL( sigma_), &vd_oldic), oldcdata, "Laplace-Beltrami on old iface");
        accus.push_back( &old_lb_accu);
        accus.push_back_acquire( make_wind_dependent_vectorP1_accu<LocalInterfaceConvectionP1CL>( &vd_oldres, &vd_oldic,  oldcdata,  make_P2Eval( MG_, Bnd_v_, oldv_), "convection on old iface"));
        accus.push_back_acquire( make_wind_dependent_vectorP1_accu<LocalInterfaceMassDivP1CL>   ( &vd_oldres, &vd_oldic,  oldcdata,  make_P2Eval( MG_, Bnd_v_, oldv_), "mass-div on old iface"));

        accumulate( accus, MG_, idx_c.TriangLevel(), idx_c.GetBndInfo());
        return VectorCL( theta_*vd_timeder.Data + (1. - theta_)*vd_oldtimeder.Data
                         + dt_*(theta_*vd_load.Data + (1. - theta_)*(vd_oldload.Data - vd_oldres.Data)));*/
    }

    VectorCL CahnHilliardcGP1CL::InitStep4 (double new_t)
    {
        // ScopeTimerCL timer( "SurfactantcGP1CL::InitStep");
        // std::cout << "SurfactantcGP1CL::InitStep:\n";

        imu.t= new_t;
        //dt_= new_t - oldt_;
//        idx_mu.CreateNumbering( oldidx_mu_.TriangLevel(), MG_, &lset_vd_, &lsetbnd_); // InitTimeStep deletes oldidx_ and swaps idx and oldidx_.
        imu.SetIdx( &idx);


        VecDescCL well_potential( &oldidx_),//double-well potential
                //vd_oldload( &idx_mu),
                vd_load( &idx),
                vd( &idx),
                vd_well( &idx),
        vd_oldic( &oldidx_);  // the initial data.
        vd_oldic.Data= oldic_;
        vd_oldic.t= oldt_;
        /* vd_oldtimeder( &idx_mu),
           vd_oldres( &idx_mu),
           vd_oldimu( &oldidx_mu_);  // the initial data.*/


        TetraAccumulatorTupleCL accus;
        InterfaceCommonDataP1CL cdata( lset_vd_, lsetbnd_);
        accus.push_back( &cdata);

        for (int i=0; i<well_potential.Data.size();i++)
        {
            well_potential.Data[i]=Potential_prime_function(oldic_[i]);//BDF1
        }

        InterfaceVectorAccuCL<LocalMatVecP1CL<LocalInterfaceMassP1CL>, InterfaceCommonDataP1CL> mass_accu( &vd,
                                                                                                           LocalMatVecP1CL<LocalInterfaceMassP1CL>( LocalInterfaceMassP1CL(), &vd_oldic), cdata, "mass");
        accus.push_back( &mass_accu);

        InterfaceVectorAccuCL<LocalMatVecP1CL<LocalInterfaceMassP1CL>, InterfaceCommonDataP1CL> well_accu( &vd_well,
                                                                                                           LocalMatVecP1CL<LocalInterfaceMassP1CL>( LocalInterfaceMassP1CL(), &well_potential), cdata, "chemical_potential");
        accus.push_back( &well_accu);


        if (rhs_fun4_)
            accus.push_back_acquire( new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>( &vd_load, LocalVectorP1CL( rhs_fun4_, new_t), cdata, "load"));

        if (theta_ == 1.0) {
            accumulate( accus, MG_, idx.TriangLevel(), idx.GetBndInfo());
            return VectorCL( (-1.)*vd_well.Data + (S_)*vd.Data + theta_*(vd_load.Data));
        }

       /* InterfaceCommonDataP1CL oldcdata( oldls_, lsetbnd_);
        accus.push_back( &oldcdata);

        if (rhs_fun4_)
            accus.push_back_acquire( new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>( &vd_oldload, LocalVectorP1CL( rhs_fun4_, oldt_), oldcdata, "load on old iface"));
*/
//        InterfaceVectorAccuCL<LocalMatVecP1CL<LocalInterfaceMassP1CL>, InterfaceCommonDataP1CL> old_mass_accu( &vd_oldtimeder,
//                                                                                                               LocalMatVecP1CL<LocalInterfaceMassP1CL>( LocalInterfaceMassP1CL(), &vd_oldimu), oldcdata, "mixed-mass on old iface");
//        accus.push_back( &old_mass_accu);
//        InterfaceVectorAccuCL<LocalMatVecP1CL<LocalLaplaceBeltramiP1CL>, InterfaceCommonDataP1CL> old_lb_accu( &vd_oldres,
//                                                                                                               LocalMatVecP1CL<LocalLaplaceBeltramiP1CL>( LocalLaplaceBeltramiP1CL( sigma_), &vd_oldimu), oldcdata, "Laplace-Beltrami on old iface");
//        accus.push_back( &old_lb_accu);
//        accus.push_back_acquire( make_wind_dependent_vectorP1_accu<LocalInterfaceConvectionP1CL>( &vd_oldres, &vd_oldimu,  oldcdata,  make_P2Eval( MG_, Bnd_v_, oldv_), "convection on old iface"));
//        accus.push_back_acquire( make_wind_dependent_vectorP1_accu<LocalInterfaceMassDivP1CL>   ( &vd_oldres, &vd_oldimu,  oldcdata,  make_P2Eval( MG_, Bnd_v_, oldv_), "mass-div on old iface"));
        return VectorCL(theta_*vd_load.Data);//+(1-theta_)*vd_oldload.Data);
    }


    void CahnHilliardcGP1CL::DoStep(const VectorCL &rhs3, const VectorCL &rhs4) {
        Update();


        double c=1.0;//BDF1

        A_.LinComb(sigma_*dt_, LaplaceM.Data, rho_*dt_, Volume_stab.Data);
        B_.LinComb(c, Mass.Data, 0*dt_ , Massd.Data, dt_, Conv.Data);

        C_.LinComb(0.,   Laplace.Data, -epsilon_, Mass.Data);
        D_.LinComb(epsilon_*epsilon_,  Laplace.Data,
                //   -1., Gprimeprime.Data,
                  S_, Mass.Data,
                  rho_*epsilon_*epsilon_, Volume_stab.Data);

        //A_.LinComb(theta_, Mass.Data, dt_ * theta_, Laplace.Data, dt_ * theta_, Massd.Data, dt_ * theta_, Conv.Data);
        //D_.LinComb(1.0, A_, 0.0, Mass.Data);
        std::cout << "      Before solve: res3 = " << norm(A_ * imu.Data  +B_ * ic.Data- rhs3) << std::endl;
        std::cout << "      Before solve: res4 = " << norm(C_ * imu.Data  +D_ * ic.Data- rhs4) << std::endl;

        {
            ScopeTimerCL timer("CahnHilliardcGP1CL::DoStep: Solve");
            block_gm_.Solve(A_, B_, C_, D_, imu.Data, ic.Data, rhs3, rhs4, imu.RowIdx->GetEx(),ic.RowIdx->GetEx());
            /*ScopeTimerCL timer( "SurfactantcGP1CL::DoStep: Solve");
            MatrixCL L;
            L.LinComb( 1.0, Mass.Data, dt_*sigma_, Laplace.Data, 1*rho_*dt_, Volume_stab.Data, 0*dt_, Massd.Data, dt_, Conv.Data);
            gm_.Solve( L, ic.Data, rhs3, ic.RowIdx->GetEx());*/

            std::cout << "      Iterations of inner 3 and 4 solver on the last outer step"  << ": " << PCGSolver3_.GetIter() << "\t" << PCGSolver4_.GetIter() << '\n';
            std::cout << "      After  solve: res3 = " << norm(A_ * imu.Data  +B_ * ic.Data- rhs3) << std::endl;
            std::cout << "      After  solve: res4 = " << norm(C_ * imu.Data  +D_ * ic.Data- rhs4) << std::endl;
        }
        std::cout << "CahnHilliardcGP1CL::DoStep: res = " << block_gm_.GetResid() << ", iter = " << block_gm_.GetIter() << std::endl;
    }

    void CahnHilliardcGP1CL:: DoStep (double new_t)
    {
        ScopeTimerCL timer( "CahnHilliardcGP1CL::DoStep");

        /*VectorCL rhs3( InitStep3( new_t));
        VectorCL rhs4( InitStep4( new_t));*/

        VectorCL rhs3,rhs4;
        InitStep( rhs3, rhs4, new_t);
        DoStep( rhs3,rhs4);
        CommitStep();
    }

    void CahnHilliardcGP1CL:: DoStep0 (double new_t)
    {
        DoStep(new_t);
    }
    void CahnHilliardcGP1CL::CommitStep ()
    {
        return;
    }

    ///StblP1CL
    void CahnHilliardNarrowBandStblP1CL::   Update()
    {
        ScopeTimerCL timer( "CahnHilliardNarrowBandStblP1CL::Update");
        std::cout << "CahnHilliardNarrowBandStblP1CL::Update:\n";

        IdxDescCL* cidx= ic.RowIdx;
        Mass.Data.clear();
        Mass.SetIdx( cidx, cidx);
        Massrho.Data.clear();
        Massrho.SetIdx( cidx, cidx);
        Laplace.Data.clear();
        Laplace.SetIdx( cidx, cidx);
        LaplaceM.Data.clear();
        LaplaceM.SetIdx( cidx, cidx);
        LaplaceNon.Data.clear();
        LaplaceNon.SetIdx( cidx, cidx);

        Volume_stab.Data.clear();
        Volume_stab.SetIdx( cidx, cidx);
        Conv.Data.clear();
        Conv.SetIdx( cidx, cidx);
        Massd.Data.clear();
        Massd.SetIdx( cidx, cidx);

        Ident.Data.clear();
        /*std::valarray<double> ones( 1., Mass.Data.num_rows());
        Ident.Data=MatrixCL(ones);*/

        VecDescCL conc_ext( &idx);//extension of the old concentration to the current NarrowBand
        Restrict(MG_, conc_extrapol, conc_ext);

        VecDescCL sp_ext( &idx);//extension of the old concentration to the current NarrowBand
        Restrict(MG_, species_extrapol, sp_ext);
//        if (oldoldic_.size()>0)
//        conc_ext.Data =  2.0 * oldic_ - 1.0 * oldoldic_;
//        else conc_ext.Data=1.0*oldic_;

        TetraAccumulatorTupleCL accus;
        InterfaceCommonDataP1CL cdata( lset_vd_, lsetbnd_);//
        accus.push_back( &cdata);

        InterfaceMatrixAccuCL<LocalInterfaceMassP1CL, InterfaceCommonDataP1CL>
                mass_accu( &Mass, LocalInterfaceMassP1CL(), cdata, "mass");
        accus.push_back( &mass_accu);

        accus.push_back_acquire( make_concentration_dependent_matrixP1_accu<LocalInterfaceMassRhoP1CL>( &Massrho,  cdata, normal_, ic.t,  make_P1Eval( MG_, Bnd_, conc_ext), "Mass_rho"));


        InterfaceMatrixAccuCL<LocalLaplaceBeltramiP1CL, InterfaceCommonDataP1CL> lb_accu( &Laplace, LocalLaplaceBeltramiP1CL( 1.), cdata, "Laplace-Beltrami");
        accus.push_back( &lb_accu);
       /* InterfaceMatrixAccuCL<LocalLaplaceP1CL, InterfaceCommonDataP1CL> lb_accu( &Laplace, LocalLaplaceP1CL( normal_, ic.t), cdata, "Laplace");
        accus.push_back( &lb_accu);*/


        accus.push_back_acquire( make_concentration_dependent_matrixP1_accu<LocalLaplaceMobilityP1CL>( &LaplaceM,  cdata, normal_, ic.t,  make_P1Eval( MG_, Bnd_, conc_ext), "Laplace_Mobility"));

        accus.push_back_acquire( make_concentration_dependent_matrixP1_accu<LocalLaplaceNonlinearP1CL>( &LaplaceNon,  cdata, normal_, ic.t,  make_P1Eval( MG_, Bnd_, sp_ext), "Laplace_Nonlinear"));


        //WARNING
       NarrowBandCommonDataP1CL bdata( lset_vd_, lsetbnd_,width_);
       accus.push_back( &bdata);

        //InterfaceMatrixAccuCL<LocalFullGradientsP1CL, InterfaceCommonDataP1CL> sb_accu( &Volume_stab, LocalFullGradientsP1CL(normal_, oldt_), cdata, "Volume_stab");
        //InterfaceMatrixAccuCL<LocalNormalLaplaceBulkP1CL, InterfaceCommonDataP1CL> sb_accu( &Volume_stab, LocalNormalLaplaceBulkP1CL( 1., dt_, normal_, oldt_), cdata, "Normal_stab");
        //NarrowBandMatrixAccuP1CL<LocalFullGradientsP1CL> sb_accu( &Volume_stab, LocalFullGradientsP1CL(normal_, oldt_), bdata, "FullLaplaceBulk"); ///* To implement the stabilization Normal gradient
        NarrowBandMatrixAccuP1CL<LocalNormalLaplaceBulkP1CL> sb_accu( &Volume_stab, LocalNormalLaplaceBulkP1CL( 1., dt_, normal_, oldt_), bdata, "NormalLaplaceBulk"); ///* To implement the stabilization Normal gradient
        accus.push_back( &sb_accu);

       // accus.push_back_acquire( make_wind_dependent_matrixP1_accu<LocalInterfaceVelocityLaplaceP1CL>( &Volume_stab,  bdata, normal_, ic.t, make_P2Eval( MG_, Bnd_v_, *v_), "VelocityLaplace"));


        /*NarrowBandMatrixAccuP1CL<LocalInterfaceConvectionP1CL> conv_accu( &Conv, LocalInterfaceConvectionP1CL( 1., dt_, normal_, oldt_), bdata, "ConvBulk");
        accus.push_back( &conv_accu);*/
        accus.push_back_acquire( make_wind_dependent_matrixP1_accu<LocalInterfaceConvectionP1CL>( &Conv,  cdata,  make_P2Eval( MG_, Bnd_v_, *v_), "convection"));
        //accus.push_back_acquire( make_wind_dependent_matrixP1_accu<LocalInterfaceMassDivP1CL>( &Massd, cdata,  make_P2Eval( MG_, Bnd_v_, *v_), "massdiv"));


        {
            ScopeTimerCL timer( "CahnHillairdNarrowBandStblP1CL::setup-Matrix");
            accumulate( accus, MG_, cidx->TriangLevel(), cidx->GetBndInfo());
        }

        //    WriteToFile( Mass.Data, "chartranspM.txt", "mass");
        //     WriteToFile( Laplace.Data, "chartranspA.txt", "Laplace-Beltrami");
        //     WriteToFile( Volume_stab.Data, "Stabilization.txt", "Stab-matrix");
        //     WriteToFile( Md.Data,"chartranspMd.txt","mass-div");
        //         WriteToFile( vd_load.Data,"chartranspload.txt","load");
        //    WriteToFile( C.Data,"convection.txt","convection");

        // std::cout << "SurfactantCharTransportP1CL::Update: Finished\n";
        MatrixCL Precond3, Precond4;
        /*std::valarray<double> ones( 1., Mass.Data.num_rows());
        MatrixCL Ident(ones);*/
//        std::valarray<double> ones( 1., Mass.Data.num_rows());
//        Ident.Data=MatrixCL(ones);

//        VectorCL id,random;
//        id.resize(Mass.Data.num_rows(), 1);
//        random.resize(Mass.Data.num_rows(), 0);
//        for (int i=0; i<Mass.Data.num_rows(); i++)
//        {
//            random[i]=  - 5 + (rand() % 10)   ;
//        }
//        std::cout << "CHECKS:" << dot(Conv.Data *random,random) << "\t" <<  dot(Mass.Data*random,random)<< std::endl;
//        //Precond3.LinComb(1.0, Mass.Data, 1.0*rho_, Volume_stab.Data, 0.0, Ident.Data);
        //Precond4.LinComb(1.0, Mass.Data, 1.0*rho_, Volume_stab.Data, 0.0, Ident.Data);

        Precond3.LinComb(1.0, Mass.Data, 1.0*rho_, Volume_stab.Data, 1.0*sigma_, LaplaceM.Data);
        Precond4.LinComb(1.0, Mass.Data, 0.0, Conv.Data, 1.0*rho_*epsilon_*epsilon_, Volume_stab.Data, 1.0*epsilon_*epsilon_, Laplace.Data);
//        Precond3.LinComb(1.0, Mass.Data, sigma_*dt_, Laplace.Data, dt_*rho_, Volume_stab.Data);
//        Precond4.LinComb(1.0, Mass.Data, epsilon_*epsilon_, Laplace.Data, epsilon_*epsilon_*rho_, Volume_stab.Data);
        block_pc_.GetPC1().Reset(Precond3);
        block_pc_.GetPC2().Reset(Precond4);
    }


// Implement Implicit Euler method
    void CahnHilliardNarrowBandStblP1CL::InitStep1 (VectorCL& rhs3, VectorCL& rhs4, VectorCL& rhs5, double new_t)
    {
        // ScopeTimerCL timer( "CahnHilliardNarrowBandStblP1CL::InitStep1");
        std::cout << "CahnHilliardNarrowBandStblP1CL::InitStep1:\n";

        ic.t= new_t;
        imu.t= new_t;
        is.t= new_t;
        ienergy.t= new_t;

        dt_= new_t - oldt_;
        idx.CreateNumbering( oldidx_.TriangLevel(), MG_, &lset_vd_, &lsetbnd_, width_); //WARNING// InitTimeStep deletes oldidx_ and swaps idx and oldidx_.
        std::cout << "new NumUnknowns: " << idx.NumUnknowns() << std::endl;
        ic.SetIdx( &idx);
        imu.SetIdx( &idx);
        is.SetIdx( &idx);
        iface.SetIdx( &idx);
        ienergy.SetIdx( &idx);

        //NarrowBand tetrahedras in VTK
        for (int i=0; i<idx.NumUnknowns();i++)
        {
            iface.Data[i]=1.;
        }

        full_idx.CreateNumbering( idx.TriangLevel(), MG_);
        std::cout << " full NumUnknowns: " << full_idx.NumUnknowns() << std::endl;
        ext1.SetIdx(&full_idx);
        ext3.SetIdx(&full_idx);

        conc_extrapol.SetIdx(&full_idx);
        species_extrapol.SetIdx(&full_idx);


        VecDescCL vd_timeder( &idx),    // right-hand sides from integrals over the old/new interface
                vd_oldtimeder( &idx),
                vd_sp_timeder( &idx),
                vd_load3( &idx),
                vd_load4( &idx),
                vd_load5( &idx),
                vd_oldres( &idx),
                vd_oldload3( &idx),
                vd_oldload4( &idx),
                well_potential( &idx),//double-well potential
                vd_well( &idx),
                vd_oldic( &oldidx_),  // the initial data.
                vd_oldis( &oldidx_);  // the initial data.


        vd_oldic.Data= oldic_;
        vd_oldic.t= oldt_;

        vd_oldis.Data= oldis_;
        vd_oldis.t= oldt_;

        DROPS::ExtendP2( MG_, vd_oldic, ext1);
        VecDescCL old_conc( &idx);//extension of the old concentration to the current NarrowBand
        Restrict( MG_, ext1, old_conc);
        conc_extrapol.Data=1.0*ext1.Data;

        DROPS::ExtendP2( MG_, vd_oldis, ext3);
        VecDescCL old_species( &idx);//extension of the old species to the current NarrowBand
        Restrict( MG_, ext3, old_species);
        species_extrapol.Data=1.0*ext3.Data;

        //new
        ic.Data=old_conc.Data;
        is.Data=old_species.Data;

        TetraAccumulatorTupleCL accus;
        InterfaceCommonDataP1CL cdata( lset_vd_, lsetbnd_);
        accus.push_back( &cdata);

        InterfaceVectorAccuCL<LocalMatVecP1CL<LocalInterfaceMassP1CL>, InterfaceCommonDataP1CL>
                mass_accu( &vd_timeder, LocalMatVecP1CL<LocalInterfaceMassP1CL>( LocalInterfaceMassP1CL(), &old_conc),
                           cdata, "concentration");
        accus.push_back( &mass_accu);

        for (int i=0; i<well_potential.Data.size();i++)
        {
            well_potential.Data[i]=Potential_prime_function(old_conc.Data[i]);//BDF1
        }

        InterfaceVectorAccuCL<LocalMatVecP1CL<LocalInterfaceMassP1CL>, InterfaceCommonDataP1CL>
                well_accu( &vd_well, LocalMatVecP1CL<LocalInterfaceMassP1CL>( LocalInterfaceMassP1CL(), &well_potential),
                           cdata, "chemical_potential");
        accus.push_back( &well_accu);

        InterfaceVectorAccuCL<LocalMatVecP1CL<LocalInterfaceMassP1CL>, InterfaceCommonDataP1CL>
                sp_accu( &vd_sp_timeder, LocalMatVecP1CL<LocalInterfaceMassP1CL>( LocalInterfaceMassP1CL(), &old_species),
                           cdata, "species");
        accus.push_back( &sp_accu);

        if (rhs_fun3_)
            accus.push_back_acquire( new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>( &vd_load3,
                                     LocalVectorP1CL( rhs_fun3_, new_t), cdata, "load3"));

        if (rhs_fun4_)
            accus.push_back_acquire( new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>( &vd_load4,
                                     LocalVectorP1CL( rhs_fun4_, new_t), cdata, "load4"));

        if (rhs_fun5_)
            accus.push_back_acquire( new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>( &vd_load5,
                                    LocalVectorP1CL( rhs_fun5_, new_t), cdata, "load5"));

        accumulate( accus, MG_, idx.TriangLevel(), idx.GetBndInfo());

        rhs3 = VectorCL( (1.0/1.0)*vd_timeder.Data + (1.0*dt_)*vd_load3.Data);
        //rhs3 = VectorCL( (1.0/dt_)*vd_timeder.Data + 1.0*vd_load3.Data);
        rhs4 = VectorCL( (-1.)*vd_well.Data + (S_)*vd_timeder.Data + vd_load4.Data);
        rhs5 = VectorCL( (1.0/1.0)*vd_sp_timeder.Data + (1.0*dt_)*vd_load5.Data);
        //        rhs3 = VectorCL( 1.0*(vd_timeder.Data + dt_*vd_load3.Data));
//        rhs4 = VectorCL( (-1.)*vd_well.Data + (S_)*vd_timeder.Data + vd_load4.Data);



       /* WriteToFile( vd_oldic.Data, "vd_oldic.txt", "vd_oldic");
        WriteToFile( rhsext1.Data, "rhsext1.txt", "rhsext1");


        WriteToFile( vd_well.Data, "vd_well.txt", "vd_well");*/

        return;

        /*if(new_t/dt_==1)//used in InitStep 3.
        {
            rhsext1.Reset();
            rhsext1.SetIdx(&full_idx);
            rhsext1.Data=(std::sqrt(dt_))*rhsext.Data;// Since we choose a smaller time step dt*0.1
        }*/
    }

// Implement BDF2 method
//In the first step, we use the backward Euler method
    void CahnHilliardNarrowBandStblP1CL::InitStep2 (VectorCL& rhs3, VectorCL& rhs4, double new_t) {

        ScopeTimerCL timer("CahnHilliardNarrowBandStblP1CL::InitStep2");
        std::cout << "CahnHilliardNarrowBandStblP1CL::InitStep2:\n";

        ic.t = new_t;
        imu.t= new_t;
        dt_ = new_t - oldt_;
        idx.CreateNumbering(oldidx_.TriangLevel(), MG_, &lset_vd_, &lsetbnd_ ,width_); // InitTimeStep deletes oldidx_ and swaps idx and oldidx_.
        std::cout << "new NumUnknowns: " << idx.NumUnknowns() << std::endl;
        ic.SetIdx(&idx);
        imu.SetIdx(&idx);
        iface.SetIdx(&idx);//NarrowBand tetrahedras in VTK
        for (int i = 0; i < idx.NumUnknowns(); i++) {
            iface.Data[i] = 1.;
        }

        VecDescCL vd_timeder(&idx),    // right-hand sides from integrals over the old/new interface
                vd_oldtimeder(&idx),
                vd_load3(&idx),
                vd_load4(&idx),
                vd_oldres(&idx),
                vd_oldload3(&idx),
                vd_oldload4(&idx),
                well_potential(&idx),//double-well potential
                vd_well(&idx),
                vd_conc_extr(&idx),
                vd_oldic(&oldidx_),  // the previous timestep data.
                vd_oldoldic(&oldoldidx_);  // the data previous to the previous data.


        full_idx.CreateNumbering(idx.TriangLevel(), MG_);
        std::cout << " full NumUnknowns: " << full_idx.NumUnknowns() << std::endl;
        ext1.SetIdx(&full_idx);
        ext2.SetIdx(&full_idx);
        conc_extrapol.SetIdx(&full_idx);


        VecDescCL old_conc(&idx);//extension of the old concentration to the current NarrowBand
        VecDescCL oldold_conc(&idx);//extension of the old old concentration to the current NarrowBand

        vd_oldic.Data = oldic_;
        vd_oldic.t = oldt_;
        DROPS::ExtendP2(MG_, vd_oldic, ext1);

        Restrict(MG_, ext1, old_conc);

        vd_oldoldic.Data = oldoldic_;
        vd_oldoldic.t = oldt_;
        DROPS::ExtendP2(MG_, vd_oldoldic, ext2);
        Restrict(MG_, ext2, oldold_conc);

        //new
        ic.Data=old_conc.Data;

        //warning
      /*  old_conc.Data=oldic_;
        oldold_conc.Data=oldoldic_;*/

        VecDescCL bdf2(&idx);
        VecDescCL conc_ext(&idx);
         {
            bdf2.Data = 2.0 * old_conc.Data - 0.5 * oldold_conc.Data;
            conc_ext.Data =  2.0 * old_conc.Data - 1.0 * oldold_conc.Data;
            conc_extrapol.Data=2.0*ext1.Data-1.0*ext2.Data;

            for (int i=0; i<well_potential.Data.size();i++)
            {
                well_potential.Data[i]=2.0*Potential_prime_function(old_conc.Data[i]) - 1.0*Potential_prime_function(oldold_conc.Data[i]);//interpolation
            }
        }
        TetraAccumulatorTupleCL accus;
        InterfaceCommonDataP1CL cdata( lset_vd_, lsetbnd_);
        accus.push_back( &cdata);


        InterfaceVectorAccuCL<LocalMatVecP1CL<LocalInterfaceMassP1CL>, InterfaceCommonDataP1CL>
                deriv_accu( &vd_timeder, LocalMatVecP1CL<LocalInterfaceMassP1CL>( LocalInterfaceMassP1CL(), &bdf2),
                           cdata, "time-der");
        accus.push_back( &deriv_accu);

        InterfaceVectorAccuCL<LocalMatVecP1CL<LocalInterfaceMassP1CL>, InterfaceCommonDataP1CL>
                mass_accu( &vd_conc_extr, LocalMatVecP1CL<LocalInterfaceMassP1CL>( LocalInterfaceMassP1CL(), &conc_ext),
                           cdata, "conc_extrap");
        accus.push_back( &mass_accu);

        InterfaceVectorAccuCL<LocalMatVecP1CL<LocalInterfaceMassP1CL>, InterfaceCommonDataP1CL>
                well_accu( &vd_well, LocalMatVecP1CL<LocalInterfaceMassP1CL>( LocalInterfaceMassP1CL(), &well_potential),
                           cdata, "chemical_potential");
        accus.push_back( &well_accu);

        if (rhs_fun3_)
            accus.push_back_acquire( new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>( &vd_load3,
                                                                                                          LocalVectorP1CL( rhs_fun3_, new_t), cdata, "load3"));

        if (rhs_fun4_)
            accus.push_back_acquire( new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>( &vd_load4,
                                                                                                          LocalVectorP1CL( rhs_fun4_, new_t), cdata, "load4"));

        accumulate( accus, MG_, idx.TriangLevel(), idx.GetBndInfo());

        rhs3 = VectorCL( 1.0*(vd_timeder.Data + dt_*vd_load3.Data));
        rhs4 = VectorCL( (-1.)*vd_well.Data + (S_)*vd_conc_extr.Data + vd_load4.Data);



        ////////////////////////////from surfactant////////
       /* if(new_t/dt_==1)
        {
            std::cout<<"test---1"<<std::endl;
            VecDescCL rhs1( &idx_c);
            //Restrict( MG_, rhsext, rhs1);
            MassMultiply(MG_,rhsext,rhs1,lset_vd_, lsetbnd_);
            //std::cout<<"test---1"<<std::endl;
            rhs1_.resize( rhs1.Data.size());
            rhs1_= rhs1.Data;
        }
        else
        {
            VecDescCL rhs1( &idx);//present time step
            //Restrict( MG_, rhsext, rhs1);
            MassMultiply(MG_,rhsext,rhs1,lset_vd_, lsetbnd_);

            //    fulltransport_->DoStep( rhsext1.Data, make_P2Eval( MG_, Bnd_v_, *v_));
            VecDescCL rhs2( &idx);//last time step
            //Restrict( MG_, rhsext1, rhs2);
            MassMultiply(MG_,rhsext1,rhs2,lset_vd_, lsetbnd_);

            rhs1_.resize( rhs1.Data.size());
            rhs1_= rhs1.Data*2.0-rhs2.Data*0.5;
        }

        rhsext1.Reset();
        rhsext1.SetIdx(&full_idx);
        rhsext1.Data=rhsext.Data;

        //  WriteToFile( rhsext1.Data, "Extended01.txt", "mass");
*/
    }

// Implement BDF2 method
// In the first step, we solve the problem with much smaller time step
    void CahnHilliardNarrowBandStblP1CL::InitStep3 (double new_t)
    {
        /*
        // ScopeTimerCL timer( "SurfactantNarrowBandStblP1CL::InitStep");
        std::cout << "SurfactantNarrowBandStblP1CL::InitStep:\n";

        ic.t= new_t;
        dt_= ic.t - oldt_;

        // idx.GetXidx().SetBound( width_); //transfer the width_ to CreatNumbering
        idx.CreateNumbering( oldidx_.TriangLevel(), MG_, &lset_vd_, &lsetbnd_,width_); // InitTimeStep deletes oldidx_ and swaps idx and oldidx_.
        std::cout << "new NumUnknowns: " << idx.NumUnknowns();

        full_idx.CreateNumbering( idx.TriangLevel(), MG_);
        std::cout << " full NumUnknowns: " << full_idx.NumUnknowns() << std::endl;
        DROPS::VecDescCL rhsext( &full_idx);


        VecDescCL rhs( &oldidx_);
        rhs.Data= (1./dt_)*oldic_;
        DROPS::ExtendP2( MG_, rhs, rhsext);
        //  WriteToFile( rhsext.Data, "Extended00.txt", "mass");



        if(new_t/dt_==1)
        {
            std::cout<<"Error occur: No solution provided on the first step, please use Initstep 2 instead"<<std::endl;
        }
        else
        {
            VecDescCL rhs1( &idx);//present time step
            //Restrict( MG_, rhsext, rhs1);
            MassMultiply(MG_,rhsext,rhs1,lset_vd_, lsetbnd_);

            //    fulltransport_->DoStep( rhsext1.Data, make_P2Eval( MG_, Bnd_v_, *v_));
            VecDescCL rhs2( &idx);//last time step
            //Restrict( MG_, rhsext1, rhs2);
            MassMultiply(MG_,rhsext1,rhs2,lset_vd_, lsetbnd_);

            rhs1_.resize( rhs1.Data.size());
            rhs1_= rhs1.Data*2.0-rhs2.Data*0.5;
        }
        rhsext1.Reset();
        rhsext1.SetIdx(&full_idx);
        rhsext1.Data=rhsext.Data;

        //  WriteToFile( rhsext1.Data, "Extended01.txt", "mass");
         */

    }
//Implicit Euler method
    void CahnHilliardNarrowBandStblP1CL::DoStep1 (VectorCL& rhs3, VectorCL& rhs4, VectorCL& rhs5)
    {
        Update();

        A_.LinComb(0.0, Mass.Data, 1.0*sigma_*dt_, LaplaceM.Data,  1.0*rho_*dt_, Volume_stab.Data);
        B_.LinComb(1.0/(1.), Mass.Data, 1.0*dt_, Conv.Data);
        C_.LinComb(-1.0*epsilon_, Massrho.Data, 0.,   Mass.Data);
        D_.LinComb(1.0*S_, Mass.Data, 1.0*epsilon_*epsilon_,  Laplace.Data,  1.0*rho_*epsilon_*epsilon_, Volume_stab.Data);


        K_.LinComb(1.0, Mass.Data, 0.1*dt_, Laplace.Data, 1.0*rho_*dt_, Volume_stab.Data);

//        double c=1.0;//BDF1
//
//        A_.LinComb(sigma_*dt_, LaplaceM.Data, rho_*dt_, Volume_stab.Data);
//        B_.LinComb(c, Mass.Data, dt_, Conv.Data);//0*dt_ , Massd.Data,
//
//                C_.LinComb(0.,   Laplace.Data, -1.0, Mass.Data);
//        D_.LinComb(epsilon_*epsilon_,  Laplace.Data,
//                //   -1., Gprimeprime.Data,
//                   S_, Mass.Data,
//                   rho_*epsilon_*epsilon_, Volume_stab.Data);

        //WriteToFile( L_, "Matrix.txt", "system");
        //A_.LinComb(theta_, Mass.Data, dt_ * theta_, Laplace.Data, dt_ * theta_, Massd.Data, dt_ * theta_, Conv.Data);
        //D_.LinComb(1.0, A_, 0.0, Mass.Data);
        std::cout << "      Before solve: res3 = " << norm(A_ * imu.Data  +B_ * ic.Data- rhs3) << std::endl;
        std::cout << "      Before solve: res4 = " << norm(C_ * imu.Data  +D_ * ic.Data- rhs4) << std::endl;
        {
            ScopeTimerCL timer("CahnHilliardNarrowBandStblP1CL::DoStep1: Solve");
            block_gm_.Solve(A_, B_, C_, D_, imu.Data, ic.Data, rhs3, rhs4, imu.RowIdx->GetEx(),ic.RowIdx->GetEx());

            rhs5 += (0.1*dt_)*(LaplaceNon.Data*ic.Data);
            PCGSolver3_.Solve(K_, is.Data, rhs5, is.RowIdx->GetEx());

            for (int i=0; i<ienergy.Data.size();i++)
            {
                ienergy.Data[i]=Potential_function(ic.Data[i]);
            }

            /*ScopeTimerCL timer( "SurfactantcGP1CL::DoStep: Solve");
            MatrixCL L;
            L.LinComb( 1.0, Mass.Data, dt_*sigma_, Laplace.Data, 1*rho_*dt_, Volume_stab.Data, 0*dt_, Massd.Data, dt_, Conv.Data);
            gm_.Solve( L, ic.Data, rhs3, ic.RowIdx->GetEx());*/

            std::cout << "      Iterations of inner 3 and 4 solver on the last outer step"  << ": " << PCGSolver3_.GetIter() << "\t" << PCGSolver4_.GetIter() << '\n';
            std::cout << "      After  solve: res3 = " << norm(A_ * imu.Data  +B_ * ic.Data- rhs3) << std::endl;
            std::cout << "      After  solve: res4 = " << norm(C_ * imu.Data  +D_ * ic.Data- rhs4) << std::endl;
        }
        //imu.Data = (epsilon_)*imu.Data;

        std::cout << "CahnHilliardNarrowBandStblP1CL::DoStep1: res = " << block_gm_.GetResid() << ", iter = " << block_gm_.GetIter() << std::endl;
    }

//BDF2 method
    void CahnHilliardNarrowBandStblP1CL::DoStep2 (VectorCL& rhs3, VectorCL& rhs4)
    {
        Update();

        /*A_.LinComb(sigma_, LaplaceM.Data, rho_, Volume_stab.Data);
        B_.LinComb(3./(2.*dt_), Mass.Data, 0., Conv.Data);
        //0*dt_ , Massd.Data,

        C_.LinComb(0.,   Laplace.Data, -epsilon_, Mass.Data);
        D_.LinComb(epsilon_*epsilon_,  Laplace.Data,
                //   -1., Gprimeprime.Data,
                   S_, Mass.Data,
                   rho_*epsilon_*epsilon_, Volume_stab.Data);
        rhs3=(1./dt_)*rhs3;*/

        A_.LinComb(sigma_*dt_, LaplaceM.Data, rho_*dt_, Volume_stab.Data);
        B_.LinComb(3./2., Mass.Data, dt_, Conv.Data);
                //0*dt_ , Massd.Data,

        C_.LinComb(0.,   Laplace.Data, -1.0, Mass.Data);
        D_.LinComb(epsilon_*epsilon_,  Laplace.Data,
                //   -1., Gprimeprime.Data,
                   S_, Mass.Data,
                   rho_*epsilon_*epsilon_, Volume_stab.Data);

        //A_.LinComb(theta_, Mass.Data, dt_ * theta_, Laplace.Data, dt_ * theta_, Massd.Data, dt_ * theta_, Conv.Data);
        //D_.LinComb(1.0, A_, 0.0, Mass.Data);
        std::cout << "      Before solve: res3 = " << norm(A_ * imu.Data  +B_ * ic.Data- rhs3) << std::endl;
        std::cout << "      Before solve: res4 = " << norm(C_ * imu.Data  +D_ * ic.Data- rhs4) << std::endl;

        {
            ScopeTimerCL timer("CahnHilliardP1BaseCL::DoStep: Solve");
            block_gm_.Solve(A_, B_, C_, D_, imu.Data, ic.Data, rhs3, rhs4, imu.RowIdx->GetEx(),ic.RowIdx->GetEx());

            std::cout << "      Iterations of inner 3 and 4 solver on the last outer step"  << ": " << PCGSolver3_.GetIter() << "\t" << PCGSolver4_.GetIter() << '\n';
            std::cout << "      After  solve: res3 = " << norm(A_ * imu.Data  +B_ * ic.Data- rhs3) << std::endl;
            std::cout << "      After  solve: res4 = " << norm(C_ * imu.Data  +D_ * ic.Data- rhs4) << std::endl;
        }
        //imu.Data = (epsilon_)*imu.Data;
        std::cout << "CahnHilliardP1BaseCL::DoStep: res = " << block_gm_.GetResid() << ", iter = " << block_gm_.GetIter() << std::endl;
    }

    void CahnHilliardNarrowBandStblP1CL::CommitStep ()
    {
        full_idx.DeleteNumbering( MG_);
    }


//only back forward Euler
    void CahnHilliardNarrowBandStblP1CL::DoStep0 (double new_t)
    {
        ScopeTimerCL timer( "CahnHilliardNarrowBandStblP1CL::::DoStep0");
        VectorCL rhs3, rhs4,rhs5;
        InitStep1( rhs3, rhs4, rhs5, new_t);
        DoStep1(rhs3, rhs4, rhs5 );/**/
        CommitStep();
        //   WriteToFile( temp_ic.Data, "Extended.txt", "mass");

    }

//BDF2 method
    void CahnHilliardNarrowBandStblP1CL::DoStep (double new_t)
    {
         ScopeTimerCL timer( "CahnHilliardNarrowBandP1CL::::DoStep");
         VectorCL rhs3, rhs4,rhs5;

        /*InitStep1( rhs3, rhs4, new_t);
        std::cout<<"BDF1 DoStep"<<std::endl;//first step~
        DoStep1(rhs3, rhs4);*/

         if(new_t/(new_t-oldt_)==1)
         {
             InitStep1( rhs3, rhs4, rhs5, new_t);
             std::cout<<"BDF1 DoStep"<<std::endl;//first step~
             DoStep1(rhs3, rhs4,rhs5);
         }
         else {
             InitStep2(rhs3, rhs4, new_t);
             std::cout<<"BDF2 DoStep"<<std::endl;//other steps
             DoStep2(rhs3, rhs4);
         }

        CommitStep();

//        InitStep3(new_t);
//        DoStep2();
//        CommitStep();/*In the first step, solve the equation by some smaller timestep size-- */

        //   WriteToFile( temp_ic.Data, "Extended.txt", "mass");

    }


    /////////////////////////SurfactantCGP1CL
    void SurfactantcGP1CL::Update()
{
    // ScopeTimerCL timer("SurfactantcGP1CL::Update");
    // std::cout << "SurfactantcGP1CL::Update:\n";

    IdxDescCL* cidx= ic.RowIdx;
    M.Data.clear();
    M.SetIdx(cidx, cidx);
    A.Data.clear();
    A.SetIdx(cidx, cidx);
    C.Data.clear();
    C.SetIdx(cidx, cidx);
    Md.Data.clear();
    Md.SetIdx(cidx, cidx);

    TetraAccumulatorTupleCL accus;
    InterfaceCommonDataP1CL cdata(lset_vd_, lsetbnd_);
    accus.push_back(&cdata);
    InterfaceMatrixAccuCL<LocalInterfaceMassP1CL, InterfaceCommonDataP1CL> mass_accu(&M, LocalInterfaceMassP1CL(), cdata, "mass");
    accus.push_back(&mass_accu);
    InterfaceMatrixAccuCL<LocalLaplaceBeltramiP1CL, InterfaceCommonDataP1CL> lb_accu(&A, LocalLaplaceBeltramiP1CL(D_), cdata, "Laplace-Beltrami");
    accus.push_back(&lb_accu);
    accus.push_back_acquire(make_wind_dependent_matrixP1_accu<LocalInterfaceConvectionP1CL>(&C,  cdata,  make_P2Eval(MG_, Bnd_v_, *v_), "convection"));
    accus.push_back_acquire(make_wind_dependent_matrixP1_accu<LocalInterfaceMassDivP1CL>   (&Md, cdata,  make_P2Eval(MG_, Bnd_v_, *v_), "massdiv"));

    if (theta_ != 1.0) {
        M2.Data.clear();
        M2.SetIdx(cidx, cidx);
        InterfaceCommonDataP1CL* oldcdata= new InterfaceCommonDataP1CL(oldls_, lsetbnd_);
        accus.push_back_acquire(oldcdata);
        accus.push_back_acquire(new InterfaceMatrixAccuCL<LocalInterfaceMassP1CL, InterfaceCommonDataP1CL>(&M2, LocalInterfaceMassP1CL(), *oldcdata, "old mass"));
    }
    accumulate(accus, MG_, cidx->TriangLevel(), cidx->GetBndInfo());

//     WriteToFile(M.Data, "cGcGM.txt", "mass");
//     WriteToFile(A.Data, "cGcGA.txt", "Laplace-Beltrami");
//     WriteToFile(C.Data, "cGcGC.txt", "material derivative");
//     WriteToFile(Md.Data,"cGcGMd.txt","mass-div");

    // std::cout << "SurfactantP1CL::Update: Finished\n";
}

VectorCL SurfactantcGP1CL::InitStep (double new_t)
{
    // ScopeTimerCL timer("SurfactantcGP1CL::InitStep");
    // std::cout << "SurfactantcGP1CL::InitStep:\n";

    ic.t= new_t;
    dt_= new_t - oldt_;
    idx.CreateNumbering(oldidx_.TriangLevel(), MG_, &lset_vd_, &lsetbnd_); // InitTimeStep deletes oldidx_ and swaps idx and oldidx_.
    std::cout << "new NumUnknowns: " << idx.NumUnknowns() << std::endl;
    ic.SetIdx(&idx);

    VecDescCL vd_timeder(&idx),    // right-hand sides from integrals over the old/new interface
              vd_oldtimeder(&idx),
              vd_load(&idx),
              vd_oldres(&idx),
              vd_oldload(&idx),
              vd_oldic(&oldidx_);  // the initial data.
    vd_oldic.Data= oldic_;
    vd_oldic.t= oldt_;

    TetraAccumulatorTupleCL accus;
    InterfaceCommonDataP1CL cdata(lset_vd_, lsetbnd_);
    accus.push_back(&cdata);
    InterfaceVectorAccuCL<LocalMatVecP1CL<LocalInterfaceMassP1CL>, InterfaceCommonDataP1CL> mass_accu(&vd_timeder,
        LocalMatVecP1CL<LocalInterfaceMassP1CL>(LocalInterfaceMassP1CL(), &vd_oldic), cdata, "mixed-mass");
    accus.push_back(&mass_accu);

    if (rhs_fun_)
        accus.push_back_acquire(new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>(&vd_load, LocalVectorP1CL(rhs_fun_, new_t), cdata, "load"));

    if (theta_ == 1.0) {
        accumulate(accus, MG_, idx.TriangLevel(), idx.GetBndInfo());
        return VectorCL(theta_*(vd_timeder.Data + dt_*vd_load.Data));
    }

    InterfaceCommonDataP1CL oldcdata(oldls_, lsetbnd_);
    accus.push_back(&oldcdata);

    if (rhs_fun_)
        accus.push_back_acquire(new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>(&vd_oldload, LocalVectorP1CL(rhs_fun_, oldt_), oldcdata, "load on old iface"));

    InterfaceVectorAccuCL<LocalMatVecP1CL<LocalInterfaceMassP1CL>, InterfaceCommonDataP1CL> old_mass_accu(&vd_oldtimeder,
        LocalMatVecP1CL<LocalInterfaceMassP1CL>(LocalInterfaceMassP1CL(), &vd_oldic), oldcdata, "mixed-mass on old iface");
    accus.push_back(&old_mass_accu);
    InterfaceVectorAccuCL<LocalMatVecP1CL<LocalLaplaceBeltramiP1CL>, InterfaceCommonDataP1CL> old_lb_accu(&vd_oldres,
        LocalMatVecP1CL<LocalLaplaceBeltramiP1CL>(LocalLaplaceBeltramiP1CL(D_), &vd_oldic), oldcdata, "Laplace-Beltrami on old iface");
    accus.push_back(&old_lb_accu);
    accus.push_back_acquire(make_wind_dependent_vectorP1_accu<LocalInterfaceConvectionP1CL>(&vd_oldres, &vd_oldic,  oldcdata,  make_P2Eval(MG_, Bnd_v_, oldv_), "convection on old iface"));
    accus.push_back_acquire(make_wind_dependent_vectorP1_accu<LocalInterfaceMassDivP1CL>   (&vd_oldres, &vd_oldic,  oldcdata,  make_P2Eval(MG_, Bnd_v_, oldv_), "mass-div on old iface"));

    accumulate(accus, MG_, idx.TriangLevel(), idx.GetBndInfo());
    return VectorCL(theta_*vd_timeder.Data + (1. - theta_)*vd_oldtimeder.Data
                   + dt_*(theta_*vd_load.Data + (1. - theta_)*(vd_oldload.Data - vd_oldres.Data)));
}

void SurfactantcGP1CL::DoStep (const VectorCL& rhs)
{
    Update();

    if (theta_ == 1.)
        L_.LinComb(theta_, M.Data, dt_*theta_, A.Data, dt_*theta_, Md.Data, dt_*theta_, C.Data);
    else {
        MatrixCL m;
        m.LinComb(theta_, M.Data, dt_*theta_, A.Data, dt_*theta_, Md.Data, dt_*theta_, C.Data);
        L_.LinComb(1., m, 1. - theta_, M2.Data);
    }
    std::cout << "Before solve: res = " << norm(L_*ic.Data - rhs) << std::endl;
    {
        ScopeTimerCL timer("SurfactantcGP1CL::DoStep: Solve");
        gm_.Solve(L_, ic.Data, rhs, ic.RowIdx->GetEx());
    }
    std::cout << "SurfactantcGP1CL::DoStep: res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
}

void SurfactantcGP1CL::CommitStep ()
{
    return;
}

void SurfactantcGP1CL::DoStep (double new_t)
{
    ScopeTimerCL timer("SurfactantcGP1CL::DoStep");

    VectorCL rhs(InitStep(new_t));
    DoStep(rhs);
    CommitStep();
}

void
InterfaceP1RepairCL::pre_refine ()
{
    DROPS::NoBndDataCL<> dummy;
    p1repair_= std::unique_ptr<RepairP1CL<double, NoBndDataCL>::type >(
        new RepairP1CL<double, NoBndDataCL>::type(mg_, fullu_, dummy));
}

void
InterfaceP1RepairCL::post_refine ()
{
    VecDescCL loc_u;
    IdxDescCL loc_idx(P1_FE);

    loc_idx.CreateNumbering(fullp1idx_.TriangLevel(), mg_);
//    loc_idx.CreateNumbering(mg_.GetLastLevel(), mg_);
    loc_u.SetIdx(&loc_idx);
    DROPS::NoBndDataCL<> dummy;

//    P1EvalCL<double, DROPS::NoBndDataCL<>, const VecDescCL> oldsol(&fullu_, &dummy, &mg_);
//    P1EvalCL<double, DROPS::NoBndDataCL<>,       VecDescCL>    sol(&loc_u , &dummy, &mg_);
//    Interpolate(sol, oldsol);

    p1repair_->repair(loc_u);

    fullu_.Clear(fullu_.t);
    fullp1idx_.swap(loc_idx);
    loc_idx.DeleteNumbering(mg_);
    fullu_.SetIdx(&fullp1idx_);
    fullu_.Data= loc_u.Data;
}

void
InterfaceP1RepairCL::pre_refine_sequence ()
{
    fullp1idx_.CreateNumbering(u_.RowIdx->TriangLevel(), mg_);
    fullu_.SetIdx(&fullp1idx_);
    Extend(mg_, u_, fullu_);
}

void
InterfaceP1RepairCL::post_refine_sequence ()
{
    u_.RowIdx->DeleteNumbering( mg_);

    if (width_>0)
         u_.RowIdx->CreateNumbering( fullp1idx_.TriangLevel(), mg_, &lset_vd_, &lset_bnd_, width_);
    else u_.RowIdx->CreateNumbering( fullp1idx_.TriangLevel(), mg_, &lset_vd_, &lset_bnd_);

    u_.SetIdx( u_.RowIdx);

    Restrict(mg_, fullu_, u_);

    fullp1idx_.DeleteNumbering(mg_);
    fullu_.Clear(u_.t);
}
#ifndef _PAR
void
Ensight6IfaceScalarCL::put (Ensight6OutCL& cf) const
{
    IdxDescCL p1idx;
    p1idx.CreateNumbering(u_.RowIdx->TriangLevel(), mg_);
    VecDescCL p1u(&p1idx);
    Extend(mg_, u_, p1u);
    BndDataCL<> bnd(0);

    cf.putScalar(make_P1Eval(mg_, bnd, p1u), varName());

    p1idx.DeleteNumbering(mg_);
}
#endif

void
VTKIfaceScalarCL::put (VTKOutCL& cf) const
{
    IdxDescCL fullidx(P1_FE);
    if (u_.RowIdx->GetFE() == P2IF_FE)
        fullidx.SetFE(P2_FE);
    fullidx.CreateNumbering(u_.RowIdx->TriangLevel(), mg_);
    std::cout<<"vtk scalar unknowns: "<<fullidx.NumUnknowns()<<std::endl;
    VecDescCL uext(&fullidx);
    Extend(mg_, u_, uext);

    if (u_.RowIdx->GetFE() == P1IF_FE)
        cf.PutScalar(make_P1Eval(mg_, BndData_, uext), varName());
    else if (u_.RowIdx->GetFE() == P2IF_FE)
        cf.PutScalar(make_P2Eval(mg_, BndData_, uext), varName());

    fullidx.DeleteNumbering(mg_);
}

void VTKIfaceVectorCL::put (VTKOutCL& cf) const
{
    IdxDescCL pidx;
    if(P2_)
      pidx.SetFE(vecP2_FE);
    else
      pidx.SetFE(vecP1_FE);
    pidx.CreateNumbering(u_.RowIdx->TriangLevel(), mg_);
    std::cout<<"vtk vector unknowns: "<<pidx.NumUnknowns()<<std::endl;
    VecDescCL pu(&pidx);
    Extend(mg_, u_, pu);

    P2EvalCL<SVectorCL<3>, const BndDataCL<Point3DCL>, const VecDescBaseCL<VectorCL> > eval2(&pu, &BndData_, &mg_);
    P1EvalCL<SVectorCL<3>, const BndDataCL<Point3DCL>, const VecDescBaseCL<VectorCL> > eval1(&pu, &BndData_, &mg_);
    if(P2_)
      cf.PutVector(eval2, varName());
    else
      cf.PutVector(eval1, varName());

    pidx.DeleteNumbering(mg_);
}


/// ==Space-Time-methods==
/// \todo Maybe introduce a new file

///\brief Bilinear space-time FE-function on a single TetraPrismCL.
LocalSTP1P1CL<> STP1P1DiscCL::ref_val[8];
LocalSTP1CL<Point4DCL> STP1P1DiscCL::ref_grad[8];
LocalSTP1CL<Point3DCL> STP1P1DiscCL::ref_gradx[8];

void STP1P1DiscCL::StaticInit()
{
    for (Uint i= 0; i < 4; ++i) {
        ref_val[i    ].at_t0()[i]= 1.;
        ref_val[i + 4].at_t1()[i]= 1.;

        const Point3DCL p1grad(FE_P1CL::DHRef(i));
        const Point4DCL& p1grad4d= MakePoint4D(p1grad[0], p1grad[1], p1grad[2], 0.);
        for (Uint j= 0; j < 4; ++j)
            ref_grad[i][j]= p1grad4d;
        ref_grad[i][i][3]= -1.;
        ref_grad[i][4][3]= i == 0 ? -1. : 0.;
        ref_grad[i + 4][i][3]=  1.;
        ref_grad[i + 4][4]= p1grad4d;
        ref_grad[i + 4][4][3]= i == 0 ? 1. : 0.;

        for (Uint j= 0; j < 4; ++j)
            ref_gradx[i][j]= p1grad;
        ref_gradx[i + 4][4]= p1grad;
    }
}

namespace {
StaticInitializerCL<STP1P1DiscCL> the_STP1P1DiscCL_initializer_;
} // end of unnamed namespace



void STP1P1IdxDescCL::CreateNumbering (Uint level, MultiGridCL& mg, const VecDescCL& oldls, const VecDescCL& newls, const BndDataCL<>& lsetbnd, double t0, double t1)
{
    TriangLevel_= level;

    IdxDescCL idx_ini(P1IF_FE);
    idx_ini.CreateNumbering(level, mg, &oldls, &lsetbnd);
    NumIniUnknowns_= idx_ini.NumUnknowns();

    IdxDescCL idx_fini(P1IF_FE);
    idx_fini.CreateNumbering(level, mg, &newls, &lsetbnd);
    NumFiniUnknowns_= idx_fini.NumUnknowns();

    NumUnknowns_= NumIniUnknowns() + NumFiniUnknowns();
    IdxT inicounter= 0, finicounter= 0;

    const TetraPrismLatticeCL& lat= TetraPrismLatticeCL::instance(2, 1);
    std::valarray<double> ls_loc(lat.vertex_size());
    LocalP2CL<> oldlocp2_ls, locp2_ls;
    LocalSTP2P1CL<> local_st_lset;
    SPatchCL<4> patch;
    QuadDomainCodim1CL<4> qdom;
    std::valarray<double> shape; // STP1P1-shape-function as integrand
    DROPS_FOR_TRIANG_TETRA (mg, level, it) {
        local_st_lset.assign(*it, oldls, newls, lsetbnd);
        evaluate_on_vertexes(local_st_lset, lat, Addr(ls_loc));
        if (!isInCutMesh(ls_loc))
            continue;

        patch.make_patch<MergeCutPolicyCL>(lat, ls_loc);
        make_CompositeQuad5DomainSTCodim1(qdom, patch, TetraPrismCL(*it, t0, t1));
        shape.resize(qdom.vertex_size());
        for (Uint i= 0; i < 8; ++i) {
            const IdxDescCL& idx= i < 4 ? idx0_ : idx1_;
            const IdxDescCL& spatial_idx= i < 4 ? idx_ini : idx_fini;
            UnknownHandleCL& unknowns= const_cast<VertexCL*>(it->GetVertex(i%4))->Unknowns;
            if (unknowns.Exist(idx.GetIdx()))
                continue;

            evaluate_on_vertexes(STP1P1DiscCL::ref_val[i], qdom, Addr(shape));
            if (quad_codim1(shape*shape, qdom) > 0.) {
                unknowns.Prepare(idx.GetIdx());
                if (unknowns.Exist(spatial_idx.GetIdx())) {
                    unknowns(idx.GetIdx())= unknowns(spatial_idx.GetIdx()) + (i < 4 ? 0 : NumIniUnknowns());
                    ++(i < 4 ? inicounter : finicounter);
                }
                else
                    unknowns(idx.GetIdx())= NumUnknowns_++;
            }
        }
    }
    if (inicounter != NumIniUnknowns())
        throw DROPSErrCL("STP1P1IdxDescCL::CreateNumbering: Wrong count of the unknowns on the old interface.\n");
    if (finicounter != NumFiniUnknowns())
        throw DROPSErrCL("STP1P1IdxDescCL::CreateNumbering: Wrong count of the unknowns on the new interface.\n");
    idx_ini.DeleteNumbering(mg);
    idx_fini.DeleteNumbering(mg);
}


void
resize_and_scatter_piecewise_spatial_normal (const SPatchCL<4>& surf, const QuadDomainCodim1CL<4>& qdom, std::valarray<Point3DCL>& spatial_normal)
{
    spatial_normal.resize(qdom.vertex_size());
    if (spatial_normal.size() == 0)
        return;
    if (surf.normal_empty()) // As qdom has vertexes, the must be facets, i.e. normals.
        throw DROPSErrCL("resize_and_scatter_piecewise_spatial_normal: normals were not precomputed.\n");

    const Uint NodesPerFacet= qdom.vertex_size()/surf.facet_size();
    if (qdom.vertex_size()%surf.facet_size() != 0)
        throw DROPSErrCL("resize_and_scatter_piecewise_spatial_normal: qdom.vertex_size is not a multiple of surf.facet_size.\n");

    const SPatchCL<4>::const_normal_iterator n= surf.normal_begin();
    for (Uint i= 0; i < surf.facet_size(); ++i) {
        const Point3DCL& tmp= MakePoint3D(n[i][0], n[i][1], n[i][2]);
        std::fill_n(&spatial_normal[i*NodesPerFacet], NodesPerFacet, tmp/tmp.norm());
    }
}


void SurfactantSTP1CL::InitStep (double new_t)
{
    // std::cout << "SurfactantSTP1CL::InitStep:\n";
    ic.t= new_t;
    dt_= ic.t - oldt_;

    st_idx_.CreateNumbering(oldidx_.TriangLevel(), MG_, oldls_, lset_vd_, lsetbnd_, oldt_, new_t);
    dim= st_idx_.NumUnknowns() - (cG_in_t_ ? st_idx_.NumIniUnknowns() : 0);
    std::cout << "SurfactantSTP1CL::InitStep: space-time Unknowns: " << st_idx_.NumUnknowns()
              << " ini: " << st_idx_.NumIniUnknowns() << " fini: " << st_idx_.NumFiniUnknowns()
              << " interior: " << st_idx_.NumUnknowns() - st_idx_.NumIniUnknowns() - st_idx_.NumFiniUnknowns()
              << " dimension of the linear system: " << dim << std::endl;
    st_ic_.resize(dim);

    st_oldic_.resize(st_idx_.NumUnknowns());
    // Copy dofs on the old interface from  old solution into st_oldic_.
    DROPS_FOR_TRIANG_VERTEX(MG_, oldidx_.TriangLevel(), it) {
        if (it->Unknowns.Exist(st_idx_.GetIdx(0))) {
            const IdxT dof= it->Unknowns(st_idx_.GetIdx(0));
            if (dof < st_idx_.NumIniUnknowns())
                st_oldic_[dof]= oldic_[dof];
        }
    }
}

void SurfactantSTP1CL::Update_cG()
{
    TetraAccumulatorTupleCL accus;
    STInterfaceCommonDataCL cdata(oldt_, ic.t,  oldls_, lset_vd_, lsetbnd_);
    accus.push_back(&cdata);
    InterfaceMatrixSTP1AccuCL<LocalLaplaceBeltramiSTP1P1CL> lb_accu(&A, &cpl_A_, &st_idx_,
        LocalLaplaceBeltramiSTP1P1CL(D_), cdata, &st_oldic_, "Laplace-Beltrami on ST-iface");
    accus.push_back(&lb_accu);

    cpl_A_.resize(dim);
    cpl_der_.resize(dim);

    VectorCL cpl_new_dummy;
    InterfaceCommonDataP1CL oldspatialcdata(oldls_, lsetbnd_);
    InterfaceMatrixSTP1AccuCL<LocalSpatialInterfaceMassSTP1P1CL> oldmass_accu(&Mold, &cpl_old_, &st_idx_,
        LocalSpatialInterfaceMassSTP1P1CL(oldspatialcdata), cdata, &st_oldic_, "mixed-mass on old iface");
    InterfaceCommonDataP1CL newspatialcdata(lset_vd_, lsetbnd_);
    InterfaceMatrixSTP1AccuCL<LocalSpatialInterfaceMassSTP1P1CL> newmass_accu(&Mnew, /*dummy*/ &cpl_new_dummy, &st_idx_,
        LocalSpatialInterfaceMassSTP1P1CL(newspatialcdata, false), cdata, &st_oldic_, "mixed-mass on new iface");

    if (use_mass_div_) {
        cpl_div_.resize(dim);
        accus.push_back_acquire(make_wind_dependent_matrixSTP1P0_1_accu<LocalMassdivSTP1P1CL>(&Mdiv, &cpl_div_, &st_idx_, cdata, &st_oldic_, make_STP2P1Eval(MG_, Bnd_v_, oldv_, *v_), "mass-div on ST-iface"));
        accus.push_back_acquire(make_wind_dependent_matrixSTP1P0_1_accu<LocalMaterialDerivativeSTP1P1CL>(&Mder, &cpl_der_, &st_idx_, cdata, &st_oldic_, make_STP2P1Eval(MG_, Bnd_v_, oldv_, *v_), "material derivative on ST-iface"));
    }
    else {
        cpl_old_.resize(dim);
        cpl_new_dummy.resize(dim);
        accus.push_back_acquire(make_wind_dependent_local_transpose_matrixSTP1P0_1_accu<LocalMaterialDerivativeSTP1P1CL>(&Mder, &cpl_der_, &st_idx_, cdata, &st_oldic_, make_STP2P1Eval(MG_, Bnd_v_, oldv_, *v_), "material derivative on ST-iface"));
        accus.push_back(&oldspatialcdata);
        accus.push_back(&oldmass_accu);
        accus.push_back(&newspatialcdata);
        accus.push_back(&newmass_accu);
    }

    if (rhs_fun_) {
        load.resize(dim);
        accus.push_back_acquire(new InterfaceVectorSTP1AccuCL<LocalVectorSTP1P1CL>(&load, &st_idx_, LocalVectorSTP1P1CL(rhs_fun_), cdata, /* cG_in_t*/ cG_in_t_, "load on ST-iface"));
    }

    accumulate(accus, MG_, st_idx_.TriangLevel(), idx.GetBndInfo());

    // WriteToFile(cpl_new_dummy, "cpl_new_dummy.txt", "coupling on new interface -- always zero.");
}

void SurfactantSTP1CL::Update_dG()
{
    TetraAccumulatorTupleCL accus;
    InterfaceCommonDataP1CL oldspatialcdata(oldls_, lsetbnd_);
    accus.push_back(&oldspatialcdata);
    STInterfaceCommonDataCL cdata(oldt_, ic.t,  oldls_, lset_vd_, lsetbnd_);
    accus.push_back(&cdata);
    InterfaceMatrixSTP1AccuCL<LocalSpatialInterfaceMassSTP1P1CL> oldmass_accu(&Mold, &st_idx_,
        LocalSpatialInterfaceMassSTP1P1CL(oldspatialcdata), cdata, "mixed-mass on old iface");
    accus.push_back(&oldmass_accu);
    InterfaceMatrixSTP1AccuCL<LocalLaplaceBeltramiSTP1P1CL> lb_accu(&A, &st_idx_,
        LocalLaplaceBeltramiSTP1P1CL(D_), cdata, "Laplace-Beltrami on ST-iface");
    accus.push_back(&lb_accu);

    InterfaceCommonDataP1CL newspatialcdata(lset_vd_, lsetbnd_);
    InterfaceMatrixSTP1AccuCL<LocalSpatialInterfaceMassSTP1P1CL> newmass_accu(&Mnew, &st_idx_,
        LocalSpatialInterfaceMassSTP1P1CL(newspatialcdata, false), cdata, "mixed-mass on new iface");
    if (use_mass_div_) {
        accus.push_back_acquire(make_wind_dependent_matrixSTP1P1_accu<LocalMaterialDerivativeSTP1P1CL>(&Mder, &st_idx_, cdata,  make_STP2P1Eval(MG_, Bnd_v_, oldv_, *v_), "material derivative on ST-iface"));
        accus.push_back_acquire(make_wind_dependent_matrixSTP1P1_accu<LocalMassdivSTP1P1CL>(&Mdiv, &st_idx_, cdata,  make_STP2P1Eval(MG_, Bnd_v_, oldv_, *v_), "mass-div on ST-iface"));
    }
    else {
        accus.push_back_acquire(make_wind_dependent_local_transpose_matrixSTP1P1_accu<LocalMaterialDerivativeSTP1P1CL>(&Mder, &st_idx_, cdata,  make_STP2P1Eval(MG_, Bnd_v_, oldv_, *v_), "material derivative on ST-iface"));
        accus.push_back(&newspatialcdata);
        accus.push_back(&newmass_accu);
    }

    if (rhs_fun_) {
        load.resize(dim);
        accus.push_back_acquire(new InterfaceVectorSTP1AccuCL<LocalVectorSTP1P1CL>(&load, &st_idx_, LocalVectorSTP1P1CL(rhs_fun_), cdata, /* cG_in_t*/ cG_in_t_, "load on ST-iface"));
    }

    accumulate(accus, MG_, st_idx_.TriangLevel(), idx.GetBndInfo());
}

void SurfactantSTP1CL::Update()
{
    ScopeTimerCL timer("SurfactantSTP1CL::Update");
    // std::cout << "SurfactantSTP1CL::Update:\n";
    if (cG_in_t_)
        Update_cG();
    else
        Update_dG();

//     WriteToFile(Mold, "Mold.txt", "mass on old iface");
//     WriteToFile(Mnew, "Mnew.txt", "mass on new iface");
//     WriteToFile(A,    "A.txt",    "Laplace-Beltrami on ST-iface");
//     WriteToFile(Mder, "Mder.txt", "material derivative on ST-iface");
//     WriteToFile(Mdiv, "Mdiv.txt", "mass-div on ST-iface");
//     WriteToFile(load, "load.txt", "load on ST-iface");
//     WriteToFile(cpl_A_,   "cpl_A.txt",   "coupling for Laplace-Beltrami on ST-iface");
//     WriteToFile(cpl_der_, "cpl_der.txt", "coupling for material derivative on ST-iface");
//     WriteToFile(cpl_div_, "cpl_div.txt", "coupling for mass-div on ST-iface");
//     WriteToFile(cpl_old_, "cpl_old.txt", "coupling ini-values on old iface");

    // std::cout << "SurfactantSTP1CL::Update: Finished\n";
}

void SurfactantSTP1CL::DoStep ()
{
    Update();

    MatrixCL L;
    VectorCL rhs(dim);
    if (cG_in_t_) {
        if (use_mass_div_) {
            L.LinComb(1., Mder, 1., Mdiv, 1., A);
            rhs= cpl_der_ + cpl_div_ + cpl_A_;
        }
        else {
            L.LinComb(-1., Mder, 1., A, 1., Mnew);
            rhs= -cpl_der_ + cpl_A_ - cpl_old_;
        }
    }
    else {
        if (use_mass_div_) {
            L.LinComb(1., Mder, 1., Mdiv, 1., A, 1., Mold);
            rhs= Mold*st_oldic_;
        }
        else {
            L.LinComb(-1., Mder, 1., A, 1., Mnew);
            rhs= Mold*st_oldic_;
        }
    }
    if (rhs_fun_ != 0)
        rhs+= load;

    std::cout << "Before solve: res = " << norm(L*st_ic_ - rhs) << std::endl;

    {
        ScopeTimerCL timer("SurfactantSTP1CL::DoStep: Solve");
        gm_.Solve(L, st_ic_, rhs, idx.GetEx());
    }
    std::cout << "SurfactantSTP1CL::DoStep: res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
}

void SurfactantSTP1CL::CommitStep ()
{
    idx.CreateNumbering(oldidx_.TriangLevel(), MG_, &lset_vd_, &lsetbnd_);
    std::cout << "new NumUnknowns at t1: " << idx.NumUnknowns() << std::endl;

    ic.SetIdx(&idx);
    // Copy dofs on the new interface from space-time-solution into ic.
    DROPS_FOR_TRIANG_VERTEX(MG_, oldidx_.TriangLevel(), it) {
        if (it->Unknowns.Exist(st_idx_.GetIdx(1))) {
            const IdxT dof= it->Unknowns(st_idx_.GetIdx(1));
            if (dof >= st_idx_.NumIniUnknowns() && dof < st_idx_.NumIniUnknowns() + st_idx_.NumFiniUnknowns())
                ic.Data[dof - st_idx_.NumIniUnknowns()]= st_ic_[dof - (cG_in_t_ ? st_idx_.NumIniUnknowns() : 0)];
        }
    }

    st_ic_.resize(0);
    st_oldic_.resize(0);
    st_idx_.DeleteNumbering(MG_);
}

void SurfactantSTP1CL::DoStep (double new_t)
{
    ScopeTimerCL timer("SurfactantSTP1CL::DoStep");

    InitStep(new_t);
    DoStep();
    CommitStep();
}


/// =Methods with transport on the domain=

class TransportP1FunctionCL
/// Solve D/Dt u = 0, u(t^0) given, with SDFEM
{
  public:
    typedef BndDataCL<>    BndDataT;
    typedef P1EvalCL<double, const BndDataT, VecDescCL>       DiscSolCL;
    typedef P1EvalCL<double, const BndDataT, const VecDescCL> const_DiscSolCL;

    const IdxDescCL&    p1idx;

  private:
    MultiGridCL&        MG_;
    double              SD_,    ///< streamline diffusion
                        theta_,
                        dt_;
    MatrixCL            L_;
    BndDataT            Bnd_;
    GSPcCL              pc_;
    GMResSolverCL<GSPcCL>  gm_;

  public:
    MatrixCL E_old,
             E_,
             H_old,
             H_;

    template<class DiscVelSolT>
    TransportP1FunctionCL (MultiGridCL& mg, const DiscVelSolT& v_old, const IdxDescCL& thep1idx, double dt, double theta= 0.5, double SD= 0., int iter= 1000, double tol= 1e-7)
    : p1idx(thep1idx), MG_(mg), SD_(SD),
        theta_(theta), dt_(dt), Bnd_(BndDataT(mg.GetBnd().GetNumBndSeg())),
        gm_(pc_, 500, iter, tol, true) {
        if (theta_ != 1.)
            SetupSystem(v_old, E_old, H_old);
    }

    /// \remarks call SetupSystem \em before calling SetTimeStep!
    template<class DiscVelSolT>
    void SetupSystem (const DiscVelSolT&, MatrixCL& E, MatrixCL& H);
    /// perform one time step
    template <class DiscVelSolT>
    void DoStep (VectorCL& u, const DiscVelSolT& /* new velocity*/);

};

template<class DiscVelSolT>
void TransportP1FunctionCL::SetupSystem (const DiscVelSolT& vel, MatrixCL& E, MatrixCL& H)
// Sets up the stiffness matrices:
// E is of mass matrix type:    E_ij = (v_j       , v_i + SD * u grad v_i)
// H describes the convection:  H_ij = (u grad v_j, v_i + SD * u grad v_i)
// where v_i, v_j denote the ansatz functions.
{
    const IdxT num_unks= p1idx.NumUnknowns();
    const Uint lvl= p1idx.TriangLevel();

    SparseMatBuilderCL<double> bE(&E, num_unks, num_unks),
                               bH(&H, num_unks, num_unks);
    IdxT Numb[4];

    std::cout << "entering TransportP1Function::SetupSystem: " << num_unks << "  unknowns. ";
    // std::cout << "SD_: " << SD_ << " dt_: " << dt_ << " theta_ : " << theta_ << "\n";

    // fill value part of matrices
    Quad5CL<Point3DCL>  u_loc;
    Point3DCL Grad[4];
    Quad5CL<> u_Grad[4], // fuer u grad v_i
              p1[4];
    double det, absdet, h_T;

    LocalP1CL<> p1dummy;
    for (int i= 0; i < 4; ++i) {
        p1dummy[i]= 1.0;
        p1[i].assign(p1dummy);
        p1dummy[i]= 0.0;
    }

    DROPS_FOR_TRIANG_CONST_TETRA(const_cast<const MultiGridCL&>(MG_), lvl, sit) {
        P1DiscCL::GetGradients(Grad, det, *sit);
        absdet= std::fabs(det);
        h_T= std::pow(absdet, 1./3.);

        GetLocalNumbP1NoBnd(Numb, *sit, p1idx);

        u_loc.assign(*sit, vel);
        for(int i=0; i<4; ++i)
            u_Grad[i]= dot(Grad[i], u_loc);

        /// \todo fixed limit for maxV (maxV_limit), any better idea?
        double maxV = 1.; // scaling of SD parameter
        for(int i= 0; i < 4; ++i)    // assemble row Numb[i]
            for(int j= 0; j < 4; ++j) {
                // E is of mass matrix type:    E_ij = (v_j       , v_i + SD * u grad v_i)
               bE(Numb[i], Numb[j])+= P1DiscCL::GetMass(i,j) * absdet
                                       + Quad5CL<>(u_Grad[i]*p1[j]).quad(absdet)*SD_/maxV*h_T;

               // H describes the convection:  H_ij = (u grad v_j, v_i + SD * u grad v_i)
               bH(Numb[i], Numb[j])+= Quad5CL<>(u_Grad[j]*p1[i]).quad(absdet)
                                       + Quad5CL<>(u_Grad[i]*u_Grad[j]).quad(absdet) * SD_/maxV*h_T;
            }
    }
    bE.Build();
    bH.Build();
    std::cout << E.num_nonzeros() << " nonzeros in E, "
              << H.num_nonzeros() << " nonzeros in H! " << std::endl;
}

template <class DiscVelSolT>
void TransportP1FunctionCL::DoStep (VectorCL& u, const DiscVelSolT& vel)
{
    SetupSystem(vel, E_, H_);
    L_.clear();
    L_.LinComb(1./dt_, E_, theta_, H_);

    VectorCL rhs((1./dt_)*u);
    if (theta_ != 1.) {
        GMResSolverCL<GSPcCL> gm(gm_);
        VectorCL tmp(rhs.size());
        gm.Solve(E_old, tmp, VectorCL(H_old*u), p1idx.GetEx());
        std::cout << "TransportP1FunctionCL::DoStep rhs: res = " << gm.GetResid() << ", iter = " << gm.GetIter() << std::endl;
        rhs-= (1. - theta_)*tmp;
    }
    gm_.Solve(L_, u, VectorCL(E_*rhs), p1idx.GetEx());

    std::cout << "TransportP1FunctionCL::DoStep: res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
}


void SurfactantCharTransportP1CL::Update()
{
    // ScopeTimerCL timer("SurfactantCharTransportP1CL::Update");
    // std::cout << "SurfactantCharTransportP1CL::Update:\n";

    IdxDescCL* cidx= ic.RowIdx;
    M.Data.clear();
    M.SetIdx(cidx, cidx);
    A.Data.clear();
    A.SetIdx(cidx, cidx);
    Md.Data.clear();
    Md.SetIdx(cidx, cidx);
    VecDescCL vd_load(&idx);

    TetraAccumulatorTupleCL accus;
    InterfaceCommonDataP1CL cdata(lset_vd_, lsetbnd_);
    accus.push_back(&cdata);
    InterfaceMatrixAccuCL<LocalInterfaceMassP1CL, InterfaceCommonDataP1CL> mass_accu(&M, LocalInterfaceMassP1CL(), cdata, "mass");
    accus.push_back(&mass_accu);
    InterfaceMatrixAccuCL<LocalLaplaceBeltramiP1CL, InterfaceCommonDataP1CL> lb_accu(&A, LocalLaplaceBeltramiP1CL(D_), cdata, "Laplace-Beltrami");
    accus.push_back(&lb_accu);
    accus.push_back_acquire(make_wind_dependent_matrixP1_accu<LocalInterfaceMassDivP1CL>(&Md, cdata,  make_P2Eval(MG_, Bnd_v_, *v_), "massdiv"));
    if (rhs_fun_)
        accus.push_back_acquire(new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>(&vd_load, LocalVectorP1CL(rhs_fun_, ic.t), cdata, "load"));

    accumulate(accus, MG_, cidx->TriangLevel(), cidx->GetBndInfo());

    load.resize(idx.NumUnknowns());
    load= vd_load.Data;

//     WriteToFile(M.Data, "chartranspM.txt", "mass");
//     WriteToFile(A.Data, "chartranspA.txt", "Laplace-Beltrami");
//     WriteToFile(Md.Data,"chartranspMd.txt","mass-div");
//     WriteToFile(vd_load.Data,"chartranspload.txt","load");

    // std::cout << "SurfactantCharTransportP1CL::Update: Finished\n";
}

void SurfactantCharTransportP1CL::InitStep (double new_t)
{
    // ScopeTimerCL timer("SurfactantcGP1CL::InitStep");
    std::cout << "SurfactantCharTransportP1CL::InitStep:\n";

    ic.t= new_t;
    dt_= ic.t - oldt_;
    idx.CreateNumbering(oldidx_.TriangLevel(), MG_, &lset_vd_, &lsetbnd_); // InitTimeStep deletes oldidx_ and swaps idx and oldidx_.
    std::cout << "new NumUnknowns: " << idx.NumUnknowns();
    full_idx.CreateNumbering(idx.TriangLevel(), MG_);
    std::cout << " full NumUnknowns: " << full_idx.NumUnknowns() << std::endl;

    fulltransport_= new TransportP1FunctionCL(MG_, make_P2Eval(MG_, Bnd_v_, oldv_), full_idx, dt_,  /*theta=*/theta_, /*SD=*/ 0.1, /*iter=*/ 2000, /*tol=*/ 0.1*gm_.GetTol());

    VecDescCL rhs(&idx);
    rhs.Data= (1./dt_)*ic.Data;
//     if (theta_ != 1.) {
//         GMResSolverCL<GSPcCL> gm(gm_);
//         VectorCL tmp(rhs.Data.size());
//         gm.Solve(M.Data, tmp, VectorCL(A.Data*ic.Data + Md.Data*ic.Data));
//         std::cout << "SurfactantP1CL::InitStep: rhs: res = " << gm.GetResid() << ", iter = " << gm.GetIter() << std::endl;
//         rhs.Data-= (1. - theta_)*tmp;
//     }
    DROPS::VecDescCL rhsext(&full_idx);
    DROPS::Extend(MG_, rhs, rhsext);
    rhs_.resize(rhsext.Data.size());
    rhs_= rhsext.Data;
}

void SurfactantCharTransportP1CL::DoStep ()
{
    VecDescCL transp_rhs(&idx),
              transp_rhsext(&full_idx);
    transp_rhsext.Data= rhs_;
    fulltransport_->DoStep(transp_rhsext.Data, make_P2Eval(MG_, Bnd_v_, *v_));
    Restrict(MG_, transp_rhsext, transp_rhs);

    ic.SetIdx(&idx);
    Update();
    // L_.LinComb(1./dt_, M.Data, theta_, A.Data, theta_, Md.Data);
    L_.LinComb(1./dt_, M.Data, 1., A.Data, 1., Md.Data);
    const VectorCL therhs(M.Data*transp_rhs.Data + load);
    std::cout << "Before solve: res = " << norm(L_*ic.Data - therhs) << std::endl;
    {
        ScopeTimerCL timer("SurfactantCharTransportP1CL::DoStep: Solve");
        gm_.Solve(L_, ic.Data, therhs, idx.GetEx());
    }
    std::cout << "SurfactantCharTransportP1CL::DoStep: res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
}

void SurfactantCharTransportP1CL::CommitStep ()
{
    full_idx.DeleteNumbering(MG_);
    delete fulltransport_;
}

void SurfactantCharTransportP1CL::DoStep (double new_t)
{
    ScopeTimerCL timer("SurfactantCharTransportP1CL::::DoStep");

    InitStep(new_t);
    DoStep();
    CommitStep();
}

class LocalCahnHilliardCL
{
  private:
    const PrincipalLatticeCL& lat;
    LocalP1CL<> P1Hat[4];
    LocalP2CL<> P2Hat[10];
    Point3DCL P1Grad[4];
    LocalP1CL<Point3DCL>  P2GradRef[10], P2Grad[10];
    GridFunctionCL<Point3DCL> qnormal;
    GridFunctionCL<Point3DCL> qSurfP1Grad[4];
    GridFunctionCL<Point3DCL> q3DsurfP1grad[4];

    GridFunctionCL<> mobility2D,well_potential_second_derivative[4], qP1Hat[4], qnormal_comp[3];


    QuadDomain2DCL  q2Ddomain;
    std::valarray<double> ls_loc;
    SurfacePatchCL spatch;

    QuadDomainCL q3Ddomain;
    GridFunctionCL<Point3DCL> q3Dnormal;
    GridFunctionCL<Point3DCL> q3DP1Grad[4],q2DP1Grad[4];

    void qP2ls(const LocalP2CL<>& ls, LocalP1CL<Point3DCL>&);

  public:
    LocalCahnHilliardCL ()
        : lat(PrincipalLatticeCL::instance(2)), ls_loc(lat.vertex_size())
    {
        P1DiscCL::GetP1Basis(P1Hat);
        P2DiscCL::GetP2Basis(P2Hat);
        P2DiscCL::GetGradientsOnRef(P2GradRef);
     }

    void calcIntegrands(const SMatrixCL<3,3>& T, const LocalP2CL<>& ls, const LocalP1CL<Point3DCL>& v,const LocalP1CL<>& chi , const TetraCL& tet); ///< has to be called before any setup method!
    void calc3DIntegrands(const SMatrixCL<3,3>& T, const LocalP2CL<>& ls, const TetraCL& tet); ///< has to be called after calcIntegrands!

    void setupM_P1 (double M_P1[4][4]);
    void setupNormalStab_P1 (double NormalStab_P1[4][4],double absdet);
    void setupTangentStab_P1 (double TangentStab_P1[4][4],double absdet);
    void setupVolumeStab_P1 (double VolumeStab_P1[4][4],double absdet);


    void setupL_P1 (double L_P1[4][4]);
    void setupLM_P1 (double LM_P1[4][4]);
    void setupGprimeprime_P1 (double Gprimeprime_P1[4][4]);


};

void LocalCahnHilliardCL::calc3DIntegrands(const SMatrixCL<3,3>& T, const LocalP2CL<>& ls, const TetraCL& tet)
{
    make_SimpleQuadDomain<Quad5DataCL> (q3Ddomain, AllTetraC);
    LocalP1CL<Point3DCL> Normals;
    qP2ls(ls, Normals);

    resize_and_evaluate_on_vertexes (Normals, q3Ddomain, q3Dnormal);

    // Scale Normals accordingly to the Euclidean Norm (only consider the ones which make a contribution in the sense of them being big enough... otherwise one has to expect problems with division through small numbers)
    for(Uint i=0; i<q3Dnormal.size(); ++i) {
         //if(q3Dnormal[i].norm()> 1e-8)
         q3Dnormal[i]= q3Dnormal[i]/q3Dnormal[i].norm();
    }

    for(int j=0; j<4; ++j) {
        q3DP1Grad[j].resize(q3Ddomain.vertex_size());
        q3DP1Grad[j]= P1Grad[j];
    }
}

double inverse_square_root(double x) {return(1./std::sqrt(x));}

double Mobility_function(double x, double t)
{
    double scaling=1.0; //1./Density_function(x, t); //std::exp(-1000*t);
	return( scaling*std::sqrt((1.-x)*(x)*(1.-x)*(x)));
    //return(1);
}

double Diffusion_function(double x, double t)
    {
        double scaling=1.; //std::exp(-1000*t);
        return( x);
        //return(1);
    }

double Density_function(double x, double t) {
    double density_ratio= 1.; //std::exp(-1000*t);
    return ( (1.-x) +  density_ratio*(x));
    //return(1);
}


double Potential_function(const double x)
    {
        return((1.-x)*(x)*(1.-x)*(x)*0.25);
    }

double Potential_prime_function(const double x)
{
    //return(0);
    return((x-1.)*(x)*(x-0.5));
}

    double Potential_prime_convex_function(const double x)
    {
      return(x/2. + x*x*x);
    }

    double Potential_prime_concave_function(const double x)
    {
        return(-(3./2.)*x*x);
    }

    double Potential_prime_prime_function(const double x)
    {
        return((x-1.)*(x)+2.*(x-0.5)*(x-0.5));
        //return (-(0.)*3.*x+(1.)*1./2. + (0.)*3.*x*x);
    }

    double Potential_prime_prime_convex_function(const double x)
    {
        return(1./2.+3.*x*x);
        //return (-(0.)*3.*x+(1.)*1./2. + (0.)*3.*x*x);
    }

void LocalCahnHilliardCL::calcIntegrands(const SMatrixCL<3,3>& T, const LocalP2CL<>& ls, const LocalP1CL<Point3DCL>& v ,const LocalP1CL<>& chi , const TetraCL& tet)
{

    P2DiscCL::GetGradients(P2Grad, P2GradRef, T);
    P1DiscCL::GetGradients(P1Grad, T);
    evaluate_on_vertexes(ls, lat, Addr(ls_loc));
    spatch.make_patch<MergeCutPolicyCL>(lat, ls_loc);

    // The routine takes the information about the tetrahedra and the cutting surface and generates a two-dimensional triangulation of the cut, including the necessary point-positions and weights for the quadrature
    make_CompositeQuad5Domain2D (q2Ddomain, spatch, tet);
    LocalP1CL<Point3DCL> Normals;
    qP2ls(ls, Normals);

    // Resize and evaluate Normals at all points which are needed for the two-dimensional quadrature-rule
    resize_and_evaluate_on_vertexes (Normals, q2Ddomain, qnormal);

    for(int j=0; j<4; ++j) {
    	q2DP1Grad[j].resize(q2Ddomain.vertex_size());
        q2DP1Grad[j]= P1Grad[j];
    }

    // Scale Normals accordingly to the Euclidean Norm (only consider the ones which make a contribution in the sense of them being big enough... otherwise one has to expect problems with division through small numbers)
    for(Uint i=0; i<qnormal.size(); ++i) {
         //if(qnormal[i].norm()> 1e-8)
         qnormal[i]= qnormal[i]/qnormal[i].norm();
    }

    // Provide all components of the normals
    for (int k=0; k<3; ++k) {
        qnormal_comp[k].resize(q2Ddomain.vertex_size());
        ExtractComponent(qnormal, qnormal_comp[k], k);
    }

    for(int j=0; j<4 ;++j) {
        resize_and_evaluate_on_vertexes(P1Hat[j], q2Ddomain, qP1Hat[j]);
        qSurfP1Grad[j].resize(q2Ddomain.vertex_size());
        qSurfP1Grad[j]= P1Grad[j];
        qSurfP1Grad[j]-= dot(qSurfP1Grad[j], qnormal)*qnormal;
    }

    for(int j=0; j<4 ;++j) {
        q3DsurfP1grad[j].resize(q3Ddomain.vertex_size());
        q3DsurfP1grad[j]= P1Grad[j];
        q3DsurfP1grad[j]-= dot(q3DsurfP1grad[j], q3Dnormal)*q3Dnormal;
    }

    //physical velocity, tetra P1 function
    LocalP1CL<Point3DCL> velocity;
    for(int i=0; i<4 ; ++i)
    {
    	velocity += v[i]*P1Hat[i];
    }

    //volume fraction , tetra P1 function
      LocalP1CL<> mobility;
      for(int i=0; i<4 ; ++i)
      {
    	mobility += Mobility_function(chi[i])*P1Hat[i];
      }
      resize_and_evaluate_on_vertexes (mobility, q2Ddomain, mobility2D);

    //LocalP1CL<> well_potential_second_derivative[4];
    for(int i=0; i<4 ; ++i)
    {
        well_potential_second_derivative[i] = Potential_prime_prime_function(chi[i])*qP1Hat[i];
    }
    //resize_and_evaluate_on_vertexes (well_potential_second_derivative, q2Ddomain, well_potential_second_derivative_2D);


}

// The P2 levelset-function is used to compute the normals which are needed for the (improved) projection onto the interface, GradLP1 has to be set before
void LocalCahnHilliardCL::qP2ls(const LocalP2CL<>& ls, LocalP1CL<Point3DCL>& Normals)
{
    for(int i=0; i<10 ; ++i)
    {
        Normals+=ls[i]*P2Grad[i];
    }

}

void LocalCahnHilliardCL::setupM_P1 (double M_P1[4][4])
{
    // Do all combinations for (i,j) i,j=4 x 4 and corresponding quadrature
    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
             M_P1[i][j]= quad_2D(qP1Hat[i]*qP1Hat[j], q2Ddomain);
        }
    }
}

    void LocalCahnHilliardCL::setupNormalStab_P1 (double NormalStab_P1[4][4], double absdet)
    {
        // Do all combinations for (i,j) i,j=4 x 4 and corresponding quadrature
        for (int i=0; i<4; ++i) {
            for (int j=0; j<4; ++j) {
                NormalStab_P1[i][j]= quad(dot(q3Dnormal, q3DP1Grad[i])*dot(q3Dnormal, q3DP1Grad[j]), absdet, q3Ddomain, AllTetraC);

            }
        }
    }


    void LocalCahnHilliardCL::setupTangentStab_P1 (double TangentStab_P1[4][4], double absdet)
    {
        // Do all combinations for (i,j) i,j=4 x 4 and corresponding quadrature
        for (int i=0; i<4; ++i) {
            for (int j=0; j<4; ++j) {
                TangentStab_P1[i][j]= quad(dot(q3DsurfP1grad[i], q3DsurfP1grad[j]), absdet, q3Ddomain, AllTetraC);

            }
        }
    }

    void LocalCahnHilliardCL::setupVolumeStab_P1 (double VolumeStab_P1[4][4], double absdet)
    {
        // Do all combinations for (i,j) i,j=4 x 4 and corresponding quadrature
        for (int i=0; i<4; ++i) {
            for (int j=0; j<4; ++j) {
                VolumeStab_P1[i][j]= quad(dot(q3DP1Grad[i], q3DP1Grad[j]), absdet, q3Ddomain, AllTetraC);

            }
        }
    }

void LocalCahnHilliardCL::setupL_P1 (double L_P1[4][4])
{
    // Do all combinations for (i,j) i,j=4 x 4 and corresponding quadrature
    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
             L_P1[i][j]= quad_2D(dot(qSurfP1Grad[i], qSurfP1Grad[j]), q2Ddomain);
        }
    }
}

void LocalCahnHilliardCL::setupLM_P1 (double LM_P1[4][4])
{
    // Do all combinations for (i,j) i,j=4 x 4 and corresponding quadrature
    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
             LM_P1[i][j]= quad_2D(mobility2D*dot(qSurfP1Grad[i], qSurfP1Grad[j]), q2Ddomain);
        }
    }
}

void LocalCahnHilliardCL::setupGprimeprime_P1 (double Gprimeprime_P1[4][4])
    {
        // Do all combinations for (i,j) i,j=4 x 4 and corresponding quadrature
        for (int i=0; i<4; ++i) {
            for (int j=0; j<4; ++j) {
                Gprimeprime_P1[i][j]= quad_2D(well_potential_second_derivative[i]*qP1Hat[j], q2Ddomain);
            }
        }
    }

/// \brief Basis Class for Accumulator to set up the matrices for interface Stokes.
class CahnHilliardIFAccumulator_P1P1CL : public TetraAccumulatorCL
{
  protected:
    const VecDescCL& lset;
    const VecDescCL& velocity;
    const VecDescCL& volume_fraction;
    const LsetBndDataCL& lset_bnd;
    const BndDataCL<Point3DCL>& velocity_bnd;
    const BndDataCL<>& volume_fraction_bnd;
    IdxDescCL &P1Idx_, &ScalarP1Idx_;


    MatrixCL &M_P1_, &NormalStab_P1_, &TangentStab_P1_,&VolumeStab_P1_, &L_P1P1_, &LM_P1P1_, &Gprimeprime_P1P1_;
    MatrixBuilderCL *mM_P1_,*mNormalStab_P1_, *mTangentStab_P1_, *mVolumeStab_P1_, *mL_P1P1_, *mLM_P1P1_, *mGprimeprime_P1P1_;
    double locM_P1[4][4], locNormalStab_P1[4][4], locTangentStab_P1[4][4],  locVolumeStab_P1[4][4], locL_P1P1[4][4],
            locLM_P1P1[4][4], locGprimeprime_P1P1[4][4];

    LocalCahnHilliardCL localCahnHilliard_;

    SMatrixCL<3,3> T;
    double det, absdet;
    LocalP2CL<> ls_loc;
    LocalP1CL<Point3DCL> v_loc;
    LocalP1CL<> chi_loc;
    LocalNumbP1CL n, nScalar;

    ///\brief Computes the mapping from local to global data "n", the local matrices in loc and.
    void local_setup (const TetraCL& tet);
    ///\brief Update the global system->
    void update_global_system ();

  public:
    CahnHilliardIFAccumulator_P1P1CL (const VecDescCL& ls, const VecDescCL& v,const VecDescCL& chi, const LsetBndDataCL& ls_bnd, const BndDataCL<Point3DCL>& v_bnd,const BndDataCL<>& chi_bnd, IdxDescCL& P1FE, IdxDescCL& ScalarP1FE,  MatrixCL& M_P1, MatrixCL& NormalStab_P1,MatrixCL& TangentStab_P1, MatrixCL& VolumeStab_P1,  MatrixCL& L_P1P1,  MatrixCL& LM_P1P1 ,  MatrixCL& Gprimeprime_P1P1);

    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    void visit (const TetraCL& sit);

    TetraAccumulatorCL* clone (int /*tid*/) { return new CahnHilliardIFAccumulator_P1P1CL (*this); }

};

CahnHilliardIFAccumulator_P1P1CL::CahnHilliardIFAccumulator_P1P1CL(const VecDescCL& ls, const VecDescCL& v, const VecDescCL& chi, const LsetBndDataCL& ls_bnd, const BndDataCL<Point3DCL>& v_bnd,const BndDataCL<>& chi_bnd, IdxDescCL& P1FE, IdxDescCL& ScalarP1FE, MatrixCL& M_P1, MatrixCL& NormalStab_P1, MatrixCL& TangentStab_P1,  MatrixCL& VolumeStab_P1,  MatrixCL& L_P1P1, MatrixCL& LM_P1P1, MatrixCL& Gprimeprime_P1P1)
 : lset(ls), velocity(v), volume_fraction(chi), lset_bnd(ls_bnd), velocity_bnd(v_bnd),volume_fraction_bnd(chi_bnd), P1Idx_(P1FE), ScalarP1Idx_(ScalarP1FE), M_P1_(M_P1), NormalStab_P1_(NormalStab_P1), TangentStab_P1_(TangentStab_P1), VolumeStab_P1_(VolumeStab_P1), L_P1P1_(L_P1P1), LM_P1P1_(LM_P1P1), Gprimeprime_P1P1_(Gprimeprime_P1P1), localCahnHilliard_()
{}

void CahnHilliardIFAccumulator_P1P1CL::begin_accumulation ()
{
    std::cout << "entering CahnHilliardIF: \n";
    const size_t num_unks_p1 = P1Idx_.NumUnknowns(), num_unks_scalarp1 = ScalarP1Idx_.NumUnknowns();

    mM_P1_= new MatrixBuilderCL(&M_P1_, num_unks_scalarp1, num_unks_scalarp1);
    mNormalStab_P1_= new MatrixBuilderCL(&NormalStab_P1_, num_unks_scalarp1, num_unks_scalarp1);
    mTangentStab_P1_= new MatrixBuilderCL(&TangentStab_P1_, num_unks_scalarp1, num_unks_scalarp1);
    mVolumeStab_P1_= new MatrixBuilderCL(&VolumeStab_P1_, num_unks_scalarp1, num_unks_scalarp1);
    mL_P1P1_= new MatrixBuilderCL(&L_P1P1_, num_unks_scalarp1, num_unks_scalarp1);
    mLM_P1P1_= new MatrixBuilderCL(&LM_P1P1_, num_unks_scalarp1, num_unks_scalarp1);
    mGprimeprime_P1P1_= new MatrixBuilderCL(&Gprimeprime_P1P1_, num_unks_scalarp1, num_unks_scalarp1);



}

void CahnHilliardIFAccumulator_P1P1CL::finalize_accumulation ()
{


    mL_P1P1_->Build();
    delete mL_P1P1_;

    mLM_P1P1_->Build();
    delete mLM_P1P1_;

    mGprimeprime_P1P1_->Build();
    delete mGprimeprime_P1P1_;

    mM_P1_->Build();
    delete mM_P1_;

    mNormalStab_P1_->Build();
    delete mNormalStab_P1_;

    mTangentStab_P1_->Build();
    delete mTangentStab_P1_;

    mVolumeStab_P1_->Build();
    delete mVolumeStab_P1_;

#ifndef _PAR
    std::cout << "CahnHilliardIF_P1P1:\t"
    		  << M_P1_.num_nonzeros() 	<< " nonzeros in Mass, "
             << NormalStab_P1_.num_nonzeros() 	<< " nonzeros in Normal, "
            << TangentStab_P1_.num_nonzeros() 	<< " nonzeros in Tang, "
            << VolumeStab_P1_.num_nonzeros() 	<< " nonzeros in Vol, "
            << L_P1P1_.num_nonzeros() << " nonzeros in Laplacian, "
			  << LM_P1P1_.num_nonzeros() << " nonzeros in LaplacianM, "
            << Gprimeprime_P1P1_.num_nonzeros() << " nonzeros in LaplacianM, "  << std::endl;

#endif
    // std::cout << '\n';
}

void CahnHilliardIFAccumulator_P1P1CL::visit (const TetraCL& tet)
{
    ls_loc.assign(tet, lset, lset_bnd);
    v_loc.assign(tet, velocity, velocity_bnd);
    chi_loc.assign(tet, volume_fraction, volume_fraction_bnd);

    if (isInCutMesh(ls_loc))
    {
        local_setup(tet);
        update_global_system();
    }

}

void CahnHilliardIFAccumulator_P1P1CL::local_setup (const TetraCL& tet)
{
    GetTrafoTr(T, det, tet);
    absdet= std::fabs(det);

    n.assign(tet, P1Idx_, P1Idx_.GetBndInfo());
    nScalar.assign(tet, ScalarP1Idx_, ScalarP1Idx_.GetBndInfo());


    localCahnHilliard_.calcIntegrands(T, ls_loc, v_loc, chi_loc, tet);
    localCahnHilliard_.calc3DIntegrands(T, ls_loc, tet);

    localCahnHilliard_.setupM_P1(locM_P1);
    localCahnHilliard_.setupNormalStab_P1(locNormalStab_P1, absdet);
    localCahnHilliard_.setupTangentStab_P1(locTangentStab_P1, absdet);
    localCahnHilliard_.setupVolumeStab_P1(locVolumeStab_P1, absdet);


    localCahnHilliard_.setupL_P1(locL_P1P1);
    localCahnHilliard_.setupLM_P1(locLM_P1P1);

    localCahnHilliard_.setupGprimeprime_P1(locGprimeprime_P1P1);


}

void CahnHilliardIFAccumulator_P1P1CL::update_global_system ()
{
    MatrixBuilderCL& mM_P1 = *mM_P1_;
    MatrixBuilderCL& mNormalStab_P1 = *mNormalStab_P1_;
    MatrixBuilderCL& mTangentStab_P1 = *mTangentStab_P1_;
    MatrixBuilderCL& mVolumeStab_P1 = *mVolumeStab_P1_;

    MatrixBuilderCL& mL_P1P1 = *mL_P1P1_;
    MatrixBuilderCL& mLM_P1P1 = *mLM_P1P1_;
    MatrixBuilderCL& mGprimeprime_P1P1 = *mGprimeprime_P1P1_;


    int count =0;

    for(int i= 0; i < 4; ++i) {
        const IdxT ii= nScalar.num[i];
        if (ii==NoIdx || !(nScalar.WithUnknowns(i))) continue;
        for(int j=0; j <4; ++j) {
            const IdxT jj= nScalar.num[j];
            if (jj==NoIdx || !(nScalar.WithUnknowns(j))) continue;
            mM_P1(ii, jj) += locM_P1[i][j];
            mNormalStab_P1(ii, jj) += locNormalStab_P1[i][j];
            mTangentStab_P1(ii, jj) += locTangentStab_P1[i][j];
            mVolumeStab_P1(ii, jj) += locVolumeStab_P1[i][j];
            mL_P1P1(ii, jj) += locL_P1P1[i][j];
            mLM_P1P1(ii, jj) += locLM_P1P1[i][j];
            mGprimeprime_P1P1(ii, jj) += locGprimeprime_P1P1[j][i];

        }
    }
}

void SetupCahnHilliardIF_P1P1(const MultiGridCL& MG_,  MatDescCL* M_P1, MatDescCL* NormalStab_P1, MatDescCL* TangentStab_P1,MatDescCL* VolumeStab_P1, MatDescCL* L_P1P1 , MatDescCL* LM_P1P1,  MatDescCL* Gprimeprime_P1P1, const VecDescCL& lset, const LsetBndDataCL& lset_bnd, const VecDescCL& velocity,  const BndDataCL<Point3DCL>& velocity_bnd, const VecDescCL& chi,const BndDataCL<>& chi_bnd)
{
  ScopeTimerCL scope("SetupCahnHilliardIF_P1P1");
  CahnHilliardIFAccumulator_P1P1CL accu(lset, velocity, chi, lset_bnd, velocity_bnd, chi_bnd, *(M_P1->RowIdx), *(L_P1P1->RowIdx),  M_P1->Data, NormalStab_P1->Data, TangentStab_P1->Data, VolumeStab_P1->Data,  L_P1P1->Data,  LM_P1P1->Data,  Gprimeprime_P1P1->Data);
  TetraAccumulatorTupleCL accus;
  //    MaybeAddProgressBar(MG_, "LapBeltr(P2) Setup", accus, RowIdx.TriangLevel());
  accus.push_back(&accu);
  accumulate(accus, MG_, M_P1->GetRowLevel(), /*M_P1->RowIdx->GetMatchingFunction(),*/ M_P1->RowIdx->GetBndInfo());

}

} // end of namespace DROPS