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
#include <cstring>
#include <cmath>


namespace DROPS {

void Extend (const MultiGridCL& mg, const VecDescCL& x, VecDescCL& xext)
{
    const Uint xidx( x.RowIdx->GetIdx()),
        xextidx( xext.RowIdx->GetIdx()),
        lvl( x.RowIdx->TriangLevel());
    xext.Data= 0.;
    const int stride= x.RowIdx->NumUnknownsVertex();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
        if (it->Unknowns.Exist( xidx) && it->Unknowns.Exist( xextidx))
            for (int k=0; k<stride; ++k)
                xext.Data[it->Unknowns( xextidx)+k]= x.Data[it->Unknowns( xidx)+k];
    }
    if (x.RowIdx->GetFE() == P1IF_FE)
        return;

    // For P2IF_FE, also fixup the edge-dofs.
    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it) {
        if (it->Unknowns.Exist( xidx) && it->Unknowns.Exist( xextidx))
            for (int k=0; k<stride; ++k)
                xext.Data[it->Unknowns( xextidx)+k]= x.Data[it->Unknowns( xidx)+k];
    }
}

///\todo extend to vector-valued FE
void Restrict (const MultiGridCL& mg, const VecDescCL& xext, VecDescCL& x)
{
    const Uint xidx( x.RowIdx->GetIdx()),
        xextidx( xext.RowIdx->GetIdx()),
        lvl( x.RowIdx->TriangLevel());

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
        if (it->Unknowns.Exist( xidx) && it->Unknowns.Exist( xextidx))
            x.Data[it->Unknowns( xidx)]= xext.Data[it->Unknowns( xextidx)];
    }
    if (x.RowIdx->GetFE() == P1IF_FE)
        return;

    // For P2IF_FE, also fixup the edge-dofs.
    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it) {
        if (it->Unknowns.Exist( xidx) && it->Unknowns.Exist( xextidx))
            x.Data[it->Unknowns( xidx)]= xext.Data[it->Unknowns( xextidx)];
    }
}

void GetLocalNumbInterface(IdxT* Numb, const TetraCL& s, const IdxDescCL& idx)
{
    const Uint sys= idx.GetIdx();
    for (Uint i= 0; i < 4; ++i)
        Numb[i]= s.GetVertex( i)->Unknowns.Exist( sys) ? s.GetVertex( i)->Unknowns( sys) : NoIdx;
    if (idx.GetFE() == P1IF_FE)
        return;

    for(Uint i= 0; i < 6; ++i)
        Numb[i+4]= s.GetEdge( i)->Unknowns.Exist( sys) ? s.GetEdge( i)->Unknowns( sys) : NoIdx;
}


InterfaceCommonDataP2CL::InterfaceCommonDataP2CL (const VecDescCL& ls_arg, const BndDataCL<>& lsetbnd_arg,
    const QuaQuaMapperCL& quaquaarg, const PrincipalLatticeCL& lat_arg)
    : ls( &ls_arg), lsetbnd( &lsetbnd_arg), lat( &lat_arg),  compute_quaddomains_( true), ls_loc( lat->vertex_size()),
      quaqua( quaquaarg)
{
    qdom_projected.compute_absdets( true);
    P2DiscCL::GetGradientsOnRef( gradrefp2);
    for (Uint i= 0; i < 10 ; ++i)
        p2[i][i]= 1.; // P2-Basis-Functions
}

InterfaceCommonDataDeformP2CL::InterfaceCommonDataDeformP2CL (const VecDescCL& ls_arg, const BndDataCL<>& lsetbnd_arg,
        VecDescCL& Psi_vdarg, const PrincipalLatticeCL& lat_arg)
    : ls( &ls_arg), lsetbnd( &lsetbnd_arg), lat( &lat_arg), Psi_vd (&Psi_vdarg), ls_loc( lat->vertex_size()),
      Phi (this)
{
    P2DiscCL::GetGradientsOnRef( gradrefp2);
    for (Uint i= 0; i < 10 ; ++i)
        p2[i][i]= 1.; // P2-Basis-Functions
}

void SetupInterfaceMassP1 (const MultiGridCL& mg, MatDescCL* matM, const VecDescCL& ls, const BndDataCL<>& lsetbnd, double alpha)
{
    //ScopeTimerCL timer( "SetupInterfaceMassP1");

    TetraAccumulatorTupleCL accus;
    InterfaceCommonDataP1CL cdata( ls, lsetbnd);
    accus.push_back( &cdata);
    InterfaceMatrixAccuCL<LocalInterfaceMassP1CL, InterfaceCommonDataP1CL> accu( matM, LocalInterfaceMassP1CL( alpha), cdata);
    accus.push_back( &accu);
    const IdxDescCL* RowIdx= matM->RowIdx;
    accumulate( accus, mg, RowIdx->TriangLevel(), RowIdx->GetBndInfo());

    // WriteToFile( matM->Data, "mass.txt", "mass");
}

void SetupInterfaceMassP1OnTriangle (const LocalP1CL<> p1[4], Quad5_2DCL<> q[4],
        const BaryCoordCL triangle[3], double det, double coup[4][4])
{
    for (int i= 0; i < 4; ++i)
        q[i].assign( p1[i], triangle);

    Quad5_2DCL<> m;
    for (int i= 0; i < 4; ++i) {
        m= q[i]*q[i];
        coup[i][i]+= m.quad( det);
        for(int j= 0; j < i; ++j) {
            m= (q[j]*q[i]);
            coup[i][j]+= m.quad( det);
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
    MatrixBuilderCL M( &mat->Data, num_rows, num_cols);

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
    DROPS_FOR_TRIANG_CONST_TETRA( MG, lvl, it) {
        // setup local matrix
        triangle.Init( *it, ls, lsetbnd);
        if (triangle.Intersects()) { // We are at the phase boundary.
            std::memset( coup, 0, 4*4*sizeof( double));

            for (int ch= 0; ch < 8; ++ch) {
                if (!triangle.ComputeForChild( ch)) // no patch for this child
                    continue;

                for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri)
                    SetupInterfaceMassP1OnTriangle( p1, q, &triangle.GetBary( tri), triangle.GetAbsDet( tri), coup);
            }

            // update global matrix

            GetLocalNumbP1NoBnd( numr, *it, *mat->RowIdx);
            GetLocalNumbP1NoBnd( numc, *it, *mat->ColIdx);
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
                    M( numr[i],   numc[j])+= Kboth*coup[j][i];
                    if (xnumc[j] == NoIdx) continue;
                    M( numr[i],   xnumc[j])+= K[j]*coup[j][i];
                }
                if (xnumr[i] == NoIdx) continue;
                for(int j= 0; j < 4; ++j) {
                    if (numc[j] == NoIdx) continue;
                    M( xnumr[i],   numc[j])+= K[i]*coup[j][i];
                    if (xnumc[j] == NoIdx || is_pos[i]!= is_pos[j]) continue;
                    M( xnumr[i],   xnumc[j])+= std::abs(K[i])*coup[j][i];
                }
            }
        }
    }

    M.Build();
}

void SetupInterfaceSorptionP1X (const MultiGridCL& MG, const VecDescCL& ls, const BndDataCL<>& lsetbnd,
        MatDescCL* R, MatDescCL* C, MatDescCL* R_i, MatDescCL* C_i, const IdxDescCL* mass_idx, const IdxDescCL* surf_idx, const double k_a[2], const double k_d[2])
{
    R->SetIdx  ( mass_idx, mass_idx);
    C->SetIdx  ( mass_idx, surf_idx);
    R_i->SetIdx( surf_idx, surf_idx);
    C_i->SetIdx( surf_idx, mass_idx);

    SetupInterfaceMassP1( MG, R_i, ls, lsetbnd, k_d[0]+k_d[1]);
    SetupInterfaceMassP1X( MG, R, ls, lsetbnd, k_a);

    const double mk_a[2]= { -k_a[0], -k_a[1] },
                 mk_d[2]= { -k_d[0], -k_d[1] };
    SetupInterfaceMassP1X( MG, C, ls, lsetbnd, mk_d);
    SetupInterfaceMassP1X( MG, C_i, ls, lsetbnd, mk_a);
}

void SetupMixedMassP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd)
{
    SetupInterfaceMassP1( mg, mat, ls, lsetbnd);
}

void SetupLBP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd, double D)
{
    //ScopeTimerCL timer( "SetuLBP1");

    TetraAccumulatorTupleCL accus;
    InterfaceCommonDataP1CL cdata( ls, lsetbnd);
    accus.push_back( &cdata);
    InterfaceMatrixAccuCL<LocalLaplaceBeltramiP1CL, InterfaceCommonDataP1CL> accu( mat, LocalLaplaceBeltramiP1CL( D), cdata);
    accus.push_back( &accu);
    const IdxDescCL* RowIdx= mat->RowIdx;
    accumulate( accus, mg, RowIdx->TriangLevel(), RowIdx->GetBndInfo());

    // WriteToFile( mat->Data, "lb.txt", "lb");
}

 /// P1
void SetupInterfaceRhsP1OnTriangle (const LocalP1CL<> p1[4],
    Quad5_2DCL<> q[4],VectorCL& v, const IdxT Numb[4],
    const TetraCL& t, const BaryCoordCL triangle[3], double det,
    instat_scalar_fun_ptr f,double time=0)
{
    for (int i= 0; i < 4; ++i)
        q[i].assign( p1[i], triangle);
    Quad5_2DCL<> qf( t, triangle, f,time), r;

    for (int i= 0; i < 4; ++i) {
        if (Numb[i] == NoIdx) continue;
        r= qf*q[i];
        v[Numb[i]]+= r.quad( det);
    }
}


void SetupInterfaceRhsP1 (const MultiGridCL& mg, VecDescCL* v,
    const VecDescCL& ls, const BndDataCL<>& lsetbnd, instat_scalar_fun_ptr f, double t)
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

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        triangle.Init( *it, ls, lsetbnd);
        if (triangle.Intersects()) { // We are at the phase boundary.
            GetLocalNumbP1NoBnd( num, *it, *v->RowIdx);

            for (int ch= 0; ch < 8; ++ch) {
                triangle.ComputeForChild( ch);
                for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri)
                    SetupInterfaceRhsP1OnTriangle( p1, q, v->Data, num,
                        *it, &triangle.GetBary( tri), triangle.GetAbsDet( tri), f,t);
            }
        }
    }
    std::cout << " Rhs set up." << std::endl;
}

// void SetupInterfaceRhsP1 (const MultiGridCL& mg, VecDescCL* v,
//     const VecDescCL& ls, const BndDataCL<>& lsetbnd, instat_scalar_fun_ptr f, double t)
// {
//     TetraAccumulatorTupleCL accus;
//     InterfaceCommonDataP1CL cdata( ls, lsetbnd);
//     accus.push_back( &cdata);
//     InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL> loadaccu( v, LocalVectorP1CL( f, v->t), cdata);
//     accus.push_back( &loadaccu);
//     accumulate( accus, mg, v->RowIdx->TriangLevel(), v->RowIdx->GetBndInfo());

//     // WriteToFile( v->Data, "rhs.txt", "Rhs");
// }

  ////////////////////////////////////////////////////////////////////////////

  /// SetupStokes
/// \brief Setup of the local Stokes system on a tetra intersected by the dividing surface.
class LocalStokesCL
{
  private:
    const PrincipalLatticeCL& lat;
    LocalP1CL<> P1Hat[4];
    LocalP2CL<> P2Hat[10];
    LocalP1CL<Point3DCL> P2GradRef[10], P2Grad[10];
    Point3DCL P1Grad[4];
    GridFunctionCL<Point3DCL> qnormal, qrotP1[12], qvP1, qvProjP1, qProj[3];
    GridFunctionCL<Point3DCL> qsurfP1grad[4], qsurfP2grad[10];
    GridFunctionCL<> qP1Hat[4], qP2Hat[10], qconvP1[4], qdivP1, qnormal_comp[3], qvProjP1_comp[3], qvP1_comp[3];


    QuadDomain2DCL  q2Ddomain;
    std::valarray<double> ls_loc;
    SurfacePatchCL spatch;

    QuadDomainCL q3Ddomain;
    GridFunctionCL<Point3DCL> q3Dnormal, q3DexactNormal, q3DP2Grad[10];
    GridFunctionCL<Point3DCL> q3DP1Grad[4], q2DP1Grad[4];

    bool fullGrad;

    void Get_Normals(const LocalP2CL<>& ls, LocalP1CL<Point3DCL>&);

  public:
    LocalStokesCL ( bool fullGradient)
        : lat( PrincipalLatticeCL::instance( 2)), ls_loc( lat.vertex_size()), fullGrad(fullGradient)
    {
        P1DiscCL::GetP1Basis( P1Hat);
        P2DiscCL::GetP2Basis( P2Hat);
        P2DiscCL::GetGradientsOnRef( P2GradRef);
        std::cout<<"full gradient="<<fullGrad<<std::endl;
    }

    void calcIntegrands(const SMatrixCL<3,3>& T, const LocalP2CL<>& ls, const TetraCL& tet); ///< has to be called before any setup method!
    void calcIntegrands(const SMatrixCL<3,3>& T, const LocalP2CL<>& ls, const LocalP1CL<Point3DCL>& v, const TetraCL& tet); // for conv term
    void calc3DIntegrands(const SMatrixCL<3,3>& T, const LocalP2CL<>& ls, const TetraCL& tet); ///< has to be called after calcIntegrands!
    void setupA_P2 (double A_P2[30][30]);
    void setupA_P2_stab (double A_P2_stab[10][10], double absdet);
    void setupA_P1 (double A_P1[12][12]);
    void setupA_P1_stab (double A_P1_stab[4][4], double absdet);
    void setupB_P1P2 (double B_P1P2[4][30]);
    void setupB_P1P1 (double B_P1P1[4][12]);
    void setupB_P2P2 (double B_P2P2[10][30]);
    void setupB_P2P1 (double B_P2P1[10][12]);
    void setupM_P2 (double M_P2[10][10]);
    void setupM_P1 (double M_P1[4][4]);
    void setupS_P2 (double S_P2[30][30]);
    void setupS_P1 (double S_P1[12][12]);
    void setupL_P1P1 (double L_P1P1[4][12]);
    void setupL_P1P1_stab (double L_P1P1_stab[4][12], double absdet);
    void setupL_P1P2 (double L_P1P2[4][30]);
    void setupL_P1P2_stab (double L_P1P2_stab[4][30], double absdet);
    void setupL_P2P1 (double L_P2P1[10][12]);
    void setupL_P2P1_stab (double L_P2P1_stab[10][12], double absdet);
    void setupL_P2P2 (double L_P2P2[10][30]);
    void setupL_P2P2_stab (double L_P2P2_stab[10][30], double absdet);
    // Navier-Stokes:
    void setupN_P1 (double N_P1[12][12]);
    void setupNT_P1 (double NT_P1[12][12]);
    void setupOmega_P1P1 (double Omega_P1P1[4][12]);
    void setupD_P1 (double D_P1[4][4]);
    void setupM_P1_stab (double M_P1_stab[4][4], double absdet);
};

DROPS::Point3DCL Normal_sphere2 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL v(p[0]/(std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])),p[1]/(std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])),p[2]/(std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])));
    //DROPS::Point3DCL w(p[0], 0., 0.);
    return v;
}

// The P2 levelset-function is used to compute the normals which are needed for the (improved) projection onto the interface, GradLP1 has to be set before
void LocalStokesCL::Get_Normals(const LocalP2CL<>& ls, LocalP1CL<Point3DCL>& Normals)
{
    for(int i=0; i<10 ; ++i)
    {
        Normals+=ls[i]*P2Grad[i];
    }

}

void LocalStokesCL::calc3DIntegrands(const SMatrixCL<3,3>& T, const LocalP2CL<>& ls, const TetraCL& tet)
{
    make_SimpleQuadDomain<Quad5DataCL> (q3Ddomain, AllTetraC);
    LocalP1CL<Point3DCL> Normals;
    Get_Normals(ls, Normals);

    resize_and_evaluate_on_vertexes (Normals, q3Ddomain, q3Dnormal);
    // Scale Normals accordingly to the Euclidean Norm (only consider the ones which make a contribution in the sense of them being big enough... otherwise one has to expect problems with division through small numbers)
    for(Uint i=0; i<q3Dnormal.size(); ++i) {
         //if(qnormal[i].norm()> 1e-8)
         q3Dnormal[i]= q3Dnormal[i]/q3Dnormal[i].norm();
    }

    resize_and_evaluate_on_vertexes (&Normal_sphere2, tet, q3Ddomain, 0., q3DexactNormal);

    for(int j=0; j<10; ++j) {
        resize_and_evaluate_on_vertexes(P2Grad[j], q3Ddomain, q3DP2Grad[j]);
    }

    for(int j=0; j<4; ++j) {
        q3DP1Grad[j].resize( q3Ddomain.vertex_size());
        q3DP1Grad[j]= P1Grad[j];
    }
}

void LocalStokesCL::calcIntegrands(const SMatrixCL<3,3>& T, const LocalP2CL<>& ls, const TetraCL& tet)
{
    P2DiscCL::GetGradients( P2Grad, P2GradRef, T);
    P1DiscCL::GetGradients( P1Grad, T);
    evaluate_on_vertexes( ls, lat, Addr( ls_loc));
    spatch.make_patch<MergeCutPolicyCL>( lat, ls_loc);
    // The routine takes the information about the tetrahedra and the cutting surface and generates a two-dimensional triangulation of the cut, including the necessary point-positions and weights for the quadrature
    make_CompositeQuad5Domain2D ( q2Ddomain, spatch, tet);
    LocalP1CL<Point3DCL> Normals;
    Get_Normals(ls, Normals);
    // Resize and evaluate Normals at all points which are needed for the two-dimensional quadrature-rule
    resize_and_evaluate_on_vertexes (Normals, q2Ddomain, qnormal);
    // Scale Normals accordingly to the Euclidean Norm (only consider the ones which make a contribution in the sense of them being big enough... otherwise one has to expect problems with division through small numbers)
    for(Uint i=0; i<qnormal.size(); ++i) {
         //if(qnormal[i].norm()> 1e-8)
         qnormal[i]= qnormal[i]/qnormal[i].norm();
    }
    // Provide all components of the normals
    for (int k=0; k<3; ++k) {
        qnormal_comp[k].resize( q2Ddomain.vertex_size());
        ExtractComponent( qnormal, qnormal_comp[k], k);
    }

    // Resize and evaluate P2 basis functions
    for(int j=0; j<10 ;++j) {
        resize_and_evaluate_on_vertexes( P2Hat[j], q2Ddomain, qP2Hat[j]);
        resize_and_evaluate_on_vertexes( P2Grad[j], q2Ddomain, qsurfP2grad[j]);
//        qsurfP2grad[j].resize( q2Ddomain.vertex_size());
//        qsurfP2grad[j]= P2Grad[j];
        qsurfP2grad[j]-= dot( qsurfP2grad[j], qnormal)*qnormal;
    }
    // Resize and evaluate of all the 4 P1 Gradient Functions and apply pointwise projections   (P grad \xi_j  for j=1..4)
    for(int j=0; j<4 ;++j) {
//        resize_and_evaluate_on_vertexes( P1Grad[j], q2Ddomain, qsurfP1grad[j]);
        resize_and_evaluate_on_vertexes( P1Hat[j], q2Ddomain, qP1Hat[j]);
        qsurfP1grad[j].resize( q2Ddomain.vertex_size());
        qsurfP1grad[j]= P1Grad[j];
        qsurfP1grad[j]-= dot( qsurfP1grad[j], qnormal)*qnormal;
    }
    // Apply pointwise projection to the 3 std basis vectors e_k
    for(int k=0; k<3 ;++k) {
        //resize_and_evaluate_on_vertexes( P1Grad[j], q2Ddomain, qsurfP1grad[j]);
        qProj[k].resize( q2Ddomain.vertex_size());
        Point3DCL e_k= DROPS::std_basis<3>( k+1);
        qProj[k]= e_k;
        qProj[k]-= dot( e_k, qnormal)*qnormal;
    }
}

// for convective term
void LocalStokesCL::calcIntegrands(const SMatrixCL<3,3>& T, const LocalP2CL<>& ls, const LocalP1CL<Point3DCL>& v , const TetraCL& tet)
{

    P2DiscCL::GetGradients( P2Grad, P2GradRef, T);
    P1DiscCL::GetGradients( P1Grad, T);
    evaluate_on_vertexes( ls, lat, Addr( ls_loc));
    spatch.make_patch<MergeCutPolicyCL>( lat, ls_loc);
    // The routine takes the information about the tetrahedra and the cutting surface and generates a two-dimensional triangulation of the cut, including the necessary point-positions and weights for the quadrature
    make_CompositeQuad5Domain2D ( q2Ddomain, spatch, tet);
    LocalP1CL<Point3DCL> Normals;
    Get_Normals(ls, Normals);
    // Resize and evaluate Normals at all points which are needed for the two-dimensional quadrature-rule

    /////////////////Y
    //resize_and_evaluate_on_vertexes (&Normal_sphere2, tet, q2Ddomain,0., qnormal);
    resize_and_evaluate_on_vertexes (Normals, q2Ddomain, qnormal);
    /////////////////

    for(int j=0; j<4; ++j) {
                q2DP1Grad[j].resize( q2Ddomain.vertex_size());
                q2DP1Grad[j]= P1Grad[j];
            }


    // Scale Normals accordingly to the Euclidean Norm (only consider the ones which make a contribution in the sense of them being big enough... otherwise one has to expect problems with division through small numbers)
    for(Uint i=0; i<qnormal.size(); ++i) {
         //if(qnormal[i].norm()> 1e-8)
         qnormal[i]= qnormal[i]/qnormal[i].norm();
    }
    // Provide all components of the normals
    for (int k=0; k<3; ++k) {
        qnormal_comp[k].resize( q2Ddomain.vertex_size());
        ExtractComponent( qnormal, qnormal_comp[k], k);
    }

    // Resize and evaluate P2 basis functions
    for(int j=0; j<10 ;++j) {
        resize_and_evaluate_on_vertexes( P2Hat[j], q2Ddomain, qP2Hat[j]);
        resize_and_evaluate_on_vertexes( P2Grad[j], q2Ddomain, qsurfP2grad[j]);
//        qsurfP2grad[j].resize( q2Ddomain.vertex_size());
//        qsurfP2grad[j]= P2Grad[j];
        qsurfP2grad[j]-= dot( qsurfP2grad[j], qnormal)*qnormal;
    }
    // Resize and evaluate of all the 4 P1 Gradient Functions and apply pointwise projections   (P grad \xi_j  for j=1..4)

    for(int j=0; j<4 ;++j) {
//        resize_and_evaluate_on_vertexes( P1Grad[j], q2Ddomain, qsurfP1grad[j]);
        resize_and_evaluate_on_vertexes( P1Hat[j], q2Ddomain, qP1Hat[j]);
        qsurfP1grad[j].resize( q2Ddomain.vertex_size());
        qsurfP1grad[j]= P1Grad[j];
        qsurfP1grad[j]-= dot( qsurfP1grad[j], qnormal)*qnormal;

    }
    // Apply pointwise projection to the 3 std basis vectors e_k
    for(int k=0; k<3 ;++k) {
        //resize_and_evaluate_on_vertexes( P1Grad[j], q2Ddomain, qsurfP1grad[j]);
        qProj[k].resize( q2Ddomain.vertex_size());
        Point3DCL e_k= DROPS::std_basis<3>( k+1);
        qProj[k]= e_k;

        qProj[k]-= dot( e_k, qnormal)*qnormal;

    }

    //physical velocity, tetra P1 function
    LocalP1CL<Point3DCL> velocity;
    for(int i=0; i<4 ; ++i)
    {
    	velocity += v[i]*P1Hat[i];
    }

    //restrict to a tri(quad) element or "trace" it
    resize_and_evaluate_on_vertexes (velocity, q2Ddomain, qvP1);


    //get tangential part of tri P1 function
    qvProjP1.resize( q2Ddomain.vertex_size());
    qvProjP1 = qvP1;
    qvProjP1 -= dot(qvP1, qnormal)*qnormal;


    //P1 function Projected discrete velocity
    GridFunctionCL<Point3DCL> qvdiscP1;
    LocalP1CL<> v_comp;

 	for(int n=0; n<3 ; ++n)
 	{
 		ExtractComponent( v, v_comp, n);
 		qvdiscP1 += v_comp*qProj[n];

 	}
    resize_and_evaluate_on_vertexes (velocity, q2Ddomain, qvdiscP1);


    qdivP1.resize( q2Ddomain.vertex_size());
    for(int m=0; m<4 ; ++m)
    {
    	qdivP1 += dot(qvdiscP1[m],qsurfP1grad[m]);
    }



    //extract components of tangential velocity as tri P1 functions
    for (int k=0; k<3; ++k) {
    	 qvProjP1_comp[k].resize( q2Ddomain.vertex_size());
          ExtractComponent( qvProjP1, qvProjP1_comp[k], k);
      }

   /* std::cout << "min|v|     : " <<std::sqrt( dot(v,v).min() ) << std::endl;

    std::cout << "max|v|     : " <<std::sqrt( dot(v,v).max() ) << std::endl;


    std::cout << "|velocity| : " <<std::sqrt( dot(velocity,velocity).max() ) << std::endl;
*/
   /* std::cout << "|qvP1|     : " <<std::sqrt( dot(qvP1,qvP1).max() ) << std::endl;
    std::cout << "|qvProjP1|     : " <<std::sqrt( dot(qvProjP1,qvProjP1).max() ) << std::endl;
    std::cout << "|qvProjP1_comp|     : " <<std::sqrt((qvProjP1_comp[0]*qvProjP1_comp[0]+qvProjP1_comp[1]*qvProjP1_comp[1]+qvProjP1_comp[2]*qvProjP1_comp[2]).max() ) << std::endl;
*/

    //for convective nonlinear term
    for(int i=0; i<4 ;++i) {
    	qconvP1[i].resize( q2Ddomain.vertex_size());
    	qconvP1[i] = dot( qsurfP1grad[i], qvP1)  ;
   }


    //for calculation curl
   for (int i=0; i<4; i++)
    {
    	for (int k=0; k<3; k++)
    	{
    		GridFunctionCL<Point3DCL>  qlocal;
    		//LocalP1CL<Point3DCL> qlocal;
    		qlocal.resize(q2Ddomain.vertex_size());

    		for (int index=0; index<3; index++)
    		{
    			if (k==index) continue;
    			int right_index;
    			int multiplier;
    			if (index==0)
    			{
    				if (k==1)
    				{
    					right_index=2;
    					multiplier=-1;
    				}
    				else if (k==2)
    				{
    					right_index=1;
    					multiplier=1;
    				}
    			}
    			else if (index==1)
    			{
    				if (k==0)
    				{
    					right_index=2;
    					multiplier=1;
    				}
    				else if (k==2)
    				{
    					right_index=0;
    					multiplier=-1;
    				}
    			}
    			else if (index==2)
    			{
    				if (k==0)
    				{
    					right_index=1;
    					multiplier=-1;
    				}
    				else if (k==1)
    				{
    					right_index=0;
    					multiplier=1;
    				}
    			}


    			GridFunctionCL<double>  ql;
    			ql.resize( q2Ddomain.vertex_size());

    			ExtractComponent( q2DP1Grad[i], ql, right_index);


    			Point3DCL e= DROPS::std_basis<3>( index+1);
    			GridFunctionCL<Point3DCL> ee ;
				ee.resize( q2Ddomain.vertex_size());
				ee= multiplier*e;
    			qlocal += ql*ee;

    		}

    		qrotP1[3*i+k].resize( q2Ddomain.vertex_size());
    		//resize_and_evaluate_on_vertexes (qlocal, q2Ddomain, qrotP1[3*i+k]);
    		qrotP1[3*i+k]=qlocal;
   		 //std::cout << "|qrotP1|     : " <<std::sqrt( dot(qrotP1[3*i+k],qrotP1[3*i+k]).max() ) << std::endl;

    	}
    }
}

// TODO: Den Fall fullGrad testen!
void LocalStokesCL::setupA_P2 (double A_P2[30][30])
{
    // Do all combinations for (i,j) i,j=30 x 30 and corresponding quadrature
    for (int i=0; i < 10; ++i) {
        for (int j=0; j<10; ++j) {
            for (int k=0; k<3; ++k) {
                Point3DCL e_k= DROPS::std_basis<3>( k+1);
                for (int l=0; l<3; ++l) {
                    Point3DCL e_l= DROPS::std_basis<3>( l+1);
                    if (!fullGrad) {
                        A_P2[3*i+k][3*j+l]= quad_2D( dot(e_k, qsurfP2grad[j])*dot(e_l, qsurfP2grad[i]) + dot(qsurfP2grad[i], qsurfP2grad[j])*dot(qProj[k], qProj[l]), q2Ddomain);
                    }
                    else {
                        if(k==l) {
                            A_P2[3*i+k][3*j+l]= quad_2D( dot(e_k, qsurfP2grad[j])*dot(e_l, qsurfP2grad[i]) + 0.5*dot(qsurfP2grad[i], qsurfP2grad[j]) + 0.5*dot(P2Grad[i], P2Grad[j])*dot(qProj[k], qProj[l]), q2Ddomain);
                        }
                        else {
                            A_P2[3*i+k][3*j+l]= quad_2D( dot(e_k, qsurfP2grad[j])*dot(e_l, qsurfP2grad[i]) + 0.5*dot(P2Grad[i], P2Grad[j])*dot(qProj[k], qProj[l]), q2Ddomain);
                        }
                    }
                }
            }
        }
    }
}

// TODO: Den Fall fullGrad testen!
void LocalStokesCL::setupA_P1 (double A_P1[12][12])
{
    // Do all combinations for (i,j) i,j=30 x 30 and corresponding quadrature
    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
            for (int k=0; k<3; ++k) {
                Point3DCL e_k= DROPS::std_basis<3>( k+1);
                for (int l=0; l<3; ++l) {
                    Point3DCL e_l= DROPS::std_basis<3>( l+1);
                    if (!fullGrad) {
                        A_P1[3*i+k][3*j+l]= quad_2D( dot(e_k, qsurfP1grad[j])*dot(e_l, qsurfP1grad[i]) + dot(qsurfP1grad[i], qsurfP1grad[j])*dot(qProj[k], qProj[l]), q2Ddomain);
                    }
                    else {
                        if(k==l) {
                            A_P1[3*i+k][3*j+l]= quad_2D( dot(e_k, qsurfP1grad[j])*dot(e_l, qsurfP1grad[i]) + 0.5*dot(qsurfP1grad[i], qsurfP1grad[j]) + 0.5*inner_prod(P1Grad[i], P1Grad[j])*dot(qProj[k], qProj[l]), q2Ddomain);
                        }
                        else {
                            A_P1[3*i+k][3*j+l]= quad_2D( dot(e_k, qsurfP1grad[j])*dot(e_l, qsurfP1grad[i]) + 0.5*inner_prod(P1Grad[i], P1Grad[j])*dot(qProj[k], qProj[l]), q2Ddomain);
                        }
                    }
                }
            }
        }
    }
}

// TODO: Den Fall fullGrad testen!
void LocalStokesCL::setupB_P1P2 (double B_P1P2[4][30])
{
    GridFunctionCL<double> qsurfgrad_k;
    qsurfgrad_k.resize( q2Ddomain.vertex_size());

    // Do all combinations for (i,j) i,j=4 x 30 and corresponding quadrature
    for (int k=0; k<3; ++k) {
        for (int i=0; i < 4; ++i) {
            ExtractComponent( qsurfP1grad[i], qsurfgrad_k, k);
            for (int j=0; j<10; ++j) {
                if (!fullGrad) {
                    B_P1P2[i][3*j+k]= quad_2D( qsurfgrad_k*qP2Hat[j], q2Ddomain);
                }
                else
                    B_P1P2[i][3*j+k]= quad_2D( P1Grad[i][k]*qP2Hat[j], q2Ddomain);
            }
        }
    }
}

// TODO: Den Fall fullGrad testen!
void LocalStokesCL::setupB_P2P2 (double B_P2P2[10][30])
{
    GridFunctionCL<double> qsurfgrad_k;
    qsurfgrad_k.resize( q2Ddomain.vertex_size());

    // Do all combinations for (i,j) i,j=4 x 30 and corresponding quadrature
    for (int k=0; k<3; ++k) {
        for (int i=0; i < 10; ++i) {
            ExtractComponent( qsurfP2grad[i], qsurfgrad_k, k);
            for (int j=0; j<10; ++j) {
                if (!fullGrad) {
                    B_P2P2[i][3*j+k]= quad_2D( qsurfgrad_k*qP2Hat[j], q2Ddomain);
                }
                else
                    std::cout << "DOES NOT WORK YET!!!" << std::endl;//B_P2P2[i][3*j+k]= quad_2D( P2Grad[i][k]*qP2Hat[j], q2Ddomain);
            }
        }
    }
}

// TODO: Den Fall fullGrad testen!
void LocalStokesCL::setupB_P1P1 (double B_P1P1[4][12])
{
    GridFunctionCL<double> qsurfgrad_k;
    qsurfgrad_k.resize( q2Ddomain.vertex_size());

    // Do all combinations for (i,j) i,j=4 x 30 and corresponding quadrature
    for (int k=0; k<3; ++k) {
        for (int i=0; i < 4; ++i) {
            ExtractComponent( qsurfP1grad[i], qsurfgrad_k, k);
            for (int j=0; j<4; ++j) {
                if (!fullGrad) {
                    B_P1P1[i][3*j+k]= quad_2D( qsurfgrad_k*qP1Hat[j], q2Ddomain);
                }
                else
                    B_P1P1[i][3*j+k]= quad_2D( P1Grad[i][k]*qP1Hat[j], q2Ddomain);
            }
        }
    }
}

// TODO: Den Fall fullGrad testen!
void LocalStokesCL::setupB_P2P1 (double B_P2P1[10][12])
{
    GridFunctionCL<double> qsurfgrad_k;
    qsurfgrad_k.resize( q2Ddomain.vertex_size());

    // Do all combinations for (i,j) i,j=4 x 30 and corresponding quadrature
    for (int k=0; k<3; ++k) {
        for (int i=0; i < 10; ++i) {
            ExtractComponent( qsurfP2grad[i], qsurfgrad_k, k);
            for (int j=0; j<4; ++j) {
                if (!fullGrad) {
                    B_P2P1[i][3*j+k]= quad_2D( qsurfgrad_k*qP1Hat[j], q2Ddomain);
                }
                else
                    std::cout << "DOES NOT WORK YET!!!" << std::endl;//B_P2P1[i][3*j+k]= quad_2D( P2Grad[i][k]*qP1Hat[j], q2Ddomain);
            }
        }
    }
}

void LocalStokesCL::setupM_P2 (double M_P2[10][10])
{
    // Do all combinations for (i,j) i,j=4 x 30 and corresponding quadrature
    for (int i=0; i < 10; ++i) {
        for (int j=0; j<10; ++j) {
             M_P2[i][j]= quad_2D( qP2Hat[i]*qP2Hat[j], q2Ddomain);
        }
    }
}

void LocalStokesCL::setupM_P1 (double M_P1[4][4])
{
    // Do all combinations for (i,j) i,j=4 x 4 and corresponding quadrature
    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
             M_P1[i][j]= quad_2D( qP1Hat[i]*qP1Hat[j], q2Ddomain);
        }
    }
}

void LocalStokesCL::setupS_P2 (double S_P2[30][30])
{
    // Do all combinations for (i,j) i,j=30 x 30 and corresponding quadrature
    for (int i=0; i < 10; ++i) {
        for (int j=0; j<10; ++j) {
            for (int k=0; k<3; ++k) {
                for (int l=0; l<3; ++l) {
                     S_P2[3*i+k][3*j+l]= quad_2D( qP2Hat[i]*qnormal_comp[k]*qP2Hat[j]*qnormal_comp[l], q2Ddomain);
                }
            }
        }
    }
}

void LocalStokesCL::setupS_P1 (double S_P1[12][12])
{
    // Do all combinations for (i,j) i,j=30 x 30 and corresponding quadrature
    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
            for (int k=0; k<3; ++k) {
                for (int l=0; l<3; ++l) {
                     S_P1[3*i+k][3*j+l]= quad_2D( qP1Hat[i]*qnormal_comp[k]*qP1Hat[j]*qnormal_comp[l], q2Ddomain);
                }
            }
        }
    }
}

void LocalStokesCL::setupL_P1P2 (double L_P1P2[4][30])
{
    //GridFunctionCL<double> qP1Hat;

    // Do all combinations for (i,j) i,j=4 x 30 and corresponding quadrature
    for (int i=0; i < 4; ++i) {
       // resize_and_evaluate_on_vertexes( P1Hat[i], q2Ddomain, qP1Hat);
        for (int j=0; j<10; ++j) {
            for (int k=0; k<3; ++k) {
                L_P1P2[i][3*j+k]= quad_2D( qP1Hat[i]*qP2Hat[j]*qnormal_comp[k], q2Ddomain);
            }
        }
    }
}

void LocalStokesCL::setupL_P1P2_stab (double L_P1P2_stab[4][30], double absdet)
{
    // Do all combinations for (i,j) i,j=4 x 30 and corresponding quadrature
    for (int i=0; i < 10; ++i) {
        for (int j=0; j<3; ++j) {
            Point3DCL e_j= DROPS::std_basis<3>( j+1);
            for (int k=0; k<4; ++k) {
                L_P1P2_stab[k][3*i+j]= quad(dot(q3Dnormal, q3DP2Grad[i])*dot(e_j, q3Dnormal)*dot(q3Dnormal, q3DP1Grad[k]), absdet, q3Ddomain, AllTetraC);
            }
        }
    }
}

void LocalStokesCL::setupL_P1P1 (double L_P1P1[4][12])
{
   // GridFunctionCL<double> qP1Hat;

    // Do all combinations for (i,j) i,j=4 x 30 and corresponding quadrature
    for (int i=0; i < 4; ++i) {
      //  resize_and_evaluate_on_vertexes( P1Hat[i], q2Ddomain, qP1Hat);
        for (int j=0; j<4; ++j) {
            for (int k=0; k<3; ++k) {
                L_P1P1[i][3*j+k]= quad_2D( qP1Hat[i]*qP1Hat[j]*qnormal_comp[k], q2Ddomain);
            }
        }
    }
}

void LocalStokesCL::setupL_P1P1_stab (double L_P1P1_stab[4][12], double absdet)
{
    // Do all combinations for (i,j) i,j=4 x 30 and corresponding quadrature
    for (int i=0; i < 4; ++i) {
        for (int j=0; j<3; ++j) {
            Point3DCL e_j= DROPS::std_basis<3>( j+1);
            for (int k=0; k<4; ++k) {
                L_P1P1_stab[k][3*i+j]= quad(dot(q3Dnormal, q3DP1Grad[i])*dot(e_j, q3Dnormal)*dot(q3Dnormal, q3DP1Grad[k]), absdet, q3Ddomain, AllTetraC);
            }
        }
    }
}

void LocalStokesCL::setupL_P2P1 (double L_P2P1[10][12])
{
    // Do all combinations for (i,j) i,j=10 x 30 and corresponding quadrature
    for (int i=0; i < 10; ++i) {
        for (int j=0; j<4; ++j) {
            for (int k=0; k<3; ++k) {
                L_P2P1[i][3*j+k]= quad_2D( qP2Hat[i]*qP1Hat[j]*qnormal_comp[k], q2Ddomain);
            }
        }
    }
}

void LocalStokesCL::setupL_P2P1_stab (double L_P2P1_stab[10][12], double absdet)
{
    // Do all combinations for (i,j) i,j=4 x 30 and corresponding quadrature
    for (int i=0; i < 4; ++i) {
        for (int j=0; j<3; ++j) {
            Point3DCL e_j= DROPS::std_basis<3>( j+1);
            for (int k=0; k<10; ++k) {
                L_P2P1_stab[k][3*i+j]= quad(dot(q3Dnormal, q3DP1Grad[i])*dot(e_j, q3Dnormal)*dot(q3Dnormal, q3DP2Grad[k]), absdet, q3Ddomain, AllTetraC);
            }
        }
    }
}

void LocalStokesCL::setupL_P2P2 (double L_P2P2[10][30])
{
    // Do all combinations for (i,j) i,j=10 x 30 and corresponding quadrature
    for (int i=0; i < 10; ++i) {
        for (int j=0; j<10; ++j) {
            for (int k=0; k<3; ++k) {
                L_P2P2[i][3*j+k]= quad_2D( qP2Hat[i]*qP2Hat[j]*qnormal_comp[k], q2Ddomain);
            }
        }
    }
}

void LocalStokesCL::setupL_P2P2_stab (double L_P2P2_stab[10][30], double absdet)
{
    // Do all combinations for (i,j) i,j=4 x 30 and corresponding quadrature
    for (int i=0; i < 10; ++i) {
        for (int j=0; j<3; ++j) {
            Point3DCL e_j= DROPS::std_basis<3>( j+1);
            for (int k=0; k<10; ++k) {
                L_P2P2_stab[k][3*i+j]= quad(dot(q3Dnormal, q3DP2Grad[i])*dot(e_j, q3Dnormal)*dot(q3Dnormal, q3DP2Grad[k]), absdet, q3Ddomain, AllTetraC);
            }
        }
    }
}

void LocalStokesCL::setupA_P2_stab (double A_P2_stab[10][10], double absdet)
{
    for (int i=0; i<10; ++i) {
        for (int j=0; j<10; ++j) {
           A_P2_stab[i][j] = quad(dot(q3Dnormal, q3DP2Grad[i])*dot(q3Dnormal, q3DP2Grad[j]), absdet, q3Ddomain, AllTetraC);
        }
    }
}

void LocalStokesCL::setupA_P1_stab (double A_P1_stab[4][4], double absdet)
{
    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
           A_P1_stab[i][j] = quad(dot(q3Dnormal, q3DP1Grad[i])*dot(q3Dnormal, q3DP1Grad[j]), absdet, q3Ddomain, AllTetraC);
        }
    }
}

// Navier-Stokes:

void LocalStokesCL::setupN_P1 (double N_P1[12][12])
{
    // Do all combinations for (i,j) i,j=12 x 12 and corresponding quadrature
    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
            for (int k=0; k<3; ++k) {
                for (int l=0; l<3; ++l) {
                    N_P1[3*i+k][3*j+l]= quad_2D( qP1Hat[j]
											   * dot(qProj[k], qProj[l])*qconvP1[i], q2Ddomain);
                }
            }
        }
    }
}

void LocalStokesCL::setupNT_P1 (double NT_P1[12][12])
{
    // Do all combinations for (i,j) i,j=12 x 12 and corresponding quadrature
    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
            for (int k=0; k<3; ++k) {
                for (int l=0; l<3; ++l) {
                	Point3DCL e_l= DROPS::std_basis<3>( l+1);

                	NT_P1[3*i+k][3*j+l] = quad_2D(qP1Hat[j]*
                			dot( e_l, qsurfP1grad[i])
							*qvProjP1_comp[k]
										   , q2Ddomain);
                }
            }
        }
    }
}

void LocalStokesCL::setupOmega_P1P1 (double Omega_P1P1[4][12])
{
	//GridFunctionCL<double> qsurfgrad_k;


    // Do all combinations for (i,j) i,j=4 x 30 and corresponding quadrature
    for (int k=0; k<3; ++k) {
        for (int i=0; i < 4; ++i) {

            for (int m=0; m<4; ++m) {
              	//Omega_P1P1[m][3*i+k]= quad_2D( dot( qrotP1[3*i+k], qrotP1[3*i+k])*qP1Hat[m], q2Ddomain);

            	Omega_P1P1[m][3*i+k]= quad_2D( dot( qnormal, qrotP1[3*i+k])*qP1Hat[m], q2Ddomain);
            }
        }
    }
}

void LocalStokesCL::setupD_P1(double D_P1[4][4])
{

    // Do all combinations for (i,j) i,j=4 x 4 and corresponding quadrature
    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
        	/*for (int k=0; k<3; ++k) {
        		for (int l=0; l<3; ++l) {
        			if (k==l)
        			{*/

        				D_P1[i][j]= quad_2D( qdivP1*
        						qP1Hat[i]*qP1Hat[j]  , q2Ddomain);


        	/*		}
        		}
        	}*/
        }
    }


}

void LocalStokesCL::setupM_P1_stab (double M_P1_stab[4][4], double absdet)
{
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

    LocalStokesCL localStokes_;

    SMatrixCL<3,3> T;
    double det, absdet;
    LocalP2CL<> ls_loc;
    LocalNumbP1CL n, nScalar;

    ///\brief Computes the mapping from local to global data "n", the local matrices in loc and.
    void local_setup (const TetraCL& tet);
    ///\brief Update the global system.
    void update_global_system ();

  public:
    StokesIFAccumulator_P1P1CL ( const VecDescCL& ls, const LsetBndDataCL& ls_bnd, IdxDescCL& P1FE, IdxDescCL& ScalarP1FE, MatrixCL& A_P1, MatrixCL& A_P1_stab, MatrixCL& B_P1P1, MatrixCL& M_P1, MatrixCL& S_P1, MatrixCL& L_P1P1, MatrixCL& L_P1P1_stab, MatrixCL& M_ScalarP1, MatrixCL& A_ScalarP1_stab, bool fullGradient);

    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    void visit (const TetraCL& sit);

    TetraAccumulatorCL* clone (int /*tid*/) { return new StokesIFAccumulator_P1P1CL ( *this); }

};

StokesIFAccumulator_P1P1CL::StokesIFAccumulator_P1P1CL( const VecDescCL& ls, const LsetBndDataCL& ls_bnd, IdxDescCL& P1FE, IdxDescCL& ScalarP1FE, MatrixCL& A_P1, MatrixCL& A_P1_stab, MatrixCL& B_P1P1, MatrixCL& M_P1, MatrixCL& S_P1, MatrixCL& L_P1P1, MatrixCL& L_P1P1_stab, MatrixCL& M_ScalarP1, MatrixCL& A_ScalarP1_stab, bool fullGradient)
 : lset(ls), lset_bnd(ls_bnd), P1Idx_(P1FE), ScalarP1Idx_(ScalarP1FE), A_P1_(A_P1), A_P1_stab_(A_P1_stab), B_P1P1_(B_P1P1), M_P1_(M_P1), S_P1_(S_P1), L_P1P1_(L_P1P1), L_P1P1_stab_(L_P1P1_stab), M_ScalarP1_(M_ScalarP1), A_ScalarP1_stab_(A_ScalarP1_stab), localStokes_( fullGradient)
{}

void StokesIFAccumulator_P1P1CL::begin_accumulation ()
{
    std::cout << "entering StokesIF: \n";
    const size_t num_unks_p1= P1Idx_.NumUnknowns(), num_unks_scalarp1= ScalarP1Idx_.NumUnknowns();
    mA_P1_= new MatrixBuilderCL( &A_P1_, num_unks_p1, num_unks_p1);
    mA_P1_stab_= new MatrixBuilderCL( &A_P1_stab_, num_unks_p1, num_unks_p1);
    mB_P1P1_= new MatrixBuilderCL( &B_P1P1_, num_unks_scalarp1, num_unks_p1);
    mM_P1_= new MatrixBuilderCL( &M_P1_, num_unks_p1, num_unks_p1);
    mS_P1_= new MatrixBuilderCL( &S_P1_, num_unks_p1, num_unks_p1);
    mL_P1P1_= new MatrixBuilderCL( &L_P1P1_, num_unks_scalarp1, num_unks_p1);
    mL_P1P1_stab_= new MatrixBuilderCL( &L_P1P1_stab_, num_unks_scalarp1, num_unks_p1);
    mM_ScalarP1_= new MatrixBuilderCL( &M_ScalarP1_, num_unks_scalarp1, num_unks_scalarp1);
    mA_ScalarP1_stab_= new MatrixBuilderCL( &A_ScalarP1_stab_, num_unks_scalarp1, num_unks_scalarp1);
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
    ls_loc.assign( tet, lset, lset_bnd);

    if (!equal_signs( ls_loc))
    {
        local_setup( tet);
        update_global_system();
    }

}

void StokesIFAccumulator_P1P1CL::local_setup (const TetraCL& tet)
{
    GetTrafoTr( T, det, tet);
    absdet= std::fabs( det);

    n.assign( tet, P1Idx_, P1Idx_.GetBndInfo());
    nScalar.assign( tet, ScalarP1Idx_, ScalarP1Idx_.GetBndInfo());

    //GetLocalNumbP1NoBnd( numP1, tet, P1Idx_);
    //GetLocalNumbP1NoBnd( numScalarP1, tet, ScalarP1Idx_);
    localStokes_.calcIntegrands( T, ls_loc, tet);
    localStokes_.calc3DIntegrands(T, ls_loc, tet);
    localStokes_.setupA_P1( locA_P1);
    localStokes_.setupA_P1_stab( locA_P1_stab, absdet);
    localStokes_.setupB_P1P1( locB_P1P1);
    localStokes_.setupM_P1( locM_P1);
    localStokes_.setupS_P1( locS_P1);
    localStokes_.setupL_P1P1( locL_P1P1);
    localStokes_.setupL_P1P1_stab( locL_P1P1_stab, absdet);
}

void StokesIFAccumulator_P1P1CL::update_global_system ()
{
    MatrixBuilderCL& mA_P1= *mA_P1_;
    MatrixBuilderCL& mA_P1_stab= *mA_P1_stab_;
    MatrixBuilderCL& mB_P1P1= *mB_P1P1_;
    MatrixBuilderCL& mM_P1= *mM_P1_;
    MatrixBuilderCL& mS_P1= *mS_P1_;
    MatrixBuilderCL& mL_P1P1= *mL_P1P1_;
    MatrixBuilderCL& mL_P1P1_stab= *mL_P1P1_stab_;
    MatrixBuilderCL& mM_ScalarP1= *mM_ScalarP1_;
    MatrixBuilderCL& mA_ScalarP1_stab= *mA_ScalarP1_stab_;

    for(int i= 0; i < 4; ++i) {
        const IdxT ii= n.num[i];
        if (ii==NoIdx || !(n.WithUnknowns( i))) continue;
        for(int j=0; j < 4; ++j) {
            const IdxT jj= n.num[j];
            if (jj==NoIdx || !(n.WithUnknowns( j))) continue;
            for (int k=0; k<3; ++k) {
                for (int l=0; l<3; ++l) {
                    mA_P1( ii+k, jj+l) += locA_P1[3*j+l][3*i+k];
                    mS_P1( ii+k, jj+l) += locS_P1[3*j+l][3*i+k];
                    if(k == l) {
                        mM_P1( ii+k, jj+l) += locM_P1[j][i];
                        mA_P1_stab( ii+k, jj+l) += locA_P1_stab[j][i];
                    } else {
                        mM_P1( ii+k, jj+l) += 0.;
                        mA_P1_stab( ii+k, jj+l) += 0.;
                    }
                }
             }
        }
        for (int j=0; j < 4; ++j) {
            if (nScalar.num[j]==NoIdx || !(nScalar.WithUnknowns( j))) continue;
            for (int k=0; k<3; ++k) {
                mB_P1P1( nScalar.num[j], ii+k) += locB_P1P1[j][3*i+k];
                mL_P1P1( nScalar.num[j], ii+k) += locL_P1P1[j][3*i+k];
                mL_P1P1_stab( nScalar.num[j], ii+k) += locL_P1P1_stab[j][3*i+k];
            }
        }
    }
    for(int i= 0; i < 4; ++i) {
        const IdxT ii= nScalar.num[i];
        if (ii==NoIdx || !(nScalar.WithUnknowns( i))) continue;
        for(int j=0; j <4; ++j) {
            const IdxT jj= nScalar.num[j];
            if (jj==NoIdx || !(nScalar.WithUnknowns( j))) continue;
            mM_ScalarP1( ii, jj) += locM_P1[i][j];
            mA_ScalarP1_stab( ii, jj) += locA_P1_stab[i][j];
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

    LocalStokesCL localStokes_;

    SMatrixCL<3,3> T;
    double det, absdet;
    LocalP2CL<> ls_loc;
    LocalP1CL<Point3DCL> v_loc;
    LocalNumbP1CL n, nScalar;

    ///\brief Computes the mapping from local to global data "n", the local matrices in loc and.
    void local_setup (const TetraCL& tet);
    ///\brief Update the global system.
    void update_global_system ();

  public:
    NavierStokesIFAccumulator_P1P1CL ( const VecDescCL& ls, const VecDescCL& v, const LsetBndDataCL& ls_bnd, const BndDataCL<Point3DCL>& v_bnd, IdxDescCL& P1FE, IdxDescCL& ScalarP1FE, MatrixCL& A_P1, MatrixCL& A_P1_stab, MatrixCL& B_P1P1,MatrixCL& Omega_P1P1, MatrixCL& N_P1, MatrixCL& NT_P1, MatrixCL& M_P1, MatrixCL& D_P1, MatrixCL& S_P1, MatrixCL& L_P1P1, MatrixCL& L_P1P1_stab, MatrixCL& M_ScalarP1, MatrixCL& A_ScalarP1_stab, MatrixCL& Schur_normalP1_stab, bool fullGradient);

    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    void visit (const TetraCL& sit);

    TetraAccumulatorCL* clone (int /*tid*/) { return new NavierStokesIFAccumulator_P1P1CL ( *this); }

};

NavierStokesIFAccumulator_P1P1CL::NavierStokesIFAccumulator_P1P1CL( const VecDescCL& ls, const VecDescCL& v, const LsetBndDataCL& ls_bnd, const BndDataCL<Point3DCL>& v_bnd, IdxDescCL& P1FE, IdxDescCL& ScalarP1FE, MatrixCL& A_P1, MatrixCL& A_P1_stab, MatrixCL& B_P1P1,MatrixCL& Omega_P1P1, MatrixCL& N_P1, MatrixCL& NT_P1, MatrixCL& M_P1,MatrixCL& D_P1, MatrixCL& S_P1, MatrixCL& L_P1P1, MatrixCL& L_P1P1_stab, MatrixCL& M_ScalarP1, MatrixCL& A_ScalarP1_stab, MatrixCL& Schur_normalP1_stab, bool fullGradient)
 : lset(ls), velocity(v), lset_bnd(ls_bnd), velocity_bnd(v_bnd), P1Idx_(P1FE), ScalarP1Idx_(ScalarP1FE), A_P1_(A_P1), A_P1_stab_(A_P1_stab), B_P1P1_(B_P1P1),Omega_P1P1_(Omega_P1P1),N_P1_(N_P1), NT_P1_(NT_P1),M_P1_(M_P1),D_P1_(D_P1), S_P1_(S_P1), L_P1P1_(L_P1P1), L_P1P1_stab_(L_P1P1_stab), M_ScalarP1_(M_ScalarP1), A_ScalarP1_stab_(A_ScalarP1_stab),Schur_normalP1_stab_(Schur_normalP1_stab), localStokes_( fullGradient)
{}

void NavierStokesIFAccumulator_P1P1CL::begin_accumulation ()
{
    std::cout << "entering NavierStokesIF: \n";
    const size_t num_unks_p1= P1Idx_.NumUnknowns(), num_unks_scalarp1= ScalarP1Idx_.NumUnknowns();
    mA_P1_= new MatrixBuilderCL( &A_P1_, num_unks_p1, num_unks_p1);
    mN_P1_= new MatrixBuilderCL( &N_P1_, num_unks_p1, num_unks_p1);
    mNT_P1_= new MatrixBuilderCL( &NT_P1_, num_unks_p1, num_unks_p1);
    mA_P1_stab_= new MatrixBuilderCL( &A_P1_stab_, num_unks_p1, num_unks_p1);
    mB_P1P1_= new MatrixBuilderCL( &B_P1P1_, num_unks_scalarp1, num_unks_p1);
    mOmega_P1P1_= new MatrixBuilderCL( &Omega_P1P1_, num_unks_scalarp1, num_unks_p1);
    mM_P1_= new MatrixBuilderCL( &M_P1_, num_unks_p1, num_unks_p1);
    mD_P1_= new MatrixBuilderCL( &D_P1_, num_unks_p1, num_unks_p1);
    mS_P1_= new MatrixBuilderCL( &S_P1_, num_unks_p1, num_unks_p1);
    mL_P1P1_= new MatrixBuilderCL( &L_P1P1_, num_unks_scalarp1, num_unks_p1);
    mL_P1P1_stab_= new MatrixBuilderCL( &L_P1P1_stab_, num_unks_scalarp1, num_unks_p1);
    mM_ScalarP1_= new MatrixBuilderCL( &M_ScalarP1_, num_unks_scalarp1, num_unks_scalarp1);
    mA_ScalarP1_stab_= new MatrixBuilderCL( &A_ScalarP1_stab_, num_unks_scalarp1, num_unks_scalarp1);
    mSchur_normalP1_stab_= new MatrixBuilderCL( &Schur_normalP1_stab_, num_unks_scalarp1, num_unks_scalarp1);
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
    ls_loc.assign( tet, lset, lset_bnd);
    v_loc.assign( tet, velocity, velocity_bnd);

    if (!equal_signs( ls_loc))
    {
        local_setup( tet);
        update_global_system();
    }

}

void NavierStokesIFAccumulator_P1P1CL::local_setup (const TetraCL& tet)
{
    GetTrafoTr( T, det, tet);
    absdet= std::fabs( det);

    n.assign( tet, P1Idx_, P1Idx_.GetBndInfo());
    nScalar.assign( tet, ScalarP1Idx_, ScalarP1Idx_.GetBndInfo());

    //GetLocalNumbP1NoBnd( numP1, tet, P1Idx_);
    //GetLocalNumbP1NoBnd( numScalarP1, tet, ScalarP1Idx_);


    localStokes_.calcIntegrands( T, ls_loc, v_loc, tet);
    localStokes_.calc3DIntegrands(T, ls_loc, tet);
    localStokes_.setupA_P1( locA_P1);
    localStokes_.setupN_P1( locN_P1);
    localStokes_.setupNT_P1( locNT_P1);
    localStokes_.setupA_P1_stab( locA_P1_stab, absdet);
    localStokes_.setupB_P1P1( locB_P1P1);
    localStokes_.setupOmega_P1P1( locOmega_P1P1);
    localStokes_.setupM_P1( locM_P1);

    localStokes_.setupD_P1( locD_P1);
    localStokes_.setupM_P1_stab( locM_P1_stab, absdet);
    localStokes_.setupS_P1( locS_P1);
    localStokes_.setupL_P1P1( locL_P1P1);
    localStokes_.setupL_P1P1_stab( locL_P1P1_stab, absdet);

}

void NavierStokesIFAccumulator_P1P1CL::update_global_system ()
{
    MatrixBuilderCL& mA_P1= *mA_P1_;
    MatrixBuilderCL& mA_P1_stab= *mA_P1_stab_;
    MatrixBuilderCL& mB_P1P1= *mB_P1P1_;
    MatrixBuilderCL& mOmega_P1P1= *mOmega_P1P1_;
    MatrixBuilderCL& mN_P1= *mN_P1_;
    MatrixBuilderCL& mNT_P1= *mNT_P1_;
    MatrixBuilderCL& mM_P1= *mM_P1_;
    MatrixBuilderCL& mD_P1= *mD_P1_;
    MatrixBuilderCL& mS_P1= *mS_P1_;
    MatrixBuilderCL& mL_P1P1= *mL_P1P1_;
    MatrixBuilderCL& mL_P1P1_stab= *mL_P1P1_stab_;
    MatrixBuilderCL& mM_ScalarP1= *mM_ScalarP1_;
    MatrixBuilderCL& mA_ScalarP1_stab= *mA_ScalarP1_stab_;
    MatrixBuilderCL& mSchur_normalP1_stab= *mSchur_normalP1_stab_;

int count =0;
    for(int i= 0; i < 4; ++i) {
        const IdxT ii= n.num[i];
        if (ii==NoIdx || !(n.WithUnknowns( i))) continue;
        for(int j=0; j < 4; ++j) {
            const IdxT jj= n.num[j];
            if (jj==NoIdx || !(n.WithUnknowns( j))) continue;
            for (int k=0; k<3; ++k) {
                for (int l=0; l<3; ++l) {
                    mA_P1( ii+k, jj+l) += locA_P1[3*j+l][3*i+k];
                    mN_P1( ii+k, jj+l) += locN_P1[3*j+l][3*i+k];
                    mNT_P1( ii+k, jj+l) += locNT_P1[3*j+l][3*i+k];
                    mS_P1( ii+k, jj+l) += locS_P1[3*j+l][3*i+k];

                    if(k == l) {
                    	mM_P1( ii+k, jj+l) += locM_P1[j][i];

                    	mD_P1( ii+k, jj+l) += locD_P1[j][i];

                    	count++;
                    	mA_P1_stab( ii+k, jj+l) += locA_P1_stab[j][i];
                    } else {
                    	mM_P1( ii+k, jj+l) += 0.;
                    	mD_P1( ii+k, jj+l) += 0.;
                    	mA_P1_stab( ii+k, jj+l) += 0.;
                    }

                }
             }
        }
        for (int j=0; j < 4; ++j) {
            if (nScalar.num[j]==NoIdx || !(nScalar.WithUnknowns( j))) continue;
            for (int k=0; k<3; ++k) {
                mB_P1P1( nScalar.num[j], ii+k) += locB_P1P1[j][3*i+k];
                mOmega_P1P1( nScalar.num[j], ii+k) += locOmega_P1P1[j][3*i+k];
                mL_P1P1( nScalar.num[j], ii+k) += locL_P1P1[j][3*i+k];
                mL_P1P1_stab( nScalar.num[j], ii+k) += locL_P1P1_stab[j][3*i+k];
            }
        }
    }
    for(int i= 0; i < 4; ++i) {
        const IdxT ii= nScalar.num[i];
        if (ii==NoIdx || !(nScalar.WithUnknowns( i))) continue;
        for(int j=0; j <4; ++j) {
            const IdxT jj= nScalar.num[j];
            if (jj==NoIdx || !(nScalar.WithUnknowns( j))) continue;
            mM_ScalarP1( ii, jj) += locM_P1[i][j];
            mA_ScalarP1_stab( ii, jj) += locM_P1_stab[i][j];
            mSchur_normalP1_stab( ii, jj) += locA_P1_stab[i][j];
        }
    }
}

void SetupStokesIF_P1P1( const MultiGridCL& MG_, MatDescCL* A_P1, MatDescCL* A_P1_stab, MatDescCL* B_P1P1, MatDescCL* M_P1, MatDescCL* S_P1, MatDescCL* L_P1P1, MatDescCL* L_P1P1_stab, MatDescCL* M_ScalarP1, MatDescCL* A_ScalarP1_stab, const VecDescCL& lset, const LsetBndDataCL& lset_bnd, bool fullgrad)
{
  ScopeTimerCL scope("SetupStokesIF_P1P1");
  StokesIFAccumulator_P1P1CL accu( lset, lset_bnd, *(A_P1->RowIdx), *(L_P1P1->RowIdx), A_P1->Data, A_P1_stab->Data, B_P1P1->Data, M_P1->Data, S_P1->Data, L_P1P1->Data, L_P1P1_stab->Data, M_ScalarP1->Data, A_ScalarP1_stab->Data, fullgrad);
  TetraAccumulatorTupleCL accus;
  //    MaybeAddProgressBar(MG_, "LapBeltr(P2) Setup", accus, RowIdx.TriangLevel());
  accus.push_back( &accu);
  accumulate( accus, MG_, A_P1->GetRowLevel(), A_P1->RowIdx->GetBndInfo());
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
        MatDescCL* Schur_normalP1_stab, const VecDescCL& lset, const LsetBndDataCL& lset_bnd, const VecDescCL& velocity, const BndDataCL<Point3DCL>& velocity_bnd, bool fullgrad)
{
  ScopeTimerCL scope("SetupNavierStokesIF_P1P1");
  NavierStokesIFAccumulator_P1P1CL accu( lset, velocity, lset_bnd, velocity_bnd, *(A_P1->RowIdx), *(L_P1P1->RowIdx), A_P1->Data, A_P1_stab->Data, B_P1P1->Data, Omega_P1P1->Data,  N_P1->Data, NT_P1->Data, M_P1->Data, D_P1->Data, S_P1->Data, L_P1P1->Data, L_P1P1_stab->Data, M_ScalarP1->Data, A_ScalarP1_stab->Data, Schur_normalP1_stab->Data, fullgrad);
  TetraAccumulatorTupleCL accus;
  //    MaybeAddProgressBar(MG_, "LapBeltr(P2) Setup", accus, RowIdx.TriangLevel());
  accus.push_back( &accu);
  accumulate( accus, MG_, A_P1->GetRowLevel(), /*A_P1->RowIdx->GetMatchingFunction(),*/ A_P1->RowIdx->GetBndInfo());
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

    LocalStokesCL localStokes_;

    SMatrixCL<3,3> T;
    double det, absdet;
    LocalP2CL<> ls_loc;
    LocalNumbP1CL n;
    LocalNumbP2CL nScalar;

    ///\brief Computes the mapping from local to global data "n", the local matrices in loc and.
    void local_setup (const TetraCL& tet);
    ///\brief Update the global system.
    void update_global_system ();

  public:
    StokesIFAccumulator_P1P2CL ( const VecDescCL& ls, const LsetBndDataCL& ls_bnd, IdxDescCL& P1FE, IdxDescCL& ScalarP2FE, MatrixCL& A_P1, MatrixCL& A_P1_stab, MatrixCL& B_P2P1, MatrixCL& M_P1, MatrixCL& S_P1, MatrixCL& L_P2P1, MatrixCL& L_P2P1_stab, MatrixCL& M_ScalarP2, MatrixCL& A_ScalarP2_stab, bool fullGradient);

    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    void visit (const TetraCL& sit);

    TetraAccumulatorCL* clone (int /*tid*/) { return new StokesIFAccumulator_P1P2CL ( *this); }

};

StokesIFAccumulator_P1P2CL::StokesIFAccumulator_P1P2CL( const VecDescCL& ls, const LsetBndDataCL& ls_bnd, IdxDescCL& P1FE, IdxDescCL& ScalarP2FE, MatrixCL& A_P1, MatrixCL& A_P1_stab, MatrixCL& B_P2P1, MatrixCL& M_P1, MatrixCL& S_P1, MatrixCL& L_P2P1, MatrixCL& L_P2P1_stab, MatrixCL& M_ScalarP2, MatrixCL& A_ScalarP2_stab, bool fullGradient)
 : lset(ls), lset_bnd(ls_bnd), P1Idx_(P1FE), ScalarP2Idx_(ScalarP2FE), A_P1_(A_P1), A_P1_stab_(A_P1_stab), B_P2P1_(B_P2P1), M_P1_(M_P1), S_P1_(S_P1), L_P2P1_(L_P2P1), L_P2P1_stab_(L_P2P1_stab), M_ScalarP2_(M_ScalarP2), A_ScalarP2_stab_(A_ScalarP2_stab), localStokes_( fullGradient)
{}

void StokesIFAccumulator_P1P2CL::begin_accumulation ()
{
    std::cout << "entering StokesIF: \n";
    const size_t num_unks_p1= P1Idx_.NumUnknowns(), num_unks_scalarp2= ScalarP2Idx_.NumUnknowns();
    mA_P1_= new MatrixBuilderCL( &A_P1_, num_unks_p1, num_unks_p1);
    mA_P1_stab_= new MatrixBuilderCL( &A_P1_stab_, num_unks_p1, num_unks_p1);
    mB_P2P1_= new MatrixBuilderCL( &B_P2P1_, num_unks_scalarp2, num_unks_p1);
    mM_P1_= new MatrixBuilderCL( &M_P1_, num_unks_p1, num_unks_p1);
    mS_P1_= new MatrixBuilderCL( &S_P1_, num_unks_p1, num_unks_p1);
    mL_P2P1_= new MatrixBuilderCL( &L_P2P1_, num_unks_scalarp2, num_unks_p1);
    mL_P2P1_stab_= new MatrixBuilderCL( &L_P2P1_stab_, num_unks_scalarp2, num_unks_p1);
    mM_ScalarP2_= new MatrixBuilderCL( &M_ScalarP2_, num_unks_scalarp2, num_unks_scalarp2);
    mA_ScalarP2_stab_= new MatrixBuilderCL( &A_ScalarP2_stab_, num_unks_scalarp2, num_unks_scalarp2);
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
    ls_loc.assign( tet, lset, lset_bnd);

    if (!equal_signs( ls_loc))
    {
        local_setup( tet);
        update_global_system();
    }

}

void StokesIFAccumulator_P1P2CL::local_setup (const TetraCL& tet)
{
    GetTrafoTr( T, det, tet);
    absdet= std::fabs( det);

    n.assign( tet, P1Idx_, P1Idx_.GetBndInfo());
    nScalar.assign( tet, ScalarP2Idx_, ScalarP2Idx_.GetBndInfo());

    //GetLocalNumbP1NoBnd( numP1, tet, P1Idx_);
    //GetLocalNumbP1NoBnd( numScalarP1, tet, ScalarP1Idx_);
    localStokes_.calcIntegrands( T, ls_loc, tet);
    localStokes_.calc3DIntegrands(T, ls_loc, tet);
    localStokes_.setupA_P1( locA_P1);
    localStokes_.setupA_P1_stab( locA_P1_stab, absdet);
    localStokes_.setupA_P2_stab( locA_P2_stab, absdet);
    localStokes_.setupB_P2P1( locB_P2P1);
    localStokes_.setupM_P1( locM_P1);
    localStokes_.setupM_P2( locM_P2);
    localStokes_.setupS_P1( locS_P1);
    localStokes_.setupL_P2P1( locL_P2P1);
    localStokes_.setupL_P2P1_stab( locL_P2P1_stab, absdet);
}

void StokesIFAccumulator_P1P2CL::update_global_system ()
{
    MatrixBuilderCL& mA_P1= *mA_P1_;
    MatrixBuilderCL& mA_P1_stab= *mA_P1_stab_;
    MatrixBuilderCL& mB_P2P1= *mB_P2P1_;
    MatrixBuilderCL& mM_P1= *mM_P1_;
    MatrixBuilderCL& mS_P1= *mS_P1_;
    MatrixBuilderCL& mL_P2P1= *mL_P2P1_;
    MatrixBuilderCL& mL_P2P1_stab= *mL_P2P1_stab_;
    MatrixBuilderCL& mM_ScalarP2= *mM_ScalarP2_;
    MatrixBuilderCL& mA_ScalarP2_stab= *mA_ScalarP2_stab_;

    for(int i= 0; i < 4; ++i) {
        const IdxT ii= n.num[i];
        if (ii==NoIdx || !(n.WithUnknowns( i))) continue;
        for(int j=0; j < 4; ++j) {
            const IdxT jj= n.num[j];
            if (jj==NoIdx || !(n.WithUnknowns( j))) continue;
            for (int k=0; k<3; ++k) {
                for (int l=0; l<3; ++l) {
                    mA_P1( ii+k, jj+l) += locA_P1[3*j+l][3*i+k];
                    mS_P1( ii+k, jj+l) += locS_P1[3*j+l][3*i+k];
                    if(k == l) {
                        mM_P1( ii+k, jj+l) += locM_P1[j][i];
                        mA_P1_stab( ii+k, jj+l) += locA_P1_stab[j][i];
                    } else {
                        mM_P1( ii+k, jj+l) += 0.;
                        mA_P1_stab( ii+k, jj+l) += 0.;
                    }
                }
             }
        }
        for (int j=0; j < 10; ++j) {
            if (nScalar.num[j]==NoIdx || !(nScalar.WithUnknowns( j))) continue;
            for (int k=0; k<3; ++k) {
                mB_P2P1( nScalar.num[j], ii+k) += locB_P2P1[j][3*i+k];
                mL_P2P1( nScalar.num[j], ii+k) += locL_P2P1[j][3*i+k];
                mL_P2P1_stab( nScalar.num[j], ii+k) += locL_P2P1_stab[j][3*i+k];
            }
        }
    }
    for(int i= 0; i < 10; ++i) {
        const IdxT ii= nScalar.num[i];
        if (ii==NoIdx || !(nScalar.WithUnknowns( i))) continue;
        for(int j=0; j < 10; ++j) {
            const IdxT jj= nScalar.num[j];
            if (jj==NoIdx || !(nScalar.WithUnknowns( j))) continue;
            mM_ScalarP2( ii, jj) += locM_P2[i][j];
            mA_ScalarP2_stab( ii, jj) += locA_P2_stab[i][j];
        }
    }
}

void SetupStokesIF_P1P2( const MultiGridCL& MG_, MatDescCL* A_P1, MatDescCL* A_P1_stab, MatDescCL* B_P2P1, MatDescCL* M_P1, MatDescCL* S_P1, MatDescCL* L_P2P1, MatDescCL* L_P2P1_stab, MatDescCL* M_ScalarP2, MatDescCL* A_ScalarP2_stab, const VecDescCL& lset, const LsetBndDataCL& lset_bnd, bool fullgrad)
{
  ScopeTimerCL scope("SetupStokesIF_P1P2");
  StokesIFAccumulator_P1P2CL accu( lset, lset_bnd, *(A_P1->RowIdx), *(L_P2P1->RowIdx), A_P1->Data, A_P1_stab->Data, B_P2P1->Data, M_P1->Data, S_P1->Data, L_P2P1->Data, L_P2P1_stab->Data, M_ScalarP2->Data, A_ScalarP2_stab->Data, fullgrad);
  TetraAccumulatorTupleCL accus;
  //    MaybeAddProgressBar(MG_, "LapBeltr(P2) Setup", accus, RowIdx.TriangLevel());
  accus.push_back( &accu);
  accumulate( accus, MG_, A_P1->GetRowLevel(), A_P1->RowIdx->GetBndInfo());
}

/// \brief Accumulator to set up the matrices for interface Stokes.
class StokesIFAccumulator_P2P1CL : public TetraAccumulatorCL
{
  private:
    const VecDescCL& lset;
    const LsetBndDataCL& lset_bnd;

    IdxDescCL& P2Idx_;
    IdxDescCL& ScalarP1Idx_;
    IdxT numP2[10], numScalarP1[4];

    MatrixCL &A_P2_, &A_P2_stab_, &B_P1P2_, &M_P2_, &S_P2_, &L_P1P2_, &L_P1P2_stab_, &M_ScalarP1_, &A_ScalarP1_stab_;
    MatrixBuilderCL *mA_P2_, *mA_P2_stab_, *mB_P1P2_, *mM_P2_, *mS_P2_, *mL_P1P2_, *mL_P1P2_stab_, *mM_ScalarP1_, *mA_ScalarP1_stab_ ;

    double locA_P2[30][30], locA_P2_stab[10][10], locB_P1P2[4][30], locM_P2[10][10], locS_P2[30][30], locL_P1P2[4][30], locL_P1P2_stab[4][30], locM_ScalarP1[4][4], locA_ScalarP1_stab[4][4];

    LocalStokesCL localStokes_;

    SMatrixCL<3,3> T;
    double det, absdet;
    LocalP2CL<> ls_loc;

    ///\brief Computes the mapping from local to global data "n", the local matrices in loc and.
    void local_setup (const TetraCL& tet);
    ///\brief Update the global system.
    void update_global_system ();

  public:
    StokesIFAccumulator_P2P1CL ( const VecDescCL& ls, const LsetBndDataCL& ls_bnd, IdxDescCL& P2FE, IdxDescCL& ScalarP1FE, MatrixCL& A_P2, MatrixCL& A_P2_stab, MatrixCL& B_P1P2, MatrixCL& M_P2, MatrixCL& S_P2, MatrixCL& L_P1P2, MatrixCL& L_P1P2_stab, MatrixCL& M_ScalarP1, MatrixCL& A_ScalarP1_stab, bool fullGradient);

    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    void visit (const TetraCL& sit);

    TetraAccumulatorCL* clone (int /*tid*/) { return new StokesIFAccumulator_P2P1CL ( *this); }
};

StokesIFAccumulator_P2P1CL::StokesIFAccumulator_P2P1CL( const VecDescCL& ls, const LsetBndDataCL& ls_bnd, IdxDescCL& P2FE, IdxDescCL& ScalarP1FE, MatrixCL& A_P2, MatrixCL& A_P2_stab, MatrixCL& B_P1P2, MatrixCL& M_P2, MatrixCL& S_P2, MatrixCL& L_P1P2, MatrixCL& L_P1P2_stab, MatrixCL& M_ScalarP1, MatrixCL& A_ScalarP1_stab, bool fullGradient)
 : lset(ls), lset_bnd(ls_bnd), P2Idx_(P2FE), ScalarP1Idx_(ScalarP1FE), A_P2_(A_P2), A_P2_stab_(A_P2_stab), B_P1P2_(B_P1P2), M_P2_(M_P2), S_P2_(S_P2), L_P1P2_(L_P1P2), L_P1P2_stab_(L_P1P2_stab), M_ScalarP1_(M_ScalarP1), A_ScalarP1_stab_(A_ScalarP1_stab), localStokes_( fullGradient)
{}

void StokesIFAccumulator_P2P1CL::begin_accumulation ()
{
    std::cout << "entering StokesIF: \n";
    const size_t num_unks_scalarp1= ScalarP1Idx_.NumUnknowns(), num_unks_p2= P2Idx_.NumUnknowns();
    mA_P2_= new MatrixBuilderCL( &A_P2_, num_unks_p2, num_unks_p2);
    mA_P2_stab_= new MatrixBuilderCL( &A_P2_stab_, num_unks_p2, num_unks_p2);
    mB_P1P2_= new MatrixBuilderCL( &B_P1P2_, num_unks_scalarp1, num_unks_p2);
    mM_P2_= new MatrixBuilderCL( &M_P2_, num_unks_p2, num_unks_p2);
    mS_P2_= new MatrixBuilderCL( &S_P2_, num_unks_p2, num_unks_p2);
    mL_P1P2_= new MatrixBuilderCL( &L_P1P2_, num_unks_scalarp1, num_unks_p2);
    mL_P1P2_stab_= new MatrixBuilderCL( &L_P1P2_stab_, num_unks_scalarp1, num_unks_p2);
    mM_ScalarP1_= new MatrixBuilderCL( &M_ScalarP1_, num_unks_scalarp1, num_unks_scalarp1);
    mA_ScalarP1_stab_= new MatrixBuilderCL( &A_ScalarP1_stab_, num_unks_scalarp1, num_unks_scalarp1);
}

void StokesIFAccumulator_P2P1CL::finalize_accumulation ()
{
    mA_P2_->Build();
    delete mA_P2_;
    mA_P2_stab_->Build();
    delete mA_P2_stab_;
    mB_P1P2_->Build();
    delete mB_P1P2_;
    mM_P2_->Build();
    delete mM_P2_;
    mS_P2_->Build();
    delete mS_P2_;
    mL_P1P2_->Build();
    delete mL_P1P2_;
    mL_P1P2_stab_->Build();
    delete mL_P1P2_stab_;
    mM_ScalarP1_->Build();
    delete mM_ScalarP1_;
    mA_ScalarP1_stab_->Build();
    delete mA_ScalarP1_stab_;

#ifndef _PAR
    std::cout << "StokesIF_P2P1:\t" << A_P2_.num_nonzeros() << " nonzeros in A, " << A_P2_stab_.num_nonzeros() << " nonzeros in A_stab, "
              << B_P1P2_.num_nonzeros() << " nonzeros in B, " << M_P2_.num_nonzeros() << " nonzeros in M, " << S_P2_.num_nonzeros() << " nonzeros in S, "
              << L_P1P2_.num_nonzeros() << " nonzeros in L, " << L_P1P2_stab_.num_nonzeros() << " nonzeros in L_stab, "
              << M_ScalarP1_.num_nonzeros() << " nonzeros in M_ScalarP1, " << A_ScalarP1_stab_.num_nonzeros() << " nonzeros in A_ScalarP1_stab\n";
#endif
    // std::cout << '\n';
}

void StokesIFAccumulator_P2P1CL::visit (const TetraCL& tet)
{
    ls_loc.assign( tet, lset, lset_bnd);

    if (!equal_signs( ls_loc))
    {
        local_setup( tet);
        update_global_system();
    }

}

void StokesIFAccumulator_P2P1CL::local_setup (const TetraCL& tet)
{
    GetTrafoTr( T, det, tet);
    absdet= std::fabs( det);

    GetLocalNumbP2NoBnd( numP2, tet, P2Idx_);
    GetLocalNumbP1NoBnd( numScalarP1, tet, ScalarP1Idx_);
    localStokes_.calcIntegrands( T, ls_loc, tet);
    localStokes_.calc3DIntegrands(T, ls_loc, tet);
    localStokes_.setupA_P2( locA_P2);
    localStokes_.setupA_P2_stab( locA_P2_stab, absdet);
    localStokes_.setupB_P1P2( locB_P1P2);
    localStokes_.setupM_P2( locM_P2);
    localStokes_.setupS_P2( locS_P2);
    localStokes_.setupL_P1P2( locL_P1P2);
    localStokes_.setupL_P1P2_stab( locL_P1P2_stab, absdet);
    localStokes_.setupM_P1( locM_ScalarP1);
    localStokes_.setupA_P1_stab( locA_ScalarP1_stab, absdet);
}

void StokesIFAccumulator_P2P1CL::update_global_system ()
{
    MatrixBuilderCL& mA_P2= *mA_P2_;
    MatrixBuilderCL& mA_P2_stab= *mA_P2_stab_;
    MatrixBuilderCL& mB_P1P2= *mB_P1P2_;
    MatrixBuilderCL& mM_P2= *mM_P2_;
    MatrixBuilderCL& mS_P2= *mS_P2_;
    MatrixBuilderCL& mL_P1P2= *mL_P1P2_;
    MatrixBuilderCL& mL_P1P2_stab= *mL_P1P2_stab_;
    MatrixBuilderCL& mM_ScalarP1= *mM_ScalarP1_;
    MatrixBuilderCL& mA_ScalarP1_stab= *mA_ScalarP1_stab_;

    for(int i= 0; i < 10; ++i) {
        const IdxT ii= numP2[i];
        if (ii==NoIdx) continue;
        for(int j=0; j <10; ++j) {
            const IdxT jj= numP2[j];
            if (jj==NoIdx) continue;
            for (int k=0; k<3; ++k) {
                for (int l=0; l<3; ++l) {
                    mA_P2( ii+k, jj+l) += locA_P2[3*j+l][3*i+k];
                    mS_P2( ii+k, jj+l) += locS_P2[3*j+l][3*i+k];
                    if(k == l) {
                        mM_P2( ii+k, jj+l) += locM_P2[j][i];
                        mA_P2_stab( ii+k, jj+l) += locA_P2_stab[j][i];
                    } else {
                        mM_P2( ii+k, jj+l) += 0.;
                        mA_P2_stab( ii+k, jj+l) += 0.;
                    }
                }
             }
        }
        for (int j=0; j < 4; ++j) {
            if (numScalarP1[j]==NoIdx) continue;
            for (int k=0; k<3; ++k) {
                mB_P1P2( numScalarP1[j], ii+k) += locB_P1P2[j][3*i+k];
                mL_P1P2( numScalarP1[j], ii+k) += locL_P1P2[j][3*i+k];
                mL_P1P2_stab( numScalarP1[j], ii+k) += locL_P1P2_stab[j][3*i+k];
            }
        }
    }
    for(int i= 0; i < 4; ++i) {
        const IdxT ii= numScalarP1[i];
        if (ii==NoIdx) continue;
        for(int j=0; j <4; ++j) {
            const IdxT jj= numScalarP1[j];
            if (jj==NoIdx) continue;
            mM_ScalarP1( ii, jj) += locM_ScalarP1[i][j];
            mA_ScalarP1_stab( ii, jj) += locA_ScalarP1_stab[i][j];
        }
    }
}

void SetupStokesIF_P2P1(
        const MultiGridCL& MG_,
        MatDescCL* A_P2,
        MatDescCL* A_P2_stab,
        MatDescCL* B_P1P2,
        MatDescCL* M_P2,
        MatDescCL* S_P2,
        MatDescCL* L_P1P2,
        MatDescCL* L_P1P2_stab,
        MatDescCL* M_ScalarP1,
        MatDescCL* A_ScalarP1_stab,
        const VecDescCL& lset,
        const LsetBndDataCL& lset_bnd,
        bool fullgrad
) {
  ScopeTimerCL scope("SetupStokesIF");
  StokesIFAccumulator_P2P1CL accu( lset, lset_bnd, *(A_P2->RowIdx), *(B_P1P2->RowIdx), A_P2->Data, A_P2_stab->Data, B_P1P2->Data, M_P2->Data, S_P2->Data, L_P1P2->Data, L_P1P2_stab->Data, M_ScalarP1->Data, A_ScalarP1_stab->Data, fullgrad);
  TetraAccumulatorTupleCL accus;
  //    MaybeAddProgressBar(MG_, "LapBeltr(P2) Setup", accus, RowIdx.TriangLevel());
  accus.push_back( &accu);
  accumulate( accus, MG_, A_P2->GetRowLevel(), A_P2->RowIdx->GetBndInfo());
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

    LocalStokesCL localStokes_;

    SMatrixCL<3,3> T;
    double det, absdet;
    LocalP2CL<> ls_loc;

    ///\brief Computes the mapping from local to global data "n", the local matrices in loc and.
    void local_setup (const TetraCL& tet);
    ///\brief Update the global system.
    void update_global_system ();

  public:
    StokesIFAccumulator_P2P2CL ( const VecDescCL& ls, const LsetBndDataCL& ls_bnd, IdxDescCL& P2FE, IdxDescCL& ScalarP2FE, MatrixCL& A_P2, MatrixCL& A_P2_stab, MatrixCL& B_P2P2, MatrixCL& M_P2, MatrixCL& S_P2, MatrixCL& L_P2P2, MatrixCL& L_P2P2_stab, MatrixCL& M_ScalarP2, MatrixCL& A_ScalarP2_stab, bool fullGradient);

    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    void visit (const TetraCL& sit);

    TetraAccumulatorCL* clone (int /*tid*/) { return new StokesIFAccumulator_P2P2CL ( *this); }
};

StokesIFAccumulator_P2P2CL::StokesIFAccumulator_P2P2CL( const VecDescCL& ls, const LsetBndDataCL& ls_bnd, IdxDescCL& P2FE, IdxDescCL& ScalarP2FE, MatrixCL& A_P2, MatrixCL& A_P2_stab, MatrixCL& B_P2P2, MatrixCL& M_P2, MatrixCL& S_P2, MatrixCL& L_P2P2, MatrixCL& L_P2P2_stab, MatrixCL& M_ScalarP2, MatrixCL& A_ScalarP2_stab, bool fullGradient)
 : lset(ls), lset_bnd(ls_bnd), P2Idx_(P2FE), ScalarP2Idx_(ScalarP2FE), A_P2_(A_P2), A_P2_stab_(A_P2_stab), B_P2P2_(B_P2P2), M_P2_(M_P2), S_P2_(S_P2), L_P2P2_(L_P2P2), L_P2P2_stab_(L_P2P2_stab), M_ScalarP2_(M_ScalarP2), A_ScalarP2_stab_(A_ScalarP2_stab), localStokes_( fullGradient)
{}

void StokesIFAccumulator_P2P2CL::begin_accumulation ()
{
    std::cout << "entering StokesIF: \n";
    const size_t num_unks_scalarp2= ScalarP2Idx_.NumUnknowns(), num_unks_p2= P2Idx_.NumUnknowns();
    mA_P2_= new MatrixBuilderCL( &A_P2_, num_unks_p2, num_unks_p2);
    mA_P2_stab_= new MatrixBuilderCL( &A_P2_stab_, num_unks_p2, num_unks_p2);
    mB_P2P2_= new MatrixBuilderCL( &B_P2P2_, num_unks_scalarp2, num_unks_p2);
    mM_P2_= new MatrixBuilderCL( &M_P2_, num_unks_p2, num_unks_p2);
    mS_P2_= new MatrixBuilderCL( &S_P2_, num_unks_p2, num_unks_p2);
    mL_P2P2_= new MatrixBuilderCL( &L_P2P2_, num_unks_scalarp2, num_unks_p2);
    mL_P2P2_stab_= new MatrixBuilderCL( &L_P2P2_stab_, num_unks_scalarp2, num_unks_p2);
    mM_ScalarP2_= new MatrixBuilderCL( &M_ScalarP2_, num_unks_scalarp2, num_unks_scalarp2);
    mA_ScalarP2_stab_= new MatrixBuilderCL( &A_ScalarP2_stab_, num_unks_scalarp2, num_unks_scalarp2);
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
    ls_loc.assign( tet, lset, lset_bnd);

    if (!equal_signs( ls_loc))
    {
        local_setup( tet);
        update_global_system();
    }

}

void StokesIFAccumulator_P2P2CL::local_setup (const TetraCL& tet)
{
    GetTrafoTr( T, det, tet);
    absdet= std::fabs( det);

    GetLocalNumbP2NoBnd( numP2, tet, P2Idx_);
    GetLocalNumbP2NoBnd( numScalarP2, tet, ScalarP2Idx_);
    localStokes_.calcIntegrands( T, ls_loc, tet);
    localStokes_.calc3DIntegrands(T, ls_loc, tet);
    localStokes_.setupA_P2( locA_P2);
    localStokes_.setupA_P2_stab( locA_P2_stab, absdet);
    localStokes_.setupB_P2P2( locB_P2P2);
    localStokes_.setupM_P2( locM_P2);
    localStokes_.setupS_P2( locS_P2);
    localStokes_.setupL_P2P2( locL_P2P2);
    localStokes_.setupL_P2P2_stab( locL_P2P2_stab, absdet);
    localStokes_.setupM_P2( locM_ScalarP2);
    localStokes_.setupA_P2_stab( locA_ScalarP2_stab, absdet);
}

void StokesIFAccumulator_P2P2CL::update_global_system ()
{
    MatrixBuilderCL& mA_P2= *mA_P2_;
    MatrixBuilderCL& mA_P2_stab= *mA_P2_stab_;
    MatrixBuilderCL& mB_P2P2= *mB_P2P2_;
    MatrixBuilderCL& mM_P2= *mM_P2_;
    MatrixBuilderCL& mS_P2= *mS_P2_;
    MatrixBuilderCL& mL_P2P2= *mL_P2P2_;
    MatrixBuilderCL& mL_P2P2_stab= *mL_P2P2_stab_;
    MatrixBuilderCL& mM_ScalarP2= *mM_ScalarP2_;
    MatrixBuilderCL& mA_ScalarP2_stab= *mA_ScalarP2_stab_;

    for(int i= 0; i < 10; ++i) {
        const IdxT ii= numP2[i];
        if (ii==NoIdx) continue;
        for(int j=0; j <10; ++j) {
            const IdxT jj= numP2[j];
            if (jj==NoIdx) continue;
            for (int k=0; k<3; ++k) {
                for (int l=0; l<3; ++l) {
                    mA_P2( ii+k, jj+l) += locA_P2[3*j+l][3*i+k];
                    mS_P2( ii+k, jj+l) += locS_P2[3*j+l][3*i+k];
                    if(k == l) {
                        mM_P2( ii+k, jj+l) += locM_P2[j][i];
                        mA_P2_stab( ii+k, jj+l) += locA_P2_stab[j][i];
                    } else {
                        mM_P2( ii+k, jj+l) += 0.;
                        mA_P2_stab( ii+k, jj+l) += 0.;
                    }
                }
             }
        }
        for (int j=0; j < 10; ++j) {
            if (numScalarP2[j]==NoIdx) continue;
            for (int k=0; k<3; ++k) {
                mB_P2P2( numScalarP2[j], ii+k) += locB_P2P2[j][3*i+k];
                mL_P2P2( numScalarP2[j], ii+k) += locL_P2P2[j][3*i+k];
                mL_P2P2_stab( numScalarP2[j], ii+k) += locL_P2P2_stab[j][3*i+k];
            }
        }
    }
    for(int i= 0; i < 10; ++i) {
        const IdxT ii= numScalarP2[i];
        if (ii==NoIdx) continue;
        for(int j=0; j < 10; ++j) {
            const IdxT jj= numScalarP2[j];
            if (jj==NoIdx) continue;
            mM_ScalarP2( ii, jj) += locM_ScalarP2[i][j];
            mA_ScalarP2_stab( ii, jj) += locA_ScalarP2_stab[i][j];
        }
    }
}

void SetupStokesIF_P2P2( const MultiGridCL& MG_, MatDescCL* A_P2, MatDescCL* A_P2_stab, MatDescCL* B_P2P2, MatDescCL* M_P2, MatDescCL* S_P2, MatDescCL* L_P2P2, MatDescCL* L_P2P2_stab, MatDescCL* M_ScalarP2, MatDescCL* A_ScalarP2_stab, const VecDescCL& lset, const LsetBndDataCL& lset_bnd, bool fullgrad)
{
  ScopeTimerCL scope("SetupStokesIF");
  StokesIFAccumulator_P2P2CL accu( lset, lset_bnd, *(A_P2->RowIdx), *(B_P2P2->RowIdx), A_P2->Data, A_P2_stab->Data, B_P2P2->Data, M_P2->Data, S_P2->Data, L_P2P2->Data, L_P2P2_stab->Data, M_ScalarP2->Data, A_ScalarP2_stab->Data, fullgrad);
  TetraAccumulatorTupleCL accus;
  //    MaybeAddProgressBar(MG_, "LapBeltr(P2) Setup", accus, RowIdx.TriangLevel());
  accus.push_back( &accu);
  accumulate( accus, MG_, A_P2->GetRowLevel(), A_P2->RowIdx->GetBndInfo());
}

void SetupInterfaceVectorRhsP1 (const MultiGridCL& mg, VecDescCL* v,
    const VecDescCL& ls, const BndDataCL<>& lsetbnd, instat_vector_fun_ptr f, double t)
{
    const IdxT num_unks= v->RowIdx->NumUnknowns();
    const Uint lvl = v->GetLevel();

    double totalsurfarea(0.);

    v->Clear(0);

    //IdxT num[10];

    LocalNumbP1CL n;

    std::cout << "entering SetupInterfaceVectorRhsP1: " << num_unks << " dof... ";

    LocalP1CL<> p1[4];
    P1DiscCL::GetP1Basis( p1);

//    LocalP2CL<> p2[10];
//    for(int i=0; i<10; i++)  // P2-Basis-Functions
//      p2[i][i] = 1.;

    Quad5_2DCL<double> q1[4]; //basis function

    // double det, absdet;

    InterfaceTriangleCL triangle;

    Point3DCL e_k;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        triangle.Init( *it, ls, lsetbnd);
        if (triangle.Intersects()) { // We are at the phase boundary.
          n.assign( *it, *v->RowIdx, lsetbnd);
          //GetLocalNumbP2NoBnd( num, *it, *v->RowIdx);

            for (int ch= 0; ch < 8; ++ch) {
                triangle.ComputeForChild( ch);
                for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri) {
                  //P2DiscCL::GetP2Basis( q2, &triangle.GetBary( tri));
                    for (int i= 0; i < 4; ++i)
                      q1[i].assign( p1[i], &triangle.GetBary( tri));

                    Quad5_2DCL<Point3DCL> qf( *it, &triangle.GetBary( tri), f,t);
                    Quad5_2DCL<double> surfarea(1.0), rhs;
                    totalsurfarea += surfarea.quad( triangle.GetAbsDet( tri));

                    for (int i= 0; i < 4; ++i) {
                        if (n.num[i] == NoIdx) continue;
                        //rhs = qf*q2[i];
                        for (int k=0; k<3; ++k) {
                            e_k= DROPS::std_basis<3>( k+1);
                            rhs = dot(e_k, qf)*q1[i];
                            v->Data[n.num[i] + k] += rhs.quad( triangle.GetAbsDet( tri));
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
    const VecDescCL& ls, const BndDataCL<>& lsetbnd, instat_vector_fun_ptr f)
{
    const IdxT num_unks= v->RowIdx->NumUnknowns();
    const Uint lvl = v->GetLevel();

    double totalsurfarea(0.);

    v->Clear(0);

    //IdxT num[10];

    LocalNumbP2CL n;

    std::cout << "entering SetupInterfaceVectorRhsP2: " << num_unks << " dof... ";

    LocalP2CL<> p2[10];
    P2DiscCL::GetP2Basis( p2);

//    LocalP2CL<> p2[10];
//    for(int i=0; i<10; i++)  // P2-Basis-Functions
//      p2[i][i] = 1.;

    Quad5_2DCL<double> q2[10]; //basis function

    // double det, absdet;

    InterfaceTriangleCL triangle;

    Point3DCL e_k;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        triangle.Init( *it, ls, lsetbnd);
        if (triangle.Intersects()) { // We are at the phase boundary.
          n.assign( *it, *v->RowIdx, lsetbnd);
          //GetLocalNumbP2NoBnd( num, *it, *v->RowIdx);

            for (int ch= 0; ch < 8; ++ch) {
                triangle.ComputeForChild( ch);
                for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri) {
                  //P2DiscCL::GetP2Basis( q2, &triangle.GetBary( tri));
                    for (int i= 0; i < 10; ++i)
                      q2[i].assign( p2[i], &triangle.GetBary( tri));

                    Quad5_2DCL<Point3DCL> qf( *it, &triangle.GetBary( tri), f);
                    Quad5_2DCL<double> surfarea(1.0), rhs;
                    totalsurfarea += surfarea.quad( triangle.GetAbsDet( tri));

                    for (int i= 0; i < 10; ++i) {
                        if (n.num[i] == NoIdx) continue;
                        //rhs = qf*q2[i];
                        for (int k=0; k<3; ++k) {
                            e_k= DROPS::std_basis<3>( k+1);
                            rhs = dot(e_k, qf)*q2[i];
                            v->Data[n.num[i] + k] += rhs.quad( triangle.GetAbsDet( tri));
                        }
                    }
                }
            }
        }
    }
    std::cout << " VectorRhs set up." << std::endl;
    std::cout << "Total surface area: "<<totalsurfarea<<std::endl;
}


void P1Init (instat_scalar_fun_ptr icf, VecDescCL& ic, const MultiGridCL& mg, double t)
{
    const Uint lvl= ic.GetLevel(),
               idx= ic.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
        if (it->Unknowns.Exist( idx))
            ic.Data[it->Unknowns( idx)]= icf( it->GetCoord(), t);
    }
    ic.t= t;
}

InterfaceDebugP2CL::InterfaceDebugP2CL (const InterfaceCommonDataP2CL& cdata)
    : cdata_( cdata), to_iface( 0), ref_dp( 0), ref_abs_det( 0), true_area( -1.)
{}

void InterfaceDebugP2CL::begin_accumulation ()
{
    max_dph_err= 0;
    surfacemeasP1= 0.;
    surfacemeasP2= 0.;
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
            cdata.quaqua.set_point( &t, std_basis<4>( i + 1))
                        .base_point();
            const TetraBaryPairT& tb= cdata.quaqua.get_base_point();
            Point3DCL offset= t.GetVertex( i)->GetCoord() - GetWorldCoord( *tb.first, tb.second);
            const size_t dof= t.GetVertex( i)->Unknowns( sys);
            std::copy( Addr( offset), Addr( offset) + 3, &to_iface->Data[dof]);
        }
    }

    for (SurfacePatchCL::const_vertex_iterator it= cdata.surf.vertex_begin(); it != cdata.surf.vertex_end(); ++it) {
        cdata.quaqua.set_point( &t, *it)
             .jacobian();
//         const TetraBaryPairT& tb= cdata.quaqua.get_base_point();
        const Point3DCL& x= GetWorldCoord( t, *it);
//         const Point3DCL& xb= GetWorldCoord( *tb.first, tb.second);
//         std::cout  << "    |x-xb|: " << (x - xb).norm();

        const SMatrixCL<3,3>& dph= cdata.quaqua.get_jacobian();
        SMatrixCL<3,3> diff_dp= ref_dp != 0 ? dph - ref_dp( x, 0.) : SMatrixCL<3,3>();
        const double dph_err= std::sqrt( frobenius_norm_sq( diff_dp));
        max_dph_err= std::max( max_dph_err, dph_err);
//         std::cout  << " |dph -dp|_F: " << dph_err;

        const double absdet= abs_det( t, *it, cdata.quaqua, cdata.surf),
                     absdet_err= ref_abs_det != 0 ? std::abs( absdet - ref_abs_det( t, *it, cdata.surf)) : 0.;
        max_absdet_err= std::max( max_absdet_err, absdet_err);
//         std::cout  << " |\\mu - \\mu^s|: " << absdet_err << std::endl;


        Point3DCL n=x/x.norm();
        SMatrixCL<3,3> diff_dp2= diff_dp - outer_product( n, transp_mul( diff_dp, n));
        Point3DCL nh;
        Point3DCL v1= GetWorldCoord( t, cdata.surf.vertex_begin()[1]) - GetWorldCoord( t, cdata.surf.vertex_begin()[0]),
                  v2= GetWorldCoord( t, cdata.surf.vertex_begin()[2]) - GetWorldCoord( t, cdata.surf.vertex_begin()[0]);
        cross_product( nh, v1, v2);
        nh/= nh.norm();
        diff_dp2= diff_dp2 - outer_product( diff_dp2*nh, nh);
        const double dph2_err= std::sqrt( frobenius_norm_sq( diff_dp2));
//         std::cout  << " |P(dph -dp)\\hat P|_F: " << dph2_err << std::endl;
        max_dph2_err= std::max( max_dph2_err, dph2_err);
    }
    surfacemeasP1+= quad_2D( std::valarray<double>( 1., cdata.qdom.vertex_size()), cdata.qdom);
    surfacemeasP2+= quad_2D( cdata.qdom_projected.absdets(), cdata.qdom);
}

// Works only if the nodes of the quadrature rule are not shared across facets in QuadDomain2DCL.
void ProjectedQuadDomain2DCL::assign (const SurfacePatchCL& p, const QuadDomain2DCL& qdomarg, const QuaQuaMapperCL& quaqua)
{
    qdom= &qdomarg;
    resize( qdom->vertex_size());
    QuadDomain2DCL::const_vertex_iterator qv= qdom->vertex_begin();

    if (!compute_absdets_) {
        for (Uint i= 0; i < qdom->vertex_size(); ++i, ++qv) {
            quaqua.set_point( quaqua.get_tetra(), *qv)
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
    const Bary2WorldCoordCL b2w( *quaqua.get_tetra());
    const SurfacePatchCL::const_vertex_iterator pv= p.vertex_begin();

    const Uint nodes_per_facet= qdom->vertex_size()/p.facet_size();
    for (Uint i= 0; i < qdom->vertex_size(); ++i, ++qv) {
        if (i % nodes_per_facet == 0) {
            const SurfacePatchCL::FacetT& facet= p.facet_begin()[i/nodes_per_facet];
            M.col( 0, b2w( pv[facet[1]]) - b2w( pv[facet[0]]));
            M.col( 1, b2w( pv[facet[2]]) - b2w( pv[facet[0]]));
            qr.prepare_solve();
            for (Uint j= 0; j < 2; ++j) {
                tmp= std_basis<3>( j + 1);
                qr.apply_Q( tmp);
                U.col( j, tmp);
            }
        }
        quaqua.set_point( quaqua.get_tetra(), *qv)
              .jacobian();
        vertexes_[i]= quaqua.get_base_point();
        Gram= GramMatrix( quaqua.get_jacobian()*U);
        absdets_[i]= std::sqrt( Gram(0,0)*Gram(1,1) - Gram(0,1)*Gram(1,0));
    }
}

void gradient_trafo (const TetraCL& tet, const BaryCoordCL& xb, const QuaQuaMapperCL& quaqua, const SurfacePatchCL& p,
                     SMatrixCL<3,3>& W)
{
    // Compute the basepoint b.
    quaqua.set_point( &tet, xb)
          .jacobian();
    const TetraBaryPairT& tb= quaqua.get_base_point();

    // nl(x)
    if (p.normal_empty())
        p.compute_normals( tet);
    Point3DCL nl= p.normal_begin()[0];
    // Evaluate the normal to the interface in b, n(y).
    Point3DCL n= quaqua.local_ls_grad( *tb.first, tb.second);
    n/= n.norm();

    // Dp_h(x)^T
    const SMatrixCL<3,3>& dph= quaqua.get_jacobian();
    SMatrixCL<3,3> dphT;
    assign_transpose( dphT, dph);

    W= dphT + outer_product( nl, n/inner_prod( nl, n) - dph*nl);
}

void SurfactantP1BaseCL::SetInitialValue (instat_scalar_fun_ptr icf, double t)
{
    P1Init ( icf, ic, MG_, t);
}

void SurfactantP1BaseCL::SetRhs (instat_scalar_fun_ptr rhs)
{
    rhs_fun_= rhs;
}

void SurfactantP1BaseCL::SetTheta (double theta)
{
    if (theta >= 0. && theta <= 1.)
        theta_= theta;
}

void SurfactantP1BaseCL::InitTimeStep ()
{
    if (oldidx_.NumUnknowns() > 0)
        oldidx_.DeleteNumbering( MG_);
    oldidx_.swap( idx);
    oldic_.resize( ic.Data.size());
    oldic_= ic.Data;
    oldt_= ic.t;

    oldls_.RowIdx= lset_vd_.RowIdx;
    oldls_.Data.resize( lset_vd_.Data.size());
    oldls_.Data= lset_vd_.Data;
    oldls_.t= lset_vd_.t;

    oldv_.SetIdx( v_->RowIdx);
    oldv_.Data= v_->Data;
    oldv_.t= v_->t;
}

void SurfactantcGP1CL::Update()
{
    // ScopeTimerCL timer( "SurfactantcGP1CL::Update");
    // std::cout << "SurfactantcGP1CL::Update:\n";

    IdxDescCL* cidx= ic.RowIdx;
    M.Data.clear();
    M.SetIdx( cidx, cidx);
    A.Data.clear();
    A.SetIdx( cidx, cidx);
    C.Data.clear();
    C.SetIdx( cidx, cidx);
    Md.Data.clear();
    Md.SetIdx( cidx, cidx);

    TetraAccumulatorTupleCL accus;
    InterfaceCommonDataP1CL cdata( lset_vd_, lsetbnd_);
    accus.push_back( &cdata);
    InterfaceMatrixAccuCL<LocalInterfaceMassP1CL, InterfaceCommonDataP1CL> mass_accu( &M, LocalInterfaceMassP1CL(), cdata, "mass");
    accus.push_back( &mass_accu);
    InterfaceMatrixAccuCL<LocalLaplaceBeltramiP1CL, InterfaceCommonDataP1CL> lb_accu( &A, LocalLaplaceBeltramiP1CL( D_), cdata, "Laplace-Beltrami");
    accus.push_back( &lb_accu);
    accus.push_back_acquire( make_wind_dependent_matrixP1_accu<LocalInterfaceConvectionP1CL>( &C,  cdata,  make_P2Eval( MG_, Bnd_v_, *v_), "convection"));
    accus.push_back_acquire( make_wind_dependent_matrixP1_accu<LocalInterfaceMassDivP1CL>   ( &Md, cdata,  make_P2Eval( MG_, Bnd_v_, *v_), "massdiv"));

    if (theta_ != 1.0) {
        M2.Data.clear();
        M2.SetIdx( cidx, cidx);
        InterfaceCommonDataP1CL* oldcdata= new InterfaceCommonDataP1CL( oldls_, lsetbnd_);
        accus.push_back_acquire( oldcdata);
        accus.push_back_acquire( new InterfaceMatrixAccuCL<LocalInterfaceMassP1CL, InterfaceCommonDataP1CL>( &M2, LocalInterfaceMassP1CL(), *oldcdata, "old mass"));
    }
    accumulate( accus, MG_, cidx->TriangLevel(), cidx->GetBndInfo());

//     WriteToFile( M.Data, "cGcGM.txt", "mass");
//     WriteToFile( A.Data, "cGcGA.txt", "Laplace-Beltrami");
//     WriteToFile( C.Data, "cGcGC.txt", "material derivative");
//     WriteToFile( Md.Data,"cGcGMd.txt","mass-div");

    // std::cout << "SurfactantP1CL::Update: Finished\n";
}

VectorCL SurfactantcGP1CL::InitStep (double new_t)
{
    // ScopeTimerCL timer( "SurfactantcGP1CL::InitStep");
    // std::cout << "SurfactantcGP1CL::InitStep:\n";

    ic.t= new_t;
    dt_= new_t - oldt_;
    idx.CreateNumbering( oldidx_.TriangLevel(), MG_, &lset_vd_, &lsetbnd_); // InitTimeStep deletes oldidx_ and swaps idx and oldidx_.
    std::cout << "new NumUnknowns: " << idx.NumUnknowns() << std::endl;
    ic.SetIdx( &idx);

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

    if (rhs_fun_)
        accus.push_back_acquire( new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>( &vd_load, LocalVectorP1CL( rhs_fun_, new_t), cdata, "load"));

    if (theta_ == 1.0) {
        accumulate( accus, MG_, idx.TriangLevel(), idx.GetBndInfo());
        return VectorCL( theta_*(vd_timeder.Data + dt_*vd_load.Data));
    }

    InterfaceCommonDataP1CL oldcdata( oldls_, lsetbnd_);
    accus.push_back( &oldcdata);

    if (rhs_fun_)
        accus.push_back_acquire( new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>( &vd_oldload, LocalVectorP1CL( rhs_fun_, oldt_), oldcdata, "load on old iface"));

    InterfaceVectorAccuCL<LocalMatVecP1CL<LocalInterfaceMassP1CL>, InterfaceCommonDataP1CL> old_mass_accu( &vd_oldtimeder,
        LocalMatVecP1CL<LocalInterfaceMassP1CL>( LocalInterfaceMassP1CL(), &vd_oldic), oldcdata, "mixed-mass on old iface");
    accus.push_back( &old_mass_accu);
    InterfaceVectorAccuCL<LocalMatVecP1CL<LocalLaplaceBeltramiP1CL>, InterfaceCommonDataP1CL> old_lb_accu( &vd_oldres,
        LocalMatVecP1CL<LocalLaplaceBeltramiP1CL>( LocalLaplaceBeltramiP1CL( D_), &vd_oldic), oldcdata, "Laplace-Beltrami on old iface");
    accus.push_back( &old_lb_accu);
    accus.push_back_acquire( make_wind_dependent_vectorP1_accu<LocalInterfaceConvectionP1CL>( &vd_oldres, &vd_oldic,  oldcdata,  make_P2Eval( MG_, Bnd_v_, oldv_), "convection on old iface"));
    accus.push_back_acquire( make_wind_dependent_vectorP1_accu<LocalInterfaceMassDivP1CL>   ( &vd_oldres, &vd_oldic,  oldcdata,  make_P2Eval( MG_, Bnd_v_, oldv_), "mass-div on old iface"));

    accumulate( accus, MG_, idx.TriangLevel(), idx.GetBndInfo());
    return VectorCL( theta_*vd_timeder.Data + (1. - theta_)*vd_oldtimeder.Data
                   + dt_*(theta_*vd_load.Data + (1. - theta_)*(vd_oldload.Data - vd_oldres.Data)));
}

void SurfactantcGP1CL::DoStep (const VectorCL& rhs)
{
    Update();

    if (theta_ == 1.)
        L_.LinComb( theta_, M.Data, dt_*theta_, A.Data, dt_*theta_, Md.Data, dt_*theta_, C.Data);
    else {
        MatrixCL m;
        m.LinComb( theta_, M.Data, dt_*theta_, A.Data, dt_*theta_, Md.Data, dt_*theta_, C.Data);
        L_.LinComb( 1., m, 1. - theta_, M2.Data);
    }
    std::cout << "Before solve: res = " << norm( L_*ic.Data - rhs) << std::endl;
    {
        ScopeTimerCL timer( "SurfactantcGP1CL::DoStep: Solve");
        gm_.Solve( L_, ic.Data, rhs, ic.RowIdx->GetEx());
    }
    std::cout << "SurfactantcGP1CL::DoStep: res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
}

void SurfactantcGP1CL::CommitStep ()
{
    return;
}

void SurfactantcGP1CL::DoStep (double new_t)
{
    ScopeTimerCL timer( "SurfactantcGP1CL::DoStep");

    VectorCL rhs( InitStep( new_t));
    DoStep( rhs);
    CommitStep();
}

void
InterfaceP1RepairCL::pre_refine ()
{
    DROPS::NoBndDataCL<> dummy;
    p1repair_= std::unique_ptr<RepairP1CL<double, NoBndDataCL>::type >(
        new RepairP1CL<double, NoBndDataCL>::type( mg_, fullu_, dummy));
}

void
InterfaceP1RepairCL::post_refine ()
{
    VecDescCL loc_u;
    IdxDescCL loc_idx( P1_FE);

    loc_idx.CreateNumbering( fullp1idx_.TriangLevel(), mg_);
//    loc_idx.CreateNumbering( mg_.GetLastLevel(), mg_);
    loc_u.SetIdx( &loc_idx);
    DROPS::NoBndDataCL<> dummy;

//    P1EvalCL<double, DROPS::NoBndDataCL<>, const VecDescCL> oldsol( &fullu_, &dummy, &mg_);
//    P1EvalCL<double, DROPS::NoBndDataCL<>,       VecDescCL>    sol( &loc_u , &dummy, &mg_);
//    Interpolate( sol, oldsol);

    p1repair_->repair( loc_u);

    fullu_.Clear( fullu_.t);
    fullp1idx_.swap( loc_idx);
    loc_idx.DeleteNumbering( mg_);
    fullu_.SetIdx( &fullp1idx_);
    fullu_.Data= loc_u.Data;
}

void
InterfaceP1RepairCL::pre_refine_sequence ()
{
    fullp1idx_.CreateNumbering( u_.RowIdx->TriangLevel(), mg_);
    fullu_.SetIdx( &fullp1idx_);
    Extend( mg_, u_, fullu_);
}

void
InterfaceP1RepairCL::post_refine_sequence ()
{
    u_.RowIdx->DeleteNumbering( mg_);
    u_.RowIdx->CreateNumbering( fullp1idx_.TriangLevel(), mg_, &lset_vd_, &lset_bnd_);
    u_.SetIdx( u_.RowIdx);

    Restrict( mg_, fullu_, u_);

    fullp1idx_.DeleteNumbering( mg_);
    fullu_.Clear( u_.t);
}
#ifndef _PAR
void
Ensight6IfaceScalarCL::put (Ensight6OutCL& cf) const
{
    IdxDescCL p1idx;
    p1idx.CreateNumbering( u_.RowIdx->TriangLevel(), mg_);
    VecDescCL p1u( &p1idx);
    Extend( mg_, u_, p1u);
    BndDataCL<> bnd( 0);

    cf.putScalar( make_P1Eval( mg_, bnd, p1u), varName());

    p1idx.DeleteNumbering( mg_);
}
#endif

void
VTKIfaceScalarCL::put (VTKOutCL& cf) const
{
    IdxDescCL fullidx( P1_FE);
    if (u_.RowIdx->GetFE() == P2IF_FE)
        fullidx.SetFE( P2_FE);
    fullidx.CreateNumbering( u_.RowIdx->TriangLevel(), mg_);
    std::cout<<"vtk scalar unknowns: "<<fullidx.NumUnknowns()<<std::endl;
    VecDescCL uext( &fullidx);
    Extend( mg_, u_, uext);

    if (u_.RowIdx->GetFE() == P1IF_FE)
        cf.PutScalar( make_P1Eval( mg_, BndData_, uext), varName());
    else if (u_.RowIdx->GetFE() == P2IF_FE)
        cf.PutScalar( make_P2Eval( mg_, BndData_, uext), varName());

    fullidx.DeleteNumbering( mg_);
}

void VTKIfaceVectorCL::put (VTKOutCL& cf) const
{
    IdxDescCL pidx;
    if(P2_)
      pidx.SetFE(vecP2_FE);
    else
      pidx.SetFE(vecP1_FE);
    pidx.CreateNumbering( u_.RowIdx->TriangLevel(), mg_);
    std::cout<<"vtk vector unknowns: "<<pidx.NumUnknowns()<<std::endl;
    VecDescCL pu( &pidx);
    Extend( mg_, u_, pu);

    P2EvalCL<SVectorCL<3>, const BndDataCL<Point3DCL>, const VecDescBaseCL<VectorCL> > eval2(&pu, &BndData_, &mg_);
    P1EvalCL<SVectorCL<3>, const BndDataCL<Point3DCL>, const VecDescBaseCL<VectorCL> > eval1(&pu, &BndData_, &mg_);
    if(P2_)
      cf.PutVector( eval2, varName());
    else
      cf.PutVector( eval1, varName());

    pidx.DeleteNumbering( mg_);
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

        const Point3DCL p1grad( FE_P1CL::DHRef( i));
        const Point4DCL& p1grad4d= MakePoint4D( p1grad[0], p1grad[1], p1grad[2], 0.);
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

    IdxDescCL idx_ini( P1IF_FE);
    idx_ini.CreateNumbering( level, mg, &oldls, &lsetbnd);
    NumIniUnknowns_= idx_ini.NumUnknowns();

    IdxDescCL idx_fini( P1IF_FE);
    idx_fini.CreateNumbering( level, mg, &newls, &lsetbnd);
    NumFiniUnknowns_= idx_fini.NumUnknowns();

    NumUnknowns_= NumIniUnknowns() + NumFiniUnknowns();
    IdxT inicounter= 0, finicounter= 0;

    const TetraPrismLatticeCL& lat= TetraPrismLatticeCL::instance( 2, 1);
    std::valarray<double> ls_loc( lat.vertex_size());
    LocalP2CL<> oldlocp2_ls, locp2_ls;
    LocalSTP2P1CL<> local_st_lset;
    SPatchCL<4> patch;
    QuadDomainCodim1CL<4> qdom;
    std::valarray<double> shape; // STP1P1-shape-function as integrand
    DROPS_FOR_TRIANG_TETRA ( mg, level, it) {
        local_st_lset.assign( *it, oldls, newls, lsetbnd);
        evaluate_on_vertexes( local_st_lset, lat, Addr( ls_loc));
        if (equal_signs( ls_loc))
            continue;

        patch.make_patch<MergeCutPolicyCL>( lat, ls_loc);
        make_CompositeQuad5DomainSTCodim1( qdom, patch, TetraPrismCL( *it, t0, t1));
        shape.resize( qdom.vertex_size());
        for (Uint i= 0; i < 8; ++i) {
            const IdxDescCL& idx= i < 4 ? idx0_ : idx1_;
            const IdxDescCL& spatial_idx= i < 4 ? idx_ini : idx_fini;
            UnknownHandleCL& unknowns= const_cast<VertexCL*>( it->GetVertex( i%4))->Unknowns;
            if (unknowns.Exist( idx.GetIdx()))
                continue;

            evaluate_on_vertexes( STP1P1DiscCL::ref_val[i], qdom, Addr( shape));
            if (quad_codim1( shape*shape, qdom) > 0.) {
                unknowns.Prepare( idx.GetIdx());
                if (unknowns.Exist( spatial_idx.GetIdx())) {
                    unknowns( idx.GetIdx())= unknowns( spatial_idx.GetIdx()) + (i < 4 ? 0 : NumIniUnknowns());
                    ++(i < 4 ? inicounter : finicounter);
                }
                else
                    unknowns( idx.GetIdx())= NumUnknowns_++;
            }
        }
    }
    if (inicounter != NumIniUnknowns())
        throw DROPSErrCL( "STP1P1IdxDescCL::CreateNumbering: Wrong count of the unknowns on the old interface.\n");
    if (finicounter != NumFiniUnknowns())
        throw DROPSErrCL( "STP1P1IdxDescCL::CreateNumbering: Wrong count of the unknowns on the new interface.\n");
    idx_ini.DeleteNumbering( mg);
    idx_fini.DeleteNumbering( mg);
}


void
resize_and_scatter_piecewise_spatial_normal (const SPatchCL<4>& surf, const QuadDomainCodim1CL<4>& qdom, std::valarray<Point3DCL>& spatial_normal)
{
    spatial_normal.resize( qdom.vertex_size());
    if (spatial_normal.size() == 0)
        return;
    if (surf.normal_empty()) // As qdom has vertexes, the must be facets, i.e. normals.
        throw DROPSErrCL( "resize_and_scatter_piecewise_spatial_normal: normals were not precomputed.\n");

    const Uint NodesPerFacet= qdom.vertex_size()/surf.facet_size();
    if (qdom.vertex_size()%surf.facet_size() != 0)
        throw DROPSErrCL( "resize_and_scatter_piecewise_spatial_normal: qdom.vertex_size is not a multiple of surf.facet_size.\n");

    const SPatchCL<4>::const_normal_iterator n= surf.normal_begin();
    for (Uint i= 0; i < surf.facet_size(); ++i) {
        const Point3DCL& tmp= MakePoint3D( n[i][0], n[i][1], n[i][2]);
        std::fill_n( &spatial_normal[i*NodesPerFacet], NodesPerFacet, tmp/tmp.norm());
    }
}


void SurfactantSTP1CL::InitStep (double new_t)
{
    // std::cout << "SurfactantSTP1CL::InitStep:\n";
    ic.t= new_t;
    dt_= ic.t - oldt_;

    st_idx_.CreateNumbering( oldidx_.TriangLevel(), MG_, oldls_, lset_vd_, lsetbnd_, oldt_, new_t);
    dim= st_idx_.NumUnknowns() - (cG_in_t_ ? st_idx_.NumIniUnknowns() : 0);
    std::cout << "SurfactantSTP1CL::InitStep: space-time Unknowns: " << st_idx_.NumUnknowns()
              << " ini: " << st_idx_.NumIniUnknowns() << " fini: " << st_idx_.NumFiniUnknowns()
              << " interior: " << st_idx_.NumUnknowns() - st_idx_.NumIniUnknowns() - st_idx_.NumFiniUnknowns()
              << " dimension of the linear system: " << dim << std::endl;
    st_ic_.resize( dim);

    st_oldic_.resize( st_idx_.NumUnknowns());
    // Copy dofs on the old interface from  old solution into st_oldic_.
    DROPS_FOR_TRIANG_VERTEX( MG_, oldidx_.TriangLevel(), it) {
        if (it->Unknowns.Exist( st_idx_.GetIdx( 0))) {
            const IdxT dof= it->Unknowns( st_idx_.GetIdx( 0));
            if (dof < st_idx_.NumIniUnknowns())
                st_oldic_[dof]= oldic_[dof];
        }
    }
}

void SurfactantSTP1CL::Update_cG()
{
    TetraAccumulatorTupleCL accus;
    STInterfaceCommonDataCL cdata( oldt_, ic.t,  oldls_, lset_vd_, lsetbnd_);
    accus.push_back( &cdata);
    InterfaceMatrixSTP1AccuCL<LocalLaplaceBeltramiSTP1P1CL> lb_accu( &A, &cpl_A_, &st_idx_,
        LocalLaplaceBeltramiSTP1P1CL( D_), cdata, &st_oldic_, "Laplace-Beltrami on ST-iface");
    accus.push_back( &lb_accu);

    cpl_A_.resize( dim);
    cpl_der_.resize( dim);

    VectorCL cpl_new_dummy;
    InterfaceCommonDataP1CL oldspatialcdata( oldls_, lsetbnd_);
    InterfaceMatrixSTP1AccuCL<LocalSpatialInterfaceMassSTP1P1CL> oldmass_accu( &Mold, &cpl_old_, &st_idx_,
        LocalSpatialInterfaceMassSTP1P1CL( oldspatialcdata), cdata, &st_oldic_, "mixed-mass on old iface");
    InterfaceCommonDataP1CL newspatialcdata( lset_vd_, lsetbnd_);
    InterfaceMatrixSTP1AccuCL<LocalSpatialInterfaceMassSTP1P1CL> newmass_accu( &Mnew, /*dummy*/ &cpl_new_dummy, &st_idx_,
        LocalSpatialInterfaceMassSTP1P1CL( newspatialcdata, false), cdata, &st_oldic_, "mixed-mass on new iface");

    if (use_mass_div_) {
        cpl_div_.resize( dim);
        accus.push_back_acquire( make_wind_dependent_matrixSTP1P0_1_accu<LocalMassdivSTP1P1CL>( &Mdiv, &cpl_div_, &st_idx_, cdata, &st_oldic_, make_STP2P1Eval( MG_, Bnd_v_, oldv_, *v_), "mass-div on ST-iface"));
        accus.push_back_acquire( make_wind_dependent_matrixSTP1P0_1_accu<LocalMaterialDerivativeSTP1P1CL>( &Mder, &cpl_der_, &st_idx_, cdata, &st_oldic_, make_STP2P1Eval( MG_, Bnd_v_, oldv_, *v_), "material derivative on ST-iface"));
    }
    else {
        cpl_old_.resize( dim);
        cpl_new_dummy.resize( dim);
        accus.push_back_acquire( make_wind_dependent_local_transpose_matrixSTP1P0_1_accu<LocalMaterialDerivativeSTP1P1CL>( &Mder, &cpl_der_, &st_idx_, cdata, &st_oldic_, make_STP2P1Eval( MG_, Bnd_v_, oldv_, *v_), "material derivative on ST-iface"));
        accus.push_back( &oldspatialcdata);
        accus.push_back( &oldmass_accu);
        accus.push_back( &newspatialcdata);
        accus.push_back( &newmass_accu);
    }

    if (rhs_fun_) {
        load.resize( dim);
        accus.push_back_acquire( new InterfaceVectorSTP1AccuCL<LocalVectorSTP1P1CL>( &load, &st_idx_, LocalVectorSTP1P1CL( rhs_fun_), cdata, /* cG_in_t*/ cG_in_t_, "load on ST-iface"));
    }

    accumulate( accus, MG_, st_idx_.TriangLevel(), idx.GetBndInfo());

    // WriteToFile( cpl_new_dummy, "cpl_new_dummy.txt", "coupling on new interface -- always zero.");
}

void SurfactantSTP1CL::Update_dG()
{
    TetraAccumulatorTupleCL accus;
    InterfaceCommonDataP1CL oldspatialcdata( oldls_, lsetbnd_);
    accus.push_back( &oldspatialcdata);
    STInterfaceCommonDataCL cdata( oldt_, ic.t,  oldls_, lset_vd_, lsetbnd_);
    accus.push_back( &cdata);
    InterfaceMatrixSTP1AccuCL<LocalSpatialInterfaceMassSTP1P1CL> oldmass_accu( &Mold, &st_idx_,
        LocalSpatialInterfaceMassSTP1P1CL(oldspatialcdata), cdata, "mixed-mass on old iface");
    accus.push_back( &oldmass_accu);
    InterfaceMatrixSTP1AccuCL<LocalLaplaceBeltramiSTP1P1CL> lb_accu( &A, &st_idx_,
        LocalLaplaceBeltramiSTP1P1CL( D_), cdata, "Laplace-Beltrami on ST-iface");
    accus.push_back( &lb_accu);

    InterfaceCommonDataP1CL newspatialcdata( lset_vd_, lsetbnd_);
    InterfaceMatrixSTP1AccuCL<LocalSpatialInterfaceMassSTP1P1CL> newmass_accu( &Mnew, &st_idx_,
        LocalSpatialInterfaceMassSTP1P1CL( newspatialcdata, false), cdata, "mixed-mass on new iface");
    if (use_mass_div_) {
        accus.push_back_acquire( make_wind_dependent_matrixSTP1P1_accu<LocalMaterialDerivativeSTP1P1CL>( &Mder, &st_idx_, cdata,  make_STP2P1Eval( MG_, Bnd_v_, oldv_, *v_), "material derivative on ST-iface"));
        accus.push_back_acquire( make_wind_dependent_matrixSTP1P1_accu<LocalMassdivSTP1P1CL>( &Mdiv, &st_idx_, cdata,  make_STP2P1Eval( MG_, Bnd_v_, oldv_, *v_), "mass-div on ST-iface"));
    }
    else {
        accus.push_back_acquire( make_wind_dependent_local_transpose_matrixSTP1P1_accu<LocalMaterialDerivativeSTP1P1CL>( &Mder, &st_idx_, cdata,  make_STP2P1Eval( MG_, Bnd_v_, oldv_, *v_), "material derivative on ST-iface"));
        accus.push_back( &newspatialcdata);
        accus.push_back( &newmass_accu);
    }

    if (rhs_fun_) {
        load.resize( dim);
        accus.push_back_acquire( new InterfaceVectorSTP1AccuCL<LocalVectorSTP1P1CL>( &load, &st_idx_, LocalVectorSTP1P1CL( rhs_fun_), cdata, /* cG_in_t*/ cG_in_t_, "load on ST-iface"));
    }

    accumulate( accus, MG_, st_idx_.TriangLevel(), idx.GetBndInfo());
}

void SurfactantSTP1CL::Update()
{
    ScopeTimerCL timer( "SurfactantSTP1CL::Update");
    // std::cout << "SurfactantSTP1CL::Update:\n";
    if (cG_in_t_)
        Update_cG();
    else
        Update_dG();

//     WriteToFile( Mold, "Mold.txt", "mass on old iface");
//     WriteToFile( Mnew, "Mnew.txt", "mass on new iface");
//     WriteToFile( A,    "A.txt",    "Laplace-Beltrami on ST-iface");
//     WriteToFile( Mder, "Mder.txt", "material derivative on ST-iface");
//     WriteToFile( Mdiv, "Mdiv.txt", "mass-div on ST-iface");
//     WriteToFile( load, "load.txt", "load on ST-iface");
//     WriteToFile( cpl_A_,   "cpl_A.txt",   "coupling for Laplace-Beltrami on ST-iface");
//     WriteToFile( cpl_der_, "cpl_der.txt", "coupling for material derivative on ST-iface");
//     WriteToFile( cpl_div_, "cpl_div.txt", "coupling for mass-div on ST-iface");
//     WriteToFile( cpl_old_, "cpl_old.txt", "coupling ini-values on old iface");

    // std::cout << "SurfactantSTP1CL::Update: Finished\n";
}

void SurfactantSTP1CL::DoStep ()
{
    Update();

    MatrixCL L;
    VectorCL rhs( dim);
    if (cG_in_t_) {
        if (use_mass_div_) {
            L.LinComb( 1., Mder, 1., Mdiv, 1., A);
            rhs= cpl_der_ + cpl_div_ + cpl_A_;
        }
        else {
            L.LinComb( -1., Mder, 1., A, 1., Mnew);
            rhs= -cpl_der_ + cpl_A_ - cpl_old_;
        }
    }
    else {
        if (use_mass_div_) {
            L.LinComb( 1., Mder, 1., Mdiv, 1., A, 1., Mold);
            rhs= Mold*st_oldic_;
        }
        else {
            L.LinComb( -1., Mder, 1., A, 1., Mnew);
            rhs= Mold*st_oldic_;
        }
    }
    if (rhs_fun_ != 0)
        rhs+= load;

    std::cout << "Before solve: res = " << norm( L*st_ic_ - rhs) << std::endl;

    {
        ScopeTimerCL timer( "SurfactantSTP1CL::DoStep: Solve");
        gm_.Solve( L, st_ic_, rhs, idx.GetEx());
    }
    std::cout << "SurfactantSTP1CL::DoStep: res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
}

void SurfactantSTP1CL::CommitStep ()
{
    idx.CreateNumbering( oldidx_.TriangLevel(), MG_, &lset_vd_, &lsetbnd_);
    std::cout << "new NumUnknowns at t1: " << idx.NumUnknowns() << std::endl;

    ic.SetIdx( &idx);
    // Copy dofs on the new interface from space-time-solution into ic.
    DROPS_FOR_TRIANG_VERTEX( MG_, oldidx_.TriangLevel(), it) {
        if (it->Unknowns.Exist( st_idx_.GetIdx( 1))) {
            const IdxT dof= it->Unknowns( st_idx_.GetIdx( 1));
            if (dof >= st_idx_.NumIniUnknowns() && dof < st_idx_.NumIniUnknowns() + st_idx_.NumFiniUnknowns())
                ic.Data[dof - st_idx_.NumIniUnknowns()]= st_ic_[dof - (cG_in_t_ ? st_idx_.NumIniUnknowns() : 0)];
        }
    }

    st_ic_.resize( 0);
    st_oldic_.resize( 0);
    st_idx_.DeleteNumbering( MG_);
}

void SurfactantSTP1CL::DoStep (double new_t)
{
    ScopeTimerCL timer( "SurfactantSTP1CL::DoStep");

    InitStep( new_t);
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
void TransportP1FunctionCL::SetupSystem (const DiscVelSolT& vel, MatrixCL& E, MatrixCL& H)
// Sets up the stiffness matrices:
// E is of mass matrix type:    E_ij = ( v_j       , v_i + SD * u grad v_i )
// H describes the convection:  H_ij = ( u grad v_j, v_i + SD * u grad v_i )
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
        p1[i].assign( p1dummy);
        p1dummy[i]= 0.0;
    }

    DROPS_FOR_TRIANG_CONST_TETRA( const_cast<const MultiGridCL&>( MG_), lvl, sit) {
        P1DiscCL::GetGradients( Grad, det, *sit);
        absdet= std::fabs( det);
        h_T= std::pow( absdet, 1./3.);

        GetLocalNumbP1NoBnd( Numb, *sit, p1idx);

        u_loc.assign( *sit, vel);
        for(int i=0; i<4; ++i)
            u_Grad[i]= dot( Grad[i], u_loc);

        /// \todo fixed limit for maxV (maxV_limit), any better idea?
        double maxV = 1.; // scaling of SD parameter
        for(int i= 0; i < 4; ++i)    // assemble row Numb[i]
            for(int j= 0; j < 4; ++j) {
                // E is of mass matrix type:    E_ij = ( v_j       , v_i + SD * u grad v_i )
               bE( Numb[i], Numb[j])+= P1DiscCL::GetMass(i,j) * absdet
                                       + Quad5CL<>( u_Grad[i]*p1[j]).quad( absdet)*SD_/maxV*h_T;

               // H describes the convection:  H_ij = ( u grad v_j, v_i + SD * u grad v_i )
               bH( Numb[i], Numb[j])+= Quad5CL<>( u_Grad[j]*p1[i]).quad( absdet)
                                       + Quad5CL<>(u_Grad[i]*u_Grad[j]).quad( absdet) * SD_/maxV*h_T;
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
    SetupSystem( vel, E_, H_);
    L_.clear();
    L_.LinComb( 1./dt_, E_, theta_, H_);

    VectorCL rhs( (1./dt_)*u);
    if (theta_ != 1.) {
        GMResSolverCL<GSPcCL> gm( gm_);
        VectorCL tmp( rhs.size());
        gm.Solve( E_old, tmp, VectorCL( H_old*u), p1idx.GetEx());
        std::cout << "TransportP1FunctionCL::DoStep rhs: res = " << gm.GetResid() << ", iter = " << gm.GetIter() << std::endl;
        rhs-= (1. - theta_)*tmp;
    }
    gm_.Solve( L_, u, VectorCL(E_*rhs), p1idx.GetEx());

    std::cout << "TransportP1FunctionCL::DoStep: res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
}


void SurfactantCharTransportP1CL::Update()
{
    // ScopeTimerCL timer( "SurfactantCharTransportP1CL::Update");
    // std::cout << "SurfactantCharTransportP1CL::Update:\n";

    IdxDescCL* cidx= ic.RowIdx;
    M.Data.clear();
    M.SetIdx( cidx, cidx);
    A.Data.clear();
    A.SetIdx( cidx, cidx);
    Md.Data.clear();
    Md.SetIdx( cidx, cidx);
    VecDescCL vd_load( &idx);

    TetraAccumulatorTupleCL accus;
    InterfaceCommonDataP1CL cdata( lset_vd_, lsetbnd_);
    accus.push_back( &cdata);
    InterfaceMatrixAccuCL<LocalInterfaceMassP1CL, InterfaceCommonDataP1CL> mass_accu( &M, LocalInterfaceMassP1CL(), cdata, "mass");
    accus.push_back( &mass_accu);
    InterfaceMatrixAccuCL<LocalLaplaceBeltramiP1CL, InterfaceCommonDataP1CL> lb_accu( &A, LocalLaplaceBeltramiP1CL( D_), cdata, "Laplace-Beltrami");
    accus.push_back( &lb_accu);
    accus.push_back_acquire( make_wind_dependent_matrixP1_accu<LocalInterfaceMassDivP1CL>( &Md, cdata,  make_P2Eval( MG_, Bnd_v_, *v_), "massdiv"));
    if (rhs_fun_)
        accus.push_back_acquire( new InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL>( &vd_load, LocalVectorP1CL( rhs_fun_, ic.t), cdata, "load"));

    accumulate( accus, MG_, cidx->TriangLevel(), cidx->GetBndInfo());

    load.resize( idx.NumUnknowns());
    load= vd_load.Data;

//     WriteToFile( M.Data, "chartranspM.txt", "mass");
//     WriteToFile( A.Data, "chartranspA.txt", "Laplace-Beltrami");
//     WriteToFile( Md.Data,"chartranspMd.txt","mass-div");
//     WriteToFile( vd_load.Data,"chartranspload.txt","load");

    // std::cout << "SurfactantCharTransportP1CL::Update: Finished\n";
}

void SurfactantCharTransportP1CL::InitStep (double new_t)
{
    // ScopeTimerCL timer( "SurfactantcGP1CL::InitStep");
    std::cout << "SurfactantCharTransportP1CL::InitStep:\n";

    ic.t= new_t;
    dt_= ic.t - oldt_;
    idx.CreateNumbering( oldidx_.TriangLevel(), MG_, &lset_vd_, &lsetbnd_); // InitTimeStep deletes oldidx_ and swaps idx and oldidx_.
    std::cout << "new NumUnknowns: " << idx.NumUnknowns();
    full_idx.CreateNumbering( idx.TriangLevel(), MG_);
    std::cout << " full NumUnknowns: " << full_idx.NumUnknowns() << std::endl;

    fulltransport_= new TransportP1FunctionCL( MG_, make_P2Eval( MG_, Bnd_v_, oldv_), full_idx, dt_,  /*theta=*/theta_, /*SD=*/ 0.1, /*iter=*/ 2000, /*tol=*/ 0.1*gm_.GetTol());

    VecDescCL rhs( &idx);
    rhs.Data= (1./dt_)*ic.Data;
//     if (theta_ != 1.) {
//         GMResSolverCL<GSPcCL> gm( gm_);
//         VectorCL tmp( rhs.Data.size());
//         gm.Solve( M.Data, tmp, VectorCL( A.Data*ic.Data + Md.Data*ic.Data));
//         std::cout << "SurfactantP1CL::InitStep: rhs: res = " << gm.GetResid() << ", iter = " << gm.GetIter() << std::endl;
//         rhs.Data-= (1. - theta_)*tmp;
//     }
    DROPS::VecDescCL rhsext( &full_idx);
    DROPS::Extend( MG_, rhs, rhsext);
    rhs_.resize( rhsext.Data.size());
    rhs_= rhsext.Data;
}

void SurfactantCharTransportP1CL::DoStep ()
{
    VecDescCL transp_rhs( &idx),
              transp_rhsext( &full_idx);
    transp_rhsext.Data= rhs_;
    fulltransport_->DoStep( transp_rhsext.Data, make_P2Eval( MG_, Bnd_v_, *v_));
    Restrict( MG_, transp_rhsext, transp_rhs);

    ic.SetIdx( &idx);
    Update();
    // L_.LinComb( 1./dt_, M.Data, theta_, A.Data, theta_, Md.Data);
    L_.LinComb( 1./dt_, M.Data, 1., A.Data, 1., Md.Data);
    const VectorCL therhs( M.Data*transp_rhs.Data + load);
    std::cout << "Before solve: res = " << norm( L_*ic.Data - therhs) << std::endl;
    {
        ScopeTimerCL timer( "SurfactantCharTransportP1CL::DoStep: Solve");
        gm_.Solve( L_, ic.Data, therhs, idx.GetEx());
    }
    std::cout << "SurfactantCharTransportP1CL::DoStep: res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
}

void SurfactantCharTransportP1CL::CommitStep ()
{
    full_idx.DeleteNumbering( MG_);
    delete fulltransport_;
}

void SurfactantCharTransportP1CL::DoStep (double new_t)
{
    ScopeTimerCL timer( "SurfactantCharTransportP1CL::::DoStep");

    InitStep( new_t);
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
    GridFunctionCL<Point3DCL> qsurfP1grad[4];
    GridFunctionCL<Point3DCL> q3DsurfP1grad[4];

    GridFunctionCL<> mobility2D,well_potential_second_derivative[4], qP1Hat[4], qnormal_comp[3];


    QuadDomain2DCL  q2Ddomain;
    std::valarray<double> ls_loc;
    SurfacePatchCL spatch;

    QuadDomainCL q3Ddomain;
    GridFunctionCL<Point3DCL> q3Dnormal;
    GridFunctionCL<Point3DCL> q3DP1Grad[4],q2DP1Grad[4];

    void Get_Normals(const LocalP2CL<>& ls, LocalP1CL<Point3DCL>&);

  public:
    LocalCahnHilliardCL ( )
        : lat( PrincipalLatticeCL::instance( 2)), ls_loc( lat.vertex_size())
    {
        P1DiscCL::GetP1Basis( P1Hat);
        P2DiscCL::GetP2Basis( P2Hat);
        P2DiscCL::GetGradientsOnRef( P2GradRef);
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
    Get_Normals(ls, Normals);

    resize_and_evaluate_on_vertexes (Normals, q3Ddomain, q3Dnormal);

    // Scale Normals accordingly to the Euclidean Norm (only consider the ones which make a contribution in the sense of them being big enough... otherwise one has to expect problems with division through small numbers)
    for(Uint i=0; i<q3Dnormal.size(); ++i) {
         //if(q3Dnormal[i].norm()> 1e-8)
         q3Dnormal[i]= q3Dnormal[i]/q3Dnormal[i].norm();
    }

    for(int j=0; j<4; ++j) {
        q3DP1Grad[j].resize( q3Ddomain.vertex_size());
        q3DP1Grad[j]= P1Grad[j];
    }
}

double Mobility_function(double x)
{
	return(std::sqrt((1.-x)*(x)*(1.-x)*(x)));
	//return(1);
}

double Potential_function(const double x)
    {
        return((1.-x)*(x)*(1.-x)*(x)*0.25);
    }

double Potential_prime_function(const double x)
{
    //return(-x);
    return((x-1.)*(x)*(x-0.5));
    //return(-(0.)*(3./2.)*x*x + (1.)*x/2. + (0.)*x*x*x);
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

    P2DiscCL::GetGradients( P2Grad, P2GradRef, T);
    P1DiscCL::GetGradients( P1Grad, T);
    evaluate_on_vertexes( ls, lat, Addr( ls_loc));
    spatch.make_patch<MergeCutPolicyCL>( lat, ls_loc);

    // The routine takes the information about the tetrahedra and the cutting surface and generates a two-dimensional triangulation of the cut, including the necessary point-positions and weights for the quadrature
    make_CompositeQuad5Domain2D ( q2Ddomain, spatch, tet);
    LocalP1CL<Point3DCL> Normals;
    Get_Normals(ls, Normals);

    // Resize and evaluate Normals at all points which are needed for the two-dimensional quadrature-rule
    resize_and_evaluate_on_vertexes (Normals, q2Ddomain, qnormal);

    for(int j=0; j<4; ++j) {
    	q2DP1Grad[j].resize( q2Ddomain.vertex_size());
        q2DP1Grad[j]= P1Grad[j];
    }

    // Scale Normals accordingly to the Euclidean Norm (only consider the ones which make a contribution in the sense of them being big enough... otherwise one has to expect problems with division through small numbers)
    for(Uint i=0; i<qnormal.size(); ++i) {
         //if(qnormal[i].norm()> 1e-8)
         qnormal[i]= qnormal[i]/qnormal[i].norm();
    }

    // Provide all components of the normals
    for (int k=0; k<3; ++k) {
        qnormal_comp[k].resize( q2Ddomain.vertex_size());
        ExtractComponent( qnormal, qnormal_comp[k], k);
    }

    for(int j=0; j<4 ;++j) {
        resize_and_evaluate_on_vertexes( P1Hat[j], q2Ddomain, qP1Hat[j]);
        qsurfP1grad[j].resize( q2Ddomain.vertex_size());
        qsurfP1grad[j]= P1Grad[j];
        qsurfP1grad[j]-= dot( qsurfP1grad[j], qnormal)*qnormal;
    }

    for(int j=0; j<4 ;++j) {
        q3DsurfP1grad[j].resize( q3Ddomain.vertex_size());
        q3DsurfP1grad[j]= P1Grad[j];
        q3DsurfP1grad[j]-= dot( q3DsurfP1grad[j], q3Dnormal)*q3Dnormal;
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
void LocalCahnHilliardCL::Get_Normals(const LocalP2CL<>& ls, LocalP1CL<Point3DCL>& Normals)
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
             M_P1[i][j]= quad_2D( qP1Hat[i]*qP1Hat[j], q2Ddomain);
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
             L_P1[i][j]= quad_2D( dot(qsurfP1grad[i], qsurfP1grad[j]), q2Ddomain);
        }
    }
}

void LocalCahnHilliardCL::setupLM_P1 (double LM_P1[4][4])
{
    // Do all combinations for (i,j) i,j=4 x 4 and corresponding quadrature
    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
             LM_P1[i][j]= quad_2D( mobility2D*dot(qsurfP1grad[i], qsurfP1grad[j]), q2Ddomain);
        }
    }
}

void LocalCahnHilliardCL::setupGprimeprime_P1 (double Gprimeprime_P1[4][4])
    {
        // Do all combinations for (i,j) i,j=4 x 4 and corresponding quadrature
        for (int i=0; i<4; ++i) {
            for (int j=0; j<4; ++j) {
                Gprimeprime_P1[i][j]= quad_2D( well_potential_second_derivative[i]*qP1Hat[j], q2Ddomain);
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
    ///\brief Update the global system.
    void update_global_system ();

  public:
    CahnHilliardIFAccumulator_P1P1CL ( const VecDescCL& ls, const VecDescCL& v,const VecDescCL& chi, const LsetBndDataCL& ls_bnd, const BndDataCL<Point3DCL>& v_bnd,const BndDataCL<>& chi_bnd, IdxDescCL& P1FE, IdxDescCL& ScalarP1FE,  MatrixCL& M_P1, MatrixCL& NormalStab_P1,MatrixCL& TangentStab_P1, MatrixCL& VolumeStab_P1,  MatrixCL& L_P1P1,  MatrixCL& LM_P1P1 ,  MatrixCL& Gprimeprime_P1P1);

    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    void visit (const TetraCL& sit);

    TetraAccumulatorCL* clone (int /*tid*/) { return new CahnHilliardIFAccumulator_P1P1CL ( *this); }

};

CahnHilliardIFAccumulator_P1P1CL::CahnHilliardIFAccumulator_P1P1CL( const VecDescCL& ls, const VecDescCL& v, const VecDescCL& chi, const LsetBndDataCL& ls_bnd, const BndDataCL<Point3DCL>& v_bnd,const BndDataCL<>& chi_bnd, IdxDescCL& P1FE, IdxDescCL& ScalarP1FE, MatrixCL& M_P1, MatrixCL& NormalStab_P1, MatrixCL& TangentStab_P1,  MatrixCL& VolumeStab_P1,  MatrixCL& L_P1P1, MatrixCL& LM_P1P1, MatrixCL& Gprimeprime_P1P1 )
 : lset(ls), velocity(v), volume_fraction(chi), lset_bnd(ls_bnd), velocity_bnd(v_bnd),volume_fraction_bnd(chi_bnd), P1Idx_(P1FE), ScalarP1Idx_(ScalarP1FE), M_P1_(M_P1), NormalStab_P1_(NormalStab_P1), TangentStab_P1_(TangentStab_P1), VolumeStab_P1_(VolumeStab_P1), L_P1P1_(L_P1P1), LM_P1P1_(LM_P1P1), Gprimeprime_P1P1_(Gprimeprime_P1P1), localCahnHilliard_()
{}

void CahnHilliardIFAccumulator_P1P1CL::begin_accumulation ()
{
    std::cout << "entering CahnHilliardIF: \n";
    const size_t num_unks_p1= P1Idx_.NumUnknowns(), num_unks_scalarp1= ScalarP1Idx_.NumUnknowns();

    mM_P1_= new MatrixBuilderCL( &M_P1_, num_unks_scalarp1, num_unks_scalarp1);
    mNormalStab_P1_= new MatrixBuilderCL( &NormalStab_P1_, num_unks_scalarp1, num_unks_scalarp1);
    mTangentStab_P1_= new MatrixBuilderCL( &TangentStab_P1_, num_unks_scalarp1, num_unks_scalarp1);
    mVolumeStab_P1_= new MatrixBuilderCL( &VolumeStab_P1_, num_unks_scalarp1, num_unks_scalarp1);
    mL_P1P1_= new MatrixBuilderCL( &L_P1P1_, num_unks_scalarp1, num_unks_scalarp1);
    mLM_P1P1_= new MatrixBuilderCL( &LM_P1P1_, num_unks_scalarp1, num_unks_scalarp1);
    mGprimeprime_P1P1_= new MatrixBuilderCL( &Gprimeprime_P1P1_, num_unks_scalarp1, num_unks_scalarp1);



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
    ls_loc.assign( tet, lset, lset_bnd);
    v_loc.assign( tet, velocity, velocity_bnd);
    chi_loc.assign( tet, volume_fraction, volume_fraction_bnd);

    if (!equal_signs( ls_loc))
    {
        local_setup( tet);
        update_global_system();
    }

}

void CahnHilliardIFAccumulator_P1P1CL::local_setup (const TetraCL& tet)
{
    GetTrafoTr( T, det, tet);
    absdet= std::fabs( det);

    n.assign( tet, P1Idx_, P1Idx_.GetBndInfo());
    nScalar.assign( tet, ScalarP1Idx_, ScalarP1Idx_.GetBndInfo());


    localCahnHilliard_.calcIntegrands( T, ls_loc, v_loc, chi_loc, tet);
    localCahnHilliard_.calc3DIntegrands(T, ls_loc, tet);

    localCahnHilliard_.setupM_P1( locM_P1);
    localCahnHilliard_.setupNormalStab_P1( locNormalStab_P1, absdet);
    localCahnHilliard_.setupTangentStab_P1( locTangentStab_P1, absdet);
    localCahnHilliard_.setupVolumeStab_P1( locVolumeStab_P1, absdet);


    localCahnHilliard_.setupL_P1( locL_P1P1);
    localCahnHilliard_.setupLM_P1( locLM_P1P1);

    localCahnHilliard_.setupGprimeprime_P1( locGprimeprime_P1P1);


}

void CahnHilliardIFAccumulator_P1P1CL::update_global_system ()
{
    MatrixBuilderCL& mM_P1= *mM_P1_;
    MatrixBuilderCL& mNormalStab_P1= *mNormalStab_P1_;
    MatrixBuilderCL& mTangentStab_P1= *mTangentStab_P1_;
    MatrixBuilderCL& mVolumeStab_P1= *mVolumeStab_P1_;

    MatrixBuilderCL& mL_P1P1= *mL_P1P1_;
    MatrixBuilderCL& mLM_P1P1= *mLM_P1P1_;
    MatrixBuilderCL& mGprimeprime_P1P1= *mGprimeprime_P1P1_;


    int count =0;

    for(int i= 0; i < 4; ++i) {
        const IdxT ii= nScalar.num[i];
        if (ii==NoIdx || !(nScalar.WithUnknowns( i))) continue;
        for(int j=0; j <4; ++j) {
            const IdxT jj= nScalar.num[j];
            if (jj==NoIdx || !(nScalar.WithUnknowns( j))) continue;
            mM_P1( ii, jj) += locM_P1[i][j];
            mNormalStab_P1( ii, jj) += locNormalStab_P1[i][j];
            mTangentStab_P1( ii, jj) += locTangentStab_P1[i][j];
            mVolumeStab_P1( ii, jj) += locVolumeStab_P1[i][j];
            mL_P1P1( ii, jj) += locL_P1P1[i][j];
            mLM_P1P1( ii, jj) += locLM_P1P1[i][j];
            mGprimeprime_P1P1( ii, jj) += locGprimeprime_P1P1[j][i];

        }
    }
}

void SetupCahnHilliardIF_P1P1( const MultiGridCL& MG_,  MatDescCL* M_P1, MatDescCL* NormalStab_P1, MatDescCL* TangentStab_P1,MatDescCL* VolumeStab_P1, MatDescCL* L_P1P1 , MatDescCL* LM_P1P1,  MatDescCL* Gprimeprime_P1P1, const VecDescCL& lset, const LsetBndDataCL& lset_bnd, const VecDescCL& velocity,  const BndDataCL<Point3DCL>& velocity_bnd, const VecDescCL& chi,const BndDataCL<>& chi_bnd)
{
  ScopeTimerCL scope("SetupCahnHilliardIF_P1P1");
  CahnHilliardIFAccumulator_P1P1CL accu( lset, velocity, chi, lset_bnd, velocity_bnd, chi_bnd, *(M_P1->RowIdx), *(L_P1P1->RowIdx),  M_P1->Data, NormalStab_P1->Data, TangentStab_P1->Data, VolumeStab_P1->Data,  L_P1P1->Data,  LM_P1P1->Data,  Gprimeprime_P1P1->Data);
  TetraAccumulatorTupleCL accus;
  //    MaybeAddProgressBar(MG_, "LapBeltr(P2) Setup", accus, RowIdx.TriangLevel());
  accus.push_back( &accu);
  accumulate( accus, MG_, M_P1->GetRowLevel(), /*M_P1->RowIdx->GetMatchingFunction(),*/ M_P1->RowIdx->GetBndInfo());

}

} // end of namespace DROPS