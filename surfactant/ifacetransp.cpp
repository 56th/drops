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

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
        if (it->Unknowns.Exist( xidx) && it->Unknowns.Exist( xextidx))
            xext.Data[it->Unknowns( xextidx)]= x.Data[it->Unknowns( xidx)];
    }
}

void Restrict (const MultiGridCL& mg, const VecDescCL& xext, VecDescCL& x)
{
    const Uint xidx( x.RowIdx->GetIdx()),
        xextidx( xext.RowIdx->GetIdx()),
        lvl( x.RowIdx->TriangLevel());

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
        if (it->Unknowns.Exist( xidx) && it->Unknowns.Exist( xextidx))
            x.Data[it->Unknowns( xidx)]= xext.Data[it->Unknowns( xextidx)];
    }
}

void update_global_matrix_P1 (MatrixBuilderCL& M, const double coup[4][4], const IdxT numr[4], const IdxT numc[4])
{
    for (int i= 0; i < 4; ++i)
        if (numr[i] != NoIdx)
            for (int j= 0; j < 4; ++j)
                if (numc[j] != NoIdx)
                    M( numr[i], numc[j])+= coup[i][j];
}

void InterfaceMassAccuP1CL::setup_local_matrix (const TetraCL& t)
{
    surf.make_patch<MergeCutPolicyCL>( lat, ls_loc);
    make_CompositeQuad2Domain2D ( qdom, surf, t);

    for (int i= 0; i < 4; ++i)
        resize_and_evaluate_on_vertexes (p1[i], qdom, q[i]);

    for (int i= 0; i < 4; ++i) {
        coup[i][i]= quad_2D( q[i]*q[i], qdom);
        for(int j= 0; j < i; ++j)
            coup[i][j]= coup[j][i]= quad_2D( q[j]*q[i], qdom);
    }
}

void InterfaceMassAccuP1CL::begin_accumulation ()
{
    const IdxT num_rows= mat_->RowIdx->NumUnknowns();
    const IdxT num_cols= mat_->ColIdx->NumUnknowns();
    std::cout << "entering InterfaceMassAccuP1CL::begin_accumulation: " << num_rows << " rows, " << num_cols << " cols, ";
    lvl = mat_->GetRowLevel();
    M= new MatrixBuilderCL( &mat_->Data, num_rows, num_cols);
}

void InterfaceMassAccuP1CL::finalize_accumulation ()
{
    M->Build();
    delete M;
    M= 0;
    std::cout << mat_->Data.num_nonzeros() << " nonzeros." << std::endl;
}

void InterfaceMassAccuP1CL::visit (const TetraCL& t)
{
    locp2_ls.assign( t, ls_, lsetbnd_);
    evaluate_on_vertexes( locp2_ls, lat, Addr( ls_loc));
    if (equal_signs( ls_loc))
        return;

    setup_local_matrix ( t);

    GetLocalNumbP1NoBnd( numr, t, *mat_->RowIdx);
    GetLocalNumbP1NoBnd( numc, t, *mat_->ColIdx);
    update_global_matrix_P1( *M, coup, numr, numc);
}

void SetupInterfaceMassP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd)
{
    ScopeTimerCL timer( "SetupInterfaceMassP1");

    InterfaceMassAccuP1CL accu( mat, ls, lsetbnd);
    TetraAccumulatorTupleCL accus;
    accus.push_back( &accu);
    const IdxDescCL* RowIdx= mat->RowIdx;
    accumulate( accus, mg, RowIdx->TriangLevel(), RowIdx->GetMatchingFunction(), RowIdx->GetBndInfo());

}

void SetupMixedMassP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd)
{
    SetupInterfaceMassP1( mg, mat, ls, lsetbnd);
}

void LaplaceBeltramiAccuP1CL::setup_local_matrix (const TetraCL& t)
{
    surf.make_patch<MergeCutPolicyCL>( lat, ls_loc);
    resize_and_evaluate_piecewise_normal( surf, t, n, &absdet);
    P1DiscCL::GetGradients( grad, dummy, t);
    for(int i= 0; i < 4; ++i) {
        q[i].resize( surf.triangle_size());
        q[i]= grad[i] - dot( grad[i], n)*n;
    }

    for (int i= 0; i < 4; ++i) {
        coup[i][i]= /*area of reference triangle*/ 0.5*(dot(q[i], q[i])*absdet).sum();
        for(int j= 0; j < i; ++j)
            coup[i][j]= coup[j][i]= /*area of reference triangle*/ 0.5*(dot(q[i], q[j])*absdet).sum();
    }
}

void LaplaceBeltramiAccuP1CL::begin_accumulation ()
{
    const IdxT num_rows= mat_->RowIdx->NumUnknowns();
    const IdxT num_cols= mat_->ColIdx->NumUnknowns();
    std::cout << "entering LaplaceBeltramiAccuP1CL::begin_accumulation: " << num_rows << " rows, " << num_cols << " cols, ";
    lvl = mat_->GetRowLevel();
    A= new MatrixBuilderCL( &mat_->Data, num_rows, num_cols);
}

void LaplaceBeltramiAccuP1CL::finalize_accumulation ()
{
    A->Build();
    delete A;
    A= 0;
    mat_->Data*= D_; // diffusion coefficient
    std::cout << mat_->Data.num_nonzeros() << " nonzeros." << std::endl;
}

void LaplaceBeltramiAccuP1CL::visit (const TetraCL& t)
{
    locp2_ls.assign( t, ls_, lsetbnd_);
    if (equal_signs( locp2_ls))
        return;

    evaluate_on_vertexes( locp2_ls, lat, Addr( ls_loc));

    setup_local_matrix( t);

    GetLocalNumbP1NoBnd( numr, t, *mat_->RowIdx);
    GetLocalNumbP1NoBnd( numc, t, *mat_->ColIdx);
    update_global_matrix_P1( *A, coup, numr, numc);
}

void SetupLBP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd, double D)
{
    ScopeTimerCL timer( "SetuLBP1");

    LaplaceBeltramiAccuP1CL accu( mat, ls, lsetbnd, D);
    TetraAccumulatorTupleCL accus;
    accus.push_back( &accu);
    const IdxDescCL* RowIdx= mat->RowIdx;
    accumulate( accus, mg, RowIdx->TriangLevel(), RowIdx->GetMatchingFunction(), RowIdx->GetBndInfo());
}

void SetupInterfaceRhsP1OnTriangle (const LocalP1CL<> p1[4],
    Quad5_2DCL<> q[4],VectorCL& v, const IdxT Numb[4],
    const TetraCL& t, const BaryCoordCL triangle[3], double det,
    instat_scalar_fun_ptr f, double time)
{
    for (int i= 0; i < 4; ++i)
        q[i].assign( p1[i], triangle);
    Quad5_2DCL<> qf( t, triangle, f, time), r;

    for (int i= 0; i < 4; ++i) {
        if (Numb[i] == NoIdx) continue;
        r= qf*q[i];
        v[Numb[i]]+= r.quad( det);
    }
}

void SetupInterfaceRhsP1 (const MultiGridCL& mg, VecDescCL* v,
    const VecDescCL& ls, const BndDataCL<>& lsetbnd, instat_scalar_fun_ptr f)
{
    const IdxT num_unks= v->RowIdx->NumUnknowns();
    const Uint lvl = v->GetLevel();
    const double t= v->t;

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
                        *it, &triangle.GetBary( tri), triangle.GetAbsDet( tri), f, t);
            }
        }
    }
    std::cout << " Rhs set up." << std::endl;
}

/// \todo This should be a generic function somewhere in num or misc.
void P1Init (instat_scalar_fun_ptr icf, VecDescCL& ic, MultiGridCL& mg, double t)
{
    const Uint lvl= ic.GetLevel(),
               idx= ic.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_VERTEX( mg, lvl, it) {
        if (it->Unknowns.Exist( idx))
            ic.Data[it->Unknowns( idx)]= icf( it->GetCoord(), t);
    }
    ic.t= t;
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

    oldv_.SetIdx( v_->RowIdx);
    oldv_.Data= v_->Data;
}

void SurfactantcGP1CL::Update()
{
    // std::cout << "SurfactantcGP1CL::Update:\n";
    IdxDescCL* cidx= ic.RowIdx;

    M.Data.clear();
    M.SetIdx( cidx, cidx);
    DROPS::SetupInterfaceMassP1( MG_, &M, lset_vd_, lsetbnd_);
    // std::cout << "M is set up.\n";
    A.Data.clear();
    A.SetIdx( cidx, cidx);
    DROPS::SetupLBP1( MG_, &A, lset_vd_, lsetbnd_, D_);
    // std::cout << "A is set up.\n";
    C.Data.clear();
    C.SetIdx( cidx, cidx);
    DROPS::SetupConvectionP1( MG_, &C, lset_vd_, lsetbnd_, make_P2Eval( MG_, Bnd_v_, *v_));
    // std::cout << "C is set up.\n";
    Md.Data.clear();
    Md.SetIdx( cidx, cidx);
    DROPS::SetupMassDivP1( MG_, &Md, lset_vd_, lsetbnd_, make_P2Eval( MG_, Bnd_v_, *v_));
    // std::cout << "Md is set up.\n";

    if (theta_ != 1.0) {
        M2.Data.clear();
        M2.SetIdx( cidx, cidx);
        DROPS::SetupInterfaceMassP1( MG_, &M2, oldls_, lsetbnd_);
        // std::cout << "M2 is set up.\n";
    }
    std::cout << "SurfactantP1CL::Update: Finished\n";
}

VectorCL SurfactantcGP1CL::InitStep (double new_t)
{

    // std::cout << "SurfactantcGP1CL::InitStep:\n";
    ic.t= new_t;
    dt_= ic.t - oldt_;
    idx.CreateNumbering( oldidx_.TriangLevel(), MG_, &lset_vd_, &lsetbnd_); // InitOld deletes oldidx_ and swaps idx and oldidx_.
    std::cout << "new NumUnknowns: " << idx.NumUnknowns() << std::endl;
    ic.SetIdx( &idx);

    MatDescCL m( &idx, &oldidx_);
    DROPS::SetupMixedMassP1( MG_, &m, lset_vd_, lsetbnd_);
    // std::cout << "mixed M on new interface is set up.\n";
    VectorCL rhs( theta_*(m.Data*oldic_));

    if (rhs_fun_) {
        VecDescCL load( &idx);
        load.t= new_t;
        DROPS::SetupInterfaceRhsP1 (MG_, &load, lset_vd_, lsetbnd_, rhs_fun_);
        rhs+= theta_*dt_*load.Data;
    }

    if (theta_ == 1.0) return rhs;

    m.Data.clear();
    DROPS::SetupMixedMassP1( MG_, &m, oldls_, lsetbnd_);
    // std::cout << "mixed M on old interface is set up.\n";
    rhs+= (1. - theta_)*(m.Data*oldic_);

    m.Data.clear();
    DROPS::SetupLBP1( MG_, &m, oldls_, lsetbnd_, D_);
    // std::cout << "mixed A on old interface is set up.\n";
    VectorCL rhs2( m.Data*oldic_);
    m.Data.clear();
    DROPS::SetupConvectionP1( MG_, &m, oldls_, lsetbnd_, make_P2Eval( MG_, Bnd_v_, oldv_));
    // std::cout << "mixed C on old interface is set up.\n";
    rhs2+= m.Data*oldic_;
    m.Data.clear();
    DROPS::SetupMassDivP1( MG_, &m, oldls_, lsetbnd_, make_P2Eval( MG_, Bnd_v_, oldv_));
    // std::cout << "mixed Md on old interface is set up.\n";
    rhs2+= m.Data*oldic_;
    if (rhs_fun_) {
        VecDescCL load( &idx);
        load.t= new_t - dt_;
        DROPS::SetupInterfaceRhsP1( MG_, &load, oldls_, lsetbnd_, rhs_fun_);
        rhs+= (1. - theta_)*dt_*load.Data;
    }

    return VectorCL( rhs - ((1. - theta_)*dt_)*rhs2);
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
    gm_.Solve( L_, ic.Data, rhs);
    std::cout << "res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
}

void SurfactantcGP1CL::CommitStep ()
{
    return;
}

void SurfactantcGP1CL::DoStep (double new_t)
{
    VectorCL rhs( InitStep( new_t));
    DoStep( rhs);
    CommitStep();
}

void
InterfaceP1RepairCL::pre_refine ()
{
    DROPS::NoBndDataCL<> dummy;
    p1repair_= std::auto_ptr<RepairP1CL<double, NoBndDataCL>::type >(
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
    IdxDescCL p1idx;
    p1idx.CreateNumbering( u_.RowIdx->TriangLevel(), mg_);
    VecDescCL p1u( &p1idx);
    Extend( mg_, u_, p1u);
    BndDataCL<> bnd( 0);

    cf.PutScalar( make_P1Eval( mg_, bnd, p1u), varName());

    p1idx.DeleteNumbering( mg_);
}

} // end of namespace DROPS
