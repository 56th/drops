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

void SetupInterfaceMassP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd)
{
    //ScopeTimerCL timer( "SetupInterfaceMassP1");

    TetraAccumulatorTupleCL accus;
    InterfaceCommonDataP1CL cdata( ls, lsetbnd);
    accus.push_back( &cdata);
    InterfaceMatrixAccuP1CL<LocalInterfaceMassP1CL> accu( mat, LocalInterfaceMassP1CL(), cdata);
    accus.push_back( &accu);
    const IdxDescCL* RowIdx= mat->RowIdx;
    accumulate( accus, mg, RowIdx->TriangLevel(), RowIdx->GetMatchingFunction(), RowIdx->GetBndInfo());

    // WriteToFile( mat->Data, "mass.txt", "mass");
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
    InterfaceMatrixAccuP1CL<LocalLaplaceBeltramiP1CL> accu( mat, LocalLaplaceBeltramiP1CL( D), cdata);
    accus.push_back( &accu);
    const IdxDescCL* RowIdx= mat->RowIdx;
    accumulate( accus, mg, RowIdx->TriangLevel(), RowIdx->GetMatchingFunction(), RowIdx->GetBndInfo());

    // WriteToFile( mat->Data, "lb.txt", "lb");
}

void SetupInterfaceRhsP1 (const MultiGridCL& mg, VecDescCL* v,
    const VecDescCL& ls, const BndDataCL<>& lsetbnd, instat_scalar_fun_ptr f)
{
    std::cout << "entering SetupInterfaceRhsP1: " << v->RowIdx->NumUnknowns() << " dof... ";

    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( 2);
    LocalP2CL<> locp2_ls;
    std::valarray<double> ls_loc( lat.vertex_size());
    SurfacePatchCL surf;
    QuadDomain2DCL qdom;
    IdxT num[4];

    LocalP1CL<> p1[4];
    p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; // P1-Basis-Functions
    std::valarray<double> qp1,
                          qf;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, v->GetLevel(), it) {
        locp2_ls.assign( *it, ls, lsetbnd);
        evaluate_on_vertexes( locp2_ls, lat, Addr( ls_loc));
        if (equal_signs( ls_loc))
            continue;

        GetLocalNumbP1NoBnd( num, *it, *v->RowIdx);
        surf.make_patch<MergeCutPolicyCL>( lat, ls_loc);
        make_CompositeQuad5Domain2D ( qdom, surf, *it);
        resize_and_evaluate_on_vertexes( f, *it, qdom, v->t, qf);
        qp1.resize( qdom.vertex_size());
        for (Uint i= 0; i < 4; ++i)
            if (num[i] != NoIdx) {
                evaluate_on_vertexes( p1[i], qdom, Addr( qp1));
                v->Data[num[i]]+= quad_2D( qf*qp1, qdom);
            }
    }
    std::cout << " Rhs set up." << std::endl;

    // WriteToFile( v->Data, "rhs.txt", "Rhs");
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


/// ==Space-Time-methods==
/// \todo Maybe introduce a new file
///\brief Bilinear space-time FE-function on a single TetraPrismCL.

LocalSTP1P1CL<> STP1P1DiscCL::ref_val[8];
LocalSTP1P1CL<Point4DCL> STP1P1DiscCL::ref_grad[8];
LocalSTP1P1CL<Point3DCL> STP1P1DiscCL::ref_gradx[8];
// LocalSTP1P1CL< SMatrixCL<4,4> > STP1P1DiscCL::ref_Hess[8];

void STP1P1DiscCL::StaticInit()
{
    for (Uint i= 0; i < 4; ++i) {
        ref_val[i    ].at_t0()[i]= 1.;
        ref_val[i + 4].at_t1()[i]= 1.;

        const Point3DCL p1grad( FE_P1CL::DHRef( i));
        for (Uint j= 0; j < 4; ++j) {
            ref_grad[i    ].at_t0()[j]= MakePoint4D( p1grad[0], p1grad[1], p1grad[2], j == i ? -1. : 0.);
            ref_grad[i + 4].at_t1()[j]= MakePoint4D( p1grad[0], p1grad[1], p1grad[2], j == i ?  1. : 0.);
            ref_gradx[i    ].at_t0()[j]= p1grad;
            ref_gradx[i + 4].at_t1()[j]= p1grad;
        }
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
    LocalSTP2P1ProxyCL<> local_st_lset( oldlocp2_ls, locp2_ls);
    SPatchCL<4> patch;
    QuadDomainCodim1CL<4> qdom;
    std::valarray<double> shape; // STP1P1-shape-function as integrand
    DROPS_FOR_TRIANG_TETRA ( mg, level, it) {
        oldlocp2_ls.assign( *it, oldls, lsetbnd);
        locp2_ls.assign(    *it, newls, lsetbnd);
        evaluate_on_vertexes( local_st_lset, lat, Addr( ls_loc));
        if (equal_signs( ls_loc))
            continue;

        patch.make_patch<MergeCutPolicyCL>( lat, ls_loc);
        make_CompositeQuad5DomainSTCodim1( qdom, patch, TetraPrismCL( *it, t0, t1));
        for (Uint i= 0; i < 8; ++i) {
            const IdxDescCL& idx= i < 4 ? idx0_ : idx1_;
            const IdxDescCL& spatial_idx= i < 4 ? idx_ini : idx_fini;
            UnknownHandleCL& unknowns= const_cast<VertexCL*>( it->GetVertex( i%4))->Unknowns;
            if (unknowns.Exist( idx.GetIdx()))
                continue;

            resize_and_evaluate_on_vertexes( STP1P1DiscCL::ref_val[i], qdom, shape);
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


VectorCL SurfactantcGdGP1CL::InitStep (double new_t)
{
    std::cout << "SurfactantcGdGP1CL::InitStep:\n";
    ic.t= new_t;
    dt_= ic.t - oldt_;

    st_idx_.CreateNumbering( oldidx_.TriangLevel(), MG_, oldls_, lset_vd_, lsetbnd_, oldt_, new_t);
    std::cout << "space-time Unknowns: " << st_idx_.NumUnknowns()
              << " ini: " << st_idx_.NumIniUnknowns() << " fini: " << st_idx_.NumFiniUnknowns()
              << " interior: " << st_idx_.NumUnknowns() - st_idx_.NumIniUnknowns() - st_idx_.NumFiniUnknowns() << std::endl;
    st_ic_.resize( 0);
    st_ic_.resize( st_idx_.NumUnknowns());

    st_oldic_.resize( 0);
    st_oldic_.resize( st_idx_.NumUnknowns());
    // Copy dofs on the old interface from  old solution into st_oldic_.
    DROPS_FOR_TRIANG_VERTEX( MG_, oldidx_.TriangLevel(), it) {
        if (it->Unknowns.Exist( st_idx_.GetIdx( 0))) {
            const IdxT dof= it->Unknowns( st_idx_.GetIdx( 0));
            if (dof < st_idx_.NumIniUnknowns())
                st_oldic_[dof]= oldic_[dof];
        }
    }


    return VectorCL( st_idx_.NumUnknowns());
}

void SurfactantcGdGP1CL::Update()
{
    std::cout << "SurfactantcGP1CL::Update:\n";

    std::cout << "SurfactantcGdGP1CL::Update: Finished\n";
}

void SurfactantcGdGP1CL::DoStep (const VectorCL& rhs)
{
    Update();

    // L_.LinComb( 1., Mder, 1., Mdiv, 1., A, 1., Mold);
    std::cout << "Before solve: res = " << norm( L_*st_ic_ - rhs) << std::endl;
    //gm_.Solve( L_, ic_st_, rhs);
    std::cout << "res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
}

void SurfactantcGdGP1CL::CommitStep ()
{
    idx.CreateNumbering( oldidx_.TriangLevel(), MG_, &lset_vd_, &lsetbnd_);
    std::cout << "new NumUnknowns at t1: " << idx.NumUnknowns() << std::endl;

    ic.SetIdx( &idx);
    // Copy dofs on the new interface from space-time-solution into ic.
    DROPS_FOR_TRIANG_VERTEX( MG_, oldidx_.TriangLevel(), it) {
        if (it->Unknowns.Exist( st_idx_.GetIdx( 1))) {
            const IdxT dof= it->Unknowns( st_idx_.GetIdx( 1));
            if (dof >= st_idx_.NumIniUnknowns() && dof < st_idx_.NumIniUnknowns() + st_idx_.NumFiniUnknowns())
                ic.Data[dof - st_idx_.NumIniUnknowns()]= st_ic_[dof];
        }
    }

    L_.clear();
    st_ic_.resize( 0);
    st_oldic_.resize( 0);
    st_idx_.DeleteNumbering( MG_);
}

void SurfactantcGdGP1CL::DoStep (double new_t)
{
    VectorCL rhs( InitStep( new_t));
    DoStep( rhs);
    CommitStep();
}

} // end of namespace DROPS
