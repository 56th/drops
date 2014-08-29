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
    if (x.RowIdx->GetFE() == P1IF_FE)
        return;

    // For P2IF_FE, also fixup the edge-dofs.
    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it) {
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
    : ls( &ls_arg), lsetbnd( &lsetbnd_arg), lat( &lat_arg), compute_absdet_( true), ls_loc( lat->vertex_size()), quaqua( quaquaarg)
{
    P2DiscCL::GetGradientsOnRef( gradrefp2);
    for (Uint i= 0; i < 10 ; ++i)
        p2[i][i]= 1.; // P2-Basis-Functions
}

void SetupInterfaceMassP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd)
{
    //ScopeTimerCL timer( "SetupInterfaceMassP1");

    TetraAccumulatorTupleCL accus;
    InterfaceCommonDataP1CL cdata( ls, lsetbnd);
    accus.push_back( &cdata);
    InterfaceMatrixAccuCL<LocalInterfaceMassP1CL, InterfaceCommonDataP1CL> accu( mat, LocalInterfaceMassP1CL(), cdata);
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
    InterfaceMatrixAccuCL<LocalLaplaceBeltramiP1CL, InterfaceCommonDataP1CL> accu( mat, LocalLaplaceBeltramiP1CL( D), cdata);
    accus.push_back( &accu);
    const IdxDescCL* RowIdx= mat->RowIdx;
    accumulate( accus, mg, RowIdx->TriangLevel(), RowIdx->GetMatchingFunction(), RowIdx->GetBndInfo());

    // WriteToFile( mat->Data, "lb.txt", "lb");
}

void SetupInterfaceRhsP1 (const MultiGridCL& mg, VecDescCL* v,
    const VecDescCL& ls, const BndDataCL<>& lsetbnd, instat_scalar_fun_ptr f)
{
    TetraAccumulatorTupleCL accus;
    InterfaceCommonDataP1CL cdata( ls, lsetbnd);
    accus.push_back( &cdata);
    InterfaceVectorAccuCL<LocalVectorP1CL, InterfaceCommonDataP1CL> loadaccu( v, LocalVectorP1CL( f, v->t), cdata);
    accus.push_back( &loadaccu);
    accumulate( accus, mg, v->RowIdx->TriangLevel(), v->RowIdx->GetMatchingFunction(), v->RowIdx->GetBndInfo());

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

    const TetraCL* tet;
    BaryCoordCL b;

    if (to_iface != 0) {
        const Uint sys= to_iface->RowIdx->GetIdx();
        for (Uint i= 0; i < 4; ++i) {
            tet= &t;
            b= std_basis<4>( i + 1);
            cdata.quaqua.base_point( tet, b);
            Point3DCL offset= t.GetVertex( i)->GetCoord() - GetWorldCoord( *tet, b);
            const size_t dof= t.GetVertex( i)->Unknowns( sys);
            std::copy( Addr( offset), Addr( offset) + 3, &to_iface->Data[dof]);
        }
    }

    for (SurfacePatchCL::const_vertex_iterator it= cdata.surf.vertex_begin(); it != cdata.surf.vertex_end(); ++it) {
        tet= &t;
        b= *it;
        cdata.quaqua.base_point( tet, b);
        const Point3DCL& x= GetWorldCoord( t, *it);
//         const Point3DCL& xb= GetWorldCoord( *tet, b);
//         std::cout  << "    |x-xb|: " << (x - xb).norm();

        SMatrixCL<3,3> dph;
        cdata.quaqua.jacobian( t, *it, *tet, b, dph);
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
    surfacemeasP2+= quad_2D( cdata.absdet, cdata.qdom);
}

void gradient_trafo (const TetraCL& tet, const BaryCoordCL& xb, const QuaQuaMapperCL& quaqua, const SurfacePatchCL& p,
                     SMatrixCL<3,3>& W)
{
    // Compute the basepoint b.
    const TetraCL* btet= &tet;
    BaryCoordCL b= xb;
    quaqua.base_point( btet, b);

    // nl(x)
    if (p.normal_empty())
        p.compute_normals( tet);
    Point3DCL nl= p.normal_begin()[0];
    // Evaluate the normal to the interface in b, n(y).
    Point3DCL n= quaqua.local_ls_grad( *btet, b);
    n/= n.norm();

    // Dp_h(x)^T
    SMatrixCL<3,3> dph, dphT;
    quaqua.jacobian( tet, xb, *btet, b, dph);
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
    accumulate( accus, MG_, cidx->TriangLevel(), cidx->GetMatchingFunction(), cidx->GetBndInfo());

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
        accumulate( accus, MG_, idx.TriangLevel(), idx.GetMatchingFunction(), idx.GetBndInfo());
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

    accumulate( accus, MG_, idx.TriangLevel(), idx.GetMatchingFunction(), idx.GetBndInfo());
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
    IdxDescCL fullidx( P1_FE);
    if (u_.RowIdx->GetFE() == P2IF_FE)
        fullidx.SetFE( P2_FE);
    fullidx.CreateNumbering( u_.RowIdx->TriangLevel(), mg_);
    VecDescCL uext( &fullidx);
    Extend( mg_, u_, uext);
    BndDataCL<> bnd( 0);

    if (u_.RowIdx->GetFE() == P1IF_FE)
        cf.PutScalar( make_P1Eval( mg_, bnd, uext), varName());
    else if (u_.RowIdx->GetFE() == P2IF_FE)
        cf.PutScalar( make_P2Eval( mg_, bnd, uext), varName());

    fullidx.DeleteNumbering( mg_);
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

    accumulate( accus, MG_, st_idx_.TriangLevel(), idx.GetMatchingFunction(), idx.GetBndInfo());

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

    accumulate( accus, MG_, st_idx_.TriangLevel(), idx.GetMatchingFunction(), idx.GetBndInfo());
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

    accumulate( accus, MG_, cidx->TriangLevel(), cidx->GetMatchingFunction(), cidx->GetBndInfo());

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

} // end of namespace DROPS
