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

#include "levelset/levelset.h"
#include "misc/scopetimer.h"
#include <cstring>

namespace DROPS {

template <class ResultContainerT>
  void
  resize_and_evaluate_piecewise_normal (const SurfacePatchCL& p, const TetraCL& t, ResultContainerT& n, std::valarray<double>* absdet)
{
    n.resize( p.facet_size());
    if (absdet)
        absdet->resize( p.facet_size());

    const typename SurfacePatchCL::const_vertex_iterator verts= p.vertex_begin();
    typename SequenceTraitCL<ResultContainerT>::iterator n_it= sequence_begin( n);
    double* a_it= absdet ? Addr( *absdet) : 0;

    Point3DCL tmp( Uninitialized);
    double tmp_norm;
    for (SurfacePatchCL::const_facet_iterator it= p.facet_begin(); it != p.facet_end(); ++it) {
        const Point3DCL& v0= GetWorldCoord( t, verts[(*it)[0]]);
        cross_product( tmp, GetWorldCoord( t, verts[(*it)[1]]) - v0,
                            GetWorldCoord( t, verts[(*it)[2]]) - v0);
        tmp_norm= tmp.norm();
        *n_it++= tmp/tmp_norm;
        if (absdet)
            *a_it++= tmp_norm;
    }
}

template <class DiscVelSolT>
void InterfaceConvectionAccuP1CL<DiscVelSolT>::setup_local_matrix (const TetraCL& t)
{
    surf.make_patch<MergeCutPolicyCL>( lat, ls_loc);
    make_CompositeQuad5Domain2D( qdom, surf, t);

    w_loc.assign( t, w_);
    resize_and_evaluate_on_vertexes( w_loc, qdom, qw);

    P1DiscCL::GetGradients( grad, dummy, t);
    for (int i= 0; i < 4; ++i)
        resize_and_evaluate_on_vertexes( p1[i], qdom, q[i]);

    for (int i= 0; i < 4; ++i)
        for(int j= 0; j < 4; ++j)
            coup[i][j]= quad_2D( dot( grad[j], qw)*q[i], qdom);
}

template <class DiscVelSolT>
void InterfaceConvectionAccuP1CL<DiscVelSolT>::begin_accumulation ()
{
    const IdxT num_rows= mat_->RowIdx->NumUnknowns();
    const IdxT num_cols= mat_->ColIdx->NumUnknowns();
    std::cout << "entering InterfaceConvectionAccuP1CL::begin_accumulation: " << num_rows << " rows, " << num_cols << " cols, ";
    lvl = mat_->GetRowLevel();
    M= new MatrixBuilderCL( &mat_->Data, num_rows, num_cols);
}

template <class DiscVelSolT>
void InterfaceConvectionAccuP1CL<DiscVelSolT>::finalize_accumulation ()
{
    M->Build();
    delete M;
    M= 0;
    std::cout << mat_->Data.num_nonzeros() << " nonzeros." << std::endl;
}

template <class DiscVelSolT>
void InterfaceConvectionAccuP1CL<DiscVelSolT>::visit (const TetraCL& t)
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

template <class DiscVelSolT>
void SetupConvectionP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd, const DiscVelSolT& w)
{
    //ScopeTimerCL timer( "SetupConvectionP1");

    InterfaceConvectionAccuP1CL<DiscVelSolT> accu( mat, ls, lsetbnd, w);
    TetraAccumulatorTupleCL accus;
    accus.push_back( &accu);
    const IdxDescCL* RowIdx= mat->RowIdx;
    accumulate( accus, mg, RowIdx->TriangLevel(), RowIdx->GetMatchingFunction(), RowIdx->GetBndInfo());
}


template <class DiscVelSolT>
void InterfaceMassDivAccuP1CL<DiscVelSolT>::setup_local_matrix (const TetraCL& t)
{
    surf.make_patch<MergeCutPolicyCL>( lat, ls_loc);
    make_CompositeQuad5Domain2D ( qdom, surf, t);

    resize_and_evaluate_piecewise_normal( surf, t, n_tri);
    n.resize( qdom.vertex_size());
    for (Uint i= 0; i < surf.facet_size(); ++i)
        n[std::slice(i*7, 7, 1)]= n_tri[i];

    GetTrafoTr( T, dummy, t);
    P2DiscCL::GetGradients( gradp2, gradrefp2, T);
    w_loc.assign( t, w_);
    qgradp2i.resize( qdom.vertex_size());
    qdivgamma_w.resize( qdom.vertex_size());
    qdivgamma_w= 0.;
    for (int i= 0; i < 10; ++i) {
        evaluate_on_vertexes( gradp2[i], qdom, Addr( qgradp2i));
        qdivgamma_w+= dot(w_loc[i], qgradp2i) - dot( w_loc[i], n)*dot( n, qgradp2i);
    }

    for (int i= 0; i < 4; ++i)
        resize_and_evaluate_on_vertexes (p1[i], qdom, q[i]);

    for (int i= 0; i < 4; ++i) {
        coup[i][i]= quad_2D( qdivgamma_w*q[i]*q[i], qdom);
        for(int j= 0; j < i; ++j)
            coup[i][j]= coup[j][i]= quad_2D( qdivgamma_w*q[j]*q[i], qdom);
    }
}

template <class DiscVelSolT>
void InterfaceMassDivAccuP1CL<DiscVelSolT>::begin_accumulation ()
{
    const IdxT num_rows= mat_->RowIdx->NumUnknowns();
    const IdxT num_cols= mat_->ColIdx->NumUnknowns();
    std::cout << "entering InterfaceMassDivAccuP1CL::begin_accumulation: " << num_rows << " rows, " << num_cols << " cols, ";
    lvl = mat_->GetRowLevel();
    M= new MatrixBuilderCL( &mat_->Data, num_rows, num_cols);
}

template <class DiscVelSolT>
void InterfaceMassDivAccuP1CL<DiscVelSolT>::finalize_accumulation ()
{
    M->Build();
    delete M;
    M= 0;
    std::cout << mat_->Data.num_nonzeros() << " nonzeros." << std::endl;
}

template <class DiscVelSolT>
void InterfaceMassDivAccuP1CL<DiscVelSolT>::visit (const TetraCL& t)
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

template <class DiscVelSolT>
void SetupMassDivP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd, const DiscVelSolT& w)
{
    //ScopeTimerCL timer( "SetupMassDivP1");

    InterfaceMassDivAccuP1CL<DiscVelSolT> accu( mat, ls, lsetbnd, w);
    TetraAccumulatorTupleCL accus;
    accus.push_back( &accu);
    const IdxDescCL* RowIdx= mat->RowIdx;
    accumulate( accus, mg, RowIdx->TriangLevel(), RowIdx->GetMatchingFunction(), RowIdx->GetBndInfo());
}

} // end of namespace DROPS
