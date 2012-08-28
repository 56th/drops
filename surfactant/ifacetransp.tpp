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

template <class ResultContainerT>
  void
  resize_and_evaluate_piecewise_normal (const SPatchCL<4>& p, const TetraPrismCL& prism, ResultContainerT& n, std::valarray<double>* absdet)
{
    n.resize( p.facet_size());
    if (absdet)
        absdet->resize( p.facet_size());

    const typename SPatchCL<4>::const_vertex_iterator verts= p.vertex_begin();
    typename SequenceTraitCL<ResultContainerT>::iterator n_it= sequence_begin( n);
    double* a_it= absdet ? Addr( *absdet) : 0;

    QRDecompCL<4,3> qr;
    SMatrixCL<4,3>& M= qr.GetMatrix();
    Point4DCL tmp( Uninitialized);
    double tmp_norm;
    for (SPatchCL<4>::const_facet_iterator it= p.facet_begin(); it != p.facet_end(); ++it) {
        const Point4DCL& v0= GetWorldCoord( prism, verts[(*it)[0]]);
        for (Uint i= 1; i < 4; ++i)
            M.col( i - 1, GetWorldCoord( prism, verts[(*it)[i]]) - v0);
        const bool is_rank_deficient= qr.prepare_solve( /*assume_full_rank*/ false);
        tmp= 0.;
        if (is_rank_deficient) {
            *n_it++= tmp;
            if (absdet)
                a_it++= 0.;
        }
        else {
            tmp[3]= 1.;
            qr.apply_Q( tmp);
            tmp_norm= tmp.norm();
            *n_it++= tmp/tmp_norm;
            if (absdet)
                *a_it++= tmp_norm;
        }
    }
}

template <class DiscVelSolT>
void LocalInterfaceConvectionP1CL<DiscVelSolT>::setup (const TetraCL& t, const InterfaceCommonDataP1CL& cdata)
{
    make_CompositeQuad5Domain2D( qdom, cdata.surf, t);

    w_loc.assign( t, w_);
    resize_and_evaluate_on_vertexes( w_loc, qdom, qw);

    P1DiscCL::GetGradients( grad, dummy, t);
    for (int i= 0; i < 4; ++i)
        resize_and_evaluate_on_vertexes( cdata.p1[i], qdom, q[i]);

    for (int i= 0; i < 4; ++i)
        for(int j= 0; j < 4; ++j)
            coup[i][j]= quad_2D( dot( grad[j], qw)*q[i], qdom);
}


template <class DiscVelSolT>
void SetupConvectionP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd, const DiscVelSolT& w)
{
    //ScopeTimerCL timer( "SetupConvectionP1");

    TetraAccumulatorTupleCL accus;
    InterfaceCommonDataP1CL cdata( ls, lsetbnd);
    accus.push_back( &cdata);
    InterfaceMatrixAccuP1CL<LocalInterfaceConvectionP1CL<DiscVelSolT> > accu( mat, LocalInterfaceConvectionP1CL<DiscVelSolT>( w), cdata);
    accus.push_back( &accu);
    const IdxDescCL* RowIdx= mat->RowIdx;
    accumulate( accus, mg, RowIdx->TriangLevel(), RowIdx->GetMatchingFunction(), RowIdx->GetBndInfo());

    // WriteToFile( mat->Data, "convection.txt", "convection");
}

template <class DiscVelSolT>
  inline InterfaceMatrixAccuP1CL<LocalInterfaceConvectionP1CL<DiscVelSolT> >*
  make_convectionP1_accu (MatDescCL* mat, const InterfaceCommonDataP1CL& cdata, const DiscVelSolT& w, std::string name= std::string())
{
    return new InterfaceMatrixAccuP1CL<LocalInterfaceConvectionP1CL<DiscVelSolT> >( mat,
                   LocalInterfaceConvectionP1CL<DiscVelSolT>( w), cdata, name);
}


template <class DiscVelSolT>
void LocalInterfaceMassDivP1CL<DiscVelSolT>::setup (const TetraCL& t, const InterfaceCommonDataP1CL& cdata)
{
    make_CompositeQuad5Domain2D ( qdom, cdata.surf, t);

    resize_and_evaluate_piecewise_normal( cdata.surf, t, n_tri);
    n.resize( qdom.vertex_size());
    for (Uint i= 0; i < cdata.surf.facet_size(); ++i)
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
        resize_and_evaluate_on_vertexes (cdata.p1[i], qdom, q[i]);

    for (int i= 0; i < 4; ++i) {
        coup[i][i]= quad_2D( qdivgamma_w*q[i]*q[i], qdom);
        for(int j= 0; j < i; ++j)
            coup[i][j]= coup[j][i]= quad_2D( qdivgamma_w*q[j]*q[i], qdom);
    }
}

template <class DiscVelSolT>
void SetupMassDivP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd, const DiscVelSolT& w)
{
    //ScopeTimerCL timer( "SetupMassDivP1");

    TetraAccumulatorTupleCL accus;
    InterfaceCommonDataP1CL cdata( ls, lsetbnd);
    accus.push_back( &cdata);
    InterfaceMatrixAccuP1CL<LocalInterfaceMassDivP1CL<DiscVelSolT> > accu( mat, LocalInterfaceMassDivP1CL<DiscVelSolT>( w), cdata);
    accus.push_back( &accu);
    const IdxDescCL* RowIdx= mat->RowIdx;
    accumulate( accus, mg, RowIdx->TriangLevel(), RowIdx->GetMatchingFunction(), RowIdx->GetBndInfo());

    // WriteToFile( mat->Data, "massdiv.txt", "massdiv");
}

template <class DiscVelSolT>
  inline InterfaceMatrixAccuP1CL<LocalInterfaceMassDivP1CL<DiscVelSolT> >*
  make_massdivP1_accu (MatDescCL* mat, const InterfaceCommonDataP1CL& cdata, const DiscVelSolT& w, std::string name= std::string())
{
    return new InterfaceMatrixAccuP1CL<LocalInterfaceMassDivP1CL<DiscVelSolT> >( mat,
                   LocalInterfaceMassDivP1CL<DiscVelSolT>( w), cdata, name);
}

} // end of namespace DROPS
