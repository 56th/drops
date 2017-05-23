/// \file ifacetransp.tpp
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

template <class LocalMatrixT>
  void
  update_global_matrix (MatrixBuilderCL& M, const LocalMatrixT& loc, const IdxT* numr, const IdxT* numc)
{
    const int num_row_dofs= LocalMatrixT::row_fe_type == P1IF_FE ? 4 : 10,
              num_col_dofs= LocalMatrixT::col_fe_type == P1IF_FE ? 4 : 10;

    for (int i= 0; i < num_row_dofs; ++i)
        if (numr[i] != NoIdx)
            for (int j= 0; j < num_col_dofs; ++j)
                if (numc[j] != NoIdx)
                    M( numr[i], numc[j])+= loc.coup[i][j];
}

template <Uint Dim>
  void
  resize_and_scatter_piecewise_normal (const SPatchCL<Dim>& surf, const QuadDomainCodim1CL<Dim>& qdom, std::valarray<typename SPatchCL<Dim>::WorldVertexT>& normal)
{
    normal.resize( qdom.vertex_size());
    if (normal.size() == 0)
        return;
    if (surf.normal_empty()) // As qdom has vertexes, the must be facets, i.e. normals.
        throw DROPSErrCL( "resize_and_scatter_piecewise_normal: normals were not precomputed.\n");

    const Uint NodesPerFacet= qdom.vertex_size()/surf.facet_size();
    if (qdom.vertex_size()%surf.facet_size() != 0)
        throw DROPSErrCL( "resize_and_scatter_piecewise_normal: qdom.vertex_size is not a multiple of surf.facet_size.\n");

    for (Uint i= 0; i < surf.facet_size(); ++i)
        std::fill_n( &normal[i*NodesPerFacet], NodesPerFacet, surf.normal_begin()[i]);
}

template <class T, class ResultIterT>
  inline ResultIterT
  evaluate_on_vertexes (T (*f)(const Point3DCL&, double), const TetraBaryPairVectorT& pos, double t, ResultIterT result_iterator)
{
    BaryEvalCL<T> eval;
    eval.set( f);
    eval.set_time( t);
    const TetraCL* prev_tetra= 0;
    for (Uint i= 0; i < pos.size(); ++i) {
        if (prev_tetra != pos[i].first) {
            prev_tetra= pos[i].first;
            eval.set( *pos[i].first);
        }
        *result_iterator++= eval( pos[i].second);
    }
    return result_iterator;
}

template <class T, class ResultContT>
  inline ResultContT&
  resize_and_evaluate_on_vertexes (T (*f)(const Point3DCL&, double), const TetraBaryPairVectorT& pos, double t, ResultContT& result_container)
{
    result_container.resize( pos.size());
    evaluate_on_vertexes( f, pos, t, sequence_begin( result_container));
    return result_container;
}

template <class PEvalT, class ResultIterT>
  inline ResultIterT
  evaluate_on_vertexes (const PEvalT& f, const TetraBaryPairVectorT& pos, ResultIterT result_iterator)
{
    typename PEvalT::LocalFET loc_f;
    const TetraCL* prev_tetra= 0;
    for (Uint i= 0; i < pos.size(); ++i) {
        if (prev_tetra != pos[i].first) {
            prev_tetra= pos[i].first;
            loc_f.assign( *pos[i].first, f);
        }
        *result_iterator++= loc_f( pos[i].second);
    }
    return result_iterator;
}

template <class PEvalT, class ResultContT>
  inline ResultContT&
  resize_and_evaluate_on_vertexes (const PEvalT& f, const TetraBaryPairVectorT& pos, ResultContT& result_container)
{
    result_container.resize( pos.size());
    evaluate_on_vertexes( f, pos, sequence_begin( result_container));
    return result_container;
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
    InterfaceMatrixAccuCL<LocalInterfaceConvectionP1CL<DiscVelSolT>, InterfaceCommonDataP1CL> accu( mat, LocalInterfaceConvectionP1CL<DiscVelSolT>( w), cdata);
    accus.push_back( &accu);
    const IdxDescCL* RowIdx= mat->RowIdx;
    accumulate( accus, mg, RowIdx->TriangLevel(), RowIdx->GetBndInfo());

    // WriteToFile( mat->Data, "convection.txt", "convection");
}


template <class DiscVelSolT>
void LocalInterfaceMassDivP1CL<DiscVelSolT>::setup (const TetraCL& t, const InterfaceCommonDataP1CL& cdata)
{
    make_CompositeQuad5Domain2D ( qdom, cdata.surf, t);
    resize_and_scatter_piecewise_normal( cdata.surf, qdom, n);

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
    InterfaceMatrixAccuCL<LocalInterfaceMassDivP1CL<DiscVelSolT>, InterfaceCommonDataP1CL> accu( mat, LocalInterfaceMassDivP1CL<DiscVelSolT>( w), cdata);
    accus.push_back( &accu);
    const IdxDescCL* RowIdx= mat->RowIdx;
    accumulate( accus, mg, RowIdx->TriangLevel(), RowIdx->GetBndInfo());

    // WriteToFile( mat->Data, "massdiv.txt", "massdiv");
}

} // end of namespace DROPS
