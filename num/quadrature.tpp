/// \file quadrature.tpp
/// \brief numerical integration at the interface
/// \author LNM RWTH Aachen: Joerg Grande; SC RWTH Aachen:

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
 * Copyright 2011 LNM/SC RWTH Aachen, Germany
 */

#include <memory>
#include "num/discretize.h"
#include "num/lattice-eval.h"

namespace DROPS {

template <class GridFunT>
  typename ValueHelperCL<GridFunT>::value_type
  quad_impl (typename CompositeQuadratureTypesNS::const_weight_iterator w_iter, const GridFunT& f, Uint begin, Uint end)
{
    typedef typename ValueHelperCL<GridFunT>::value_type value_type;
    value_type sum= value_type();
    while (begin != end)
       sum+= (*w_iter++)*f[begin++];
    return sum;
}

template <class GridFunT, class QuadDataT>
  inline typename ValueHelperCL<GridFunT>::value_type
  quad (const GridFunT& f, double absdet, const QuadDataT&)
{
    return quad_impl( QuadDataT::Weight, f, 0, QuadDataT::NumNodesC)*absdet;
}

template <class GridFunT, class QuadDataT, class WeightSelectorT>
  inline typename ValueHelperCL<GridFunT>::value_type
  quad (const GridFunT& f, double absdet, const QuadDataT& q, const WeightSelectorT& weightsel)
{
    return quad_impl( q.weights( weightsel), f, 0, QuadDataT::NumNodesC)*absdet;;
}

template <class GridFunT>
  inline typename ValueHelperCL<GridFunT>::value_type
  quad (const GridFunT& f, double absdet, const QuadDomainCL& dom, TetraSignEnum s)
{
    return quad_impl( dom.weight_begin( s), f, dom.dof_begin( s), dom.dof_end( s))*absdet;
}

template <class GridFunT>
  inline void
  quad (const GridFunT& f, double absdet, const QuadDomainCL& dom,
    typename ValueHelperCL<GridFunT>::value_type& neg_int,
    typename ValueHelperCL<GridFunT>::value_type& pos_int)
{
    neg_int= quad( f, absdet, dom, NegTetraC);
    pos_int= quad( f, absdet, dom, PosTetraC);
}

///\brief Helper to quad_{neg,pos}_integrand
/// Integrate a integrand, that is defined only on either the negative or the positive tetras. It does not work for standard integrands.
template <class GridFunT>
  inline typename ValueHelperCL<GridFunT>::value_type
  quad_single_domain_integrand (const GridFunT& f, double absdet, const QuadDomainCL& dom, TetraSignEnum s)
{
    return quad_impl( dom.weight_begin( s), f, 0, dom.dof_end( s) - dom.dof_begin( s))*absdet;
}

template <class GridFunT>
  inline typename ValueHelperCL<GridFunT>::value_type
  quad_neg_part_integrand (const GridFunT& f, double absdet, const QuadDomainCL& dom)
{
    return quad_single_domain_integrand( f, absdet, dom, NegTetraC);
}

template <class GridFunT>
  inline typename ValueHelperCL<GridFunT>::value_type
  quad_pos_part_integrand (const GridFunT& f, double absdet, const QuadDomainCL& dom)
{
    return quad_single_domain_integrand( f, absdet, dom, PosTetraC);
}

template <class GridFunT, Uint Dim>
  inline typename ValueHelperCL<GridFunT>::value_type
  quad_codim1 (const GridFunT& f, const QuadDomainCodim1CL<Dim>& dom)
{
    return quad_impl( dom.weight_begin(), f, dom.dof_begin(), dom.dof_end());
}

template <class GridFunT>
  inline typename ValueHelperCL<GridFunT>::value_type
  quad_2D (const GridFunT& f, const QuadDomainCodim1CL<3>& dom)
{
    // return quad_impl( dom.weight_begin(), f, dom.dof_begin(), dom.dof_end());
    return quad_codim1( f, dom);
}

inline const Quad2DataCL&
make_Quad2Data ()
{
    static const Quad2DataCL quad2data;
    return quad2data;
}

template <class QuadDataT>
  const QuadDomainCL&
  make_SimpleQuadDomain (QuadDomainCL& q, const TetraSignEnum& s)
{
    static const QuadDataT quaddata; // ensure that the fields of QuadDataT are initialized.

    q.vertexes_.resize( QuadDataT::NumNodesC);
    std::copy( QuadDataT::Node, QuadDataT::Node + QuadDataT::NumNodesC, q.vertexes_.begin());
    q.pos_begin_= q.neg_end_= (s == NegTetraC ? QuadDataT::NumNodesC : 0);

    q.weights_.resize( QuadDataT::NumNodesC);
    q.weights_= QuadDomainCL::WeightContT( QuadDataT::Weight, QuadDataT::NumNodesC);
    q.pos_weights_begin_= (s == NegTetraC ? QuadDataT::NumNodesC : 0);
    q.all_weights_begin_= 0;

    return q;
}

template <class QuadDataT>
  const QuadDomainCL&
  make_CompositeQuadDomain (QuadDomainCL& q, const TetraPartitionCL& p)
{
    const Uint num_nodes= QuadDataT::NumNodesC;

    q.vertexes_.resize( 0);
    q.vertexes_.reserve( num_nodes*p.tetra_size());
    q.pos_begin_= q.neg_end_= num_nodes*p.tetra_size( NegTetraC);
    q.weights_.resize( num_nodes*p.tetra_size());
    q.all_weights_begin_= 0;
    q.pos_weights_begin_= q.pos_begin_;

    const typename TetraPartitionCL::const_vertex_iterator partition_vertexes= p.vertex_begin();
    const typename QuadDomainCL::WeightContT tetra_weights( QuadDataT::Weight, num_nodes);
    Uint w_begin= 0;
    SMatrixCL<4,4> T;
    double absdet;
    for (typename TetraPartitionCL::const_tetra_iterator it= p.tetra_begin(); it != p.tetra_end();
        ++it, w_begin+= num_nodes) {
        for (int i= 0; i < 4; ++i)
            T.col( i, partition_vertexes[(*it)[i]]);
        for (Uint i= 0; i < num_nodes; ++i)
            q.vertexes_.push_back( T*QuadDataT::Node[i]);
        absdet= std::fabs( VolFrac(T));
        q.weights_[std::slice( w_begin, num_nodes, 1)]= absdet*tetra_weights;
    }
    return q;
}

inline const QuadDomainCL&
make_CompositeQuad5Domain (QuadDomainCL& q, const TetraPartitionCL& p)
{
    return make_CompositeQuadDomain<Quad5DataCL>( q, p);
}

inline const QuadDomainCL&
make_CompositeQuad3Domain (QuadDomainCL& q, const TetraPartitionCL& p)
{
    return make_CompositeQuadDomain<Quad3DataCL>( q, p);
}

template <class SubdivisionT>
  ExtrapolationToZeroCL::ExtrapolationToZeroCL (Uint num_level, const SubdivisionT& s)
    : lattices_( num_level), f0_( num_level)
{
    if (num_level == 0)
        throw DROPSErrCL( "ExtrapolationToZeroCL: At least one level is needed.");

    std::vector<VecT> f( num_level);
    VecT x( num_level);
    for (Uint i= 0; i < num_level; ++i) {
        lattices_[i]= &PrincipalLatticeCL::instance( s( i));
        x[i]= 1./num_intervals( i);
        f[i].resize( num_level);
        f[i][i]= 1.;
    }
    compute_divided_differences( x, f);
    VecT der0( num_level);
    evaluate_newton_polynomial_and_derivative( x, f, 0., f0_, der0);
    eliminate_linear_term( x, f0_, der0);
//    std::cerr.precision(12);
//    for (Uint i= 0; i < num_level; ++i)
//        std::cerr << weights()[i] << ' ';
//    std::cerr << std::endl;
}

/// \brief Multiply the weight for each level with the extrapolation factor and copy it to weights.
void
copy_weights (const std::vector<CompositeQuadratureTypesNS::WeightContT>& w_vec, const std::vector<Uint>& w_pos_begin,
    const std::valarray<double>& w_factor, CompositeQuadratureTypesNS::WeightContT& weights);

template <class QuadDataT, class LocalFET>
  const QuadDomainCL&
  make_ExtrapolatedQuadDomain (QuadDomainCL& q, const LocalFET& ls, const ExtrapolationToZeroCL& extra)
{
    q.vertexes_.resize( 0);
    q.weights_.resize( 0);

    typename QuadDomainCL::VertexContT pos_vertexes; // temporary container for the positive vertexes
    std::vector<QuadDomainCL::WeightContT> w_vec; // the weights for each level
    w_vec.reserve( extra.num_level());
    std::vector<Uint> w_pos_begin; // begin of the positive weights on each level
    w_pos_begin.reserve( extra.num_level());

    TetraPartitionCL partition;
    QuadDomainCL qdom;
    std::valarray<double> ls_val; // values of the level-set function in the lattice-vertexes
    // Accumulate quadrature-points and weights for each level
    for (Uint i= extra.num_level() - 1; i < extra.num_level(); --i) {
        const PrincipalLatticeCL& lat= extra.lattice( i);
        resize_and_evaluate_on_vertexes( ls, lat, ls_val);
        partition.make_partition<SortedVertexPolicyCL, MergeCutPolicyCL>( lat, ls_val);
        make_CompositeQuadDomain<QuadDataT>( qdom, partition);
        // if (i == extra.num_level() - 1 && lat.tetra_size() == partition.tetra_size()) // No interface cut; no extrapolation; this makes sense, if the order through extrapolation is not higher than the order of the base quadrature in the uncut case.
        //     return make_SimpleQuadDomain<QuadDataT>( q, qdom.vertex_size( NegTetraC) > 0 ? NegTetraC : PosTetraC);

        std::copy( qdom.vertex_begin( NegTetraC), qdom.vertex_end( NegTetraC), std::back_inserter( q.vertexes_));
        std::copy( qdom.vertex_begin( PosTetraC), qdom.vertex_end( PosTetraC), std::back_inserter( pos_vertexes));
        w_vec.push_back( QuadDomainCL::WeightContT( qdom.weight_begin(), qdom.vertex_size()));
        w_pos_begin.push_back( qdom.vertex_size( NegTetraC));
    }
    // Setup the data for the quadrature points
    q.pos_begin_= q.neg_end_= q.vertexes_.size();
    q.vertexes_.resize( q.vertexes_.size() + pos_vertexes.size());
    std::copy( pos_vertexes.begin(), pos_vertexes.end(), q.vertexes_.begin() + q.pos_begin_);

    // Compute the extrapolated weights
    q.pos_weights_begin_= q.pos_begin_;
    q.all_weights_begin_= 0;
    copy_weights( w_vec, w_pos_begin, extra.weights(), q.weights_);
    return q;
}

template <class LocalFET>
  const QuadDomainCL&
  make_ExtrapolatedQuad5Domain (QuadDomainCL& q, const LocalFET& ls, const ExtrapolationToZeroCL& extra)
{
    return make_ExtrapolatedQuadDomain<Quad5DataCL>( q, ls, extra);
}


template <Uint Dim>
QuadDomainCodim1CL<Dim>&
QuadDomainCodim1CL<Dim>::operator= (const QuadDomainCodim1CL<Dim>& q)
{
    if ( &q == this)
        return *this;

    vertexes_= q.vertexes_;
    weights_.resize( q.weights_.size());
    weights_= q.weights_;

    return *this;
}


inline double
codim1_absdet (const TetraCL& t, const BaryCoordCL face_vert[3])
{
    const Point3DCL& f0= GetWorldCoord( t, face_vert[0]);
    return FuncDet2D( GetWorldCoord( t, face_vert[1]) - f0,
                      GetWorldCoord( t, face_vert[2]) - f0);
}


inline Point4DCL GetWorldCoord (const TetraPrismCL& p, const STCoordCL& c)
{
    const Point3DCL& x= GetWorldCoord( p.t, c.x_bary);
    return MakePoint4D( x[0], x[1], x[2], (1. - c.t_ref)*p.t0 + c.t_ref*p.t1);
}

inline double
codim1_absdet (const TetraPrismCL& prism, const STCoordCL face_vert[4])
{
    QRDecompCL<4,3> qr;
    SMatrixCL<4,3>& M= qr.GetMatrix();
    const Point4DCL& v0= GetWorldCoord( prism, face_vert[0]);
    for (Uint i= 1; i < 4; ++i)
        M.col( i - 1, GetWorldCoord( prism, face_vert[i]) - v0);

    const bool is_rank_deficient= qr.prepare_solve( /*assume_full_rank*/ false);
    return !is_rank_deficient ? std::fabs( qr.Determinant_R()) : 0.;
}

inline double
codim1_spatial_absdet (const TetraPrismCL& prism, const STCoordCL face_vert[4])
{
    QRDecompCL<4,3> qr;
    SMatrixCL<4,3>& M= qr.GetMatrix();
    const Point4DCL& v0= GetWorldCoord( prism, face_vert[0]);
    for (Uint i= 1; i < 4; ++i)
        M.col( i - 1, GetWorldCoord( prism, face_vert[i]) - v0);

    const bool is_rank_deficient= qr.prepare_solve( /*assume_full_rank*/ false);
    if (is_rank_deficient)
        return 0.;

    Point4DCL tmp;
    tmp[3]= 1.;
    qr.apply_Q( tmp);
    const double ttt= MakePoint3D( tmp[0], tmp[1], tmp[2]).norm();
    return ttt*std::fabs( qr.Determinant_R());
}

template <Uint Dim>
struct TrivialAbsdetCL
{
    static double codim1_absdet (const typename DimensionTraitsCL<Dim>::WorldBodyT&,
                                 const typename DimensionTraitsCL<Dim>::VertexT[Dim])
    { return 1.; }
};

template <Uint Dim>
struct Codim1AbsdetCL
{
    static double codim1_absdet (const typename DimensionTraitsCL<Dim>::WorldBodyT& body,
                                 const typename DimensionTraitsCL<Dim>::VertexT     facet_vertexes[Dim])
    { return DROPS::codim1_absdet( body, facet_vertexes); }
};

template <Uint Dim>
struct SpatialAbsdetCL
{
    static double codim1_absdet (const typename DimensionTraitsCL<Dim>::WorldBodyT& body,
                                 const typename DimensionTraitsCL<Dim>::VertexT     facet_vertexes[Dim])
    { return DROPS::codim1_spatial_absdet( body, facet_vertexes); }
};

template <class QuadDataT, template <Uint> class AbsdetPolicyT, Uint Dim>
  const QuadDomainCodim1CL<Dim>&
  make_CompositeQuadDomainCodim1 (QuadDomainCodim1CL<Dim>& q,
                                  const SPatchCL<Dim>& p,
                                  const typename DimensionTraitsCL<Dim>::WorldBodyT& t)
{
    const Uint num_nodes= QuadDataT::NumNodesC;

    q.vertexes_.clear(); // Zero init needed due to upddating with += below.
    q.vertexes_.resize( num_nodes*p.facet_size());
    q.weights_.resize(  num_nodes*p.facet_size());

    const typename SPatchCL<Dim>::const_vertex_iterator partition_vertexes= p.vertex_begin();
    const typename QuadDomainCL::WeightContT facet_weights( QuadDataT::Weight, num_nodes);
    Uint beg= 0;
    typename DimensionTraitsCL<Dim>::VertexT fac_vert[Dim];
    for (typename SPatchCL<Dim>::const_facet_iterator it= p.facet_begin(); it != p.facet_end();
        ++it, beg+= num_nodes) {
        for (Uint i= 0; i < Dim; ++i)
            fac_vert[i]= partition_vertexes[(*it)[i]];
        SetInterface<Dim>( fac_vert, num_nodes, q.vertexes_.begin() + beg, QuadDataT::Node);
        const double absdet= (p.is_boundary_facet( it) ? 0.5 : 1.)*AbsdetPolicyT<Dim>::codim1_absdet( t, fac_vert);
        q.weights_[std::slice( beg, num_nodes, 1)]= absdet*facet_weights;
    }

    return q;
}

inline const QuadDomain2DCL&
make_CompositeQuad1Domain2D (QuadDomain2DCL& q, const SurfacePatchCL& p, const TetraCL& t)
{
    return make_CompositeQuadDomainCodim1<Quad1_2DDataCL, Codim1AbsdetCL, 3>( q, p, t);
}

inline const QuadDomain2DCL&
make_CompositeQuad2Domain2D (QuadDomain2DCL& q, const SurfacePatchCL& p, const TetraCL& t)
{
    return make_CompositeQuadDomainCodim1<Quad2_2DDataCL, Codim1AbsdetCL, 3>( q, p, t);
}

inline const QuadDomain2DCL&
make_CompositeQuad5Domain2D (QuadDomain2DCL& q, const SurfacePatchCL& p, const TetraCL& t)
{
    return make_CompositeQuadDomainCodim1<Quad5_2DDataCL, Codim1AbsdetCL, 3>( q, p, t);
}

inline const QuadDomainCodim1CL<4>&
make_CompositeQuad2DomainSTCodim1 (QuadDomainCodim1CL<4>& q, const SPatchCL<4>& p, const TetraPrismCL& t)
{
    return make_CompositeQuadDomainCodim1<Quad2DataCL, Codim1AbsdetCL, 4>( q, p, t);
}

inline const QuadDomainCodim1CL<4>&
make_CompositeQuad5DomainSTCodim1 (QuadDomainCodim1CL<4>& q, const SPatchCL<4>& p, const TetraPrismCL& t)
{
    return make_CompositeQuadDomainCodim1<Quad5DataCL, Codim1AbsdetCL, 4>( q, p, t);
}

inline const QuadDomainCodim1CL<4>&
make_CompositeQuad5DomainSTCodim1WithoutAbsdet (QuadDomainCodim1CL<4>& q, const SPatchCL<4>& p, const TetraPrismCL& t)
{
    return make_CompositeQuadDomainCodim1<Quad5DataCL, TrivialAbsdetCL, 4>( q, p, t);
}

inline const QuadDomainCodim1CL<4>&
make_CompositeQuad2DomainSTCodim1WithoutAbsdet (QuadDomainCodim1CL<4>& q, const SPatchCL<4>& p, const TetraPrismCL& t)
{
    return make_CompositeQuadDomainCodim1<Quad2DataCL, TrivialAbsdetCL, 4>( q, p, t);
}

inline const QuadDomainCodim1CL<4>&
make_CompositeQuad2DomainSTCodim1SpatialAbsdet (QuadDomainCodim1CL<4>& q, const SPatchCL<4>& p, const TetraPrismCL& t)
{
    return make_CompositeQuadDomainCodim1<Quad2DataCL, SpatialAbsdetCL, 4>( q, p, t);
}

/// \brief Multiply the weight for each level with the extrapolation factor and copy it to weights.
void
copy_weights_surface (const std::vector<CompositeQuadratureTypesNS::WeightContT>& w_vec,
    const std::valarray<double >& w_factor, CompositeQuadratureTypesNS::WeightContT& weights);

template <class QuadDataT, class LocalFET>
  const QuadDomain2DCL&
  make_ExtrapolatedQuadDomain2D (QuadDomain2DCL& q, const LocalFET& ls, const TetraCL& t, const ExtrapolationToZeroCL& extra)
{
    q.vertexes_.resize( 0);

    std::vector<QuadDomain2DCL::WeightContT> w_vec; // the weights for each level
    w_vec.reserve( extra.num_level());

    SurfacePatchCL partition;
    QuadDomain2DCL qdom;
    std::valarray<double> ls_val; // values of the level-set function in the lattice-vertexes
    // Accumulate quadrature-points and weights for each level
    for (Uint i= extra.num_level()-1; i < extra.num_level(); --i) {
        const PrincipalLatticeCL& lat= extra.lattice( i);
        resize_and_evaluate_on_vertexes( ls, lat, ls_val);
        partition.make_patch<MergeCutPolicyCL>( lat, ls_val);
        make_CompositeQuadDomainCodim1<QuadDataT, 3>( qdom, partition, t);
        std::copy( qdom.vertex_begin(), qdom.vertex_end(), std::back_inserter( q.vertexes_));
        w_vec.push_back( QuadDomain2DCL::WeightContT( qdom.weight_begin(), qdom.vertex_size()));
    }

    // Compute the extrapolated weights
    copy_weights_surface( w_vec, extra.weights(), q.weights_);
    return q;

}

template <class LocalFET>
  inline const QuadDomain2DCL&
  make_ExtrapolatedQuad5Domain2D (QuadDomain2DCL& q, const LocalFET& ls, const TetraCL& t, const ExtrapolationToZeroCL& extra)
{
    return make_ExtrapolatedQuadDomain2D<Quad5_2DDataCL>( q, ls, t, extra);
}

} // end of namespace DROPS
