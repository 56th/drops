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

inline const Quad5DataCL&
make_Quad5Data ()
{
    static const Quad5DataCL quad5data;
    return quad5data;
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


template <class QuadDataT, AbsdetPolicyEnum AbsdetPolicy, Uint Dim>
  const QuadDomainCodim1CL<Dim>&
  make_CompositeQuadDomainCodim1 (QuadDomainCodim1CL<Dim>& q,
                                  const SPatchCL<Dim>& p,
                                  const typename DimensionTraitsCL<Dim>::WorldBodyT& wb)
{
    q.clear();
    if (p.empty())
        return q;

    const Uint num_nodes= QuadDataT::NumNodesC;
    q.vertexes_.resize( num_nodes*p.facet_size());
    q.weights_.resize(  num_nodes*p.facet_size());

    if (AbsdetPolicy == SpaceProjectedCodim1Absdet && p.normal_empty())
        p.compute_normals( wb);
    else if (AbsdetPolicy == Codim1Absdet && p.absdets_empty())
        p.compute_absdets( wb);

    const typename SPatchCL<Dim>::const_vertex_iterator p_vertexes= p.vertex_begin();
    const typename SPatchCL<Dim>::const_normal_iterator p_normals= p.normal_begin();
    const typename SPatchCL<Dim>::const_absdet_iterator p_absdets= p.absdet_begin();

    const typename QuadDomainCL::WeightContT facet_weights( QuadDataT::Weight, num_nodes);
    double absdet;
    Uint f= 0;
    for (typename SPatchCL<Dim>::const_facet_iterator it= p.facet_begin(); it != p.facet_end(); ++it, ++f) {
        const typename SPatchCL<Dim>::FacetT& facet= *it;
        const typename QuadDomainCodim1CL<Dim>::VertexContT::iterator NodeInBody= q.vertexes_.begin() + f*num_nodes;
        for (Uint i= 0; i < num_nodes; ++i) { // NodeInBody[i]= facet*Node[i]
            NodeInBody[i]= p_vertexes[facet[0]]*QuadDataT::Node[i][0];
            for (Uint j= 1; j < Dim; ++j)
                NodeInBody[i]+= p_vertexes[facet[j]]*QuadDataT::Node[i][j];
        }
        absdet= p.is_boundary_facet( it) ? 0.5 : 1.;
        switch (AbsdetPolicy) {
          case TrivialAbsdet: break;
          case SpaceProjectedCodim1Absdet:
            absdet*= std::sqrt( (1. - p_normals[f][Dim-1])*(1. + p_normals[f][Dim-1]));
            // fall through
          case Codim1Absdet:
            absdet*= p_absdets[f];
            break;
          default: throw DROPSErrCL( "make_CompositeQuadDomainCodim1: Unknown AbsdetPolicy.\n");
        }
        q.weights_[std::slice( f*num_nodes, num_nodes, 1)]= absdet*facet_weights;
    }

    return q;
}

inline const QuadDomain2DCL&
make_CompositeQuad1Domain2D (QuadDomain2DCL& q, const SurfacePatchCL& p, const TetraCL& t)
{
    return make_CompositeQuadDomainCodim1<Quad1_2DDataCL, Codim1Absdet, 3>( q, p, t);
}

inline const QuadDomain2DCL&
make_CompositeQuad2Domain2D (QuadDomain2DCL& q, const SurfacePatchCL& p, const TetraCL& t)
{
    return make_CompositeQuadDomainCodim1<Quad2_2DDataCL, Codim1Absdet, 3>( q, p, t);
}

inline const QuadDomain2DCL&
make_CompositeQuad5Domain2D (QuadDomain2DCL& q, const SurfacePatchCL& p, const TetraCL& t)
{
    return make_CompositeQuadDomainCodim1<Quad5_2DDataCL, Codim1Absdet, 3>( q, p, t);
}

inline const QuadDomainCodim1CL<4>&
make_CompositeQuad2DomainSTCodim1 (QuadDomainCodim1CL<4>& q, const SPatchCL<4>& p, const TetraPrismCL& t)
{
    return make_CompositeQuadDomainCodim1<Quad2DataCL, Codim1Absdet, 4>( q, p, t);
}

inline const QuadDomainCodim1CL<4>&
make_CompositeQuad5DomainSTCodim1 (QuadDomainCodim1CL<4>& q, const SPatchCL<4>& p, const TetraPrismCL& t)
{
    return make_CompositeQuadDomainCodim1<Quad5DataCL, Codim1Absdet, 4>( q, p, t);
}

inline const QuadDomainCodim1CL<4>&
make_CompositeQuad5DomainSTCodim1WithoutAbsdet (QuadDomainCodim1CL<4>& q, const SPatchCL<4>& p, const TetraPrismCL& t)
{
    return make_CompositeQuadDomainCodim1<Quad5DataCL, TrivialAbsdet, 4>( q, p, t);
}

inline const QuadDomainCodim1CL<4>&
make_CompositeQuad2DomainSTCodim1WithoutAbsdet (QuadDomainCodim1CL<4>& q, const SPatchCL<4>& p, const TetraPrismCL& t)
{
    return make_CompositeQuadDomainCodim1<Quad2DataCL, TrivialAbsdet, 4>( q, p, t);
}

inline const QuadDomainCodim1CL<4>&
make_CompositeQuad2DomainSTCodim1SpatialAbsdet (QuadDomainCodim1CL<4>& q, const SPatchCL<4>& p, const TetraPrismCL& t)
{
    return make_CompositeQuadDomainCodim1<Quad2DataCL, SpaceProjectedCodim1Absdet, 4>( q, p, t);
}

inline const QuadDomainCodim1CL<4>&
make_CompositeQuad5DomainSTCodim1SpatialAbsdet (QuadDomainCodim1CL<4>& q, const SPatchCL<4>& p, const TetraPrismCL& t)
{
    return make_CompositeQuadDomainCodim1<Quad5DataCL, SpaceProjectedCodim1Absdet, 4>( q, p, t);
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
        make_CompositeQuadDomainCodim1<QuadDataT, Codim1Absdet, 3>( qdom, partition, t);
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


template <class QuadDataT>
  const QuadDomainCL&
  make_CompositeQuadBndDomain2D (QuadDomainCL& q, const BndTriangPartitionCL& p, const TetraCL& t)
{

    const Uint num_nodes= QuadDataT::NumNodesC;

    q.vertexes_.resize( num_nodes*p.triangle_size());

    q.pos_begin_= q.neg_end_= num_nodes*p.triangle_size( NegTetraC); // will be added later after BndTriangPartitionCL includes triangle_size( NegTetraC)

    q.weights_.resize( num_nodes*p.triangle_size());
    q.all_weights_begin_= 0;
    q.pos_weights_begin_= q.pos_begin_;



    const typename BndTriangPartitionCL::const_vertex_iterator partition_vertexes= p.vertex_begin();
    const typename QuadDomainCL::WeightContT triangle_weights( QuadDataT::Weight, num_nodes);

    Uint beg= 0;
    BaryCoordCL tri_bary[3];
    Point3DCL   tri[3];

    for (BndTriangPartitionCL::const_triangle_iterator it= p.triangle_begin(); it != p.triangle_end();
        ++it, beg+= num_nodes) {
        for (int i= 0; i < 3; ++i) {
            tri_bary[i]= partition_vertexes[(*it)[i]];
            tri[i]= GetWorldCoord( t, tri_bary[i]);
        }
        QuadDataT::SetInterface( tri_bary, q.vertexes_.begin() + beg);
        const double absdet= FuncDet2D( tri[1] - tri[0], tri[2] - tri[0]);
        q.weights_[std::slice( beg, num_nodes, 1)]= absdet*triangle_weights;
    }

    return q;

}

inline const QuadDomainCL&
make_CompositeQuad5BndDomain2D  (QuadDomainCL& q, const BndTriangPartitionCL& p, const TetraCL& t)
{
    return make_CompositeQuadBndDomain2D<Quad5_2DDataCL>( q, p, t);
}

} // end of namespace DROPS
