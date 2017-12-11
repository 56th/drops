/// \file quadrature.h
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

#ifndef DROPS_QUADRATURE_H
#define DROPS_QUADRATURE_H

#include "misc/container.h"
#include "geom/subtriangulation.h"

#include <valarray>

namespace DROPS {

class QuadDomainCL;   ///< forward declaration for quad
template <Uint Dim>
class QuadDomainCodim1CL; ///< forward declaration for quad_codim1

typedef QuadDomainCodim1CL<3> QuadDomain2DCL;


///\brief Integrate on a tetra using the QuadDataCL-rules from num/discretize.h
/// Uses QuadDataT::Weight as weights.
template <class GridFunT, class QuadDataT>
  inline typename ValueHelperCL<GridFunT>::value_type
  quad (const GridFunT& f, double absdet, const QuadDataT&);

///\brief Integrate on a tetra using the QuadDataCL-rules from num/discretize.h with special weights
/// Calls q.weights( weightsel) to determine the weights.
template <class GridFunT, class QuadDataT, class WeightSelectorT>
  inline typename ValueHelperCL<GridFunT>::value_type
  quad (const GridFunT& f, double absdet, const QuadDataT& q, const WeightSelectorT& weightsel);

/// Integration of full-sized integrands, which have size( AllTetraC) components.
///@{
/// \brief Integrate on the negative, the positive part or on the whole tetra.
template <class GridFunT>
  inline typename ValueHelperCL<GridFunT>::value_type
  quad (const GridFunT& f, double absdet, const QuadDomainCL& dom, TetraSignEnum s= AllTetraC);

/// \brief Integrate on the negative and the positive part of a tetra.
template <class GridFunT>
  inline void
  quad (const GridFunT& f, double absdet, const QuadDomainCL& dom,
    typename ValueHelperCL<GridFunT>::value_type& neg_int,
    typename ValueHelperCL<GridFunT>::value_type& pos_int);
///@}

/// Integration of small integrands, which have size( NegTetraC) or size( PosTetraC) components
///@{
/// \brief Integrate an integrand, that is defined only on the negative tetras. It does not work for full-sized integrands. Use quad for the latter.
template <class GridFunT>
  inline typename ValueHelperCL<GridFunT>::value_type
  quad_neg_part_integrand (const GridFunT& f, double absdet, const QuadDomainCL& dom);

/// \brief Integrate an integrand, that is defined only on the positive tetras. It does not work for standard integrands. Use quad for the latter.
template <class GridFunT>
  inline typename ValueHelperCL<GridFunT>::value_type
  quad_pos_part_integrand (const GridFunT& f, double absdet, const QuadDomainCL& dom);
///@}

/// \brief Integrate on a surface-patch of codimension 1.
///@{
template <class GridFunT, Uint Dim>
  inline typename ValueHelperCL<GridFunT>::value_type
  quad_codim1 (const GridFunT& f, const QuadDomainCodim1CL<Dim>& dom);

template <class GridFunT>
  inline typename ValueHelperCL<GridFunT>::value_type
  quad_2D (const GridFunT& f, const QuadDomain2DCL& dom);
///@}

namespace CompositeQuadratureTypesNS {

typedef std::valarray<double> WeightContT;
typedef const double* const_weight_iterator;
typedef       double*       weight_iterator;

} // end of namespace DROPS::CompositeQudratureTypesNS


class ExtrapolationToZeroCL; ///< forward declaration for the factory-method
class Quad2DataCL; ///< forward declaration for the factory method

/// \brief returns a Quad2DataCL-object.
/// Used as a type selector for the quad()-function, that integrates on a tetra with the base-quadrature-rule
inline const Quad2DataCL&
make_Quad2Data ();

class Quad5DataCL; ///< forward declaration for the factory method
inline const Quad5DataCL&
make_Quad5Data ();

/// \brief Create a quadrature rule that equals a QuadDataCL-type rule.
/// This lets one use a QuadDataCL-rule as QuadDomainCL. Generally, one should use the appropriate quad()-function for QuadDataCL-types directly and spare the copying. For obscure extrapolation rules, this approach is not possible, as a QuadDomainCL is required.
/// The template-parameter QuadDataT must be given explicitly.
template <class QuadDataT>
  const QuadDomainCL&
  make_SimpleQuadDomain (QuadDomainCL& q, const TetraSignEnum& s);

/// \brief Create a composite quadrature rule.
/// No sharing of quadrature points is performed. The sequence of weights for the whole tetra is the concatenation of the sequences of weights for the negative and positive dof.
/// The template-parameter QuadDataT must be given explicitly.
/// Helpers for common QuadDataCL are given below.
template <class QuadDataT>
  const QuadDomainCL&
  make_CompositeQuadDomain (QuadDomainCL& q,const TetraPartitionCL& p);

template <class QuadDataT>
  const QuadDomainCL&
  make_CompositeQuadDomainBnd2D (QuadDomainCL& q, const BndTriangPartitionCL& p);


///\brief Initialize q as a composite Quad3DataCL-quadrature-rule.
inline const QuadDomainCL&
make_CompositeQuad3Domain (QuadDomainCL& q, const TetraPartitionCL& p);

///\brief Initialize q as a composite Quad5DataCL-quadrature-rule.
inline const QuadDomainCL&
make_CompositeQuad5Domain (QuadDomainCL& q, const TetraPartitionCL& p);

/// \brief Create a composite quadrature rule of degree 2 with sharing of dof.
/// The vertices (which are all quadrature points) are shared by all adjacent tetras. Thei qudrature points are: [negative barycenters..., ...negative vertexes..., ...zero vertexes..., positive vertexes..., positive barycenters). The sequence of weights for the whole tetra is an appropriately interleaved sum of the sequences of weights for the negative and positive dof. It starts at all_weights_begin_.
const QuadDomainCL&
make_CompositeQuad2Domain (QuadDomainCL& qdom, const TetraPartitionCL& p);

/// \brief Create an extrapolated quadrature rule.
/// No sharing of quadrature points is performed. The sequence of weights for the whole tetra is the concatenation of the sequences of weights for the negative and positive dof.
/// The extrapolation method is determined by extra.
/// ls can be anything that has the interface of e.g. LocalP2CL for evaluation on a tetra.
/// The template-parameter QuadDataT must be given explicitly.
template <class QuadDataT, class LocalFET>
  const QuadDomainCL&
  make_ExtrapolatedQuadDomain (QuadDomainCL& q, const LocalFET& ls, const ExtrapolationToZeroCL& extra);

///\brief Initialize q as an extrapolated Quad5DataCL-quadrature-rule.
/// The extrapolation method is determined by extra.
/// ls can be anything that has the interface of e.g. LocalP2CL for evaluation on a tetra.
template <class LocalFET>
  inline const QuadDomainCL&
  make_ExtrapolatedQuad5Domain (QuadDomainCL& q, const LocalFET& ls, const ExtrapolationToZeroCL& extra);

/// \brief Create a composite quadrature rule for a boundary-triangle cut by moving contact lines.
/// No sharing of quadrature points is performed.
/// The template-parameter QuadDataT must be given explicitly.
/// Helpers for common QuadData_2DCL are given below.
template <class QuadDataT>
  const QuadDomainCL&
  make_CompositeQuadBndDomain2D (QuadDomainCL& q, const BndTriangPartitionCL& p, const TetraCL& t);
///\brief Initialize q as a composite Quad5_2DDataCL-quadrature-rule for boundary triangles
inline const QuadDomainCL&
make_CompositeQuad5BndDomain2D (QuadDomainCL& q, const BndTriangPartitionCL& p, const TetraCL& t);

/// \brief General quadrature-domain
/// A quadrature rule is defined (and implemented) as a collection of quadrature points and a corresponding collection of weights.
class QuadDomainCL
{
  public:
     /// \brief Container for barycentric coordinates of quadrature points.
    typedef LatticePartitionTypesNS::VertexContT           VertexContT;
    typedef LatticePartitionTypesNS::const_vertex_iterator const_vertex_iterator;

     ///\brief Container for the quadrature weights
    typedef CompositeQuadratureTypesNS::WeightContT           WeightContT;
    typedef CompositeQuadratureTypesNS::const_weight_iterator const_weight_iterator;
    typedef CompositeQuadratureTypesNS::      weight_iterator       weight_iterator;

    /// Friend declaration for the factory methods; if their number becomes to big, a more elaborate factory-design is in order.
    ///@{
    template <class QuadDataT>
      friend const QuadDomainCL&
      make_SimpleQuadDomain (QuadDomainCL&, const TetraSignEnum& s);

    template <class QuadDataT>
      friend const QuadDomainCL&
      make_CompositeQuadDomain (QuadDomainCL&, const TetraPartitionCL&);

    friend const QuadDomainCL&
    make_CompositeQuad2Domain (QuadDomainCL&, const TetraPartitionCL&);

    template <class QuadDataT, class LocalFET>
      friend const QuadDomainCL&
      make_ExtrapolatedQuadDomain (QuadDomainCL&, const LocalFET&, const ExtrapolationToZeroCL&);

    template <class QuadDataT>
      friend const QuadDomainCL&
      make_CompositeQuadBndDomain2D (QuadDomainCL& q, const BndTriangPartitionCL& p, const TetraCL& t);
    ///@}

  private:
    VertexContT vertexes_;  ///< sequence of all vertexes; some may be used for both, the positive and the negative domain
    Uint        pos_begin_; ///< begin of the subsequence of vertexes of positive tetras
    Uint        neg_end_;   ///< end of the subsequence of vertexes of negative tetras

    WeightContT weights_; ///< sequence of all weights; if there are vertexes on the interface, which are used for both domains, the weights for the whole domain are appended and all_weights_begin_ > 0.
    Uint        pos_weights_begin_;
    Uint        all_weights_begin_;

  public:
    QuadDomainCL () ///< empty default constructor
        : pos_begin_( 0), neg_end_( 0), weights_( 0), pos_weights_begin_( 0), all_weights_begin_( 0) {}

    /// The default copy-constructor does the right thing
    /// \brief copy assignment: resize the valarray for weights to make it behave like a container
    QuadDomainCL& operator= (const QuadDomainCL&);

    /// \brief sequence of the indices of the vertexes (quadrature points) for the given domain
    ///@{
    Uint dof_begin (TetraSignEnum s= AllTetraC) const
        { return s == PosTetraC ? pos_begin_ : 0; }
    Uint dof_end   (TetraSignEnum s= AllTetraC) const
        { return s == NegTetraC ? neg_end_ : vertexes_.size(); }
    ///@}

    Uint vertex_size (TetraSignEnum s= AllTetraC) const ///< Number of quadrature points in the given domain
        { return dof_end( s) - dof_begin( s); }

    /// \brief Begin of the sequence of weights for integration on the given domain
    const_weight_iterator weight_begin (TetraSignEnum s= AllTetraC) const {
        return Addr( weights_) + (s == NegTetraC ? 0
            : (s == PosTetraC ? pos_weights_begin_
                              : all_weights_begin_));
    }
    weight_iterator weight_begin (TetraSignEnum s= AllTetraC) {
        return Addr( weights_) + (s == NegTetraC ? 0
            : (s == PosTetraC ? pos_weights_begin_
                              : all_weights_begin_));
    }

    /// \brief sequence of quadrature points in the given domain.
    ///@{
    const_vertex_iterator vertex_begin (TetraSignEnum s= AllTetraC) const
        { return vertexes_.begin() + (s == PosTetraC ? pos_begin_ : 0); }
    const_vertex_iterator vertex_end   (TetraSignEnum s= AllTetraC) const
        { return s == NegTetraC ? vertexes_.begin() + neg_end_ : vertexes_.end(); }
    ///@}
};


enum AbsdetPolicyEnum { TrivialAbsdet, Codim1Absdet, SpaceProjectedCodim1Absdet };

/// \brief Create a composite quadrature rule for a surface-patch.
/// No sharing of quadrature points is performed.
/// The template-parameter QuadDataT must be given explicitly.
/// The Policy-Template AbsdetPolicyT depends on the dimension Dim; choices are Codim1Absdet for the standard absdet and TrivialAbsdet for absdet==1, SpaceProjectedCodim1Absdet for the standard absdet divided by \sqrt(1 + (w*n)^2). The latter is actually useful, because on the spacetime, the iterated integral \int_t \int_{\Gamma(t)} equals (approximately) the spacetime-integral \int_\mathcal{G}, where one uses SpatialAbsdet.
/// Helpers for common QuadData_2DCL are given below.
template <class QuadDataT,  AbsdetPolicyEnum AbsdetPolicy, Uint Dim>
  const QuadDomainCodim1CL<Dim>&
  make_CompositeQuadDomainCodim1 (QuadDomainCodim1CL<Dim>& q, const SPatchCL<Dim>& p, const  typename DimensionTraitsCL<Dim>::WorldBodyT& t);

///\brief Initialize q as a composite Quad0_2DDataCL-quadrature-rule.
inline const QuadDomain2DCL&
make_CompositeQuad1Domain2D (QuadDomain2DCL& q, const SurfacePatchCL& p, const TetraCL& t);

///\brief Initialize q as a composite Quad2_2DDataCL-quadrature-rule.
inline const QuadDomain2DCL&
make_CompositeQuad2Domain2D (QuadDomain2DCL& q, const SurfacePatchCL& p, const TetraCL& t);

///\brief Initialize q as a composite Quad5_2DDataCL-quadrature-rule.
inline const QuadDomain2DCL&
make_CompositeQuad5Domain2D (QuadDomain2DCL& q, const SurfacePatchCL& p, const TetraCL& t);

inline const QuadDomainCodim1CL<4>&
make_CompositeQuad2DomainSTCodim1 (QuadDomainCodim1CL<4>& q, const SPatchCL<4>& p, const TetraPrismCL& t);

inline const QuadDomainCodim1CL<4>&
make_CompositeQuad5DomainSTCodim1 (QuadDomainCodim1CL<4>& q, const SPatchCL<4>& p, const TetraPrismCL& t);

inline const QuadDomainCodim1CL<4>&
make_CompositeQuad5DomainSTCodim1WithoutAbsdet (QuadDomainCodim1CL<4>& q, const SPatchCL<4>& p, const TetraPrismCL& t);

inline const QuadDomainCodim1CL<4>&
make_CompositeQuad2DomainSTCodim1WithoutAbsdet (QuadDomainCodim1CL<4>& q, const SPatchCL<4>& p, const TetraPrismCL& t);

inline const QuadDomainCodim1CL<4>&
make_CompositeQuad2DomainSTCodim1SpatialAbsdet (QuadDomainCodim1CL<4>& q, const SPatchCL<4>& p, const TetraPrismCL& t);

inline const QuadDomainCodim1CL<4>&
make_CompositeQuad5DomainSTCodim1SpatialAbsdet (QuadDomainCodim1CL<4>& q, const SPatchCL<4>& p, const TetraPrismCL& t);

/// \brief Create an extrapolated quadrature rule.
/// No sharing of quadrature points is performed.
/// The extrapolation method is determined by extra.
/// ls can be anything that has the interface of e.g. LocalP2CL for evaluation on a tetra.
/// The template-parameter QuadDataT must be given explicitly.
template <class QuadDataT, class LocalFET>
  const QuadDomain2DCL&
  make_ExtrapolatedQuadDomain2D (QuadDomain2DCL&, const LocalFET&, const TetraCL&, const ExtrapolationToZeroCL&);

///\brief Initialize q as an extrapolated Quad5_2DDataCL-quadrature-rule.
/// The extrapolation method is determined by extra.
/// ls can be anything that has the interface of e.g. LocalP2CL for evaluation on a tetra.
template <class LocalFET>
  inline const QuadDomain2DCL&
  make_ExtrapolatedQuad5Domain2D (QuadDomain2DCL&, const LocalFET&, const TetraCL&, const ExtrapolationToZeroCL&);

/// \brief General 2D-quadrature-domain
/// A quadrature rule is defined (and implemented) as a collection of quadrature points and a corresponding collection of weights.
template <Uint Dim>
class QuadDomainCodim1CL
{
  public:
    typedef DimensionTraitsCL<Dim> DimTraitsT;

     /// \brief Container for barycentric coordinates of quadrature points.
    typedef typename DimTraitsT::VertexT               VertexT;
    typedef typename DimTraitsT::VertexContT           VertexContT;
    typedef typename DimTraitsT::const_vertex_iterator const_vertex_iterator;
    typedef typename DimTraitsT::      vertex_iterator       vertex_iterator;

     ///\brief Container for the quadrature weights
    typedef CompositeQuadratureTypesNS::WeightContT           WeightContT;
    typedef CompositeQuadratureTypesNS::const_weight_iterator const_weight_iterator;
    typedef CompositeQuadratureTypesNS::      weight_iterator       weight_iterator;

    /// Friend declaration for the factory methods; if their number becomes to big, a more elaborate factory-design is in order.
    ///@{
    template <class QuadDataT,  AbsdetPolicyEnum AbsdetPolicy, Uint D>
      friend const QuadDomainCodim1CL<D>&
      make_CompositeQuadDomainCodim1 (QuadDomainCodim1CL<D>& q,
                                  const SPatchCL<D>& p,
                                  const  typename DimensionTraitsCL<D>::WorldBodyT& t);

    template <class QuadDataT, class LocalFET>
      friend const QuadDomain2DCL&
      make_ExtrapolatedQuadDomain2D (QuadDomain2DCL&, const LocalFET&, const TetraCL&, const ExtrapolationToZeroCL&);
    ///@}

  private:
    VertexContT vertexes_;  ///< sequence of all vertexes
    WeightContT weights_;

  public:
    QuadDomainCodim1CL () ///< empty default constructor
        {}

    /// The default copy-constructor and (since C++11) the assignment operator do the right thing.

    /// \brief XXX: Hack to incrementally construct mapped QuadDomains with QuaQuaQuadDomainMapperAccuCL.
    void push_back_quad_node (const VertexT& v, double w) {
        vertexes_.push_back( v);
        std::valarray<double> wc( weights_.size() + 1);
        wc[std::slice( 0, weights_.size(), 1)]= weights_;
        wc[weights_.size()]= w;
        weights_.resize( wc.size());
        weights_= wc;
    }

    bool empty () const { return vertexes_.empty(); } ///< True, iff there are no quadrature-points
    void clear () { vertexes_.clear(); weights_.resize( 0); } ///< reset to empty default-state

    /// \brief sequence of the indices of the vertexes (quadrature points)
    ///@{
    Uint dof_begin () const { return 0; }
    Uint dof_end   () const { return vertexes_.size(); }
    ///@}

    Uint vertex_size () const ///< Number of quadrature points
        { return dof_end() - dof_begin(); }

    /// \brief Begin of the sequence of weights for integration
    const_weight_iterator weight_begin () const { return Addr( weights_); }
          weight_iterator weight_begin ()       { return Addr( weights_); }

    /// \brief sequence of quadrature points
    ///@{
    const_vertex_iterator vertex_begin () const { return vertexes_.begin(); }
    const_vertex_iterator vertex_end   () const { return vertexes_.end(); }
    vertex_iterator vertex_begin () { return vertexes_.begin(); }
    vertex_iterator vertex_end   () { return vertexes_.end(); }
    ///@}
};

/// Determine, how many subdivisions of the tetra-edges are required for extrapolation on level i.
///@{
///\brief The step size is halved for each additional level
struct RombergSubdivisionCL {
    Uint operator() (Uint i) const { return 1 << i; }
};

///\brief The step size is 1/(i+1) on level i.
struct HarmonicSubdivisionCL {
    Uint operator() (Uint i) const { return i + 1; }
};
///@}

///\brief Computes the extrapolation weights for an extrapolated quadrature rule on num_level() levels.
/// The step sizes are taken from the subdivisionT-argument to the constructor.
class ExtrapolationToZeroCL
{
  private:
    typedef std::valarray<double> VecT;

    std::vector<const PrincipalLatticeCL*> lattices_; ///< The lattice used to subdivide the tetra on each level.
    VecT f0_; ///< Vector of weighting factors for the level-wise quadrature weights due to extrapolation to zero of the error and due to elimination of the linear error term.

    /// \brief Compute the coefficients of the polynomial interpolating (x_i, w_i) in the Newton-basis.
    /// The computation is performed in-place.
    void compute_divided_differences (const VecT& x, std::vector<VecT>& w);
    /// \brief Evaluate the first derivative of a polynomial and the polynomial itself given in the Newton-basis.
    /// This is nested multiplication combined with the product rule.
    void evaluate_newton_polynomial_and_derivative (const VecT& x, const std::vector<VecT>& c, double p, VecT& f, VecT& der);
    /// \brief Raise the extrapolation order by one by exploiting, that there is no linear term in the error.
    /// The equation 'linear-coefficient = 0' is used to compute the next Newton-basis-coefficient of the extrapolation-polynomial.
    void eliminate_linear_term (const VecT& x, VecT& val0, const VecT& der0);

  public:
    template <class SubdivisionT> ///< Compute the extrapolation weights.
      ExtrapolationToZeroCL (Uint num_level, const SubdivisionT& s);

    ///\brief returns the vector of extrapolation weights. These must be multiplied with the weights of the composite quadrature rule on the corresponding level.
    const std::valarray<double>& weights () const { return f0_; }
    ///\brief return the number of extrapolation-levels.
    Uint num_level () const { return f0_.size(); }
    ///\brief return the number of subdivisions of the tetra's edges on level i
    Uint num_intervals (Uint i) const { return lattices_[i]->num_intervals(); }
    ///\brief return the PrincipalLatticeCL on level i
    const PrincipalLatticeCL& lattice (Uint i) const { return *lattices_[i]; }
};

} // end of namespace DROPS

#include "num/quadrature.tpp"

#endif
