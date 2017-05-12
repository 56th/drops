/// \file levelsetmapper.h
/// \brief Implements mappings from a neighborhood of the levelset to the level set.
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
 * Copyright 2014--2017 LNM RWTH Aachen, Germany
*/

#ifndef DROPS_LEVELSETMAPPER_H
#define DROPS_LEVELSETMAPPER_H

#include "geom/locator.h"
#include "geom/subtriangulation.h"
#include "num/fe.h"
#include "num/discretize.h"
#include "num/newton.h"
#include "num/oswald_projection.h"
#include "misc/problem.h"

#include <tr1/unordered_map>
#include <tr1/unordered_set>

namespace DROPS
{

typedef std::tr1::unordered_set<const TetraCL*>            TetraSetT;
typedef std::tr1::unordered_map<const TetraCL*, TetraSetT> TetraToTetrasT;

typedef std::pair<const TetraCL*, BaryCoordCL> TetraBaryPairT;

// Forward declaration of helper classes of QuaQuaMapperCL.
class base_point_newton_cacheCL;
class QuaQuaMapperFunctionCL;
class QuaQuaMapperLineSearchFunctionCL;

class QuaQuaMapperCL
{
  private:
    friend QuaQuaMapperFunctionCL;
    friend QuaQuaMapperLineSearchFunctionCL;

    mutable NewtonSolverCL newton_solver_,
                           inner_newton_solver_;
    bool use_line_search_;

    // The level set function.
    NoBndDataCL<> nobnddata;
    P2EvalCL<double, const NoBndDataCL<>, const VecDescCL> ls;

    // The recovered gradient of ls.
    NoBndDataCL<Point3DCL> nobnddata_vec;
    P2EvalCL<Point3DCL, const NoBndDataCL<Point3DCL>, const VecDescCL> ls_grad_rec;

    LocalP1CL<Point3DCL> gradrefp2[10];

    // The neighborhoods around each tetra in which base points are searched for.
    TetraToTetrasT* neighborhoods_= nullptr;
    mutable MyLocatorCL locator_;

    mutable std::unique_ptr<base_point_newton_cacheCL> cache_;
    std::unique_ptr<QuaQuaMapperFunctionCL> f_;
    std::unique_ptr<QuaQuaMapperLineSearchFunctionCL> f_line_search_;

    mutable const TetraCL* tet= nullptr;
    mutable BaryCoordCL xb;
    mutable const TetraCL* btet= nullptr;
    mutable BaryCoordCL bxb;
    mutable double dh;
    mutable SMatrixCL<3,3> dph;
    mutable bool have_dph= false;

    void locate_new_point (const Point3DCL& x, const Point3DCL& dx, const TetraCL*& tet, BaryCoordCL& xb, double& d) const;
    bool line_search (Point3DCL& v, const Point3DCL& nx, const TetraCL*& tetra, BaryCoordCL& bary) const;
    void base_point_with_line_search () const;
    void base_point_newton () const;

  public:
    QuaQuaMapperCL (const MultiGridCL& mg, VecDescCL& lsarg, const VecDescCL& ls_grad_recarg,
        TetraToTetrasT* neighborhoods= 0, int maxiter= 100, double tol= 1e-7,
        bool use_line_search= true, double armijo_c= 1e-4, Uint max_damping_steps= 8);
    QuaQuaMapperCL (const QuaQuaMapperCL&);
    ~QuaQuaMapperCL ();

    void set_inner_iter_tol (Uint i, double t) {
        inner_newton_solver_.SetMaxIter (i);
        inner_newton_solver_.SetTol (t);
        num_inner_iter.resize (i + 1);
    }

    void set_tetra_neighborhoods (TetraToTetrasT& neigborhoods) { neighborhoods_= &neigborhoods; }
    MyLocatorCL& get_locator () { return locator_; }

    const QuaQuaMapperCL& set_point (const TetraCL* tetarg, const BaryCoordCL& xbarg) const;
    const QuaQuaMapperCL& base_point () const;
    const QuaQuaMapperCL& jacobian   () const;

    const TetraCL*         get_tetra () const { return tet; }
    const BaryCoordCL&     get_bary  () const { return xb; }

    TetraBaryPairT         get_base_point () const { return std::make_pair( btet, bxb); }
    const TetraCL*         get_base_tetra () const { return btet; }
    const BaryCoordCL&     get_base_bary  () const { return bxb; }
    double                 get_dh ()         const { return dh; }
    const SMatrixCL<3, 3>& get_jacobian   () const { return dph; }

    /// Return the local level set function and its gradient on tet; only for convenience. @{
    LocalP2CL<> local_ls      (const TetraCL& tet) const { return LocalP2CL<>( tet, ls); }
    Point3DCL   local_ls_grad (const TetraCL& tet, const BaryCoordCL& xb) const;
    ///@}

    // Count number of iterations iter->#computations with iter iterations.
    mutable std::vector<size_t> num_outer_iter;
    mutable std::vector<size_t> num_inner_iter;

    mutable double base_point_time= 0.,
                   locate_new_point_time= 0.;
    mutable Uint cur_num_outer_iter= 0,
                 min_outer_iter= -1u,
                 max_outer_iter= 0,
                 total_outer_iter= 0,
                 total_inner_iter= 0,
                 total_damping_iter= 0,
                 total_base_point_calls= 0,
                 total_locate_new_point_calls= 0;
};


void compute_tetra_neighborhoods (const DROPS::MultiGridCL& mg, const VecDescCL& lsetPhi, const BndDataCL<>& lsetbnd, const PrincipalLatticeCL& lat, TetraToTetrasT& tetra_neighborhoods);

/// Deals only with patches that have 1 facet! (Use ProjectedQuadDomain2DCL for the general case)
double abs_det (const TetraCL& tet, const BaryCoordCL& xb, const QuaQuaMapperCL& quaqua, const SurfacePatchCL& p);


///\brief Assumes that \sum_i b_i = 1.
inline bool is_in_ref_tetra (const BaryCoordCL& b, double eps= 1e-11)
{
    for (int i= 0; i < 4; ++i)
        if (b[i] < -eps)
            return false;
    return true; // Implies b[i] < 1 + 3eps
}


// forward declaration of helper of LocalQuaMapperCL.
class LocalQuaMapperFunctionCL;

class LocalQuaMapperCL
{
  public:
    enum DeformationMethodE {MAP_LOCAL_LEVEL_SETS, MAP_ZERO_LEVEL_SETS};

  private:
    mutable NewtonSolverCL newton_solver_;

    // The level set function.
    NoBndDataCL<> nobnddata;
    P2EvalCL<double, const NoBndDataCL<>, const VecDescCL> ls;

    mutable LocalP1CL<> loclsp1;
    mutable Point3DCL gp1;

    LocalP1CL<Point3DCL> gradrefp2[10];

    mutable World2BaryCoordCL w2b;
    mutable Bary2WorldCoordCL b2w;

    mutable const TetraCL* tet= nullptr;
    mutable BaryCoordCL xb; // point to be projected
    mutable BaryCoordCL bxb; // projection of xb
    mutable double dh;
    mutable bool have_base_point= false,
                 base_in_trust_region= false;
    double lower_bary_for_trust_region= -1.;

    mutable DeformationMethodE deformation_method= MAP_ZERO_LEVEL_SETS;
    mutable bool have_deformation= false;
    mutable Point3DCL deformation;

    mutable std::unique_ptr<LocalQuaMapperFunctionCL> localF;

    void base_point_newton () const;

  public:
    LocalQuaMapperCL (const MultiGridCL& mg, VecDescCL& lsarg, int maxiter= 100, double tol= 1e-7, double armijo_c= 1e-4, Uint max_damping_steps= 8);
    LocalQuaMapperCL (const LocalQuaMapperCL&);
    ~LocalQuaMapperCL ();

    const LocalQuaMapperCL& set_point (const BaryCoordCL& xbarg) const;
    const LocalQuaMapperCL& set_tetra (const TetraCL* tetarg) const;
    const LocalQuaMapperCL& base_point () const;
    const LocalQuaMapperCL& set_deformation_method (DeformationMethodE m) const {
        deformation_method= m;
        return *this;
    }
    const LocalQuaMapperCL& set_trust_region (double lb);

    const LocalQuaMapperCL& compute_deformation () const;

    const TetraCL*         get_tetra () const { return tet; }
    const BaryCoordCL&     get_bary  () const { return xb; }

    TetraBaryPairT         get_base_point () const { return std::make_pair( tet, bxb); }
    const TetraCL*         get_base_tetra () const { return tet; }
    const BaryCoordCL&     get_base_bary  () const { return bxb; }
    double                 get_dh ()         const { return dh; }
    bool                   base_in_trust_region_p () const { return base_in_trust_region; }

    Point3DCL get_deformation () const { return deformation; }

    // Count number of iterations iter->#computations with iter iterations.
    mutable std::vector<size_t> num_outer_iter;
    mutable std::vector<size_t> num_inner_iter;

    mutable double base_point_time= 0.,
                   locate_new_point_time= 0.;
    mutable Uint cur_num_outer_iter= 0,
                 min_outer_iter= -1u,
                 max_outer_iter= 0,
                 total_outer_iter= 0,
                 total_damping_iter= 0,
                 total_base_point_calls= 0;
};


// Compute the average of a LocalQuaMapperCL in all P2-dofs.
template <class TraitT>
class LocalQuaMapperP2CL
{
  public:
    using value_type= typename TraitT::value_type;
    static const int num_components= TraitT::num_components;

  private:
    LocalQuaMapperCL f_;

    value_type loc_[10];
    bool valid_[10];

  public:
    LocalQuaMapperP2CL (const LocalQuaMapperCL& f)
        : f_( f) {}

    void set_tetra (const TetraCL* t) {
        f_.set_tetra (t);
        for (Uint i= 0; i < 10; ++i) {
            loc_[i]= TraitT::get (FE_P2CL::bary_coord[i], f_);
            valid_[i]= f_.base_in_trust_region_p();
        }
    }
    value_type&       operator[] (size_t i)       { return loc_[i]; }
    const value_type& operator[] (size_t i) const { return loc_[i]; }
    bool invalid_p (size_t i) const { return !valid_[i]; }
    void finalize_accumulation () {
        std::cout << "LocalQuaMapperP2CL::Distribution of outer iterations:\n";
        seq_out( f_.num_outer_iter.begin(), f_.num_outer_iter.end(), std::cout);
    }
};

struct LocalQuaMapperP2DistanceTraitsCL
{
    using value_type= double;
    static const int num_components= 1;
    static value_type get (const BaryCoordCL& b, const LocalQuaMapperCL& f) {
        return f.set_point (b)
                .base_point ()
                .get_dh ();
    }
};

using LocalQuaMapperDistanceP2CL= LocalQuaMapperP2CL<LocalQuaMapperP2DistanceTraitsCL>;

struct LocalQuaMapperP2DeformationTraitsCL
{
    using value_type= Point3DCL;
    static const int num_components= 3;
    static value_type get (const BaryCoordCL& b, const LocalQuaMapperCL& f) {
        return f.set_point (b)
                .compute_deformation ()
                .get_deformation ();
    }
};

using LocalQuaMapperDeformationP2CL= LocalQuaMapperP2CL<LocalQuaMapperP2DeformationTraitsCL>;

} // end of namespace DROPS

#endif
