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
 * Copyright 2014 LNM RWTH Aachen, Germany
*/

#ifndef DROPS_LEVELSETMAPPER_H
#define DROPS_LEVELSETMAPPER_H

#include "geom/locator.h"
#include "geom/subtriangulation.h"
#include "num/fe.h"
#include "num/discretize.h"
#include "misc/problem.h"

#include <tr1/unordered_map>
#include <tr1/unordered_set>

namespace DROPS
{

typedef std::tr1::unordered_set<const TetraCL*>            TetraSetT;
typedef std::tr1::unordered_map<const TetraCL*, TetraSetT> TetraToTetrasT;

typedef std::pair<const TetraCL*, BaryCoordCL> TetraBaryPairT;

class base_point_newton_cacheCL
{
  private:
    const TetraCL* tet;

    const P2EvalCL<double, const NoBndDataCL<>, const VecDescCL>&             ls_;
    const P2EvalCL<Point3DCL, const NoBndDataCL<Point3DCL>, const VecDescCL>& ls_grad_rec_;

    const LocalP1CL<Point3DCL> (& gradrefp2_)[10];

    LocalP2CL<>          locls_;
    LocalP2CL<Point3DCL> loc_gh_;
    LocalP1CL<Point3DCL> gradp2_[10];
    World2BaryCoordCL    w2b_;

    SMatrixCL<3,3> hessp2_[10];
    double h_;

  public:
    base_point_newton_cacheCL (const P2EvalCL<double, const NoBndDataCL<>, const VecDescCL>& ls,
                               const P2EvalCL<Point3DCL, const NoBndDataCL<Point3DCL>, const VecDescCL>& ls_grad_rec,
                               const LocalP1CL<Point3DCL> (& gradrefp2)[10])
        : tet( 0), ls_( ls), ls_grad_rec_( ls_grad_rec), gradrefp2_( gradrefp2) {}

    void set_tetra (const TetraCL* newtet);

    const LocalP2CL<>&          locls  () const { return locls_; }
    const LocalP2CL<Point3DCL>& loc_gh () const { return loc_gh_; }
    const LocalP1CL<Point3DCL>& gradp2 (Uint i) const { return gradp2_[i]; }
    const SMatrixCL<3,3>&       hessp2 (Uint i) const { return hessp2_[i]; }
    const World2BaryCoordCL&    w2b    () const { return w2b_; }
    double                      get_h  () const { return h_; }
};


class QuaQuaMapperCL
{
  private:
    int maxiter_;
    double tol_;
    int maxinneriter_;
    const double innertol_;
    bool use_line_search_;
    double armijo_c_;

    // The level set function.
    NoBndDataCL<> nobnddata;
    P2EvalCL<double, const NoBndDataCL<>, const VecDescCL> ls;

    // The recovered gradient of ls.
    NoBndDataCL<Point3DCL> nobnddata_vec;
    P2EvalCL<Point3DCL, const NoBndDataCL<Point3DCL>, const VecDescCL> ls_grad_rec;

    LocalP1CL<Point3DCL> gradrefp2[10];

    // The neighborhoods around each tetra in which base points are searched for.
    TetraToTetrasT* neighborhoods_;
    mutable MyLocatorCL locator_;

    mutable base_point_newton_cacheCL cache_;

    mutable const TetraCL* tet;
    mutable BaryCoordCL xb;
    mutable const TetraCL* btet;
    mutable BaryCoordCL bxb;
    mutable double dh;
    mutable SMatrixCL<3,3> dph;
    mutable bool have_dph;

    void locate_new_point (const Point3DCL& x, const Point3DCL& dx, const TetraCL*& tet, BaryCoordCL& xb, double& d) const;
    bool line_search (Point3DCL& v, const Point3DCL& nx, const TetraCL*& tetra, BaryCoordCL& bary) const;
    void base_point_with_line_search () const;
    void base_point_newton () const;

  public:
    QuaQuaMapperCL (const MultiGridCL& mg, VecDescCL& lsarg, const VecDescCL& ls_grad_recarg, TetraToTetrasT* neighborhoods= 0, int maxiter= 100, double tol= 1e-7, bool use_line_search= true, double armijo_c= 1e-4)
        : maxiter_( maxiter), tol_( tol), maxinneriter_( 100), innertol_( 5e-9),
          use_line_search_( use_line_search), armijo_c_( armijo_c), 
          ls( &lsarg, &nobnddata, &mg), ls_grad_rec( &ls_grad_recarg, &nobnddata_vec, &mg),
          neighborhoods_( neighborhoods), locator_( mg, lsarg.GetLevel(), /*greedy*/ false),
          cache_( ls, ls_grad_rec, gradrefp2), tet( 0), btet( 0), have_dph( false),
          num_outer_iter( maxiter + 1), num_inner_iter( maxinneriter_ + 1)
    { P2DiscCL::GetGradientsOnRef( gradrefp2); }

    void set_tetra_neighborhoods (TetraToTetrasT& neigborhoods) { neighborhoods_= &neigborhoods; }

    const QuaQuaMapperCL& set_point (const TetraCL* tetarg, const BaryCoordCL& xbarg) const;
    const QuaQuaMapperCL& base_point () const;
    const QuaQuaMapperCL& jacobian   () const;

    double                 get_dh ()         const { return dh; }
    TetraBaryPairT         get_base_point () const { return std::make_pair( btet, bxb); }
    const TetraCL*         get_base_tetra () const { return btet; }
    const BaryCoordCL&     get_base_bary  () const { return bxb; }
    const SMatrixCL<3, 3>& get_jacobian   () const { return dph; }

    /// Return the local level set function and its gradient on tet; only for convenience. @{
    LocalP2CL<> local_ls      (const TetraCL& tet) const { return LocalP2CL<>( tet, ls); }
    Point3DCL   local_ls_grad (const TetraCL& tet, const BaryCoordCL& xb) const;
    ///@}

    // Count number of iterations iter->#computations with iter iterations.
    mutable std::vector<size_t> num_outer_iter;
    mutable std::vector<size_t> num_inner_iter;
};

void compute_tetra_neighborhoods (const DROPS::MultiGridCL& mg, const VecDescCL& lsetPhi, const BndDataCL<>& lsetbnd, const PrincipalLatticeCL& lat, TetraToTetrasT& tetra_neighborhoods);

/// XXX Refactor to deal with patches that have more than 1 facet!
double abs_det (const TetraCL& tet, const BaryCoordCL& xb, const QuaQuaMapperCL& quaqua, const SurfacePatchCL& p);

inline bool is_in_ref_tetra (const BaryCoordCL& b, double eps= 1e-10)
{
    for (int i= 0; i < 4; ++i)
        if (b[i] < -eps || b[i] > 1. + eps)
            return false;
    return true;
}


} // end of namespace DROPS

#endif
