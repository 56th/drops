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
#include "num/oswald_projection.h"
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
    double h_;

    bool compute_gradp2_;

  public:
    base_point_newton_cacheCL (const P2EvalCL<double, const NoBndDataCL<>, const VecDescCL>& ls,
                               const P2EvalCL<Point3DCL, const NoBndDataCL<Point3DCL>, const VecDescCL>& ls_grad_rec,
                               const LocalP1CL<Point3DCL> (& gradrefp2)[10])
        : tet( 0), ls_( ls), ls_grad_rec_( ls_grad_rec), gradrefp2_( gradrefp2), compute_gradp2_( true)
    {}

    void set_tetra (const TetraCL* newtet);

    void set_compute_gradp2 (bool b);

    const LocalP2CL<>&          locls  () const { return locls_; }
    const LocalP2CL<Point3DCL>& loc_gh () const { return loc_gh_; }
    const LocalP1CL<Point3DCL>& gradp2 (Uint i) const { return gradp2_[i]; }
    const World2BaryCoordCL&    w2b    () const { return w2b_; }
    double                      get_h  () const { return h_; }
};


class QuaQuaMapperCL
{
  private:
    int maxiter_;
    double tol_;
    int maxinneriter_;
    double innertol_;
    bool use_line_search_;
    double armijo_c_;
    Uint max_damping_steps_;

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
    QuaQuaMapperCL (const MultiGridCL& mg, VecDescCL& lsarg, const VecDescCL& ls_grad_recarg, TetraToTetrasT* neighborhoods= 0, int maxiter= 100, double tol= 1e-7, bool use_line_search= true, double armijo_c= 1e-4, Uint max_damping_steps= 8)
        : maxiter_( maxiter), tol_( tol), maxinneriter_( 100), innertol_( 5e-9),
          use_line_search_( use_line_search), armijo_c_( armijo_c), max_damping_steps_( max_damping_steps),
          ls( &lsarg, &nobnddata, &mg), ls_grad_rec( &ls_grad_recarg, &nobnddata_vec, &mg),
          neighborhoods_( neighborhoods), locator_( mg, lsarg.GetLevel(), /*greedy*/ false),
          cache_( ls, ls_grad_rec, gradrefp2), tet( 0), btet( 0), have_dph( false),
          num_outer_iter( maxiter + 1), num_inner_iter( maxinneriter_ + 1),
          base_point_time( 0.), locate_new_point_time( 0.), cur_num_outer_iter( 0), min_outer_iter(-1u), max_outer_iter( 0),
          total_outer_iter( 0), total_inner_iter( 0), total_damping_iter( 0), total_base_point_calls( 0), total_locate_new_point_calls( 0)
    { P2DiscCL::GetGradientsOnRef( gradrefp2); }

    void set_inner_iter_tol (Uint i, double t) {
        maxinneriter_= i;
        innertol_= t;
        num_inner_iter.resize( i + 1);
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

    mutable double base_point_time,
                   locate_new_point_time;
    mutable Uint cur_num_outer_iter,
                 min_outer_iter,
                 max_outer_iter,
                 total_outer_iter,
                 total_inner_iter,
                 total_damping_iter,
                 total_base_point_calls,
                 total_locate_new_point_calls;
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



// Take LocalP2CL<> ls, check whether there is a linear interface. If yes, map to piecewise quadratic interface and return LocalP2CL<Point3DCL> with displacement for each P2-dof.
// Two modes of operation: 2nd: For each dof, map the dof to the quadratic level set of the linear ls value in dof (Christophs method)
// //                         1st: For each dof map base point on linear interface to base point on quadratic interface.
// class LocalLinearToQuadraticIfaceCL
// {
//   private:
//     enum map_type {MAP_BASE_POINTS, MAP_LEVEL_SETS};
//     map_type t_;
//     bool intersection_p_;
//     LocalP2CL<Point3DCL> map_;
// 
//   public:
//     LocalLinearToQuadraticIfaceCL () : t_(MAP_BASE_POINTS), intersection_p_ (false) {}
//     LocalLinearToQuadraticIfaceCL& assign (const LocalP2CL<>& ls) { return *this; }
// 
//     bool intersection_p () const { return intersection_p_; }
//     const LocalP2CL<Point3DCL>& get_map () const { return map_; }
// };


///\brief returns the L2-projector which turns a LocalP2CL into a LocalP1CL (M_1^{-1} M_{interpolation} M_2).
inline SMatrixCL<4, 10> local_p2_to_p1_L2_projection ()
{
    QRDecompCL<4,4> qrp1;
    SMatrixCL<4,4>& Mp1= qrp1.GetMatrix ();
    for (Uint i= 0; i < 4; ++i) {
        Mp1( i, i)= P1DiscCL::GetMass (i, i);
        for (Uint j= 0; j < i; ++j)
            Mp1 (i, j)= Mp1 (j, i)= P1DiscCL::GetMass (i, j);
    }
    qrp1.prepare_solve ();
    SMatrixCL<10,10> Mp2;
    for (Uint i= 0; i < 10; ++i) {
        Mp2( i, i)= P2DiscCL::GetMass (i, i);
        for (Uint j= 0; j < i; ++j)
            Mp2 (i, j)= Mp2 (j, i)= P2DiscCL::GetMass (i, j);
    }
    SMatrixCL<4,10> Mp1p2;
    LocalP1CL<> p1;
    LocalP2CL<> p2;
    for (Uint i= 0; i < 4; ++i) {
        p1= 0.;
        p1[i]= 1.;
        p2.assign (p1); // Interpolation
        for (Uint j= 0; j < 10; ++j)
            Mp1p2 (i, j)= p2[j];
    }
    qrp1.Solve (Mp1p2);
    return Mp1p2*Mp2;
}


// The function of which we search a root is F(s) = \phi(M\inv*(p - s*v, 1)) =  = \phi(a - s*b) with a=M\inv*(p,1), b= M\inv*(v, 0). Here, \phi is the (local) level set function and M\inv*(x,1) for some matrix M is the affine function mapping world coordinates x to barycentric coordinates.
// Its Jacobian is the dF(s) = -d\phi(p - s*v)*v.
class LocalQuaLineSearchFunctionCL
{
  public:
    typedef double value_type;
    typedef double derivative_type;

  private:
    Point3DCL p, // Point to be projected.
              v; // search direction
    LocalP2CL<> locls;
    LocalP2CL<Point3DCL> locls_grad;
    BaryCoordCL a,
                b;
    double h; // local mesh width at tet.
    World2BaryCoordCL w2b;

    bool xF_p, xdF_p;  // Remember if xF, xdF have been initialized.
    value_type xF, xdF;  // Points at which F, dF are set up.
    BaryCoordCL bxF, bxdF;       // Barycentric coordinates of the spatial part of xF, xdF, xdFinv.
    Point3DCL g_ls_xF, g_ls_xdF; // Gradient of ls at xF, xdF, xdFinv.

    value_type F;           // value at xF;
    derivative_type  dF;    // Jacobian at xdFinv.

    void compute_F (const value_type& x);     // Compute F and set xF.
    void compute_dF (const value_type& x);    // Compute dF and set xdF.

  public:
    LocalQuaLineSearchFunctionCL (
        const LocalP2CL<>& loclsarg,
        const LocalP1CL<Point3DCL>& locls_gradarg,
        double harg,
        const World2BaryCoordCL& w2barg)
        : locls (loclsarg), locls_grad (locls_gradarg),
          h (harg), w2b (w2barg), xF_p (false), xdF_p (false) {}
    LocalQuaLineSearchFunctionCL () : xF_p (false), xdF_p (false) {}

    LocalQuaLineSearchFunctionCL& set_tetra (
        const LocalP2CL<>& loclsarg,
        const LocalP1CL<Point3DCL>& locls_gradarg,
        double harg,
        const World2BaryCoordCL& w2barg) {
        locls= loclsarg;
        locls_grad= locls_gradarg;
        h= harg;
        w2b= w2barg;
        xF_p= false;
        xdF_p= false;
        return *this;
    }

    LocalQuaLineSearchFunctionCL& set_point_and_direction (const Point3DCL& parg, const Point3DCL& varg) {
        p= parg;
        v= varg;
        a= w2b (p);
        b= w2b.map_direction (v);
        return *this;
    }

    value_type value (const value_type& x);
    value_type apply_derivative (const value_type& x, const value_type& v);
    value_type apply_derivative_inverse (const value_type& x, const value_type& v);
    value_type apply_derivative_transpose (const value_type& x, const value_type& v);

    double initial_damping_factor (const value_type& x, const value_type& dx, const value_type& F);
};


// The function of which we search a root is F(x, s) = ( p - x - s*gh(x), -(ls(x) - level_value ) ).
// Its Jacobian is the blockmatrix dF(x, s) = (-I - s H_ls(x) | - g_ls(x), -g_ls(x)^T | 0).
class LocalQuaMapperFunctionCL
{
  public:
    typedef SVectorCL<4> value_type;

  private:
    Point3DCL p;          // Point to be projected.
    LocalP2CL<> locls;
    double level_value;
    LocalP1CL<Point3DCL> locls_grad;
    SMatrixCL<3, 3> locls_H; // Hessian of ls on tet.
    double h;             // local mesh width at tet.
    World2BaryCoordCL w2b;
    double lower_bary;   // default -0.5; lower bound for the barycentric coordinates in initial_damping_factor.

    bool xF_p, xdF_p, xdFinv_p;  // Remember if xF, xdF, xdFinv have been initialized.
    value_type xF, xdF, xdFinv;  // Points at which F, dF, dFinv are set up.
    BaryCoordCL bxF, bxdF;       // Barycentric coordinates of the spatial part of xF, xdF, xdFinv.
    Point3DCL g_ls_xF, g_ls_xdF; // Gradient of ls at xF, xdF, xdFinv.

    value_type F;           // value at xF;
    SMatrixCL<4, 4>  dF;    // Jacobian at xdFinv.
    QRDecompCL<4, 4> dFinv; // Solver for Jacobian at xdFinv.

    void compute_F (const value_type& x);     // Compute F and set xF.
    void compute_dF (const value_type& x);    // Compute dF and set xdF.
    void compute_dFinv (const value_type& x); // Compute dFinv and set xdFinv.

  public:
    LocalQuaMapperFunctionCL (
        const LocalP2CL<>& loclsarg,
        const LocalP1CL<Point3DCL>& locls_gradarg,
        const SMatrixCL<3, 3>& locls_Harg,
        double harg,
        const World2BaryCoordCL& w2barg)
        : locls (loclsarg), level_value (0.), locls_grad (locls_gradarg), locls_H (locls_Harg),
          h (harg), w2b (w2barg), lower_bary (-0.5), xF_p (false), xdF_p (false), xdFinv_p (false) {}
    LocalQuaMapperFunctionCL () : level_value (0.), lower_bary (-0.5), xF_p (false), xdF_p (false), xdFinv_p (false) {}

    LocalQuaMapperFunctionCL& set_tetra (
        const LocalP2CL<>& loclsarg,
        const LocalP1CL<Point3DCL>& locls_gradarg,
        const SMatrixCL<3, 3>& locls_Harg,
        double harg,
        const World2BaryCoordCL& w2barg) {
        locls= loclsarg;
        locls_grad= locls_gradarg;
        locls_H= locls_Harg;
        h= harg;
        w2b= w2barg;
        xF_p= false;
        xdF_p= false;
        xdFinv_p= false;
        return *this;
    }

    LocalQuaMapperFunctionCL& set_point (const Point3DCL& x) {
        p= x;
        return *this;
    }

    LocalQuaMapperFunctionCL& set_level_value (double lsval) {
        level_value= lsval;
        return *this;
    }

    LocalQuaMapperFunctionCL& set_lower_bound_bary (double lb) {
        lower_bary= lb;
        return *this;
    }

    value_type value (const value_type& x);
    value_type apply_derivative (const value_type& x, const value_type& v);
    value_type apply_derivative_inverse (const value_type& x, const value_type& v);
    value_type apply_derivative_transpose (const value_type& x, const value_type& v);

    double initial_damping_factor (const value_type& x, const value_type& dx, const value_type& F);

    const LocalP1CL<Point3DCL>& get_locls_grad () const { return locls_grad; }
};

class LocalQuaMapperCL
{
  public:
    enum DeformationMethodE {MAP_LOCAL_LEVEL_SETS, MAP_ZERO_LEVEL_SETS};

  private:
    int maxiter_;
    double tol_;
    double armijo_c_;
    Uint max_damping_steps_;

    // The level set function.
    NoBndDataCL<> nobnddata;
    P2EvalCL<double, const NoBndDataCL<>, const VecDescCL> ls;

// //     mutable SMatrixCL<4,10> p2top1;
    mutable LocalP1CL<> loclsp1;
    mutable Point3DCL gp1;
//     mutable double c_lin_dist;

    LocalP1CL<Point3DCL> gradrefp2[10];

    mutable World2BaryCoordCL w2b;
    mutable Bary2WorldCoordCL b2w;

    mutable const TetraCL* tet;
    mutable BaryCoordCL xb; // point to be projected
    mutable BaryCoordCL bxb; // projection of xb
    mutable double dh;
    mutable bool have_base_point,
                 base_in_trust_region;

    mutable DeformationMethodE deformation_method;
    mutable bool have_deformation;
    mutable Point3DCL deformation;

    mutable LocalQuaMapperFunctionCL localF;

    void base_point_newton () const;

  public:
    LocalQuaMapperCL (const MultiGridCL& mg, VecDescCL& lsarg, int maxiter= 100, double tol= 1e-7, double armijo_c= 1e-4, Uint max_damping_steps= 8)
        : maxiter_( maxiter), tol_( tol),
          armijo_c_( armijo_c), max_damping_steps_( max_damping_steps),
          ls( &lsarg, &nobnddata, &mg),
          tet( 0), have_base_point (false), deformation_method (MAP_ZERO_LEVEL_SETS), have_deformation (false),
          num_outer_iter( maxiter + 1),
          base_point_time( 0.), locate_new_point_time( 0.), cur_num_outer_iter( 0), min_outer_iter(-1u), max_outer_iter( 0),
          total_outer_iter( 0), total_damping_iter( 0), total_base_point_calls( 0) {
        P2DiscCL::GetGradientsOnRef( gradrefp2);
//         p2top1= local_p2_to_p1_L2_projection ();
    }

    const LocalQuaMapperCL& set_point (const BaryCoordCL& xbarg) const;
    const LocalQuaMapperCL& set_tetra (const TetraCL* tetarg) const;
    const LocalQuaMapperCL& base_point () const;
    const LocalQuaMapperCL& set_trust_region (double lb) {
        localF.set_lower_bound_bary (-lb);
        return *this;
    }
    const LocalQuaMapperCL& set_deformation_method (DeformationMethodE m) const {
        deformation_method= m;
        return *this;
    }

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

    mutable double base_point_time,
                   locate_new_point_time;
    mutable Uint cur_num_outer_iter,
                 min_outer_iter,
                 max_outer_iter,
                 total_outer_iter,
                 total_damping_iter,
                 total_base_point_calls;
};

// Compute the average of a LocalQuaMapperCL in all P2-dofs.
class LocalQuaMapperP2CL
{
  private:
    LocalQuaMapperCL f_;

    double loc_[10];

  public:
    typedef double value_type;
    static const int num_components= 1;

    LocalQuaMapperP2CL (const LocalQuaMapperCL& f)
        : f_( f) {}

    void set_tetra (const TetraCL* t) {
        f_.set_tetra (t);
        for (Uint i= 0; i < 10; ++i) {
            f_.set_point(FE_P2CL::bary_coord[i])
              .base_point ();
            loc_[i]= f_.base_in_trust_region_p () ? f_.get_dh () : std::numeric_limits<double>::max ();
        }
    }
    value_type&       operator[] (size_t i)       { return loc_[i]; }
    const value_type& operator[] (size_t i) const { return loc_[i]; }
    bool invalid_p (size_t i) const { return loc_[i] == std::numeric_limits<double>::max (); }
    void finalize_accumulation () {
        std::cout << "LocalQuaMapperP2CL::Distribution of outer iterations:\n";
        seq_out( f_.num_outer_iter.begin(), f_.num_outer_iter.end(), std::cout);
    }
};

// Compute the average of the mesh deformation in LocalQuaMapperCL in all P2-dofs.
class LocalQuaMapperDeformationP2CL
{
  private:
    LocalQuaMapperCL f_;

    Point3DCL loc_[10];

  public:
    typedef Point3DCL value_type;
    static const int num_components= 3;

    LocalQuaMapperDeformationP2CL (const LocalQuaMapperCL& f)
        : f_( f) {}

    void set_tetra (const TetraCL* t) {
        f_.set_tetra (t);
        for (Uint i= 0; i < 10; ++i)
            loc_[i]= f_.set_point(FE_P2CL::bary_coord[i])
                       .compute_deformation ()
                       .get_deformation ();
    }
    value_type&       operator[] (size_t i)       { return loc_[i]; }
    const value_type& operator[] (size_t i) const { return loc_[i]; }
    bool invalid_p (size_t i) const { return loc_[i][0] == std::numeric_limits<double>::max (); }
    void finalize_accumulation () {
        std::cout << "LocalQuaMapperDeformationP2CL::Distribution of outer iterations:\n";
        seq_out( f_.num_outer_iter.begin(), f_.num_outer_iter.end(), std::cout);
    }
};

} // end of namespace DROPS

#endif
