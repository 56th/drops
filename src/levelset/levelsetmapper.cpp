/// \file levelsetmapper.cpp
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

#include "levelset/levelsetmapper.h"
#include "misc/scopetimer.h"
#include "num/lattice-eval.h"

namespace DROPS
{

class base_point_newton_cacheCL
{
  private:
    const TetraCL* tet= 0;

    const P2EvalCL<double, const NoBndDataCL<>, const VecDescCL>&             ls_;
    const P2EvalCL<Point3DCL, const NoBndDataCL<Point3DCL>, const VecDescCL>& ls_grad_rec_;

    const LocalP1CL<Point3DCL> (& gradrefp2_)[10];

    LocalP2CL<>          locls_;
    LocalP2CL<Point3DCL> loc_gh_;
    LocalP1CL<Point3DCL> gradp2_[10];
    World2BaryCoordCL    w2b_;
    double h_;

    bool compute_gradp2_= true;

  public:
    base_point_newton_cacheCL (const P2EvalCL<double, const NoBndDataCL<>, const VecDescCL>& ls,
                               const P2EvalCL<Point3DCL, const NoBndDataCL<Point3DCL>, const VecDescCL>& ls_grad_rec,
                               const LocalP1CL<Point3DCL> (& gradrefp2)[10])
        : ls_( ls), ls_grad_rec_( ls_grad_rec), gradrefp2_( gradrefp2)
    {}

    void           set_tetra (const TetraCL* newtet);
    const TetraCL* get_tetra () const { return tet; }

    void set_compute_gradp2 (bool b);
    bool get_compute_gradp2 () const { return compute_gradp2_; }

    const LocalP2CL<>&          locls  () const { return locls_; }
    const LocalP2CL<Point3DCL>& loc_gh () const { return loc_gh_; }
    Point3DCL                   loc_grad_ls (const BaryCoordCL& bary) const {
        Point3DCL tmp;
        for (Uint i= 0; i < 10; ++i)
            tmp+= locls_[i]*gradp2_[i](bary);
        return tmp;
    }
    const LocalP1CL<Point3DCL>& gradp2 (Uint i) const { return gradp2_[i]; }
    const World2BaryCoordCL&    w2b    () const { return w2b_; }
    double                      get_h  () const { return h_; }
};

void base_point_newton_cacheCL::set_tetra (const TetraCL* newtet)
{
    if (tet == newtet)
        return;

    tet= newtet;

    locls_.assign( *tet, ls_);
    loc_gh_.assign( *tet, ls_grad_rec_);

    if (compute_gradp2_) {
        SMatrixCL<3,3> T( Uninitialized);
        double dummy;
        GetTrafoTr( T, dummy, *tet);
        P2DiscCL::GetGradients( gradp2_, gradrefp2_, T);

        h_= ::cbrt( std::abs( dummy));
    }
    else {
        h_= 6.*tet->GetVolume();
    }

    w2b_.assign( *tet);
}

void base_point_newton_cacheCL::set_compute_gradp2 (bool b)
{
    if (tet != 0 && b && !compute_gradp2_) {
        SMatrixCL<3,3> T( Uninitialized);
        double dummy;
        GetTrafoTr( T, dummy, *tet);
        P2DiscCL::GetGradients( gradp2_, gradrefp2_, T);
    }
    compute_gradp2_= b;
}

class QuaQuaMapperFunctionCL
{
  public:
    typedef SVectorCL<4> value_type;

  private:
    mutable QuaQuaMapperCL* quaqua_= 0;

    bool dF_p=    false,
         dFinv_p= false;  // Remember if dF, dFinv have been initialized.
    value_type xcur;        // Point at which F, dF, dFinv are set up.
    const TetraCL* btet= 0; // "
    BaryCoordCL bxb;        // "
    Point3DCL x0; // World coordinates of initial point.
    Point3DCL gh; // Recovered gradient at xcur.

    value_type F;           // value at xcur;
    SMatrixCL<4, 4>  dF;    // Jacobian at xcur.
    QRDecompCL<4, 4> dFinv; // Solver for Jacobian at xcur.

    void compute_F ();     // Compute F at xcur.
    void compute_dF ();    // Compute dF at xcur, requires compute_F.
    void compute_dFinv (); // Compute dFinv at xcur, requires compute_F.

  public:
    QuaQuaMapperFunctionCL (QuaQuaMapperCL* quaqua)
        : quaqua_ (quaqua) {}

    // Set initial point x0, btet, and bxb; compute xcur gh, and F.
    QuaQuaMapperFunctionCL& set_initial_point (const TetraCL* tet, const BaryCoordCL& xb);

    value_type      get_point () const { return xcur; }
    BaryCoordCL     get_bary  () const { return bxb; }
    const TetraCL*  get_tetra () const { return btet; }
    Point3DCL       get_gh    () const { return gh; }

    bool set_point (const value_type& x);
    value_type value ();
    value_type apply_derivative (const value_type& v);
    value_type apply_derivative_inverse ( const value_type& v);
    value_type apply_derivative_transpose (const value_type& v);

    double initial_damping_factor (const value_type& dx, const value_type& F);
};

QuaQuaMapperFunctionCL&
QuaQuaMapperFunctionCL::set_initial_point (const TetraCL* tet, const BaryCoordCL& xb)
{
    btet= tet;
    bxb=  xb;
    x0= GetWorldCoord (*btet, bxb);
    std::copy (x0.begin(), x0.end(), xcur.begin());
    xcur[3]= 0.;

    quaqua_->cache_->set_tetra( btet);
    gh= quaqua_->cache_->loc_gh() (bxb);
    F[0]= F[1]= F[2]= 0.;
    F[3]= -quaqua_->cache_->locls() (bxb);

    dF_p=    false;
    dFinv_p= false;

    return *this;
}

bool
QuaQuaMapperFunctionCL::set_point (const value_type& x)
{
    if (xcur == x)
        return true;

    const Point3DCL xx    (x.begin (),    x.begin  ()   + 3),
                    xcurx (xcur.begin (), xcur.begin () + 3);
    double l= 1.;
    BaryCoordCL newbxb= quaqua_->cache_->w2b()( xx);
    const TetraCL* newbtet= btet;
    try {
        quaqua_->locate_new_point( xcurx, xx - xcurx, newbtet, newbxb, l);
    } catch (DROPSErrCL) {
        return false;
    }
    if (l != 1.) {
        std::cerr << " QuaQuaMapperFunctionCL::set_point: l: " << l << std::endl;
        return false;
    }

    dF_p=    false;
    dFinv_p= false;
    xcur= x;
    btet= newbtet;
    bxb= newbxb;
    quaqua_->cache_->set_tetra (btet);
    gh= quaqua_->cache_->loc_gh() (bxb);
    compute_F ();

    return true;
}

void
QuaQuaMapperFunctionCL::compute_F ()
{
    // Setup Fnew= (p - xnew - snew*gh, -locls( bxbnew)).
    for (Uint i= 0; i < 3; ++i)
        F[i]= x0[i] - xcur[i] - xcur[3]*gh[i];
    F[3]= -quaqua_->cache_->locls() ( bxb);
}

QuaQuaMapperFunctionCL::value_type
QuaQuaMapperFunctionCL::value ()
{
    return F;
}

void
QuaQuaMapperFunctionCL::compute_dF ()
{
    // Evaluate the gradient of locls in bxb: g_ls.
    // Evaluate the Jacobian of gh in bxb: dgh.
    Point3DCL g_ls;
    SMatrixCL<3,3> dgh;
    Point3DCL tmp;
    for (Uint i= 0; i < 10; ++i) {
        tmp= quaqua_->cache_->gradp2 (i) (bxb);
        g_ls+= quaqua_->cache_->locls()[i]*tmp;
        dgh+= outer_product (quaqua_->cache_->loc_gh()[i], tmp);
    }

    // Setup the blockmatrix M= (-I - s dgh | - gh, -g_ls^T | 0).
    for (Uint i= 0; i < 3; ++i) {
        for (Uint j= 0; j < 3; ++j)
            dF (i, j)= -xcur[3]*dgh (i, j);
        dF (i, i)-= 1.;
        dF (i, 3)= -gh[i];
        dF (3, i)= -g_ls[i];
    }
    dF (3, 3)= 0.;

    dF_p= true;
}

QuaQuaMapperFunctionCL::value_type
QuaQuaMapperFunctionCL::apply_derivative (const value_type& v)
{
    if (!dF_p)
        compute_dF();
    return dF*v;
}

QuaQuaMapperFunctionCL::value_type
QuaQuaMapperFunctionCL::apply_derivative_transpose (const value_type& v)
{
    if (!dF_p)
        compute_dF ();
    return transp_mul(dF, v);
}

void
QuaQuaMapperFunctionCL::compute_dFinv ()
{
    if (!dF_p)
        compute_dF ();
    dFinv.GetMatrix ()= dF;
    dFinv.prepare_solve ();
    dFinv_p= true;
}

QuaQuaMapperFunctionCL::value_type
QuaQuaMapperFunctionCL::apply_derivative_inverse (const value_type& v)
{
    if (!dFinv_p)
        compute_dFinv ();
    value_type ret (v);
    dFinv.Solve (ret);
    return ret;
}

double
QuaQuaMapperFunctionCL::initial_damping_factor (const value_type& dx, const value_type& /*F*/)
{
    return 0.5*quaqua_->cache_->get_h()/MakePoint3D( dx[0], dx[1], dx[2]).norm();
}


class QuaQuaMapperLineSearchFunctionCL
{
  public:
    typedef double value_type;

  private:
    mutable QuaQuaMapperCL* quaqua_= 0;

    bool dF_p= false;  // Remember if dF has been initialized.
    value_type xcur;        // Point at which F, dF, dFinv are set up.
    const TetraCL* btet= 0; // "
    BaryCoordCL bxb;        // "
    Point3DCL x0,    // World coordinates of initial point.
              nline; // Direction defining a line through x0.
    Point3DCL gh; // (Recovered) gradient at xcur.
    // If use_recovered_gradient_in_place_of_gradient_= false (the default) this is the method from Joerg Grande "Analysis of Highly Accurate Finite Element Based Algorithms For Computing Distances To Level Sets", SINUM. Otherwise it is the ad hoc method from Arnold Reusken "A finite element level set redistancing method based on gradient recovery", SINUM, with a bit slower (and unproven) convergence.
    bool use_recovered_gradient_in_place_of_gradient_= false;

    value_type F; // value at xcur;
    double  dF;   // Jacobian at xcur.

    void compute_F () {  // Compute F at xcur.
        F= quaqua_->cache_->locls() ( bxb);
    }

    void compute_dF () { // Compute dF at xcur.
        dF= -inner_prod( gh, nline);
        dF_p= true;
    }

  public:
    QuaQuaMapperLineSearchFunctionCL (QuaQuaMapperCL* quaqua)
        : quaqua_ (quaqua) {}

    // Set initial point x0, direction of the line nline, btet, and bxb; compute xcur gh, and F.
    QuaQuaMapperLineSearchFunctionCL& set_initial_point (const TetraCL* tet, const BaryCoordCL& xb, const Point3DCL& x0arg, const Point3DCL& n);

    value_type      get_point () const { return xcur; }
    BaryCoordCL     get_bary  () const { return bxb; }
    const TetraCL*  get_tetra () const { return btet; }

    bool set_point (const value_type& x);
    value_type value () { return F; };
    value_type apply_derivative (const value_type& v) {
        if (!dF_p)
            compute_dF();
        return dF*v;
    }
    value_type apply_derivative_inverse ( const value_type& v) {
        if (!dF_p)
            compute_dF ();
        return v/dF;
    }
    value_type apply_derivative_transpose (const value_type& v) {
        if (!dF_p)
            compute_dF ();
        return dF*v;
    }

    double initial_damping_factor (const value_type& dx, const value_type& /*F*/) {
        return 0.5*quaqua_->cache_->get_h()/std::fabs( dx);
    }
};

QuaQuaMapperLineSearchFunctionCL&
QuaQuaMapperLineSearchFunctionCL::set_initial_point (const TetraCL* tet, const BaryCoordCL& xb, const Point3DCL& x0arg, const Point3DCL& n)
{
    btet= tet;
    bxb=  xb;
    x0= x0arg;
    xcur= 0.;
    nline= n;

    quaqua_->cache_->set_tetra( btet);
    gh= use_recovered_gradient_in_place_of_gradient_ ?
        quaqua_->cache_->loc_gh() (bxb) : quaqua_->cache_->loc_grad_ls (bxb);
    F= quaqua_->cache_->locls() (bxb);

    dF_p= false;

    return *this;
}

bool
QuaQuaMapperLineSearchFunctionCL::set_point (const value_type& x)
{
    if (xcur == x)
        return true;

    double l= 1.;
    BaryCoordCL newbxb= quaqua_->cache_->w2b()( x0 - x*nline);
    const TetraCL* newbtet= btet;
    try {
        quaqua_->locate_new_point( x0 - xcur*nline, (xcur - x)*nline, newbtet, newbxb, l);
    } catch (DROPSErrCL) {
        return false;
    }
    if (l != 1.) {
        std::cerr << " QuaQuaMapperFunctionCL::set_point: l: " << l << std::endl;
        return false;
    }

    dF_p= false;
    xcur= x;
    btet= newbtet;
    bxb= newbxb;
    quaqua_->cache_->set_tetra (btet);
    gh= use_recovered_gradient_in_place_of_gradient_ ?
        quaqua_->cache_->loc_gh() (bxb) : quaqua_->cache_->loc_grad_ls(bxb);
    compute_F ();

    return true;
}


// Return a tetra from neighborhood that contains v up to precision eps in barycentric coordinates.
// Returns 0 on failure.
void enclosing_tetra (const Point3DCL& v, const TetraSetT& neighborhood, double eps, const TetraCL*& tetra, BaryCoordCL& bary)
{
    World2BaryCoordCL w2b;
    for (TetraSetT::const_iterator tit = neighborhood.begin(); tit != neighborhood.end(); ++tit) {
        w2b.assign( **tit);
        bary= w2b( v);
        if (is_in_ref_tetra( bary, eps)) {
            tetra= *tit;
            return;
        }
    }
    tetra= 0;
}

QuaQuaMapperCL::QuaQuaMapperCL (const MultiGridCL& mg, VecDescCL& lsarg, const VecDescCL& ls_grad_recarg,
    TetraToTetrasT* neighborhoods, int maxiter, double tol,
    bool use_line_search, double armijo_c, Uint max_damping_steps)
    : newton_solver_ (maxiter, tol, false, nullptr, armijo_c, max_damping_steps),
      inner_newton_solver_ ( 100, 5e-9, false, nullptr, armijo_c, max_damping_steps),
      use_line_search_( use_line_search),
      ls( &lsarg, &nobnddata, &mg), ls_grad_rec( &ls_grad_recarg, &nobnddata_vec, &mg),
      neighborhoods_( neighborhoods), locator_( mg, lsarg.GetLevel(), /*greedy*/ false),
      cache_( new base_point_newton_cacheCL (ls, ls_grad_rec, gradrefp2)), f_ (new QuaQuaMapperFunctionCL (this)), f_line_search_ (new QuaQuaMapperLineSearchFunctionCL (this)),
      num_outer_iter( maxiter + 1), num_inner_iter( inner_newton_solver_.GetMaxIter() + 1)
{
    P2DiscCL::GetGradientsOnRef( gradrefp2);
}

QuaQuaMapperCL::QuaQuaMapperCL (const QuaQuaMapperCL& q)
    : newton_solver_( q.newton_solver_), inner_newton_solver_( q.inner_newton_solver_),
      use_line_search_( q.use_line_search_),
      ls( q.ls), ls_grad_rec( q.ls_grad_rec),
      neighborhoods_( q.neighborhoods_), locator_( q.locator_),
      cache_( new base_point_newton_cacheCL (ls, ls_grad_rec, gradrefp2)), f_ (new QuaQuaMapperFunctionCL (this)), f_line_search_ (new QuaQuaMapperLineSearchFunctionCL (this)), tet( q.tet), btet( q.btet), have_dph( q.have_dph),
      num_outer_iter( q.num_outer_iter), num_inner_iter(q.num_inner_iter),
      base_point_time( q.base_point_time), locate_new_point_time( q.locate_new_point_time), cur_num_outer_iter( q.cur_num_outer_iter), min_outer_iter(q.min_outer_iter), max_outer_iter( q.max_outer_iter),
      total_outer_iter( q.total_outer_iter), total_inner_iter( q.total_inner_iter), total_damping_iter( q.total_damping_iter), total_base_point_calls( q.total_base_point_calls), total_locate_new_point_calls( q.total_locate_new_point_calls)
{
    ls.SetBndData (&nobnddata);
    ls_grad_rec.SetBndData (&nobnddata_vec);
    P2DiscCL::GetGradientsOnRef( gradrefp2);
    cache_->set_compute_gradp2 (q.cache_->get_compute_gradp2());
    cache_->set_tetra (q.cache_->get_tetra());
}


QuaQuaMapperCL::~QuaQuaMapperCL ()
{}

// Find a tetra in neighborhoods_[tet] enclosing x + d*dx; computes the barycentric coordinates in xb.
// On entry, xb must contain the barycentric coordinates of x + d*dx.
void QuaQuaMapperCL::locate_new_point (const Point3DCL& x, const Point3DCL& dx, const TetraCL*& tet, BaryCoordCL& xb, double& d) const
{
    TimerCL timer;
    ++total_locate_new_point_calls;

    const double eps= 1e-10;
    if (!is_in_ref_tetra( xb, eps)) {
        if (neighborhoods_ != 0) {
            const int max_damping= 10;
            const TetraSetT& neighborhood= neighborhoods_[0][tet];

            for (int k= 0; tet == 0 && k < max_damping; ++k, d*= 0.5)
                enclosing_tetra( x + d*dx, neighborhood, eps, tet, xb);
            if (tet == 0) {
                std::cout << "x: " << x << "\tdx: " << dx << "\td: " << d << "\tx + d*dx: "<< x + d*dx << std::endl;
                throw DROPSErrCL("QuaQuaMapperCL::locate_new_point: Coord not in given tetra set.\n");
            }
        }
        else {
//             locator_.set_epsilon( eps);
            locator_.locate( x + d*dx);
            tet= &locator_.get_tetra();
            xb= locator_.get_bary();
        }
    }
    timer.Stop();
    locate_new_point_time+= timer.GetTime();
}

const QuaQuaMapperCL& QuaQuaMapperCL::set_point (const TetraCL* tetarg, const BaryCoordCL& xbarg) const
{
    if (tet != tetarg || (xb - xbarg).norm() >= 1e-15) {
        tet= tetarg;
        xb= xbarg;
        btet= 0;
        have_dph= false;
    }
    return *this;
}

bool QuaQuaMapperCL::line_search (Point3DCL& x, const Point3DCL& nx, const TetraCL*& tetra, BaryCoordCL& bary) const
{
    f_line_search_->set_initial_point (btet, bxb, x, nx);
    double alpha= f_line_search_->get_point();
    inner_newton_solver_.Solve (*f_line_search_, alpha);
    tetra= f_line_search_->get_tetra();
    bary=  f_line_search_->get_bary();
    x-= alpha*nx;

    total_damping_iter+= inner_newton_solver_.get_num_damping_steps();
    total_inner_iter+= inner_newton_solver_.GetIter();
    ++num_inner_iter[inner_newton_solver_.GetIter()];

    return inner_newton_solver_.GetResid() < inner_newton_solver_.GetTol();
}

void QuaQuaMapperCL::base_point_with_line_search () const
{
    btet= tet;
    bxb= xb;

    Point3DCL x, xold, x0; // World coordinates of current and previous xb and of the inital point.
    x0= x= GetWorldCoord( *btet, bxb);

    Point3DCL n; // Current search direction.

//     cache_->set_compute_gradp2( false);

    bool found_zero_level= false;
    int iter;
    for (iter= 1; iter < newton_solver_.GetMaxIter(); ++iter) {
        xold= x;
        cache_->set_tetra( btet);
        n=  cache_->loc_gh()( bxb);

        found_zero_level= line_search( x, n, btet, bxb);

        if ((xold - x).norm() < newton_solver_.GetTol() && found_zero_level)
            break;
    }
    if (iter >= newton_solver_.GetMaxIter()) {
        std::cout << "QuaQuaMapperCL::base_point_with_line_search: max iteration number exceeded; |x - xold|: " << norm( xold - x) << "\tx: " << x << "\t n: " << n << "\t ls(x): " << ls.val( *btet, bxb) << "found_zero_level: " << found_zero_level << std::endl;
    }
    cur_num_outer_iter= iter;
    ++num_outer_iter[iter];
    dh=  inner_prod( n/n.norm(), x0 - x);
}

void QuaQuaMapperCL::base_point_newton () const
{
    // The function of which we search a root, ( x0 - x - s*gh(x), -ls(x) ).
    f_->set_initial_point (tet, xb);
    SVectorCL<4> x= f_->get_point();
    newton_solver_.Solve (*f_, x);
    btet= f_->get_tetra();
    bxb= f_->get_bary();
    // Compute the quasi-distance dh:
    dh= x[3]*f_->get_gh().norm();

    cur_num_outer_iter= newton_solver_.GetIter();
    ++num_outer_iter[newton_solver_.GetIter()];
    total_damping_iter+= newton_solver_.get_num_damping_steps();
}


const QuaQuaMapperCL& QuaQuaMapperCL::base_point () const
// tet and xb specify the point which is projected to the zero level.
// On return tet and xb are the resulting base point.
{
    ScopeTimerCL scopetimer( "QuaQuaMapperCL::base_point");
    TimerCL timer;

    if (btet == 0) {
        if (use_line_search_ == true)
            base_point_with_line_search();
        else
            base_point_newton();
    }
    timer.Stop();
    base_point_time+= timer.GetTime();
    min_outer_iter= std::min( min_outer_iter, cur_num_outer_iter);
    max_outer_iter= std::max( max_outer_iter, cur_num_outer_iter);
    total_outer_iter+= cur_num_outer_iter;
    ++total_base_point_calls;
    return *this;
}

const QuaQuaMapperCL& QuaQuaMapperCL::jacobian () const
{
    if (have_dph)
        return *this;

    if (btet == 0)
        base_point();

//     cache_->set_compute_gradp2( true);
    cache_->set_tetra( btet);

    // Evaluate the quasi normal field in bxb.
    const LocalP2CL<Point3DCL>& loc_gh= cache_->loc_gh();
    const Point3DCL gh= loc_gh( bxb);
    const Point3DCL q_n= gh/gh.norm();

    // Evaluate the normal to the interface in bxb.
    Point3DCL n;
    const LocalP2CL<>& locls= cache_->locls();
    Point3DCL gradp2_b[10];
    for (Uint i= 0; i < 10; ++i) {
        gradp2_b[i]= cache_->gradp2( i)( bxb);
        n+= locls[i]*gradp2_b[i];
    }
    n/= n.norm();

    // Evaluate the Jacobian of the quasi-normal field in bxb: dn_h= 1/|G_h| P dG_h.
    SMatrixCL<3,3> dgh;
    for (Uint i= 0; i < 10; ++i)
        dgh+= outer_product( loc_gh[i], gradp2_b[i]);
    const SMatrixCL<3,3> dn= 1./gh.norm()*(dgh - outer_product( q_n, transp_mul( dgh, q_n)));

    // Compute Q.
    const SMatrixCL<3,3> Q= eye<3,3>() - outer_product( q_n/inner_prod( q_n, n), n);

    // Compute the Jacobian of p_h.
    QRDecompCL<3,3> qr;
    qr.GetMatrix()= eye<3,3>() + dh*Q*dn;
    qr.prepare_solve();
    Point3DCL tmp;
    for (Uint i= 0; i < 3; ++i) {
        tmp= Q.col( i);
        qr.Solve( tmp);
        dph.col( i, tmp);
    }

    have_dph= true;
    return *this;
}


Point3DCL QuaQuaMapperCL::local_ls_grad (const TetraCL& tet, const BaryCoordCL& xb) const
{
    LocalP2CL<> locls( tet, ls);
    LocalP1CL<Point3DCL> gradp2[10];
    SMatrixCL<3,3> T;
    double dummy;
    GetTrafoTr( T, dummy, tet);
    P2DiscCL::GetGradients( gradp2, gradrefp2, T);
    Point3DCL grad;
    for (Uint i= 0; i < 10; ++i)
        grad+= locls[i]*gradp2[i]( xb);
    return grad;
}

void compute_tetra_neighborhoods (const MultiGridCL& mg, const VecDescCL& lsetPhi, const BndDataCL<>& lsetbnd, const PrincipalLatticeCL& lat, TetraToTetrasT& tetra_neighborhoods)
{
    std::tr1::unordered_set<const VertexCL*> iface_vertices;
    LocalP2CL<> locp2_ls;
    std::valarray<double> ls_loc( lat.vertex_size());
    DROPS_FOR_TRIANG_CONST_TETRA( mg, lsetPhi.GetLevel(), it) {
        locp2_ls.assign( *it, lsetPhi, lsetbnd);
        evaluate_on_vertexes( locp2_ls, lat, Addr( ls_loc));
        if (equal_signs( ls_loc))
            continue;

        tetra_neighborhoods[&*it]; // insert it with an empty set of neighbors.
        for (Uint i= 0; i < 4; ++i)
            iface_vertices.insert( it->GetVertex( i));
    }

    typedef std::tr1::unordered_map<const VertexCL*, TetraSetT> VertexToTetrasT;
    VertexToTetrasT vertex_neighborhoods;
    DROPS_FOR_TRIANG_CONST_TETRA( mg, lsetPhi.GetLevel(), it) {
        for (Uint i= 0; i < 4; ++i)
            if (iface_vertices.count( it->GetVertex( i)) == 1)
                vertex_neighborhoods[it->GetVertex( i)].insert( &*it);
    }
    for (TetraToTetrasT::iterator it= tetra_neighborhoods.begin(); it != tetra_neighborhoods.end(); ++it) {
        const TetraCL& tet= *it->first;
        for (Uint v= 0; v < 4; ++v) {
            const TetraSetT& tset= vertex_neighborhoods[tet.GetVertex( v)];
            it->second.insert( tset.begin(), tset.end());
        }
    }
}

double abs_det (const TetraCL& tet, const BaryCoordCL& xb, const QuaQuaMapperCL& quaqua, const SurfacePatchCL& p)
{
    if (p.empty())
        return 0.;

    // Compute the jacobian of p_h.
    quaqua.set_point( &tet, xb)
          .jacobian();
    const SMatrixCL<3,3>& dph= quaqua.get_jacobian();

    const Bary2WorldCoordCL b2w( *quaqua.get_tetra());
    QRDecompCL<3,2> qr;
    SMatrixCL<3,2>& M= qr.GetMatrix();
    M.col(0, b2w( p.vertex_begin()[1]) - b2w( p.vertex_begin()[0]));
    M.col(1, b2w( p.vertex_begin()[2]) - b2w( p.vertex_begin()[0]));
    qr.prepare_solve();
    SMatrixCL<3,2> U;
    Point3DCL tmp;
    for (Uint i= 0; i < 2; ++i) {
        tmp= std_basis<3>( i + 1);
        qr.apply_Q( tmp);
        U.col( i, tmp);
    }
    const SMatrixCL<2,2> Gram= GramMatrix( dph*U);
    return std::sqrt( Gram(0,0)*Gram(1,1) - Gram(0,1)*Gram(1,0));
}


// The function of which we search a root is F(x, s) = ( p - x - s*gh(x), -(ls(x) - level_value ) ).
// Its Jacobian is the blockmatrix dF(x, s) = (-I - s H_ls(x) | - g_ls(x), -g_ls(x)^T | 0).
class LocalQuaMapperFunctionCL
{
  public:
    typedef SVectorCL<4> value_type;

  private:
    Point3DCL p;          // Point to be projected.
    LocalP2CL<> locls;
    double level_value= 0.;
    LocalP1CL<Point3DCL> locls_grad;
    SMatrixCL<3, 3> locls_H; // Hessian of ls on tet.
    double h;             // local mesh width at tet.
    World2BaryCoordCL w2b;
    double max_bary_step= 0.75;   // default 0.7; max. step size in barycentric coordinates in initial_damping_factor.

    bool dF_p=    false,
         dFinv_p= false; // Remember if dF, dFinv have been initialized.
    value_type xcur;  // Point at which F, dF, dFinv are set up.
    BaryCoordCL bcur;       // Barycentric coordinates of the spatial part of xcur.
    Point3DCL g_ls_cur; // Gradient of ls at xcur.

    value_type F;           // value at xF;
    SMatrixCL<4, 4>  dF;    // Jacobian at xdFinv.
    QRDecompCL<4, 4> dFinv; // Solver for Jacobian at xdFinv.

    void compute_F ();     // Compute F and set xF.
    void compute_dF ();    // Compute dF and set xdF.
    void compute_dFinv (); // Compute dFinv and set xdFinv.

  public:
    LocalQuaMapperFunctionCL (
        const LocalP2CL<>& loclsarg,
        const LocalP1CL<Point3DCL>& locls_gradarg,
        const SMatrixCL<3, 3>& locls_Harg,
        double harg,
        const World2BaryCoordCL& w2barg)
        : locls (loclsarg), locls_grad (locls_gradarg), locls_H (locls_Harg),
          h (harg), w2b (w2barg) {}
    LocalQuaMapperFunctionCL (const LocalQuaMapperFunctionCL&)= default;
    LocalQuaMapperFunctionCL () {}

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
        dF_p= false;
        dFinv_p= false;
        return *this;
    }

    LocalQuaMapperFunctionCL& set_initial_point (const Point3DCL& x) {
        p= x;
        return *this;
    }

    LocalQuaMapperFunctionCL& set_level_value (double lsval) {
        level_value= lsval;
        return *this;
    }

    LocalQuaMapperFunctionCL& set_max_bary_step (double l) {
        max_bary_step= l;
        return *this;
    }

    bool set_point (const value_type& x);
    value_type value ();
    value_type apply_derivative (const value_type& v);
    value_type apply_derivative_inverse (const value_type& v);
    value_type apply_derivative_transpose (const value_type& v);

    double initial_damping_factor (const value_type& dx, const value_type& F);

    const LocalP1CL<Point3DCL>& get_locls_grad () const { return locls_grad; }
};

bool
LocalQuaMapperFunctionCL::set_point (const value_type& x)
{
    if (xcur == x)
        return true;

    xcur= x;
    bcur= w2b (Point3DCL (x.begin (), x.begin () + 3));
    g_ls_cur= locls_grad( bcur);
    dF_p= dFinv_p= false;
    compute_F ();

    return true;
}

void
LocalQuaMapperFunctionCL::compute_F ()
{
    for (Uint i= 0; i < 3; ++i)
        F[i]= p[i] - xcur[i] - xcur[3]*g_ls_cur[i];
    F[3]= -(locls( bcur) - level_value);
}

LocalQuaMapperFunctionCL::value_type
LocalQuaMapperFunctionCL::value ()
{
    return F;
}

void
LocalQuaMapperFunctionCL::compute_dF ()
{
    dF_p= true;

    for (Uint i= 0; i < 3; ++i) {
        for (Uint j= 0; j < 3; ++j)
            dF (i,j)= -xcur[3]*locls_H (i,j);
        dF (i,i)-= 1.;
        dF (i, 3)= -g_ls_cur[i];
        dF (3, i)= -g_ls_cur[i];
    }
    dF( 3,3)= 0.;
}

LocalQuaMapperFunctionCL::value_type
LocalQuaMapperFunctionCL::apply_derivative (const value_type& v)
{
    if (!dF_p)
        compute_dF ();
    return dF*v;
}

LocalQuaMapperFunctionCL::value_type
LocalQuaMapperFunctionCL::apply_derivative_transpose (const value_type& v)
{
    if (!dF_p)
        compute_dF ();
    return transp_mul (dF, v);
}

void
LocalQuaMapperFunctionCL::compute_dFinv ()
{
    dFinv_p= true;

    compute_dF ();
    dFinv.GetMatrix ()= dF;
    dFinv.prepare_solve ();
}

LocalQuaMapperFunctionCL::value_type
LocalQuaMapperFunctionCL::apply_derivative_inverse (const value_type& v)
{
    if (!dFinv_p)
        compute_dFinv ();
    value_type ret (v);
    dFinv.Solve (ret);
    return ret;
}

double
LocalQuaMapperFunctionCL::initial_damping_factor (const value_type& dx, const value_type& /*F*/)
{
    // Choose initial damping factor so small that the step size in barycentric coordinates is at most max_bary_step.
    const BaryCoordCL bdir= w2b.map_direction (Point3DCL( dx.begin (), dx.begin () + 3));
    return bdir.norm() < 1e-12 ? 1. : max_bary_step/bdir.norm();
}


LocalQuaMapperCL::LocalQuaMapperCL (const MultiGridCL& mg, VecDescCL& lsarg, int maxiter, double tol, double armijo_c, Uint max_damping_steps)
    : newton_solver_(maxiter, tol, false, nullptr, armijo_c, max_damping_steps),
      ls( &lsarg, &nobnddata, &mg),
      localF (new LocalQuaMapperFunctionCL()),
      num_outer_iter( newton_solver_.GetMaxIter() + 1)
{
    P2DiscCL::GetGradientsOnRef( gradrefp2);
}

LocalQuaMapperCL::LocalQuaMapperCL (const LocalQuaMapperCL& q)
    : newton_solver_ (q.newton_solver_),
      ls( q.ls), loclsp1( q.loclsp1), gp1( q.gp1),
      w2b( q.w2b), b2w( q.b2w), tet( q.tet), xb( q.xb), bxb( q.bxb),
      dh( q.dh), have_base_point( q.have_base_point),
      base_in_trust_region( q.base_in_trust_region), lower_bary_for_trust_region( q.lower_bary_for_trust_region),
      deformation_method( q.deformation_method), have_deformation( q.have_deformation), deformation( q.deformation),
      localF( new LocalQuaMapperFunctionCL( *q.localF)),
      num_outer_iter(q.num_outer_iter), num_inner_iter( q.num_inner_iter),
      base_point_time( q.base_point_time), locate_new_point_time( q.locate_new_point_time),
      cur_num_outer_iter( q.cur_num_outer_iter), min_outer_iter( q.min_outer_iter), max_outer_iter( q.max_outer_iter), total_outer_iter( q.total_outer_iter), total_damping_iter( q.total_damping_iter), total_base_point_calls( q.total_base_point_calls)
{
    ls.SetBndData (&nobnddata);
    P2DiscCL::GetGradientsOnRef( gradrefp2);
}

LocalQuaMapperCL::~LocalQuaMapperCL () {}

const LocalQuaMapperCL& LocalQuaMapperCL::set_point (const BaryCoordCL& xbarg) const
{
    xb= xbarg;
    have_base_point= false;
    have_deformation= false;
    return *this;
}

const LocalQuaMapperCL& LocalQuaMapperCL::set_tetra (const TetraCL* tetarg) const
{
    if (tet != tetarg) {
        tet= tetarg;
        have_base_point= false;
        have_deformation= false;
        const LocalP2CL<> locls ( *tet, ls);

        SMatrixCL<3,3> T( Uninitialized);
        double dummy;
        GetTrafoTr( T, dummy, *tet);
        LocalP1CL<Point3DCL> gradp2[10];
        P2DiscCL::GetGradients( gradp2, gradrefp2, T);
        LocalP1CL<Point3DCL> loc_ls_grad;
        for (Uint i= 0; i < 10; ++i)
            loc_ls_grad+= locls[i]*gradp2[i];
        const double h= ::cbrt( std::abs( dummy));
        SMatrixCL<3,3> H[10];
        P2DiscCL::GetHessians (H, T);
        SMatrixCL<3,3> H_ls;
        for (Uint i= 0; i < 10; ++i)
            H_ls+= locls[i]*H[i];
        w2b.assign (*tet);
        b2w.assign (*tet);

        std::copy (&locls[0], &locls[4], &loclsp1[0]);
        Point3DCL gradp1[4];
        P1DiscCL::GetGradients( gradp1, T);
        gp1= Point3DCL();
        for (Uint i= 0; i < 4; ++i)
            gp1+= loclsp1[i]*gradp1[i];

        localF->set_tetra (locls, loc_ls_grad, H_ls, h, w2b);
    }
    return *this;
}

const LocalQuaMapperCL& LocalQuaMapperCL::set_trust_region (double lb)
{
    lower_bary_for_trust_region= -lb;
    localF->set_max_bary_step (.75);
    return *this;
}

const LocalQuaMapperCL& LocalQuaMapperCL::base_point () const
// tet and xb specify the point which is projected to the zero level.
// On return tet and xb are the resulting base point.
{
    ScopeTimerCL scopetimer( "LocalQuaMapperCL::base_point");
    TimerCL timer;

    if (!have_base_point) {
        const Point3DCL p= b2w (xb);
        localF->set_initial_point (p);
        // Setup initial value.
        LocalQuaMapperFunctionCL::value_type x;
        std::copy (p.begin (), p.end (), x.begin ());
        x[3]= 0.;
        newton_solver_.Solve (*localF, x);
        const Point3DCL x_spatial (x.begin (), x.begin () + 3);
        bxb= w2b (x_spatial);
        dh= (p - x_spatial).norm () * (x[3] > 0. ? 1. : -1.); // x[3]*loc_ls_grad (bxb).norm ();
        base_in_trust_region= newton_solver_.GetResid() < newton_solver_.GetTol()
            && *std::min_element (bxb.begin(), bxb.end()) > lower_bary_for_trust_region; // We trust the base point, if the newton_solver converged and the base point is close enough to the tetra.
        have_base_point= true;

        min_outer_iter= std::min( min_outer_iter, (Uint) newton_solver_.GetIter());
        max_outer_iter= std::max( max_outer_iter, (Uint) newton_solver_.GetIter());
        total_outer_iter+= newton_solver_.GetIter();
        total_damping_iter+= newton_solver_.get_num_damping_steps();
        ++num_outer_iter[newton_solver_.GetIter()];
    }
    timer.Stop();
    base_point_time+= timer.GetTime();
    ++total_base_point_calls;
    return *this;
}


const LocalQuaMapperCL& LocalQuaMapperCL::compute_deformation () const
{
    ScopeTimerCL scopetimer( "LocalQuaMapperCL::deformation");

    if (have_deformation)
        return *this;

    if (deformation_method == MAP_LOCAL_LEVEL_SETS)
        localF->set_level_value (loclsp1 (xb));
    else
        localF->set_level_value (0.);

    have_base_point= false;
    base_point ();
    if (!base_in_trust_region_p ()) {
        have_deformation= false;
        return *this;
    }

    if (deformation_method == MAP_LOCAL_LEVEL_SETS)
        deformation= b2w (bxb) - b2w (xb);
    else {
        Point3DCL n (localF->get_locls_grad ()(bxb));
        n/= n.norm ();
        deformation= (loclsp1 (xb) - dh)*n;
    }

    have_deformation= true;
    return *this;
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

void
LocalQuaLineSearchFunctionCL::compute_F (const value_type& x)
{
    xF_p= true;
    xF= x;
    bxF= a - x*b;
    F= locls (bxF);
}

LocalQuaLineSearchFunctionCL::value_type
LocalQuaLineSearchFunctionCL::value (const value_type& x)
{
    if (!xF_p || !(xF == x))
        compute_F (x);
    return F;
}

void
LocalQuaLineSearchFunctionCL::compute_dF (const value_type& x)
{
    xdF_p= true;
    xdF= x;

    bxdF= a - x*b;
    g_ls_xdF= locls_grad( bxdF);
    dF= -inner_prod (g_ls_xdF, v);
}

LocalQuaLineSearchFunctionCL::value_type
LocalQuaLineSearchFunctionCL::apply_derivative (const value_type& x, const value_type& v)
{
    if (!xdF_p || !(xdF == x))
        compute_dF (x);
    return dF*v;
}

LocalQuaLineSearchFunctionCL::value_type
LocalQuaLineSearchFunctionCL::apply_derivative_transpose (const value_type& x, const value_type& v)
{
    if (!xdF_p || !(xdF == x))
        compute_dF (x);
    return dF*v; // dF is scalar, hence symmetric.
}

LocalQuaLineSearchFunctionCL::value_type
LocalQuaLineSearchFunctionCL::apply_derivative_inverse (const value_type& x, const value_type& v)
{
    if (!xdF_p || !(xdF == x))
        compute_dF (x);
    return v/dF; // dF is scalar.
}

double
LocalQuaLineSearchFunctionCL::initial_damping_factor (const value_type& x, const value_type& dx, const value_type& /*F*/)
{
    const Point3DCL g_ls= xF_p && x == xF ? g_ls_xF : locls_grad (a - x*b);

    return 0.5*h/std::max(1e-3*h, dx*g_ls.norm());
}


} // end of namespace DROPS
