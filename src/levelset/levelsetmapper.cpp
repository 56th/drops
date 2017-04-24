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
 * Copyright 2014 LNM RWTH Aachen, Germany
*/

#include "levelset/levelsetmapper.h"
#include "misc/scopetimer.h"
#include "num/lattice-eval.h"
#include "num/newton.h"

namespace DROPS
{

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

    void           set_tetra (const TetraCL* newtet);
    const TetraCL* get_tetra () const { return tet; }

    void set_compute_gradp2 (bool b);
    bool get_compute_gradp2 () const { return compute_gradp2_; }

    const LocalP2CL<>&          locls  () const { return locls_; }
    const LocalP2CL<Point3DCL>& loc_gh () const { return loc_gh_; }
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

    bool maybe_change_current_point (const value_type& x);

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
    : maxiter_( maxiter), tol_( tol), maxinneriter_( 100), innertol_( 5e-9),
      use_line_search_( use_line_search), armijo_c_( armijo_c), max_damping_steps_( max_damping_steps),
      ls( &lsarg, &nobnddata, &mg), ls_grad_rec( &ls_grad_recarg, &nobnddata_vec, &mg),
      neighborhoods_( neighborhoods), locator_( mg, lsarg.GetLevel(), /*greedy*/ false),
      cache_( new base_point_newton_cacheCL (ls, ls_grad_rec, gradrefp2)), f_ (new QuaQuaMapperFunctionCL (this)), tet( 0), btet( 0), have_dph( false),
      num_outer_iter( maxiter + 1), num_inner_iter( maxinneriter_ + 1),
      base_point_time( 0.), locate_new_point_time( 0.), cur_num_outer_iter( 0), min_outer_iter(-1u), max_outer_iter( 0),
      total_outer_iter( 0), total_inner_iter( 0), total_damping_iter( 0), total_base_point_calls( 0), total_locate_new_point_calls( 0)
{
    P2DiscCL::GetGradientsOnRef( gradrefp2);
}

QuaQuaMapperCL::QuaQuaMapperCL (const QuaQuaMapperCL& q)
    : maxiter_( q.maxiter_), tol_( q.tol_), maxinneriter_( q.maxinneriter_), innertol_( q.innertol_),
      use_line_search_( q.use_line_search_), armijo_c_( q.armijo_c_), max_damping_steps_( q.max_damping_steps_),
      ls( q.ls), ls_grad_rec( q.ls_grad_rec),
      neighborhoods_( q.neighborhoods_), locator_( q.locator_),
      cache_( new base_point_newton_cacheCL (ls, ls_grad_rec, gradrefp2)), f_ (new QuaQuaMapperFunctionCL (this)), tet( q.tet), btet( q.btet), have_dph( q.have_dph),
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
    bool found_newtet;

    double dalpha= 0.;
    double l;
    int inneriter= 0;

    Point3DCL xnew;
    const TetraCL* newtet= tetra;
    BaryCoordCL newbary;
    double Fnew= std::numeric_limits<double>::max(); // Silence false warning.

    const LocalP2CL<>&           locls= cache_->locls();
    const LocalP2CL<Point3DCL>& loc_gh= cache_->loc_gh();

    cache_->set_tetra( tetra);
    bary= cache_->w2b()( x);
    double F= locls( bary);
    for (; inneriter < maxinneriter_; ++inneriter) {
        if (std::abs( F) < innertol_)
            break;

        // Compute undamped Newton correction dalpha.
        const Point3DCL& gradval= loc_gh( bary);
        const double slope= inner_prod( gradval, nx);
        if (std::abs( slope) < 1.0e-8)
            std::cout << "g_phi: " << gradval << "\tgy: " << nx << std::endl;
        dalpha= F/slope;

        l= std::min( 1., 0.5*cache_->get_h()/std::abs( dalpha));
        Uint j;
        found_newtet= false;
        for (j= 0; j < max_damping_steps_; ++j, l*= 0.5) {
            if (l < 1e-7) {
                std::cerr << "QuaQuaMapperCL::line_search: Too much damping. inneriter: " << inneriter <<  " x: " << x << " dalpha: " << dalpha << " ls(x): " << locls( bary) << " l: " << l << " slope: " << slope << std::endl;
            }
            xnew= x - (l*dalpha)*nx;
            cache_->set_tetra( tetra);
            newbary= cache_->w2b()( xnew);
            try {
                locate_new_point( x, -dalpha*nx, newtet, newbary, l);
            }
            catch (DROPSErrCL e) {
                continue;
            }
            found_newtet= true;
            cache_->set_tetra( newtet);
            // Setup Fnew= locls( newbary)).
            Fnew= locls( newbary);

            // Armijo-rule
            if (std::abs( Fnew) < std::abs( F) + armijo_c_*F*slope*dalpha/std::abs( F)*l) {
                break;
            }
        }
        total_damping_iter+= j;
        if (!found_newtet)
            throw DROPSErrCL( "QuaQuaMapperCL::line_search: Could not find the right tetra.\n");
        tetra= newtet;
        x= xnew;
        bary= newbary;
        F= Fnew;
    }

//     if (inneriter >= maxinneriter_)
//         std::cout <<"QuaQuaMapperCL::line_search: Warning: max inner iteration number at x : " << x << " exceeded; ls.val: " << ls.val( *tetra, bary) << std::endl;
    total_inner_iter+= inneriter;
    ++num_inner_iter[inneriter];

    return inneriter < maxinneriter_;
}

void QuaQuaMapperCL::base_point_with_line_search () const
{
    btet= tet;
    bxb= xb;

    double alpha= 0.;
    Point3DCL x, xold, x0; // World coordinates of current and previous xb and of the inital point.
    x0= x= GetWorldCoord( *btet, bxb);

    Point3DCL n; // Current search direction.

//     cache_->set_compute_gradp2( false);
    cache_->set_tetra( btet); // To make get_h() well-defined.

    bool found_zero_level= false;
    int iter;
    for (iter= 1; iter < maxiter_; ++iter) {
        xold= x;
        n=  cache_->loc_gh()( bxb);
        n/= norm( n);
        x= x0 - alpha*n;
        found_zero_level= line_search( x, n, btet, bxb);
        alpha= inner_prod( n, x0 - x);

        if ((xold - x).norm() < tol_ && found_zero_level)
            break;
    }
    if (iter >= maxiter_) {
        std::cout << "QuaQuaMapperCL::base_point_with_line_search: max iteration number exceeded; |x - xold|: " << norm( xold - x) << "\tx: " << x << "\t n: " << n << "\t ls(x): " << ls.val( *btet, bxb) << "found_zero_level: " << found_zero_level << std::endl;
    }
    cur_num_outer_iter= iter;
    ++num_outer_iter[iter];
    dh=  inner_prod( n, x0 - x);
}

void QuaQuaMapperCL::base_point_newton () const
{
    // The function of which we search a root, ( x0 - x - s*gh(x), -ls(x) ).
    f_->set_initial_point (tet, xb);
    size_t iter= maxiter_;
    double tol= tol_;
    size_t max_damping_steps= max_damping_steps_;
    SVectorCL<4> x= f_->get_point();
    newton_solve_1 (*f_, x, iter, tol, max_damping_steps, armijo_c_);
    btet= f_->get_tetra();
    bxb= f_->get_bary();
    // Compute the quasi-distance dh:
    dh= x[3]*f_->get_gh().norm();

    cur_num_outer_iter= iter;
    ++num_outer_iter[iter];
    total_damping_iter+= max_damping_steps;
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


void
LocalQuaMapperFunctionCL::compute_F (const value_type& x)
{
    xF_p= true;
    xF= x;
    bxF= w2b (Point3DCL (x.begin (), x.begin () + 3));
    g_ls_xF= locls_grad( bxF);
    for (Uint i= 0; i < 3; ++i)
        F[i]= p[i] - x[i] - x[3]*g_ls_xF[i];
    F[3]= -(locls( bxF) - level_value);
}

LocalQuaMapperFunctionCL::value_type
LocalQuaMapperFunctionCL::value (const value_type& x)
{
    if (!xF_p || !(xF == x))
        compute_F (x);
    return F;
}

void
LocalQuaMapperFunctionCL::compute_dF (const value_type& x)
{
    xdF_p= true;
    xdF= x;

    if (xF_p && xF == xdF) {
        bxdF= bxF;
        g_ls_xdF= g_ls_xF;
    }
    else {
        bxdF= w2b (Point3DCL (x.begin (), x.begin () + 3));
        g_ls_xdF= locls_grad( bxdF);
    }

    for (Uint i= 0; i < 3; ++i) {
        for (Uint j= 0; j < 3; ++j)
            dF (i,j)= -x[3]*locls_H (i,j);
        dF (i,i)-= 1.;
        dF (i, 3)= -g_ls_xdF[i];
        dF (3, i)= -g_ls_xdF[i];
    }
    dF( 3,3)= 0.;
}

LocalQuaMapperFunctionCL::value_type
LocalQuaMapperFunctionCL::apply_derivative (const value_type& x, const value_type& v)
{
    if (!xdF_p || !(xdF == x))
        compute_dF (x);
// std::cerr << "LocalQuaMapperFunctionCL::apply_derivative: xdF: " << w2b (Point3DCL(x.begin (), x.begin () + 3)) << std::endl;
    return dF*v;
}

LocalQuaMapperFunctionCL::value_type
LocalQuaMapperFunctionCL::apply_derivative_transpose (const value_type& x, const value_type& v)
{
    if (!xdF_p || !(xdF == x))
        compute_dF (x);
    return transp_mul(dF, v);
}

void
LocalQuaMapperFunctionCL::compute_dFinv (const value_type& x)
{
    xdFinv_p= true;
    xdFinv= x;

    compute_dF (x);
    SMatrixCL<4, 4>& M= dFinv.GetMatrix ();
    M= dF;
    dFinv.prepare_solve ();
}

LocalQuaMapperFunctionCL::value_type
LocalQuaMapperFunctionCL::apply_derivative_inverse (const value_type& x, const value_type& v)
{
    if (!xdFinv_p || !(xdFinv == x))
        compute_dFinv (x);
    value_type ret (v);
    dFinv.Solve (ret);
    return ret;
}

double
LocalQuaMapperFunctionCL::initial_damping_factor (const value_type& /*x*/, const value_type& dx, const value_type& /*F*/)
{
    // Choose l so small that all barycentric coordinates of the new point are >= lower_bary.
    const BaryCoordCL bdir= w2b.map_direction (Point3DCL( dx.begin (), dx.begin () + 3));
    double l= 1.;
    for (Uint i= 0; i < 4; ++i)
        if (std::abs (bdir[i]) > 1e-12)
            l= std::min (l, std::abs ((bxdF[i] - lower_bary)/bdir[i]));
    return l;
}


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

        localF.set_tetra (locls, loc_ls_grad, H_ls, h, w2b);
    }
    return *this;
}


const LocalQuaMapperCL& LocalQuaMapperCL::base_point () const
// tet and xb specify the point which is projected to the zero level.
// On return tet and xb are the resulting base point.
{
    ScopeTimerCL scopetimer( "LocalQuaMapperCL::base_point");
    TimerCL timer;

    if (!have_base_point) {
        size_t maxiter= maxiter_;
        double tol= tol_;
        size_t max_damping_steps= max_damping_steps_;
        const Point3DCL p= b2w (xb);
        localF.set_point (p);
        // Setup initial value.
        LocalQuaMapperFunctionCL::value_type x;
        std::copy (p.begin (), p.end (), x.begin ());
        x[3]= 0.;
        newton_solve (localF, x, maxiter, tol, max_damping_steps, armijo_c_);
        const Point3DCL x_spatial (x.begin (), x.begin () + 3);
        bxb= w2b (x_spatial);
        dh= (p - x_spatial).norm () * (x[3] > 0. ? 1. : -1.); // x[3]*loc_ls_grad (bxb).norm ();
//         if (maxiter >= (size_t) maxiter_) {
//             dh= std::numeric_limits<double>::max();
//             size_t maxiter= maxiter_;
//             double tol= tol_;
//             size_t max_damping_steps= max_damping_steps_;
//             std::copy (p.begin (), p.end (), x.begin ());
//             x[3]= 0.;
//             newton_solve (localF, x, maxiter, tol, max_damping_steps, armijo_c_, &std::cerr);
//             std::cerr << "Press enter to continue.";
//             char tmp;
//             std::cin >> tmp;
//         }
        base_in_trust_region= maxiter < (size_t) maxiter_ && tol < tol_; // Otherwise, the iteration was terminated prematurely because the boundary of the trust region was reached.
        have_base_point= true;

        min_outer_iter= std::min( min_outer_iter, (Uint) maxiter);
        max_outer_iter= std::max( max_outer_iter, (Uint) maxiter);
        total_outer_iter+= maxiter;
        total_damping_iter+= max_damping_steps;
        ++num_outer_iter[maxiter];
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
        localF.set_level_value (loclsp1 (xb));
    else
        localF.set_level_value (0.);

    have_base_point= false;
    base_point ();
    if (!base_in_trust_region_p ()) {
        deformation= MakePoint3D (std::numeric_limits<double>::max (), 0., 0.);
        have_deformation= false;
        return *this;
    }

    if (deformation_method == MAP_LOCAL_LEVEL_SETS)
        deformation= b2w (bxb) - b2w (xb);
    else {
        Point3DCL n (localF.get_locls_grad ()(bxb));
        n/= n.norm ();
        deformation= (loclsp1 (xb) - dh)*n;
    }

    have_deformation= true;
    return *this;
}

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
