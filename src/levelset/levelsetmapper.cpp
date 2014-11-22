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

namespace DROPS
{

void base_point_newton_cacheCL::set_tetra (const TetraCL* newtet)
{
    if (tet == newtet)
        return;

    tet= newtet;

    locls_.assign( *tet, ls_);
    loc_gh_.assign( *tet, ls_grad_rec_);

    SMatrixCL<3,3> T( Uninitialized);
    double dummy;
    GetTrafoTr( T, dummy, *tet);
    P2DiscCL::GetGradients( gradp2_, gradrefp2_, T);
    P2DiscCL::GetHessians( hessp2_, T);

    w2b_.assign( *tet);

    h_= ::cbrt( std::abs( dummy));
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


// Find a tetra in neighborhoods_[tet] enclosing x + d*dx; computes the barycentric coordinates in xb.
// On entry, xb must contain the barycentric coordinates of x + d*dx.
void QuaQuaMapperCL::locate_new_point (const Point3DCL& x, const Point3DCL& dx, const TetraCL*& tet, BaryCoordCL& xb, double& d) const
{
    const double eps= 1e-10;

    if (is_in_ref_tetra( xb, eps))
        return;

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
//         locator_.set_epsilon( eps);
        locator_.locate( x + d*dx);
        tet= &locator_.get_tetra();
        xb= locator_.get_bary();
    }
}

bool QuaQuaMapperCL::line_search (Point3DCL& x, const Point3DCL& nx, const TetraCL*& tetra, BaryCoordCL& bary) const
{
    double dalpha= 0.;
    int inneriter= 0;
    World2BaryCoordCL w2b;
    for (; inneriter < maxinneriter_; ++inneriter) {
        // Evaluate ls in x and check convergence.
        w2b.assign( *tetra);
        bary= w2b( x);
        const double lsval= ls.val( *tetra, bary);

        if (std::abs( lsval) < innertol_)
            break;

        // Compute undamped Newton correction dalpha.
        const Point3DCL& gradval= ls_grad_rec.val( *tetra, bary);
        const double slope= inner_prod( gradval, nx);
        if (std::abs( slope) < 1.0e-8)
            std::cout << "g_phi: " << gradval << "\tgy: " << nx << std::endl;
        dalpha= lsval/slope;

        bary= w2b( x - dalpha*nx);
        double damp= 1.;
        locate_new_point( x, -dalpha*nx, tetra, bary, damp); // Find tetra enclosing x - dalpha*nx; depending on the method, the step must be damped with 0 < damp < 1.
        x-= (damp*dalpha)*nx;
    }

    if (inneriter >= maxinneriter_)
        std::cout <<"QuaQuaMapperCL::line_search: Warning: max inner iteration number at x : " << x << " exceeded; ls.val: " << ls.val( *tetra, bary) << std::endl;
    ++num_inner_iter[inneriter];

    return inneriter < maxinneriter_;
}

double QuaQuaMapperCL::base_point_with_line_search (const TetraCL*& tet, BaryCoordCL& xb) const
// tet and xb specify the point which is projected to the zero level.
// On return tet and xb are the resulting base point.
{
    ScopeTimerCL timer( "QuaQuaMapperCL::base_point_with_line_search");

    Point3DCL x, xold, x0; // World coordinates of current and previous xb and of the inital point.
    x0= x= GetWorldCoord( *tet, xb);

    Point3DCL n; // Current search direction.

    int iter;
    for (iter= 0; iter < maxiter_; ++iter) {
        xold= x;
        n=  ls_grad_rec.val( *tet, xb);
        n/= norm( n);
        const bool found_zero_level= line_search( x, n, tet, xb);

        if (norm( xold - x) < tol_ && found_zero_level)
            break;
    }
    if (iter >= maxiter_) {
        std::cout << "QuaQuaMapperCL::base_point_with_line_search: max iteration number exceeded; |x - xold|: " << norm( xold - x) << "\tx: " << x << "\t n: " << n << "\t ls(x): " << ls.val( *tet, xb) << std::endl;
    }
    ++num_outer_iter[iter];
    return inner_prod( n, x0 - x);
}

double QuaQuaMapperCL::base_point_newton (const TetraCL*& tet, BaryCoordCL& xb) const
// tet and xb specify the point which is projected to the zero level.
// On return tet and xb are the resulting base point.
{
    ScopeTimerCL timer( "QuaQuaMapperCL::base_point_newton");

    bool need_damping= true,
         have_mup= false;

    Point3DCL x0, x, dx; // World coordinates of initial and current xb and the coordinate part of fdx.
    x0= x= GetWorldCoord( *tet, xb);
    SVectorCL<4> F, Fold, // The function of which we search a root, ( x0 - x - s*gh(x), -ls(x) ).
        fdx; // Newton correction.
    double s= 0., ds;  // scaled quasi-distance and increment part of fdx
    double l= 1e-4; // Damping factor for trust region.
    double mup, oldmup= -1.;

    QRDecompCL<4,4> qr;
    SMatrixCL<4,4>& M= qr.GetMatrix();

    cache_.set_tetra( tet);
    const LocalP2CL<>& locls= cache_.locls();
    const LocalP2CL<Point3DCL>& loc_gh= cache_.loc_gh();
    Point3DCL gradp2_xb[10];

    Point3DCL g_ls, // gradient of the level set function.
              gh;   // recovered gradient
    double ghnorm; // norm of the recovered gradient gh.
    SMatrixCL<3,3> dgh; // Jacobian of the recovered gradient.

    int iter;
    for (iter= 0; iter < maxiter_; ++iter) {
        gh= loc_gh( xb);
        ghnorm= gh.norm();

        // Setup F= (x0 - x - s*gh, -locls( xb)); compute Newton update.
        for (Uint i= 0; i < 3; ++i)
            F[i]= x0[i] - x[i] - s*gh[i];
        F[3]= -locls( xb);

        if (F.norm() < tol_)
            break;

        g_ls= Point3DCL();
        for (Uint i= 0; i < 10; ++i) {
            gradp2_xb[i]= cache_.gradp2( i)( xb);
            g_ls+= locls[i]*gradp2_xb[i];
        }

        // Evaluate the Jacobian of gh in xb: dgh.
        dgh= SMatrixCL<3,3>();
        for (Uint i= 0; i < 10; ++i)
            dgh+= outer_product( loc_gh[i], gradp2_xb[i]);

        // Setup the blockmatrix M= (-I - s dgh | - gh, -g_ls^T | 0).
        for (Uint i= 0; i < 3; ++i) {
            for (Uint j= 0; j < 3; ++j) {
                M( i,j)= -s*dgh( i,j);
            }
            M( i,i)-= 1.;
            M( i, 3)= -gh[i];
            M( 3, i)= -g_ls[i];
        }
        M( 3,3)= 0.;
        qr.prepare_solve();

        fdx= F;
        qr.Solve( fdx);
        dx= MakePoint3D( fdx[0], fdx[1], fdx[2]);
        ds= fdx[3];

        if (need_damping && have_mup) {
            const double mu= Fold.norm()/F.norm() * mup;
            l= std::min( 1., mu);
        }
        if (have_mup)
            oldmup=mup;

        for (Uint j= 0; need_damping && j< 10; ++j) {
            if (l < 1e-6) {
                std::cerr << "QuaQuaMapperCL::base_point_newton: Too much damping, giving up. iter: " << iter << " x0: " << x0 << " x: " << x << " dx: " << dx << " ls(x): " << ls.val( *tet, xb) << " l: " << l << " s: " << s << " ghnorm: " << ghnorm << " g_ls: " << g_ls << std::endl;
                l=1e-3;
                need_damping= false;
                break;
            }
            Point3DCL xnew= x - l*dx;
            double snew= s - l*ds;
            cache_.set_tetra( tet);
            BaryCoordCL newxb= cache_.w2b()( xnew);
            const TetraCL* newtet= tet;
            try {
                locate_new_point( x, -dx, newtet, newxb, l); // XXX Do not use neighborhoods (modifies l).
            }
            catch (DROPSErrCL e) {
//                 std::cerr << "Hallo Welt.\n";
                l*= 0.5;
                continue;
            }
            cache_.set_tetra( newtet);

            gh= loc_gh( newxb);
            // Setup Fnew= (x0 - xnew - snew*gh, -locls( newxb)); compute Newton update.
            SVectorCL<4> Fnew;
            for (Uint i= 0; i < 3; ++i)
                Fnew[i]= x0[i] - xnew[i] - snew*gh[i];
            Fnew[3]= -locls( newxb);

            const double theta= Fnew.norm()/F.norm(); // Error monitoring quantities Deuflhard "Newton Methods for Nonlinear Problems" NLEQ-ERR
            mup= 0.5 * F.norm()/(Fnew - (1. - l)*F).norm() * l*l;
            have_mup= true;
// std::cout << "l: " << l << " theta: " << theta << " mup: " << mup << ".\n";
            if (theta >= 1.) { // more damping required.
                l= std::min( 0.5*l, mup);
                continue;
            }
            const double lp= std::min( 1., mup);
            if (l >= 1. && lp >= 1. && theta < 0.7) {
                l= 1.;
                need_damping= false;
                break;
            }
            else if (lp > 4.*l) {
//                 l= lp;
                l*= 4.;
            }
            else {
                break;
            }
        }
        cache_.set_tetra( tet);
        for (Uint j= 0; j < 10; ++j) {
            try {
                const TetraCL* newtet= tet;
                BaryCoordCL newxb= cache_.w2b()( x - l*dx);
                locate_new_point( x, -dx, newtet, newxb, l);
                cache_.set_tetra( newtet);
                xb= newxb;
                break;
            }
            catch (DROPSErrCL e) {
                e.what( std::cerr);
                l*= 0.5;
                if (have_mup)
                    mup=oldmup;
                continue;
            }
        }
        x-= l*dx;
        s-= l*ds;
        Fold= F;
    }
    ++num_outer_iter[iter];
    // Compute the quasi-distance dh:
    gh= loc_gh( xb);
    const double dh= s*gh.norm();
    if (iter >= maxiter_) {
        std::cout << "QuaQuaMapperCL::base_point_newton: max iteration number exceeded; x0: " << x0 << "\tx: " << x << "\t dx: " << dx << "\t s: " << s << "\t ls(x): " << ls.val( *tet, xb) << "\tl: " << l << "\tdh: " << dh << "\tghnorm: " << ghnorm << "\tg_ls: " << g_ls << std::endl;
    }
    return dh;
}

/// Error based Newton method
// double QuaQuaMapperCL::base_point_newton (const TetraCL*& tet, BaryCoordCL& xb) const
// // tet and xb specify the point which is projected to the zero level.
// // On return tet and xb are the resulting base point.
// {
//     ScopeTimerCL timer( "QuaQuaMapperCL::base_point_newton");
// 
//     bool need_damping= true,
//          solution_found= false;
// 
//     Point3DCL x0, x, dx; // World coordinates of initial and current xb and the coordinate part of fdx.
//     x0= x= GetWorldCoord( *tet, xb);
//     SVectorCL<4> F, // The function of which we search a root, ( x0 - x - s*gh(x), -ls(x) ).
//         fdxold, fdx, fdxbar; // (old) Newton correction and simplified Newton correction.
//     double s= 0.;  // scaled quasi-distance
//     double l= 1e-4; // Damping factor for trust region.
// 
//     QRDecompCL<4,4> qr;
//     SMatrixCL<4,4>& M= qr.GetMatrix();
// 
//     cache_.set_tetra( tet);
//     const LocalP2CL<>& locls= cache_.locls();
//     const LocalP2CL<Point3DCL>& loc_gh= cache_.loc_gh();
//     Point3DCL gradp2_xb[10];
// 
//     Point3DCL g_ls, // gradient of the level set function.
//               gh;   // recovered gradient
//     double ghnorm; // norm of the recovered gradient gh.
//     SMatrixCL<3,3> dgh; // Jacobian of the recovered gradient.
// 
//     int iter;
//     for (iter= 0; iter < maxiter_ && !solution_found; ++iter) {
//         gh= loc_gh( xb);
//         ghnorm= gh.norm();
// 
//         g_ls= Point3DCL();
//         for (Uint i= 0; i < 10; ++i) {
//             gradp2_xb[i]= cache_.gradp2( i)( xb);
//             g_ls+= locls[i]*gradp2_xb[i];
//         }
// 
//         // Evaluate the Jacobian of gh in xb: dg_h.
//         dgh= SMatrixCL<3,3>();
//         for (Uint i= 0; i < 10; ++i)
//             dgh+= outer_product( loc_gh[i], gradp2_xb[i]);
// 
//         // Setup the blockmatrix M= (-I - s dgh | - gh, -g_ls^T | 0).
//         for (Uint i= 0; i < 3; ++i) {
//             for (Uint j= 0; j < 3; ++j) {
//                 M( i,j)= -s*dgh( i,j);
//             }
//             M( i,i)-= 1.;
//             M( i, 3)= -gh[i];
//             M( 3, i)= -g_ls[i];
//         }
//         M( 3,3)= 0.;
//         qr.prepare_solve();
// 
//         // Setup F= (x0 - x - s*gh, -locls( xb)); compute Newton update.
//         for (Uint i= 0; i < 3; ++i)
//             F[i]= x0[i] - x[i] - s*gh[i];
//         F[3]= -locls( xb);
//         fdx= F;
//         qr.Solve( fdx);
//         dx= MakePoint3D( fdx[0], fdx[1], fdx[2]);
// 
//         if (fdx.norm() < tol_) { // Solution found.
//             double dummy= 1.;
//             locate_new_point( x, -dx, tet, xb, dummy);
//             s-= fdx[3];
//             cache_.set_tetra( tet);
//             break;
//         }
//         if (need_damping && iter > 0) {
//             const double mu= l * fdxold.norm()/(fdxbar - fdxold).norm() * fdxbar.norm()/fdx.norm();
//             l= std::min( 1., mu);
//         }
// 
//         while (need_damping && !solution_found) {
//             if (l < 1e-6) {
//                 std::cerr << "QuaQuaMapperCL::base_point_newton: Too much damping, giving up. iter: " << iter << " x0: " << x0 << " x: " << x << " dx: " << dx << " ls(x): " << ls.val( *tet, xb) << " l: " << l << " s: " << s << " ghnorm: " << ghnorm << " g_ls: " << g_ls << std::endl;
// l=1e-4;
// need_damping= false;
//                 break;
//             }
//             cache_.set_tetra( tet);
//             BaryCoordCL newxb= cache_.w2b()( x - l*dx);
//             const TetraCL* newtet= tet;
//             locate_new_point( x, -dx, newtet, newxb, l); // XXX Do not use neighborhoods (modifies l).
//             Point3DCL xnew= x - l*dx;
//             double snew= s - l*fdx[3];
//             cache_.set_tetra( newtet);
// 
//             gh= loc_gh( newxb);
//             // Setup F= (x0 - xnew - snew*gh, -locls( newxb)); compute Newton update.
//             for (Uint i= 0; i < 3; ++i)
//                 F[i]= x0[i] - xnew[i] - snew*gh[i];
//             F[3]= -locls( newxb);
//             fdxbar= F;
//             qr.Solve( fdxbar);
// 
//             const double theta= fdxbar.norm()/fdx.norm(); // Error monitoring quantities Deuflhard "Newton Methods for Nonlinear Problems" NLEQ-ERR
//             const double mup= 0.5*fdx.norm()*l*l/(fdxbar - (1. - l)*fdx).norm();
// std::cout << "theta: " << theta << " mup: " << mup << ".\n";
//             if (theta >= 1.) { // more damping required.
//                 l= std::min( 0.5*l, mup);
//                 continue;
//             }
//             const double lp= std::min( 1., mup);
//             if (l > 0.999 && lp > 0.999) {
//                 if (fdxbar.norm() < tol_) {
//                     l= 1.;
//                     fdx= fdxbar;
//                     dx= MakePoint3D( fdx[0], fdx[1], fdx[2]);
//                     solution_found= true;
//                     locate_new_point( x, -dx, tet, xb, l);
//                     break;
//                 }
//                 if (theta < 0.5) {
//                     l= 1.;
//                     need_damping= false;
//                 }
//             }
//             else if (lp > 4.*l) {
//                 l= lp;
//             }
//             else {
//                 break;
//             }
//         }
//         cache_.set_tetra( tet);
//         xb= cache_.w2b()( x - l*dx);
//         locate_new_point( x, -dx, tet, xb, l);
//         x-= l*dx;
//         s-= l*fdx[3];
//         cache_.set_tetra( tet);
//         fdxold= fdx;
//     }
//     ++num_outer_iter[iter];
//     // Compute the quasi-distance dh:
//     gh= loc_gh( xb);// The recovered gradient.
//     const double dh= s*gh.norm();
//     if (iter >= maxiter_) {
//         std::cout << "QuaQuaMapperCL::base_point_newton: max iteration number exceeded; x0: " << x0 << "\tx: " << x << "\t dx: " << dx << "\t s: " << s << "\t ls(x): " << ls.val( *tet, xb) << "\tl: " << l << "\tdh: " << dh << "\tghnorm: " << ghnorm << "\tg_ls: " << g_ls << std::endl;
//     }
// // std::cout << "F: " << F << std::endl;
// // char c;
// // std::cin >> c;
// // if (c == 'q')
// //     exit( 0);
// //     else {
// //         std::cout << "QuaQuaMapperCL::base_point_newton: Ok! x0: " << x0 << "\tx: " << x << "\t dx: " << dx << "\t ls(x): " << ls.val( *tet, xb) << "\tdamping: " << damp << "\tdh: " << dh << "\titer: " << iter << std::endl;
// //     }
//     return dh;
// }

double QuaQuaMapperCL::base_point (const TetraCL*& tet, BaryCoordCL& xb) const
// tet and xb specify the point which is projected to the zero level.
// On return tet and xb are the resulting base point.
{
    return use_line_search_ == true
           ? base_point_with_line_search( tet, xb)
           : base_point_newton( tet, xb);
}

void QuaQuaMapperCL::jacobian (const TetraCL& tet, const BaryCoordCL& xb, const TetraCL& btet, const BaryCoordCL& b, SMatrixCL<3,3>& dph) const
{
    cache_.set_tetra( &btet);

    // Evaluate the quasi normal field in b.
    const LocalP2CL<Point3DCL>& loc_gh= cache_.loc_gh();
    const Point3DCL gh= loc_gh( b);
    const Point3DCL q_n= gh/gh.norm();

    // Evaluate the normal to the interface in b.
    Point3DCL n;
    const LocalP2CL<>& locls= cache_.locls();
    Point3DCL gradp2_b[10];
    for (Uint i= 0; i < 10; ++i) {
        gradp2_b[i]= cache_.gradp2( i)( b);
        n+= locls[i]*gradp2_b[i];
    }
    n/= n.norm();

    // Evaluate the Jacobian of the quasi-normal field in b: dn_h= 1/|G_h| P dG_h.
    SMatrixCL<3,3> dgh;
    for (Uint i= 0; i < 10; ++i)
        dgh+= outer_product( loc_gh[i], gradp2_b[i]);
    const SMatrixCL<3,3> dn= 1./gh.norm()*(dgh - outer_product( q_n, transp_mul( dgh, q_n)));

    // Compute Q.
    const SMatrixCL<3,3> Q= eye<3,3>() - outer_product( q_n/inner_prod( q_n, n), n);

    // Compute d_h(x).
    const Point3DCL x= GetWorldCoord( tet, xb),
                    y= GetWorldCoord( btet, b);
    const double dhx= inner_prod( q_n, x - y);

    // Compute the Jacobian of p_h.
    QRDecompCL<3,3> qr;
    qr.GetMatrix()= eye<3,3>() + dhx*Q*dn;
    qr.prepare_solve();
    Point3DCL tmp;
    for (Uint i= 0; i < 3; ++i) {
        tmp= Q.col( i);
        qr.Solve( tmp);
        dph.col( i, tmp);
    }
}

void QuaQuaMapperCL::jacobian (const TetraCL& tet, const BaryCoordCL& xb, SMatrixCL<3,3>& dph) const
{
    // Compute the basepoint b.
    const TetraCL* btet= &tet;
    BaryCoordCL b= xb;
    base_point( btet, b);

    jacobian( tet, xb, *btet, b, dph);
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

double abs_det (const TetraCL& tet, const BaryCoordCL& xb, const TetraCL& btet, const BaryCoordCL& b, const QuaQuaMapperCL& quaqua, const SurfacePatchCL& p)
{
    if (p.empty())
        return 0.;

    // Compute the jacobian of p_h.
    SMatrixCL<3,3> dph;
    quaqua.jacobian( tet, xb, btet, b, dph);

    const Bary2WorldCoordCL b2w( tet);
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

double abs_det (const TetraCL& tet, const BaryCoordCL& xb, const QuaQuaMapperCL& quaqua, const SurfacePatchCL& p)
{
    if (p.empty())
        return 0.;

    // Compute the basepoint b.
    const TetraCL* btet= &tet;
    BaryCoordCL b= xb;
    quaqua.base_point( btet, b);

    return abs_det( tet, xb, *btet, b, quaqua, p);
}

} // end of namespace DROPS
