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
#include "num/lattice-eval.h"

namespace DROPS
{

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

bool QuaQuaMapperCL::line_search (const Point3DCL& v, const Point3DCL& nx, const TetraCL*& tetra, BaryCoordCL& bary, const TetraSetT& neighborhood) const
{
    const int max_inneriter= 100;
    const int max_damping= 10;
    const double eps= 1e-10;
    const double inner_tol= 5e-9;

    double alpha= 0.,
           dalpha= 0.;
    int inneriter= 0;
    World2BaryCoordCL w2b;
    for (; inneriter < max_inneriter; ++inneriter) {
        // Evaluate ls in v - alpha*nx and check convergence.
        w2b.assign( *tetra);
        bary= w2b( v - alpha*nx);
        const double lsval= ls.val( *tetra, bary);

        if (std::abs( lsval) < inner_tol)
            break;

        // Compute undamped Newton correction dalpha.
        const Point3DCL& gradval= ls_grad_rec.val( *tetra, bary);
        const double slope= inner_prod( gradval, nx);
        if (std::abs( slope) < 1.0e-8)
            std::cout << "g_phi: " << gradval << "\tgy: " << nx << std::endl;
        dalpha= lsval/slope;

        // Apply damping to dalpha until v - (alpha + dalpha)nx is in the neighborhood of tet.
        bary= w2b( v - (alpha + dalpha)*nx);
        if (!is_in_ref_tetra( bary, eps))
            tetra= 0;
        for (int k= 0; tetra == 0 && k < max_damping; ++k, dalpha*= 0.5)
            enclosing_tetra( v - (alpha + dalpha)*nx, neighborhood, eps, tetra, bary);
        if (tetra == 0) {
            std::cout << "v: " << v << "\talpha: " << alpha << "\tnx: " << nx <<  "\tv - alpha*nx: "<< v - alpha*nx << std::endl;
            throw DROPSErrCL("QuaQuaMapperCL::line_search: Coord not in given tetra set.\n");
        }
        alpha+= dalpha;
    }

    if (inneriter >= max_inneriter)
        std::cout <<"QuaQuaMapperCL::line_search: Warning: max inner iteration number at v : " << v << " exceeded; ls.val: " << ls.val( *tetra, bary) << std::endl;

    return inneriter < max_inneriter;
}

void QuaQuaMapperCL::base_point (const TetraCL*& tet, BaryCoordCL& xb) const
// tet and xb specify the point which is projected to the zero level.
// On return tet and xb are the resulting base point.
{
    Point3DCL x, xold; // World coordinates of current and previous xb.
    x= GetWorldCoord( *tet, xb);

    Point3DCL n; // Current search direction.

    int iter;
    for (iter= 0; iter < maxiter_; ++iter) {
        xold= x;
        n=  ls_grad_rec.val( *tet, xb);
        n/= norm( n);
        const bool found_zero_level= line_search( x, n, tet, xb, neighborhoods_[tet]);
        x= GetWorldCoord( *tet, xb);

        if (norm( xold - x) < tol_ && found_zero_level)
            break;
    }
    if (iter >= maxiter_) {
        std::cout << "QuaQuaMapperCL::base_point: max iteration number exceeded; |x - xold|: " << norm( xold - x) << "\tx: " << x << "\t n: " << n << "\t ls(x): " << ls.val( *tet, xb) << std::endl;
    }
}

void QuaQuaMapperCL::jacobian (const TetraCL& tet, const BaryCoordCL& xb, SMatrixCL<3,3>& dph) const
{
    // Compute the basepoint b.
    const TetraCL* btet= &tet;
    BaryCoordCL b= xb;
    base_point( btet, b);

    // Evaluate the quasi normal field in b.
    LocalP2CL<Point3DCL> loc_gh( *btet, ls_grad_rec);
    const Point3DCL gh= loc_gh( b);
    const Point3DCL q_n= gh/gh.norm();

    // Evaluate the normal to the interface in b.
    Point3DCL n;
    LocalP2CL<> locls( *btet, ls);
    LocalP1CL<Point3DCL> gradp2[10];
    SMatrixCL<3,3> T;
    double dummy;
    GetTrafoTr( T, dummy, *btet);
    P2DiscCL::GetGradients( gradp2, gradrefp2, T);
    for (Uint i= 0; i < 10; ++i)
        n+= locls[i]*gradp2[i]( b);
    n/= n.norm();

    // Evaluate the Jacobian of the quasi-normal field in b: dn_h= 1/|G_h| P dG_h.
    SMatrixCL<3,3> dgh;
    for (Uint i= 0; i < 10; ++i)
        dgh+= outer_product( loc_gh[i], gradp2[i]( b));
    const SMatrixCL<3,3> dn= 1/gh.norm()*(dgh - outer_product( q_n, transp_mul( dgh, q_n)));

    // Compute Q.
    const SMatrixCL<3,3> Q= eye<3,3>() - outer_product( q_n/inner_prod( q_n, n), n);

    // Compute d_h(x).
    const Point3DCL x= GetWorldCoord( tet, xb),
                    y= GetWorldCoord( *btet, b);
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
    SMatrixCL<3,3> dph;
    quaqua.jacobian( tet, xb, dph);

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


} // end of namespace DROPS
