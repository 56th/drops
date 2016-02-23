/// \file newton.h
/// \brief The damped Newton method as class template.
/// \author Joerg Grande

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
 * Copyright 2016 Joerg Grande, Aachen, Germany
*/


#include "misc/problem.h"

namespace DROPS
{

// class QuaQuaBasePointFunctionCL
// {
//   public:
//     typedef SVectorCL<4>   value_type;
//     typedef SMatrixCL<4,4> derivative_type;
// 
//     bool found_newtet;
//     btet= tet;
//     bxb= xb;
// 
//     double    snew;
//     BaryCoordCL bxbnew;
//     const TetraCL* newtet;
// 
//         // Evaluate the gradient of locls in bxb: g_ls.
//         g_ls= Point3DCL();
//         for (Uint i= 0; i < 10; ++i) {
//             gradp2_xb[i]= cache_.gradp2( i)( bxb);
//             g_ls+= locls[i]*gradp2_xb[i];
//         }
// 
//         // Evaluate the Jacobian of gh in bxb: dgh.
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
//         Msave= M;
//         qr.prepare_solve();
// 
// 
//         dx= MakePoint3D( fdx[0], fdx[1], fdx[2]);
//         ds= fdx[3];
// 
//         l= std::min( 1., 0.5*cache_.get_h()/dx.norm());
// 
// };
// 
//     double s= 0., // scaled quasi-distance
//            ds,    // increment part of fdx
// 
//     Point3DCL x0, x, dx; // World coordinates of initial and current bxb and the coordinate part of fdx.
//     x0= x= GetWorldCoord( *btet, bxb);
// 
// Fnew.norm() < F.norm() + armijo_c_*inner_prod( transp_mul( Msave, F), fdx)/F.norm()*l)
// 
//         if (!found_newtet)
//             throw DROPSErrCL( "QuaQuaMapperCL::base_point_newton: Could not find the right tetra.\n");


namespace NewtonImplNS {

template <typename T>
struct dotCL
{};

template <>
struct dotCL<double>
{
    static double dot (double x, double y) { return x*y; }
};

template <Uint Size>
struct dotCL< SVectorCL<Size> >
{
    static double dot (const SVectorCL<Size>& x, const SVectorCL<Size>& y) { return inner_prod (x, y); }
};

template <typename T>
double dot (T x, T y)
{
    return dotCL<T>::dot (x, y);
}

} // end of namespace NewtonImplNS


/// FunctionT must define: * type_def for value_type,
///     * value (x),
///     * apply_derivative (x, v),
///     * apply_derivative_inverse (x, v),
///     * initial_damping_factor (x, dx, F) (must return a nonnegative value which is further limited to 1 in the algorithm).
template <typename FunctionT>
void newton_solve (FunctionT& fun, typename FunctionT::value_type& x, size_t& maxiter, double& tol, size_t& max_damping_steps, double armijo_c, std::ostream* os= 0)
{
    const double min_step_length= 1e-7; // Minimum value for step-length factor l.

    typedef typename FunctionT::value_type value_type;
    value_type dx, // Newton correction.
               F,  // The function of which we search a root.
               Fnew; // Helper in the step-size computation.
    double l; // Damping factor for line search.
    size_t total_damping_iter= 0; // Count all rejected steps sizes.
    double normF;

    int iter;
    for (iter= 0; true; ++iter) {
        F= fun.value (x);
        normF= std::sqrt (NewtonImplNS::dot (F, F));
        if (iter >= maxiter || normF < tol)
            break;
        dx= fun.apply_derivative_inverse (x, F);
//         if (iter == 0) // Initial gradient descent step
//             dx= fun.apply_derivative_transpose (x, F/normF);
        const double armijo_slope= armijo_c*inner_prod( F, fun.apply_derivative (x, dx))/normF;
        // Armijo-rule with backtracking
        l= std::min( 1., fun.initial_damping_factor (x, dx, F));
        Uint j;
        for (j= 0; j < max_damping_steps && l >= min_step_length; ++j, l*= 0.5) {
            Fnew= fun.value (x - l*dx);
            if (std::sqrt (NewtonImplNS::dot (Fnew, Fnew)) < normF + armijo_slope*l)
                break;
        }
        total_damping_iter+= j;
        if (l < min_step_length) {
            std::cerr << "newton_solve: Too much damping. iter: " << iter << " x: " << x << " dx: " << dx << " l: " << l << " F: " << F << std::endl;
            break;
        }
        if (os)
            (*os) << "iter: " << iter << " x: " << x << " dx: " << dx << " l: " << l << " F: " << F << std::endl;
        x-= l*dx;
    }
    if (iter >= maxiter) {
        std::cout << "newton_solve: max iteration number exceeded; \tx: " << x << "\t dx: " << dx << "\tl: " << l << "\t F: " << F << std::endl;
    }
    maxiter= iter;
    tol= normF;
    max_damping_steps= total_damping_iter;
}

} // end of namespace DROPS
