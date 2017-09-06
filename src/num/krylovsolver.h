/// \file krylovsolver.h
/// \brief iterative solvers
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Helmut Jarausch, Volker Reichelt; SC RWTH Aachen:

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
 * Copyright 2012 LNM/SC RWTH Aachen, Germany
*/

#ifndef DROPS_SOLVER_H
#define DROPS_SOLVER_H

#include <vector>
#include "misc/container.h"
#include "num/spmat.h"
#include "num/spblockmat.h"

namespace DROPS
{

//*****************************************************************************
//
//  Iterative solvers: CG, PCG, PCGNE, GMRES, PMINRES, MINRES, BiCGStab, GCR,
//                     GMRESR, IDR(s), QMR
//
//*****************************************************************************

//=============================================================================
//  Implementation of the methods
//=============================================================================

/// \brief Parallel CG-Algorithm
template <typename Mat, typename Vec, typename ExCL>
bool CG(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, int& max_iter,
        double& tol, bool measure_relative_tol, std::ostream* output=0)
    /// \param[in]     A        local distributed coefficients-matrix of the linear equation system
    /// \param[in,out] x_acc    start vector and the solution in accumulated form
    /// \param[in]     b        rhs of the linear equation system (distributed form)
    /// \param[in]     ExX      ExchangeCL corresponding to the RowIdx of x and the ColIdx of A
    /// \param[in,out] max_iter IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol      IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol  true: stop if |b - Ax|/|b| <= tol, false: stop if |b - Ax| <= tol.
    /// \return                 convergence within max_iter iterations
{
    Vec r( A*x_acc - b ), r_acc(r);
    double normb= ExX.Norm( b, false),
           res,
           resid=ExX.Norm_sq( r, false, &r_acc);
    Vec d_acc( -r_acc);

    if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;

    if ((res= std::sqrt( resid)/normb) <= tol)
    {
        tol= res;
        max_iter= 0;
        return true;
    }

    for (int i= 1; i <= max_iter; ++i)
    {
        const Vec    Ad= A*d_acc;
        const double delta= ExX.ParDot( Ad, false, d_acc, true);
        const double alpha= resid/delta;
        double       beta= resid;

        axpy(alpha, d_acc, x_acc);  // x+= alpha*d;
        axpy(alpha, Ad, r);         // r+= alpha*Ad;

        resid= ExX.Norm_sq( r, false, &r_acc);

        if ( output){
            (*output) << "CG: " << i << " resid " << std::sqrt(resid) << std::endl;
        }

        if ((res= std::sqrt( resid)/normb) <= tol)
        {
            tol= res;
            max_iter= i;
            return true;
        }
        beta= resid / beta;
        d_acc= beta*d_acc-r_acc;
    }
    tol= res;
    return false;
}


/// \brief Parallel preconditioned CG-Algorithm
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool PCG(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX,
            PreCon& M, int& max_iter, double& tol, bool measure_relative_tol=false,
            std::ostream* output=0)
    /// \param[in]     A                    local distributed coefficients-matrix of the linear equation system
    /// \param[in,out] x_acc                start vector and the solution in accumulated form
    /// \param[in]     b                    rhs of the linear equation system (distributed form)
    /// \param[in]     ExX                  ExchangeCL corresponding to the RowIdx of x and the ColIdx of A
    /// \param[in,out] M                    Preconditioner
    /// \param[in,out] max_iter             IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol                  IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol measure relative residual
    /// \return                             convergence within max_iter iterations
{
    // Check if preconditioner needs diagonal of matrix. The preconditioner
    // only computes the diagonal new, if the matrix has changed
    if (M.NeedDiag())
        M.SetDiag(A, ExX);

    const size_t n= x_acc.size();
    Vec p_acc(n), z(n), z_acc(n), q(n), q_acc(n), r( b - A*x_acc), r_acc(r);

    double rho,
           rho_1= 0.0,
           resid= ExX.Norm( r, false, &r_acc),
           normb= ExX.Norm( b, false);

    if (normb == 0.0 || measure_relative_tol == false)
        normb= 1.0;
    resid = resid/normb;


    if (resid<=tol){
        tol= resid;
        max_iter= 0;
        return true;
    }

    M.Apply(A, z, r, ExX);
    rho = ExX.ParDot( z, M.RetAcc(), r_acc, true, &z_acc);
    p_acc= z_acc;

    for (int i=1; i<=max_iter; ++i)
    {
        q= A*p_acc;

        const double lambda = ExX.ParDot( q, false, p_acc, true, &q_acc);
        const double alpha  = rho/lambda;

        axpy( alpha, p_acc, x_acc);        // x+= alpha*p;
        axpy( -alpha, q, r);               // r-= alpha*q;
        axpy( -alpha, q_acc, r_acc);       // r-= alpha*q;

        resid= ExX.Norm( r_acc, true) / normb;

        if ( output){
            (*output) << "PCG: " << i << " resid " << resid << std::endl;
        }

        if (resid<=tol){
            tol= resid;
            max_iter= i;
            return true;
        }

        M.Apply(A, z, r, ExX);
        rho_1= rho;
        rho = ExX.ParDot( z, M.RetAcc(), r_acc, true, &z_acc);
        z_xpay( p_acc, z_acc, rho/rho_1, p_acc); // p= z + (rho/rho_1)*p;
    }
    tol= resid;
    return false;
}


/// \brief PCGNE: Preconditioned CG for the normal equations (error-minimization)
///
/// Solve A*A^T x = b with left preconditioner M. This is more stable than PCG with
/// a CompositeMatrixCL.
///
/// The return value indicates convergence within max_iter (input)
/// iterations (true), or no convergence within max_iter iterations (false).
/// Upon successful return, output arguments have the following values:
///
/// \param A - matrix (not necessarily quadratic)
/// \param b - right hand side
/// \param u - approximate solution to A*A^T u = b
/// \param M - preconditioner
/// \param max_iter - number of iterations performed before tolerance was reached
/// \param tol - 2-norm of the (relative, see below) residual after the final iteration
/// \param measure_relative_tol - If true, stop if |b - A*A^T u|/|b| <= tol,
///        if false, stop if |b - A*A^T u| <= tol.
template <typename Mat, typename Vec, typename PreCon, typename ExACL, typename ExATranspCL>
bool
PCGNE(const Mat& A, Vec& u, const Vec& b, const ExACL& ExAX,  const ExATranspCL& ExATranspX, PreCon& M,
    int& max_iter, double& tol, bool measure_relative_tol=false, std::ostream* output=0)
{
    M.SetDiag(A, ExATranspX, ExAX);
    Vec Atranspu( transp_mul( A, u));
    ExAX.Accumulate( Atranspu);
    Vec r( b - A*Atranspu);

    double normb= ExATranspX.Norm( b, false);
    if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;

    double resid= ExATranspX.Norm( r, false)/normb;
    if (output)
        (*output) << "PCGNE: iter: 0 resid: " << resid <<'\n';
    if (resid <= tol) {
        tol= resid;
        max_iter= 0;
        return true;
    }

    const size_t n= A.num_rows();
    const size_t num_cols= A.num_cols();

    Vec z( n);
    M.Apply( A, z, r, ExATranspX, ExAX);
    Vec pt( num_cols);
    Vec qtacc( z), ptacc( pt), racc(r), zacc(z);


    double rho= ExATranspX.ParDot( z, M.RetAcc(), r, false, &qtacc), rho_1;

    for (int i= 1; i <= max_iter; ++i) {
        pt= transp_mul( A, qtacc);
        const double alpha= rho/ExAX.Norm_sq( pt, false, &ptacc);
        u+= alpha*qtacc;
        r-= alpha*(A*ptacc);
        M.Apply( A, z, r, ExATranspX, ExAX);

        resid= ExATranspX.Norm( r, false, &racc)/normb;
        if ( output && i%10 == 0) (*output) << "PCGNE: iter: " << i << " resid: " << resid <<'\n';
        if (resid <= tol) {
            tol= resid;
            max_iter= i;
            return true;
        }
        rho_1= rho;
        rho= ExATranspX.ParDot( z, M.RetAcc(), racc, true, &zacc);
        qtacc= zacc + (rho/rho_1)*qtacc;
    }
    tol= resid;
    return false;
}

//-----------------------------------------------------------------------------
// GMRES:
//-----------------------------------------------------------------------------

inline void GMRES_GeneratePlaneRotation(double &dx, double &dy, double &cs, double &sn)
{
    if (dy == 0.0)
    {
        cs = 1.0;
        sn = 0.0;
    }
    else
    {
        const double r=std::sqrt(dx*dx+dy*dy);
        cs = dx/r;
        sn = dy/r;
    }
}


inline void GMRES_ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn)
{
    const double tmp = cs*dx + sn*dy;
    dy = -sn*dx + cs*dy;
    dx = tmp;
}


template <typename Mat, typename Vec>
void GMRES_Update(Vec &x, int k, const Mat &H, const Vec &s, const std::vector<Vec> &v)
{
    Vec y(s);

    // Backsolve:
    for ( int i=k; i>=0; --i )
    {
        y[i] /= H(i,i);
        for ( int j=i-1; j>=0; --j )
            y[j] -= H(j,i) * y[i];
    }

    for ( int i=0; i<=k; ++i )
        x += y[i] * v[i];
}
enum PreMethGMRES { RightPreconditioning, LeftPreconditioning};

#ifndef _PAR
//-----------------------------------------------------------------------------
// GMRES: The return value indicates convergence within max_iter (input)
// iterations (true), or no convergence within max_iter iterations (false).
//
// Upon successful return, output arguments have the following values:
//
// x - approximate solution to Ax = b
// max_iter - number of iterations performed before tolerance was reached
// tol - 2-norm of the (relative, see below) preconditioned residual after the
//     final iteration
//
// measure_relative_tol - If true, stop if |M^(-1)( b - Ax)|/|M^(-1)b| <= tol,
//     if false, stop if |M^(-1)( b - Ax)| <= tol.
// calculate2norm - If true, the unpreconditioned (absolute) residual is
//     calculated in every step for debugging. This is rather expensive.
// method - If RightPreconditioning, solve the right-preconditioned GMRES: A*M^(-1)*u = b, u=Mx,
//     if LeftPreconditioning, solve the left-preconditioned GMRES:   M^(-1)*A*x= M^(-1)*b,
//TODO: bug in RightPreconditioning?
//-----------------------------------------------------------------------------
template <typename Mat, typename Vec, typename PreCon, typename ExT>
bool
GMRES(const Mat& A, Vec& x, const Vec& b, const ExT& ex, PreCon& M,
      int /*restart parameter*/ m, int& max_iter, double& tol,
      bool measure_relative_tol= true, bool calculate2norm= false, PreMethGMRES method = LeftPreconditioning)
{
    m= (m <= max_iter) ? m : max_iter; // m > max_iter only wastes memory.

    DMatrixCL<double> H( m, m);
    Vec               s( m), cs( m), sn( m), w( b.size()), r( b.size());
    std::vector<Vec>  v( m);
    double            beta, normb, resid;
    Vec z(x.size()), t(x.size());
    for (int i= 0; i < m; ++i)
        v[i].resize( b.size());

     if (method == RightPreconditioning)
     {
          r= b - A*x;
          beta= norm( r);
          normb= norm(b);
     }
     else
     {
          M.Apply( A, r, Vec( b - A*x), ex);
          beta= norm( r);
          M.Apply( A, w, b, ex);
          normb= norm( w);
     }
    if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;

    resid = beta/normb;
    if (resid <= tol) {
        tol= resid;
        max_iter= 0;
        return true;
    }

    int j= 1;
    while (j <= max_iter) {
        v[0]= r*(1.0/beta);
        s= 0.0;
        s[0]= beta;

        int i;
        for (i= 0; i < m - 1 && j <= max_iter; ++i, ++j) {
            if (method == RightPreconditioning)
            {
                M.Apply( A, w, v[i], ex);
                w=A*w;
            }
            else M.Apply( A, w, A*v[i], ex);
            for (int k= 0; k <= i; ++k ) {
                H( k, i)= dot( w, v[k]);
                w-= H( k, i)*v[k];
            }

            H( i + 1, i)= norm( w);
            v[i + 1]= w*(1.0/H( i + 1, i));

            for (int k= 0; k < i; ++k)
                GMRES_ApplyPlaneRotation( H(k,i), H(k + 1, i), cs[k], sn[k]);

            GMRES_GeneratePlaneRotation( H(i,i), H(i+1,i), cs[i], sn[i]);
            GMRES_ApplyPlaneRotation( H(i,i), H(i+1,i), cs[i], sn[i]);
            GMRES_ApplyPlaneRotation( s[i], s[i+1], cs[i], sn[i]);

            resid= std::abs( s[i+1])/normb;
            if (calculate2norm == true) { // debugging aid
                Vec y( x);
                if (method == RightPreconditioning)
                {
                    z=0.;
                    GMRES_Update( z, i, H, s, v);
                    M.Apply( A, t, z, ex);
                    y+=t;
                }
                else GMRES_Update( y, i, H, s, v);
                double resid2= norm( Vec( b - A*y));
                std::cout << "GMRES: absolute residual 2-norm: " << resid2
                          << "\tabsolute preconditioned residual 2-norm: "
                          << std::fabs( s[i+1]) << '\n';
            }
            if (resid <= tol) {
                if (method == RightPreconditioning)
                {
                    z=0.;
                    GMRES_Update( z, i, H, s, v);
                    M.Apply( A, t, z, ex);
                    x+=t;
                }
                else GMRES_Update( x, i, H, s, v);
                tol= resid;
                max_iter= j;
                return true;
            }
        }

        if (method == RightPreconditioning)
        {
            z=0.;
            GMRES_Update( z, i - 1, H, s, v);
            M.Apply( A, t, z, ex);
            x+=t;
            r= b - A*x;
        }
        else
        {
            GMRES_Update( x, i - 1, H, s, v);
            M.Apply( A, r, Vec( b - A*x), ex);
        }
        beta=norm(r);
        resid= beta/normb;
        if (resid <= tol) {
            tol= resid;
            max_iter= j;
            return true;
        }
    }
    tol= resid;
    return false;
}

#else
/// \brief Computes an orthogonal vector on i vectors by the standard Gram-Schmidt method
template <typename Vec, typename ExCL>
void StandardGramSchmidt(DMatrixCL<double>& H,
                          Vec& w, bool acc_w, const std::vector<Vec>& v, bool acc_v,
                          int i, const ExCL& ex, Vec& tmpHCol)
/// This routine uses only one synchronization point.
/// \param H       Hessenberg matrix
/// \param w       orthogonal vector
/// \param acc_w   is w accumulated
/// \param v       vectors used to orthogonalize w
/// \param acc_v   is v accumulated
/// \param i       number of vectors v
/// \param ex      class to accumulate a vector
/// \param tmpHCol temporary vector to store a column of H (as parameter so no mem is allocated in this routine)
{
    if (acc_v&&acc_w){                          // both vectors are accumulated
        for (int k=0; k<=i; ++k)
            tmpHCol[k] = ex.LocalDot( w, true, v[k], true);
    }
    else if (!acc_v&&!acc_w){                   // one of the both vectors have to be accumulated: really bad!
        for (int k=0; k<=i; ++k)
            tmpHCol[k] = ex.LocalDot( w, false, v[k], false);
    }
    else        // update of w do only works on same types
        throw DROPSErrCL("StandardGramSchmidt: Cannot do Gram Schmidt on that kind of vectors!");

    // Syncpoint!
    ProcCL::GlobalSum(Addr(tmpHCol), H.GetCol(i), i+1);
    for (int k=0; k<=i; ++k)
        w -= H(k,i)*v[k];
}


/// \brief Computes an orthogonal vector on i vectors by the modified Gram-Schmidt method
template <typename Vec, typename ExCL>
void ModifiedGramSchmidt(DMatrixCL<double>& H, Vec& w, bool acc_w, const std::vector<Vec>& v, bool acc_v, int i, const ExCL& ex)
/// \param H     Hessenberg matrix
/// \param w     orthogonal vector
/// \param acc_w is w accumulated
/// \param v     vectors used to orthogonalize w
/// \param acc_v is v accumulated
/// \param i     number of vectors v
/// \param ex    class to accumulate a vector
{
    if (acc_v&&acc_w){
        for (int k=0; k<=i; ++k){
            H( k, i)= ex.ParDot( w, true, v[k], true);
            w-= H( k, i)*v[k];
        }
    }
    else if (!acc_v&&!acc_w){
        for (int k=0; k<=i; ++k){
            H( k, i)= ex.ParDot( w, false, v[k], false);
            w-= H( k, i)*v[k];
        }
    }
    else
        throw DROPSErrCL("StandardGramSchmidt: Cannot do Gram Schmidt on that kind of vectors!");
}

/// \brief Parallel GMRES-method.
///
/// This method is the same algorithm as the serial algorithm. For performance issues take the ParModGMRES
/// procedure.
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
  bool GMRES(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M,
                int m, int& max_iter, double& tol,
                bool measure_relative_tol=true, PreMethGMRES method=LeftPreconditioning, std::ostream* output=0)
    /// \param[in]     A                    local distributed coefficients-matrix of the linear equation system
    /// \param[in,out] x_acc                start vector and the solution in accumulated form
    /// \param[in]     b                    rhs of the linear equation system (distributed form)
    /// \param[in]     ExX                  ExchangeCL corresponding to the RowIdx of x and the ColIdx of A
    /// \param[in,out] M                    Preconditioner
    /// \param[in]     m                    number of steps after a restart is performed
    /// \param[in,out] max_iter             IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol                  IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol if true stop if |M^(-1)(b-Ax)|/|M^(-1)b| <= tol, else stop if |M^(-1)(b-Ax)|<=tol
    /// \param[in]     method               left or right preconditioning (see solver.h for definition and declaration)
    /// \return  convergence within max_iter iterations
    /// \pre     the preconditioner should be able to handle an accumulated b
{
    // Check if preconditioner needs diagonal of matrix. The preconditioner
    // only computes the diagonal new, if the matrix has changed
    if (M.NeedDiag())
        M.SetDiag(A, ExX);

    m= (m <= max_iter) ? m : max_iter; // m > max_iter only wastes memory.

    // The preconditioner can have one of the following properties by calling M.Apply(A,x,b) with b distributed:
    // 1. x has distributed form    (like Dummy, Jacobi, BlockSSOR)
    // 2. x has accumulated form    (like SolverAsPreCL)
    // Unfortunately there are a lot of differences, so we have to implement two strategies

    if (!M.RetAcc()){                                   // case 1
        if (method==LeftPreconditioning)
            throw DROPSErrCL("ParGMRES: Left preconditioning does not work!");
        /// \todo (of) Left preconditioning funktioniert nicht!
        DMatrixCL<double> H( m, m);
        Vec               s( m), cs( m), sn( m),
                          w( b.size()), r( b.size());
        std::vector<Vec>  v( m);
        double            beta, normb, resid;
        Vec z(x_acc.size()), t(x_acc.size());
        for (int i= 0; i < m; ++i)
            v[i].resize( b.size());

        if (method == RightPreconditioning)
        {
            r= b - A*x_acc;
            beta = ExX.Norm(r, false);
            normb= ExX.Norm(b, false);
        }
        else
        {
            M.Apply( A, r, Vec( b - A*x_acc), ExX);
            beta= ExX.Norm(r, false);
            M.Apply( A, w, b, ExX);
            normb= ExX.Norm(w, false);
        }
        if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;

        resid = beta/normb;
        if (resid <= tol) {
            tol= resid;
            max_iter= 0;
            return true;
        }

        int j= 1;
        while (j <= max_iter) {
            v[0]= r*(1.0/beta);
            s= 0.0;
            s[0]= beta;

            int i;
            for (i= 0; i<m-1 && j<=max_iter; ++i, ++j) {
                if (method == RightPreconditioning)
                {
                    M.Apply( A, w, v[i], ExX);
                    w=A*ExX.GetAccumulate(w);
                }
                else
                    M.Apply( A, w, A*ExX.GetAccumulate(v[i]), ExX);
                for (int k= 0; k <= i; ++k ) {
                    H( k, i)= ExX.ParDot( w, false, v[k], false);
                    w-= H( k, i)*v[k];
                }

                if (i == m - 1) break;

                H( i + 1, i)= ExX.Norm(w, false);
                v[i + 1]= w*(1.0/H( i + 1, i));

                for (int k= 0; k < i; ++k)
                    GMRES_ApplyPlaneRotation( H(k,i), H(k + 1, i), cs[k], sn[k]);

                GMRES_GeneratePlaneRotation( H(i,i), H(i+1,i), cs[i], sn[i]);
                GMRES_ApplyPlaneRotation( H(i,i), H(i+1,i), cs[i], sn[i]);
                GMRES_ApplyPlaneRotation( s[i], s[i+1], cs[i], sn[i]);

                resid= std::abs( s[i+1])/normb;
                if (output && j%10 == 0)
                    (*output) << "ParGMRES: " << j << " resid " << resid << std::endl;
                if (resid <= tol) {
                    if (method == RightPreconditioning)
                    {
                        z=0.;
                        GMRES_Update( z, i, H, s, v);
                        M.Apply( A, t, z, ExX);
                        x_acc+=ExX.GetAccumulate(t);
                    }
                    else{
                        z=0.;
                        GMRES_Update( z, i, H, s, ExX.GetAccumulate(v));
                        x_acc += z;
                    }
                    tol= resid;
                    max_iter= j;
                    return true;
                }
            }

            if (method == RightPreconditioning)
            {
                z=0.;
                GMRES_Update( z, i-1, H, s, v);
                M.Apply( A, t, z, ExX);
                x_acc+=ExX.GetAccumulate(t);
                r= b - A*x_acc;
            }
            else
            {
                GMRES_Update( x_acc, i-1, H, s, ExX.GetAccumulate(v));
                M.Apply( A, r, Vec( b - A*x_acc), ExX);
            }
            beta=ExX.Norm(r, false);
            resid= beta/normb;
            if (resid <= tol) {
                tol= resid;
                max_iter= j;
                return true;
            }
        }
        tol= resid;
        return false;
    }
    else                                                   // case 2: Apply returns accumulated vector
    {
        DMatrixCL<double> H( m, m);
        Vec               s( m), cs( m), sn( m),
                          w( b.size()), r( b.size());
        std::vector<Vec>  v( m);
        double            beta, normb, resid;
        Vec z(x_acc.size()), t(x_acc.size());
        for (int i= 0; i < m; ++i)
            v[i].resize( b.size());

        if (method == RightPreconditioning)
        {
            Vec r_tmp(b - A*x_acc);
            beta = ExX.Norm(r_tmp, false, &r);
            normb= ExX.Norm(b, false);
        }
        else
        {
            M.Apply( A, r, Vec( b - A*x_acc), ExX);
            beta = ExX.Norm(r, true);
            M.Apply( A, w, b, ExX);
            normb= ExX.Norm(w, true);
        }
        if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;

        resid = beta/normb;
        if (resid <= tol) {
            tol= resid;
            max_iter= 0;
            return true;
        }

        int j= 1;
        while (j <= max_iter) {
            v[0]= r*(1.0/beta);                             // v has accumulated form because r is accumulated
            s= 0.0;
            s[0]= beta;

            int i;
            for (i= 0; i<m-1 && j<=max_iter; ++i, ++j) {
                if (method == RightPreconditioning)
                {
                    M.Apply( A, w, v[i], ExX);                   // hopefully, preconditioner do right things with accumulated v[i]
                    w=A*w;
                    ExX.Accumulate(w);
                }
                else
                    M.Apply( A, w, A*v[i], ExX);
                for (int k= 0; k <= i; ++k ) {
                    H( k, i)= ExX.ParDot(w, true, v[k], true);
                    w-= H( k, i)*v[k];
                }

                if (i == m - 1) break;

                H( i + 1, i)= ExX.Norm(w, true);
                v[i + 1]= w*(1.0/H( i + 1, i));

                for (int k= 0; k < i; ++k)
                    GMRES_ApplyPlaneRotation( H(k,i), H(k + 1, i), cs[k], sn[k]);

                GMRES_GeneratePlaneRotation( H(i,i), H(i+1,i), cs[i], sn[i]);
                GMRES_ApplyPlaneRotation( H(i,i), H(i+1,i), cs[i], sn[i]);
                GMRES_ApplyPlaneRotation( s[i], s[i+1], cs[i], sn[i]);

                resid= std::abs( s[i+1])/normb;
                if (output)
                    (*output) << "ParGMRES: " << j << " resid " << resid << std::endl;

                if (resid <= tol) {
                    if (method == RightPreconditioning)
                    {
                        z=0.;
                        GMRES_Update( z, i, H, s, v);
                        M.Apply( A, t, z, ExX);
                        x_acc+=t;
                    }
                    else
                        GMRES_Update( x_acc, i, H, s, v);
                    tol= resid;
                    max_iter= j;
                    return true;
                }
            }

            if (method == RightPreconditioning)
            {
                z=0.;
                GMRES_Update( z, i-1, H, s, v);
                M.Apply( A, t, z, ExX);
                x_acc+=t;
                r= ExX.GetAccumulate(Vec(b - A*x_acc));
            }
            else
            {
                GMRES_Update( x_acc, i-1, H, s, v);
                M.Apply( A, r, Vec( b - A*x_acc), ExX);
            }
            beta=ExX.Norm(r, true);
            resid= beta/normb;
            if (resid <= tol) {
                tol= resid;
                max_iter= j;
                return true;
            }
        }
        tol= resid;
        return false;
    }
}


/// \brief Parallel GMRES-method with modifications for better scalability.
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
  bool ModGMRES(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M,
                   int m, int& max_iter, double& tol,
                   bool measure_relative_tol=true, bool useMGS=false, PreMethGMRES method=LeftPreconditioning, std::ostream* output=0)
    /// \param[in]     A                    local distributed coefficients-matrix of the linear equation system
    /// \param[in,out] x_acc                start vector and the solution in accumulated form
    /// \param[in]     b                    rhs of the linear equation system (distributed form)
    /// \param[in]     ExX                  ExchangeCL corresponding to the RowIdx of x and the ColIdx of A
    /// \param[in,out] M                    Preconditioner
    /// \param[in]     m                    number of steps after a restart is performed
    /// \param[in,out] max_iter             IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol                  IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol if true stop if |M^(-1)(b-Ax)|/|M^(-1)b| <= tol, else stop if |M^(-1)(b-Ax)|<=tol
    /// \param[in]     useMGS               use modified Gram-Schmidt orthogonalization (many more sync-points exists!)
    /// \param[in]     method               left or right preconditioning (see solver.h for definition and declaration)
    /// \return  convergence within max_iter iterations
    /// \pre     the preconditioner should be able to handle a accumulated b
{
    // Check if preconditioner needs diagonal of matrix. The preconditioner
    // only recomputes the diagonal, if the matrix has changed
    if (M.NeedDiag())
        M.SetDiag(A, ExX);

    const size_t n = x_acc.size();      // dimension

    DMatrixCL<double> H(m,m);           // upper Hessenberg-matrix
    VectorCL tmpHCol(m);
    double beta, normb, resid;

    Vec r(n), w(n), w_acc(n), r_acc(n), z_acc(n), t_acc(n);
    Vec c(m), s(m), gamma(m);

    std::vector<Vec> v_acc(m);    // basis of the krylov-subspaces
    for (int i=0; i<m; ++i)
        v_acc[i].resize(n);

    if (method == RightPreconditioning){
        r    = b - A*x_acc;
        beta = ExX.Norm(r, false, &r_acc);
        normb= ExX.Norm(b, false);
    }
    else{
        M.Apply(A, r, VectorCL( b-A*x_acc), ExX);
        beta = ExX.Norm(r, M.RetAcc(), &r_acc);
        M.Apply(A, w, b, ExX);
        normb = ExX.Norm(w, M.RetAcc(), &w_acc);
    }

    if (normb == 0. || measure_relative_tol==false) normb=1.0;

    resid = beta/normb;
    if (resid<=tol){                        // finished
        tol = resid; max_iter = 0; return true;
    }

    int j=1;                                // number of steps
    while (j<=max_iter)
    {
        v_acc[0] = r_acc * (1./beta);
        gamma    = 0.;
        gamma[0] = beta;

        int i;
        for (i=0; i<m-1 && j<=max_iter; ++i, ++j)
        {
            if (method == RightPreconditioning){
                M.Apply(A, w_acc, v_acc[i], ExX);                // hopefully M does the right thing
                w = A*w_acc;
                w_acc = ExX.GetAccumulate(w);
            }
            else{
                M.Apply(A, w, A*v_acc[i], ExX);
                if (!M.RetAcc())
                    w_acc = ExX.GetAccumulate(w);
                else
                    w_acc = w;
            }

            // Gram-Schmidt orthogonalization without  update of w!
            if (!useMGS)
                StandardGramSchmidt(H, w_acc, true, v_acc, true, i, ExX, tmpHCol);
            else
                ModifiedGramSchmidt(H, w_acc, true, v_acc, true, i, ExX);

            H(i+1,i) = ExX.Norm(w_acc, true);
            v_acc[i+1] = w_acc * (1.0 / H(i+1,i));

            for (int k=0; k<i; ++k)
                GMRES_ApplyPlaneRotation(H(k,i),H(k+1,i), c[k], s[k]);

            GMRES_GeneratePlaneRotation(H(i,i), H(i+1,i), c[i], s[i]);
            GMRES_ApplyPlaneRotation(H(i,i), H(i+1,i), c[i], s[i]);
            GMRES_ApplyPlaneRotation(gamma[i], gamma[i+1], c[i], s[i]);

            resid = std::abs(gamma[i+1])/normb;
            if (output)
                (*output) << "ParModGMRES: " << j << " resid " << resid << std::endl;

            if (resid<=tol){            // finished
                if (method == RightPreconditioning){
                    z_acc=0.;
                    GMRES_Update( z_acc, i, H, gamma, v_acc);
                    M.Apply( A, t_acc, z_acc, ExX);              // hopefully M does the right thing
                    x_acc+=t_acc;
                }
                else
                    GMRES_Update(x_acc, i, H, gamma, v_acc);
                tol = resid; max_iter = j; return true;
            }
        }
        if (method == RightPreconditioning){
            z_acc=0.;
            GMRES_Update( z_acc, i-1, H, gamma, v_acc);
            M.Apply( A, t_acc, z_acc, ExX);                      // hopefully M does the right thing
            x_acc += t_acc;
            r      = b-A*x_acc;
            beta = ExX.Norm(r, false, &r_acc);
        }
        else{
            GMRES_Update(x_acc, i-1, H, gamma, v_acc);
            M.Apply(A, r, static_cast<Vec>( b-A*x_acc), ExX);
            beta = ExX.Norm(r, M.RetAcc(), &r_acc);
        }


        resid = beta/normb;
        if (resid<=tol){                // finished
            tol = resid; max_iter=j; return true;
        }
    }
    tol = std::fabs(gamma[0]);
    return false;
}
#endif

/** One recursive step of Lanzcos' algorithm for computing an ONB (q1, q2, q3,...)
 * of the Krylovspace of A for a given starting vector r.
 *
 * This is a three term recursion, computing the next q_i from the two previous ones.
 * See Arnold Reusken, "Numerical methods for elliptic partial differential equations", p. 148.
 * Returns false for 'lucky breakdown' (see below), true in the generic case.
 */
template <typename Mat, typename Vec>
bool LanczosStep(const Mat& A,
            const Vec& q0, const Vec& q1, Vec& q2,
            double& a1,
            const double b0, double& b1)
{
    q2= A*q1 - b0*q0;
    a1= dot( q2, q1);
    q2-= a1*q1;
    b1= norm( q2);
    // Lucky breakdown; the Krylov-space K up to q1 is A-invariant. Thus,
    // the correction dx needed to solve A(x0+dx)=b is in this space and
    // the Minres-algo will terminate with the exact solution in the
    // following step.
    if (b1 < 1e-15) return false;
    q2/= b1;
    return true;
}

template <typename Vec>
class LanczosONBCL
{
  private:
    bool nobreakdown_;
    double norm_r0_;

  public:
    SBufferCL<Vec, 3> q;
    double a0;
    SBufferCL<double, 2> b;

    /// Sets up initial values and computes q0.
    template <typename Mat, typename ExT>
    void new_basis(const Mat& A, const Vec& r0, const ExT&)
    {
        q[-1].resize( r0.size(), 0.);
        norm_r0_= norm( r0);
        q[0].resize( r0.size(), 0.); q[0]= r0/norm_r0_;
        q[1].resize( r0.size(), 0.);
        b[-1]= 0.;
        nobreakdown_= LanczosStep( A, q[-1], q[0], q[1], a0, b[-1], b[0]);
    }

    double norm_r0()   const { return norm_r0_; }
    bool   breakdown() const { return !nobreakdown_; }

    /// Computes new q_i, a_i, b_1, q_{i+1} in q0, a0, b0, q1 and moves old values to qm1, bm1.
    template <typename Mat, typename ExT>
    bool next( const Mat& A, const ExT&)
    {
        q.rotate(); b.rotate();
        return (nobreakdown_= LanczosStep( A, q[-1], q[0], q[1], a0, b[-1], b[0]));
    }
};


/** One recursive step of the preconditioned Lanzcos' algorithm for computing an ONB (q1, q2, q3,...)
 * of the Krylovspace of A for a given starting vector r.
 *
 * This is a three term recursion, computing the next q_i from the two previous ones.
 * See Arnold Reusken, "Numerical methods for elliptic partial differential equations", p. 153.
 * Returns false for 'lucky breakdown' (see below), true in the generic case.
 */
template <typename Mat, typename Vec, typename PreCon, typename ExT>
bool PLanczosStep(const Mat& A,
             const PreCon& M,
             const Vec& q1, Vec& q2,
             const Vec& t0, const Vec& t1, Vec& t2,
             double& a1,
             const double b0, double& b1, const ExT& ex)
{
    t2= A*q1 - b0*t0;
    a1= dot( t2, q1);
    t2-= a1*t1;
    M.Apply( A, q2, t2, ex);
    const double b1sq= dot( q2, t2);
    Assert( b1sq >= 0.0, "PLanczosStep: b1sq is negative!\n", DebugNumericC);
    b1= std::sqrt( b1sq);
    if (b1 < 1e-15) return false;
    t2*= 1./b1;
    q2*= 1./b1;
    return true;
}

template <typename Vec, typename PreCon>
class PLanczosONBCL
{
  private:
    bool nobreakdown_;
    double norm_r0_;

  public:
    const PreCon& M;
    SBufferCL<Vec, 2> q;
    SBufferCL<Vec, 3> t;
    double a0;
    SBufferCL<double, 2> b;

    // Sets up initial values and computes q0.
    PLanczosONBCL(const PreCon& M_)
      : M( M_) {}

    /// Sets up initial values and computes q0.
    template <typename Mat, typename ExT>
    void new_basis(const Mat& A, const Vec& r0, const ExT& ex)
    {
        t[-1].resize( r0.size(), 0.);
        q[-1].resize( r0.size(), 0.); M.Apply( A, q[-1], r0, ex);
        norm_r0_= std::sqrt( dot( q[-1], r0));
        t[0].resize( r0.size(), 0.); t[0]= r0/norm_r0_;
        q[0].resize( r0.size(), 0.); q[0]= q[-1]/norm_r0_;
        t[1].resize( r0.size(), 0.);
        b[-1]= 0.;
        nobreakdown_= PLanczosStep( A, M, q[0], q[1], t[-1], t[0], t[1], a0, b[-1], b[0], ex);
    }

    double norm_r0()   const { return norm_r0_; }
    bool   breakdown() const { return !nobreakdown_; }

    /// Computes new q_i, t_i, a_i, b_1, q_{i+1} in q0, t_0, a0, b0, q1 and moves old values to qm1, tm1, bm1.
    template <typename Mat, typename ExT>
    bool next(const Mat& A, const ExT& ex)
    {
        q.rotate(); t.rotate(); b.rotate();
        return (nobreakdown_= PLanczosStep( A, M, q[0], q[1], t[-1], t[0], t[1], a0, b[-1], b[0], ex));
    }
};

//-----------------------------------------------------------------------------
// PMINRES: The return value indicates convergence within max_iter (input)
// iterations (true), or no convergence within max_iter iterations (false).
// See Arnold Reusken, "Numerical methods for elliptic partial differential
// equations", pp. 149 -- 154
//
// Upon successful return, output arguments have the following values:
//
//        x - approximate solution to Ax = rhs
// max_iter - number of iterations performed before tolerance was reached
//      tol - (relative, see below) residual b - Ax measured in the (M^-1 ., .)-
//     inner-product-norm.
// measure_relative_tol - If true, stop if (M^(-1)( b - Ax), b - Ax)/(M^(-1)b, b) <= tol,
//     if false, stop if (M^(-1)( b - Ax), b - Ax) <= tol.
//-----------------------------------------------------------------------------
template <typename Mat, typename Vec, typename Lanczos, typename ExT>
bool PMINRES(const Mat& A, Vec& x, const Vec&, const ExT& ex, Lanczos& q, int& max_iter, double& tol, bool measure_relative_tol= false)
{
    Vec dx( x.size());
    const double norm_r0= q.norm_r0();
    double normb= std::fabs( norm_r0);
    double res= norm_r0;
    bool lucky= q.breakdown();
    SBufferCL<double, 3> c;
    SBufferCL<double, 3> s;
    SBufferCL<SVectorCL<3>, 3> r;
    SBufferCL<Vec, 3> p;
    p[0].resize( x.size()); p[1].resize( x.size()); p[2].resize( x.size());
    SBufferCL<SVectorCL<2>, 2> b;

    if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;

    if ((res= norm_r0/normb) <= tol) {
        tol= res;
        max_iter= 0;
        return true;
    }
    std::cout << "PMINRES: k: 0\tresidual: " << res << std::endl;
    for (int k= 1; k <= max_iter; ++k) {
        switch (k) {
          case 1:
            // Compute r1
            GMRES_GeneratePlaneRotation( q.a0, q.b[0], c[0], s[0]);
            r[0][0]= std::sqrt( q.a0*q.a0 + q.b[0]*q.b[0]);
            // Compute p1
            // p[0]= q.q[0]/r[0][0];
            p[0]= q.q[0]/r[0][0];
            // Compute b11
            b[0][0]= 1.; b[0][1]= 0.;
            GMRES_ApplyPlaneRotation(b[0][0], b[0][1], c[0], s[0]);
            break;
          case 2:
            // Compute r2
            r[0][0]= q.b[-1]; r[0][1]= q.a0; r[0][2]= q.b[0];
            GMRES_ApplyPlaneRotation( r[0][0], r[0][1], c[-1], s[-1]);
            GMRES_GeneratePlaneRotation( r[0][1], r[0][2], c[0], s[0]);
            GMRES_ApplyPlaneRotation( r[0][1], r[0][2], c[0], s[0]);
            // Compute p2
            // p[0]= (q.q[0] - r[0][0]*p[-1])/r[0][1];
            p[0]= (q.q[0] - r[0][0]*p[-1])/r[0][1];
            // Compute b22
            b[0][0]= b[-1][1]; b[0][1]= 0.;
            GMRES_ApplyPlaneRotation( b[0][0], b[0][1], c[0], s[0]);
            break;
          default:
            r[0][0]= 0.; r[0][1]= q.b[-1]; r[0][2]= q.a0;
            double tmp= q.b[0];
            GMRES_ApplyPlaneRotation( r[0][0], r[0][1], c[-2], s[-2]);
            GMRES_ApplyPlaneRotation( r[0][1], r[0][2], c[-1], s[-1]);
            GMRES_GeneratePlaneRotation( r[0][2], tmp, c[0], s[0]);
            GMRES_ApplyPlaneRotation( r[0][2], tmp, c[0], s[0]);
            // p[0]= (q.q[0] - r[0][0]*p[-2] -r[0][1]*p[-1])/r[0][2];
            p[0]= (q.q[0] - r[0][0]*p[-2] -r[0][1]*p[-1])*(1/r[0][2]);
            b[0][0]= b[-1][1]; b[0][1]= 0.;
            GMRES_ApplyPlaneRotation( b[0][0], b[0][1], c[0], s[0]);
        }
        dx= norm_r0*b[0][0]*p[0];
        x+= dx;

        res= std::fabs( norm_r0*b[0][1])/normb;
        if (k%10==0) std::cout << "PMINRES: k: " << k << "\tresidual: " << res << std::endl;
        if (res<= tol || lucky==true) {
            tol= res;
            max_iter= k;
            return true;
        }
        q.next( A, ex);
        if (q.breakdown()) {
            lucky= true;
            std::cout << "PMINRES: lucky breakdown" << std::endl;
        }
        c.rotate(); s.rotate(); r.rotate(); p.rotate(); b.rotate();
    }
    tol= res;
    return false;
}


template <typename Mat, typename Vec, typename ExT>
bool MINRES(const Mat& A, Vec& x, const Vec& rhs, const ExT& ex, int& max_iter, double& tol, bool measure_relative_tol= false)
{
    LanczosONBCL<Vec> q;
    q.new_basis( A, Vec( rhs - A*x), ex);
    return PMINRES( A,  x, rhs, ex, q, max_iter, tol, measure_relative_tol);
}


//*****************************************************************
// BiCGSTAB
//
// BiCGSTAB solves the non-symmetric linear system Ax = b
// using the Preconditioned BiConjugate Gradient Stabilized method.
//
// BiCGSTAB follows the algorithm described on p. 27 of the
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (true), or no convergence or breakdown within
// max_iter iterations (false). In cases of breakdown a message is printed
// std::cout and the iteration returns the approximate solution found.
//
// Upon successful return, output arguments have the following values:
//
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
// measure_relative_tol - If true, stop if |b - Ax|/|b| <= tol,
//     if false, stop if |b - Ax| <= tol. ( |.| is the euclidean norm.)
//
//*****************************************************************
template <class Mat, class Vec, class Preconditioner, class ExT>
bool BICGSTAB( const Mat& A, Vec& x, const Vec& b, const ExT& ex,
         Preconditioner& M, int& max_iter, double& tol, bool measure_relative_tol= true)
{
    double rho_1= 0.0, rho_2= 0.0, alpha= 0.0, beta= 0.0, omega= 0.0;
    Vec p( x.size()), phat( x.size()), s( x.size()), shat( x.size()),
        t( x.size()), v( x.size());

    double normb= norm( b);
    Vec r( b - A*x);
    Vec rtilde= r;

    if (normb == 0.0 || measure_relative_tol == false) normb = 1.0;

    double resid= norm( r)/normb;
    if (resid <= tol) {
        tol= resid;
        max_iter= 0;
        return true;
    }

    for (int i= 1; i <= max_iter; ++i) {
        rho_1= dot( rtilde, r);
        if (rho_1 == 0.0) {
            tol = norm( r)/normb;
            max_iter= i;
            std::cout << "BiCGSTAB: Breakdown with rho_1 = 0.\n";
            return false;
        }
        if (i == 1)
            p= r;
        else {
            beta= (rho_1/rho_2)*(alpha/omega);
            p= r + beta*(p - omega*v);
        }
        M.Apply( A, phat, p, ex);
        v= A*phat;
        alpha= rho_1/dot( rtilde, v);
        s= r - alpha*v;
        if ((resid= norm( s)/normb) < tol) {
            x+= alpha*phat;
            tol= resid;
            max_iter= i;
            return true;
        }
        M.Apply( A, shat, s, ex);
        t= A*shat;
        omega= dot( t, s)/dot( t, t);
        x+= alpha*phat + omega*shat;
        r= s - omega*t;

        rho_2= rho_1;
        if ((resid= norm( r)/normb) < tol) {
            tol= resid;
            max_iter= i;
            return true;
        }
        if (omega == 0.0) {
            tol= norm( r)/normb;
            max_iter= i;
            std::cout << "BiCGSTAB: Breakdown with omega_ = 0.\n";
            return false;
        }
    }
    tol= resid;
    return false;
}


/// \brief Preconditioned BiCGStab-Method with modifications to minimize global communications
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
  bool ModBICGSTAB(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX,
                   PreCon& M, int& max_iter, double& tol, bool measure_relative_tol=true)
    /// \param[in]     A                    local distributed coefficients-matrix of the linear equation system
    /// \param[in,out] x_acc                start vector and the solution in accumulated form
    /// \param[in]     b                    rhs of the linear equation system (distributed form)
    /// \param[in]     ExX                  ExchangeCL corresponding to the RowIdx of x and the ColIdx of A
    /// \param[in]     M                    Preconditioner
    /// \param[in,out] max_iter             IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol                  IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol if true stop if |M^(-1)(b-Ax)|/|M^(-1)b| <= tol, else stop if |M^(-1)(b-Ax)|<=tol
    /// \pre                                the preconditioner must fulfill the following condition: M.Apply(A,x,b): b is accumulated => x is accumulated
    /// \return                             convergence within max_iter iterations
{
    // Check if preconditioner needs diagonal of matrix. The preconditioner
    // only computes the diagonal new, if the matrix has changed
    if (M.NeedDiag())
        M.SetDiag(A, ExX);

    const size_t n = x_acc.size();
    Vec r(b-A*x_acc), r_acc(n), r0hat_acc(n),
        p_acc(n), phat_acc(n),
        v(n), v_acc(n),
        s_acc(n), shat_acc(n),
        t(n), t_acc(n), that(n), that_acc(n);

    double resid, rho=1, alpha=1, beta, omega=1, sigma;

    VectorCL dots(4), glob_dots(4);

    double normb = ExX.Norm(b, false);
    if ( normb==0. || !measure_relative_tol) normb= 1.;

    sigma = ExX.Norm(r, false, &r_acc);
    resid = std::sqrt(sigma) / normb;

    M.Apply(A, r0hat_acc, r, ExX);
    if (!M.RetAcc())
        r0hat_acc= ExX.GetAccumulate(r0hat_acc);

    for (int i=0; i<max_iter; ++i)
    {
        beta = (sigma*alpha)/(rho*omega);
        rho  = sigma;
        if (i>0)
            p_acc    = r_acc + beta*p_acc - (beta*omega)*v_acc;
        else
            p_acc    = r_acc;

        M.Apply(A, phat_acc, p_acc, ExX);

        v= A*phat_acc;

        dots[0]= ExX.LocalDot( v, false, r0hat_acc, true, &v_acc);
        dots[1]= ExX.LocalNorm_sq( r_acc, true);
#ifdef _PAR
        ProcCL::GlobalSum(Addr(dots), Addr(glob_dots), 2);
#else
        glob_dots = dots;
#endif

        resid = std::sqrt(glob_dots[1])/normb;
        if (resid<tol){
            tol= resid;
            max_iter=i;
            return true;
        }

        sigma = glob_dots[0];

        if (glob_dots[0]==0.)
        {
            std::cout << ">>>>> BREAKDOWN of BiCGStab!" <<std::endl;
            tol = resid;
            max_iter=i;
            return false;
        }

        alpha = rho/sigma;
        s_acc = r_acc -alpha*v_acc;

        M.Apply(A, shat_acc, s_acc, ExX);
        t = A*shat_acc;

        if (!M.RetAcc()){
            M.Apply(A,that,t, ExX);
            dots[0]= ExX.LocalDot( that, false, shat_acc, true, &that_acc);
        }
        else{
            M.Apply(A,that_acc,t, ExX);
            dots[0]= ExX.LocalDot(that_acc, true, shat_acc, true);
        }

        dots[1]= ExX.LocalDot(that_acc, true, that_acc, true);
        dots[2]= ExX.LocalDot(r0hat_acc, true, s_acc, true);
        dots[3]= ExX.LocalDot( t, false, r0hat_acc, true, &t_acc);

#ifdef _PAR
        ProcCL::GlobalSum(Addr(dots), Addr(glob_dots), 4);
#else
        glob_dots = dots;
#endif

        if (glob_dots[1]==0.)
        {
            std::cout << ">>>>> BREAKDOWN of BiCGStab!" <<std::endl;
            tol = resid;
            max_iter=i;
            return false;
        }

        omega = glob_dots[0]/glob_dots[1];
        sigma = glob_dots[2] - omega * glob_dots[3];

        r_acc =  s_acc - omega * t_acc;
        x_acc += alpha * phat_acc + omega * shat_acc;
    }

    tol = resid;
    return false;
}


//*****************************************************************
// GCR
//
// GCR solves the non-symmetric linear system Ax = b
// using the Preconditioned Generalized Conjugate Residuals method;
//
// The return value indicates convergence within max_iter (input)
// iterations (true), or no convergence within
// max_iter iterations (false).
//
// Upon successful return, output arguments have the following values:
//
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
//
// m -- truncation parameter; only m >= 1 residual vectors are kept; the
//     one with the smallest a (in modulus) is overwritten by the
//     new vector sn (min-alpha strategy).
// measure_relative_tol -- If true, stop if |b - Ax|/|b| <= tol,
//     if false, stop if |b - Ax| <= tol. ( |.| is the euclidean norm.)
//
//*****************************************************************
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool GCR(const Mat& A, Vec& x, const Vec& b, const ExCL& ExX, PreCon& M,
    int m, int& max_iter, double& tol, bool measure_relative_tol=true, std::ostream* output=0)
{
    if (M.NeedDiag())
        M.SetDiag(A, ExX);

    m= (m <= max_iter) ? m : max_iter; // m > max_iter only wastes memory.

    Vec r( b - A*x);
    Vec racc( r.size());
    Vec sn( b.size()), vn( b.size());
    Vec vnacc( b.size());
    std::vector<Vec> s, v;
    std::vector<Vec> vacc;         // memory usage vs communication overhead
    std::vector<double> a( m);

    double normb= ExX.Norm( b, false);
    if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;
    double resid= ExX.Norm( r, false, &racc)/normb;
    for (int k= 0; k < max_iter; ++k) {
        if (k%10==0 && output)
            (*output) << "GCR: k: " << k << "\tresidual: " << resid << std::endl;
        if (resid < tol) {
            tol= resid;
            max_iter= k;
            return true;
        }
        M.Apply( A, sn, r, ExX);
        if (!M.RetAcc())
            vn=A*ExX.GetAccumulate(sn);
        else
            vn= A*sn;
        vnacc = ExX.GetAccumulate(vn);
        for (int i= 0; i < k && i < m; ++i) {
//            const double alpha= ExX.ParDot( vn, false, v[i], false);
            const double alpha= ExX.ParDot( vnacc, true, vacc[i], true);
            a[i]= alpha;
            vn-= alpha*v[i];
            vnacc-= alpha*vacc[i];  // memory usage vs communication overhead
            sn-= alpha*s[i];
        }
//        const double beta= ExX.Norm( vn, false, &vnacc);
        const double beta= ExX.Norm( vnacc, true);
        vn/= beta;
        vnacc/= beta;
        sn/= beta;
        const double gamma= ExX.ParDot( racc, true, vnacc, true);
        if (!M.RetAcc())
            x+= gamma*ExX.GetAccumulate(sn);
        else
            x+= gamma*sn;
        r-= gamma*vn;
        racc -= gamma*vnacc;
        resid= ExX.Norm( racc, true)/normb;
        if (k < m) {
            s.push_back( sn);
            v.push_back( vn);
            vacc.push_back( vnacc);   // memory usage vs communication overhead
        }
        else {
#ifdef _PAR
            throw DROPSErrCL("GCR: Sorry, parallel truncation not implemented");
#endif
            int min_idx= 0;
            double a_min= std::fabs( a[0]); // m >= 1, thus this access is valid.
            for (int i= 1; i < k && i < m; ++i)
                if ( std::fabs( a[i]) < a_min) {
                    min_idx= i;
                    a_min= std::fabs( a[i]);
                }
            s[min_idx]= sn;
            v[min_idx]= vn;
        }
    }
    tol= resid;
    return false;
}


//*****************************************************************
// GMRESR
//
// GMRESR solves the non-symmetric linear system Ax = b
// using the GMRES Recursive method;
//
// The return value indicates convergence within max_iter (input)
// iterations (true), or no convergence within
// max_iter iterations (false).
//
// Upon successful return, output arguments have the following values:
//
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
// measure_relative_tol - If true, stop if |b - Ax|/|b| <= tol,
//     if false, stop if |b - Ax| <= tol. ( |.| is the euclidean norm.)
//
//*****************************************************************
template <class Mat, class Vec, class Preconditioner, class ExT>
bool GMRESR( const Mat& A, Vec& x, const Vec& b, const ExT& ex, Preconditioner& M,
    int /*restart parameter m*/ m, int& max_iter, int& inner_max_iter, double& tol, double& inner_tol,
    bool measure_relative_tol= true, PreMethGMRES method = RightPreconditioning)
{
    Vec r( b - A*x);
    std::vector<Vec> u(1), c(1); // Positions u[0], c[0] are unused below.
    double normb= norm( b);
    if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;
    double resid= -1.0;

    for (int k= 0; k < max_iter; ++k) {
        if ((resid= norm( r)/normb) < tol) {
            tol= resid;
            max_iter= k;
            return true;
        }
        std::cout << "GMRESR: k: " << k << "\tresidual: " << resid << std::endl;
        u.push_back( Vec( b.size()));
        u[k+1]=0;
        inner_tol=0.0;
        double in_tol = inner_tol;
        int in_max_iter = inner_max_iter;
#ifdef _PAR
        GMRES(A, u[k+1], r, ex, M, m, in_max_iter, in_tol, true, method);
#else
        GMRES(A, u[k+1], r, ex, M, m, in_max_iter, in_tol, true, false, method);
#endif
        std::cout << "norm of u_k_0: "<<norm(u[k+1])<<"\n";
        std::cout << "inner iteration:  " << in_max_iter << " GMRES iteration(s),\tresidual: " << in_tol << std::endl;
        if (norm(A*u[k+1]-r)>0.999*norm(r) && norm(u[k+1]) < 1e-3)
        {
            u[k+1] = transp_mul(A, r);
            std::cout<<"LSQR switch!\n";
        }
        c.push_back( A*u[k+1]);
        for (int i= 1; i <= k; ++i) {
            const double alpha= dot( c[k+1], c[i]);
            c[k+1]-= alpha*c[i];
            u[k+1]-= alpha*u[i];
        }
        const double beta= norm( c[k+1]);
        c[k+1]/= beta;
        u[k+1]/= beta;
        const double gamma= dot( r, c[k+1]);
        x+= gamma*u[k+1];
        r-= gamma*c[k+1];
    }
    tol= resid;
    return false;
}


//*****************************************************************
// IDR(s)
//
// IDR(s) solves the non-symmetric linear system Ax = b
// using the method by Martin B. van Gijzen / Peter Sonneveld
// An Elegant IDR(s) Variant that Efficiently Exploits
// Bi-Orthogonality Properties
// ISSN 1389-6520
// Department of Applied Mathematical Analysis,  Report 10-16
// Delft University of Technology, 2010
//
// The return value indicates convergence within max_iter (input)
// iterations (true), or no convergence within
// max_iter iterations (false).
//
// Upon successful return, output arguments have the following values:
//
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
// measure_relative_tol - If true, stop if |b - Ax|/|b| <= tol,
//     if false, stop if |b - Ax| <= tol. ( |.| is the Euclidean norm.)
//
//*****************************************************************
template <class Mat, class Vec, class PC, class ExT>
bool IDRS( const Mat& A, Vec& x, const Vec& rhs, const ExT& ex, PC& pc, int& max_iter, double& tol,
        bool measure_relative_tol= false, const int s=4, typename Vec::value_type omega_bound=0.7)
{
    typedef typename Vec::value_type ElementTyp;
    int n= x.size(),  it;

    ElementTyp omega = 1.0;
    Vec v(n), t(n), f(s);

    double normb = norm(rhs);
    if (normb == 0.0 || measure_relative_tol == false)
    normb= 1.0;

    Vec resid( rhs - A*x);
    double normres = norm (resid);
    if ( normres/normb <= tol ) {
        tol= normres;
        max_iter= 0;
        return true;
    }

    std::vector<Vec> P(s, Vec( x.size()));
    srand ( time(NULL) );
    for (int i=0; i<s; ++i){
        for ( size_t j=0; j< x.size(); ++j)
            P[i][j] = (double) (2.0 *rand()) / RAND_MAX - 1.0;  // random in [-1,1)
    }

    // orthonormalize P   (mod. Gram-Schmidt)
    for ( int k= 0; k < s; k++ ) {
        for ( int j= 0; j < k; j++ ) {
            ElementTyp sm= dot( P[j], P[k]);
            P[k]-= sm*P[j];
        }
        P[k]/= norm(P[k]);
    }

    std::vector<Vec> G(s, Vec( x.size())), U(s, Vec( x.size()));
    DMatrixCL<ElementTyp> M(s,s);
    for (int i = 0; i <s; ++i)
        M(i,i) = 1.0;
    it= 0;
    while ( normres/normb > tol && it < max_iter ) {
        for (int i=0; i < s; i++)  f[i]= dot(P[i], resid);
        for (int k=0; k < s; k++) {
//            if (it % 10 == 0) std::cout << "IDR(s) iter: " << it << "\tres: " << normres << std::endl;
            // solve the lower tridiagonal system  M[k:s,k:s] c = f[k:s]
            Vec c(s-k);
            for (int j=0; j < s-k; j++) {
                double cs = f[k+j];
                for (int l=0; l < j; l++) cs-= M(k+j,k+l)*c[l];
                c[j]= cs/M(k+j,k+j);
            }
            v= resid;
            for (int j=0; j < s-k; j++) v-= c[j]*G[k+j];
            pc.Apply( A, v, v, ex);

            // Compute new U(:,k) and G(:,k), G(:,k) is in space G_j
            U[k]= c[0]*U[k] + omega*v;
            for (int j=1; j < s-k; j++) U[k]+= c[j]*U[k+j];
            G[k]= A * U[k];
            // Bi-Orthogonalize the new basis vectors
            for (int i= 0; i < k; i++) {
                ElementTyp alpha= dot (P[i], G[k])/M(i,i);
                G[k]-= alpha*G[i];
                U[k]-= alpha*U[i];
            }
            // compute new column of M (first k-1 entries are zero)
            for (int j=0; j < s-k; j++) M(k+j,k)= dot (P[k+j], G[k]);
            if ( M(k,k) == 0 )
            {
                throw DROPSErrCL( "IDR(s): M(k,k) ==0");
            }

            //    make  R orthogonal to  G
            ElementTyp beta = f[k] / M(k,k);
            resid-= beta*G[k];
            x+= beta*U[k];
            normres= norm(resid);
            it++;
            if ( normres/normb <= tol)   break;
            if ( k+1 < s ) {
                for (int j=1; j < s-k; j++) f[k+j]-= beta*M(k+j,k);
            }
        }  //  end of  k = 0 .. (s-1)

        if ( normres/normb <= tol)   break;
        // Entering  G+
        pc.Apply( A, v, resid, ex);
        t= A*v;
        double tn = norm(t), tr = dot(t, resid);
        omega= tr/(tn*tn);
        ElementTyp rho= std::abs(tr)/(tn*normres);
        if ( rho < omega_bound )  omega*= omega_bound/rho;
        if ( omega == 0 ) {
            throw DROPSErrCL( "IDR(s): omega ==0");
        }
        resid-= omega*t;  x+= omega*v;
        normres= norm(resid);
        it++;
    }
    if (tol > normres/normb) {
        max_iter = it;
        tol = normres;
        return false;
    }
    max_iter = it;
    tol = normres;

    return true;
}

#ifdef _PAR
/// \brief Preconditioned QMR-Method
template <typename Mat, typename Vec, typename Lanczos, typename ExCL>
bool QMR(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, Lanczos lan,
            int& max_iter, double& tol, bool measure_relative_tol)
    /// \param[in]     A                    local distributed coefficients-matrix of the linear equation system
    /// \param[in,out] x_acc                start vector and the solution in accumulated form
    /// \param[in]     b                    rhs of the linear equation system (distributed form)
    /// \param[in]     ExX                  ExchangeCL corresponding to the RowIdx of x and the ColIdx of A
    /// \param[in]     lan                  Lanczos-Class for computing a bi-orthogonal basis of Krylov subspaces
    /// \param[in,out] max_iter             IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol                  IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol if true stop if |M^(-1)(b-Ax)|/|M^(-1)b| <= tol, else stop if |M^(-1)(b-Ax)|<=tol
    /// \return                             convergence within max_iter iterations
{
    tol*=tol;

    Vec d_acc(x_acc.size()), s(d_acc), r(b-A*x_acc), r_acc(r);

    double normb = ExX.Norm_sq(b, false);
    if (normb==0.0 || measure_relative_tol==false) normb=1.0;
    double norm_r = ExX.Norm( r, false, &r_acc);

    if (norm_r/normb<std::sqrt(tol))
    {
        tol = std::fabs(norm_r)/std::sqrt(normb);
        max_iter=1;
        return true;
    }

    lan.Init(A,static_cast<Vec>(r*(1./norm_r)),static_cast<Vec>(r*(1./norm_r)),ExX,1);
    double tetha,
    kappa=-1,
    lambda=1,
    beta,
    rho = norm_r,
    rho_last = norm_r;
//  rho_last=lan.GetRho();
    double quot;

    for (int j=1; j<=max_iter; ++j)
    {
        if (!lan.Step()){
            tol = std::sqrt(norm_r/normb);
            max_iter=j-1;
            return false;
        }

        beta    = lan.GetBeta();
        rho_last= rho;
        rho     = lan.GetRho();
        quot    = (lambda*beta*beta+rho*rho);
        tetha   = beta*beta * (1-lambda)/quot;
        kappa   = -rho_last*beta*kappa/quot;
        lambda  = lambda*beta*beta/quot;

        d_acc   = tetha * d_acc + kappa*lan.GetAccP();
        s       = tetha * s     + kappa*lan.GetAp();
        x_acc  += d_acc;
        r      -= s;

        norm_r  = ExX.Norm(r, false, &r_acc);

        if (norm_r<0)
            std::cout << "["<<ProcCL::MyRank()<<"]==> negative squared norm of residual in QMR because of accumulation!" << std::endl;
        if (norm_r/normb<tol)
        {
            tol = std::sqrt(std::fabs(norm_r)/normb);
            max_iter=j;
            return true;
        }
    }
    tol = std::sqrt(norm_r);
    return false;
}
#endif

//=============================================================================
//  Drivers
//=============================================================================

/// What every iterative solver should have
class SolverBaseCL
{
  protected:
    int            maxiter_;
    mutable int    iter_;
    double         tol_;
    mutable double res_;
    bool           rel_;

    mutable std::ostream* output_;

    SolverBaseCL (int maxiter, double tol, bool rel= false, std::ostream* output= 0)
        : maxiter_( maxiter), iter_( -1), tol_( tol), res_( -1.),
          rel_( rel), output_( output)  {}
    virtual ~SolverBaseCL() {}

  public:
    virtual void   SetTol     (double tol) { tol_= tol; }
    virtual void   SetMaxIter (int iter)   { maxiter_= iter; }
    virtual void   SetRelError(bool rel)   { rel_= rel; }

    virtual double GetTol     () const { return tol_; }
    virtual int    GetMaxIter () const { return maxiter_; }
    virtual double GetResid   () const { return res_; }
    virtual int    GetIter    () const { return iter_; }
    virtual bool   GetRelError() const { return rel_; }

    virtual void   SetOutput( std::ostream* os) { output_=os; }
};

/// Bare CG solver
class CGSolverCL : public SolverBaseCL
{
  public:
    CGSolverCL(int maxiter, double tol, bool rel= false, std::ostream* output = 0)
        : SolverBaseCL(maxiter, tol, rel, output) {}

    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex)
    {
        res_=  tol_;
        iter_= maxiter_;
        CG(A, x, b, ex, iter_, res_, rel_, output_);
    }
    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex, int& numIter, double& resid) const
    {
        resid=   tol_;
        numIter= maxiter_;
        CG(A, x, b, ex, numIter, resid, rel_, output_);
    }
};

/// CG with preconditioner
template <typename PC>
class PCGSolverCL : public SolverBaseCL
{
  private:
    PC& _pc;

  public:
    PCGSolverCL(PC& pc, int maxiter, double tol, bool rel= false, std::ostream* output = 0)
        : SolverBaseCL(maxiter, tol, rel, output), _pc(pc) {}

    PC&       GetPc ()       { return _pc; }
    const PC& GetPc () const { return _pc; }

    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex)
    {
        res_=  tol_;
        iter_= maxiter_;
        PCG(A, x, b, ex, _pc, iter_, res_, rel_, output_);
    }
    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex, int& numIter, double& resid) const
    {
        resid=   tol_;
        numIter= maxiter_;
        PCG(A, x, b, ex, _pc, numIter, resid, rel_, output_);
    }
};

///\brief Solver for A*A^Tx=b with Craig's method and left preconditioning.
///
/// A preconditioned CG version for matrices of the form A*A^T. Note that *A* must be
/// supplied, not A*A^T, to the Solve-method.
template <typename PC>
class PCGNESolverCL : public SolverBaseCL
{
  private:
    PC& pc_;

  public:
    PCGNESolverCL(PC& pc, int maxiter, double tol, bool rel= false, std::ostream* output=0)
        : SolverBaseCL( maxiter, tol, rel, output), pc_( pc) {}

    PC&       GetPc ()       { return pc_; }
    const PC& GetPc () const { return pc_; }

    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex, const ExT& ex_transp)
    {
        res_=  tol_;
        iter_= maxiter_;
        PCGNE( A, x, b, ex, ex_transp, pc_, iter_, res_, rel_, output_);
    }
    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex, const ExT& ex_transp, int& numIter, double& resid) const
    {
        resid=   tol_;
        numIter= maxiter_;
        PCGNE(A, x, b, ex, ex_transp, pc_, numIter, resid, rel_, output_);
    }
};


/// Bare MINRES solver
class MResSolverCL : public SolverBaseCL
{
  public:
    MResSolverCL(int maxiter, double tol, bool rel= false)
        : SolverBaseCL( maxiter, tol, rel) {}

    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex)
    {
        res_=  tol_;
        iter_= maxiter_;
        MINRES( A, x, b, ex, iter_, res_, rel_);
    }
    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex, int& numIter, double& resid) const
    {
        resid=   tol_;
        numIter= maxiter_;
        MINRES( A, x, b, ex, numIter, resid, rel_);
    }
};

/// Preconditioned MINRES solver
template <typename Lanczos>
class PMResSolverCL : public SolverBaseCL
{
  private:
    Lanczos& q_;

  public:
    PMResSolverCL(Lanczos& q, int maxiter, double tol, bool rel= false)
      :SolverBaseCL( maxiter, tol, rel), q_( q) {}

    Lanczos&       GetONB ()       { return q_; }
    const Lanczos& GetONB () const { return q_; }

    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex)
    {
        res_=  tol_;
        iter_= maxiter_;
        q_.new_basis( A, Vec( b - A*x), ex);
        PMINRES( A, x, b, ex, q_, iter_, res_, rel_);
    }
    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex, int& numIter, double& resid) const
    {
        resid=   tol_;
        numIter= maxiter_;
        q_.new_basis( A, Vec( b - A*x), ex);
        PMINRES( A, x, b, ex, q_, numIter, resid, rel_);
    }
};

/// GMRES
template <typename PC>
class GMResSolverCL : public SolverBaseCL
{
  private:
    PC&          pc_;
    int          restart_;
    bool         calculate2norm_;
    PreMethGMRES method_;
    bool         mod_;
    bool         useModGS_;

  public:
    GMResSolverCL( PC& pc, int restart, int maxiter, double tol,
        bool relative= true, bool calculate2norm= false, PreMethGMRES method= LeftPreconditioning, bool mod = true, bool useModGS = false, std::ostream* output=0)
        : SolverBaseCL( maxiter, tol, relative, output), pc_(pc), restart_(restart),
          calculate2norm_(calculate2norm), method_(method), mod_(mod), useModGS_(useModGS){}

    PC&       GetPc      ()       { return pc_; }
    const PC& GetPc      () const { return pc_; }
    int       GetRestart () const { return restart_; }

    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex)
    {
        res_=  tol_;
        iter_= maxiter_;
#ifndef _PAR
        GMRES(A, x, b, ex, pc_, restart_, iter_, res_, rel_, calculate2norm_, method_);
#else
        if (mod_)
            ModGMRES(A, x, b, ex, pc_, restart_, iter_, res_, rel_, useModGS_, method_, output_);
        else
            GMRES(A, x, b, ex, pc_, restart_, iter_, res_, rel_, method_, output_);
#endif
    }
    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex, int& numIter, double& resid) const
    {
        resid=   tol_;
        numIter= maxiter_;
#ifndef _PAR
        GMRES(A, x, b, ex, pc_, restart_, numIter, resid, rel_, calculate2norm_, method_);
#else
        if (mod_)
            ModGMRES(A, x, b, ex, pc_, restart_, numIter, resid, rel_, useModGS_, method_, output_);
        else
            GMRES(A, x, b, ex, pc_, restart_, numIter, resid, rel_, method_, output_);
#endif
    }
};


/// BiCGStab
template <typename PC>
class BiCGStabSolverCL : public SolverBaseCL
{
  private:
    PC& pc_;
    bool mod_;

  public:
    BiCGStabSolverCL( PC& pc, int maxiter, double tol, bool relative= true, bool mod=true)
        : SolverBaseCL( maxiter, tol, relative), pc_( pc), mod_(mod){}

          PC& GetPc ()       { return pc_; }
    const PC& GetPc () const { return pc_; }

    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex)
    {
        res_=  tol_;
        iter_= maxiter_;
#ifdef _PAR
        ModBICGSTAB( A, x, b, ex, pc_, iter_, res_, rel_);
#else
        BICGSTAB( A, x, b, ex, pc_, iter_, res_, rel_);
#endif
    }
    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex, int& numIter, double& resid) const
    {
        resid=   tol_;
        numIter= maxiter_;
#ifdef _PAR
        ModBICGSTAB(A, x, b, ex, pc_, numIter, resid, rel_);
#else
        BICGSTAB(A, x, b, ex, pc_, numIter, resid, rel_);
#endif
    }
};

/// GCR
template <typename PC>
class GCRSolverCL : public SolverBaseCL
{
  private:
    PC& pc_;
    int truncate_; // no effect for the parallel version

  public:
    GCRSolverCL( PC& pc, int truncate, int maxiter, double tol,
        bool relative= true, std::ostream* output= 0)
        : SolverBaseCL( maxiter, tol, relative, output), pc_( pc),
          truncate_( truncate) {}

    PC&       GetPc      ()       { return pc_; }
    const PC& GetPc      () const { return pc_; }
    int       GetTruncate() const { return truncate_; }

    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex)
    {
        res_=  tol_;
        iter_= maxiter_;
        GCR( A, x, b, ex, pc_, truncate_, iter_, res_, rel_, output_);
        if (output_ != 0)
            *output_ << "GCRSolverCL: iterations: " << GetIter()
                     << "\tresidual: " << GetResid() << std::endl;
    }
    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex, int& numIter, double& resid) const
    {
        resid=   tol_;
        numIter= maxiter_;
        GCR(A, x, b, ex, pc_, truncate_, numIter, resid, rel_, output_);
        if (output_ != 0)
            *output_ << "GCRSolverCL: iterations: " << GetIter()
                    << "\tresidual: " << GetResid() << std::endl;
    }
};

/// GMRESR
template <typename PC>
class GMResRSolverCL : public SolverBaseCL
{
  private:
    PC&          pc_;
    int          restart_;
    int          inner_maxiter_;
    double       inner_tol_;
    PreMethGMRES method_;
  public:
    GMResRSolverCL( PC& pc, int restart, int maxiter, int  inner_maxiter,
        double tol, double  inner_tol, bool relative= true, PreMethGMRES method = RightPreconditioning, std::ostream* output= 0)
        : SolverBaseCL( maxiter, tol, relative, output), pc_(pc), restart_(restart),
         inner_maxiter_( inner_maxiter), inner_tol_(inner_tol), method_(method){}

    PC&       GetPc      ()       { return pc_; }
    const PC& GetPc      () const { return pc_; }
    int       GetRestart () const { return restart_; }
    void   SetInnerTol     (double tol) { inner_tol_= tol; }
    void   SetInnerMaxIter (int iter)   { inner_maxiter_= iter; }

    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex)
    {
        res_=  tol_;
        iter_= maxiter_;
        GMRESR(A, x, b, ex, pc_, restart_, iter_, inner_maxiter_, res_, inner_tol_, rel_, method_);
        if (output_ != 0)
            *output_ << "GmresRSolverCL: iterations: " << GetIter()
                     << "\tresidual: " << GetResid() << std::endl;
            std::cout << "GmresRSolverCL: iterations: " << GetIter()
                     << "\tresidual: " << GetResid() << std::endl;
    }
    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex, int& numIter, double& resid) const
    {
        resid=   tol_;
        numIter= maxiter_;
        GMRESR(A, x, b, ex, pc_, restart_, numIter, inner_maxiter_, resid, inner_tol_, rel_, method_);
        if (output_ != 0)
            *output_ << "GmresRSolverCL: iterations: " << GetIter()
                     << "\tresidual: " << GetResid() << std::endl;
            std::cout << "GmresRSolverCL: iterations: " << GetIter()
                     << "\tresidual: " << GetResid() << std::endl;
    }
};

/// IDR(s)
template <typename PC>
class IDRsSolverCL : public SolverBaseCL
{
  private:
    PC&          pc_;
    const int    s_;
    const double omega_bound_;

  public:
    IDRsSolverCL( PC& pc, int maxiter, double tol, bool relative= true, const int s = 4, const double omega_bound = 0.7, std::ostream* output= 0)
        : SolverBaseCL( maxiter, tol, relative, output), pc_(pc), s_(s), omega_bound_( omega_bound) {}

    PC&       GetPc      ()       { return pc_; }
    const PC& GetPc      () const { return pc_; }

    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex)
    {
        res_=  tol_;
        iter_= maxiter_;
        IDRS(A, x, b, ex, pc_, iter_, res_, rel_, s_, omega_bound_);
        if (output_ != 0)
            *output_ << "IDRsSolverCL: iterations: " << GetIter()
                     << "\tresidual: " << GetResid() << std::endl;
    }
    template <typename Mat, typename Vec, typename ExT>
    void Solve(const Mat& A, Vec& x, const Vec& b, const ExT& ex, int& numIter, double& resid) const
    {
        resid=   tol_;
        numIter= maxiter_;
        IDRS(A, x, b, ex, pc_, iter_, res_, rel_, s_, omega_bound_);
        if (output_ != 0)
            *output_ << "IDRsSolverCL: iterations: " << GetIter()
                     << "\tresidual: " << GetResid() << std::endl;
    }
};

#ifdef _PAR

// ***************************************************************************
/// \brief Parallel QMR-Solver class
// ***************************************************************************
template<typename Lanczos>
class ParQMRSolverCL : public SolverBaseCL
{
  private:
    Lanczos *lan_;
    typedef SolverBaseCL base;

  public:
    /// \brief Constructor of the parallel preconditioned QMR-Solver
    /** Tries to solve a linear equation system within \a maxiter steps with
        accuracy \a tol. The ExCL \a ex is used to do parallel inner products. The
        Lanczos-Algorithm to compute the bi-orthogonal-basis is given by a Lanczos class \a lan.
        If \a measure_relative_tol is given, the residual is computed relative.*/
    ParQMRSolverCL(int maxiter, double tol, Lanczos &lan, bool rel=true, std::ostream* output=0) :
        base(maxiter, tol, rel, output), lan_(&lan) {}

    /// \brief Solve a linear equation system with Quasi Minimal Residual Method
    template <typename Mat, typename Vec, typename ExT>
      void Solve(const Mat& A, Vec& x, const Vec &b, const ExT& ex)
    /// Solve the linear equation system with coefficient matrix \a A and rhs \a b iterative with
    /// QMR algorithm, uses \a x as start-vector and result vector.
    /// \post x has accumulated form
    {
        res_  = tol_;
        iter_ = maxiter_;
        QMR(A, x, b, ex, *lan_, iter_, res_, rel_);
    }
};
#endif

} // end of namespace DROPS

#endif
