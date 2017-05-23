/// \file pardisosolver.h
/// \brief interface to pardiso solver (MKL 11.0)
/// \author LNM RWTH Aachen: Patrick Esser

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

#ifndef PARDISOSOLVER_H_
#define PARDISOSOLVER_H_

#include "mkl_pardiso.h"
#include "mkl_types.h"

namespace DROPS
{
/*******************************************************************
*   P A R D I S O S O L V E R   C L                                *
*******************************************************************/
/// \brief Use direct solvers for Ax=b with sparse matrix A.
/** This class offers an interface to the PARDISO package from
    INTEL MKL 11.0                                                 */
/*******************************************************************
*   P A R D I S O S O L V E R   C L                                *
*******************************************************************/
class PardisoSolverCL {

private:
    /* Matrix data. */
    MKL_INT n;
    MKL_INT* ia;         /* row_beg */
    MKL_INT* ja;         /* col_ind */
    double* a;           /* non zeros */
    size_t mat_version_;

    MKL_INT mtype;       /* 11 Real unsymmetric matrix */
    MKL_INT nrhs;        /* Number of right hand sides. */

    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void *pt[64];
    /* Pardiso control parameters. */
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    double ddum;          /* Double dummy */
    MKL_INT idum;         /* Integer dummy. */

public:

/// store and factorize the matrix A
    PardisoSolverCL(const MatrixCL& A)
  : ia(0), ja(0), a(0), mat_version_(0), mtype(11), nrhs(1), phase(0)
{
    /* -------------------------------------------------------------------- */
    /* .. Setup Pardiso control parameters. */
    /* -------------------------------------------------------------------- */
      for (int i = 0; i < 64; i++)
        {
          iparm[i] = 0;
        }
      iparm[0] = 1;         /* No solver default */
      iparm[1] = 2;         /* Fill-in reordering from METIS */
      /* Numbers of processors, value of OMP_NUM_THREADS */
      iparm[2] = 1;
      iparm[3] = 0;         /* No iterative-direct algorithm */
      iparm[4] = 0;         /* No user fill-in reducing permutation */
      iparm[5] = 0;         /* Write solution into x */
      iparm[6] = 0;         /* Not in use */
      iparm[7] = 2;         /* Max numbers of iterative refinement steps */
      iparm[8] = 0;         /* Not in use */
      iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
      iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
      iparm[11] = 0;        /* Conjugate transposed/transpose solve */
      iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
      iparm[13] = 0;        /* Output: Number of perturbed pivots */
      iparm[14] = 0;        /* Not in use */
      iparm[15] = 0;        /* Not in use */
      iparm[16] = 0;        /* Not in use */
      iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
      iparm[18] = -1;       /* Output: Mflops for LU factorization */
      iparm[19] = 0;        /* Output: Numbers of CG Iterations */
      maxfct = 1;           /* Maximum number of numerical factorizations. */
      mnum = 1;             /* Which factorization to use. */
      msglvl = getenv("PARDISOMSG") ? 1 : 0;           /* Print statistical information in file */
      error = 0;            /* Initialize error flag */
    /* -------------------------------------------------------------------- */
    /* .. Initialize the internal solver memory pointer. This is only */
    /* necessary for the FIRST call of the PARDISO solver. */
    /* -------------------------------------------------------------------- */
      for (int i = 0; i < 64; i++)
        {
          pt[i] = 0;
        }

    Update(A);
}

~PardisoSolverCL()
{
    /* -------------------------------------------------------------------- */
    /* .. Termination and release of memory. */
    /* -------------------------------------------------------------------- */
    phase = -1;           /* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);
    delete[] ia;
    delete[] ja;
    delete[] a;
}
/// delete stored matrix and store/factorize the matrix A
void Update(const MatrixCL& A)
{
    if (phase != 0) {
    /* -------------------------------------------------------------------- */
    /* .. Termination and release of memory. */
    /* -------------------------------------------------------------------- */
        phase = -1;           /* Release internal memory. */
        PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, &ddum, ia, ja, &idum, &nrhs,
                 iparm, &msglvl, &ddum, &ddum, &error);
    }
    delete[] ia;
    delete[] ja;
    delete[] a;

    /* Matrix data. */
    n = A.num_rows();
    ia = new MKL_INT[n+1];
    ja = new MKL_INT[A.num_nonzeros()];
    a  = new double [A.num_nonzeros()];

    for (size_t row=0; row<A.num_rows()+1; ++row)
        ia[row] = A.row_beg(row)+1;

    for (size_t nz=0; nz<A.num_nonzeros(); ++nz){
        ja[nz] = A.col_ind(nz)+1;
        a[nz]  = A.val(nz);
    }
    mat_version_ = A.Version();

  /* -------------------------------------------------------------------- */
  /* .. Reordering and Symbolic Factorization. This step also allocates */
  /* all memory that is necessary for the factorization. */
  /* -------------------------------------------------------------------- */
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
         &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0)
      {
        printf ("\nERROR during symbolic factorization: %lld", error);
        exit (1);
      }
    printf ("\nReordering completed ... ");
    printf ("\nNumber of nonzeros in factors = %lld", iparm[17]);
    printf ("\nNumber of factorization MFLOPS = %lld", iparm[18]);
  /* -------------------------------------------------------------------- */
  /* .. Numerical factorization. */
  /* -------------------------------------------------------------------- */
    phase = 22;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
         &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0)
      {
        printf ("\nERROR during numerical factorization: %lld", error);
        exit (2);
      }
    printf ("\nFactorization completed ... ");
}

void Solve(const MatrixCL& A, VectorCL &xx, const VectorCL& bb)
{
    if (A.Version() != mat_version_)
        Update(A);
    /* RHS and solution vectors. */
    double *b, *x;
    b = new double[bb.size()];
    x = new double[xx.size()];

    /* Set right hand side to one. */
    for (int i = 0; i < n; i++)
        b[i] = bb[i];

    /* -------------------------------------------------------------------- */
    /* .. Back substitution and iterative refinement. */
    /* -------------------------------------------------------------------- */
    phase = 33;
    printf ("\n\nSolvingsystem...\n");
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
    if (error != 0)
    {
        printf ("\nERROR during solution: %lld", error);
        exit (3);
    }

    for (size_t i= 0; i< xx.size(); ++i)
        xx[i] = x[i];
    delete[] x;
    delete[] b;
}
};
} // end of namespace


#endif /* PARDISOSOLVER_H_ */
