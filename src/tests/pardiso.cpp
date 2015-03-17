/// \file pardiso.cpp
/// \brief tests pardiso interface
/// \author LNM RWTH Aachen: Patrick Esser (Christoph Lehrenfeld);

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
 * Copyright 2013 LNM/SC RWTH Aachen, Germany
*/

#include "misc/container.h"
#include "num/spmat.h"
#ifdef DROPS_PARDISO
#include "num/pardisosolver.h"
#endif
#include <iostream>
#include <math.h>



int PardisoTest()
{
    std::cout << "PARDISO: 500x500:\n" << std::endl;
    DROPS::MatrixCL A;
    const int N = 500;
    DROPS::MatrixBuilderCL AB(&A, N*N, N*N);
    DROPS::VectorCL b(N*N);
    DROPS::VectorCL x(N*N);
    const double h = 1.0/N;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            AB(i*N+j,i*N+j) += 4.0/(h*h);
            if (i>0)
                AB(i*N+j,(i-1)*N+j) -= 1.0/(h*h) + 1.0/h;
            if (i<N-1)
                AB(i*N+j,(i+1)*N+j) -= 1.0/(h*h) + 1.0/h;
            if (j>0)
                AB(i*N+j,i*N+j-1) -= 1.0/(h*h) - 1.0/h;
            if (j<N-1)
                AB(i*N+j,i*N+j+1) -= 1.0/(h*h) - 1.0/h;
            b[i*N+j] += 1.0;
            
        }
    }
    AB.Build();
    
#ifdef DROPS_PARDISO
    DROPS::PardisoSolverCL SolveA( A);
    SolveA.Solve(A, x, b);
#endif
    DROPS::VectorCL r(A*x - b);
    std::cout << "PARDISO 1 Residual:  "<< norm(r) << std::endl;
    return 0;
}


int main (int, char**)
{
#ifdef DROPS_PARDISO
  try {
      return PardisoTest();
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
#endif
}
