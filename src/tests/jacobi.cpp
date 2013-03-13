/// \file jacobi.cpp
/// \brief Tests Jacobi eigenvalue-solver for small dense matrices.
/// \author LNM RWTH Aachen: Joerg Grande; SC RWTH Aachen:

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
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#include "misc/container.h"

const double dat[]= { 1., 1., 1., 1.,
                      1., 2., 3., 4.,
                      1., 3., 6., 10.,
                      1., 4., 10., 20. };

using namespace DROPS;

int Test1()
{
    SMatrixCL<4,4> M(dat+0);
    std::cout<<"Data: Matrix M: " << M << std::endl;
    cyclic_jacobi( M, 1e-10);
    std::cout <<"Eigenvalues: "<< M << std::endl;

    SMatrixCL<4,4> M2(dat+0), Q;
    std::cout<<"Data: Matrix M2: " << M2 << std::endl;
    cyclic_jacobi( M2, 1e-10, &Q);
    std::cout <<"Eigenvalues: "<< M2 << std::endl;
    SMatrixCL<4,4> M3(dat+0);
    std::cout << "Eigenvectors: " << Q << "\n"
              << "Orthogonality: " << GramMatrix( Q) << "\n"
              << "Eigenvector-property: " << M3*Q - Q*M2 << "\n";

    return 0;
}

int Test2()
{
    DMatrixCL<double> D( 4, 4);
    std::copy( dat+0, dat+16, D.GetCol( 0));
    std::cout<<"Data: Matrix D: ";
    seq_out ( D.GetCol( 0), D.GetCol( 0) + 16, std::cout);
    cyclic_jacobi( D, 1e-10);
    std::cout <<"Eigenvalues: ";
    seq_out ( D.GetCol( 0), D.GetCol( 0) + 16, std::cout);

    return 0;
}


int main()
{
    return Test1()+Test2();
}

