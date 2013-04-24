/// \file minres.cpp
/// \brief tests implementation of Minres and Lanczos algorithms
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande; SC RWTH Aachen:

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

#include "num/krylovsolver.h"
#include "num/precond.h"
#include <iostream>

class TrivialPreCL
{
  public:
    TrivialPreCL() {}

    template <typename Mat, typename Vec>
    void
    Apply(const Mat& K, Vec& x, const Vec& b) const {
        x= b; }
};

int TestLanczos()
{
    std::cout << "TEST1 - METHOD: Lanczos" << std::endl;
    DROPS::MatrixCL A;
    DROPS::MatrixBuilderCL AB(&A, 2, 2);
    AB( 0, 1)= 1.; AB( 1, 0)= 1.;
    AB.Build();
    DROPS::VectorCL r0( 1., 2);
    r0[0]= 1.; r0[1]= 2.;

    std::cout << "A = " << A << '\n' << "r0 = " << r0 << std::endl;
    DROPS::LanczosONBCL<DROPS::VectorCL> onb;
    onb.new_basis(A, r0, DROPS::DummyExchangeCL());
    std::vector<DROPS::VectorCL> basis;
    do {
        basis.push_back( onb.q[0]);
        std::cout << onb.q[0] << std::endl;
        std::cout << "a " << onb.a0 << " b " << onb.b[0] << std::endl;
    } while (onb.next( A, DROPS::DummyExchangeCL()));
    basis.push_back( onb.q[0]);
    std::cout << onb.q[0] << std::endl;
    std::cout << "a " << onb.a0 << " b " << onb.b[0] << std::endl;
    std::cout << "lucky breakdown after " << basis.size() << " vectors."
              <<std::endl;
    for (unsigned int i= 0; i<basis.size(); ++i) {
        for (unsigned int j= 0; j<basis.size(); ++j) {
            std::cout << dot( basis[i], basis[j]) << '\t';
        }
        std::cout << std::endl;
    }
    std::cout  << '\n';
    return 0;
}

int TestMinres()
{
    std::cout << '\n' << "TEST2 - METHOD: Minres, 2x2 example" << std::endl;
    DROPS::MatrixCL A;
    DROPS::MatrixBuilderCL AB(&A, 2, 2);
    AB( 0, 1)= 1.; AB( 1, 0)= 1.;
    AB.Build();
    DROPS::VectorCL b( 1., 2);
    b[0]= 1.; b[1]= 2.;
    DROPS::VectorCL x( 0., 2);

    std::cout << "A =" << A  << '\n' << "b = " << b << std::endl;
    int mi= 10;
    double tol= 1e-10;
    std::cout << "Minres: " ;
    MINRES( A, x, b, DROPS::DummyExchangeCL(), mi, tol);
    std::cout << "TEST2 x: " << x << "TEST2 Ax-b: " << DROPS::VectorCL( A*x - b) << '\n'  << "TEST2 Iterationnumber: " << mi
              << '\n'<< "TEST2 Residual: "  << tol << '\n'  << std::endl;
    return 0;
}

int TestMinres2()
{
    std::cout << '\n' << "TEST3 - METHOD: Minres, 4x4 example" << std::endl;
    DROPS::MatrixCL A;
    DROPS::MatrixBuilderCL AB(&A, 4, 4);
    AB( 0, 0)= -249.;AB( 1, 0)= -453.;AB( 2, 0)= -397.;AB( 3, 0)= -52.;
    AB( 0, 1)= -453.;AB( 1, 1)= -731.;AB( 2, 1)= -601.;AB( 3, 1)= -78.;
    AB( 0, 2)= -397.;AB( 1, 2)= -601.;AB( 2, 2)= -648.;AB( 3, 2)= -91.;
    AB( 0, 3)=  -52.;AB( 1, 3)=  -78.;AB( 2, 3)=  -91.;AB( 3, 3)= -13.;
    AB.Build();
    DROPS::VectorCL b( 0., 4);
    b[0]= -1277./2/1.;
    b[1]= -2015./2.;
    b[2]= -3907./4.;
    b[3]= -533./4.;
    DROPS::VectorCL x( 0., 4);

    std::cout << "A = " << A  << '\n' << "b = " << b << std::endl;
    int mi= 10;
    double tol= 1e-10;
    std::cout << "Minres2: "; 
    MINRES( A, x, b, DROPS::DummyExchangeCL(), mi, tol);
    std::cout << "TEST3 x: " << x << "TEST3 Ax-b: " << DROPS::VectorCL( A*x - b) << '\n'  << "TEST3 Iterationnumber: " << mi
              << '\n'<< "TEST3 Residual: "  << tol << '\n'  << std::endl;
    return 0;
}

int TestPMinres()
{
    std::cout << '\n' << "TEST4 - METHOD: PMinres, 4x4 example" << std::endl;
    DROPS::MatrixCL A;
    DROPS::MatrixBuilderCL AB(&A, 4, 4);
    AB( 0, 0)= -249.;AB( 1, 0)= -453.;AB( 2, 0)= -397.;AB( 3, 0)= -52.;
    AB( 0, 1)= -453.;AB( 1, 1)= -731.;AB( 2, 1)= -601.;AB( 3, 1)= -78.;
    AB( 0, 2)= -397.;AB( 1, 2)= -601.;AB( 2, 2)= -648.;AB( 3, 2)= -91.;
    AB( 0, 3)=  -52.;AB( 1, 3)=  -78.;AB( 2, 3)=  -91.;AB( 3, 3)= -13.;
    AB.Build();
    DROPS::VectorCL b( 0., 4);
    b[0]= -1277./2/1.;
    b[1]= -2015./2.;
    b[2]= -3907./4.;
    b[3]= -533./4.;
    DROPS::VectorCL x( 0., 4);

    std::cout << "A = " << A  << '\n' << "b = " << b << std::endl;
    int mi= 10;
    double tol= 1e-10;
    DROPS::LanczosONBCL<DROPS::VectorCL> q;
    q.new_basis( A, b, DROPS::DummyExchangeCL());
    DROPS::PMResSolverCL<DROPS::LanczosONBCL<DROPS::VectorCL> > pmr( q, mi, tol);
    std::cout << "PMinres: "; 
    pmr.Solve( A, x, b, DROPS::DummyExchangeCL());
    std::cout << "TEST4 x: " << x << "TEST4 Ax-b: " << DROPS::VectorCL( A*x - b) << '\n' << "TEST4 Iterationnumber: " << pmr.GetIter() 
              << '\n' << "TEST4 Residual: " << pmr.GetResid() << '\n' << std::endl;
    return 0;
}

int TestPMinres2()
{
    std::cout << '\n' << "TEST5 - METHOD: PMinres2, 4x4 example" << std::endl;
    DROPS::MatrixCL A;
    DROPS::MatrixBuilderCL AB(&A, 4, 4);
    AB( 0, 0)= -249.;AB( 1, 0)= -453.;AB( 2, 0)= -397.;AB( 3, 0)= -52.;
    AB( 0, 1)= -453.;AB( 1, 1)= -731.;AB( 2, 1)= -601.;AB( 3, 1)= -78.;
    AB( 0, 2)= -397.;AB( 1, 2)= -601.;AB( 2, 2)= -648.;AB( 3, 2)= -91.;
    AB( 0, 3)=  -52.;AB( 1, 3)=  -78.;AB( 2, 3)=  -91.;AB( 3, 3)= -13.;
    AB.Build();
    DROPS::VectorCL b( 0., 4);
    b[0]= -1277./2/1.;
    b[1]= -2015./2.;
    b[2]= -3907./4.;
    b[3]= -533./4.;
    DROPS::VectorCL x( 0., 4);

    std::cout << "A = " << '\n' << A << "b = " << b << std::endl;
    int mi= 10;
    double tol= 1e-10;
    DROPS::DummyPcCL pc;
    DROPS::PLanczosONBCL<DROPS::VectorCL, DROPS::DummyPcCL> q( pc);
    DROPS::PMResSolverCL<DROPS::PLanczosONBCL<DROPS::VectorCL, DROPS::DummyPcCL> > pmr( q, mi, tol);
    std::cout << "PMinres2: "; 
    pmr.Solve( A, x, b, DROPS::DummyExchangeCL());
    std::cout << "TEST5 x: " << x << "TEST5 Ax-b: " << DROPS::VectorCL( A*x - b) << '\n' << "TEST5 Iterationnumber: " << pmr.GetIter() 
              << '\n' << "TEST5 Residual: " << pmr.GetResid() << '\n' << std::endl;
    return 0;
}

int main (int, char**)
{
  try {
    return TestLanczos() + TestMinres() + TestMinres2()
           + TestPMinres() + TestPMinres2();
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
