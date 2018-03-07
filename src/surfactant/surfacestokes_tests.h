/// \file surfacestokes_tests.h
/// \brief Some test routines for surfacestokes.cpp
/// \author LNM RWTH Aachen: Thomas Jankuhn

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

#ifndef DROPS_SURFACESTOKES_TESTS_H
#define DROPS_SURFACESTOKES_TESTS_H

#include "surfactant/ifacetransp.h"
#include "surfactant/surfacestokes_funcs.h"
#include "surfactant/surfacestokes_utils.h"

using namespace DROPS;

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Matrix/Vector Tests ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

// Test matrix B with levelset function sphere_2 (??)
void TestB(const MultiGridCL& MG, DROPS::MatDescCL& B, IdxDescCL* ifaceidx, IdxDescCL* ifaceidy)
{
    DROPS::VecDescCL p(ifaceidy);
    DROPS::VecDescCL v(ifaceidx);

    std::cout << "TestB: " << std::endl;

    std::cout << "Test 1: p(x,y,z) = x" << std::endl;
    InitScalar(MG, p, TestB_PressureFun1);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt p^T*B*v_1 = " << dot(p.Data, B.Data*v.Data) << " (8.377580410)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt p^T*B*v_2 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt p^T*B*v_3 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt p^T*B*v_4 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt p^T*B*v_5 = " << dot(p.Data, B.Data*v.Data) << " (-.8377580410)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt p^T*B*v_6 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;

    std::cout << "Test 1: p(x,y,z) = yz" << std::endl;
    InitScalar(MG, p, TestB_PressureFun2);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt p^T*B*v_1 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt p^T*B*v_2 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt p^T*B*v_3 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt p^T*B*v_4 = " << dot(p.Data, B.Data*v.Data) << " (2.513274123)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt p^T*B*v_5 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt p^T*B*v_6 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;

    std::cout << "Test 1: p(x,y,z) = x^2" << std::endl;
    InitScalar(MG, p, TestB_PressureFun3);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt p^T*B*v_1 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt p^T*B*v_2 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt p^T*B*v_3 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt p^T*B*v_4 = " << dot(p.Data, B.Data*v.Data) << " (3.351032164)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt p^T*B*v_5 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt p^T*B*v_6 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;

    std::cout << "Test 1: p(x,y,z) = x(z-y)" << std::endl;
    InitScalar(MG, p, TestB_PressureFun4);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt p^T*B*v_1 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt p^T*B*v_2 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt p^T*B*v_3 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt p^T*B*v_4 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt p^T*B*v_5 = " << dot(p.Data, B.Data*v.Data) << " (2.513274123)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt p^T*B*v_6 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
}

//void TestM(const MultiGridCL& MG, DROPS::MatDescCL& M, IdxDescCL* ifaceidx)
//{
//    DROPS::VecDescCL w(ifaceidx);
//    DROPS::VecDescCL v(ifaceidx);

//    std::cout << "TestM: " << std::endl;

//    std::cout << "Test 1: w(x,y,z) = [1, 0, 0]^T" << std::endl;
//    InitVector(MG, &w, TestRhsVectorFun_v1);
//    InitVector(MG, &v, TestRhsVectorFun_v1);
//    std::cout << "Mit v_1 = [1, 0, 0]^T gilt w^T*M*v_1 = " << dot(w.Data, M.Data*v.Data) << " (8.377580410)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v2);
//    std::cout << "Mit v_2 = [0, 1, 0]^T gilt w^T*M*v_2 = " << dot(w.Data, M.Data*v.Data) << " (0)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v3);
//    std::cout << "Mit v_3 = [0, 0, 1]^T gilt w^T*M*v_3 = " << dot(w.Data, M.Data*v.Data) << " (0)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v4);
//    std::cout << "Mit v_4 = [x, z, 0]^T gilt w^T*M*v_4 = " << dot(w.Data, M.Data*v.Data) << " (0)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v5);
//    std::cout << "Mit v_5 = [0, x*y, x]^T gilt w^T*M*v_5 = " << dot(w.Data, M.Data*v.Data) << " (-.8377580410)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v6);
//    std::cout << "Mit v_6 = [y, z*z, x]^T gilt w^T*M*v_6 = " << dot(w.Data, M.Data*v.Data) << " (0)" << std::endl;

//    std::cout << "Test 2: w(x,y,z) = [x, z, 0]^T" << std::endl;
//    P2Init(TestB_PressureFun2, p, MG, 0.);
//    InitVector(MG, &v, TestRhsVectorFun_v1);
//    std::cout << "Mit v_1 = [1, 0, 0]^T gilt p^T*B*v_1 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v2);
//    std::cout << "Mit v_2 = [0, 1, 0]^T gilt p^T*B*v_2 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v3);
//    std::cout << "Mit v_3 = [0, 0, 1]^T gilt p^T*B*v_3 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v4);
//    std::cout << "Mit v_4 = [x, z, 0]^T gilt p^T*B*v_4 = " << dot(p.Data, B.Data*v.Data) << " (2.513274123)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v5);
//    std::cout << "Mit v_5 = [0, x*y, x]^T gilt p^T*B*v_5 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v6);
//    std::cout << "Mit v_6 = [y, z*z, x]^T gilt p^T*B*v_6 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;

//    std::cout << "Test 3: w(x,y,z) = [0, x*y, x]^T" << std::endl;
//    P2Init(TestB_PressureFun3, p, MG, 0.);
//    InitVector(MG, &v, TestRhsVectorFun_v1);
//    std::cout << "Mit v_1 = [1, 0, 0]^T gilt p^T*B*v_1 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v2);
//    std::cout << "Mit v_2 = [0, 1, 0]^T gilt p^T*B*v_2 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v3);
//    std::cout << "Mit v_3 = [0, 0, 1]^T gilt p^T*B*v_3 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v4);
//    std::cout << "Mit v_4 = [x, z, 0]^T gilt p^T*B*v_4 = " << dot(p.Data, B.Data*v.Data) << " (3.351032164)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v5);
//    std::cout << "Mit v_5 = [0, x*y, x]^T gilt p^T*B*v_5 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v6);
//    std::cout << "Mit v_6 = [y, z*z, x]^T gilt p^T*B*v_6 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;

//    std::cout << "Test 4: w(x,y,z) = [y, z*z, x]^T" << std::endl;
//    P2Init(TestB_PressureFun4, p, MG, 0.);
//    InitVector(MG, &v, TestRhsVectorFun_v1);
//    std::cout << "Mit v_1 = [1, 0, 0]^T gilt p^T*B*v_1 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v2);
//    std::cout << "Mit v_2 = [0, 1, 0]^T gilt p^T*B*v_2 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v3);
//    std::cout << "Mit v_3 = [0, 0, 1]^T gilt p^T*B*v_3 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v4);
//    std::cout << "Mit v_4 = [x, z, 0]^T gilt p^T*B*v_4 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v5);
//    std::cout << "Mit v_5 = [0, x*y, x]^T gilt p^T*B*v_5 = " << dot(p.Data, B.Data*v.Data) << " (2.513274123)" << std::endl;
//    InitVector(MG, &v, TestRhsVectorFun_v6);
//    std::cout << "Mit v_6 = [y, z*z, x]^T gilt p^T*B*v_6 = " << dot(p.Data, B.Data*v.Data) << " (0)" << std::endl;
//}

void TestL(const MultiGridCL& MG, DROPS::MatDescCL& L, IdxDescCL* ifaceidx, IdxDescCL* ifaceidy)
{
    //for the sphere_2

    DROPS::VecDescCL q(ifaceidy);
    DROPS::VecDescCL v(ifaceidx);

    std::cout << "TestL: " << std::endl;

    std::cout << "Test 1: q(x,y,z) = x" << std::endl;
    InitScalar(MG, q, TestL_LagrangeFun1);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt q^T*L*v_1 = " << dot(q.Data, L.Data*v.Data) << " (4.188790205)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt q^T*L*v_2 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt q^T*L*v_3 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt q^T*L*v_4 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt q^T*L*v_5 = " << dot(q.Data, L.Data*v.Data) << " (0.8377580410)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt q^T*L*v_6 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;

    std::cout << "Test 1: q(x,y,z) = yz" << std::endl;
    InitScalar(MG, q, TestL_LagrangeFun2);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt q^T*L*v_1 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt q^T*L*v_2 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt q^T*L*v_3 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt q^T*L*v_4 = " << dot(q.Data, L.Data*v.Data) << " (0.8377580410)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt q^T*L*v_5 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt q^T*L*v_6 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;

    std::cout << "Test 1: q(x,y,z) = x^2" << std::endl;
    InitScalar(MG, q, TestL_LagrangeFun3);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt q^T*L*v_1 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt q^T*L*v_2 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt q^T*L*v_3 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt q^T*L*v_4 = " << dot(q.Data, L.Data*v.Data) << " (2.513274123)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt q^T*L*v_5 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt q^T*L*v_6 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;

    std::cout << "Test 1: q(x,y,z) = x(z-y)" << std::endl;
    InitScalar(MG, q, TestL_LagrangeFun4);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt q^T*L*v_1 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt q^T*L*v_2 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt q^T*L*v_3 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt q^T*L*v_4 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt q^T*L*v_5 = " << dot(q.Data, L.Data*v.Data) << " (0.8377580410)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt q^T*L*v_6 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
}

void TestL2(const MultiGridCL& MG, DROPS::MatDescCL& L, IdxDescCL* ifaceidx, IdxDescCL* ifaceidy)
{
    //for the xy_plane using No_Bnd_Conditions and z = 0.1

    DROPS::VecDescCL q(ifaceidy);
    DROPS::VecDescCL v(ifaceidx);

    std::cout << "TestL: " << std::endl;

    std::cout << "Test 1: q(x,y,z) = x" << std::endl;
    InitScalar(MG, q, TestL_LagrangeFun21);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt q^T*L*v_1 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt q^T*L*v_2 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt q^T*L*v_3 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt q^T*L*v_4 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt q^T*L*v_5 = " << dot(q.Data, L.Data*v.Data) << " (21.33333333)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt q^T*L*v_6 = " << dot(q.Data, L.Data*v.Data) << " (21.33333333)" << std::endl;

    std::cout << "Test 1: q(x,y,z) = y" << std::endl;
    InitScalar(MG, q, TestL_LagrangeFun22);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt q^T*L*v_1 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt q^T*L*v_2 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt q^T*L*v_3 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt q^T*L*v_4 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt q^T*L*v_5 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt q^T*L*v_6 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;

    std::cout << "Test 1: q(x,y,z) = z" << std::endl;
    InitScalar(MG, q, TestL_LagrangeFun23);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt q^T*L*v_1 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt q^T*L*v_2 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt q^T*L*v_3 = " << dot(q.Data, L.Data*v.Data) << " (1.6)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt q^T*L*v_4 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt q^T*L*v_5 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt q^T*L*v_6 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;

    std::cout << "Test 1: q(x,y,z) = 1" << std::endl;
    InitScalar(MG, q, TestL_LagrangeFun24);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt q^T*L*v_1 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt q^T*L*v_2 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt q^T*L*v_3 = " << dot(q.Data, L.Data*v.Data) << " (16)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt q^T*L*v_4 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt q^T*L*v_5 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt q^T*L*v_6 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
}

void TestL_stab(const MultiGridCL& MG, DROPS::MatDescCL& L, IdxDescCL* ifaceidx, IdxDescCL* ifaceidy)
{
    //for the xy_plane using No_Bnd_Conditions integrated on [-2,2]x[-2,2]x[0,1].

    DROPS::VecDescCL q(ifaceidy);
    DROPS::VecDescCL v(ifaceidx);

    std::cout << "TestL_stab: " << std::endl;

    std::cout << "Test 1: q(x,y,z) = z" << std::endl;
    InitScalar(MG, q, TestL_stab_LagrangeFun1);
    InitVector(MG, v, TestL_stabVectorFun_v1);
    std::cout << "Mit v_1 = [0, 0, z]^T gilt q^T*L*v_1 = " << dot(q.Data, L.Data*v.Data) << " (16)" << std::endl;
    InitVector(MG, v, TestL_stabVectorFun_v2);
    std::cout << "Mit v_2 = [0, 0, x*z]^T gilt q^T*L*v_2 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestL_stabVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, y*z]^T gilt q^T*L*v_3 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestL_stabVectorFun_v4);
    std::cout << "Mit v_4 = [0, 0, z*z]^T gilt q^T*L*v_4 = " << dot(q.Data, L.Data*v.Data) << " (16)" << std::endl;

    std::cout << "Test 2: q(x,y,z) = x*z" << std::endl;
    InitScalar(MG, q, TestL_stab_LagrangeFun2);
    InitVector(MG, v, TestL_stabVectorFun_v1);
    std::cout << "Mit v_1 = [0, 0, z]^T gilt q^T*L*v_1 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestL_stabVectorFun_v2);
    std::cout << "Mit v_2 = [0, 0, x*z]^T gilt q^T*L*v_2 = " << dot(q.Data, L.Data*v.Data) << " (21.33333333)" << std::endl;
    InitVector(MG, v, TestL_stabVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, y*z]^T gilt q^T*L*v_3 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestL_stabVectorFun_v4);
    std::cout << "Mit v_4 = [0, 0, z*z]^T gilt q^T*L*v_4 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;

    std::cout << "Test 3: q(x,y,z) = y*z" << std::endl;
    InitScalar(MG, q, TestL_stab_LagrangeFun3);
    InitVector(MG, v, TestL_stabVectorFun_v1);
    std::cout << "Mit v_1 = [0, 0, z]^T gilt q^T*L*v_1 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestL_stabVectorFun_v2);
    std::cout << "Mit v_2 = [0, 0, x*z]^T gilt q^T*L*v_2 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestL_stabVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, y*z]^T gilt q^T*L*v_3 = " << dot(q.Data, L.Data*v.Data) << " (21.33333333)" << std::endl;
    InitVector(MG, v, TestL_stabVectorFun_v4);
    std::cout << "Mit v_4 = [0, 0, z*z]^T gilt q^T*L*v_4 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;

    std::cout << "Test 4: q(x,y,z) = z*z" << std::endl;
    InitScalar(MG, q, TestL_stab_LagrangeFun4);
    InitVector(MG, v, TestL_stabVectorFun_v1);
    std::cout << "Mit v_1 = [0, 0, z]^T gilt q^T*L*v_1 = " << dot(q.Data, L.Data*v.Data) << " (16)" << std::endl;
    InitVector(MG, v, TestL_stabVectorFun_v2);
    std::cout << "Mit v_2 = [0, 0, x*z]^T gilt q^T*L*v_2 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestL_stabVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, y*z]^T gilt q^T*L*v_3 = " << dot(q.Data, L.Data*v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestL_stabVectorFun_v4);
    std::cout << "Mit v_4 = [0, 0, z*z]^T gilt q^T*L*v_4 = " << dot(q.Data, L.Data*v.Data) << " (21.33333333)" << std::endl;
}

// Test Setup Rhs with levelset function sphere_2 (??)
void TestRhs(const MultiGridCL& MG, LevelsetP2CL& lset, IdxDescCL *ifaceidx)
{
    DROPS::VecDescCL v(ifaceidx);
    DROPS::VecDescCL rhs(ifaceidx);

    std::cout << "TestRhs: " << std::endl;

    std::cout << "Test 1: f(x,y,z) = [-y, x, 0]^T" << std::endl;
    DROPS::SetupInterfaceVectorRhsP1(MG, &rhs, lset.Phi, lset.GetBndData(), TestRhsVectorFun1);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt rhs^T*v_1 = " << dot(rhs.Data, v.Data) << " (3.456944444*10^(-16))" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt rhs^T*v_2 = " << dot(rhs.Data, v.Data) << " (6.495035630*10^(-13))" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt rhs^T*v_3 = " << dot(rhs.Data, v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt rhs^T*v_4 = " << dot(rhs.Data, v.Data) << " (-8.958868921*10^(-17))" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt rhs^T*v_5 = " << dot(rhs.Data, v.Data) << " (1.325906936*10^(-16))" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt rhs^T*v_6 = " << dot(rhs.Data, v.Data) << "(-4.188790205)" << std::endl;

    std::cout << "Test 2: f(x,y,z) = P [-z, x, 0]^T" << std::endl;
    DROPS::SetupInterfaceVectorRhsP1(MG, &rhs, lset.Phi, lset.GetBndData(), TestRhsVectorFun2);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt rhs^T*v_1 = " << dot(rhs.Data, v.Data) << " (4.371503159*10^(-16))" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt rhs^T*v_2 = " << dot(rhs.Data, v.Data) << " (6.512264686*10^(-13))" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt rhs^T*v_3 = " << dot(rhs.Data, v.Data) << " (1.628003288*10^(-13))" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt rhs^T*v_4 = " << dot(rhs.Data, v.Data) << " (4.163336342*10^(-17))" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt rhs^T*v_5 = " << dot(rhs.Data, v.Data) << " (0.8377580410)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt rhs^T*v_6 = " << dot(rhs.Data, v.Data) << " (2.733212962*10^(-13))" << std::endl;

    std::cout << "Test 3: f(x,y,z) = [x*z, -y, z]^T" << std::endl;
    DROPS::SetupInterfaceVectorRhsP1(MG, &rhs, lset.Phi, lset.GetBndData(), TestRhsVectorFun3);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt rhs^T*v_1 = " << dot(rhs.Data, v.Data) << " (-2.037029639*10^(-29))" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt rhs^T*v_2 = " << dot(rhs.Data, v.Data) << " (3.456944444*10^(-16))" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt rhs^T*v_3 = " << dot(rhs.Data, v.Data) << " (-1.933333537*10^(-16))" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt rhs^T*v_4 = " << dot(rhs.Data, v.Data) << " (-2.567390744*10^(-16))" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt rhs^T*v_5 = " << dot(rhs.Data, v.Data) << " (3.752423719*10^(-16))" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt rhs^T*v_6 = " << dot(rhs.Data, v.Data) << " (-1.179611964*10^(-16))" << std::endl;

    std::cout << "Test 4: f(x,y,z) = P [0, y*z, x]^T" << std::endl;
    DROPS::SetupInterfaceVectorRhsP1(MG, &rhs, lset.Phi, lset.GetBndData(), TestRhsVectorFun4);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt rhs^T*v_1 = " << dot(rhs.Data, v.Data) << " (2.205757161*10^(-17))" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt rhs^T*v_2 = " << dot(rhs.Data, v.Data) << " (-5.540273101*10^(-17))" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt rhs^T*v_3 = " << dot(rhs.Data, v.Data) << " (-.8377580410)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt rhs^T*v_4 = " << dot(rhs.Data, v.Data) << " (-4.033232082*10^(-17))" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt rhs^T*v_5 = " << dot(rhs.Data, v.Data) << " (3.351032164)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt rhs^T*v_6 = " << dot(rhs.Data, v.Data) << " (3.351032164)" << std::endl;
}

void TestA(const MultiGridCL& MG, DROPS::MatDescCL& A, IdxDescCL *ifaceidx)
{
    DROPS::VecDescCL u(ifaceidx);
    DROPS::VecDescCL v(ifaceidx);

    std::cout << "TestA: " << std::endl;

    std::cout << "Test 1: f(x,y,z) = [, , ]^T" << std::endl;
    InitVector(MG, u, TestA_VectorFun1);
    InitVector(MG, v, TestA_VectorFun_v1);
    std::cout << "Mit v_1 = [, , ]^T gilt v_1^T*A*u = " << dot(A.Data*u.Data, v.Data) << " ()" << std::endl;
}

void TestRhsP2(const MultiGridCL& MG, LevelsetP2CL& lset, IdxDescCL *ifaceidx)
{
    //xy-plane mit z = Pi/30
    DROPS::VecDescCL v(ifaceidx);
    DROPS::VecDescCL rhs(ifaceidx);

    std::cout << "TestRhsP2: " << std::endl;

    std::cout << "Test 1: f(x,y,z) = [-y, x, 0]^T" << std::endl;
    DROPS::SetupInterfaceVectorRhsP2(MG, &rhs, lset.Phi, lset.GetBndData(), TestRhsVectorFun1);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt rhs^T*v_1 = " << dot(rhs.Data, v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt rhs^T*v_2 = " << dot(rhs.Data, v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt rhs^T*v_3 = " << dot(rhs.Data, v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt rhs^T*v_4 = " << dot(rhs.Data, v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt rhs^T*v_5 = " << dot(rhs.Data, v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt rhs^T*v_6 = " << dot(rhs.Data, v.Data) << "(−21.333333333)" << std::endl;

    std::cout << "Test 2: f(x,y,z) = P [-z, x, 0]^T" << std::endl;
    DROPS::SetupInterfaceVectorRhsP2(MG, &rhs, lset.Phi, lset.GetBndData(), TestRhsVectorFun2);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt rhs^T*v_1 = " << dot(rhs.Data, v.Data) << " (~)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt rhs^T*v_2 = " << dot(rhs.Data, v.Data) << " (~)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt rhs^T*v_3 = " << dot(rhs.Data, v.Data) << " (~)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt rhs^T*v_4 = " << dot(rhs.Data, v.Data) << " (~)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt rhs^T*v_5 = " << dot(rhs.Data, v.Data) << " (~)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt rhs^T*v_6 = " << dot(rhs.Data, v.Data) << " (~)" << std::endl;

    std::cout << "Test 3: f(x,y,z) = [x*z, -y, z]^T" << std::endl;
    DROPS::SetupInterfaceVectorRhsP2(MG, &rhs, lset.Phi, lset.GetBndData(), TestRhsVectorFun3);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt rhs^T*v_1 = " << dot(rhs.Data, v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt rhs^T*v_2 = " << dot(rhs.Data, v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt rhs^T*v_3 = " << dot(rhs.Data, v.Data) << " (1.675516082)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt rhs^T*v_4 = " << dot(rhs.Data, v.Data) << " (2.234021443)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt rhs^T*v_5 = " << dot(rhs.Data, v.Data) << " (0)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt rhs^T*v_6 = " << dot(rhs.Data, v.Data) << " (0)" << std::endl;

    std::cout << "Test 4: f(x,y,z) = P [0, y*z, x]^T" << std::endl;
    DROPS::SetupInterfaceVectorRhsP2(MG, &rhs, lset.Phi, lset.GetBndData(), TestRhsVectorFun4);
    InitVector(MG, v, TestRhsVectorFun_v1);
    std::cout << "Mit v_1 = [1, 0, 0]^T gilt rhs^T*v_1 = " << dot(rhs.Data, v.Data) << " (~)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v2);
    std::cout << "Mit v_2 = [0, 1, 0]^T gilt rhs^T*v_2 = " << dot(rhs.Data, v.Data) << " (~)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v3);
    std::cout << "Mit v_3 = [0, 0, 1]^T gilt rhs^T*v_3 = " << dot(rhs.Data, v.Data) << " (~)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v4);
    std::cout << "Mit v_4 = [x, z, 0]^T gilt rhs^T*v_4 = " << dot(rhs.Data, v.Data) << " (~)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v5);
    std::cout << "Mit v_5 = [0, x*y, x]^T gilt rhs^T*v_5 = " << dot(rhs.Data, v.Data) << " (~)" << std::endl;
    InitVector(MG, v, TestRhsVectorFun_v6);
    std::cout << "Mit v_6 = [y, z*z, x]^T gilt rhs^T*v_6 = " << dot(rhs.Data, v.Data) << " (~)" << std::endl;
}

bool CheckValuesSym(const MatrixCL& AP1, const MatrixCL& AP2, const VectorCL& aP1, const VectorCL& bP1, const VectorCL& aP2, const VectorCL& bP2)
{
    double intP1 = dot(AP1*bP1, aP1);
    double intP2 = dot(AP2*bP2, aP2);

//    std::cout << "intP1 is: " << intP1 << ",  " << "intP2 is: " << intP2 << std::endl;

    return std::abs(intP1-intP2) < 1e-14;
}

bool CheckValuesAsym(const MatrixCL& AP1, const MatrixCL& AP2, const VectorCL& aP1, const VectorCL& cP1, const VectorCL& aP2)
{
    bool equal = 0;
    double intP1 = dot(AP1*aP1, cP1);
    double intP2 = dot(AP2*aP2, cP1);

    if( std::abs(intP1-intP2) < 1e-14) {
        equal = 1;
    }
    return equal;
}

void TestAllP2MatricesWithP1Matrices(const MultiGridCL& MG, LevelsetP2CL& lset, IdxDescCL& ifaceVecP2idx, IdxDescCL& ifaceVecP1idx, IdxDescCL& ifaceP1idx, bool fullgrad)
{
    DROPS::MatDescCL A_P1, A_P1_stab, B_P1P1, M_P1, S_P1, L_P1P1, L_P1P1_stab, Schur_P1, Schur_P1_stab;

    A_P1.SetIdx( &ifaceVecP1idx, &ifaceVecP1idx);
    A_P1_stab.SetIdx( &ifaceVecP1idx, &ifaceVecP1idx);
    B_P1P1.SetIdx( &ifaceP1idx, &ifaceVecP1idx);
    M_P1.SetIdx( &ifaceVecP1idx, &ifaceVecP1idx);
    S_P1.SetIdx( &ifaceVecP1idx, &ifaceVecP1idx);
    L_P1P1.SetIdx( &ifaceP1idx, &ifaceVecP1idx);
    L_P1P1_stab.SetIdx( &ifaceP1idx, &ifaceVecP1idx);
    Schur_P1.SetIdx(&ifaceP1idx, &ifaceP1idx);
    Schur_P1_stab.SetIdx(&ifaceP1idx, &ifaceP1idx);

    SetupStokesIF_P1P1(MG, &A_P1, &A_P1_stab, &B_P1P1, &M_P1, &S_P1, &L_P1P1, &L_P1P1_stab, &Schur_P1, &Schur_P1_stab, lset.Phi, lset.GetBndData(), fullgrad);

    DROPS::MatDescCL A_P2, A_P2_stab, B_P1P2, M_P2, S_P2, L_P1P2, L_P1P2_stab, Schur_P2, Schur_P2_stab;

    A_P2.SetIdx( &ifaceVecP2idx, &ifaceVecP2idx);
    A_P2_stab.SetIdx( &ifaceVecP2idx, &ifaceVecP2idx);
    B_P1P2.SetIdx( &ifaceP1idx, &ifaceVecP2idx);
    M_P2.SetIdx( &ifaceVecP2idx, &ifaceVecP2idx);
    S_P2.SetIdx( &ifaceVecP2idx, &ifaceVecP2idx);
    L_P1P2.SetIdx( &ifaceP1idx, &ifaceVecP2idx);
    L_P1P2_stab.SetIdx( &ifaceP1idx, &ifaceVecP2idx);
    Schur_P2.SetIdx(&ifaceP1idx, &ifaceP1idx);
    Schur_P2_stab.SetIdx(&ifaceP1idx, &ifaceP1idx);

    SetupStokesIF_P2P1(MG, &A_P2, &A_P2_stab, &B_P1P2, &M_P2, &S_P2, &L_P1P2, &L_P1P2_stab, &Schur_P2, &Schur_P2_stab, lset.Phi, lset.GetBndData(), fullgrad);

    DROPS::VecDescCL a, b, c, a_P2, b_P2;
    a.SetIdx( &ifaceVecP1idx);
    size_t n = a.Data.size();
    std:: cout << "The outer loop size is: " << n << std::endl;
    c.SetIdx( &ifaceP1idx);
    size_t m = c.Data.size();
    std:: cout << "The second inner loop size is: " << m << std::endl;


    for (size_t i=0; i<n; i++) {
        std::cout << "The outer loop counter is: " << i << std::endl;
        a.SetIdx( &ifaceVecP1idx);
        a.Data[i] = 1.;
        std::cout << "The norm of a is: " << norm(a.Data) << std::endl;
        a_P2.SetIdx( &ifaceVecP2idx);
        InterpolateP2Vec(MG, a, a_P2);
        std::cout << "The norm of a_P2 is: " << norm(a_P2.Data) << std::endl;
        for (size_t j=0; j<n; j++) {
            b.SetIdx( &ifaceVecP1idx);
            b.Data[j] = 1.;
            b_P2.SetIdx( &ifaceVecP2idx);
            InterpolateP2Vec(MG, b, b_P2);

            if(!CheckValuesSym(A_P1.Data, A_P2.Data, a.Data, b.Data, a_P2.Data, b_P2.Data)) {
                std::cout << "Something is wrong with A!" << std::endl;
            }
            if(!CheckValuesSym(A_P1_stab.Data, A_P2_stab.Data, a.Data, b.Data, a_P2.Data, b_P2.Data)) {
                std::cout << "Something is wrong with A_stab!" << std::endl;
            }
            if(!CheckValuesSym(M_P1.Data, M_P2.Data, a.Data, b.Data, a_P2.Data, b_P2.Data)) {
                std::cout << "Something is wrong with M!" << std::endl;
            }
            if(!CheckValuesSym(S_P1.Data, S_P2.Data, a.Data, b.Data, a_P2.Data, b_P2.Data)) {
                std::cout << "Something is wrong with S!" << std::endl;
            }
        }
        for (size_t k=0; k<m; k++) {
            c.SetIdx( &ifaceP1idx);
            c.Data[k] = 1.;

            if(!CheckValuesAsym(L_P1P1.Data, L_P1P2.Data, a.Data, c.Data, a_P2.Data)) {
               std::cout << "Something is wrong with L!" << std::endl;
            }
            if(!CheckValuesAsym(L_P1P1_stab.Data, L_P1P2_stab.Data, a.Data, c.Data, a_P2.Data)) {
               std::cout << "Something is wrong with L_stab!" << std::endl;
            }
        }
    }
}

void TestP2Matrices(const MultiGridCL& MG, LevelsetP2CL& lset, IdxDescCL &ifaceVecP2idx, DROPS::IdxDescCL &vecP2idx, BndDataCL<Point3DCL> bndvec, const MatrixCL A, const MatrixCL M)
{
    //For matrix A
    DROPS::VecDescCL v_A, rhs_A, v_Axtent, res_A;
    v_A.SetIdx( &ifaceVecP2idx);
    rhs_A.SetIdx( &ifaceVecP2idx);
    InitVector(MG, v_A, TestP2Matrices_A_Vector_Fun_1);
    DROPS::SetupInterfaceVectorRhsP2(MG, &rhs_A, lset.Phi, lset.GetBndData(), TestP2Matrices_A_Rhs_1);
    std::cout << "norm(A*v_A - rhs_A) is: " << norm(A*v_A.Data - rhs_A.Data) << std::endl;

    res_A.SetIdx( &ifaceVecP2idx);
    v_Axtent.SetIdx( &vecP2idx);
    res_A.Data = A*v_A.Data - rhs_A.Data;
    Extend(MG, res_A, v_Axtent);
    std::cout << "The L2-Norm of A*v_A - rhs_A is: " << L2_Vector_error(MG, lset.Phi, lset.GetBndData(), make_P2Eval(MG, bndvec, v_Axtent), ZeroVectorFun) << std::endl;

    //For matrix M
    DROPS::VecDescCL v_M, rhs_M, v_Mxtent, res_M;
    v_M.SetIdx( &ifaceVecP2idx);
    rhs_M.SetIdx( &ifaceVecP2idx);
    InitVector(MG, v_M, TestP2Matrices_A_Vector_Fun_1);
    DROPS::SetupInterfaceVectorRhsP2(MG, &rhs_M, lset.Phi, lset.GetBndData(), TestP2Matrices_A_Vector_Fun_1);
    std::cout << "norm(M*v_M - rhs_M) is: " << norm(M*v_M.Data - rhs_M.Data) << std::endl;

    res_M.SetIdx( &ifaceVecP2idx);
    v_Mxtent.SetIdx( &vecP2idx);
    res_M.Data = M*v_M.Data - rhs_M.Data;
    Extend(MG, res_M, v_Mxtent);
    std::cout << "The L2-Norm of M*v_M - rhs_M is: " << L2_Vector_error(MG, lset.Phi, lset.GetBndData(), make_P2Eval(MG, bndvec, v_Mxtent), ZeroVectorFun) << std::endl;

/*    //For matrix S
    DROPS::VecDescCL v_S, rhs_S;
    v_S.SetIdx( &ifaceVecP2idx);
    rhs_S.SetIdx( &ifaceVecP2idx);
    InitVector(MG, v_S, TestP2Matrices_S_Vector_Fun_1);
    DROPS::SetupInterfaceVectorRhsP2(MG, &rhs_S, lset.Phi, lset.GetBndData(), TestP2Matrices_S_Rhs_1);
    std::cout << "norm(S*v_S - rhs_S) is: " << norm(S*v_S.Data - rhs_S.Data) << std::endl;

    //For matrix A_stab
    DROPS::VecDescCL v_A_stab, rhs_A_stab;
    v_A_stab.SetIdx( &ifaceVecP2idx);
    rhs_A_stab.SetIdx( &ifaceVecP2idx);
    InitVector(MG, v_A_stab, TestP2Matrices_A_stab_Vector_Fun_1);
    DROPS::SetupInterfaceVectorRhsP2(MG, &rhs_A_stab, lset.Phi, lset.GetBndData(), TestP2Matrices_A_stab_Rhs_1);
    std::cout << "norm(A_stab*v_A_stab - rhs_A_stab) is: " << norm(A_stab*v_A_stab.Data - rhs_A_stab.Data) << std::endl;
*/

}

#endif
