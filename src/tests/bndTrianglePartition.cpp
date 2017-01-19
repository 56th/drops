/// \file bndTrianglePartition.cpp
/// \brief tests for refTrianglePartionCL, BndTrainglePartitionCL, and related quadrature class and functions
/// \author LNM RWTH Aachen: Liang Zhang; SC RWTH Aachen:

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

#include "geom/principallattice.h"
#include "geom/subtriangulation.h"
#include "num/quadrature.h"
#include "misc/container.h"
#include "num/discretize.h"
#include "num/lattice-eval.h"
#include "geom/multigrid.h"

#include <iostream>
#include <sstream>


void test_triangle_partition()
{
    std::cout<<"=========================TrianglePartitionCL test: \n"
             <<"all 26 level-set sign patterns (except 0, 0, 0) will be prescribed, then the corresponding reference triangle partition info will be given"<<std::endl;
    DROPS::byte VertexNum =0;
    bool consistency= true;
    DROPS::byte ls[4];
    DROPS::GridFunctionCL<> ls_value(4);
    ls[VertexNum] = 0;
    ls_value[VertexNum] = 0;
    int c= 0;
    DROPS::BndTriangPartitionCL BndTri;
    for (int i= -1; i <= 1; ++i)
      for (int j= -1; j <= 1; ++j)
        for (int k= -1; k <= 1; ++k) {
              if (i == 0 && j == 0 && k == 0) continue;
              ls[ (VertexNum +1)%4 ]= i; ls[(VertexNum +2)%4]= j; ls[(VertexNum +3)%4]= k;
              ls_value[ (VertexNum +1)%4 ]= i; ls_value[(VertexNum +2)%4]= j; ls_value[(VertexNum +3)%4]= k;
              DROPS::RefTrianglePartitionCL tri (ls, VertexNum);
              std::cout <<"RefTrianglePartitionCL====================================================="<<std::endl;
              std::cout << "case: " << c << " ls: " << int(ls[0]) << ' ' << int(ls[1]) << ' ' << int(ls[2]) << ' ' << int(ls[3]) << std::endl;
              std::cout << "The triangle is cut to " << tri.size() << " Sub-triangles;"<<std::endl;
              for (DROPS::RefTrianglePartitionCL::const_triangle_iterator it= tri.triangle_begin(), end= tri.triangle_end(); it != end; ++it)
              std::cout << "Sign of the triangle: " <<tri.sign(it)<< " Indices of vertices: "<<int((*it)[0])<<" "<< int((*it)[1])<< " "<<int((*it)[2]) <<std::endl;
              std::cout <<"TrianglePartitionCL------------------------------------------------------"<<std::endl; 
              BndTri.make_partition2D<DROPS::SortedVertexPolicyCL, DROPS::MergeCutPolicyCL> ( DROPS::PrincipalLatticeCL::instance(1), int(VertexNum), ls_value);
              std::cout<< "Size of negative triangle(s): "<< BndTri.triangle_size(DROPS::NegTetraC)<<std::endl;
              for (DROPS::BndTriangPartitionCL::const_triangle_iterator it= BndTri.triangle_begin(), end= BndTri.triangle_end(); it != end; ++it)
              std::cout << " Indices of vertices: "<<int((*it)[0])<<" "<< int((*it)[1])<< " "<<int((*it)[2]) <<std::endl;
              if ( BndTri.vertex_size() != (4 + tri.size() -1) || BndTri.triangle_size() != tri.size() )
                  consistency=false;
              c++;
          }
    std::cout <<"TrianglePartitionCL------------------------------------------------------"<<std::endl;
    if(!consistency)
        std::cout<<"The vertex size and the triangle size of BndTrianglePartitionCL is NOT consistent with Ref..."<<std::endl;
    else
        std::cout<<"The vertex size and the triangle size of BndTrianglePartitionCL is consistent with Ref..."<<std::endl;
}

inline double cylinder_instat (const DROPS::Point3DCL& p, double)
{
    return p[0]*p[0] + p[2]*p[2] -0.5*0.5;
}

inline double ball_instat (const DROPS::Point3DCL& p, double)
{
    return p.norm() - 0.5;
}

void test_bnd_integral()
{
    std::cout<<"=========================cut boundary integral test: \n"
             <<"The cut boundary x=0 is integrated separately by positive part and negative part"<<std::endl;
    DROPS::Uint num_sub = 32;
    DROPS::Uint num_sub_lattice = 2;
    DROPS::Point3DCL orig;
    orig[0] = -2;
    orig[2] = -1;
    orig[1] = -1;
    // [-2, 0] x[-1, 1] x[-1, 1] brick
    DROPS::BrickBuilderCL brick(orig, 2.*DROPS::std_basis<3>(1), 2.*DROPS::std_basis<3>(2), 2.*DROPS::std_basis<3>(3), num_sub, num_sub, num_sub);
    DROPS::MultiGridCL mg( brick);
    const DROPS::PrincipalLatticeCL& lat= DROPS::PrincipalLatticeCL::instance( num_sub_lattice);
    DROPS::GridFunctionCL<> ls( lat.vertex_size());
    DROPS::BndTriangPartitionCL BndTri;
    double area_neg= 0.;
    double area_pos= 0.;
    DROPS::QuadDomainCL qdom;
    bool onbnd = false;
    DROPS_FOR_TRIANG_TETRA( mg, 0, it) {
        evaluate_on_vertexes( ball_instat, *it, lat, 0., Addr( ls));
        DROPS::Uint facenum =3;
        onbnd = (*it).IsBndSeg(facenum);   //it seems segment 3 is the face number for all bnd faces on x=0
        if(onbnd)
        {
            DROPS::Point3DCL normal;
            (*it).GetOuterNormal(facenum, normal);
            if(normal[0]==1){
                BndTri.make_partition2D<DROPS::PartitionedVertexPolicyCL, DROPS::MergeCutPolicyCL>( lat, facenum, ls);
                DROPS::make_CompositeQuad5BndDomain2D(qdom, BndTri, *it); 
                DROPS::GridFunctionCL<> integrand( 1., qdom.vertex_size()); // Gridfunction with constant 1 everywhere
                //double tmp_neg, tmp_pos;
                area_neg +=quad( integrand, 1., qdom, DROPS::NegTetraC);
                area_pos +=quad( integrand, 1., qdom, DROPS::PosTetraC);
            }
        }
    }
    //analytical solution of negative area is 0.78539815
    const double area_neg_ref= 0.78539815;
    std::cout << "Area of the negative part: " << area_neg << ", area of the positive part: " << area_pos << std::endl;
    if (fabs(area_neg - area_neg_ref) > 1e-9 || fabs(area_neg + area_pos - 4.) > 1e-9)
        std::cout << "This is NOT consistent with the reference solution" << std::endl;
    else
        std::cout << "This is consistent with the reference solution" << std::endl;

}

int main()
{
    try {
        test_triangle_partition();
        test_bnd_integral();
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
    return 0;
}
