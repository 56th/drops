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


void ref_triangle_partition()
{
	std::cout<<"=========================RefTrianglePartitionCL test: \n"
	         <<"all 26 level-set sign patterns (except 0, 0, 0) will be prescribed, then the corresponding reference triangle partition info will be given"<<std::endl;
	DROPS::byte VertexNum =1;
	DROPS::byte ls[4];
	ls[VertexNum] = 0;
    int c= 0;
    for (int i= -1; i <= 1; ++i)
      for (int j= -1; j <= 1; ++j)
        for (int k= -1; k <= 1; ++k) {
              if (i == 0 && j == 0 && k == 0) continue;
               ls[ (VertexNum +1)%4 ]= i; ls[(VertexNum +2)%4]= j; ls[(VertexNum +3)%4]= k;
			  DROPS::RefTrianglePartitionCL tri (ls, VertexNum);
			  std::cout <<"====================================================="<<std::endl;
              std::cout << "case: " << c << " ls: " << int(ls[0]) << ' ' << int(ls[1]) << ' ' << int(ls[2]) << ' ' << int(ls[3]) << std::endl;
              std::cout << "The triangle is cut to " << tri.size() << " Sub-triangles;"<<std::endl;
			  for (DROPS::RefTrianglePartitionCL::const_triangle_iterator it= tri.triangle_begin(), end= tri.triangle_end(); it != end; ++it)
			  std::cout << "Sign of the triangle: " <<tri.sign(it)<< " Indices of vertices: "<<int((*it)[0])<<" "<< int((*it)[1])<< " "<<int((*it)[2]) <<std::endl;
			  c++;
          }
}

int main()
{
    try {
        ref_triangle_partition();
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
    return 0;
}
