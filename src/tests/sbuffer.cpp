/// \file sbuffer.cpp
/// \brief tests a buffer or ring-buffer of fixed size with wrap-around indices
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
#include <iostream>


int TestPushBack()
{
    DROPS::SBufferCL<int, 3> b;
    std::cout << "TestPushBack Row1: " << b[0] << " " << b[1] << " " << b[2] << std::endl;
    b[0]= 0; b[1]= 1; b[2]= 2;
    std::cout << "TestPushBack Row2: " << b[0] << " " << b[1] << " " << b[2] << std::endl;
    b.push_back( 3);
    std::cout << "TestPushBack Row3: " << b[0] << " " << b[1] << " " << b[2] << std::endl;
    std::cout << "TestPushBack Row4: " << b[-1]<< " " << b[-2]<< " " << b[-3]<< std::endl;
    std::cout << "TestPushBack Row5: " << b[3] << " " << b[4] << " " << b[5] << std::endl;
    std::cout << std::endl;
    return 0;
}

int TestRotate()
{
    DROPS::SBufferCL<int, 3> b;
    std::cout << "TestRotate Row1: " << b[0] << " " << b[1] << " " << b[2] << std::endl;
    b[0]= 0; b[1]= 1; b[2]= 2;
    std::cout << "TestRotate Row2: " << b[0] << " " << b[1] << " " << b[2] << std::endl;
    b.rotate( 2);
    std::cout << "TestRotate Row3: " << b[0] << " " << b[1] << " " << b[2] << std::endl;
    b.rotate( -1);
    b.rotate( -1);
    std::cout << "TestRotate Row4: " << b[0] << " " << b[1] << " " << b[2] << std::endl;
    return 0;
}

int main (int, char**)
{
  try {
    return TestPushBack() + TestRotate();
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
