/// \file  dynload.cpp
/// \brief testint dyn. loading
/// \author LNM RWTH Aachen: Christoph Lehrenfeld

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

// include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>

// include function container
#include "misc/funcmap.h"
#include "misc/dynamicload.h"

int main (int , char** )
{    try
    {
        std::vector<std::string> dls(2);
        dls[0] = "misc/libmisc-scalarFunctions";
        dls[1] = "misc/libmisc-vectorFunctions";

        DROPS::dynamicLoad("../", dls);

        //create geometry

        /// TEST BEGIN
        DROPS::BaryCoordCL a;
        DROPS::Point3DCL b;
        
        DROPS::InVecMap & invec = DROPS::InVecMap::getInstance();
        DROPS::InScaMap & insca = DROPS::InScaMap::getInstance();
        DROPS::ScaTetMap & scatet = DROPS::ScaTetMap::getInstance();

        std::cout << " invec.PrintAll(): " << std::endl;
        invec.PrintAll();
        std::cout << " insca.PrintAll(): " << std::endl;
        insca.PrintAll();
        std::cout << " scatet.PrintAll(): " << std::endl;
        scatet.PrintAll();

        std::cout << " &invec = " << &invec << std::endl;
        std::cout << " &insca = " << &insca << std::endl;
        std::cout << " &scatet = " << &scatet << std::endl;

        DROPS::instat_scalar_fun_ptr test2(insca["Zero"]);
        std::cout << " function pointer address of \"Zero\" in InScaMap: " 
                  << test2 << std::endl;

        std::cout << " q2(a) is " << test2( b,0) << std::endl;
        std::cout << " q2(a) should be 0" << std::endl;
        std::cout << " q1(a) is " << insca["One"]( b,0) << std::endl;
        std::cout << " q1(a) should be 1 " << std::endl;

    
        DROPS::scalar_tetra_function test(scatet["Zero"]);

        std::cout << " function pointer address of \"Zero\" in ScaTetMap: " 
                  << test << std::endl;

        DROPS::TetraCL tet;
        std::cout << " p2(a) is " << test( tet,a,0) << std::endl;
        std::cout << " p2(a) should be 0" << std::endl;
        std::cout << " p1(a) is " << scatet["One"]( tet,a,0) << std::endl;
        std::cout << " p1(a) should be 1" << std::endl;

        return 0;
    }
    catch (DROPS::DROPSErrCL& err) { err.handle(); }

}

