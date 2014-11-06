/// \file spacetime_decomp.cpp
/// \brief tests implementation of the decomposition of space time prism-4
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
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#include "num/spacetime_geom.h"
#include "num/spacetime_quad.h"
#include "num/discretize.h"
#include "misc/scopetimer.h"

using namespace DROPS;


double test_levelset1(const Point3DCL& p, __UNUSED__ double t)
{
    return p[0]*p[0]+p[1]*p[1]+p[2]*p[2] - 0.0625;
    //return sin(p[0]*100.0 + p[1]*20.0)+cos(p[2]+1.7); //chaotic function... for testing robustness...
}

double test_f1(const Point3DCL& p, __UNUSED__ double t)
{
    // return 1.0;
    return std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
}

double test_one(__UNUSED__ const Point3DCL& p, __UNUSED__ double t)
{
    return 1.0;
}

typedef double    (*instat_scalar_fun_ptr)(const DROPS::Point3DCL&, double);


// Decompose reference prism4 into pentatopes. 
// Decompose reference hypertrig into pentatopes. 
// Cut reference pentatope with the help of a specific lset function
// Decompose pentatope accordingly and measure the volume on the separate
// parts.
// Consistency check checks if measure adds up to the measure fo the pentatope
void TestSpaceTimeDecompositionsOnAPentatope()
{
    Point4DContainer pcontref;
    GeneralizedPrism4CL a(pcontref);
    // std::cout << " reference prism4 " << std::endl;
    std::vector<PentatopeCL> pentas;
    a.decompose_add_to_pentas(pentas);
    // std::cout << "pentas.size() = " << pentas.size() << std::endl;
    // for (Uint i = 0; i < pentas.size(); i++)
    //     std::cout << pentas[i] << std::endl;

    HyperTrigCL b(pcontref);
    // std::cout << " reference hypertrig " << std::endl;
    // std::cout << b << std::endl;
    b.decompose_to_pentas(pentas);
    // std::cout << " ------------------------------ " << std::endl;
    // for (Uint i = 0; i < pentas.size(); i++)
    //     std::cout << pentas[i] << std::endl;

    PentatopeCL p(pcontref);
    std::vector<PentatopeCL> pospentas;
    std::vector<PentatopeCL> negpentas;
    std::vector<Tetra4DCL> iftets;
    // std::cout << " reference pentatope " << std::endl;
    // std::cout << p << std::endl;
    SArrayCL<double,5> lsetvals;
    lsetvals[0] = -3.0;
    lsetvals[1] = 4.0;
    lsetvals[2] = -2.0;
    lsetvals[3] = -1.0;
    lsetvals[4] = 1.0;

    p.decompose_add_signed_pentas_4dtets(lsetvals,negpentas,pospentas,iftets);

    double posmeas = 0, negmeas = 0;
    // std::cout << " neg pentas: " << std::endl;

    for (Uint i = 0; i < negpentas.size(); i++)
    {
        // std::cout << negpentas[i] << std::endl;
        negmeas += negpentas[i].Measure();
    }
    
    // std::cout << " pos pentas: " << std::endl;

    for (Uint i = 0; i < pospentas.size(); i++){
        // std::cout << pospentas[i] << std::endl;
        posmeas += pospentas[i].Measure();
    }

    std::cout << " total neg measure = " << negmeas << std::endl;
    std::cout << " total pos measure = " << posmeas << std::endl;
    std::cout << " total     measure = " << posmeas+negmeas << std::endl;

    if (std::fabs(posmeas+negmeas-1.0/24.0) > 1e-12)
        std::cout << " |||||||||||||||||||| WARNING : measures do not add up correctly |||||||||||||||||||| " << std::endl;

}

void TestSpaceTimeDecompositionsOnAPrism4()
{
    std::cout << "\n/// --- Begin of TestSpaceTimeDecompositionsOnAPrism4 --- \\\\\\" << std::endl;

    //make the reference tetraeder 
    Point3DCL p0(0.), p1(0.), p2(0.), p3(0.);
    p1[0]= p2[1]= p3[2]= 1;
    TetraBuilderCL tet (0,p0, p1, p2, p3);
    MultiGridCL MG( tet);
    TetraCL & reftet = *MG.GetTetrasBegin();
    //make reference time interval
    TimeInterval refinterv(0.0,1.0);


    {
        ScopeTimer timingmain("The main test case");

        // std::cout <<    "|||-------------------------------------------------------|||\n";
        std::cout <<    " lset and f not depending on time\n";
        // std::cout <<    "----------------------------------------------------;
        std::cout <<    "principal lattice intervals per edge = 16\n";
        //create decomposition and composite quadrature rule according to test_levelset
        CompositeSTQuadCL<QuadRule> cstquad(reftet,refinterv, test_levelset1, /* ints_per_space_edge */ 16, 1);
        cstquad.Report(std::cout,"","");
        // evaluate test_f at quadrature points of negative part
        GridFunctionCL<double> f = cstquad.EvalOnPart( test_f1, /* posPart */ false);
        GridFunctionCL<double> one = cstquad.EvalOnInterface( test_one);
        
        __UNUSED__ const GridFunctionCL<Point4DCL> & normals = cstquad.GetNormalsOnInterface( );
        __UNUSED__ const GridFunctionCL<double> & nu = cstquad.GetNuOnInterface( );
        
        // do quadrature of test_f on pos part...
        double res = cstquad.QuadOnPart( f, /* posPart */ false);
        // do quadrature of test_one on interface...
        double surfmeas = cstquad.QuadOnInterface( one );
        // compare to exact solution:
        const double ref = M_PI/2048.0; // ref. solution (volume)
        const double ref2 = M_PI/32.0; // ref. solution (surface)
        std::cout << "int_h = " << res << "(=!=" <<  ref  << " = int )" << std::endl;
        std::cout << "rel_error = " << std::fabs(res-ref)/ref << std::endl;
        std::cout << "correct implementation yields: total neg measure = 0.00763306 " << std::endl;
        std::cout << "correct implementation yields: total pos measure = 0.159034 " << std::endl;
        std::cout << "correct implementation yields: total measure = 0.166667 " << std::endl;
        std::cout << "correct implementation yields: int_h = 0.00140492 " << std::endl;
        std::cout << "surface measure int_if_h = " << surfmeas << "(=!=" << ref2 << " = int_if )" << std::endl;
        std::cout << "rel_error = " << std::fabs(surfmeas-ref2)/ref2 << std::endl;
        std::cout << "correct implementation yields: int_if_h = 0.0955709" << std::endl;
    }

    std::cout << "End of TestSpaceTimeDecompositionsOnAPrism4\n" << std::endl;
}

int main ()
{
  try 
    {
      // std::cout << "------------------------------------------------------------------------------------------" << std::endl;
      TestSpaceTimeDecompositionsOnAPentatope();
      std::cout << "Tested space time decomposition + quadrature on Pentatope" << std::endl;
      std::cout << "------------------------------------------------------------------------------------------" << std::endl;
      TestSpaceTimeDecompositionsOnAPrism4();
      std::cout << "Tested space time decomposition + quadrature on RefPrism4" << std::endl;
      // std::cout << "------------------------------------------------------------------------------------------" << std::endl;
      return 0;
    }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
