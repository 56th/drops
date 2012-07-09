/// \file poissonCoeff.cpp
/// \brief boundary and source functions for the poisson-type problems
/// \author LNM RWTH Aachen: Christoph Lehrenfeld, Joerg Grande, Thorolf Schulte, Liang Zhang

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

#include "misc/bndmap.h"
#include "misc/params.h"
#include "poisson/poissonCoeff.h"
#define PI 3.14159265

//======================================================================================================================
//                                  static members declarations of poissonCoeff class
//======================================================================================================================

namespace DROPS{
ParamCL PoissonCoeffCL::C_;
double PoissonCoeffCL::dx_;
double PoissonCoeffCL::dy_;
int PoissonCoeffCL::nx_;
int PoissonCoeffCL::ny_;
double PoissonCoeffCL::dt_;
scalar_tetra_function PoissonCoeffCL::q;
double PoissonCoeffCL::alpha;
scalar_tetra_function PoissonCoeffCL::f;
instat_scalar_fun_ptr PoissonCoeffCL::Solution;
instat_scalar_fun_ptr PoissonCoeffCL::InitialCondition;
vector_tetra_function PoissonCoeffCL::Vel;
instat_scalar_fun_ptr PoissonCoeffCL::interface;
}

//======================================================================================================================
//                                  special Functions for scalar.cpp
//======================================================================================================================

extern DROPS::ParamCL P;

double Heat(const DROPS::Point3DCL&, double)
{
    extern DROPS::ParamCL P;
    return P.get<double>("Exp.Heat")/P.get<double>("Exp.Lambda")*1e-3;
}

/// boundary description of a neumann problem
// uses constant function f = (-1)^seg *4.0
template<int sel>
double NeuConst( const DROPS::Point3DCL& , double ) { return std::pow(-1.,sel)*4.0; }

/// boundary description of a neumann problem
// uses exp-function f =  e^(t)* e^(px + py + pz)
template<int sel>
double NeuExp( const DROPS::Point3DCL& p, double t) { return std::pow(-1.,sel)*std::exp(t)*std::exp(p[0]+p[1]+p[2]); }

/// boundary description of a neumann problem
// uses polynomial function
double NeuPoly( const DROPS::Point3DCL& p, double ) { return -64.0*p[0]*p[1]*(1.0-p[0])*(1.0-p[1]);}

/// \brief Nusselt velocity profile for flat film
DROPS::Point3DCL Nusselt(const DROPS::TetraCL& tet, const DROPS::BaryCoordCL& b, double)
{
    static bool first = true;
    static double dx, dy;
    static double Rho, Mu;   //density, viscosity
    //dirty hack
    if (first) {
        std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
        size_t idx_;
        while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
            mesh[idx_]= ' ';
        std::istringstream brick_info( mesh);
        brick_info >> dx >> dy;
        Rho = P.get<double>("Exp.Rho");
        Mu  = P.get<double>("Exp.Mu");
        first = false;
    }

    DROPS::Point3DCL ret;
    const double d= DROPS::GetWorldCoord(tet,b)[1]/dy,
        U= Rho*9.81*dy*dy/2/Mu;  //U=gh^2/(2*nu)
    ret[0]= U*(2-d)*d;
    ret[1]=0.;
    ret[2]=0.;

    return ret;
}

static DROPS::RegisterScalarFunction regscaheat("Heat", Heat);
static DROPS::RegisterScalarFunction regscaconstpos("NeuConstPos", NeuConst<0>);
static DROPS::RegisterScalarFunction regscaconstneg("NeuConstNeg", NeuConst<1>);
static DROPS::RegisterScalarFunction regscaexppos("NeuExpPos", NeuExp<0>);
static DROPS::RegisterScalarFunction regscaexpneg("NeuExpNeg", NeuExp<1>);
static DROPS::RegisterScalarFunction regscapoly("NeuPoly", NeuPoly);
static DROPS::RegisterVectorFunction regvecnus("Nusselt", Nusselt);

//======================================================================================================================
//
//               Registered functions in function container which are used in corresponding *Ex.json file
//
//======================================================================================================================

/****************************************
 *  Example 1:                        \n*
 *   - stationary setup               \n*
 *   - constant diffusion             \n*
 *   - no convection                  \n*
 *   - no reaction                    \n*
 *   - given solution (adapted r.h.s) \n*
 ****************************************/
 namespace statpoissonExample {
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 0.0;
    }
    /// \brief Convection: no convection
    DROPS::Point3DCL Flowfield(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return DROPS::Point3DCL(0.);
    }
    /// \brief Right-hand side
    double Source(const DROPS::TetraCL& tet, const DROPS::BaryCoordCL& b, double) {
        const DROPS::Point3DCL p= DROPS::GetWorldCoord( tet, b);
        return 128*(p[1]*p[2]*(1.-p[1])*(1.-p[2])
                + p[0]*p[2]*(1.-p[0])*(1.-p[2])
                + p[0]*p[1]*(1.-p[0])*(1.-p[1]));
    }
    /// \brief Diffusion
    double Diffusion(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 1.;
    }
    ///Neumann boundary condition at x=1
    double Neumann(const DROPS::Point3DCL& p, double) {
        return -64.*p[0]*p[1]*p[2]*(1-p[1])*(1-p[2]);
    }
    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double) {
        return 1 + 64.*p[0]*p[1]*p[2]*(1-p[0])*(1-p[1])*(1-p[2]);
    }
    /// \brief Initial value
    double InitialValue( const DROPS::Point3DCL& , double) {
        return 1.;
    }

    static DROPS::RegisterScalarFunction regscaq("Example1_Reaction",     Reaction    );
    static DROPS::RegisterScalarFunction regscaf("Example1_Source",       Source      );
    static DROPS::RegisterScalarFunction regscas("Example1_Solution",     Solution    );
    static DROPS::RegisterScalarFunction regscaa("Example1_Diffusion",    Diffusion   );
    static DROPS::RegisterScalarFunction regscan("Example1_Neumann",      Neumann     );
    static DROPS::RegisterVectorFunction regscav("Example1_Flowfield",    Flowfield   );
    static DROPS::RegisterScalarFunction regscai("Example1_InitialValue", InitialValue);

} //end of namespace


/****************************************
 *  Example 2:                        \n*
 *   - instationary setup             \n*
 *   - constant diffusion             \n*
 *   - no convection                  \n*
 *   - no reaction                    \n*
 *   - given solution (adapted r.h.s) \n*
 ****************************************/
 namespace instatpoissonExample {
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 0.0;
    }
    /// \brief Convection: no convection
    DROPS::Point3DCL Flowfield(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return DROPS::Point3DCL(0.);
    }
    /// \brief Right-hand side
    double Source(const DROPS::TetraCL& tet, const DROPS::BaryCoordCL& b, double t) {
        const DROPS::Point3DCL p= DROPS::GetWorldCoord( tet, b);
        return (-2.0*std::exp(t)*std::exp(p[0]+p[1]+p[2]));
    }
    /// \brief Diffusion
    double Diffusion(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 1.;
    }
    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double t) {
        return (std::exp(t)*std::exp(p[0]+p[1]+p[2]));
    }

    static DROPS::RegisterScalarFunction regscaq("Example2_Reaction",     Reaction    );
    static DROPS::RegisterScalarFunction regscaf("Example2_Source",       Source      );
    static DROPS::RegisterScalarFunction regscas("Example2_Solution",     Solution    );
    static DROPS::RegisterScalarFunction regscaa("Example2_Diffusion",    Diffusion   );
    static DROPS::RegisterVectorFunction regscav("Example2_Flowfield",    Flowfield   );

} //end of namespace


/****************************************
 *  Example 3:                        \n*
 *   - stationary setup               \n*
 *   - quasi-1D setup                 \n*
 *   - constant diffusion             \n*
 *   - convection                     \n*
 *   - no reaction                    \n*
 *   - source = 1                     \n*
 ****************************************/
 namespace convdiffExample {

    /// \brief Reaction: no reaction
    double Reaction(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 0.0;
    }
    /// \brief Convection: constant flow in x direction
    DROPS::Point3DCL Flowfield(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        DROPS::Point3DCL v(0.);
        v[0] = 1.0;
        return v;
    }
    /// \brief Right-hand side
    double Source(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 1.0;
    }
    /// \brief Diffusion
    double Diffusion(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 1.0;
    }
    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double) {
        return 1.0 + p[0] + (1-exp(p[0]))/exp(1.0);
    }
    /// \brief Initial value
    double InitialValue( const DROPS::Point3DCL& , double) {
        return 1.0;
    }

    static DROPS::RegisterScalarFunction regscaq("Example3_Reaction",     Reaction    );
    static DROPS::RegisterScalarFunction regscaf("Example3_Source",       Source      );
    static DROPS::RegisterScalarFunction regscas("Example3_Solution",     Solution    );
    static DROPS::RegisterScalarFunction regscaa("Example3_Diffusion",    Diffusion   );
    static DROPS::RegisterVectorFunction regscav("Example3_Flowfield",    Flowfield   );
    static DROPS::RegisterScalarFunction regscai("Example3_InitialValue", InitialValue);

} //end of namespace
/****************************************
 *  Example 4:                        \n*
 *   - stationary setup               \n*
 *   - quasi-1D setup                 \n*
 *   - constant diffusion             \n*
 *   - convection                     \n*
 *   - no reaction                    \n*
 *   - source = 1                     \n*
 ****************************************/
 namespace statSUPGExample {
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 0.0;
    }
    /// \brief Convection: constant flow in x direction
    DROPS::Point3DCL Flowfield(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        DROPS::Point3DCL v(0.);
        v[0] = 1.0;
        return v;
    }
    /// \brief Right-hand side
    double Source(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 1.0;
    }
    /// \brief Diffusion
    double Diffusion(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 1.0;
    }
    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double)
    {
        double D = P.get<double>("PoissonCoeff.Diffusion");
        return p[0] - (1 - exp(p[0]/D))/(1 - exp(1./D));
    }
    static DROPS::RegisterScalarFunction regscaq("SUPG_Reaction",     Reaction    );
    static DROPS::RegisterScalarFunction regscaf("SUPG_Source",       Source      );
    static DROPS::RegisterScalarFunction regscas("SUPG_Solution",     Solution    );
    static DROPS::RegisterScalarFunction regscaa("SUPG_Diffusion",    Diffusion   );
    static DROPS::RegisterVectorFunction regscav("SUPG_Flowfield",    Flowfield   );
} //end of namespace

/****************************************
 *  Example 5:                        \n*
 *   - instationary setup             \n*
 *   - constant diffusion             \n*
 *   - convection                     \n*
 *   - no reaction                    \n*
 ****************************************/
 namespace instatSUPGExample {
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 0.0;
    }
    /// \brief Convection: constant flow in x direction
    DROPS::Point3DCL Flowfield(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
    //(1, 0, 0)^T
        DROPS::Point3DCL v(0.);
        v[0] = 1.0;
        return v;
    }
    /// \brief Right-hand side
    double Source(const DROPS::TetraCL& tet, const DROPS::BaryCoordCL& b, double t) {
    // 1 - (1 + D)e^(-t-y)
        const DROPS::Point3DCL p= DROPS::GetWorldCoord( tet, b);
        static bool first = true;
        //alpha to control diffusion parameter
        static double alpha;
        if (first) {
            alpha= P.get<double>("PoissonCoeff.Diffusion");
            first= false;
        }
        double ret=0.;
        ret = 1. - (1 + alpha) * exp(-t-p[1]);
        return ret;
    }
    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double t)
    {
        static bool first = true;
        //alpha to control diffusion parameter
        static double alpha;
        if (first) {
            alpha = P.get<double>("PoissonCoeff.Diffusion");
            first=false;
        }
        return exp(-t-p[1]) + p[0] - (1 - exp(p[0]/alpha))/(1 - exp(1./alpha));
    }
    static DROPS::RegisterScalarFunction regscaq("instatSUPG_Reaction",     Reaction    );
    static DROPS::RegisterScalarFunction regscaf("instatSUPG_Source",       Source      );
    static DROPS::RegisterScalarFunction regscas("instatSUPG_Solution",     Solution    );
    static DROPS::RegisterVectorFunction regscav("instatSUPG_Flowfield",    Flowfield   );
} //end of namespace

/****************************************
 *  Example 6:                        \n*
 *   - instationary setup             \n*
 *   - constant diffusion             \n*
 *   - convection                     \n*
 *   - no reaction                    \n*
 ****************************************/
 namespace AdjointExample {
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 0.0;
    }
    /// \brief Convection: constant flow in x direction
    DROPS::Point3DCL Flowfield(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        DROPS::Point3DCL v(0.);
        v[0] = 1.;
        return v;
    }
    /// \brief Right-hand side
    double Source(const DROPS::TetraCL& tet, const DROPS::BaryCoordCL& b, double t) {
        const DROPS::Point3DCL p= DROPS::GetWorldCoord( tet, b);
        return -exp(-t)*exp(-p[0]);
    }
    /// \brief Diffusion
    double Diffusion(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 1.0;
    }
    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double t)
    {
        return exp(-t)*exp(-p[0]);
    }
    static DROPS::RegisterScalarFunction regscaq("Adjoint_Reaction",     Reaction    );
    static DROPS::RegisterScalarFunction regscaf("Adjoint_Source",       Source      );
    static DROPS::RegisterScalarFunction regscas("Adjoint_Solution",     Solution    );
    static DROPS::RegisterScalarFunction regscaa("Adjoint_Diffusion",    Diffusion   );
    static DROPS::RegisterVectorFunction regscav("Adjoint_Flowfield",    Flowfield   );
} //end of namespace

/****************************************
 *  Example 7:                        \n*
 *   - instationary setup             \n*
 *   - constant diffusion             \n*
 *   - convection                     \n*
 *   - no reaction                    \n*
 ****************************************/
 namespace ALEExample1 {
    //refH, Mag, paraX and paraT are used to change the free surface functions
    double refH   = 0.2;
    double Mag    = 0.25;
    double paraX  = 6. * PI;
    double paraT  = 10. * PI;
    //Free surface
    double Interface( const DROPS::Point3DCL& p, double t)
    {

        double h= refH + refH * Mag * sin ( paraX * p[0]  + paraT * t );
        return h;
    }
    //Transform the physical to coordinates to reference coordinates
    DROPS::Point3DCL TransBack(const DROPS::Point3DCL& p, double t)
    {
        DROPS::Point3DCL ref(p);
        ref[1]*= refH/Interface(p, t);
        return ref;
    }
    //Gradx, Grady, Grad1 and Grad2 are used for source term
    double Gradx( const DROPS::Point3DCL& ref, double t)   //a=h_x
    {
        return ref[1]* paraX * Mag * cos(paraX * ref[0]  + paraT * t);
    }
    double Grady( const DROPS::Point3DCL& ref, double t)   //b=h_y
    {
        return 1. + Mag * sin(paraX * ref[0]  + paraT * t);
    }
    double Grad1( const DROPS::Point3DCL& ref, double t)   //\nabla_y(a/b)
    {
        double ay=paraX * Mag * cos(paraX * ref[0]  + paraT * t);
        return ay/Grady(ref, t);
    }
    double Grad2( const DROPS::Point3DCL& ref, double t)    //\nabla_x(a/b)
    {
        double ax     = -ref[1]* paraX * paraX * Mag * sin(paraX * ref[0]  + paraT * t);
        double bx     = paraX * Mag * cos(paraX * ref[0]  + paraT * t);
        return (ax*Grady(ref, t) - Gradx(ref, t)*bx)/(Grady(ref,t)*Grady(ref, t));
    }
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 0.0;
    }
    DROPS::Point3DCL Flowfield(const DROPS::TetraCL& tet, const DROPS::BaryCoordCL& b, double t) {
        DROPS::Point3DCL ref = TransBack(DROPS::GetWorldCoord( tet, b), t);
        DROPS::Point3DCL v(0.);
        v[0] = 1.0;
        v[1] = ref[1] * Mag * paraT * cos(paraX * ref[0]  + paraT * t);
        return v;
    }
    DROPS::Point3DCL ALEFlowfield(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        DROPS::Point3DCL v(0.);
        v[0] = 1.;
        v[1] = 0.;
        return v;
    }

    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double t)
    {
        DROPS::Point3DCL ref = TransBack(p, t);
        return exp(t)*exp(ref[0] + ref[1] + ref[2]);
    }
    /// \brief Right-hand side
    double Source(const DROPS::TetraCL& tet, const DROPS::BaryCoordCL& bary, double t) {
        static bool first = true;
        //alpha to control diffusion parameter, you could change a small number to make problem convection-dominated
        static double alpha;
        if (first) {
            alpha= P.get<double>("PoissonCoeff.Diffusion");
            first= false;
        }
        const DROPS::Point3DCL p= DROPS::GetWorldCoord( tet, bary);
        DROPS::Point3DCL ref = TransBack( p, t);
        double a= Gradx(ref,t);
        double b= Grady(ref,t);
        double c= Grad1(ref,t);
        double d= Grad2(ref,t);
        double sol= Solution(p, t);
        double timederiv = sol;
        double conv      = (1. - a/b)*sol;
        double diff1     = -alpha*(1.-d - 2.*a/b + a/b*c + a*a/b/b)*sol;
        double diff2     = -alpha/(b*b)*sol;
        double diff3     = -alpha*sol;
        return timederiv + conv + diff1 +diff2 +diff3;
    }
    static DROPS::RegisterScalarFunction regscaq("TestALE_Reaction",     Reaction    );
    static DROPS::RegisterScalarFunction regscaf("TestALE_Source",       Source      );
    static DROPS::RegisterScalarFunction regscaint("TestALE_Interface",  Interface   );
    static DROPS::RegisterScalarFunction regscas("TestALE_Solution",     Solution    );
    static DROPS::RegisterVectorFunction regscav("TestOrigin_Velocity",  Flowfield   );
    static DROPS::RegisterVectorFunction regscaalev("TestALE_Velocity",  ALEFlowfield);
} //end of namespace

/****************************************
 *  Example 8:                        \n*
 *   - instationary setup             \n*
 *   - constant diffusion             \n*
 *   - convection                     \n*
 *   - no reaction                    \n*
 ****************************************/
 namespace ALEExample2{
    //refH, Mag, paraX and paraT are used to change the free surface functions
    double refH   = 0.2;
    double Mag    = 0.25;
    double paraT  = 10. * PI;
    //Free surface
    double Interface( const DROPS::Point3DCL&, double t)
    {

        double h= refH + refH * Mag * sin (paraT * t );
        return h;
    }
    //Transform the physical to coordinates to reference coordinates
    DROPS::Point3DCL TransBack(const DROPS::Point3DCL &p, double t)
    {
        DROPS::Point3DCL ref(0.);
        ref[0] = p[0];
        ref[1] = refH * p[1]/Interface(p, t);
        ref[2] = p[2];        
        return ref;
    }
    //Grady
    double Grady( const DROPS::Point3DCL&, double t)   //b=h_y
    {
        return 1. + Mag * sin( paraT * t);
    } 
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 0.0;
    }
    DROPS::Point3DCL Flowfield(const DROPS::TetraCL& tet, const DROPS::BaryCoordCL& b, double t) {
        DROPS::Point3DCL ref = TransBack(DROPS::GetWorldCoord( tet, b), t);
        DROPS::Point3DCL v(0.);
        v[0] = 1.0;
        v[1] = ref[1] * Mag * paraT * cos( paraT * t);
        return v; 
    }

    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double t)
    {
        DROPS::Point3DCL ref = TransBack(p, t);
        return 100. * exp(10*t) * (ref[0] - ref[0] * ref [0] ) * (0.2 * ref[1] - ref[1]* ref[1]) * (0.2 * ref[2] - ref[2] * ref[2]);
    }
    ///ALL partial derivatives are based on ALE coordinates
    double Solx( const DROPS::Point3DCL& ref, double t)
    {
        return 100. * exp(10*t) * (1. -2. * ref[0]) * (0.2 * ref[1] - ref[1] * ref[1]) * (0.2 * ref[2] - ref[2]* ref[2]);
    }
    double Solxx( const DROPS::Point3DCL& ref, double t)
    {
        return 100. * exp(10*t) * (-2.) * (0.2 * ref[1] - ref[1] * ref[1]) * (0.2 * ref[2] - ref[2]* ref[2]);
    }
    double Solyy( const DROPS::Point3DCL& ref, double t)
    {
        return 100. * exp(10*t) * (ref[0] - ref[0]* ref [0]) * ( -2.) * (0.2 * ref[2] - ref[2]* ref[2]);
    }
    double Solzz( const DROPS::Point3DCL& ref, double t)
    {
        return 100. * exp(10*t) * (ref[0] - ref[0]* ref [0]) * (0.2 * ref[1] - ref[1]* ref[1]) * ( -2.);
    }
    /// \brief Right-hand side
    double Source(const DROPS::TetraCL& tet, const DROPS::BaryCoordCL& bary, double t) {
        static bool first = true;
        //alpha to control diffusion parameter, you could change a small number to make problem convection-dominated
        static double alpha;
        if (first) {
            alpha= P.get<double>("PoissonCoeff.Diffusion");
            first= false;
        }
        const DROPS::Point3DCL p= DROPS::GetWorldCoord( tet, bary);
        DROPS::Point3DCL ref = TransBack(p, t);
        double b= Grady(ref,t);
        double sol= Solution(p, t);
        double timederiv = 10. * sol;
        double conv      =  Solx(ref, t); 
        double diff1     = -alpha* Solxx(ref, t);
        double diff2     = -alpha/(b*b)*Solyy(ref, t);
        double diff3     = -alpha*Solzz(ref, t);
        return timederiv + conv + diff1 +diff2 +diff3; 
    }
    static DROPS::RegisterScalarFunction regscaq("ALEEx2_Reaction",     Reaction    );
    static DROPS::RegisterScalarFunction regscaf("ALEEx2_Source",       Source      );
    static DROPS::RegisterScalarFunction regscaint("ALEEx2_Interface",  Interface   );
    static DROPS::RegisterScalarFunction regscas("ALEEx2_Solution",     Solution   );
    static DROPS::RegisterVectorFunction regscav("ALEEx2_Velocity",    Flowfield   );
}//end of namespace

/****************************************
 *  Example 9:                        \n*
 *   - instationary setup             \n*
 *   - constant diffusion             \n*
 *   - convection                     \n*
 *   - no reaction                    \n*
 ****************************************/
 namespace ALEExample3{
    //refH, Mag, paraX and paraT are used to change the free surface functions
    double refH   = 0.2;
    double Mag    = 0.25;
    double paraX  = 6. * PI;
    double paraT  = 10. * PI;
    //Free surface
    double Interface( const DROPS::Point3DCL& p, double t)
    {

        double h= refH + refH * Mag * sin (paraX * p[0] - paraT * t );
        return h;
    }
    //Transform the physical to coordinates to reference coordinates
    DROPS::Point3DCL TransBack(const DROPS::Point3DCL &p, double t)
    {
        DROPS::Point3DCL ref(0.);
        ref[0] = p[0];
        ref[1] = refH * p[1]/Interface(p, t);
        ref[2] = p[2];        
        return ref;
    }
    //Gradx, Grady, Grad1 and Grad2 are used for source term
    double Gradx( const DROPS::Point3DCL& ref, double t)   //a=h_x
    {
        return ref[1]* paraX * Mag * cos(paraX * ref[0]  + paraT * t);
    }
    double Grady( const DROPS::Point3DCL& ref, double t)   //b=h_y
    {
        return 1. + Mag * sin(paraX * ref[0]  + paraT * t);
    } 
    double Grad1( const DROPS::Point3DCL& ref, double t)   //\nabla_y(a/b)
    {
        double ay=paraX * Mag * cos(paraX * ref[0]  + paraT * t);
        return ay/Grady(ref, t);
    }
    double Grad2( const DROPS::Point3DCL& ref, double t)    //\nabla_x(a/b)
    {
        double ax     = -ref[1]* paraX * paraX * Mag * sin(paraX * ref[0]  + paraT * t);
        double bx     = paraX * Mag * cos(paraX * ref[0]  + paraT * t);
        return (ax*Grady(ref, t) - Gradx(ref, t)*bx)/(Grady(ref,t)*Grady(ref, t));
    }
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 0.0;
    }
    DROPS::Point3DCL Flowfield(const DROPS::TetraCL& tet, const DROPS::BaryCoordCL& b, double t) {
        DROPS::Point3DCL ref = TransBack(DROPS::GetWorldCoord( tet, b), t);
        DROPS::Point3DCL v(0.);
        v[0] = 1.;
        v[1] = ref[1] * Mag * paraT * cos( paraX * ref[0]  + paraT * t);
        return v; 
    }
    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double t)
    {
        DROPS::Point3DCL ref = TransBack(p, t);
        return 100. * exp(10*t) * (ref[0] - ref[0] * ref [0] ) * (0.2 * ref[1] - ref[1]* ref[1]) * (0.2 * ref[2] - ref[2] * ref[2]);
    }
    ///ALL partial derivatives are based on ALE coordinates
    double Solx( const DROPS::Point3DCL& ref, double t)
    {
        return 100. * exp(10*t) * (1. -2. * ref[0]) * (0.2 * ref[1] - ref[1] * ref[1]) * (0.2 * ref[2] - ref[2]* ref[2]);
    }
    double Solxx( const DROPS::Point3DCL& ref, double t)
    {
        return 100. * exp(10*t) * (-2.) * (0.2 * ref[1] - ref[1] * ref[1]) * (0.2 * ref[2] - ref[2]* ref[2]);
    }
    double Soly( const DROPS::Point3DCL& ref, double t)
    {
        return 100. * exp(10*t) * (ref[0] - ref[0]* ref [0]) * (0.2 - 2.* ref[1]) * (0.2 * ref[2] - ref[2]* ref[2]);
    }
    double Solyy( const DROPS::Point3DCL& ref, double t)
    {
        return 100. * exp(10*t) * (ref[0] - ref[0]* ref [0]) * ( -2.) * (0.2 * ref[2] - ref[2]* ref[2]);
    }
    double Solzz( const DROPS::Point3DCL& ref, double t)
    {
        return 100. * exp(10*t) * (ref[0] - ref[0]* ref [0]) * (0.2 * ref[1] - ref[1]* ref[1]) * ( -2.);
    }
    double Solxy( const DROPS::Point3DCL& ref, double t)
    {
        return 100. * exp(10*t) * (1. -2. * ref[0]) * (0.2  - 2. * ref[1]) * (0.2 * ref[2] - ref[2]* ref[2]);
    }
    /// \brief Right-hand side
    double Source(const DROPS::TetraCL& tet, const DROPS::BaryCoordCL& bary, double t) {
        static bool first = true;
        //alpha to control diffusion parameter, you could change a small number to make problem convection-dominated
        static double alpha;
        if(first){
        alpha = P.get<double>("PoissonCoeff.Diffusion");
        first=false;                                    
        }
        const DROPS::Point3DCL p= DROPS::GetWorldCoord( tet, bary);
        DROPS::Point3DCL ref = TransBack(p, t);
        double a= Gradx(ref,t);
        double b= Grady(ref,t);
        double c= Grad1(ref,t);
        double d= Grad2(ref,t);        
        double sol= Solution(p, t);
        double timederiv = 10. * sol;
        double conv      =  1. * (Solx(ref, t) - a/b * Soly(ref, t)); 
        double diff1     = -alpha* (Solxx(ref, t) - 2 * a/b * Solxy(ref,t ) - d * Soly(ref, t) + a/b * c * Soly(ref, t) + a * a /b/b *Solyy(ref, t));
        double diff2     = -alpha/(b*b)*Solyy(ref, t);
        double diff3     = -alpha*Solzz(ref, t);
        return timederiv + conv + diff1 +diff2 +diff3; 
    }
    static DROPS::RegisterScalarFunction regscaq("ALEEx3_Reaction",     Reaction    );
    static DROPS::RegisterScalarFunction regscaf("ALEEx3_Source",       Source      );
    static DROPS::RegisterScalarFunction regscaint("ALEEx3_Interface",  Interface   );
    static DROPS::RegisterScalarFunction regscas("ALEEx3_Solution",     Solution   );
    static DROPS::RegisterVectorFunction regscav("ALEEx3_Velocity",    Flowfield   );
}//end of namespace

/****************************************
 *  Example 10:                        \n*
 *   - instationary setup             \n*
 *   - constant diffusion             \n*
 *   - convection                     \n*
 *   - no reaction                    \n*
 ****************************************/
 namespace ALEExample4{
    //refH, Mag, paraX and paraT are used to change the free surface functions
    double refH   = 0.2;
    double Mag    = 0.25;
    double paraT  = 10. * PI;
    //Free surface
    double Interface( const DROPS::Point3DCL&, double t)
    {

        double h= refH + refH * Mag * sin (paraT * t );
        return h;
    }
    //Transform the physical to coordinates to reference coordinates
    DROPS::Point3DCL TransBack(const DROPS::Point3DCL &p, double t)
    {
        DROPS::Point3DCL ref(0.);
        ref[0] = p[0];
        ref[1] = refH * p[1]/Interface(p, t);
        ref[2] = p[2];        
        return ref;
    }
    //Grady
    double Grady( const DROPS::Point3DCL&, double t)   //b=h_y
    {
        return 1. + Mag * sin( paraT * t);
    } 
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 0.0;
    }
    DROPS::Point3DCL Flowfield(const DROPS::TetraCL& tet, const DROPS::BaryCoordCL& b, double t) {
        DROPS::Point3DCL ref = TransBack(DROPS::GetWorldCoord( tet, b), t);
        DROPS::Point3DCL v(0.);
        v[0] = exp(5.* t);
        v[1] = ref[1] * Mag * paraT * cos( paraT * t);
        return v; 
    }
    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double t)
    {
        DROPS::Point3DCL ref = TransBack(p, t);
        return 100. * exp(10*t) * (ref[0] - ref[0] * ref [0] ) * (0.2 * ref[1] - ref[1]* ref[1]) * (0.2 * ref[2] - ref[2] * ref[2]);
    }
    ///ALL partial derivatives are based on ALE coordinates
    double Solx( const DROPS::Point3DCL& ref, double t)
    {
        return 100. * exp(10*t) * (1. -2. * ref[0]) * (0.2 * ref[1] - ref[1] * ref[1]) * (0.2 * ref[2] - ref[2]* ref[2]);
    }
    double Solxx( const DROPS::Point3DCL& ref, double t)
    {
        return 100. * exp(10*t) * (-2.) * (0.2 * ref[1] - ref[1] * ref[1]) * (0.2 * ref[2] - ref[2]* ref[2]);
    }
    double Solyy( const DROPS::Point3DCL& ref, double t)
    {
        return 100. * exp(10*t) * (ref[0] - ref[0]* ref [0]) * ( -2.) * (0.2 * ref[2] - ref[2]* ref[2]);
    }
    double Solzz( const DROPS::Point3DCL& ref, double t)
    {
        return 100. * exp(10*t) * (ref[0] - ref[0]* ref [0]) * (0.2 * ref[1] - ref[1]* ref[1]) * ( -2.);
    }
    /// \brief Right-hand side
    double Source(const DROPS::TetraCL& tet, const DROPS::BaryCoordCL& bary, double t) {
        static bool first = true;
        //alpha to control diffusion parameter, you could change a small number to make problem convection-dominated
        static double alpha;
        if(first){
        alpha = P.get<double>("PoissonCoeff.Diffusion");
        first=false;                                    
        }
        const DROPS::Point3DCL p= DROPS::GetWorldCoord( tet, bary);
        DROPS::Point3DCL ref = TransBack(p, t);
        double b= Grady(ref,t);
        double sol= Solution(p, t);
        double timederiv = 10. * sol;
        double conv      =  exp(5.* t) * Solx(ref, t); 
        double diff1     = -alpha* Solxx(ref, t);
        double diff2     = -alpha/(b*b)*Solyy(ref, t);
        double diff3     = -alpha*Solzz(ref, t);
        return timederiv + conv + diff1 +diff2 +diff3; 
    }
    static DROPS::RegisterScalarFunction regscaq("ALEEx4_Reaction",     Reaction    );
    static DROPS::RegisterScalarFunction regscaf("ALEEx4_Source",       Source      );
    static DROPS::RegisterScalarFunction regscaint("ALEEx4_Interface",  Interface   );
    static DROPS::RegisterScalarFunction regscas("ALEEx4_Solution",     Solution   );
    static DROPS::RegisterVectorFunction regscav("ALEEx4_Velocity",    Flowfield   );
}//end of namespace
 namespace ALEBala {
    //refH, Mag, paraX and paraT are used to change the free surface functions
    double refH   = 0.2;
    double Mag    = 0.25;
    double paraX  = 0.2 * PI;
    double paraT  = 10. * PI;
    //Free surface
    double Interface( const DROPS::Point3DCL& p, double t)
    {
        double h= refH + refH * Mag * sin ( paraX * p[0]  + paraT * t );
        return h;
    }
    double Inter( const DROPS::Point3DCL&, double)
    {
        double ret =0.01;
        return ret;
    }
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 0.0;
    }
    double QC( const DROPS::Point3DCL& , double t)
    {
        double ret= pow(refH + refH * Mag * sin ( paraT * t ), 3) +  refH * Mag * paraT *sin ( paraT * t )/paraX ;
        return ret;
    }
    double flux(const DROPS::Point3DCL& p, double t)
    {
        double ret= QC(p, t) - refH * Mag * paraT *sin ( paraX * p[0] + paraT * t )/paraX;
        return ret;
    }

    double Hx(const DROPS::Point3DCL& p, double t)
    {
        double ret= refH * Mag * paraX* cos ( paraX * p[0]  + paraT * t );
        return ret;
    }
    double Qx(const DROPS::Point3DCL& p, double t)
    {
        double ret = - refH * Mag * paraT * cos( paraX * p[0] + paraT * t );
        return ret;
    }
    DROPS::Point3DCL Flowfield(const DROPS::TetraCL& tet, const DROPS::BaryCoordCL& b, double t) {
        const DROPS::Point3DCL p= DROPS::GetWorldCoord( tet, b);
        DROPS::Point3DCL v(0.);
        double h = Interface(p,t);
        double q = flux(p, t);
        double hx =Hx(p,t);
        double qx =Qx(p,t);
        v[0] = -0.5*(3.*h*p[1]*p[1] - pow(p[1], 3))*(qx/pow(h, 3) + (-3.)*q*hx/pow(h, 4));
        v[1] = 3.*q/pow(h, 3)*(h*p[1] - p[1]*p[1]/2.);
        return v;
    }
    /// \brief Right-hand side
    double Source(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 0.;
    }
    static DROPS::RegisterScalarFunction regscaq("ALEBala_Reaction",     Reaction    );
    static DROPS::RegisterScalarFunction regscaf("ALEBala_Source",       Source      );
    static DROPS::RegisterScalarFunction regscaint("ALEBala_Interface",  Interface   );
    static DROPS::RegisterScalarFunction regscainter("ALEBala_Inter",    Inter       );
    static DROPS::RegisterVectorFunction regscav("ALEBala_Velocity",     Flowfield   );
} //end of namespace

 namespace ALE {
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 0.0;
    }
    double Interface( const DROPS::Point3DCL& p, double t)
    {
        double h;
        h = 0.2 + 0.05 * sin ( p[0] * 0.1 * PI - 20. * PI * t );
        return h;
    }
    DROPS::Point3DCL TransBack(const DROPS::Point3DCL &p, double t)
    {
        DROPS::Point3DCL ref(p);
        ref[1]*= 0.2/Interface(p, t);
        return ref;
    }
    /// \brief Convection: constant flow in x direction
    DROPS::Point3DCL Flowfield(const DROPS::TetraCL& tet, const DROPS::BaryCoordCL& b, double t) {
        DROPS::Point3DCL ref = TransBack( DROPS::GetWorldCoord( tet, b), t);
        DROPS::Point3DCL ret;
        const double d= ref[1]/0.2,
            U= 200.;  //U=gh^2/(2*nu)
        ret[0]= U*(2-d)*d;
        ret[1]=0.;
        ret[2]=0.;
        return ret;
    }
    /// \brief Right-hand side
    double Source(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) {
        return 0.0;
    }
    double BInter(const DROPS::Point3DCL&, double) {
        return 0.01;
    }
    /// \brief Solution ====================================careful

    static DROPS::RegisterScalarFunction regscaq("ALE_Reaction",     Reaction    );
    static DROPS::RegisterScalarFunction regscaf("ALE_Source",       Source      );
    static DROPS::RegisterScalarFunction regscaint("ALE_Interface",  Interface   );
    static DROPS::RegisterVectorFunction regscav("ALE_Velocity",     Flowfield   );
    static DROPS::RegisterScalarFunction regscas("ALE_Inter",        BInter      );
} //end of namespace
