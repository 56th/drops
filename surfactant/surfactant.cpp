/// \file
/// \brief Solve a non-stationary convection-diffusion-equation on a moving interface
/// \author LNM RWTH Aachen: Joerg Grande

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

#include "surfactant/ifacetransp.h"
#include "misc/params.h"
#include "geom/builder.h"
#include "levelset/levelset.h"
#include "levelset/adaptriang.h"
#include "levelset/surfacetension.h"
#include "out/output.h"
#include "out/vtkOut.h"
#include "misc/bndmap.h"

#include <fstream>
#include <string>

using namespace DROPS;

DROPS::ParamCL P;

std::auto_ptr<VTKOutCL> vtkwriter;

DROPS::InVecMap& invecmap= DROPS::InVecMap::getInstance();
DROPS::InScaMap& inscamap= DROPS::InScaMap::getInstance();

instat_vector_fun_ptr the_wind_fun;
instat_scalar_fun_ptr the_lset_fun;
instat_scalar_fun_ptr the_rhs_fun;
instat_scalar_fun_ptr the_sol_fun;

typedef DROPS::Point3DCL (*bnd_val_fun) (const DROPS::Point3DCL&, double);

DROPS::BndCondT bc_wind[6]= { DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC };
bnd_val_fun bf_wind[6];

instat_scalar_fun_ptr sigma( 0);
SurfaceTensionCL sf( sigma, 0);
DROPS::LsetBndDataCL lsbnd( 6);


DROPS::MGBuilderCL* make_MGBuilder (const DROPS::ParamCL& P)
{
    const std::string type= P.get<std::string>( "Type");
    if (type != std::string("BrickBuilder"))
        throw DROPS::DROPSErrCL(std::string( "make_MGBuilder: Unknown Domain: ") + type + std::string("\n"));

    const Point3DCL orig( P.get<Point3DCL>( "Origin")),
                      e1( P.get<Point3DCL>( "E1")),
                      e2( P.get<Point3DCL>( "E2")),
                      e3( P.get<Point3DCL>( "E3"));
    const Uint n1= P.get<Uint>( "N1"),
               n2= P.get<Uint>( "N2"),
               n3= P.get<Uint>( "N3");
    return new BrickBuilderCL( orig, e1, e2, e3, n1, n2, n3);
}

DROPS::Point3DCL WindVelocity;
DROPS::Point3DCL constant_wind (const DROPS::Point3DCL&, double)
{
    return WindVelocity;
}
static RegisterVectorFunction regvec_constant_wind( "ConstantWind", constant_wind);

DROPS::Point3DCL RadDrop;
DROPS::Point3DCL PosDrop;
double ellipsoid (const DROPS::Point3DCL& p, double)
{
    const DROPS::Point3DCL x= (p - PosDrop)/RadDrop;
    return x.norm_sq() - 1.;
}
static RegisterScalarFunction regsca_ellipsoid_lset( "Ellipsoid", ellipsoid);

// ==non-stationary test case "HeatConduction"
// "ConstantWind" 0
// "Ellipsoid" pos:[0,0,0], rad:[1,1,1]

double heat_conduction_u0 (const Point3DCL& p)
{
    return p[0]*p[1];
}

double heat_conduction_rhs (const Point3DCL&, double)
{
    return 0.;
}
static RegisterScalarFunction regsca_heat_conduction_rhs( "HeatConductionRhs", heat_conduction_rhs);

double heat_conduction_sol (const Point3DCL& p, double t)
{
    return std::exp( -6.*t)*heat_conduction_u0( p);
}
static RegisterScalarFunction regsca_heat_conduction_sol( "HeatConductionSol", heat_conduction_sol);


// ==non-stationary test case "AxisScaling"==

double axis_scaling (double t)
{
    return 1. + 0.25*std::sin( t);
}

DROPS::Point3DCL axis_scaling_wind (const DROPS::Point3DCL&, double t)
{
    return MakePoint3D( 0.25*std::cos( t)/axis_scaling( t), 0., 0.);
}
static RegisterVectorFunction regvec_axis_scaling_wind( "AxisScalingWind", axis_scaling_wind);

double axis_scaling_lset (const Point3DCL& p, double t)
{
    return std::pow( p[0]/axis_scaling( t), 2) + std::pow( p[1], 2) + std::pow( p[2], 2) - 1.;
}
static RegisterScalarFunction regsca_axis_scaling_lset( "AxisScalingLset", axis_scaling_lset);

double axis_scaling_sol (const Point3DCL& p, double t)
{
    return std::exp( -6.*t)*heat_conduction_u0( p); // The same as for heat conduction test case.
}
static RegisterScalarFunction regsca_axis_scaling_sol( "AxisScalingSol", axis_scaling_sol);

double axis_scaling_rhs (const Point3DCL& p, double t)
{
    const double a= axis_scaling( t);
    const double bf4= (p/MakePoint3D( std::pow( a, 2), 1., 1.)).norm_sq();
    const double bf6= (p/MakePoint3D( std::pow( a, 3), 1., 1.)).norm_sq();

    const double mat_der=  0.25*std::cos( t)/a - 6.;
    const double reaction= 0.25*std::cos( t)/a/bf4*( std::pow( p[1], 2) + std::pow( p[2], 2)
        - 2.*std::pow( p[0]/a, 2)*(1. + 1./std::pow( a, 2) - bf6/bf4));
    const double diffusion= (1. + 1./std::pow( a, 2))/bf4*(1. + 2./std::pow( a, 2) + 2.*bf6/bf4);

    return (mat_der + reaction - diffusion)*axis_scaling_sol( p, t);
}
static RegisterScalarFunction regsca_axis_scaling_rhs( "AxisScalingRhs", axis_scaling_rhs);


double sphere_dist (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL x( p - PosDrop);
    return x.norm() - RadDrop[0];
}
static RegisterScalarFunction regsca_sphere_dist_lset( "SphereDist", sphere_dist);


typedef double (*dist_funT) (const DROPS::Point3DCL&, double);

double sphere_2move (const DROPS::Point3DCL& p, double t)
{
    DROPS::Point3DCL x( p - (PosDrop + t*constant_wind(p, t)));
    return x.norm() - RadDrop[0];
}

// ==stationary test case "LaplaceBeltrami0"==
// Sphere around 0, RadDrop 1, wind == 0
// A right hand side from C.J. Heine...
const double a( -13./8.*std::sqrt( 35./M_PI));
double laplace_beltrami_0_rhs (const DROPS::Point3DCL& p, double)
{
    return a*(3.*p[0]*p[0]*p[1] - p[1]*p[1]*p[1]);
}
static RegisterScalarFunction regsca_laplace_beltrami_0_rhs( "LaplaceBeltrami0Rhs", laplace_beltrami_0_rhs);

// ...and the corresponding solution (extended)
double laplace_beltrami_0_sol (const DROPS::Point3DCL& p, double)
{
    return p.norm_sq()/(12. + p.norm_sq())*laplace_beltrami_0_rhs( p, 0.);
//    return 1. + p.norm_sq()/(12. + p.norm_sq())*laplace_beltrami_0_rhs( p, 0.);
}
static RegisterScalarFunction regsca_laplace_beltrami_0_sol( "LaplaceBeltrami0Sol", laplace_beltrami_0_sol);


double sol0t (const DROPS::Point3DCL& p, double t)
{
    const Point3DCL q( p - (PosDrop + t*constant_wind(p, t)));
    const double val( a*(3.*q[0]*q[0]*q[1] - q[1]*q[1]*q[1]));

//    return q.norm_sq()/(12. + q.norm_sq())*val;
    return 1. + q.norm_sq()/(12. + q.norm_sq())*val;
}


template<class DiscP1FunType>
double L2_error (const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
    const DiscP1FunType& discsol, DROPS::instat_scalar_fun_ptr extsol)
{
    double d( 0.);
    const DROPS::Uint lvl = ls.GetLevel();
    const double t= discsol.GetTime();
    const DROPS::MultiGridCL& mg= discsol.GetMG();
    DROPS::InterfaceTriangleCL triangle;
    DROPS::Quad5_2DCL<> qsol, qdiscsol;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        triangle.Init( *it, ls, lsbnd);
        if (triangle.Intersects()) { // We are at the phase boundary.
            for (int ch= 0; ch < 8; ++ch) {
                triangle.ComputeForChild( ch);
                for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri) {
                    qsol.assign( *it, &triangle.GetBary( tri), extsol, t);
                    qdiscsol.assign(  *it, &triangle.GetBary( tri), discsol);
                    d+= DROPS::Quad5_2DCL<>( std::pow( qdiscsol - qsol, 2)).quad( triangle.GetAbsDet( tri));
                }
            }
        }
    }
    return std::sqrt( d);
}

double L2_norm (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
    DROPS::instat_scalar_fun_ptr extsol)
{
    double d( 0.);
    const DROPS::Uint lvl= ls.GetLevel();
    const double        t= ls.t;
    DROPS::InterfaceTriangleCL triangle;
    DROPS::Quad5_2DCL<> qsol;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        triangle.Init( *it, ls, lsbnd);
        if (triangle.Intersects()) { // We are at the phase boundary.
            for (int ch= 0; ch < 8; ++ch) {
                triangle.ComputeForChild( ch);
                for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri) {
                    qsol.assign( *it, &triangle.GetBary( tri), extsol, t);
                    d+= DROPS::Quad5_2DCL<>( qsol*qsol).quad( triangle.GetAbsDet( tri));
                }
            }
        }
    }
    return std::sqrt( d);
}

void LinearLSInit (const DROPS::MultiGridCL& mg, DROPS::VecDescCL& ls, DROPS::instat_scalar_fun_ptr d, double t= 0.)
{
    const DROPS::Uint lvl= ls.GetLevel(),
                      idx= ls.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= d( it->GetCoord(), t);

    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= ls.Data[it->Unknowns( idx)]=
            0.5*(ls.Data[it->GetVertex( 0)->Unknowns( idx)] + ls.Data[it->GetVertex( 1)->Unknowns( idx)]);
    ls.t= t;
}

void LSInit (const DROPS::MultiGridCL& mg, DROPS::VecDescCL& ls, dist_funT d, double t= 0.)
{
    const DROPS::Uint lvl= ls.GetLevel(),
                      idx= ls.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= d( it->GetCoord(), t);

    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= d( 0.5*(it->GetVertex( 0)->GetCoord() + it->GetVertex( 1)->GetCoord()), t);
    ls.t= t;
}

void InitVel ( const MultiGridCL& mg, VecDescCL* vec, BndDataCL<Point3DCL>& Bnd, instat_vector_fun_ptr LsgVel, double t= 0.)
{
    VectorCL& lsgvel= vec->Data;
    const Uint lvl  = vec->GetLevel(),
               vidx = vec->RowIdx->GetIdx();

   DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, sit) {
        if (!Bnd.IsOnDirBnd( *sit))
            DoFHelperCL<Point3DCL, VectorCL>::set( lsgvel, sit->Unknowns( vidx),
                LsgVel(sit->GetCoord(), t));
    }
    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, sit) {
        if (!Bnd.IsOnDirBnd( *sit))
            DoFHelperCL<Point3DCL, VectorCL>::set( lsgvel, sit->Unknowns( vidx),
                LsgVel( (sit->GetVertex(0)->GetCoord() + sit->GetVertex(1)->GetCoord())/2., t));
    }
    vec->t= t;
}

void Strategy (DROPS::MultiGridCL& mg, DROPS::AdapTriangCL& adap, DROPS::LevelsetP2CL& lset)
{
    using namespace DROPS;

    adap.MakeInitialTriang( the_lset_fun);

    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
    // LinearLSInit( mg, lset.Phi, the_lset_fun, 0.);
    LSInit( mg, lset.Phi, the_lset_fun, 0.);

    //DROPS::LevelsetP2CL lset2( mg, lsbnd, sf, P.get<double>("Levelset.Theta"), P.get<double>("Levelset.SD")); // Only for output
    //lset2.idx.CreateNumbering( mg.GetLastLevel(), mg);
    //lset2.Phi.SetIdx( &lset2.idx);
    //LSInit( mg, lset2.Phi, the_lset_fun, 0.);

    const double Vol= lset.GetVolume();
    std::cout << "droplet volume: " << Vol << std::endl;

    BndDataCL<Point3DCL> Bnd_v( 6, bc_wind, bf_wind);
    IdxDescCL vidx( vecP2_FE);
    vidx.CreateNumbering( mg.GetLastLevel(), mg, Bnd_v);
    VecDescCL v( &vidx);
    InitVel( mg, &v, Bnd_v, the_wind_fun, 0.);

    //lset2.SetupSystem( make_P2Eval( mg, Bnd_v, v), P.get<double>("Time.StepSize"));

    SurfactantcGP1CL timedisc( mg,
        P.get<double>("SurfTransp.Theta"), P.get<double>("SurfTransp.Visc"),
        &v, Bnd_v, lset.Phi, lset.GetBndData(),
        P.get<int>("SurfTransp.Iter"), P.get<double>("SurfTransp.Tol"),
        P.get<double>("SurfTransp.OmitBound"));

    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);
    InterfaceP1RepairCL ic_repair( mg, lset.Phi, lset.GetBndData(), timedisc.ic);
    adap.push_back( &ic_repair);
    //LevelsetRepairCL lset2repair( lset2);
    //adap.push_back( &lset2repair);

    // Init Interface-Sol
    timedisc.idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
    std::cout << "NumUnknowns: " << timedisc.idx.NumUnknowns() << std::endl;
    timedisc.ic.SetIdx( &timedisc.idx);
    timedisc.SetInitialValue( the_sol_fun, 0.);

    ScalarFunAsP2EvalCL the_sol_eval( the_sol_fun, 0., &mg);
    if (vtkwriter.get() != 0) {
        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "Levelset") );
        vtkwriter->Register( make_VTKIfaceScalar( mg, timedisc.ic,  "InterfaceSol"));
        vtkwriter->Register( make_VTKVector(      make_P2Eval( mg, Bnd_v, v),      "Velocity"));
        //vtkwriter->Register( make_VTKScalar(      lset2.GetSolution(),             "Levelset2"));
        vtkwriter->Register( make_VTKScalar( the_sol_eval, "TrueSol"));
        //vtkwriter->Register( make_VTKScalar( ScalarFunAsP2EvalCL( the_sol_fun, 0., &mg), "TrueSol"));
        vtkwriter->Write( 0.);
    }
    //if (P.get<int>( "SolutionOutput.Freq") > 0)
    //    DROPS::WriteFEToFile( timedisc.ic, mg, P.get<std::string>( "SolutionOutput.Path"), P.get<bool>( "SolutionOutput.Binary"));

    double L_2x_err= L2_error( lset.Phi, lset.GetBndData(), timedisc.GetSolution(), the_sol_fun);
    std::cout << "L_2x-error: " << L_2x_err
              << "\nnorm of true solution: " << L2_norm( mg, lset.Phi, lset.GetBndData(), the_sol_fun)
              << std::endl;
    double L_inftL_2x_err= L_2x_err;
    std::cout << "L_inftL_2x-error: " <<  L_inftL_2x_err << std::endl;
    BndDataCL<> ifbnd( 0);
    std::cerr << "initial surfactant on \\Gamma: " << Integral_Gamma( mg, lset.Phi, lset.GetBndData(), make_P1Eval(  mg, ifbnd, timedisc.ic)) << '\n';

    const double dt= P.get<double>("Time.StepSize");
    for (int step= 1; step <= P.get<int>("Time.NumSteps"); ++step) {
        std::cout << "======================================================== step " << step << ":\n";
        const double cur_time= step*dt;
        // Assumes (as the rest of Drops), that the current triangulation is acceptable to perform the time-step.
        // If dt is large and AdapRef.Width is small, this may not be true.
        // Watch for large differences in numbers of old and new dof.
        timedisc.InitOld();
        LSInit( mg, lset.Phi, the_lset_fun, cur_time);
        InitVel( mg, &v, Bnd_v, the_wind_fun, cur_time);
        timedisc.DoStep( cur_time);
        std::cout << "surfactant on \\Gamma: " << Integral_Gamma( mg, lset.Phi, lset.GetBndData(), make_P1Eval(  mg, ifbnd, timedisc.ic)) << '\n';
        L_2x_err= L2_error( lset.Phi, lset.GetBndData(), timedisc.GetSolution(), the_sol_fun);
        std::cout << "L_2x-error: " << L_2x_err
                  << "\nnorm of true solution: " << L2_norm( mg, lset.Phi, lset.GetBndData(), the_sol_fun)
                  << std::endl;
        L_inftL_2x_err= std::max( L_inftL_2x_err, L_2x_err);
        std::cout << "L_inftL_2x-error: " << L_inftL_2x_err << std::endl;
        if (vtkwriter.get() != 0 && step % P.get<int>( "VTK.VTKOut") == 0) {
            the_sol_eval.SetTime( cur_time);
            vtkwriter->Write( cur_time);
        }
        if (P.get<int>( "SolutionOutput.Freq") > 0 && step % P.get<int>( "SolutionOutput.Freq") == 0)
            DROPS::WriteFEToFile( timedisc.ic, mg, P.get<std::string>( "SolutionOutput.Path"), P.get<bool>( "SolutionOutput.Binary"));

//        lset2.DoStep();
//        VectorCL rhs( lset2.Phi.Data.size());
//        lset2.ComputeRhs( rhs);
//        lset2.SetupSystem( make_P2Eval( mg, Bnd_v, v, cur_time));
//        lset2.SetTimeStep( P.get<double>("Time.StepSize"));
//        lset2.DoStep( rhs);

//         std::cout << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
//         if (P.get("Levelset.VolCorr", 0)) {
//             double dphi= lset.AdjustVolume( Vol, 1e-9);
//             std::cout << "volume correction is " << dphi << std::endl;
//             lset.Phi.Data+= dphi;
//             std::cout << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
//         }
        //if (C.rpm_Freq && step%C.rpm_Freq==0) { // reparam levelset function
            // lset.ReparamFastMarching( C.rpm_Method);
        const bool doGridMod= P.get<int>("AdaptRef.Freq") && step%P.get<int>("AdaptRef.Freq") == 0;
        const bool gridChanged= doGridMod ? adap.UpdateTriang( lset) : false;
        if (gridChanged) {
            std::cout << "Triangulation changed.\n";
            vidx.DeleteNumbering( mg);
            vidx.CreateNumbering( mg.GetLastLevel(), mg, Bnd_v);
            v.SetIdx( &vidx);
            InitVel( mg, &v, Bnd_v, the_wind_fun, cur_time);
            LSInit( mg, lset.Phi, the_lset_fun, cur_time);
            // timedisc.Update(); // Called unconditionally in DoStep.

            //lset2.SetupSystem( make_P2Eval( mg, Bnd_v, v), P.get<double>("Time.StepSize"));

            std::cout << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            if (P.get<int>( "Levelset.VolCorr")) {
                double dphi= lset.AdjustVolume( Vol, 1e-9);
                std::cout << "volume correction is " << dphi << std::endl;
                lset.Phi.Data+= dphi;
                std::cout << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            }
        }
    }
    std::cout << std::endl;
}


int main (int argc, char* argv[])
{
  try {
    std::ifstream param;
    if (argc!=2) {
        std::cout << "Using default parameter file: surfactant.json\n";
        param.open( "surfactant.json");
    }
    else
        param.open( argv[1]);
    if (!param) {
        std::cout << "error while opening parameter file\n";
        return 1;
    }
    param >> P;
    param.close();
    std::cout << P << std::endl;
 
    WindVelocity= P.get<DROPS::Point3DCL>("Exp.Velocity");
    RadDrop=      P.get<DROPS::Point3DCL>("Exp.RadDrop");
    PosDrop=      P.get<DROPS::Point3DCL>("Exp.PosDrop");
    the_wind_fun= invecmap[P.get<std::string>("Exp.Wind")];
    the_lset_fun= inscamap[P.get<std::string>("Exp.Levelset")];
    the_rhs_fun=  inscamap[P.get<std::string>("Exp.Rhs")];
    the_sol_fun=  inscamap[P.get<std::string>("Exp.Solution")];
    for (Uint i= 0; i < 6; ++i)
        bf_wind[i]= the_wind_fun;

    std::cout << "Setting up domain:\n";
    std::auto_ptr<MGBuilderCL> builder( make_MGBuilder( P.get_child( "Domain")));
    DROPS::MultiGridCL mg( *builder);
    DROPS::AdapTriangCL adap( mg, P.get<double>("AdaptRef.Width"),
                              P.get<int>("AdaptRef.CoarsestLevel"), P.get<int>("AdaptRef.FinestLevel"));

    DROPS::LevelsetP2CL lset( mg, lsbnd, sf);
 
   if (P.get<int>("VTK.VTKOut",0))
        vtkwriter= std::auto_ptr<VTKOutCL>( new VTKOutCL(
            adap.GetMG(),
            "DROPS data",
            P.get<int>("Time.NumSteps")/P.get<int>("VTK.VTKOut") + 1,
            P.get<std::string>("VTK.VTKDir"),
            P.get<std::string>("VTK.VTKName"),
            P.get<std::string>("VTK.TimeFileName"),
            P.get<int>("VTK.Binary"), 
            P.get<bool>("VTK.UseOnlyP1"),
            -1,  /* <- level */
            P.get<bool>("VTK.ReUseTimeFile")));

    if (P.get<bool>("Exp.StationaryPDE") == false) { // Time dependent tests in Strategy
        Strategy( mg, adap, lset);
        return 0;
    }
    adap.MakeInitialTriang( sphere_dist);

    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
    // LinearLSInit( mg, lset.Phi, &sphere_dist);
    LSInit( mg, lset.Phi, sphere_dist, 0.);

    DROPS::IdxDescCL ifaceidx( P1IF_FE);
    ifaceidx.GetXidx().SetBound( P.get<double>("SurfTransp.OmitBound"));
    ifaceidx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
    std::cout << "NumUnknowns: " << ifaceidx.NumUnknowns() << std::endl;

    DROPS::MatDescCL M( &ifaceidx, &ifaceidx);
    DROPS::SetupInterfaceMassP1( mg, &M, lset.Phi, lset.GetBndData());
    std::cout << "M is set up.\n";
    DROPS::MatDescCL A( &ifaceidx, &ifaceidx);
    DROPS::SetupLBP1( mg, &A, lset.Phi, lset.GetBndData(), P.get<double>("SurfTransp.Visc"));
    DROPS::MatrixCL L;
    L.LinComb( 1.0, A.Data, 1.0, M.Data);
    DROPS::VecDescCL b( &ifaceidx);
    DROPS::SetupInterfaceRhsP1( mg, &b, lset.Phi, lset.GetBndData(), laplace_beltrami_0_rhs);

    //DROPS::WriteToFile( M.Data, "m_iface.txt", "M");
    //DROPS::WriteToFile( A.Data, "a_iface.txt", "A");
    //DROPS::WriteFEToFile( b, mg, "rhs_iface.txt", /*binary=*/ true);

    typedef DROPS::SSORPcCL SurfPcT;
    SurfPcT surfpc;
    typedef DROPS::PCGSolverCL<SurfPcT> SurfSolverT;
    SurfSolverT surfsolver( surfpc, P.get<int>("SurfTransp.Iter"), P.get<double>("SurfTransp.Tol"), true);

    DROPS::VecDescCL x( &ifaceidx);
    surfsolver.Solve( L, x.Data, b.Data);
    std::cout << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';

    if (P.get<int>( "SolutionOutput.Freq") > 0)
        DROPS::WriteFEToFile( x, mg, P.get<std::string>( "SolutionOutput.Path"), P.get<bool>( "SolutionOutput.Binary"));

    DROPS::IdxDescCL ifacefullidx( DROPS::P1_FE);
    ifacefullidx.CreateNumbering( mg.GetLastLevel(), mg);
    DROPS::VecDescCL xext( &ifacefullidx);
    DROPS::Extend( mg, x, xext);
    DROPS::NoBndDataCL<> nobnd;
    if (vtkwriter.get() != 0) {
        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "Levelset") );
        vtkwriter->Register( make_VTKIfaceScalar( mg, x, "InterfaceSol"));
        vtkwriter->Write( 0.);
    }

    double L2_err( L2_error( lset.Phi, lset.GetBndData(), make_P1Eval( mg, nobnd, xext), &laplace_beltrami_0_sol));
    std::cout << "L_2-error: " << L2_err << std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
