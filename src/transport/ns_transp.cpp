/// \file ns_transp.cpp
/// \brief flow in measurement cell or brick,
///        basically a copy of twophasedrops.cpp with extensions for mass transport with NitscheXFEM
/// \author LNM RWTH Aachen: Hieu Nguyen, Patrick Esser, Joerg Grande, Sven Gross, Martin Horsky, Christoph Lehrenfeld; SC RWTH Aachen: Oliver Fortmeier

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

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "navstokes/instatnavstokes2phase.h"
//#include "stokes/integrTime.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "out/vtkOut.h"
#include "levelset/coupling.h"
#include "misc/params.h"
#include "levelset/marking_strategy.h"
#include "levelset/adaptriang.h"
#include "levelset/mzelle_hdr.h"
#include "levelset/surfacetension.h"
#include "transport/transportNitsche.h"
#include "surfactant/ifacetransp.h"
#include "levelset/twophaseutils.h"
#include "misc/funcmap.h"
#include "misc/progressaccu.h"
#include "misc/dynamicload.h"

#include "num/stokessolverfactory.h"
#ifdef _PAR
#include "parallel/loadbal.h"
#include "parallel/parmultigrid.h"
#endif
#include <fstream>
#include <sstream>

/// \todo serialization of mass transp. data
/// \todo boundary values for concentration should get more elegant (like poisson perhaps?)
/// \todo solutiononpart-output for ensight AND vtk
/// \todo surfacetension, varSurfaceTension ... flags, output and cases!
DROPS::ParamCL P;
// rho*du/dt - mu*laplace u + Dp = f + rho*g - okn
//                        -div u = 0
//                             u = u0, t=t0


typedef DROPS::BndDataCL<> cBndDataCL;
/// \name Initial data for transport equation
//@{
typedef DROPS::BndDataCL<> cBndDataCL;
typedef cBndDataCL::bnd_val_fun  c_bnd_val_fun;


void TransformToDropOrig( DROPS::Point3DCL& p)
{
    static DROPS::Point3DCL dropPos= P.get<DROPS::Point3DCL>("Levelset.PosDrop");
    p-= dropPos;
}

/// \name Initial data and rhs for coupled mass/surfactant transport
//@{
// TestCase: Sphere around 0, RadDrop 1, vel == 0
double surf_sol (const DROPS::Point3DCL& q, double)
{
    DROPS::Point3DCL p(q);
    TransformToDropOrig( p);
    return 3.*p[0]*p[0]*p[1] - p[1]*p[1]*p[1];
}
// ...and the corresponding right-hand side
double surf_rhs (const DROPS::Point3DCL& q, double)
{
    return surf_sol( q, 0.)*12.;
}
// some nice looking initial data
double surf_ini (const DROPS::Point3DCL& p, double)
{
    return 1. + std::sin( atan2( p[0] -P.get<DROPS::Point3DCL>("Levelset.PosDrop")[0], p[2] -P.get<DROPS::Point3DCL>("Levelset.PosDrop")[2]));
}


double mass_sol_p (const DROPS::Point3DCL& q, double)
{
    DROPS::Point3DCL p(q);
    TransformToDropOrig( p);
    return surf_sol( q, 0.)*exp(1. - p.norm_sq());
}
double mass_sol_n (const DROPS::Point3DCL& q, double)
{
    static double mu_ratio= P.get<double>("Transp.DiffPos")/P.get<double>("Transp.DiffNeg");
    return mass_sol_p( q, 0.)*mu_ratio;
}
double mass_rhs (const DROPS::Point3DCL& q, double)
{
    DROPS::Point3DCL p(q);
    TransformToDropOrig( p);
    return mass_sol_p( q, 0.)*(-4*p.norm_sq() + 18.);
}

DROPS::Point3DCL mass_grad_p (const DROPS::Point3DCL& q, double)
{
    DROPS::Point3DCL p(q), g;
    TransformToDropOrig( p);
    // g= [ 6xy, 3x^2 - 3y^2, 0 ]^T
    g[0]= 6*p[0]*p[1];
    g[1]= 3*(p[0]*p[0] - p[1]*p[1]);
    return exp(1. - p.norm_sq())*(-2*surf_sol( q, 0.)*p + g);
}
DROPS::Point3DCL mass_grad_n (const DROPS::Point3DCL& q, double)
{
    static double mu_ratio= P.get<double>("Transp.DiffPos")/P.get<double>("Transp.DiffNeg");
    return mu_ratio*mass_grad_p( q, 0.);
}

//@}

DROPS::Point3DCL RotateVelXZ (const DROPS::Point3DCL& q, double)
{
    DROPS::Point3DCL p(q);
    TransformToDropOrig( p);
    // return [z, 0, -x]^T
    double h= p[0];
    p[0]= p[2];
    p[1]= 0.;
    p[2]= -h;
    return 0.1*p;
}

static DROPS::RegisterScalarFunction regsca_surfrhs("rhsSurfNoVel", surf_rhs);
static DROPS::RegisterScalarFunction regsca_surfsol("solSurf", surf_sol);
static DROPS::RegisterScalarFunction regsca_massrhs("rhsMassNoVel", mass_rhs);
static DROPS::RegisterScalarFunction regsca_masssolp("solMassPos", mass_sol_p);
static DROPS::RegisterScalarFunction regsca_masssoln("solMassNeg", mass_sol_n);
static DROPS::RegisterVectorFunction regvec_gradp("gradMassPos", mass_grad_p);
static DROPS::RegisterVectorFunction regvec_gradn("gradMassNeg", mass_grad_n);
static DROPS::RegisterVectorFunction regvec_rotxz("RotatingFlowfield", RotateVelXZ);



/// L2 error on interface
template <typename DiscP1FunT>
double L2_error_interface (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const DROPS::BndDataCL<>& lsbnd,
    const DiscP1FunT& discsol, DROPS::instat_scalar_fun_ptr sol, double t)
{
    double d( 0.);
    const DROPS::Uint lvl = ls.GetLevel();
    DROPS::InterfaceTriangleCL triangle;
    DROPS::Quad5_2DCL<> qdiscsol, qsol;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        DROPS_FOR_TETRA_INTERFACE_BEGIN( *it, ls, lsbnd, triangle, tri) {
            qdiscsol.assign(  *it, &triangle.GetBary( tri), discsol);
            qsol.assign( *it, &triangle.GetBary( tri), sol, t);
            qdiscsol-= qsol;
            d+= DROPS::Quad5_2DCL<>( qdiscsol*qdiscsol).quad( triangle.GetAbsDet( tri));
        }
        DROPS_FOR_TETRA_INTERFACE_END
    }
    return std::sqrt( d);
}

namespace DROPS // for Strategy
{

void InitVel (const MultiGridCL& MG, VecDescCL& v, instat_vector_fun_ptr vf)
{
    const Uint lvl= v.GetLevel(),
               idx= v.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( MG, lvl, it)
        if (it->Unknowns.Exist( idx))
            DoFHelperCL<Point3DCL, VectorCL>::set( v.Data, it->Unknowns( idx),
                vf( it->GetCoord(), 0.));

    DROPS_FOR_TRIANG_CONST_EDGE( MG, lvl, it)
        if (it->Unknowns.Exist( idx))
            DoFHelperCL<Point3DCL, VectorCL>::set( v.Data, it->Unknowns( idx),
                vf( BaryCenter( it->GetVertex( 0)->GetCoord(), it->GetVertex( 1)->GetCoord()), 0.));
}


void  StatMassSurfTransportStrategy( MultiGridCL& MG, InstatNavierStokes2PhaseP2P1CL& Stokes, LsetBndDataCL& lsetbnddata, AdapTriangCL& adap)
// stationary coupled mass/surfactant transport
{
    std::cout << "Entering Strategy for stationary mass/surfactant transport...\n";
    InVecMap & tdvectormap = InVecMap::getInstance();
    InScaMap & tdscalarmap = InScaMap::getInstance();
    instat_vector_fun_ptr Flowfield = tdvectormap[P.get<std::string>("Transp.Flow")];
    instat_scalar_fun_ptr Reaction = tdscalarmap["Zero"]; // tdscalarmap["ReactionFct"];
    instat_scalar_fun_ptr massRhs = tdscalarmap[P.get<std::string>("Transp.Rhs")];
    instat_scalar_fun_ptr distance = tdscalarmap[P.get<std::string>("Transp.Levelset")];
    instat_scalar_fun_ptr massSol_p= tdscalarmap[P.get<std::string>("Transp.SolPos")];
    instat_scalar_fun_ptr massSol_n= tdscalarmap[P.get<std::string>("Transp.SolNeg")];
    /* unused
    instat_vector_fun_ptr massGrad_p= tdvectormap[P.get<std::string>("Transp.GradSolPos")];
    instat_vector_fun_ptr massGrad_n= tdvectormap[P.get<std::string>("Transp.GradSolNeg")];
    */

    instat_scalar_fun_ptr surfRhs = tdscalarmap[P.get<std::string>("SurfTransp.Rhs")];
    instat_scalar_fun_ptr surfSol = tdscalarmap[P.get<std::string>("SurfTransp.Sol")];


    cBndDataCL pBnd_c, pBnd_ct;
    read_BndData( pBnd_c, MG, P.get_child( "Transp.BoundaryData"));
    read_BndData( pBnd_ct,MG, P.get_child( "Transp.BoundaryDataTransformed"));
    cBndDataCL & Bnd_c(pBnd_c);
    cBndDataCL & Bnd_ct(pBnd_ct);

    DROPS::instat_scalar_fun_ptr sigmap = 0;
    SurfaceTensionCL sf( sigmap, Bnd_c);

    LevelsetP2CL & lset( * LevelsetP2CL::Create( MG, lsetbnddata, sf, false, 1, -1) ) ;

    // levelset wrt the previous time step:
    LevelsetP2CL & oldlset( * LevelsetP2CL::Create( MG, lsetbnddata, sf, false, 1, -1) ) ;
    //Prolongate and Restrict solution vector levelset from old mesh to new mesh after mesh adaptation:
    //always act on the same grid with possibly different interface position
    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);
    LevelsetRepairCL oldlsetrepair( oldlset);
    adap.push_back( &oldlsetrepair);
    lset.CreateNumbering( MG.GetLastLevel());
    oldlset.CreateNumbering( MG.GetLastLevel()); // interface at previous time step
    lset.Init( distance );
    oldlset.Init( distance);
    std::cout << "initial volume(abs value): " << lset.GetVolume() << std::endl;

    // grid refinement
    typedef UniformMarkingStrategyCL MarkerT;
    MarkerT marker( P.get<int>("Mesh.AdaptRef.FinestLevel") );
    adap.set_marking_strategy( &marker );
    adap.UpdateTriang();

    lset.Init( distance );
    oldlset.Init( distance);

    //VelocityContainer vel(Stokes.v,Stokes.GetBndData().Vel,MG);
    VelocityContainer vel(Flowfield);

    TransportP1XCL massTransp( MG, Bnd_c, Bnd_ct, vel, lsetbnddata, lset.Phi, oldlset.Phi,P,0,Reaction,massRhs);
    TransportXRepairCL transprepair(massTransp, MG.GetLastLevel());

    // index of the concentration wrt the interface at actual time step:
    MLIdxDescCL* cidx= &massTransp.idx;

    // index of the concentration wrt the interface at previous time step:
    MLIdxDescCL* cidx_old= &massTransp.oldidx;

    //This following Vector c_out is responsable for the communication from concentration to surface tension.
    //Before a new surface tension is computed c_out should be updated (via GetSolutionOnPart)
    //Important: This vector has to be kept in memory as long as the surface tension is computed!
    VecDescCL c_out;
    //c_in: s. c_out but c_in is only used for visualization until now
    VecDescCL c_in;
    IdxDescCL p1idx(P1_FE,Bnd_c,0);
    p1idx.CreateNumbering( MG.GetLastLevel(), MG, Bnd_c, &lset.Phi,&lsetbnddata);
    c_in.SetIdx( &p1idx);
    c_out.SetIdx( &p1idx);

    {
        adap.push_back(&transprepair);
        massTransp.CreateNumbering( MG.GetLastLevel(), cidx, cidx_old, lset.Phi, oldlset.Phi);
        massTransp.SetNewIdx();  // associate all matrices and vectors in massTransp with the FE space
        massTransp.ct.SetIdx( cidx);
        massTransp.c.SetIdx( cidx);
        massTransp.oldct.SetIdx( cidx_old);
        std::cout << massTransp.ct.Data.size() << " concentration unknowns,\n";
        massTransp.Init( massSol_n, massSol_p, 0);
    }

    Stokes.CreateNumberingVel( MG.GetLastLevel(), &Stokes.vel_idx);
    Stokes.v.SetIdx( &Stokes.vel_idx);
    InitVel( MG, Stokes.v, Flowfield);
    const double dt= P.get<double>("Time.FinalTime")/P.get<int>("Time.NumSteps");
    SurfactantcGP1CL surfTransp( MG, Stokes.GetBndData().Vel, P.get<double>("Time.Theta"), P.get<double>("SurfTransp.Visc"), &Stokes.v, *lset.PhiC, lset.GetBndData(),
                                 dt, P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"), P.get<double>("SurfTransp.XFEMReduced"));
    InterfaceP1RepairCL surf_repair( MG, *lset.PhiC, lset.GetBndData(), surfTransp.ic);
    {
        adap.push_back( &surf_repair);
        surfTransp.idx.CreateNumbering( MG.GetLastLevel(), MG, lset.PhiC, &lset.GetBndData());
        surfTransp.ic.SetIdx( &surfTransp.idx);
        surfTransp.Init( surfSol);
        std::cout << surfTransp.ic.Data.size() << " interface concentration unknowns" << std::endl;
    }
    const VectorCL surfExSol= surfTransp.ic.Data;

    // writer for vtk-format
    VTKOutCL * vtkwriter = nullptr;
    if (P.get<int>("VTK.Freq")){
        vtkwriter = new VTKOutCL(adap.GetMG(), "DROPS data", /*num time steps*/ 1,
                                 P.get<std::string>("VTK.VTKDir"), P.get<std::string>("VTK.VTKName"),
                                 P.get<std::string>("VTK.TimeFileName"),
                                 P.get<int>("VTK.Binary"),
                                 P.get<int>("VTK.UseOnlyP1"),
                                 false, /* P2DG */
                                 -1,  /* <- level */
                                 P.get<int>("VTK.ReUseTimeFile") );

        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "level-set") );

        vtkwriter->Register( make_VTKScalar( massTransp.GetSolution( massTransp.ct,false), "TransConcentration") );
        vtkwriter->Register( make_VTKScalar( massTransp.GetSolution( c_out,false), "XConcentrationPos") );
        vtkwriter->Register( make_VTKScalar( massTransp.GetSolution( c_in,false), "XConcentrationNeg") );
        vtkwriter->Register( make_VTKIfaceScalar( MG, surfTransp.ic, "SurfConcentration") );
    }

    // setup system
    massTransp.SetupInstatSystem( massTransp.A, massTransp.cplA, massTransp.M, massTransp.cplM, massTransp.C, massTransp.cplC,
            massTransp.b, 0.);

    surfTransp.Update(); // associate all matrices in surfTransp with the FE space and setup matrices
    VecDescCL surf_rhs( &surfTransp.idx);
    SetupInterfaceRhsP1( MG, &surf_rhs, lset.Phi, lset.GetBndData(), surfRhs);
    MatDescCL R, C, R_i, C_i;
    const double k_a[2]= { P.get<double>("SurfTransp.AdsorpPos"), P.get<double>("SurfTransp.AdsorpNeg")};
    const double k_d[2]= { P.get<double>("SurfTransp.DesorpPos"), P.get<double>("SurfTransp.DesorpNeg")};
    SetupInterfaceSorptionP1X( MG, lset.Phi, lset.GetBndData(), &R, &C, &R_i, &C_i, &massTransp.idx.GetFinest(), &surfTransp.idx, k_a, k_d);

    MatrixCL A, A_i;
    A.LinComb  ( 1., massTransp.A.Data.GetFinest(), 1., massTransp.C.Data.GetFinest(), 1., R.Data);
    massTransp.b.Data+= massTransp.cplA.Data + massTransp.cplC.Data + massTransp.C.Data.GetFinest()*massTransp.ct.Data;
    A_i.LinComb( 1., surfTransp.A.Data,             1., surfTransp.C.Data,             1., R_i.Data);
    surf_rhs.Data+= surfTransp.C.Data*surfTransp.ic.Data;

    // free some memory
    massTransp.M.Data.clear();
    massTransp.A.Data.clear();
    massTransp.C.Data.clear();

    // create block system
    BlockMatrixCL bmat( &A, MUL, &C.Data, MUL, &C_i.Data, MUL, &A_i, MUL);
    VectorCL &u= massTransp.ct.Data;
    VectorCL &v= surfTransp.ic.Data;
    size_t nmass= A.num_rows(), nsurf= A_i.num_rows();
    VectorCL rhs( nmass+nsurf), x( nmass+nsurf);
    rhs[std::slice( 0, nmass, 1)]= massTransp.b.Data;
    rhs[std::slice( nmass, nsurf, 1)]= surf_rhs.Data;
    x[std::slice( 0, nmass, 1)]= u;
    x[std::slice( nmass, nsurf, 1)]= v;
    DummyExchangeBlockCL exBlock;
    exBlock.AttachTo( massTransp.idx.GetEx());
    exBlock.AttachTo( surfTransp.idx.GetEx());

    // solve system
    typedef SSORPcCL PcT;
    typedef BlockDiagPreCL<PcT,PcT> SPPcT;
    typedef GCRSolverCL<SPPcT> SPSolverT;
    PcT pcM, pcS;
    SPPcT spPc( pcM, pcS);
    SPSolverT solver(spPc, P.get<int>("Transp.Solver.Restart",100), P.get<int>("Transp.Solver.Iter"), P.get<double>("Transp.Solver.Tol"), /*rel*/ false, &std::cout);

    std::cout << "Solving coupled system..." << std::endl;
    solver.Solve( bmat, x, rhs, exBlock);
    std::cout << "resid = " << solver.GetResid() << " after " << solver.GetIter() << " iterations.\n";

    u= x[std::slice( 0, nmass, 1)];
    v= x[std::slice( nmass, nsurf, 1)];

    // compute errors
    std::cout << "=== bulk error ===\n";
    /* unused
    const double L2_mass= massTransp.CheckSolution( massSol_n, massSol_p, massGrad_n, massGrad_p, 0);
    */
    const double c_mean = massTransp.MeanDropConcentration();
    std::cout << "Mean concentration in drop: " << c_mean <<"\n";
    const double L2_surf= L2_error_interface( MG, lset.Phi, lsetbnddata, surfTransp.GetSolution(), surfSol, 0);
    const VectorCL surfDiff(v - surfExSol);
    const double H1_surf= std::sqrt( L2_surf*L2_surf + dot( surfTransp.A.Data*surfDiff, surfDiff));
    std::cout << "=== surface error ===\nL2 = " << L2_surf << "\nH1 = " << H1_surf << std::endl;

    // VTK output
    massTransp.GetSolutionOnPart( c_out, true , false);
    massTransp.GetSolutionOnPart( c_in, false , false);
    vtkwriter->Write(0);

    adap.set_marking_strategy(0);
    std::cout << std::endl;

    delete vtkwriter;
    delete &lset;
    delete &oldlset;
}

void  OnlyTransportStrategy( MultiGridCL& MG, LsetBndDataCL& lsetbnddata, AdapTriangCL& adap)    // do just the transport stuff
{
    InVecMap & tdvectormap = InVecMap::getInstance();
    InScaMap & tdscalarmap = InScaMap::getInstance();
    instat_vector_fun_ptr Flowfield = tdvectormap[P.get<std::string>("Transp.Flow")];
    instat_scalar_fun_ptr Reaction = tdscalarmap["ReactionFct"];
    instat_scalar_fun_ptr Rhs = tdscalarmap[P.get<std::string>("Transp.Rhs")];
    instat_scalar_fun_ptr Initialcneg = tdscalarmap[P.get<std::string>("Transp.InitialConcNeg")];
    instat_scalar_fun_ptr Initialcpos = tdscalarmap[P.get<std::string>("Transp.InitialConcPos")];
    instat_scalar_fun_ptr distance = tdscalarmap[P.get<std::string>("Transp.Levelset")];

    cBndDataCL pBnd_c, pBnd_ct;
    read_BndData( pBnd_c, MG, P.get_child( "Transp.BoundaryData"));
    read_BndData( pBnd_ct,MG, P.get_child( "Transp.BoundaryDataTransformed"));
    cBndDataCL & Bnd_c(pBnd_c);
    cBndDataCL & Bnd_ct(pBnd_ct);

    DROPS::instat_scalar_fun_ptr sigmap = 0;
    SurfaceTensionCL sf( sigmap, Bnd_c);

    LevelsetP2CL & lset( * LevelsetP2CL::Create( MG, lsetbnddata, sf, false, 1, -1) ) ;

    // levelset wrt the previous time step:
    LevelsetP2CL & oldlset( * LevelsetP2CL::Create( MG, lsetbnddata, sf, false, 1, -1) ) ;
    //Prolongate and Restrict solution vector levelset from old mesh to new mesh after mesh adaptation:
    //always act on the same grid with possibly different interface position
    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);
    LevelsetRepairCL oldlsetrepair( oldlset);
    adap.push_back( &oldlsetrepair);
    lset.CreateNumbering( MG.GetLastLevel());
    oldlset.CreateNumbering( MG.GetLastLevel()); // interface at previous time step
    SetInitialLevelsetConditions( lset, MG, P);
    SetInitialLevelsetConditions( oldlset, MG, P);
    lset.Init( distance );
    oldlset.Init( distance);
    DisplayDetailedGeom( MG);
    //const double Vol= lset.GetVolume(); //0.5 * 0.125 * M_PI; //EllipsoidCL::GetVolume();
    std::cout << "initial volume(abs value): " << lset.GetVolume() << std::endl;

    //VelocityContainer vel(Stokes.v,Stokes.GetBndData().Vel,MG);
    VelocityContainer vel(Flowfield);

    TransportP1XCL massTransp( MG, Bnd_c, Bnd_ct, vel, lsetbnddata, lset.Phi, oldlset.Phi,P,0,Reaction,Rhs);
    TransportXRepairCL transprepair(massTransp, MG.GetLastLevel());

    // index of the concentration wrt the interface at actual time step:
    MLIdxDescCL* cidx= &massTransp.idx;

    // index of the concentration wrt the interface at previous time step:
    MLIdxDescCL* cidx_old= &massTransp.oldidx;

    //This following Vector c_out is responsable for the communication from concentration to surface tension.
    //Before a new surface tension is computed c_out should be updated (via GetSolutionOnPart)
    //Important: This vector has to be kept in memory as long as the surface tension is computed!
    VecDescCL c_out;
    //c_in: s. c_out but c_in is only used for visualization until now
    VecDescCL c_in;
    IdxDescCL p1idx(P1_FE,Bnd_c,0);
    p1idx.CreateNumbering( MG.GetLastLevel(), MG, Bnd_c, &lset.Phi,&lsetbnddata);
    c_in.SetIdx( &p1idx);
    c_out.SetIdx( &p1idx);

    {
        adap.push_back(&transprepair);
        massTransp.CreateNumbering( MG.GetLastLevel(), cidx, cidx_old, lset.Phi, oldlset.Phi);
        massTransp.ct.SetIdx( cidx);
        massTransp.c.SetIdx( cidx);
        massTransp.oldct.SetIdx( cidx_old);
        std::cout << massTransp.ct.Data.size() << " concentration unknowns,\n";

        if (!ReadInitialConditionFromFile(P))
          massTransp.Init( Initialcneg, Initialcpos, 0);
        else
          ReadFEFromFile( massTransp.ct, MG, P.get<std::string>("Restart.InputData")+"concentrationTransf");

        massTransp.GetSolutionOnPart( c_out, true , false);
        massTransp.GetSolutionOnPart( c_in, false , false);
        double c_mean = massTransp.MeanDropConcentration();
        std::cout << "START:: Mean concentration in drop: " << std::setprecision(12) << c_mean <<"\n";
    }

    // for serialization of geometry and numerical data
    if (P.get("Transp.Enable", 0))
      std::cout << "WARNING: mass transport data is not serialized, yet!" << std::endl;

    // Initialize Ensight6 output
    //Update c from ct
    //massTransp.TransformWithScaling(massTransp.ct, massTransp.c, 1.0/massTransp.GetHenry(true), 1.0/massTransp.GetHenry(false));
    // Output-Registrations:
    Ensight6OutCL* ensight = nullptr;
    if (P.get<int>("Ensight.Freq",0)){
        std::string ensf( P.get<std::string>("Ensight.EnsDir") + "/" + P.get<std::string>("Ensight.EnsCase"));
        ensight = new Ensight6OutCL( P.get<std::string>("Ensight.EnsCase") + ".case",
                                     P.get<int>("Time.NumSteps")/P.get("Ensight.Freq", 0)+1,
                                     P.get<int>("Ensight.Binary"));
        ensight->Register( make_Ensight6Geom      ( MG, MG.GetLastLevel(),
                                                    P.get<std::string>("Ensight.GeomName"),
                                                    ensf + ".geo", true));
        ensight->Register( make_Ensight6Scalar    ( lset.GetSolution(),      "Levelset",
                                                    ensf + ".scl", true));
        ensight->Register( make_Ensight6Scalar( massTransp.GetSolution( massTransp.ct,true),
                                                "TransConc",     ensf + ".ct",  true));
        ensight->Register( make_Ensight6P1XScalar( MG, lset.Phi, massTransp.ct, "XTransConcentration",
                                                   ensf + ".xconc", true));
        ensight->Write(0);
    }

    // writer for vtk-format
    VTKOutCL * vtkwriter = nullptr;
    if (P.get<int>("VTK.Freq")){
        vtkwriter = new VTKOutCL(adap.GetMG(), "DROPS data",
                                 P.get<int>("Time.NumSteps")/P.get<int>("VTK.Freq")+1,
                                 P.get<std::string>("VTK.VTKDir"), P.get<std::string>("VTK.VTKName"),
                                 P.get<std::string>("VTK.TimeFileName"),
                                 P.get<int>("VTK.Binary"),
                                 P.get<int>("VTK.UseOnlyP1"),
                                 false, /* P2DG */
                                 -1,  /* <- level */
                                 P.get<int>("VTK.ReUseTimeFile") );

        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "level-set") );

        vtkwriter->Register( make_VTKScalar( massTransp.GetSolution( massTransp.ct,false), "TransConcentration") );
        vtkwriter->Register( make_VTKScalar( massTransp.GetSolution( c_out,false), "XConcentrationPos") );
        vtkwriter->Register( make_VTKScalar( massTransp.GetSolution( c_in,false), "XConcentrationNeg") );
        vtkwriter->Write(0);
    }

    // massTransp.CheckSolution(Solutioncneg,Solutioncpos,0);
    double cmean_old = massTransp.MeanDropConcentration();

    typedef DistMarkingStrategyCL MarkerT;
    MarkerT marker( lset, P.get<double>("Mesh.AdaptRef.Width"),
                          P.get<int>("Mesh.AdaptRef.CoarsestLevel" ),
                          P.get<int>("Mesh.AdaptRef.FinestLevel") );
    adap.set_marking_strategy( &marker );

    const double dt= P.get<double>("Time.FinalTime")/P.get<int>("Time.NumSteps");
    for (int step= 1; step<=P.get<int>("Time.NumSteps"); ++step)
    {
        std::cout << "============================================================ step " << std::setw(8) << step << "  /  " << std::setw(8) << P.get<int>("Time.NumSteps") << std::endl;
        double c_mean = massTransp.MeanDropConcentration();
        double t= dt * step;
        std::cout << "Mean concentration in drop: " << c_mean <<"\n";
        if(step > 5 && std::abs(cmean_old - c_mean)/(std::abs(cmean_old) * dt) < 1e-5 ){
          std::cout << "I think I found a stationary solution! " << std::endl;
          break;
        }
        cmean_old = c_mean;

        // grid modification
        bool doGridMod= P.get<int>("Mesh.AdaptRef.Freq") && step%P.get<int>("Mesh.AdaptRef.Freq") == 0;
        if (doGridMod) {
            adap.UpdateTriang();
        }

        massTransp.DoStep( t);

        massTransp.GetSolutionOnPart( c_out, true , false);
        massTransp.GetSolutionOnPart( c_in, false , false);

        bool ensightoutnow = ensight && step%P.get<int>("Ensight.Freq")==0;
        bool vtkoutnow = vtkwriter && (step%P.get<int>("VTK.Freq")==0 || step < 20);
        if (ensightoutnow)
            ensight->Write(t);
        if (vtkoutnow)
            vtkwriter->Write(t);
    }
    adap.set_marking_strategy(0);
    std::cout << std::endl;

    delete vtkwriter;
    delete ensight;
    delete &lset;
    delete &oldlset;
}

void Strategy( InstatNavierStokes2PhaseP2P1CL& Stokes,  LsetBndDataCL& lsetbnddata, AdapTriangCL& adap)
// flow control
{
    MultiGridCL& MG= Stokes.GetMG();

    InVecMap & tdvectormap = InVecMap::getInstance();
    InScaMap & tdscalarmap = InScaMap::getInstance();

    instat_vector_fun_ptr Flowfield = tdvectormap["ZeroVel"];
    instat_scalar_fun_ptr Reaction = tdscalarmap["ReactionFct"];
    instat_scalar_fun_ptr Rhs = tdscalarmap[P.get<std::string>("Transp.Rhs")];
    instat_scalar_fun_ptr Initialcneg = tdscalarmap[P.get<std::string>("Transp.InitialConcNeg")];
    instat_scalar_fun_ptr Initialcpos = tdscalarmap[P.get<std::string>("Transp.InitialConcPos")];

    instat_scalar_fun_ptr distance = tdscalarmap[P.get<std::string>("Transp.Levelset")];

    instat_scalar_fun_ptr sigmap = tdscalarmap[P.get<std::string>("NavStokes.Coeff.SurfTens.VarTensionFunc")];

    cBndDataCL pBnd_c, pBnd_ct;
    read_BndData( pBnd_c, MG, P.get_child( "Transp.BoundaryData"));
    read_BndData( pBnd_ct,MG, P.get_child( "Transp.BoundaryDataTransformed"));
    cBndDataCL & Bnd_c(pBnd_c);
    cBndDataCL & Bnd_ct(pBnd_ct);


    // initialization of surface tension
    SurfaceTensionCL * sf;
    sf = new SurfaceTensionCL( sigmap, Bnd_c);

    double cp=0., coeffC[5];
    //coefficients for ansatz of var. surface tension
    //coeffC[0]= 1.625; coeffC[1]= 0.0; coeffC[2]= 0.0; coeffC[3]= coeffC[4]= 0.;
    coeffC[0]= 1.625; coeffC[1]= -28.07768; coeffC[2]= 222.7858; coeffC[3]= coeffC[4]= 0.;

    if (!P.get<bool>("NavStokes.Coeff.SurfTens.ConcentrationDep"))
        sf->SetInputMethod(Sigma_X);
    else{
        sf->SetInputMethod(Sigma_C);
        sf->SetCoeff(coeffC, cp);
    }

    LevelsetP2CL & lset( * LevelsetP2CL::Create( MG, lsetbnddata, *sf, P.get_child("Levelset")) );

    // levelset wrt the previous time step:
    LevelsetP2CL & oldlset( * LevelsetP2CL::Create( MG, lsetbnddata, *sf, P.get_child("Levelset")) );
    //Prolongate and Restrict solution vector levelset from old mesh to new mesh after mesh adaptation:
    //always act on the same grid with possibly different interface position
    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);
    LevelsetRepairCL oldlsetrepair( oldlset);
    adap.push_back( &oldlsetrepair);
    //Prolongate and Restrict solution vector for the velocity from old mesh to new mesh after mesh adaptation:
    //always act on the same grid with possibly different interface position
    VelocityRepairCL velrepair( Stokes);
    adap.push_back( &velrepair);
    PressureRepairCL prrepair( Stokes, lset);
    adap.push_back( &prrepair);
    MLIdxDescCL* vidx= &Stokes.vel_idx;
    IdxDescCL old_vidx(vecP2_FE);
    MLIdxDescCL* pidx= &Stokes.pr_idx;

    lset.CreateNumbering( MG.GetLastLevel());
    VecDescCL old_v(&old_vidx);

    //Prolongate and Restrict solution vector old_v from old mesh to new mesh after mesh adaptation:
    DROPS::VelTranspRepairCL old_vrepair(old_v, MG, Stokes.GetBndData().Vel, old_vidx, 0.);
    adap.push_back( &old_vrepair);
    oldlset.CreateNumbering( MG.GetLastLevel()); // interface at previous time step

    lset.SetSurfaceForce( SF_ImprovedLBVar);

    ParamCL PSolver= P.get_child("CouplingSolver.NavStokesSolver.OseenSolver");
    if ( StokesSolverFactoryHelperCL().VelMGUsed(PSolver)){
        Stokes.SetNumVelLvl ( Stokes.GetMG().GetNumLevel());
        lset.SetNumLvl(Stokes.GetMG().GetNumLevel());
    }
    if ( StokesSolverFactoryHelperCL().PrMGUsed(PSolver)){
        Stokes.SetNumPrLvl  ( Stokes.GetMG().GetNumLevel());
        lset.SetNumLvl(Stokes.GetMG().GetNumLevel());
    }
    lset.CreateNumbering( MG.GetLastLevel());

    SetInitialLevelsetConditions( lset, MG, P);
    SetInitialLevelsetConditions( oldlset, MG, P);
    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, &lset);
    old_vidx.CreateNumbering( MG.GetLastLevel(), MG, Stokes.GetBndData().Vel);
    old_v.SetIdx  ( &old_vidx);

    Stokes.SetIdx();
    Stokes.v.SetIdx  ( vidx);
    Stokes.p.SetIdx  ( pidx);
    Stokes.InitVel( &Stokes.v, Flowfield);
    SetInitialConditions( Stokes, lset, MG, P);
    InitVel(MG, old_v, Flowfield);

    lset.Init( distance);
    oldlset.Init( distance);

    DisplayDetailedGeom( MG);
    DisplayUnks(Stokes, lset, MG);

    const double Vol= EllipsoidCL::GetVolume();
    lset.InitVolume( Vol);

    std::cout << "initial volume: " << lset.GetVolume()/Vol << std::endl;
    lset.AdjustVolume();
    std::cout << "initial lset volume adjustment:\n";
    lset.GetVolumeAdjuster()->DebugOutput( std::cout);

    VelocityContainer vel(Stokes.v,Stokes.GetBndData().Vel,MG);
    //VelocityContainer vel(Flowfield);
    TransportP1XCL massTransp( MG, Bnd_c, Bnd_ct, vel, lsetbnddata, lset.Phi, oldlset.Phi,P,0,Reaction,Rhs);
    TransportXRepairCL transprepair(massTransp, MG.GetLastLevel());

    // index of the concentration wrt the interface at actual time step:
    MLIdxDescCL* cidx= &massTransp.idx;

    // index of the concentration wrt the interface at previous time step:
    MLIdxDescCL* cidx_old= &massTransp.oldidx;

    //This following Vector c_out is responsible for the communication from concentration to surface tension.
    //Before a new surface tension is computed c_out should be updated (via GetSolutionOnPart)
    //Important: This vector has to be kept in memory as long as the surface tension is computed!
    VecDescCL c_out;
    //c_in: s. c_out but c_in is only used for visualization until now
    VecDescCL c_in;
    IdxDescCL p1idx(P1_FE,Bnd_c,0);
    /// \todo for periodic stuff: matching function here
    p1idx.CreateNumbering( MG.GetLastLevel(), MG, Bnd_c, 0, &lset.Phi,&lsetbnddata);
    c_in.SetIdx( &p1idx);
    c_out.SetIdx( &p1idx);
    if (P.get<int>("Transp.Enable")) {
        adap.push_back(&transprepair);
        massTransp.CreateNumbering( MG.GetLastLevel(), cidx, cidx_old, lset.Phi, oldlset.Phi);
        massTransp.ct.SetIdx( cidx);
        massTransp.c.SetIdx( cidx);
        massTransp.oldct.SetIdx( cidx_old);
        if (!ReadInitialConditionFromFile(P))
          massTransp.Init( Initialcpos, Initialcpos, true);
        else
          ReadFEFromFile( massTransp.ct, MG, P.get<std::string>("Restart.InputData")+"concentrationTransf");

        std::cout << massTransp.ct.Data.size() << " concentration unknowns,\n";

        if (P.get<bool>("NavStokes.Coeff.SurfTens.ConcentrationDep")){
            massTransp.GetSolutionOnPart( c_out, true , false);
            massTransp.GetSolutionOnPart( c_in, false , false);
//            P1XtoP1 (*massTransp.c.RowIdx, massTransp.c.Data, p1idx, c_out.Data, c_in.Data, lset.Phi, MG);
            sf->SetConcentration(&c_out);
            sf->SetInputMethod(Sigma_C);
            sf->SetTime(0.);
        }

        if (!ReadInitialConditionFromFile(P))
          massTransp.Init( Initialcneg, Initialcpos, true);
        else
          ReadFEFromFile( massTransp.ct, MG, P.get<std::string>("Restart.InputData")+"concentrationTransf");
        massTransp.GetSolutionOnPart( c_out, true , false);
        massTransp.GetSolutionOnPart( c_in, false , false);
//        P1XtoP1 (*massTransp.c.RowIdx, massTransp.c.Data, p1idx, c_out.Data, c_in.Data, lset.Phi, MG);

        double c_mean = massTransp.MeanDropConcentration();
        std::cout << "START:: Mean concentration in drop: " << std::setprecision(12) << c_mean <<"\n";
    }

    const double dt= P.get<double>("Time.FinalTime")/P.get<int>("Time.NumSteps");
    /// \todo rhs beruecksichtigen
    SurfactantcGP1CL surfTransp( MG, Stokes.GetBndData().Vel, P.get<double>("Time.Theta"), P.get<double>("SurfTransp.Visc"), &Stokes.v, lset.Phi, lset.GetBndData(), dt, P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"), P.get<double>("SurfTransp.XFEMReduced"));
    InterfaceP1RepairCL surf_repair( MG, lset.Phi, lset.GetBndData(), surfTransp.ic);
    if (P.get<int>("SurfTransp.Enable"))
    {
        adap.push_back( &surf_repair);
        surfTransp.idx.CreateNumbering( MG.GetLastLevel(), MG, &lset.Phi, &lset.GetBndData());
        std::cout << "Surfactant transport: NumUnknowns: " << surfTransp.idx.NumUnknowns() << std::endl;
        surfTransp.ic.SetIdx( &surfTransp.idx);
        surfTransp.Init( &surf_sol);
    }

    // Stokes-Solver
    StokesSolverFactoryCL<InstatNavierStokes2PhaseP2P1CL> stokessolverfactory(Stokes, PSolver, P.get_child("Time"));
    StokesSolverBaseCL* stokessolver = stokessolverfactory.CreateStokesSolver();
//     StokesSolverAsPreCL pc (*stokessolver1, 1);
//     GCRSolverCL<StokesSolverAsPreCL> gcr(pc, C.stk_OuterIter, C.stk_OuterIter, C.stk_OuterTol, /*rel*/ false);
//     BlockMatrixSolverCL<GCRSolverCL<StokesSolverAsPreCL> >* stokessolver =
//             new BlockMatrixSolverCL<GCRSolverCL<StokesSolverAsPreCL> > (gcr);

    // Navier-Stokes-Solver
    NSSolverBaseCL<InstatNavierStokes2PhaseP2P1CL>* navstokessolver = 0;
    if (P.get<double>("CouplingSolver.NavStokesSolver.Nonlinear")==0.0)
        navstokessolver = new NSSolverBaseCL<InstatNavierStokes2PhaseP2P1CL>(Stokes, *stokessolver);
    else
        navstokessolver = new AdaptFixedPtDefectCorrCL<InstatNavierStokes2PhaseP2P1CL>(Stokes, *stokessolver, P.get<int>("CouplingSolver.NavStokesSolver.Iter"), P.get<double>("CouplingSolver.NavStokesSolver.Tol"), P.get<double>("CouplingSolver.NavStokesSolver.Reduction"));
    // Level-Set-Solver
#ifndef _PAR
    SSORPcCL ssorpc;
    GMResSolverCL<SSORPcCL>* gm = new GMResSolverCL<SSORPcCL>( ssorpc, 100, P.get<int>("CouplingSolver.LevelsetSolver.Iter"), P.get<double>("CouplingSolver.LevelsetSolver.Tol"));
#else
    ParJac0CL jacparpc( *lidx);
    ParPreGMResSolverCL<ParJac0CL>* gm = new ParPreGMResSolverCL<ParJac0CL>
           (/*restart*/100, P.get<int>("CouplingSolver.LevelsetSolver.Iter"), P.get<double>("CouplingSolver.LevelsetSolver.Tol"), *lidx, jacparpc,/*rel*/true, /*acc*/ true, /*modGS*/false, LeftPreconditioning, /*parmod*/true);
#endif

    LevelsetModifyCL lsetmod( P.get<int>("Levelset.Reparam.Freq"), P.get<int>("Levelset.Reparam.Method"), P.get<double>("Levelset.Reparam.MaxGrad"), P.get<double>("Levelset.Reparam.MinGrad"));

    UpdateProlongationCL<Point3DCL> PVel( Stokes.GetMG(), stokessolverfactory.GetPVel(), &Stokes.vel_idx, &Stokes.vel_idx);
    adap.push_back( &PVel);
    UpdateProlongationCL<double> PPr ( Stokes.GetMG(), stokessolverfactory.GetPPr(), &Stokes.pr_idx, &Stokes.pr_idx);
    adap.push_back( &PPr);
    UpdateProlongationCL<double> PLset( lset.GetMG(), lset.GetProlongation(), lset.idxC, lset.idxC);
    adap.push_back( &PLset);
    Stokes.P_ = stokessolverfactory.GetPVel();

    // For a two-level MG-solver: P2P1 -- P2P1X;
//     MakeP1P1XProlongation ( Stokes.vel_idx.NumUnknowns(), Stokes.pr_idx.NumUnknowns(),
//         Stokes.pr_idx.GetFinest().GetXidx().GetNumUnknownsStdFE(),
//         stokessolverfactory.GetPVel()->GetFinest(), stokessolverfactory.GetPPr()->GetFinest());

    // Time discretisation + coupling
    TimeDisc2PhaseCL* timedisc= CreateTimeDisc(Stokes, lset, navstokessolver, gm, P, lsetmod);

    if (P.get<int>("Time.NumSteps") != 0){
        timedisc->SetTimeStep( dt);
        timedisc->SetSchurPrePtr( stokessolverfactory.GetSchurPrePtr() );
    }

    if (P.get<double>("CouplingSolver.NavStokesSolver.Nonlinear")!=0.0 || P.get<int>("Time.NumSteps") == 0) {
        stokessolverfactory.SetMatrixA( &navstokessolver->GetAN()->GetFinest());
            //for Stokes-MGM
        stokessolverfactory.SetMatrices( navstokessolver->GetAN(), &Stokes.B.Data,
                                         &Stokes.M.Data, &Stokes.prM.Data, &Stokes.pr_idx);
    }
    else {
        stokessolverfactory.SetMatrixA( &timedisc->GetUpperLeftBlock()->GetFinest());
            //for Stokes-MGM
        stokessolverfactory.SetMatrices( timedisc->GetUpperLeftBlock(), &Stokes.B.Data,
                                         &Stokes.M.Data, &Stokes.prM.Data, &Stokes.pr_idx);
    }

    std::ofstream* infofile = 0;
    IF_MASTER {
        infofile = new std::ofstream ((P.get<std::string>("VTK.VTKName")+".info").c_str());
    }
    IFInfo.Init(infofile);
    IFInfo.WriteHeader();

    if (P.get<int>("Time.NumSteps") == 0)
        SolveStatProblem( Stokes, lset, *navstokessolver);

    // for serialization of geometry and numerical data
    if (P.get<int>("Transp.Enable"))
      std::cout << "WARNING: mass transport data is not serialized, yet!" << std::endl;
    TwoPhaseStoreCL<InstatNavierStokes2PhaseP2P1CL> ser(MG, Stokes, lset, nullptr,
                                                        P.get<std::string>("Restart.OutputData"),
                                                        P.get<std::string>("Restart.OutputGrid"),
                                                        P.get<int>("Restart.OutputOverwrite"),
                                                        P.get<int>("Restart.Binary"));




    // Output-Registrations:
    Ensight6OutCL* ensight = nullptr;
    if (P.get<int>("Ensight.Freq")){
        // Initialize Ensight6 output
        std::string ensf( P.get<std::string>("Ensight.EnsDir") + "/" + P.get<std::string>("Ensight.EnsCase"));
        ensight = new Ensight6OutCL( P.get<std::string>("Ensight.EnsCase") + ".case",
                                     P.get<int>("Time.NumSteps")/P.get<int>("Ensight.Freq")+1,
                                     P.get<int>("Ensight.Binary"));
        ensight->Register( make_Ensight6Geom      ( MG, MG.GetLastLevel(), P.get<std::string>("Ensight.GeomName"),
                                                    ensf + ".geo", true));
        ensight->Register( make_Ensight6Scalar    ( lset.GetSolution(),      "Levelset",      ensf + ".scl", true));
        ensight->Register( make_Ensight6Scalar    ( Stokes.GetPrSolution(),  "Pressure",      ensf + ".pr",  true));
        ensight->Register( make_Ensight6Vector    ( Stokes.GetVelSolution(), "Velocity",      ensf + ".vel", true));
        ensight->Register( make_Ensight6Scalar    ( ScalarFunAsP2EvalCL( sigmap, 0., &MG, MG.GetLastLevel()),
                                                    "Surfaceforce",  ensf + ".sf",  true));

        if (P.get<int>("Transp.Enable")) {
            ensight->Register( make_Ensight6Scalar( massTransp.GetSolution(),"Concentration", ensf + ".c",   true));
            ensight->Register( make_Ensight6P1XScalar( MG, lset.Phi, massTransp.ct,
                                                      "XTransConcentration",   ensf + ".xconc", true));

        }
        if (P.get<int>("SurfTransp.Enable")) {
            ensight->Register( make_Ensight6IfaceScalar( MG, surfTransp.ic,  "InterfaceSol",  ensf + ".sur", true));
        }

#ifndef _PAR
        if (Stokes.UsesXFEM())
            ensight->Register( make_Ensight6P1XScalar( MG, lset.Phi, Stokes.p, "XPressure",   ensf + ".pr", true));
#endif
        ensight->Write( Stokes.v.t);
    }

    // writer for vtk-format
    VTKOutCL * vtkwriter = nullptr;
    if (P.get<int>("VTK.Freq")){
        vtkwriter = new VTKOutCL(adap.GetMG(), "DROPS data",
                                 P.get<int>("Time.NumSteps")/P.get<int>("VTK.Freq")+1,
                                 P.get<std::string>("VTK.VTKDir"), P.get<std::string>("VTK.VTKName"),
                                 P.get<std::string>("VTK.TimeFileName"),
                                 P.get<int>("VTK.Binary"),
                                 P.get<int>("VTK.UseOnlyP1"),
                                 false, /* P2DG */
                                 -1,  /* <- level */
                                 P.get<int>("VTK.ReUseTimeFile") );
        vtkwriter->Register( make_VTKVector( Stokes.GetVelSolution(), "velocity") );
        vtkwriter->Register( make_VTKScalar( Stokes.GetPrSolution(), "pressure") );
        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "level-set") );

        if (P.get<int>("Transp.Enable")) {
            vtkwriter->Register( make_VTKScalar( massTransp.GetSolution( massTransp.ct,false), "TransConcentration") );
            vtkwriter->Register( make_VTKScalar( massTransp.GetSolution( c_out,false), "XConcentrationPos") );
            vtkwriter->Register( make_VTKScalar( massTransp.GetSolution( c_in,false), "XConcentrationNeg") );

        }

        if (P.get<int>("SurfTransp.Enable")) {
            vtkwriter->Register( make_VTKIfaceScalar( MG, surfTransp.ic,  "InterfaceSol"));
        }
        vtkwriter->Write(Stokes.v.t);
    }

    typedef DistMarkingStrategyCL MarkerT;
    MarkerT marker( lset, P.get<double>("Mesh.AdaptRef.Width"),
                          P.get<int>("Mesh.AdaptRef.CoarsestLevel" ),
                          P.get<int>("Mesh.AdaptRef.FinestLevel") );
    adap.set_marking_strategy( &marker );

    for (int step= 1; step<=P.get<int>("Time.NumSteps"); ++step)
    {
        std::cout << "============================================================ step " << step << std::endl;
        double c_mean = massTransp.MeanDropConcentration();
        double t= Stokes.v.t;
        std::cout << "Mean concentration in drop: " << c_mean <<"\n";

        IFInfo.Update( lset, Stokes.GetVelSolution());
//         IFInfo.Write(Stokes.t, c_mean);
        IFInfo.Write(Stokes.v.t);

        if (P.get<int>("SurfTransp.Enable")) surfTransp.InitOld();
        timedisc->DoStep( P.get<int>("CouplingSolver.Iter"));

//         if (P.get("Transp.DoTransp", 0)) massTransp.DoStep( step*dt);
        if (P.get<int>("SurfTransp.Enable")) {
            surfTransp.DoStep( step*dt);
            BndDataCL<> ifbnd( 0);
            std::cout << "surfactant on \\Gamma: " << Integral_Gamma( MG, lset.Phi, lset.GetBndData(), make_P1Eval(  MG, ifbnd, surfTransp.ic)) << '\n';
        }

        // WriteMatrices( Stokes, step);

        // grid modification
        bool doGridMod= P.get<int>("Mesh.AdaptRef.Freq") && step%P.get<int>("Mesh.AdaptRef.Freq") == 0;
        if (doGridMod) {
            adap.UpdateTriang();
        }

        if (P.get<int>("Transp.Enable")) {
            massTransp.DoStep( t);
            old_vrepair.SetTime(t);

        //}

            if (P.get<bool>("NavStokes.Coeff.SurfTens.ConcentrationDep")){
                sf->SetConcentration(&c_out);
                sf->SetInputMethod(Sigma_C);
                sf->SetTime(t);
            }
        }

        timedisc->Update();

        //Update c from ct
//        massTransp.TransformWithScaling(massTransp.ct, massTransp.c, 1.0/massTransp.GetHenry(true), 1.0/massTransp.GetHenry(false));

		massTransp.GetSolutionOnPart( c_out, true , false);
		massTransp.GetSolutionOnPart( c_in, false , false);

/// \todo for periodic stuff: matching function here
//        p1idx.CreateNumbering( MG.GetLastLevel(), MG, Bnd_c, 0, &lset.Phi,&lsetbnddata);
//        P1XtoP1 (*massTransp.c.RowIdx, massTransp.c.Data, p1idx, c_out.Data, c_in.Data, lset.Phi, MG);

        if (ensight && step%P.get<int>("Ensight.Freq")==0)
            ensight->Write( Stokes.v.t);
        if (vtkwriter && step%P.get<int>("VTK.Freq")==0)
            vtkwriter->Write(Stokes.v.t);
//        if (C.rst_Serialization && step%C.rst_Serialization==0)
//            ser.Write();
    }
    IFInfo.Update( lset, Stokes.GetVelSolution());
    IFInfo.Write(Stokes.v.t);
    adap.set_marking_strategy(0);
    std::cout << std::endl;
    delete timedisc;
    delete navstokessolver;
    delete stokessolver;
    delete gm;
    if (infofile) delete infofile;
//     delete stokessolver1;
    if (vtkwriter) delete vtkwriter;
    if (ensight) delete ensight;
    delete &lset;
    delete &oldlset;
    if (sf) delete sf;

}

} // end of namespace DROPS

int main (int argc, char** argv)
{
#ifdef _PAR
    DROPS::ProcInitCL procinit(&argc, &argv);
#endif
  try
  {
#ifdef _PAR
    DROPS::ParMultiGridInitCL pmginit;
#endif
    DROPS::read_parameter_file_from_cmdline( P, argc, argv, "../../param/transport/ns_transp/risingbutanoldroplet.json");
    P.put_if_unset<std::string>("VTK.TimeFileName",P.get<std::string>("VTK.VTKName"));
    std::cout << P << std::endl;

    DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"), P.get<std::vector<std::string> >("General.DynamicLibs") );
    if (P.get<int>("General.ProgressBar"))
        DROPS::ProgressBarTetraAccumulatorCL::Activate();

    DROPS::MultiGridCL* mg= 0;
    DROPS::StokesVelBndDataCL velbnddata;
    DROPS::StokesPrBndDataCL  prbnddata;
    DROPS::LsetBndDataCL lsetbnddata;

    std::unique_ptr<DROPS::MGBuilderCL> builder( DROPS::make_MGBuilder( P));
    mg = new DROPS::MultiGridCL( *builder);
    if (P.exists("Mesh.PeriodicBnd"))
        DROPS::read_PeriodicBoundaries( *mg, P.get_child("Mesh.PeriodicBnd"));
    read_BndData( velbnddata, *mg, P.get_child( "NavStokes.BoundaryData.Velocity"));
    read_BndData( prbnddata,  *mg, P.get_child( "NavStokes.BoundaryData.Pressure"));
    read_BndData( lsetbnddata,*mg, P.get_child( "Levelset.BoundaryData"));
    DROPS::StokesBndDataCL bnddata( velbnddata, prbnddata);

    std::cout << "Generated MG of " << mg->GetLastLevel() << " levels." << std::endl;

    using DROPS::EllipsoidCL;
    using DROPS::AdapTriangCL;
    using DROPS::DistMarkingStrategyCL;

    EllipsoidCL::Init(P.get<DROPS::Point3DCL>("Levelset.PosDrop"), P.get<DROPS::Point3DCL>("Levelset.RadDrop"));
    AdapTriangCL adap( *mg, 0,
                      (!ReadInitialConditionFromFile(P) ?  P.get<int>("Mesh.AdaptRef.LoadBalStrategy") :
                                                          -P.get<int>("Mesh.AdaptRef.LoadBalStrategy")));


    // If we read the Multigrid, it shouldn't be modified;
    // otherwise the pde-solutions from the ensight files might not fit.
    if (!P.get<int>("Transp.UseNSSol",1) || !ReadInitialConditionFromFile(P))
    {
        DROPS::InScaMap & tdscalarmap = DROPS::InScaMap::getInstance();
        DROPS::instat_scalar_fun_ptr distance = tdscalarmap[P.get("Transp.Levelset", std::string("Ellipsoid"))];

        DistMarkingStrategyCL marker( distance,
                                      P.get<double>("Mesh.AdaptRef.Width"),
                                      P.get<int>("Mesh.AdaptRef.CoarsestLevel"),
                                      P.get<int>("Mesh.AdaptRef.FinestLevel") );
        adap.set_marking_strategy( &marker );
        adap.MakeInitialTriang();
        adap.set_marking_strategy(0);
    }

    std::cout << DROPS::SanityMGOutCL(*mg) << std::endl;
#ifdef _PAR
    adap.GetLb().GetLB().SetWeightFnct(3);
    if (DROPS::ProcCL::Check( CheckParMultiGrid( adap.GetPMG())))
        std::cout << "As far as I can tell the ParMultigridCl is sane\n";
#endif

    if (P.get<int>("Transp.UseNSSol",1)){
        DROPS::InstatNavierStokes2PhaseP2P1CL prob( *mg, DROPS::TwoPhaseFlowCoeffCL(P), bnddata, P.get<double>("NavStokes.XFEMReduced")<0 ? DROPS::P1_FE : DROPS::P1X_FE, P.get<double>("NavStokes.XFEMReduced"));
        Strategy( prob, lsetbnddata, adap);    // do all the stuff
    }
    else if (!P.get<int>("SurfTransp.Enable",0))
    {
        OnlyTransportStrategy( *mg, lsetbnddata, adap);    // do just the transport stuff
    }
    else {
        DROPS::InstatNavierStokes2PhaseP2P1CL prob( *mg, DROPS::TwoPhaseFlowCoeffCL(0,0,0,0,0,DROPS::Point3DCL()), bnddata, DROPS::P1_FE);
        StatMassSurfTransportStrategy( *mg, prob, lsetbnddata, adap); // stationary mass/surfactant transport
    }


    delete mg;
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
