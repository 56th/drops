/// \file prJump.cpp
/// \brief test FE spaces for pressure
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen:

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

#include <vector>
#include <memory>
#include <limits>

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "stokes/instatstokes2phase.h"
#include "num/krylovsolver.h"
#include "num/precond.h"
#include "out/output.h"
#include "out/vtkOut.h"
#include "levelset/coupling.h"
#include "levelset/adaptriang.h"
#include "levelset/marking_strategy.h"
#include "levelset/mzelle_hdr.h"
#include "levelset/surfacetension.h"
#include "num/stokessolverfactory.h"
#include "misc/funcmap.h"
#include "misc/dynamicload.h"

using namespace DROPS;

// Some useful abbreviations.
typedef StokesVelBndDataCL::bnd_val_fun VelBndFunT;
typedef LsetBndDataCL::bnd_val_fun     LsetBndFunT;
// This has already been declared somewhere in the headers.
// typedef InstatStokes2PhaseP2P1CL StokesT;
typedef DistMarkingStrategyCL MarkerT;

ParamCL P;

void Assemble( StokesT& Stokes, LevelsetP2CL& lset, VelVecDescCL& curv );
void Solve( StokesT& Stokes, LevelsetP2CL& lset, VelVecDescCL& curv );

// Functions for the computation of the pressure error.
double L2ErrorPr( VecDescCL& p, const MatrixCL& M, double delta_p,
                  const MultiGridCL& MG, const FiniteElementT prFE,
                  const ExtIdxDescCL& Xidx );
void ExactPr( VecDescCL& p, double delta_p, const MultiGridCL& mg,
              const FiniteElementT prFE, const ExtIdxDescCL& Xidx );
void NormalisePr( VectorCL& p, const MatrixCL& M,
                  const FiniteElementT prFE, const ExtIdxDescCL& Xidx );

// Computation of the LBB constant.
void LBB_constant( const StokesT &Stokes );

// VTK output.
void output( StokesT &Stokes, LevelsetP2CL& lset, BndDataCL<>& lsetbnd );

//double DistanceFct( const Point3DCL& p )
//{
//    // Plane perpendicular to n=PosDrop with distance Radius[0] from origin.
//    return inner_prod( C.Mitte/norm(C.Mitte), p) - C.Radius[0];
//}

int main ( int argc, char** argv )
{
    try
    {
        // Read input.
        read_parameter_file_from_cmdline( P, argc, argv, "prJump.json" );
        //std::cout << P << std::endl;

        dynamicLoad( P.get<std::string>("General.DynamicLibsPrefix"),
                     P.get< std::vector<std::string> >("General.DynamicLibs") );

        // Define Boundary conditions.
        const BndCondT bc[6] = { WallBC, WallBC, WallBC,
                                 WallBC, WallBC, WallBC };

        VelBndFunT ZeroVel = InVecMap::getInstance().find("ZeroVel")->second;
        const VelBndFunT bnd_fun[6] = { ZeroVel, ZeroVel, ZeroVel,
                                        ZeroVel, ZeroVel, ZeroVel };

        const BndCondT bcls[6] = { NoBC, NoBC, NoBC,
                                   NoBC, NoBC, NoBC };
        const LsetBndFunT bfunls[6] = { 0, 0, 0, 0, 0, 0 };
        LsetBndDataCL lsbnd( 6, bcls, bfunls );


        // Create Initial Grid builder.
        const double L = 1; // Vol = 8 L³
        Point3DCL orig(-L), e1, e2, e3;
        e1[0]= e2[1]= e3[2]= 2*L;

        const int n = P.get<int>("DomainCond.MeshFile");
        BrickBuilderCL builder( orig, e1, e2, e3, n, n, n );

        // Create Stokes Problem
        FiniteElementT prFE = P.get<double>("Stokes.XFEMStab") < 0 ? P1_FE : P1X_FE;
        StokesT prob( builder, TwoPhaseFlowCoeffCL(P),
                      StokesBndDataCL( 6, bc, bnd_fun ), prFE,
                      P.get<double>("Stokes.XFEMStab") );

        MultiGridCL& MG = prob.GetMG();
   
        // Refine the mesh around the interface.
        EllipsoidCL::Init( P.get<Point3DCL>("Exp.PosDrop"),
                           P.get<Point3DCL>("Exp.RadDrop") );

        MarkerT marker( EllipsoidCL::DistanceFct,
                        P.get<double>("AdaptRef.Width"),
                        P.get<int>("AdaptRef.CoarsestLevel"),
                        P.get<int>("AdaptRef.FinestLevel") ); 
        AdapTriangCL adap( MG, &marker );
        adap.MakeInitialTriang();

        // Create and initialise the levelset function.
        sigma = P.get<double>( "SurfTens.SurfTension" );
        SurfaceTensionCL sf( sigmaf, 0 );

        std::auto_ptr<LevelsetP2CL> lset_ptr
        (
            LevelsetP2CL::Create( MG, lsbnd, sf,
                                  P.get<double>("Levelset.SD"),
                                  P.get<double>("Levelset.CurvDiff") )
        );
        LevelsetP2CL& lset = *lset_ptr;
        lset.SetNumLvl( MG.GetNumLevel() );
        lset.CreateNumbering( MG.GetLastLevel(), &lset.idx );
        lset.Phi.SetIdx( &lset.idx );
        lset.Init( EllipsoidCL::DistanceFct );

        std::cout << SanityMGOutCL(MG) << std::endl;

        VelVecDescCL rhs;
        Assemble( prob, lset, rhs );

        LBB_constant( prob );

        Solve( prob, lset, rhs );

        // Create an unscaled version of the stokes version.
        // Its mass and stiffness matrices allow us to easily compute the
        // L2- and H1- errors of the solution.
        P.put<double>( "Mat.ViscDrop", 1 );
        P.put<double>( "Mat.DensDrop", 1 );
        P.put<double>( "Mat.ViscFluid", 1 );
        P.put<double>( "Mat.DensFluid", 1 );
        StokesT unscal_prob( MG, TwoPhaseFlowCoeffCL(P),
                             StokesBndDataCL( 6, bc, bnd_fun ), prFE,
                             P.get<double>("Stokes.XFEMStab") );
        Assemble( unscal_prob, lset, rhs );

        const double radius = P.get<DROPS::Point3DCL>("Exp.RadDrop")[0];
        const double delta_p = sigma*2.0/radius;
        double unscal_error = L2ErrorPr( prob.p, unscal_prob.prM.Data.GetFinest(), delta_p,
                                         MG, prFE, prob.GetXidx() );
        std::cout << "||e_p||_L2 = " << unscal_error << std::endl;

        double l2_vel_error = dot( unscal_prob.M.Data.GetFinest()*prob.v.Data, prob.v.Data );
        double h1_vel_error = dot( unscal_prob.A.Data.GetFinest()*prob.v.Data, prob.v.Data );
               h1_vel_error += l2_vel_error;
        l2_vel_error = std::sqrt( l2_vel_error );
        h1_vel_error = std::sqrt( h1_vel_error );
        std::cout << "||e_v||_L2 = " << l2_vel_error << std::endl;
        std::cout << "||e_v||_H1 = " << h1_vel_error << std::endl;
 
        output( prob, lset, lsbnd );
    }
    catch ( DROPSErrCL err )
    {
        err.handle();
    }

    return 0;
}

void Assemble( StokesT& Stokes, LevelsetP2CL& lset, VelVecDescCL& curv )
{
    // Initialise the pressure, velocity and level-set vectors. 
    MultiGridCL& MG = Stokes.GetMG();
    MLIdxDescCL* vidx = &Stokes.vel_idx;
    MLIdxDescCL* pidx = &Stokes.pr_idx;
    if ( StokesSolverFactoryHelperCL().VelMGUsed(P) )
    {
        Stokes.SetNumVelLvl( MG.GetNumLevel() );
        lset.SetNumLvl( MG.GetNumLevel() );
    }

    if ( StokesSolverFactoryHelperCL().PrMGUsed(P) )
    {
        Stokes.SetNumPrLvl( MG.GetNumLevel() );
        lset.SetNumLvl( MG.GetNumLevel() );
    }

    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx );
    Stokes.CreateNumberingPr( MG.GetLastLevel(), pidx, NULL, &lset );

    Stokes.SetIdx();
    Stokes.v.SetIdx( vidx );
    Stokes.p.SetIdx( pidx );
    Stokes.InitVel( &Stokes.v, InVecMap::getInstance().find("ZeroVel")->second );

    // Initialise prolongation matrices
    StokesSolverFactoryCL<InstatStokes2PhaseP2P1CL> factory( Stokes, P );
    UpdateProlongationCL<Point3DCL> PVel( MG, factory.GetPVel(), vidx, vidx );
    UpdateProlongationCL<double> PPr ( MG, factory.GetPPr(), pidx, pidx );
    UpdateProlongationCL<double> PLset( MG, lset.GetProlongation(),
                                        lset.idxC, lset.idxC );
    lset.UpdateMLPhi();

    // Output system information.
    MG.SizeInfo( std::cout );
    std::cout << Stokes.p.Data.size()  << " pressure unknowns,\n";
    std::cout << Stokes.v.Data.size()  << " velocity unknowns,\n";
    std::cout << lset.Phi.Data.size()  << " levelset unknowns.\n";


    // Assemble mass and stiffness matrices for the pressure space.
    TimerCL time;
    time.Reset();
    Stokes.SetupPrMass ( &Stokes.prM, lset );
    Stokes.SetupPrStiff( &Stokes.prA, lset ); // makes no sense for P0
    time.Stop();
    std::cout << "Assembling the pressure space's mass- and stiffness matrices "
              << "took: " << time.GetTime() << " seconds." << std::endl;


    // Assemble the stabilisation matrix
    time.Reset();
    Stokes.SetupC( &Stokes.C, lset, P.get<double>("Stokes.epsP") );
    time.Stop();
    std::cout << "Assembling the stabilisation matrix took: " << time.GetTime()
              << " seconds." << std::endl;

    // Assemble system matrices.
    curv = VelVecDescCL( vidx );

    time.Reset();
    Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b,
                         &curv, lset, Stokes.v.t );
    Stokes.SetupSystem2( &Stokes.B, &Stokes.c, lset, Stokes.v.t );

    time.Stop();
    std::cout << "Assembling the system matrices took: " << time.GetTime()
              << " seconds." << std::endl;


    // Compute right hand side.
    time.Reset();
    curv.Clear( Stokes.v.t );

    lset.SetSurfaceForce( SF_ImprovedLB );
//  lset.SetSurfaceForce( SF_LB );
//  lset.SetSurfaceForce( SF_Const );
    lset.AccumulateBndIntegral( curv );

    time.Stop();
    std::cout << "Assembling the right hand side took: " << time.GetTime()
              << " seconds." << std::endl;
}

void Solve( StokesT& Stokes, LevelsetP2CL& lset, VelVecDescCL& curv )
{
    const MultiGridCL& MG = Stokes.GetMG();
    MLIdxDescCL* vidx = &Stokes.vel_idx;
    MLIdxDescCL* pidx = &Stokes.pr_idx;

    // Initialise prolongation matrices
    StokesSolverFactoryCL<InstatStokes2PhaseP2P1CL> factory( Stokes, P );
    UpdateProlongationCL<Point3DCL> PVel( MG, factory.GetPVel(), vidx, vidx );
    UpdateProlongationCL<double> PPr ( MG, factory.GetPPr(), pidx, pidx );
    UpdateProlongationCL<double> PLset( MG, lset.GetProlongation(),
                                        lset.idxC, lset.idxC );

    // Create a solver...
    factory.SetMatrixA ( &Stokes.A.Data.GetFinest() );
    factory.SetMatrices( &Stokes.A.Data, &Stokes.B.Data,
                         &Stokes.M.Data, &Stokes.prM.Data, pidx );
    std::auto_ptr<StokesSolverBaseCL> solver( factory.CreateStokesSolver() );

    // ...and solve.
    TimerCL time;
    time.Reset();

    solver->Solve( Stokes.A.Data, Stokes.B.Data, Stokes.C.Data,
                   Stokes.v.Data, Stokes.p.Data,
                   curv.Data, Stokes.c.Data,
                   Stokes.v.RowIdx->GetEx(), Stokes.p.RowIdx->GetEx() );
    time.Stop();

    std::cout << "iter: "  << solver->GetIter()  << '\t'
              << "resid: " << solver->GetResid() << std::endl;
    std::cout << "Solving the system took: " << time.GetTime()
              << " seconds." << std::endl;

}

double L2ErrorPr( VecDescCL& p, const MatrixCL& M, double delta_p,
                  const MultiGridCL& MG, const FiniteElementT prFE,
                  const ExtIdxDescCL& Xidx )
    //This if clause filters out too large/small pressure values due to xfem
{
    NormalisePr( p.Data, M, prFE, Xidx );

    VecDescCL p_exact( p.RowIdx );
    ExactPr( p_exact, delta_p, MG, prFE, Xidx );
    NormalisePr( p_exact.Data, M, prFE, Xidx );

    VectorCL error = p_exact.Data;
    error -= p.Data;
    return std::sqrt( dot( M*error, error ) );
}

void ExactPr( VecDescCL& p, double delta_p, const MultiGridCL& mg,
              const FiniteElementT prFE, const ExtIdxDescCL& Xidx )
{
    const Uint lvl    = p.RowIdx->TriangLevel();
    const Uint idxnum = p.RowIdx->GetIdx();

    typedef MultiGridCL::const_TriangTetraIteratorCL TetraIterT;
    const TetraIterT tbegin = mg.GetTriangTetraBegin( lvl );
    const TetraIterT tend   = mg.GetTriangTetraEnd( lvl );

    typedef MultiGridCL::const_TriangVertexIteratorCL VertexIterT;
    const VertexIterT vbegin = mg.GetTriangVertexBegin( lvl );
    const VertexIterT vend   = mg.GetTriangVertexEnd( lvl );

    delta_p /= 2;
    switch ( prFE )
    {
    case P0_FE:
        for( TetraIterT it = tbegin; it != tend; ++it )
        {
            const double dist = EllipsoidCL::DistanceFct( GetBaryCenter(*it), 0 );
            p.Data[ it->Unknowns(idxnum) ]= dist > 0 ? -delta_p : delta_p;
        }
        break;

    case P1X_FE:
        for( VertexIterT it = vbegin; it != vend; ++it )
        {
            const IdxT idx = it->Unknowns(idxnum);
            if ( Xidx[ idx ] != NoIdx )
            {
                p.Data[ Xidx[ idx ] ]= -2*delta_p; // jump height
            }
        }
        // No break!

    case P1_FE: // and P1X_FE
        for( VertexIterT it = vbegin; it != vend; ++it )
        {
            const double dist = EllipsoidCL::DistanceFct( it->GetCoord(), 0. );
            p.Data[ it->Unknowns(idxnum) ] = InterfacePatchCL::Sign( dist ) == 1 ? -delta_p : delta_p;
        }
        break;
    P.put_if_unset<int>("Levelset.Discontinuous", 0);

    default:
        throw DROPSErrCL( "ExactPr not implemented for this FE type!" );
    }
	
    P.put_if_unset<std::string>("Exp.Solution_Vel", "None");
    P.put_if_unset<std::string>("Exp.Solution_GradPr", "None");
    P.put_if_unset<std::string>("Exp.Solution_Pr", "None");
}

void NormalisePr( VectorCL& p, const MatrixCL& M,
                  const FiniteElementT prFE, const ExtIdxDescCL& Xidx )
{
    VectorCL ones( 1.0, p.size() );
    if ( prFE == P1X_FE )
    {
        for ( Uint i = Xidx.GetNumUnknownsStdFE(); i < ones.size(); ++i )
        {
            ones[i]= 0;
        }
    }
    const double p_avg = dot( M*p, ones ) / dot( M*ones, ones );
    p -= p_avg*ones;
}

void LBB_constant( const StokesT &Stokes )
{
    std::cout << "The computation of the LBB constant takes a long time.\n"
              << "Do you wish to proceed? (y/n)" << std::endl;

    char answer;
    std::cin >> answer;
    if ( answer != 'y' )
    {
        std::cout << "Skipping the computation of the LBB constant." << std::endl;
        return;
    }
    std::cout << "Computing the value of the LBB-constant..." << std::endl;

    // Inverse power-iteration on M^-1 B A^-1 B^T
    const MatrixCL &A(Stokes.A.Data.GetFinest());
    const MatrixCL &B(Stokes.B.Data.GetFinest());
    const MatrixCL &C(Stokes.C.Data.GetFinest());
    const MatrixCL &M(Stokes.prM.Data.GetFinest());
    const MatrixCL &Mvel(Stokes.M.Data.GetFinest());

    // Create a preconditioner for A. A PCG
    typedef SSORPcCL SymmPcPcT;
    typedef PCGSolverCL<SymmPcPcT> PCGSolverT;
    typedef SolverAsPreCL<PCGSolverT> APcT;

    SymmPcPcT symmpcpc;
    PCGSolverT Asolver( symmpcpc, 500, 0.01, true );
    APcT Apc(Asolver);

    // Approximate Schur.
    typedef ApproximateSchurComplMatrixCL<APcT, MatrixCL, DummyExchangeCL > SchurCL;
    SchurCL schur( A, Apc, B, C, DummyExchangeCL() );
    

    // Schur-Inverter...
    int tmp = 5000;
    typedef GCRSolverCL<ISBBT_Stab_PreCL> SolverT;
    ISBBT_Stab_PreCL Sprecond( &B, &C, &M, &Mvel, Stokes.pr_idx.GetFinest() );
    SolverT gcr( Sprecond, tmp, 5000, 1e-10, true );
 
    VectorCL y_k( Stokes.p.Data.size() ), y_k1( Stokes.p.Data.size() );
    double lambda_prev = 2, lambda_prev2 = 3, lambda = 1;
    double q, rel_error;
    rel_error = std::numeric_limits<double>::max();

    const ExtIdxDescCL& Xidx = Stokes.GetXidx();
    VectorCL ones( 1.0, Stokes.p.Data.size() );
    if ( Stokes.GetPrFE() == P1X_FE )
        for ( size_t i = Xidx.GetNumUnknownsStdFE(); i < ones.size(); ++i )
            ones[i]= 0;

    y_k1[0] = 1;

    int i = 0;
    while ( true )
    {
    if ( i++ > 50 && rel_error < 1e-3 )
        break;

    y_k = y_k1;
    y_k -= dot(M*y_k,ones)/dot(M*ones,ones)*ones;
    y_k /= std::sqrt( dot(y_k,y_k) );

    VectorCL rhs = M*y_k;
    y_k1 = std::abs((1.0/lambda))*y_k;
    gcr.Solve( schur, y_k1, rhs, DummyExchangeCL() );
    std::cout << "GCR-Residual: " << gcr.GetResid()
              << ".\t GCR-Iterations: " << gcr.GetIter() << ".\n";

    lambda_prev2 = lambda_prev;
    lambda_prev = lambda;
    lambda = 1.0/dot( y_k, y_k1 );

    q = ( lambda - lambda_prev ) / ( lambda_prev - lambda_prev2 );
    double error  = std::abs((q/(1.0-q))*(lambda-lambda_prev));
    rel_error = error/std::abs(lambda);
    std::cout << "Eigenvalue estimate: " << lambda << "\tError estimate: "
              << rel_error * 100 << "%\tLBB-Estimate: " << std::sqrt(lambda) << '\n';
    std::cout.flush();
    }      
    std::cout << "The value of the LBB-constant is: " << std::sqrt(lambda)
              << std::endl;
}

void output( StokesT &Stokes, LevelsetP2CL& lset, BndDataCL<>& lsetbnd )
{
    using std::string;

    if ( ! P.get<bool>("VTK.vtkOut") ) return;

    MultiGridCL& MG = Stokes.GetMG();
    string casename = P.get<string>("VTK.vtkCase");
    string casedir  = P.get<string>("VTK.vtkDir");
    bool   binary   = P.get<bool>("VTK.Binary");
    
    VTKOutCL vtk( MG, casename, 1, casedir, casename, casename, binary );

    vtk.Register( make_VTKScalar( Stokes.GetPrSolution(), "pressure" ) );

    FiniteElementT prFE = P.get<double>("Stokes.XFEMStab") < 0 ? P1_FE : P1X_FE;
    if ( prFE == P1X_FE )
    {
        vtk.Register( make_VTKP1XScalar( MG, lset.Phi, Stokes.p, lsetbnd, "Xpressure" ) );
    }
    vtk.Register( make_VTKScalar( make_P2Eval( MG, lsetbnd, lset.Phi ), "levelset" ) );

    vtk.Register( make_VTKVector( Stokes.GetVelSolution() , "velocity" ) );
    vtk.Write( 0 );
}

