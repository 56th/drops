/// \file prJumpKirchhart.cpp
/// \brief test FE spaces for pressure
/// \author LNM RWTH Aachen: Matthias Kirchhart, Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen:

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

#include <cfloat>

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "stokes/instatstokes2phase.h"
#include "num/krylovsolver.h"
#include "num/precond.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "levelset/coupling.h"
#include "levelset/adaptriang.h"
#include "levelset/mzelle_hdr.h"
#include "levelset/surfacetension.h"
#include "num/stokessolverfactory.h"
#include "misc/funcmap.h"
#include <fstream>
#include <iomanip>
#include <vector>
#include "misc/dynamicload.h"

// Some abbrevations.
using namespace DROPS;
typedef InstatStokes2PhaseP2P1CL        StokesCL;
typedef StokesVelBndDataCL::bnd_val_fun bnd_vel_fun; 
typedef LsetBndDataCL::bnd_val_fun      bnd_lst_fun;



// The global parameter variable.
ParamCL P;

void read_params( int argc, char **argv );
void setup_problem( StokesCL*& Stokes, LevelsetP2CL*& lset, AdapTriangCL*& adap, MarkingStrategyCL*& marker,
                    VelVecDescCL*& curv, StokesSolverFactoryCL<StokesCL>*& factory );

void solve( StokesCL& Stokes, VelVecDescCL& curv,
            StokesSolverFactoryCL<StokesCL>& factory );

void output_matrices( StokesCL& Stokes );
void output( StokesCL& Stokes, LevelsetP2CL& lset );

double LBB_constant( const StokesCL &Stokes );

// Set p to the exact pressure solution when using the XFEM space.
void exact_pressure( VecDescCL& p, double delta_p, const MultiGridCL& mg,
                     const FiniteElementT prFE, const ExtIdxDescCL& Xidx );

// Compute the error of the numerical solution.
void L2ErrorPr( const VecDescCL& p, const MatrixCL& prM, double delta_p,
                const MultiGridCL& mg, const FiniteElementT prFE, const ExtIdxDescCL& Xidx );

// Create additional fields for output in Ensight files.
void PostProcessPr( const VecDescCL& p, VecDescCL& new_p, const MultiGridCL& mg );

double my_abs( double x );


int main (int argc, char** argv)
{
  try
  {
    read_params( argc, argv );

    DROPS::dynamicLoad( P.get<std::string>("General.DynamicLibsPrefix"),
                        P.get< std::vector<std::string> >("General.DynamicLibs") );

	StokesCL* stokes = 0;
	LevelsetP2CL* lset = 0;
	AdapTriangCL* adap = 0;
        MarkingStrategyCL* marker = 0;
	VelVecDescCL* curv = 0;
	StokesSolverFactoryCL<StokesCL>* factory = 0;

	setup_problem( stokes, lset, adap, marker, curv, factory );
	output_matrices( *stokes );
	solve( *stokes, *curv, *factory );
	output( *stokes, *lset );

	LBB_constant( *stokes );

	delete factory; factory = 0;
	delete curv; curv = 0;
	delete adap; adap = 0;
        delete marker;
	delete lset; lset = 0;
	delete stokes; stokes = 0;

    return 0;
  }
  catch ( DROPSErrCL err ) { err.handle(); }
}

/*!
 * \brief This function reads in the parameters from a file.
 *
 * When the program is started with a parameter, it is assumed to be the
 * filename where the parameters are stored. Otherwise the parameters are
 * read from the default location "prJumpKirchhart.json".
 */
void read_params( int argc, char **argv )
{
	std::string filename( "prJumpKirchhart.json" );
	if ( argc > 2 )
	{
		std::cerr << "Expecting at most one parameter, the input-filename."
		          << std::endl;
		std::abort();
	}
	else if ( argc == 2 )
	{
	    filename = argv[1];
	}

	std::cout << "Using \"" << filename << "\" as input-file." << std::endl;
	std::ifstream param( filename.c_str() );
	if ( ! param )
	{
		std::cerr << "Could not open the input-file." << std::endl;
		std::abort();
	}

	if ( ! ( param >> P ) )
	{
		std::cerr << "An error occured while reading the input file."
		          << std::endl;
		std::abort();
	}

	try
	{
		P.get<double>("Stokes.epsP");
	}
	catch ( ... )
	{
		P.put<double>( "Stokes.epsP", 0.0 );
	}

	try
	{
		P.get<bool>( "Stokes.PrintMatrices" );
	}
	catch ( ... )
	{
		P.put<bool>( "Stokes.PrintMatrices", false );
	}
}

void setup_problem( StokesCL*& Stokes, LevelsetP2CL*& lset, AdapTriangCL*& adap, MarkingStrategyCL*& marker,
                    VelVecDescCL*& curv, StokesSolverFactoryCL<StokesCL>*& factory )
{
	// Create the four corner nodes of the brick.
	const double L = 1; 
	Point3DCL orig(-L), e1, e2, e3;
	e1[0] = e2[1]= e3[2] = 2*L;

	// Create the mesh-creator.
	// n is the number of tetras in each dimension.
	const int n = P.get<int>("DomainCond.MeshFile");
	BrickBuilderCL builder( orig, e1, e2, e3, n, n, n );

	// Now specify the boundary conditions.
	const bnd_vel_fun ZeroVel = InVecMap::getInstance().find("ZeroVel")->second;
	const BndCondT bc[6] = { WallBC, WallBC, WallBC, WallBC, WallBC, WallBC };
	const bnd_vel_fun bnd_fun[6] = { ZeroVel, ZeroVel, ZeroVel, ZeroVel, ZeroVel, ZeroVel};

	Stokes = new StokesCL( builder, TwoPhaseFlowCoeffCL( P ),
                               StokesBndDataCL( 6, bc, bnd_fun ),
                               P.get<double>("Stokes.XFEMStab") < 0 ? P1_FE : P1X_FE,
                               P.get<double>("Stokes.XFEMStab") );

	// Now define the interface and refine the mesh around it.
	const BndCondT bcls[6] = { NoBC, NoBC, NoBC, NoBC, NoBC, NoBC };
	const bnd_lst_fun bfunls[6] = {0,0,0,0,0,0};
	LsetBndDataCL lsbnd( 6, bcls, bfunls );

	EllipsoidCL::Init( P.get<Point3DCL>("Exp.PosDrop"),
	                   P.get<Point3DCL>("Exp.RadDrop") );

        typedef DistMarkingStrategyCL MarkerT;
        marker = new MarkerT( EllipsoidCL::DistanceFct,
                              P.get<double>("AdaptRef.Width"),
                              P.get<int>("AdaptRef.CoarsestLevel"),
                              P.get<int>("AdaptRef.FinestLevel") );
	adap = new AdapTriangCL( Stokes->GetMG(), marker );
	adap->MakeInitialTriang();

	// Set up the surface tension functional settings.
	MultiGridCL& MG = Stokes->GetMG();
	sigma = P.get<double>("SurfTens.SurfTension");
	SurfaceTensionCL sf( sigmaf, 0 );
	lset = LevelsetP2CL::Create( MG, lsbnd, sf,
	                             P.get<double>("Levelset.SD"),
	                             P.get<double>("Levelset.CurvDiff") );
	lset->SetSurfaceForce( SF_ImprovedLB );


	// Create the numberings and set the indeces.
	MLIdxDescCL *const lidx = &(lset->idx);
	MLIdxDescCL *const vidx = &(Stokes->vel_idx);
	MLIdxDescCL *const pidx = &(Stokes->pr_idx);
	if ( StokesSolverFactoryHelperCL().VelMGUsed(P) )
	{
		Stokes->SetNumVelLvl( Stokes->GetMG().GetNumLevel() );
		lset->SetNumLvl( Stokes->GetMG().GetNumLevel() );
	}
	if ( StokesSolverFactoryHelperCL().PrMGUsed(P) )
	{
		Stokes->SetNumPrLvl( Stokes->GetMG().GetNumLevel() );
		lset->SetNumLvl( Stokes->GetMG().GetNumLevel() );
	}

	lset->CreateNumbering( MG.GetLastLevel(), lidx );
	lset->Phi.SetIdx( lidx );
	lset->Init( EllipsoidCL::DistanceFct );

	Stokes->CreateNumberingVel( MG.GetLastLevel(), vidx );
	Stokes->CreateNumberingPr ( MG.GetLastLevel(), pidx, NULL, lset );

	Stokes->SetIdx();
	Stokes->v.SetIdx(vidx);
	Stokes->p.SetIdx(pidx);
	std::cout << Stokes->p.Data.size() << " pressure unknowns,\n";
	std::cout << Stokes->v.Data.size() << " velocity unknowns,\n";
	std::cout << lset->Phi.Data.size() << " levelset unknowns.\n";


	// Set initial guess for velocity solutiom.
	Stokes->InitVel( &Stokes->v, InVecMap::getInstance().find("ZeroVel")->second );


	// Create prolongation matrices
	factory = new StokesSolverFactoryCL<StokesCL>( *Stokes, P );
	UpdateProlongationCL<Point3DCL> PVel( MG, factory->GetPVel(),
	                                      vidx, vidx );
	adap->push_back( &PVel );

	UpdateProlongationCL<double> PPr ( MG, factory->GetPPr(),
	                                   pidx, pidx );
	adap->push_back( &PPr );

	UpdateProlongationCL<double> PLset( MG, lset->GetProlongation(),
	                                    lidx, lidx ); 
	adap->push_back( &PLset );
	lset->UpdateMLPhi();

	factory->SetMatrixA(  &Stokes->A.Data.GetFinest() );
	factory->SetMatrices( &Stokes->A.Data, &Stokes->B.Data,
	                      &Stokes->M.Data, &Stokes->prM.Data, &Stokes->pr_idx );

	// Assemble the system matrices.
	TimerCL time;
	time.Reset();
	Stokes->SetupPrMass(  &Stokes->prM, *lset );
	Stokes->SetupPrStiff( &Stokes->prA, *lset );
	time.Stop();
	std::cout << "Setup of pressure mass and stiffnes matrices took "
	          << time.GetTime() << " seconds.\n";

	time.Reset();
	Stokes->SetupPrJ( &Stokes->prJ, *lset );
	Stokes->prJ.Data.GetFinest() *= -P.get<double>( "Stokes.epsP" );
	time.Stop();
	std::cout << "Setup of pressure stabilisation matrix took "
	          << time.GetTime() << " seconds.\n";
	std::cout.flush();

	curv = new VelVecDescCL( vidx );
	time.Reset();
	Stokes->SetupSystem1( &Stokes->A, &Stokes->M, &Stokes->b, &Stokes->b, curv, *lset, Stokes->v.t );
	Stokes->SetupSystem2( &Stokes->B, &Stokes->c, *lset, Stokes->v.t );
	time.Stop();
	std::cout << "Assembly of system matrices took: "
	          << time.GetTime() << " seconds.\n";

	// Compute surface tension functional.
	time.Reset();
	curv->Clear( Stokes->v.t );
	lset->AccumulateBndIntegral( *curv );
	time.Stop();
	std::cout << "Discretisation of surface tension functional took "
	          << time.GetTime() << " seconds.\n";
}

void solve( StokesCL& Stokes, VelVecDescCL& curv,
            StokesSolverFactoryCL<StokesCL>& factory )
{
	StokesSolverBaseCL* solver = factory.CreateStokesSolver();

	TimerCL time;
	time.Reset();
	solver->Solve( Stokes.A.Data, Stokes.B.Data, Stokes.prJ.Data,
	               Stokes.v.Data, Stokes.p.Data, curv.Data, Stokes.c.Data,
	               Stokes.v.RowIdx->GetEx(), Stokes.p.RowIdx->GetEx() );
	time.Stop();
	std::cout << "iter: " << solver->GetIter()
	          << "\tresid: " << solver->GetResid() << std::endl;
	std::cout << "Solving Stokes took "<<time.GetTime()<<" seconds.\n";

	delete solver;
}


void output_matrices( StokesCL& Stokes )
{
	if ( P.get<bool>("Stokes.PrintMatrices") == true )
	{
		std::ofstream Afile( "A.txt" );
		std::ofstream Bfile( "B.txt" );
		std::ofstream Cfile( "C.txt" );
		std::ofstream Mfile( "M.txt" );

		Afile << std::setprecision(16);
		Afile << Stokes.A.Data.GetFinest();
		Bfile << std::setprecision(16);
		Bfile << Stokes.B.Data.GetFinest();
		Cfile << std::setprecision(16);
		Cfile << Stokes.prJ.Data.GetFinest();
		Mfile << std::setprecision(16);
		Mfile << Stokes.prM.Data.GetFinest();
	}
}

void output( StokesCL& Stokes, LevelsetP2CL& lset )
{
	const double radius = P.get<DROPS::Point3DCL>("Exp.RadDrop")[0];
	const double prJump = sigma*2.0/radius;
	VecDescCL new_pr; 
	new_pr.SetIdx( &lset.idx );

	MultiGridCL& MG = Stokes.GetMG();
	const VectorCL& u = Stokes.v.Data;
	std::cout << "\n----------------\n || u ||_oo = " << supnorm(u)
	          << "\n || u ||_L2  = " << std::sqrt( dot( Stokes.M.Data*u, u ) )
	          << "\n || u ||_H1  = " << std::sqrt( dot( Stokes.M.Data*u, u ) + dot( Stokes.A.Data*u, u ) )
	          << "\n----------------\n";

	L2ErrorPr( Stokes.p, Stokes.prM.Data.GetFinest(), prJump, Stokes.GetMG(), Stokes.GetPrFE(), Stokes.GetXidx() );
	PostProcessPr( Stokes.p, new_pr, MG );

	// Initialize Ensight6 output
	std::string ensf( P.get<std::string>("Ensight.EnsDir") + "/" + P.get<std::string>("Ensight.EnsCase"));
	Ensight6OutCL ensight( P.get<std::string>("Ensight.EnsCase") + ".case", P.get<int>("Time.NumSteps") + 1, P.get<int>("Ensight.Binary") );
	ensight.Register( make_Ensight6Geom  ( MG, MG.GetLastLevel(), P.get<std::string>("Ensight.GeomName"),           ensf + ".geo"));
	ensight.Register( make_Ensight6Scalar( lset.GetSolution(),             "Levelset",  ensf + ".scl"));
	ensight.Register( make_Ensight6Scalar( Stokes.GetPrSolution( new_pr),  "Pressure",  ensf + ".pr"));
	ensight.Register( make_Ensight6Vector( Stokes.GetVelSolution(),        "Velocity",  ensf + ".vel"));

	if (Stokes.UsesXFEM())
	    ensight.Register( make_Ensight6P1XScalar( MG, lset.Phi, Stokes.p,  "XPressure", ensf + ".pr"));

    	if ( P.get<int>("Ensight.EnsightOut") )
		ensight.Write();

	std::cout << std::endl;
}

void exact_pressure( VecDescCL& p, double delta_p, const MultiGridCL& mg,
                     const FiniteElementT prFE, const ExtIdxDescCL& Xidx )
{
	const Uint lvl = p.RowIdx->TriangLevel();
	const Uint idxnum = p.RowIdx->GetIdx();

	delta_p /= 2;
	switch ( prFE )
	{
	case P0_FE:
		for ( MultiGridCL::const_TriangTetraIteratorCL it = mg.GetTriangTetraBegin(lvl),
		     end = mg.GetTriangTetraEnd(lvl); it != end; ++it)
		{
			const double dist = EllipsoidCL::DistanceFct( GetBaryCenter( *it), 0);
			p.Data[it->Unknowns(idxnum)]= dist > 0 ? -delta_p : delta_p;
		}
		break;
	case P1X_FE:
		for( MultiGridCL::const_TriangVertexIteratorCL it = mg.GetTriangVertexBegin(lvl),
		     end = mg.GetTriangVertexEnd(lvl); it != end; ++it)
		{
			const IdxT idx = it->Unknowns(idxnum);
			if (Xidx[idx]==NoIdx) continue;
			p.Data[Xidx[idx]]= -2*delta_p; // jump height
		}

	case P1_FE: // and P1X_FE
		for( MultiGridCL::const_TriangVertexIteratorCL it= mg.GetTriangVertexBegin(lvl),
		    end = mg.GetTriangVertexEnd(lvl); it != end; ++it)
		{
			const double dist= EllipsoidCL::DistanceFct( it->GetCoord(), 0.);
			p.Data[it->Unknowns(idxnum)]= InterfacePatchCL::Sign(dist)==1 ? -delta_p : delta_p;
		}
		break;
      default:
        std::cout << "exact_pressure not implemented yet for this FE type!\n";
    }
}

void L2ErrorPr( const VecDescCL& p,
                const MatrixCL& prM, double delta_p, const MultiGridCL& mg,
                const FiniteElementT prFE, const ExtIdxDescCL& Xidx )
{
    // Shift pressure values such that the average over the total volume is 0.
    VectorCL ones( 1.0, p.Data.size() );
    if (prFE==P1X_FE)
        for (int i=Xidx.GetNumUnknownsStdFE(), n = ones.size(); i<n; ++i)
            ones[i]= 0;

    const double Vol   = dot( prM*ones, ones );
    const double p_avg = dot( prM*p.Data, ones ) / Vol;
    VectorCL shifted_pressure(p.Data - p_avg*ones);

    VecDescCL p_exactunshifted( p.RowIdx );
    exact_pressure( p_exactunshifted, delta_p, mg, prFE, Xidx );
    const double p_exactunshifted_avg = dot( prM*p_exactunshifted.Data, ones ) / Vol;
    VectorCL p_exact( p_exactunshifted.Data - p_exactunshifted_avg*ones );
    
    VectorCL error = VectorCL( shifted_pressure - p_exact );

    const double L2 = std::sqrt( dot(prM*error, error) );
    std::cout << "||e_p||_oo = " << supnorm(error) << std::endl;
    std::cout << "||e_p||_L2 = " << L2 << std::endl; 
}

void PostProcessPr( const VecDescCL& p, VecDescCL& new_p, const MultiGridCL& mg )
// as Ensight output is for cont. data only, put values of P0 DoFs in tetras;;
// into P1 DoFs in vertices (s.t. the average over each tetra gives
// the former P0 value)
{
    VectorCL num( new_p.Data.size()), sum( new_p.Data.size());

    const Uint lvl= p.RowIdx->TriangLevel(),
        idxnum= p.RowIdx->GetIdx(),
        idxnum2= new_p.RowIdx->GetIdx();

    if (p.RowIdx->NumUnknownsTetra())
        for( MultiGridCL::const_TriangTetraIteratorCL it= mg.GetTriangTetraBegin(lvl),
            end= mg.GetTriangTetraEnd(lvl); it != end; ++it)
        {
            const double val= p.Data[it->Unknowns(idxnum)];
            for (int i=0; i<4; ++i)
            {
                const IdxT nr2= it->GetVertex(i)->Unknowns(idxnum2);
                sum[nr2]+= val;
                num[nr2]+= 1;
            }
        }

    for( MultiGridCL::const_TriangVertexIteratorCL it= mg.GetTriangVertexBegin(lvl),
        end= mg.GetTriangVertexEnd(lvl); it != end; ++it)
    {
        const IdxT nr= it->Unknowns(idxnum),
            nr2= it->Unknowns(idxnum2);
        double val= p.Data[nr];
        new_p.Data[nr2]= p.RowIdx->NumUnknownsTetra() ? val + sum[nr2]/num[nr2]
                                                      : val;
    }
}

double LBB_constant( const StokesCL &Stokes )
{
    std::cout << "Computing the value of the LBB-constant..." << std::endl;

    // Inverse power-iteration on M^-1 B A^-1 B^T
    const MatrixCL &A(Stokes.A.Data.GetFinest());
    const MatrixCL &B(Stokes.B.Data.GetFinest());
    const MatrixCL &C( Stokes.prJ.Data.GetFinest() );
    const MatrixCL &M(Stokes.prM.Data.GetFinest());
    const MatrixCL &Mvel(Stokes.prM.Data.GetFinest());
    const VectorCL  Dvec( M.GetDiag() );
    const VectorAsDiagMatrixCL  D( &Dvec );

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
    rel_error = DBL_MAX;

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
        //y_k -= dot(D*y_k,ones)/dot(D*ones,ones)*ones;
        y_k /= std::sqrt( dot(y_k,y_k) );

	VectorCL rhs = M*y_k;
	//VectorCL rhs = D*y_k;
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
    return std::sqrt(lambda);
}

inline
double my_abs( double x )
{
    return std::abs(x);
}

