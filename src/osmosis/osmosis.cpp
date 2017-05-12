/// \file osmosis.cpp
/// \brief osmosis model
///        basically a copy of twophasedrops.cpp with extensions for osmosis
/// 	   Problem:
/// 		du/dt - mu*laplace u = 0    on Omega
///     	    u*V_n + mu*du/dn = 0    on dOmega
///			alpha*kappa + beta*u = V_n  on dOmega
/// \author LNM RWTH Aachen: Christoph Lehrenfeld, Thorolf Schulte

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
 * Copyright 2012 LNM/SC RWTH Aachen, Germany
*/

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "navstokes/instatnavstokes2phase.h"
#include "stokes/integrTime.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "out/vtkOut.h"
#include "levelset/coupling.h"
#include "misc/params.h"
#include "levelset/adaptriang.h"
#include "levelset/marking_strategy.h"
#include "levelset/mzelle_hdr.h"
#include "levelset/surfacetension.h"
#include "osmosis/osmosisSetup.h"
#include "surfactant/ifacetransp.h"
#include "levelset/twophaseutils.h"
#include "misc/funcmap.h"
#include "misc/scopetimer.h"


#include "num/stokessolverfactory.h"
#ifndef _PAR
#include "num/oseensolver.h"
#else
#include "num/parstokessolver.h"
#include "parallel/loadbal.h"
#include "parallel/parmultigrid.h"
#endif
#include <fstream>
#include <sstream>

#include "misc/progressaccu.h"
#include "misc/dynamicload.h"


DROPS::ParamCL P;

// Problem:
// du/dt - mu*laplace u = 0 on Omega
//     u*V_n + mu*du/dn = 0 on dOmega
//                  V_n = alpha*kappa + beta*u on dOmega

typedef DROPS::BndDataCL<> cBndDataCL;
typedef cBndDataCL::bnd_val_fun  c_bnd_val_fun;


namespace DROPS // for Strategy
{

void  OnlyOsmosisStrategy( MultiGridCL& MG, LsetBndDataCL& lsetbnddata, AdapTriangCL& adap)
{
    InVecMap & tdvectormap = InVecMap::getInstance();
    InScaMap & tdscalarmap = InScaMap::getInstance();
    InScaMap & scalarmap = InScaMap::getInstance();
    instat_vector_fun_ptr Flowfield = tdvectormap[P.get("Transp.Flow", std::string("ZeroVel"))];
    instat_scalar_fun_ptr Reaction = tdscalarmap["ReactionFct"];
    instat_scalar_fun_ptr Rhs = tdscalarmap["Rhs"];
    instat_scalar_fun_ptr Initialcneg = tdscalarmap[P.get("Osmosis.InitialConc", std::string("One"))];
    instat_scalar_fun_ptr Zero = tdscalarmap["ZeroFct"];
    instat_vector_fun_ptr ZeroVel = tdvectormap["ZeroFct"];
    instat_scalar_fun_ptr distance = scalarmap[P.get("Osmosis.Levelset", std::string("Ellipsoid"))];

    const c_bnd_val_fun c_bfun[6]= {Zero,Zero,Zero,Zero,Zero,Zero};
    const instat_vector_fun_ptr v_bfun[6]= {ZeroVel,ZeroVel,ZeroVel,ZeroVel,ZeroVel,ZeroVel};
    const c_bnd_val_fun c_bfunt[6]= {Zero,Zero,Zero,Zero,Zero,Zero};

    const DROPS::BndCondT c_bc[6]= {
        DROPS::OutflowBC, DROPS::OutflowBC, DROPS::OutflowBC,
        DROPS::OutflowBC, DROPS::OutflowBC, DROPS::OutflowBC
    };

    const DROPS::BndCondT v_bc[6]= {
        DROPS::DirBC, DROPS::DirBC, DROPS::DirBC,
        DROPS::DirBC, DROPS::DirBC, DROPS::DirBC
    };

    cBndDataCL Bnd_c( 6, c_bc, c_bfun);
    cBndDataCL Bnd_ct( 6, c_bc, c_bfunt);
    bool vdirvals = P.get("ZeroVelAtBnd.Active", 0);
    BndDataCL<Point3DCL> VelBnd(6, vdirvals? v_bc : c_bc , v_bfun); // dir. boundary conditions for extendedvel.?
    DROPS::instat_scalar_fun_ptr sigmap = 0;
    SurfaceTensionCL sf( sigmap, Bnd_c);

    LevelsetP2CL & lset( * LevelsetP2CL::Create( MG, lsetbnddata, sf, P.get_child("Levelset")) );
    // levelset wrt the previous time step:
    LevelsetP2CL & oldlset( * LevelsetP2CL::Create( MG, lsetbnddata, sf, P.get_child("Levelset")) );
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

	VelocityContainer vel(Flowfield);

    OsmosisP1CL osmosis( MG, Bnd_c, VelBnd, /*TODO: remove*/  vel, lsetbnddata, lset, oldlset, P, 0, Reaction, Rhs);
    OsmosisRepairCL transprepair(osmosis, MG.GetLastLevel());

    // index of the concentration wrt the interface at actual time step:
    MLIdxDescCL* cidx= &osmosis.idx;
    // index of the global velocity field at actual time step:
    IdxDescCL* velidx= &osmosis.Velidx;
    IdxDescCL* oldvelidx= &osmosis.oldVelidx;

    // index of the concentration wrt the interface at previous time step:
    MLIdxDescCL* cidx_old= &osmosis.oldidx;

    {
        adap.push_back(&transprepair);
        osmosis.CreateNumbering( MG.GetLastLevel(), cidx, cidx_old, velidx, oldvelidx, lset.Phi, oldlset.Phi);
        osmosis.conc.SetIdx( cidx);
        osmosis.Vn_.SetIdx( velidx);
        osmosis.oldconc.SetIdx( cidx_old);
        osmosis.oldVn_.SetIdx( oldvelidx);
        std::cout << osmosis.conc.Data.size() << " concentration unknowns,\n";
        std::cout << osmosis.Vn_.Data.size() << " velocity unknowns,\n";

        if (P.get<int>("NavStokes.InitialValue") != -1)
        {
            if (P.get<int>("Osmosis.ScaleInitialConc") == 1)
                osmosis.InitWithScaling( Initialcneg, P.get<double>("Osmosis.TotalConcentration"), 0);
            else
                osmosis.Init( Initialcneg, 0);
        }
        else
          ReadFEFromFile( osmosis.conc, MG, P.get<std::string>("Restart.InitialData")+"concentrationTransf");

	}

    //LevelsetModifyCL lsetmod( P.get<int>("Reparam.Freq"), P.get<int>("Reparam.Method"), P.get<double>("Reparam.MaxGrad"), P.get<double>("Reparam.MinGrad"));

    // GNUPLOT
    GNUPlotCL gnu(P.get("GNUPlot.Plotname", (std::string)"").c_str(), &osmosis, &OsmosisP1CL::SolutionErrorCase1);
	if (P.get("GNUPlot.Out", 0))
		gnu.Write(0.0);

    // writer for vtk-format
    VTKOutCL vtkwriter(adap.GetMG(), "DROPS data",
                       (P.get("VTK.Freq", 0) ?
                        P.get<int>("Time.NumSteps")/P.get("VTK.Freq", 0)+1 : 0),
                       P.get<std::string>("VTK.VTKDir"), P.get<std::string>("VTK.VTKName"),
                       P.get<std::string>("VTK.TimeFileName"),
                       P.get<int>("VTK.Binary"),
                       P.get<int>("VTK.UseOnlyP1"),
                       false, /* -< P2DG */
                       -1,  /* <- level */
                       P.get<int>("VTK.ReUseTimeFile") );

    vtkwriter.Register( make_VTKScalar( lset.GetSolution(), "level-set") );

    vtkwriter.Register( make_VTKScalar( osmosis.GetSolution( osmosis.conc), "Concentration") );
    vtkwriter.Register( make_VTKVector( osmosis.GetVelP1Solution( osmosis.Vn_), "V") );


    if (P.get("VTK.Freq", 0))
        vtkwriter.Write(0);


    const double dt = P.get<int>("Time.NumSteps")!=0 ? P.get<double>("Time.FinalTime")/P.get<int>("Time.NumSteps") : 0;

    // Create the marking strategy for the adaptive mesh refinement.
    typedef DistMarkingStrategyCL MarkerT;
    MarkerT marker( lset, P.get<double>("Mesh.AdaptRef.Width"),
                          P.get<int>("Mesh.AdaptRef.CoarsestLevel"),
                          P.get<int>("Mesh.AdaptRef.FinestLevel") );
    adap.set_marking_strategy( &marker );

    for (int step= 1; step<=P.get<int>("Time.NumSteps"); ++step)
    {
        std::cout << "============================================================ step " << std::setw(8) << step << "  /  " << std::setw(8) << P.get<int>("Time.NumSteps") << std::endl;
        double t= dt * step;

        // grid modification
        bool doGridMod= P.get<int>("Mesh.AdaptRef.Freq") && step%P.get<int>("Mesh.AdaptRef.Freq") == 0;
        if (doGridMod) {
            adap.UpdateTriang();
        }

		std::cout << osmosis.conc.Data.size() << " concentration unknowns\n";
		std::cout << osmosis.Vn_.Data.size() << " interface unknowns\n";

        osmosis.DoStep( t);

        static bool gnuoutnow = P.get("GNUPlot.Out", 0);
        if (gnuoutnow)
			gnu.Write(t);

        bool vtkoutnow = P.get("VTK.Freq", 0) && (step%P.get("VTK.Freq", 0)==0);// || step < 20);
        if (vtkoutnow){
            ScopeTimer vtktime("VTK-output");
            vtkwriter.Write(t);
        }
    }

    adap.set_marking_strategy(0);

    gnu.Close();

    std::cout << std::endl;
    delete &lset;
    delete &oldlset;
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

    DROPS::read_parameter_file_from_cmdline( P, argc, argv, "../../param/osmosis/osmosis/osmosis.json");
    P.put_if_unset<std::string>("VTK.TimeFileName",P.get<std::string>("VTK.VTKName"));
    std::cout << P << std::endl;

    DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"),
                       P.get<std::vector<std::string> >("General.DynamicLibs") );

    if (P.get<int>("General.ProgressBar"))
        DROPS::ProgressBarTetraAccumulatorCL::Activate();

    DROPS::MultiGridCL* mg= 0;
    DROPS::LsetBndDataCL lsetbnddata;
    std::auto_ptr<DROPS::MGBuilderCL> builder( DROPS::make_MGBuilder( P));
    mg = new DROPS::MultiGridCL( *builder);
    read_BndData( lsetbnddata,*mg, P.get_child( "Levelset.BoundaryData"));

    std::cout << "Generated MG of " << mg->GetLastLevel() << " levels." << std::endl;

    DROPS::EllipsoidCL::Init(P.get<DROPS::Point3DCL>("Levelset.PosDrop"), P.get<DROPS::Point3DCL>("Levelset.RadDrop"));
    DROPS::AdapTriangCL adap( *mg );
    // If we read the Multigrid, it shouldn't be modified;
    // otherwise the pde-solutions from the ensight files might not fit.
    if (!P.get<int>("Transp.UseNSSol") || (P.get<std::string>("Restart.InputData","") == "")){
        DROPS::InScaMap & scalarmap = DROPS::InScaMap::getInstance();
        DROPS::instat_scalar_fun_ptr distance = scalarmap[P.get("Osmosis.Levelset", std::string("Ellipsoid"))];

	typedef DROPS::DistMarkingStrategyCL InitMarkerT;
	InitMarkerT initmarker( distance, P.get<double>("Mesh.AdaptRef.Width"),
                                          P.get<int>( "Mesh.AdaptRef.CoarsestLevel" ),
                                          P.get<int>( "Mesh.AdaptRef.FinestLevel" ) );
        adap.set_marking_strategy( &initmarker );
        adap.MakeInitialTriang();
        adap.set_marking_strategy(0);
    }

    std::cout << DROPS::SanityMGOutCL(*mg) << std::endl;
#ifdef _PAR
    adap.GetLb().GetLB().SetWeightFnct(3);
    if (DROPS::ProcCL::Check( CheckParMultiGrid( adap.GetPMG())))
        std::cout << "As far as I can tell the ParMultigridCl is sane\n";
#endif


    // Osmosis without Navier Stokes etc.
    OnlyOsmosisStrategy( *mg, lsetbnddata, adap);


    delete mg;
    return 0;
  }

  catch (DROPS::DROPSErrCL& err) { err.handle(); }
  return -1;
}
