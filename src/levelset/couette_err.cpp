/// \file couette_err.cpp
/// \brief This file is modified from twophasedrops.cpp. 
///        Read solution (velocity, pressure, levelset) in reference mesh from files through serialization code. 
///        The solutions in coarse meshes are also read from files, which are prolongated with prolongation matrix. 
///        The error (discrepancy) between reference mesh and coarse meshes are computed. 
/// \author LNM RWTH Aachen: Liang Zhang

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
 * Copyright 2015 LNM/SC RWTH Aachen, Germany
*/

//multigrid
#include "geom/multigrid.h"
#include "geom/builder.h"
//time integration
#include "navstokes/instatnavstokes2phase.h"
#include "stokes/integrTime.h"
//output
#include "out/output.h"
#ifndef _PAR
#include "out/ensightOut.h"
#endif
#include "out/vtkOut.h"
//levelset
#include "levelset/coupling.h"
#include "levelset/marking_strategy.h"
#include "levelset/adaptriang.h"
#include "levelset/mzelle_hdr.h"
#include "levelset/twophaseutils.h"
//surfactants
#include "surfactant/ifacetransp.h"
//function map
#include "misc/funcmap.h"
//solver factory for stokes
#include "num/stokessolverfactory.h"
#include "num/oseensolver.h"
#include "num/prolongation.h"
#include "num/stokespardiso.h" 
#ifdef _PAR
#include "parallel/loadbal.h"
#include "parallel/parmultigrid.h"
#endif
//general: streams
#include <fstream>
#include <sstream>

#include "misc/progressaccu.h"
#include "misc/dynamicload.h"

#include <sys/resource.h>

DROPS::ParamCL P;

// rho*du/dt - mu*laplace u + Dp = f + rho*g - okn
//                        -div u = 0
//                             u = u0, t=t0

namespace DROPS // for Strategy
{
double GetTimeOffset(){
    double timeoffset = 0.0;
    const std::string restartfilename = P.get<std::string>("DomainCond.InitialFile");
    if (P.get<int>("DomainCond.InitialCond") == -1){
        const std::string timefilename = restartfilename + "time";
        std::ifstream f_(timefilename.c_str());
        f_ >> timeoffset;
        std::cout << "used time offset file is " << timefilename << std::endl;
        std::cout << "time offset is " << timeoffset << std::endl;
    }
    return timeoffset;
}


void Strategy( InstatNavierStokes2PhaseP2P1CL& Stokes, LsetBndDataCL& lsetbnddata, AdapTriangCL& adap)
// flow control
{
DROPS::InScaMap & inscamap = DROPS::InScaMap::getInstance();
    DROPS::InVecMap & invecmap = DROPS::InVecMap::getInstance();
    DROPS::MatchMap & matchmap = DROPS::MatchMap::getInstance();

    instat_scalar_fun_ptr the_Young_angle;
    instat_vector_fun_ptr the_Bnd_outnormal;

    the_Young_angle= inscamap[P.get<std::string>("SpeBnd.CtAngle")];
    the_Bnd_outnormal= invecmap[P.get<std::string>("SpeBnd.BndOutNormal")];

    bool is_periodic = P.get<std::string>("DomainCond.PeriodicMatching", "none") != "none";
    match_fun periodic_match = is_periodic ? matchmap[P.get("DomainCond.PeriodicMatching", std::string("periodicx"))] : 0;

    MultiGridCL& MG= Stokes.GetMG();

    // initialization of surface tension
    sigma= Stokes.GetCoeff().SurfTens;
    eps= P.get<double>("SurfTens.JumpWidth");    lambda= P.get<double>("SurfTens.RelPos");    sigma_dirt_fac= P.get<double>("SurfTens.DirtFactor");
    instat_scalar_fun_ptr sigmap  = 0;
    if (P.get<double>("SurfTens.VarTension"))
    {
        sigmap  = &sigma_step;
    }
    else
    {
        sigmap  = &sigmaf;
    }
    SurfaceTensionCL sf( sigmap);

    // Creates new Levelset-Object, has to be cleaned manually
    LevelsetP2CL & lset( * LevelsetP2CL::Create( MG, lsetbnddata, sf, P.get_child("Levelset")) );

    lset.SetYoungAngle(the_Young_angle);//set Young's Contact angle on the solid boundary
    lset.SetBndOutNormal(the_Bnd_outnormal);//set outnormal of the domain boundary
     if (is_periodic) //CL: Anyone a better idea? perDirection from ParameterFile?
    {
        DROPS::Point3DCL dx;
        //hack:
        std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
        size_t idx_;
        while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
            mesh[idx_]= ' ';
        std::istringstream brick_info( mesh);
        brick_info >> dx[0] >> dx[1] >> dx[2] ;
        int n = 0;
        if (P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicx" || P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicy" || P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicz")
            n = 1;
        if (P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicxy" || P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicxz" || P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicyz")
            n = 2;
        LevelsetP2CL::perDirSetT pdir(n);
        if (P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicx") pdir[0][0] = dx[0];
        if (P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicy") pdir[0][1] = dx[1];
        if (P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicz") pdir[0][2] = dx[2];
        if (P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicxy") {pdir[0][0] = dx[0]; pdir[1][1] = dx[1];}
        if (P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicxz") {pdir[0][0] = dx[0]; pdir[1][2] = dx[2];}
        if (P.get("DomainCond.PeriodicMatching", std::string("periodicx")) == "periodicyz") {pdir[0][1] = dx[1]; pdir[1][2] = dx[2];}
        if (P.get("DomainCond.PeriodicMatching", std::string("periodicx")) != "periodicx" && P.get("DomainCond.PeriodicMatching", std::string("periodicx")) != "periodicy" && P.get("DomainCond.PeriodicMatching", std::string("periodicx")) != "periodicz" &&
          P.get("DomainCond.PeriodicMatching", std::string("periodicx")) != "periodicxy" && P.get("DomainCond.PeriodicMatching", std::string("periodicx")) != "periodicxz" && P.get("DomainCond.PeriodicMatching", std::string("periodicx")) != "periodicyz"){
            std::cout << "WARNING: could not set periodic directions! Reparametrization can not work correctly now!" << std::endl;
            std::cout << "Press any key to continue" << std::endl; getchar();
        }
        lset.SetPeriodicDirections(&pdir);
    }

    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);
    VelocityRepairCL velrepair( Stokes);
    adap.push_back( &velrepair);
    PressureRepairCL prrepair( Stokes, lset);
    adap.push_back( &prrepair);

    MLIdxDescCL* lidx= &lset.idx;
    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;

    if ( StokesSolverFactoryHelperCL().VelMGUsed(P)){
        Stokes.SetNumVelLvl ( Stokes.GetMG().GetNumLevel());
        lset.SetNumLvl(Stokes.GetMG().GetNumLevel());
    }
    if ( StokesSolverFactoryHelperCL().PrMGUsed(P)){
        Stokes.SetNumPrLvl  ( Stokes.GetMG().GetNumLevel());
        lset.SetNumLvl(Stokes.GetMG().GetNumLevel());
    }
    lset.CreateNumbering( MG.GetLastLevel(), lidx, periodic_match);
    lset.Phi.SetIdx( lidx);

    if (lset.IsDiscontinuous())
    {
        LevelsetP2DiscontCL& lsetD (dynamic_cast<LevelsetP2DiscontCL&>(lset));
        MLIdxDescCL* lidxc = lsetD.idxC;
        lsetD.CreateNumbering( MG.GetLastLevel(), lidxc, periodic_match);
        lsetD.PhiContinuous.SetIdx( lidxc);
    }

    PermutationT lset_downwind;
    if (P.get<double>("SurfTens.VarTension"))
        lset.SetSurfaceForce( SF_ImprovedLBVar);
    else
        lset.SetSurfaceForce( SF_ImprovedLB);

    if ( StokesSolverFactoryHelperCL().VelMGUsed(P))
        Stokes.SetNumVelLvl ( Stokes.GetMG().GetNumLevel());
    if ( StokesSolverFactoryHelperCL().PrMGUsed(P))
        Stokes.SetNumPrLvl  ( Stokes.GetMG().GetNumLevel());

    SetInitialLevelsetConditions( lset, MG, P);

    double Vol = 0;
    //Vol = lset.GetVolume();
    if ((P.get("Exp.InitialLSet", std::string("Ellipsoid")) == "Ellipsoid" || P.get("Exp.InitialLSet", std::string("Ellipsoid")) == "Cylinder" 
        || P.get("Exp.InitialLSet", std::string("Ellipsoid")) == "ContactDroplet" || P.get("Exp.InitialLSet", std::string("Ellipsoid")) == "HalfEllipsoid" || P.get("Exp.InitialLSet", std::string("Ellipsoid")) == "TaylorFlowDistance") 
        && P.get<int>("Levelset.VolCorrection") != 0)
    {  
        if (P.get<double>("Exp.InitialVolume",-1.0) > 0 )
            Vol = P.get<double>("Exp.InitialVolume");      
        std::string InitialLSet= P.get("Exp.InitialLSet", std::string("Ellipsoid"));
        if (InitialLSet == "Ellipsoid")
            Vol = EllipsoidCL::GetVolume();
        if (InitialLSet == "HalfEllipsoid")
            Vol = HalfEllipsoidCL::GetVolume();
        if (InitialLSet == "ContactDroplet")
            Vol = ContactDropletCL::GetVolume();
        if (InitialLSet.find("Cylinder")==0)
            Vol = CylinderCL::GetVolume();
        std::cout << "initial rel. volume: " << lset.GetVolume()/Vol << std::endl;
        double dphi= lset.AdjustVolume( Vol, 1e-9);
        std::cout << "initial lset offset for correction is " << dphi << std::endl;
        lset.Phi.Data+= dphi;
        std::cout << "new initial rel. volume: " << lset.GetVolume()/Vol << std::endl;
    }else{
        Vol = lset.GetVolume();
    }

    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx, periodic_match);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, periodic_match, &lset);
    PermutationT vel_downwind;

    StokesVelBndDataCL::bnd_val_fun ZeroVel = InVecMap::getInstance().find("ZeroVel")->second;
    Stokes.SetIdx();
    Stokes.v.SetIdx  ( vidx);
    Stokes.p.SetIdx  ( pidx);
    Stokes.InitVel( &Stokes.v, ZeroVel);

    IteratedDownwindCL navstokes_downwind( P.get_child( "NavStokes.Downwind"));
    if (P.get<int>( "NavStokes.Downwind.Frequency") > 0) {
        if (StokesSolverFactoryHelperCL().VelMGUsed( P))
            throw DROPSErrCL( "Strategy: Multigrid-solver and downwind-numbering cannot be used together. Sorry.\n");
        vel_downwind= Stokes.downwind_numbering( lset, navstokes_downwind);
    }
    IteratedDownwindCL levelset_downwind( P.get_child( "Levelset.Downwind"));
    if (P.get<int>( "Levelset.Downwind.Frequency") > 0)
        lset_downwind= lset.downwind_numbering( Stokes.GetVelSolution(), levelset_downwind);

    DisplayDetailedGeom( MG);
    DisplayUnks(Stokes, lset, MG);

    // Stokes-Solver
    StokesSolverFactoryCL<InstatNavierStokes2PhaseP2P1CL> stokessolverfactory(Stokes, P);
    
    StokesSolverBaseCL* stokessolver = NULL;
    if (! P.get("Stokes.DirectSolve", 0))
        stokessolver = stokessolverfactory.CreateStokesSolver();
    else
        stokessolver = new StokesPardisoSolverCL(); 

    // Navier-Stokes-Solver
    NSSolverBaseCL<InstatNavierStokes2PhaseP2P1CL>* navstokessolver = 0;
    if (P.get<double>("NavStokes.Nonlinear")==0.0)
        navstokessolver = new NSSolverBaseCL<InstatNavierStokes2PhaseP2P1CL>(Stokes, *stokessolver);
    else
        navstokessolver = new AdaptFixedPtDefectCorrCL<InstatNavierStokes2PhaseP2P1CL>(Stokes, *stokessolver, P.get<int>("NavStokes.Iter"), P.get<double>("NavStokes.Tol"), P.get<double>("NavStokes.Reduction"));

    // Level-Set-Solver
#ifndef _PAR
    typedef GSPcCL  LsetPcT;
#else
    typedef JACPcCL LsetPcT;
#endif
    LsetPcT lset_pc;
    GMResSolverCL<LsetPcT>* gm = new GMResSolverCL<LsetPcT>( lset_pc, 200, P.get<int>("Levelset.Iter"), P.get<double>("Levelset.Tol"));

    LevelsetModifyCL lsetmod( P.get<int>("Reparam.Freq"), P.get<int>("Reparam.Method"), P.get<double>("Reparam.MaxGrad"), P.get<double>("Reparam.MinGrad"), P.get<int>("Levelset.VolCorrection"), Vol, is_periodic);

    UpdateProlongationCL<Point3DCL> PVel( Stokes.GetMG(), stokessolverfactory.GetPVel(), &Stokes.vel_idx, &Stokes.vel_idx);
    adap.push_back( &PVel);
    UpdateProlongationCL<double> PPr ( Stokes.GetMG(), stokessolverfactory.GetPPr(), &Stokes.pr_idx, &Stokes.pr_idx);
    adap.push_back( &PPr);
    UpdateProlongationCL<double> PLset( lset.GetMG(), lset.GetProlongation(), lset.idxC, lset.idxC);
    adap.push_back( &PLset);
    Stokes.P_ = stokessolverfactory.GetPVel();

//The velocity and pressure are read from files
    SetInitialConditions( Stokes, lset, MG, P);
        
    //read in velocity and pressure in different levels of meshes.
    VecDescCL VelMesh[4]; 
    VecDescCL PreMesh[4];
    int MeshLev = P.get<int>("Exp.MeshLevel");
    MLIdxDescCL  pr_idx(P1X_FE, MG.GetNumLevel(), Stokes.BndData_.Pr, 0, 0.1);
    pr_idx.CreateNumbering(MG.GetLastLevel(), MG, Stokes.BndData_.Pr, 0, &lset.Phi, &lset.GetBndData());
    pr_idx.UpdateXNumbering(MG, *lset.PhiC, lset.GetBndData());
    //std::cout<<"Size of vidx: "<<vidx->size()<<std::endl;
    //std::cout<<"Size of pidx: "<<pr_idx.size()<<std::endl;
    MLIdxDescCL::iterator it = vidx->begin();
    MLIdxDescCL::iterator itp = pr_idx.begin();
    for (int i=0; i< MeshLev; i++)
    {
        std::stringstream s;
        std::string filename = "Exp.InitialFile";
        s << filename << i;
        VelMesh[i].SetIdx(&(*it));
        PreMesh[i].SetIdx(&(*itp));
        std::cout<<"Number of unknows of velocity of level "<<i<<" is "<<(*it).NumUnknowns()<<std::endl;
        std::cout<<"Number of unknows of pressure of level "<<i<<" is "<<(*itp).NumUnknowns()<<std::endl;
        ReadFEFromFile( VelMesh[i], MG, P.get<std::string>(s.str())+"velocity", P.get<int>("Restart.Binary"));
        ReadFEFromFile( PreMesh[i], MG, P.get<std::string>(s.str())+"pressure", P.get<int>("Restart.Binary"), lset.PhiC);
        ++it;
        ++itp;
    }
    
    //prolongate the original matrices;
    VectorCL VelVecProlong[4];    //used to store the prolongated data of solutions 
    VecDescCL VelProlongMesh[4];  //the vector description objects with the data are the VelVecProlong;  
    MLDataCL<ProlongationCL<Point3DCL> > prolongVel;   //store the prolongation matrices
    
    //VectorCL PreVecProlong[4];    //used to store the prolongated data of solutions 
    //VecDescCL PreProlongMesh[4];  //the vector description objects with the data are the PreVecProlong;  
    //MLDataCL<ProlongationCL<Point3DCL> > prolongPre;   //store the prolongation matrices
    
    for(int i=0; i<MeshLev; i++)
    {
        VelVecProlong[i]= VelMesh[i].Data;
        //PreVecProlong[i]= PreMesh[i].Data;
    }
    SetupProlongationMatrix(MG, prolongVel, &Stokes.vel_idx , &Stokes.vel_idx);     //set up prolongation matrices
    //SetupProlongationMatrix(MG, prolongPre, &pr_idx , &pr_idx);     //set up prolongation matrices, not implemented yet
    int counter=1;
    for (MLDataCL<ProlongationCL<Point3DCL> >::iterator proVelPtr=(++prolongVel.begin()); proVelPtr != prolongVel.end(); ++proVelPtr)
    {
        for(int i=0; i<counter; i++)
        {
            VelVecProlong[i] = (*proVelPtr)* VelVecProlong[i];   // the original matricies are prolongated
        }
        if (counter < MeshLev)
            counter++;
    }
    /*int counter=1;
    for (MLDataCL<ProlongationCL<Point3DCL> >::iterator proPrePtr=(++prolongPre.begin()); proPrePtr != prolongPre.end(); ++proPrePtr)
    {
        for(int i=0; i<counter; i++)
        {
            PreVecProlong[i] = (*proPrePtr)* PreVecProlong[i];   // the original matricies are prolongated
        }
        if (counter < MeshLev)
            counter++;
    }*/
    for (int i=0; i < MeshLev; i++)
    {
        VelProlongMesh[i].SetIdx(vidx->GetFinestPtr());
        //PreProlongMesh[i].SetIdx(pr_idx.GetFinestPtr());
        VelProlongMesh[i].Data = VelVecProlong[i];     // use the vector description class to store the data; 
        //PreProlongMesh[i].Data = PreVecProlong[i];
        std::cout<<"The error of mesh 4"<<i<< " is :" <<std::endl;
        Stokes.CheckTwoPhaseSolution( &Stokes.v, &Stokes.p, lset, &VelProlongMesh[i], &Stokes.p);
    }
    std::cout<<"The L2 norm of Vel is :" <<std::endl;
    Stokes.CheckTwoPhaseSolution( &Stokes.v, &Stokes.p, lset, Stokes.Coeff_.RefVel, Stokes.Coeff_.RefPr);
    //Stokes.v.t += GetTimeOffset();
    std::cout << std::endl;
    delete navstokessolver;
    delete stokessolver;
    delete gm;
    delete &lset;
}

} // end of namespace DROPS


/// \brief Set Default parameters here s.t. they are initialized.
/// The result can be checked when Param-list is written to the output.
void SetMissingParameters(DROPS::ParamCL& P){
    P.put_if_unset<int>("VTK.ReUseTimeFile",0);
    P.put_if_unset<int>("VTK.UseDeformation",0);
    P.put_if_unset<int>("VTK.UseOnlyP1",0);
    P.put_if_unset<int>("VTK.AddP1XPressure",0);
    P.put_if_unset<int>("VTK.AddDGOutput",0);
    P.put_if_unset<int>("Transp.DoTransp",0);
    P.put_if_unset<std::string>("Restart.Inputfile","none");
    P.put_if_unset<int>("NavStokes.Downwind.Frequency", 0);
    P.put_if_unset<double>("NavStokes.Downwind.MaxRelComponentSize", 0.05);
    P.put_if_unset<double>("NavStokes.Downwind.WeakEdgeRatio", 0.2);
    P.put_if_unset<double>("NavStokes.Downwind.CrosswindLimit", std::cos( M_PI/6.));
    P.put_if_unset<int>("Levelset.Discontinuous", 0);
    P.put_if_unset<int>("Levelset.Downwind.Frequency", 0);
    P.put_if_unset<double>("Levelset.Downwind.MaxRelComponentSize", 0.05);
    P.put_if_unset<double>("Levelset.Downwind.WeakEdgeRatio", 0.2);
    P.put_if_unset<double>("Levelset.Downwind.CrosswindLimit", std::cos( M_PI/6.));

    P.put_if_unset<std::string>("Exp.VolForce", "ZeroVel");
    P.put_if_unset<double>("Mat.DensDrop", 0.0);
    P.put_if_unset<double>("Mat.ShearVisco", 0.0);
    P.put_if_unset<double>("Mat.DilatationalVisco", 0.0);
    P.put_if_unset<double>("SurfTens.ShearVisco", 0.0);
    P.put_if_unset<double>("SurfTens.DilatationalVisco", 0.0);

    P.put_if_unset<int>("General.ProgressBar", 0);
    P.put_if_unset<std::string>("General.DynamicLibsPrefix", "../");
	//contactangle problem--------------------------------------------
	P.put_if_unset<double>("SpeBnd.alpha", 0.0);
    P.put_if_unset<double>("SpeBnd.beta1", 0.0);
    P.put_if_unset<double>("SpeBnd.beta2", 0.0);
	P.put_if_unset<double>("SpeBnd.SmoothZone", 0.0);
	P.put_if_unset<std::string>("SpeBnd.CtAngle", "ConstantAngle");
	P.put_if_unset<double>("SpeBnd.contactangle", 0.0);
	P.put_if_unset<std::string>("SpeBnd.BndOutNormal", "OutNormalBrick");
	P.put_if_unset<std::string>("SpeBnd.posDrop", "[0.5, 0, 0.5 ]");

	P.put_if_unset<std::string>("Exp.Solution_Vel", "None");
	P.put_if_unset<std::string>("Exp.Solution_GradPr", "None");
	P.put_if_unset<std::string>("Exp.Solution_Pr", "None");
	P.put_if_unset<int>("Exp.OutputInfo",1);
	//---------------------------------------------------------------
   P.put_if_unset<double>("Exp.SimuType", 0.0);
   P.put_if_unset<double>("Stokes.epsP", 0.0);
   P.put_if_unset<double>("Stokes.DirectSolve", 0);
}

int main (int argc, char** argv)
{
#ifdef _PAR
    DROPS::ProcCL::Instance(&argc, &argv);
#endif
  try
  {
    std::cout << "Boost version: " << BOOST_LIB_VERSION << std::endl;

    DROPS::read_parameter_file_from_cmdline( P, argc, argv, "risingdroplet.json");
    SetMissingParameters(P);
    std::cout << P << std::endl;

    DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"), P.get<std::vector<std::string> >("General.DynamicLibs") );

    if (P.get<int>("General.ProgressBar"))
        DROPS::ProgressBarTetraAccumulatorCL::Activate();

    // check parameter file
    if (P.get<double>("SurfTens.DilatationalVisco")< P.get<double>("SurfTens.ShearVisco"))
    {
        throw DROPS::DROPSErrCL("Parameter error : Dilatational viscosity must be larger than surface shear viscosity");
    }

    DROPS::MatchMap & matchmap = DROPS::MatchMap::getInstance();
    bool is_periodic = P.get<std::string>("DomainCond.PeriodicMatching", "none") != "none";
    DROPS::match_fun periodic_match = is_periodic ? matchmap[P.get<std::string>("DomainCond.PeriodicMatching", "periodicx")] : 0;

    DROPS::MultiGridCL* mg= 0;
    typedef DROPS::BndDataCL<DROPS::Point3DCL> VelBndDataCL;
    typedef DROPS::BndDataCL<double>    PrBndDataCL;
    VelBndDataCL *velbnddata = 0;
    PrBndDataCL *prbnddata = 0;
    DROPS::LsetBndDataCL* lsetbnddata= 0;

    //you cannot pass a double& per P.get, so you need to use this indirect way
    double ExpRadInlet = P.get<double>("Exp.RadInlet");

    try
    {
        std::auto_ptr<DROPS::MGBuilderCL> builder( DROPS::make_MGBuilder( P));
        mg = new DROPS::MultiGridCL( *builder);
    }
    catch (DROPS::DROPSParamErrCL& e)
    {
        std::cout << "\n"
                  << "  /----------------------------------------------------------------\\ \n"
                  << "  | WARNING: It seems you are using the old domain descriptions    | \n"
                  << "  |          or your \"Domain\" section is not correct.              | \n"
                  << "  |          Please adapt your json-file to the new description.   | \n"
                  <<"  \\----------------------------------------------------------------/ \n"
                  << std::endl;
        DROPS::BuildDomain( mg, P.get<std::string>("DomainCond.MeshFile"), P.get<int>("DomainCond.GeomType"), P.get<std::string>("Restart.Inputfile"), ExpRadInlet);
    }




    P.put("Exp.RadInlet", ExpRadInlet);

    std::cout << "Generated MG of " << mg->GetLastLevel() << " levels." << std::endl;

    std::string perbndtypestr;
    std::string zerobndfun;
    for( size_t i= 1; i<=mg->GetBnd().GetNumBndSeg(); ++i) {
        zerobndfun += "Zero";
        if (i!=mg->GetBnd().GetNumBndSeg())
          zerobndfun += "!";
    }
    DROPS::BuildBoundaryData( mg, velbnddata, P.get<std::string>("DomainCond.BoundaryType"), P.get<std::string>("DomainCond.BoundaryFncs"), periodic_match, &perbndtypestr);
    std::cout << "Generated boundary conditions for velocity, ";
    DROPS::BuildBoundaryData( mg, prbnddata, perbndtypestr, zerobndfun, periodic_match);
    std::cout << "pressure, ";
    //DROPS::BuildBoundaryData( mg, lsetbnddata,  perbndtypestr, zerobndfun, periodic_match);
    //DROPS::BuildBoundaryData( mg, lsetbnddata,P.get<std::string>("DomainCond.BoundaryType"), zerobndfun, periodic_match);
    DROPS::BuildlsetBoundaryData( mg, lsetbnddata, perbndtypestr, P.get<std::string>("DomainCond.BoundaryType"), zerobndfun, periodic_match);
    //Hope this will not affect solving levelset equation!?
    std::cout << "and levelset." << std::endl;
    DROPS::StokesBndDataCL bnddata(*velbnddata,*prbnddata);

    std::string InitialLSet= P.get("Exp.InitialLSet", std::string("Ellipsoid"));
    if (InitialLSet == "Ellipsoid")
        DROPS::EllipsoidCL::Init( P.get<DROPS::Point3DCL>("Exp.PosDrop"), P.get<DROPS::Point3DCL>("Exp.RadDrop"));
    if (InitialLSet == "HalfEllipsoid")
        DROPS::HalfEllipsoidCL::Init( P.get<DROPS::Point3DCL>("Exp.PosDrop"), P.get<DROPS::Point3DCL>("Exp.RadDrop"));
    if (InitialLSet == "ContactDroplet")
        DROPS::ContactDropletCL::Init( P.get<DROPS::Point3DCL>("Exp.PosDrop"), P.get<DROPS::Point3DCL>("Exp.RadDrop"), P.get<double>("Exp.AngleDrop"));
    if  (InitialLSet == "TwoEllipsoid")
        DROPS::TwoEllipsoidCL::Init( P.get<DROPS::Point3DCL>("Exp.PosDrop"), P.get<DROPS::Point3DCL>("Exp.RadDrop"), P.get<DROPS::Point3DCL>("Exp.PosDrop2"), P.get<DROPS::Point3DCL>("Exp.RadDrop2"));
    if  (InitialLSet.find("Layer")==0){
      	DROPS::LayerCL::Init( P.get<DROPS::Point3DCL>("Exp.PosDrop"), P.get<DROPS::Point3DCL>("Exp.RadDrop"), InitialLSet[5]-'X');
    	P.put("Exp.InitialLSet", InitialLSet= "Layer");
    }
    if (InitialLSet.find("Cylinder")==0){
        DROPS::CylinderCL::Init( P.get<DROPS::Point3DCL>("Exp.PosDrop"), P.get<DROPS::Point3DCL>("Exp.RadDrop"), InitialLSet[8]-'X');
        P.put("Exp.InitialLSet", InitialLSet= "Cylinder");
    }
    typedef DROPS::DistMarkingStrategyCL MarkerT;
    MarkerT InitialMarker( DROPS::InScaMap::getInstance()[InitialLSet],
                           P.get<double>("AdaptRef.Width"),
                           P.get<double>("AdaptRef.CoarsestLevel"), P.get<double>("AdaptRef.FinestLevel") );

    DROPS::AdapTriangCL adap( *mg, &InitialMarker,
                              ((P.get<std::string>("Restart.Inputfile") == "none") ? P.get<int>("AdaptRef.LoadBalStrategy") : -P.get<int>("AdaptRef.LoadBalStrategy")));
      // If we read the Multigrid, it shouldn't be modified;
    // otherwise the pde-solutions from the ensight files might not fit.
    if (P.get("Restart.Inputfile", std::string("none")) == "none")
    {
        adap.MakeInitialTriang();
    }

    std::cout << DROPS::SanityMGOutCL(*mg) << std::endl;
#ifdef _PAR
    if ( DROPS::ProcCL::Check( DROPS::DiST::InfoCL::Instance().IsSane( std::cerr)))
        std::cout << " DiST-module seems to be alright!" << std::endl;
    else
        std::cout << " DiST-module seems to be broken!" << std::endl;
    if ( DROPS::CheckParMultiGrid())
        std::cout << "As far as I can tell the ParMultigridCL is sane\n";
#endif

    DROPS::InstatNavierStokes2PhaseP2P1CL prob( *mg, DROPS::TwoPhaseFlowCoeffCL(P), bnddata, P.get<double>("Stokes.XFEMStab")<0 ? DROPS::P1_FE : DROPS::P1X_FE, P.get<double>("Stokes.XFEMStab"));

    //DROPS::TimerCL time;
	//time.Reset();
    Strategy( prob, *lsetbnddata, adap);    // do all the stuff
    //time.Stop();
	//std::cout<<"In strategy function it took "<<time.GetTime()<<std::endl;
    delete mg;
    delete velbnddata;
    delete prbnddata;
    delete lsetbnddata;

    rusage usage;
    getrusage( RUSAGE_SELF, &usage);

#ifdef _PAR
    printf( "[%i]: ru_maxrss: %li kB.\n", DROPS::ProcCL::MyRank(), usage.ru_maxrss);
#else
    printf( "ru_maxrss: %li kB.\n", usage.ru_maxrss);
#endif
    std::cout << " twophasedrops finished regularly" << std::endl;
    return 0;
  }
  catch (DROPS::DROPSErrCL& err) { err.handle(); }
}

