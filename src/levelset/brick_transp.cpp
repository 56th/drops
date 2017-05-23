/// \file brick_transp.cpp
/// \brief two-phase flow in square pipe with transport
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt

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
#include "out/output.h"
#include "out/ensightOut.h"
#include "misc/params.h"
#include "levelset/mzelle_hdr.h"
#include "levelset/surfacetension.h"
#include "num/bndData.h"
#include "poisson/transport2phase.h"
#include "misc/funcmap.h"
#include <fstream>

#include "misc/progressaccu.h"
#include "misc/dynamicload.h"

DROPS::ParamCL P;

double Initialcneg (const DROPS::Point3DCL&, double)
{
    return 1.0;
}

double Initialcpos (const DROPS::Point3DCL&, double)
{
    return 0.3; // 0.5;
}

typedef DROPS::BndDataCL<DROPS::Point3DCL> VelBndDataCL;
typedef VelBndDataCL::bnd_val_fun  vel_bnd_val_fun;
typedef DROPS::BndDataCL<> cBndDataCL;
typedef cBndDataCL::bnd_val_fun  c_bnd_val_fun;

double D[2]= { /*pos. part*/ 10e-3, /*neg. part*/ 5e-3 };
double H= 0.5; // in the neg. part

namespace DROPS
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

typedef P2EvalCL<SVectorCL<3>, const VelBndDataCL, const VecDescCL> const_DiscVelSolCL;

void Strategy (MultiGridCL& MG, const LsetBndDataCL& lsbnd)
{
    DROPS::InScaMap & inscamap = DROPS::InScaMap::getInstance();

    const DROPS::BndCondT v_bc[6]= { DROPS::Dir0BC, DROPS::Dir0BC, DROPS::DirBC, DROPS::DirBC, DROPS::Dir0BC, DROPS::Dir0BC };
    const vel_bnd_val_fun v_bfun[6]= { 0, 0, DROPS::InVecMap::getInstance()["InflowBrickTransp"], DROPS::InVecMap::getInstance()["InflowBrickTransp"], 0, 0};
    const DROPS::BndCondT c_bc[6]= { DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC };
    const c_bnd_val_fun c_bfun[6]= { & Initialcpos,  & Initialcpos, & Initialcpos,& Initialcpos, & Initialcpos, & Initialcpos };

    // initialization of surface tension
    // choose a proper model for surface tension coefficient, see levelset/surfacetension.h
    SurfaceTensionCL * sf;
    sf = new SurfaceTensionCL( inscamap[P.get<std::string>("NavStokes.Coeff.SurfTens.VarTensionFunc")]);
    sf->SetInputMethod( Sigma_X);

    // LevelsetP2CL lset( MG, lsbnd, sf, P.get<double>("Levelset.SD"), P.get<double>("Levelset.CurvDiff"));
    // Creates new Levelset-Object, has to be cleaned manually
    LevelsetP2CL & lset( * LevelsetP2CL::Create( MG, lsbnd, *sf, P.get_child("Levelset")) );

    lset.CreateNumbering( MG.GetLastLevel());
    lset.Init( EllipsoidCL::DistanceFct);

    VelBndDataCL Bnd_v( 6, v_bc, v_bfun);
    IdxDescCL  vidx( vecP2_FE, Bnd_v);
    vidx.CreateNumbering( MG.GetLastLevel(), MG);
    VecDescCL v( &vidx);
    InitVel( MG, v, DROPS::InVecMap::getInstance()["InflowBrickTransp"]);

    cBndDataCL Bnd_c( 6, c_bc, c_bfun);
    const double dt= P.get<double>("Time.FinalTime")/P.get<int>("Time.NumSteps");

    TransportP1CL c( MG, Bnd_c, Bnd_v, /*theta*/ 0.5, D, H, &v, lset,
        dt, P.get<int>("Transp.Solver.Iter"), P.get<double>("Transp.Solver.Tol"));
    MLIdxDescCL* cidx= &c.idx;
    c.CreateNumbering( MG.GetLastLevel(), cidx);
    c.ct.SetIdx( cidx);
    c.Init( &Initialcneg, &Initialcpos);
    c.Update();

    // Initialize Ensight6 output
    std::string ensf( P.get<std::string>("Ensight.EnsDir") + "/" + P.get<std::string>("Ensight.EnsCase"));
    Ensight6OutCL ensight( P.get<std::string>("Ensight.EnsCase") + ".case", P.get<int>("Time.NumSteps") + 1);
    ensight.Register( make_Ensight6Geom  ( MG, MG.GetLastLevel(),   "Messzelle",     ensf + ".geo"));
    ensight.Register( make_Ensight6Scalar( lset.GetSolution(),      "Levelset",      ensf + ".scl", true));
    ensight.Register( make_Ensight6Vector( const_DiscVelSolCL( &v, &Bnd_v, &MG),
                                                                    "Velocity",      ensf + ".vel", true));
    ensight.Register( make_Ensight6Scalar( c.GetSolution(),         "Concentration", ensf + ".c",   true));
    ensight.Register( make_Ensight6Scalar( c.GetSolution( c.ct),    "TransConc",     ensf + ".ct",  true));

    MG.SizeInfo( std::cout);
    std::cout << c.c.Data.size() << " concentration unknowns,\n";
    std::cout << v.Data.size() << " velocity unknowns,\n";
    std::cout << lset.Phi.Data.size() << " levelset unknowns.\n";

    const double Vol= EllipsoidCL::GetVolume();
    std::cout << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;

    if (P.get<int>("Ensight.EnsightOut"))
        ensight.Write();

    c.SetTimeStep( dt);
    for (int step= 1; step <= P.get<int>("Time.NumSteps"); ++step) {
        std::cout << "======================================================== Schritt " << step << ":\n";
        c.DoStep( step*dt);
        if (P.get<int>("Ensight.EnsightOut"))
            ensight.Write( step*dt);
    }
    std::cout << std::endl;
    delete &lset;
    if(sf) delete sf;
}

} // end of namespace DROPS


int main (int argc, char** argv)
{
  try {
    DROPS::read_parameter_file_from_cmdline( P, argc, argv, "../../param/levelset/brick_transp/brick_transp.json");
    std::cout << P << std::endl;

    DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"), P.get<std::vector<std::string> >("General.DynamicLibs") );

    if (P.get<int>("General.ProgressBar"))
        DROPS::ProgressBarTetraAccumulatorCL::Activate();

    DROPS::Point3DCL e1(0.), e2(0.), e3(0.), orig;
    e1[0]=2*P.get<double>("Inflow.RadInlet"); e2[1]=1.0; e3[2]= 2*P.get<double>("Inflow.RadInlet");
    DROPS::BrickBuilderCL builder( orig, e1, e2, e3, 20, 20, 20);
    DROPS::MultiGridCL mg( builder);

    const DROPS::BndCondT bcls[6]= { DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC };
    const DROPS::LsetBndDataCL::bnd_val_fun bfunls[6]= { 0,0,0,0,0,0};
    DROPS::LsetBndDataCL lsbnd( 6, bcls, bfunls);

    std::cout << DROPS::SanityMGOutCL( mg) << std::endl;
    DROPS::EllipsoidCL::Init( P.get<DROPS::Point3DCL>("Levelset.PosDrop"), P.get<DROPS::Point3DCL>("Levelset.RadDrop"));

    Strategy( mg, lsbnd);    // do all the stuff

    std::cout << " brick_transp finished regularly" << std::endl;
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}

