/// \file mzelle_hdr.h
/// \brief header for all main programs simulating the measuring cell
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Yuanjun Zhang; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_MZELLE_HDR_H
#define DROPS_MZELLE_HDR_H

#include "geom/multigrid.h"
#include "levelset/levelset.h"
#include "num/discretize.h"
#include "stokes/instatstokes2phase.h"
#include "poisson/transport2phase.h"
#include "geom/geomselect.h"
#include "num/nssolver.h"
#include "levelset/coupling.h"
#include "spacetimetransp/spacetime_sol.h"
#include "misc/funcmap.h"
#include "misc/params.h"

extern DROPS::ParamCL P;
namespace DROPS
{


/// The level set has two disjoint singular points on the p[1]-axis; contained in \f$(-0.3,0.3)\times(-1,1)\times(-0.3,0.3)\f$.
double zitrus (const Point3DCL& p, double)
{
    return 0.2*(std::pow( p[0], 2) + std::pow( p[2], 2)) + std::pow( p[1], 3)*std::pow( p[1] - 1., 3);
}
static DROPS::RegisterScalarFunction regsca_zitrus2("zitrus", zitrus);

/// The level set is a torus, wrapped around p[2]-axis, that touches itself in the origin; contained in \f$(-1,1)\times(-1,1)\times(-0.5,0.5)\f$.
double dullo (const Point3DCL& p, double)
{
    return std::pow( p.norm_sq(), 2) - (std::pow( p[0], 2) + std::pow( p[1], 2));
}
/// The level set is a heart with an inward and an outward cusp on the p[2]-axis; contained in \f$(-1.2,1.2)\times(-0.7,0.7)\times(-1.1,1.3)\f$.
double suess (const Point3DCL& p, double)
{
    return std::pow( std::pow( p[0], 2) + 2.25*std::pow( p[1], 2) + std::pow( p[2], 2) - 1., 3) - std::pow( p[0], 2)*std::pow( p[2], 3) - 9./80.*std::pow( p[1], 2)*std::pow( p[2], 3);
}

class EllipsoidCL
{
  private:
    static Point3DCL Mitte_;
    static Point3DCL Radius_;

  public:
    EllipsoidCL( const Point3DCL& Mitte, const Point3DCL& Radius)
    { Init( Mitte, Radius); }
    static void Init( const Point3DCL& Mitte, const Point3DCL& Radius)
    { Mitte_= Mitte;    Radius_= Radius; }
    static double DistanceFct( const Point3DCL& p, double)
    {
        Point3DCL d= p - Mitte_;
        const double avgRad= cbrt(Radius_[0]*Radius_[1]*Radius_[2]);
        d/= Radius_;
        return std::abs( avgRad)*d.norm() - avgRad;
    }
    static double GetVolume() { return 4./3.*M_PI*Radius_[0]*Radius_[1]*Radius_[2]; }
    static Point3DCL& GetCenter() { return Mitte_; }
    static Point3DCL& GetRadius() { return Radius_; }
};

Point3DCL EllipsoidCL::Mitte_;
Point3DCL EllipsoidCL::Radius_;

static DROPS::RegisterScalarFunction regsca_ellipsoid("Ellipsoid", DROPS::EllipsoidCL::DistanceFct);

class HalfEllipsoidCL
{
  private:
    static Point3DCL Mitte_;
    static Point3DCL Radius_;

  public:
    HalfEllipsoidCL( const Point3DCL& Mitte, const Point3DCL& Radius)
    { Init( Mitte, Radius); }
    static void Init( const Point3DCL& Mitte, const Point3DCL& Radius)
    { Mitte_= Mitte;    Radius_= Radius; }
    static double DistanceFct( const Point3DCL& p, double)
    {
        Point3DCL d= p - Mitte_;
        const double avgRad= cbrt(Radius_[0]*Radius_[1]*Radius_[2]);
        d/= Radius_;
        return std::abs( avgRad)*d.norm() - avgRad;
    }
    static double GetVolume() { return 2./3.*M_PI*Radius_[0]*Radius_[1]*Radius_[2]; }
    static Point3DCL& GetCenter() { return Mitte_; }
    static Point3DCL& GetRadius() { return Radius_; }
};

Point3DCL HalfEllipsoidCL::Mitte_;
Point3DCL HalfEllipsoidCL::Radius_;

static DROPS::RegisterScalarFunction regsca_halfellipsoid("HalfEllipsoid", DROPS::HalfEllipsoidCL::DistanceFct);

class ContactDropletCL
{
  private:
    static Point3DCL Mitte_;
    static Point3DCL Radius_;
    static double Angle_;

  public:
    ContactDropletCL( const Point3DCL& Mitte, const Point3DCL& Radius, const double Angle)
    { Init( Mitte, Radius, Angle); }
    static void Init( const Point3DCL& Mitte, const Point3DCL& Radius, const double Angle)
    { Mitte_= Mitte;    Radius_= Radius;  Angle_= Angle;}
    static double DistanceFct( const Point3DCL& p, double)
    {
        Point3DCL d= p - Mitte_;
        const double avgRad= cbrt(Radius_[0]*Radius_[1]*Radius_[2]);
        d/= Radius_;
        return std::abs( avgRad)*d.norm() - avgRad;
    }
    static double GetVolume() 
    {
        double R_Angle_ = Angle_/180 * M_PI;
        return 1./3.*M_PI*Radius_[0]*Radius_[1]*Radius_[2] * (2.- 2.*std::cos(R_Angle_)- std::sin(R_Angle_)*std::sin(R_Angle_)*std::cos(R_Angle_)); 
    }
    static Point3DCL& GetCenter() { return Mitte_; }
    static Point3DCL& GetRadius() { return Radius_; }
    static double GetAngle() { return Angle_; }
};

Point3DCL ContactDropletCL::Mitte_;
Point3DCL ContactDropletCL::Radius_;
double ContactDropletCL::Angle_;

static DROPS::RegisterScalarFunction regsca_contactdroplet("ContactDroplet", DROPS::ContactDropletCL::DistanceFct);

/// \brief Represents a cylinder with ellipsoidal cross section
class CylinderCL
{
  private:
    static Point3DCL Mitte_;
    static Point3DCL Radius_;  ///< stores radii of ellipse and axis length (latter stored in entry given by axisDir_)
    static Usint     axisDir_; ///< direction of axis (one of 0,1,2)

  public:
    CylinderCL( const Point3DCL& Mitte, const Point3DCL& RadiusLength, Usint axisDir)
    { Init( Mitte, RadiusLength, axisDir); }
    static void Init( const Point3DCL& Mitte, const Point3DCL& RadiusLength, Usint axisDir)
    { Mitte_= Mitte;    Radius_= RadiusLength;    axisDir_= axisDir; }
    static double DistanceFct( const Point3DCL& p, double)
    {
        Point3DCL d= p - Mitte_;
        const double avgRad= std::sqrt(Radius_[0]*Radius_[1]*Radius_[2]/Radius_[axisDir_]);
        d/= Radius_;
        d[axisDir_]= 0;
        return std::abs( avgRad)*d.norm() - avgRad;
    }
    static double GetVolume() { return M_PI*Radius_[0]*Radius_[1]*Radius_[2]; }
    static Point3DCL& GetCenter() { return Mitte_; }
};

Point3DCL CylinderCL::Mitte_;
Point3DCL CylinderCL::Radius_;
Usint     CylinderCL::axisDir_;

static DROPS::RegisterScalarFunction regsca_cylinder("Cylinder", DROPS::CylinderCL::DistanceFct);


/// \brief Represents a thin layer bounded by two parallel planes perpendicular to one axis
class LayerCL
{
  private:
    static Point3DCL Mitte_;
    static Point3DCL Radius_;  ///< stores the width of the layer and axis lengths (former stored in entry given by axisDir_)
    static Usint     axisDir_; ///< direction of axis (one of 0,1,2)

  public:
    LayerCL( const Point3DCL& Mitte, const Point3DCL& RadiusLength, Usint axisDir)
    { Init( Mitte, RadiusLength, axisDir); }
    static void Init( const Point3DCL& Mitte, const Point3DCL& RadiusLength, Usint axisDir)
    { Mitte_= Mitte;    Radius_= RadiusLength;    axisDir_= axisDir; }
    static double DistanceFct( const Point3DCL& p, double)
    {
        double d= p[axisDir_] - Mitte_[axisDir_];

        return std::abs( d)  - Radius_[axisDir_]/2;
    }
    static double GetVolume() { return Radius_[0]*Radius_[1]*Radius_[2]; }
    static Point3DCL& GetCenter() { return Mitte_; }
};

Point3DCL LayerCL::Mitte_;
Point3DCL LayerCL::Radius_;
Usint     LayerCL::axisDir_;

static DROPS::RegisterScalarFunction regsca_layer("Layer", DROPS::LayerCL::DistanceFct);

// collision setting (rising butanol droplet in water)
//  RadDrop1 =  1.50e-3  1.500e-3  1.50e-3
//  PosDrop1 =  6.00e-3  3.000e-3  6.00e-3
//  RadDrop2 =  0.75e-3  0.750e-3  0.75e-3
//  PosDrop2 =  6.00e-3  5.625e-3  6.00e-3
//  MeshFile =  12e-3x30e-3x12e-3@4x10x4

class TwoEllipsoidCL
{
  private:
    static Point3DCL Mitte1_,  Mitte2_;
    static Point3DCL Radius1_, Radius2_;

  public:
    TwoEllipsoidCL( const Point3DCL& Mitte1, const Point3DCL& Radius1,
                    const Point3DCL& Mitte2, const Point3DCL& Radius2)
    { Init( Mitte1, Radius1, Mitte2, Radius2); }
    static void Init( const Point3DCL& Mitte1, const Point3DCL& Radius1,
                      const Point3DCL& Mitte2, const Point3DCL& Radius2)
    { Mitte1_= Mitte1;    Radius1_= Radius1;
      Mitte2_= Mitte2;    Radius2_= Radius2;}
    static double DistanceFct( const Point3DCL& p, double)
    {
        Point3DCL d1= p - Mitte1_;
        const double avgRad1= cbrt(Radius1_[0]*Radius1_[1]*Radius1_[2]);
        d1/= Radius1_;
        Point3DCL d2= p - Mitte2_;
        const double avgRad2= cbrt(Radius2_[0]*Radius2_[1]*Radius2_[2]);
        d2/= Radius2_;
        return std::min(std::abs( avgRad1)*d1.norm() - avgRad1, std::abs( avgRad2)*d2.norm() - avgRad2);
    }
    static double GetVolume() { return 4./3.*M_PI*Radius1_[0]*Radius1_[1]*Radius1_[2]
        + 4./3.*M_PI*Radius2_[0]*Radius2_[1]*Radius2_[2]; }
};

Point3DCL TwoEllipsoidCL::Mitte1_;
Point3DCL TwoEllipsoidCL::Radius1_;
Point3DCL TwoEllipsoidCL::Mitte2_;
Point3DCL TwoEllipsoidCL::Radius2_;

static DROPS::RegisterScalarFunction regsca_twoellipsoid("TwoEllipsoid", DROPS::TwoEllipsoidCL::DistanceFct);

class TwoPhaseInfoBaseCL
{
  protected:
    std::ofstream* file_;    ///< write information, to this file

    void Update (const LevelsetP2CL& ls) {
        std::pair<double, double> h= h_interface( ls.GetMG().GetTriangEdgeBegin( ls.PhiC->RowIdx->TriangLevel()), ls.GetMG().GetTriangEdgeEnd( ls.PhiC->RowIdx->TriangLevel()), *ls.PhiC);
        h_min= h.first; h_max= h.second;
    }

  public:
    double h_min, h_max;

    TwoPhaseInfoBaseCL() : file_(nullptr), h_min(0), h_max(0) {}

    virtual void WriteHeader() = 0;
    virtual void Write (double time) = 0;
    /// \brief Set file for writing
    void Init(std::ofstream* file) { file_= file; }
};

class InterfaceInfoCL: public TwoPhaseInfoBaseCL
{
    typedef TwoPhaseInfoBaseCL BaseT;
  public:
    Point3DCL bary, vel, min, max;
    double maxGrad, Vol, surfArea, sphericity;


    template<class DiscVelSolT>
    void Update (const LevelsetP2CL& ls, const DiscVelSolT& u) {
        BaseT::Update( ls);
        ls.GetInfo( maxGrad, Vol, bary, vel, u, min, max, surfArea);
        // sphericity is the ratio of surface area of a sphere of same volume and surface area of the approximative interface
        sphericity= std::pow(6*Vol, 2./3.)*std::pow(M_PI, 1./3.)/surfArea;
    }
    void WriteHeader() {
        if (file_)
          (*file_) << "# time maxGradPhi volume bary_drop min_drop max_drop vel_drop h_min h_max surfArea sphericity\n";
    }
    void Write (double time) {
        if (file_)
          (*file_) << time << " " << maxGrad << " " << Vol << " " << bary << " " << min << " " << max << " " << vel << " " << h_min << " " << h_max << " " << surfArea << " " << sphericity << std::endl;
    }
} IFInfo;

class FilmInfoCL: public TwoPhaseInfoBaseCL
{
    typedef TwoPhaseInfoBaseCL BaseT;
  private:
    double x_, z_;
  public:
    Point3DCL vel;
    double maxGrad, Vol, point_h;

    template<class DiscVelSolT>
    void Update (const LevelsetP2CL& ls, const DiscVelSolT& u) {
        BaseT::Update( ls);
        ls.GetFilmInfo( maxGrad, Vol, vel, u, x_, z_, point_h);
    }
    void WriteHeader() {
        if (file_) {
            (*file_) << "Film Reynolds number Re_f = "
                     << P.get<double>("Mat.DensFluid")*P.get<double>("Mat.DensFluid")*P.get<DROPS::Point3DCL>("Exp.Gravity")[0]*std::pow(P.get<double>("Exp.Thickness"),3)/P.get<double>("Mat.ViscFluid")/P.get<double>("Mat.ViscFluid")/3 << '\n';
            (*file_) << "max. inflow velocity at film surface = "
                     << P.get<double>("Mat.DensFluid")*P.get<DROPS::Point3DCL>("Exp.Gravity")[0]*P.get<double>("Exp.Thickness")*P.get<double>("Exp.Thickness")/P.get<double>("Mat.ViscFluid")/2 << '\n';
            (*file_) << "# time maxGradPhi volume vel_drop h_min h_max point_x point_z point_h\n";
        }
    }
    void Write (double time) {
        if (file_)
          (*file_) << time << " " << maxGrad << " " << Vol << " " << vel << " " << h_min << " " << h_max << " " << x_ << " " << z_<< " " << point_h << std::endl;
    }
    /// \brief Set file for writing and film length (x) and width (z). Note: film height (y) is determined from level set.
    void Init(std::ofstream* file, double x, double z) { file_= file; x_= x; z_= z; }
} FilmInfo;

/// \brief factory for the time discretization schemes
template<class LevelSetSolverT>
TimeDisc2PhaseCL* CreateTimeDisc( InstatNavierStokes2PhaseP2P1CL& Stokes, LevelsetP2CL& lset,
    NSSolverBaseCL<InstatNavierStokes2PhaseP2P1CL>* stokessolver, LevelSetSolverT* lsetsolver, ParamCL& P, LevelsetModifyCL& lsetmod)
{
    int numSteps        = P.get<int>("Time.NumSteps");

    if( numSteps == 0)
        return 0;

    double theta        = P.get<double>("Time.Theta");
    double lsTheta      = theta;
    double tEnd         = P.get<double>("Time.FinalTime");
    double stepSize     = tEnd / numSteps;
    double implLB       = P.get<double>("CouplingSolver.ImplLB");
    double couplingTol  = P.get<double>("CouplingSolver.Tol");
    double couplingProj = P.get<double>("CouplingSolver.Projection");
    double nonlinear    = P.get<double>("CouplingSolver.NavStokesSolver.Nonlinear");

    switch (P.get<int>("Time.Scheme"))
    {
        case 1 :
            return (new LinThetaScheme2PhaseCL<LevelSetSolverT>
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, stepSize, theta, lsTheta, nonlinear, implLB));
        break;
        case 3 :
            std::cout << "[WARNING] use of ThetaScheme2PhaseCL is deprecated using RecThetaScheme2PhaseCL instead\n";
        case 2 :
            return (new RecThetaScheme2PhaseCL<LevelSetSolverT >
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, stepSize, couplingTol, theta, lsTheta, nonlinear, couplingProj, implLB));
        break;
        case 4 :
            return (new OperatorSplitting2PhaseCL<LevelSetSolverT>
                        (Stokes, lset, stokessolver->GetStokesSolver(), *lsetsolver, lsetmod, stepSize, P.get<int>("CouplingSolver.NavStokesSolver.OseenSolver.Iter"), P.get<double>("CouplingSolver.NavStokesSolver.OseenSolver.Tol"), nonlinear));
        break;
        case 6 :
            return (new SpaceTimeDiscTheta2PhaseCL<LevelSetSolverT>
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, stepSize, couplingTol, theta, lsTheta, nonlinear, couplingProj, implLB, false));
        break;
        case 7 :
            return (new SpaceTimeDiscTheta2PhaseCL< LevelSetSolverT>
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, stepSize, couplingTol, theta, lsTheta, nonlinear, couplingProj, implLB, true));
        break;
        case 8 :
            return (new EulerBackwardScheme2PhaseCL<LevelSetSolverT>
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, stepSize, couplingTol, nonlinear, couplingProj, implLB));
        break;
        case 9 :
            return (new CrankNicolsonScheme2PhaseCL<RecThetaScheme2PhaseCL, LevelSetSolverT>
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, stepSize, couplingTol, nonlinear, couplingProj, implLB));
        break;
        case 10 :
            return (new CrankNicolsonScheme2PhaseCL<SpaceTimeDiscTheta2PhaseCL, LevelSetSolverT>
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, stepSize, couplingTol, nonlinear, couplingProj, implLB));
        break;
        case 11 :
            return (new FracStepScheme2PhaseCL<RecThetaScheme2PhaseCL, LevelSetSolverT >
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, stepSize, couplingTol, nonlinear, couplingProj, implLB));
        break;
        case 12 :
            return (new FracStepScheme2PhaseCL<SpaceTimeDiscTheta2PhaseCL, LevelSetSolverT >
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, stepSize, couplingTol, nonlinear, couplingProj, implLB));
        break;
        case 13 :
            return (new Frac2StepScheme2PhaseCL<RecThetaScheme2PhaseCL, LevelSetSolverT >
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, stepSize, couplingTol, nonlinear, couplingProj, implLB));
        break;
        case 14 :
            return (new Frac2StepScheme2PhaseCL<SpaceTimeDiscTheta2PhaseCL, LevelSetSolverT >
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, stepSize, couplingTol, nonlinear, couplingProj, implLB));
        break;
        default : throw DROPSErrCL("Unknown TimeDiscMethod");
    }
}

template <class StokesT>
void SolveStatProblem( StokesT& Stokes, LevelsetP2CL& lset, NSSolverBaseCL<StokesT >& solver)
{
#ifndef _PAR
    TimerCL time;
#else
    ParTimerCL time;
#endif
    double duration;
    time.Reset();
    VelVecDescCL cplM, cplN;
    VecDescCL curv;
    cplM.SetIdx( &Stokes.vel_idx);
    cplN.SetIdx( &Stokes.vel_idx);
    curv.SetIdx( &Stokes.vel_idx);
    Stokes.SetIdx();
    Stokes.SetLevelSet( lset);
    lset.UpdateMLPhi();
    lset.AccumulateBndIntegral( curv);
    Stokes.AccumulateYoungForce (lset, curv);
    Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &cplM, lset, Stokes.v.t);
    Stokes.SetupPrStiff( &Stokes.prA, lset);
    Stokes.SetupPrMass ( &Stokes.prM, lset);
    Stokes.SetupSystem2( &Stokes.B, &Stokes.C, &Stokes.c, lset, Stokes.v.t);
    time.Stop();
    duration = time.GetTime();
    std::cout << "Discretizing took "<< duration << " sec.\n";
    time.Reset();
    Stokes.b.Data += curv.Data;
    solver.Solve( Stokes.A.Data, Stokes.B.Data, Stokes.C.Data, Stokes.v, Stokes.p.Data, Stokes.b.Data, cplN, Stokes.c.Data, Stokes.vel_idx.GetEx(), Stokes.pr_idx.GetEx(), 1.0);
    time.Stop();
    duration = time.GetTime();
    std::cout << "Solving (Navier-)Stokes took "<<  duration << " sec.\n";
    std::cout << "iter: " << solver.GetIter() << "\tresid: " << solver.GetResid() << std::endl;
}

bool ReadInitialConditionFromFile( const ParamCL &P )
{
    if( P.exists( std::string("Restart.InputData")) )
    {
        std::string inputDat = P.get<std::string>("Restart.InputData");
        if( inputDat.empty() || inputDat == "none" )
            return false;

        return true;
    }
    return false;
}

void SetInitialLevelsetConditions( LevelsetP2CL& lset, MultiGridCL& MG, ParamCL& P)
{
    if( ReadInitialConditionFromFile( P ) )
    {
        ReadFEFromFile( lset.Phi, MG, P.get<std::string>("Restart.InputData")+"levelset",
                        P.get<int>("Restart.Binary") );
        return;
    }

    switch ( P.get<int>("NavStokes.InitialValue") )
    {
#ifndef _PAR
      case -10: // read from ensight-file [deprecated]
      {
        std::cout << "[DEPRECATED] read from ensight-file [DEPRECATED]\n";
        //        ReadEnsightP2SolCL reader( MG);
        //        reader.ReadScalar( P.get<std::string>("DomainCond.InitialFile")+".scl", lset.Phi, lset.GetBndData());
      } break;
#endif
      case -1: // read from file
      {
        throw DROPSErrCL("Specify \"Restart.InputData\" to read initial condition from file");
        //ReadFEFromFile( lset.Phi, MG, P.get<std::string>("Restart.InputData")+"levelset", P.get<int>("Restart.Binary"));
      } break;
      case 0: case 1: case 10:
          //lset.Init( EllipsoidCL::DistanceFct);
          lset.Init( DROPS::InScaMap::getInstance()[P.get("Levelset.InitialValue", std::string("Ellipsoid"))]);
        break;
      case  2: //flow without droplet
          lset.Init( DROPS::InScaMap::getInstance()["One"]);
      break;
      default : throw DROPSErrCL("Unknown initial condition");
    }
}

template <typename StokesT>
void SetInitialConditions(StokesT& Stokes, LevelsetP2CL& lset, MultiGridCL& MG, const ParamCL& P)
{
    MLIdxDescCL* pidx= &Stokes.pr_idx;

    if( ReadInitialConditionFromFile( P ) )
    {
        ReadFEFromFile( Stokes.v, MG, P.get<std::string>("Restart.InputData")+"velocity",
                        P.get<int>("Restart.Binary"));
        Stokes.UpdateXNumbering( pidx, lset);
        Stokes.p.SetIdx( pidx);
        // pass also level set, as p may be extended
        ReadFEFromFile( Stokes.p, MG, P.get<std::string>("Restart.InputData")+"pressure",
                        P.get<int>("Restart.Binary"), lset.PhiC);
    }

    const int inival= P.get<int>("NavStokes.InitialValue");
    switch ( inival )
    {
#ifndef _PAR
      case -10: // read from ensight-file [deprecated]
      {
        std::cout << "[DEPRECATED] read from ensight-file [DEPRECATED]\n";
        /*        ReadEnsightP2SolCL reader( MG);
        reader.ReadVector( P.get<std::string>("DomainCond.InitialFile")+".vel", Stokes.v, Stokes.GetBndData().Vel);
        Stokes.UpdateXNumbering( pidx, lset);
        Stokes.p.SetIdx( pidx);
        if (Stokes.UsesXFEM()) {
            VecDescCL pneg( pidx), ppos( pidx);
            reader.ReadScalar( P.get<std::string>("DomainCond.InitialFile")+".prNeg", pneg, Stokes.GetBndData().Pr);
            reader.ReadScalar( P.get<std::string>("DomainCond.InitialFile")+".prPos", ppos, Stokes.GetBndData().Pr);
            P1toP1X ( pidx->GetFinest(), Stokes.p.Data, pidx->GetFinest(), ppos.Data, pneg.Data, lset.Phi, MG);
        }
        else
        reader.ReadScalar( P.get<std::string>("DomainCond.InitialFile")+".pr", Stokes.p, Stokes.GetBndData().Pr);*/
      } break;
#endif
      case -1: // read from file
      {
        throw DROPSErrCL("Specify \"Restart.InputData\" to read initial condition from file");
//        ReadFEFromFile( Stokes.v, MG, P.get<std::string>("Restart.InputData")+"velocity", P.get<int>("Restart.Binary"));
//        Stokes.UpdateXNumbering( pidx, lset);
//        Stokes.p.SetIdx( pidx);
//        ReadFEFromFile( Stokes.p, MG, P.get<std::string>("Restart.InputData")+"pressure", P.get<int>("Restart.Binary"), lset.PhiC); // pass also level set, as p may be extended
      } break;
      case 0: // zero initial condition
          Stokes.UpdateXNumbering( pidx, lset);
          Stokes.p.SetIdx( pidx);
        break;
      case 1: case 10: // stationary flow with/without surface tension
      {
        Stokes.UpdateXNumbering( pidx, lset);
        Stokes.p.SetIdx( pidx);
#ifdef _PAR
        typedef JACPcCL PcT;
#else
        typedef SSORPcCL PcT;
#endif
        PcT pc;
        PCGSolverCL<PcT> PCGsolver( pc, P.get<int>("CouplingSolver.NavStokesSolver.OseenSolver.PcAIter"), P.get<double>("CouplingSolver.NavStokesSolver.OseenSolver.PcATol"), true);
        typedef SolverAsPreCL<PCGSolverCL<PcT> > PCGPcT;
        PCGPcT apc( PCGsolver);
        SchurPreBaseCL *schurpc=0;
        const int iter= P.get<int>("CouplingSolver.NavStokesSolver.OseenSolver.Iter");
        const double tol= P.get<double>("CouplingSolver.NavStokesSolver.OseenSolver.Tol"),
                PcSTol= P.get<double>("CouplingSolver.NavStokesSolver.OseenSolver.PcSTol");
        if( Stokes.usesGhostPen() )
            schurpc= new ISGhPenPreCL( &Stokes.prA.Data.GetFinest(), &Stokes.prM.Data.GetFinest(), &Stokes.C.Data.GetFinest(), 0., 1., PcSTol, PcSTol, 200);
        else
            schurpc= new ISBBTPreCL( &Stokes.B.Data.GetFinest(), &Stokes.prM.Data.GetFinest(), &Stokes.M.Data.GetFinest(), Stokes.pr_idx.GetFinest(), 0.0, 1.0, PcSTol, PcSTol);
        typedef BlockPreCL<ExpensivePreBaseCL, SchurPreBaseCL, LowerBlockPreCL>    LowerBlockPcT;
        LowerBlockPcT blockPc(apc, *schurpc);
        GCRSolverCL< LowerBlockPcT > gcrsolver( blockPc, iter, iter, tol, false);
        BlockMatrixSolverCL< GCRSolverCL<LowerBlockPcT> > ssolver( gcrsolver );

        if (inival==1) { // with surface tension
            NSSolverBaseCL<StokesT> stokessolver( Stokes, ssolver);
            SolveStatProblem( Stokes, lset, stokessolver);
        } else { // without surface tension
            TimerCL time;
            VelVecDescCL curv( &Stokes.vel_idx);
            time.Reset();
            Stokes.SetupPrMass(  &Stokes.prM, lset/*, P.get<double>("Mat.ViscFluid"), C.mat_ViscGas*/);
            Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &curv, lset, Stokes.v.t);
            Stokes.SetupSystem2( &Stokes.B, &Stokes.C, &Stokes.c, lset, Stokes.v.t);
            curv.Clear( Stokes.v.t);
            lset.AccumulateBndIntegral( curv);
            time.Stop();
            std::cout << "Discretizing Stokes/Curv for initial velocities took "<<time.GetTime()<<" sec.\n";

            time.Reset();
            ssolver.Solve( Stokes.A.Data, Stokes.B.Data, Stokes.C.Data,
                Stokes.v.Data, Stokes.p.Data, Stokes.b.Data, Stokes.c.Data, Stokes.vel_idx.GetEx(), Stokes.pr_idx.GetEx());
            time.Stop();
            std::cout << "Solving Stokes for initial velocities took "<<time.GetTime()<<" sec.\n";
        }
        delete schurpc;
      } break;

      case  2: //Stationary flow without droplet
          Stokes.UpdateXNumbering( pidx, lset);
          Stokes.p.SetIdx( pidx);
      break;
      default : throw DROPSErrCL("Unknown initial condition");
    }
}

/// \brief Class for serializing a two-phase flow problem, i.e., storing
///    the multigrid and the numerical data
/// \todo  Storing of transport data!
template <typename StokesT>
class TwoPhaseStoreCL
{
  private:
    MultiGridCL&         mg_;
    const StokesT&       Stokes_;
    const LevelsetP2CL&  lset_;
    const TransportP1CL* transp_;//old
    const SpaceTimeXSolutionCL* st_transp_;//new
    std::string          path_;
    std::string          mgPath_;
    Uint                 numRecoverySteps_;
    Uint                 recoveryStep_;
    bool                 binary_;
    const PermutationT&  vel_downwind_;
    const PermutationT&  lset_downwind_;

    /// \brief Write time info
    void WriteTime( std::string filename)
    {
        std::ofstream file( filename.c_str());
        if (!file) throw DROPSErrCL("TwoPhaseStoreCL::WriteTime: Cannot open file "+filename+" for writing");
        file << Stokes_.v.t << "\n";
        file.close();
    }

  public:
      /// \brief Construct a class for storing a two-phase flow problem in files
      /** This class generates multiple files, all with prefix path, for storing
       *  the geometric as well as the numerical data.
       *  \param recoverySteps number of backup steps before overwriting files
       *  \param mg Multigrid
       *  \param Stokes Stokes flow field
       *  \param lset Level Set field
       *  \param transp mass transport concentration field
       *  \param path location for storing output
       *  \param binary save output  binary?
       *  */
    TwoPhaseStoreCL(MultiGridCL& mg, const StokesT& Stokes, const LevelsetP2CL& lset, const TransportP1CL* transp,
                    const std::string& solPath, const std::string& mgPath, Uint recoverySteps=2, bool binary= false, const PermutationT& vel_downwind= PermutationT(), const PermutationT& lset_downwind= PermutationT())
      : mg_(mg), Stokes_(Stokes), lset_(lset), transp_(transp), st_transp_(NULL), path_(solPath), mgPath_(mgPath), numRecoverySteps_(recoverySteps),
        recoveryStep_(0), binary_( binary), vel_downwind_( vel_downwind), lset_downwind_( lset_downwind){}

//new cstr with spacetimetransport-object (as reference to distinguish from std-cstr... - could be cleaned..)
    TwoPhaseStoreCL(MultiGridCL& mg, const StokesT& Stokes, const LevelsetP2CL& lset, const SpaceTimeXSolutionCL& st_transp,
                    const std::string& solPath, const std::string& mgPath, Uint recoverySteps=2, bool binary= false,
                    const PermutationT& vel_downwind= PermutationT(), const PermutationT& lset_downwind= PermutationT())
        : mg_(mg), Stokes_(Stokes), lset_(lset), transp_(NULL), st_transp_(&st_transp), path_(solPath), mgPath_(mgPath), numRecoverySteps_(recoverySteps),
        recoveryStep_(0), binary_( binary), vel_downwind_( vel_downwind), lset_downwind_( lset_downwind){}

    /// \brief Write all information in a file
    void Write()
    {
        if( path_.empty() && mgPath_.empty() ) return;

        // Create filename
        std::stringstream filename;
        std::stringstream mgFileName;
        const size_t postfix= numRecoverySteps_==0 ? recoveryStep_++ : (recoveryStep_++)%numRecoverySteps_;
        if( path_.empty() )
            filename << mgPath_ << postfix;
        else
            filename << path_ << postfix;
        mgFileName << mgPath_ << postfix;

        // first master writes time info
        IF_MASTER
            WriteTime( filename.str() + "time");

        // write multigrid
        if( !mgPath_.empty() )
        {
            MGSerializationCL ser( mg_, mgFileName.str() );
            ser.WriteMG();
        }

        // write numerical data if path is specified
        if( path_.empty() ) return;
        VecDescCL vel= Stokes_.v;
        permute_Vector( vel.Data, invert_permutation( vel_downwind_), 3);
        WriteFEToFile( vel, mg_, filename.str() + "velocity", binary_);
        VecDescCL ls= lset_.Phi;
        permute_Vector( ls.Data, invert_permutation( lset_downwind_));
        WriteFEToFile( ls, mg_, filename.str() + "levelset", binary_);
        WriteFEToFile(Stokes_.p, mg_, filename.str() + "pressure", binary_, &lset_.Phi); // pass also level set, as p may be extended
        if (transp_) WriteFEToFile(transp_->ct, mg_, filename.str() + "concentrationTransf", binary_);
        if (st_transp_) WriteFEToFile(st_transp_->GetFutureTrace_Neg(), mg_, filename.str() + "concentration_neg", binary_);
        if (st_transp_) WriteFEToFile(st_transp_->GetFutureTrace_Pos(), mg_, filename.str() + "concentration_pos", binary_);
    }
};

}   // end of namespace DROPS

namespace SurfTens
{
//==============================================================================
//          Functions for twophasedrops-executable (surface tension coefficient)
//==============================================================================
// tau is constant
double sigmaf (const DROPS::Point3DCL&, double) {
    static double sigma = P.get<double>("NavStokes.Coeff.SurfTens.SurfTension");
    return sigma;
}

DROPS::Point3DCL gsigma (const DROPS::Point3DCL&, double) { return DROPS::Point3DCL(); }

//linear decrease of surface tension coefficient in y direction
double lin_in_y(const DROPS::Point3DCL& p)  // linear function in y-direction, zero in the middle of the domain
{
    static double sigma = P.get<double>("NavStokes.Coeff.SurfTens.SurfTension",1.0); // surface tension coefficient sigma : if tau is a variable, sigma is the constant part of tau.
    static double Ly = P.get<DROPS::Point3DCL>("Mesh.E2")[1];  // domain size in y
    static double grad_tau = P.get<double>("SurfTens.GradTau", 0.);  // gradient of tau

    return sigma + (p[1] - Ly*0.5)*grad_tau;  // grad_tau is the gradient of tau in y-direction
}

double sigma_const_grad_tau_y(const DROPS::Point3DCL& p, double) { return lin_in_y( p); }

// Registration of functions
static DROPS::RegisterScalarFunction regconstsurftens("ConstTau", sigmaf);
static DROPS::RegisterScalarFunction regconstgradtau_y("ConstGradTau_y", sigma_const_grad_tau_y);

}// end of namespace SurfTens

#endif

