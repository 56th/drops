/// \file instatstokes2phase.h
/// \brief classes that constitute the 2-phase Stokes problem
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt, Thomas Ludescher, Liang Zhang; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_INSTATSTOKES2PHASE_H
#define DROPS_INSTATSTOKES2PHASE_H

#include <memory>

#include "stokes/stokes.h"
#include "levelset/levelset.h"
#include "levelset/mgobserve.h"
#include "misc/params.h"
#include "num/MGsolver.h"
#include "num/fe_repair.h"
#include "misc/funcmap.h"
extern DROPS::ParamCL P;
namespace DROPS
{

/// \brief Repair a P1X-vector if grid changes occur
///
/// Create such an object with the variable to be repaired before any grid-modifications.
/// Repair the linear part however you like and call the operator() to repair the extended part.
class P1XRepairCL
{
  private:
    bool UsesXFEM_;
    MultiGridCL& mg_;
    IdxDescCL idx_;
    VecDescCL ext_;
    VecDescCL& p_;

  public:
    P1XRepairCL( MultiGridCL& mg, VecDescCL& p);

    VecDescCL* GetExt() { return &ext_; }

    void operator() ();
};

/// \brief Compute the main diagonal of the unscaled \f$L_2(\Omega)\f$-mass-matrix.
void SetupMassDiag_P1 (const MultiGridCL& MG, VectorCL& M, const IdxDescCL& RowIdx,
                       const BndCondCL& bnd= BndCondCL( 0));

/// \brief Compute the main diagonal of the unscaled \f$L_2(\Omega)\f$-mass-matrix.
void SetupMassDiag_P1X (const MultiGridCL& MG, VectorCL& M, const IdxDescCL& RowIdx, const VecDescCL& lset,
                       const BndDataCL<>& lsetbnd, const BndCondCL& bnd= BndCondCL( 0));

/// \brief Compute the main diagonal of the unscaled \f$L_2(\Omega)\f$-mass-matrix.
void SetupMassDiag_vecP2 (const MultiGridCL& MG, VectorCL& M, const IdxDescCL& RowIdx,
                          const BndCondCL& bnd= BndCondCL( 0));
/// \brief Compute the main diagonal of the unscaled \f$L_2(\Omega)\f$-mass-matrix.
void SetupMassDiag (const MultiGridCL& MG, VectorCL& M, const IdxDescCL& RowIdx,
                    const BndCondCL& bnd= BndCondCL( 0), const VecDescCL* lsetp=0, const BndDataCL<>* lsetbnd=0);

/// \brief Compute the unscaled lumped \f$L_2(\Omega)\f$-mass-matrix, i.e., M = diag( \f$\int_\Omega v_i dx\f$).
void SetupLumpedMass (const MultiGridCL& MG, VectorCL& M, const IdxDescCL& RowIdx,
                    const BndCondCL& bnd= BndCondCL( 0), const VecDescCL* lsetp=0, const BndDataCL<>* lsetbnd=0);


// rho*du/dt - mu*laplace u + Dp = f + rho*g - okn
//                        -div u = 0
//                             u = u0, t=t0

/// \brief Parameter class describing a twophase flow
class TwoPhaseFlowCoeffCL
{
// \Omega_1 = Tropfen,    \Omega_2 = umgebendes Fluid

  private:
    bool film; // TL: delete?
    bool ns_shiftframe;
    double surfTens;
    double rho_koeff1, rho_koeff2, mu_koeff1, mu_koeff2;
    double beta_coeff1, beta_coeff2; ///< slip bnd coeff for fluid 1/2

  public:
    DROPS::instat_vector_fun_ptr volforce;
    DROPS::instat_vector_fun_ptr BndOutNormal;
    const SmoothedJumpCL rho, mu;
    const SmoothedJumpCL beta;      ///< Coefficient for SlipBC
    const double SurfTens, DilVisco, ShearVisco;
    const double betaL,             ///< contact line coefficient
                 alpha;             ///< Nitsche coeff for slip bnd
    const Point3DCL g;
    const Point3DCL framevel;
    DROPS::instat_scalar_fun_ptr var_tau_fncs;

    TwoPhaseFlowCoeffCL( ParamCL& P, bool dimless = false)
      //big question: film or measurecell? 1: measure, 2: film
            //If we merge film.cpp and twophasedrops.cpp, we don't need film flag anymore.
        : film( false ), /// \todo change to something meaningful
        ns_shiftframe( (P.get<int>("NavStokes.ShiftFrame", 0) == 1) ),
        surfTens( P.get<double>("NavStokes.Coeff.SurfTens.SurfTension") ),
        rho_koeff1( P.get<double>("NavStokes.Coeff.DensPos") ),
        rho_koeff2( P.get<double>("NavStokes.Coeff.DensNeg") ),
        mu_koeff1( P.get<double>("NavStokes.Coeff.ViscPos") ),
        mu_koeff2( P.get<double>("NavStokes.Coeff.ViscNeg") ),
        beta_coeff1(P.get<double>("NavStokes.BoundaryData.SlipBnd.Beta1")),
        beta_coeff2(P.get<double>("NavStokes.BoundaryData.SlipBnd.Beta2")),

        rho( dimless ? JumpCL( 1., rho_koeff1/rho_koeff2)
          : JumpCL( rho_koeff2, rho_koeff1), H_sm, P.get<double>("NavStokes.Coeff.SmoothZone")),
        mu( dimless ? JumpCL( 1., mu_koeff1/mu_koeff2)
          : JumpCL( mu_koeff2, mu_koeff1), H_sm, P.get<double>("NavStokes.Coeff.SmoothZone")),
        beta(dimless ? JumpCL( 1., beta_coeff2/beta_coeff1)
                     : JumpCL(beta_coeff2,beta_coeff1), H_sm, 0),
        SurfTens (dimless ? surfTens/rho_koeff2 : surfTens),
        DilVisco( P.get<double>("NavStokes.Coeff.SurfTens.DilatationalVisco") ),
        ShearVisco( P.get<double>("NavStokes.Coeff.SurfTens.ShearVisco") ),
        betaL(P.get<double>("NavStokes.BoundaryData.SlipBnd.BetaL")), alpha(P.get<double>("NavStokes.BoundaryData.SlipBnd.NitschePenalty")),
        g( P.get<DROPS::Point3DCL>("NavStokes.Coeff.Gravity")),
        framevel( ns_shiftframe ? P.get<DROPS::Point3DCL>("NavStokes.FrameVel", DROPS::Point3DCL(0.0)) : DROPS::Point3DCL(0.0) )
        {
        volforce = InVecMap::getInstance()[P.get<std::string>("NavStokes.Coeff.VolForce")];
        var_tau_fncs = InScaMap::getInstance()[P.get<std::string>("NavStokes.Coeff.SurfTens.VarTensionFunc")];
        if( P.get<std::string>("NavStokes.BoundaryData.SlipBnd.BndOuterNormal").compare("None")!=0)
            BndOutNormal = InVecMap::getInstance()[P.get<std::string>("NavStokes.BoundaryData.SlipBnd.BndOuterNormal")];
        else
            BndOutNormal = nullptr;
    }

    TwoPhaseFlowCoeffCL( double rho1, double rho2, double mu1, double mu2, double surftension, Point3DCL gravity, Point3DCL framevelocity = Point3DCL(0.0), bool dimless = false, double dilatationalvisco = 0.0, double shearvisco = 0.0, 
                         double betaL_=0, double alpha_ = 1.0, double beta1 = 0.0, double beta2 =0.0)
      : rho( dimless ? JumpCL( 1., rho2/rho1)
                     : JumpCL( rho1, rho2), H_sm, 0),
        mu(  dimless ? JumpCL( 1., mu2/mu1)
                     : JumpCL( mu1, mu2), H_sm, 0),
        beta( dimless? JumpCL( 1., beta2/beta1)
                     : JumpCL(beta2,beta1), H_sm, 0),
        SurfTens( dimless ? surftension/rho1 : surftension),
        DilVisco( dilatationalvisco),
        ShearVisco( shearvisco),
        betaL(betaL_),
        alpha(alpha_), 
        g( gravity),
        framevel( framevelocity)
    {
        volforce = InVecMap::getInstance()["ZeroVel"];
        var_tau_fncs = InScaMap::getInstance()["ConstTau"];
        BndOutNormal = InVecMap::getInstance()["ZeroVel"];
    }
};

/// problem class for instationary two-pase Stokes flow
class InstatStokes2PhaseP2P1CL : public ProblemCL<TwoPhaseFlowCoeffCL, StokesBndDataCL>
{
private:
    void computeGhostPenaltyKernel(MultiGridCL& mg, const LevelsetP2CL& lset , const IdxDescCL &prIdx) const;
    double epsP;           ///< ghost penalty stabilization parameter

  public:
    typedef ProblemCL<TwoPhaseFlowCoeffCL, StokesBndDataCL>       base_;
    typedef base_::BndDataCL                                      BndDataCL;
    using base_::MG_;
    using base_::Coeff_;
    using base_::BndData_;
    using base_::GetBndData;
    using base_::GetMG;

    typedef P1EvalCL<double, const StokesPrBndDataCL, VecDescCL>   DiscPrSolCL;
    typedef P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, VelVecDescCL> DiscVelSolCL;
    typedef P1EvalCL<double, const StokesPrBndDataCL, const VecDescCL>   const_DiscPrSolCL;
    typedef P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, const VelVecDescCL> const_DiscVelSolCL;

  public:
    MLIdxDescCL  vel_idx;  ///< for velocity unknowns
    MLIdxDescCL  pr_idx;   ///< for pressure unknowns
    VelVecDescCL v;        ///< velocity
    VecDescCL    p;        ///< pressure
    VelVecDescCL b;
    VecDescCL    c;
    MLMatDescCL  A,
                 B,
                 C,
                 M,
                 prA,
                 prM,
                 prMhat;
    mutable VectorBaseCL<VectorCL> cKernel;
    SurfaceTensionCL*     SurfTension_;
    instat_scalar_fun_ptr CtAngleFnc_;
    instat_vector_fun_ptr BndOutNormal_;

  public:
    InstatStokes2PhaseP2P1CL( const MGBuilderCL& mgb, const TwoPhaseFlowCoeffCL& coeff, const BndDataCL& bdata, FiniteElementT prFE= P1_FE, double XFEMstab=0.1, FiniteElementT velFE= vecP2_FE, double EpsP = 0.0 )
        : base_(mgb, coeff, bdata), epsP(EpsP), vel_idx(velFE, 1, bdata.Vel, XFEMstab), pr_idx(prFE, 1, bdata.Pr, XFEMstab), cKernel(0), SurfTension_(nullptr) { }
    InstatStokes2PhaseP2P1CL( MultiGridCL& mg, const TwoPhaseFlowCoeffCL& coeff, const BndDataCL& bdata, FiniteElementT prFE= P1_FE, double XFEMstab=0.1, FiniteElementT velFE= vecP2_FE, double EpsP = 0.0)
        : base_(mg, coeff, bdata), epsP(EpsP),  vel_idx(velFE, 1, bdata.Vel, XFEMstab), pr_idx(prFE, 1, bdata.Pr, XFEMstab), cKernel(0), SurfTension_(nullptr) { }

    /// \name Numbering
    //@{
    /// Create/delete numbering of unknowns
    void CreateNumberingVel( Uint level, MLIdxDescCL* idx, const LevelsetP2CL* lsetp= 0);
    void CreateNumberingPr ( Uint level, MLIdxDescCL* idx, const LevelsetP2CL* lsetp= 0);
    /// \brief Only used for XFEM
    void UpdateXNumbering( MLIdxDescCL* idx, const LevelsetP2CL& lset)
        {
            if (UsesXFEM()) idx->UpdateXNumbering( MG_, *lset.PhiC, lset.GetBndData());
        }
    /// \brief Only used for XFEM
    void UpdatePressure( VecDescCL* p)
        {
            if (UsesXFEM()) p->RowIdx->GetXidx().Old2New( p);
        }
    void DeleteNumbering( MLIdxDescCL* idx)
        { idx->DeleteNumbering( MG_); }
    //@}

    /// \name Discretization
    //@{
    /// Returns whether extended FEM are used for pressure
    bool UsesXFEM() const { return pr_idx.GetFinest().IsExtended(); }
    /// Set up matrices A, M and rhs b (depending on phase bnd)
    void SetupSystem1( MLMatDescCL* A, MLMatDescCL* M, VecDescCL* b, VecDescCL* cplA, VecDescCL* cplM, const LevelsetP2CL& lset, double t) const;
    MLTetraAccumulatorTupleCL& system1_accu (MLTetraAccumulatorTupleCL& accus, MLMatDescCL* A, MLMatDescCL* M, VecDescCL* b, VecDescCL* cplA, VecDescCL* cplM, const LevelsetP2CL& lset, double t) const;
    /// Set up rhs b (depending on phase bnd)
    void SetupRhs1( VecDescCL* b, const LevelsetP2CL& lset, double t) const;
    /// Set up coupling terms for M matrix at given time t for time integration
    void SetupCplM( VecDescCL *cplM, const LevelsetP2CL &lset, double t) const;
    /// Set up matrix vector product v = A^n*u^n - cplA at given time t
    void SetupAdotU( VecDescCL *AdotU, const VecDescCL &un, const LevelsetP2CL &lset, double t) const;
    /// Set up the Laplace-Beltrami-Operator
    void SetupLB( MLMatDescCL* A, VecDescCL* cplA, const LevelsetP2CL& lset, double t) const;
    /// Set up the Boussinesq-Scriven Law of surface stress
    void SetupBS( MLMatDescCL* A, VecDescCL* cplA, const LevelsetP2CL& lset, double t) const;
    /// Set up matrix B and rhs c
    void SetupSystem2( MLMatDescCL* B, MLMatDescCL *C, VecDescCL* c, const LevelsetP2CL& lset, double t) const;
    MLTetraAccumulatorTupleCL& system2_accu (MLTetraAccumulatorTupleCL& accus, MLMatDescCL* B, VecDescCL* c, const LevelsetP2CL& lset, double t) const;
    /// Set up rhs c
    void SetupRhs2( VecDescCL* c, const LevelsetP2CL& lset, double t) const;
    /// Set up the time-derivative of B times velocity
    void SetupBdotv (VecDescCL* Bdotv, const VelVecDescCL* vel, const LevelsetP2CL& lset, double t) const;
    /// Set up the mass matrix for the pressure, scaled by \f$\mu^{-1}\f$.
    void SetupPrMass( MLMatDescCL* prM, const LevelsetP2CL& lset) const;
    /// Set up the overlapping mass matrix for the pressure, scaled by \f$\mu^{-1}\f$.
    void SetupPrMassHat( MLMatDescCL* prMhat, const LevelsetP2CL& lset) const;
    /// Set up the stabilisation matrix for the pressure.
    void SetupC( MLMatDescCL* matC, const LevelsetP2CL& lset, double eps_p ) const;
    /// Set up the stiffness matrix for the pressure, scaled by \f$\rho^{-1}\f$.
    void SetupPrStiff(MLMatDescCL* prA, const LevelsetP2CL& lset, double lambda=1.0) const;
    //@}

    /// Initialize velocity field
    void InitVel( VelVecDescCL*, instat_vector_fun_ptr, double t0= 0.) const;
    /// Smooth velocity field
    void SmoothVel( VelVecDescCL*, int num= 1, double tau=0.5);
    /// Clear all matrices, should be called after grid change to avoid reuse of matrix pattern
    void ClearMat() { A.Data.clear(); B.Data.clear(); C.Data.clear(); M.Data.clear(); prA.Data.clear(); prM.Data.clear(); }
    /// Set all indices
    void SetIdx();
    /// Set number of used levels
    void SetNumVelLvl( size_t n);
    void SetNumPrLvl ( size_t n);
    /// Get FE type for velocity space
    FiniteElementT GetVelFE() const { return vel_idx.GetFinest().GetFE(); }
    /// Get FE type for pressure space
    FiniteElementT GetPrFE() const { return pr_idx.GetFinest().GetFE(); }
    /// \name Get extended index (only makes sense for P1X_FE)
    //@{
    const ExtIdxDescCL& GetXidx() const { return pr_idx.GetFinest().GetXidx(); }
    ExtIdxDescCL&       GetXidx()       { return pr_idx.GetFinest().GetXidx(); }
    //@}
    /// Get pressure solution on inner/outer part (especially for P1X_FE)
    void GetPrOnPart( VecDescCL& p_part, const LevelsetP2CL& lset, bool posPart= true); // false = inner = Phi<0, true = outer = Phi>0
    /// Get CFL restriction for explicit time stepping
    double GetCFLTimeRestriction( LevelsetP2CL& lset);
    /// check whether Ghost Penalty is used or not
    bool usesGhostPen(){ return epsP > 0.0; }
    /// get Ghost Penalty stabilization factor
    double getGhPenStab(){ return epsP; }
    /// set Ghost Penalty stabilization factor
    void setGhPenStab( double EpsP ){ epsP = EpsP; }
    
    /// Set Equilibrium Contact Angle
    void SetYoungAngle(instat_scalar_fun_ptr CtAngleFnc) { CtAngleFnc_= CtAngleFnc; }
    /// Set out normal function of the slip boundary
    void SetBndOutNormal(instat_vector_fun_ptr outnormal) { BndOutNormal_= outnormal; }
    /// Set the surface force type and the surface tension
    void SetSurfTension(SurfaceTensionCL* Sf) { SurfTension_= Sf; }
    /// Discretize Young Force on the three-phase contact line
    void AccumulateYoungForce(const LevelsetP2CL& lset, VecDescCL& f) const;

    /// \name Evaluate Solution
    //@{
    /// Get solution as FE-function for evaluation
    const_DiscPrSolCL GetPrSolution() const
        { return const_DiscPrSolCL( &p, &GetBndData().Pr, &GetMG()); }
    const_DiscVelSolCL GetVelSolution() const
        { return const_DiscVelSolCL( &v, &GetBndData().Vel, &GetMG()); }

    const_DiscPrSolCL GetPrSolution( const VecDescCL& pr) const
        { return const_DiscPrSolCL( &pr, &GetBndData().Pr, &GetMG()); }
    const_DiscVelSolCL GetVelSolution( const VelVecDescCL& vel) const
        { return const_DiscVelSolCL( &vel, &GetBndData().Vel, &GetMG()); }
    //@}
};

/// \brief Observes the MultiGridCL-changes by AdapTriangCL to repair the Function stokes_.v.
///
/// The actual work is done in post_refine().
class VelocityRepairCL : public MGObserverCL
{
  private:
    InstatStokes2PhaseP2P1CL& stokes_;
    std::unique_ptr<RepairP2CL<Point3DCL>::type > p2repair_;

  public:
    VelocityRepairCL (InstatStokes2PhaseP2P1CL& stokes)
        : stokes_( stokes) {}
    void pre_refine  ();
    void post_refine ();
    void pre_refine_sequence  () {}
    void post_refine_sequence () {}
    const IdxDescCL* GetIdxDesc() const { return stokes_.v.RowIdx; }
#ifdef _PAR
    const VectorCL*  GetVector()  const { return &stokes_.v.Data; }
    void swap( IdxDescCL& idx, VectorCL& v) { stokes_.v.RowIdx->swap(idx); stokes_.v.Data.swap(v); }
#endif
};

/// \brief Observes the MultiGridCL-changes by AdapTriangCL to repair the Function stokes_.pr.
///
/// For the P1-part, the actual work is done in post_refine().
/// For the P1X-part, a P1XRepairCL is created in pre_refine_sequence() and used in
/// post_refine_sequence(). Holding the P1XRepairCL* in an auto_ptr simplifies the use
/// of heap-memory: No memory is lost, even if successive calls of pre_refine_sequence()
/// occur without interleaved post_refine_sequence()-calls.
class PressureRepairCL : public MGObserverCL
{
  private:
    InstatStokes2PhaseP2P1CL& stokes_;
    std::unique_ptr<P1XRepairCL> p1xrepair_;
    std::unique_ptr<RepairP1CL<double>::type> p1repair_;
    const LevelsetP2CL& ls_;

  public:
    PressureRepairCL ( InstatStokes2PhaseP2P1CL& stokes, const LevelsetP2CL& ls)
        : stokes_( stokes), ls_( ls) {}
    void pre_refine  ();
    void post_refine ();
    void pre_refine_sequence  ();
    void post_refine_sequence ();
    const IdxDescCL* GetIdxDesc() const { return stokes_.p.RowIdx; }
#ifdef _PAR
    const VectorCL*  GetVector()  const { return &stokes_.p.Data; }
    void swap( IdxDescCL& idx, VectorCL& v) { stokes_.p.RowIdx->swap(idx); stokes_.p.Data.swap(v); }
#endif
};

} // end of namespace DROPS

#include "stokes/instatstokes2phase.tpp"

#endif

