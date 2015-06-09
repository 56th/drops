/// \file instatstokes2phase.h
/// \brief classes that constitute the 2-phase Stokes problem
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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
    bool film;
    double surfTens;
    double rho_koeff1, rho_koeff2, mu_koeff1, mu_koeff2;
    double beta_coeff1, beta_coeff2;

  public:
    DROPS::instat_vector_fun_ptr volforce;
    DROPS::instat_vector_fun_ptr RefVel;
    DROPS::instat_vector_fun_ptr RefGradPr;
    DROPS::instat_scalar_fun_ptr RefPr;
    DROPS::instat_vector_fun_ptr Bndoutnormal;
    const SmoothedJumpCL rho, mu;
    const SmoothedJumpCL beta;    //slip length
    const double SurfTens, DilVisco, ShearVisco;
	const double alpha;
    const Point3DCL g;

    TwoPhaseFlowCoeffCL( ParamCL& P, bool dimless = false)
      //big question: film or measurecell? 1: measure, 2: film
        : film( (P.get<double>("Mat.DensDrop") == 0.0) ),
        surfTens( film ? P.get<double>("Mat.SurfTension") : P.get<double>("SurfTens.SurfTension")),
        rho_koeff1( film ? P.get<double>("Mat.DensGas") : P.get<double>("Mat.DensFluid")),
        rho_koeff2( film ? P.get<double>("Mat.DensFluid") : P.get<double>("Mat.DensDrop")),
        mu_koeff1( film ? P.get<double>("Mat.ViscGas") : P.get<double>("Mat.ViscFluid")),
        mu_koeff2( film ? P.get<double>("Mat.ViscFluid") : P.get<double>("Mat.ViscDrop")),
        beta_coeff1(P.get<double>("SpeBnd.beta1")),
        beta_coeff2(P.get<double>("SpeBnd.beta2")),
		
        rho( dimless ? JumpCL( 1., rho_koeff1/rho_koeff2)
          : JumpCL( rho_koeff2, rho_koeff1), H_sm, P.get<double>("Mat.SmoothZone")),
        mu( dimless ? JumpCL( 1., mu_koeff1/mu_koeff2)
          : JumpCL( mu_koeff2, mu_koeff1), H_sm, P.get<double>("Mat.SmoothZone")),
        beta(dimless ? JumpCL( 1., beta_coeff2/beta_coeff1)
                     :JumpCL(beta_coeff2,beta_coeff1), H_sm, P.get<double>("SpeBnd.SmoothZone")),
        SurfTens (dimless ? surfTens/rho_koeff2 : surfTens),
        DilVisco( film ? P.get<double>("Mat.DilatationalVisco") : P.get<double>("SurfTens.DilatationalVisco")),
        ShearVisco( film ? P.get<double>("Mat.ShearVisco") : P.get<double>("SurfTens.ShearVisco")),
		alpha(P.get<double>("SpeBnd.alpha")),
        g( P.get<DROPS::Point3DCL>("Exp.Gravity"))
        {
			volforce = InVecMap::getInstance()[P.get<std::string>("Exp.VolForce")];
			if( P.get<std::string>("Exp.Solution_Vel").compare("None")!=0)
				RefVel = InVecMap::getInstance()[P.get<std::string>("Exp.Solution_Vel")];
			else
				RefVel = NULL;
			if( P.get<std::string>("Exp.Solution_GradPr").compare("None")!=0)
				RefGradPr = InVecMap::getInstance()[P.get<std::string>("Exp.Solution_GradPr")];
			else
				RefGradPr = NULL;
			if( P.get<std::string>("Exp.Solution_Pr").compare("None")!=0)
				RefPr = InScaMap::getInstance()[P.get<std::string>("Exp.Solution_Pr")];
			else
				RefPr = NULL;
          if( P.get<std::string>("SpeBnd.BndOutNormal").compare("None")!=0)
				Bndoutnormal = InVecMap::getInstance()[P.get<std::string>("SpeBnd.BndOutNormal")];
			else
				Bndoutnormal = NULL;
    }

    TwoPhaseFlowCoeffCL( double rho1, double rho2, double mu1, double mu2, double surftension, Point3DCL gravity, bool dimless = false, double dilatationalvisco = 0.0, double shearvisco = 0.0, double alpha_ = 1.0, double beta1 = 0.0, double beta2 =0.0
    							)//double sl1=1, double sl2=1
      : rho( dimless ? JumpCL( 1., rho2/rho1)
                     : JumpCL( rho1, rho2), H_sm, 0),
        mu(  dimless ? JumpCL( 1., mu2/mu1)
                     : JumpCL( mu1, mu2), H_sm, 0),
        beta( dimless ? JumpCL( 1., beta2/beta1)
                     :JumpCL(beta2,beta1), H_sm, 0),
        SurfTens( dimless ? surftension/rho1 : surftension),
        DilVisco( dilatationalvisco),
        ShearVisco( shearvisco),
		alpha(alpha_), 
        g( gravity)
		{
          volforce   = InVecMap::getInstance()["ZeroVel"];
		  RefVel     = InVecMap::getInstance()["ZeroVel"];
		  RefGradPr  = InVecMap::getInstance()["ZeroVel"];
		  RefPr      = InScaMap::getInstance()["Zero"];
          Bndoutnormal = InVecMap::getInstance()["ZeroVel"];
        }
};

/// problem class for instationary two-pase Stokes flow


class InstatStokes2PhaseP2P1CL : public ProblemCL<TwoPhaseFlowCoeffCL, StokesBndDataCL>
{
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
                 prM;
    //VelVecDescCL dv;        ///< velocity error, not used now
  
  public:
    InstatStokes2PhaseP2P1CL( const MGBuilderCL& mgb, const TwoPhaseFlowCoeffCL& coeff, const BndDataCL& bdata, FiniteElementT prFE= P1_FE, double XFEMstab=0.1, FiniteElementT velFE= vecP2_FE)
        : base_(mgb, coeff, bdata), vel_idx(velFE, 1, bdata.Vel, 0, XFEMstab), pr_idx(prFE, 1, bdata.Pr, 0, XFEMstab) {}
    InstatStokes2PhaseP2P1CL( MultiGridCL& mg, const TwoPhaseFlowCoeffCL& coeff, const BndDataCL& bdata, FiniteElementT prFE= P1_FE, double XFEMstab=0.1, FiniteElementT velFE= vecP2_FE)
        : base_(mg, coeff, bdata),  vel_idx(velFE, 1, bdata.Vel, 0, XFEMstab), pr_idx(prFE, 1, bdata.Pr, 0, XFEMstab) {}

    /// \name Numbering
    //@{
    /// Create/delete numbering of unknowns
    void CreateNumberingVel( Uint level, MLIdxDescCL* idx, match_fun match= 0, const LevelsetP2CL* lsetp= 0);
    void CreateNumberingPr ( Uint level, MLIdxDescCL* idx, match_fun match= 0, const LevelsetP2CL* lsetp= 0);
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
    /// Set up the Laplace-Beltrami-Operator
    void SetupLB( MLMatDescCL* A, VecDescCL* cplA, const LevelsetP2CL& lset, double t) const;
    /// Set up the Boussinesq-Scriven Law of surface stress
    void SetupBS( MLMatDescCL* A, VecDescCL* cplA, const LevelsetP2CL& lset, double t) const;
    /// Set up matrix B and rhs c
    void SetupSystem2( MLMatDescCL* B, VecDescCL* c, const LevelsetP2CL& lset, double t) const;
    MLTetraAccumulatorTupleCL& system2_accu (MLTetraAccumulatorTupleCL& accus, MLMatDescCL* B, VecDescCL* c, const LevelsetP2CL& lset, double t) const;
    /// Set up rhs c
    void SetupRhs2( VecDescCL* c, const LevelsetP2CL& lset, double t) const;
    /// Set up the time-derivative of B times velocity
    void SetupBdotv (VecDescCL* Bdotv, const VelVecDescCL* vel, const LevelsetP2CL& lset, double t) const;
    /// Set up the mass matrix for the pressure, scaled by \f$\mu^{-1}\f$.
    void SetupPrMass( MLMatDescCL* prM, const LevelsetP2CL& lset) const;
    /// Set up the stabilisation matrix for the pressure.
    void SetupC( MLMatDescCL* matC, const LevelsetP2CL& lset, double eps_p ) const;
    /// Set up the stiffness matrix for the pressure, scaled by \f$\rho^{-1}\f$.
    void SetupPrStiff(MLMatDescCL* prA, const LevelsetP2CL& lset) const;
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

	//This function can only be applied for one phase simulation, because the function doesn't consider any discontinuity in solutions;
	//It checks the L2 norm of discretized error of velocity and gradient of pressure
    void CheckOnePhaseSolution(const VelVecDescCL* DescVel, const VecDescCL* DescPr, const instat_vector_fun_ptr RefVel, const instat_vector_fun_ptr RefGradPr , const instat_scalar_fun_ptr RefPr) const;
    void CheckTwoPhaseSolution(const VelVecDescCL* DescVel, const VecDescCL* DescPr, const LevelsetP2CL& lset, const VelVecDescCL* RefVel, const VecDescCL* RefPr);
    void CheckTwoPhaseSolution(const VelVecDescCL* DescVel, const VecDescCL* DescPr, const LevelsetP2CL& lset, const instat_vector_fun_ptr RefVel, const instat_scalar_fun_ptr RefPr);

    //Get the total kinetic energy
    double GetKineticEnergy(const LevelsetP2CL& lset) const;
	// To setup u-u_h
	//void SetupVelError(const instat_vector_fun_ptr RefVel);
    /// \name Evaluate Solution
    //@{
    /// Get solution as FE-function for evaluation
    const_DiscPrSolCL GetPrSolution() const
        { return const_DiscPrSolCL( &p, &GetBndData().Pr, &GetMG()); }
    const_DiscVelSolCL GetVelSolution() const
        { return const_DiscVelSolCL( &v, &GetBndData().Vel, &GetMG()); }
	// to get error in velocity	
    //const_DiscVelSolCL GetVelError(StokesBndDataCL::VelBndDataCL& bnd) const
    //    { return const_DiscVelSolCL( &dv, &bnd, &GetMG()); }		

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
    std::auto_ptr<RepairP2CL<Point3DCL>::type > p2repair_;

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
    std::auto_ptr<P1XRepairCL> p1xrepair_;
    std::auto_ptr<RepairP1CL<double>::type> p1repair_;
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


//\brief Handle the slip or symmetric boundary condition
class SpecialBndHandleSystem2_P2P1XCL
{
 private:
    BndTriangPartitionCL 	  bndpartition_;
    QuadDomainCL              bndq5dom_;

	const StokesBndDataCL& BndData_;
	Point3DCL normal;
	Uint unknownIdx[6];
  public:
	SpecialBndHandleSystem2_P2P1XCL(const StokesBndDataCL& BndData): BndData_(BndData) {}
	void setupB(SMatrixCL<1, 3> loc_b[4][10], const TetraCL& tet, int ls_sign[4], const PrincipalLatticeCL& lat,
                     const std::valarray<double>& ls_loc_, IdxT prNumb[4], const ExtIdxDescCL* Xidx);
};

} // end of namespace DROPS

#include "stokes/instatstokes2phase.tpp"

#endif

