/// \file levelset.h
/// \brief levelset equation for two phase flow problems
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

#ifndef DROPS_LEVELSET_H
#define DROPS_LEVELSET_H

#include "num/spmat.h"
#include "num/discretize.h"
#include "num/bndData.h"
#include "num/fe.h"
#include "num/fe_repair.h"
#include "levelset/mgobserve.h"
#include "levelset/surfacetension.h"
#include "levelset/volume_adjustment.h"
#include "num/interfacePatch.h"
#include "num/renumber.h"
#include "num/prolongation.h"
#include <vector>

#ifdef _PAR
# include "parallel/exchange.h"
# include "misc/container.h"
#endif

namespace DROPS
{

typedef BndDataCL<double>    LsetBndDataCL;

enum SurfaceForceT
/// different types of surface forces
{
    SF_LB=0,             ///< Laplace-Beltrami discretization: \f$\mathcal O(h^{1/2})\f$
    SF_Const=2,          ///< surface force with constant curvature
    SF_ImprovedLBVar=3   ///< improved Laplace-Beltrami discretization with variable surface tension
};

/// not used at the moment
class LevelsetCoeffCL {

};

/// forward declaration
class VolumeAdjustmentCL;

class LevelsetP2CL : public ProblemCL< LevelsetCoeffCL, LsetBndDataCL>
/// abstract base class for continuous and discontinuous P2 level set discretization
/// P2-discretization and solution of the level set equation for two phase flow problems. Bnd_ will be used to impose boundary data on the inflow boundary.
/// At the moment setting all boundary conditions to NoBC is the only valid case.
{
  public:
    typedef P2EvalCL<double, const LsetBndDataCL, VecDescCL>       DiscSolCL;
    typedef P2EvalCL<double, const LsetBndDataCL, const VecDescCL> const_DiscSolCL;
    typedef std::vector<Point3DCL> perDirSetT;
    MLIdxDescCL  idx;
    MLIdxDescCL* idxC;
    MLVecDescCL  MLPhi;     ///< level set function on all level (needed in parallel for multigrid)
    VecDescCL    Phi;       ///< level set function
    VecDescCL *  PhiC;      ///< continuous version of level set function

    typedef ProblemCL<LevelsetCoeffCL, LsetBndDataCL> base_;
    using base_::MG_;
    using base_::Coeff_;
    using base_::BndData_;
    using base_::GetBndData;
    using base_::GetMG;

  protected:
    double              curvDiff_, ///< amount of diffusion in curvature calculation
                        SD_;       ///< streamline diffusion
    SurfaceForceT       SF_;

    SurfaceTensionCL&   sf_;      ///< data for surface tension
    void SetupSmoothSystem ( MatrixCL&, MatrixCL&)               const;
    void SmoothPhi( VectorCL& SmPhi, double diff)                const;
    double GetVolume_Composite( double translation, int l)       const;
    double GetVolume_Extrapolation( double translation, int l)   const;
    perDirSetT* perDirections;    ///< periodic directions
    typedef MLDataCL<ProlongationCL<double> > ProlongationT;
    ProlongationT P_;
    mutable std::unique_ptr<VolumeAdjustmentCL> volume_adjuster_;

    bool IsDG;

  public:
    MatrixCL            E, H;  ///< E: mass matrix, H: convection matrix
    VecDescCL           rhs;  ///< rhs due to boundary conditions

    bool IsDiscontinuous(){ return IsDG; }

LevelsetP2CL( MultiGridCL& mg, const LsetBndDataCL& bnd, SurfaceTensionCL& sf, FiniteElementT fetype, double SD= 0, double curvDiff= -1)
    : base_( mg, LevelsetCoeffCL(), bnd), idx(fetype), idxC(NULL), MLPhi( &idx), PhiC(NULL), curvDiff_( curvDiff), SD_( SD),
        SF_(SF_ImprovedLBVar), sf_(sf), perDirections(NULL), IsDG(false)
    {}

    virtual ~LevelsetP2CL(){
      if (perDirections) delete perDirections;
    }

    static LevelsetP2CL * Create(  MultiGridCL& mg, const LsetBndDataCL& bnd, SurfaceTensionCL& sf, const ParamCL & P);
    static LevelsetP2CL * Create(  MultiGridCL& mg, const LsetBndDataCL& bnd, SurfaceTensionCL& sf,
                                   bool discontinuous = false, double SD = 0, double curvdiff = -1);


    /// Update PhiC (do nothing if continuous anyway)
    virtual void UpdateContinuous( ) = 0;
    /// Update Phi (do nothing if continuous anyway)
    virtual void UpdateDiscontinuous( ) = 0;

    /// \name Numbering
    ///@{
    void CreateNumbering(Uint level);
    void CreateNumbering( Uint level, MLIdxDescCL* idx);
    void DeleteNumbering( MLIdxDescCL* idx)
        { idx->DeleteNumbering( MG_); }
    ///@}

    /// initialize level set function
    virtual void Init( instat_scalar_fun_ptr, double t = 0) = 0;

    /// \remarks call SetupSystem \em before calling SetTimeStep!
    template<class DiscVelSolT>
    void SetupSystem( const DiscVelSolT&, const double);
    /// Reparametrization of the level set function.
    void Reparam( int method=03, bool Periodic= false);

    /// \brief Perform downwind numbering
    template <class DiscVelSolT>
    PermutationT downwind_numbering (const DiscVelSolT& vel, IteratedDownwindCL dw);

    /// returns information about level set function and interface.
    template<class DiscVelSolT>
    void   GetInfo( double& maxGradPhi, double& Volume, Point3DCL& bary, Point3DCL& vel, const DiscVelSolT& vel_sol, Point3DCL& minCoord, Point3DCL& maxCoord, double& surfArea) const;
    /// returns information about level set function and interface for film.
    template<class DiscVelSolT>
    void   GetFilmInfo( double& maxGradPhi, double& Volume, Point3DCL& vel, const DiscVelSolT& vel_sol,  double& x, double& z, double& h) const;
    /// returns the maximum and minimum of the gradient of phi
    void   GetMaxMinGradPhi(double& maxGradPhi, double& minGradPhi) const;
    /// returns approximate volume of domain where level set function is negative. For l > 0 the level set function is evaluated as a linear FE-function on the principal lattice of order l.
    /// l = 1 : integration on the tetra itself. l = 2 integration on the regular refinement.
    /// l < 0 : extrapolation from current level lvl to lvl - l - 1
    double GetVolume( double translation= 0, int l= 2) const;
    /// Volume correction to ensure no loss or gain of mass.
    void AdjustVolume() const;
    void InitVolume (double vol);
    const VolumeAdjustmentCL* GetVolumeAdjuster() const { return volume_adjuster_.get(); }
    VolumeAdjustmentCL* GetVolumeAdjuster() { return volume_adjuster_.get(); }
    /// Apply smoothing to \a SmPhi, if curvDiff_ > 0
    void MaybeSmooth( VectorCL& SmPhi) const { if (curvDiff_>0) SmoothPhi( SmPhi, curvDiff_); }
    /// Set type of surface force.
    void SetSurfaceForce( SurfaceForceT SF) { SF_= SF; }
    /// Get type of surface force.
    SurfaceForceT GetSurfaceForce() const { return SF_; }

    ///returns the area of the two-phase flow interface(\phi=0)
    double GetInterfaceArea() const;
    ///returns the area of the solid-liquid(phi<0) interface
    double GetWetArea() const;
    /// Discretize surface force
    void   AccumulateBndIntegral( VecDescCL& f) const;
    /// Clear all matrices, should be called after grid change to avoid reuse of matrix pattern
    void   ClearMat() { E.clear(); H.clear(); }
    /// \name Evaluate Solution
    ///@{
    const_DiscSolCL GetSolution() const
        { return const_DiscSolCL( PhiC, &BndData_, &MG_); }
    const_DiscSolCL GetSolution( const VecDescCL& MyPhi) const
        { return const_DiscSolCL( &MyPhi, &BndData_, &MG_); }
    ///@}
    /// Set number of used levels
    void SetNumLvl( size_t n);

    ///Set PeriodicDirections
    void SetPeriodicDirections( perDirSetT* pperDirections)
    {
        if (perDirections != NULL) delete perDirections;
        if (pperDirections != NULL) perDirections = new perDirSetT(pperDirections->size());
        *perDirections = *pperDirections;
    }

    void UpdateMLPhi()
    {
        if (idxC->size() == 1) return;
        MLPhi.SetIdx(idxC);
        ProlongationT::const_iterator prolong = P_.GetFinestIter();
        MLDataCL<VecDescCL>::iterator mlphi = MLPhi.GetFinestIter();
        mlphi->Data = PhiC->Data;
        for (size_t lvl = MG_.GetLastLevel(); lvl >= 1; --lvl){
            const VectorCL tmp (transp_mul(*prolong, mlphi->Data)); // tmp is distributed
            (--mlphi)->Data = tmp;
#ifdef _PAR
            mlphi->RowIdx->GetEx().Accumulate(mlphi->Data);
#endif
            --prolong;
        }
    }
    ProlongationT* GetProlongation() { return &P_; }
};


/// levelset class for continuous P2FE
///differs from discontinuous P2FE by Init and SetUpSystem
class LevelsetP2ContCL: public LevelsetP2CL
{
  public:
    typedef P2EvalCL<double, const LsetBndDataCL, VecDescCL>       DiscSolCL;
    typedef P2EvalCL<double, const LsetBndDataCL, const VecDescCL> const_DiscSolCL;
    typedef std::vector<Point3DCL> perDirSetT;
    typedef LevelsetP2CL base_;
    using base_::MG_;
    using base_::Coeff_;
    using base_::BndData_;
    using base_::GetBndData;
    using base_::GetMG;
    using base_::idx;
    using base_::idxC;
    using base_::Phi;
    using base_::PhiC;
  protected:
    using base_::IsDG;

    public:
    LevelsetP2ContCL( MultiGridCL& mg, const LsetBndDataCL& bnd, SurfaceTensionCL& sf, double SD= 0, double curvDiff= -1)
        : base_( mg, bnd, sf, P2_FE, SD, curvDiff)
    {
        PhiC = &Phi;
        idxC = &idx;
    }

    /// Update PhiC (do nothing)
    virtual void UpdateContinuous( );
    /// Update Phi (do nothing)
    virtual void UpdateDiscontinuous( );

    void Init( instat_scalar_fun_ptr, double t = 0); //void Init( instat_scalar_fun_ptr, double);

    template<class DiscVelSolT>
    void SetupSystem( const DiscVelSolT&, const double);
};


/// levelset class for discontinuous P2FE
///differs from continuous P2FE by Init and SetUpSystem
class LevelsetP2DiscontCL: public LevelsetP2CL
{
  public:
    typedef P2EvalCL<double, const LsetBndDataCL, VecDescCL>       DiscSolCL;
    typedef P2EvalCL<double, const LsetBndDataCL, const VecDescCL> const_DiscSolCL;
    typedef std::vector<Point3DCL> perDirSetT;
    typedef LevelsetP2CL base_;
    using base_::MG_;
    using base_::Coeff_;
    using base_::BndData_;
    using base_::GetBndData;
    using base_::GetMG;
    using base_::idx;
    using base_::idxC;
    using base_::Phi;
    using base_::PhiC;
  protected:
    using base_::IsDG;
  public:
    MLIdxDescCL             idxContinuous;
    VecDescCL             PhiContinuous;        ///< level set function


    public:
    LevelsetP2DiscontCL( MultiGridCL& mg, const LsetBndDataCL& bnd, SurfaceTensionCL& sf, double SD= 0, double curvDiff= -1)
        : base_( mg, bnd, sf, P2D_FE, SD, curvDiff), idxContinuous(P2_FE, mg.GetNumLevel())
    {
        IsDG = true;
        PhiC = &PhiContinuous;
        idxC = &idxContinuous;
    }

    /// \name Numbering
    ///@{
    /* virtual void CreateNumbering( Uint level, IdxDescCL* idx, match_fun match= 0); */
    /* virtual void DeleteNumbering( IdxDescCL* idx); */
    ///@}

    const_DiscSolCL GetDSolution() const
        { return const_DiscSolCL( &Phi, &BndData_, &MG_); }


    /// Update PhiC (Clement-Call...)
    virtual void UpdateContinuous( );
    /// Update Phi (Prolongation...)
    virtual void UpdateDiscontinuous( );

    void InitProjection( instat_scalar_fun_ptr, double t = 0);
    void Init( instat_scalar_fun_ptr, double t = 0);

    void ApplyZeroOrderClementInterpolation();
    void ApplyClementInterpolation();
    void ProjectContinuousToDiscontinuous();

    template<class DiscVelSolT>
    void SetupSystem( const DiscVelSolT&, const double);
};


/// \brief Observes the MultiGridCL-changes by AdapTriangCL to repair the Function ls.Phi.
///
/// Sequential: The actual work is done in post_refine().<br>
/// Parallel:
/// - In pre_refine_sequence the parallel multigrid is informed about
///  the DOF of the level-set function in order to handle them during the
///  refinement and the load-migration.
/// - In post_refine_sequence the actual work is done.
class LevelsetRepairCL : public MGObserverCL
{
  private:
    LevelsetP2CL& ls_;
    std::unique_ptr<RepairP2CL<double>::type > p2repair_;

  public:
    /// \brief Construct a levelset repair class
    LevelsetRepairCL (LevelsetP2CL& ls)
        : ls_( ls) {}

    void pre_refine  ();
    void post_refine ();

    void pre_refine_sequence  ();
    void post_refine_sequence ();
    const IdxDescCL* GetIdxDesc() const { return ls_.Phi.RowIdx; }
#ifdef _PAR
    const VectorCL*  GetVector()  const { return &ls_.Phi.Data; }
    void swap( IdxDescCL& idx, VectorCL& v) { ls_.Phi.RowIdx->swap(idx); ls_.Phi.Data.swap(v); }
#endif
};

/// \brief volume correction and reparametrization of level set function phi
class LevelsetModifyCL
{
private:
    int    rpm_Freq_,
           rpm_Method_;
    double rpm_MaxGrad_,
           rpm_MinGrad_;

    int    step_;
    bool   per_;

public:
    LevelsetModifyCL( int rpm_Freq, int rpm_Method, double rpm_MaxGrad, double rpm_MinGrad, bool periodic=false) :
        rpm_Freq_( rpm_Freq), rpm_Method_( rpm_Method), rpm_MaxGrad_( rpm_MaxGrad),
        rpm_MinGrad_( rpm_MinGrad), step_( 0), per_(periodic) {}


    void maybeDoReparam( LevelsetP2CL& lset);
    void maybeDoVolCorr( LevelsetP2CL& lset);

    void init() {
        step_++;
    }
};


/// marks all tetrahedra in the band |\p DistFct(x)| < \p width for refinement
bool MarkInterface (instat_scalar_fun_ptr DistFct, double width, MultiGridCL&, Uint f_level=(Uint)(-1), Uint c_level=(Uint)(-1), double t=0.);
/// marks all tetrahedra in the band |\p lset(x)| < \p width for refinement
void MarkInterface ( const LevelsetP2CL::const_DiscSolCL& lset, double width, MultiGridCL& mg);

} // end of namespace DROPS

#include "levelset/levelset.tpp"

#endif

