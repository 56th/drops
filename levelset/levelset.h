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
#include "num/solver.h"
#include "num/bndData.h"
#include "num/fe.h"
#include "levelset/mgobserve.h"
#include "levelset/surfacetension.h"
#include "num/interfacePatch.h"
#include <vector>

#ifdef _PAR
# include "parallel/exchange.h"
# include "num/parsolver.h"
# include "num/parprecond.h"
# include "misc/container.h"
#endif

namespace DROPS
{

typedef BndDataCL<double>    LsetBndDataCL;

enum SurfaceForceT
/// different types of surface forces
{
    SF_LB=0,             ///< Laplace-Beltrami discretization: \f$\mathcal O(h^{1/2})\f$
    SF_ImprovedLB=1,     ///< improved Laplace-Beltrami discretization: \f$\mathcal O(h)\f$
    SF_Const=2,          ///< surface force with constant curvature
    SF_ImprovedLBVar=3   ///< improved Laplace-Beltrami discretization with variable surface tension
};

/// not used at the moment
class LevelsetCoeffCL {

};

class LevelsetP2CL : public ProblemCL< LevelsetCoeffCL, LsetBndDataCL>
/// P2-discretization and solution of the level set equation for two phase flow problems. Bnd_ will be used to impose boundary data on the inflow boundary.
/// At the moment setting all boundary conditions to NoBC is the only valid case.
{
  public:
    typedef P2EvalCL<double, const LsetBndDataCL, VecDescCL>       DiscSolCL;
    typedef P2EvalCL<double, const LsetBndDataCL, const VecDescCL> const_DiscSolCL;
    typedef std::vector<Point3DCL> perDirSetT;
    IdxDescCL             idx;
    VecDescCL             Phi;        ///< level set function

    typedef ProblemCL<LevelsetCoeffCL, LsetBndDataCL> base_;
    using base_::MG_;
    using base_::Coeff_;
    using base_::BndData_;
    using base_::GetBndData;
    using base_::GetMG;

  private:
    double              curvDiff_, ///< amount of diffusion in curvature calculation
                        SD_;       ///< streamline diffusion
    SurfaceForceT       SF_;

    SurfaceTensionCL&   sf_;      ///< data for surface tension
    void SetupSmoothSystem ( MatrixCL&, MatrixCL&)               const;
    void SmoothPhi( VectorCL& SmPhi, double diff)                const;
    double GetVolume_Composite( double translation, int l)    const;
    double GetVolume_Extrapolation( double translation, int l) const;
    perDirSetT* perDirections;    ///< periodic directions

  public:
    MatrixCL            E, H;

    LevelsetP2CL( MultiGridCL& mg, const LsetBndDataCL& bnd, SurfaceTensionCL& sf, double SD= 0, double curvDiff= -1)
    : base_( mg, LevelsetCoeffCL(), bnd), idx( P2_FE), curvDiff_( curvDiff), SD_( SD),
        SF_(SF_ImprovedLB), sf_(sf), perDirections(NULL)
    {}

    ~LevelsetP2CL(){
      if (perDirections) delete perDirections;
    }
    /// \name Numbering
    ///@{
    void CreateNumbering( Uint level, IdxDescCL* idx, match_fun match= 0);
    void DeleteNumbering( IdxDescCL* idx)
        { idx->DeleteNumbering( MG_); }
    ///@}

    /// initialize level set function
    void Init( scalar_fun_ptr);

    /// \remarks call SetupSystem \em before calling SetTimeStep!
    template<class DiscVelSolT>
    void SetupSystem( const DiscVelSolT&, const double);
    /// Reparametrization of the level set function.
    void Reparam( int method=03, bool Periodic= false);

    /// returns information about level set function and interface.
    template<class DiscVelSolT>
    void   GetInfo( double& maxGradPhi, double& Volume, Point3DCL& bary, Point3DCL& vel, const DiscVelSolT& vel_sol, Point3DCL& minCoord, Point3DCL& maxCoord, double& surfArea) const;
    /// returns the maximum and minimum of the gradient of phi
    void   GetMaxMinGradPhi(double& maxGradPhi, double& minGradPhi) const;
    /// returns approximate volume of domain where level set function is negative. For l > 0 the level set function is evaluated as a linear FE-function on the principal lattice of order l. 
    /// l = 1 : integration on the tetra itself. l = 2 integration on the regular refinement.
    /// l < 0 : extrapolation from current level lvl to lvl - l - 1   
    double GetVolume( double translation= 0, int l= 2) const;
    /// volume correction to ensure no loss or gain of mass. The parameter l is passed to GetVolume().
    double AdjustVolume( double vol, double tol, double surf= 0., int l= 2) const;
    /// Set type of surface force.
    void   SetSurfaceForce( SurfaceForceT SF) { SF_= SF; }
    /// Get type of surface force.
    SurfaceForceT GetSurfaceForce() const { return SF_; }
    /// Discretize surface force
    void   AccumulateBndIntegral( VecDescCL& f) const;
    /// Clear all matrices, should be called after grid change to avoid reuse of matrix pattern
    void   ClearMat() { E.clear(); H.clear(); }
    /// \name Evaluate Solution
    ///@{
    const_DiscSolCL GetSolution() const
        { return const_DiscSolCL( &Phi, &BndData_, &MG_); }
    const_DiscSolCL GetSolution( const VecDescCL& MyPhi) const
        { return const_DiscSolCL( &MyPhi, &BndData_, &MG_); }
    ///@}
    
    ///Set PeriodicDirections
    void SetPeriodicDirections( perDirSetT* pperDirections)
    {
        if (perDirections != NULL) delete perDirections;
        if (pperDirections != NULL) perDirections = new perDirSetT(pperDirections->size());
        *perDirections = *pperDirections;
    }
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

  public:
    /// \brief Construct a levelset repair class
    LevelsetRepairCL (LevelsetP2CL& ls)
        : ls_( ls) {}

    void pre_refine  ();
    void post_refine ();

    void pre_refine_sequence  () {}
    void post_refine_sequence () {}
    const IdxDescCL* GetIdxDesc() const { return ls_.Phi.RowIdx; }
};

/// \brief volume correction and reparametrization of level set function phi
class LevelsetModifyCL
{
private:
    int    rpm_Freq_,
           rpm_Method_;
    double rpm_MaxGrad_,
           rpm_MinGrad_;
    int    lvs_VolCorrection_;

    const double Vol_;

    int    step_;
    bool   per_;

public:
    LevelsetModifyCL( int rpm_Freq, int rpm_Method, double rpm_MaxGrad, double rpm_MinGrad, int lvs_VolCorrection, double Vol, bool periodic=false) :
        rpm_Freq_( rpm_Freq), rpm_Method_( rpm_Method), rpm_MaxGrad_( rpm_MaxGrad),
        rpm_MinGrad_( rpm_MinGrad), lvs_VolCorrection_( lvs_VolCorrection), Vol_( Vol), step_( 0), per_(periodic) {}


    void maybeDoReparam( LevelsetP2CL& lset) {
        bool doReparam= rpm_Freq_ && step_%rpm_Freq_ == 0;
        bool doVolCorr= lvs_VolCorrection_ && step_%lvs_VolCorrection_ == 0;

        double lsetmaxGradPhi, lsetminGradPhi;

        if (doReparam) {
            lset.GetMaxMinGradPhi( lsetmaxGradPhi, lsetminGradPhi);
            doReparam = (lsetmaxGradPhi > rpm_MaxGrad_ || lsetminGradPhi < rpm_MinGrad_);
        }

        // reparam levelset function
        if (doReparam) {
            std::cout << "before reparametrization: minGradPhi " << lsetminGradPhi << "\tmaxGradPhi " << lsetmaxGradPhi << '\n';
            lset.Reparam( rpm_Method_, per_);
            lset.GetMaxMinGradPhi( lsetmaxGradPhi, lsetminGradPhi);
            std::cout << "after  reparametrization: minGradPhi " << lsetminGradPhi << "\tmaxGradPhi " << lsetmaxGradPhi << '\n';
            // volume correction after reparametrization
            if (doVolCorr) {
                double dphi= lset.AdjustVolume( Vol_, 1e-9);
                std::cout << "volume correction is " << dphi << std::endl;
                lset.Phi.Data+= dphi;
                std::cout << "new rel. volume: " << lset.GetVolume()/Vol_ << std::endl;
            }
        }
    }

    double maybeDoVolCorr( LevelsetP2CL& lset) {
        bool doVolCorr= lvs_VolCorrection_ && step_%lvs_VolCorrection_ == 0;
        double dphi = 0.0;

        if (!doVolCorr) return dphi;

        if (doVolCorr) {
            dphi= lset.AdjustVolume( Vol_, 1e-9);
            std::cout << "volume correction is " << dphi << std::endl;
            lset.Phi.Data+= dphi;
            std::cout << "new rel. volume: " << lset.GetVolume()/Vol_ << std::endl;
        }
        return dphi;
    }

    void init() {
        step_++;
    }
};



/// marks all tetrahedra in the band |\p DistFct(x)| < \p width for refinement
bool MarkInterface (scalar_fun_ptr DistFct, double width, MultiGridCL&, Uint f_level=(Uint)(-1), Uint c_level=(Uint)(-1));
/// marks all tetrahedra in the band |\p lset(x)| < \p width for refinement
void MarkInterface ( const LevelsetP2CL::const_DiscSolCL& lset, double width, MultiGridCL& mg);

} // end of namespace DROPS

#include "levelset/levelset.tpp"

#endif

