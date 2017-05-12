/// \file instatstokes2phase.tpp
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

namespace DROPS
{

/// \name Routines for SetupSystem2
//@{
/// \brief P2 / P0 FEs for vel/pr
void SetupSystem2_P2P0( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData,
                        MatrixCL* B, VecDescCL* c, IdxDescCL* RowIdx, IdxDescCL* ColIdx, double t);

/// \brief P2 / P1 FEs (Taylor-Hood) for vel/pr
void SetupSystem2_P2P1( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData,
                        MatrixCL* B, VecDescCL* c, IdxDescCL* RowIdx, IdxDescCL* ColIdx, double t);

/// \brief P2 / P1X FEs (X=extended) for vel/pr
void SetupSystem2_P2P1X( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData,
                        MatrixCL* B, VecDescCL* c, const LevelsetP2CL& lset, IdxDescCL* RowIdx, IdxDescCL* ColIdx, double t);

/// \brief P2X / P1X FEs (X=extended) for vel/pr
void SetupSystem2_P2RP1X( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData,
                        MatrixCL* B, VecDescCL* c, const LevelsetP2CL& lset, IdxDescCL* RowIdx, IdxDescCL* ColIdx, double t);

/// \brief P2X / P1 FEs (X=extended) for vel/pr
void SetupSystem2_P2RP1( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData,
                        MatrixCL* B, VecDescCL* c, const LevelsetP2CL& lset, IdxDescCL* RowIdx, IdxDescCL* ColIdx, double t);

/// \brief P2 / P1D FEs for vel/pr
void SetupSystem2_P2P1D( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData,
                        MatrixCL* B, VecDescCL* c, IdxDescCL* RowIdx, IdxDescCL* ColIdx, double t);
//@}


/// \name Routines for SetupRhs2
//@{
/// \brief P2 / P0 FEs for vel/pr
void SetupRhs2_P2P0( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData, VecDescCL* c, double t);

/// \brief P2 / P1 FEs (Taylor-Hood) for vel/pr
void SetupRhs2_P2P1( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData, VecDescCL* c, double t);

/// \brief P2 / P1X FEs (X=extended) for vel/pr
void SetupRhs2_P2P1X( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData, VecDescCL* c, const LevelsetP2CL& lset, double t);

/// \brief P2 / P1D FEs for vel/pr
void SetupRhs2_P2P1D( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData, VecDescCL* c, double t);
//@}


/// \name Routines for SetupPrMass
//@{
/// \brief P0 FEs for pr
void SetupPrMass_P0(const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& matM, IdxDescCL& RowIdx, const LevelsetP2CL& lset);

/// \brief P1 FEs for pr
void SetupPrMass_P1(const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& matM, IdxDescCL& RowIdx, const LevelsetP2CL& lset);

/// \brief P1X FEs for pr
void SetupPrMass_P1X(const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& matM, IdxDescCL& RowIdx, const LevelsetP2CL& lset);

/// \brief PD FEs for pr
void SetupPrMass_P1D(const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& matM, IdxDescCL& RowIdx, const LevelsetP2CL& lset);

/// \brief P1X FEs for pr hat
void SetupPrMassHat_P1X(const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& matM, IdxDescCL& RowIdx, const LevelsetP2CL& lset);
//@}


/// \name Routines for SetupPrSiff
//@{
/// \brief P1 FEs for pr
void SetupPrStiff_P1( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& A_pr, IdxDescCL& RowIdx, IdxDescCL& ColIdx, const LevelsetP2CL& lset);

/// \brief P1X FEs for pr
/// \todo: As in SetupPrMass_P1X, replace the smoothed density-function with integration
///        over the inner and outer part.
void SetupPrStiff_P1X(const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& A_pr, IdxDescCL& RowIdx, IdxDescCL& ColIdx, const LevelsetP2CL& lset, const StokesVelBndDataCL &velbnd, double lambda);

/// \brief P1D FEs for pr
void SetupPrStiff_P1D( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& A_pr, IdxDescCL& RowIdx, IdxDescCL& ColIdx, const LevelsetP2CL& lset);
//@}

//helper routine for calculating P1P1 products like in mass matrix
void computeLocalP2_pipj(LocalP2CL<> (&pipj)[4][4] );


//*****************************************************************************
//                               VelocityRepairCL
//*****************************************************************************
inline void VelocityRepairCL::pre_refine()
{
    p2repair_= std::unique_ptr<RepairP2CL<Point3DCL>::type >(
        new RepairP2CL<Point3DCL>::type( stokes_.GetMG(), stokes_.v, stokes_.GetBndData().Vel));
}

inline void
  VelocityRepairCL::post_refine ()
{
    VelVecDescCL loc_v;
    VelVecDescCL& v= stokes_.v;
    Uint LastLevel= stokes_.GetMG().GetLastLevel();
    MLIdxDescCL loc_vidx( vecP2_FE, stokes_.vel_idx.size());

    loc_vidx.CreateNumbering( LastLevel, stokes_.GetMG(), stokes_.GetBndData().Vel);
    /*
    if (LastLevel != v.RowIdx->TriangLevel()) {
        std::cout << "LastLevel: " << LastLevel
                  << " old v->TriangLevel(): " << v.RowIdx->TriangLevel() << std::endl;
        throw DROPSErrCL( "VelocityRepairCL::post_refine: Sorry, not yet implemented.");
    }
    */
    loc_v.SetIdx( &loc_vidx);

    p2repair_->repair( loc_v);

    v.Clear( v.t);
    stokes_.vel_idx.DeleteNumbering( stokes_.GetMG());

    stokes_.vel_idx.swap( loc_vidx);
    v.SetIdx( &stokes_.vel_idx);
    v.Data= loc_v.Data;
}


//*****************************************************************************
//                               PressureRepairCL
//*****************************************************************************

inline void PressureRepairCL::pre_refine()
{
    p1repair_= std::unique_ptr<RepairP1CL<double>::type >(
        new RepairP1CL<double>::type( stokes_.GetMG(), stokes_.p, stokes_.GetBndData().Pr));
}

inline void
  PressureRepairCL::post_refine ()
{
    VecDescCL loc_p;
    MLIdxDescCL loc_pidx( stokes_.GetPrFE(), stokes_.pr_idx.size());
    VecDescCL& p= stokes_.p;

    loc_pidx.CreateNumbering( stokes_.GetMG().GetLastLevel(), stokes_.GetMG(), stokes_.GetBndData().Pr, ls_.PhiC, &ls_.GetBndData());
    loc_p.SetIdx( &loc_pidx);

    p1repair_->repair( loc_p);

    p.Clear( p.t);
    stokes_.pr_idx.DeleteNumbering( stokes_.GetMG());
    stokes_.pr_idx.swap( loc_pidx);
    p.SetIdx( &stokes_.pr_idx);
    p.Data= loc_p.Data;
}

inline void
  PressureRepairCL::pre_refine_sequence ()
{
    p1xrepair_= std::unique_ptr<P1XRepairCL>( new P1XRepairCL( stokes_.GetMG(), stokes_.p));
}

inline void
  PressureRepairCL::post_refine_sequence ()
{
    (*p1xrepair_)();
    p1xrepair_.reset();
}
} // end of namespace DROPS
