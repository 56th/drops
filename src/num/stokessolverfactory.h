/// \file stokessolverfactory.h
/// \brief creates several standard Stokes-solver
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross; SC RWTH Aachen:

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


#include "num/oseensolver.h"
#include "num/prolongation.h"
#include "misc/params.h"

namespace DROPS {

/// codes for Oseen solvers
enum OseenSolverE {
    GCR_OS= 1, iUzawa_OS= 2, MinRes_OS= 3, GMRes_OS= 4, GMResR_OS= 5, IDRs_OS= 7, StokesMGM_OS= 30
};

/// codes for velocity preconditioners (also including smoothers for the StokesMGM_OS)
enum APcE {
    MG_APC= 1, MGsymm_APC= 2, PCG_APC= 3, GMRes_APC= 4, BiCGStab_APC= 5, VankaBlock_APC= 6, IDRs_APC=7, GS_GMRes_APC= 8, AMG_APC= 20, // preconditioners
    PVanka_SM= 30, BraessSarazin_SM= 31 // smoothers, nevertheless listed here
};

/// codes for the pressure Schur complement preconditioners
enum SPcE {
    ISBBT_SPC= 1, DummyPre_SPC=13, ISGhPenKernel_SPC = 10, ISGhPen_SPC =11, MinComm_SPC= 2, ISPre_SPC= 3, ISMG_SPC= 7, BDinvBT_SPC= 5, SIMPLER_SPC=8, MSIMPLER_SPC=9, VankaSchur_SPC= 4, VankaBlock_SPC=6, ISNonlinear_SPC=12
};

/// collects some information on the different Oseen solvers and preconditioners
struct StokesSolverInfoCL
{
    static std::string GetOseenSolverName( int solver) {
        switch(solver) {
            case GCR_OS:       return "GCR";
            case iUzawa_OS:    return "inexact Uzawa";
            case MinRes_OS:    return "PMinRes";
            case GMRes_OS:     return "GMRes";
            case GMResR_OS:    return "GMResR";
            case StokesMGM_OS: return "Stokes MG";
            case IDRs_OS:      return "IDR(s)";
            default:           return "unknown";
        }
    }
    static std::string GetVelPreName( int pre) {
        switch(pre) {
            case MG_APC:           return "multigrid V-cycle";
            case MGsymm_APC:       return "symm. multigrid V-cycle";
            case PCG_APC:          return "PCG iterations";
            case GMRes_APC:        return "Jacobi-GMRes iterations";
            case BiCGStab_APC:     return "BiCGStab iterations";
            case AMG_APC:          return "algebraic multigrid";
            case VankaBlock_APC:   return "block Vanka";
            case PVanka_SM:        return "Vanka smoother";
            case BraessSarazin_SM: return "Braess-Sarazin smoother";
            case IDRs_APC:         return "IDR(s) iterations";
            case GS_GMRes_APC:     return "Gauss-Seidel-GMRes iterations";
            default:               return "unknown";
        }
    }
    static std::string GetSchurPreName( int pre) {
        switch(pre) {
            case ISBBT_SPC:        return "ISBBT (modified Cahouet-Chabard)";            
            case MinComm_SPC:      return "MinComm (minimal commutator)";
            case ISPre_SPC:        return "ISPre (Cahouet-Chabard)";
            case ISGhPenKernel_SPC:return "ISGhPenKernelPre (Cahouet-Chabard for ghost penalty) with kernel info";
            case ISGhPen_SPC:      return "ISGhPenPre (Cahouet-Chabard for ghost penalty)";
            case ISNonlinear_SPC:  return "ISNonlinearPreCL (Cahouet-Chabard)";
            case ISMG_SPC:         return "ISMGPre (multigrid Cahouet-Chabard)";
            case BDinvBT_SPC:      return "B D^-1 B^T";
            case SIMPLER_SPC:      return "SIMPLER";
            case MSIMPLER_SPC:     return "MSIMPLER";
            case VankaSchur_SPC:   return "Vanka Schur";
            case VankaBlock_SPC:   return "block Vanka";
            case PVanka_SM:        return "Vanka smoother";
            case BraessSarazin_SM: return "Braess-Sarazin smoother";
            case DummyPre_SPC:     return "No Preconditioner";
            default:               return "unknown";
        }
    }
    static bool IsBlockPre( APcE pre) { return pre==VankaBlock_APC; }
    static bool IsBlockPre( SPcE pre) { return pre==SIMPLER_SPC || pre==MSIMPLER_SPC; }
    static bool IsSmoother( int pre)  { return pre==PVanka_SM || pre==BraessSarazin_SM; }
    static bool EqualStokesMGSmoother (APcE preA, SPcE preS)
        { return (preA == PVanka_SM && int(preS) == int(PVanka_SM)) || (preA == BraessSarazin_SM && int(preS) == int(BraessSarazin_SM)); }
};

/*******************************************************************
*   S t o k e s S o l v e r F a c t o r y B a s e   C L            *
*******************************************************************/
/// \brief Creates a StokesSolverCL* and manages the preconditioner.
/// Interface for all Stokes solver factories.
/// Construction of an Oseen solver, e.g. Inexact Uzawa with GMRes and BBT preconditioner: 2*10000 + 4*100 + 2. Note: Not all combinations are implemented!
/**
    <table border="3">
    <tr><th> no </th><th> Oseen-Solver-Type </th><th> Type of Preconditioner for A-Block </th><th> Type of Preconditioner for S </th></tr>
    <tr><td>  1 </td><td> GCR               </td><td> MultiGrid V-cycle                  </td><td> ISBBTPreCL                   </td></tr>
    <tr><td>  2 </td><td> Inexact Uzawa     </td><td> symm. Multigrid V-cycle            </td><td> MinCommPreCL                 </td></tr>
    <tr><td>  3 </td><td> MinRes            </td><td> PCG                                </td><td> ISPreCL                      </td></tr>
    <tr><td>  4 </td><td> GMRes             </td><td> Jacobi-GMRes                       </td><td> VankaSchurPreCL              </td></tr>
    <tr><td>  5 </td><td> GMResR            </td><td> BiCGStab                           </td><td> BD^{-1}BT                    </td></tr>
    <tr><td>  6 </td><td>                   </td><td> VankaPreCL                         </td><td> VankaPreCL                   </td></tr>
    <tr><td>  7 </td><td> IDR(s)            </td><td> IDR(s)                             </td><td> ISMGPreCL                    </td></tr>
    <tr><td>  8 </td><td>                   </td><td> Gauss-Seidel-GMRes                 </td><td> SIMPLER                      </td></tr>
    <tr><td>  9 </td><td>                   </td><td>                                    </td><td> MSIMPLER                     </td></tr>
    <tr><td> 10 </td><td>                   </td><td>                                    </td><td> ISGhPenKernelPreCL           </td></tr>
    <tr><td> 11 </td><td>                   </td><td>                                    </td><td> ISGhPenPreCL                 </td></tr>
    <tr><td> 12 </td><td>                   </td><td>                                    </td><td> ISNonlinearPreCL             </td></tr>
    <tr><td> 13 </td><td>                   </td><td>                                    </td><td> DummyPreCL                   </td></tr>
    <tr><td> 20 </td><td>                   </td><td> HYPRE-AMG                          </td><td>                              </td></tr>
    <tr><td> 30 </td><td> StokesMGM         </td><td> PVankaSmootherCL                   </td><td> PVankaSmootherCL             </td></tr>
    <tr><td> 31 </td><td>                   </td><td> BSSmootherCL                       </td><td> BSSmootherCL                 </td></tr>
    </table>*/
template <class StokesT, class ProlongationVelT= MLDataCL<ProlongationCL<Point3DCL> >, class ProlongationPT= MLDataCL<ProlongationCL<double> > >
class StokesSolverFactoryBaseCL
{
  protected:
    StokesT& Stokes_;           ///< Stokes problem
    ParamCL& P_;                ///< Parameter for tolerances, iteration number, type of solver, ...
    SPcE     SPc_;              ///< type of preconditioner for S
    APcE     APc_;              ///< type of preconditioner for A-block
    int      OseenSolver_;      ///< type of Oseen solver

  public:
    StokesSolverFactoryBaseCL( StokesT& Stokes, ParamCL& P) : Stokes_( Stokes), P_( P),
                               SPc_( static_cast<SPcE>( P_.get<int>("Solver") % 100)),
                               APc_( static_cast<APcE>(  (P_.get<int>("Solver") / 100) % 100)),
                               OseenSolver_( (P_.get<int>("Solver") /10000) % 100) {}
    virtual ~StokesSolverFactoryBaseCL() {}

    /// print some infos about solver combination
    void PrintSolverInfo( std::ostream&) const;
    /// Set the A-block in the minimal commutator
    virtual void       SetMatrixA ( const MatrixCL*) = 0;
    /// Set all matrices in Schur complement preconditioner
    virtual void       SetMatrices( const MLMatrixCL*, const MLMatrixCL*, const MLMatrixCL*, const MLMatrixCL*, const MLIdxDescCL* pr_idx) = 0;
    /// Returns pointer to prolongation for velocity
    virtual ProlongationVelT* GetPVel() = 0;
    /// Returns pointer to prolongation for pressure
    virtual ProlongationPT*   GetPPr() = 0;
    /// Returns a stokes solver with specifications from ParamsT C
    virtual StokesSolverBaseCL* CreateStokesSolver() = 0;
};

/*******************************************************************
*   S t o k e s S o l v e r F a c t o r y H e l p e r  C L         *
********************************************************************/
class StokesSolverFactoryHelperCL
{
  public:
    int GetOseenSolver( const ParamCL& P) const { return (P.get<int>("Solver") / 10000) % 100; }
    int GetAPc( const ParamCL& P) const { return (P.get<int>("Solver") / 100) % 100; }
    int GetSPc( const ParamCL& P) const { return P.get<int>("Solver") % 100; }
    bool VelMGUsed ( const ParamCL& P) const
    {
        const int APc = GetAPc( P);
        return (( APc == MG_APC) || (APc == MGsymm_APC) || (APc == PVanka_SM) || (APc == BraessSarazin_SM));
    }
    bool PrMGUsed  ( const ParamCL& P) const
    {
        const int APc = GetAPc( P),
            SPc = GetSPc( P);
        return (( APc == PVanka_SM) || ( APc == BraessSarazin_SM) || (SPc == ISMG_SPC));
    }
};

/*******************************************************************
*   S t o k e s S o l v e r F a c t o r y  C L                     *
********************************************************************/
template <class StokesT, class ProlongationVelT= MLDataCL<ProlongationCL<Point3DCL> >, class ProlongationPT= MLDataCL<ProlongationCL<double> > >
class StokesSolverFactoryCL : public StokesSolverFactoryBaseCL<StokesT, ProlongationVelT, ProlongationPT>
{
  private:
    typedef StokesSolverFactoryBaseCL<StokesT, ProlongationVelT, ProlongationPT> base_;
    using base_::Stokes_;
    using base_::P_;
    using base_::OseenSolver_;
    using base_::APc_;
    using base_::SPc_;
    using base_::PrintSolverInfo;
    double kA_, kM_;

// generic preconditioners
    JACPcCL  JACPc_;
    GSPcCL   GSPc_;
#ifndef _PAR
    typedef SSORPcCL      SymmPcPcT;
#else
    //typedef ChebyshevPcCL SymmPcPcT;
    typedef JACPcCL  SymmPcPcT;
#endif
    SymmPcPcT symmPcPc_;

// PC for instat. Schur complement
    SchurPreBaseCL  *spc_;
    DummyPreCL      nopc_;
    ISBBTPreCL      bbtispc_;
    MinCommPreCL    mincommispc_;
    BDinvBTPreCL    bdinvbtispc_;
    VankaSchurPreCL vankaschurpc_;
    ISPreCL         isprepc_;
    ISGhPenKernelPreCL  isgpkernpc_;
    ISGhPenPreCL        isgppc_;
    ISMGPreCL<ProlongationPT> ismgpre_;
    typedef PCGSolverCL<SymmPcPcT> PCGSolverT;
    PCGSolverT isnonlinearprepc1_, isnonlinearprepc2_;
    ISNonlinearPreCL<PCGSolverT> isnonlinearpc_;

// PC for A-block
    ExpensivePreBaseCL *apc_;
    // MultiGrid symm.
#ifdef _PAR
    typedef MLSmootherCL<ChebyshevsmoothCL> SmootherT;
    //typedef MLSmootherCL<JORsmoothCL> SmootherT;
#else
    typedef MLSmootherCL<SSORsmoothCL> SmootherT;
#endif

    SmootherT smoother_;
    PCGSolverT   coarsesolversymm_;
    MGSolverCL<SmootherT, PCGSolverT, ProlongationVelT> MGSolversymm_;
    typedef SolverAsPreCL<MGSolverCL<SmootherT, PCGSolverT, ProlongationVelT> > MGsymmPcT;
    MGsymmPcT MGPcsymm_;

    // Multigrid nonsymm.
    GMResSolverCL<JACPcCL> coarsesolver_;
    MGSolverCL<SmootherT, GMResSolverCL<JACPcCL>, ProlongationVelT > MGSolver_;
    typedef SolverAsPreCL<MGSolverCL<SmootherT, GMResSolverCL<JACPcCL>, ProlongationVelT> > MGPcT;
    MGPcT MGPc_;

    //JAC-GMRes
    typedef GMResSolverCL<JACPcCL> GMResSolverT;
    GMResSolverT GMResSolver_;
    typedef SolverAsPreCL<GMResSolverT> GMResPcT;
    GMResPcT GMResPc_;

    //GS-GMRes
    typedef GMResSolverCL<GSPcCL> GS_GMResSolverT;
    GS_GMResSolverT GS_GMResSolver_;
    typedef SolverAsPreCL<GS_GMResSolverT> GS_GMResPcT;
    GS_GMResPcT GS_GMResPc_;

    //JAC-BiCGStab
    typedef BiCGStabSolverCL<JACPcCL> BiCGStabSolverT;
    BiCGStabSolverT BiCGStabSolver_;
    typedef SolverAsPreCL<BiCGStabSolverT> BiCGStabPcT;
    BiCGStabPcT BiCGStabPc_;

    //PCG
    PCGSolverT PCGSolver_;
    typedef SolverAsPreCL<PCGSolverT> PCGPcT;
    PCGPcT PCGPc_;

    //IDR(s)
    typedef IDRsSolverCL<JACPcCL> IDRsSolverT;
    IDRsSolverT IDRsSolver_;
    typedef SolverAsPreCL<IDRsSolverT> IDRsPcT;
    IDRsPcT IDRsPc_;

// Block PC for Oseen problem
    typedef BlockPreCL<ExpensivePreBaseCL, SchurPreBaseCL, DiagSpdBlockPreCL>  DiagBlockPcT;
    typedef BlockPreCL<ExpensivePreBaseCL, SchurPreBaseCL, LowerBlockPreCL>    LowerBlockPcT;
    typedef BlockPreCL<ExpensivePreBaseCL, BDinvBTPreCL, SIMPLERBlockPreCL>    SIMPLERBlockPcT;
    typedef BlockPreCL<ExpensivePreBaseCL, SchurPreBaseCL, UpperBlockPreCL>    UpperBlockPcT;

    DiagBlockPcT    *DBlock_;
    LowerBlockPcT   *LBlock_;
    UpperBlockPcT   *UBlock_;
    SIMPLERBlockPcT *SBlock_;
    VankaPreCL      vankapc_;

//GCR solver
    typedef GCRSolverCL<LowerBlockPcT>   GCR_LBlockT;
    typedef GCRSolverCL<SIMPLERBlockPcT> GCR_SBlockT;
    typedef GCRSolverCL<VankaPreCL>      GCR_VankaT;

    GCR_LBlockT *GCRLBlock_;
    GCR_SBlockT *GCRSBlock_;
    GCR_VankaT  *GCRVanka_;

//GMRes solver
    typedef GMResSolverCL<LowerBlockPcT> GMRes_LBlockT;
    typedef GMResSolverCL<VankaPreCL>    GMRes_VankaT;

    GMRes_LBlockT *GMResLBlock_;
    GMRes_VankaT  *GMResVanka_;

// GMResR solver
    typedef GMResRSolverCL<LowerBlockPcT> GMResR_LBlockT;
    typedef GMResRSolverCL<VankaPreCL>    GMResR_VankaT;

    GMResR_LBlockT *GMResRLBlock_;
    GMResR_VankaT  *GMResRVanka_;

// Lanczos
    typedef PLanczosONBCL<VectorCL, DiagBlockPcT> LanczosT;

    LanczosT *lanczos_;

// MinRes solver
    typedef PMResSolverCL<LanczosT> MinResT;
    MinResT *MinRes_;

// IDR(s) solver
    typedef IDRsSolverCL<LowerBlockPcT> IDRs_LBlockT;
    typedef IDRsSolverCL<VankaPreCL>    IDRs_VankaT;

    IDRs_LBlockT *IDRsLBlock_;
    IDRs_VankaT  *IDRsVanka_;

// coarse grid solver
    DiagBlockPcT DiagPCGBBTOseenPc_, DiagGMResMinCommPc_;
    LanczosT lanczosPCGBBT_;
    MinResT minressolver_;
    BlockMatrixSolverCL<MinResT> coarse_blockminressolver_;

    //GCR solver
    GCRSolverCL<DiagBlockPcT> gcrsolver_;
    BlockMatrixSolverCL<GCRSolverCL<DiagBlockPcT> > coarse_blockgcrsolver_;

    PVankaSmootherCL vankasmoother_;
    BSSmootherCL bssmoother_;

    //StokesMGSolver
    StokesMGSolverCL<PVankaSmootherCL, ProlongationVelT, ProlongationPT>* mgvankasolver_;
    StokesMGSolverCL<BSSmootherCL,     ProlongationVelT, ProlongationPT>* mgbssolver_;

    ExpensivePreBaseCL* CreateAPc();
    SchurPreBaseCL*     CreateSPc();

  public:
    StokesSolverFactoryCL(StokesT& Stokes, ParamCL& PSolver, const ParamCL& PTime);
    ~StokesSolverFactoryCL();

    // checks, whether the combination of Oseen solver, A and S preconditioner is valid.
    bool ValidSolverCombination( std::ostream* os= 0) const;
    /// Set the A-block in the minimal commutator
    void       SetMatrixA ( const MatrixCL* A) { mincommispc_.SetMatrixA(A); bdinvbtispc_.SetMatrixA(A); }
    /// Set all matrices in Schur complement preconditioner (only for StokesMGM)
    void       SetMatrices( const MLMatrixCL* A, const MLMatrixCL* B, const MLMatrixCL* Mvel, const MLMatrixCL* M, const MLIdxDescCL* pr_idx);

    /// Set all matrices in Schur complement preconditioner (only for StokesMGM)
    void       SetMatrices( const MLMatrixCL* A, const MLMatrixCL* B, const MLMatrixCL* C, const MLMatrixCL* Mvel, const MLMatrixCL* M, const MLIdxDescCL* pr_idx);

    /// Returns pointer to prolongation for velocity
    ProlongationVelT* GetPVel();
    /// Returns pointer to prolongation for pressure
    ProlongationPT*   GetPPr();
    /// Returns a stokes solver with specifications from ParamsT C
    StokesSolverBaseCL* CreateStokesSolver();

    /// Returns a pointer to the used velocity preconditioner for the upper left block (aka A block)
    PreBaseCL*      GetVelPrePtr()   { return apc_; }
    /// Returns a pointer to the used pressure preconditioner for the schur complement
    SchurPreBaseCL* GetSchurPrePtr() { return spc_; }

    PVankaSmootherCL&      GetVankaSmoother () { return vankasmoother_; }
    VankaSchurPreCL&       GetVankaSchurPc ()  { return vankaschurpc_; }
};

template <class StokesT, class ProlongationVelT, class ProlongationPT>
StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::
    StokesSolverFactoryCL(StokesT& Stokes, ParamCL& PSolver, const ParamCL& PTime)
    : base_(Stokes, PSolver),
        kA_(PTime.get<int>("NumSteps") != 0 ? PTime.get<double>("NumSteps")/PTime.get<double>("FinalTime") : 0.0), // PTime.get<int>("NumSteps") == 0: stat. problem
        kM_(PTime.get<double>("Theta")),
        // schur complement preconditioner        
        nopc_ (kA_,kM_),
        bbtispc_    ( &Stokes_.B.Data.GetFinest(), &Stokes_.prM.Data.GetFinest(), &Stokes_.M.Data.GetFinest(), Stokes_.pr_idx.GetFinest(), kA_, kM_, PSolver.get<double>("PcSTol"), PSolver.get<double>("PcSTol") /* enable regularization: , 0.707*/),
        mincommispc_( 0, &Stokes_.B.Data.GetFinest(), &Stokes_.M.Data.GetFinest(), &Stokes_.prM.Data.GetFinest(),Stokes_.pr_idx.GetFinest(), PSolver.get<double>("PcSTol") /* enable regularization: , 0.707*/),
        bdinvbtispc_( 0, &Stokes_.B.Data.GetFinest(), &Stokes_.C.Data.GetFinest(), &Stokes_.M.Data.GetFinest(), &Stokes_.prM.Data.GetFinest(),Stokes_.pr_idx.GetFinest(), PSolver.get<double>("PcSTol") /* enable regularization: , 0.707*/),
        vankaschurpc_( &Stokes.pr_idx), isprepc_( Stokes.prA.Data, Stokes.prM.Data, kA_, kM_),
        isgpkernpc_( &Stokes_.prA.Data.GetFinest(), &Stokes_.prM.Data.GetFinest(), &Stokes_.C.Data.GetFinest(), Stokes_.cKernel, kA_, kM_, PSolver.get<double>("PcSTol"), PSolver.get<double>("PcSTol"), PSolver.get<int>("PcSIter",150), &std::cout),
        isgppc_( &Stokes_.prA.Data.GetFinest(), &Stokes_.prM.Data.GetFinest(), &Stokes_.C.Data.GetFinest(), kA_, kM_, PSolver.get<double>("PcSTol"), PSolver.get<double>("PcSTol"), PSolver.get<int>("PcSIter",150) ),
        ismgpre_( Stokes.prA.Data, Stokes.prM.Data, kA_, kM_, Stokes.pr_idx),
        isnonlinearprepc1_( symmPcPc_, 100, PSolver.get<double>("PcSTol"), true),
        isnonlinearprepc2_( symmPcPc_, 100, PSolver.get<double>("PcSTol"), true),
        isnonlinearpc_( isnonlinearprepc1_, isnonlinearprepc2_, Stokes_.prA.Data.GetFinest(), Stokes_.prM.Data.GetFinest(), kA_, kM_),        
        // preconditioner for A
        smoother_( 1.0), coarsesolversymm_( symmPcPc_, 500, 1e-6, true),
        MGSolversymm_ ( smoother_, coarsesolversymm_, PSolver.get<int>("PcAIter"), PSolver.get<double>("PcATol"), Stokes.vel_idx, false),
        MGPcsymm_( MGSolversymm_),
        coarsesolver_( JACPc_, 500, 500, 1e-6, true),        
        MGSolver_ ( smoother_, coarsesolver_, PSolver.get<int>("PcAIter"), PSolver.get<double>("PcATol"), Stokes.vel_idx, false), MGPc_( MGSolver_),
        GMResSolver_( JACPc_, PSolver.get<int>("PcAIter"), /*restart*/ 100, PSolver.get<double>("PcATol"), /*rel*/ true), GMResPc_( GMResSolver_),
        GS_GMResSolver_( GSPc_, PSolver.get<int>("PcAIter"), /*restart*/ 100, PSolver.get<double>("PcATol"), /*rel*/ true), GS_GMResPc_( GS_GMResSolver_),
        BiCGStabSolver_( JACPc_, PSolver.get<int>("PcAIter"), PSolver.get<double>("PcATol"), /*rel*/ true),BiCGStabPc_( BiCGStabSolver_),
        PCGSolver_( symmPcPc_, PSolver.get<int>("PcAIter"), PSolver.get<double>("PcATol"), true), PCGPc_( PCGSolver_),
        IDRsSolver_( JACPc_, PSolver.get<int>("PcAIter"), PSolver.get<double>("PcATol"), true), IDRsPc_( IDRsSolver_),
        // block precondtioner
        DBlock_(0), LBlock_(0), SBlock_(0),
        vankapc_( &Stokes.pr_idx),
        // GCR solver
        GCRLBlock_(0), GCRSBlock_(0), GCRVanka_(0),
        // GMRes solver
        GMResLBlock_(0),  GMResVanka_(0),
        GMResRLBlock_(0), GMResRVanka_(0),
        // lanczos objects
        lanczos_ (0),
        // PMinRes solver
        MinRes_(0),
        // IDRs solver
        IDRsLBlock_(0),
        IDRsVanka_(0),
        // coarse grid/direct solver for StokesMGM
        DiagPCGBBTOseenPc_( PCGPc_, bbtispc_), DiagGMResMinCommPc_( GMResPc_, mincommispc_), lanczosPCGBBT_ (DiagPCGBBTOseenPc_),
        minressolver_( lanczosPCGBBT_, 500, 1e-6, true), coarse_blockminressolver_(minressolver_),        
        gcrsolver_( DiagGMResMinCommPc_, 500, 500, 1e-6, true), coarse_blockgcrsolver_(gcrsolver_),
        vankasmoother_( 0, 0.8, &Stokes.pr_idx)
{
    apc_= CreateAPc();
    spc_= CreateSPc();
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::
    ~StokesSolverFactoryCL()
{
    delete MinRes_; delete lanczos_;
    delete GMResRVanka_; delete GMResRLBlock_;
    delete GMResVanka_; delete GMResLBlock_;
    delete GCRVanka_; delete GCRLBlock_; delete GCRSBlock_;
    delete SBlock_; delete LBlock_; delete DBlock_;
    delete IDRsVanka_; delete IDRsLBlock_;
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
bool StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::ValidSolverCombination( std::ostream* os) const
{
    std::string msg;
    bool ok= false;

    if (StokesSolverInfoCL::GetOseenSolverName(OseenSolver_)=="unknown")
        msg= "unknown Oseen solver";
    else if (StokesSolverInfoCL::GetVelPreName(APc_)=="unknown")
        msg= "unknown vel preconditioner";
    else if (StokesSolverInfoCL::GetSchurPreName(SPc_)=="unknown")
        msg= "unknown pr preconditioner";
    else if (OseenSolver_==StokesMGM_OS && (!StokesSolverInfoCL::EqualStokesMGSmoother( APc_, SPc_) || !StokesSolverInfoCL::IsSmoother(APc_) ))
        msg= "Stokes multigrid method requires smoother";
    else if ((StokesSolverInfoCL::IsSmoother(APc_) || StokesSolverInfoCL::IsSmoother(SPc_)) && OseenSolver_!=StokesMGM_OS)
        msg= "smoother makes no sense without multigrid solver";
    else if (OseenSolver_==iUzawa_OS && (StokesSolverInfoCL::IsBlockPre(APc_) || StokesSolverInfoCL::IsBlockPre(SPc_) ))
        msg= "block preconditioner not allowed for inexact Uzawa";
    else if (OseenSolver_==MinRes_OS && (StokesSolverInfoCL::IsBlockPre(APc_) || StokesSolverInfoCL::IsBlockPre(SPc_) ))
        msg= "MinRes requires diagonal block preconditioner";
    else if ((StokesSolverInfoCL::IsBlockPre(APc_) || StokesSolverInfoCL::IsBlockPre(SPc_)) && !StokesSolverInfoCL::EqualStokesMGSmoother( APc_, SPc_) && SPc_!=SIMPLER_SPC && SPc_!=MSIMPLER_SPC)
        msg= "block preconditioner should be the same for vel and pr part";
#ifdef _PAR
    else if (APc_ == GS_GMRes_APC)
        msg= "Gauss-Seidel is not available in parallel, yet";
#endif
    else if (Stokes_.UsesXFEM())
    { // check whether solver is well-defined for XFEM
        ok = true;
        if (OseenSolver_ == StokesMGM_OS)
        {
            msg = "StokesMGM not implemented for P1X-elements";
            ok  = false;
        }
        if (SPc_ == ISMG_SPC)
        {
            msg = "ISMGPreCL not implemented for P1X-elements";
            ok  = false;
        }
/*
        if( Stokes_.usesGhostPen() )
        {
            if( !( SPc_ == ISGhPen_SPC || SPc_ == ISGhPenKernel_SPC || SPc_ == DummyPre_SPC ) )
            {
                msg = StokesSolverInfoCL::GetSchurPreName( SPc_ ) + " not suited for Ghost Penalty stabilization";
                ok = false;
            }
        }*/
    }
    else // all tests passed successfully
        ok= true;

    if (os && !ok)
        (*os) << "invalid solver combination in StokesSolverFactoryCL:\t" << msg << std::endl;
    return ok;
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
void StokesSolverFactoryBaseCL<StokesT, ProlongationVelT, ProlongationPT>::PrintSolverInfo( std::ostream& os) const
{
    os << "Oseen solver info:\t" << StokesSolverInfoCL::GetOseenSolverName( OseenSolver_)
       << "\n + vel precond.  :\t" << StokesSolverInfoCL::GetVelPreName( APc_)
       << "\n + pr  precond.  :\t" << StokesSolverInfoCL::GetSchurPreName( SPc_) << std::endl;
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
ExpensivePreBaseCL* StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::CreateAPc()
{
    switch (APc_) {
        case MG_APC:       return &MGPc_;
        case MGsymm_APC:   return &MGPcsymm_;
        case PCG_APC:      return &PCGPc_;
        case GMRes_APC:    return &GMResPc_;
        case GS_GMRes_APC: return &GS_GMResPc_;
        case BiCGStab_APC: return &BiCGStabPc_;
        case IDRs_APC:     return &IDRsPc_;
        default:           return 0;
    }
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
SchurPreBaseCL* StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::CreateSPc()
{
    switch (SPc_) {
        case ISBBT_SPC:         return &bbtispc_;
        case MinComm_SPC:       return &mincommispc_;
        case ISPre_SPC:         return &isprepc_;
        case ISGhPenKernel_SPC: return &isgpkernpc_;
        case ISGhPen_SPC:       return &isgppc_;
        case ISMG_SPC:          return &ismgpre_;
        case SIMPLER_SPC:
        case MSIMPLER_SPC:
        case BDinvBT_SPC:       return &bdinvbtispc_;
        case VankaSchur_SPC:    return &vankaschurpc_;
        case ISNonlinear_SPC:   return &isnonlinearpc_;
        case DummyPre_SPC:      return &nopc_;
        default:                return 0;
    }
}


template <class StokesT, class ProlongationVelT, class ProlongationPT>
StokesSolverBaseCL* StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::CreateStokesSolver()
{
    PrintSolverInfo( std::cout);
    if (!ValidSolverCombination( &std::cout))
        throw DROPSErrCL("StokesSolverFactoryCL::CreateStokesSolver(): Invalid solver combination");

    StokesSolverBaseCL* stokessolver = 0;

    switch (OseenSolver_) {
        case iUzawa_OS: {
            if (APc_==MGsymm_APC) // symmetric A preconditioner -> use more efficient version of inexact Uzawa
                stokessolver= new InexactUzawaCL<ExpensivePreBaseCL, SchurPreBaseCL, APC_SYM>  ( *apc_, *spc_, P_.template get<int>("Iter"), P_.template get<double>("Tol"), P_.template get<double>("UzawaInnerTol"), P_.template get<int>("UzawaInnerIter"));
            else
                stokessolver= new InexactUzawaCL<ExpensivePreBaseCL, SchurPreBaseCL, APC_OTHER>( *apc_, *spc_, P_.template get<int>("Iter"), P_.template get<double>("Tol"), P_.template get<double>("UzawaInnerTol"), P_.template get<int>("UzawaInnerIter"));
        }
        break;

        case MinRes_OS: { // MinRes requires symmetric block preconditioner, hence we can only use diagonal block preconditioners
            DBlock_= new DiagBlockPcT( *apc_, *spc_);
            lanczos_= new LanczosT( *DBlock_);
            MinRes_= new MinResT( *lanczos_,  P_.template get<int>("Iter"), P_.template get<double>("Tol"), /*relative*/ P_.template get<bool>("RelTol",false));
            stokessolver= new BlockMatrixSolverCL<MinResT>( *MinRes_);
        }
        break;

        case GCR_OS: {
            if (APc_==VankaBlock_APC) {
                GCRVanka_= new GCR_VankaT( vankapc_,  P_.template get<int>("Iter"), P_.template get<int>("Iter"), P_.template get<double>("Tol"), /*rel*/ false);
                stokessolver= new BlockMatrixSolverCL<GCR_VankaT> ( *GCRVanka_);
            } else if (SPc_==SIMPLER_SPC || SPc_==MSIMPLER_SPC) {
                bdinvbtispc_.SetMassLumping( SPc_==MSIMPLER_SPC);
                SBlock_= new SIMPLERBlockPcT( *apc_, bdinvbtispc_);
                GCRSBlock_= new GCR_SBlockT( *SBlock_,  P_.template get<int>("Iter"), P_.template get<int>("Iter"), P_.template get<double>("Tol"), /*rel*/ false);
                stokessolver= new BlockMatrixSolverCL<GCR_SBlockT>( *GCRSBlock_);
            } else {
                LBlock_= new LowerBlockPcT( *apc_, *spc_);
                //GCRLBlock_= new GCR_LBlockT( *LBlock_,  P_.template get<int>("Iter"), P_.template get<int>("Iter"), P_.template get<double>("Tol"), /*rel*/ false);
                GCRLBlock_= new GCR_LBlockT( *LBlock_,  P_.template get<int>("Iter"), P_.template get<int>("Iter"), P_.template get<double>("Tol"), /*rel*/  P_.template get<bool>("RelTol",false));
                stokessolver= new BlockMatrixSolverCL<GCR_LBlockT>( *GCRLBlock_);
            }
        }
        break;

        case GMRes_OS: {
            if (APc_==VankaBlock_APC) {
                GMResVanka_= new GMRes_VankaT( vankapc_,  P_.template get<int>("Iter"), P_.template get<int>("Iter"), P_.template get<double>("Tol"), /*rel*/ false, false, RightPreconditioning);
                stokessolver= new BlockMatrixSolverCL<GMRes_VankaT> ( *GMResVanka_);
            } else {
                LBlock_= new LowerBlockPcT( *apc_, *spc_);
                GMResLBlock_= new GMRes_LBlockT( *LBlock_,  P_.template get<int>("Iter"), P_.template get<int>("Iter"), P_.template get<double>("Tol"), /*rel*/ false, false, RightPreconditioning);
                stokessolver= new BlockMatrixSolverCL<GMRes_LBlockT>( *GMResLBlock_);
            }
        }
        break;

        case GMResR_OS: {
            if (APc_==VankaBlock_APC) {
                GMResRVanka_= new GMResR_VankaT( vankapc_,  P_.template get<int>("Iter"), P_.template get<int>("Iter"), P_.template get<int>("GMResRInnerIter"), P_.template get<double>("Tol"), P_.template get<double>("GMResRInnerTol"), /*rel*/ false);
                stokessolver= new BlockMatrixSolverCL<GMResR_VankaT> ( *GMResRVanka_);
            } else {
                LBlock_= new LowerBlockPcT( *apc_, *spc_);
                GMResRLBlock_= new GMResR_LBlockT( *LBlock_,  P_.template get<int>("Iter"), P_.template get<int>("Iter"), P_.template get<int>("GMResRInnerIter"), P_.template get<double>("Tol"), P_.template get<double>("GMResRInnerTol"), /*rel*/ false);
                stokessolver= new BlockMatrixSolverCL<GMResR_LBlockT>( *GMResRLBlock_);
            }
        }
        break;

        case StokesMGM_OS: {
            if (APc_==PVanka_SM) {
                if (P_.template get<double>("NavierStokes.Nonlinear", 0.0)==0.0) // Stokes
                    mgvankasolver_ = new StokesMGSolverCL<PVankaSmootherCL, ProlongationVelT, ProlongationPT>
                               ( Stokes_.prM.Data, vankasmoother_, coarse_blockminressolver_, P_.template get<int>("Iter"), P_.template get<double>("Tol"), false, 2);
                else
                    mgvankasolver_ = new StokesMGSolverCL<PVankaSmootherCL, ProlongationVelT, ProlongationPT>
                               ( Stokes_.prM.Data, vankasmoother_, coarse_blockgcrsolver_, P_.template get<int>("Iter"), P_.template get<double>("Tol"), false, 2);
                stokessolver = mgvankasolver_;
            }
            else if (APc_==BraessSarazin_SM) {
                if (P_.template get<double>("NavierStokes.Nonlinear", 0.0) ==0.0) // Stokes
                    mgbssolver_ = new StokesMGSolverCL<BSSmootherCL, ProlongationVelT, ProlongationPT>
                               ( Stokes_.prM.Data, bssmoother_, coarse_blockminressolver_, P_.template get<int>("Iter"), P_.template get<double>("Tol"), false, 2);
                else
                    mgbssolver_ = new StokesMGSolverCL<BSSmootherCL, ProlongationVelT, ProlongationPT>
                               ( Stokes_.prM.Data, bssmoother_, coarse_blockgcrsolver_, P_.template get<int>("Iter"), P_.template get<double>("Tol"), false, 2);
                stokessolver = mgbssolver_;
            }
        }
        break;

        case IDRs_OS: {
            if (APc_==VankaBlock_APC) {
                IDRsVanka_= new IDRs_VankaT( vankapc_,  P_.template get<int>("Iter"), P_.template get<double>("Tol"), /*rel*/ false, false, RightPreconditioning);
                stokessolver= new BlockMatrixSolverCL<IDRs_VankaT> ( *IDRsVanka_);
            } else {
                LBlock_= new LowerBlockPcT( *apc_, *spc_);
                IDRsLBlock_= new IDRs_LBlockT( *LBlock_,  P_.template get<int>("Iter"), P_.template get<double>("Tol"), /*rel*/ false, false, RightPreconditioning);
                stokessolver= new BlockMatrixSolverCL<IDRs_LBlockT>( *IDRsLBlock_);
            }
        }
        break;

        default: throw DROPSErrCL("StokesSolverFactoryCL: Unknown Oseen solver");
    }
    if (stokessolver==0)
        throw DROPSErrCL("StokesSolverFactoryCL: Sorry, this solver combination is not implemented, yet");
    stokessolver->SetOutput(&std::cout);
    return stokessolver;
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
void StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::
    SetMatrices( const MLMatrixCL* A, const MLMatrixCL* B, const MLMatrixCL* Mvel, const MLMatrixCL* M, const MLIdxDescCL* pr_idx) {
    if ( APc_ == PVanka_SM || APc_ == BraessSarazin_SM) { //  Vanka or Braess Sarazin smoother
        bbtispc_.SetMatrices(B->GetCoarsestPtr(), Mvel->GetCoarsestPtr(), M->GetCoarsestPtr(), pr_idx->GetCoarsestPtr());
        return;
    }
    if ( SPc_ == VankaSchur_SPC) {              // VankaSchur
        vankaschurpc_.SetAB(A->GetCoarsestPtr(), B->GetCoarsestPtr());
    }
    else {
        mincommispc_.SetMatrices(A->GetFinestPtr(), B->GetFinestPtr(), Mvel->GetFinestPtr(), M->GetFinestPtr(), pr_idx->GetFinestPtr());
        bdinvbtispc_.SetMatrices(A->GetFinestPtr(), B->GetFinestPtr(), Mvel->GetFinestPtr(), M->GetFinestPtr(), pr_idx->GetFinestPtr());
        bbtispc_.SetMatrices(B->GetFinestPtr(), Mvel->GetFinestPtr(), M->GetFinestPtr(), pr_idx->GetFinestPtr());
    }
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
void StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::
    SetMatrices( const MLMatrixCL* A, const MLMatrixCL* B, const MLMatrixCL* C, const MLMatrixCL* Mvel, const MLMatrixCL* M, const MLIdxDescCL* pr_idx) {
    if ( APc_ == PVanka_SM || APc_ == BraessSarazin_SM) { //  Vanka or Braess Sarazin smoother
        bbtispc_.SetMatrices(B->GetCoarsestPtr(), Mvel->GetCoarsestPtr(), M->GetCoarsestPtr(), pr_idx->GetCoarsestPtr());
        return;
    }
    if ( SPc_ == VankaSchur_SPC) {              // VankaSchur
        vankaschurpc_.SetAB(A->GetCoarsestPtr(), B->GetCoarsestPtr());
    }
    else {
        mincommispc_.SetMatrices(A->GetFinestPtr(), B->GetFinestPtr(), Mvel->GetFinestPtr(), M->GetFinestPtr(), pr_idx->GetFinestPtr());
        bdinvbtispc_.SetMatrices(A->GetFinestPtr(), B->GetFinestPtr(), Mvel->GetFinestPtr(), M->GetFinestPtr(), pr_idx->GetFinestPtr());
        bbtispc_.SetMatrices(B->GetFinestPtr(), Mvel->GetFinestPtr(), M->GetFinestPtr(), pr_idx->GetFinestPtr());
    }
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
ProlongationVelT* StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::GetPVel()
{
    switch ( APc_) {
        case MG_APC           : return MGSolver_.GetProlongation();     break;  // general MG
        case MGsymm_APC       : return MGSolversymm_.GetProlongation(); break;  // symm. MG
        case PVanka_SM        : return mgvankasolver_->GetPVel(); break;
        case BraessSarazin_SM : return mgbssolver_->GetPVel();    break;
        default: return 0;
    }
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
ProlongationPT* StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::GetPPr()
{
    switch ( APc_) {
        case PVanka_SM        : return mgvankasolver_->GetPPr(); break;
        case BraessSarazin_SM : return mgbssolver_->GetPPr();    break;
        default: /*silence warning*/;
    }
    if (SPc_ == ISMG_SPC )
        return ismgpre_.GetProlongation();
    return 0;
}

} // end of namespace DROPS

