/// \file Poissonsolverfactory.h
/// \brief Solver for Poisson problems functions
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Eva Loch, Volker Reichelt, Marcus Soemers, Yuanjun Zhang, Thorolf Schulte; SC RWTH Aachen: Oliver Fortmeier

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

/** We solve \f$ -\Delta u = f\;\mbox{in}\; \Omega:=[0,1]^3 \f$ for the given
    solution \f$ u(x,y,z):= 64 \cdot xyz (1-x) (1-y) (1-z) \f$, i.e. homogeneous
    Dirichlet conditions are used. A uniform tetrahedral grid is applied as
    a triangulation of \f$ \Omega \f$. GMRES is used as a linear solver for the
    discretized linear equation system. Note, that CG-type methods can be used
    as well because the resulting linear equation system is s.p.d. However,
    since this program acts as a base performance test, GMRES is used here.
*/

#ifndef POISSONSOLVERFACTORY_H_
#define POISSONSOLVERFACTORY_H_

#include "num/krylovsolver.h"
#include "num/precond.h"
#include "misc/params.h"
#ifdef _HYPRE
#include "num/hypre.h"
#endif
#include "num/prolongation.h"

namespace DROPS
{


class PoissonSolverFactoryHelperCL
{
  public:
    bool MGUsed ( ParamCL& P)
    {
        const int PM = P.get<int>(std::string("Poisson.Solver.Solver"));
        return ( PM / 100 == 1 || PM % 10 == 1);
    }
};


class PoissonSolverBaseCL : public SolverBaseCL
{
  public:
    PoissonSolverBaseCL (int maxiter, double tol, bool rel= false, std::ostream* output= 0)
        : SolverBaseCL(maxiter, tol, rel, output){}
#ifdef _PAR
    virtual void Solve( const MatrixCL& A, VectorCL& x, const VectorCL& b, const ExchangeCL& ex) = 0;
    virtual void Solve( const MLMatrixCL& A, VectorCL& x, const VectorCL& b, const ExchangeCL& ex) = 0;
#endif
    virtual void Solve( const MatrixCL& A, VectorCL& x, const VectorCL& b, const DummyExchangeCL& ex) = 0;
    virtual void Solve( const MLMatrixCL& A, VectorCL& x, const VectorCL& b, const DummyExchangeCL& ex) = 0;

};

template <class SolverT>
class PoissonSolverCL : public PoissonSolverBaseCL
{
  private:
    SolverT& solver_;

  public:
    PoissonSolverCL( SolverT& solver) : PoissonSolverBaseCL( -1, -1.0), solver_( solver) {}
#ifdef _PAR
    void Solve( const MatrixCL& A, VectorCL& x, const VectorCL& b, const ExchangeCL& ex) { solver_.Solve( A, x, b, ex); }
    void Solve( const MLMatrixCL& A, VectorCL& x, const VectorCL& b, const ExchangeCL& ex) { solver_.Solve( A, x, b, ex); }
#endif
    void Solve( const MatrixCL& A, VectorCL& x, const VectorCL& b, const DummyExchangeCL& ex) { solver_.Solve( A, x, b, ex); }
    void Solve( const MLMatrixCL& A, VectorCL& x, const VectorCL& b, const DummyExchangeCL& ex) { solver_.Solve( A, x, b, ex); }
	// We overwrite these functions.
    void   SetTol     (double tol) { solver_.SetTol( tol); }
    void   SetMaxIter (int iter)   { solver_.SetMaxIter( iter); }
    void   SetRelError(bool rel)   { solver_.SetRelError( rel); }

    double GetTol     () const { return solver_.GetTol(); }
    int    GetMaxIter () const { return solver_.GetMaxIter(); }
    double GetResid   () const { return solver_.GetResid(); }
    int    GetIter    () const { return solver_.GetIter(); }
    bool   GetRelError() const { return solver_.GetRelError(); }

    void SetOutput( std::ostream* os) { solver_.SetOutput(os); }
};

/// \brief Create a Poisson-Solver (Design-Pattern: Factory class)
/// Construction of a Poisson solver, e.g. CG and SSOR preconditioner: 2*100 + 3. Note: Not all combinations are implemented!
/**
    <table border="3">
    <tr><th> no </th><th> Poisson-Solver-Type </th><th> Type of PC/smoother  </th></tr>
    <tr><td>  0 </td><td>                     </td><td> DummyPc              </td></tr>
    <tr><td>  1 </td><td> MultiGrid V-cycle   </td><td>                      </td></tr>
    <tr><td>  2 </td><td> Preconditioned CG   </td><td> JOR                  </td></tr>
    <tr><td>  3 </td><td> GMRes               </td><td> SSOR                 </td></tr>
    <tr><td>  4 </td><td> Hypre-AMG           </td><td> GS                   </td></tr>
    <tr><td>  5 </td><td>                     </td><td> SGS                  </td></tr>
    <tr><td>  6 </td><td>                     </td><td> SOR                  </td></tr>
    <tr><td>  7 </td><td>                     </td><td> Chebychev            </td></tr>
    </table>*/

template <class ProlongationT= MLDataCL<ProlongationCL<double> > >
class PoissonSolverFactoryCL
{
  private:
    ParamCL P_;
    MLIdxDescCL& idx_;
    ProlongationT* prolongptr_;

// generic preconditioners
#ifdef _PAR
    //typedef ChebyshevPcCL CoarseSolverPcT;
    typedef JACPcCL CoarseSolverPcT;
#else
    typedef SSORPcCL      CoarseSolverPcT;
#endif
    CoarseSolverPcT genpc_;
    JACPcCL         JACPc_;
    SSORPcCL        SSORPc_;
    ChebyshevPcCL   ChebyPc_;

    // MultiGrid symm.
    MLSmootherCL<JORsmoothCL>  jorsmoother_;         // Jacobi
    MLSmootherCL<GSsmoothCL>   gssmoother_;          // Gauss-Seidel
    MLSmootherCL<SGSsmoothCL>  sgssmoother_;         // symmetric Gauss-Seidel
    MLSmootherCL<SORsmoothCL>  sorsmoother_;         // Gauss-Seidel with over-relaxation
    MLSmootherCL<SSORsmoothCL> ssorsmoother_;        // symmetric Gauss-Seidel with over-relaxation
    MLSmootherCL<ChebyshevsmoothCL> chebysmoother_;  // Chebychev-polynomial based smoother
    typedef PCGSolverCL<CoarseSolverPcT> CoarseSolverT;

    CoarseSolverT   coarsesolversymm_;
    typedef MGSolverCL<MLSmootherCL<JORsmoothCL>, CoarseSolverT, ProlongationT> MGSolversymmJORT;
    MGSolversymmJORT MGSolversymmJOR_;
    typedef MGSolverCL<MLSmootherCL<GSsmoothCL>, CoarseSolverT, ProlongationT> MGSolversymmGST;
    MGSolversymmGST MGSolversymmGS_;
    typedef MGSolverCL<MLSmootherCL<SGSsmoothCL>, CoarseSolverT, ProlongationT> MGSolversymmSGST;
    MGSolversymmSGST MGSolversymmSGS_;
    typedef MGSolverCL<MLSmootherCL<SORsmoothCL>, CoarseSolverT, ProlongationT> MGSolversymmSORT;
    MGSolversymmSORT MGSolversymmSOR_;
    typedef MGSolverCL<MLSmootherCL<SSORsmoothCL>, CoarseSolverT, ProlongationT> MGSolversymmSSORT;
    MGSolversymmSSORT MGSolversymmSSOR_;
    typedef MGSolverCL<MLSmootherCL<ChebyshevsmoothCL>, CoarseSolverT, ProlongationT> MGSolversymmChebychevT;
    MGSolversymmChebychevT MGSolversymmChebychev_;

    //JAC-GMRes
    typedef GMResSolverCL<JACPcCL> GMResSolverT;
    GMResSolverT GMResSolver_;
    typedef GMResSolverCL<SSORPcCL> GMResSolverSSORT;
    GMResSolverSSORT GMResSolverSSOR_;

    //PCG
    typedef PCGSolverCL<JACPcCL>       PCGSolverJACT;
    PCGSolverJACT PCGSolverJAC_;
    typedef PCGSolverCL<SSORPcCL>      PCGSolverSSORT;
    PCGSolverSSORT PCGSolverSSOR_;
    typedef PCGSolverCL<ChebyshevPcCL> PCGSolverChebychevT;
    PCGSolverChebychevT PCGSolverChebychev_;

  public:
    PoissonSolverFactoryCL( const ParamCL& P, MLIdxDescCL& idx);
    ~PoissonSolverFactoryCL() {}

    /// Returns pointer to prolongation for velocity
    ProlongationT* GetProlongation();
    PoissonSolverBaseCL* CreatePoissonSolver();

};


template <class ProlongationT>
PoissonSolverFactoryCL<ProlongationT>::
    PoissonSolverFactoryCL(const ParamCL& P, MLIdxDescCL& idx)
    : P_(P), idx_(idx), prolongptr_( 0), JACPc_( P.get<double>("Relax")), SSORPc_( P.get<double>("Relax")), ChebyPc_( P.get<double>("Relax")),
        jorsmoother_( P.get<double>("Relax")), gssmoother_( P.get<double>("Relax")), sgssmoother_( P.get<double>("Relax")), sorsmoother_( P.get<double>("Relax")), ssorsmoother_( P.get<double>("Relax")), chebysmoother_( P.get<double>("Relax")),
        coarsesolversymm_( genpc_, 500, 1e-14, true),
        MGSolversymmJOR_( jorsmoother_, coarsesolversymm_, P.get<int>("Iter"), P.get<double>("Tol"), idx, P.get<double>("useRelTol"), P.get<int>("MG.SmoothingSteps"), P.get<int>("MG.NumLvl")),
        MGSolversymmGS_( gssmoother_, coarsesolversymm_, P.get<int>("Iter"), P.get<double>("Tol"), idx, P.get<double>("useRelTol"), P.get<int>("MG.SmoothingSteps"), P.get<int>("MG.NumLvl")),
        MGSolversymmSGS_( sgssmoother_, coarsesolversymm_, P.get<int>("Iter"), P.get<double>("Tol"), idx, P.get<double>("useRelTol"), P.get<int>("MG.SmoothingSteps"), P.get<int>("MG.NumLvl")),
        MGSolversymmSOR_( sorsmoother_, coarsesolversymm_, P.get<int>("Iter"), P.get<double>("Tol"), idx, P.get<double>("useRelTol"), P.get<int>("MG.SmoothingSteps"), P.get<int>("MG.NumLvl")),
        MGSolversymmSSOR_( ssorsmoother_, coarsesolversymm_, P.get<int>("Iter"), P.get<double>("Tol"), idx, P.get<double>("useRelTol"), P.get<int>("MG.SmoothingSteps"), P.get<int>("MG.NumLvl")),
        MGSolversymmChebychev_( chebysmoother_, coarsesolversymm_, P.get<int>("Iter"), P.get<double>("Tol"), idx, P.get<double>("useRelTol"), P.get<int>("MG.SmoothingSteps"), P.get<int>("MG.NumLvl")),
        GMResSolver_( JACPc_, P.get<int>("Restart"), P.get<int>("Iter"), P.get<double>("Tol"), P.get<double>("useRelTol")),
        GMResSolverSSOR_( SSORPc_, P.get<int>("Restart"), P.get<int>("Iter"), P.get<double>("Tol"), P.get<double>("useRelTol")),
        PCGSolverJAC_( JACPc_, P.get<int>("Iter"), P.get<double>("Tol"), P.get<double>("useRelTol")),
        PCGSolverSSOR_( SSORPc_, P.get<int>("Iter"), P.get<double>("Tol"), P.get<double>("useRelTol")),
        PCGSolverChebychev_( ChebyPc_, P.get<int>("Iter"), P.get<double>("Tol"), P.get<double>("useRelTol"))
        {}

template <class ProlongationT>
PoissonSolverBaseCL* PoissonSolverFactoryCL<ProlongationT>::CreatePoissonSolver()
{
    PoissonSolverBaseCL* Poissonsolver = 0;
    switch (P_.get<int>("Solver"))
    {
        case  102 : {
            Poissonsolver = new PoissonSolverCL<MGSolversymmJORT>( MGSolversymmJOR_);
            prolongptr_ = MGSolversymmJOR_.GetProlongation();
        } break;
        case  103 : {
            Poissonsolver = new PoissonSolverCL<MGSolversymmSSORT>( MGSolversymmSSOR_);
            prolongptr_ = MGSolversymmSSOR_.GetProlongation();
        } break;
        case  104 : {
            Poissonsolver = new PoissonSolverCL<MGSolversymmGST>( MGSolversymmGS_);
            prolongptr_ = MGSolversymmGS_.GetProlongation();
        } break;
        case  105 : {
            Poissonsolver = new PoissonSolverCL<MGSolversymmSGST>( MGSolversymmSGS_);
            prolongptr_ = MGSolversymmSGS_.GetProlongation();
        } break;
        case  106 : {
            Poissonsolver = new PoissonSolverCL<MGSolversymmSORT>( MGSolversymmSOR_);
            prolongptr_ = MGSolversymmSOR_.GetProlongation();
        } break;
        case  107 : {
            Poissonsolver = new PoissonSolverCL<MGSolversymmChebychevT>( MGSolversymmChebychev_);
            prolongptr_ = MGSolversymmChebychev_.GetProlongation();
        } break;
        case  302 : Poissonsolver = new PoissonSolverCL<GMResSolverT>( GMResSolver_);  break;
        case  303 : Poissonsolver = new PoissonSolverCL<GMResSolverSSORT>( GMResSolverSSOR_);  break;
        case  202 : Poissonsolver = new PoissonSolverCL<PCGSolverJACT>( PCGSolverJAC_); break;
        case  203 : Poissonsolver = new PoissonSolverCL<PCGSolverSSORT>( PCGSolverSSOR_); break;
        case  207 : Poissonsolver = new PoissonSolverCL<PCGSolverChebychevT>( PCGSolverChebychev_); break;
        default: throw DROPSErrCL("PoissonSolverFactoryCL: Unknown Poisson solver");
    }
    return Poissonsolver;
}

template <class ProlongationT>
ProlongationT* PoissonSolverFactoryCL<ProlongationT>::GetProlongation()
{
    return prolongptr_;
}

} //end of namespace DROPS

#endif /* PoissonSOLVERFACTORY_H_ */
