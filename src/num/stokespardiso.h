#ifndef DROPS_STOKESPARDISOSOLVER_H
#define DROPS_STOKESPARDISOSOLVER_H

#include "num/oseensolver.h"
#ifdef _PAR
#  include "parallel/exchange.h"
#endif
namespace DROPS
{


//=============================================================================
//  The Stokes solvers solve systems of the form
//    A v + BT p = b
//    B v        = c
//=============================================================================

// What every iterative stokes-solver should have
class StokesPardisoSolverCL: public StokesSolverBaseCL
{
  public:
    StokesPardisoSolverCL (std::ostream* output= 0)
        : StokesSolverBaseCL(1, 0, false, output){}

#ifdef _PAR
    virtual void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                        const VectorCL& b, const VectorCL& c, const ExchangeCL& vel_ex, const ExchangeCL& pr_ex){
        throw DROPSErrCL("No parallel pardiso...");
    }
    virtual void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                        const VectorCL& b, const VectorCL& c, const ExchangeCL& vel_ex, const ExchangeCL& pr_ex){
        throw DROPSErrCL("No parallel pardiso...");
    }
#endif
    virtual void Solve( const MatrixCL& A, const MatrixCL& B, const MatrixCL* C, VectorCL& v, VectorCL& p,
                        const VectorCL& b, const VectorCL& c, const DummyExchangeCL& vel_ex, const DummyExchangeCL& pr_ex);

    virtual void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                        const VectorCL& b, const VectorCL& c, const DummyExchangeCL& vel_ex, const DummyExchangeCL& pr_ex){
        Solve(A,B,NULL,v,p,b,c,vel_ex,pr_ex);       
    }

    virtual void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                        const VectorCL& b, const VectorCL& c, const DummyExchangeCL& vel_ex, const DummyExchangeCL& pr_ex){
         Solve(A.GetFinest(), B.GetFinest(), NULL, v, p, b, c, vel_ex, pr_ex);
    }

    virtual void Solve( const MatrixCL& A, const MatrixCL& B, const MatrixCL& C, VectorCL& v, VectorCL& p,
                        const VectorCL& b, const VectorCL& c, const DummyExchangeCL& vel_ex, const DummyExchangeCL& pr_ex)
    {
         Solve(A,B,&C,v,p,b,c,vel_ex,pr_ex);       
    }


    virtual void Solve( const MLMatrixCL& A, const MLMatrixCL& B, const MLMatrixCL& C, VectorCL& a, VectorCL& b,
                        const VectorCL& c, const VectorCL& d, const DummyExchangeCL& da, const DummyExchangeCL&db)
    {
        Solve(A.GetFinest(), B.GetFinest(), &(C.GetFinest()), a, b, c, d, da, db);        
    }


};


} //end of namespace
#endif
