#include "num/stokespardiso.h"
#ifdef DROPS_PARDISO
#include "num/pardisosolver.h"
#endif
namespace DROPS
{

//=============================================================================
//  The Stokes solvers solve systems of the form
//    A v + BT p = b
//    B v        = c
//=============================================================================

void StokesPardisoSolverCL::Solve( const MatrixCL& A, const MatrixCL& B, const MatrixCL& C, VectorCL& v, VectorCL& p,
                                   const VectorCL& b, const VectorCL& c, const DummyExchangeCL&, const DummyExchangeCL&){

    size_t ndof_v = v.size();
    size_t ndof_p = p.size();
    size_t ndof = ndof_v + ndof_p;
    VectorCL rhs(ndof);
    for (size_t i = 0; i < ndof_v; i++)
        rhs[i] = b[i];
    for (size_t i = 0; i < ndof_p; i++)
        rhs[ndof_v + i] = c[i];
    VectorCL sol(ndof);

    BlockMatrixCL BlockM( &A, MUL, &B, TRANSP_MUL, &B, MUL, &C, MUL);

    MatrixCL M(BuildMatrix( BlockM ));

#ifdef DROPS_PARDISO
    {
        DROPS::PardisoSolverCL SolveA( M);
        SolveA.Solve(M, sol, rhs);
        DROPS::VectorCL r(M * sol - rhs);
        const double normr = norm(r);
        std::cout << "PARDISO 1 Residual absolute:  "<< normr << std::endl;
        std::cout << "PARDISO 1 Residual relative:  "<< normr/norm(rhs) << std::endl;
    }
#else
    throw DROPSErrCL("StokesPardisoSolverCL called, but PARDISO-Flag is not active");
#endif

    for (size_t i = 0; i < ndof_v; i++)
        v[i] = sol[i];
    for (size_t i = 0; i < ndof_p; i++)
        p[i] = sol[ndof_v +i];
}


} //end of namespace
