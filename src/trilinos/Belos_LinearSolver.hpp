//
// Created by Alexander Zhiliakov alex@math.uh.edu on 4/8/21.
//

#ifndef DROPS_BELOS_LINEARSOLVER_HPP
#define DROPS_BELOS_LINEARSOLVER_HPP

#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "SingletonLogger.hpp"

namespace DROPS {

    struct Belos_LinearSolver {
        using SV = Epetra_Vector;
        using MV = Epetra_MultiVector;
        using OP = Epetra_Operator;
        using MT = Epetra_CrsMatrix;
        using ST = double;
        struct {
            Teuchos::RCP<SV> lhs, rhs;
            Teuchos::RCP<OP> mtx, pre;
        } system;
        Teuchos::RCP<Belos::SolverManager<ST, MV, OP>> solver;
        std::function<void()> updatePreconditioner;
        bool zeroIniGuess = false;
        bool mute = false;
        struct Stats {
            int problemSize = -1;
            Belos::ReturnType result = Belos::Converged;
            struct { double r_0 = 0., r_i = 0., b = 0.; } norm;
            struct { double solve = 0., updatePreconditioner = 0.; } time;
        };
        const Stats stats;
        void solve() {
            using namespace Teuchos;
            using namespace Belos;
            std::string funcName = __func__;
            if (!system.mtx || !system.lhs || !system.rhs) throw std::invalid_argument(funcName + ": setup system mtx, lhs, and rhs");
            if (system.lhs->GlobalLength() != system.rhs->GlobalLength()) throw std::invalid_argument(funcName + ": lhs and rhs sizes are inconsistent");
            auto& logger = SingletonLogger::instance();
            std::swap(logger.mute, mute);
            auto& s = const_cast<Stats&>(stats);
            auto updated = false;
            if (updatePreconditioner && s.problemSize != system.lhs->GlobalLength()) {
                if (s.problemSize == -1) logger.beg("build initial preconditioner");
                else logger.beg("update preconditioner (problem size changed)");
                    updatePreconditioner();
                s.time.updatePreconditioner = logger.end();
                s.problemSize = system.lhs->GlobalLength();
                updated = true;
            }
            else s.time.updatePreconditioner = 0.;
            logger.beg("linear solve");
                if (zeroIniGuess) system.lhs->PutScalar(0.);
                system.rhs->Norm2(&s.norm.b);
                auto residualNorm = [&]() {
                    SV r(system.rhs->Map(), true);
                    system.mtx->Apply(*system.lhs, r);
                    r.Update(-1., *system.rhs, 1.);
                    double nrm;
                    r.Norm2(&nrm);
                    return nrm;
                };
                s.norm.r_0 = residualNorm();
                auto params = rcp_const_cast<ParameterList>(solver->getCurrentParameters());
                auto tol = params->get<double>("Convergence Tolerance");
                params->set("Convergence Tolerance", s.norm.r_0 ? tol * s.norm.b / s.norm.r_0 : 0.);
                solver->setParameters(params);
                logger.buf
                    << "b norm: " << s.norm.b << '\n'
                    << "r_0 norm: " << s.norm.r_0 << '\n'
                    << "tol: " << tol;
                logger.log();
                LinearProblem<ST, MV, OP> problem(system.mtx, system.lhs, system.rhs);
                if (system.pre) problem.setRightPrec(system.pre);
                problem.setProblem();
                solver->setProblem(rcpFromRef(problem));
                s.result = solver->solve();
            s.time.solve = logger.end();
            if (s.result != Converged) {
                logger.wrn("belos did not converge");
                if (updatePreconditioner && !updated) {
                    logger.beg("update preconditioner");
                        updatePreconditioner();
                    s.time.updatePreconditioner += logger.end();
                    logger.beg("linear solve w/ updated preconditioner");
                        if (zeroIniGuess) system.lhs->PutScalar(0.);
                        else {
                            auto r_0 = residualNorm();
                            params->set("Convergence Tolerance", r_0 ? tol * s.norm.b / r_0 : 0.);
                            solver->setParameters(params);
                        }
                        LinearProblem<ST, MV, OP> problem(system.mtx, system.lhs, system.rhs);
                        if (system.pre) problem.setRightPrec(system.pre);
                        problem.setProblem();
                        solver->setProblem(rcpFromRef(problem));
                        s.result = solver->solve();
                        if (s.result != Converged) logger.wrn("belos did not converge");
                    s.time.solve += logger.end();
                }
            }
            s.norm.r_i = residualNorm();
            params->set("Convergence Tolerance", tol);
            solver->setParameters(params);
            std::swap(logger.mute, mute);
        }
    };

}
#endif //DROPS_BELOS_LINEARSOLVER_HPP
