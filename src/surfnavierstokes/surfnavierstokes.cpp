/// \file surfnavierstokes.cpp
/// \brief Trace FEM discretization of a surface Navier-Stokes problem
/// \author Alexander Zhiliakov alex@math.uh.edu

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
 * GNU Lesser General Public License for more dtau_uils.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#include <fstream>
#include "misc/params.h"
#include "geom/builder.h"
#include "levelset/levelset.h"
#include "levelset/adaptriang.h"
#include "levelset/surfacetension.h"
#include "misc/dynamicload.h"
#include "num/bndData.h"
#include "surfactant/ifacetransp.h"
#include "out/VTKWriter.hpp"
#include "SurfNavierStokesData.hpp"
#include "SingletonLogger.hpp"
// belos (iterative solvers)
#include "BelosSolverFactory.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
// #include "BelosStatusTestImpResNorm.hpp"
// amesos2 (sparse-direct solvers)
#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"
// epetra (vectors and operators / matrices)
#include "epetra/Epetra_OperatorApply.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#ifdef HAVE_MPI
    #include "Epetra_MpiComm.h"
#else
    #include "Epetra_SerialComm.h"
#endif

using namespace DROPS;

int main(int argc, char* argv[]) {
    auto& logger = SingletonLogger::instance();
    try {
        logger.beg("mpi init");
            #ifdef HAVE_MPI
                MPI_Init(&argc, &argv);
                Epetra_MpiComm comm(MPI_COMM_WORLD);
                logger.log("using Epetra_MpiComm");
            #else
                Epetra_SerialComm comm;
                logger.log("using Epetra_SerialComm");
            #endif
            auto myRank = comm.MyPID();
            if (myRank == 0) {
                logger.buf << "num proc = " << comm.NumProc();
                logger.log();
            }
        logger.end();
        std::cout << std::scientific;
        logger.beg("read input .json");
            ParamCL inpJSON;
            read_parameter_file_from_cmdline(inpJSON, argc, argv, "../../param/surfnavierstokes/No_Bnd_Condition.json");
            dynamicLoad(inpJSON.get<std::string>("General.DynamicLibsPrefix"), inpJSON.get<std::vector<std::string>>("General.DynamicLibs"));
            auto dirName = inpJSON.get<std::string>("Output.Directory");
            system((
                "mkdir -p " + dirName + " ; rm -r -f " + dirName + "/* ; "
                "mkdir " + dirName + "/stats ; mkdir " + dirName + "/matrices ; mkdir " + dirName + "/vectors ; mkdir " + dirName + "/vtk"
            ).c_str());
            {
                std::ofstream input(dirName + "/input.json");
                input << inpJSON;
                logger.buf << inpJSON;
                logger.log();
            }
            auto h = inpJSON.get<Point3DCL>("Mesh.E1")[0] / inpJSON.get<double>("Mesh.N1") * std::pow(2., -inpJSON.get<double>("Mesh.AdaptRef.FinestLevel"));
            auto tau_u = inpJSON.get<double>("SurfNavierStokes.NormalPenalty.Scaling") * pow(h, inpJSON.get<double>("SurfNavierStokes.NormalPenalty.Power")); // constant for normal penalty
            auto rho_u = inpJSON.get<double>("SurfNavierStokes.VelocityStab.Scaling")  * pow(h, inpJSON.get<double>("SurfNavierStokes.VelocityStab.Power")); // constant for velocity stabilisation
            auto rho_p = inpJSON.get<double>("SurfNavierStokes.PressureStab.Scaling")  * pow(h, inpJSON.get<double>("SurfNavierStokes.PressureStab.Power")); // constant for pressure stabilisation
            auto numSteps = inpJSON.get<size_t>("Time.NumbOfSteps");
            auto finalTime = inpJSON.get<double>("Time.FinalTime");
            auto stepSize = finalTime / numSteps;
            auto everyStep = inpJSON.get<int>("Output.EveryStep");
            auto exportMatrices = inpJSON.get<bool>("Output.Matrices");
            auto binary = inpJSON.get<bool>("Output.Binary");
            auto mtxFormat = binary ? ".mat" : ".mtx";
            auto mtxExpFunc = binary ? &MatrixCL::exportMAT : &MatrixCL::exportMTX;
            auto exportVectors = inpJSON.get<bool>("Output.Vectors");
            auto surface = surfaceFactory(inpJSON);
            auto testName = inpJSON.get<std::string>("SurfNavierStokes.IC.Name");
            auto surfNavierStokesData = surfNavierStokesDataFactory(*surface, testName, inpJSON);
            FESystem surfNSystem;
            auto& gamma = surfNSystem.params.surfNavierStokesParams.gamma;
            gamma = inpJSON.get<double>("SurfNavierStokes.gamma");
            auto& nu = surfNSystem.params.surfNavierStokesParams.nu;
            nu = inpJSON.get<double>("SurfNavierStokes.nu");
            if (nu <= 0.) throw std::invalid_argument("viscosity must be non-negative");
            auto& t = surfNSystem.params.t;
            auto& w_T = surfNSystem.params.surfNavierStokesParams.w_T;
            auto& Pe = surfNSystem.params.surfNavierStokesParams.Pe;
            auto& u_N = surfNSystem.params.surfNavierStokesParams.u_N;
            auto& u_N_max = surfNSystem.params.surfNavierStokesParams.u_N_max;
            auto& levelSet = surfNSystem.params.levelSet;
            surfNSystem.params.numbOfVirtualSubEdges = inpJSON.get<size_t>("SurfNavierStokes.NumbOfVirtualSubEdges");
            surfNSystem.params.surfNavierStokesParams.m_g = surfNavierStokesData.m_g;
            surfNSystem.params.surfNavierStokesParams.f_T = surfNavierStokesData.f_T;
            auto formulation = inpJSON.get<std::string>("SurfNavierStokes.Formulation");
            auto stab = inpJSON.get<std::string>("SurfNavierStokes.PressureStab.Type");
            auto useTangMassMat = inpJSON.get<bool>("SurfNavierStokes.UseTangentialMassMatrix");
            auto useInnerIters = inpJSON.get<bool>("Solver.Inner.Use");
            auto usePrevGuess = inpJSON.get<bool>("Solver.UsePreviousFrameAsInitialGuess");
            auto BDF = inpJSON.get<size_t>("SurfNavierStokes.BDF");
            if (BDF != 1 && BDF != 2) throw std::invalid_argument("use BDF = 1 or 2");
            auto tol = inpJSON.get<double>("Solver.Outer.RelResTol");
            logger.buf
                << surface->description() << '\n'
                << surfNavierStokesData.description << '\n'
                << "exact soln: " << (surfNavierStokesData.exact ? "yes" : "no") << '\n'
                << "delta t = " << stepSize;
            logger.log();
        logger.end();
        logger.beg("build mesh");
            logger.beg("build initial bulk mesh");
                std::auto_ptr<MGBuilderCL> builder(make_MGBuilder(inpJSON));
                MultiGridCL mg(*builder);
                const ParamCL::ptree_type* ch= 0;
                try {
                    ch = &inpJSON.get_child("Mesh.Periodicity");
                } catch (DROPSParamErrCL) {}
                if (ch) read_PeriodicBoundaries(mg, *ch);
                // levelset shift
                auto shift = inpJSON.get<Point3DCL>("Levelset.ShiftDir", Point3DCL(0., 0., 0.));
                shift /= shift.norm();
                auto shiftNorm = fabs(inpJSON.get<double>("Levelset.ShiftNorm", 0.));
                shift *= shiftNorm;
                mg.Transform([&](Point3DCL const & p) { return p - shift; });
                logger.buf << "surface shift: " << shift;
                logger.log();
            logger.end();
            logger.beg("refine towards the surface");
                AdapTriangCL adap(mg);
                DistMarkingStrategyCL markerLset([&](Point3DCL const & x, double) { return surface->phi(x); }, inpJSON.get<double>("Mesh.AdaptRef.Width"), inpJSON.get<int>("Mesh.AdaptRef.CoarsestLevel"), inpJSON.get<int>("Mesh.AdaptRef.FinestLevel"));
                adap.set_marking_strategy(&markerLset);
                adap.MakeInitialTriang();
                adap.set_marking_strategy(nullptr);
                auto numbOfTetras = mg.GetNumTriangTetra();
                logger.buf
                    << "h = " << h << '\n'
                    << "numb of tetras = " << numbOfTetras;
                logger.log();
            logger.end();
            BndDataCL<double> lsetBnd(0); read_BndData(lsetBnd, mg, inpJSON.get_child("Levelset.BndData"));
            BndDataCL<Point3DCL> vecBnd(0); read_BndData(vecBnd, mg, inpJSON.get_child("Stokes.VelocityBndData"));
            BndDataCL<double> scaBnd(0); read_BndData(scaBnd, mg, inpJSON.get_child("Stokes.PressureBndData"));
        logger.end();
        logger.beg("interpolate level-set");
            MLIdxDescCL lstIdx(P2_FE);
            lstIdx.CreateNumbering(mg.GetLastLevel(), mg, lsetBnd);
            levelSet.SetIdx(&lstIdx);
            levelSet.Interpolate(mg, [&](Point3DCL const & x) { return surface->phi(x); });
        logger.end();
        logger.beg("set up FE spaces");
            IdxDescCL velIdx(vecP2IF_FE, vecBnd); {
                velIdx.GetXidx().SetBound(inpJSON.get<double>("SurfTransp.OmitBound"));
                auto numbOfActiveTetras = velIdx.CreateNumbering(mg.GetLastLevel(), mg, &levelSet, &lsetBnd);
                logger.buf
                    << "numb of active (cut) tetras is: " << numbOfActiveTetras << " ("
                    << (100. * numbOfActiveTetras) / numbOfTetras << "%)\n";
                logger.log();
            }
            IdxDescCL preIdx(P1IF_FE, scaBnd); {
                preIdx.GetXidx().SetBound(inpJSON.get<double>("SurfTransp.OmitBound"));
                preIdx.CreateNumbering(mg.GetLastLevel(), mg, &levelSet, &lsetBnd);
            }
            IdxDescCL speedIdx(P2IF_FE, scaBnd); {
                speedIdx.GetXidx().SetBound(inpJSON.get<double>("SurfTransp.OmitBound"));
                speedIdx.CreateNumbering(mg.GetLastLevel(), mg, &levelSet, &lsetBnd);
            }
            size_t n = velIdx.NumUnknowns();
            auto n_i = n / 3;
            size_t m = preIdx.NumUnknowns();
            logger.beg("set up FE system");
                u_N.SetIdx(&speedIdx);
                w_T.SetIdx(&velIdx);
                FEMatDescCL
                    M_u(&velIdx, &velIdx, useTangMassMat ? &LocalAssembler::M_t_vecP2vecP2 : &LocalAssembler::M_vecP2vecP2),
                    AL_u(&velIdx, &velIdx, &LocalAssembler::AL_vecP2vecP2),
                    N_u(&velIdx, &velIdx, &LocalAssembler::N_vecP2vecP2),
                    H_u(&velIdx, &velIdx, &LocalAssembler::H_vecP2vecP2),
                    A_u(&velIdx, &velIdx, formulation == "Consistent" ? &LocalAssembler::A_consistent_vecP2vecP2 : &LocalAssembler::A_vecP2vecP2),
                    S_u(&velIdx, &velIdx, &LocalAssembler::S_vecP2vecP2),
                    C_u(&velIdx, &velIdx, &LocalAssembler::C_n_vecP2vecP2),
                    B_pu(&preIdx, &velIdx, &LocalAssembler::B_P1vecP2),
                    Q_pu(&preIdx, &velIdx, &LocalAssembler::Q_P1vecP2),
                    M_p(&preIdx, &preIdx, &LocalAssembler::M_P1P1),
                    C_p(&preIdx, &preIdx, stab == "Normal" ? &LocalAssembler::C_n_P1P1 : &LocalAssembler::C_full_P1P1),
                    A_p(&preIdx, &preIdx, &LocalAssembler::A_P1P1);
                FEVecDescCL
                    F_u(&velIdx, &LocalAssembler::F_momentum_vecP2),
                    G_p(&preIdx, &LocalAssembler::F_continuity_P1);
                VecDescCL
                    u_star(&velIdx), u(&velIdx), u_prev(&velIdx), surf_curl_u(&preIdx),
                    p_star(&preIdx), p(&preIdx);
            logger.end();
        logger.end();
        logger.beg("set up vtk");
            VTKWriter vtkWriter(dirName + "/vtk/" + testName + inpJSON.get<std::string>("Surface.Name"), mg, binary);
            vtkWriter.add({ "level-set", levelSet });
            if (everyStep > 0) {
                if (inpJSON.get<bool>("Output.VTK.Velocity")) {
                    vtkWriter.add({ "u_h", u });
                    if (surfNavierStokesData.exact) vtkWriter.add({"u_*", u_star });
                }
                if (inpJSON.get<bool>("Output.VTK.SurfSpeed")) vtkWriter.add({ "u_N", u_N });
                if (inpJSON.get<bool>("Output.VTK.Vorticity")) vtkWriter.add({ "w_h", surf_curl_u });
                if (inpJSON.get<bool>("Output.VTK.Pressure")) {
                    vtkWriter.add({ "p_h", p });
                    if (surfNavierStokesData.exact) vtkWriter.add({"p_*", p_star });
                }
            }
            auto exportVTK = vtkWriter.numVars() > 0;
        logger.end();
        logger.beg("set up linear solver");
            using SV = Epetra_Vector;
            using MV = Epetra_MultiVector;
            using OP = Epetra_Operator;
            using MT = Epetra_CrsMatrix;
            using ST = double;
            using namespace Teuchos;
            using namespace Belos;
            auto belosParams = parameterList();
            belosParams->set("Num Blocks", inpJSON.get<int>("Solver.Outer.KrylovSubspaceSize"));
            belosParams->set("Maximum Iterations", inpJSON.get<int>("Solver.Outer.MaxIter"));
            belosParams->set("Output Frequency", inpJSON.get<int>("Solver.Outer.OutputFrequency"));
            belosParams->set("Verbosity", Errors + Warnings + StatusTestDetails + TimingDetails + FinalSummary + IterationDetails);
            belosParams->set<int>("Output Style", Brief);
            SolverFactory<ST, MV, OP> belosFactory;
            logger.buf << "available iterations: " << belosFactory.supportedSolverNames();
            logger.log();
            auto belosSolver = belosFactory.create(inpJSON.get<std::string>("Solver.Outer.Iteration"), belosParams);
            logger.buf << "outer solver: " << belosSolver->description();
            logger.log();
            Epetra_Map mapVelocity(static_cast<int>(n), 0, comm), mapVelocityComp(static_cast<int>(n_i), 0, comm), mapVelocityPressure(static_cast<int>(n + m), 0, comm), mapPressure(static_cast<int>(m), 0, comm);
            SV belosLHS(mapVelocityPressure), belosRHS(mapVelocityPressure);
            RCP<MT> A, B, C;
            RCP<OP> belosMTX = rcp(new Epetra_OperatorApply([&](MV const & X, MV& Y) {
                double* view;
                X(0)->ExtractView(&view);
                SV x1(Epetra_DataAccess::View, mapVelocity, view);
                SV x2(Epetra_DataAccess::View, mapPressure, view + n);
                Y(0)->ExtractView(&view);
                SV y11(Epetra_DataAccess::View, mapVelocity, view), y12(mapVelocity, true);
                SV y21(Epetra_DataAccess::View, mapPressure, view + n), y22(mapPressure, true);
                // y1
                A->Multiply(false, x1, y11);
                B->Multiply(true, x2, y12);
                y11.Update(1., y12, 1.);
                // y2
                B->Multiply(false, x1, y21);
                C->Multiply(false, x2, y22);
                y21.Update(1., y22, 1.);
            }));
            auto belosRES = [&]() {
                SV residual(mapVelocityPressure);
                belosMTX->Apply(belosLHS, residual);
                residual.Update(-1., belosRHS, 1.);
                double nrm;
                residual.Norm2(&nrm);
                return nrm;
            };
            logger.beg("set up preconditioners");
                logger.beg("diffusion-convection-reaction block");
                    auto iterationA = inpJSON.get<std::string>("Solver.Inner.A.Iteration");
                    if (!Amesos2::query(iterationA)) throw std::invalid_argument("solver " + iterationA + " is not available for Amesos2");
                    logger.log("iteration: " + iterationA);
                    auto precTypeA = inpJSON.get<std::string>("Solver.Inner.A.Type");
                    logger.log("preconditioner type: " + precTypeA);
                    size_t numBlocksA = precTypeA == "Full" ? 1 : 3;
                    auto& mapVelocityBlock = numBlocksA == 1 ? mapVelocity : mapVelocityComp;
                    if (numBlocksA == 1) logger.log("factorizing full velocity matrix");
                    else logger.log("factorizing diagonal blocks A_{ii} of velocity matrix");
                    std::vector<std::vector<RCP<MT>>> A_block(numBlocksA, std::vector<RCP<MT>>(numBlocksA));
                    std::vector<RCP<Amesos2::Solver<MT, MV>>> amesosSolver(numBlocksA);
                    auto runFactorization = [&]() {
                        for (size_t i = 0; i < numBlocksA; ++i)
                            amesosSolver[i] = Amesos2::create<MT, MV>(iterationA, A_block[i][i]);
                        logger.buf
                            << "solver: " << amesosSolver[0]->name() << '\n'
                            << "description: " << amesosSolver[0]->description();
                        logger.log();
                        #pragma omp parallel for
                        for (size_t i = 0; i < numBlocksA; ++i) {
                            amesosSolver[i]->symbolicFactorization();
                            amesosSolver[i]->numericFactorization();
                        }
                        auto amesosStatus = amesosSolver[0]->getStatus();
                        logger.buf << "numb of nonzeros in L + U = " << amesosStatus.getNnzLU() << " ("
                                   << (100. * amesosStatus.getNnzLU()) / (static_cast<double>(A_block[0][0]->NumGlobalRows()) * A_block[0][0]->NumGlobalCols()) << "%)";
                        logger.log();
                    };
                    size_t numItersA = 0;
                    Epetra_OperatorApply::ApplyType invA = [&](MV const &X, MV &Y) {
                        double *viewX, *viewY;
                        X(0)->ExtractView(&viewX);
                        Y(0)->ExtractView(&viewY);
                        std::vector<SV> x, y;
                        x.reserve(numBlocksA);
                        y.reserve(numBlocksA);
                        for (size_t i = 0; i < numBlocksA; ++i) {
                            x.emplace_back(Epetra_DataAccess::View, mapVelocityBlock, viewX + i * n_i);
                            y.emplace_back(Epetra_DataAccess::View, mapVelocityBlock, viewY + i * n_i);
                        }
                        #pragma omp parallel for
                        for (size_t i = 0; i < numBlocksA; ++i) {
                            amesosSolver[i]->setB(rcpFromRef(x[i]));
                            amesosSolver[i]->setX(rcpFromRef(y[i]));
                            amesosSolver[i]->solve();
                        }
                        numItersA++;
                    };
                    if (precTypeA == "BlockTriangular")
                        invA = [&](MV const &X, MV &Y) {
                            double *viewX, *viewY;
                            X(0)->ExtractView(&viewX);
                            Y(0)->ExtractView(&viewY);
                            std::vector<SV> x, y;
                            x.reserve(numBlocksA);
                            y.reserve(numBlocksA);
                            for (size_t i = 0; i < numBlocksA; ++i) {
                                x.emplace_back(Epetra_DataAccess::View, mapVelocityBlock, viewX + i * n_i);
                                y.emplace_back(Epetra_DataAccess::View, mapVelocityBlock, viewY + i * n_i);
                            }
                            for (int i = numBlocksA - 1; i >= 0; --i) { // backward substitution
                                SV rhs(Epetra_DataAccess::Copy, x[i], 0);
                                for (size_t j = i + 1; j < numBlocksA; ++j) {
                                    SV A_y(mapVelocityBlock, true);
                                    A_block[i][j]->Multiply(false, y[j], A_y);
                                    rhs.Update(-1., A_y, 1.);
                                }
                                amesosSolver[i]->setB(rcpFromRef(rhs));
                                amesosSolver[i]->setX(rcpFromRef(y[i]));
                                amesosSolver[i]->solve();
                            }
                            numItersA++;
                        };
                logger.end();
                logger.beg("schur complement block");
                    auto belosParamsS = parameterList();
                    belosParamsS->set("Maximum Iterations", inpJSON.get<int>("Solver.Inner.S.MaxIter"));
                    belosParamsS->set("Convergence Tolerance", inpJSON.get<double>("Solver.Inner.S.RelResTol"));
                    auto belosSolverS_M = belosFactory.create(inpJSON.get<std::string>("Solver.Inner.S.Iteration"), belosParamsS);
                    auto belosSolverS_L = belosFactory.create(inpJSON.get<std::string>("Solver.Inner.S.Iteration"), belosParamsS);
                    logger.buf << "inner solver: " << belosSolverS_M->description();
                    logger.log();
                    RCP<MT> S_L, S_M;
                    size_t numItersS_M = 0, numItersS_L = 0;
                    auto invS = [&](MV const & X, MV& Y) {
                        Y.PutScalar(0.); // ini guess
                        auto& Y_M = Y;
                        auto  Y_L = Y;
                        auto X_nrm = X; { // normalized rhs
                            double mean;
                            X_nrm.MeanValue(&mean);
                            for (size_t i = 0; i < m; ++i) (*X_nrm(0))[i] -= mean;
                        }
                        // Y_M
                        LinearProblem<ST, MV, OP> belosProblemS_M(S_M, rcpFromRef(Y_M), rcpFromRef(X_nrm));
                        belosProblemS_M.setProblem();
                        belosSolverS_M->setProblem(rcpFromRef(belosProblemS_M));
                        belosSolverS_M->solve();
                        numItersS_M += belosSolverS_M->getNumIters();
                        // Y_L (TODO: make parallel)
                        if (inpJSON.get<bool>("Solver.Inner.S.S_L")) {
                            LinearProblem<ST, MV, OP> belosProblemS_L(S_L, rcpFromRef(Y_L), rcpFromRef(X_nrm));
                            belosProblemS_L.setProblem();
                            belosSolverS_L->setProblem(rcpFromRef(belosProblemS_L));
                            belosSolverS_L->solve();
                            numItersS_L += belosSolverS_L->getNumIters();
                        }
                        // result
                        Y_M.Update(1., Y_L, 1.);
                    };
                logger.end();
                auto precType = inpJSON.get<std::string>("Solver.Inner.Type");
                logger.log("preconditioner type: " + precType);
                auto belosPRE = rcp(new Epetra_OperatorApply([&](MV const & X, MV& Y) {
                    double* view;
                    X(0)->ExtractView(&view);
                    SV x1(Epetra_DataAccess::View, mapVelocity, view);
                    SV x2(Epetra_DataAccess::View, mapPressure, view + n);
                    Y(0)->ExtractView(&view);
                    SV y1(Epetra_DataAccess::View, mapVelocity, view);
                    SV y2(Epetra_DataAccess::View, mapPressure, view + n);
                    #pragma omp parallel
                    {
                        #pragma omp single
                        {
                            #pragma omp task
                            invA(x1, y1); // y1
                            #pragma omp task
                            invS(x2, y2); // y2
                        }
                    }
                }));
                if (precType == "BlockTriangular")
                    belosPRE = rcp(new Epetra_OperatorApply([&](MV const & X, MV& Y) {
                        double* view;
                        X(0)->ExtractView(&view);
                        SV x1(Epetra_DataAccess::View, mapVelocity, view);
                        SV x2(Epetra_DataAccess::View, mapPressure, view + n);
                        Y(0)->ExtractView(&view);
                        SV y1(Epetra_DataAccess::View, mapVelocity, view);
                        SV y2(Epetra_DataAccess::View, mapPressure, view + n);
                        // y2
                        invS(x2, y2);
                        // y1
                        SV rhs(mapVelocity, true);
                        B->Multiply(true, y2, rhs);
                        rhs.Update(1., x1, -1.);
                        invA(rhs, y1);
                    }));
            logger.end();
            LinearProblem<ST, MV, OP> belosProblem(belosMTX, rcpFromRef(belosLHS), rcpFromRef(belosRHS));
            if (useInnerIters) belosProblem.setRightPrec(belosPRE);
        logger.end();
        logger.beg("t = t_0 = 0");
            t = 0.;
            logger.beg("assemble");
                surfNSystem.matrices = { &M_u, &A_u, &AL_u, &S_u, &C_u, &M_p, &C_p, &A_p, &B_pu, &Q_pu };
                setupFESystem(mg, surfNSystem);
                C_p.Data *= -rho_p;
                MatrixCL A_sum;
                VectorCL I_p(1., m);
            auto assembleTime = logger.end();
            auto factorizationTime = 0.;
            logger.beg("interpolate initial data");
                u_N.Interpolate(mg, [&](Point3DCL const & x) { return surface->u_N(x, t); });
                u.Interpolate(mg, [&](Point3DCL const & x) { return surfNavierStokesData.u_T(x, t); });
                p.Interpolate(mg, [&](Point3DCL const & x) { return surfNavierStokesData.p(x, t); });
                p.Data -= dot(M_p.Data * p.Data, I_p) / dot(M_p.Data * I_p, I_p) * I_p;
                u_prev = u;
            auto solveTime = logger.end();
            auto solveWastedTime = 0.;
            logger.beg("project surface vorticity");
                auto belosParamsW = parameterList();
                belosParamsW->set("Maximum Iterations", 1000);
                belosParamsW->set("Convergence Tolerance", 1e-8);
                belosParamsW->set("Output Frequency", 10);
                belosParamsW->set("Verbosity", Errors + Warnings + StatusTestDetails);
                belosParamsW->set<int>("Output Style", Brief);
                auto belosSolverW = belosFactory.create("CG", belosParamsW);
                ReturnType belosSolverResultW;
                auto projectVorticity = [&]() {
                    auto mtx = static_cast<RCP<MT>>(MatrixCL(1., M_p.Data, -1., C_p.Data));
                    auto rhs = static_cast<SV>(Q_pu.Data * u.Data);
                    auto sln = static_cast<SV>(surf_curl_u.Data);
                    sln.PutScalar(0.);
                    LinearProblem<ST, MV, OP> belosProblemW(mtx, rcpFromRef(sln), rcpFromRef(rhs));
                    belosProblemW.setProblem();
                    belosSolverW->setProblem(rcpFromRef(belosProblemW));
                    belosSolverResultW = belosSolverW->solve();
                    if (belosSolverResultW == Converged) logger.log("belos converged");
                    else logger.wrn("belos did not converge");
                    for (size_t i = 0; i < m; ++i) surf_curl_u.Data[i] = sln[i];
                };
                projectVorticity();
            auto projectTime = logger.end();
            logger.beg("output");
                double alpha, b_norm, r_0_norm;
                ReturnType belosSolverResult;
                auto exportStats = [&](size_t i) {
                    std::ofstream stats(dirName + "/stats/t_" + std::to_string(i) + ".json");
                    ParamCL tJSON;
                    tJSON.put("t", t);
                    tJSON.put("h", h);
                    tJSON.put("nu", nu);
                    tJSON.put("gamma", gamma);
                    tJSON.put("DOF.Velocity", n);
                    tJSON.put("DOF.Pressure", m);
                    tJSON.put("MeshDepParams.rho_p", rho_p);
                    tJSON.put("MeshDepParams.rho_u", rho_u);
                    tJSON.put("MeshDepParams.tau_u", tau_u);
                    tJSON.put("ElapsedTime.Assemble", assembleTime);
                    tJSON.put("ElapsedTime.Factorization", factorizationTime);
                    tJSON.put("ElapsedTime.LinearSolve", solveTime);
                    tJSON.put("ElapsedTime.LinearSolveWasted", solveWastedTime);
                    tJSON.put("ElapsedTime.ProjectVorticity", projectTime);
                    tJSON.put("Solver.ProjectVorticity.TotalIters", belosSolverW->getNumIters());
                    tJSON.put("Solver.ProjectVorticity.Converged", belosSolverResultW == Converged);
                    tJSON.put("Solver.ProjectVorticity.ResidualNormBelosRelative", belosSolverW->achievedTol());
                    auto surfArea = dot(I_p, M_p.Data * I_p);
                    tJSON.put("Integral.SurfacaAreaP1", surfArea);
                    tJSON.put("Integral.FESolution.PressureMean", dot(I_p, M_p.Data * p.Data) / surfArea);
                    tJSON.put("Integral.FESolution.PressureL2", sqrt(dot(p.Data, M_p.Data * p.Data)));
                    auto velL2Sq = dot(u.Data, M_u.Data * u.Data);
                    tJSON.put("Integral.FESolution.KineticEnergy", .5 * velL2Sq);
                    tJSON.put("Integral.FESolution.VelocityL2", sqrt(velL2Sq));
                    tJSON.put("Integral.FESolution.VelocityNormalL2", sqrt(dot(u.Data, S_u.Data * u.Data)));
                    tJSON.put("Integral.FESolution.VelocitySurfaceDivergenceL2", sqrt(dot(u.Data, AL_u.Data * u.Data)));
                    auto vorL2Sq = dot(surf_curl_u.Data, M_p.Data * surf_curl_u.Data);
                    tJSON.put("Integral.FESolution.Enstrophy", .5 * vorL2Sq);
                    tJSON.put("Integral.FESolution.SurfaceVorticityL2", sqrt(vorL2Sq));
                    auto vorH1Sq = dot(surf_curl_u.Data, A_p.Data * surf_curl_u.Data);
                    tJSON.put("Integral.FESolution.Palinstrophy", .5 * vorH1Sq);
                    tJSON.put("Integral.FESolution.SurfaceVorticityH1", sqrt(vorH1Sq));
                    if (surfNavierStokesData.exact) {
                        u_star.Interpolate(mg, [&](Point3DCL const & x) { return surfNavierStokesData.u_T(x, t); });
                        p_star.Interpolate(mg, [&](Point3DCL const & x) { return surfNavierStokesData.p(x, t); });
                        p_star.Data -= dot(M_p.Data * p_star.Data, I_p) / dot(M_p.Data * I_p, I_p) * I_p;
                        tJSON.put("Integral.ExactSolution.PressureMean", dot(I_p, M_p.Data * p_star.Data) / surfArea);
                        tJSON.put("Integral.ExactSolution.PressureL2", sqrt(dot(p_star.Data, M_p.Data * p_star.Data)));
                        auto velL2Sq = dot(u_star.Data, M_u.Data * u_star.Data);
                        tJSON.put("Integral.ExactSolution.KineticEnergy", .5 * velL2Sq);
                        tJSON.put("Integral.ExactSolution.VelocityL2", sqrt(velL2Sq));
                        tJSON.put("Integral.ExactSolution.VelocityNormalL2", sqrt(dot(u_star.Data, S_u.Data * u_star.Data)));
                        tJSON.put("Integral.ExactSolution.VelocitySurfaceDivergenceL2", sqrt(dot(u_star.Data, AL_u.Data * u_star.Data)));
                        auto u_diff = u_star.Data - u.Data;
                        tJSON.put("Integral.Error.VelocityL2", sqrt(dot(u_diff, M_u.Data * u_diff)));
                        tJSON.put("Integral.Error.VelocityH1", sqrt(dot(u_diff, A_u.Data * u_diff)));
                        auto p_diff = p_star.Data - p.Data;
                        tJSON.put("Integral.Error.PressureL2", sqrt(dot(p_diff, M_p.Data * p_diff)));
                    }
                    if (i > 0) {
                        auto r_i_norm = belosRES();
                        tJSON.put("Solver.Outer.ResidualNorm.r_i", r_i_norm);
                        tJSON.put("Solver.Outer.ResidualNorm.b", b_norm);
                        tJSON.put("Solver.Outer.ResidualNorm.r_0", r_0_norm);
                        tJSON.put("Solver.Outer.ResidualNorm.r_i/b", r_i_norm / b_norm);
                        tJSON.put("Solver.Outer.ResidualNorm.r_i/r_0", r_i_norm / r_0_norm);
                        tJSON.put("Solver.Outer.ResidualNorm.BelosRelative", belosSolver->achievedTol());
                        tJSON.put("Solver.Outer.TotalIters", belosSolver->getNumIters());
                        tJSON.put("Solver.Outer.Converged", belosSolverResult == Converged);
                        tJSON.put("Solver.Inner.A.TotalIters", numItersA);
                        tJSON.put("Solver.Inner.A.MeanIters", static_cast<double>(numItersA) / belosSolver->getNumIters());
                        tJSON.put("Solver.Inner.S.TotalIters.S_M", numItersS_M);
                        tJSON.put("Solver.Inner.S.TotalIters.S_L", numItersS_L);
                        tJSON.put("Solver.Inner.S.MeanIters.S_M", static_cast<double>(numItersS_M) / belosSolver->getNumIters());
                        tJSON.put("Solver.Inner.S.MeanIters.S_L", static_cast<double>(numItersS_L) / belosSolver->getNumIters());
                        tJSON.put("Peclet", Pe);
                        tJSON.put("MaxSurfaceSpeed", u_N_max);
                        tJSON.put("MassMatrixCoef", alpha);
                    }
                    auto exportFiles = everyStep > 0 && (i % everyStep == 0 || i == numSteps);
                    if (exportVTK && exportFiles) {
                        logger.beg("write vtk");
                            vtkWriter.write(t);
                        auto vtkTime = logger.end();
                        tJSON.put("ElapsedTime.VTK", vtkTime);
                    }
                    if (exportVectors && exportFiles) {
                        logger.beg("export vectors");
                            auto expVec = [&](VecDescCL const & v, std::string const name) {
                                logger.beg("export " + name + " vec");
                                    logger.buf << "size: " << v.Data.size();
                                    logger.log();
                                    auto path = dirName + "/vectors/" + name + '_' + std::to_string(i) + ".dat";
                                    std::ofstream file(path);
                                    v.Write(file, binary);
                                    tJSON.put("Vectors." + name, path);
                                logger.end();
                            };
                            expVec(u, "u");
                            expVec(p, "p");
                        auto vecTime = logger.end();
                        tJSON.put("ElapsedTime.VectorExport", vecTime);
                    }
                    if (exportMatrices && exportFiles) {
                        auto expMat = [&](MatrixCL &A, std::string const name) {
                            logger.beg("export " + name + " mtx");
                                logger.buf << "size: " << A.num_rows() << 'x' << A.num_cols();
                                logger.log();
                                auto path = dirName + "/matrices/" + name + '_' + std::to_string(i) + mtxFormat;
                                (A.*mtxExpFunc)(path);
                                tJSON.put("Matrices." + name, path);
                            logger.end();
                        };
                        logger.beg("export matrices");
                            expMat(M_u.Data, "M_u");
                            if (i > 0) {
                                expMat(A_sum, "A_u");
                                for (size_t i = 0; i < numBlocksA; ++i)
                                    for (size_t j = 0; j < numBlocksA; ++j)
                                        expMat(A_sum.Split(numBlocksA, numBlocksA)[i][j], "A_u_" + std::to_string(i + 1) + std::to_string(j + 1));
                            }
                            // ...
                        auto mtxTime = logger.end();
                        tJSON.put("ElapsedTime.MatrixExport", mtxTime);
                    }
                    stats << tJSON;
                    logger.buf << tJSON;
                    logger.log();;
                };
                exportStats(0);
            logger.end();
        logger.end();
        for (size_t i = 1; i <= numSteps; ++i) {
            t = i * stepSize;
            std::stringstream header; header << "t = t_" << i << " = " << t << " (" << (100. * i) / numSteps << "%)";
            logger.beg(header.str());
                logger.beg("initialize surface speed u_N");
                    u_N.Interpolate(mg, [&](Point3DCL const & x) { return surface->u_N(x, t); });
                logger.end();
                logger.beg("assemble");
                    surfNSystem.matrices = {};
                    if (surfNavierStokesData.w_T) w_T.Interpolate(mg, [&](Point3DCL const & x) { return surfNavierStokesData.w_T(x, t); }); // Oseen (and Stokes) case
                    else w_T.Data = 2. * u.Data - u_prev.Data; // Navier-Stokes case
                    Pe = supnorm(w_T.Data) / nu;
                    if (Pe) {
                        logger.log("assembling convection mtx");
                        surfNSystem.matrices.push_back(&N_u);
                    }
                    u_N_max = supnorm(u_N.Data);
                    if (u_N_max) {
                        logger.log("assembling surf speed mtx");
                        surfNSystem.matrices.push_back(&H_u);
                    }
                    alpha = i == 1 || BDF == 1 ? 1. / stepSize : 1.5 / stepSize;
                    logger.buf
                        << "alpha = " << alpha << '\n'
                        << "Pe = " << Pe << '\n'
                        << "max |u_N| = " << u_N_max;
                    logger.log();
                    surfNSystem.vectors = { &F_u, &G_p };
                    setupFESystem(mg, surfNSystem);
                    // system mtx
                    A_sum.LinComb(alpha, M_u.Data, gamma, AL_u.Data, nu, A_u.Data, tau_u, S_u.Data, rho_u, C_u.Data);
                    if (Pe) A_sum.LinComb(1., MatrixCL(A_sum), 1., N_u.Data);
                    if (u_N_max) A_sum.LinComb(1., MatrixCL(A_sum), 1., H_u.Data);
                    // system rhs
                    if (i == 1 || BDF == 1) F_u.Data += (1. / stepSize) * (M_u.Data * u.Data);
                    else F_u.Data += (2. / stepSize) * (M_u.Data * u.Data) - (.5 / stepSize) * (M_u.Data * u_prev.Data);
                    G_p.Data -= (dot(G_p.Data, I_p) / dot(I_p, I_p)) * I_p;
                    // for the next step
                    u_prev = u;
                assembleTime = logger.end();
                logger.beg("convert to Epetra");
                    A = static_cast<RCP<MT>>(A_sum);
                    logCRS(*A, "A");
                    {
                        auto A_sum_blocks = A_sum.Split(numBlocksA, numBlocksA);
                        for (size_t i = 0; i < numBlocksA; ++i)
                            for (size_t j = 0; j < numBlocksA; ++j)
                                A_block[i][j] = static_cast<RCP<MT>>(A_sum_blocks[i][j]);
                        if (numBlocksA != 1)
                            logCRS(*A_block[0][0], "A_{ij}");
                    }
                    B = static_cast<RCP<MT>>(B_pu.Data);
                    logCRS(*B, "B");
                    C = static_cast<RCP<MT>>(C_p.Data);
                    logCRS(*C, "C := -rho_p (pressure volume stabilization mtx)");
                    S_M = static_cast<RCP<MT>>(MatrixCL(1. / (gamma + nu), M_p.Data, -1., C_p.Data));
                    logCRS(*S_M, "S_M := (nu + gamma)^{-1} M_p - C");
                    S_L = static_cast<RCP<MT>>(MatrixCL(1. / alpha, A_p.Data, -1., C_p.Data));
                    logCRS(*S_L, "S_L := \\alpha^{-1} L_p - C");
                    belosRHS = static_cast<SV>(F_u.Data.append(G_p.Data));
                logger.end();
                factorizationTime = 0.;
                if (useInnerIters && i == 1) {
                    logger.beg("build initial factorization");
                        runFactorization();
                    factorizationTime = logger.end();
                }
                solveWastedTime = 0.;
                logger.beg("linear solve");
                    numItersA = numItersS_M = numItersS_L = 0;
                    belosLHS.PutScalar(0.);
                    b_norm = r_0_norm = belosRES();
                    if (usePrevGuess) {
                        belosLHS = static_cast<SV>(u.Data.append(p.Data));
                        r_0_norm = belosRES();
                    }
                    belosSolver->setParameters(rcpFromRef(belosParams->set("Convergence Tolerance", tol * b_norm / r_0_norm)));
                    belosProblem.setProblem();
                    belosSolver->setProblem(rcpFromRef(belosProblem));
                    belosSolverResult = belosSolver->solve();
                solveTime = logger.end();
                if (belosSolverResult == Converged) logger.log("belos converged");
                else {
                    logger.wrn("belos did not converge");
                    if (useInnerIters && i > 1) {
                        solveWastedTime = solveTime;
                        numItersA = numItersS_M = numItersS_L = 0;
                        logger.beg("build new factorization");
                            runFactorization();
                        factorizationTime = logger.end();
                        logger.beg("linear solve w/ new factorization");
                            if (usePrevGuess) belosLHS = static_cast<SV>(u.Data.append(p.Data));
                            else belosLHS.PutScalar(0.);
                            belosProblem.setProblem();
                            belosSolver->setProblem(rcpFromRef(belosProblem));
                            belosSolverResult = belosSolver->solve();
                            if (belosSolverResult == Converged) logger.log("belos converged");
                            else logger.wrn("belos did not converge");
                        solveTime = logger.end();
                    }
                }
                logger.beg("convert from Epetra");
                    for (size_t i = 0; i < n; ++i) u.Data[i] = belosLHS[i];
                    for (size_t i = 0; i < m; ++i) p.Data[i] = belosLHS[n + i];
                    p.Data -= dot(M_p.Data * p.Data, I_p) / dot(M_p.Data * I_p, I_p) * I_p;
                logger.end();
                logger.beg("project surface vorticity");
                    projectVorticity();
                projectTime = logger.end();
                logger.beg("output");
                    exportStats(i);
                logger.end();
            logger.end();
        }
        #ifdef HAVE_MPI
            MPI_Finalize() ;
        #endif
        return 0;
    } catch (std::exception const & e) {
        logger.err(e.what());
    } catch (DROPSErrCL const & e) {
        e.what(logger.buf);
        logger.err(logger.buf.str());
    } catch (...) {
        logger.err("unknown error");
    }
    return 1;
}