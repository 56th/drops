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
#include "levelset/adaptriang.h"
#include "misc/dynamicload.h"
#include "surfactant/ifacetransp.h"
#include "out/VTKWriter.hpp"
#include "SurfNavierStokesData.hpp"
#include "SingletonLogger.hpp"
// belos (iterative solvers)
#include "BelosSolverFactory.hpp"
#include "trilinos/Belos_LinearSolver.hpp"
// amesos2 (sparse-direct solvers)
#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"
// trilinos (vectors and operators / matrices)
#include "trilinos/Epetra_OperatorApply.hpp"
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
            auto BDF = inpJSON.get<size_t>("Time.BDF");
            if (BDF != 1 && BDF != 2) throw std::invalid_argument("use BDF = 1 or 2");
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
            InstatScalarFunction phi = [&](Point3DCL const & x, double t) { return surface->phi(x, t); };
            if (inpJSON.get<bool>("Surface.UseExactDistanceFunc")) phi = [&](Point3DCL const & x, double t) { return surface->dist(x, t); };
            auto velName = inpJSON.get<std::string>("SurfNavierStokes.IC.Velocity.Name");
            auto preName = inpJSON.get<std::string>("SurfNavierStokes.IC.Pressure.Name");
            auto surfNavierStokesData = surfNavierStokesDataFactory(*surface, velName, preName, inpJSON);
            auto u_N_min = inpJSON.get<double>("Surface.MinSurfaceSpeed");
            auto NarrowBandWidthScaling = inpJSON.get<double>("Surface.NarrowBandWidthScaling");
            FESystem surfNSystem;
            auto& gamma = surfNSystem.params.surfNavierStokesParams.gamma;
            gamma = inpJSON.get<double>("SurfNavierStokes.gamma");
            auto nu = inpJSON.get<double>("SurfNavierStokes.nu");
            if (nu <= 0.) throw std::invalid_argument("viscosity must be non-negative");
            surfNSystem.params.surfNavierStokesParams.nu = { nu, nu };
            auto& t = surfNSystem.params.t;
            auto& levelSet = surfNSystem.params.levelSet;
            surfNSystem.params.numbOfVirtualSubEdges = inpJSON.get<size_t>("SurfNavierStokes.NumbOfVirtualSubEdges");
            surfNSystem.params.surfNavierStokesParams.m_g = surfNavierStokesData.m_g;
            surfNSystem.params.surfNavierStokesParams.f_T = surfNavierStokesData.f_T;
            auto formulation = inpJSON.get<std::string>("SurfNavierStokes.Formulation");
            auto stab = inpJSON.get<std::string>("SurfNavierStokes.PressureStab.Type");
            auto useTangMassMat = inpJSON.get<bool>("SurfNavierStokes.UseTangentialMassMatrix");
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
                auto shift = inpJSON.get<Point3DCL>("Surface.Shift.Dir", Point3DCL(0., 0., 0.)); {
                    if (shift.norm()) shift /= shift.norm();
                    shift *= std::abs(inpJSON.get<double>("Surface.Shift.Norm", 0.));
                }
                mg.Transform([&](Point3DCL const & p) { return p - shift; });
                logger.buf << "surface shift: " << shift;
                logger.log();
            logger.end();
            logger.beg("refine towards the surface");
                AdapTriangCL adap(mg);
                DistMarkingStrategyCL markerLset([&](Point3DCL const & x, double) { return phi(x, 0.); }, inpJSON.get<double>("Mesh.AdaptRef.Width"), inpJSON.get<int>("Mesh.AdaptRef.CoarsestLevel"), inpJSON.get<int>("Mesh.AdaptRef.FinestLevel"));
                adap.set_marking_strategy(&markerLset);
                adap.MakeInitialTriang();
                adap.set_marking_strategy(nullptr);
                logger.buf
                    << "h = " << h << '\n'
                    << "numb of tetras = " << mg.GetNumTriangTetra();
                logger.log();
            logger.end();
        logger.end();
        logger.beg("define FE space for level-set");
            struct { size_t vel, pre, speed, lst; } numActiveTetras;
            IdxDescCL lstIdx(P2_FE);
            numActiveTetras.lst = lstIdx.DistributeDOFs(mg.GetLastLevel(), mg);
            levelSet.SetIdx(&lstIdx);
            VecDescCL distFunc(&lstIdx);
            logger.buf << "numb of active tetras for levelset: " << numActiveTetras.lst << " (" << (100. * numActiveTetras.lst) / mg.GetNumTriangTetra() << "%)";
            logger.log();
        logger.end();
        logger.beg("define (bi)linear forms");
            IdxDescCL velIdx(vecP2IF_FE), preIdx(P1IF_FE), speedIdx(P2IF_FE);
            FEMatDescCL
                M_u(&velIdx, &velIdx, useTangMassMat ? &LocalAssembler::M_t_vecP2vecP2 : &LocalAssembler::M_vecP2vecP2),
                AL_u(&velIdx, &velIdx, &LocalAssembler::AL_vecP2vecP2),
                N_u(&velIdx, &velIdx, &LocalAssembler::rho_N_vecP2vecP2),
                H_u(&velIdx, &velIdx, &LocalAssembler::rho_H_vecP2vecP2),
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
                u_star(&velIdx), u(&velIdx), u_prev(&velIdx), surf_curl_u(&preIdx), w_T(&velIdx), u_N(&speedIdx),
                p_star(&preIdx), p(&preIdx);
        logger.end();
        logger.beg("set up vtk");
            VTKWriter::VTKParams vtkParams; {
                vtkParams.path = dirName + "/vtk/" + velName + inpJSON.get<std::string>("Surface.Name");
                vtkParams.mg = &mg;
                vtkParams.binary = binary;
                vtkParams.staticGrid = surface->isStationary();
                vtkParams.builder = inpJSON.get<bool>("Output.VTK.ExportFullGrid") ? &VTKWriter::buildFullGrid : &VTKWriter::buildInterfaceGrid;
            }
            VTKWriter vtkWriter(vtkParams);
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
            SolverFactory<ST, MV, OP> belosFactory;
            logger.buf << "available iterations: " << belosFactory.supportedSolverNames();
            logger.log();
            Belos_LinearSolver linearSolver; {
                auto belosParams = parameterList();
                belosParams->set("Convergence Tolerance", inpJSON.get<double>("Solver.Outer.RelResTol"));
                belosParams->set("Num Blocks", inpJSON.get<int>("Solver.Outer.KrylovSubspaceSize"));
                belosParams->set("Maximum Iterations", inpJSON.get<int>("Solver.Outer.MaxIter"));
                belosParams->set("Output Frequency", inpJSON.get<int>("Solver.Outer.OutputFrequency"));
                belosParams->set("Verbosity", Errors + Warnings + StatusTestDetails/* + TimingDetails + FinalSummary + IterationDetails*/);
                belosParams->set<int>("Output Style", Brief);
                linearSolver.solver = belosFactory.create(inpJSON.get<std::string>("Solver.Outer.Iteration"), belosParams);
                logger.buf << "outer solver: " << linearSolver.solver->description();
                logger.log();
            }
            linearSolver.zeroIniGuess = !inpJSON.get<bool>("Solver.UsePreviousFrameAsInitialGuess");
            RCP<Epetra_Map> mapVelocity, mapVelocityComp, mapPressure;
            RCP<MT> A, B, C;
            linearSolver.system.mtx = rcp(new Epetra_OperatorApply([&](MV const & X, MV& Y) {
                double* view;
                X(0)->ExtractView(&view);
                SV x1(Epetra_DataAccess::View, *mapVelocity, view);
                SV x2(Epetra_DataAccess::View, *mapPressure, view + velIdx.NumUnknowns());
                Y(0)->ExtractView(&view);
                SV y11(Epetra_DataAccess::View, *mapVelocity, view), y12(*mapVelocity, true);
                SV y21(Epetra_DataAccess::View, *mapPressure, view + velIdx.NumUnknowns()), y22(*mapPressure, true);
                // y1
                A->Multiply(false, x1, y11);
                B->Multiply(true, x2, y12);
                y11.Update(1., y12, 1.);
                // y2
                B->Multiply(false, x1, y21);
                C->Multiply(false, x2, y22);
                y21.Update(1., y22, 1.);
            }));
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
                    linearSolver.updatePreconditioner = [&]() {
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
                    Epetra_OperatorApply::ApplyType invA = [&](MV const &X, MV &Y) {
                        double *viewX, *viewY;
                        X(0)->ExtractView(&viewX);
                        Y(0)->ExtractView(&viewY);
                        std::vector<SV> x, y;
                        x.reserve(numBlocksA);
                        y.reserve(numBlocksA);
                        auto n_i = velIdx.NumUnknowns() / 3;
                        for (size_t i = 0; i < numBlocksA; ++i) {
                            x.emplace_back(Epetra_DataAccess::View, *mapVelocityBlock, viewX + i * n_i);
                            y.emplace_back(Epetra_DataAccess::View, *mapVelocityBlock, viewY + i * n_i);
                        }
                        #pragma omp parallel for
                        for (size_t i = 0; i < numBlocksA; ++i) {
                            amesosSolver[i]->setB(rcpFromRef(x[i]));
                            amesosSolver[i]->setX(rcpFromRef(y[i]));
                            amesosSolver[i]->solve();
                        }
                    };
                    if (precTypeA == "BlockTriangular")
                        invA = [&](MV const &X, MV &Y) {
                            double *viewX, *viewY;
                            X(0)->ExtractView(&viewX);
                            Y(0)->ExtractView(&viewY);
                            std::vector<SV> x, y;
                            x.reserve(numBlocksA);
                            y.reserve(numBlocksA);
                            auto n_i = velIdx.NumUnknowns() / 3;
                            for (size_t i = 0; i < numBlocksA; ++i) {
                                x.emplace_back(Epetra_DataAccess::View, *mapVelocityBlock, viewX + i * n_i);
                                y.emplace_back(Epetra_DataAccess::View, *mapVelocityBlock, viewY + i * n_i);
                            }
                            for (int i = numBlocksA - 1; i >= 0; --i) { // backward substitution
                                SV rhs(Epetra_DataAccess::Copy, x[i], 0);
                                for (size_t j = i + 1; j < numBlocksA; ++j) {
                                    SV A_y(*mapVelocityBlock, true);
                                    A_block[i][j]->Multiply(false, y[j], A_y);
                                    rhs.Update(-1., A_y, 1.);
                                }
                                amesosSolver[i]->setB(rcpFromRef(rhs));
                                amesosSolver[i]->setX(rcpFromRef(y[i]));
                                amesosSolver[i]->solve();
                            }
                        };
                logger.end();
                logger.beg("schur complement block");
                    Belos_LinearSolver linearSolverS_M, linearSolverS_L; {
                        auto belosParams = parameterList();
                        belosParams->set("Maximum Iterations", inpJSON.get<int>("Solver.Inner.S.MaxIter"));
                        belosParams->set("Convergence Tolerance", inpJSON.get<double>("Solver.Inner.S.RelResTol"));
                        linearSolverS_M.solver = belosFactory.create(inpJSON.get<std::string>("Solver.Inner.S.Iteration"), belosParams);
                        linearSolverS_L.solver = belosFactory.create(inpJSON.get<std::string>("Solver.Inner.S.Iteration"), belosParams);
                        logger.buf << "inner solver: " << linearSolverS_M.solver->description();
                        logger.log();
                    }
                    linearSolverS_M.mute = linearSolverS_L.mute = true;
                    auto invS = [&](MV const & X, MV& Y) {
                        Y.PutScalar(0.); // ini guess
                        auto& Y_M = Y;
                        auto  Y_L = Y;
                        auto X_nrm = X; { // normalized rhs
                            double mean;
                            X_nrm.MeanValue(&mean);
                            for (size_t i = 0; i < preIdx.NumUnknowns(); ++i) (*X_nrm(0))[i] -= mean;
                        }
                        // Y_M
                        linearSolverS_M.system.lhs = rcpFromRef(*Y_M(0));
                        linearSolverS_M.system.rhs = rcpFromRef(*X_nrm(0));
                        linearSolverS_M.solve();
                        // Y_L
                        if (inpJSON.get<bool>("Solver.Inner.S.S_L")) {
                            linearSolverS_L.system.lhs = rcpFromRef(*Y_L(0));
                            linearSolverS_L.system.rhs = rcpFromRef(*X_nrm(0));
                            linearSolverS_L.solve();
                        }
                        // result
                        Y_M.Update(1., Y_L, 1.);
                    };
                logger.end();
                auto precType = inpJSON.get<std::string>("Solver.Inner.Type");
                logger.log("preconditioner type: " + precType);
                linearSolver.system.pre = precType == "BlockDiagonal"
                    ? rcp(new Epetra_OperatorApply([&](MV const & X, MV& Y) {
                        double* view;
                        X(0)->ExtractView(&view);
                        SV x1(Epetra_DataAccess::View, *mapVelocity, view);
                        SV x2(Epetra_DataAccess::View, *mapPressure, view + velIdx.NumUnknowns());
                        Y(0)->ExtractView(&view);
                        SV y1(Epetra_DataAccess::View, *mapVelocity, view);
                        SV y2(Epetra_DataAccess::View, *mapPressure, view + velIdx.NumUnknowns());
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
                    }))
                    : rcp(new Epetra_OperatorApply([&](MV const & X, MV& Y) {
                        double* view;
                        X(0)->ExtractView(&view);
                        SV x1(Epetra_DataAccess::View, *mapVelocity, view);
                        SV x2(Epetra_DataAccess::View, *mapPressure, view + velIdx.NumUnknowns());
                        Y(0)->ExtractView(&view);
                        SV y1(Epetra_DataAccess::View, *mapVelocity, view);
                        SV y2(Epetra_DataAccess::View, *mapPressure, view + velIdx.NumUnknowns());
                        // y2
                        invS(x2, y2);
                        // y1
                        SV rhs(*mapVelocity, true);
                        B->Multiply(true, y2, rhs);
                        rhs.Update(1., x1, -1.);
                        invA(rhs, y1);
                    }));
            logger.end();
        logger.end();
        logger.beg("t = t_0 = 0");
            t = 0.;
            logger.beg("interpolate level-set");
                levelSet.Interpolate(mg, [&](Point3DCL const & x) { return phi(x, t); });
                distFunc.Interpolate(mg, [&](Point3DCL const & x) { return surface->dist(x, t); });
            logger.end();
            logger.beg("set up FE spaces");
                logger.beg("set up pressure space");
                    numActiveTetras.pre = preIdx.DistributeDOFs(mg.GetLastLevel(), mg, &distFunc);
                    logger.buf << "numb of active tetras for pressure: " << numActiveTetras.pre << " (" << (100. * numActiveTetras.pre) / mg.GetNumTriangTetra() << "%)";
                    logger.log();
                    mapPressure = rcp(new Epetra_Map(static_cast<int>(preIdx.NumUnknowns()), 0, comm)); // for vorticity linear solve
                logger.end();
                logger.beg("set up velocity space");
                    logger.beg("initialize surface speed u_N");
                        speedIdx.DistributeDOFs(mg.GetLastLevel(), mg, &distFunc);
                        u_N.Interpolate(mg, [&](Point3DCL const & x) { return surface->u_N(x, t); });
                        auto u_N_max = supnorm(u_N.Data);
                        if (u_N_max) u_N_max = std::max(u_N_min, u_N_max);
                        auto narrowBandWidth = NarrowBandWidthScaling * BDF * u_N_max * stepSize;
                        logger.buf
                            << "max |u_N| = " << u_N_max << '\n'
                            << "narrow band width = " << narrowBandWidth;
                        logger.log();
                    logger.end();
                    numActiveTetras.vel = velIdx.DistributeDOFs(mg.GetLastLevel(), mg, &distFunc, narrowBandWidth);
                    logger.buf << "numb of active tetras for velocity: " << numActiveTetras.vel << " (" << (100. * numActiveTetras.vel) / mg.GetNumTriangTetra() << "%)";
                    logger.log();
                logger.end();
            auto feTime = logger.end();
            logger.beg("assemble");
                surfNSystem.matrices = { &M_u, &A_u, &AL_u, &S_u, &C_u, &M_p, &C_p, &A_p, &B_pu, &Q_pu };
                setupFESystem(mg, surfNSystem);
                MatrixCL A_sum;
                VectorCL I_p(1., preIdx.NumUnknowns());
            auto assembleTime = logger.end();
            logger.beg("interpolate initial data");
                u.Interpolate(mg, [&](Point3DCL const & x) { return surfNavierStokesData.u_T(x, t); });
                p.Interpolate(mg, [&](Point3DCL const & x) { return surfNavierStokesData.p(x, t); });
                p.Data -= dot(M_p.Data * p.Data, I_p) / dot(M_p.Data * I_p, I_p) * I_p;
                u_prev = u;
            logger.end();
            logger.beg("project surface vorticity");
                Belos_LinearSolver linearSolverVorticity; {
                    auto belosParams = parameterList();
                    belosParams->set("Maximum Iterations", 1000);
                    belosParams->set("Convergence Tolerance", 1e-8);
                    belosParams->set("Output Frequency", 10);
                    belosParams->set("Verbosity", Errors + Warnings + StatusTestDetails);
                    belosParams->set<int>("Output Style", Brief);
                    linearSolverVorticity.solver = belosFactory.create("CG", belosParams);
                }
                auto projectVorticity = [&]() {
                    linearSolverVorticity.system.mtx = static_cast<RCP<MT>>(MatrixCL(1., M_p.Data, rho_p, C_p.Data));
                    linearSolverVorticity.system.rhs = static_cast<RCP<SV>>(Q_pu.Data * u.Data);
                    linearSolverVorticity.system.lhs = rcp(new SV(*mapPressure, true));
                    linearSolverVorticity.solve();
                    for (size_t i = 0; i < preIdx.NumUnknowns(); ++i) surf_curl_u.Data[i] = (*linearSolverVorticity.system.lhs)[i];
                };
                projectVorticity();
            logger.end();
            logger.beg("output");
                double Pe, alpha;
                auto exportStats = [&](size_t i) {
                    std::ofstream stats(dirName + "/stats/t_" + std::to_string(i) + ".json");
                    ParamCL tJSON;
                    tJSON.put("t", t);
                    tJSON.put("h", h);
                    tJSON.put("nu", nu);
                    tJSON.put("gamma", gamma);
                    tJSON.put("DOF.Velocity", velIdx.NumUnknowns());
                    tJSON.put("DOF.Pressure", preIdx.NumUnknowns());
                    tJSON.put("NumTetras.Bulk", mg.GetNumTriangTetra());
                    tJSON.put("NumTetras.Pressure", numActiveTetras.pre);
                    tJSON.put("NumTetras.Velocity", numActiveTetras.vel);
                    tJSON.put("NarrowBandWidth", narrowBandWidth);
                    tJSON.put("MaxSurfaceSpeed", u_N_max);
                    tJSON.put("MeshDepParams.rho_p", rho_p);
                    tJSON.put("MeshDepParams.rho_u", rho_u);
                    tJSON.put("MeshDepParams.tau_u", tau_u);
                    tJSON.put("ElapsedTime.RemapAndDistributeDOF", feTime);
                    tJSON.put("ElapsedTime.Assemble", assembleTime);
                    tJSON.put("ElapsedTime.ProjectVorticity", linearSolverVorticity.stats.time.solve);
                    tJSON.put("Solver.ProjectVorticity.TotalIters", linearSolverVorticity.solver->getNumIters());
                    tJSON.put("Solver.ProjectVorticity.Converged", linearSolverVorticity.stats.result == Converged);
                    tJSON.put("Solver.ProjectVorticity.ResidualNormBelosRelative", linearSolverVorticity.solver->achievedTol());
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
                        tJSON.put("ElapsedTime.Factorization", linearSolver.stats.time.updatePreconditioner);
                        tJSON.put("ElapsedTime.LinearSolve", linearSolver.stats.time.solve);
                        tJSON.put("Solver.Outer.ResidualNorm.r_i", linearSolver.stats.norm.r_i);
                        tJSON.put("Solver.Outer.ResidualNorm.b", linearSolver.stats.norm.b);
                        tJSON.put("Solver.Outer.ResidualNorm.r_0", linearSolver.stats.norm.r_0);
                        tJSON.put("Solver.Outer.ResidualNorm.r_i/b", linearSolver.stats.norm.r_i / linearSolver.stats.norm.b);
                        tJSON.put("Solver.Outer.ResidualNorm.r_i/r_0", linearSolver.stats.norm.r_i / linearSolver.stats.norm.r_0);
                        tJSON.put("Solver.Outer.ResidualNorm.BelosRelative", linearSolver.solver->achievedTol());
                        tJSON.put("Solver.Outer.TotalIters", linearSolver.solver->getNumIters());
                        tJSON.put("Solver.Outer.Converged", linearSolver.stats.result == Converged);
                        tJSON.put("Solver.Inner.S.MeanIters.S_M", linearSolverS_M.solver->getNumIters());
                        tJSON.put("Solver.Inner.S.MeanIters.S_L", linearSolverS_L.solver->getNumIters());
                        tJSON.put("Peclet", Pe);
                        tJSON.put("MassMatrixCoef", alpha);
                    }
                    auto exportFiles = everyStep > 0 && (i % everyStep == 0 || i == numSteps);
                    if (exportVTK && exportFiles) {
                        logger.beg("write vtk");
                            vtkWriter.write(t);
                        auto vtkTime = logger.end();
                        tJSON.put("ElapsedTime.VTK", vtkTime);
                    }
                    if (exportVectors && exportFiles || i > numSteps - BDF) {
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
                logger.beg("interpolate level-set");
                    levelSet.Interpolate(mg, [&](Point3DCL const & x) { return phi(x, t); });
                    distFunc.Interpolate(mg, [&](Point3DCL const & x) { return surface->dist(x, t); });
                logger.end();
                logger.beg("set up FE spaces");
                    logger.beg("set up pressure space");
                        {
                            auto n = preIdx.DistributeDOFs(mg.GetLastLevel(), mg, &distFunc);
                            logger.buf << "numb of active tetras for pressure: " << numActiveTetras.pre << " -> " << n << " (" << (100. * n) / mg.GetNumTriangTetra() << "%)";
                            logger.log();
                            std::swap(n, numActiveTetras.pre);
                        }
                        mapPressure = rcp(new Epetra_Map(static_cast<int>(preIdx.NumUnknowns()), 0, comm));
                    logger.end();
                    logger.beg("set up velocity space");
                        logger.beg("initialize surface speed u_N");
                            speedIdx.DistributeDOFs(mg.GetLastLevel(), mg, &distFunc);
                            u_N.Interpolate(mg, [&](Point3DCL const & x) { return surface->u_N(x, t); });
                            u_N_max = supnorm(u_N.Data);
                            if (u_N_max) u_N_max = std::max(u_N_min, u_N_max);
                            narrowBandWidth = NarrowBandWidthScaling * BDF * u_N_max * stepSize;
                            logger.buf
                                << "max |u_N| = " << u_N_max << '\n'
                                << "narrow band width = " << narrowBandWidth;
                            logger.log();
                        logger.end();
                        {
                            auto n = velIdx.DistributeDOFs(mg.GetLastLevel(), mg, &distFunc, narrowBandWidth);
                            logger.buf << "numb of active tetras for velocity: " << numActiveTetras.vel << " -> " << n << " (" << (100. * n) / mg.GetNumTriangTetra() << "%)";
                            logger.log();
                            std::swap(n, numActiveTetras.vel);
                        }
                        mapVelocity = rcp(new Epetra_Map(static_cast<int>(velIdx.NumUnknowns()), 0, comm));
                        mapVelocityComp = rcp(new Epetra_Map(static_cast<int>(velIdx.NumUnknowns() / 3), 0, comm));
                    logger.end();
                feTime = logger.end();
                logger.beg("assemble");
                    alpha = i == 1 || BDF == 1 ? 1. / stepSize : 1.5 / stepSize;
                    if (surfNavierStokesData.w_T) w_T.Interpolate(mg, [&](Point3DCL const & x) { return surfNavierStokesData.w_T(x, t); }); // Oseen (and Stokes) case
                    else { // Navier-Stokes case
                        if (BDF == 1) w_T.Data = u.Data;
                        else w_T.Data = 2. * u.Data - u_prev.Data;
                    }
                    Pe = supnorm(w_T.Data) / nu;
                    logger.buf
                        << "alpha = " << alpha << '\n'
                        << "Pe = " << Pe;
                    logger.log();
                    surfNSystem.params.surfNavierStokesParams.w_T = Pe ? &w_T : nullptr;
                    surfNSystem.params.surfNavierStokesParams.u_N = u_N_max ? &u_N : nullptr;
                    if (u_N_max) {
                        logger.log("assembling all matrices due to nnz surf speed");
                        surfNSystem.matrices = { &N_u, &H_u, &M_u, &A_u, &AL_u, &S_u, &C_u, &M_p, &C_p, &A_p, &B_pu, &Q_pu };
                    }
                    else if (Pe) {
                        logger.log("assembling convection mtx");
                        surfNSystem.matrices = { &N_u };
                    }
                    surfNSystem.vectors = { &F_u, &G_p };
                    setupFESystem(mg, surfNSystem);
                    // system mtx
                    A_sum.LinComb(alpha, M_u.Data, gamma, AL_u.Data, nu, A_u.Data, tau_u, S_u.Data, rho_u, C_u.Data);
                    if (Pe || u_N_max) A_sum.LinComb(1., MatrixCL(A_sum), 1., N_u.Data);
                    if (u_N_max) A_sum.LinComb(1., MatrixCL(A_sum), 1., H_u.Data);
                    // system rhs
                    if (i == 1 || BDF == 1) F_u.Data += (1. / stepSize) * (M_u.Data * u.Data);
                    else F_u.Data += (2. / stepSize) * (M_u.Data * u.Data) - (.5 / stepSize) * (M_u.Data * u_prev.Data);
                    I_p.resize(preIdx.NumUnknowns(), 1.);
                    G_p.Data -= (dot(G_p.Data, I_p) / dot(I_p, I_p)) * I_p;
                    // for the next step
                    u_prev = u;
                assembleTime = logger.end();
                logger.beg("convert to Epetra");
                    A = static_cast<RCP<MT>>(A_sum);
                    logCRS(A, "A");
                    {
                        auto A_sum_blocks = A_sum.Split(numBlocksA, numBlocksA);
                        for (size_t i = 0; i < numBlocksA; ++i)
                            for (size_t j = 0; j < numBlocksA; ++j)
                                A_block[i][j] = static_cast<RCP<MT>>(A_sum_blocks[i][j]);
                        if (numBlocksA != 1)
                            logCRS(A_block[0][0], "A_{ij}");
                    }
                    B = static_cast<RCP<MT>>(B_pu.Data);
                    logCRS(B, "B");
                    C = static_cast<RCP<MT>>(MatrixCL(-rho_p, C_p.Data));
                    logCRS(C, "C := -rho_p (pressure volume stabilization mtx)");
                    linearSolverS_M.system.mtx = static_cast<RCP<MT>>(MatrixCL(1. / (gamma + nu), M_p.Data, rho_p, C_p.Data));
                    logCRS(linearSolverS_M.system.mtx, "S_M := (nu + gamma)^{-1} M_p - C");
                    linearSolverS_L.system.mtx = static_cast<RCP<MT>>(MatrixCL(1. / alpha, A_p.Data, rho_p, C_p.Data));
                    logCRS(linearSolverS_L.system.mtx, "S_L := alpha^{-1} L_p - C");
                    linearSolver.system.rhs = static_cast<RCP<SV>>(F_u.Data.append(G_p.Data));
                    linearSolver.system.lhs = static_cast<RCP<SV>>(u.Data.append(p.Data));
                logger.end();
                linearSolver.solve();
                logger.beg("convert from Epetra");
                    for (size_t i = 0; i < velIdx.NumUnknowns(); ++i) u.Data[i] = (*linearSolver.system.lhs)[i];
                    for (size_t i = 0; i < preIdx.NumUnknowns(); ++i) p.Data[i] = (*linearSolver.system.lhs)[velIdx.NumUnknowns() + i];
                    p.Data -= dot(M_p.Data * p.Data, I_p) / dot(M_p.Data * I_p, I_p) * I_p;
                logger.end();
                logger.beg("project surface vorticity");
                    projectVorticity();
                logger.end();
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