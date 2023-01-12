/// \file SurfCahnHilliard.cpp
/// \brief Trace FEM discretization of a surface Cahn-Hilliard problem
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
#include "hermite_cubic/hermite_cubic.hpp"
#include "SurfCahnHilliardData.hpp"
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
            read_parameter_file_from_cmdline(inpJSON, argc, argv, "../../param/savch/No_Bnd_Condition.json");
            dynamicLoad(inpJSON.get<std::string>("General.DynamicLibsPrefix"), inpJSON.get<std::vector<std::string>>("General.DynamicLibs"));
            auto dirName = inpJSON.get<std::string>("Output.Directory");
            system((
                "mkdir -p " + dirName + " ; rm -r -f " + dirName + "/* ; "
                "mkdir " + dirName + "/stats ; mkdir " + dirName + "/vectors ; mkdir " + dirName + "/vtk"
            ).c_str());
            {
                std::ofstream input(dirName + "/input.json");
                input << inpJSON;
                logger.buf << inpJSON;
                logger.log();
            }
            auto h = inpJSON.get<Point3DCL>("Mesh.E1")[0] / inpJSON.get<double>("Mesh.N1") * std::pow(2., -inpJSON.get<double>("Mesh.AdaptRef.FinestLevel"));
            auto eps = inpJSON.get<double>("SurfCahnHilliard.Epsilon");
            auto rho_vol = inpJSON.get<double>("SurfCahnHilliard.VolumeStab.Scaling") * pow(h, inpJSON.get<double>("SurfCahnHilliard.VolumeStab.Power")); // constant for volume stabilization
            auto finalTime = inpJSON.get<double>("Time.FinalTime");
            auto everyStep = inpJSON.get<int>("Output.EveryStep");
            auto binary = inpJSON.get<bool>("Output.Binary");
            auto exportVectors = inpJSON.get<bool>("Output.Vectors");
            auto surface = surfaceFactory(inpJSON);
            auto testName = inpJSON.get<std::string>("SurfCahnHilliard.IC.Name");
            auto surfCahnHilliardData = surfCahnHilliardDataFactory(*surface, testName, inpJSON);
            auto mobilityScaling = inpJSON.get<double>("SurfCahnHilliard.MobilityScaling");
            auto beta_s = inpJSON.get<double>("SurfCahnHilliard.Beta_s");
            FESystem surfCHSystem;
            auto& useDegenerateMobility = surfCHSystem.params.surfCahnHilliardParams.useDegenerateMobility;
            useDegenerateMobility = inpJSON.get<bool>("SurfCahnHilliard.UseDegenerateMobility");
            auto& t = surfCHSystem.params.t;
            auto& levelSet = surfCHSystem.params.levelSet;
            surfCHSystem.params.numbOfVirtualSubEdges = inpJSON.get<size_t>("SurfCahnHilliard.NumbOfVirtualSubEdges");
            surfCHSystem.params.surfCahnHilliardParams.f = surfCahnHilliardData.f;
            auto meshCoarseLevel = std::max(0, inpJSON.get<int>("Mesh.AdaptRef.CoarsestLevel"));
            auto meshFineLevel = std::max(0, inpJSON.get<int>("Mesh.AdaptRef.FinestLevel"));
            if (meshFineLevel < meshCoarseLevel) meshFineLevel = meshCoarseLevel;
            auto prolongationLevel = inpJSON.get<int>("SurfCahnHilliard.ProlongateFromLevelNo");
            bool import_chi_from_file = inpJSON.get<bool>("SurfCahnHilliard.IC.ImportChi");
            if (prolongationLevel < meshCoarseLevel) prolongationLevel = meshCoarseLevel;
            if (prolongationLevel > meshFineLevel) prolongationLevel = meshFineLevel;
            auto F_rho = inpJSON.get<double>("Time.Adaptive.rho");
            auto F_tol = inpJSON.get<double>("Time.Adaptive.Tol");
            auto F_min = inpJSON.get<double>("Time.Adaptive.MinStepSize", 0.);
            auto F_max = inpJSON.get<double>("Time.Adaptive.MaxStepSize", std::numeric_limits<double>::max());
            auto F = [=](double e, double dt) {
                if (!e) return inpJSON.get<double>("Time.StepSize");
                auto f = F_rho * std::sqrt(F_tol / e) * dt;
                if (f > F_max) return F_max;
                if (f < F_min) return F_min;
                return f;
            };
            auto r_sav = 0.0;
            auto E_1 = 0.0;
            logger.buf
                << surface->description() << '\n'
                << surfCahnHilliardData.description << '\n'
                << "exact soln: " << (surfCahnHilliardData.exact ? "yes" : "no") << '\n'
                << "prolongate CH IC from lvl " << prolongationLevel;
            logger.log();
        logger.end();
        logger.beg("set up chemical potential");
            auto xi = inpJSON.get<double>("SurfCahnHilliard.ChemicalPotentialScaling");
            auto c0 = inpJSON.get<double>("SurfCahnHilliard.c_0");
            auto c0_l = std::min(c0, 1. - c0) / sqrt(3.);
            auto f_0 = [&](VectorCL const & c) {
                auto res = c;
                for (auto& el : res) el = xi * std::pow(el * (1. - el), 2.)+1.001;
                return res;
            };
        auto chemicalPotential = [&](double c) { // f'_0
            if (c < 0.) return xi * (1. - c0) * c;
            if (c > 1.) return xi * c0 * (c - 1.);
            double x[1], f[1], d[1], s[1], t[1];
            x[0] = c;
            if      (c < c0 - c0_l) hermite_cubic_value(0., 0., 1. - c0, c0 - c0_l, 1. / (12. * sqrt(3.)), 0., 1, x, f, d, s, t);
            else if (c < c0)        hermite_cubic_value(c0 - c0_l, 1. / (12. * sqrt(3.)), 0., c0, 0., -std::max(c0 * c0, (1. - c0) * (1. - c0)), 1, x, f, d, s, t);
            else if (c < c0 + c0_l) hermite_cubic_value(c0, 0., -std::max(c0 * c0, (1. - c0) * (1. - c0)), c0 + c0_l, -1. / (12. * sqrt(3.)), 0., 1, x, f, d, s, t);
            else                    hermite_cubic_value(c0 + c0_l, -1. / (12. * sqrt(3.)), 0., 1., 0., c0, 1, x, f, d, s, t);
            return xi * f[0];
        };
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
                DistMarkingStrategyCL markerLset([&](Point3DCL const & x, double) { return surface->dist(x); }, inpJSON.get<double>("Mesh.AdaptRef.Width"), meshCoarseLevel, meshFineLevel);
                adap.set_marking_strategy(&markerLset);
                adap.MakeInitialTriang();
                adap.set_marking_strategy(nullptr);
                auto numbOfTetras = mg.GetNumTriangTetra();
                logger.buf
                    << "h = " << h << '\n'
                    << "numb of tetras = " << numbOfTetras;
                logger.log();
            logger.end();
        logger.end();
        logger.beg("interpolate level-set");
            IdxDescCL lstIdx(P2_FE); {
                auto n = lstIdx.DistributeDOFs(mg.GetLastLevel(), mg);
                logger.buf << "numb of active tetras for levelset: " << n << " (" << (100. * n) / mg.GetNumTriangTetra() << "%)";
                logger.log();
            }
            levelSet.SetIdx(&lstIdx);
            levelSet.Interpolate(mg, [&](Point3DCL const & x) { return surface->dist(x); });
        logger.end();
        logger.beg("set up FE spaces");
            IdxDescCL chiIdx(P1IF_FE); {
                auto n = chiIdx.DistributeDOFs(mg.GetLastLevel(), mg, &levelSet);
                logger.buf << "numb of active tetras for concentration: " << n << " (" << (100. * n) / mg.GetNumTriangTetra() << "%)";
                logger.log();
            }
            IdxDescCL velIdx(vecP2IF_FE); {
                auto n = velIdx.DistributeDOFs(mg.GetLastLevel(), mg, &levelSet);
                logger.buf << "numb of active tetras for velocity: " << n << " (" << (100. * n) / mg.GetNumTriangTetra() << "%)";
                logger.log();
            }
            size_t m = chiIdx.NumUnknowns();
            logger.beg("set up FE system");
                FEMatDescCL
                    M(&chiIdx, &chiIdx, &LocalAssembler::M_P1P1),
                    C(&chiIdx, &chiIdx, &LocalAssembler::C_n_P1P1),
                    N(&chiIdx, &chiIdx, &LocalAssembler::N_P1P1),
                    A_one(&chiIdx, &chiIdx, &LocalAssembler::A_P1P1),
                    A_deg(&chiIdx, &chiIdx, &LocalAssembler::LaplaceM_P1P1);
                FEVecDescCL f(&chiIdx, &LocalAssembler::F_concentration_P1);
                VecDescCL
                    chi_star(&chiIdx), chi(&chiIdx), chi_BDF1(&chiIdx), chi_BDF2(&chiIdx), chi_extrap(&chiIdx), omega(&chiIdx), wind(&velIdx),
                    F_chi(&chiIdx), F_omega(&chiIdx), chePot(&chiIdx);
            logger.end();
        logger.end();
        logger.beg("set up vtk");
            VTKWriter::VTKParams vtkParams; {
                vtkParams.path = dirName + "/vtk/" + testName + inpJSON.get<std::string>("Surface.Name");
                vtkParams.mg = &mg;
                vtkParams.binary = binary;
                vtkParams.staticGrid = surface->isStationary();
                vtkParams.builder = inpJSON.get<bool>("Output.VTK.ExportFullGrid") ? &VTKWriter::buildFullGrid : &VTKWriter::buildInterfaceGrid;
            }
            VTKWriter vtkWriter(vtkParams);
            vtkWriter.add({ "level-set", levelSet });
            if (everyStep > 0) {
                if (inpJSON.get<bool>("Output.VTK.Concentration")) {
                    vtkWriter.add({"c_h", chi});
                    if (surfCahnHilliardData.exact) vtkWriter.add({"c_*", chi_star});
                }
                if (inpJSON.get<bool>("Output.VTK.Wind")) vtkWriter.add({"wind", wind});
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
                belosParams->set("Verbosity", Errors + Warnings + StatusTestDetails);
                belosParams->set<int>("Output Style", Brief);
                linearSolver.solver = belosFactory.create(inpJSON.get<std::string>("Solver.Outer.Iteration"), belosParams);
                logger.buf << "outer solver: " << linearSolver.solver->description();
                logger.log();
            }
            linearSolver.zeroIniGuess = !inpJSON.get<bool>("Solver.UsePreviousFrameAsInitialGuess");
            logger.beg("set up preconditioner");
                Teuchos::RCP<OP> preMTX;
                RCP<Amesos2::Solver<MT, MV>> amesosSolver;
                linearSolver.system.pre = rcp(new Epetra_OperatorApply([&](MV const &X, MV &Y) {
                    amesosSolver->setB(rcpFromRef(X));
                    amesosSolver->setX(rcpFromRef(Y));
                    amesosSolver->solve();
                }));
                linearSolver.updatePreconditioner = [&]() {
                    auto mtx = rcp_dynamic_cast<MT>(preMTX);
                    amesosSolver = Amesos2::create<MT, MV>(inpJSON.get<std::string>("Solver.Inner.Iteration"), mtx);
                    amesosSolver->symbolicFactorization();
                    amesosSolver->numericFactorization();
                    auto amesosStatus = amesosSolver->getStatus();
                    logger.buf << "numb of nonzeros in L + U = " << amesosStatus.getNnzLU() << " (" << (100. * amesosStatus.getNnzLU()) / (static_cast<double>(mtx->NumGlobalRows()) * mtx->NumGlobalCols()) << "%)";
                    logger.log();
                };
            logger.end();
        logger.end();
        logger.beg("t = t_0 = 0");
            size_t i = 0;
            t = 0.;
            logger.beg("assemble");
                surfCHSystem.matrices = { &M, &C, &A_one };
                setupFESystem(mg, surfCHSystem);
                VectorCL I_p(1., m);
            auto assembleTime = logger.end();
            logger.beg("interpolate initial data");
                if(import_chi_from_file){
                    logger.beg("Reading initial values of chi from file");
                    auto path = inpJSON.get<std::string>("SurfCahnHilliard.IC.ChiPath");
                    std::ifstream file(path);
                    chi.Read(file, binary);
                    logger.end();
                }
                else {
                    IdxDescCL chiIdxExtCoarse(P1_FE), chiIdxExt(P1_FE);
                    chiIdxExtCoarse.CreateNumbering(prolongationLevel, mg);
                    chiIdxExt.CreateNumbering(meshFineLevel, mg);
                    VecDescCL chiExtCoarse(&chiIdxExtCoarse), chiExt(&chiIdxExt);
                    std::vector<ProlongationCL<double>> P; {
                        for (size_t l = prolongationLevel + 1; l <= meshFineLevel; ++l) {
                            P.emplace_back(ProlongationCL<double>(mg));
                            IdxDescCL chiIdxExtFrom(P1_FE), chiIdxExtTo(P1_FE);
                            chiIdxExtFrom.CreateNumbering(l - 1, mg);
                            chiIdxExtTo.CreateNumbering(l, mg);
                            P.back().Create(&chiIdxExtFrom, &chiIdxExtTo);
                        }
                    }
                    double raftFraction, raftFractionError;
                    do {
                        chiExtCoarse.Interpolate(mg , [&](Point3DCL const & x) { return surfCahnHilliardData.chi(x, t); });
                        chiExt.Data = chiExtCoarse.Data;
                        for (auto const & p : P)
                            chiExt.Data = p * chiExt.Data;
                        Restrict(mg, chiExt, chi);
                        if (surfCahnHilliardData.raftRatio > 0.) {
                            raftFraction = dot(M.Data * chi.Data, I_p) / dot(M.Data * I_p, I_p);
                            raftFractionError = std::fabs(surfCahnHilliardData.raftRatio - raftFraction) / surfCahnHilliardData.raftRatio;
                            logger.buf << "raft ratio error (%) = " << raftFractionError;
                            logger.log();
                        }
                    } while (surfCahnHilliardData.raftRatio > 0. && raftFractionError > .001);
                }
                auto chi_prev = chi;
                wind.Interpolate(mg, [&](Point3DCL const & x) { return surfCahnHilliardData.wind(x, t); });
                omega.Data = 0.;
                linearSolver.system.lhs = static_cast<RCP<SV>>(chi.Data.append(omega.Data));
            logger.end();
            logger.beg("output");
                size_t numTries;
                double e = 0., dt, alpha;
                dt = inpJSON.get<double>("Time.StepSize");
                auto exportStats = [&](size_t i) {
                    std::ofstream stats(dirName + "/stats/t_" + std::to_string(i) + ".json");
                    ParamCL tJSON;
                    tJSON.put("t", t);
                    tJSON.put("h", h);
                    tJSON.put("c_h.Max", chi.Data.max());
                    tJSON.put("c_h.Min", chi.Data.min());
                    tJSON.put("DOF.c_h", m);
                    tJSON.put("DOF.mu_h", m);
                    tJSON.put("MeshDepParams.rho_vol", rho_vol);
                    tJSON.put("ElapsedTime.Assemble", assembleTime);
                    auto surfArea = dot(I_p, M.Data * I_p);
                    tJSON.put("Integral.SurfacaAreaP1", surfArea);
                    tJSON.put("Integral.FESolution.RaftFraction", dot(M.Data * chi.Data, I_p) / surfArea);
                    tJSON.put("Integral.FESolution.PerimeterEstimate", eps * dot(A_one.Data * chi.Data, chi.Data));
                    tJSON.put("Integral.FESolution.LyapunovEnergy", dot(M.Data * I_p, f_0(chi.Data)) + .5 * eps * eps * dot(A_one.Data * chi.Data, chi.Data));
                    tJSON.put("Integral.FESolution.ModifiedSavEnergy", r_sav*r_sav + .5 * eps * eps * dot(A_one.Data * chi.Data, chi.Data) + 0.5* eps * eps *rho_vol * dot(C.Data * chi.Data, chi.Data));
                    tJSON.put("Integral.FESolution.r_SAV", r_sav);
                    tJSON.put("Integral.FESolution.r_SAV_E_1", E_1);
                    if (surfCahnHilliardData.exact) {
                        chi_star.Interpolate(mg, [&](Point3DCL const & x) { return surfCahnHilliardData.chi(x, t); });
                        tJSON.put("Integral.ExactSolution.RaftFraction", dot(M.Data * chi_star.Data, I_p) / surfArea);
                        tJSON.put("Integral.ExactSolution.PerimeterEstimate", eps * dot(A_one.Data * chi_star.Data, chi_star.Data));
                        tJSON.put("Integral.ExactSolution.LyapunovEnergy", dot(M.Data * I_p, f_0(chi_star.Data)) + .5 * eps * eps * dot(A_one.Data * chi_star.Data, chi_star.Data));
                        auto chi_diff = chi_star.Data - chi.Data;
                        tJSON.put("Integral.Error.ConcentrationL2", sqrt(dot(chi_diff, M.Data * chi_diff)));
                        tJSON.put("Integral.Error.ConcentrationH1", sqrt(dot(chi_diff, A_one.Data * chi_diff)));
                    }
                    if (i > 0) {
                        tJSON.put("ElapsedTime.LinearSolve", linearSolver.stats.time.solve);
                        tJSON.put("ElapsedTime.Factorization", linearSolver.stats.time.updatePreconditioner);
                        tJSON.put("dt", dt);
                        tJSON.put("Solver.ResidualNorm.r_i", linearSolver.stats.norm.r_i);
                        tJSON.put("Solver.ResidualNorm.b", linearSolver.stats.norm.b);
                        tJSON.put("Solver.ResidualNorm.r_0", linearSolver.stats.norm.r_0);
                        tJSON.put("Solver.ResidualNorm.r_i/b", linearSolver.stats.norm.r_i / linearSolver.stats.norm.b);
                        tJSON.put("Solver.ResidualNorm.r_i/r_0", linearSolver.stats.norm.r_i / linearSolver.stats.norm.r_0);
                        tJSON.put("Solver.ResidualNorm.BelosRelative", linearSolver.solver->achievedTol());
                        tJSON.put("Solver.TotalIters", linearSolver.solver->getNumIters());
                        tJSON.put("Solver.Converged", linearSolver.stats.result == Converged);
                        tJSON.put("MassMatrixCoef", alpha);
                    }
                    if (i > 1) {
                        tJSON.put("ApaptiveTimeStep.NumbOfTries", numTries);
                        tJSON.put("ApaptiveTimeStep.RelativeError", e);
                    }
                    auto exportFiles = everyStep > 0 && (i % everyStep == 0 || t >= finalTime);
                    if (exportVTK && exportFiles) {
                        logger.beg("write vtk");
                            vtkWriter.write(t);
                        auto vtkTime = logger.end();
                        tJSON.put("ElapsedTime.VTK", vtkTime);
                    }
                    if (exportVectors && exportFiles || t >= finalTime) {
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
                            expVec(chi, "chi");
                        auto vecTime = logger.end();
                        tJSON.put("ElapsedTime.VectorExport", vecTime);
                    }
                    stats << tJSON;
                    logger.buf << tJSON;
                    logger.log();;
                };
                exportStats(i);
            logger.end();
        logger.end();
        while (t < finalTime) {
            logger.beg("t = t_" + std::to_string(++i));
                auto t_prev = t;
                auto dt_prev = dt;
                logger.beg("update time");
                    t = t_prev + dt;
                    logger.buf
                        << "t  = " << t << '\n'
                        << "dt = " << dt;
                    logger.log();
                logger.end();
                logger.beg("BDF1 step");
                    logger.beg("assemble");
                        surfCHSystem.matrices = {};
                        surfCHSystem.params.surfCahnHilliardParams.chi = &chi;
                        if (useDegenerateMobility) surfCHSystem.matrices.push_back(&A_deg);
                        surfCHSystem.params.surfCahnHilliardParams.u_T = &wind;
                        //wind.Interpolate(mg, [&](Point3DCL const & x) { return surfCahnHilliardData.wind(x, t); });
                        //auto wind_max = supnorm(wind.Data);
                        //if (wind_max) surfCHSystem.matrices.push_back(&N);
                        surfCHSystem.vectors = { &f };
                        setupFESystem(mg, surfCHSystem);
                        if(i == 1){chi_extrap.Data = chi.Data;}
                        else {chi_extrap.Data = 2 * chi.Data -  chi_prev.Data;}
                        F_chi.Data = f.Data + (1. / dt) * (M.Data * chi.Data);
                        for (size_t i = 0; i < m; ++i) F_omega.Data[i] = chemicalPotential(chi.Data[i]);
                        auto omega_rhs_scaling = dot(M.Data * F_omega.Data,chi.Data);
                        E_1 = dot(M.Data * f_0(chi.Data),I_p);
                        if(i==1) r_sav = std::sqrt(E_1);
                        omega_rhs_scaling = omega_rhs_scaling / (2.0 * E_1) - r_sav / std::sqrt(E_1);
                        chePot.Data = M.Data * F_omega.Data;
                        F_omega.Data =  omega_rhs_scaling * chePot.Data;
                        MatrixCL M_sav;
                        M_sav.OuterProduct(chePot.Data,chePot.Data);
                    assembleTime = logger.end();
                    logger.beg("convert to Epetra");
                        alpha = 1. / dt;
                        logger.buf << "alpha = " << alpha;
                        logger.log();
                        linearSolver.system.mtx = static_cast<RCP<MT>>(MatrixCL(
                            MatrixCL(alpha, M.Data), MatrixCL(mobilityScaling, useDegenerateMobility ? A_deg.Data : A_one.Data, rho_vol, C.Data),
                            MatrixCL(eps * eps, A_one.Data, beta_s, M.Data, 1/(2*E_1), M_sav, eps * eps * rho_vol, C.Data), MatrixCL(-1., M.Data)
                        ));
                        preMTX = static_cast<RCP<MT>>(MatrixCL(
                                MatrixCL(alpha, M.Data), MatrixCL(mobilityScaling, useDegenerateMobility ? A_deg.Data : A_one.Data, rho_vol, C.Data),
                                MatrixCL(eps * eps, A_one.Data,  1, M.Data, eps * eps * rho_vol, C.Data) , MatrixCL(-1., M.Data)
                        ));
                        logCRS(linearSolver.system.mtx, "{A, B; C, D} block mtx");
                        linearSolver.system.rhs = static_cast<RCP<SV>>(F_chi.Data.append(F_omega.Data));
                    logger.end();
                    linearSolver.solve();
                    logger.beg("convert from Epetra");
                        for (size_t i = 0; i < m; ++i) chi_BDF1.Data[i] = (*linearSolver.system.lhs)[i];
                    logger.end();
                logger.end();
                logger.beg("save BDF1 soln");
                    for (size_t i = 0; i < m; ++i) omega.Data[i] = (*linearSolver.system.lhs)[i + m];
                    chi_prev = chi;
                    chi = chi_BDF1;
                    //for (size_t i = 0; i < m; ++i) chePot.Data[i] = chemicalPotential(chi_prev.Data[i]);
                    //E_1 = dot(M.Data * f_0(chi_prev.Data),I_p);
                    r_sav =r_sav + dot(chePot.Data, chi.Data - chi_prev.Data) / (2.0 * std::sqrt(E_1));
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