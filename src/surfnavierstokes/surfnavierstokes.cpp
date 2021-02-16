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
#include <unordered_map>

#include "surfactant/ifacetransp.h"
#include "misc/params.h"
#include "geom/builder.h"
#include "levelset/levelset.h"
#include "levelset/adaptriang.h"
#include "levelset/surfacetension.h"
#include "misc/dynamicload.h"
#include "num/bndData.h"
#include "surfactant/ifacetransp.h"

#include "VTKWriter.hpp"
#include "SurfNavierStokesData.hpp"
#include "surfnavierstokes_utils.h"
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
#include "Epetra_OperatorApply.hpp"
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
            system(("mkdir -p " + dirName).c_str());
            system(("rm -r -f " + dirName + "/*").c_str());
            system((
                "mkdir " + dirName + "/matrices ; "
                "mkdir " + dirName + "/vtk ; "
                "mkdir " + dirName + "/stats"
            ).c_str());
            {
                std::ofstream input(dirName + "/input.json");
                input << inpJSON;
                logger.buf << inpJSON;
                logger.log();
            }
            auto h = inpJSON.get<Point3DCL>("Mesh.E1")[0] / inpJSON.get<double>("Mesh.N1") * std::pow(2., -inpJSON.get<double>("Mesh.AdaptRef.FinestLevel"));
            auto tau_u_order  = inpJSON.get<double>("SurfNavStokes.NormalPenaltyPower");
            auto tau_u_factor = inpJSON.get<double>("SurfNavStokes.NormalPenaltyFactor");
            auto tau_u 		  = tau_u_factor * pow(h, tau_u_order); // constant for normal penalty
            auto rho_u_order  = inpJSON.get<double>("SurfNavStokes.VelocityVolumestabPower");
            auto rho_u_factor = inpJSON.get<double>("SurfNavStokes.VelocityVolumestabFactor");
            auto rho_u        = rho_u_factor  * pow(h, rho_u_order); // constant for velocity stabilisation
            auto rho_p_order  = inpJSON.get<double>("SurfNavStokes.PressureVolumestabPower");
            auto rho_p_factor = inpJSON.get<double>("SurfNavStokes.PressureVolumestabFactor");
            auto rho_p        = rho_p_factor  * pow(h, rho_p_order); // constant for pressure stabilisation
            auto numSteps     = inpJSON.get<size_t>("Time.NumbOfSteps");
            auto finalTime    = inpJSON.get<double>("Time.FinalTime");
            auto stepSize     = finalTime / numSteps;
            auto everyStep = inpJSON.get<int>("Output.EveryStep");
            auto testName = inpJSON.get<std::string>("SurfNavStokes.TestName");
            auto surfNavierStokesData = SurfNavierStokesDataFactory(testName, inpJSON);
            FESystem surfOseenSystem;
            auto& gamma = surfOseenSystem.params.surfOseenParams.gamma;
            auto& nu = surfOseenSystem.params.surfOseenParams.nu;
            auto& t = surfOseenSystem.params.t;
            auto& w_T = surfOseenSystem.params.surfOseenParams.w_T;
            auto& Pe = surfOseenSystem.params.surfOseenParams.Pe;
            auto& u_N = surfOseenSystem.params.surfOseenParams.u_N;
            auto& u_N_max = surfOseenSystem.params.surfOseenParams.u_N_max;
            auto& levelSet = surfOseenSystem.params.levelSet;
            nu = inpJSON.get<double>("SurfNavStokes.nu");
            if (nu <= 0.) throw std::invalid_argument("viscosity must be non-negative");
            gamma = inpJSON.get<double>("SurfNavStokes.gamma");
            surfOseenSystem.params.numbOfVirtualSubEdges = inpJSON.get<size_t>("SurfNavStokes.NumbOfVirtualSubEdges");
            surfOseenSystem.params.surfOseenParams.m_g = surfNavierStokesData.m_g;
            surfOseenSystem.params.surfOseenParams.f_T = surfNavierStokesData.f_T;
            auto formulation = inpJSON.get<std::string>("SurfNavStokes.Formulation") == "Consistent" ? "Consistent" : "Inonsistent";
            auto stab = inpJSON.get<std::string>("SurfNavStokes.PressureVolumestabType") == "Normal" ? "Normal" : "Full";
            auto useTangMassMat = inpJSON.get<bool>("SurfNavStokes.UseTangentialMassMatrix");
            auto useInnerIters = inpJSON.get<bool>("Solver.Inner.Use");
            auto usePrevGuess = inpJSON.get<bool>("Solver.UsePreviousFrameAsInitialGuess");
            auto BDF = inpJSON.get<size_t>("SurfNavStokes.BDF");
            if (BDF != 1 && BDF != 2) throw std::invalid_argument("use BDF = 1 or 2");
            logger.buf
                << surfNavierStokesData.description
                << "$\\Delta t$ = " << stepSize << '\n'
                << "BDF: " << BDF << '\n'
                << "$\\nu$ = " << nu << '\n'
                << "$\\gamma$ = " << gamma << " (AL / grad-div stabilization constant)\n"
                << "penalty formulation type: " << formulation << '\n'
                << "pressure volume stabilization type: " << stab << '\n';
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
                DistMarkingStrategyCL markerLset(surfNavierStokesData.surface.phi, inpJSON.get<double>("Mesh.AdaptRef.Width"), inpJSON.get<int>("Mesh.AdaptRef.CoarsestLevel"), inpJSON.get<int>("Mesh.AdaptRef.FinestLevel"));
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
            InitScalar(mg, levelSet, surfNavierStokesData.surface.phi, 0.);
        logger.end();
        logger.beg("set up FE spaces");
            IdxDescCL velExtIdx(vecP2_FE, vecBnd); {
                velExtIdx.GetXidx().SetBound(inpJSON.get<double>("SurfTransp.OmitBound"));
                velExtIdx.CreateNumbering(mg.GetLastLevel(), mg);
            }
            IdxDescCL velIdx(vecP2IF_FE, vecBnd); {
                velIdx.extIdx = &velExtIdx;
                velIdx.GetXidx().SetBound(inpJSON.get<double>("SurfTransp.OmitBound"));
                auto numbOfActiveTetras = velIdx.CreateNumbering(mg.GetLastLevel(), mg, &levelSet, &lsetBnd);
                logger.buf
                    << "numb of active (cut) tetras is: " << numbOfActiveTetras << " ("
                    << (100. * numbOfActiveTetras) / numbOfTetras << "%)\n";
                logger.log();
            }
            IdxDescCL preExtIdx(P1_FE, scaBnd); {
                preExtIdx.GetXidx().SetBound(inpJSON.get<double>("SurfTransp.OmitBound"));
                preExtIdx.CreateNumbering(mg.GetLastLevel(), mg);
            }
            IdxDescCL preIdx(P1IF_FE, scaBnd); {
                preIdx.extIdx = &preExtIdx;
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
                    mM(&velIdx, &velIdx, useTangMassMat ? &LocalAssembler::M_t_vecP2P2 : &LocalAssembler::M_vecP2P2),
                    mAL(&velIdx, &velIdx, &LocalAssembler::AL_vecP2P2),
                    mN(&velIdx, &velIdx, &LocalAssembler::N_vecP2P2),
                    mH(&velIdx, &velIdx, &LocalAssembler::H_vecP2P2),
                    mA(&velIdx, &velIdx, formulation == "Consistent" ? &LocalAssembler::A_consistent_vecP2P2 : &LocalAssembler::A_vecP2P2),
                    mS(&velIdx, &velIdx, &LocalAssembler::S_vecP2P2),
                    mC(&velIdx, &velIdx, &LocalAssembler::C_n_vecP2P2),
                    mB(&preIdx, &velIdx, &LocalAssembler::B_P1vecP2),
                    mQ(&preIdx, &velIdx, &LocalAssembler::Q_P1vecP2),
                    mM_p(&preIdx, &preIdx, &LocalAssembler::M_P1P1),
                    mC_p(&preIdx, &preIdx, stab == "Normal" ? &LocalAssembler::C_n_P1P1 : &LocalAssembler::C_full_P1P1),
                    mA_p(&preIdx, &preIdx, &LocalAssembler::A_P1P1);
                FEVecDescCL
                    vF(&velIdx, &LocalAssembler::F_vecP2),
                    vG(&preIdx, &LocalAssembler::G_P1);
                VecDescCL
                    u_star(&velIdx), u(&velIdx), u_prev(&velIdx), surf_curl_u(&preIdx),
                    p_star(&preIdx), p(&preIdx);
            logger.end();
        logger.end();
        logger.beg("set up vtk");
            VTKWriter vtkWriter(dirName + "/vtk/" + testName, mg, inpJSON.get<bool>("Output.Binary"));
            vtkWriter.add({ "level-set", levelSet });
            if (everyStep > 0) {
                if (inpJSON.get<bool>("Output.Velocity")) {
                    vtkWriter.add({ "u_h", u });
                    if (surfNavierStokesData.exactSoln) vtkWriter.add({ "u_*", u_star });
                }
                if (inpJSON.get<bool>("Output.Vorticity")) vtkWriter.add({ "w_h", surf_curl_u });
                if (inpJSON.get<bool>("Output.Pressure")) {
                    vtkWriter.add({ "p_h", p });
                    if (surfNavierStokesData.exactSoln) vtkWriter.add({ "p_*", p_star });
                }
            }
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
            belosParams->set("Convergence Tolerance", inpJSON.get<double>("Solver.Outer.RelResTol"));
            belosParams->set( "Output Frequency", inpJSON.get<int>("Solver.Outer.OutputFrequency"));
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
            MT A(Epetra_DataAccess::Copy, mapVelocity, 0, true), B(Epetra_DataAccess::Copy, mapPressure, 0, true), C(Epetra_DataAccess::Copy, mapPressure, 0, true), S_L(Epetra_DataAccess::Copy, mapPressure, 0, true), S_M(Epetra_DataAccess::Copy, mapPressure, 0, true);
            RCP<OP> belosMTX = rcp(new Epetra_OperatorApply([&](MV const & X, MV& Y) {
                double* view;
                X(0)->ExtractView(&view);
                SV x1(Epetra_DataAccess::View, mapVelocity, view);
                SV x2(Epetra_DataAccess::View, mapPressure, view + n);
                Y(0)->ExtractView(&view);
                SV y11(Epetra_DataAccess::View, mapVelocity, view);
                auto y12 = y11;
                SV y21(Epetra_DataAccess::View, mapPressure, view + n);
                auto y22 = y21;
                // y1
                A.Multiply(false, x1, y11);
                B.Multiply(true, x2, y12);
                y11.Update(1., y12, 1.);
                // y2
                B.Multiply(false, x1, y21);
                C.Multiply(false, x2, y22);
                y21.Update(1., y22, 1.);
            }));
            logger.beg("set up preconditioners");
                auto identity = [](MV const & X, MV& Y) {
                    double* view;
                    Y(0)->ExtractView(&view);
                    X(0)->ExtractCopy(view);
                };
                logger.beg("diffusion-convection-reaction block");
                    size_t numItersA = 0;
                    Epetra_OperatorApply::ApplyType invA = identity;
                    logger.log(inpJSON.get<std::string>("Solver.Inner.A.Comment"));
                    auto iterationA = inpJSON.get<std::string>("Solver.Inner.A.Iteration");
                    if (!Amesos2::query(iterationA)) throw std::invalid_argument("solver " + iterationA + " is not available for Amesos2");
                    logger.log("factorizing full velocity matrix");
                    auto amesosSolver = Amesos2::create<MT, MV>(iterationA, rcpFromRef(A));
                    logger.buf
                        << "solver:      " << amesosSolver->name() << '\n'
                        << "description: " << amesosSolver->description();
                    logger.log();
                    auto factorizationTime = 0.;
                    auto runFactorization = [&]() {
                        logger.beg("factorization");
                            logger.beg("symbolic factorization");
                                amesosSolver->symbolicFactorization();
                            logger.end();
                            logger.beg("numeric factorization");
                                amesosSolver->numericFactorization();
                                auto amesosStatus = amesosSolver->getStatus();
                                logger.buf
                                    << "numb of nonzeros in L + U = " << amesosStatus.getNnzLU() << " ("
                                    << (100. * amesosStatus.getNnzLU()) / (static_cast<double>(n) * n) << "%)";
                                logger.log();
                            logger.end();
                        factorizationTime = logger.end();
                    };
                    invA = [&](MV const &X, MV &Y) {
                        amesosSolver->setB(rcpFromRef(X));
                        amesosSolver->setX(rcpFromRef(Y));
                        amesosSolver->solve();
                        numItersA++;
                    };
                logger.end();
                logger.beg("schur complement block");
                    size_t numItersS_M = 0, numItersS_L = 0;
                    Epetra_OperatorApply::ApplyType invS = identity;
                    auto belosParamsS = parameterList();
                    belosParamsS->set("Maximum Iterations", inpJSON.get<int>("Solver.Inner.S.MaxIter"));
                    belosParamsS->set("Convergence Tolerance", inpJSON.get<double>("Solver.Inner.S.RelResTol"));
                    auto belosSolverS_M = belosFactory.create(inpJSON.get<std::string>("Solver.Inner.S.Iteration"), belosParamsS);
                    auto belosSolverS_L = belosFactory.create(inpJSON.get<std::string>("Solver.Inner.S.Iteration"), belosParamsS);
                    logger.buf << "inner solver: " << belosSolverS_M->description();
                    logger.log();
                    invS = [&](MV const & X, MV& Y) {
                        // ini guess
                        Y.PutScalar(0.);
                        auto& Y_M = Y;
                        auto  Y_L = Y;
                        // normalized rhs
                        auto X_nrm = X; {
                            double mean;
                            X_nrm.MeanValue(&mean);
                            for (size_t i = 0; i < m; ++i) (*X_nrm(0))[i] -= mean;
                        }
                        // Y_M
                        LinearProblem<ST, MV, OP> belosProblemS_M(rcpFromRef(S_M), rcpFromRef(Y_M), rcpFromRef(X_nrm));
                        belosProblemS_M.setProblem();
                        belosSolverS_M->setProblem(rcpFromRef(belosProblemS_M));
                        belosSolverS_M->solve();
                        numItersS_M += belosSolverS_M->getNumIters();
                        // Y_L (TODO: make parallel)
                        if (inpJSON.get<bool>("Solver.Inner.S.S_L")) {
                            LinearProblem<ST, MV, OP> belosProblemS_L(rcpFromRef(S_L), rcpFromRef(Y_L), rcpFromRef(X_nrm));
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
                RCP<OP> belosPRE;
                if (precType == "BlockDiagonal") {
                    logger.log("using block-diagonal preconditioner");
                    belosPRE = rcp(new Epetra_OperatorApply([&](MV const & X, MV& Y) {
                        double* view;
                        X(0)->ExtractView(&view);
                        SV x1(Epetra_DataAccess::View, mapVelocity, view);
                        SV x2(Epetra_DataAccess::View, mapPressure, view + n);
                        Y(0)->ExtractView(&view);
                        SV y1(Epetra_DataAccess::View, mapVelocity, view);
                        SV y2(Epetra_DataAccess::View, mapPressure, view + n);
                        #ifdef _OPENMP
                            #pragma omp parallel num_threads(2)
                            {
                                #pragma omp single
                                {
                                    #pragma omp task
                                        invA(x1, y1); // y1
                                    #pragma omp task
                                        invS(x2, y2); // y2
                                }
                            }
                        #else
                            invA(x1, y1);
                            invS(x2, y2);
                        #endif
                    }));
                }
                else {
                    logger.log("using block-triangular preconditioner");
                    belosPRE = rcp(new Epetra_OperatorApply([&](MV const & X, MV& Y) {
                        double* view;
                        X(0)->ExtractView(&view);
                        SV x1(Epetra_DataAccess::View, mapVelocity, view);
                        SV x2(Epetra_DataAccess::View, mapPressure, view + n);
                        Y(0)->ExtractView(&view);
                        SV y1(Epetra_DataAccess::View, mapVelocity, view);
                        auto y1_tmp = y1;
                        SV y2(Epetra_DataAccess::View, mapPressure, view + n);
                        // y2
                        invS(x2, y2);
                        // y1
                        B.Multiply(true, y2, y1_tmp);
                        y1_tmp.Update(1., x1, -1.);
                        invA(y1_tmp, y1);
                    }));
                }
            logger.end();
            LinearProblem<ST, MV, OP> belosProblem(belosMTX, rcpFromRef(belosLHS), rcpFromRef(belosRHS));
            if (useInnerIters) belosProblem.setRightPrec(belosPRE);
        logger.end();
        t = 0.;
        logger.beg("t = t_0 = 0");
            logger.beg("assemble");
                surfOseenSystem.matrices = { &mM, &mA, &mAL, &mS, &mC, &mM_p, &mC_p, &mA_p, &mB, &mQ };
                setupFESystem(mg, surfOseenSystem);
                mC_p.Data *= -rho_p;
                MatrixCL B_T, A_sum;
                transpose(mB.Data, B_T);
                VectorCL I_p(1., m);
            auto assembleTime = logger.end();
            logger.beg("interpolate initial data");
                InitVector(mg, u, surfNavierStokesData.u_T, t);
                InitScalar(mg, p, surfNavierStokesData.p, t);
                p.Data -= dot(mM_p.Data * p.Data, I_p) / dot(mM_p.Data * I_p, I_p) * I_p;
                u_prev = u;
            logger.end();
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
                    auto mtx = static_cast<MT>(MatrixCL(1., mM_p.Data, -1., mC_p.Data));
                    auto rhs = static_cast<SV>(mQ.Data * u.Data);
                    auto sln = static_cast<SV>(surf_curl_u.Data);
                    sln.PutScalar(0.);
                    LinearProblem<ST, MV, OP> belosProblemW(rcpFromRef(mtx), rcpFromRef(sln), rcpFromRef(rhs));
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
                double alpha, b_norm, r0_norm, solveTime, solveWastedTime;
                ReturnType belosSolverResult;
                auto residual = [&](VectorCL const & u, VectorCL const & p) {
                    auto velResSq = norm_sq(A_sum * u + B_T * p - vF.Data);
                    auto preResSq = norm_sq(mB.Data * u + mC_p.Data * p - vG.Data);
                    return std::tuple<double, double, double>(sqrt(velResSq), sqrt(preResSq), sqrt(velResSq + preResSq));
                };
                auto exportStats = [&](size_t i) {
                    std::ofstream stats(dirName + "/stats/t_" + std::to_string(i) + ".json");
                    ParamCL tJSON;
                    tJSON.put("Time", t);
                    tJSON.put("h", h);
                    tJSON.put("nu", nu);
                    tJSON.put("gamma", gamma);
                    tJSON.put("MeshDepParams.rho_p", rho_p);
                    tJSON.put("MeshDepParams.rho_u", rho_u);
                    tJSON.put("MeshDepParams.tau_u", tau_u);
                    tJSON.put("ElapsedTime.Assemble", assembleTime);
                    tJSON.put("ElapsedTime.ProjectVorticity", projectTime);
                    tJSON.put("Solver.ProjectVorticity.TotalIters", belosSolverW->getNumIters());
                    tJSON.put("Solver.ProjectVorticity.Converged", belosSolverResultW == Converged);
                    tJSON.put("Solver.ProjectVorticity.ResidualNormRelative", belosSolverW->achievedTol());
                    auto surfArea = dot(I_p, mM_p.Data * I_p);
                    tJSON.put("Integral.SurfacaArea", surfArea);
                    tJSON.put("Integral.PressureMean", dot(I_p, mM_p.Data * p.Data) / surfArea);
                    tJSON.put("Integral.PressureL2", sqrt(dot(p.Data, mM_p.Data * p.Data)));
                    auto velL2Sq = dot(u.Data, mM.Data * u.Data);
                    tJSON.put("Integral.KineticEnergy", .5 * velL2Sq);
                    tJSON.put("Integral.VelocityL2", sqrt(velL2Sq));
                    tJSON.put("Integral.VelocitySurfaceDivergenceL2", sqrt(dot(u.Data, mAL.Data * u.Data)));
                    tJSON.put("Integral.VelocityNormalL2", sqrt(dot(u.Data, mS.Data * u.Data)));
                    auto vorL2Sq = dot(surf_curl_u.Data, mM_p.Data * surf_curl_u.Data);
                    tJSON.put("Integral.Enstrophy", .5 * vorL2Sq);
                    tJSON.put("Integral.SurfaceVorticityL2", sqrt(vorL2Sq));
                    auto vorH1Sq = dot(surf_curl_u.Data, mA_p.Data * surf_curl_u.Data);
                    tJSON.put("Integral.Palinstrophy", .5 * vorH1Sq);
                    tJSON.put("Integral.SurfaceVorticityH1", sqrt(vorH1Sq));
                    if (surfNavierStokesData.exactSoln) {
                        InitVector(mg, u_star, surfNavierStokesData.u_T, t);
                        InitScalar(mg, p_star, surfNavierStokesData.p, t);
                        p_star.Data -= dot(mM_p.Data * p_star.Data, I_p) / dot(mM_p.Data * I_p, I_p) * I_p;
                        tJSON.put("Integral.Exact.PressureMean", dot(I_p, mM_p.Data * p_star.Data) / surfArea);
                        tJSON.put("Integral.Exact.PressureL2", sqrt(dot(p_star.Data, mM_p.Data * p_star.Data)));
                        auto velL2Sq = dot(u_star.Data, mM.Data * u_star.Data);
                        tJSON.put("Integral.Exact.KineticEnergy", .5 * velL2Sq);
                        tJSON.put("Integral.Exact.VelocityL2", sqrt(velL2Sq));
                        tJSON.put("Integral.Exact.VelocitySurfaceDivergenceL2", sqrt(dot(u_star.Data, mAL.Data * u_star.Data)));
                        tJSON.put("Integral.Exact.VelocityNormalL2", sqrt(dot(u_star.Data, mS.Data * u_star.Data)));
                        auto u_diff = u_star.Data - u.Data;
                        tJSON.put("Integral.Error.VelocityL2", sqrt(dot(u_diff, mM.Data * u_diff)));
                        tJSON.put("Integral.Error.VelocityH1", sqrt(dot(u_diff, mA.Data * u_diff)));
                        auto p_diff = p_star.Data - p.Data;
                        tJSON.put("Integral.Error.PressureL2", sqrt(dot(p_diff, mM_p.Data * p_diff)));
                        if (i > 0) {
                            auto && [velRes, preRes, fullRes] = residual(u_star.Data, p_star.Data);
                            tJSON.put("Solver.Outer.ResidualNorm.ExactSolnAbsolute.Velocity", velRes);
                            tJSON.put("Solver.Outer.ResidualNorm.ExactSolnAbsolute.Pressure", preRes);
                            tJSON.put("Solver.Outer.ResidualNorm.ExactSolnAbsolute.Full", fullRes);
                        }
                    }
                    if (i > 0) {
                        auto && [velRes, preRes, fullRes] = residual(u.Data, p.Data);
                        tJSON.put("Solver.Outer.ResidualNorm.TrueAbsolute.Velocity", velRes);
                        tJSON.put("Solver.Outer.ResidualNorm.TrueAbsolute.Pressure", preRes);
                        tJSON.put("Solver.Outer.ResidualNorm.TrueAbsolute.Full", fullRes);
                        tJSON.put("Solver.Outer.ResidualNorm.b", b_norm);
                        tJSON.put("Solver.Outer.ResidualNorm.r_0", r0_norm);
                        tJSON.put("Solver.Outer.ResidualNorm.r_i/b", fullRes / b_norm);
                        tJSON.put("Solver.Outer.ResidualNorm.r_i/r_0", fullRes / r0_norm);
                        tJSON.put("Solver.Outer.ResidualNorm.SolverRelative", belosSolver->achievedTol());
                        tJSON.put("Solver.Outer.TotalIters", belosSolver->getNumIters());
                        tJSON.put("Solver.Outer.DOF", n + m);
                        tJSON.put("Solver.Outer.Converged", belosSolverResult == Converged);
                        tJSON.put("ElapsedTime.LinearSolve", solveTime);
                        tJSON.put("ElapsedTime.LinearSolveWasted", solveWastedTime);
                        tJSON.put("ElapsedTime.Factorization", factorizationTime);
                        tJSON.put("Solver.Inner.A.TotalIters", numItersA);
                        tJSON.put("Solver.Inner.A.MeanIters", static_cast<double>(numItersA) / belosSolver->getNumIters());
                        tJSON.put("Solver.Inner.A.DOF", n);
                        tJSON.put("Solver.Inner.S.TotalIters.S_M", numItersS_M);
                        tJSON.put("Solver.Inner.S.TotalIters.S_L", numItersS_L);
                        tJSON.put("Solver.Inner.S.MeanIters.S_M", static_cast<double>(numItersS_M) / belosSolver->getNumIters());
                        tJSON.put("Solver.Inner.S.MeanIters.S_L", static_cast<double>(numItersS_L) / belosSolver->getNumIters());
                        tJSON.put("Solver.Inner.S.DOF", m);
                        tJSON.put("Peclet", Pe);
                        tJSON.put("MaxSurfaceSpeed", u_N_max);
                        tJSON.put("MassMatrixCoef", alpha);
                        if (inpJSON.get<bool>("SurfNavStokes.ExportMatrices")) {
                            std::string format = inpJSON.get<std::string>("SurfNavStokes.ExportMatricesFormat") == "mtx" ? ".mtx" : ".mat";
                            auto expFunc = format == ".mtx" ? &MatrixCL::exportMTX : &MatrixCL::exportMAT;
                            auto expMat = [&](MatrixCL &A, std::string const a, std::string const &b) {
                                logger.beg(a);
                                logger.buf << "size: " << A.num_rows() << "x" << A.num_cols();
                                logger.log();
                                (A.*expFunc)(dirName + "/matrices/" + b + format);
                                tJSON.put("Matrices." + a, "../matrices/" + b + format);
                                logger.end();
                            };
                            expMat(mA_p.Data, "PressureStiffness", "A_p");
                            logger.log();
                        }
                    }
                    auto vtkTime = 0.;
                    if (everyStep > 0 && (i == 0 || i == numSteps || (i - 1) % everyStep == 0)) {
                        logger.beg("write vtk");
                            vtkWriter.write(t);
                        vtkTime = logger.end();
                    }
                    tJSON.put("ElapsedTime.VTK", vtkTime);
                    stats << tJSON;
                    logger.buf << tJSON;
                    logger.log();
                    numItersA = numItersS_M = numItersS_L = 0;
                    factorizationTime = solveWastedTime = 0.;
                };
                exportStats(0);
            logger.end();
        logger.end();
        for (size_t i = 1; i <= numSteps; ++i) {
            t = i * stepSize;
            std::stringstream header; header << "t = t_" << i << " = " << t << " (" << (100. * i) / numSteps << "%)";
            logger.beg(header.str());
                logger.beg("initialize surface speed $u_N$");
                    InitScalar(mg, u_N, surfNavierStokesData.surface.u_N, t);
                logger.end();
                logger.beg("assemble");
                    surfOseenSystem.matrices = {};
                    if (surfNavierStokesData.w_T) InitVector(mg, w_T, surfNavierStokesData.w_T, t); // Oseen (and Stokes) case
                    else w_T.Data = 2. * u.Data - u_prev.Data; // Navier-Stokes case
                    Pe = supnorm(w_T.Data) / nu;
                    if (Pe) {
                        logger.log("assembling convection mtx");
                        surfOseenSystem.matrices.push_back(&mN);
                    }
                    u_N_max = supnorm(u_N.Data);
                    if (u_N_max) {
                        logger.log("assembling surf speed mtx");
                        surfOseenSystem.matrices.push_back(&mH);
                    }
                    alpha = i == 1 || BDF == 1 ? 1. / stepSize : 1.5 / stepSize;
                    logger.buf
                        << "$\\alpha$ = " << alpha << '\n'
                        << "Pe = " << Pe << '\n'
                        << "max |u_N| = " << u_N_max;
                    logger.log();
                    surfOseenSystem.vectors = { &vF, &vG };
                    setupFESystem(mg, surfOseenSystem);
                    // system mtx
                    A_sum.LinComb(alpha, mM.Data, gamma, mAL.Data, nu, mA.Data, tau_u, mS.Data, rho_u, mC.Data);
                    if (Pe) A_sum.LinComb(1., MatrixCL(A_sum), 1., mN.Data);
                    if (u_N_max) A_sum.LinComb(1., MatrixCL(A_sum), 1., mH.Data);
                    // system rhs
                    if (i == 1 || BDF == 1) vF.Data += (1. / stepSize) * (mM.Data * u.Data);
                    else vF.Data += (2. / stepSize) * (mM.Data * u.Data) - (.5 / stepSize) * (mM.Data * u_prev.Data);
                    vG.Data -= (dot(vG.Data, I_p) / dot(I_p, I_p)) * I_p;
                    // for the next step
                    u_prev = u;
                assembleTime = logger.end();
                logger.beg("convert to Epetra");
                    A = static_cast<MT>(A_sum);
                    auto printStat = [&](std::string const & name, MT const & A) {
                        logger.buf << name << ": " << A.NumGlobalRows() << 'x' << A.NumGlobalCols() << ", " << A.NumGlobalNonzeros() << " nonzeros";
                        if (A.NumGlobalNonzeros()) logger.buf << " (" << (100. * A.NumGlobalNonzeros()) / (static_cast<double>(A.NumGlobalRows()) * A.NumGlobalCols()) << "%)";
                        logger.log();
                    };
                    printStat("A", A);
                    B = static_cast<MT>(mB.Data);
                    printStat("B", B);
                    C = static_cast<MT>(mC_p.Data);
                    printStat("C := -rho_p (pressure volume stabilization mtx)", C);
                    S_M = static_cast<MT>(MatrixCL(1. / (gamma + nu), mM_p.Data, -1., mC_p.Data));
                    printStat("S_M := (\\nu + \\gamma)^{-1} M_p - C", S_M);
                    S_L = static_cast<MT>(MatrixCL(1. / alpha, mA_p.Data, -1., mC_p.Data));
                    printStat("S_L := \\alpha^{-1} L_p - C", S_L);
                    belosRHS = static_cast<SV>(vF.Data.append(vG.Data));
                logger.end();
                if (useInnerIters && i == 1) runFactorization();
                logger.beg("linear solve");
                    b_norm = sqrt(norm_sq(vF.Data) + norm_sq(vG.Data));
                    if (usePrevGuess) {
                        r0_norm = std::get<2>(residual(u.Data, p.Data));
                        belosLHS = static_cast<SV>(u.Data.append(p.Data));
                        belosParams->set("Convergence Tolerance", inpJSON.get<double>("Solver.Outer.RelResTol") * b_norm / r0_norm);
                        belosSolver = belosFactory.create(inpJSON.get<std::string>("Solver.Outer.Iteration"), belosParams);
                    } else {
                        r0_norm = b_norm;
                        belosLHS.PutScalar(0.);
                    }
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
                        runFactorization();
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
                    p.Data -= dot(mM_p.Data * p.Data, I_p) / dot(mM_p.Data * I_p, I_p) * I_p;
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