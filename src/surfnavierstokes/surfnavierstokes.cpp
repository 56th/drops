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
// amesos (sparse-direct solvers)
#include "Amesos.h"
#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
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
        logger.end();
        logger.beg("set up test case");
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
            auto numbOfSteps  = inpJSON.get<size_t>("Time.NumbOfSteps");
            auto finalTime    = inpJSON.get<double>("Time.FinalTime");
            auto stepSize     = finalTime / numbOfSteps;
            auto everyStep = inpJSON.get<int>("Output.EveryStep");
            logger.buf << "$\\Delta t$ = " << stepSize << '\n';
            // parse FE types and some other parameters from json file
            auto FE = inpJSON.get<std::string>("SurfNavStokes.FE");
            if (FE != "P2P1")
                throw std::invalid_argument("TODO: refactor P1P1 assembly");
            auto velFE = FE.substr(0,2);
            auto preFE = FE.substr(2,2);
            logger.buf
                << "velocity FE = " << velFE << '\n'
                << "pressure FE = " << preFE << '\n';
            auto testName = inpJSON.get<std::string>("SurfNavStokes.TestName");
            auto nu = inpJSON.get<double>("SurfNavStokes.nu");
            auto gamma = inpJSON.get<double>("SurfNavStokes.gamma");
            auto surfNavierStokesData = SurfNavierStokesDataFactory(testName, nu, inpJSON);
            logger.buf << surfNavierStokesData.description;
            logger.buf << "$\\nu$ = " << nu << '\n';
            logger.buf << "$\\gamma$ = " << gamma << " (AL / grad-div stabilization constant)\n";
            SurfOseenParam param;
            param.input.f = inpJSON.get<bool>("SurfNavStokes.IntegrateRHS") ? surfNavierStokesData.f_T : nullptr;
            param.input.g = inpJSON.get<bool>("SurfNavStokes.IntegrateRHS") ? surfNavierStokesData.m_g : nullptr;
            param.input.exactNormal = inpJSON.get<bool>("SurfNavStokes.ComputeNormalErr") ? surfNavierStokesData.surface.n : nullptr;
            param.input.exactShape = inpJSON.get<bool>("SurfNavStokes.ComputeShapeErr") ? surfNavierStokesData.surface.H : nullptr;
            param.input.numbOfVirtualSubEdges = inpJSON.get<size_t>("SurfNavStokes.NumbOfVirtualSubEdges");
            param.input.usePatchNormal = inpJSON.get<bool>("SurfNavStokes.UsePatchNormals");
            if (param.input.usePatchNormal)  logger.log("using patch normals in surf integrands");
            auto formulation = inpJSON.get<std::string>("SurfNavStokes.Formulation");
            if (formulation == "consistent") {
                param.input.formulation = SurfOseenParam::Formulation::consistent;
                logger.buf << "using consistent penalty formulation\n";
            }
            else {
                param.input.formulation = SurfOseenParam::Formulation::inconsistent;
                logger.buf << "using inconsistent penalty formulation\n";
            }
            SurfOseenSystem surfOseenSystem;
            auto stab = inpJSON.get<std::string>("SurfNavStokes.PressureVolumestabType");
            if (stab == "full") {
                param.input.stab = SurfOseenParam::PressureVolumestab::full;
                logger.buf << "using full pressure volume stabilization\n";
            }
            else {
                param.input.stab = SurfOseenParam::PressureVolumestab::normal;
                logger.buf << "using normal pressure volume stabilization\n";
            }
            if (!inpJSON.get<bool>("SurfNavStokes.ExportMatrices")) {
                surfOseenSystem.LB.assemble = false;
                surfOseenSystem.LB_stab.assemble = false;
            }
            if (gamma == 0.)
                surfOseenSystem.AL.assemble = false;
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
                mg.Transform([&](Point3DCL const & p) {
                    return p - shift;
                });
                logger.buf << "surface shift = " << shift;
                logger.log();
            logger.end();
            logger.beg("refine towards the surface");
                AdapTriangCL adap(mg);
                // adaptive mesh refinement based on level set function
                DistMarkingStrategyCL initmarker(surfNavierStokesData.surface.phi, inpJSON.get<double>("Mesh.AdaptRef.Width"), inpJSON.get<int>("Mesh.AdaptRef.CoarsestLevel"), inpJSON.get<int>("Mesh.AdaptRef.FinestLevel"));
                adap.set_marking_strategy(&initmarker);
                adap.MakeInitialTriang();
                adap.set_marking_strategy(0);
                auto numbOfTetras = mg.GetNumTriangTetra();
                logger.buf
                    << "h = " << h << '\n'
                    << "numb of tetras = " << numbOfTetras;
                logger.log();
            logger.end();
        logger.end();
        logger.beg("interpolate level-set");
            InstatScalarFunction sigma(0);
            SurfaceTensionCL sf(sigma, 0);
            BndDataCL<double> lsbnd(0);
            read_BndData(lsbnd, mg, inpJSON.get_child("Levelset.BndData"));
            std::shared_ptr<LevelsetP2CL> lsetPtr(LevelsetP2CL::Create(mg, lsbnd, sf));
            auto& lset = *lsetPtr;
            lset.CreateNumbering(mg.GetLastLevel(), &lset.idx);
            lset.Phi.SetIdx(&lset.idx);
            lset.Init(surfNavierStokesData.surface.phi);
        logger.end();
        logger.beg("build FE spaces");
            BndDataCL<Point3DCL> vbnd(0);
            read_BndData(vbnd, mg, inpJSON.get_child("Stokes.VelocityBndData"));
            BndDataCL<double> pbnd(0);
            read_BndData(pbnd, mg, inpJSON.get_child("Stokes.PressureBndData"));
            IdxDescCL ifaceVecP2idx(vecP2IF_FE, vbnd);
            IdxDescCL ifaceVecP1idx(vecP1IF_FE, vbnd);
            IdxDescCL ifaceP1idx(P1IF_FE, pbnd);
            IdxDescCL ifaceP2idx(P2IF_FE, pbnd);
            IdxDescCL vecP2idx(vecP2_FE, vbnd);
            IdxDescCL vecP1idx(vecP1_FE, vbnd);
            IdxDescCL P1FEidx(P1_FE, pbnd);
            IdxDescCL P2FEidx(P2_FE, pbnd);
            ifaceVecP2idx.GetXidx().SetBound(inpJSON.get<double>("SurfTransp.OmitBound"));
            auto numbOfActiveTetrasP2 = ifaceVecP2idx.CreateNumbering(mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
            ifaceVecP1idx.GetXidx().SetBound(inpJSON.get<double>("SurfTransp.OmitBound"));
            auto numbOfActiveTetrasP1 = ifaceVecP1idx.CreateNumbering(mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
            if (numbOfActiveTetrasP2 != numbOfActiveTetrasP1) {
                std::stringstream err;
                err << "inconsistent numb of cut tetras for P2 and P1: " << numbOfActiveTetrasP2 << " vs. " << numbOfActiveTetrasP1;
                throw std::logic_error(err.str());
            }
            logger.buf << "numb of active (cut) tetras is: " << numbOfActiveTetrasP1 << " (" << (100. * numbOfActiveTetrasP1) / numbOfTetras << "%)\n";
            logger.log();
            ifaceP1idx.GetXidx().SetBound(inpJSON.get<double>("SurfTransp.OmitBound"));
            ifaceP1idx.CreateNumbering(mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
            ifaceP2idx.GetXidx().SetBound(inpJSON.get<double>("SurfTransp.OmitBound"));
            ifaceP2idx.CreateNumbering(mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
            vecP2idx.GetXidx().SetBound(inpJSON.get<double>("SurfTransp.OmitBound"));
            vecP2idx.CreateNumbering(mg.GetLastLevel(), mg);
            vecP1idx.GetXidx().SetBound(inpJSON.get<double>("SurfTransp.OmitBound"));
            vecP1idx.CreateNumbering(mg.GetLastLevel(), mg);
            P1FEidx.GetXidx().SetBound(inpJSON.get<double>("SurfTransp.OmitBound"));
            P1FEidx.CreateNumbering(mg.GetLastLevel(), mg);
            P2FEidx.GetXidx().SetBound(inpJSON.get<double>("SurfTransp.OmitBound"));
            P2FEidx.CreateNumbering(mg.GetLastLevel(), mg);
            logger.buf
                    << "numb of d.o.f. vector interface P2: " << ifaceVecP2idx.NumUnknowns() << '\n'
                    << "numb of d.o.f. vector interface P1: " << ifaceVecP1idx.NumUnknowns() << '\n'
                    << "numb of d.o.f. scalar interface P2: " << ifaceP2idx.NumUnknowns() << '\n'
                    << "numb of d.o.f. scalar interface P1: " << ifaceP1idx.NumUnknowns() << '\n'
                    << "numb of d.o.f. vector P2: " << vecP2idx.NumUnknowns() << '\n'
                    << "numb of d.o.f. vector P1: " << vecP1idx.NumUnknowns() << '\n'
                    << "numb of d.o.f. scalar P2: " << P2FEidx.NumUnknowns() << '\n'
                    << "numb of d.o.f. scalar P1: " << P1FEidx.NumUnknowns();
            logger.log();
            VecDescCL u_star, u_star_ext, u, u_ext, u_prev, u_prev_prev, p_star, p_star_ext, p, p_ext;
            size_t n, m;
            if (FE == "P2P1") {
                n = ifaceVecP2idx.NumUnknowns();
                m = ifaceP1idx.NumUnknowns();
                u.SetIdx(&ifaceVecP2idx);
                u_ext.SetIdx(&vecP2idx);
                u_star.SetIdx(&ifaceVecP2idx);
                u_star_ext.SetIdx(&vecP2idx);
                u_prev.SetIdx(&ifaceVecP2idx);
                u_prev_prev.SetIdx(&ifaceVecP2idx);
                surfOseenSystem.w.SetIdx(&ifaceVecP2idx);
                surfOseenSystem.fRHS.SetIdx(&ifaceVecP2idx);
                surfOseenSystem.alRHS.SetIdx(&ifaceVecP2idx);
                surfOseenSystem.A.SetIdx(&ifaceVecP2idx, &ifaceVecP2idx);
                surfOseenSystem.AL.SetIdx(&ifaceVecP2idx, &ifaceVecP2idx);
                surfOseenSystem.N.SetIdx(&ifaceVecP2idx, &ifaceVecP2idx);
                surfOseenSystem.A_stab.SetIdx(&ifaceVecP2idx, &ifaceVecP2idx);
                surfOseenSystem.B.SetIdx(&ifaceP1idx, &ifaceVecP2idx);
                surfOseenSystem.M.SetIdx(&ifaceVecP2idx, &ifaceVecP2idx);
                surfOseenSystem.S.SetIdx(&ifaceVecP2idx, &ifaceVecP2idx);
            }
            else throw std::invalid_argument(FE + ": unknown FE pair");
            p.SetIdx(&ifaceP1idx);
            p_ext.SetIdx(&P1FEidx);
            p_star.SetIdx(&ifaceP1idx);
            p_star_ext.SetIdx(&P1FEidx);
            surfOseenSystem.gRHS.SetIdx(&ifaceP1idx);
            surfOseenSystem.M_p.SetIdx(&ifaceP1idx, &ifaceP1idx);
            surfOseenSystem.C.SetIdx(&ifaceP1idx, &ifaceP1idx);
            logger.beg("interpolate initial data");
                InitVector(mg, u_star, surfNavierStokesData.u_T, 0.);
                InitScalar(mg, p_star, surfNavierStokesData.p, 0.);
                u = u_star;
                p = p_star;
            logger.end();
        logger.end();
        logger.beg("set up vtk");
            VTKWriter vtkWriter(dirName + "/vtk/" + testName, mg, inpJSON.get<bool>("Output.Binary"));
            auto writeVTK = [&](double t) {
                Extend(mg, u_star, u_star_ext);
                Extend(mg, u, u_ext);
                Extend(mg, p_star, p_star_ext);
                Extend(mg, p, p_ext);
                vtkWriter.write(t);
            };
            if (everyStep > 0) {
                VTKWriter::VTKVar vtkLevelSet;
                vtkLevelSet.name = "level-set";
                vtkLevelSet.value = &lset.Phi.Data;
                vtkLevelSet.type = VTKWriter::VTKVar::Type::P2;
                vtkWriter
                    .add(vtkLevelSet)
                    .add(VTKWriter::VTKVar({ "u_h", &u_ext.Data, VTKWriter::VTKVar::Type::vecP2 }))
                    .add(VTKWriter::VTKVar({ "p_h", &p_ext.Data, VTKWriter::VTKVar::Type::P1 }));
                if (surfNavierStokesData.exactSoln)
                    vtkWriter
                        .add(VTKWriter::VTKVar({ "u_*", &u_star_ext.Data, VTKWriter::VTKVar::Type::vecP2 }))
                        .add(VTKWriter::VTKVar({ "p_*", &p_star_ext.Data, VTKWriter::VTKVar::Type::P1 }));
                logger.beg("write initial condition to vtk");
                    writeVTK(0.);
                logger.end();
            }
        logger.end();
        logger.beg("set up linear solver");
            using MV = Epetra_MultiVector;
            using OP = Epetra_Operator;
            using ST = double;
            using namespace Teuchos;
            using namespace Belos;
            RCP<ParameterList> belosParams = parameterList();
            belosParams->set("Num Blocks", inpJSON.get<int>("Solver.Outer.KrylovSubspaceSize"));
            belosParams->set("Maximum Iterations", inpJSON.get<int>("Solver.Outer.MaxIter"));
            belosParams->set("Convergence Tolerance", inpJSON.get<double>("Solver.Outer.RelResTol"));
            belosParams->set( "Output Frequency", inpJSON.get<int>("Solver.Outer.OutputFrequency"));
            belosParams->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails + Belos::TimingDetails + Belos::FinalSummary + Belos::IterationDetails);
            belosParams->set<int>("Output Style", Belos::Brief);
            SolverFactory<ST, MV, OP> belosFactory;
            logger.buf << "available iterations: " << belosFactory.supportedSolverNames();
            logger.log();
            auto belosSolver = belosFactory.create(inpJSON.get<std::string>("Solver.Outer.Iteration"), belosParams);
            logger.buf << "outer solver: " << belosSolver->description();
            logger.log();
        logger.end();
        logger.beg("t = t_1");
            logger.beg("assemble");
                InitVector(mg, u_prev, surfNavierStokesData.u_T, 0.);
                param.input.t = stepSize;
                auto setWind = [&](VectorCL const & wind) {
                    if (surfNavierStokesData.w_T) // Oseen (and Stokes) case
                        InitVector(mg, surfOseenSystem.w, surfNavierStokesData.w_T, param.input.t);
                    else surfOseenSystem.w.Data = wind; // Navier-Stokes case
                };
                setWind(u_prev.Data);
                auto Pe = supnorm(surfOseenSystem.w.Data) / nu;
                logger.buf << "Pe = " << Pe;
                logger.log();
                SetupSurfOseen_P2P1(mg, lset, &surfOseenSystem, &param);
                VectorCL I_p;
                I_p.resize(m, 1.);
                surfOseenSystem.gRHS.Data -= (dot(surfOseenSystem.gRHS.Data, I_p) / dot(I_p, I_p)) * I_p;
                surfOseenSystem.fRHS.Data += (1. / stepSize) * (surfOseenSystem.M.Data * u_prev.Data) - gamma * surfOseenSystem.alRHS.Data;
                surfOseenSystem.sumA.Data.LinComb(1. / stepSize, surfOseenSystem.M.Data, gamma, surfOseenSystem.AL.Data, 1., surfOseenSystem.N.Data, nu, surfOseenSystem.A.Data, tau_u, surfOseenSystem.S.Data, rho_u, surfOseenSystem.A_stab.Data);
                surfOseenSystem.C.Data *= -rho_p;
                surfOseenSystem.Schur.Data.LinComb(-1., surfOseenSystem.C.Data, 1., surfOseenSystem.M_p.Data);
                surfOseenSystem.Schur.Data *= 1. / (nu + gamma);
            logger.end();
            logger.beg("convert to Epetra");
                logger.beg("cast matrices");
                    auto A = static_cast<Epetra_CrsMatrix>(surfOseenSystem.sumA.Data);
                    auto printStat = [&](std::string const & name, Epetra_CrsMatrix const & A) {
                        logger.buf << name << ": " << A.NumGlobalRows() << 'x' << A.NumGlobalCols() << ", " << A.NumGlobalNonzeros() << " nonzeros (" << (100. * A.NumGlobalNonzeros()) / (A.NumGlobalRows() * A.NumGlobalCols()) << "%)";
                        logger.log();
                    };
                    printStat("A", A);
                    auto B = static_cast<Epetra_CrsMatrix>(surfOseenSystem.B.Data);
                    printStat("B", B);
                    auto C = static_cast<Epetra_CrsMatrix>(surfOseenSystem.C.Data);
                    printStat("C := -rho_p (pressure volume stabilization mtx)", C);
                    auto S = static_cast<Epetra_CrsMatrix>(surfOseenSystem.Schur.Data);
                    printStat("S := 1/(nu + gamma) (M_p - C)", S);
                logger.end();
                logger.beg("cast rhs vectors");
                    auto fRHS = static_cast<Epetra_Vector>(surfOseenSystem.fRHS.Data);
                    auto gRHS = static_cast<Epetra_Vector>(surfOseenSystem.gRHS.Data);
                    Epetra_Map mapVelocity(static_cast<int>(n), 0, comm), mapPressure(static_cast<int>(m), 0, comm), mapVelocityPressure(static_cast<int>(n + m), 0, comm);
                    Epetra_Vector belosLHS(mapVelocityPressure), belosRHS(mapVelocityPressure);
                    auto joinEpetraVectors = [&](Epetra_Vector const & a, Epetra_Vector const & b, MV& ab) {
                        for (size_t i = 0; i < n; ++i) (*ab(0))[i] = a[i];
                        for (size_t i = 0; i < m; ++i) (*ab(0))[i + n] = b[i];
                    };
                    joinEpetraVectors(fRHS, gRHS, belosRHS);
                logger.end();
                RCP<OP> belosMTX = rcp(new Epetra_OperatorApply([&](MV const & X, MV& Y) {
                    double* view;
                    X(0)->ExtractView(&view);
                    Epetra_Vector x1(Epetra_DataAccess::View, mapVelocity, view);
                    Epetra_Vector x2(Epetra_DataAccess::View, mapPressure, view + n);
                    Y(0)->ExtractView(&view);
                    Epetra_Vector y11(Epetra_DataAccess::View, mapVelocity, view);
                    auto y12 = y11;
                    Epetra_Vector y21(Epetra_DataAccess::View, mapPressure, view + n);
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
            logger.end();
            logger.beg("build preconditioners");
                auto identity  = [](MV const & X, MV& Y) {
                    double* view;
                    Y(0)->ExtractView(&view);
                    X(0)->ExtractCopy(view);
                };
                logger.beg("diffusion-convection-reaction block");
                    Epetra_OperatorApply::ApplyType invA = identity;
                    size_t numItersA = 0;
                    logger.log(inpJSON.get<std::string>("Solver.Inner.A.Comment"));
                    auto iterationA = inpJSON.get<std::string>("Solver.Inner.A.Iteration");
                    // amesos
                    Epetra_LinearProblem amesosProblem;
                    RCP<Amesos_BaseSolver> amesosSolver;
                    Amesos amesosFactory;
                    // belos
                    auto belosParamsA = parameterList();
                    decltype(belosSolver) belosSolverA;
                    // choose
                    if (iterationA.find("Amesos") != std::string::npos) { // Amesos
                        logger.log("using Amesos");
                        if (!amesosFactory.Query(iterationA))
                            throw std::invalid_argument("solver " + iterationA + " is not available");
                        logger.log("solver: " + iterationA);
                        amesosProblem.SetOperator(&A);
                        amesosSolver = rcp(amesosFactory.Create(iterationA, amesosProblem));
                        logger.beg("symbolic factorization");
                            AMESOS_CHK_ERR(amesosSolver->SymbolicFactorization());
                        logger.end();
                        logger.beg("numeric factorization");
                            AMESOS_CHK_ERR(amesosSolver->NumericFactorization());
                        logger.end();
                        invA = [&](MV const &X, MV &Y) {
                            amesosProblem.SetLHS(&Y);
                            amesosProblem.SetRHS(const_cast<MV*>(&X));
                            AMESOS_CHK_ERR(amesosSolver->Solve());
                            numItersA++;
                        };
                    } else { // Belos
                        logger.log("using Belos");
                        belosParamsA->set("Maximum Iterations", inpJSON.get<int>("Solver.Inner.A.MaxIter"));
                        belosParamsA->set("Convergence Tolerance", inpJSON.get<double>("Solver.Inner.A.RelResTol"));
                        belosSolverA = belosFactory.create(iterationA, belosParamsA);
                        logger.buf << "inner solver: " << belosSolverA->description();
                        logger.log();
                        invA = [&](MV const &X, MV &Y) {
                            Y.PutScalar(0.);
                            auto belosProblemA = LinearProblem<ST, MV, OP>(rcpFromRef(A), rcpFromRef(Y), rcpFromRef(X));
                            belosProblemA.setProblem();
                            belosSolverA->setProblem(rcpFromRef(belosProblemA));
                            belosSolverA->solve();
                            numItersA += belosSolverA->getNumIters();
                        };
                    }
                logger.end();
                logger.beg("schur complement block");
                    Epetra_OperatorApply::ApplyType invS = identity;
                    size_t numItersS = 0;
                    RCP<ParameterList> belosParamsS = parameterList();
                    belosParamsS->set("Maximum Iterations", inpJSON.get<int>("Solver.Inner.S.MaxIter"));
                    belosParamsS->set("Convergence Tolerance", inpJSON.get<double>("Solver.Inner.S.RelResTol"));
                    auto belosSolverS = belosFactory.create(inpJSON.get<std::string>("Solver.Inner.S.Iteration"), belosParamsS);
                    logger.buf << "inner solver: " << belosSolverS->description();
                    logger.log();
                    invS = [&](MV const & X, MV& Y) {
                        Y.PutScalar(0.);
                        auto belosProblemS = LinearProblem<ST, MV, OP>(rcpFromRef(S), rcpFromRef(Y), rcpFromRef(X));
                        belosProblemS.setProblem();
                        belosSolverS->setProblem(rcpFromRef(belosProblemS));
                        belosSolverS->solve();
                        numItersS += belosSolverS->getNumIters();
                    };
                logger.end();
                auto precType = inpJSON.get<std::string>("Solver.Inner.Type");
                RCP<OP> belosPRE;
                if (precType == "BlockDiagonal") {
                    logger.log("using block-diagonal preconditioner");
                    belosPRE = rcp(new Epetra_OperatorApply([&](MV const & X, MV& Y) {
                        double* view;
                        X(0)->ExtractView(&view);
                        Epetra_Vector x1(Epetra_DataAccess::View, mapVelocity, view);
                        Epetra_Vector x2(Epetra_DataAccess::View, mapPressure, view + n);
                        Y(0)->ExtractView(&view);
                        Epetra_Vector y1(Epetra_DataAccess::View, mapVelocity, view);
                        Epetra_Vector y2(Epetra_DataAccess::View, mapPressure, view + n);
                        #ifdef _OPENMP
                            #pragma omp parallel num_threads(2)
                            {
                                // y1
                                #pragma omp single
                                invA(x1, y1);
                                // y2
                                #pragma omp single
                                invS(x2, y2);
                            }
                        #else
                            // y1
                            invA(x1, y1);
                            // y2
                            invS(x2, y2);
                        #endif
                    }));
                }
                else {
                    logger.log("using block-triangular preconditioner");
                    belosPRE = rcp(new Epetra_OperatorApply([&](MV const & X, MV& Y) {
                        double* view;
                        X(0)->ExtractView(&view);
                        Epetra_Vector x1(Epetra_DataAccess::View, mapVelocity, view);
                        Epetra_Vector x2(Epetra_DataAccess::View, mapPressure, view + n);
                        Y(0)->ExtractView(&view);
                        Epetra_Vector y1(Epetra_DataAccess::View, mapVelocity, view);
                        auto y1_tmp = y1;
                        Epetra_Vector y2(Epetra_DataAccess::View, mapPressure, view + n);
                        // y2
                        invS(x2, y2);
                        // y1
                        B.Multiply(true, y2, y1_tmp);
                        y1_tmp.Update(1., x1, -1.);
                        invA(y1_tmp, y1);
                    }));
                }
            logger.end();
            logger.beg("linear solve");
            if (inpJSON.get<bool>("Solver.UsePreviousFrameAsInitialGuess")) {
                for (size_t i = 0; i < n; ++i) belosLHS[i] = u.Data[i];
                for (size_t i = 0; i < m; ++i) belosLHS[n + i] = p.Data[i];
            } else belosLHS.PutScalar(0.);
            auto belosProblem = LinearProblem<ST, MV, OP>(belosMTX, rcpFromRef(belosLHS), rcpFromRef(belosRHS));
                if (inpJSON.get<bool>("Solver.Inner.Use")) belosProblem.setRightPrec(belosPRE);
                belosProblem.setProblem();
                belosSolver->setProblem(rcpFromRef(belosProblem));
                std::cout << std::scientific;
                auto belosSolverResult = belosSolver->solve();
                if (myRank == 0) {
                    if (belosSolverResult == Belos::Converged)
                        logger.log("belos converged");
                    else
                        logger.wrn("belos did not converge");
                }
            auto solveTime = logger.end();
            logger.beg("convert from Epetra");
                for (size_t i = 0; i < n; ++i) u.Data[i] = belosLHS[i];
                for (size_t i = 0; i < m; ++i) p.Data[i] = belosLHS[n + i];
                p.Data -= dot(surfOseenSystem.M_p.Data * p.Data, I_p) / dot(surfOseenSystem.M_p.Data * I_p, I_p) * I_p;
            logger.end();
            logger.beg("output");
                MatrixCL B_T;
                auto residual = [&](VectorCL const & u, VectorCL const & p) {
                    auto velResSq = norm_sq(surfOseenSystem.sumA.Data * u + B_T * p - surfOseenSystem.fRHS.Data);
                    auto preResSq = norm_sq(surfOseenSystem.B.Data * u + surfOseenSystem.C.Data * p - surfOseenSystem.gRHS.Data);
                    return std::tuple<double, double, double>(sqrt(velResSq), sqrt(preResSq), sqrt(velResSq + preResSq));
                };
                auto exportStats = [&](size_t i) {
                    std::ofstream stats(dirName + "/stats/t_" + std::to_string(i) + ".json");
                    ParamCL tJSON;
                    transpose(surfOseenSystem.B.Data, B_T);
                    auto res = residual(u.Data, p.Data);
                    tJSON.put("Solver.Outer.ResidualNorm.TrueAbsolute.Velocity", std::get<0>(res));
                    tJSON.put("Solver.Outer.ResidualNorm.TrueAbsolute.Pressure", std::get<1>(res));
                    tJSON.put("Solver.Outer.ResidualNorm.TrueAbsolute.Full", std::get<2>(res));
                    tJSON.put("Solver.Outer.ResidualNorm.SolverRelative", belosSolver->achievedTol());
                    tJSON.put("Solver.Outer.TotalIters", belosSolver->getNumIters());
                    tJSON.put("Solver.Outer.DOF", n + m);
                    tJSON.put("Solver.Outer.Converged", belosSolverResult == Belos::Converged);
                    tJSON.put("ElapsedTime.LinearSolve", solveTime);
                    tJSON.put("Solver.Inner.A.TotalIters", numItersA);
                    tJSON.put("Solver.Inner.A.MeanIters", static_cast<double>(numItersA) / belosSolver->getNumIters());
                    tJSON.put("Solver.Inner.A.DOF", n);
                    tJSON.put("Solver.Inner.S.TotalIters", numItersS);
                    tJSON.put("Solver.Inner.S.MeanIters", static_cast<double>(numItersS) / belosSolver->getNumIters());
                    tJSON.put("Solver.Inner.S.DOF", m);
                    auto t = i * stepSize;
                    tJSON.put("Time", t);
                    tJSON.put("h", h);
                    tJSON.put("Peclet", Pe);
                    tJSON.put("MeshDepParams.rho_p", rho_p);
                    tJSON.put("MeshDepParams.rho_u", rho_u);
                    tJSON.put("MeshDepParams.tau_u", tau_u);
                    auto velL2 = sqrt(dot(u.Data, surfOseenSystem.M.Data * u.Data));
                    auto preL2 = sqrt(dot(p.Data, surfOseenSystem.M_p.Data * p.Data));
                    tJSON.put("Integral.PressureL2", preL2);
                    tJSON.put("Integral.VelocityL2", velL2);
                    tJSON.put("Integral.KineticEnergy", .5 * velL2 * velL2);
                    if (surfNavierStokesData.exactSoln) {
                        InitVector(mg, u_star, surfNavierStokesData.u_T, t);
                        InitScalar(mg, p_star, surfNavierStokesData.p, t);
                        res = residual(u_star.Data, p_star.Data);
                        tJSON.put("Solver.Outer.ResidualNorm.ExactSolnAbsolute.Velocity", std::get<0>(res));
                        tJSON.put("Solver.Outer.ResidualNorm.ExactSolnAbsolute.Pressure", std::get<1>(res));
                        tJSON.put("Solver.Outer.ResidualNorm.ExactSolnAbsolute.Full", std::get<2>(res));
                        VectorCL u_star_minus_u = u_star.Data - u.Data, p_star_minus_p = p_star.Data - p.Data;
                        auto velL2err = dot(u_star_minus_u, surfOseenSystem.M.Data * u_star_minus_u);
                        auto velNormalL2 = dot(u.Data, surfOseenSystem.S.Data * u.Data);
                        auto velTangenL2 = sqrt(velL2err - velNormalL2);
                        velL2err = sqrt(velL2err);
                        velNormalL2 = sqrt(velNormalL2);
                        auto velH1err = sqrt(dot(u_star_minus_u, surfOseenSystem.A.Data * u_star_minus_u));
                        auto preL2err = sqrt(dot(p_star_minus_p, surfOseenSystem.M_p.Data * p_star_minus_p));
                        tJSON.put("Integral.Error.VelocityL2", velL2err);
                        tJSON.put("Integral.Error.VelocityTangentialL2", velTangenL2);
                        tJSON.put("Integral.Error.VelocityNormalL2", velNormalL2);
                        tJSON.put("Integral.Error.VelocityH1", velH1err);
                        tJSON.put("Integral.Error.PressureL2", preL2err);
                    }
                    if (inpJSON.get<bool>("SurfNavStokes.ExportMatrices")) {
                        std::string format = inpJSON.get<std::string>("SurfNavStokes.ExportMatricesFormat") == "mtx" ? ".mtx" : ".mat";
                        auto expFunc = format == ".mtx" ? &MatrixCL::exportMTX : &MatrixCL::exportMAT;
                        auto expMat = [&](MatrixCL& A, std::string const a, std::string const & b) {
                            logger.beg(a);
                                (A.*expFunc)(dirName + "/matrices/" + b + format);
                                tJSON.put("Matrices." + a, "../matrices/" + b + format);
                            logger.end();
                        };
                        expMat(surfOseenSystem.sumA.Data, "DiffusionConvectionReaction", "A");
                        expMat(surfOseenSystem.B.Data, "Divergence", "B");
                        expMat(surfOseenSystem.C.Data, "PressureVolumeStab", "C");
                        expMat(surfOseenSystem.LB.Data, "VelocityScalarLaplaceBeltrami", "LB");
                        expMat(surfOseenSystem.LB_stab.Data, "VelocityScalarLaplaceBeltramiNormalStab", "LB_stab");
                        logger.log();
                    }
                    if (everyStep > 0 && (i-1) % everyStep == 0) {
                        logger.beg("write vtk");
                            writeVTK(t);
                        auto vtkTime = logger.end();
                        tJSON.put("ElapsedTime.VTK", vtkTime);
                    }
                    stats << tJSON;
                    logger.buf << tJSON;
                    logger.log();
                };
                exportStats(1);
            logger.end();
        logger.end();
        // no need to re-assemble these matrices:
        surfOseenSystem.A.assemble = false;
        surfOseenSystem.A_stab.assemble = false;
        surfOseenSystem.B.assemble = false;
        surfOseenSystem.M.assemble = false;
        surfOseenSystem.S.assemble = false;
        surfOseenSystem.M_p.assemble = false;
        surfOseenSystem.C.assemble = false;
        surfOseenSystem.LB.assemble = false;
        surfOseenSystem.LB_stab.assemble = false;
        for (size_t i = 2; i <= numbOfSteps; ++i) {
            numItersA = 0;
            numItersS = 0;
            std::stringstream header;
            header << "t = t_" << i << " (" << (100. * i) / numbOfSteps << "%)";
            logger.beg(header.str());
                logger.beg("assemble");
                    u_prev_prev = u_prev;
                    u_prev = u;
                    param.input.t = i * stepSize;
                    auto wind = u_prev.Data;
                    wind *= 2.;
                    wind -= u_prev_prev.Data;
                    setWind(wind);
                    Pe = supnorm(surfOseenSystem.w.Data) / nu;
                    logger.buf << "Pe = " << Pe;
                    logger.log();
                    SetupSurfOseen_P2P1(mg, lset, &surfOseenSystem, &param);
                    surfOseenSystem.gRHS.Data -= (dot(surfOseenSystem.gRHS.Data, I_p) / dot(I_p, I_p)) * I_p;
                    surfOseenSystem.fRHS.Data += (2. / stepSize) * (surfOseenSystem.M.Data * u_prev.Data) - (.5 / stepSize) * (surfOseenSystem.M.Data * u_prev_prev.Data) - gamma * surfOseenSystem.alRHS.Data;
                    surfOseenSystem.sumA.Data.LinComb(1.5 / stepSize, surfOseenSystem.M.Data, gamma, surfOseenSystem.AL.Data, 1., surfOseenSystem.N.Data, nu, surfOseenSystem.A.Data, tau_u, surfOseenSystem.S.Data, rho_u, surfOseenSystem.A_stab.Data);
                    surfOseenSystem.C.Data *= -rho_p;
                    surfOseenSystem.Schur.Data.LinComb(-1., surfOseenSystem.C.Data, 1., surfOseenSystem.M_p.Data);
                    surfOseenSystem.Schur.Data *= 1. / (nu + gamma);
                logger.end();
                logger.beg("convert to Epetra");
                    logger.beg("cast matrices");
                        A = static_cast<Epetra_CrsMatrix>(surfOseenSystem.sumA.Data);
                        printStat("A", A);
                        B = static_cast<Epetra_CrsMatrix>(surfOseenSystem.B.Data);
                        printStat("B", B);
                        C = static_cast<Epetra_CrsMatrix>(surfOseenSystem.C.Data);
                        printStat("C := -rho_p (pressure volume stabilization mtx)", C);
                        S = static_cast<Epetra_CrsMatrix>(surfOseenSystem.Schur.Data);
                        printStat("S := 1/(nu + gamma) (M_p - C)", S);
                    logger.end();
                    logger.beg("cast rhs vectors");
                        fRHS = static_cast<Epetra_Vector>(surfOseenSystem.fRHS.Data);
                        gRHS = static_cast<Epetra_Vector>(surfOseenSystem.gRHS.Data);
                        joinEpetraVectors(fRHS, gRHS, belosRHS);
                    logger.end();
                logger.end();
                if (iterationA.find("Amesos") != std::string::npos) {
                    logger.beg("symbolic factorization");
                        AMESOS_CHK_ERR(amesosSolver->SymbolicFactorization());
                    logger.end();
                    logger.beg("numeric factorization");
                        AMESOS_CHK_ERR(amesosSolver->NumericFactorization());
                    logger.end();
                }
                logger.beg("linear solve");
                    if (!inpJSON.get<bool>("Solver.UsePreviousFrameAsInitialGuess")) belosLHS.PutScalar(0.);
                    belosProblem = LinearProblem<ST, MV, OP>(belosMTX, rcpFromRef(belosLHS), rcpFromRef(belosRHS));
                    if (inpJSON.get<bool>("Solver.Inner.Use")) belosProblem.setRightPrec(belosPRE);
                    belosProblem.setProblem();
                    belosSolver->setProblem(rcpFromRef(belosProblem));
                    std::cout << std::scientific;
                    auto belosSolverResult = belosSolver->solve();
                    if (myRank == 0) {
                        if (belosSolverResult == Belos::Converged)
                            logger.log("belos converged");
                        else
                            logger.wrn("belos did not converge");
                    }
                solveTime = logger.end();
                logger.beg("convert from Epetra");
                    for (size_t i = 0; i < n; ++i) u.Data[i] = belosLHS[i];
                    for (size_t i = 0; i < m; ++i) p.Data[i] = belosLHS[n + i];
                    p.Data -= dot(surfOseenSystem.M_p.Data * p.Data, I_p) / dot(surfOseenSystem.M_p.Data * I_p, I_p) * I_p;
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