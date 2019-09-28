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
#include <surfactant/ifacetransp.h>

#include "SingletonLogger.hpp"

#include "surfactant/ifacetransp.h"
#include "misc/params.h"
#include "geom/builder.h"
#include "levelset/levelset.h"
#include "levelset/adaptriang.h"
#include "levelset/surfacetension.h"
#include "misc/dynamicload.h"
#include "out/vtkOut.h"
#include "num/bndData.h"
#include "num/precond.h"

#include "num/oseensolver.h"
#include "surfactant/ifacetransp.h"

#include "SurfNavierStokesData.hpp"

// #include "surfnavierstokes/surfnavierstokes_funcs.h"
#include "surfnavierstokes/surfnavierstokes_utils.h"
// #include "surfnavierstokes/surfnavierstokes_tests.h"

using namespace DROPS;

DROPS::Point3DCL shift;
DROPS::Point3DCL shiftTransform(const Point3DCL& p) {
    return p - shift;
}
bool fpEqual(double a, double b) {
    return std::fabs(a - b) < 1e-10;
}

int main(int argc, char* argv[]) {
    auto& logger = SingletonLogger::instance();
    try {
        logger.beg("read input .json");
            ParamCL P;
            read_parameter_file_from_cmdline(P, argc, argv, "../../param/surfnavierstokes/No_Bnd_Condition.json");
            dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"), P.get<std::vector<std::string>>("General.DynamicLibs"));
            logger.buf << P;
            logger.log();
        logger.end();
        logger.beg("set up test case");
            // parse FE types and some other parameters from json file
            auto FE = P.get<std::string>("SurfNavStokes.FE");
            auto velFE = FE.substr(0,2);
            auto preFE = FE.substr(2,2);
            logger.buf
                << "velocity FE = " << velFE << '\n'
                << "pressure FE = " << preFE;
            logger.log();
            auto testName = P.get<std::string>("SurfNavStokes.TestName");
            auto rho = P.get<double>("SurfNavStokes.rho");
            auto mu = P.get<double>("SurfNavStokes.mu");
            auto surfNavierStokesData = SurfNavierStokesDataFactory(testName, rho, mu);
            logger.wrn("TODO: set mu and rho in the solution");
            logger.buf
                << "rho = " << rho << '\n'
                << "mu  = " << mu;
            logger.log();
            LocalStokesParam param;
            param.input.f = P.get<bool>("SurfNavStokes.IntegrateRHS") ? surfNavierStokesData.f_T : nullptr;
            param.input.g = P.get<bool>("SurfNavStokes.IntegrateRHS") ? surfNavierStokesData.m_g : nullptr;
            param.input.exactNormal = P.get<bool>("SurfNavStokes.ComputeNormalErr") ? surfNavierStokesData.surface.n : nullptr;
            param.input.exactShape = P.get<bool>("SurfNavStokes.ComputeShapeErr")  ? surfNavierStokesData.surface.H : nullptr;
            param.input.computeMatrices = P.get<bool>("SurfNavStokes.ComputeMatrices");
            param.input.numbOfVirtualSubEdges = P.get<size_t>("SurfNavStokes.NumbOfVirtualSubEdges");
            param.input.usePatchNormal = P.get<bool>("SurfNavStokes.UsePatchNormals");
            if (param.input.usePatchNormal)  logger.log("using patch normals in surf integrands");
            auto formulation = P.get<std::string>("SurfNavStokes.formulation");
            if (formulation == "consistent") {
                param.input.formulation = LocalStokesParam::Formulation::consistent;
                logger.log("using consistent penalty formulation");
            }
            else {
                param.input.formulation = LocalStokesParam::Formulation::inconsistent;
                logger.log("using inconsistent penalty formulation");
            }
            StokesSystem stokesSystem;
            if (!P.get<bool>("SurfNavStokes.ExportMatrices")) {
                stokesSystem.Schur_full_stab.assemble = P.get<std::string>("SurfNavStokes.stab") == "full";
                stokesSystem.Schur_normal_stab.assemble = P.get<std::string>("SurfNavStokes.stab") == "normal";
                stokesSystem.LB.assemble = false;
                stokesSystem.LB_stab.assemble = false;
            }
            auto h = P.get<DROPS::Point3DCL>("Mesh.E1")[0] / P.get<double>("Mesh.N1") * std::pow(2., -P.get<double>("Mesh.AdaptRef.FinestLevel"));
            auto tau_u_order  = P.get<double>("SurfNavStokes.normal_penalty_pow");
            auto tau_u_factor = P.get<double>("SurfNavStokes.normal_penalty_fac");
            auto tau_u 		  = tau_u_factor * pow(h, tau_u_order); // constant for normal penalty
            auto rho_u_order  = P.get<double>("SurfNavStokes.vel_volumestab_pow");
            auto rho_u_factor = P.get<double>("SurfNavStokes.vel_volumestab_fac");
            auto rho_u        = rho_u_factor  * pow(h, rho_u_order); // constant for velocity stabilisation
            auto rho_p_order  = P.get<double>("SurfNavStokes.pre_volumestab_pow");
            auto rho_p_factor = P.get<double>("SurfNavStokes.pre_volumestab_fac");
            auto rho_p        = rho_p_factor  * pow(h, rho_p_order); // constant for pressure stabilisation
            auto numbOfSteps  = P.get<size_t>("Time.NumbOfSteps");
            auto finalTime    = P.get<double>("Time.FinalTime");
            auto stepSize     = finalTime / numbOfSteps;
            logger.buf << "$\\Delta t$ = " << stepSize;
            logger.log();
            auto dirName = P.get<std::string>("Output.Directory");
            system(("mkdir -p " + dirName).c_str());
            std::ofstream output(dirName + "/output.json");
        logger.end();
        logger.beg("build mesh");
            logger.beg("build initial bulk mesh");
                std::auto_ptr<DROPS::MGBuilderCL> builder(DROPS::make_MGBuilder(P));
                DROPS::MultiGridCL mg(*builder);
                const DROPS::ParamCL::ptree_type* ch= 0;
                try {
                    ch = &P.get_child("Mesh.Periodicity");
                } catch (DROPS::DROPSParamErrCL) {}
                if (ch) read_PeriodicBoundaries(mg, *ch);
                // levelset shift
                shift = P.get<DROPS::Point3DCL>("Levelset.ShiftDir", DROPS::Point3DCL(0., 0., 0.));
                shift /= shift.norm();
                auto shiftNorm = fabs(P.get<double>("Levelset.ShiftNorm", 0.));
                shift *= shiftNorm;
                mg.Transform(shiftTransform);
                logger.buf << "surface shift = " << shift;
                logger.log();
            logger.end();
            logger.beg("refine towards the surface");
                DROPS::AdapTriangCL adap(mg);
                // adaptive mesh refinement based on level set function
                DROPS::DistMarkingStrategyCL initmarker(surfNavierStokesData.surface.phi, P.get<double>("Mesh.AdaptRef.Width"), P.get<int>("Mesh.AdaptRef.CoarsestLevel"), P.get<int>("Mesh.AdaptRef.FinestLevel"));
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
            instat_scalar_fun_ptr sigma(0);
            SurfaceTensionCL sf(sigma, 0);
            BndDataCL<double> lsbnd(0);
            read_BndData(lsbnd, mg, P.get_child("Levelset.BndData"));
            DROPS::LevelsetP2CL & lset(*DROPS::LevelsetP2CL::Create(mg, lsbnd, sf));
            lset.CreateNumbering(mg.GetLastLevel(), &lset.idx);
            lset.Phi.SetIdx(&lset.idx);
            lset.Init(surfNavierStokesData.surface.phi);
        logger.end();
        logger.beg("build FE spaces");
            BndDataCL<Point3DCL> vbnd(0);
            read_BndData(vbnd, mg, P.get_child("Stokes.VelocityBndData"));
            BndDataCL<double> pbnd(0);
            read_BndData(pbnd, mg, P.get_child("Stokes.PressureBndData"));
            DROPS::IdxDescCL ifaceVecP2idx(vecP2IF_FE, vbnd);
            DROPS::IdxDescCL ifaceVecP1idx(vecP1IF_FE, vbnd);
            DROPS::IdxDescCL ifaceP1idx(P1IF_FE, pbnd);
            DROPS::IdxDescCL ifaceP2idx(P2IF_FE, pbnd);
            DROPS::IdxDescCL vecP2idx(vecP2_FE, vbnd);
            DROPS::IdxDescCL vecP1idx(vecP1_FE, vbnd);
            DROPS::IdxDescCL P1FEidx(P1_FE, pbnd);
            DROPS::IdxDescCL P2FEidx(P2_FE, pbnd);
            ifaceVecP2idx.GetXidx().SetBound(P.get<double>("SurfTransp.OmitBound"));
            auto numbOfActiveTetrasP2 = ifaceVecP2idx.CreateNumbering(mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
            ifaceVecP1idx.GetXidx().SetBound(P.get<double>("SurfTransp.OmitBound"));
            auto numbOfActiveTetrasP1 = ifaceVecP1idx.CreateNumbering(mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
            if (numbOfActiveTetrasP2 != numbOfActiveTetrasP1) {
                std::stringstream err;
                err << "inconsistent numb of cut tetras for P2 and P1: " << numbOfActiveTetrasP2 << " vs. " << numbOfActiveTetrasP1;
                throw std::logic_error(err.str());
            }
            logger.buf << "numb of active (cut) tetras is: " << numbOfActiveTetrasP1 << " (" << (100. * numbOfActiveTetrasP1) / numbOfTetras << "%)\n";
            logger.log();
            ifaceP1idx.GetXidx().SetBound(P.get<double>("SurfTransp.OmitBound"));
            ifaceP1idx.CreateNumbering(mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
            ifaceP2idx.GetXidx().SetBound(P.get<double>("SurfTransp.OmitBound"));
            ifaceP2idx.CreateNumbering(mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
            vecP2idx.GetXidx().SetBound(P.get<double>("SurfTransp.OmitBound"));
            vecP2idx.CreateNumbering(mg.GetLastLevel(), mg);
            vecP1idx.GetXidx().SetBound(P.get<double>("SurfTransp.OmitBound"));
            vecP1idx.CreateNumbering(mg.GetLastLevel(), mg);
            P1FEidx.GetXidx().SetBound(P.get<double>("SurfTransp.OmitBound"));
            P1FEidx.CreateNumbering(mg.GetLastLevel(), mg);
            P2FEidx.GetXidx().SetBound(P.get<double>("SurfTransp.OmitBound"));
            P2FEidx.CreateNumbering(mg.GetLastLevel(), mg);
            logger.buf
                    << "numb of d.o.f. vector IFP2: " << ifaceVecP2idx.NumUnknowns() << '\n'
                    << "numb of d.o.f. vector IFP1: " << ifaceVecP1idx.NumUnknowns() << '\n'
                    << "numb of d.o.f. scalar IFP2: " << ifaceP2idx.NumUnknowns() << '\n'
                    << "numb of d.o.f. scalar IFP1: " << ifaceP1idx.NumUnknowns() << '\n'
                    << "numb of d.o.f. vector P2: " << vecP2idx.NumUnknowns() << '\n'
                    << "numb of d.o.f. vector P1: " << vecP1idx.NumUnknowns() << '\n'
                    << "numb of d.o.f. scalar P2: " << P2FEidx.NumUnknowns() << '\n'
                    << "numb of d.o.f. scalar P1: " << P1FEidx.NumUnknowns();
            logger.log();
            DROPS::VecDescCL u_0, u, u_prev, u_prev_prev, p_0, p;
            if (FE == "P2P1") {
                u.SetIdx(&ifaceVecP2idx);
                u_0.SetIdx(&ifaceVecP2idx);
                u_prev.SetIdx(&ifaceVecP2idx);
                u_prev_prev.SetIdx(&ifaceVecP2idx);
                stokesSystem.fRHS.SetIdx(&ifaceVecP2idx);
                stokesSystem.A.SetIdx(&ifaceVecP2idx, &ifaceVecP2idx);
                stokesSystem.N.SetIdx(&ifaceVecP2idx, &ifaceVecP2idx);
                stokesSystem.A_stab.SetIdx(&ifaceVecP2idx, &ifaceVecP2idx);
                stokesSystem.B.SetIdx(&ifaceP1idx, &ifaceVecP2idx);
                stokesSystem.M.SetIdx(&ifaceVecP2idx, &ifaceVecP2idx);
                stokesSystem.S.SetIdx(&ifaceVecP2idx, &ifaceVecP2idx);
            }
            else if(FE == "P1P1") {
                u.SetIdx(&ifaceVecP1idx);
                u_0.SetIdx(&ifaceVecP1idx);
                u_prev.SetIdx(&ifaceVecP1idx);
                u_prev_prev.SetIdx(&ifaceVecP1idx);
                stokesSystem.fRHS.SetIdx(&ifaceVecP1idx);
                stokesSystem.gRHS.SetIdx(&ifaceP1idx);
                stokesSystem.A.SetIdx(&ifaceVecP1idx, &ifaceVecP1idx);
                stokesSystem.N.SetIdx(&ifaceVecP1idx, &ifaceVecP1idx);
                stokesSystem.A_stab.SetIdx(&ifaceVecP1idx, &ifaceVecP1idx);
                stokesSystem.B.SetIdx(&ifaceP1idx, &ifaceVecP1idx);
                stokesSystem.M.SetIdx(&ifaceVecP1idx, &ifaceVecP1idx);
                stokesSystem.S.SetIdx(&ifaceVecP1idx, &ifaceVecP1idx);
            }
            else throw std::invalid_argument("unknown FE pair");
            p.SetIdx(&ifaceP1idx);
            p_0.SetIdx(&ifaceP1idx);
            stokesSystem.gRHS.SetIdx(&ifaceP1idx);
            stokesSystem.Schur.SetIdx(&ifaceP1idx, &ifaceP1idx);
            stokesSystem.Schur_full_stab.SetIdx(&ifaceP1idx, &ifaceP1idx);
            stokesSystem.Schur_normal_stab.SetIdx(&ifaceP1idx, &ifaceP1idx);
        logger.end();
        logger.beg("interpolate initial condition");
            InitVector(mg, u_prev, surfNavierStokesData.u_T);
            // InitScalar(mg, p_0, surfNavierStokesData.p);
        logger.end();
        logger.beg("assemble");
            MatrixCL Schur_hat;
            if (FE == "P2P1")
                SetupStokesIF_P2P1(mg, lset, &stokesSystem, &param);
            else if (FE == "P1P1") {
                // ...
                throw std::invalid_argument("TODO: refactor P1P1 assembly");
            }
            logger.buf
                << "numb of cut tetras is         " << param.output.numbOfCutTetras << '\n'
                << "stiffness mtx is              " << stokesSystem.A.Data.num_rows() << " * " << stokesSystem.A.Data.num_cols() << '\n'
                << "velocity soln size is         " << u.Data.size() << '\n'
                << "exact velocity interp size is " << u_0.Data.size() << '\n'
                << "pre mass mtx is               " << stokesSystem.Schur.Data.num_rows() << " * " << stokesSystem.Schur.Data.num_cols() << '\n'
                << "pressure soln size is         " << p.Data.size() << '\n'
                << "exact pressure interp size is " << p_0.Data.size() << '\n'
                << "f size is                     " << stokesSystem.fRHS.Data.size() << '\n'
                << "g size is                     " << stokesSystem.gRHS.Data.size();
        logger.end();
        if (param.input.computeMatrices && P.get<bool>("SurfNavStokes.ExportMatrices")) {
            logger.beg("export matrices to " + dirName);
                MatrixCL A_final;
                A_final.LinComb(1., stokesSystem.A.Data, 1., stokesSystem.M.Data, tau_u, stokesSystem.S.Data, rho_u, stokesSystem.A_stab.Data);
                auto C_full = stokesSystem.Schur_full_stab.Data;
                C_full *= rho_p;
                auto C_n = stokesSystem.Schur_normal_stab.Data;
                C_n *= rho_p;
                auto format = P.get<std::string>("SurfNavStokes.ExportMatricesFormat") == "mtx" ? ".mtx" : ".mat";
                auto expFunc = format == ".mtx" ? &MatrixCL::exportMTX : &MatrixCL::exportMAT;
                (A_final.*expFunc)(dirName + "/A" + format);
                (stokesSystem.B.Data.*expFunc)(dirName + "/B" + format);
                (stokesSystem.Schur.Data.*expFunc)(dirName + "/M" + format);
                (C_full.*expFunc)(dirName + "/C_full" + format);
                (C_n.*expFunc)(dirName + "/C_n" + format);
                (stokesSystem.LB.Data.*expFunc)(dirName + "/LB" + format);
                (stokesSystem.LB_stab.Data.*expFunc)(dirName + "/LB_stab" + format);
            logger.end();
        }
        if (P.get<bool>("SurfNavStokes.ComputeNormalErr")) {
            P.put("Result.LevelsetNormalErr", sqrt(param.output.normalErrSq.lvset));
            P.put("Result.PatchNormalErr", sqrt(param.output.normalErrSq.patch));
        }
        if (P.get<bool>("SurfNavStokes.ComputeShapeErr"))
            P.put("Result.ShapeErr", sqrt(param.output.shapeErrSq));
        if (!P.get<bool>("SurfNavStokes.Solve")) {
            output << P;
            logger.buf << P;
            logger.log();
            return 0;
        }
        logger.beg("build preconditioners");
            if (P.get<std::string>("SurfNavStokes.stab") == "full")
                Schur_hat.LinComb(1., stokesSystem.Schur.Data, rho_p, stokesSystem.Schur_full_stab.Data);
            else if (P.get<std::string>("SurfNavStokes.stab") == "normal")
                Schur_hat.LinComb(1., stokesSystem.Schur.Data, rho_p, stokesSystem.Schur_normal_stab.Data);
            else {
                logger.wrn("no stabilization is used! Schur precond = I");
                Schur_hat = std::valarray<double>(1., stokesSystem.Schur.Data.num_rows());
            }
            typedef SSORPcCL SymmPcPcT;
            typedef PCGSolverCL<SymmPcPcT> PCGSolverT;
            typedef SolverAsPreCL<PCGSolverT> PCGPcT;
            SymmPcPcT symmPcPc_;
            std::stringstream Astream;
            PCGSolverT PCGSolver_(symmPcPc_, P.get<int>("Solver.PcAIter"), P.get<double>("Solver.PcATol"), true, &Astream);
            PCGPcT PCGPc_(PCGSolver_);//velocity preconditioner
            std::stringstream Schurstream;
            PCGSolverT SchurPCGSolver (symmPcPc_, P.get<int>("Solver.PcBIter"), P.get<double>("Solver.PcBTol"), true, &Schurstream);
            SchurPreBaseCL *spc_ = new SurfaceLaplacePreCL<PCGSolverT>(Schur_hat, SchurPCGSolver);//pressure preconditioner
            //construct symmteric block iterative solver
            typedef BlockPreCL<ExpensivePreBaseCL, SchurPreBaseCL, DiagSpdBlockPreCL>  DiagBlockPcT;
            typedef PLanczosONBCL<VectorCL, DiagBlockPcT> LanczosT;
            typedef PMResSolverCL<LanczosT> MinResT;
            DiagBlockPcT    *DBlock_ = new DiagBlockPcT(PCGPc_, *spc_);
            LanczosT *lanczos_ = new LanczosT(*DBlock_);
            MinResT *MinRes_ = new MinResT(*lanczos_,  P.get<int>("Solver.Iter"), P.get<double>("Solver.Tol"), /*relative*/ false);
            BlockMatrixSolverCL<MinResT> *symStokesSolver= new BlockMatrixSolverCL<MinResT>(*MinRes_);
            // construct preconditioners for nonsymmetric case
            typedef /*JACPcCL*/ SSORPcCL NonSymmPcPcT;
            typedef GMResSolverCL<NonSymmPcPcT> GMResSolverT;
            typedef SolverAsPreCL<GMResSolverT> GMResPcT;
            NonSymmPcPcT nonsymmPcPc_;
            std::stringstream PCstream;
            GMResSolverT GMResSolver_(nonsymmPcPc_, P.get<int>("Solver.PcAIter"), P.get<int>("Solver.PcAIter"), P.get<double>("Solver.PcATol"),
                     /*bool relative=*/ true, /*bool calculate2norm=*/ false, /*PreMethGMRES method=*/ LeftPreconditioning,
                         /*bool mod =*/ true,      /*bool useModGS =*/ false, &PCstream);
            GMResPcT GMResPc_nonsym(GMResSolver_);
            //setup iterative block solver for nonsymmetric case
            typedef GMResSolverCL<DiagBlockPcT> GMResT;
            DiagBlockPcT    *DBlock_nonsym = new DiagBlockPcT(GMResPc_nonsym, *spc_);
            GMResT *GMRes_ = new GMResT(*DBlock_nonsym, P.get<int>("Solver.Iter"), P.get<int>("Solver.Iter"), P.get<double>("Solver.Tol"),
                     /*bool relative=*/ false, /*bool calculate2norm=*/ false, /*PreMethGMRES method=*/ LeftPreconditioning,
                          /*bool mod =*/ true,      /*bool useModGS =*/ false);
            BlockMatrixSolverCL<GMResT> *nonsymStokesSolver = new BlockMatrixSolverCL<GMResT>(*GMRes_);
            StokesSolverBaseCL* stokesSolver;
        logger.end();
        logger.beg("solve");
            //set up discrete vectors and matrices
            VectorCL id, id2 ;
            id.resize(stokesSystem.fRHS.Data.size(), 1.);
            id2.resize(stokesSystem.gRHS.Data.size(),1.);
            MatrixCL A_dyn, C, B_T;
//            MatrixCL A_dyn, Ahat, B, Chat, Mhat;
//            MatrixCL B_T, NTranspose;
//            DROPS::VecDescCL instantrhs;
//            DROPS::VecDescCL vxtent, pxtent;
//            //set up output
//            std::string filename = "test" + testcase + "_";
//            std::cout << "dirName: " << dirName << '\n';
//            //set up VTK output
            VTKOutCL* vtkwriter = nullptr;
            vtkwriter = new VTKOutCL(mg, "DROPS data", (int)P.get<double>("Time.NumSteps")/P.get<int>("Output.EveryStep"), dirName , testName + "_", testName, 0);
            vtkwriter->Register(make_VTKScalar(lset.GetSolution(), "level-set"));
            vtkwriter->Register(make_VTKIfaceVector(mg, u_0, "u_*", velFE, vbnd));
            vtkwriter->Register(make_VTKIfaceVector(mg, stokesSystem.fRHS, "f_T", velFE, vbnd));
            vtkwriter->Register(make_VTKIfaceVector(mg, u, "u_h", velFE, vbnd));
            vtkwriter->Register(make_VTKIfaceScalar(mg, p, "p_h", /*preFE,*/ pbnd));
            vtkwriter->Register(make_VTKIfaceScalar(mg, stokesSystem.gRHS, "-g", /*preFE,*/ pbnd));
            vtkwriter->Register(make_VTKIfaceScalar(mg, p_0, "p_*", /*preFE,*/ pbnd));

            std::unordered_map<std::string, std::vector<double>> stat;
//
//            //set up a txt file for custom time output
//            std::ofstream log_solo(dirName +"/"
//                                  + "tau_u="+ std::to_string(float(tau_u))+ "_"
//                                  + "l="  + P.get<std::string>("Mesh.AdaptRef.FinestLevel")
//                                  + "_nu=" + P.get<std::string>("SurfNavStokes.kinematic_viscosity")
//                                  + "_Plot_" + filename
//                                  + P.get<std::string>("SurfNavStokes.nonlinear_term") + "_"
//                                  + P.get<std::string>("SurfNavStokes.instationary") + "="
//                                  + std::to_string(float(tau)) + ".txt");
//            log_solo << "Time" <<  "\tKinetic\t" << "\tMomentum\t" << "\tWork" <<  '\n';
//
//            //set up a txt file for error time output
//            std::ofstream log_error(dirName +"/"
//                    + "tau_u="+ std::to_string(float(tau_u))+ "_"
//                    + "l="  + P.get<std::string>("Mesh.AdaptRef.FinestLevel")
//                    + "_nu=" + P.get<std::string>("SurfNavStokes.kinematic_viscosity")
//                    + "_Error_" + filename
//                    + P.get<std::string>("SurfNavStokes.nonlinear_term") + "_"
//                    + P.get<std::string>("SurfNavStokes.instationary") + "="
//                    + std::to_string(float(tau)) + ".txt");
//            log_error << "Time" << "\tL_2(uT-ext_uT)\t" <<  "\tadvH_1(u-ext_u)\t" <<  "\tsurfH_1(u-ext_u)\t" <<  "\tH_1(u-ext_u)\t" << "\tL_2(uN)\t" << "\tL_2(p-ext_p)" << '\n';
//            double mu = P.get<double>("SurfNavStokes.kinematic_viscosity");
//            std::cout << "viscosity: " << mu << '\n';
            transpose(stokesSystem.B.Data, B_T);
            if (P.get<std::string>("SurfNavStokes.stab") == "full")
                C.LinComb(0, stokesSystem.Schur.Data, -rho_p, stokesSystem.Schur_full_stab.Data);
            else if (P.get<std::string>("SurfNavStokes.stab") == "normal")
                C.LinComb(0, stokesSystem.Schur.Data, -rho_p, stokesSystem.Schur_normal_stab.Data);
            else
                C.LinComb(0, stokesSystem.Schur.Data, 0, stokesSystem.Schur.Data);

            logger.beg("t = t_1");
                logger.beg("linear solve");
                    stokesSystem.fRHS.Data += (1. / stepSize) * u_prev.Data;
                    stokesSystem.gRHS.Data -= (dot(stokesSystem.gRHS.Data, id2) / dot(id2, id2)) * id2; // renormalize
                    auto wind = rho * u_prev.Data;
                    auto windRises = fpEqual(norm(wind), 0.);
                    if (windRises) logger.log("using unsym solver");
                    else logger.log("using sym solver");
                    stokesSolver = windRises ? nonsymStokesSolver : symStokesSolver;
                    A_dyn.LinComb(mu, stokesSystem.A.Data, rho / stepSize, stokesSystem.M.Data, tau_u, stokesSystem.S.Data, rho_u, stokesSystem.A_stab.Data);
                    stokesSolver->Solve(A_dyn, stokesSystem.B.Data, C, u.Data, p.Data, stokesSystem.fRHS.Data, stokesSystem.gRHS.Data, u.RowIdx->GetEx(), p.RowIdx->GetEx());
                logger.end();
                if (P.get<int>("Output.EveryStep") > 0) {
                    logger.beg("write vtk");
                        vtkwriter->Write(1);
                    logger.end();
                }
                logger.beg("export stat");
                    auto velResSq = norm_sq(A_dyn * u.Data + B_T * p.Data - stokesSystem.fRHS.Data);
                    auto preResSq = norm_sq(stokesSystem.B.Data * u.Data + C * p.Data - stokesSystem.gRHS.Data);
                    auto residual = sqrt(velResSq + preResSq);
                    std::cout << "residual: " << residual << '\n';
                    p.Data -=  dot(stokesSystem.Schur.Data * p.Data, id2) / dot(stokesSystem.Schur.Data*id2, id2) * id2;
                    std::stringstream logss;
                    logss << "norm(A_dyn * v + B^T  * p -  rhs): " << sqrt(velResSq) << '\n'
                        << "norm(B    * v + C * p - stokesSystem.gRHS): " << sqrt(preResSq) << '\n';
                    // send to logfiles

                    logss << "Time is: " << std::to_string((float)1) << '\n';
                    logss << "h is: " << h << '\n';
                    logss << "rho_p is: "   << rho_p << '\n';
                    logss << "rho_u is: " << rho_u << '\n';
                    logss << "tau_u is: "     << tau_u << '\n';
                    logss << "Total iterations: " << Solver->GetIter() << '\n';
                    logss	<< "Final MINRES residual: " << Solver->GetResid() << '\n';
                    VectorCL u_exactMinusV = u_exact.Data - v.Data, p_exactMinusP = p_exact.Data - p.Data;
                    auto velL2          = sqrt(dot(v.Data, M.Data * v.Data));
                    auto velL2err       = dot(u_exactMinusV, M.Data * u_exactMinusV);
                    auto velNormalL2    = dot(v.Data, S.Data * v.Data);
                    auto velTangenL2    = sqrt(velL2err - velNormalL2);
                    velL2err = sqrt(velL2err);
                    velNormalL2 = sqrt(velNormalL2);
                    auto velH1err       = sqrt(dot(u_exactMinusV, A.Data * u_exactMinusV));
                    auto preL2          = sqrt(dot(p.Data, Schur.Data * p.Data));
                    auto preL2err       = sqrt(dot(p_exactMinusP, Schur.Data * p_exactMinusP));
                    logss << "The L2-Norm of v - u_exact is: " << velL2err << '\n';
                    logss << "The L2-Norm of v  is: " << velL2 << '\n';


                    log_solo << std::to_string((float)(velL2*velL2*0.5)) << '\n';

                    logss << "The H1-Norm of v - u_exact is: " << velH1err << '\n';
                    logss << "The L2-Norm of v * n is: " << velNormalL2 << '\n';
                    logss << "The L2-Norm of p - p_exact is: " << preL2err << '\n';
                    logss << "The L2-Norm of p is: " << preL2 << '\n';
                    logss << "The L2-Norm of v_T - u_exact is: " << velTangenL2 << '\n';
                    logss << "Actual residual is: " << residual << '\n';
                logger.end();
            logger.end();
            logger.beg("t = t_2, ..., t_n");
                for (size_t i = 2; i <= numbOfSteps; ++i) {
                    logger.wrn("TODO!");
                    logger.pro(i, numbOfSteps);
                }
            logger.end();
        logger.end();
//        // error
//        double Aaverage(0.), Avariation(0.), Schuraverage(0.), Schurvariation(0.);
//        double PCaverage(0.), PCvariation(0.);
//        std::stringstream Astreamcopy(Astream.str());
//        std::stringstream PCstreamcopy(PCstream.str());
//        std::stringstream Schurstreamcopy(Schurstream.str());
//        RightComputeAverageIterations(Astream, Aaverage);
//        ComputeAverageIterations(Schurstream, Schuraverage);
//        RightComputeAverageIterations(PCstream, PCaverage);
//        RightComputeVariationFromAverageIterations(Astreamcopy, Aaverage, Avariation);
//        RightComputeVariationFromAverageIterations(PCstreamcopy, PCaverage, PCvariation);
//        ComputeVariationFromAverageIterations(Schurstreamcopy, Schuraverage, Schurvariation);
//        std::cout << "The average iterationsnumber of the A-preconditioner is: " <<     Aaverage     << '\n' << " ...with a variation of: " << Avariation << '\n';
//        std::cout << "The average iterationsnumber of the nonsymmetric A-preconditioner is: " <<     PCaverage     << '\n' << " ...with a variation of: " << PCvariation << '\n';
//        std::cout << "The average iterationsnumber of the Schur-preconditioner is: " << Schuraverage << '\n' << " ...with a variation of: " << Schurvariation << '\n';
//        delete &lset;
        return 0;
    } catch (std::exception const & e) {
        logger.err(e.what());
        return 1;
    }
}