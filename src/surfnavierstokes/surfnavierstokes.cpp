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
#include "surfnavierstokes/surfnavierstokes_utils.h"

using namespace DROPS;

DROPS::Point3DCL shift;
DROPS::Point3DCL shiftTransform(const Point3DCL& p) {
    return p - shift;
}

bool windRises(double nu, VecDescCL const & w, double& Pe) {
    Pe = supnorm(w.Data) / nu;
    return Pe > .1; // if the wind is more than 10% of the diffusion
}

int main(int argc, char* argv[]) {
    auto& logger = SingletonLogger::instance();
    try {
        logger.beg("read input .json");
            ParamCL inpJSON, outJSON;
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
            std::ofstream output(dirName + "/output.json");
            {
                std::ofstream input(dirName + "/input.json");
                input << inpJSON;
                logger.buf << inpJSON;
                logger.log();
            }
        logger.end();
        logger.beg("set up test case");
            auto h = inpJSON.get<DROPS::Point3DCL>("Mesh.E1")[0] / inpJSON.get<double>("Mesh.N1") * std::pow(2., -inpJSON.get<double>("Mesh.AdaptRef.FinestLevel"));
            auto tau_u_order  = inpJSON.get<double>("SurfNavStokes.normal_penalty_pow");
            auto tau_u_factor = inpJSON.get<double>("SurfNavStokes.normal_penalty_fac");
            auto tau_u 		  = tau_u_factor * pow(h, tau_u_order); // constant for normal penalty
            auto rho_u_order  = inpJSON.get<double>("SurfNavStokes.vel_volumestab_pow");
            auto rho_u_factor = inpJSON.get<double>("SurfNavStokes.vel_volumestab_fac");
            auto rho_u        = rho_u_factor  * pow(h, rho_u_order); // constant for velocity stabilisation
            auto rho_p_order  = inpJSON.get<double>("SurfNavStokes.pre_volumestab_pow");
            auto rho_p_factor = inpJSON.get<double>("SurfNavStokes.pre_volumestab_fac");
            auto rho_p        = rho_p_factor  * pow(h, rho_p_order); // constant for pressure stabilisation
            auto numbOfSteps  = inpJSON.get<size_t>("Time.NumbOfSteps");
            auto finalTime    = inpJSON.get<double>("Time.FinalTime");
            auto stepSize     = finalTime / numbOfSteps;
            auto everyStep = inpJSON.get<int>("Output.EveryStep");
            int numbOfStepsVTK = ceil(static_cast<double>(everyStep) / everyStep);
            if (numbOfStepsVTK < 0) numbOfStepsVTK = 0;
            logger.buf << "$\\Delta t$ = " << stepSize << '\n';
            // parse FE types and some other parameters from json file
            auto FE = inpJSON.get<std::string>("SurfNavStokes.FE");
            auto velFE = FE.substr(0,2);
            auto preFE = FE.substr(2,2);
            logger.buf
                << "velocity FE = " << velFE << '\n'
                << "pressure FE = " << preFE << '\n';
            auto testName = inpJSON.get<std::string>("SurfNavStokes.TestName");
            auto nu = inpJSON.get<double>("SurfNavStokes.nu");
            auto surfNavierStokesData = SurfNavierStokesDataFactory(testName, nu);
            logger.buf << "$\\nu$ = " << nu << '\n';
            LocalStokesParam param;
            param.input.t = stepSize;
            param.input.f = inpJSON.get<bool>("SurfNavStokes.IntegrateRHS") ? surfNavierStokesData.f_T : nullptr;
            param.input.g = inpJSON.get<bool>("SurfNavStokes.IntegrateRHS") ? surfNavierStokesData.m_g : nullptr;
            param.input.exactNormal = inpJSON.get<bool>("SurfNavStokes.ComputeNormalErr") ? surfNavierStokesData.surface.n : nullptr;
            param.input.exactShape = inpJSON.get<bool>("SurfNavStokes.ComputeShapeErr") ? surfNavierStokesData.surface.H : nullptr;
            param.input.numbOfVirtualSubEdges = inpJSON.get<size_t>("SurfNavStokes.NumbOfVirtualSubEdges");
            param.input.usePatchNormal = inpJSON.get<bool>("SurfNavStokes.UsePatchNormals");
            if (param.input.usePatchNormal)  logger.log("using patch normals in surf integrands");
            auto formulation = inpJSON.get<std::string>("SurfNavStokes.formulation");
            if (formulation == "consistent") {
                param.input.formulation = LocalStokesParam::Formulation::consistent;
                logger.buf << "using consistent penalty formulation\n";
            }
            else {
                param.input.formulation = LocalStokesParam::Formulation::inconsistent;
                logger.buf << "using inconsistent penalty formulation\n";
            }
            StokesSystem stokesSystem;
            auto stab = inpJSON.get<std::string>("SurfNavStokes.stab");
            if (!inpJSON.get<bool>("SurfNavStokes.ExportMatrices")) {
                stokesSystem.Schur_full_stab.assemble   = (stab == "full");
                stokesSystem.Schur_normal_stab.assemble = (stab == "normal");
                stokesSystem.LB.assemble = false;
                stokesSystem.LB_stab.assemble = false;
            }
            logger.log();
        logger.end();
        logger.beg("build mesh");
            logger.beg("build initial bulk mesh");
                std::auto_ptr<DROPS::MGBuilderCL> builder(DROPS::make_MGBuilder(inpJSON));
                DROPS::MultiGridCL mg(*builder);
                const DROPS::ParamCL::ptree_type* ch= 0;
                try {
                    ch = &inpJSON.get_child("Mesh.Periodicity");
                } catch (DROPS::DROPSParamErrCL) {}
                if (ch) read_PeriodicBoundaries(mg, *ch);
                // levelset shift
                shift = inpJSON.get<DROPS::Point3DCL>("Levelset.ShiftDir", DROPS::Point3DCL(0., 0., 0.));
                shift /= shift.norm();
                auto shiftNorm = fabs(inpJSON.get<double>("Levelset.ShiftNorm", 0.));
                shift *= shiftNorm;
                mg.Transform(shiftTransform);
                logger.buf << "surface shift = " << shift;
                logger.log();
            logger.end();
            logger.beg("refine towards the surface");
                DROPS::AdapTriangCL adap(mg);
                // adaptive mesh refinement based on level set function
                DROPS::DistMarkingStrategyCL initmarker(surfNavierStokesData.surface.phi, inpJSON.get<double>("Mesh.AdaptRef.Width"), inpJSON.get<int>("Mesh.AdaptRef.CoarsestLevel"), inpJSON.get<int>("Mesh.AdaptRef.FinestLevel"));
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
            read_BndData(lsbnd, mg, inpJSON.get_child("Levelset.BndData"));
            DROPS::LevelsetP2CL & lset(*DROPS::LevelsetP2CL::Create(mg, lsbnd, sf));
            lset.CreateNumbering(mg.GetLastLevel(), &lset.idx);
            lset.Phi.SetIdx(&lset.idx);
            lset.Init(surfNavierStokesData.surface.phi);
        logger.end();
        logger.beg("build FE spaces");
            BndDataCL<Point3DCL> vbnd(0);
            read_BndData(vbnd, mg, inpJSON.get_child("Stokes.VelocityBndData"));
            BndDataCL<double> pbnd(0);
            read_BndData(pbnd, mg, inpJSON.get_child("Stokes.PressureBndData"));
            DROPS::IdxDescCL ifaceVecP2idx(vecP2IF_FE, vbnd);
            DROPS::IdxDescCL ifaceVecP1idx(vecP1IF_FE, vbnd);
            DROPS::IdxDescCL ifaceP1idx(P1IF_FE, pbnd);
            DROPS::IdxDescCL ifaceP2idx(P2IF_FE, pbnd);
            DROPS::IdxDescCL vecP2idx(vecP2_FE, vbnd);
            DROPS::IdxDescCL vecP1idx(vecP1_FE, vbnd);
            DROPS::IdxDescCL P1FEidx(P1_FE, pbnd);
            DROPS::IdxDescCL P2FEidx(P2_FE, pbnd);
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
            VecDescCL u_star, u, u_prev, u_prev_prev, p_star, p;
            if (FE == "P2P1") {
                u.SetIdx(&ifaceVecP2idx);
                u_star.SetIdx(&ifaceVecP2idx);
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
                u_star.SetIdx(&ifaceVecP1idx);
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
            else throw std::invalid_argument(FE + ": unknown FE pair");
            p.SetIdx(&ifaceP1idx);
            p_star.SetIdx(&ifaceP1idx);
            stokesSystem.gRHS.SetIdx(&ifaceP1idx);
            stokesSystem.Schur.SetIdx(&ifaceP1idx, &ifaceP1idx);
            stokesSystem.Schur_full_stab.SetIdx(&ifaceP1idx, &ifaceP1idx);
            stokesSystem.Schur_normal_stab.SetIdx(&ifaceP1idx, &ifaceP1idx);
        logger.end();
        logger.beg("interpolate initial condition");
            InitVector(mg, u_prev, surfNavierStokesData.u_T);
            stokesSystem.w = u_prev;
            logger.wrn("TODO: add normal component to wind field");
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
                << "exact velocity interp size is " << u_star.Data.size() << '\n'
                << "pre mass mtx is               " << stokesSystem.Schur.Data.num_rows() << " * " << stokesSystem.Schur.Data.num_cols() << '\n'
                << "pressure soln size is         " << p.Data.size() << '\n'
                << "exact pressure interp size is " << p_star.Data.size() << '\n'
                << "f size is                     " << stokesSystem.fRHS.Data.size() << '\n'
                << "g size is                     " << stokesSystem.gRHS.Data.size();
            logger.log();
        logger.end();
        if (inpJSON.get<bool>("SurfNavStokes.ExportMatrices")) {
            logger.beg("export matrices to " + dirName);
                MatrixCL A_final;
                A_final.LinComb(1., stokesSystem.A.Data, 1., stokesSystem.M.Data, tau_u, stokesSystem.S.Data, rho_u, stokesSystem.A_stab.Data);
                auto C_full = stokesSystem.Schur_full_stab.Data;
                C_full *= rho_p;
                auto C_n = stokesSystem.Schur_normal_stab.Data;
                C_n *= rho_p;
                std::string format = inpJSON.get<std::string>("SurfNavStokes.ExportMatricesFormat") == "mtx" ? ".mtx" : ".mat";
                auto expFunc = format == ".mtx" ? &MatrixCL::exportMTX : &MatrixCL::exportMAT;
                auto expMat = [&](MatrixCL& A, std::string const a, std::string const & b) {
                    logger.beg(a);
                        (A.*expFunc)(dirName + "/matrices/" + b + format);
                        outJSON.put("Matrices." + a, "matrices/" + b + format);
                    logger.end();
                };
                expMat(A_final, "DiffusionReaction", "A");
                expMat(stokesSystem.N.Data, "Convection", "N");
                expMat(stokesSystem.B.Data, "Divergence", "B");
                expMat(C_full, "PressureFullStab", "C_full");
                expMat(C_n, "PressureNormalStab", "C_n");
                expMat(stokesSystem.LB.Data, "VelocityScalarLaplaceBeltrami", "LB");
                expMat(stokesSystem.LB_stab.Data, "VelocityScalarLaplaceBeltramiNormalStab", "LB_stab");
                stokesSystem.LB.assemble = false;
                stokesSystem.LB_stab.assemble = false;
            logger.buf << outJSON;
                logger.log();
            logger.end();
        }
        if (param.input.exactNormal) {
            param.input.exactNormal = nullptr;
            outJSON.put("Surface.LevelsetNormalErr", sqrt(param.output.normalErrSq.lvset));
            outJSON.put("Surface.PatchNormalErr", sqrt(param.output.normalErrSq.patch));
        }
        if (param.input.exactShape) {
            param.input.exactShape = nullptr;
            outJSON.put("Surface.ShapeErr", sqrt(param.output.shapeErrSq));
        }
        if (!inpJSON.get<bool>("SurfNavStokes.Solve")) {
            output << outJSON;
            logger.buf << outJSON;
            logger.log();
            return 0;
        }
        logger.beg("build preconditioners");
            logger.beg("Schur block");
                if (stab == "full")
                    Schur_hat.LinComb(1., stokesSystem.Schur.Data, rho_p, stokesSystem.Schur_full_stab.Data);
                else if (stab == "normal")
                    Schur_hat.LinComb(1., stokesSystem.Schur.Data, rho_p, stokesSystem.Schur_normal_stab.Data);
                else {
                    logger.wrn("no stabilization is used! Schur precond = I");
                    Schur_hat = std::valarray<double>(1., stokesSystem.Schur.Data.num_rows());
                }
            logger.end();
            logger.beg("diffusion-convection-reaction block");
                typedef SSORPcCL SymmPcPcT;
                typedef PCGSolverCL<SymmPcPcT> PCGSolverT;
                typedef SolverAsPreCL<PCGSolverT> PCGPcT;
                SymmPcPcT symmPcPc_;
                std::stringstream Astream;
                PCGSolverT PCGSolver_(symmPcPc_, inpJSON.get<int>("Solver.PcAIter"), inpJSON.get<double>("Solver.PcATol"), true, &Astream);
                PCGPcT PCGPc_(PCGSolver_);//velocity preconditioner
                std::stringstream Schurstream;
                PCGSolverT SchurPCGSolver (symmPcPc_, inpJSON.get<int>("Solver.PcBIter"), inpJSON.get<double>("Solver.PcBTol"), true, &Schurstream);
                SchurPreBaseCL *spc_ = new SurfaceLaplacePreCL<PCGSolverT>(Schur_hat, SchurPCGSolver);//pressure preconditioner
                //construct symmteric block iterative solver
                typedef BlockPreCL<ExpensivePreBaseCL, SchurPreBaseCL, DiagSpdBlockPreCL>  DiagBlockPcT;
                typedef PLanczosONBCL<VectorCL, DiagBlockPcT> LanczosT;
                typedef PMResSolverCL<LanczosT> MinResT;
                DiagBlockPcT    *DBlock_ = new DiagBlockPcT(PCGPc_, *spc_);
                LanczosT *lanczos_ = new LanczosT(*DBlock_);
                MinResT *MinRes_ = new MinResT(*lanczos_, inpJSON.get<int>("Solver.Iter"), inpJSON.get<double>("Solver.Tol"), /*relative*/ false);
                BlockMatrixSolverCL<MinResT> *symStokesSolver= new BlockMatrixSolverCL<MinResT>(*MinRes_);
                // construct preconditioners for nonsymmetric case
                typedef /*JACPcCL*/ SSORPcCL NonSymmPcPcT;
                typedef GMResSolverCL<NonSymmPcPcT> GMResSolverT;
                typedef SolverAsPreCL<GMResSolverT> GMResPcT;
                NonSymmPcPcT nonsymmPcPc_;
                std::stringstream PCstream;
                GMResSolverT GMResSolver_(nonsymmPcPc_, inpJSON.get<int>("Solver.PcAIter"), inpJSON.get<int>("Solver.PcAIter"), inpJSON.get<double>("Solver.PcATol"),
                         /*bool relative=*/ true, /*bool calculate2norm=*/ false, /*PreMethGMRES method=*/ LeftPreconditioning,
                             /*bool mod =*/ true,      /*bool useModGS =*/ false, &PCstream);
                GMResPcT GMResPc_nonsym(GMResSolver_);
                //setup iterative block solver for nonsymmetric case
                typedef GMResSolverCL<DiagBlockPcT> GMResT;
                DiagBlockPcT    *DBlock_nonsym = new DiagBlockPcT(GMResPc_nonsym, *spc_);
                GMResT *GMRes_ = new GMResT(*DBlock_nonsym, inpJSON.get<int>("Solver.Iter"), inpJSON.get<int>("Solver.Iter"), inpJSON.get<double>("Solver.Tol"),
                         /*bool relative=*/ false, /*bool calculate2norm=*/ false, /*PreMethGMRES method=*/ LeftPreconditioning,
                              /*bool mod =*/ true,      /*bool useModGS =*/ false);
                BlockMatrixSolverCL<GMResT> *nonsymStokesSolver = new BlockMatrixSolverCL<GMResT>(*GMRes_);
                StokesSolverBaseCL* stokesSolver;
            logger.end();
        logger.end();
        logger.beg("solve");
            VectorCL I_p;
            I_p.resize(stokesSystem.gRHS.Data.size(),1.);
            MatrixCL A_dyn, C, B_T;
            transpose(stokesSystem.B.Data, B_T);
            if (stab == "full")
                C.LinComb(0., stokesSystem.Schur.Data, -rho_p, stokesSystem.Schur_full_stab.Data);
            else if (stab == "normal")
                C.LinComb(0., stokesSystem.Schur.Data, -rho_p, stokesSystem.Schur_normal_stab.Data);
            else
                C.LinComb(0., stokesSystem.Schur.Data, 0., stokesSystem.Schur.Data);
            double Pe; // peclet number
            VTKOutCL* vtkWriter = nullptr;
            vtkWriter = new VTKOutCL(mg, "DROPS data", numbOfStepsVTK, dirName + "/vtk", testName + "_", testName, 0);
            vtkWriter->Register(make_VTKScalar(lset.GetSolution(), "level-set"));
            vtkWriter->Register(make_VTKIfaceVector(mg, u_star, "u_*", velFE, vbnd));
            vtkWriter->Register(make_VTKIfaceVector(mg, stokesSystem.fRHS, "f_T", velFE, vbnd));
            vtkWriter->Register(make_VTKIfaceVector(mg, u, "u_h", velFE, vbnd));
            vtkWriter->Register(make_VTKIfaceScalar(mg, p, "p_h", /*preFE,*/ pbnd));
            vtkWriter->Register(make_VTKIfaceScalar(mg, stokesSystem.gRHS, "-g", /*preFE,*/ pbnd));
            vtkWriter->Register(make_VTKIfaceScalar(mg, p_star, "p_*", /*preFE,*/ pbnd));
            if (everyStep > 0) {
                logger.beg("write initial condition to vtk");
                    u_star = u_prev;
                    u = u_prev;
                    vtkWriter->Write(0);
                logger.end();
            }
            logger.beg("t = t_1");
                logger.beg("linear solve");
                    stokesSystem.fRHS.Data += (1. / stepSize) * (stokesSystem.M.Data * u_prev.Data);
                    stokesSystem.gRHS.Data -= (dot(stokesSystem.gRHS.Data, I_p) / dot(I_p, I_p)) * I_p;
                    A_dyn.LinComb(1. / stepSize, stokesSystem.M.Data, 1., stokesSystem.N.Data, nu, stokesSystem.A.Data, tau_u, stokesSystem.S.Data, rho_u, stokesSystem.A_stab.Data);
                    stokesSolver = symStokesSolver;
                    if (windRises(nu, stokesSystem.w, Pe)) stokesSolver = nonsymStokesSolver;
                    logger.buf << "Pe = " << Pe << ", using " << (stokesSolver == symStokesSolver ? "" : "non-") << "symmetric solver";
                    logger.log();
                    stokesSolver->Solve(A_dyn, stokesSystem.B.Data, C, u.Data, p.Data, stokesSystem.fRHS.Data, stokesSystem.gRHS.Data, u.RowIdx->GetEx(), p.RowIdx->GetEx());
                logger.end();
                logger.beg("output");
                    auto exportStats = [&](size_t i) {
                        std::ofstream stats(dirName + "/stats/t_" + std::to_string(i) + ".json");
                        ParamCL tJSON;
                        auto velResSq = norm_sq(A_dyn * u.Data + B_T * p.Data - stokesSystem.fRHS.Data);
                        auto preResSq = norm_sq(stokesSystem.B.Data * u.Data + C * p.Data - stokesSystem.gRHS.Data);
                        auto residual = sqrt(velResSq + preResSq);
                        tJSON.put("Solver.ResidualNorm.True.Full", residual);
                        tJSON.put("Solver.ResidualNorm.True.Velocity", sqrt(velResSq));
                        tJSON.put("Solver.ResidualNorm.True.Pressure", sqrt(preResSq));
                        tJSON.put("Solver.ResidualNorm.GMRES", stokesSolver->GetResid());
                        tJSON.put("Solver.TotalIters", stokesSolver->GetIter());
                        tJSON.put("Solver.DOF.Velocity", stokesSystem.A.Data.num_rows());
                        tJSON.put("Solver.DOF.Pressure", stokesSystem.Schur.Data.num_rows());
                        auto t = i * stepSize;
                        tJSON.put("Time", t);
                        tJSON.put("h", h);
                        tJSON.put("Peclet", Pe);
                        tJSON.put("MeshDepParams.rho_p", rho_p);
                        tJSON.put("MeshDepParams.rho_u", rho_u);
                        tJSON.put("MeshDepParams.tau_u", tau_u);
                        InitVector(mg, u_star, surfNavierStokesData.u_T, t);
                        InitScalar(mg, p_star, surfNavierStokesData.p, t);
                        p.Data -= dot(stokesSystem.Schur.Data * p.Data, I_p) / dot(stokesSystem.Schur.Data * I_p, I_p) * I_p;
                        VectorCL u_star_minus_u = u_star.Data - u.Data, p_star_minus_p = p_star.Data - p.Data;
                        auto velL2 = sqrt(dot(u.Data, stokesSystem.M.Data * u.Data));
                        auto velL2err = dot(u_star_minus_u, stokesSystem.M.Data * u_star_minus_u);
                        auto velNormalL2 = dot(u.Data, stokesSystem.S.Data * u.Data);
                        auto velTangenL2 = sqrt(velL2err - velNormalL2);
                        velL2err = sqrt(velL2err);
                        velNormalL2 = sqrt(velNormalL2);
                        auto velH1err = sqrt(dot(u_star_minus_u, stokesSystem.A.Data * u_star_minus_u));
                        auto preL2 = sqrt(dot(p.Data, stokesSystem.Schur.Data * p.Data));
                        auto preL2err = sqrt(dot(p_star_minus_p, stokesSystem.Schur.Data * p_star_minus_p));
                        tJSON.put("Integral.Error.VelocityL2", velL2err);
                        tJSON.put("Integral.Error.VelocityTangentialL2", velTangenL2);
                        tJSON.put("Integral.Error.VelocityNormalL2", velNormalL2);
                        tJSON.put("Integral.Error.VelocityH1", velH1err);
                        tJSON.put("Integral.Error.PressureL2", preL2err);
                        tJSON.put("Integral.PressureL2", preL2);
                        tJSON.put("Integral.VelocityL2", velL2);
                        tJSON.put("Integral.KineticEnergy", .5 * velL2 * velL2);
                        stats << tJSON;
                        logger.buf << tJSON;
                        logger.log();
                        if (everyStep > 0 && (i-1) % everyStep == 0) {
                            logger.beg("write vtk");
                                vtkWriter->Write(t);
                            logger.end();
                        }
                    };
                    exportStats(1);
                logger.end();
            logger.end();
            // no need to re-assemble these matrices:
            stokesSystem.A.assemble = false;
            stokesSystem.A_stab.assemble = false;
            stokesSystem.B.assemble = false;
            stokesSystem.M.assemble = false;
            stokesSystem.S.assemble = false;
            stokesSystem.Schur.assemble = false;
            stokesSystem.Schur_full_stab.assemble = false;
            stokesSystem.Schur_normal_stab.assemble = false;
            for (size_t i = 2; i <= numbOfSteps; ++i) {
                std::stringstream header;
                header << "t = t_" << i << " (" << (100. * i) / numbOfSteps << "%)";
                logger.beg(header.str());
                    logger.beg("assemble");
                        u_prev_prev = u_prev;
                        u_prev = u;
                        param.input.t = i * stepSize;
                        stokesSystem.w.Data = 2. * u_prev.Data - u_prev_prev.Data;
                        SetupStokesIF_P2P1(mg, lset, &stokesSystem, &param);
                        stokesSystem.fRHS.Data += (2. / stepSize) * (stokesSystem.M.Data * u_prev.Data) - (.5 / stepSize) * (stokesSystem.M.Data * u_prev_prev.Data);
                        stokesSystem.gRHS.Data -= (dot(stokesSystem.gRHS.Data, I_p) / dot(I_p, I_p)) * I_p;
                        A_dyn.LinComb(1.5 / stepSize, stokesSystem.M.Data, 1., stokesSystem.N.Data, nu, stokesSystem.A.Data, tau_u, stokesSystem.S.Data, rho_u, stokesSystem.A_stab.Data);
                    logger.end();
                    logger.beg("linear solve");
                        stokesSolver = symStokesSolver;
                        if (windRises(nu, stokesSystem.w, Pe)) stokesSolver = nonsymStokesSolver;
                        logger.buf << "Pe = " << Pe << ", using " << (stokesSolver == symStokesSolver ? "" : "non-") << "symmetric solver";
                        logger.log();
                        stokesSolver->Solve(A_dyn, stokesSystem.B.Data, C, u.Data, p.Data, stokesSystem.fRHS.Data, stokesSystem.gRHS.Data, u.RowIdx->GetEx(), p.RowIdx->GetEx());
                    logger.end();
                    logger.beg("output");
                        exportStats(i);
                    logger.end();
                logger.end();
            }
        logger.end();
        // error
        double Aaverage(0.), Avariation(0.), Schuraverage(0.), Schurvariation(0.);
        double PCaverage(0.), PCvariation(0.);
        std::stringstream Astreamcopy(Astream.str());
        std::stringstream PCstreamcopy(PCstream.str());
        std::stringstream Schurstreamcopy(Schurstream.str());
        RightComputeAverageIterations(Astream, Aaverage);
        ComputeAverageIterations(Schurstream, Schuraverage);
        RightComputeAverageIterations(PCstream, PCaverage);
        RightComputeVariationFromAverageIterations(Astreamcopy, Aaverage, Avariation);
        RightComputeVariationFromAverageIterations(PCstreamcopy, PCaverage, PCvariation);
        ComputeVariationFromAverageIterations(Schurstreamcopy, Schuraverage, Schurvariation);
        logger.buf << "The average iterationsnumber of the A-preconditioner is: " <<     Aaverage     << '\n' << " ...with a variation of: " << Avariation << '\n';
        logger.buf << "The average iterationsnumber of the nonsymmetric A-preconditioner is: " <<     PCaverage     << '\n' << " ...with a variation of: " << PCvariation << '\n';
        logger.buf << "The average iterationsnumber of the Schur-preconditioner is: " << Schuraverage << '\n' << " ...with a variation of: " << Schurvariation << '\n';
        logger.log();
        delete &lset;
        return 0;
    } catch (std::exception const & e) {
        logger.err(e.what());
    } catch (DROPS::DROPSErrCL const & e) {
        e.what(logger.buf);
        logger.err(logger.buf.str());
    } catch (...) {
        logger.err("unknown error");
    }
    return 1;
}