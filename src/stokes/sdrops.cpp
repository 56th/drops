/// \file sdrops.cpp
/// \brief stokes problem
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt, Thorolf Schulte; SC RWTH Aachen: Oliver Fortmeier

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
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "stokes/stokes.h"
#include "num/oseensolver.h"
#include <fstream>

#include "misc/params.h"

DROPS::ParamCL P;
inline DROPS::SVectorCL<3> LsgVel(const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret;
    ret[0]=    std::sin(p[0])*std::sin(p[1])*std::sin(p[2]);
    ret[1]=  - std::cos(p[0])*std::cos(p[1])*std::sin(p[2]);
    ret[2]= 2.*std::cos(p[0])*std::sin(p[1])*std::cos(p[2]);
    return ret/3.;
}

inline double LsgPr(const DROPS::Point3DCL& p, double)
{
    return std::cos(p[0])*std::sin(p[1])*std::sin(p[2]);
//     return 1.;
}

// boundary value functions (in 2D-bnd-coords)
//DROPS::StokesVelBndDataCL::bnd_type bnd_val_e2e3(const Point2DCL&);
//DROPS::StokesVelBndDataCL::bnd_type bnd_val_e1e3(const Point2DCL&);
//DROPS::StokesVelBndDataCL::bnd_type bnd_val_e1e2(const Point2DCL&);


// q*u - nu*laplace u + Dp = f
//                  -div u = 0
class StokesCoeffCL
{
  public:
    static double q(const DROPS::TetraCL&, const DROPS::BaryCoordCL&, double) { return 0.0; }
    static DROPS::SVectorCL<3> f(const DROPS::TetraCL& t, const DROPS::BaryCoordCL& b, double)
    {
        const DROPS::Point3DCL p= DROPS::GetWorldCoord( t, b);
        DROPS::SVectorCL<3> ret(0.0);
        ret[2]= 3.*std::cos(p[0])*std::sin(p[1])*std::cos(p[2]);
        return ret;
    }
/*    {
        const DROPS::Point3DCL p= DROPS::GetWorldCoord( t, b);
        SVectorCL<3> ret;
        ret[0]=    std::sin(p[0])*std::sin(p[1])*std::sin(p[2]);
        ret[1]=  - std::cos(p[0])*std::cos(p[1])*std::sin(p[2]);
        ret[2]= 2.*std::cos(p[0])*std::sin(p[1])*std::cos(p[2]);
        return ret;
    }
*/
    const double nu;

    StokesCoeffCL() : nu(1.0) {}
};

    typedef DROPS::StokesP1BubbleP1CL<StokesCoeffCL>
            StokesOnBrickCL;
//    typedef DROPS::StokesP2P1CL<StokesCoeffCL>
//            StokesOnBrickCL;
    typedef StokesOnBrickCL MyStokesCL;

namespace DROPS // for Strategy
{
using ::MyStokesCL;

template<class Coeff>
void Strategy(StokesP1BubbleP1CL<Coeff>& Stokes, double omega, double inner_iter_tol, Uint maxStep, double rel_red)
// flow control
{
    MultiGridCL& MG= Stokes.GetMG();
    const typename MyStokesCL::BndDataCL::PrBndDataCL& PrBndData= Stokes.GetBndData().Pr;
    const typename MyStokesCL::BndDataCL::VelBndDataCL& VelBndData= Stokes.GetBndData().Vel;

    MLIdxDescCL  loc_vidx, loc_pidx;
    MLIdxDescCL* vidx1= &Stokes.vel_idx;
    MLIdxDescCL* pidx1= &Stokes.pr_idx;
    MLIdxDescCL* vidx2= &loc_vidx;
    MLIdxDescCL* pidx2= &loc_pidx;
    VecDescCL     loc_p;
    VelVecDescCL  loc_v;
    VelVecDescCL* v1= &Stokes.v;
    VelVecDescCL* v2= &loc_v;
    VecDescCL*    p1= &Stokes.p;
    VecDescCL*    p2= &loc_p;
    VelVecDescCL* b= &Stokes.b;
    VelVecDescCL* c= &Stokes.c;
    MLMatDescCL* A= &Stokes.A;
    MLMatDescCL* B= &Stokes.B;
    Uint step= 0;
    // measure of cube: (Pi/4)^3==0.484...
    StokesDoerflerMarkCL<typename MyStokesCL::est_fun, MyStokesCL>
        Estimator(rel_red, 0., .484473073129685, true, &MyStokesCL::ResidualErrEstimator, Stokes );
    bool new_marks;
//    double akt_glob_err;

    vidx1->SetFE( vecP1Bubble_FE);
    vidx2->SetFE( vecP1Bubble_FE);
    pidx1->SetFE( P1_FE);
    pidx2->SetFE( P1_FE);

    TimerCL time;
//    err_idx->Set(0, 0, 0, 1);
    do
    {
//        akt_glob_err= glob_err;
//        MarkAll(MG);
        MG.Refine();
        Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx1);
        Stokes.CreateNumberingPr( MG.GetLastLevel(), pidx1);
        std::cout << "altes und neues TriangLevel: " << vidx2->TriangLevel() << ", "
                  << vidx1->TriangLevel() << std::endl;
        MG.SizeInfo(std::cout);
        b->SetIdx(vidx1);
        c->SetIdx(pidx1);
        p1->SetIdx(pidx1);
        v1->SetIdx(vidx1);
        std::cout << "Anzahl der Druck-Unbekannten: " << p2->Data.size() << ", "
                  << p1->Data.size() << std::endl;
        std::cout << "Anzahl der Geschwindigkeitsunbekannten: " << v2->Data.size() << ", "
                  << v1->Data.size() << std::endl;

        if (p2->RowIdx)
        {
            const StokesBndDataCL& BndData= Stokes.GetBndData();
            P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, const VecDescCL>  pr2(p2, &BndData.Pr, &MG);
            P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, VecDescCL>        pr1(p1, &BndData.Pr, &MG);
//            P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, const VelVecDescCL> vel2(v2, &BndData.Vel, &MG);
//            P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, VelVecDescCL>       vel1(v1, &BndData.Vel, &MG);
            Interpolate(pr1, pr2);
//            Interpolate(vel1, vel2);
//            CheckSolution(v1,p1,&LsgVel,&LsgPr);
            v2->Reset();
            p2->Reset();
        }
        A->SetIdx(vidx1, vidx1);
        B->SetIdx(pidx1, vidx1);
        time.Reset();
        time.Start();
        Stokes.SetupSystem(A, b, B, c);
        time.Stop();
        std::cout << time.GetTime() << " seconds for setting up all systems!" << std::endl;
        time.Reset();
        time.Start();
        A->Data * v1->Data;
        time.Stop();
        std::cout << " A*x took " << time.GetTime() << " seconds!" << std::endl;
        time.Reset();
        time.Start();
        transp_mul( A->Data, v1->Data);
        time.Stop();
        std::cout << "AT*x took " << time.GetTime() << " seconds!" << std::endl;
/*
        { // write system in files for MatLab
            std::ofstream Adat("Amat.dat"), Bdat("Bmat.dat"), bdat("fvec.dat"), cdat("gvec.dat");
            Adat << A->Data;   Bdat << B->Data;    bdat << b->Data;    cdat << c->Data;
        }
*/
        Stokes.GetDiscError(&LsgVel, &LsgPr);
//std::cout << A->Data << std::endl << b->Data << std::endl;
/*        double half= M_PI/8;
        MultiGridCL::TriangVertexIteratorCL vert= MG.GetTriangVertexBegin(A->GetRowLevel());
        while (vert->GetCoord()[0]!=half || vert->GetCoord()[1]!=half || vert->GetCoord()[2]!=half) ++vert;
        IdxT unk= vert->Unknowns(A->RowIdx->Idx);
        std::cout << vert->GetCoord() << " has index " << unk << std::endl;
        std::cout << "A(i,i) = " << A->Data(unk,unk) <<std::endl;
        std::cout << "B(i,j) = " << B->Data(vert->Unknowns(B->RowIdx->Idx),unk) << std::endl;
*/
        int max_iter;
        double tol;
        time.Reset();

        SSORPcCL  pc(omega);
        MLMatDescCL prM;
        prM.SetIdx( pidx1, pidx1);
        Stokes.SetupPrMass( &prM);
        PreGSOwnMatCL<P_SSOR0> schur_pc(prM.Data.GetFinest());
        SSORPcCL poissonpc;
        typedef PCGSolverCL<SSORPcCL> PCG_SsorCL;
        PCG_SsorCL poissonsolver( poissonpc, 500, inner_iter_tol);
        SchurComplMatrixCL<PCG_SsorCL, MLMatrixCL, DummyExchangeCL> BABT( poissonsolver, A->Data, B->Data, DummyExchangeCL());
        double outer_tol;
        std::cout << "tol = "; std::cin >> outer_tol;
        time.Start();
        VectorCL rhs( -c->Data);
        {
            double tol= inner_iter_tol;
            int max_iter= 200;
            VectorCL tmp(vidx1->NumUnknowns());
            PCG(A->Data, tmp, b->Data, DummyExchangeCL(), pc, max_iter, tol);
            std::cout << "Iterationen: " << max_iter << "    Norm des Residuums: " << tol << std::endl;
            rhs+= B->Data*tmp;
        }
        std::cout << "rhs has been set!" << std::endl;
//            tol= 1.0e-14;
        max_iter= 200;
        tol= outer_tol;
//        PCG(A->Data, new_x->Data, b->Data, DummyExchangeCL(), pc, max_iter, tol);
        PCG(BABT, p1->Data, rhs, DummyExchangeCL(), schur_pc, max_iter, tol);
        std::cout << "Iterationen: " << max_iter << "    Norm des Residuums: " << tol << std::endl;

        tol= outer_tol;
        max_iter= 200;
        PCG(A->Data, v1->Data, VectorCL( b->Data - transp_mul(B->Data, p1->Data)), DummyExchangeCL(), pc, max_iter, tol);
        time.Stop();

        std::cout << "Das Verfahren brauchte "<<time.GetTime()<<" Sekunden.\n";
        Stokes.CheckSolution(v1, p1, &LsgVel, &LsgPr);
        typename MyStokesCL::const_DiscPrSolCL  pr(p1, &PrBndData, &MG);
        typename MyStokesCL::const_DiscVelSolCL vel(v1, &VelBndData, &MG);
        if (step==0)
        {
            Estimator.Init( pr, vel);
        }
        new_marks= Estimator.Estimate( pr, vel);
        A->Reset();
        B->Reset();
        b->Reset();
        c->Reset();
//        std::cout << "Loesung Druck: " << p1->Data << std::endl;
//        CreateNumbering(new_idx->TriangLevel, err_idx);
//        err->SetIdx(err_idx);
//        NumMarked= EstimateError(new_x, 0.2, &Estimator);
//TODO: Fehler schaetzen
//        new_marks= EstimateError(new_x, 1.5, akt_glob_err, &ResidualErrEstimator);
//        err->Reset();
//        DeleteNumbering(err_idx);
        std::swap(v2, v1);
        std::swap(p2, p1);
        std::swap(vidx2, vidx1);
        std::swap(pidx2, pidx1);
        std::cout << std::endl;
    }
    while (new_marks && ++step<maxStep);
    // we want the solution to be in Stokes.v, Stokes.pr
    if (v2 == &loc_v)
    {
        Stokes.vel_idx.swap( loc_vidx);
        Stokes.pr_idx.swap( loc_pidx);
        Stokes.v.SetIdx(&Stokes.vel_idx);
        Stokes.p.SetIdx(&Stokes.pr_idx);

        Stokes.v.Data= loc_v.Data;
        Stokes.p.Data= loc_p.Data;
    }
}

} // end of namespace DROPS


int main (int argc, char** argv)
{
  try
  {
    if (argc!=4)
    {
        std::cout << "You have to specify three parameters:\n\tsdrops <omega> <inner_iter_tol> <relative error reduction>" << std::endl;
        return 1;
    }

    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= M_PI/4.;


    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 1, 1, 1);
//    DROPS::BBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 2, 2, 2);
//    DROPS::LBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 2, 2);
    const bool IsNeumann[6]=
        {false, false, false, false, false, false};
    const DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &LsgVel, &LsgVel, &LsgVel, &LsgVel, &LsgVel, &LsgVel};

    StokesOnBrickCL prob(brick, StokesCoeffCL(), DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
    DROPS::MultiGridCL& mg = prob.GetMG();
    DROPS::RBColorMapperCL colormap;
    double omega= std::atof(argv[1]);
    double inner_iter_tol= std::atof(argv[2]);
    double rel_red= std::atof(argv[3]);
    std::cout << "Omega: " << omega << " inner iter tol: " << inner_iter_tol << " rel. error reduction: " << rel_red << std::endl;
    Strategy(prob, omega, inner_iter_tol, 8, rel_red);
    std::cout << "hallo" << std::endl;
    std::cout << DROPS::SanityMGOutCL(mg) << std::endl;
    std::ofstream fil("ttt.off");
    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    fil << DROPS::GeomSolOutCL<MyStokesCL::const_DiscPrSolCL>(mg, prob.GetPrSolution(), &colormap, -1, false, 0.0, min, max) << std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
