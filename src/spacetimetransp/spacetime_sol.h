/// \file spacetime_sol.h
/// \brief classes that handle space time solutions (vectors, trace vectors, output, etc..)
/// \author LNM RWTH Aachen: Christoph Lehrenfeld

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
 * Copyright 2012 LNM/SC RWTH Aachen, Germany
*/

#ifndef DROPS_SPACETIME_SOL_H
#define DROPS_SPACETIME_SOL_H

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "misc/problem.h"
#include "num/discretize.h"
#include "misc/params.h"
#include "levelset/levelset.h"

// #include <list>
// #include <vector>
// #include <algorithm>
// #include <iostream>
// #include <iterator>
// #include "misc/utils.h"
// #include "geom/boundary.h"
// #include "geom/topo.h"
// #include "geom/simplex.h"
// #include "geom/principallattice.h"
// #include "num/unknowns.h"
// 


namespace DROPS
{

// fwd declarations:
// class PentatopeCL;
// class GeneralizedPrism4CL;
// class HyperTrigCL;
// class Tetra4DCL;

// typedef double    (*instat_scalar_fun_ptr)(const DROPS::Point3DCL&, double);
// typedef std::pair<double,double> TimeInterval;

typedef double    (*instat_scalar_fun_ptr)(const DROPS::Point3DCL&, double);
typedef std::pair<double,double> TimeInterval;

class SpaceTimeXSolutionCL
{
private:
    MultiGridCL& mg_;
    const BndDataCL<> & Bndneg_;
    const BndDataCL<> & Bndpos_;

    const double stxfem_omit_param_;
    const double weight_pos_;
    const double weight_neg_;

    IdxDescCL st_x_idx_; // space time p1p1x index
    IdxDescCL fut_x_idx_; // p1x index according to new levelset
    IdxDescCL past_x_idx_; // p1x index according to old levelset
    /* IdxDescCL neg_space_idx_; // std. p1 index */
    /* IdxDescCL pos_space_idx_; // std. p1 index */

    
    /* double lasttraces_time_past; */
    /* double lasttraces_time_future; */

    const LevelsetP2CL * lsetp2old;
    const LevelsetP2CL * lsetp2new;

    VecDescCL sol;
    VecDescCL futurevec, pastvec;

    IdxDescCL p1idxn;
    IdxDescCL p1idxp;
    VecDescCL future_vec_pos;
    VecDescCL future_vec_neg;
    VecDescCL past_vec_pos;
    VecDescCL past_vec_neg;

public:

    typedef P2EvalCL<double, const BndDataCL<>, const VecDescCL> const_DiscSolCL;

    SpaceTimeXSolutionCL( MultiGridCL& mg, const BndDataCL<>& Bndneg, const BndDataCL<>& Bndpos, ParamCL & P);
    void UpdateTimeSlab( const LevelsetP2CL & lsetold, const LevelsetP2CL & lsetnew, double vmax = 1.0, instat_scalar_fun_ptr lset_fpt = NULL);

    void EvalFutureTrace();
    void EvalPastTrace();
    void EvalTraces(){ 
        EvalFutureTrace(); 
        EvalPastTrace();
    }

    IdxDescCL & GetIdx()
    { return st_x_idx_; }

    const VecDescCL & GetFutureTrace() const
    { return futurevec; }
    const VecDescCL & GetPastTrace() const
    { return pastvec; }
    VecDescCL & GetFutureTrace()
    { return futurevec; }
    VecDescCL & GetPastTrace()
    { return pastvec; }

    const VecDescCL & GetFutureTrace_Pos() const
    { return future_vec_pos; }
    const VecDescCL & GetFutureTrace_Neg() const
    { return future_vec_neg; }
    const VecDescCL & GetPastTrace_Pos() const
    { return past_vec_pos; }
    const VecDescCL & GetPastTrace_Neg() const
    { return past_vec_neg; }
    VecDescCL & GetFutureTrace_Pos()
    { return future_vec_pos; }
    VecDescCL & GetFutureTrace_Neg()
    { return future_vec_neg; }
    VecDescCL & GetPastTrace_Pos()
    { return past_vec_pos; }
    VecDescCL & GetPastTrace_Neg()
    { return past_vec_neg; }

    const VecDescCL & GetSolution() const
    { return sol; }
    VecDescCL & GetSolution()
    { return sol; }

    virtual ~SpaceTimeXSolutionCL();



    inline double CheckFutureSolution(instat_scalar_fun_ptr Lsgn, instat_scalar_fun_ptr Lsgp, double time)
    {
        IdxDescCL idx;
        VecDescCL cn(&idx);
        VecDescCL cp(&idx);
        if (idx.NumUnknowns() != 0)
            idx.DeleteNumbering( mg_);
        idx.CreateNumbering( futurevec.RowIdx->TriangLevel(), mg_, *futurevec.RowIdx);
        P1XtoP1 ( *futurevec.RowIdx, futurevec.Data, idx, cp.Data, cn.Data, lsetp2new->Phi, mg_);
        cp.Data *= 1.0/weight_pos_;
        cn.Data *= 1.0/weight_neg_;
        cp.t = lsetp2new->Phi.t;
        cn.t = lsetp2new->Phi.t;

        double errl2p = 0, errl2n = 0;
        double errl1p = 0, errl1n = 0;
        double massn = 0, massp = 0;
        double voln = 0, volp = 0;
        InterfaceTetraCL patch;

        LocalP2CL<double> ones( 1.);
        std::cout << "Difference to exact solution:" << std::endl;

        /* const Uint lvl= ct.GetLevel(); */
        DROPS_FOR_TRIANG_TETRA( mg_, -1, it)
        {
            LocalP2CL<double> lp2_soln (*it, Lsgn, time);
            LocalP2CL<double> lp2_solp (*it, Lsgp, time);
            LocalP1CL<double> lp1_p (*it, cp, Bndpos_);
            LocalP2CL<double> lp2_p (lp1_p);
            LocalP1CL<double> lp1_n (*it, cn, Bndneg_);
            LocalP2CL<double> lp2_n (lp1_n);
            SMatrixCL<3,3> M;
            double det;
            GetTrafoTr(M,det,*it);
            double absdet= std::fabs(det);
            patch.Init( *it,  lsetp2new->Phi,0.);
            if (patch.Intersects()){
                patch.ComputeSubTets();
                Uint NumTets=patch.GetNumTetra(); //# of subtetras
                for (Uint k=0; k<NumTets; ++k)
                {
                    bool pPart = (k>=patch.GetNumNegTetra());
                    const SArrayCL<BaryCoordCL,4>& T =patch.GetTetra(k);
                    BaryCoordCL* nodes = Quad3CL<>::TransformNodes(T);
                    Quad3CL<> q3_sol = Quad3CL<>(pPart?lp2_solp : lp2_soln, nodes);
                    Quad3CL<> q3_dsol = Quad3CL<>(pPart? lp2_p : lp2_n, nodes);
                    Quad3CL<> q3_diff; q3_diff = q3_dsol - q3_sol;
                    Quad3CL<> q3_diff2; q3_diff2 = q3_diff * q3_diff;
                    Quad3CL<> q3_diffabs; q3_diffabs = std::abs(q3_dsol - q3_sol);

                    double Vol = absdet*VolFrac(T);
                    if (pPart)
                        errl2p += q3_diff2.quad(Vol);
                    else
                        errl2n += q3_diff2.quad(Vol);
                    if (pPart)
                        errl1p += q3_diffabs.quad(Vol);
                    else
                        errl1n += q3_diffabs.quad(Vol);

                    if (pPart)
                        massp += q3_dsol.quad(Vol);
                    else
                        massn += q3_dsol.quad(Vol);

                    if (pPart)
                        volp += 1.0/6.0 * Vol;
                    else
                        voln += 1.0/6.0 * Vol;
                    delete[] nodes;
                }

            }
            else
            {
                bool pPart= (patch.GetSign( 0) == 1);
                Quad3CL<> q3_sol = Quad3CL<>(pPart?lp2_solp : lp2_soln);
                Quad3CL<> q3_dsol = Quad3CL<>(pPart? lp2_p : lp2_n);
                Quad3CL<> q3_diff; q3_diff = q3_dsol - q3_sol;
                Quad3CL<> q3_diff2; q3_diff2 = q3_diff * q3_diff;
                Quad3CL<> q3_diffabs; q3_diffabs = std::abs(q3_dsol - q3_sol);
                if (pPart)
                    errl2p += q3_diff2.quad(absdet);
                else
                    errl2n += q3_diff2.quad(absdet);
                if (pPart)
                    errl1p += q3_diffabs.quad(absdet);
                else
                    errl1n += q3_diffabs.quad(absdet);
                
                if (pPart)
                    massp += q3_dsol.quad(absdet);
                else
                    massn += q3_dsol.quad(absdet);

                if (pPart)
                    volp += absdet / 6.0;
                else
                    voln += absdet / 6.0;
            }
        }



        double maxval = -1e99;
        double dmaxval = -1e99;
        double lsetv_maxdiff = -1e99;
        double maxdiff = -1e99;
        
        Point3DCL maxcoord;
        Point3DCL dmaxcoord;
        Point3DCL dmaxdiff;
        double dmaxdiff_val = 0, dmaxdiff_dval = 0;
        Uint idxn = idx.GetIdx();
        Uint lset_idxn = lsetp2new->Phi.RowIdx->GetIdx();
//        Uint idxn = idx.GetIdx();
        DROPS_FOR_TRIANG_VERTEX( mg_, -1, it)
        {
            IdxT unkn = (*it).Unknowns.Exist(idxn) ? (*it).Unknowns(idxn) : NoIdx;            
            IdxT lset_unkn = (*it).Unknowns.Exist(lset_idxn) ? (*it).Unknowns(lset_idxn) : NoIdx;            
            if (unkn == NoIdx)
                continue;
            if ((lset_unkn == NoIdx) || (lsetp2new->Phi.Data[lset_unkn] < 0))
                continue;
            const double val = Lsgp((*it).GetCoord(),time);

            if (val > maxval)
            {
                maxval = val;
                maxcoord = (*it).GetCoord();
            }

            const double dval = cp.Data[unkn];

            if (dval > dmaxval)
            {
                dmaxval = dval;
                dmaxcoord = (*it).GetCoord();
            }

            const double diff = std::abs(dval-val);
            if (diff > maxdiff)
            {
                maxdiff = diff;
                lsetv_maxdiff = lsetp2new->Phi.Data[lset_unkn];
                dmaxdiff = (*it).GetCoord();
                dmaxdiff_val = val;
                dmaxdiff_dval = dval;
            }
            
            /* LocalP2CL<double> lp2_soln (*it, Lsgn, time); */
            /* LocalP2CL<double> lp2_solp (*it, Lsgp, time); */
            /* LocalP1CL<double> lp1_p (*it, cp, Bndpos_); */
            /* LocalP2CL<double> lp2_p (lp1_p); */
            /* LocalP1CL<double> lp1_n (*it, cn, Bndneg_); */
            /* LocalP2CL<double> lp2_n (lp1_n); */
        }

        std::cout << " time = " << time << std::endl;
        std::cout << " dmaxval = " << dmaxval;
        std::cout << " at " << dmaxcoord << std::endl;
        std::cout << " maxval = " << maxval;
        std::cout << " at " << maxcoord << std::endl;
        std::cout << " maxdiff  = " << maxdiff;
        std::cout << " with val =" << dmaxdiff_val;
        std::cout << " and dval = " << dmaxdiff_dval;
        std::cout << " at " << dmaxdiff;
        std::cout << " with lset " << lsetv_maxdiff << std::endl;
        double errl2 = std::sqrt(errl2n + errl2p);
        double errl1 = (errl1n + errl1p);
        std::cout << "errl2p = " << std::sqrt(errl2p) << "\t";
        std::cout << "errl2n = " << std::sqrt(errl2n) << "\t";
        std::cout << "errl2  = " << errl2 << std::endl;
        std::cout << "errl1p = " << errl1p << "\t";
        std::cout << "errl1n = " << errl1n << "\t";
        std::cout << "errl1  = " << errl1 << std::endl;

        std::cout << "mass in neg = " << massn << std::endl;
        std::cout << "mass in pos = " << massp << std::endl;
        std::cout << "total mass = " << massn + massp << std::endl;
        std::cout << "vol in neg = " << voln << std::endl;
        std::cout << "vol in pos = " << volp << std::endl;

        return errl2;
    }



    // hack for reference solution calculations
    inline double CalcL2NormOfFuture()
    {
        VecDescCL &cn(GetFutureTrace_Neg());
        VecDescCL &cp(GetFutureTrace_Pos());
        /* IdxDescCL & idx = *(cn.RowIdx); */

        double errl2p = 0, errl2n = 0;
        InterfaceTetraCL patch;

        LocalP2CL<double> ones( 1.);
        std::cout << "Difference to reference solution:" << std::endl;

        /* const Uint lvl= ct.GetLevel(); */
        DROPS_FOR_TRIANG_TETRA( mg_, -1, it)
        {
            LocalP1CL<double> lp1_p (*it, cp, Bndpos_);
            LocalP2CL<double> lp2_p (lp1_p);
            LocalP1CL<double> lp1_n (*it, cn, Bndneg_);
            LocalP2CL<double> lp2_n (lp1_n);
            SMatrixCL<3,3> M;
            double det;
            GetTrafoTr(M,det,*it);
            double absdet= std::fabs(det);
            patch.Init( *it,  lsetp2new->Phi,0.);
            if (patch.Intersects()){
                patch.ComputeSubTets();
                Uint NumTets=patch.GetNumTetra(); //# of subtetras
                for (Uint k=0; k<NumTets; ++k)
                {
                    bool pPart = (k>=patch.GetNumNegTetra());
                    const SArrayCL<BaryCoordCL,4>& T =patch.GetTetra(k);
                    BaryCoordCL* nodes = Quad3CL<>::TransformNodes(T);
                    Quad3CL<> q3_dsol = Quad3CL<>(pPart? lp2_p : lp2_n, nodes);
                    Quad3CL<> q3_dsol2; q3_dsol2 = q3_dsol * q3_dsol;

                    double Vol = absdet*VolFrac(T);
                    if (pPart)
                        errl2p += q3_dsol2.quad(Vol);
                    else
                        errl2n += q3_dsol2.quad(Vol);
                    /* if (pPart) */
                    /*     errl1p += q3_diffabs.quad(Vol); */
                    /* else */
                    /*     errl1n += q3_diffabs.quad(Vol); */
                    delete[] nodes;
                }

            }
            else
            {
                bool pPart= (patch.GetSign( 0) == 1);
                Quad3CL<> q3_dsol = Quad3CL<>(pPart? lp2_p : lp2_n);
                Quad3CL<> q3_dsol2; q3_dsol2 = q3_dsol * q3_dsol;

                if (pPart)
                    errl2p += q3_dsol2.quad(absdet);
                else
                    errl2n += q3_dsol2.quad(absdet);
                /* if (pPart) */
                /*     errl1p += q3_diffabs.quad(absdet); */
                /* else */
                /*     errl1n += q3_diffabs.quad(absdet); */
                
            }
        }

        double errl2 = std::sqrt(errl2n + errl2p);
        /* double errl1 = (errl1n + errl1p); */
        std::cout << "errl2p = " << std::sqrt(errl2p) << "\t";
        std::cout << "errl2n = " << std::sqrt(errl2n) << "\t";
        std::cout << "errl2  = " << errl2 << std::endl;

        return errl2;
    }



};







/// \brief Observes the MultiGridCL-changes by AdapTriangCL to repair the Function ...
///
class MassTranspRepairCL : public MGObserverCL
{
  private:
    MultiGridCL & MG_;
    VecDescCL & c_;
    BndDataCL<> & cBnd_;
    std::unique_ptr<RepairP1CL<double>::type> p1repair_;
    const LevelsetP2CL& ls_;
  public:
    MassTranspRepairCL ( MultiGridCL & MG, VecDescCL & c, BndDataCL<> & cBnd, const LevelsetP2CL& ls)
        : MG_(MG), c_(c), cBnd_(cBnd), ls_( ls) {}
    void pre_refine  ();
    void post_refine ();
    void pre_refine_sequence  ();
    void post_refine_sequence ();
    const IdxDescCL* GetIdxDesc() const { return c_.RowIdx; }
#ifdef _PAR
    const VectorCL*  GetVector()  const { return &c_.Data; }
    void swap( IdxDescCL& idx, VectorCL& v) { c_.RowIdx->swap(idx); c_.Data.swap(v); }
#endif
};

} // end of namespace DROPS

#endif
