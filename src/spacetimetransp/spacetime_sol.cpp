/// \file spacetime_sol.cpp
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

#include "spacetimetransp/spacetime_sol.h"
#include "spacetimetransp/stxfem.h"
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

SpaceTimeXSolutionCL::SpaceTimeXSolutionCL(MultiGridCL& mg, const BndDataCL<>& Bndneg, const BndDataCL<>& Bndpos,
                                           ParamCL & P)
    :
    mg_(mg),
    Bndneg_(Bndneg),
    Bndpos_(Bndpos),
    stxfem_omit_param_(P.get<double>("Transp.XFEMReduced",0.0)),
    weight_pos_(P.get<double>("Transp.HPos",1.0)),
    weight_neg_(P.get<double>("Transp.HNeg",1.0)),
    st_x_idx_(P1SP1TX_FE,Bndneg,Bndpos,stxfem_omit_param_),
    fut_x_idx_(P1X_FE,Bndneg,Bndpos,0.0),
    past_x_idx_(P1X_FE,Bndneg,Bndpos,0.0),
    // neg_space_idx_(P1_FE,Bndneg,0.0),
    // pos_space_idx_(P1_FE,Bndpos,0.0),
    lsetp2old(0), lsetp2new(0),
    p1idxn(P1_FE),
    p1idxp(P1_FE),
    future_vec_pos(&p1idxp),
    future_vec_neg(&p1idxn),
    past_vec_pos(&p1idxp),
    past_vec_neg(&p1idxn)
{
    sol.SetIdx( &st_x_idx_);
    futurevec.SetIdx( &fut_x_idx_);
    pastvec.SetIdx( &past_x_idx_);

}


void SpaceTimeXSolutionCL::UpdateTimeSlab(const LevelsetP2CL & lsetold, const LevelsetP2CL & lsetnew,
                                          double vmax, instat_scalar_fun_ptr lset_fpt)
{
    lsetp2old = &lsetold;
    lsetp2new = &lsetnew;
    if (st_x_idx_.NumUnknowns() != 0)
        st_x_idx_.DeleteNumbering( mg_);
    st_x_idx_.CreateNumbering( mg_.GetLastLevel(), mg_, Bndneg_, Bndpos_, &lsetold.Phi, &lsetold.GetBndData());
    const double told = lsetold.Phi.t;
    const double tnew = lsetnew.Phi.t;
    STXFEM::UpdateSTXNumbering( &st_x_idx_, mg_, lsetold.Phi, lsetnew.Phi, lset_fpt, TimeInterval(told,tnew), lsetold.GetBndData(), true, vmax);
    sol.SetIdx( &st_x_idx_);
}

void SpaceTimeXSolutionCL::EvalFutureTrace()
{
    if (fut_x_idx_.NumUnknowns() != 0)
        fut_x_idx_.DeleteNumbering( mg_);
    fut_x_idx_.CreateNumbering( mg_.GetLastLevel(), mg_, Bndneg_, Bndpos_, &(lsetp2new->Phi), &(lsetp2new->GetBndData()));
    futurevec.SetIdx(&fut_x_idx_);
    STXFEM::GetFutureTrace(sol, futurevec, mg_);
    futurevec.t = lsetp2new->Phi.t;

    p1idxn.CreateNumbering( mg_.GetLastLevel(), mg_, fut_x_idx_);
    p1idxp.CreateNumbering( mg_.GetLastLevel(), mg_, fut_x_idx_);
    future_vec_neg.RowIdx = &p1idxn;
    future_vec_pos.RowIdx = &p1idxp;

    P1XtoP1 ( fut_x_idx_, futurevec.Data, p1idxn, future_vec_pos.Data, 
              future_vec_neg.Data, lsetp2new->Phi, mg_);
    future_vec_neg.Data *= 1.0/weight_neg_;
    future_vec_pos.Data *= 1.0/weight_pos_;
    future_vec_neg.t = lsetp2new->Phi.t;
    future_vec_pos.t = lsetp2new->Phi.t;
    

}

void SpaceTimeXSolutionCL::EvalPastTrace()
{
    if (past_x_idx_.NumUnknowns() != 0)
        past_x_idx_.DeleteNumbering( mg_);
    past_x_idx_.CreateNumbering( mg_.GetLastLevel(), mg_, Bndneg_, Bndpos_, &(lsetp2old->Phi), &(lsetp2old->GetBndData()));
    pastvec.SetIdx(&past_x_idx_);
    STXFEM::GetPastTrace(sol, pastvec, mg_);
    pastvec.t = lsetp2old->Phi.t;

    p1idxn.CreateNumbering( mg_.GetLastLevel(), mg_, past_x_idx_);
    p1idxp.CreateNumbering( mg_.GetLastLevel(), mg_, past_x_idx_);
    P1XtoP1 ( past_x_idx_, pastvec.Data, p1idxn, past_vec_pos.Data, 
              past_vec_neg.Data, lsetp2old->Phi, mg_);
    past_vec_neg.Data *= 1.0/weight_neg_;
    past_vec_pos.Data *= 1.0/weight_pos_;
    past_vec_neg.t = lsetp2old->Phi.t;
    past_vec_pos.t = lsetp2old->Phi.t;

}

SpaceTimeXSolutionCL::~SpaceTimeXSolutionCL()
{
}


//*****************************************************************************
//                               ConcentrationRepairCL
//*****************************************************************************

inline void MassTranspRepairCL::pre_refine()
{
    // std::cout << " before " << std::endl;
    // std::cout << " c_.Data = " << c_.Data << std::endl;
    p1repair_= std::unique_ptr<RepairP1CL<double>::type >(
        new RepairP1CL<double>::type( MG_, c_, cBnd_));
}

inline void
  MassTranspRepairCL::post_refine ()
{
    VecDescCL loc_c;
    IdxDescCL loc_cidx( P1_FE);

    loc_cidx.CreateNumbering( MG_.GetLastLevel(), MG_, cBnd_, &ls_.Phi, &ls_.GetBndData());
    loc_c.SetIdx( &loc_cidx);

    p1repair_->repair( loc_c);
    c_.Clear( c_.t);
    c_.RowIdx->DeleteNumbering( MG_);
    c_.RowIdx->swap( loc_cidx);
    c_.SetIdx( c_.RowIdx);
    // loc_c.Data = 1.0;
    c_.Data= loc_c.Data;

    // std::cout << " after " << std::endl;
    // std::cout << " c_.Data = " << c_.Data << std::endl;

}

inline void
  MassTranspRepairCL::pre_refine_sequence ()
{
}

inline void
  MassTranspRepairCL::post_refine_sequence ()
{
}


} // end of namespace DROPS
