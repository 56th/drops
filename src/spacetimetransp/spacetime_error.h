/// \file spacetime_error.h
/// \brief setup routines for calculating errors
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

#ifndef DROPS_SPACETIME_ERROR_H
#define DROPS_SPACETIME_ERROR_H

#include "misc/problem.h"
#include "misc/params.h"
#include "num/spacetime_geom.h"
#include "num/spacetime_quad.h"
#include "spacetimetransp/stxfem.h"
#include "spacetimetransp/spacetime_setup.h"
#include "num/interfacePatch.h"
#include "num/discretize.h"
#include "num/accumulator.h"
#include "num/quadrature.h"
#include "levelset/levelset.h"

#include <iostream>
#include <fstream>

namespace DROPS
{


// Accumulator for space time integrals - measure error in jump [ beta  u ], different versions..
class EnergyNormErrorAccumulatorCL : public TetraAccumulatorCL
{
    protected:
    const MultiGridCL& MG_;

    const LevelsetP2CL * lsetp2old;
    const LevelsetP2CL * lsetp2new;
    instat_scalar_fun_ptr lset_fpt;

    instat_scalar_fun_ptr sol_neg;
    instat_scalar_fun_ptr sol_pos;

    instat_scalar_fun_ptr sol_dt_neg;
    instat_scalar_fun_ptr sol_dt_pos;

    instat_vector_fun_ptr sol_grad_neg;
    instat_vector_fun_ptr sol_grad_pos;

    double det;
    double absdet;
    Point3DCL G[4];
    SMatrixCL<3,3> T;
    bool sign[8]; //sign of levelset at vertices (4(space) x 2(time))

    const Uint ints_per_space_edge;
    const Uint subtimeintervals;

    const double told;
    const double tnew;

    double iscut;
    double * h10norm2; // ||.||_{1,0}^2
    double * h01norm2; // ||.||_{0,1}^2
    double * l2norm2;  // ||.||_{0,0}^2

    LocalP2CL<double> phi[4];
    Quad5CL<double> q_phi[4];
    LocalP2CL<double> phi_i_phi_j[4][4];

    const VecDescCL & old_sol_neg;
    const VecDescCL & old_sol_pos;
    const VecDescCL & new_sol_neg;
    const VecDescCL & new_sol_pos;
    const BndDataCL<> & Bnd_neg;
    const BndDataCL<> & Bnd_pos;

    const double beta_neg;
    const double beta_pos;
    const double alpha_neg;
    const double alpha_pos;
    const double lambda_stab;
    const double vmax;

public:
    EnergyNormErrorAccumulatorCL (const MultiGridCL& MG, const LevelsetP2CL * lsetp2old_in,
                                const LevelsetP2CL * lsetp2new_in,
                                instat_scalar_fun_ptr lset_fpt, const double t1, const double t2,
                                const VecDescCL & oldsol_neg_in, const VecDescCL & oldsol_pos_in,
                                const VecDescCL & newsol_neg_in, const VecDescCL & newsol_pos_in,
                                const BndDataCL<> & Bnd_neg_in,
                                const BndDataCL<> & Bnd_pos_in,
                                const ParamCL::ptree_type & P);
    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();
    virtual void visit (const TetraCL&);
    virtual TetraAccumulatorCL* clone (int /*tid*/) { return new EnergyNormErrorAccumulatorCL( *this); }
};


}//end of namespace DROPS

#endif
