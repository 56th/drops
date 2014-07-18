/// \file levelsetmapper.h
/// \brief Implements mappings from a neighborhood of the levelset to the level set.
/// \author LNM RWTH Aachen: Joerg Grande

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
 * Copyright 2014 LNM RWTH Aachen, Germany
*/

#ifndef DROPS_LEVELSETMAPPER_H
#define DROPS_LEVELSETMAPPER_H

#include "geom/subtriangulation.h"
#include "num/fe.h"
#include "num/discretize.h"
#include "misc/problem.h"

#include <tr1/unordered_map>
#include <tr1/unordered_set>

namespace DROPS
{

typedef std::tr1::unordered_set<const TetraCL*>            TetraSetT;
typedef std::tr1::unordered_map<const TetraCL*, TetraSetT> TetraToTetrasT;


class QuaQuaMapperCL
{
  private:
    int maxiter_;
    double tol_;
    bool use_line_search_;

    // The level set function.
    NoBndDataCL<> nobnddata;
    P2EvalCL<double, const NoBndDataCL<>, const VecDescCL> ls;

    // The recovered gradient of ls.
    NoBndDataCL<Point3DCL> nobnddata_vec;
    P2EvalCL<Point3DCL, const NoBndDataCL<Point3DCL>, const VecDescCL> ls_grad_rec;

    LocalP1CL<Point3DCL> gradrefp2[10];

    // The neighborhoods around each tetra in which base points are searched for.
    TetraToTetrasT& neighborhoods_;

    bool line_search (const Point3DCL& v, const Point3DCL& nx, const TetraCL*& tetra, BaryCoordCL& bary, const TetraSetT& neighborhood) const;
    void base_point_with_line_search (const TetraCL*& tet, BaryCoordCL& xb) const;
    void base_point_newton (const TetraCL*& tet, BaryCoordCL& xb) const;

  public:
    QuaQuaMapperCL (const MultiGridCL& mg, VecDescCL& lsarg, const VecDescCL& ls_grad_recarg, TetraToTetrasT& neigborhoods, int maxiter= 100, double tol= 1e-7, bool use_line_search= true)
        : maxiter_( maxiter), tol_( tol), use_line_search_( use_line_search),
          ls( &lsarg, &nobnddata, &mg), ls_grad_rec( &ls_grad_recarg, &nobnddata_vec, &mg), neighborhoods_( neigborhoods)
    { P2DiscCL::GetGradientsOnRef( gradrefp2); }


    void base_point (const TetraCL*& tet, BaryCoordCL& xb) const;
    void jacobian (const TetraCL& tet, const BaryCoordCL& xb, SMatrixCL<3,3>& dph) const;

    /// Return the local level set function and its gradient on tet; only for convenience. @{
    LocalP2CL<> local_ls      (const TetraCL& tet) const { return LocalP2CL<>( tet, ls); }
    Point3DCL   local_ls_grad (const TetraCL& tet, const BaryCoordCL& xb) const;
    ///@}
};

void compute_tetra_neighborhoods (const DROPS::MultiGridCL& mg, const VecDescCL& lsetPhi, const BndDataCL<>& lsetbnd, const PrincipalLatticeCL& lat, TetraToTetrasT& tetra_neighborhoods);

double abs_det (const TetraCL& tet, const BaryCoordCL& xb, const QuaQuaMapperCL& quaqua, const SurfacePatchCL& p);

inline bool is_in_ref_tetra (const BaryCoordCL& b, double eps= 1e-10)
{
    for (int i= 0; i < 4; ++i)
        if (b[i] < -eps || b[i] > 1. + eps)
            return false;
    return true;
}


} // end of namespace DROPS

#endif
