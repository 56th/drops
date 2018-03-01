/// \file gradient_recovery.cpp
/// \brief Compute a Lipschitz approximation of the gradient of a finite element function.
/// \author Joerg Grande

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
 * Copyright 2014, 2015, 2016 Joerg Grande, Aachen, Germany
*/

#include "num/gradient_recovery.h"
#include "geom/topo.h"
#include "geom/multigrid.h"
#include "num/bndData.h"
#include "num/fe.h"
#include "num/discretize.h"
#include "num/accumulator.h"
#include "num/oswald_projection.h"

namespace DROPS
{

// Compute the gradient of a P2 function in all P2-dofs.
class LocalP2GradientCL
{
  private:
    const VecDescCL& f_;
    const BndDataCL<>& fbnd_;

    SMatrixCL<3,3> M;
    LocalP1CL<Point3DCL> GradRefLP1[10],
                         GradLP1[10],
                         p1grad;
    LocalP2CL<> p2;
    LocalP2CL<Point3DCL> p2grad;

  public:
    typedef Point3DCL value_type;
    static const int num_components= 3;

    LocalP2GradientCL (const VecDescCL& f, const BndDataCL<>& fbnd)
        : f_( f), fbnd_( fbnd)  { P2DiscCL::GetGradientsOnRef( GradRefLP1); }

    void set_tetra (const TetraCL* t);
    value_type&       operator[] (size_t i)       { return p2grad[i]; }
    const value_type& operator[] (size_t i) const { return p2grad[i]; }
    bool invalid_p (size_t /*i*/) const { return false; }
    void finalize_accumulation () const {}
};

void LocalP2GradientCL::set_tetra (const TetraCL* t)
{
    double det; // dummy
    GetTrafoTr( M, det, *t);
    P2DiscCL::GetGradients( GradLP1, GradRefLP1, M);
    p2.assign( *t, f_, fbnd_);
    p1grad= Point3DCL();
    for (Uint i= 0; i < 10; ++i)
        p1grad+= p2[i]*GradLP1[i];
    for (Uint i= 0; i < 4; ++i)
        p2grad[i]= p1grad[i];
    for (Uint i= 0; i < 6; ++i)
        p2grad[i+4]= 0.5*(p1grad[VertOfEdge( i, 0)] + p1grad[VertOfEdge( i, 1)]);

}


void averaging_P2_gradient_recovery (const MultiGridCL& mg, const VecDescCL& f, const BndDataCL<>& fbnd, VecDescCL& grad)
{
    LocalP2GradientCL loc( f, fbnd);
    OswaldProjectionP2AccuCL<LocalP2GradientCL> accu( loc, grad);
    TetraAccumulatorTupleCL accus;
    accus.push_back( &accu);
    accumulate( accus, mg, f.RowIdx->TriangLevel(), f.RowIdx->GetBndInfo());
}

} // end of namespace DROPS

