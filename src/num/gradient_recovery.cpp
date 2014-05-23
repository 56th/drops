/// \file gradient_recovery.cpp
/// \brief Compute a Lipschitz approximation of the gradient of a finite element function.
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

#include "num/gradient_recovery.h"
#include "geom/topo.h"
#include "geom/multigrid.h"
#include "num/bndData.h"
#include "num/fe.h"
#include "num/discretize.h"
#include "num/accumulator.h"

namespace DROPS
{

class AveragingP2GradientRecoveryAccuCL : public TetraAccumulatorCL
{
  private:
    std::valarray<double>* n_;
    const VecDescCL& f_;
    const BndDataCL<>& fbnd_;
    VecDescCL& grad_;

    SMatrixCL<3,3> M;
    LocalP1CL<Point3DCL> GradRefLP1[10],
                         GradLP1[10],
                         p1grad;
    LocalP2CL<> p2;
    LocalP2CL<Point3DCL> p2grad;

    LocalNumbP2CL numg;

    void set_n (std::valarray<double>* n) { n_= n; } // The clones must refer to the n_ of thread 0.

  public:
    AveragingP2GradientRecoveryAccuCL (const VecDescCL& f, const BndDataCL<>& fbnd, VecDescCL& grad)
        : n_( 0), f_( f), fbnd_( fbnd), grad_( grad) { P2DiscCL::GetGradientsOnRef( GradRefLP1); }

    virtual void begin_accumulation   () {
        n_= new std::valarray<double>( grad_.Data.size()/3);
    }
    virtual void finalize_accumulation() {
        delete n_;
    }

    virtual void visit (const TetraCL& t) {
        double det; // dummy
        GetTrafoTr( M, det, t);
        P2DiscCL::GetGradients( GradLP1, GradRefLP1, M);
        p2.assign( t, f_, fbnd_);
        p1grad= Point3DCL();
        for (Uint i= 0; i < 10; ++i)
            p1grad+= p2[i]*GradLP1[i];
        for (Uint i= 0; i < 4; ++i)
            p2grad[i]= p1grad[i];
        for (Uint i= 0; i < 6; ++i)
            p2grad[i+4]= 0.5*(p1grad[VertOfEdge( i, 0)] + p1grad[VertOfEdge( i, 1)]);
        numg.assign_indices_only( t, *grad_.RowIdx);
        for (Uint i= 0; i < 10; ++i) {
            if (!numg.WithUnknowns( i))
                continue;
            const IdxT dof= numg.num[i];
            double& n= n_[0][dof/3]; // This assumes that the local gradient is in the components dof/3..dof/3 + 2.
            n+= 1.;
            const Point3DCL& oldgrad= DoFHelperCL<Point3DCL, VectorCL>::get( grad_.Data, dof);
            DoFHelperCL<Point3DCL, VectorCL>::set( grad_.Data, dof, ((n - 1.)/n)*oldgrad + (1./n)*p2grad[i]);
        }
    }

    virtual TetraAccumulatorCL* clone (int /*clone_id*/) {
        AveragingP2GradientRecoveryAccuCL* p= new AveragingP2GradientRecoveryAccuCL( f_, fbnd_, grad_);
        p->set_n( n_);
        return p;
    }
};

void averaging_P2_gradient_recovery (const MultiGridCL& mg, const VecDescCL& f, const BndDataCL<>& fbnd, VecDescCL& grad)
{
    AveragingP2GradientRecoveryAccuCL accu( f, fbnd, grad);
    TetraAccumulatorTupleCL accus;
    accus.push_back( &accu);
    accumulate( accus, mg, f.RowIdx->TriangLevel(), f.RowIdx->GetMatchingFunction(), f.RowIdx->GetBndInfo());
}

} // end of namespace DROPS

