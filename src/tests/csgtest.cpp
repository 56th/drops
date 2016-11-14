/// \file csgtest.cpp
/// \brief Test CSG for level sets. Test refinement strategies for the interface.
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
 * Copyright 2012 LNM/SC RWTH Aachen, Germany
*/

#include <fstream>

#include "geom/csg.h"
#include "misc/funcmap.h"
#include "misc/utils.h"
#include "geom/multigrid.h"
#include "geom/principallattice.h"
#include "geom/subtriangulation.h"
#include "levelset/surfacetension.h"
#include "levelset/levelset.h"
#include "geom/builder.h"
#include "num/fe.h"
#include "num/lattice-eval.h"
#include "misc/problem.h"
#include "out/vtkOut.h"


using namespace DROPS;

static void GetP2Hessian( SMatrixCL<3,3> G[10], const SMatrixCL<3,3>& M)
// G[d]= M*Href[d]*M^T
{
    for (Uint d= 0; d< 10; ++d) {
        std::memset( &G[d], 0, 3*3*sizeof( double));
        for (Uint i= 0; i < 3; ++i)
            for (Uint j= 0; j < 3; ++j)
                for (Uint k= 0; k < 3; ++k)
                    for (Uint l= 0; l < 3; ++l)
                        G[d](i,j)+= M(i, k)*M(j, l)*FE_P2CL::D2HRef( d, k, l);
    }
}

/// \brief Accumulate different error measures for the approximation of the level
/// set function $\varphi$ by the piecewise linear approximation $\varphi_h$.
class InterfaceApproxErrorAccuCL : public TetraAccumulatorCL
{
  private:
    VecDescCL* yH_, ///< error-term for H^2-smooth \varphi
             * yG_, ///< error-term for H^1 smooth \varphi
             * yQ_; ///< yH_/yG_ (clipped at 1e16)

    IdxT numry[4];
    double vec[4];

    const BaryCoordCL bc;
    SMatrixCL<3,3> M;
    double det;

    SMatrixCL<3,3> Hp2[10];
    LocalP1CL<Point3DCL> Grefp2[10],
                         Gp2[10];

    const PrincipalLatticeCL& lat;

    const VecDescCL*   ls;      // P2-level-set
    const BndDataCL<>* lsetbnd; // boundary data for the level set function
    LocalP2CL<> locp2_ls;
    std::valarray<double> ls_loc;

    bool setrefmarks_;

  public:
    InterfaceApproxErrorAccuCL (const LevelsetP2CL& lset,
        VecDescCL* yh, VecDescCL* yg, VecDescCL* yq, bool setrefmarks= false)
        : yH_( yh), yG_( yg), yQ_( yq), bc( 0.25), lat( PrincipalLatticeCL::instance( 2)), ls( &lset.Phi), lsetbnd( &lset.GetBndData()), ls_loc( 10), setrefmarks_( setrefmarks) { P2DiscCL::GetGradientsOnRef( Grefp2); }
    virtual ~InterfaceApproxErrorAccuCL () {}

    virtual void begin_accumulation () {
        std::cout << "#InterfaceApproxErrorAccuCL::begin_accumulation"
                     ": " << yH_->RowIdx->NumUnknowns() << " rows.\n";
    }

    virtual void finalize_accumulation() {}

    virtual void visit (const TetraCL& t) {
        locp2_ls.assign( t, *ls, *lsetbnd);
        evaluate_on_vertexes( locp2_ls, lat, Addr( ls_loc));
        if (equal_signs( ls_loc))
            return;

        GetTrafoTr( M, det, t);
        GetP2Hessian( Hp2, M);
        SMatrixCL<3,3> H;
        P2DiscCL::GetGradients( Gp2, Grefp2, M);
        Point3DCL G;
        for (Uint d= 0; d < 10; ++d) {
            H+= Hp2[d]*locp2_ls[d];
            G+= Gp2[d]( bc)*locp2_ls[d];
        }
        const double h= ::cbrt( std::abs( det));

        /// For $C^2$-regular $\varphi$, errQ$\in\mathcal{O}(h)$, $h\to 0$.
        /// One could declare the curvature a resolved for $errQ \lesssim 1$. This
        /// corresponds to the curvature being less than $h^{-1}$ (for
        /// distance-like $\varphi$ with $\|D\varphi\|\approx 1$).
        ///
        /// For kinks in $\varphi$, errQ will remain constant (around
        /// $[D\varphi]/{D\varphi}$). This is actually a feature as one can
        /// reliably detect kinks by $1\lesssim errQ$.
        const double errH= h*h*std::sqrt( trace( GramMatrix( H))),
                     errG= h*G.norm(),
                     errQ= errG < 1e-16 ? 1e16 : errH/errG;

        if (setrefmarks_ && errQ > 1.)
            const_cast<TetraCL&>( t).SetRegRefMark();

        GetLocalNumbP1NoBnd( numry, t, *yH_->RowIdx);
        for (int i= 0; i < 4; ++i) {
            yH_->Data[numry[i]]=  std::max( yH_->Data[numry[i]], errH);
            yG_->Data[numry[i]]= std::max( yG_->Data[numry[i]], errG);
            yQ_->Data[numry[i]]= std::max( yQ_->Data[numry[i]], errQ);
        }
    }

    virtual InterfaceApproxErrorAccuCL* clone (int /*clone_id*/) { return new InterfaceApproxErrorAccuCL( *this); }
};


/// The level set is a heart with an inward and an outward cusp on the p[2]-axis; contained in \f$(-1.2,1.2)\times(-0.7,0.7)\times(-1.1,1.3)\f$.
double suess (const Point3DCL& p, double)
{
    return std::pow( std::pow( p[0], 2) + 2.25*std::pow( p[1], 2) + std::pow( p[2], 2) - 1., 3) - std::pow( p[0], 2)*std::pow( p[2], 3) - 9./80.*std::pow( p[1], 2)*std::pow( p[2], 3);
}
RegisterScalarFunction reg_suess( "suess", &suess);

//dummy
double sigmaf (const Point3DCL&, double) { return 0.; }

const CSG::BodyCL* thebody;

inline double csg_fun (const Point3DCL& x, double t)
{
    return (*thebody)( x, t);
}

int TestExamples (MultiGridCL& mg)
{
    SurfaceTensionCL sf( sigmaf);   // dummy class
    LsetBndDataCL lsbnd( 6);
    LevelsetP2CL & lset( * LevelsetP2CL::Create( mg, lsbnd, sf) );

    lset.CreateNumbering( mg.GetLastLevel());

    IdxDescCL p1idx;
    p1idx.CreateNumbering( mg.GetLastLevel(), mg, lsbnd);
    VecDescCL ierr( &p1idx),
              ierrg( &p1idx),
              ierrq( &p1idx);

    std::ifstream jsonfile( "../../param/geom/unspecified/csg-examples.json");

    ParamCL p;
    jsonfile >> p;
    const size_t num= 2*std::distance( p.begin(), p.end());

   // writer for vtk-format
    VTKOutCL vtkwriter( mg, "DROPS data", num,
                        ".", "csg-examples",
                        "csg-examples", /* <- time file name */
                        true, /* <- binary */
                        false, /* <- onlyp1 */
                        false, /* <- p2dg */
                        -1, /* <- level */
                        0, 0);
    vtkwriter.Register( make_VTKScalar( lset.GetSolution(), "level-set") );
    vtkwriter.Register( make_VTKScalar( make_P1Eval( mg, lsbnd, ierr), "interpolation-errorH") );
    vtkwriter.Register( make_VTKScalar( make_P1Eval( mg, lsbnd, ierrg), "interpolation-errorG") );
    vtkwriter.Register( make_VTKScalar( make_P1Eval( mg, lsbnd, ierrq), "error-quotient") );

    // InterfaceApproxErrAccuCL accu( lset, &ierr, &ierrg, &ierrq);
    InterfaceApproxErrorAccuCL accu( lset, &ierr, &ierrg, &ierrq, true);
    TetraAccumulatorTupleCL accus;
    accus.push_back( &accu);

    size_t i= 0;
    for (ParamCL::ptree_const_iterator_type it= p.begin(); it != p.end(); ++it, ++i) {
        std::cout << "\n\n#Processing example \"" << it->first <<"\"." << std::endl;
        const bool is_cake= it->first == std::string( "Birthday cake (7 connected components)");

        thebody= CSG::body_builder( it->second);
        lset.Init( csg_fun);
        ierr.Data= 0.;
        ierrg.Data= 0.;
        ierrq.Data= 0.;
        accus( mg.GetTriangTetraBegin(), mg.GetTriangTetraEnd());
        vtkwriter.Write( 2*i);

        lset.idx.DeleteNumbering( mg);
        p1idx.DeleteNumbering( mg);
        mg.Refine();
        lset.CreateNumbering( mg.GetLastLevel());
        lset.Init( csg_fun);
        p1idx.CreateNumbering( mg.GetLastLevel(), mg, lsbnd);
        ierr.SetIdx( &p1idx);
        ierrg.SetIdx( &p1idx);
        ierrq.SetIdx( &p1idx);
        accus( mg.GetTriangTetraBegin(), mg.GetTriangTetraEnd());

        if (is_cake) { // needs an additional refinement to look good
            lset.idx.DeleteNumbering( mg);
            p1idx.DeleteNumbering( mg);
            mg.Refine();
            lset.CreateNumbering( mg.GetLastLevel());
            lset.Init( csg_fun);
            p1idx.CreateNumbering( mg.GetLastLevel(), mg, lsbnd);
            ierr.SetIdx( &p1idx);
            ierrg.SetIdx( &p1idx);
            ierrq.SetIdx( &p1idx);
            accus( mg.GetTriangTetraBegin(), mg.GetTriangTetraEnd());
        }

        vtkwriter.Write( 2*i + 1);

        lset.idx.DeleteNumbering( mg);
        p1idx.DeleteNumbering( mg);
        DROPS_FOR_TRIANG_TETRA( mg, mg.GetLastLevel(), it)
            it->SetRemoveMark(); // level 0 will not be removed.
        mg.Refine();
        if (is_cake) {
            DROPS_FOR_TRIANG_TETRA( mg, mg.GetLastLevel(), it)
                it->SetRemoveMark(); // level 0 will not be removed.
            mg.Refine();
        }
        lset.CreateNumbering( mg.GetLastLevel());
        lset.Init( csg_fun);
        p1idx.CreateNumbering( mg.GetLastLevel(), mg, lsbnd);
        ierr.SetIdx( &p1idx);
        ierrg.SetIdx( &p1idx);
        ierrq.SetIdx( &p1idx);

        // seq_out( Addr( ierr.Data), Addr( ierr.Data) + ierr.Data.size(), std::cout);
        delete thebody;
    }
    lset.idx.DeleteNumbering( mg);
    p1idx.DeleteNumbering( mg);
    std::cout << "Successfully proccessed '../geom/csg-examples.json'.\n";
    delete &lset;
    return 0;
}

int main ()
{
  try {
    BrickBuilderCL mgb( Point3DCL( -1.), 2*std_basis<3>( 1), 2*std_basis<3>( 2), 2*std_basis<3>( 3), 16, 16, 16);
    MultiGridCL mg( mgb);
    return  TestExamples( mg);
  }
  catch (DROPSErrCL err) { err.handle(); }
}
