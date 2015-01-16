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
#include "misc/omp_variable.h"
#include "geom/multigrid.h"
#include "geom/principallattice.h"
#include "geom/subtriangulation.h"
#include "levelset/levelset.h"
#include "levelset/adaptriang.h"
#include "levelset/levelset.h"
#include "levelset/levelsetmapper.h"
#include "levelset/marking_strategy.h"
#include "geom/builder.h"
#include "num/fe.h"
#include "num/gradient_recovery.h"
#include "num/lattice-eval.h"
#include "misc/problem.h"
#include "misc/auto_diff.h"
#include "out/vtkOut.h"
#include "surfactant/ifacetransp.h"


using namespace DROPS;

ParamCL P;

/// \brief Accumulate different error measures for the approximation of the level
/// set function $\varphi$ by the piecewise linear approximation $\varphi_h$.
class InterfaceApproxErrorAccuCL : public TetraAccumulatorCL
{
  private:
    VecDescCL* yH_, ///< error-term for H^2-smooth \varphi
             * yG_, ///< error-term for H^1 smooth \varphi
             * yQ_, ///< yH_/yG_ (clipped at 1e16)
             * ydist_; ///< distance per QuaQuamapper

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

    QuaQuaMapperCL* mapper_;

  public:
    InterfaceApproxErrorAccuCL (const LevelsetP2CL& lset,
        VecDescCL* yh, VecDescCL* yg, VecDescCL* yq, VecDescCL* ydist, bool setrefmarks= false)
        : yH_( yh), yG_( yg), yQ_( yq), ydist_( ydist), bc( 0.25), lat( PrincipalLatticeCL::instance( 2)), ls( &lset.Phi), lsetbnd( &lset.GetBndData()), ls_loc( 10), setrefmarks_( setrefmarks), mapper_( 0) { P2DiscCL::GetGradientsOnRef( Grefp2); }
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
        P2DiscCL::GetHessians( Hp2, M);
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
        SurfacePatchCL p;
        if (mapper_) {
            BaryCoordCL xb;
            p.make_patch<MergeCutPolicyCL>( lat, ls_loc);
            for (Uint tri= 0; tri < p.facet_size(); ++tri) {
                const SurfacePatchCL::FacetT& f= p.facet_begin()[tri];
                xb= (p.vertex_begin()[f[0]] + p.vertex_begin()[f[1]] + p.vertex_begin()[f[2]])/3.;
                mapper_->set_point( &t, xb)
                        .base_point();
            }

//             for (Uint i= 0; i < 4; ++i) {
//                 if (ydist_->Data[numry[i]] != 0.)
//                     continue;
//                 mapper_->set_point( &t,std_basis<4>( i + 1))
//                            .base_point();
//                 ydist_->Data[numry[i]]= mapper_->get_dh();
//             }
        }
        for (int i= 0; i < 4; ++i) {
//             yH_->Data[numry[i]]=  std::max( yH_->Data[numry[i]], errH);
//             yG_->Data[numry[i]]= std::max( yG_->Data[numry[i]], errG);
//             yQ_->Data[numry[i]]= std::max( yQ_->Data[numry[i]], errQ);
            yH_->Data[numry[i]]=  std::max( yH_->Data[numry[i]], errH/h);
            yG_->Data[numry[i]]= std::max( yG_->Data[numry[i]], errG/h);
            yQ_->Data[numry[i]]= std::max( yQ_->Data[numry[i]], errQ);
        }
    }

    virtual InterfaceApproxErrorAccuCL* clone (int /*clone_id*/) { return new InterfaceApproxErrorAccuCL( *this); }

    void set_mapper (QuaQuaMapperCL* mapper) { mapper_= mapper; }
};

/// \brief Accumulate L2-norms and errors on the higher order zero level.
/// Works for P1IF_FE, P2IF_FE, and C-functions. All functions are evaluated on the P2-levelset.
class InterfaceL2AccuP2CL : public TetraAccumulatorCL
{
  private:
    const InterfaceCommonDataP2CL& cdata_;
    const MultiGridCL& mg;
    std::string name_;

    instat_scalar_fun_ptr f;
    double f_time;

    bool sample_curvature_;

    VecDescCL* redist_err_vd;

  public:
    OpenMPVar_ZeroInit_Sum_CL<double> f_norm,
           f_int,
           err,
           area,
           areap1,
           time_dh,
           time_locate,
           iter_dh,
           inner_iter_dh,
           num_damping,
           invocations_dh,
           calls_locate;
    OpenMPVar_MinInit_Max_CL<double> max_curv,
                                     max_iter,
                                     max_grad;
    OpenMPVar_MaxInit_Min_CL<double> min_iter,
                                     min_grad;

    InterfaceL2AccuP2CL (const InterfaceCommonDataP2CL& cdata, const MultiGridCL& mg_arg, std::string name= std::string())
        : cdata_( cdata), mg( mg_arg), name_( name), f( 0), f_time( 0.), sample_curvature_( true), redist_err_vd( 0) {}
    virtual ~InterfaceL2AccuP2CL () {}

    void set_name (const std::string& n) { name_= n; }
    void set_function (const instat_scalar_fun_ptr farg, double f_time_arg= 0.) {
        f= farg;
        f_time= f_time_arg;
    }

    void set_sample_curvature (bool b) { sample_curvature_= b; }
    void compute_redist_err (VecDescCL* vd) { redist_err_vd= vd; }

    virtual void begin_accumulation () {
        std::cout << "InterfaceL2AccuP2CL::begin_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\".\n";

        if (redist_err_vd)
            redist_err_vd->Data= 128.; // std::numeric_limits<double>::max();

        f_norm.scatter();
        f_int.scatter();
        err.scatter();
        area.scatter();
        areap1.scatter();
        time_dh.scatter();
        time_locate.scatter();
        iter_dh.scatter();
        inner_iter_dh.scatter();
        num_damping.scatter();
        invocations_dh.scatter();
        calls_locate.scatter();
        max_curv.scatter();
        max_grad.scatter();
        min_grad.scatter();
        max_iter.scatter();
        min_iter.scatter();
    }

    virtual void finalize_accumulation() {
        std::cout << "InterfaceL2AccuP2CL::finalize_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\":";

        area.reduce();
        std::cout.precision( 16);
        std::cout << "\n\tarea: " << area.value();
        areap1.reduce();
        std::cout << "\n\tareap1: " << areap1.value();
        time_dh.reduce();
        std::cout << "\n\ttime_dh: " << time_dh.value();
        time_locate.reduce();
        std::cout << "\n\ttime_locate: " << time_locate.value();
        iter_dh.reduce();
        std::cout << "\n\titer_dh: " << iter_dh.value();
        inner_iter_dh.reduce();
        std::cout << "\n\tinner_iter_dh: " << inner_iter_dh.value() << "\t average_per_outer_iter: " << inner_iter_dh.value()/iter_dh.value();
        num_damping.reduce();
        std::cout << "\n\tnum_damping: " << num_damping.value();
        invocations_dh.reduce();
        std::cout << "\n\tinvocations_dh: " << invocations_dh.value() << "\t average_iter: " << iter_dh.value()/invocations_dh.value();
        calls_locate.reduce();
        std::cout << "\n\tcalls_locate: " << calls_locate.value() << "\t average: " << calls_locate.value()/invocations_dh.value();
        max_curv.reduce();
        std::cout << "\n\tmax_curv: " << max_curv.value();
        max_grad.reduce();
        std::cout << "\n\tmax_grad: " << max_grad.value();
        min_grad.reduce();
        std::cout << "\n\tmin_grad: " << min_grad.value();
        max_iter.reduce();
        std::cout << "\n\tmax_iter: " << max_iter.value();
        if (redist_err_vd) {
            double tmp= 0.;
            for (size_t i= 0; i < redist_err_vd->Data.size(); ++i)
                if (redist_err_vd->Data[i] != 128. /*std::numeric_limits<double>::max()*/ && std::abs( redist_err_vd->Data[i]) > tmp)
                    tmp= std::abs( redist_err_vd->Data[i]);
            std::cout << "\n\tredist_err: " << tmp;
        }
        min_iter.reduce();
        std::cout << "\n\tmin_iter: " << min_iter.value();
//         if (f != 0) {
//             f_norm.reduce();
//             f_norm.value()= std::sqrt( f_norm.value());
//             f_int.reduce();
//             std::cout << "\n\t|| f ||_L2: " << f_norm.value()
//                       << "\t integral: " << f_int.value();
//         }
        std::cout << std::endl;
    }

    virtual void visit (const TetraCL& t) {
        const InterfaceCommonDataP2CL& cdata= cdata_.get_clone();
        if (cdata.empty())
            return;

        const int tid= omp_get_thread_num();

        area.value( tid)+= quad_2D( cdata.qdom_projected.absdets(), cdata.qdom);
        areap1.value( tid)+= quad_2D( std::valarray<double>( 1., cdata.qdom.vertex_size()), cdata.qdom);

        std::valarray<double> qf;
//         if (f != 0) {
//             resize_and_evaluate_on_vertexes( f, cdata.qdom_projected.vertexes(), f_time, qf);
//             f_int.value( tid)+= quad_2D( cdata.qdom_projected.absdets()*qf, cdata.qdom);
//             f_norm.value( tid)+= quad_2D( cdata.qdom_projected.absdets()*qf*qf, cdata.qdom);
//         }

        // Sample the maximal curvature radius in all quadrature points in t.
        if (sample_curvature_) {
            SMatrixCL<3,3> M,
                           P;
            double det;
            GetTrafoTr( M, det, t);
            SMatrixCL<3,3> Hp2[10];
            LocalP1CL<Point3DCL> Gp2[10];
            P2DiscCL::GetHessians( Hp2, M);
            SMatrixCL<3,3> H;
            P2DiscCL::GetGradients( Gp2, cdata.gradrefp2, M);
            Point3DCL G;

            for (Uint d= 0; d < 10; ++d) {
                H+= Hp2[d]*cdata.locp2_ls[d];
            }
            for (Uint i= 0; i < cdata.qdom_projected.vertexes().size(); ++ i) {
                const TetraBaryPairT& tb= cdata.qdom_projected.vertexes()[i];
                if (tb.first != &t)
                    continue;
                G= Point3DCL();
                for (Uint d= 0; d < 10; ++d)
                    G+= Gp2[d]( tb.second)*cdata.locp2_ls[d];
                P= eye<3,3>() - outer_product(G/G.norm_sq(), G);
                M= P*H*P/G.norm();
                cyclic_jacobi( M, 1e-6);
                for (Uint j= 0; j < 3; ++j)
                    max_curv.value( tid)= std::max( max_curv.value( tid), M(j, j));
                max_grad.value( tid)= std::max( max_grad.value( tid), G.norm());
                min_grad.value( tid)= std::min( min_grad.value( tid), G.norm());
            }
        }

        if (redist_err_vd) {
            for (Uint i= 0; i < 10; ++i) {
                const IdxT idx= i < 4 ? t.GetVertex( i)->  Unknowns( redist_err_vd->RowIdx->GetIdx())
                                      : t.GetEdge( i - 4)->Unknowns( redist_err_vd->RowIdx->GetIdx());
                if (redist_err_vd->Data[idx] != 128. /*std::numeric_limits<double>::max()*/)
                    continue;
                const BaryCoordCL& b= i < 4 ? std_basis<4>( i + 1)
                                            : BaryCenter( std_basis<4>( VertOfEdge( i - 4, 0) + 1), std_basis<4>( VertOfEdge( i - 4, 1) + 1));
                const Point3DCL& p= i < 4 ? t.GetVertex( i)->GetCoord()
                                          : GetBaryCenter( *t.GetEdge( i - 4));
                cdata.quaqua.set_point( &t, b);
                redist_err_vd->Data[idx]= std::abs( cdata.quaqua.base_point().get_dh() - f( p, 0.));
            }
        }

        // Update QuaQua statistics
        time_dh.value( tid)= cdata.quaqua.base_point_time;
        time_locate.value( tid)= cdata.quaqua.locate_new_point_time;
        iter_dh.value( tid)= cdata.quaqua.total_outer_iter;
        inner_iter_dh.value( tid)= cdata.quaqua.total_inner_iter;
        num_damping.value( tid)= cdata.quaqua.total_damping_iter;
        invocations_dh.value( tid)= cdata.quaqua.total_base_point_calls;
        calls_locate.value( tid)= cdata.quaqua.total_locate_new_point_calls;
        min_iter.value( tid)= cdata.quaqua.min_outer_iter;
        max_iter.value( tid)= cdata.quaqua.max_outer_iter;

    }

    virtual InterfaceL2AccuP2CL* clone (int /*clone_id*/) {
        InterfaceL2AccuP2CL* tmp= new InterfaceL2AccuP2CL( *this);
        tmp->f_norm.make_reference_to( f_norm);
        tmp->f_int.make_reference_to( f_int);
        tmp->err.make_reference_to( err);
        tmp->area.make_reference_to( area);
        tmp->areap1.make_reference_to( areap1);
        tmp->time_dh.make_reference_to( time_dh);
        tmp->time_locate.make_reference_to( time_locate);
        tmp->iter_dh.make_reference_to( iter_dh);
        tmp->inner_iter_dh.make_reference_to( inner_iter_dh);
        tmp->num_damping.make_reference_to( num_damping);
        tmp->invocations_dh.make_reference_to( invocations_dh);
        tmp->calls_locate.make_reference_to( calls_locate);
        tmp->max_curv.make_reference_to( max_curv);
        tmp->max_grad.make_reference_to( max_grad);
        tmp->min_grad.make_reference_to( min_grad);
        tmp->max_iter.make_reference_to( max_iter);
        tmp->min_iter.make_reference_to( min_iter);
        return tmp;
    }
};


/// The level set is a heart with an inward and an outward cusp on the p[2]-axis; contained in \f$(-1.2,1.2)\times(-0.7,0.7)\times(-1.1,1.3)\f$.
double suess (const Point3DCL& p, double)
{
    return std::pow( std::pow( p[0], 2) + 2.25*std::pow( p[1], 2) + std::pow( p[2], 2) - 1., 3) - std::pow( p[0], 2)*std::pow( p[2], 3) - 9./80.*std::pow( p[1], 2)*std::pow( p[2], 3);
}
RegisterScalarFunction reg_suess( "suess", &suess);


double deco_cube_radius;
double deco_cube_shift;

void deco_cube_val_grad (const Point3DCL& pp, double& v, Point3DCL& g)
{
    typedef ADFwdCL<double> Adf;
    const Adf cc( deco_cube_radius*deco_cube_radius),
              shift( deco_cube_shift),
              one( 1.);
    Adf p[3],
        gi;
    for (Uint i= 0; i < 3; ++i) {
        for (Uint j= 0; j < 3; ++j)
            p[j]= Adf( pp[j]);
        p[i].seed();
        gi= (pow( p[0]*p[0] + p[1]*p[1] - cc, 2) + pow( p[2]*p[2] - one, 2))
           *(pow( p[1]*p[1] + p[2]*p[2] - cc, 2) + pow( p[0]*p[0] - one, 2))
           *(pow( p[2]*p[2] + p[0]*p[0] - cc, 2) + pow( p[1]*p[1] - one, 2))
           + shift;
        g[i]= gi.derivative();
    }
    v= gi.value();
}

double deco_cube (const Point3DCL& p, double)
{
    const double cc= deco_cube_radius*deco_cube_radius;
    return  (std::pow( p[0]*p[0] + p[1]*p[1] - cc, 2) + std::pow( p[2]*p[2] - 1, 2))
           *(std::pow( p[1]*p[1] + p[2]*p[2] - cc, 2) + std::pow( p[0]*p[0] - 1, 2))
           *(std::pow( p[2]*p[2] + p[0]*p[0] - cc, 2) + std::pow( p[1]*p[1] - 1, 2)) + deco_cube_shift;
//     const double cc= deco_cube_radius*deco_cube_radius;
//     return  (std::abs( p[0]*p[0] + p[1]*p[1] - cc) + std::abs( p[2]*p[2] - 1))
//            *(std::abs( p[1]*p[1] + p[2]*p[2] - cc) + std::abs( p[0]*p[0] - 1))
//            *(std::abs( p[2]*p[2] + p[0]*p[0] - cc) + std::abs( p[1]*p[1] - 1)) + deco_cube_shift;
//     double v;
//     Point3DCL g;
//     deco_cube_val_grad( p, v, g);
// //     std::cout << "p: " << p << "\t\t " << v << "\t " << g.norm() << "\t\t " << v/g.norm() << std::endl;
// //     return v/std::sqrt( 10. + g.norm_sq());
//     return v;
}
RegisterScalarFunction reg_deco_cube( "deco_cube", &deco_cube);

const CSG::BodyCL* torus_body; // The (unperturbed) distance function
double torus_distance (const Point3DCL& p, double)
{
    return torus_body[0]( p, 0.);
}
double perturbed_torus_alpha; // perturbation parameter
double perturbed_torus (const Point3DCL& p, double)
{
    return torus_body[0]( p, 0.)*(9. + 4.*std::cos( perturbed_torus_alpha*p[0]*p[1]/p.norm()));
}
RegisterScalarFunction reg_perturbed_torus( "perturbed_torus", &perturbed_torus);

class LevelsetReinitCL : public MGObserverCL
{
  private:
    LevelsetP2CL& ls_;
    instat_scalar_fun_ptr f_;

  public:
    LevelsetReinitCL (LevelsetP2CL& ls, instat_scalar_fun_ptr f)
        : ls_( ls), f_( f) {}

    void pre_refine  () {};
    void post_refine ();

    void pre_refine_sequence  () {}
    void post_refine_sequence () {}
    const IdxDescCL* GetIdxDesc() const { return ls_.Phi.RowIdx; }
};

void
LevelsetReinitCL::post_refine ()
{
    ls_.DeleteNumbering( &ls_.idx);
    ls_.CreateNumbering( ls_.GetMG().GetLastLevel(), &ls_.idx);
    ls_.Phi.SetIdx( &ls_.idx);
    ls_.Init( f_);
}

//dummy
double sigmaf (const Point3DCL&, double) { return 0.; }

const CSG::BodyCL* thebody;

inline double csg_fun (const Point3DCL& x, double t)
{
    return (*thebody)( x, t);
}


int TestExamples (MultiGridCL& mg, ParamCL& p)
{
    SurfaceTensionCL sf( sigmaf);   // dummy class
    LsetBndDataCL lsbnd( 6);
    std::auto_ptr<LevelsetP2CL> lset_ptr( LevelsetP2CL::Create( mg, lsbnd, sf));
    LevelsetP2CL& lset= *lset_ptr;

    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);

    IdxDescCL p1idx;
    p1idx.CreateNumbering( mg.GetLastLevel(), mg, lsbnd);
    IdxDescCL vecp2idx( vecP2_FE);
    NoBndDataCL<Point3DCL> vecp2bnd;
    vecp2idx.CreateNumbering( mg.GetLastLevel(), mg, vecp2bnd);
    VecDescCL ierr( &p1idx),
              ierrg( &p1idx),
              ierrq( &p1idx),
              idist( &p1idx),
              lsgradrec( &vecp2idx);

    const size_t num= 2*std::distance( p.begin(), p.end());

   // writer for vtk-format
    VTKOutCL vtkwriter( mg, "DROPS data", num,
                        P.get<std::string>( "VTK.VTKDir"),
                        P.get<std::string>( "VTK.VTKName"),
                        P.get<std::string>( "VTK.TimeFileName"),
                        P.get<bool>( "VTK.Binary"),
                        P.get<bool>( "VTK.UseOnlyP1"), /* <- onlyp1 */
                        false, /* <- p2dg */
                        -1, /* <- level */
                        false, /* <- reusepvd */
                        false /* <- usedeformed */
              );
    vtkwriter.Register( make_VTKScalar( lset.GetSolution(), "level-set") );
    vtkwriter.Register( make_VTKScalar( make_P1Eval( mg, lsbnd, ierr), "interpolation-errorH") );
    vtkwriter.Register( make_VTKScalar( make_P1Eval( mg, lsbnd, ierrg), "interpolation-errorG") );
    vtkwriter.Register( make_VTKScalar( make_P1Eval( mg, lsbnd, ierrq), "error-quotient") );
    vtkwriter.Register( make_VTKScalar( make_P1Eval( mg, lsbnd, idist), "dh") );
    vtkwriter.Register( make_VTKVector( make_P2Eval( mg, vecp2bnd, lsgradrec), "ls_grad_rec") );

    // InterfaceApproxErrAccuCL accu( lset, &ierr, &ierrg, &ierrq);
    InterfaceApproxErrorAccuCL accu( lset, &ierr, &ierrg, &ierrq, &idist, true);
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
        idist.Data= 0.;
//         averaging_P2_gradient_recovery( mg, lset.Phi, lset.GetBndData(), lsgradrec);
//         QuaQuaMapperCL quaqua( mg, lset.Phi, lsgradrec, /*neighborhoods*/ 0, /*maxiter*/ 100, /*tol*/ 1e-7, /*use_line_search*/ false, /*armijo_c*/ 1e-2);
//         accu.set_mapper( &quaqua);
        accus( mg.GetTriangTetraBegin(), mg.GetTriangTetraEnd());
//         accu.set_mapper( 0);

        vtkwriter.Write( 2*i);

        lset.idx.DeleteNumbering( mg);
        p1idx.DeleteNumbering( mg);
        vecp2idx.DeleteNumbering( mg);
        mg.Refine();
        lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
        lset.Phi.SetIdx( &lset.idx);
        lset.Init( csg_fun);
        p1idx.CreateNumbering( mg.GetLastLevel(), mg, lsbnd);
        ierr.SetIdx( &p1idx);
        ierrg.SetIdx( &p1idx);
        ierrq.SetIdx( &p1idx);
        idist.SetIdx( &p1idx);
        vecp2idx.CreateNumbering( mg.GetLastLevel(), mg, vecp2bnd);
        lsgradrec.SetIdx( &vecp2idx);
        averaging_P2_gradient_recovery( mg, lset.Phi, lset.GetBndData(), lsgradrec);
        QuaQuaMapperCL quaqua( mg, lset.Phi, lsgradrec, /*neighborhoods*/ 0, /*maxiter*/ 100, /*tol*/ 1e-7, /*use_line_search*/ false, /*armijo_c*/ 1e-2);
//         accu.set_mapper( &quaqua);
        accus( mg.GetTriangTetraBegin(), mg.GetTriangTetraEnd());
        accu.set_mapper( 0);

        if (is_cake) { // needs an additional refinement to look good
            lset.idx.DeleteNumbering( mg);
            p1idx.DeleteNumbering( mg);
            vecp2idx.DeleteNumbering( mg);
            mg.Refine();
            lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
            lset.Phi.SetIdx( &lset.idx);
            lset.Init( csg_fun);
            p1idx.CreateNumbering( mg.GetLastLevel(), mg, lsbnd);
            ierr.SetIdx( &p1idx);
            ierrg.SetIdx( &p1idx);
            ierrq.SetIdx( &p1idx);
            idist.SetIdx( &p1idx);
            vecp2idx.CreateNumbering( mg.GetLastLevel(), mg, vecp2bnd);
            lsgradrec.SetIdx( &vecp2idx);
            accus( mg.GetTriangTetraBegin(), mg.GetTriangTetraEnd());
        }

        vtkwriter.Write( 2*i + 1);

        lset.idx.DeleteNumbering( mg);
        p1idx.DeleteNumbering( mg);
        vecp2idx.DeleteNumbering( mg);
        DROPS_FOR_TRIANG_TETRA( mg, mg.GetLastLevel(), it)
            it->SetRemoveMark(); // level 0 will not be removed.
        mg.Refine();
        if (is_cake) {
            DROPS_FOR_TRIANG_TETRA( mg, mg.GetLastLevel(), it)
                it->SetRemoveMark(); // level 0 will not be removed.
            mg.Refine();
        }
        lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
        lset.Phi.SetIdx( &lset.idx);
        lset.Init( csg_fun);
        p1idx.CreateNumbering( mg.GetLastLevel(), mg, lsbnd);
        ierr.SetIdx( &p1idx);
        ierrg.SetIdx( &p1idx);
        ierrq.SetIdx( &p1idx);
        idist.SetIdx( &p1idx);
        vecp2idx.CreateNumbering( mg.GetLastLevel(), mg, vecp2bnd);
        lsgradrec.SetIdx( &vecp2idx);

        // seq_out( Addr( ierr.Data), Addr( ierr.Data) + ierr.Data.size(), std::cout);
        delete thebody;
    }
    lset.idx.DeleteNumbering( mg);
    p1idx.DeleteNumbering( mg);
    vecp2idx.DeleteNumbering( mg);
    std::cout << "Successfully proccessed all examples.\n";
    return 0;
}


int TestAdap (MultiGridCL& mg, ParamCL& p)
{
    SurfaceTensionCL sf( sigmaf);   // dummy class
    LsetBndDataCL lsbnd( 6);
    std::auto_ptr<LevelsetP2CL> lset_ptr( LevelsetP2CL::Create( mg, lsbnd, sf));
    LevelsetP2CL& lset= *lset_ptr;

    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);

    DistMarkingStrategyCL dist_marker( csg_fun,
        P.get<double>("AdaptRef.Width"),
        P.get<Uint>(  "AdaptRef.CoarsestLevel"),
        P.get<Uint>(  "AdaptRef.FinestLevel"));
    CurvatureMarkingStrategyCL curv_marker( lset,
        PrincipalLatticeCL::instance( 2),
        P.get<Uint>( "AdaptRef.CurvatureFinestLevel"));
    StrategyCombinerCL marker;
    marker.push_back( dist_marker);
    if (curv_marker.GetFineLevel() > dist_marker.GetFineLevel())
        marker.push_back( curv_marker);

    DROPS::AdapTriangCL adap( mg, &marker);
    LevelsetReinitCL lsetreinit( lset, &csg_fun);
    adap.push_back( &lsetreinit);

    IdxDescCL p1idx;
    IdxDescCL p2idx( P2_FE);
    IdxDescCL vecp2idx( vecP2_FE);
    NoBndDataCL<> p2bnd;
    NoBndDataCL<Point3DCL> vecp2bnd;
    VecDescCL ierr( &p1idx),
              ierrg( &p1idx),
              ierrq( &p1idx),
              idist( &p1idx),
              lsgradrec( &vecp2idx),
              redist_err( &p2idx);

    const size_t num= std::distance( p.begin(), p.end());

   // writer for vtk-format
    VTKOutCL vtkwriter( mg, "DROPS data", num,
                        P.get<std::string>( "VTK.VTKDir"),
                        P.get<std::string>( "VTK.VTKName"),
                        P.get<std::string>( "VTK.TimeFileName"),
                        P.get<bool>( "VTK.Binary"),
                        P.get<bool>( "VTK.UseOnlyP1"), /* <- onlyp1 */
                        false, /* <- p2dg */
                        -1, /* <- level */
                        false, /* <- reusepvd */
                        false /* <- usedeformed */
              );
    vtkwriter.Register( make_VTKScalar( lset.GetSolution(), "level-set") );
    vtkwriter.Register( make_VTKScalar( make_P1Eval( mg, lsbnd, ierr), "interpolation-errorH") );
    vtkwriter.Register( make_VTKScalar( make_P1Eval( mg, lsbnd, ierrg), "interpolation-errorG") );
    vtkwriter.Register( make_VTKScalar( make_P1Eval( mg, lsbnd, ierrq), "error-quotient") );
    vtkwriter.Register( make_VTKScalar( make_P1Eval( mg, lsbnd, idist), "dh") );
    vtkwriter.Register( make_VTKVector( make_P2Eval( mg, vecp2bnd, lsgradrec), "ls_grad_rec") );
    vtkwriter.Register( make_VTKScalar( make_P2Eval( mg, p2bnd, redist_err), "redist_err") );

    // InterfaceApproxErrAccuCL accu0( lset, &ierr, &ierrg, &ierrq);
//     InterfaceApproxErrorAccuCL accu0( lset, &ierr, &ierrg, &ierrq, &idist, true);
//     TetraAccumulatorTupleCL accus0;
//     accus0.push_back( &accu0);

    size_t i= 0;
    for (ParamCL::ptree_const_iterator_type it= p.begin(); it != p.end(); ++it, ++i) {
        std::cout << "\n\n#Processing example \"" << it->first <<"\"." << std::endl;

        thebody= CSG::body_builder( it->second);
        lset.Init( csg_fun);
        adap.UpdateTriang();
        p1idx.CreateNumbering( mg.GetLastLevel(), mg, lsbnd);
        p2idx.CreateNumbering( mg.GetLastLevel(), mg, p2bnd);
        vecp2idx.CreateNumbering( mg.GetLastLevel(), mg, vecp2bnd);
        ierr.SetIdx( &p1idx);
        ierrg.SetIdx( &p1idx);
        ierrq.SetIdx( &p1idx);
        idist.SetIdx( &p1idx);
        lsgradrec.SetIdx( &vecp2idx);
        redist_err.SetIdx( &p2idx);
        averaging_P2_gradient_recovery( mg, lset.Phi, lset.GetBndData(), lsgradrec);
        QuaQuaMapperCL quaqua( mg, lset.Phi, lsgradrec, /*neighborhoods*/ 0,
            /*maxiter*/ P.get<int>( "LevelsetMapper.Iter"),
            /*tol*/ P.get<double>( "LevelsetMapper.Tol"),
            /*use_line_search*/ P.get<std::string>( "LevelsetMapper.Method") == "FixedPointWithLineSearch",
            /*armijo_c*/ P.get<double>( "LevelsetMapper.ArmijoConstant"),
            /*max_damping_steps*/ P.get<Uint>( "LevelsetMapper.MaxDampingSteps"));
        quaqua.set_inner_iter_tol( P.get<int>(    "LevelsetMapper.InnerIter"),
                                   P.get<double>( "LevelsetMapper.InnerTol"));
        quaqua.get_locator().set_structured_coarse_grid( P.get_child( "Domain"));

        TetraAccumulatorTupleCL accus;
        InterfaceCommonDataP2CL cdatap2( lset.Phi, lset.GetBndData(), quaqua, PrincipalLatticeCL::instance( P.get<Uint>( "Subsampling")));
        cdatap2.compute_quaddomains( P.get<std::string>("Testcase") == "area");
        accus.push_back( &cdatap2);
        InterfaceL2AccuP2CL l2accu( cdatap2, mg, "Area-Accumulator");
        l2accu.set_sample_curvature(P.get<std::string>("Testcase") == "area" && P.get<bool>("SampleCurvature") == true);
        if (P.get<std::string>("Testcase") == "redistancing") {
            l2accu.compute_redist_err( &redist_err);
            l2accu.set_function( &torus_distance);
        }
        accus.push_back( &l2accu);
//         InterfaceDebugP2CL p2debugaccu( cdatap2);
//         accus.push_back( &p2debugaccu);
//         accumulate( accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetMatchingFunction(), ifacep2idx.GetBndInfo());
        if (P.get<std::string>("Testcase") == "area") {
            const MultiGridCL& cmg= mg; // To get the const-iterators.
            ColorClassesCL cc;
            cc.make_single_color_class( cmg.GetTriangTetraBegin(), cmg.GetTriangTetraEnd());
            accus( cc);
        }
        else {
            accumulate( accus, mg, p2idx.TriangLevel(), p2idx.GetMatchingFunction(), p2idx.GetBndInfo());
        }
        std::cout << "Distribution of outer iterations in thread 0:\n";
        seq_out( cdatap2.quaqua.num_outer_iter.begin(), cdatap2.quaqua.num_outer_iter.end(), std::cout);
        std::cout << "Distribution of inner iterations in thread 0:\n";
        seq_out( cdatap2.quaqua.num_inner_iter.begin(), cdatap2.quaqua.num_inner_iter.end(), std::cout);

//         accu0.set_mapper( &quaqua);
//         accus0( mg.GetTriangTetraBegin(), mg.GetTriangTetraEnd());
//         accu0.set_mapper( 0);
//         seq_out( quaqua.num_outer_iter.begin(), quaqua.num_outer_iter.end(), std::cout);

        if (P.get<std::string>( "VTK.VTKName") != "")
            vtkwriter.Write( i);

        lset.idx.DeleteNumbering( mg);
        p1idx.DeleteNumbering( mg);
        p2idx.DeleteNumbering( mg);
        vecp2idx.DeleteNumbering( mg);

        delete thebody;
    }
    std::cout << "Successfully proccessed all examples.\n";
    return 0;
}

/// \brief Set Default parameters here s.t. they are initialized.
/// The result can be checked when Param-list is written to the output.
void SetMissingParameters (DROPS::ParamCL& P)
{
//     P.put_if_unset<std::string>("VTK.VTKName", "csgtest");
    P.put_if_unset<std::string>("VTK.VTKDir", ".");
    P.put_if_unset<std::string>( "VTK.TimeFileName", P.get<std::string>("VTK.VTKName"));
    P.put_if_unset<int>( "VTK.Binary", 1);
    P.put_if_unset<int>( "VTK.UseOnlyP1", 0);
    P.put_if_unset<int>( "VTK.ReUseTimeFile", 0);
    P.put_if_unset<int>("VTK.UseDeformation", 0);
}

int main (int argc, char** argv)
{
    ScopeTimerCL timer( "main");
  try {
    DROPS::read_parameter_file_from_cmdline( P, argc, argv, "csgtest.json");
    SetMissingParameters(P);
    std::cout << P << std::endl;

    deco_cube_radius= P.get<double>( "deco_cube_parameters.radius");
    deco_cube_shift=  P.get<double>( "deco_cube_parameters.shift");

    if (P.get<std::string>( "Testcase") == "redistancing") {
        torus_body= CSG::body_builder( P.get_child( "perturbed_torus_parameters.Torus"));
        perturbed_torus_alpha= P.get<double>( "perturbed_torus_parameters.alpha");
    }

    std::auto_ptr<DROPS::MGBuilderCL> builder( DROPS::make_MGBuilder( P.get_child( "Domain")));
    DROPS::MultiGridCL mg( *builder);
    if (P.get<std::string>( "Testcase") == "area") {
        MarkAll( mg);
        mg.Refine();
        MarkAll( mg);
        mg.Refine();
    }

    // If the key CSGLevelsets.File exists and is not empty, read the examples from the given file; otherwise, assume that the section CSGLevelsets itself contains the examples.
    std::string examples;
    try {
        examples= P.get<std::string>( "CSGLevelsets.File");
    }
    catch (DROPSParamErrCL&) {}
    ParamCL ex;
    if (examples != "") {
        std::ifstream jsonfile( examples.c_str());
        jsonfile >> ex;
    }
    else {
        ex= P.get_child( "CSGLevelsets");
        std::cout << ex << std::endl;
    }
//     return TestExamples( mg, ex);
    TestAdap( mg, ex);

    if (P.get<std::string>( "Testcase") == "redistancing")
        delete torus_body;

    rusage usage;
    getrusage( RUSAGE_SELF, &usage);
    printf( "ru_maxrss: %li kB.\n", usage.ru_maxrss);
    return 0;
  }
  catch (DROPSErrCL err) { err.handle(); }
}
