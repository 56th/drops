/// \file
/// \brief Discretization for PDEs on an interface.
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
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#include "geom/multigrid.h"
#include "misc/problem.h"
#include "num/discretize.h"
#include "num/precond.h"
#include "num/krylovsolver.h"
#include "num/interfacePatch.h"
#include "num/accumulator.h"
#include "levelset/levelset.h"
#include "levelset/levelsetmapper.h"
#include "levelset/mgobserve.h"
#include "out/ensightOut.h"
#include "out/vtkOut.h"
#include "num/fe_repair.h"
#include "num/quadrature.h"
#include "num/oseensolver.h"


#ifndef DROPS_IFACETRANSP_H
#define DROPS_IFACETRANSP_H

namespace DROPS {

/// \brief Given a P1IF_FE/P2IF_FE-vecdesc (on the interface) x and a P1/P2-vecdesc xext (on mg), the values of x
///        are copied to xext and the remaining values of xext are set to 0.
void Extend (const MultiGridCL& mg, const VecDescCL& x, VecDescCL& xext);

/// \brief Given a P1IF_FE/P2IF_FE-vecdesc (on the interface) x and a P1/P2-vecdesc xext (on mg), the values of xext
///        are copied to x, where x has an unknown.
void Restrict (const MultiGridCL& mg, const VecDescCL& xext, VecDescCL& x);

/// \brief Helper for the accumulators: inserts the local matrix loc.coup into M.
/// Works for P1IF_FE and P2IF_FE.
template <class LocalMatrixT>
  void
  update_global_matrix (MatrixBuilderCL& M, const LocalMatrixT& loc, const IdxT* numr, const IdxT* numc);

/// \brief The routine sets up the mass-matrix in matM on the interface defined by ls.
///        It belongs to the FE induced by standard P1-elements.
void SetupInterfaceMassP1 (const MultiGridCL& MG, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd, double alpha= 1.);

/// \brief The routine sets up the mass-matrix in matM on the interface defined by ls.
///        The FE spaces associated with rows and columns can be mixed among P1 and P1X.
///        The jumping coefficient \a alpha is defined w.r.t. the positive(0) and negative(1) part of the domain. If omitted, \a alpha is assumed to be 1 on the whole domain.
void SetupInterfaceMassP1X (const MultiGridCL& MG, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd, const double alpha[2]= 0);

/// \brief The routine equips the matrices with the right FE spaces \a mass_idx and \a surf_idx and
///        sets up the ad/de-sorption terms for coupled mass/surfactant transport.
///        The jumping coefficients \a k_a, \a k_d are defined w.r.t. the positive(0) and negative(1) part of the domain.
void SetupInterfaceSorptionP1X (const MultiGridCL& MG, const VecDescCL& ls, const BndDataCL<>& lsetbnd,
        MatDescCL* R, MatDescCL* C, MatDescCL* R_i, MatDescCL* C_i, const IdxDescCL* mass_idx, const IdxDescCL* surf_idx, const double k_a[2], const double k_d[2]);

/// \brief Copies P1IF_FE-unknown-indices or P2IF_FE-indices from idx on s into Numb.
/// Non-existent dofs get NoIdx.
void GetLocalNumbInterface(IdxT* Numb, const TetraCL& s, const IdxDescCL& idx);


/// \brief Helper for the accumulators: inserts the local matrix coup into M.
    void update_global_matrix_P1 (MatrixBuilderCL& M, const double coup[4][4], const IdxT numr[4], const IdxT numc[4]);

/// \todo This should be a generic function somewhere in num or misc.
void P1Init (InstatScalarFunction icf, VecDescCL& ic, const MultiGridCL& mg, double t);

/// \brief Resize normal according to qdom and fill in surf.normal. The normal must be precomputed.
template <Uint Dim>
  void
  resize_and_scatter_piecewise_normal (const SPatchCL<Dim>& surf, const QuadDomainCodim1CL<Dim>& qdom, std::valarray<typename SPatchCL<Dim>::WorldVertexT>& normal);


/// \brief The routine sets up the mass-matrix in matM on the interface defined by ls.
///        It belongs to the FE induced by standard P1-elements.
void SetupInterfaceMassP1 (const MultiGridCL& mg, MatDescCL* matM, const VecDescCL& ls, const BndDataCL<>& lsetbnd, double alpha);

/// \brief The routine sets up the mixed mass-matrix on the interface given by ls: The rows belong
///        to the new timestep, the columns to the old timestep. It belongs to the FE induced by
///        standard P1-elements. Actually only an alias for SetupInterfaceMassP1 with RowIdx != ColIdx
void SetupMixedMassP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsbnd);


class InterfaceCommonDataP1CL : public TetraAccumulatorCL
{
  private:
    InterfaceCommonDataP1CL** the_clones;

    const VecDescCL*   ls;      // P2-level-set
    const BndDataCL<>* lsetbnd; // boundary data for the level set function
    LocalP2CL<> locp2_ls;

  public:
    const PrincipalLatticeCL& lat;
    LocalP1CL<> p1[4];

    std::valarray<double> ls_loc;
    SurfacePatchCL surf;

    const InterfaceCommonDataP1CL& get_clone () const {
        const int tid= omp_get_thread_num();
        return tid == 0 ? *this : the_clones[tid][0];
    }

    bool empty () const { return surf.empty(); }

    InterfaceCommonDataP1CL (const VecDescCL& ls_arg, const BndDataCL<>& lsetbnd_arg)
        : ls( &ls_arg), lsetbnd( &lsetbnd_arg), lat( PrincipalLatticeCL::instance( 2)), ls_loc( lat.vertex_size())
    { p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; } // P1-Basis-Functions
    InterfaceCommonDataP1CL ()
        : ls( 0), lsetbnd( 0), lat( PrincipalLatticeCL::instance( 2))
    { p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; } // P1-Basis-Functions

    virtual ~InterfaceCommonDataP1CL () {}

    virtual void begin_accumulation   () {
        the_clones= new InterfaceCommonDataP1CL*[omp_get_max_threads()];
        the_clones[0]= this;
    }
    virtual void finalize_accumulation() {
        delete[] the_clones;
    }
    virtual void visit (const TetraCL& t) {
        surf.clear();
        locp2_ls.assign( t, *ls, *lsetbnd);
        evaluate_on_vertexes( locp2_ls, lat, Addr( ls_loc));
        if (equalSigns(ls_loc))
            return;
        surf.make_patch<MergeCutPolicyCL>( lat, ls_loc);
    }
    virtual InterfaceCommonDataP1CL* clone (int clone_id) {
        return the_clones[clone_id]= new InterfaceCommonDataP1CL( *this);
    }
};

typedef std::pair<const TetraCL*, BaryCoordCL> TetraBaryPairT;
typedef std::vector<TetraBaryPairT>            TetraBaryPairVectorT;

template <class T, class ResultIterT>
  inline ResultIterT
  evaluate_on_vertexes (std::function<T(const Point3DCL&, double)> const & f, const TetraBaryPairVectorT& pos, double t, ResultIterT result_iterator);

template <class T, class ResultContT>
  inline ResultContT&
  resize_and_evaluate_on_vertexes (std::function<T(const Point3DCL&, double)> const & f, const TetraBaryPairVectorT& pos, double t, ResultContT& result_container);

template <class PEvalT, class ResultIterT>
  inline ResultIterT
  evaluate_on_vertexes (const PEvalT& f, const TetraBaryPairVectorT& pos, ResultIterT result_iterator);

template <class PEvalT, class ResultContT>
  inline ResultContT&
  resize_and_evaluate_on_vertexes (const PEvalT& f, const TetraBaryPairVectorT& pos, ResultContT& result_container);


class ProjectedQuadDomain2DCL
{
  private:
    bool compute_absdets_;
    const QuadDomain2DCL* qdom;
    TetraBaryPairVectorT  vertexes_;
    std::valarray<double> absdets_;

    void resize (size_t s) {
        vertexes_.resize( s);
        if (compute_absdets_)
            absdets_.resize( s);
    }

  public:
    ProjectedQuadDomain2DCL () : compute_absdets_( true), qdom( 0) {}
    void assign (const SurfacePatchCL& p, const QuadDomain2DCL& qdomarg, const QuaQuaMapperCL& quaqua);

    void compute_absdets (bool b) { compute_absdets_= b; }

    const TetraBaryPairVectorT&  vertexes() const { return vertexes_; }
    const std::valarray<double>& absdets () const { return absdets_; }
};


class InterfaceCommonDataP2CL : public TetraAccumulatorCL
{
  private:
    InterfaceCommonDataP2CL** the_clones;

    const VecDescCL*   ls;      // P2-level-set
    const BndDataCL<>* lsetbnd; // boundary data for the level set function

    const PrincipalLatticeCL* lat;

    bool compute_quaddomains_;

  public:
    /// common data @{
    LocalP2CL<> locp2_ls;

    LocalP2CL<>          p2[10];
    LocalP1CL<Point3DCL> gradrefp2[10];

    std::valarray<double>     ls_loc;
    SurfacePatchCL            surf;
    QuadDomain2DCL            qdom;
    ProjectedQuadDomain2DCL   qdom_projected;
    QuaQuaMapperCL            quaqua;

    const PrincipalLatticeCL& get_lattice () const { return *lat; }
    /// @}

    LocalP2CL<> get_local_p2_ls (const TetraCL& t) const { return LocalP2CL<>( t, *ls, *lsetbnd); };

    const InterfaceCommonDataP2CL& get_clone () const {
        const int tid= omp_get_thread_num();
        return tid == 0 ? *this : the_clones[tid][0];
    }

    bool empty () const { return surf.empty(); }

    void compute_absdet (bool b) { qdom_projected.compute_absdets( b); }
    void compute_quaddomains (bool b) { compute_quaddomains_= b; }

    void set_lattice (const PrincipalLatticeCL& newlat) {
        lat= &newlat;
        ls_loc.resize( lat->vertex_size());
    }

    InterfaceCommonDataP2CL (const VecDescCL& ls_arg, const BndDataCL<>& lsetbnd_arg,
        const QuaQuaMapperCL& quaquaarg, const PrincipalLatticeCL& lat_arg);
    virtual ~InterfaceCommonDataP2CL () {}

    virtual void begin_accumulation () {
        the_clones= new InterfaceCommonDataP2CL*[omp_get_max_threads()];
        the_clones[0]= this;
    }
    virtual void finalize_accumulation() { delete[] the_clones; }

    virtual void visit (const TetraCL& t) {
        surf.clear();
        locp2_ls.assign( t, *ls, *lsetbnd);
        evaluate_on_vertexes( locp2_ls, *lat, Addr( ls_loc));
        if (equalSigns(ls_loc))
            return;
        surf.make_patch<MergeCutPolicyCL>( *lat, ls_loc);
        if (surf.empty())
            return;

        if (compute_quaddomains_) {
            make_CompositeQuad5Domain2D ( qdom, surf, t);
            quaqua.set_point( &t, BaryCoordCL()); // set the current tetra.
            qdom_projected.assign( surf, qdom, quaqua);
        }
    }

    virtual InterfaceCommonDataP2CL* clone (int clone_id) {
        return the_clones[clone_id]= new InterfaceCommonDataP2CL( *this);
    }
};


class InterfaceCommonDataDeformP2CL;

class LocalMeshTransformationCL
{
  public:
    InterfaceCommonDataDeformP2CL* cdata; // Must be set before use of members.

    LocalP2CL<Point3DCL> Psi; // Phi= id + Psi
    Bary2WorldCoordCL b2w;
    World2BaryCoordCL w2b;

    Point3DCL n_lin;
    SMatrixCL<3, 2> Q;

    LocalP1CL<SMatrixCL<3, 3> > dPhi;
    SMatrixCL<3, 3> G,
                    Ginv_wwT;
    QRDecompCL<3, 3> Gqr;
    Point3DCL w;  // the pull-back of n_Gamma (scaled to unit length)

    double JPhi,
           JPhiQ;

    LocalMeshTransformationCL (InterfaceCommonDataDeformP2CL* cdataarg= 0)
        : cdata (cdataarg) {}

    void set_tetra (const TetraCL* t);

    void set_surface_patch (const BaryCoordCL verts[3], const Point3DCL& pos_pt) { // Set Q, n_lin
        QRDecompCL<3, 2> qr;
        SMatrixCL<3, 2>& M= qr.GetMatrix ();
        M.col( 0, b2w( verts[1]) - b2w( verts[0]));
        M.col( 1, b2w( verts[2]) - b2w( verts[0]));
        qr.prepare_solve();
        const SMatrixCL<3, 3> QQ= qr.get_Q ();
        Q.col(0, QQ.col(0));
        Q.col(1, QQ.col(1));
        n_lin= QQ.col (2)*sign (inner_prod (pos_pt - b2w (verts[1]), QQ.col (2))); // n_lin points out of the neg. domain.
    }

    void set_point (const BaryCoordCL& xb, bool surface_data_p) {
        SMatrixCL<3, 3> dPhix= dPhi(xb);
        G= GramMatrix (dPhix);
        Gqr.GetMatrix ()= G;
        Gqr.prepare_solve ();
        JPhi= std::sqrt (std::abs (Gqr.Determinant_R ())); // |det dPhix|
        if (!surface_data_p)
            return;

        w= n_lin;
        Gqr.Solve (w);
        double nlinT_Ginv_nlin= inner_prod (w, n_lin);
        w/= std::sqrt (nlinT_Ginv_nlin);
        Ginv_wwT= eye<3, 3> ();
        Gqr.Solve (Ginv_wwT);
        Ginv_wwT-= outer_product (w, w);

//        QRDecompCL<3, 3> qrdphi;
//        qrdphi.GetMatrix() = dPhix;
//        qrdphi.prepare_solve();
//        SMatrixCL<3, 3> dphiinv = eye<3, 3> ();
//        SMatrixCL<3, 3> dphiinvT = eye<3, 3> ();
//        qrdphi.Solve(dphiinv);
//        assign_transpose(dphiinvT, dphiinv);
//        SMatrixCL<3, 3> Proj = eye<3, 3> ();
//        Point3DCL p = b2w(xb);
//        DROPS::Point3DCL v(p[0]/(std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])),p[1]/(std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])),p[2]/(std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])));
//        Proj-= outer_product(v,v);
//        Ginv_wwT = dphiinv*Proj*dphiinvT;



        const SMatrixCL<2, 2> gr= GramMatrix (dPhix*Q);
        JPhiQ= std::sqrt(gr(0, 0)*gr(1, 1) - gr(0, 1)*gr(1, 0));
    }

    void map_QuadDomain (QuadDomainCL& qdom);

    void map_QuadDomain2D (QuadDomain2DCL& qdom, const SurfacePatchCL& p, const Point3DCL& pos_pt) {
        const Uint nodes_per_facet= qdom.vertex_size()/p.facet_size();
        for (Uint i= 0; i < qdom.vertex_size(); ++i) {
            if (i % nodes_per_facet == 0) {
                const SurfacePatchCL::FacetT& facet= p.facet_begin()[i/nodes_per_facet];
                const BaryCoordCL verts[3]= { p.vertex_begin()[facet[0]],
                                              p.vertex_begin()[facet[1]],
                                              p.vertex_begin()[facet[2]] };
                set_surface_patch (verts, pos_pt);
            }
            set_point (qdom.vertex_begin ()[i], /*surface_data_p=*/ true);
            qdom.vertex_begin ()[i]+= w2b.map_direction( Psi (qdom.vertex_begin ()[i]));
            qdom.weight_begin ()[i]*= JPhiQ;
        }
    }
};

class InterfaceCommonDataDeformP2CL : public TetraAccumulatorCL
{
  private:
    InterfaceCommonDataDeformP2CL** the_clones;

    const VecDescCL*   ls;      // P2-level-set
    const BndDataCL<>* lsetbnd; // boundary data for the level set function

    const PrincipalLatticeCL* lat;

  public:
    /// common data @{
    VecDescCL* Psi_vd;
    LocalP2CL<> locp2_ls;

    double det_T;
    SMatrixCL<3,3> T;
    LocalP2CL<>          p2[10];
    LocalP1CL<Point3DCL> gradrefp2[10];
    LocalP1CL<Point3DCL> gradp2[10];

    std::valarray<double>     ls_loc;
    SurfacePatchCL            surf;
    QuadDomain2DCL            qdom2d_only_weights, // For the surface bilinear forms: points on linear interface, mapped weights.
                              qdom2d_full; // For the right-hand side and error measurement: mapped points and weights
    QuadDomainCL              qdom; // points not mapped, weights mapped.
    Point3DCL                 pos_pt;

    mutable LocalMeshTransformationCL Phi;

    const PrincipalLatticeCL& get_lattice () const { return *lat; }
    /// @}

    InterfaceCommonDataDeformP2CL& get_clone () {
        const int tid= omp_get_thread_num();
        return tid == 0 ? *this : the_clones[tid][0];
    }
    const InterfaceCommonDataDeformP2CL& get_clone () const {
        const int tid= omp_get_thread_num();
        return tid == 0 ? *this : the_clones[tid][0];
    }

    bool empty () const { return surf.empty(); }

    void set_lattice (const PrincipalLatticeCL& newlat) {
        lat= &newlat;
        ls_loc.resize( lat->vertex_size());
    }

    InterfaceCommonDataDeformP2CL (const VecDescCL& ls_arg, const BndDataCL<>& lsetbnd_arg,
        VecDescCL& Psi_vdarg, const PrincipalLatticeCL& lat_arg);
    virtual ~InterfaceCommonDataDeformP2CL () {}

    virtual void begin_accumulation () {
        the_clones= new InterfaceCommonDataDeformP2CL*[omp_get_max_threads()];
        the_clones[0]= this;
    }
    virtual void finalize_accumulation() {
        delete[] the_clones;
    }

    virtual void visit (const TetraCL& t) {
        surf.clear();
        locp2_ls.assign( t, *ls, *lsetbnd);
        evaluate_on_vertexes( locp2_ls, *lat, Addr( ls_loc));
        if (equalSigns(ls_loc))
            return;
        surf.make_patch<MergeCutPolicyCL>( *lat, ls_loc);
        if (surf.empty())
            return;

        GetTrafoTr( T, det_T, t);
        P2DiscCL::GetGradients( gradp2, gradrefp2, T);

        Phi.set_tetra (&t);
        make_CompositeQuad5Domain2D (qdom2d_full, surf, t);
        qdom2d_only_weights= qdom2d_full;
        Uint i= 0;
        for (; ls_loc[i] <= 0.; ++i)
            ;
        if (i > 3)
            std::cerr << "InterfaceCommonDataDeformP2CL::visit: No positive vertex.\n";
        pos_pt= t.GetVertex (i)->GetCoord ();
        Phi.map_QuadDomain2D (qdom2d_full, surf, pos_pt);
        std::copy (qdom2d_full.weight_begin (), qdom2d_full.weight_begin () + qdom2d_full.vertex_size (),
                   qdom2d_only_weights.weight_begin ());

        make_SimpleQuadDomain<Quad5DataCL> (qdom, AllTetraC);
        Phi.map_QuadDomain (qdom);
    }

    virtual InterfaceCommonDataDeformP2CL* clone (int clone_id) {
        the_clones[clone_id]= new InterfaceCommonDataDeformP2CL( *this);
        the_clones[clone_id]->Phi.cdata= the_clones[clone_id];
        return the_clones[clone_id];
    }
};


inline void LocalMeshTransformationCL::set_tetra (const TetraCL* t)
{ // Set Psi, dPhi
    b2w.assign (*t);
    w2b.assign (*t);
    const NoBndDataCL<Point3DCL> nobnd;
    Psi.assign (*t, *cdata->Psi_vd, nobnd);
    dPhi= eye<3, 3> ();
    for (Uint i= 0; i < 10; ++i)
        dPhi+= outer_product (Psi[i], cdata->gradp2[i]);
}

inline void LocalMeshTransformationCL::map_QuadDomain (QuadDomainCL& qdom)
{
    for (Uint i= 0; i < qdom.vertex_size (); ++i) {
        set_point (qdom.vertex_begin ()[i], /*surface_data_p=*/ false);
        qdom.weight_begin ()[i]*= JPhi*std::abs(cdata->det_T);
    }
}


/// \brief Map the QuadDomain2DCL from InterfaceCommonDataP2CL to the quadratic levelset.
/// The result is represented as unordered map from tetras to (mapped) QuadDomain2DCL.
/// The quadrature weights are adjusted by the absdet of the mapping.
/// The mapping uses QuaQuaMapperCL.
/// XXX Not OpenMP-thread-safe.
class QuaQuaQuadDomainMapperAccuCL : public TetraAccumulatorCL
{
  private:
    const InterfaceCommonDataP2CL& cdata_;

  public:
    DROPS_STD_UNORDERED_MAP<const TetraCL*, QuadDomain2DCL> qmap;

    QuaQuaQuadDomainMapperAccuCL (const InterfaceCommonDataP2CL& cdata)
        : cdata_( cdata)
    {}

    void begin_accumulation () {
        std::cerr << "QuaQuaQuadDomainMapperAccuCL::begin_accumulation.\n";
        qmap.clear();
    }
    void finalize_accumulation() {
        std::cerr << "QuaQuaQuadDomainMapperAccuCL::finalize_accumulation.\n";
    }

    void visit (const TetraCL&) {
        const InterfaceCommonDataP2CL& cdata= cdata_.get_clone();
        if (cdata.empty())
            return;

        for (Uint i= 0; i < cdata.qdom.vertex_size(); ++i) {
            const TetraBaryPairT& p= cdata.qdom_projected.vertexes()[i];
            const double newweight= cdata.qdom.weight_begin()[i]*cdata.qdom_projected.absdets()[i];
            qmap[p.first].push_back_quad_node( p.second, newweight);
        }
    };

    TetraAccumulatorCL* clone (int /*tid*/) { return new QuaQuaQuadDomainMapperAccuCL( *this); };
};


template <class LocalMatrixT, class InterfaceCommonDataT>
class InterfaceMatrixAccuCL : public TetraAccumulatorCL
{
  private:
    const InterfaceCommonDataT& cdata_;
    std::string name_;

    MatDescCL* mat_; // the matrix
    MatrixBuilderCL* M;

    LocalMatrixT local_mat;

    Uint lvl;
    SArrayCL<IdxT, LocalMatrixT::row_fe_type == P1IF_FE ? 4 : 10> numr;
    SArrayCL<IdxT, LocalMatrixT::col_fe_type == P1IF_FE ? 4 : 10> numc;

  public:
    InterfaceMatrixAccuCL (MatDescCL* Mmat, const LocalMatrixT& loc_mat, const InterfaceCommonDataT& cdata,
                           std::string name= std::string())
        : cdata_( cdata), name_( name), mat_( Mmat), M( 0), local_mat( loc_mat) {}
    virtual ~InterfaceMatrixAccuCL () {}

    void set_name (const std::string& n) { name_= n; }

    virtual void begin_accumulation () {
        const IdxT num_rows= mat_->RowIdx->NumUnknowns();
        const IdxT num_cols= mat_->ColIdx->NumUnknowns();

        const FiniteElementT mat_row_fe= mat_->RowIdx->GetFE(),
                             mat_col_fe= mat_->ColIdx->GetFE();

        std::cout << "InterfaceMatrixAccuCL::begin_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\"";
        std::cout  << ": " << num_rows << " rows, " << num_cols << " cols.\n";
        if (mat_row_fe != LocalMatrixT::row_fe_type)
            std::cout << "Warning: LocalMatrixT and MatDescCL have different FE-type for rows.\n";
        if (mat_col_fe != LocalMatrixT::col_fe_type)
            std::cout << "Warning: LocalMatrixT and MatDescCL have different FE-type for cols.\n";

        lvl = mat_->GetRowLevel();
        M= new MatrixBuilderCL( &mat_->Data, num_rows, num_cols);
    }

    virtual void finalize_accumulation () {
        M->Build();
        delete M;
        M= 0;
       /* std::cout << "InterfaceMatrixAccuCL::finalize_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\"";
        std::cout << ": " << mat_->Data.num_nonzeros() << " nonzeros." << std::endl;*/
    }

    virtual void visit (const TetraCL& t) {
        const InterfaceCommonDataT& cdata= cdata_.get_clone();
        if (cdata.empty())
            return;
        local_mat.setup( t, cdata);
        GetLocalNumbInterface( numr.begin(), t, *mat_->RowIdx);
        GetLocalNumbInterface( numc.begin(), t, *mat_->ColIdx);
        update_global_matrix( *M, local_mat, numr.begin(), numc.begin());
    }

    virtual InterfaceMatrixAccuCL* clone (int /*clone_id*/) { return new InterfaceMatrixAccuCL( *this); }
};

/// \brief Accumulate an interface-vector.
template <class LocalVectorT, class InterfaceCommonDataT>
class InterfaceVectorAccuCL : public TetraAccumulatorCL
{
  private:
    const InterfaceCommonDataT& cdata_;
    std::string name_;

    VecDescCL* y_;
    LocalVectorT local_vec;
    SArrayCL<IdxT, LocalVectorT::row_fe_type == P1IF_FE ? 4 : 10> numry;

  public:
    InterfaceVectorAccuCL (VecDescCL* y, const LocalVectorT& loc_vec, const InterfaceCommonDataT& cdata,
                             std::string name= std::string())
        : cdata_( cdata), name_( name), y_( y), local_vec( loc_vec) {}
    virtual ~InterfaceVectorAccuCL () {}

    void set_name (const std::string& n) { name_= n; }

    virtual void begin_accumulation () {
        /*std::cout << "InterfaceVectorAccuCL::begin_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\"";
        std::cout  << ": " << y_->RowIdx->NumUnknowns() << " rows.\n";*/
        if (y_->RowIdx->GetFE() != LocalVectorT::row_fe_type)
            std::cout << "Warning: LocalVectorT and VecDescCL have different FE-type for rows.\n";
    }

    virtual void finalize_accumulation() {}

    virtual void visit (const TetraCL& t) {
        const InterfaceCommonDataT& cdata= cdata_.get_clone();
        if (cdata.empty())
            return;
        GetLocalNumbInterface( numry.begin(), t, *y_->RowIdx);
        local_vec.setup( t, cdata, numry.begin());
        const int num_row_dofs= LocalVectorT::row_fe_type == P1IF_FE ? 4 : 10;
        for (int i= 0; i < num_row_dofs; ++i)
            if (numry[i] != NoIdx)
                y_->Data[numry[i]]+= local_vec.vec[i];
    }

    virtual InterfaceVectorAccuCL* clone (int /*clone_id*/) { return new InterfaceVectorAccuCL( *this); }
};


    class NarrowBandCommonDataP1CL : public TetraAccumulatorCL
    {
    private:
        NarrowBandCommonDataP1CL** the_clones;

        const VecDescCL*   ls;      // P2-level-set
        const BndDataCL<>* lsetbnd; // boundary data for the level set function
        const double & dist;//to characterize the width of the narrow band
        LocalP2CL<> locp2_ls;
        bool inband_;

    public:
        const PrincipalLatticeCL& lat;
        LocalP1CL<> p1[4];

        std::valarray<double> ls_loc;
        SurfacePatchCL surf;

        const NarrowBandCommonDataP1CL& get_clone () const {
            const int tid= omp_get_thread_num();
            return tid == 0 ? *this : the_clones[tid][0];
        }

        bool empty () const { return surf.empty(); }
        bool in_band() const {return inband_;}

        NarrowBandCommonDataP1CL (const VecDescCL& ls_arg, const BndDataCL<>& lsetbnd_arg, const double & dst)//need check the distantce~!!!!!!
                : ls( &ls_arg), lsetbnd( &lsetbnd_arg),dist(dst) ,lat( PrincipalLatticeCL::instance( 2)), ls_loc( lat.vertex_size())
        { p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; inband_=true;} // P1-Basis-Functions
        NarrowBandCommonDataP1CL ()
                : ls( 0), lsetbnd( 0),dist(0), lat( PrincipalLatticeCL::instance( 2))
        { p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; } // P1-Basis-Functions

        virtual ~NarrowBandCommonDataP1CL () {}

        virtual void begin_accumulation   () {
            the_clones= new NarrowBandCommonDataP1CL*[omp_get_max_threads()];
            the_clones[0]= this;
        }
        virtual void finalize_accumulation() {
            delete[] the_clones;
        }
        virtual void visit (const TetraCL& t) {
            surf.clear();
            locp2_ls.assign( t, *ls, *lsetbnd);
            evaluate_on_vertexes( locp2_ls, lat, Addr( ls_loc));
            //  std::cout<<dist<<std::endl;
            inband_=true;
            if (distance( ls_loc)>dist)
                // if (equal_signs( ls_loc))
            {
                inband_=false;
                return;
            }
            surf.make_patch<MergeCutPolicyCL>( lat, ls_loc);
        }
        virtual NarrowBandCommonDataP1CL* clone (int clone_id) {
            return the_clones[clone_id]= new NarrowBandCommonDataP1CL( *this);
        }
    };

    template <class LocalMatrixT>
    class NarrowBandMatrixAccuP1CL : public TetraAccumulatorCL
    {
    private:
        const NarrowBandCommonDataP1CL& cdata_;
        std::string name_;

        MatDescCL* mat_; // the matrix
        MatrixBuilderCL* M;

        LocalMatrixT local_mat;

        Uint lvl;
        IdxT numr[4],
                numc[4];

    public:
        NarrowBandMatrixAccuP1CL (MatDescCL* Mmat, const LocalMatrixT& loc_mat, const NarrowBandCommonDataP1CL& cdata,
                                  std::string name= std::string())
                : cdata_( cdata), name_( name), mat_( Mmat), M( 0), local_mat( loc_mat) {}
        virtual ~NarrowBandMatrixAccuP1CL () {}

        void set_name (const std::string& n) { name_= n; }

        virtual void begin_accumulation () {
            const IdxT num_rows= mat_->RowIdx->NumUnknowns();
            const IdxT num_cols= mat_->ColIdx->NumUnknowns();
            std::cout << "NarrowBandMatrixAccuP1CL::begin_accumulation";
            if (name_ != std::string())
                std::cout << " for \"" << name_ << "\"";
            std::cout  << ": " << num_rows << " rows, " << num_cols << " cols.\n";
            lvl = mat_->GetRowLevel();
            M= new MatrixBuilderCL( &mat_->Data, num_rows, num_cols);
        }

        virtual void finalize_accumulation () {
            M->Build();
            delete M;
            M= 0;
            std::cout << "NarrowBandMatrixAccuP1CL::finalize_accumulation";
            if (name_ != std::string())
                std::cout << " for \"" << name_ << "\"";
            std::cout << ": " << mat_->Data.num_nonzeros() << " nonzeros." << std::endl;
        }

        virtual void visit (const TetraCL& t) {

            const NarrowBandCommonDataP1CL& cdata= cdata_.get_clone();
            if (!cdata.in_band())
                return;
            local_mat.setup( t, cdata);
            GetLocalNumbP1NoBnd( numr, t, *mat_->RowIdx);
            GetLocalNumbP1NoBnd( numc, t, *mat_->ColIdx);
            update_global_matrix_P1( *M, local_mat.coup, numr, numc);
        }

        virtual NarrowBandMatrixAccuP1CL* clone (int /*clone_id*/) { return new NarrowBandMatrixAccuP1CL( *this); }
    };

/// \brief Compute the load-vector corresponding to the function f on a single tetra.
class LocalVectorP1CL
{
  private:
    InstatScalarFunction f_;
    double time_;

    LocalP1CL<> p1[4];
    std::valarray<double> qp1,
                          qf;
    QuadDomain2DCL qdom;

  public:
    static const FiniteElementT row_fe_type= P1IF_FE;

    double vec[4];

    LocalVectorP1CL (InstatScalarFunction f, double time) : f_( f), time_( time) { p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; }

    void setup (const TetraCL& t, const InterfaceCommonDataP1CL& cdata, const IdxT numr[4]) {
        make_CompositeQuad5Domain2D ( qdom, cdata.surf, t);
        resize_and_evaluate_on_vertexes( f_, t, qdom, time_, qf);
        qp1.resize( qdom.vertex_size());
        for (Uint i= 0; i < 4; ++i) {
                if (numr[i] == NoIdx)
                    continue;
                evaluate_on_vertexes( p1[i], qdom, Addr( qp1));
                vec[i]= quad_2D( qf*qp1, qdom);
        }
    }
};

/// \brief Compute local_mat_*x_ on a single tetra.
template <class  LocalMatrixT>
class LocalMatVecP1CL
{
  private:
    LocalMatrixT local_mat_;
    IdxT numx_[4];
    const VecDescCL& x_;

  public:
    static const FiniteElementT row_fe_type= P1IF_FE;

    double vec[4];

    LocalMatVecP1CL (const LocalMatrixT& local_mat, const VecDescCL* x) : local_mat_( local_mat), x_( *x) {}

    void setup (const TetraCL& t, const InterfaceCommonDataP1CL& cdata, const IdxT numr[4]) {
        local_mat_.setup( t, cdata);
        GetLocalNumbP1NoBnd( numx_, t, *x_.RowIdx);
        for (int i= 0; i < 4; ++i) {
            vec[i]= 0.;
            if (numr[i] != NoIdx)
                for (int j= 0; j < 4; ++j)
                    if (numx_[j] != NoIdx)
                        vec[i]+= local_mat_.coup[i][j]*x_.Data[numx_[j]];
        }
    }
};


/// \brief Convenience-function to reduce the number of explicit template-parameters for the massdiv- and the convection-rhs.
template <template <class> class LocalMatrixT, class DiscVelSolT>
  inline InterfaceVectorAccuCL<LocalMatVecP1CL< LocalMatrixT<DiscVelSolT> >, InterfaceCommonDataP1CL>*
  make_wind_dependent_vectorP1_accu (VecDescCL* y, const VecDescCL* x, const InterfaceCommonDataP1CL& cdata, const DiscVelSolT& wind, std::string name= std::string())
{
    return new InterfaceVectorAccuCL<LocalMatVecP1CL< LocalMatrixT<DiscVelSolT> >, InterfaceCommonDataP1CL>( y,
        LocalMatVecP1CL< LocalMatrixT<DiscVelSolT> >( LocalMatrixT<DiscVelSolT>( wind), x), cdata, name);
}

/// \brief Convenience-function to reduce the number of explicit template-parameters for the massdiv- and the convection-matrix.
template <template <class> class LocalMatrixT, class DiscVelSolT>
  inline InterfaceMatrixAccuCL< LocalMatrixT<DiscVelSolT>, InterfaceCommonDataP1CL>*
  make_wind_dependent_matrixP1_accu (MatDescCL* mat, const InterfaceCommonDataP1CL& cdata, const DiscVelSolT& wind, std::string name= std::string())
{
    return new InterfaceMatrixAccuCL< LocalMatrixT<DiscVelSolT>, InterfaceCommonDataP1CL>( mat,
        LocalMatrixT<DiscVelSolT>( wind), cdata, name);
}

    template <template <class> class LocalMatrixT, class DiscVelSolT>
    inline InterfaceMatrixAccuCL< LocalMatrixT<DiscVelSolT>, InterfaceCommonDataP1CL>*
    make_wind_dependent_matrixP1_accu (MatDescCL* mat, const InterfaceCommonDataP1CL& cdata, const  InstatVectorFunction normal, double time, const DiscVelSolT& wind, std::string name= std::string())
    {

        return new InterfaceMatrixAccuCL< LocalMatrixT<DiscVelSolT>, InterfaceCommonDataP1CL>( mat,
                                                                         LocalMatrixT<DiscVelSolT>( wind,normal,time), cdata, name);
    }
/// \brief Convenience-function to reduce the number of explicit template-parameters for the massdiv- and the convection-matrix.
    template <template <class> class LocalMatrixT, class DiscVelSolT>
    inline NarrowBandMatrixAccuP1CL< LocalMatrixT<DiscVelSolT>>*
    make_wind_dependent_matrixP1_accu (MatDescCL* mat, const NarrowBandCommonDataP1CL& cdata, const DiscVelSolT& wind, std::string name= std::string())
    {
        return new NarrowBandMatrixAccuP1CL< LocalMatrixT<DiscVelSolT>>( mat,
                                                                         LocalMatrixT<DiscVelSolT>( wind), cdata, name);
    }
/// \brief Convenience-function to reduce the number of explicit template-parameters for the massdiv- and the convection-matrix.
    template <template <class> class LocalMatrixT, class DiscVelSolT>
    inline NarrowBandMatrixAccuP1CL< LocalMatrixT<DiscVelSolT>>*
    make_wind_dependent_matrixP1_accu (MatDescCL* mat, const NarrowBandCommonDataP1CL& cdata, const  InstatVectorFunction normal, double time, const DiscVelSolT& wind, std::string name= std::string())
    {
        return new NarrowBandMatrixAccuP1CL< LocalMatrixT<DiscVelSolT>>( mat,
                LocalMatrixT<DiscVelSolT>( wind,normal,time), cdata, name);
    }

/// \brief Convenience-function to reduce the number of explicit template-parameters for the massdiv- and the convection-matrix.
    template <template <class> class LocalMatrixT, class DiscVelSolT>
    inline InterfaceMatrixAccuCL< LocalMatrixT<DiscVelSolT>, InterfaceCommonDataP1CL>*
    make_concentration_dependent_matrixP1_accu (MatDescCL* mat, const InterfaceCommonDataP1CL& cdata, const  InstatVectorFunction normal, double time, const DiscVelSolT& wind, std::string name= std::string())
    {
        return new InterfaceMatrixAccuCL< LocalMatrixT<DiscVelSolT>, InterfaceCommonDataP1CL>( mat,
                LocalMatrixT<DiscVelSolT>( wind, normal, time), cdata, name);
    }

class LocalInterfaceMassP1CL
{
  private:
    std::valarray<double> q[4];
    QuadDomain2DCL qdom;
    double alpha_;

  public:
    static const FiniteElementT row_fe_type= P1IF_FE,
                                col_fe_type= P1IF_FE;

    double coup[4][4];

    void setup (const TetraCL& t, const InterfaceCommonDataP1CL& cdata) {
        make_CompositeQuad2Domain2D ( qdom, cdata.surf, t);
        for (int i= 0; i < 4; ++i)
            resize_and_evaluate_on_vertexes ( cdata.p1[i], qdom, q[i]);
        for (int i= 0; i < 4; ++i) {
            coup[i][i]= quad_2D( q[i]*q[i], qdom);
            for(int j= 0; j < i; ++j)
                coup[i][j]= coup[j][i]= alpha_*quad_2D( q[j]*q[i], qdom);
        }
    }

    LocalInterfaceMassP1CL (double alpha= 1.) : alpha_( alpha) {}
};

/// \brief The routine sets up the Laplace-Beltrami-matrix in mat on the interface defined by ls.
///        It belongs to the FE induced by standard P1-elements.
///
/// D is the diffusion-coefficient
void SetupLBP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsbnd, double D);

class LocalAssembler {
public:
    struct LocalAssemblerParams {
        size_t numbOfVirtualSubEdges = 2;
        double t = 0.;
        VecDescCL levelSet;
        struct {
            double nu = 1., gamma = 0., Pe = 0., u_N_max = 0., rho_max = 1., rho_min = 1., lineTension = 0.;
            VecDescCL w_T, u_N, chi, omega;
            InstatVectorFunction f_T = nullptr; // moment rhs
            InstatScalarFunction m_g = nullptr; // - continuity eqn rhs
        } surfNavierStokesParams;
        struct {
            bool useDegenerateMobility = false;
            double mobilityScaling = 1.;
            VecDescCL u_T;
            InstatScalarFunction f = nullptr;
        } surfCahnHilliardParams;
    };
    using vec = std::vector<double>;
    using mtx = std::vector<vec>;
    using LinearForm = vec(LocalAssembler::*)();
    using BilinearForm = mtx(LocalAssembler::*)();
    LocalP2CL<> levelSetTet;
private:
    struct {
        size_t P1 = 4;
        size_t P2 = 10;
        size_t vecP2 = 30;
    } const n;
    LocalAssemblerParams const & params;
    TetraCL const & tet;
    SMatrixCL<3, 3> T;
    double absDet, rho_delta;
    QuadDomain2DCL qDomain;
    QuadDomainCL q3Domain;
    GridFunctionCL<> qHatP2[10], qHatP1[4], qLsGradNorm, qSurfSpeed, qG, qF_c, qChi, qOmega, qCutOffFunc;
    GridFunctionCL<Point3DCL> qGradP2[10], qSurfGradP2[10], q3DGradP2[10], qGradP1[4], qSurfGradP1[4], q3DGradP1[4], qNormal, q3DNormal, qSurfSpeedSurfGrad, qChiSurfGrad, qOmegaSurfGrad, qF, qHatP2CrossN[30];
    struct { GridFunctionCL<Point3DCL> NS, CH; } qWind;
    GridFunctionCL<SMatrixCL<3,3>> qP, qH, qE[30];
    mtx createMtx(size_t n, size_t m, double v) { return mtx(n, vec(v, m)); }
    mtx createMtx(size_t n, double v) { return createMtx(n, n, v); }
    mtx createMtx(size_t n, size_t m) { return mtx(n, vec(m)); }
    mtx createMtx(size_t n) { return createMtx(n, n); }
    mtx scaleMtx(double a, mtx const & A) {
        auto B = A;
        for (auto& row : B)
            for (auto& val : row)
                val *= a;
        return B;
    }
    using GridFunctionBuilder = void(LocalAssembler::*)();
    template<class T>
    void require(GridFunctionCL<T> const & function, GridFunctionBuilder builder) {
        if (!function.size()) (this->*builder)();
        if (!function.size()) throw std::logic_error(__func__ + std::string(": invalid builder"));
    }
    std::pair<size_t, size_t> ind(size_t n, size_t i) {
        /*
        auto is = i / 3; // scalar shape index
        auto in = i - 3 * is; // nonzero vect component
        return { is, in };
        */
        auto in = i / n; // nonzero vect component
        auto is = i - in * n; // scalar shape index
        return { is, in };
    }
    GridFunctionCL<Point3DCL> getSurfGradP1(GridFunctionCL<> const & f) {
        require(qSurfGradP1[0], &LocalAssembler::buildSurfGradP1);
        auto res = f[0] * qSurfGradP1[0];
        for (size_t i = 1; i < n.P1; ++i) res += f[i] * qSurfGradP1[i];
        return res;
    }
    GridFunctionCL<Point3DCL> getGradP2(GridFunctionCL<> const & f) {
        require(qGradP2[0], &LocalAssembler::buildGradP2);
        auto res = f[0] * qGradP2[0];
        for (size_t i = 1; i < n.P2; ++i) res += f[i] * qGradP2[i];
        return res;
    }
    GridFunctionCL<Point3DCL> get3DGradP2(GridFunctionCL<> const & f) {
        require(q3DGradP2[0], &LocalAssembler::buildGradP2);
        auto res = f[0] * q3DGradP2[0];
        for (size_t i = 1; i < n.P2; ++i) res += f[i] * q3DGradP2[i];
        return res;
    }
    GridFunctionCL<Point3DCL> getSurfGradP2(GridFunctionCL<> const & f) {
        require(qSurfGradP2[0], &LocalAssembler::buildSurfGradP2);
        auto res = f[0] * qSurfGradP2[0];
        for (size_t i = 1; i < 10; ++i) res += f[i] * qSurfGradP2[i];
        return res;
    }
    SMatrixCL<3,3> getHessP2(GridFunctionCL<> const & f) {
        SMatrixCL<3,3> hessP2[10];
        P2DiscCL::GetHessians(hessP2, T);
        auto res = f[0] * hessP2[0];
        for (size_t i = 1; i < 10 ; ++i) res += f[i] * hessP2[i];
        return res;
    }
    // builders
    void buildHatP1() {
        LocalP1CL<> hatP1[4];
        P1DiscCL::GetP1Basis(hatP1);
        for (size_t i = 0; i < 4; ++i)
            resize_and_evaluate_on_vertexes(hatP1[i], qDomain, qHatP1[i]);
    }
    void buildGradP1() {
        Point3DCL gradP1[4];
        P1DiscCL::GetGradients(gradP1, T);
        for (size_t i = 0; i < 4; ++i) {
            qGradP1[i].resize(qDomain.vertex_size());
            qGradP1[i] = gradP1[i];
            q3DGradP1[i].resize(q3Domain.vertex_size());
            q3DGradP1[i] = gradP1[i];
        }
    }
    void buildSurfGradP1() {
        require(qNormal, &LocalAssembler::buildNormal);
        require(qGradP1[0], &LocalAssembler::buildGradP1);
        for (size_t i = 0; i < 4; ++i) qSurfGradP1[i] = qGradP1[i] - dot(qGradP1[i], qNormal) * qNormal;
    }
    void buildHatP2() {
        LocalP2CL<> hatP2[10];
        P2DiscCL::GetP2Basis(hatP2);
        for (size_t i = 0; i < 10; ++i)
            resize_and_evaluate_on_vertexes(hatP2[i], qDomain, qHatP2[i]);
    }
    void buildGradP2() {
        LocalP1CL<Point3DCL> gradRefP2[10], gradP2[10];
        P2DiscCL::GetGradientsOnRef(gradRefP2);
        P2DiscCL::GetGradients(gradP2, gradRefP2, T);
        for (size_t i = 0; i < 10; ++i) {
            resize_and_evaluate_on_vertexes(gradP2[i], qDomain, qGradP2[i]);
            resize_and_evaluate_on_vertexes(gradP2[i], q3Domain, q3DGradP2[i]);
        }
    }
    void buildSurfGradP2() {
        require(qNormal, &LocalAssembler::buildNormal);
        require(qGradP2[0], &LocalAssembler::buildGradP2);
        for (size_t i = 0; i < 10; ++i) qSurfGradP2[i] = qGradP2[i] - dot(qGradP2[i], qNormal) * qNormal;
    }
    void buildNormal() {
        qNormal = getGradP2(levelSetTet);
        qLsGradNorm = sqrt(dot(qNormal, qNormal));
        qNormal = qNormal / qLsGradNorm;
        q3DNormal = get3DGradP2(levelSetTet);
        q3DNormal = q3DNormal / sqrt(dot(q3DNormal, q3DNormal));
    }
    void buildProjector() {
        require(qNormal, &LocalAssembler::buildNormal);
        qP.resize(qDomain.vertex_size());
        qP = eye<3, 3>() - outer_product(qNormal, qNormal);
    }
    void buildShapeOp() {
        require(qLsGradNorm, &LocalAssembler::buildNormal);
        require(qP, &LocalAssembler::buildProjector);
        qH.resize(qDomain.vertex_size());
        qH = getHessP2(levelSetTet);
        qH = qP * (qH / qLsGradNorm) * qP;
    }
    void buildRateOfStrainTensor() {
        require(qGradP2[0], &LocalAssembler::buildGradP2);
        require(qP, &LocalAssembler::buildProjector);
        auto qVectGrad = [&](size_t vecShapeIndex, GridFunctionCL<Point3DCL>* qGrad) {
            GridFunctionCL<SMatrixCL<3,3>> res(SMatrixCL<3,3>(), qDomain.vertex_size());
            auto && [ scaShapeIndex, row ] = ind(n.P2, vecShapeIndex);
            for (size_t i = 0; i < res.size(); ++i) {
                SMatrixCL<3, 3> mtx(0.);
                mtx.col(row, qGrad[scaShapeIndex][i]);
                assign_transpose(res[i], mtx);
            }
            return res;
        };
        for (size_t i = 0; i < n.vecP2; ++i)
            qE[i] = qP * sym_part(qVectGrad(i, qGradP2)) * qP;
    }
    void buildHatP2CrossN() { // compute velocity shape func cross normal vector
        require(qHatP2[0], &LocalAssembler::buildHatP2);
        require(qNormal, &LocalAssembler::buildNormal);
        auto qHatCrossN = [&](size_t vecShapeIndex, GridFunctionCL<>* qHat) {
            GridFunctionCL<Point3DCL> phi(Point3DCL(0., 0., 0.), qDomain.vertex_size());
            auto && [ is, in ] = ind(n.P2, vecShapeIndex);
            for (size_t i = 0; i < phi.size(); ++i)
                phi[i][in] = qHat[is][i];
            return cross_product(qNormal, phi);
        };
        for (size_t i = 0; i < n.vecP2; ++i)
            qHatP2CrossN[i] = qHatCrossN(i, qHatP2);
    }
    void buildWindNS() {
        LocalP2CL<Point3DCL> windTet;
        auto I = params.surfNavierStokesParams.w_T.RowIdx->Loc2Glo(tet);
        for (size_t i = 0; i < I.size(); ++i) {
            auto && [is, in] = ind(n.P2, i);
            windTet[is][in] = params.surfNavierStokesParams.w_T.Data[I[i]];
        }
        // windTet.assign(tet, params.surfNavierStokesParams.w_T, BndDataCL<Point3DCL>());
        resize_and_evaluate_on_vertexes(windTet, qDomain, qWind.NS);
        if (params.surfNavierStokesParams.u_N_max) {
            require(qNormal, &LocalAssembler::buildNormal);
            require(qSurfSpeed, &LocalAssembler::buildSurfSpeed);
            qWind.NS += qSurfSpeed * qNormal;
        }
    }
    void buildWindCH() {
        LocalP2CL<Point3DCL> windTet;
        auto I = params.surfCahnHilliardParams.u_T.RowIdx->Loc2Glo(tet);
        for (size_t i = 0; i < I.size(); ++i) {
            auto && [is, in] = ind(n.P2, i);
            windTet[is][in] = params.surfCahnHilliardParams.u_T.Data[I[i]];
        }
        resize_and_evaluate_on_vertexes(windTet, qDomain, qWind.CH);
    }
    void buildSurfSpeed() {
        LocalP2CL<> uNTet;
        uNTet.assign(tet, params.surfNavierStokesParams.u_N, BndDataCL<>());
        resize_and_evaluate_on_vertexes(uNTet, qDomain, qSurfSpeed);
        qSurfSpeedSurfGrad = getSurfGradP2(uNTet);
    }
    void buildF() {
        resize_and_evaluate_on_vertexes(params.surfNavierStokesParams.f_T, tet, qDomain, params.t, qF);
    }
    void buildF_c() {
        resize_and_evaluate_on_vertexes(params.surfCahnHilliardParams.f, tet, qDomain, params.t, qF_c);
    }
    void buildG() {
        resize_and_evaluate_on_vertexes(params.surfNavierStokesParams.m_g, tet, qDomain, params.t, qG);
        if (params.surfNavierStokesParams.u_N_max) {
            require(qH, &LocalAssembler::buildShapeOp);
            require(qSurfSpeed, &LocalAssembler::buildSurfSpeed);
            qG = qG + qSurfSpeed * trace(qH);
        }
    }
    void buildConcentration() {
        LocalP1CL<> chiTet;
        chiTet.assign(tet, params.surfNavierStokesParams.chi, BndDataCL<>());
        resize_and_evaluate_on_vertexes(chiTet, qDomain, qChi);
        qChiSurfGrad = getSurfGradP1(chiTet);
    }
    void buildChemPotential() {
        LocalP1CL<> omegaTet;
        omegaTet.assign(tet, params.surfNavierStokesParams.omega, BndDataCL<>());
        resize_and_evaluate_on_vertexes(omegaTet, qDomain, qOmega);
        qOmegaSurfGrad = getSurfGradP1(omegaTet);
    }
    void buildCutOffFunc() {
        require(qChi, &LocalAssembler::buildConcentration);
        qCutOffFunc = params.surfNavierStokesParams.rho_min + rho_delta * qChi;
        for (auto& val : qCutOffFunc)
            if (val < params.surfNavierStokesParams.rho_min)
                val = params.surfNavierStokesParams.rho_min;
    }
public:
    LocalAssembler(TetraCL const & tet, LocalAssemblerParams const & params) : tet(tet), params(params) {
        levelSetTet.assign(tet, params.levelSet, BndDataCL<>());
        auto const & lattice = PrincipalLatticeCL::instance(params.numbOfVirtualSubEdges);
        GridFunctionCL<> levelSetTetLat(lattice.vertex_size());
        evaluate_on_vertexes(levelSetTet, lattice, Addr(levelSetTetLat));
        SurfacePatchCL spatch;
        spatch.make_patch<MergeCutPolicyCL>(lattice, levelSetTetLat);
        make_CompositeQuad5Domain2D(qDomain, spatch, tet);
        make_SimpleQuadDomain<Quad5DataCL>(q3Domain, AllTetraC);
        GetTrafoTr(T, absDet, tet);
        absDet = std::fabs(absDet);
        rho_delta = params.surfNavierStokesParams.rho_max - params.surfNavierStokesParams.rho_min;
    }
    mtx A_vecP2vecP2() {
        require(qE[0], &LocalAssembler::buildRateOfStrainTensor);
        auto A = createMtx(n.vecP2);
        for (size_t i = 0; i < n.vecP2; ++i)
            for (size_t j = i; j < n.vecP2; ++j) {
                A[i][j] = quad_2D(contract(qE[j], qE[i]), qDomain);
                A[j][i] = A[i][j];
            }
        return A;
    }
    mtx A_consistent_vecP2vecP2() {
        require(qHatP2[0], &LocalAssembler::buildHatP2);
        require(qNormal, &LocalAssembler::buildNormal);
        require(qH, &LocalAssembler::buildShapeOp);
        require(qE[0], &LocalAssembler::buildRateOfStrainTensor);
        auto A = createMtx(n.vecP2);
        for (size_t i = 0; i < n.vecP2; ++i) {
            auto && [ is, in ] = ind(n.P2, i);
            auto e_in = std_basis<3>(in + 1);
            for (size_t j = i; j < n.vecP2; ++j) {
                auto && [ js, jn ] = ind(n.P2, j);
                auto e_jn = std_basis<3>(jn + 1);
                A[i][j] = quad_2D(contract(qE[j] - (qHatP2[js] * dot(e_jn, qNormal)) * qH, qE[i] - (qHatP2[is] * dot(e_in, qNormal)) * qH), qDomain);
                A[j][i] = A[i][j];
            }
        }
        return A;
    }
    mtx M_vecP2vecP2() {
        require(qHatP2[0], &LocalAssembler::buildHatP2);
        auto A = createMtx(n.vecP2);
        for (size_t i = 0; i < n.vecP2; ++i) {
            auto && [ is, in ] = ind(n.P2, i);
            for (size_t j = i; j < n.vecP2; ++j) {
                auto && [ js, jn ] = ind(n.P2, j);
                A[i][j] = in == jn ? quad_2D(qHatP2[js] * qHatP2[is], qDomain) : 0.;
                A[j][i] = A[i][j];
            }
        }
        return A;
    }
    mtx rho_M_vecP2vecP2() {
        if (!rho_delta) return scaleMtx(params.surfNavierStokesParams.rho_max, M_vecP2vecP2());
        require(qCutOffFunc, &LocalAssembler::buildCutOffFunc);
        require(qHatP2[0], &LocalAssembler::buildHatP2);
        auto A = createMtx(n.vecP2);
        for (size_t i = 0; i < n.vecP2; ++i) {
            auto && [ is, in ] = ind(n.P2, i);
            for (size_t j = i; j < n.vecP2; ++j) {
                auto && [ js, jn ] = ind(n.P2, j);
                A[i][j] = in == jn ? quad_2D(qCutOffFunc * qHatP2[js] * qHatP2[is], qDomain) : 0.;
                A[j][i] = A[i][j];
            }
        }
        return A;
    }
    mtx M_t_vecP2vecP2() {
        require(qHatP2[0], &LocalAssembler::buildHatP2);
        require(qP, &LocalAssembler::buildProjector);
        auto A = createMtx(n.vecP2);
        for (size_t i = 0; i < n.vecP2; ++i) {
            auto && [ is, in ] = ind(n.P2, i);
            for (size_t j = i; j < n.vecP2; ++j) {
                auto && [ js, jn ] = ind(n.P2, j);
                A[i][j] = quad_2D(qHatP2[js] * take(qP, jn, in) * qHatP2[is], qDomain);
                A[j][i] = A[i][j];
            }
        }
        return A;
    }
    mtx rho_M_t_vecP2vecP2() {
        if (!rho_delta) return scaleMtx(params.surfNavierStokesParams.rho_max, M_t_vecP2vecP2());
        require(qCutOffFunc, &LocalAssembler::buildCutOffFunc);
        require(qHatP2[0], &LocalAssembler::buildHatP2);
        require(qP, &LocalAssembler::buildProjector);
        auto A = createMtx(n.vecP2);
        for (size_t i = 0; i < n.vecP2; ++i) {
            auto && [ is, in ] = ind(n.P2, i);
            for (size_t j = i; j < n.vecP2; ++j) {
                auto && [ js, jn ] = ind(n.P2, j);
                A[i][j] = quad_2D(qCutOffFunc * qHatP2[js] * take(qP, jn, in) * qHatP2[is], qDomain);
                A[j][i] = A[i][j];
            }
        }
        return A;
    }
    mtx N_vecP2vecP2() {
        if (!params.surfNavierStokesParams.Pe) return createMtx(n.vecP2, 0.);
        require(qHatP2[0], &LocalAssembler::buildHatP2);
        require(qGradP2[0], &LocalAssembler::buildGradP2);
        require(qP, &LocalAssembler::buildProjector);
        require(qWind.NS, &LocalAssembler::buildWindNS);
        auto A = createMtx(n.vecP2);
        for (size_t i = 0; i < n.vecP2; ++i) {
            auto && [ is, in ] = ind(n.P2, i);
            for (size_t j = 0; j < n.vecP2; ++j) {
                auto && [ js, jn ] = ind(n.P2, j);
                A[i][j] = quad_2D(dot(qWind.NS, qGradP2[js]) * take(qP, in, jn) * qHatP2[is], qDomain);
            }
        }
        return A;
    }
    mtx rho_N_vecP2vecP2() {
        if (!rho_delta) return scaleMtx(params.surfNavierStokesParams.rho_max, N_vecP2vecP2());
        require(qCutOffFunc, &LocalAssembler::buildCutOffFunc);
        require(qHatP2[0], &LocalAssembler::buildHatP2);
        require(qGradP2[0], &LocalAssembler::buildGradP2);
        require(qP, &LocalAssembler::buildProjector);
        require(qWind.NS, &LocalAssembler::buildWindNS);
        auto A = createMtx(n.vecP2);
        for (size_t i = 0; i < n.vecP2; ++i) {
            auto && [ is, in ] = ind(n.P2, i);
            for (size_t j = 0; j < n.vecP2; ++j) {
                auto && [ js, jn ] = ind(n.P2, j);
                A[i][j] = quad_2D(qCutOffFunc * dot(qWind.NS, qGradP2[js]) * take(qP, in, jn) * qHatP2[is], qDomain);
            }
        }
        return A;
    }
    mtx AL_vecP2vecP2() {
        require(qE[0], &LocalAssembler::buildRateOfStrainTensor);
        auto A = createMtx(n.vecP2);
        for (size_t i = 0; i < n.vecP2; ++i)
            for (size_t j = i; j < n.vecP2; ++j) {
                A[i][j] = quad_2D(trace(qE[j]) * trace(qE[i]), qDomain);
                A[j][i] = A[i][j];
            }
        return A;
    }
    mtx H_vecP2vecP2() {
        if (!params.surfNavierStokesParams.u_N_max) return createMtx(n.vecP2, 0.);
        require(qHatP2[0], &LocalAssembler::buildHatP2);
        require(qH, &LocalAssembler::buildShapeOp);
        require(qSurfSpeed, &LocalAssembler::buildSurfSpeed);
        auto A = createMtx(n.vecP2);
        for (size_t i = 0; i < n.vecP2; ++i) {
            auto && [ is, in ] = ind(n.P2, i);
            for (size_t j = i; j < n.vecP2; ++j) {
                auto && [ js, jn ] = ind(n.P2, j);
                A[i][j] = quad_2D(qHatP2[js] * qSurfSpeed * take(qH, jn, in) * qHatP2[is], qDomain);
                A[j][i] = A[i][j];
            }
        }
        return A;
    }
    mtx S_vecP2vecP2() {
        require(qHatP2[0], &LocalAssembler::buildHatP2);
        require(qNormal, &LocalAssembler::buildNormal);
        auto A = createMtx(n.vecP2);
        for (size_t i = 0; i < n.vecP2; ++i) {
            auto && [ is, in ] = ind(n.P2, i);
            auto e_in = std_basis<3>(in + 1);
            for (size_t j = i; j < n.vecP2; ++j) {
                auto && [ js, jn ] = ind(n.P2, j);
                auto e_jn = std_basis<3>(jn + 1);
                A[i][j] = quad_2D(qHatP2[js] * dot(e_jn, qNormal) * qHatP2[is] * dot(e_in, qNormal), qDomain);
                A[j][i] = A[i][j];
            }
        }
        return A;
    }
    mtx C_n_vecP2vecP2() {
        require(q3DGradP2[0], &LocalAssembler::buildGradP2);
        require(q3DNormal, &LocalAssembler::buildNormal);
        auto A = createMtx(n.vecP2);
        for (size_t i = 0; i < n.vecP2; ++i) {
            auto && [ is, in ] = ind(n.P2, i);
            for (size_t j = i; j < n.vecP2; ++j) {
                auto && [ js, jn ] = ind(n.P2, j);
                A[i][j] = in == jn ? quad(dot(q3DNormal, q3DGradP2[js]) * dot(q3DNormal, q3DGradP2[is]), absDet, q3Domain, AllTetraC) : 0.;
                A[j][i] = A[i][j];
            }
        }
        return A;
    }
    vec F_momentum_vecP2() {
        require(qHatP2[0], &LocalAssembler::buildHatP2);
        require(qF, &LocalAssembler::buildF);
        vec b(n.vecP2);
        if (rho_delta) {
            require(qCutOffFunc, &LocalAssembler::buildCutOffFunc);
            for (size_t i = 0; i < n.vecP2; ++i) {
                auto && [is, in] = ind(n.P2, i);
                auto e_in = std_basis<3>(in + 1);
                b[i] = quad_2D(qCutOffFunc * dot(e_in, qF) * qHatP2[is], qDomain);
            }
            if (params.surfNavierStokesParams.u_N_max) {
                require(qSurfSpeed, &LocalAssembler::buildSurfSpeed);
                require(qSurfSpeedSurfGrad, &LocalAssembler::buildSurfSpeed);
                require(qE[0], &LocalAssembler::buildRateOfStrainTensor);
                require(qH, &LocalAssembler::buildShapeOp);
                for (size_t i = 0; i < n.vecP2; ++i) {
                    auto && [is, in] = ind(n.P2, i);
                    auto e_in = std_basis<3>(in + 1);
                    b[i] += quad_2D(qCutOffFunc * qSurfSpeed * dot(e_in, qSurfSpeedSurfGrad) * qHatP2[is] - params.surfNavierStokesParams.nu * qSurfSpeed * contract(qH, qE[i]), qDomain);
                }
            }
        }
        else {
            for (size_t i = 0; i < n.vecP2; ++i) {
                auto && [is, in] = ind(n.P2, i);
                auto e_in = std_basis<3>(in + 1);
                b[i] = params.surfNavierStokesParams.rho_max * quad_2D(dot(e_in, qF) * qHatP2[is], qDomain);
            }
            if (params.surfNavierStokesParams.u_N_max) {
                require(qSurfSpeed, &LocalAssembler::buildSurfSpeed);
                require(qSurfSpeedSurfGrad, &LocalAssembler::buildSurfSpeed);
                require(qE[0], &LocalAssembler::buildRateOfStrainTensor);
                require(qH, &LocalAssembler::buildShapeOp);
                for (size_t i = 0; i < n.vecP2; ++i) {
                    auto && [is, in] = ind(n.P2, i);
                    auto e_in = std_basis<3>(in + 1);
                    b[i] += quad_2D(params.surfNavierStokesParams.rho_max * qSurfSpeed * dot(e_in, qSurfSpeedSurfGrad) * qHatP2[is] - params.surfNavierStokesParams.nu * qSurfSpeed * contract(qH, qE[i]), qDomain);
                }
            }
        }
        if (params.surfNavierStokesParams.lineTension) {
            require(qChi, &LocalAssembler::buildConcentration);
            require(qChiSurfGrad, &LocalAssembler::buildConcentration);
            require(qOmega, &LocalAssembler::buildChemPotential);
            for (size_t i = 0; i < n.vecP2; ++i) {
                auto && [is, in] = ind(n.P2, i);
                auto e_in = std_basis<3>(in + 1);
                b[i] -= params.surfNavierStokesParams.lineTension * quad_2D(qOmega * dot(e_in, qChiSurfGrad) * qHatP2[is], qDomain);
            }
        }
        if (params.surfNavierStokesParams.gamma) {
            require(qG, &LocalAssembler::buildG);
            require(qE[0], &LocalAssembler::buildRateOfStrainTensor);
            for (size_t i = 0; i < n.vecP2; ++i)
                b[i] -= params.surfNavierStokesParams.gamma * quad_2D(qG * trace(qE[i]), qDomain);
        }
        return b;
    }
    mtx B_P1vecP2() {
        require(qSurfGradP1[0], &LocalAssembler::buildSurfGradP1);
        require(qHatP2[0], &LocalAssembler::buildHatP2);
        auto A = createMtx(n.P1, n.vecP2);
        for (size_t i = 0; i < n.P1; ++i)
            for (size_t j = 0; j < n.vecP2; ++j) {
                auto && [ js, jn ] = ind(n.P2, j);
                auto e_jn = std_basis<3>(jn + 1);
                A[i][j] = quad_2D(qHatP2[js] * dot(e_jn, qSurfGradP1[i]), qDomain);
            }
        return A;
    }
    mtx Q_P1vecP2() { // rhs-curl-projection mtx
        require(qSurfGradP1[0], &LocalAssembler::buildSurfGradP1);
        require(qHatP2CrossN[0], &LocalAssembler::buildHatP2CrossN);
        auto A = createMtx(n.P1, n.vecP2);
        for (size_t i = 0; i < n.P1; ++i)
            for (size_t j = 0; j < n.vecP2; ++j)
                A[i][j] = quad_2D(dot(qHatP2CrossN[j], qSurfGradP1[i]), qDomain);
        return A;
    }
    mtx A_P1P1() {
        require(qSurfGradP1[0], &LocalAssembler::buildSurfGradP1);
        auto A = createMtx(n.P1);
        for (size_t i = 0; i < n.P1; ++i)
            for (size_t j = i; j < n.P1; ++j) {
                A[i][j] = quad_2D(dot(qSurfGradP1[j], qSurfGradP1[i]), qDomain);
                A[j][i] = A[i][j];
            }
        return A;
    }
    mtx C_n_P1P1() {
        require(q3DGradP1[0], &LocalAssembler::buildGradP1);
        require(q3DNormal, &LocalAssembler::buildNormal);
        auto A = createMtx(n.P1);
        for (size_t i = 0; i < n.P1; ++i)
            for (size_t j = i; j < n.P1; ++j) {
                A[i][j] = quad(dot(q3DNormal, q3DGradP1[j]) * dot(q3DNormal, q3DGradP1[i]), absDet, q3Domain, AllTetraC);
                A[j][i] = A[i][j];
            }
        return A;
    }
    mtx C_full_P1P1() {
        require(q3DGradP1[0], &LocalAssembler::buildGradP1);
        auto A = createMtx(n.P1);
        for (size_t i = 0; i < n.P1; ++i)
            for (size_t j = i; j < n.P1; ++j) {
                A[i][j] = quad(dot(q3DGradP1[j], q3DGradP1[i]), absDet, q3Domain, AllTetraC);
                A[j][i] = A[i][j];
            }
        return A;
    }
    mtx M_P1P1() {
        require(qHatP1[0], &LocalAssembler::buildHatP1);
        auto A = createMtx(n.P1);
        for (size_t i = 0; i < n.P1; ++i)
            for (size_t j = i; j < n.P1; ++j) {
                A[i][j] = quad_2D(qHatP1[j] * qHatP1[i], qDomain);
                A[j][i] = A[i][j];
            }
        return A;
    }
    mtx N_P1P1() {
        require(qHatP1[0], &LocalAssembler::buildHatP1);
        require(qSurfGradP1[0], &LocalAssembler::buildSurfGradP1);
        require(qWind.CH, &LocalAssembler::buildWindCH);
        auto A = createMtx(n.P1);
        for (size_t i = 0; i < n.P1; ++i)
            for (size_t j = 0; j < n.P1; ++j)
                A[i][j] = quad_2D(dot(qWind.CH, qSurfGradP1[j]) * qHatP1[i], qDomain);
        return A;
    }
    vec F_continuity_P1() {
        require(qHatP1[0], &LocalAssembler::buildHatP1);
        require(qG, &LocalAssembler::buildG);
        vec b(n.P1);
        for (size_t i = 0; i < n.P1; ++i)
            b[i] = quad_2D(qG * qHatP1[i], qDomain);
        return b;
    }
    vec F_concentration_P1(){
        require(qHatP1[0], &LocalAssembler::buildHatP1);
        require(qF_c, &LocalAssembler::buildF_c);
        vec b(n.P1);
        for (size_t i = 0; i < n.P1; ++i)
            b[i] = quad_2D(qF_c * qHatP1[i], qDomain);
        return b;
    }
};

struct FEMatDescCL : MatDescCL {
    LocalAssembler::BilinearForm form;
    FEMatDescCL(IdxDescCL* r, IdxDescCL* c, LocalAssembler::BilinearForm form) : MatDescCL(r, c), form(form) {}
};

struct FEVecDescCL : VecDescCL {
    LocalAssembler::LinearForm form;
    FEVecDescCL(IdxDescCL* r, LocalAssembler::LinearForm form) : VecDescCL(r), form(form) {}
};

struct FESystem {
    std::vector<FEMatDescCL*> matrices;
    std::vector<FEVecDescCL*> vectors;
    LocalAssembler::LocalAssemblerParams params;
};

void setupFESystem(MultiGridCL const &, FESystem&);

void SetupCahnHilliardIF_P1P1( const MultiGridCL& MG_,  MatDescCL* M_P1, MatDescCL* NormalStab_P1, MatDescCL* TangentStab_P1, MatDescCL* VolumeStab_P1, MatDescCL* L_P1P1 ,MatDescCL* LM_P1P1 ,MatDescCL* Gprimeprime_P1P1 , const VecDescCL& lset, const LsetBndDataCL& lset_bnd, const VecDescCL& velocity, const BndDataCL<Point3DCL>& velocity_bnd,const VecDescCL& volume_fraction, const BndDataCL<>& volume_fraction_bnd);

    double Mobility_function(double x,double t=0);
    double Diffusion_function(double x,double t=0);
    double Density_function(double x,double t=0);


    double inverse_square_root(double x);
    double Potential_function(double x);
    double Potential_prime_function(double x);
    double Potential_prime_convex_function(double x);
    double Potential_prime_concave_function(double x);

void SetupInterfaceVectorRhsP1 (const MultiGridCL& mg, VecDescCL* v,
                                const VecDescCL& ls, const BndDataCL<>& lsbnd, InstatVectorFunction f, double t = 0.);
void SetupInterfaceVectorRhsP2 (const MultiGridCL& mg, VecDescCL* v,
                                const VecDescCL& ls, const BndDataCL<>& lsbnd, InstatVectorFunction f);

class LocalLaplaceBeltramiP1CL
{
  private:
    double D_; // diffusion coefficient

    Point3DCL grad[4];
    double dummy;
    GridFunctionCL<Point3DCL> n,
                              q[4];
    std::valarray<double> absdet;

  public:
    static const FiniteElementT row_fe_type= P1IF_FE,
                                col_fe_type= P1IF_FE;

    double coup[4][4];

    void setup (const TetraCL& t, const InterfaceCommonDataP1CL& cdata) {
        n.resize( cdata.surf.facet_size());
        absdet.resize( cdata.surf.facet_size());
        if (cdata.surf.normal_empty())
            cdata.surf.compute_normals( t);
        std::copy( cdata.surf.normal_begin(), cdata.surf.normal_end(), sequence_begin( n));
        std::copy( cdata.surf.absdet_begin(), cdata.surf.absdet_end(), sequence_begin( absdet));
        P1DiscCL::GetGradients( grad, dummy, t);
        for(int i= 0; i < 4; ++i) {
            q[i].resize( cdata.surf.facet_size());
            q[i]= grad[i] - dot( grad[i], n)*n;
        }
        for (int i= 0; i < 4; ++i) {
            coup[i][i]= D_* 0.5*(dot(q[i], q[i])*absdet).sum();
            for(int j= 0; j < i; ++j)
                coup[i][j]= coup[j][i]= D_* 0.5*(dot(q[i], q[j])*absdet).sum();
        }
    }

    LocalLaplaceBeltramiP1CL (double D)
        :D_( D) {}
};

    template <typename DiscVelSolT>
    class LocalLaplaceMobilityP1CL
    {
    private:

        QuadDomain2DCL qdom;
        double time_;
        InstatVectorFunction normal_;

        std::valarray<double> mobility;
        std::valarray<double> qmobility;

        const DiscVelSolT concentr_;

        LocalP1CL<double> concentr_loc;
        GridFunctionCL<> qconcentr;

        LocalP1CL<double> P1Hat[4];

        Point3DCL grad[4];
        double dummy;
        GridFunctionCL<Point3DCL> n, q[4], qq[4];
        LocalP1CL<Point3DCL> Normals;

        GridFunctionCL<Point3DCL> qnormal;

        std::valarray<double> absdet;

    public:
        static const FiniteElementT row_fe_type= P1IF_FE,
                col_fe_type= P1IF_FE;

        double coup[4][4];

        void setup (const TetraCL& t, const InterfaceCommonDataP1CL& cdata) {

            make_CompositeQuad5Domain2D( qdom, cdata.surf, t);
            concentr_loc.assign( t, concentr_);
            //resize_and_evaluate_on_vertexes( concentr_loc, qdom, qconcentr);

            P1DiscCL::GetGradients( grad, dummy, t);

            //qnormal.assign(t, normal_, time_);
            resize_and_evaluate_on_vertexes( normal_, t, qdom, time_, qnormal);

            /*// Scale Normals accordingly to the Euclidean Norm (only consider the ones which make a contribution in the sense of them being big enough... otherwise one has to expect problems with division through small numbers)
            for(Uint i=0; i<qnormal.size(); ++i) {
                //if(qnormal[i].norm()> 1e-8)
                qnormal[i]= qnormal[i]/qnormal[i].norm();
            }*/

            for(int j=0; j<4 ;++j) {
                qq[j].resize( qdom.vertex_size());
                qq[j]= grad[j];
                qq[j]-= dot( qq[j], qnormal)*qnormal;
            }

            LocalP1CL<> mobility;
            for(int i=0; i<4 ; ++i)
            {
                mobility += Mobility_function(concentr_loc[i], time_)*cdata.p1[i];
            }
            resize_and_evaluate_on_vertexes (mobility, qdom, qmobility);


            for (int i= 0; i < 4; ++i)
                for(int j= 0; j < 4; ++j) {
                    coup[i][j]= quad_2D( qmobility*dot(qq[i],qq[j]), qdom);
                }

        }

        LocalLaplaceMobilityP1CL (const DiscVelSolT& conc, InstatVectorFunction normal, double t)
                :concentr_(conc), normal_(normal), time_(t) {}
    };

    template <typename DiscVelSolT>
    class LocalLaplaceNonlinearP1CL
    {
    private:

        QuadDomain2DCL qdom;
        double time_;
        InstatVectorFunction normal_;


        std::valarray<double> mobility;
        std::valarray<double> qmobility;

        const DiscVelSolT concentr_;

        LocalP1CL<double> concentr_loc;
        GridFunctionCL<> qconcentr;

        LocalP1CL<double> P1Hat[4];

        Point3DCL grad[4];
        double dummy;
        GridFunctionCL<Point3DCL> n, q[4], qq[4];
        LocalP1CL<Point3DCL> Normals;

        GridFunctionCL<Point3DCL> qnormal;

        std::valarray<double> absdet;

    public:
        static const FiniteElementT row_fe_type= P1IF_FE,
                col_fe_type= P1IF_FE;

        double coup[4][4];

        void setup (const TetraCL& t, const InterfaceCommonDataP1CL& cdata) {

            make_CompositeQuad5Domain2D( qdom, cdata.surf, t);
            concentr_loc.assign( t, concentr_);
            //resize_and_evaluate_on_vertexes( concentr_loc, qdom, qconcentr);

            P1DiscCL::GetGradients( grad, dummy, t);

            //qnormal.assign(t, normal_, time_);
            resize_and_evaluate_on_vertexes( normal_, t, qdom, time_, qnormal);

            /*// Scale Normals accordingly to the Euclidean Norm (only consider the ones which make a contribution in the sense of them being big enough... otherwise one has to expect problems with division through small numbers)
            for(Uint i=0; i<qnormal.size(); ++i) {
                //if(qnormal[i].norm()> 1e-8)
                qnormal[i]= qnormal[i]/qnormal[i].norm();
            }*/

            for(int j=0; j<4 ;++j) {
                qq[j].resize( qdom.vertex_size());
                qq[j]= grad[j];
                qq[j]-= dot( qq[j], qnormal)*qnormal;
            }

            LocalP1CL<> mobility;
            for(int i=0; i<4 ; ++i)
            {
                mobility += Diffusion_function(concentr_loc[i], time_)*cdata.p1[i];
            }
            resize_and_evaluate_on_vertexes (mobility, qdom, qmobility);


            for (int i= 0; i < 4; ++i)
                for(int j= 0; j < 4; ++j) {
                    coup[i][j]= quad_2D( qmobility*dot(qq[i],qq[j]), qdom);
                }

        }

        LocalLaplaceNonlinearP1CL (const DiscVelSolT& conc, InstatVectorFunction normal, double t)
                :concentr_(conc), normal_(normal), time_(t)  {}
    };

    class LocalFullGradientsP1CL
    {
    private:

        QuadDomain2DCL qdom;
        double time_;
        InstatVectorFunction normal_;

        std::valarray<double> mobility;
        std::valarray<double> qmobility;


        LocalP1CL<double> concentr_loc;
        GridFunctionCL<> qconcentr;

        LocalP1CL<double> P1Hat[4];
        Quad5CL<double> U_Grad[4];

        Point3DCL grad[4];
        double dummy;
        GridFunctionCL<Point3DCL> n, q[4], qq[4];
        LocalP1CL<Point3DCL> Normals;

        GridFunctionCL<Point3DCL> qnormal;

        double absdet;

    public:
        static const FiniteElementT row_fe_type= P1IF_FE,
                col_fe_type= P1IF_FE;

        double coup[4][4];
        void setup (const TetraCL& tet, const NarrowBandCommonDataP1CL& bdata) {

            P1DiscCL::GetGradients( grad, dummy, tet);

            //qnormal.assign(tet,normal_,time_,Quad5DataCL::Node);

            //for(int i=0; i<4; ++i)
            //    U_Grad[i]=dot( qnormal, Quad5CL<Point3DCL>( grad[i]));

            absdet=std::abs(dummy);

            dummy=std::pow(absdet,1./3); //of order h--meshsize of the tetra hedra~~!!!!
            //std::cout<<"dummy  "<<dummy<<std::endl;
            for (int i= 0; i < 4; ++i) {
                for(int j= 0; j <= i; ++j)
                {
                    Quad5CL<double> res3( dot(Quad5CL<Point3DCL>(grad[i]),Quad5CL<Point3DCL>(grad[j])));
                    coup[i][j]= coup[j][i]= 1.0*res3.quad(absdet/6);//dummy*
                    // 	std::cout<<coup[i][j]<<std::endl;
                }
            }
            // std::cin>>dummy;
        }
        void setup (const TetraCL& tet, const InterfaceCommonDataP1CL& cdata) {

            P1DiscCL::GetGradients( grad, dummy, tet);

            //qnormal.assign(tet,normal_,time_,Quad5DataCL::Node);

            //for(int i=0; i<4; ++i)
            //    U_Grad[i]=dot( qnormal, Quad5CL<Point3DCL>( grad[i]));

            absdet=std::abs(dummy);

            dummy=std::pow(absdet,1./3); //of order h--meshsize of the tetra hedra~~!!!!
            //std::cout<<"dummy  "<<dummy<<std::endl;
            for (int i= 0; i < 4; ++i) {
                for(int j= 0; j <= i; ++j)
                {
                    Quad5CL<double> res3( dot(Quad5CL<Point3DCL>(grad[i]),Quad5CL<Point3DCL>(grad[j])));
                    coup[i][j]= coup[j][i]= 1.0*res3.quad(absdet/6);//dummy*
                    // 	std::cout<<coup[i][j]<<std::endl;
                }
            }
            // std::cin>>dummy;
        }
        /*void setup (const TetraCL& t, const InterfaceCommonDataP1CL& cdata) {

            QuadDomainCL q3Domain;

            make_SimpleQuadDomain<Quad5DataCL> (q3Domain, AllTetraC);

            make_CompositeQuad5Domain2D( qdom, cdata.surf, t);
            //resize_and_evaluate_on_vertexes( concentr_loc, qdom, qconcentr);

            P1DiscCL::GetGradients( grad, dummy, t);

            //qnormal.assign(t, normal_, time_);
            resize_and_evaluate_on_vertexes( normal_, t, qdom, time_, qnormal);

            // Scale Normals accordingly to the Euclidean Norm (only consider the ones which make a contribution in the sense of them being big enough... otherwise one has to expect problems with division through small numbers)
            //for(Uint i=0; i<qnormal.size(); ++i) {
                //if(qnormal[i].norm()> 1e-8)
            //    qnormal[i]= qnormal[i]/qnormal[i].norm();
            //}

            for(int j=0; j<4 ;++j) {
                qq[j].resize( qdom.vertex_size());
                qq[j]= grad[j];
                qq[j]-= dot( qq[j], qnormal)*qnormal;
            }

            for (int i= 0; i < 4; ++i)
                for(int j= 0; j < 4; ++j) {
                    coup[i][j]= quad( dot(grad[i],grad[j]), q3Domain);
                }*//*
            for(int i=0; i<4; ++i)
                U_Grad[i]=dot( qnormal, Quad5CL<Point3DCL>( grad[i]));

            absdet=std::abs(dummy);

            dummy=std::pow(absdet,1./3); //of order h--meshsize of the tetra hedra~~!!!!
            //std::cout<<"dummy  "<<dummy<<std::endl;
            for (int i= 0; i < 4; ++i) {
                for(int j= 0; j <= i; ++j)
                {
                    Quad5CL<double> res3( U_Grad[i] * U_Grad[j]);
                    coup[i][j]= coup[j][i]= res3.quad(absdet/6);//dummy*
                    // 	std::cout<<coup[i][j]<<std::endl;
                }
            }

        }*/

        LocalFullGradientsP1CL ( InstatVectorFunction normal, double time)
                :normal_(normal), time_(time) {}
    };

/// \brief The routine sets up the Laplace-matrix in mat in bulk element in a narrow band near the interface defined by ls.
///        It belongs to the FE induced by standard P1-elements.
///        Notmal gradient is used for stabilization in solving surface diffusion equation
/// D is the diffusion-coefficient

    class LocalNormalLaplaceBulkP1CL
    {
    private:
        double D_; // diffusion coefficient- stabilization parameter?
        double dt_;
        InstatVectorFunction normal_;
        double time_;

        Point3DCL grad[4];
        double dummy;
        double absdet;
        //GridFunctionCL<Point3DCL> q[4];
        Quad5CL<double> U_Grad[4];
        Quad5CL<Point3DCL> qnormal;

    public:
        static const FiniteElementT row_fe_type= P1IF_FE,
                col_fe_type= P1IF_FE;
        double coup[4][4];

        void setup (const TetraCL& tet, const NarrowBandCommonDataP1CL& cdata) {

            P1DiscCL::GetGradients( grad, dummy, tet);

            qnormal.assign(tet,normal_,time_,Quad5DataCL::Node);

            for(int i=0; i<4; ++i)
                U_Grad[i]=dot( qnormal, Quad5CL<Point3DCL>( grad[i]));

            absdet=std::abs(dummy);

            dummy=std::pow(absdet,1./3); //of order h--meshsize of the tetra hedra~~!!!!
            //  std::cout<<"dummy  "<<dummy<<std::endl;
            for (int i= 0; i < 4; ++i) {
                for(int j= 0; j <= i; ++j)
                {
                    Quad5CL<double> res3( U_Grad[i] * U_Grad[j]);
                    coup[i][j]= coup[j][i]=D_*//(1.0+1.0/(dummy+dt_))*
                            res3.quad(absdet/6);//D_*(D_/dummy+dummy/dt_+0.2)*// D_*
                    //  	std::cout<<i<<" "<<j<<" : "<<coup[i][j]<<"  "<<coup[j][i]<<" ; ";
                }
                //   std::cout<<std::endl;
            }
            // std::cin>>dummy;
        }
        void setup (const TetraCL& tet, const InterfaceCommonDataP1CL& cdata) {

            P1DiscCL::GetGradients( grad, dummy, tet);

            qnormal.assign(tet,normal_,time_,Quad5DataCL::Node);

            for(int i=0; i<4; ++i)
                U_Grad[i]=dot( qnormal, Quad5CL<Point3DCL>( grad[i]));

            absdet=std::abs(dummy);

            dummy=std::pow(absdet,1./3); //of order h--meshsize of the tetra hedra~~!!!!
            //std::cout<<"dummy  "<<dummy<<std::endl;
            for (int i= 0; i < 4; ++i) {
                for(int j= 0; j <= i; ++j)
                {
                    Quad5CL<double> res3( U_Grad[i] * U_Grad[j]);
                    coup[i][j]= coup[j][i]= D_*res3.quad(absdet/6);//dummy*
                    // 	std::cout<<coup[i][j]<<std::endl;
                }
            }
            // std::cin>>dummy;
        }
        LocalNormalLaplaceBulkP1CL (double D,double dt,InstatVectorFunction normal, double t)
                :D_( D),dt_(dt),normal_(normal),time_(t) {}
    };


/// \brief The routine sets up the convection-matrix in mat on the interface defined by ls.
///        It belongs to the FE induced by standard P1-elements.
///
/// The template-parameter is only used to circumvent the exact type of the discrete
/// velocity solution in the Stokes classes.
template <class DiscVelSolT>
void SetupConvectionP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsbnd, const DiscVelSolT& v);

template <class DiscVelSolT>
class LocalInterfaceConvectionP1CL
{
  private:
    const DiscVelSolT w_; // wind

    QuadDomain2DCL qdom;
    Point3DCL grad[4];
    double dummy;
    LocalP2CL<Point3DCL> w_loc;
    std::valarray<double> q[4];
    GridFunctionCL<Point3DCL> qw;

  public:
    static const FiniteElementT row_fe_type= P1IF_FE,
                                col_fe_type= P1IF_FE;

    double coup[4][4];

    void setup (const TetraCL& t, const InterfaceCommonDataP1CL& cdata);
    void setup (const TetraCL& t, const NarrowBandCommonDataP1CL& cdata);


    LocalInterfaceConvectionP1CL (const DiscVelSolT& w)
        :  w_( w) {}
};
//used to penalize gradients in the direction of ambient velocity field
    template <class DiscVelSolT>
    class LocalInterfaceVelocityLaplaceP1CL
    {
    private:
        const DiscVelSolT w_; // wind
        InstatVectorFunction normal_;
        double time_;
        Quad5CL<Point3DCL> qnormal;
        QuadDomain2DCL qdom;
        Point3DCL grad[4];
        double dummy;
        double absdet;
        Quad5CL<Point3DCL> w_loc;
        Quad5CL<double> U_Gradw[4];
        std::valarray<double> q[4];
        GridFunctionCL<Point3DCL> qw;

    public:
        static const FiniteElementT row_fe_type= P1IF_FE,
                col_fe_type= P1IF_FE;

        double coup[4][4];

        void setup (const TetraCL& t, const InterfaceCommonDataP1CL& cdata);
        void setup (const TetraCL& t, const NarrowBandCommonDataP1CL& cdata);

        LocalInterfaceVelocityLaplaceP1CL (const DiscVelSolT& w, InstatVectorFunction normal, double time)
                :  w_( w), normal_(normal), time_(time) {}
    };


    template <class DiscVelSolT>
    void LocalInterfaceVelocityLaplaceP1CL<DiscVelSolT>::setup (const TetraCL& tet, const InterfaceCommonDataP1CL& cdata)
    {
        //make_CompositeQuad5Domain2D( qdom, bdata.surf, t);
        w_loc.assign( tet, w_);
        //resize_and_evaluate_on_vertexes( w_loc, qdom, qw);
//
//        P1DiscCL::GetGradients( grad, dummy, t);
//        for (int i= 0; i < 4; ++i)
//            resize_and_evaluate_on_vertexes( cdata.p1[i], qdom, q[i]);
//
//        for (int i= 0; i < 4; ++i)
//            for(int j= 0; j < 4; ++j)
//                coup[i][j]= quad_2D( dot( grad[j], qw)*q[i], qdom);
//
        ////////

        P1DiscCL::GetGradients( grad, dummy, tet);
        qnormal.assign(tet,normal_,time_);

        GridFunctionCL<double> inv_norm = dot(w_loc,w_loc).apply(inverse_square_root);
//        for (int i=0;i<norm.size();i++) norm.apply(sqr)

        for(int i=0; i<4; ++i) {
            U_Gradw[i] = dot(qnormal, Quad5CL<Point3DCL>(grad[i]));

        }
        absdet=std::abs(dummy);

        dummy=std::pow(absdet,1./3); //of order h--meshsize of the tetra hedra~~!!!!
        //  std::cout<<"dummy  "<<dummy<<std::endl;
        for (int i= 0; i < 4; ++i) {
            for(int j= 0; j <= i; ++j)
            {
                Quad5CL<double> res3( U_Gradw[i] * U_Gradw[j]);
                coup[i][j]= coup[j][i]=//(1.0+1.0/(dummy+dt_))*
                        res3.quad(absdet/6);//D_*(D_/dummy+dummy/dt_+0.2)*// D_*
                //  	std::cout<<i<<" "<<j<<" : "<<coup[i][j]<<"  "<<coup[j][i]<<" ; ";
            }
            //   std::cout<<std::endl;
        }
    }

    template <class DiscVelSolT>
    void LocalInterfaceVelocityLaplaceP1CL<DiscVelSolT>::setup (const TetraCL& tet, const NarrowBandCommonDataP1CL& bdata)
    {
        //make_CompositeQuad5Domain2D( qdom, bdata.surf, t);
        w_loc.assign( tet, w_);
        //resize_and_evaluate_on_vertexes( w_loc, qdom, qw);
//
//        P1DiscCL::GetGradients( grad, dummy, t);
//        for (int i= 0; i < 4; ++i)
//            resize_and_evaluate_on_vertexes( cdata.p1[i], qdom, q[i]);
//
//        for (int i= 0; i < 4; ++i)
//            for(int j= 0; j < 4; ++j)
//                coup[i][j]= quad_2D( dot( grad[j], qw)*q[i], qdom);
//
        ////////

        P1DiscCL::GetGradients( grad, dummy, tet);
        qnormal.assign(tet,normal_,time_,Quad5DataCL::Node);

        GridFunctionCL<double> inv_norm = dot(w_loc,w_loc).apply(inverse_square_root);
//        for (int i=0;i<norm.size();i++) norm.apply(sqr)

        for(int i=0; i<4; ++i) {
            U_Gradw[i] = 0.1*dot(qnormal, Quad5CL<Point3DCL>(grad[i]));
            U_Gradw[i] += inv_norm*dot(w_loc, Quad5CL<Point3DCL>(grad[i]));

        }

        absdet=std::abs(dummy);

        dummy=std::pow(absdet,1./3); //of order h--meshsize of the tetra hedra~~!!!!
        //  std::cout<<"dummy  "<<dummy<<std::endl;
        for (int i= 0; i < 4; ++i) {
            for(int j= 0; j <= i; ++j)
            {
                Quad5CL<double> res3( U_Gradw[i] * U_Gradw[j]);
                coup[i][j]= coup[j][i]=//(1.0+1.0/(dummy+dt_))*
                        res3.quad(absdet/6);//D_*(D_/dummy+dummy/dt_+0.2)*// D_*
                //  	std::cout<<i<<" "<<j<<" : "<<coup[i][j]<<"  "<<coup[j][i]<<" ; ";
            }
            //   std::cout<<std::endl;
        }
    }

    /// \brief The routine sets up the mass-matrix scaled with \f$ div_\Gamma v \f$ in mat on the interface
///        defined by ls. It belongs to the FE induced by standard P1-elements.
///
/// The template-parameter is only used to circumvent the exact type of the discrete
/// velocity solution in the Stokes classes.
template <class DiscVelSolT>
void SetupMassDivP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsbnd, const DiscVelSolT& w);

template <typename DiscVelSolT>
class LocalInterfaceMassDivP1CL
{
  private:
    const DiscVelSolT w_;

    QuadDomain2DCL qdom;
    LocalP2CL<Point3DCL> w_loc;
    std::valarray<double> q[4];
    double dummy;
    SMatrixCL<3,3> T;
    GridFunctionCL<Point3DCL> n,
                              qgradp2i;
    std::valarray<double> qdivgamma_w;
    LocalP1CL<Point3DCL> gradrefp2[10],
                         gradp2[10];


  public:
    static const FiniteElementT row_fe_type= P1IF_FE,
                                col_fe_type= P1IF_FE;

    double coup[4][4];

    void setup (const TetraCL& t, const InterfaceCommonDataP1CL& cdata);

    LocalInterfaceMassDivP1CL (const DiscVelSolT& w)
        : w_( w) { P2DiscCL::GetGradientsOnRef( gradrefp2); }
};

    template <typename DiscVelSolT>
    class LocalInterfaceMassRhoP1CL
    {
    private:
        const DiscVelSolT w_;
        double time_;
        QuadDomain2DCL qdom;
        LocalP1CL<Point3DCL> w_loc;
        std::valarray<double> q[4];
        double dummy;
        std::valarray<double> qdensity;

        //SMatrixCL<3,3> T;
        //GridFunctionCL<Point3DCL> n,
               // qgradp2i;
        const DiscVelSolT concentr_;
        LocalP1CL<double> concentr_loc;

        //LocalP1CL<Point3DCL> gradrefp2[10],
        //        gradp2[10];


    public:
        static const FiniteElementT row_fe_type= P1IF_FE,
                col_fe_type= P1IF_FE;

        double coup[4][4];

        void setup (const TetraCL& t, const InterfaceCommonDataP1CL& cdata);

        LocalInterfaceMassRhoP1CL(const DiscVelSolT& conc, InstatVectorFunction normal, double t)
        :concentr_(conc), time_(t){}

    };

/// \brief Compute the P2 load vector corresponding to the function f on a single tetra.
class LocalVectorP2CL
{
  private:
    InstatScalarFunction f_;
    double time_;

    std::valarray<double> qp2,
                          qf;

  public:
    static const FiniteElementT row_fe_type= P2IF_FE;

    double vec[10];

    LocalVectorP2CL (InstatScalarFunction f, double time) : f_( f), time_( time) {}

    void setup (const TetraCL&, const InterfaceCommonDataP2CL& cdata, const IdxT numr[10]) {
        resize_and_evaluate_on_vertexes( f_, cdata.qdom_projected.vertexes(), time_, qf);
        qp2.resize( cdata.qdom.vertex_size());
        for (Uint i= 0; i < 10; ++i) {
                if (numr[i] == NoIdx)
                    continue;
                evaluate_on_vertexes( cdata.p2[i], cdata.qdom, Addr( qp2));
                vec[i]= quad_2D( cdata.qdom_projected.absdets()*qf*qp2, cdata.qdom);
        }
    }
};


/// \brief Trafo of the interfacial gradient on the linear interface to the quadratic iface under a QuaQuaMapperCL.
/// Computes W from La. 5.1 of the high order paper. The transformation of the gradient in the discretization requires W^{-1}.
void gradient_trafo (const TetraCL& tet, const BaryCoordCL& xb, const QuaQuaMapperCL& quaqua, const SurfacePatchCL& p, SMatrixCL<3,3>& W);

class LocalLaplaceBeltramiP2CL
{
  private:
    double D_; // diffusion coefficient

    LocalP1CL<Point3DCL> gradp2[10];
    GridFunctionCL<Point3DCL> qgradp2[10];

    GridFunctionCL<Point3DCL> nl;
    GridFunctionCL<SMatrixCL<3,3> > Winv;

  public:
    static const FiniteElementT row_fe_type= P2IF_FE,
                                col_fe_type= P2IF_FE;

    double coup[10][10];

    const GridFunctionCL<Point3DCL>& get_qgradp2 (size_t i) { return qgradp2[i]; }

    void setup (const TetraCL& t, const InterfaceCommonDataP2CL& cdata) {
        if (cdata.surf.normal_empty())
            cdata.surf.compute_normals( t);
        resize_and_scatter_piecewise_normal( cdata.surf, cdata.qdom, nl);

        Winv.resize( cdata.qdom.vertex_size());
        QRDecompCL<3,3> qr;
        SVectorCL<3> tmp;
        for (Uint i= 0; i < cdata.qdom.vertex_size(); ++i) {
            gradient_trafo( t, cdata.qdom.vertex_begin()[i], cdata.quaqua, cdata.surf, qr.GetMatrix());
            qr.prepare_solve();
            for (Uint j= 0; j < 3; ++j) {
                tmp= std_basis<3>( j + 1);
                qr.Solve( tmp);
                Winv[i].col( j, tmp);
            }
        }

        double dummy;
        SMatrixCL<3,3> T;
        GetTrafoTr( T, dummy, t);
        P2DiscCL::GetGradients( gradp2, cdata.gradrefp2, T);
        for (int i= 0; i < 10; ++i) {
            resize_and_evaluate_on_vertexes ( gradp2[i], cdata.qdom, qgradp2[i]);
            for (Uint j= 0; j < qgradp2[i].size(); ++j) {
                tmp=  qgradp2[i][j] - inner_prod( nl[j], qgradp2[i][j])*nl[j];
                qgradp2[i][j]= Winv[j]*tmp;
            }
        }

        for (int i= 0; i < 10; ++i) {
            coup[i][i]= quad_2D( cdata.qdom_projected.absdets()*dot( qgradp2[i], qgradp2[i]), cdata.qdom);
            for(int j= 0; j < i; ++j)
                coup[i][j]= coup[j][i]= quad_2D( cdata.qdom_projected.absdets()*dot( qgradp2[j], qgradp2[i]), cdata.qdom);
        }
    }

    LocalLaplaceBeltramiP2CL (double D)
        :D_( D) {}
};

class LocalMassP2CL
{
  private:
    std::valarray<double> qp2[10];

  public:
    static const FiniteElementT row_fe_type= P2IF_FE,
                                col_fe_type= P2IF_FE;

    double coup[10][10];

    void setup (const TetraCL&, const InterfaceCommonDataP2CL& cdata) {
        for (int i= 0; i < 10; ++i)
            resize_and_evaluate_on_vertexes ( cdata.p2[i], cdata.qdom, qp2[i]);

        for (int i= 0; i < 10; ++i) {
            coup[i][i]= quad_2D( cdata.qdom_projected.absdets()*qp2[i]*qp2[i], cdata.qdom);
            for(int j= 0; j < i; ++j)
                coup[i][j]= coup[j][i]= quad_2D( cdata.qdom_projected.absdets()*qp2[j]*qp2[i], cdata.qdom);
        }
    }

    LocalMassP2CL () {}
};

/// \brief Compute the P2 load vector corresponding to the function f on a single tetra.
class LocalVectorDeformP2CL
{
  private:
    InstatScalarFunction f_;
    double time_;

    std::valarray<double> qp2,
                          qf;

  public:
    static const FiniteElementT row_fe_type= P2IF_FE;

    double vec[10];

    LocalVectorDeformP2CL (InstatScalarFunction f, double time) : f_( f), time_( time) {}

    void setup (const TetraCL& t, const InterfaceCommonDataDeformP2CL& cdata, const IdxT numr[10]) {
        resize_and_evaluate_on_vertexes( f_, t, cdata.qdom2d_full, time_, qf);
        qp2.resize( cdata.qdom2d_full.vertex_size());
        for (Uint i= 0; i < 10; ++i) {
                if (numr[i] == NoIdx)
                    continue;
                evaluate_on_vertexes( cdata.p2[i], cdata.qdom2d_only_weights, Addr( qp2));
                vec[i]= quad_2D( qf*qp2, cdata.qdom2d_full);
        }
    }
};

class LocalMassDeformP2CL
{
  private:
    std::valarray<double> qp2[10];

  public:
    static const FiniteElementT row_fe_type= P2IF_FE,
                                col_fe_type= P2IF_FE;

    double coup[10][10];

    void setup (const TetraCL&, const InterfaceCommonDataDeformP2CL& cdata) {
        for (int i= 0; i < 10; ++i)
            resize_and_evaluate_on_vertexes ( cdata.p2[i], cdata.qdom2d_only_weights, qp2[i]);

        for (int i= 0; i < 10; ++i) {
            coup[i][i]= quad_2D (qp2[i]*qp2[i], cdata.qdom2d_only_weights);
            for(int j= 0; j < i; ++j)
                coup[i][j]= coup[j][i]= quad_2D (qp2[j]*qp2[i], cdata.qdom2d_only_weights);
        }
    }

    LocalMassDeformP2CL () {}
};

class LocalLaplaceBeltramiDeformP2CL
{
  private:
    double D_; // diffusion coefficient

    GridFunctionCL<Point3DCL> qgradp2[10];

    GridFunctionCL<Point3DCL> nl;
    GridFunctionCL<SMatrixCL<3,3> > Ginv_wwT;

  public:
    static const FiniteElementT row_fe_type= P2IF_FE,
                                col_fe_type= P2IF_FE;

    double coup[10][10];

    void setup (const TetraCL&, const InterfaceCommonDataDeformP2CL& cdata) {
        Ginv_wwT.resize( cdata.qdom2d_only_weights.vertex_size());
        for (Uint i= 0; i < cdata.qdom2d_only_weights.vertex_size(); ++i) {
            cdata.Phi.set_point (cdata.qdom2d_only_weights.vertex_begin()[i], true);
            Ginv_wwT[i]= cdata.Phi.Ginv_wwT;
        }

        for (int i= 0; i < 10; ++i) {
            resize_and_evaluate_on_vertexes ( cdata.gradp2[i], cdata.qdom2d_only_weights, qgradp2[i]);
        }

        for (int i= 0; i < 10; ++i) {
            coup[i][i]= D_*quad_2D( dot( qgradp2[i], Ginv_wwT, qgradp2[i]), cdata.qdom2d_only_weights);
            for(int j= 0; j < i; ++j)
                coup[i][j]= coup[j][i]= D_*quad_2D( dot( qgradp2[j], Ginv_wwT, qgradp2[i]), cdata.qdom2d_only_weights);
        }
    }

    LocalLaplaceBeltramiDeformP2CL (double D)
        :D_( D) {}
};

class LocalNormalLaplaceDeformP2CL
{
  private:
    double D_; // diffusion coefficient

    GridFunctionCL<Point3DCL> qgradp2[10];
    GridFunctionCL<> qngradp2[10];

  public:
    static const FiniteElementT row_fe_type= P2IF_FE,
                                col_fe_type= P2IF_FE;

    double coup[10][10];

    void setup (const TetraCL&, const InterfaceCommonDataDeformP2CL& cdata) {
        for (int i= 0; i < 10; ++i) {
            resize_and_evaluate_on_vertexes ( cdata.gradp2[i], cdata.qdom, qgradp2[i]);
            qngradp2[i].resize (cdata.qdom.vertex_size ());
            for (Uint j= 0; j < qgradp2[i].size(); ++j) {
                cdata.Phi.set_point (cdata.qdom.vertex_begin()[j], true);
                qngradp2[i][j]= inner_prod( cdata.Phi.w, qgradp2[i][j]);
             }
        }

        for (int i= 0; i < 10; ++i) {
            coup[i][i]= D_*quad( qngradp2[i]*qngradp2[i], 1., cdata.qdom);
            for(int j= 0; j < i; ++j)
                coup[i][j]= coup[j][i]= D_*quad( qngradp2[j]*qngradp2[i], 1., cdata.qdom);
        }
    }

    LocalNormalLaplaceDeformP2CL (double D)
        :D_( D) {}
};



/// \brief The routine sets up the load-vector in v on the interface defined by ls.
///        It belongs to the FE induced by standard P1-elements.
void SetupInterfaceRhsP1 (const MultiGridCL& mg, VecDescCL* v,
    const VecDescCL& ls, const BndDataCL<>& lsbnd, InstatScalarFunction f, double t = 0.);


///\brief Initialize the QuadDomain2DCL-object qdom for quadrature with Quad5_2DCL on the lattice lat of t, given the level set in ls and bnd.
inline const QuadDomain2DCL&
make_CompositeQuad5Domain2D (QuadDomain2DCL& qdom, const TetraCL& t, const PrincipalLatticeCL& lat, const DROPS::VecDescCL& ls, const DROPS::BndDataCL<>& bnd)
{
    LocalP2CL<> locp2_ls( t, ls, bnd);
    std::valarray<double> ls_loc( lat.vertex_size());
    evaluate_on_vertexes( locp2_ls, lat, Addr( ls_loc));
    if (equalSigns(ls_loc)) {
        qdom.clear();
        return qdom;
    }
    SurfacePatchCL surf;
    surf.make_patch<MergeCutPolicyCL>( lat, ls_loc);
    return make_CompositeQuad5Domain2D ( qdom, surf, t);
}

/// \brief Short-hand for integral on the interface.
template <typename T>
double Integral_Gamma (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const DROPS::BndDataCL<>& bnd,
        const InstatFunction<T>& discsol, Uint lattice_num_intervals= 2){
    const DROPS::Uint lvl = ls.GetLevel();
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( lattice_num_intervals);

    std::valarray<double> q;
    QuadDomain2DCL qdom;
    double d( 0.);
    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        make_CompositeQuad5Domain2D( qdom, *it, lat, ls, bnd);
        resize_and_evaluate_on_vertexes( discsol, *it, qdom, 0., q);
        d+= quad_2D( q, qdom);
    }
    return d;
}

template <typename DiscP1FunT>
double Integral_Gamma (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const DROPS::BndDataCL<>& bnd,
    const DiscP1FunT& discsol, Uint lattice_num_intervals= 2)
{
    const DROPS::Uint lvl = ls.GetLevel();
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( lattice_num_intervals);

    std::valarray<double> q;
    QuadDomain2DCL qdom;
    double d( 0.);
    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        make_CompositeQuad5Domain2D( qdom, *it, lat, ls, bnd);
        resize_and_evaluate_on_vertexes( discsol, *it, qdom, q);
        d+= quad_2D( q, qdom);
    }
    return d;
}

/// \brief Short-hand for integral on the interface, h^3-version through interpolation.
template <typename DiscP1FunT>
double Integral_Gamma_Extrapolate (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const DROPS::BndDataCL<>& bnd,
    const DiscP1FunT& discsol)
{
    return (4.*Integral_Gamma( mg, ls, bnd, discsol, 2) - Integral_Gamma_Coarse( mg, ls, bnd, discsol, 1))/3.;
}

/// Compute some data with QuaQuaMapperCL and simple integrals on the higher order interface.
/// Optionally, compute the vector field p_h explicitly for visualization.
/// This is only correct for 1 OpenMP-Thread due to update races.
class InterfaceDebugP2CL : public TetraAccumulatorCL
{
  private:
    const InterfaceCommonDataP2CL& cdata_;

    VecDescCL* to_iface; // For all P2-dofs x at the interface: p_h(x) - x. Computed if to_iface != 0.

    InstatMatrixFunction ref_dp;
    double (*ref_abs_det) (const TetraCL& t, const BaryCoordCL& b, const SurfacePatchCL& surf);

    double max_dph_err,
           surfacemeasP1,
           surfacemeasP2,
           max_absdet_err,
           max_dph2_err;
    double true_area;

  public:
    void store_offsets( VecDescCL& to_ifacearg) { to_iface= &to_ifacearg; }
    void set_true_area( double a) { true_area= a; }
    void set_ref_dp   (InstatMatrixFunction rdp) { ref_dp= rdp; }
    void set_ref_abs_det   ( double (*rad) (const TetraCL& t, const BaryCoordCL& b, const SurfacePatchCL& surf)) { ref_abs_det= rad; }

    InterfaceDebugP2CL (const InterfaceCommonDataP2CL& cdata);
    virtual ~InterfaceDebugP2CL () {}

    virtual void begin_accumulation   ();
    virtual void finalize_accumulation();

    virtual void visit (const TetraCL& t);

    virtual InterfaceDebugP2CL* clone (int /*clone_id*/) { return new InterfaceDebugP2CL( *this); }
};

/// \brief P1-discretization and solution of an equation on the interface
class SurfacePDEP1BaseCL
        {
    public:
        typedef BndDataCL<Point3DCL>  VelBndDataT;
        typedef NoBndDataCL<>         BndDataT;

    protected:
        MultiGridCL&  MG_;
        //SolverBaseCL& solver_;

        double  theta_, ///< time scheme parameter
                dt_;    ///< time step size

        BndDataT            Bnd_;    ///< Dummy boundary data for interface solution

        const VelBndDataT&  Bnd_v_;  ///< Boundary condition for the velocity
        VecDescCL*          v_;      ///< velocity at current time step
        VecDescCL&          lset_vd_;///< levelset at current time step
        const BndDataCL<>&  lsetbnd_;///< level set boundary

        VecDescCL           oldls_;  ///< levelset at old time
        VecDescCL           oldv_;   ///< velocity at old time
        double              oldt_;   ///< old time

    public:
        SurfacePDEP1BaseCL(MultiGridCL& mg,
                         double theta, VecDescCL* v, const VelBndDataT& Bnd_v, VecDescCL& lset_vd, const BndDataCL<>& lsetbnd,
                         int iter= 1000, double tol= 1e-7, double omit_bound= -1.)
                : MG_( mg), theta_( theta), dt_( 0.), Bnd_v_( Bnd_v), v_( v), lset_vd_( lset_vd), lsetbnd_( lsetbnd)
        {}
        virtual ~SurfacePDEP1BaseCL () {}

        const MultiGridCL& GetMG() const { return MG_; }

        //SolverBaseCL& GetSolver() { return solver_; }

        /// set the parameter of the theta-scheme for time stepping
        void SetTheta (double theta);

        /// perform one time step to new_t.
        virtual void DoStep (double /*new_t*/) {}
    };



/// \brief P1-discretization and solution of the transport equation on the interface
    class SurfactantP1BaseCL: public SurfacePDEP1BaseCL
    {
    public:
        typedef P1EvalCL<double, const BndDataT, VecDescCL>       DiscSolCL;
        typedef P1EvalCL<double, const BndDataT, const VecDescCL> const_DiscSolCL;

        IdxDescCL idx; ///< index desctription for concentration at current time
        VecDescCL ic;  ///< concentration on the interface at current time

        VecDescCL iface;  ///< interface mesh at current time
        VecDescCL iface_old;  ///< interface mesh at current time

    protected:
        double        D_;     ///< diffusion coefficient

        InstatScalarFunction rhs_fun_; ///< function for a right-hand side

        IdxDescCL           oldidx_; ///< idx that corresponds to old time (and oldls_)
        VectorCL            oldic_;  ///< interface concentration at old time

        GSPcCL                  pc_;
        GMResSolverCL<GSPcCL>   gm_;

//    //block solver`
//    DiagBlockPcT block_pc_;
//    GMResSolverCL<DiagBlockPcT> GMRes_;
//    BlockMatrixSolverCL<GMResBlockT> block_gm_;

        double omit_bound_; ///< not used atm

    public:
        SurfactantP1BaseCL (MultiGridCL& mg,
                            double theta, double D, VecDescCL* v, const VelBndDataT& Bnd_v, VecDescCL& lset_vd, const BndDataCL<>& lsetbnd,
                            int iter= 1000, double tol= 1e-7, double omit_bound= -1.)
                : SurfacePDEP1BaseCL(mg, theta, v, Bnd_v, lset_vd, lsetbnd),
                  idx( P1IF_FE), D_( D),  rhs_fun_( 0), oldidx_( P1IF_FE), gm_( pc_, 100, iter, tol, true),
                  omit_bound_( omit_bound)
        { idx.GetXidx().SetBound( omit_bound);}

        /*SurfactantP1BaseCL (MultiGridCL& mg,
                            double theta, double D, VecDescCL* v, const VelBndDataT& Bnd_v, VecDescCL& lset_vd, const BndDataCL<>& lsetbnd,double width,
                            int iter= 1000, double tol= 1e-10, double omit_bound= -1.)//This is to implement a method in a trip near the interface
        // width is the distance between the interface and the boundary of the strip.
        // omit_bound is not used any more
                : SurfacePDEP1BaseCL(mg, theta, v, Bnd_v, lset_vd, lsetbnd),
                idx( P1IF_FE), D_( D), rhs_fun_( 0), oldidx_( P1IF_FE), gm_( pc_, 1000, iter, tol, true),
                  omit_bound_( omit_bound)
        { idx.GetXidx().SetBound( omit_bound); }//we use the idx.extIdx_.omit_bound_ to transfer the information of the width of the strip
        //This will be used in IndexDescCL::CreatNumbering()*/
        virtual ~SurfactantP1BaseCL () {}

        GMResSolverCL<GSPcCL>& GetSolver() { return gm_; }

        const_DiscSolCL GetSolution() const
        { return const_DiscSolCL( &ic, &Bnd_, &MG_); }
        const_DiscSolCL GetSolution( const VecDescCL& Myic) const
        { return const_DiscSolCL( &Myic, &Bnd_, &MG_); }

        /// initialize the interface concentration
        void SetInitialValue (InstatScalarFunction, double t= 0.);

        /// set the parameter of the theta-scheme for time stepping
        void SetRhs (InstatScalarFunction);

        /// save a copy of the old level-set and velocity; moves ic to oldic; must be called before DoStep.
        virtual void InitTimeStep ();

        virtual void DoStep0 (double /*new_t*/) {}  //used only in Class: SurfactantExtensionP1CL

        /// perform one time step to new_t.
        virtual void DoStep (double /*new_t*/) {}
    };

    class TransportP2FunctionCL;//???

///The class is based on a stablized term in a narrow band near the interface
    class SurfactantNarrowBandStblP1CL: public SurfactantP1BaseCL
    {
    public:
        IdxDescCL full_idx;
        MatDescCL Laplace,  ///< diffusion matrix,
                Volume_stab, ///< stabilization matrix,
                Mass,  ///< mass matrix
                Conv,  ///< convection matrix
                Massd; ///< mass matrix with interface-divergence of velocity

        VecDescCL rhsext1;
        VectorCL load, ///< for a load-function
                rhs1_, ///< for the extension initial data
                rhs2_; ///< for the extension initial data
        const double width_;
        const double rho_;///<stabilization parameter for Volume_stab


    private:
        MatrixCL      L_; ///< sum of matrices
        //  InstatScalarFunction lvlset_; ///< must be the signed distance function
        InstatVectorFunction normal_; ///< the level-set function
        TransportP2FunctionCL* fulltransport_;

    public:
        SurfactantNarrowBandStblP1CL (MultiGridCL& mg,
                                      double theta, double D, VecDescCL* v, const VelBndDataT& Bnd_v, VecDescCL& lset_vd, const BndDataCL<>& lsetbnd,
                                      InstatVectorFunction normal,const double width, double rho,
                                      int iter= 1000, double tol= 1e-7, double omit_bound= -1.)
                :  SurfactantP1BaseCL( mg, theta, D, v, Bnd_v, lset_vd, lsetbnd, iter, tol, omit_bound),
                //: SurfactantP1BaseCL( mg, theta, D, v, Bnd_v, lset_vd, lsetbnd, width, iter, tol, omit_bound),
        full_idx( P2_FE), rhsext1(), normal_(normal), width_(width), rho_(rho)
                   {}

        /// \remarks call SetupSystem \em before calling SetTimeStep!
        //void SetTimeStep( double dt, double theta=-1);

        /// perform one time step

        void DoStep0 (double new_t); //Backward Euler
        void DoStep (double new_t);  //BDF2 method

        const_DiscSolCL GetSolution() const
        { return const_DiscSolCL( &ic, &Bnd_, &MG_); }
        const_DiscSolCL GetSolution( const VecDescCL& Myic) const
        { return const_DiscSolCL( &Myic, &Bnd_, &MG_); }
        ///@}

        /// \name For internal use only
        /// The following member functions are added to enable an easier implementation
        /// of the coupling navstokes-levelset. They should not be called by a common user.
        /// Use DoStep() instead.
        ///@{
        void InitStep1 (double new_t);// one-step--Implicit Euler
        void InitStep2 (double new_t);// two-step--BDF2 Euler, Use Backward Euler in the first time step
        void InitStep3 (double new_t);// two-step--BDF2 Euler, In the first time step, solve the equation with smaller time step
        void DoStep1 (); //one step-- Implicit Euler
        void DoStep2 (); //two step-- BDF2 method
        void CommitStep ();
        void Update ();
        ///@}
    };

/// \brief P1-discretization and solution of the CahnHilliard equation on the interface
class CahnHilliardP1BaseCL: public SurfacePDEP1BaseCL
    {
    public:
    MatDescCL
            Laplace,  ///< diffusion matrix div_Gamma( grad_Gamma)
            LaplaceM,  ///< diffusion matrix with mobility div_Gamma(M grad_Gamma)
            Volume_stab, ///< stabilization matrix, tetra integral over normal gradients
            Ident,
            Mass,  ///< mass matrix
            Conv,  ///< convection matrix
            Massd, ///< mass matrix with interface-divergence of velocity
            Mass2; ///< mass matrix: new trial- and test- functions on old interface

        typedef P1EvalCL<double, const BndDataT, VecDescCL>       DiscSolCL;
        typedef P1EvalCL<double, const BndDataT, const VecDescCL> const_DiscSolCL;

        typedef PCGSolverCL<SSORPcCL> PCGSolverT;
        typedef BlockPreCL<SchurPreBaseCL, SchurPreBaseCL, DiagBlockPreCL>  DiagBlockPcT; //block preconditioner
        typedef GMResSolverCL<DiagBlockPcT> GMResBlockT; //block GMRES solver
        typedef BlockMatrixSolverCL<GMResBlockT> GMResBlockSolver;

        IdxDescCL idx; ///< index desctription for concentration at current time
        VecDescCL ic;  ///< concentration on the interface at current time
        VecDescCL ienergy;  ///< energy on the interface at current time

    VecDescCL conc_extrapol;
    VecDescCL species_extrapol;

    /*IdxDescCL idx_mu; ///< index desctription for chemical potential at current time*/
        VecDescCL imu;  ///< chemical potential on the interface at current time
    VecDescCL is;  ///< third species concentration on the interface at current time

    VecDescCL iface;  ///< interface mesh at current time


protected:
        double        sigma_;     ///< mobility coefficient
        double        epsilon_;  ///< epsilon coefficient

        InstatScalarFunction rhs_fun3_; ///< function for a right-hand side to concentration equation
        InstatScalarFunction rhs_fun4_; ///< function for a right-hand side chemical potential equation
        InstatScalarFunction rhs_fun5_; ///< function for a right-hand side species equation


        IdxDescCL           oldidx_; ///< idx that corresponds to old time (and oldoldls_)
        IdxDescCL           oldoldidx_; ///< idx that corresponds to old old time (and oldoldls_)

        VectorCL            oldic_;  ///< interface concentration at old time
         VectorCL           oldis_;  ///< interface species at old time

        VectorCL            oldoldic_;  ///< interface concentration at old old time

        VectorCL            oldimu_;  ///< interface chemical potential at old time



        GSPcCL                  pc_;
        GMResSolverCL<GSPcCL>   gm_;
        MatrixCL dummy_matrix3_, dummy_matrix4_;

        SSORPcCL symmPcPc_;
        PCGSolverT PCGSolver3_, PCGSolver4_;
        SurfaceLaplacePreCL<PCGSolverT> spc3_, spc4_;
        DiagBlockPcT block_pc_;
        GMResBlockT GMRes_;
        double omit_bound_; ///< not used atm

    public:

    GMResBlockSolver block_gm_;

    CahnHilliardP1BaseCL (MultiGridCL& mg,
                              double theta, double sigma, double epsilon, VecDescCL* v, const VelBndDataT& Bnd_v, VecDescCL& lset_vd, const BndDataCL<>& lsetbnd,
                              //MatrixCL Precond3 , MatrixCL Precond4,
                              int iter= 1000, double tol= 1e-7, double iterA=500, double tolA=1e-3, double iterB=500, double tolB=1e-3, double omit_bound= -1.)
                : SurfacePDEP1BaseCL(mg, theta, v, Bnd_v, lset_vd, lsetbnd), idx( P1IF_FE), /*idx_mu( P1IF_FE),*/ omit_bound_( omit_bound),
                sigma_( sigma), epsilon_( epsilon),  rhs_fun3_( 0), rhs_fun4_( 0), oldidx_( P1IF_FE),oldoldidx_( P1IF_FE),
                  conc_extrapol(),//yushutin
                PCGSolver3_(symmPcPc_, iterA, tolA, true),
                PCGSolver4_(symmPcPc_, iterB, tolB, true),
                spc3_( dummy_matrix3_, PCGSolver3_), spc4_( dummy_matrix4_, PCGSolver4_),
                block_pc_(spc3_,spc4_),
                GMRes_(block_pc_, 10, iter, tol, false, false, RightPreconditioning, true, false),
                gm_(pc_, 50, iter, tol, true),
                block_gm_(GMRes_)
        {
            idx.GetXidx().SetBound( omit_bound);
            //idx_mu.GetXidx().SetBound( omit_bound);
        }

        virtual ~CahnHilliardP1BaseCL () {}

        BlockMatrixSolverCL<GMResBlockT>& GetSolver() { return block_gm_; }

        const_DiscSolCL GetConcentr() const
        { return const_DiscSolCL( &ic, &Bnd_, &MG_); }
        const_DiscSolCL GetConcentr( const VecDescCL& Myic) const
        { return const_DiscSolCL( &Myic, &Bnd_, &MG_); }

        const_DiscSolCL GetSpecies() const
        { return const_DiscSolCL( &is, &Bnd_, &MG_); }
        const_DiscSolCL GetSpecies( const VecDescCL& Myis) const
        { return const_DiscSolCL( &Myis, &Bnd_, &MG_); }

        const_DiscSolCL GetPotential() const
        { return const_DiscSolCL( &imu, &Bnd_, &MG_); }
        const_DiscSolCL GetPotential( const VecDescCL& Myimu) const
        { return const_DiscSolCL( &Myimu, &Bnd_, &MG_); }

        /// initialize the interface concentration
        void SetInitialValue (InstatScalarFunction,InstatScalarFunction,InstatScalarFunction, double t= 0.);

        /// set the parameter of the theta-scheme for time stepping
        void SetRhs (InstatScalarFunction,InstatScalarFunction);

        /// save a copy of the old level-set and velocity; moves ic to oldic; must be called before DoStep.
        virtual void InitTimeStep ();

        /// perform one time step to new_t.
        virtual void DoStep (double /*new_t*/) {}
        virtual void DoStep0 (double /*new_t*/) {}

};

class CahnHilliardcGP1CL : public CahnHilliardP1BaseCL {
    public:
        MatDescCL
            LaplaceNon,  ///< nonlinear diffusion matrix with div_Gamma(f grad_Gamma)
            Massrho;  ///< mass matrix

    const double& width_;///< we extend only a band near the zero leve set with a width
    const double rho_;///<stabilization parameter for Volume_stab
    const double S_;//stabilization parameter for time derivative;
    private:
        MatrixCL      A_, B_, C_, D_; ///< blocks of the matrix
        InstatVectorFunction normal_; ///< the normal vector function


    public:
        CahnHilliardcGP1CL (MultiGridCL& mg, double theta, double sigma, double epsilon,
                            VecDescCL* v, const VelBndDataT& Bnd_v, VecDescCL& lset_vd, const BndDataCL<>& lsetbnd,
                            InstatVectorFunction normal,const double & width, double rho, double S,
                            int iter= 999, double tol= 1.1e-7, double iterA=499, double tolA=1.1e-3, double iterB=499, double tolB=1.1e-3,double omit_bound= -1.)
                : CahnHilliardP1BaseCL( mg, theta, sigma, epsilon, v, Bnd_v, lset_vd, lsetbnd,
                        iter, tol, iterA, tolA, iterB, tolB, omit_bound),
                  normal_(normal),width_(width), rho_(rho), S_(S)
        {}

        /// save a copy of the old level-set and velocity; moves ic to oldic; must be called before DoStep.
        // void InitTimeStep (); // as in the base

        /// perform one time step
        void DoStep (double new_t);
        void DoStep0 (double new_t); //Backward Euler

        /// \name For internal use only
        /// The following member functions are added to enable an easier implementation
        /// of the coupling navstokes-levelset. They should not be called by a common user.
        /// Use DoStep() instead.
        ///@{
        VectorCL InitStep3 (double new_t);
        VectorCL InitStep4 (double new_t);
        void InitStep(VectorCL&, VectorCL&, double new_t);

        void DoStep (const VectorCL&, const VectorCL&);
        void CommitStep ();
        void Update ();
        ///@}
    };

///The class is based on a stablized term in a narrow band near the interface
class    CahnHilliardNarrowBandStblP1CL: public CahnHilliardP1BaseCL {
    public:
        IdxDescCL full_idx;
        MatDescCL
            LaplaceNon,  ///< nonlinear diffusion matrix with div_Gamma(f grad_Gamma)
            Massrho;  ///< mass matrix

        const double S_;//stabilization parameter for time derivative;

        VecDescCL ext1, ext2, ext3;
        VectorCL load, ///< for a load-function
                rhs1_, ///< for the extension initial data
                rhs2_; ///< for the extension initial data

        const double width_;
        const double rho_;///<stabilization parameter for Volume_stab


    private:
        MatrixCL      A_, B_, C_, D_, K_; ///< blocks of the matrix
        ///< //  InstatScalarFunction lvlset_; ///< must be the signed distance function
        InstatVectorFunction normal_; ///< the level-set function
        //TransportP2FunctionCL* fulltransport_;

    public:
        CahnHilliardNarrowBandStblP1CL (MultiGridCL& mg, double theta, double sigma, double epsilon,
                VecDescCL* v, const VelBndDataT& Bnd_v, VecDescCL& lset_vd, const BndDataCL<>& lsetbnd,
                InstatVectorFunction normal,const double  width, double rho, double S,
        int iter= 999, double tol= 1.1e-7, double iterA=499, double tolA=1.1e-3, double iterB=499, double tolB=1.1e-3,double omit_bound= -1.)
        : CahnHilliardP1BaseCL( mg, theta, sigma, epsilon, v, Bnd_v, lset_vd, lsetbnd,
                iter, tol, iterA, tolA, iterB, tolB, omit_bound),
          full_idx( P2_FE), ext1(),ext2(), ext3(),
          //conc_extrapol(),
          normal_(normal),width_(width), rho_(rho), S_(S)
        {

            //ext1.SetIdx(&full_idx);

            //ext2.SetIdx(&full_idx);

        }

        /// \remarks call SetupSystem \em before calling SetTimeStep!
        //void SetTimeStep( double dt, double theta=-1);

        /// perform one time step

        void DoStep0 (double new_t); //Backward Euler
        void DoStep (double new_t);  //BDF2 method

       /* const_DiscSolCL GetSolution() const
        { return const_DiscSolCL( &ic, &Bnd_, &MG_); }
        const_DiscSolCL GetSolution( const VecDescCL& Myic) const
        { return const_DiscSolCL( &Myic, &Bnd_, &MG_); }
        ///@}*/

        /// \name For internal use only
        /// The following member functions are added to enable an easier implementation
        /// of the coupling navstokes-levelset. They should not be called by a common user.
        /// Use DoStep() instead.
        ///@{
        void InitStep1 (VectorCL&, VectorCL&, VectorCL&, double new_t);// one-step--Implicit Euler
        void InitStep2 (VectorCL&, VectorCL&, double new_t);// two-step--BDF2 Euler, Use Backward Euler in the first time step
        void InitStep3 (double new_t);// two-step--BDF2 Euler, In the first time step, solve the equation with smaller time step
        void DoStep1 (VectorCL&, VectorCL&,VectorCL&); //one step-- Implicit Euler
        void DoStep2 (VectorCL& ,VectorCL& ); //two step-- BDF2 method
        void CommitStep ();
        void Update ();
        ///@}
    };

/// \brief P1-discretization and solution of the transport equation on the interface
class SurfactantcGP1CL : public SurfactantP1BaseCL
{
  public:
    MatDescCL A,  ///< diffusion matrix
              M,  ///< mass matrix
              C,  ///< convection matrix
              Md, ///< mass matrix with interface-divergence of velocity
              M2; ///< mass matrix: new trial- and test- functions on old interface

  private:
    MatrixCL      L_; ///< sum of matrices

  public:
    SurfactantcGP1CL (MultiGridCL& mg,
        double theta, double D, VecDescCL* v, const VelBndDataT& Bnd_v, VecDescCL& lset_vd, const BndDataCL<>& lsetbnd,
        int iter= 1000, double tol= 1e-7, double omit_bound= -1.)
    : SurfactantP1BaseCL( mg, theta, D, v, Bnd_v, lset_vd, lsetbnd, iter, tol, omit_bound)
    {}

    /// save a copy of the old level-set and velocity; moves ic to oldic; must be called before DoStep.
    // void InitTimeStep (); // as in the base

    /// perform one time step
    virtual void DoStep (double new_t);

    /// \name For internal use only
    /// The following member functions are added to enable an easier implementation
    /// of the coupling navstokes-levelset. They should not be called by a common user.
    /// Use DoStep() instead.
    ///@{
    VectorCL InitStep (double new_t);
    void DoStep (const VectorCL&);
    void CommitStep ();
    void Update ();
    ///@}
};

/// \brief Observes the MultiGridCL-changes by AdapTriangCL to repair an interface p1-function.
///
/// The function is extended to a P1-function on \f$ \Omega \f$ and then repaired as such. In
/// post_refine_sequence() the new interface-index is generated and the function is restricted
/// to it.
class InterfaceP1RepairCL : public MGObserverCL
{
  private:
    MultiGridCL& mg_;
    const VecDescCL&   lset_vd_;
    const BndDataCL<>& lset_bnd_;
    VecDescCL&         u_;

    IdxDescCL          fullp1idx_;
    VecDescCL          fullu_;
    std::unique_ptr<RepairP1CL<double, NoBndDataCL>::type > p1repair_;
 double width_;
  public:
    InterfaceP1RepairCL (MultiGridCL& mg, const VecDescCL& lset_vd, const BndDataCL<>& lset_bnd,  VecDescCL& u, double width=0)
        : mg_( mg), lset_vd_( lset_vd), lset_bnd_( lset_bnd), u_( u), fullp1idx_( P1_FE), fullu_( &fullp1idx_), width_(width) {}

    void pre_refine  ();
    void post_refine ();

    void pre_refine_sequence  ();
    void post_refine_sequence ();
    const IdxDescCL* GetIdxDesc() const { return u_.RowIdx; }
#ifdef _PAR
    const VectorCL*  GetVector()  const { return &fullu_.Data; }
    void swap( IdxDescCL& idx, VectorCL& v) { fullu_.RowIdx->swap(idx); fullu_.Data.swap(v); }
#endif
};

/// \brief Observes the MultiGridCL-changes by AdapTriangCL to repair an interface p1-function.
///
/// The function is extended to a P1-function on \f$ \Omega \f$ and then repaired as such. In
/// post_refine_sequence() the new interface-index is generated and the function is restricted
/// to it.
class NarrowBandP1RepairCL : public MGObserverCL
{
  private:
    MultiGridCL& mg_;
    const VecDescCL&   lset_vd_;
    const BndDataCL<>& lset_bnd_;
    VecDescCL&         u_;

    IdxDescCL          fullp1idx_;
    VecDescCL          fullu_;
    std::unique_ptr<RepairP1CL<double, NoBndDataCL>::type > p1repair_;

  public:
    NarrowBandP1RepairCL (MultiGridCL& mg, const VecDescCL& lset_vd, const BndDataCL<>& lset_bnd, VecDescCL& u)
        : mg_( mg), lset_vd_( lset_vd), lset_bnd_( lset_bnd), u_( u), fullp1idx_( P1_FE), fullu_( &fullp1idx_) {}

    void pre_refine  ();
    void post_refine ();

    void pre_refine_sequence  ();
    void post_refine_sequence ();
    const IdxDescCL* GetIdxDesc() const { return u_.RowIdx; }
#ifdef _PAR
    const VectorCL*  GetVector()  const { return &fullu_.Data; }
    void swap( IdxDescCL& idx, VectorCL& v) { fullu_.RowIdx->swap(idx); fullu_.Data.swap(v); }
#endif
};

#ifndef _PAR
///\brief Represents a scalar P1 function on the interface as Ensight6 variable by extension to the
///       whole domain.
class Ensight6IfaceScalarCL : public Ensight6VariableCL
{
  private:
    const VecDescCL&   u_;
    MultiGridCL&       mg_;

  public:
    Ensight6IfaceScalarCL (MultiGridCL& mg, const VecDescCL& u, std::string varName, std::string fileName, bool timedep= false)
        : Ensight6VariableCL( varName, fileName, timedep), u_( u), mg_( mg) {}

    void Describe (Ensight6OutCL& cf) const { cf.DescribeVariable( this->varName(), true); }
    void put      (Ensight6OutCL& cf) const;
};

///\brief Create an Ensight6IfaceP1ScalarCL with operator new.
///
/// This is just for uniform code; the analogous functions for scalars and vectors are more useful
/// because they help to avoid template parameters in user code.
inline Ensight6IfaceScalarCL&
make_Ensight6IfaceScalar (MultiGridCL& mg, const VecDescCL& u,
    std::string varName, std::string fileName, bool timedep= false)
{
    return *new Ensight6IfaceScalarCL( mg, u, varName, fileName, timedep);
}
#endif

///\brief Represents a scalar P1/P2 function on the interface as VTK variable by extension to the
///       whole domain.
class VTKIfaceScalarCL : public VTKVariableCL
{
  private:
    const VecDescCL&   u_;
    MultiGridCL&       mg_;
    const BndDataCL<double> BndData_;

  public:
    VTKIfaceScalarCL (MultiGridCL& mg, const VecDescCL& u, std::string varName, const BndDataCL<double>& BndData = BndDataCL<double>(0))
        : VTKVariableCL( varName), u_( u), mg_( mg), BndData_(BndData) {}

    void put      (VTKOutCL& cf) const;
    Uint GetDim() const { return 1; }
};

///\brief Create an VTKIfaceP1ScalarCL with operator new.
///
/// This is just for uniform code; the analogous functions for scalars and vectors are more useful
/// because they help to avoid template parameters in user code.
inline VTKIfaceScalarCL&
make_VTKIfaceScalar (MultiGridCL& mg, const VecDescCL& u,
    std::string varName, const BndDataCL<double>& BndData = BndDataCL<double>(0))
{
    return *new VTKIfaceScalarCL( mg, u, varName, BndData);
}

class VTKIfaceVectorCL : public VTKVariableCL
{
  private:
    const VecDescCL& u_;
    MultiGridCL& mg_;
    int P2_;
    const BndDataCL<Point3DCL>& BndData_;

  public:
    VTKIfaceVectorCL (MultiGridCL& mg, const VecDescCL& u, std::string varName, int P2, const BndDataCL<Point3DCL>& BndData)
        : VTKVariableCL( varName), u_( u), mg_( mg), P2_(P2), BndData_(BndData) {}

      void put      (VTKOutCL& cf) const;
      Uint GetDim() const { return 1; }
};

inline VTKIfaceVectorCL&
make_VTKIfaceVector (MultiGridCL& mg, const VecDescCL& u,
                     std::string varName, std::string FEdegree, const BndDataCL<Point3DCL>& BndData)
{
  //std::cout<<"Element :"<< FEdegree <<std::endl;
  if(!FEdegree.compare("P2"))
    return *new VTKIfaceVectorCL( mg, u, varName, 1, BndData);
  else
    return *new VTKIfaceVectorCL( mg, u, varName, 0, BndData);
}/// ==Space-Time-methods==
/// \todo Maybe introduce a new header

///\brief Bilinear space-time FE-function on a single TetraPrismCL.
template <typename T= double>
class LocalSTP1P1CL
{
  public:
    typedef T value_type;
    typedef T (*instat_fun_ptr) (const Point3DCL&, double);
    typedef LocalP1CL<T> spatial_localfe_type;

    enum { Dim= 8 }; ///< local number of unknowns

  protected:
    typedef LocalSTP1P1CL<T> self_;

  private:
    LocalP1CL<T> v0_,
                 v1_;

  public:
    LocalSTP1P1CL() : v0_(), v1_() {}
    LocalSTP1P1CL(const value_type& t) : v0_( t), v1_( t) {}
    // Initialize from a given function
    LocalSTP1P1CL(const TetraPrismCL& prism, instat_fun_ptr f)
        : v0_( prism.t, f, prism.t0), v1_( prism.t, f, prism.t1) {}

    // These "assignment-operators" correspond to the constructors
    // with multiple arguments
    inline self_&
    assign(const TetraPrismCL& prism, instat_fun_ptr f) {
        v0_.assign( prism.t, f, prism.t0);
        v1_.assign( prism.t, f, prism.t1);
        return *this;
    }

    // pointwise evaluation in STCoordCL-coordinates
    inline value_type operator() (const STCoordCL& p) const {
        return (1. - p.t_ref)*v0_( p.x_bary) + p.t_ref*v1_( p.x_bary);
    }

    // The LocalP1-function at t0
    spatial_localfe_type& at_t0 () { return v0_; }
    const spatial_localfe_type& at_t0 () const { return v0_; }
    // The LocalP1-function at t1
    spatial_localfe_type& at_t1 () { return v1_; }
    const spatial_localfe_type& at_t1 () const { return v1_; }
};

///\brief Finear space-time FE-function on a single TetraPrismCL.
/// The values on the penta are defined as those on the tetra at t0 and the value in the vertex (v0, t1).
template <typename T= double>
class LocalSTP1CL : public GridFunctionCL<T>
{
  public:
    typedef GridFunctionCL<T> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::instat_fun_ptr instat_fun_ptr;

    enum { Dim= 5 }; ///< local number of unknowns

  protected:
    typedef LocalSTP1CL<T> self_;

  public:
    LocalSTP1CL() : base_type( value_type(), Dim) {}
    LocalSTP1CL(const value_type& t): base_type( t, Dim) {}

DROPS_DEFINE_VALARRAY_DERIVATIVE(LocalSTP1CL, T, base_type)

    // pointwise evaluation in STCoords.
    value_type operator() (const STCoordCL& p) const {
        return LinearCombinationCL<self_, value_type>::do_it( *this, p.x_bary[0] - p.t_ref,
                                                                     p.x_bary[1],
                                                                     p.x_bary[2],
                                                                     p.x_bary[3],
                                                                     p.t_ref);
    }
};

///\brief Quadratic in space and linear in time space-time FE-function on a single TetraPrismCL.
/// The class holds references to two LocalP2CL-objects.
template <typename T= double>
class LocalSTP2P1CL
{
  public:
    typedef T value_type;
    typedef T (*instat_fun_ptr) (const Point3DCL&, double);
    typedef LocalP2CL<T> spatial_localfe_type;

  protected:
    typedef LocalSTP2P1CL<T> self_;

  private:
    LocalP2CL<T> v0_,
                 v1_;

  public:
    LocalSTP2P1CL() {}
    LocalSTP2P1CL(const spatial_localfe_type& v0, const spatial_localfe_type& v1)
        : v0_( v0), v1_( v1) {}

    template<class BndDataT>
      self_&
      assign (const TetraCL& t, const VecDescCL& v0, const VecDescCL& v1, const BndDataT& bnd) {
        v0_.assign( t, v0, bnd);
        v1_.assign( t, v1, bnd);
        return *this;
    }
    template <class STP2P1FunT>
      self_&
      assign(const TetraCL& t, const STP2P1FunT& f) {
        v0_.assign( t, *f.GetSolution( 0), *f.GetBndData());
        v1_.assign( t, *f.GetSolution( 1), *f.GetBndData());
        return *this;
    }

    // pointwise evaluation in STCoordCL-coordinates
    inline value_type operator() (const STCoordCL& p) const {
        return (1. - p.t_ref)*v0_( p.x_bary) + p.t_ref*v1_( p.x_bary);
    }

    // The LocalP1-function at t0
    const spatial_localfe_type& at_t0 () const { return v0_; }
          spatial_localfe_type& at_t0 ()       { return v0_; }
    // The LocalP1-function at t1
    const spatial_localfe_type& at_t1 () const { return v1_; }
          spatial_localfe_type& at_t1 ()       { return v1_; }
};

class STWindProxyCL
{
  private:
    const LocalSTP2P1CL<Point3DCL>& w_;

  public:
    STWindProxyCL (const LocalSTP2P1CL<Point3DCL>& w) : w_( w) {}

    // pointwise evaluation in STCoordCL-coordinates
    inline Point4DCL operator() (const STCoordCL& p) const {
        const Point3DCL tmp( w_( p));
        return MakePoint4D( tmp[0], tmp[1], tmp[2], 1.);
    }
};

class STP1P1DiscCL
{
  public:
    static LocalSTP1P1CL<>        ref_val[8];
    static LocalSTP1CL<Point4DCL> ref_grad[8];
    static LocalSTP1CL<Point3DCL> ref_gradx[8];

    static void StaticInit();
    static void StaticDestruct() {}

    static void GetGradients (const TetraPrismCL& prism, LocalSTP1CL<Point4DCL> grad[8]) {
        SMatrixCL<3,3> M( Uninitialized);
        double det= 0.; // dummy
        GetTrafoTr( M, det, prism.t);
        const double dt_inv= 1./(prism.t1 - prism.t0);
        for (Uint i= 0; i < 4; ++i) {
            const Point3DCL& p1grad= M*FE_P1CL::DHRef( i);
            const Point4DCL& p1grad4d= MakePoint4D( p1grad[0], p1grad[1], p1grad[2], 0.);
            // grad[i]= Point4DCL();
            for (Uint j= 0; j < 4; ++j) {
                grad[i][j]= p1grad4d;
                grad[i][j][3]= j == i ? -dt_inv : 0.;
            }
            grad[i][4][3]= i == 0 ? -dt_inv : 0.;
            // grad[i + 4]= Point4DCL();
            grad[i + 4][i][3]=  dt_inv;
            grad[i + 4][4]= p1grad4d;
            grad[i + 4][4][3]= i == 0 ? dt_inv : 0.;

        }
    }
    static void GetSpatialGradients (const TetraPrismCL& prism, LocalSTP1CL<Point3DCL> grad[8]) {
        SMatrixCL<3,3> M( Uninitialized);
        double det= 0.; // dummy
        GetTrafoTr( M, det, prism.t);
        for (Uint i= 0; i < 4; ++i) {
            const Point3DCL& p1grad= M*FE_P1CL::DHRef( i);
            for (Uint j= 0; j < 4; ++j)
                grad[i][j]= p1grad;
            grad[i + 4][4]= p1grad;
        }
    }
};


///\brief Evaluates a function expecting world-coordinates and time as arguments in STCoord-coordinates on a tetra-prism.
template <class T= double>
class STCoordEvalCL
{
  public:
    using fun_type = std::function<T(Point3DCL const &, double)>;
    typedef T value_type;
  private:
    STCoord2WorldCoordCL mapper_;
    fun_type f_;

  public:
    STCoordEvalCL (const TetraPrismCL& prism, fun_type f)
        : mapper_( prism), f_(f) {}

    value_type operator() (const STCoordCL& b) const { return f_( mapper_.space( b), mapper_.time( b)); }
};

template <class T, class DomainT, class ResultIterT>
  inline ResultIterT
  evaluate_on_vertexes (std::function<T(const Point3DCL&, double)> const & f, const TetraPrismCL& prism, const DomainT& dom, ResultIterT result_iterator)
{
    return std::transform( dom.vertex_begin(), dom.vertex_end(), result_iterator, STCoordEvalCL<T>( prism, f));
}

template <class T, class DomainT, class ResultContT>
  inline const ResultContT&
  resize_and_evaluate_on_vertexes (std::function<T(const Point3DCL&, double)> const & f, const TetraPrismCL& prism, const DomainT& dom, ResultContT& result_container)
{
    result_container.resize( dom.vertex_size());
    evaluate_on_vertexes( f, prism, dom, sequence_begin( result_container));
    return result_container;
}


//idx_ini: oldls, oldt, IFP1_FE (0.. N_ini-1)
//idx0, idx1: ST-P1P1: 2x IFP1_FE: idxi for dof at ti; idx0 agrees with idx_ini, where the latter is defined (well-posed, see lemma 1)
class STP1P1IdxDescCL
{
  private:
    IdxDescCL idx0_,
              idx1_;

    Uint TriangLevel_;
    IdxT NumUnknowns_;
    IdxT NumIniUnknowns_,
         NumFiniUnknowns_;

  public:
    STP1P1IdxDescCL () : idx0_( P1IF_FE), idx1_( P1IF_FE) {}

    /// \brief Returns the number of the index. This can be used to access
    ///     the numbering on the simplices.
    Uint GetIdx (Uint i) const { return i == 0 ? idx0_.GetIdx() : idx1_.GetIdx(); }

    /// \brief Triangulation of the index.
    Uint TriangLevel () const { return TriangLevel_; }
    /// \brief total number of unknowns on the triangulation
    IdxT NumUnknowns () const { return NumUnknowns_; }
    /// \brief The first NumIniUnknowns shape-functions do not vanish on the old interface.
    IdxT NumIniUnknowns () const { return NumIniUnknowns_; }
    /// \brief The unknowns in [NumIniUnknowns, NumFiniUnknowns + NumIniUnknowns)  do not vanish on the new interface.
    IdxT NumFiniUnknowns () const { return NumFiniUnknowns_; }

    void CreateNumbering (Uint level, MultiGridCL& mg, const VecDescCL& oldls, const VecDescCL& newls, const BndDataCL<>& lsetbnd, double t0, double t1);

    void DeleteNumbering (MultiGridCL& mg) {
        DeleteNumbOnSimplex( idx0_.GetIdx(), mg.GetAllVertexBegin( TriangLevel()), mg.GetAllVertexEnd( TriangLevel()));
        DeleteNumbOnSimplex( idx1_.GetIdx(), mg.GetAllVertexBegin( TriangLevel()), mg.GetAllVertexEnd( TriangLevel()));
        TriangLevel_= static_cast<Uint>( -1);
        NumUnknowns_= NumIniUnknowns_= NumFiniUnknowns_= 0;
    }
};


/// \brief Collect indices of unknowns -- bilinear space-time FE with standard Lagrange-basis in space and time
///
/// This is convenient for discretisation of operators in the Setup-routines.
class LocalNumbSTP1P1CL
{
  private:
    size_t NumIniUnknowns_,
           NumFiniUnknowns_;

  public:
    enum { Dim= 8 }; ///< local number of unknowns (max.)

    /// \brief Field of unknown-indices; NoIdx, iff the degree of freedom lies
    /// on a boundary without unknowns.
    IdxT num[Dim];


    /// \brief The default constructors leaves everything uninitialized.
    LocalNumbSTP1P1CL() {}

    /// \brief Read indices only from a tetrahedron.
    LocalNumbSTP1P1CL(const TetraCL& s, const STP1P1IdxDescCL& idx)
        { assign_indices_only( s, idx); }

    /// \brief Compute the indices only.
    /// Only num is set up.
    void assign_indices_only (const TetraCL& s, const STP1P1IdxDescCL& idx) {
        const Uint sys0= idx.GetIdx( 0),
                   sys1= idx.GetIdx( 1);
        for (Uint i= 0; i < 4; ++i) {
            num[i]=   s.GetVertex( i)->Unknowns.Exist( sys0) ? s.GetVertex( i)->Unknowns( sys0) : NoIdx;
            num[i+4]= s.GetVertex( i)->Unknowns.Exist( sys1) ? s.GetVertex( i)->Unknowns( sys1) : NoIdx;
        }
        NumIniUnknowns_=  idx.NumIniUnknowns();
        NumFiniUnknowns_= idx.NumFiniUnknowns();
    }

    /// \brief True, iff index i has a dof associated with it.
    bool WithUnknowns (IdxT i) const { return num[i] != NoIdx; }

    /// \brief True, iff the shape function exists on the old interface.
    bool IsIni  (IdxT i) const { return num[i] < NumIniUnknowns_; }
    /// \brief True, iff the shape function exists on the new interface.
    bool IsFini (IdxT i) const { return num[i] >= NumIniUnknowns_ && num[i] < NumIniUnknowns_ + NumFiniUnknowns_; }

    /// \brief The first NumIniUnknowns shape-functions do not vanish on the old interface.
    IdxT NumIniUnknowns () const { return NumIniUnknowns_; }
};

template<class Data, class BndData_, class VD_>
class STP2P1EvalCL
{
public:
    typedef Data     DataT;
    typedef BndData_ BndDataCL;
    typedef VD_      VecDescT;

    typedef LocalSTP2P1CL<DataT> LocalFET;

protected:
    // numerical data
    VecDescT*          v0_;
    VecDescT*          v1_;
    // boundary-data
    BndDataCL*         bnd_;
    // the multigrid
    const MultiGridCL* MG_;

public:
    STP2P1EvalCL() :v0_( 0), v1_( 0), bnd_( 0), MG_( 0) {}
    STP2P1EvalCL(VecDescT* v0, VecDescT* v1, BndDataCL* bnd, const MultiGridCL* MG)
        : v0_( v0), v1_( v1), bnd_( bnd), MG_( MG) {}
    //default copy-ctor, dtor, assignment-op
    // copying is safe - it is a flat copy, which is fine,
    // as STP2P1EvalCL does not take possession of the pointed to.

    void // set / get the container of numerical data
    SetSolution(Uint i, VecDescT* v) { (i== 0 ? v0_ : v1_)= v; }
    const VecDescT*
    GetSolution(Uint i) const { return i == 0 ? v0_: v1_; }
    VecDescT*
    GetSolution(Uint i) { return i == 0 ? v0_: v1_; }
    void // set / get the container of boundary-data
    SetBndData(BndDataCL* bnd) { bnd_= bnd; }
    BndDataCL*
    GetBndData() { return bnd_; }
    const BndDataCL*
    GetBndData() const { return bnd_; }
    const MultiGridCL& // the multigrid we refer to
    GetMG() const { return *MG_; }
    Uint // Triangulation level of this function
    GetLevel() const { return v0_->GetLevel(); }
    // The time at which data is evaluated.
    double GetTime(Uint i) const { return i == 0 ? v0_->t : v1_->t; }

};

// Create a STP2P1EvalCL without the agonizing template-pain.
template<class BndData_, class VD_>
  STP2P1EvalCL<typename BndData_::bnd_type, BndData_, VD_>
    make_STP2P1Eval (const MultiGridCL& mg, BndData_& bnd, VD_& vd0, VD_& vd1)
{
    return STP2P1EvalCL<typename BndData_::bnd_type, BndData_, VD_>( &vd0, &vd1, &bnd, &mg);
}

enum ZeroPolicyEnum { KeepLocalZeros, RemoveExactLocalZeros };

template <ZeroPolicyEnum ZeroPolicy, typename LocalRowNumbT, typename LocalColNumbT>
  inline void
  update_global_matrix (MatrixBuilderCL& M, const double coup[LocalRowNumbT::Dim][LocalColNumbT::Dim], const LocalRowNumbT& r, const LocalColNumbT& c)
{
    for (Uint i= 0; i < LocalRowNumbT::Dim; ++i)
        if (r.num[i] != NoIdx)
            for (Uint j= 0; j < LocalColNumbT::Dim; ++j)
                if (c.num[j] != NoIdx)
                    if (ZeroPolicy == KeepLocalZeros || (ZeroPolicy == RemoveExactLocalZeros && coup[i][j] != 0.))
                        M( r.num[i], c.num[j])+= coup[i][j];
}

template <ZeroPolicyEnum ZeroPolicy, typename LocalRowNumbT, typename LocalColNumbT>
  inline void
  update_global_matrix_and_coupling (MatrixBuilderCL& M, const double coup[LocalRowNumbT::Dim][LocalColNumbT::Dim], const LocalRowNumbT& r, const LocalColNumbT& c, VectorCL* ini_cpl, const VectorCL* ini_data)
{
    for (Uint i= 0; i < LocalRowNumbT::Dim; ++i)
        if (r.WithUnknowns( i) && !r.IsIni( i))
            for (Uint j= 0; j < LocalColNumbT::Dim; ++j)
                if (c.WithUnknowns( j)) {
                    if (c.IsIni( j)) // Put initial data into coupling vector on the rhs.
                        ini_cpl[0][r.num[i] - r.NumIniUnknowns()]-= coup[i][j]*ini_data[0][c.num[j]];
                    else  if (ZeroPolicy == KeepLocalZeros || (ZeroPolicy == RemoveExactLocalZeros && coup[i][j] != 0.))
                        M( r.num[i] - r.NumIniUnknowns(), c.num[j] - c.NumIniUnknowns())+= coup[i][j];
                }
}

class STInterfaceCommonDataCL : public TetraAccumulatorCL
{
  private:
    STInterfaceCommonDataCL** the_clones;

    LocalSTP2P1CL<> st_local_ls;

    const VecDescCL*   old_ls;  // P2-level-set at t0
    const VecDescCL*   new_ls;  // P2-level-set at t1
    const BndDataCL<>* lsetbnd; // boundary data for the level set function

  public:
    const TetraPrismLatticeCL& lat;

    std::valarray<double> ls_loc;
    SPatchCL<4> surf;
    QuadDomainCodim1CL<4> q5dom;

    const double t0,
                 t1;

    const STInterfaceCommonDataCL& get_clone () const {
        const int tid= omp_get_thread_num();
        return tid == 0 ? *this : the_clones[tid][0];
    }
    bool empty () const { return surf.empty(); }

    STInterfaceCommonDataCL (double t0arg, double t1arg, const VecDescCL& oldls_arg, const VecDescCL& newls_arg, const BndDataCL<>& lsetbnd_arg)
        : old_ls( &oldls_arg), new_ls( &newls_arg), lsetbnd( &lsetbnd_arg), lat( TetraPrismLatticeCL::instance( 2, 1)), ls_loc( lat.vertex_size()), t0( t0arg), t1( t1arg)
    {}

    virtual ~STInterfaceCommonDataCL () {}

    virtual void begin_accumulation   () {
        the_clones= new STInterfaceCommonDataCL*[omp_get_max_threads()];
        the_clones[0]= this;
    }
    virtual void finalize_accumulation() {
        delete[] the_clones;
    }
    virtual void visit (const TetraCL& t) {
        surf.clear();
        q5dom.clear();
        st_local_ls.assign( t, *old_ls, *new_ls, *lsetbnd);
        evaluate_on_vertexes( st_local_ls, lat, Addr( ls_loc));
        if (equalSigns(ls_loc))
            return;
        surf.make_patch<MergeCutPolicyCL>( lat, ls_loc);
        surf.compute_normals( TetraPrismCL( t, t0, t1));
        make_CompositeQuad5DomainSTCodim1SpatialAbsdet( q5dom, surf, TetraPrismCL( t, t0, t1));
    }
    virtual STInterfaceCommonDataCL* clone (int clone_id) {
        return the_clones[clone_id]= new STInterfaceCommonDataCL( *this);
    }
};

template <class LocalMatrixT>
class InterfaceMatrixSTP1AccuCL : public TetraAccumulatorCL
{
  private:
    const STInterfaceCommonDataCL& cdata_;
    std::string name_;

    MatrixCL* mat_; // the matrix
    VectorCL* cpl_; // vector of eliminated initial data
    const STP1P1IdxDescCL* idx_;
    const VectorCL* u_; // the initial data as STP1P1-function

    MatrixBuilderCL* M;

    LocalMatrixT      local_mat;
    LocalNumbSTP1P1CL local_idx;

  public:
    ///\brief For dG in time; does not use cpl_ and u_.
    InterfaceMatrixSTP1AccuCL (MatrixCL* mat, const STP1P1IdxDescCL* idx,
                                 const LocalMatrixT& loc_mat, const STInterfaceCommonDataCL& cdata, std::string name= std::string())
        : cdata_( cdata), name_( name), mat_( mat), cpl_( 0), idx_( idx), u_( 0), M( 0), local_mat( loc_mat) {}
    ///\brief For cG in time; uses cpl_, u_ and changes the test-space.
    InterfaceMatrixSTP1AccuCL (MatrixCL* mat, VectorCL* cpl, const STP1P1IdxDescCL* idx,
                                   const LocalMatrixT& loc_mat, const STInterfaceCommonDataCL& cdata, const VectorCL* u, std::string name= std::string())
        : cdata_( cdata), name_( name), mat_( mat), cpl_( cpl), idx_( idx), u_( u), M( 0), local_mat( loc_mat) {}
    virtual ~InterfaceMatrixSTP1AccuCL () {}

    void set_name (const std::string& n) { name_= n; }

    virtual void begin_accumulation () {
        const IdxT dim= idx_->NumUnknowns() - (cpl_ ? idx_->NumIniUnknowns() : 0);
        std::cout << "InterfaceMatrixSTP1AccuCL::begin_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\"";
        std::cout  << ": " << dim << " rows and cols.\n";
        M= new MatrixBuilderCL( mat_, dim, dim);
    }

    virtual void finalize_accumulation () {
        M->Build();
        delete M;
        M= 0;
        std::cout << "InterfaceMatrixSTP1AccuCL::finalize_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\"";
        std::cout << ": " << mat_->num_nonzeros() << " nonzeros." << std::endl;
    }

    virtual void visit (const TetraCL& t) {
        const STInterfaceCommonDataCL& cdata= cdata_.get_clone();
        if (cdata.empty())
            return;
        local_mat.setup( TetraPrismCL( t, cdata.t0, cdata.t1), cdata);
        local_idx.assign_indices_only( t, *idx_);
        if (cpl_) {
            for (Uint i= 0; i < 4; ++i) // The first four test-functions remain as they are, the last four are made constant in time.
                for (Uint j= 0; j < 8; ++j)
                    local_mat.coup[i + 4][j]+= local_mat.coup[i][j];
            update_global_matrix_and_coupling<LocalMatrixT::ZeroPolicy>( *M, local_mat.coup, local_idx, local_idx, cpl_, u_);

        }
        else
            update_global_matrix<LocalMatrixT::ZeroPolicy>( *M, local_mat.coup, local_idx, local_idx);
    }

    virtual InterfaceMatrixSTP1AccuCL* clone (int /*clone_id*/) { return new InterfaceMatrixSTP1AccuCL( *this); }
};

template <class LocalVectorT>
class InterfaceVectorSTP1AccuCL : public TetraAccumulatorCL
{
  private:
    const STInterfaceCommonDataCL& cdata_;
    std::string name_;

    bool cG_in_t_;

    VectorCL* vec_; // the vector
    const STP1P1IdxDescCL* rowidx_;

    LocalVectorT local_vec;
    LocalNumbSTP1P1CL local_row;

  public:
    InterfaceVectorSTP1AccuCL (VectorCL* vec, const STP1P1IdxDescCL* rowidx,
                               const LocalVectorT& loc_vec, const STInterfaceCommonDataCL& cdata, bool cG_in_t, std::string name= std::string())
        : cdata_( cdata), name_( name), cG_in_t_( cG_in_t), vec_( vec), rowidx_( rowidx), local_vec( loc_vec) {}
    virtual ~InterfaceVectorSTP1AccuCL () {}

    void set_name (const std::string& n) { name_= n; }

    virtual void begin_accumulation () {
        const IdxT num_rows= rowidx_->NumUnknowns() - ( cG_in_t_ ? rowidx_->NumIniUnknowns() : 0);
        std::cout << "STInterfaceVectorSTP1AccuCL::begin_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\"";
        std::cout  << ": " << num_rows << " rows.\n";
    }

    virtual void finalize_accumulation () {}

    virtual void visit (const TetraCL& t) {
        const STInterfaceCommonDataCL& cdata= cdata_.get_clone();
        if (cdata.empty())
            return;
        local_row.assign_indices_only( t, *rowidx_);
        local_vec.setup( TetraPrismCL( t, cdata.t0, cdata.t1), cdata, local_row.num);
        if (cG_in_t_) {
            for (Uint i= 0; i < 4; ++i)
                if (local_row.WithUnknowns( i) && local_row.WithUnknowns( i + 4)) // The last four components are made constant in time.
                    local_vec.vec[i + 4]+= local_vec.vec[i];
            for (int i= 0; i < 8; ++i)
                if (local_row.WithUnknowns( i) && !local_row.IsIni( i))
                    vec_[0][local_row.num[i] - local_row.NumIniUnknowns()]+= local_vec.vec[i];
        }
        else
            for (int i= 0; i < 8; ++i)
                if (local_row.WithUnknowns( i))
                    vec_[0][local_row.num[i]]+= local_vec.vec[i];
    }

    virtual InterfaceVectorSTP1AccuCL* clone (int /*clone_id*/) { return new InterfaceVectorSTP1AccuCL( *this); }
};

/// \brief Compute the load-vector corresponding to the function f on a single tetra-prism.
class LocalVectorSTP1P1CL
{
  private:
    InstatScalarFunction f_;

    std::valarray<double> qp1,
                          qf;

  public:
    double vec[8];

    LocalVectorSTP1P1CL (InstatScalarFunction f) : f_( f) {}

    void setup (const TetraPrismCL& prism, const STInterfaceCommonDataCL& cdata, const IdxT numr[8]) {
        resize_and_evaluate_on_vertexes( f_, prism, cdata.q5dom, qf);
        qp1.resize( cdata.q5dom.vertex_size());
        for (Uint i= 0; i < 8; ++i) {
                if (numr[i] == NoIdx)
                    continue;
                evaluate_on_vertexes( STP1P1DiscCL::ref_val[i], cdata.q5dom, Addr( qp1));
                vec[i]= quad_codim1( qf*qp1, cdata.q5dom);
        }
    }
};

class LocalSpatialInterfaceMassSTP1P1CL
{
  private:
    const InterfaceCommonDataP1CL& spatialcdata_;
    LocalInterfaceMassP1CL spatial_mass_;
    const bool on_old_iface_;

  public:
    double coup[8][8];
    static const ZeroPolicyEnum ZeroPolicy= RemoveExactLocalZeros;

    LocalSpatialInterfaceMassSTP1P1CL (const InterfaceCommonDataP1CL& spatialcdata, bool on_old_iface= true)
        : spatialcdata_( spatialcdata), on_old_iface_( on_old_iface)
    { std::memset( coup, 0, 8*8*sizeof(double)); }

    void setup (const TetraPrismCL& p, const STInterfaceCommonDataCL&) {
        spatial_mass_.setup( p.t, spatialcdata_.get_clone());
        // on_old_iface_==true: The basis functions for t1 are all zero on the interface at t0
        //                      --> zero-init in the constructor.
        // else: The basis functions for t0 are all zero on the interface at t1.
        // For use_cG_in_time_ == true && use_massdiv_ == false: The rows [4..7) are modified, thus:
        std::memset( coup, 0, 8*8*sizeof(double));
        const Uint offset= on_old_iface_ ? 0 : 4;
        for (Uint i= 0; i < 4 ; ++i) {
            coup[i + offset][i + offset]= spatial_mass_.coup[i][i];
            for(Uint j= 0; j < i; ++j)
                coup[i + offset][j + offset]= coup[j + offset][i + offset]= spatial_mass_.coup[i][j];
        }
    }
};

void
resize_and_scatter_piecewise_spatial_normal (const SPatchCL<4>& surf, const QuadDomainCodim1CL<4>& qdom, std::valarray<Point3DCL>& spatial_normal);

class LocalLaplaceBeltramiSTP1P1CL
{
  private:
    double D_; // diffusion coefficient

    LocalSTP1CL<Point3DCL> gradx[8];
    double dummy;
    GridFunctionCL<Point3DCL> qgradx,
                              q[8],
                              spatial_n;
    QuadDomainCodim1CL<4> qdom;

  public:
    double coup[8][8];
    static const ZeroPolicyEnum ZeroPolicy= KeepLocalZeros;

    void setup (const TetraPrismCL& prism, const STInterfaceCommonDataCL& cdata) {
        make_CompositeQuad2DomainSTCodim1SpatialAbsdet( qdom, cdata.surf, prism);
        resize_and_scatter_piecewise_spatial_normal( cdata.surf, qdom, spatial_n);

        STP1P1DiscCL::GetSpatialGradients( prism, gradx);
        qgradx.resize( qdom.vertex_size());
        for (int i= 0; i < 8; ++i) {
            evaluate_on_vertexes( gradx[i], qdom, Addr( qgradx));
            q[i].resize( qdom.vertex_size());
            q[i]= qgradx - dot( qgradx, spatial_n)*spatial_n;
        }

        for (int i= 0; i < 8; ++i) {
            coup[i][i]= D_* quad_codim1( dot( q[i], q[i]), qdom);
            for(int j= 0; j < i; ++j)
                coup[i][j]= coup[j][i]= D_* quad_codim1( dot( q[i], q[j]), qdom);
        }
    }

    LocalLaplaceBeltramiSTP1P1CL (double D)
        :D_( D) {}
};

template <class DiscVelSolT>
class LocalMaterialDerivativeSTP1P1CL
{
  private:
    const DiscVelSolT w_; // wind
    LocalSTP2P1CL<Point3DCL> loc_w_;

    const bool transpose_local_matrix_;

    LocalSTP1CL<Point4DCL> grad[8];
    double dummy;
    GridFunctionCL<Point4DCL> qw,
                              qgrad;
    std::valarray<double>     q[8],
                              qwdotgrad[8];

  public:
    double coup[8][8];
    static const ZeroPolicyEnum ZeroPolicy= KeepLocalZeros;

    void setup (const TetraPrismCL& prism, const STInterfaceCommonDataCL& cdata) {
        loc_w_.assign( prism.t, w_);
        resize_and_evaluate_on_vertexes( STWindProxyCL( loc_w_), cdata.q5dom, qw);
        qgrad.resize( cdata.q5dom.vertex_size());
        STP1P1DiscCL::GetGradients( prism, grad);
        for (int i= 0; i < 8; ++i) {
            evaluate_on_vertexes( grad[i], cdata.q5dom, Addr( qgrad));
            qwdotgrad[i].resize( cdata.q5dom.vertex_size());
            qwdotgrad[i]= dot( qw, qgrad);
            resize_and_evaluate_on_vertexes( STP1P1DiscCL::ref_val[i], cdata.q5dom, q[i]);
        }

        for (int i= 0; i < 8; ++i) {
            for(int j= 0; j < 8; ++j)
                (transpose_local_matrix_ ? coup[j][i] : coup[i][j])= quad_codim1( qwdotgrad[j]*q[i], cdata.q5dom);
        }
    }

    LocalMaterialDerivativeSTP1P1CL (const DiscVelSolT& w, const bool transpose_local_matrix= false)
        :w_( w), transpose_local_matrix_( transpose_local_matrix) {}
};

template <class DiscVelSolT>
class LocalMassdivSTP1P1CL
{
  private:
    const DiscVelSolT                w_; // wind
    LocalSTP2P1CL<Point3DCL>         loc_w_;
    LocalSTP1P1CL< SMatrixCL<3,3> >  dw;
    GridFunctionCL< SMatrixCL<3,3> > qdw;
    std::valarray<double> qdivgamma_w;

    GridFunctionCL<Point3DCL> spatial_n;

    double dummy;
    SMatrixCL<3,3> T;
    LocalP1CL<Point3DCL> gradrefp2[10],
                         gradp2[10];

    std::valarray<double> q[8];

  public:
    double coup[8][8];
    static const ZeroPolicyEnum ZeroPolicy= KeepLocalZeros;

    void setup (const TetraPrismCL& prism, const STInterfaceCommonDataCL& cdata) {
        resize_and_scatter_piecewise_spatial_normal( cdata.surf, cdata.q5dom, spatial_n);
        loc_w_.assign( prism.t, w_);
        GetTrafoTr( T, dummy, prism.t);
        P2DiscCL::GetGradients( gradp2, gradrefp2, T);

        dw= LocalSTP1P1CL< SMatrixCL<3,3> >();
        for (int i= 0; i < 10; ++i) {
            dw.at_t0()+= outer_product(loc_w_.at_t0()[i], gradp2[i]);
            dw.at_t1()+= outer_product(loc_w_.at_t1()[i], gradp2[i]);
        }
        resize_and_evaluate_on_vertexes( dw, cdata.q5dom, qdw);
        qdivgamma_w.resize( cdata.q5dom.vertex_size());
        qdivgamma_w= trace( qdw) - dot( spatial_n, qdw, spatial_n);

        for (int i= 0; i < 8; ++i)
            resize_and_evaluate_on_vertexes( STP1P1DiscCL::ref_val[i], cdata.q5dom, q[i]);

        for (int i= 0; i < 8; ++i) {
            coup[i][i]= quad_codim1( qdivgamma_w*q[i]*q[i], cdata.q5dom);
            for(int j= 0; j < i; ++j)
                coup[j][i]= coup[i][j]= quad_codim1( qdivgamma_w*q[i]*q[j], cdata.q5dom);
        }
    }

    LocalMassdivSTP1P1CL (const DiscVelSolT& w)
        :w_( w) { P2DiscCL::GetGradientsOnRef( gradrefp2); }
};

/// \brief Convenience-function to reduce the number of explicit template-parameters for the spacetime-massdiv- and the -convection-matrix.
template <template <class> class LocalMatrixT, class DiscVelSolT>
  inline InterfaceMatrixSTP1AccuCL< LocalMatrixT<DiscVelSolT> >*
  make_wind_dependent_matrixSTP1P1_accu (MatrixCL* mat, const STP1P1IdxDescCL* idx,
                                         const STInterfaceCommonDataCL& cdata, const DiscVelSolT& wind, std::string name= std::string())
{
    return new InterfaceMatrixSTP1AccuCL< LocalMatrixT<DiscVelSolT> >( mat, idx,
        LocalMatrixT<DiscVelSolT>( wind), cdata, name);
}

/// \brief Convenience-function to reduce the number of explicit template-parameters for the spacetime-massdiv- and the -convection-matrix.
template <template <class> class LocalMatrixT, class DiscVelSolT>
  inline InterfaceMatrixSTP1AccuCL< LocalMatrixT<DiscVelSolT> >*
  make_wind_dependent_local_transpose_matrixSTP1P1_accu (MatrixCL* mat, const STP1P1IdxDescCL* idx,
                                         const STInterfaceCommonDataCL& cdata, const DiscVelSolT& wind, std::string name= std::string())
{
    return new InterfaceMatrixSTP1AccuCL< LocalMatrixT<DiscVelSolT> >( mat, idx,
        LocalMatrixT<DiscVelSolT>( wind, true), cdata, name);
}

/// \brief Convenience-function to reduce the number of explicit template-parameters for the spacetime-massdiv- and the -convection-matrix.
template <template <class> class LocalMatrixT, class DiscVelSolT>
  inline InterfaceMatrixSTP1AccuCL< LocalMatrixT<DiscVelSolT> >*
  make_wind_dependent_matrixSTP1P0_1_accu (MatrixCL* mat, VectorCL* cpl, const STP1P1IdxDescCL* idx,
                                           const STInterfaceCommonDataCL& cdata, const VectorCL* u_ini, const DiscVelSolT& wind, std::string name= std::string())
{
    return new InterfaceMatrixSTP1AccuCL< LocalMatrixT<DiscVelSolT> >( mat, cpl, idx,
        LocalMatrixT<DiscVelSolT>( wind), cdata, u_ini, name);
}

/// \brief Convenience-function to reduce the number of explicit template-parameters for the spacetime-massdiv- and the -convection-matrix.
template <template <class> class LocalMatrixT, class DiscVelSolT>
  inline InterfaceMatrixSTP1AccuCL< LocalMatrixT<DiscVelSolT> >*
  make_wind_dependent_local_transpose_matrixSTP1P0_1_accu (MatrixCL* mat, VectorCL* cpl, const STP1P1IdxDescCL* idx,
                                           const STInterfaceCommonDataCL& cdata, const VectorCL* u_ini, const DiscVelSolT& wind, std::string name= std::string())
{
    return new InterfaceMatrixSTP1AccuCL< LocalMatrixT<DiscVelSolT> >( mat, cpl, idx,
        LocalMatrixT<DiscVelSolT>( wind, true), cdata, u_ini, name);
}


/// \brief Space-time P1-discretization and solution of the transport equation on the interface
/// Two methods are implemented:
///     * cG_in_t_ == false: dG-coupling in time; both, the trial- and the test-space, are the full STP1P1 FE-space.
///     * cG_in_t_ == true: The initial values are eliminated; additionally, the test space is changed: for Ini-Unknowns, there remains only a test-function, which is constant in time.
class SurfactantSTP1CL : public SurfactantP1BaseCL
{
  public:
    MatrixCL A,    ///< ST-diffusion matrix
             Mder, ///< ST-material-derivative matrix
             Mdiv, ///< ST-mass-matrix with interface-divergence of velocity
             Mold, ///< mass matrix on the old spatial interface; only used for cG_in_t_ == false.
             Mnew; ///< mass matrix on the new spatial interface; only used for cG_in_t_ == false and use_mass_div_ == false.

    VectorCL load, ///< load-vector
             cpl_A_,   ///< The cpl-Vectors are used only for cG_in_t_ == true.
             cpl_der_,
             cpl_div_,
             cpl_old_;

    size_t dim; ///< Dimension of the linear system.

  private:
    bool cG_in_t_,
         use_mass_div_;

    STP1P1IdxDescCL st_idx_;
    VectorCL st_oldic_, ///< the old solution represented in the full space-time-FE-basis.
             st_ic_;    ///< the new solution represented in the space-time-FE-basis with eliminated initial values.

    void Update_cG();
    void Update_dG();

  public:
    SurfactantSTP1CL (MultiGridCL& mg,
        double theta, double D, VecDescCL* v, const VelBndDataT& Bnd_v, VecDescCL& lset_vd, const BndDataCL<>& lsetbnd, bool cG_in_t, bool use_mass_div,
        int iter= 1000, double tol= 1e-7, double omit_bound= -1.)
    : SurfactantP1BaseCL( mg, theta, D, v, Bnd_v, lset_vd, lsetbnd, iter, tol, omit_bound), cG_in_t_( cG_in_t), use_mass_div_( use_mass_div)
    {}

    /// save a copy of the old level-set and velocity; moves ic to oldic; must be called before DoStep.
    // virtual void InitTimeStep (); // as in the base

    /// perform one time step
    virtual void DoStep (double new_t);

    /// \name For internal use only
    /// The following member functions are added to enable an easier implementation
    /// of the coupling navstokes-levelset. They should not be called by a common user.
    /// Use DoStep() instead.
    ///@{
    void InitStep (double new_t);
    void DoStep ();
    void CommitStep ();
    void Update ();
    ///@}
};


/// =Methods with transport on the domain=

class TransportP1FunctionCL; ///< forward declaration


/// \brief P1-discretization and solution of the transport equation on the interface
class SurfactantCharTransportP1CL: public SurfactantP1BaseCL
{
  public:
    IdxDescCL full_idx;
    MatDescCL A,  ///< diffusion matrix
              M,  ///< mass matrix
              C,  ///< convection matrix
              Md; ///< mass matrix with interface-divergence of velocity

    VectorCL load, ///< for a load-function
             rhs_; ///< for the transported initial data

  private:
    MatrixCL      L_; ///< sum of matrices
    TransportP1FunctionCL* fulltransport_;

  public:
    SurfactantCharTransportP1CL (MultiGridCL& mg,
        double theta, double D, VecDescCL* v, const VelBndDataT& Bnd_v, VecDescCL& lset_vd, const BndDataCL<>& lsetbnd,
        int iter= 1000, double tol= 1e-7, double omit_bound= -1.)
        : SurfactantP1BaseCL( mg, theta, D, v, Bnd_v, lset_vd, lsetbnd, iter, tol, omit_bound),
          full_idx( P1_FE), fulltransport_( 0)
    {}

    /// \remarks call SetupSystem \em before calling SetTimeStep!
    void SetTimeStep( double dt, double theta=-1);

    /// perform one time step
    void DoStep (double new_t);

    const_DiscSolCL GetSolution() const
        { return const_DiscSolCL( &ic, &Bnd_, &MG_); }
    const_DiscSolCL GetSolution( const VecDescCL& Myic) const
        { return const_DiscSolCL( &Myic, &Bnd_, &MG_); }
    ///@}

    /// \name For internal use only
    /// The following member functions are added to enable an easier implementation
    /// of the coupling navstokes-levelset. They should not be called by a common user.
    /// Use DoStep() instead.
    ///@{
    void InitStep (double new_t);
    void DoStep ();
    void CommitStep ();
    void Update ();
    ///@}
};

} // end of namespace DROPS

#include "surfactant/ifacetransp.tpp"

#endif

