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
#include "levelset/levelsetmapper.h"
#include "levelset/mgobserve.h"
#include "out/ensightOut.h"
#include "out/vtkOut.h"
#include "num/fe_repair.h"
#include "num/quadrature.h"

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

/// \brief Copies P1IF_FE-unknown-indices or P2IF_FE-indices from idx on s into Numb.
/// Non-existent dofs get NoIdx.
void GetLocalNumbInterface(IdxT* Numb, const TetraCL& s, const IdxDescCL& idx);

/// \todo This should be a generic function somewhere in num or misc.
void P1Init (instat_scalar_fun_ptr icf, VecDescCL& ic, const MultiGridCL& mg, double t);

/// \brief Resize normal according to qdom and fill in surf.normal. The normal must be precomputed.
template <Uint Dim>
  void
  resize_and_scatter_piecewise_normal (const SPatchCL<Dim>& surf, const QuadDomainCodim1CL<Dim>& qdom, std::valarray<typename SPatchCL<Dim>::WorldVertexT>& normal);


/// \brief The routine sets up the mass-matrix in matM on the interface defined by ls.
///        It belongs to the FE induced by standard P1-elements.
void SetupInterfaceMassP1 (const MultiGridCL& mg, MatDescCL* matM, const VecDescCL& ls, const BndDataCL<>& lsetbnd);

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
        if (equal_signs( ls_loc))
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
  evaluate_on_vertexes (T (*f)(const Point3DCL&, double), const TetraBaryPairVectorT& pos, double t, ResultIterT result_iterator);

template <class T, class ResultContT>
  inline ResultContT&
  resize_and_evaluate_on_vertexes (T (*f)(const Point3DCL&, double), const TetraBaryPairVectorT& pos, double t, ResultContT& result_container);

template <class PEvalT, class ResultIterT>
  inline ResultIterT
  evaluate_on_vertexes (const PEvalT& f, const TetraBaryPairVectorT& pos, ResultIterT result_iterator);

template <class PEvalT, class ResultContT>
  inline ResultContT&
  resize_and_evaluate_on_vertexes (const PEvalT& f, const TetraBaryPairVectorT& pos, ResultContT& result_container);


class InterfaceCommonDataP2CL : public TetraAccumulatorCL
{
  private:
    InterfaceCommonDataP2CL** the_clones;

    const VecDescCL*   ls;      // P2-level-set
    const BndDataCL<>* lsetbnd; // boundary data for the level set function

    const PrincipalLatticeCL* lat;

  public:
    /// common data @{
    LocalP2CL<> locp2_ls;

    LocalP2CL<>          p2[10];
    LocalP1CL<Point3DCL> gradrefp2[10];

    std::valarray<double>     ls_loc;
    SurfacePatchCL            surf;
    QuadDomain2DCL            qdom;
    TetraBaryPairVectorT      qdom_projected;
    std::valarray<double>     absdet;
    QuaQuaMapperCL            quaqua;

    const PrincipalLatticeCL& get_lattice () const { return *lat; }
    /// @}

    LocalP2CL<> get_local_p2_ls (const TetraCL& t) const { return LocalP2CL<>( t, *ls, *lsetbnd); };

    const InterfaceCommonDataP2CL& get_clone () const {
        const int tid= omp_get_thread_num();
        return tid == 0 ? *this : the_clones[tid][0];
    }

    bool empty () const { return surf.empty(); }

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
        if (equal_signs( ls_loc))
            return;
        surf.make_patch<MergeCutPolicyCL>( *lat, ls_loc);
        if (surf.empty())
            return;

        make_CompositeQuad5Domain2D ( qdom, surf, t);
        qdom_projected.clear();
        qdom_projected.reserve( qdom.vertex_size());
        absdet.resize( qdom.vertex_size());
        const TetraCL* tet;
        BaryCoordCL b;
        QuadDomain2DCL::const_vertex_iterator v= qdom.vertex_begin();
        for (Uint i= 0; i < qdom.vertex_size(); ++i, ++v) {
            tet= &t;
            b= *v;
            quaqua.base_point( tet, b);
            qdom_projected.push_back( std::make_pair( tet, b));
            try { // Ignore very small triangles by setting absdet to 0.
                absdet[i]= abs_det( t, *v, *tet, b, quaqua, surf);
            } catch (const DROPSErrCL&) {
//               t.DebugInfo( std::cerr);
                absdet[i]= 0.; // Implied by absdet.resize above the loop...
            }

        }
    }

    virtual InterfaceCommonDataP2CL* clone (int clone_id) {
        return the_clones[clone_id]= new InterfaceCommonDataP2CL( *this);
    }
};

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
            const TetraBaryPairT& p= cdata.qdom_projected[i];
            const double newweight= cdata.qdom.weight_begin()[i]*cdata.absdet[i];
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
        std::cout << "InterfaceMatrixAccuCL::finalize_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\"";
        std::cout << ": " << mat_->Data.num_nonzeros() << " nonzeros." << std::endl;
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
        std::cout << "InterfaceVectorAccuCL::begin_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\"";
        std::cout  << ": " << y_->RowIdx->NumUnknowns() << " rows.\n";
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


/// \brief Compute the load-vector corresponding to the function f on a single tetra.
class LocalVectorP1CL
{
  private:
    instat_scalar_fun_ptr f_;
    double time_;

    LocalP1CL<> p1[4];
    std::valarray<double> qp1,
                          qf;
    QuadDomain2DCL qdom;

  public:
    static const FiniteElementT row_fe_type= P1IF_FE;

    double vec[4];

    LocalVectorP1CL (instat_scalar_fun_ptr f, double time) : f_( f), time_( time) { p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; }

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


class LocalInterfaceMassP1CL
{
  private:
    std::valarray<double> q[4];
    QuadDomain2DCL qdom;

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
                coup[i][j]= coup[j][i]= quad_2D( q[j]*q[i], qdom);
        }
    }
};


/// \brief The routine sets up the Laplace-Beltrami-matrix in mat on the interface defined by ls.
///        It belongs to the FE induced by standard P1-elements.
///
/// D is the diffusion-coefficient
void SetupLBP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsbnd, double D);

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
            coup[i][i]= D_* /*area of reference triangle*/ 0.5*(dot(q[i], q[i])*absdet).sum();
            for(int j= 0; j < i; ++j)
                coup[i][j]= coup[j][i]= D_* /*area of reference triangle*/ 0.5*(dot(q[i], q[j])*absdet).sum();
        }
    }

    LocalLaplaceBeltramiP1CL (double D)
        :D_( D) {}
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

    LocalInterfaceConvectionP1CL (const DiscVelSolT& w)
        :  w_( w) {}
};

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

/// \brief Compute the P2 load vector corresponding to the function f on a single tetra.
class LocalVectorP2CL
{
  private:
    instat_scalar_fun_ptr f_;
    double time_;

    std::valarray<double> qp2,
                          qf;

  public:
    static const FiniteElementT row_fe_type= P2IF_FE;

    double vec[10];

    LocalVectorP2CL (instat_scalar_fun_ptr f, double time) : f_( f), time_( time) {}

    void setup (const TetraCL&, const InterfaceCommonDataP2CL& cdata, const IdxT numr[10]) {
        resize_and_evaluate_on_vertexes( f_, cdata.qdom_projected, time_, qf);
        qp2.resize( cdata.qdom.vertex_size());
        for (Uint i= 0; i < 10; ++i) {
                if (numr[i] == NoIdx)
                    continue;
                evaluate_on_vertexes( cdata.p2[i], cdata.qdom, Addr( qp2));
                vec[i]= quad_2D( cdata.absdet*qf*qp2, cdata.qdom);
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
            coup[i][i]= quad_2D( cdata.absdet*dot( qgradp2[i], qgradp2[i]), cdata.qdom);
            for(int j= 0; j < i; ++j)
                coup[i][j]= coup[j][i]= quad_2D( cdata.absdet*dot( qgradp2[j], qgradp2[i]), cdata.qdom);
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
            coup[i][i]= quad_2D( cdata.absdet*qp2[i]*qp2[i], cdata.qdom);
            for(int j= 0; j < i; ++j)
                coup[i][j]= coup[j][i]= quad_2D( cdata.absdet*qp2[j]*qp2[i], cdata.qdom);
        }
    }

    LocalMassP2CL () {}
};

/// \brief The routine sets up the load-vector in v on the interface defined by ls.
///        It belongs to the FE induced by standard P1-elements.
void SetupInterfaceRhsP1 (const MultiGridCL& mg, VecDescCL* v,
    const VecDescCL& ls, const BndDataCL<>& lsbnd, instat_scalar_fun_ptr f);


///\brief Initialize the QuadDomain2DCL-object qdom for quadrature with Quad5_2DCL on the lattice lat of t, given the level set in ls and bnd.
inline const QuadDomain2DCL&
make_CompositeQuad5Domain2D (QuadDomain2DCL& qdom, const TetraCL& t, const PrincipalLatticeCL& lat, const DROPS::VecDescCL& ls, const DROPS::BndDataCL<>& bnd)
{
    LocalP2CL<> locp2_ls( t, ls, bnd);
    std::valarray<double> ls_loc( lat.vertex_size());
    evaluate_on_vertexes( locp2_ls, lat, Addr( ls_loc));
    if (equal_signs( ls_loc)) {
        qdom.clear();
        return qdom;
    }
    SurfacePatchCL surf;
    surf.make_patch<MergeCutPolicyCL>( lat, ls_loc);
    return make_CompositeQuad5Domain2D ( qdom, surf, t);
}

/// \brief Short-hand for integral on the interface.
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

    instat_matrix_fun_ptr ref_dp;
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
    void set_ref_dp   ( instat_matrix_fun_ptr rdp) { ref_dp= rdp; }
    void set_ref_abs_det   ( double (*rad) (const TetraCL& t, const BaryCoordCL& b, const SurfacePatchCL& surf)) { ref_abs_det= rad; }

    InterfaceDebugP2CL (const InterfaceCommonDataP2CL& cdata);
    virtual ~InterfaceDebugP2CL () {}

    virtual void begin_accumulation   ();
    virtual void finalize_accumulation();

    virtual void visit (const TetraCL& t);

    virtual InterfaceDebugP2CL* clone (int /*clone_id*/) { return new InterfaceDebugP2CL( *this); }
};


/// \brief P1-discretization and solution of the transport equation on the interface
class SurfactantP1BaseCL
{
  public:
    typedef BndDataCL<Point3DCL>                              VelBndDataT;
    typedef NoBndDataCL<>                                     BndDataT;
    typedef P1EvalCL<double, const BndDataT, VecDescCL>       DiscSolCL;
    typedef P1EvalCL<double, const BndDataT, const VecDescCL> const_DiscSolCL;

    IdxDescCL idx; ///< index desctription for concentration at current time
    VecDescCL ic;  ///< concentration on the interface at current time

  protected:
    MultiGridCL&  MG_;
    double        D_,     ///< diffusion coefficient
                  theta_, ///< time scheme parameter
                  dt_;    ///< time step size
    instat_scalar_fun_ptr rhs_fun_; ///< function for a right-hand side

    BndDataT            Bnd_;    ///< Dummy boundary data for interface solution

    const VelBndDataT&  Bnd_v_;  ///< Boundary condition for the velocity
    VecDescCL*          v_;      ///< velocity at current time step
    VecDescCL&          lset_vd_;///< levelset at current time step
    const BndDataCL<>&  lsetbnd_;///< level set boundary

    IdxDescCL           oldidx_; ///< idx that corresponds to old time (and oldls_)
    VectorCL            oldic_;  ///< interface concentration at old time
    VecDescCL           oldls_;  ///< levelset at old time
    VecDescCL           oldv_;   ///< velocity at old time
    double              oldt_;   ///< old time

    GSPcCL                  pc_;
    GMResSolverCL<GSPcCL>   gm_;
    __UNUSED__ double omit_bound_; ///< not used atm

  public:
    SurfactantP1BaseCL (MultiGridCL& mg,
        double theta, double D, VecDescCL* v, const VelBndDataT& Bnd_v, VecDescCL& lset_vd, const BndDataCL<>& lsetbnd,
        int iter= 1000, double tol= 1e-7, double omit_bound= -1.)
    : idx( P1IF_FE), MG_( mg), D_( D), theta_( theta), dt_( 0.), rhs_fun_( 0),
        Bnd_v_( Bnd_v), v_( v), lset_vd_( lset_vd), lsetbnd_( lsetbnd), oldidx_( P1IF_FE), gm_( pc_, 100, iter, tol, true),
        omit_bound_( omit_bound)
    { idx.GetXidx().SetBound( omit_bound); }
    virtual ~SurfactantP1BaseCL () {}

    const MultiGridCL& GetMG() const { return MG_; }
    GMResSolverCL<GSPcCL>& GetSolver() { return gm_; }

    const_DiscSolCL GetSolution() const
        { return const_DiscSolCL( &ic, &Bnd_, &MG_); }
    const_DiscSolCL GetSolution( const VecDescCL& Myic) const
        { return const_DiscSolCL( &Myic, &Bnd_, &MG_); }

    /// initialize the interface concentration
    void SetInitialValue (instat_scalar_fun_ptr, double t= 0.);

    /// set the parameter of the theta-scheme for time stepping
    void SetRhs (instat_scalar_fun_ptr);

    /// set the parameter of the theta-scheme for time stepping
    void SetTheta (double theta);

    /// save a copy of the old level-set and velocity; moves ic to oldic; must be called before DoStep.
    virtual void InitTimeStep ();

    /// perform one time step to new_t.
    virtual void DoStep (double /*new_t*/) {}
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
    std::auto_ptr<RepairP1CL<double, NoBndDataCL>::type > p1repair_;

  public:
    InterfaceP1RepairCL (MultiGridCL& mg, const VecDescCL& lset_vd, const BndDataCL<>& lset_bnd, VecDescCL& u)
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

///\brief Represents a scalar P1 function on the interface as VTK variable by extension to the
///       whole domain.
class VTKIfaceScalarCL : public VTKVariableCL
{
  private:
    const VecDescCL&   u_;
    MultiGridCL&       mg_;

  public:
    VTKIfaceScalarCL (MultiGridCL& mg, const VecDescCL& u, std::string varName)
        : VTKVariableCL( varName), u_( u), mg_( mg) {}

    void put      (VTKOutCL& cf) const;
    Uint GetDim() const { return 1; }
};

///\brief Create an VTKIfaceP1ScalarCL with operator new.
///
/// This is just for uniform code; the analogous functions for scalars and vectors are more useful
/// because they help to avoid template parameters in user code.
inline VTKIfaceScalarCL&
make_VTKIfaceScalar (MultiGridCL& mg, const VecDescCL& u,
    std::string varName)
{
    return *new VTKIfaceScalarCL( mg, u, varName);
}

/// ==Space-Time-methods==
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
    typedef T (*fun_type)(const Point3DCL&, double);
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
  evaluate_on_vertexes (T (*f)(const Point3DCL&, double), const TetraPrismCL& prism, const DomainT& dom, ResultIterT result_iterator)
{
    return std::transform( dom.vertex_begin(), dom.vertex_end(), result_iterator, STCoordEvalCL<T>( prism, f));
}

template <class T, class DomainT, class ResultContT>
  inline const ResultContT&
  resize_and_evaluate_on_vertexes (T (*f)(const Point3DCL&, double), const TetraPrismCL& prism, const DomainT& dom, ResultContT& result_container)
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
        if (equal_signs( ls_loc))
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
    instat_scalar_fun_ptr f_;

    std::valarray<double> qp1,
                          qf;

  public:
    double vec[8];

    LocalVectorSTP1P1CL (instat_scalar_fun_ptr f) : f_( f) {}

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

