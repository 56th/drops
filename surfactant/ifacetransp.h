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
#include "levelset/mgobserve.h"
#include "out/ensightOut.h"
#include "out/vtkOut.h"
#include "num/fe_repair.h"
#include "num/quadrature.h"

#ifndef DROPS_IFACETRANSP_H
#define DROPS_IFACETRANSP_H

namespace DROPS {

/// \brief Given a P1-vecdesc on the interface x and a P1-vecdesc xext on mg, the values of x
///        are copied to xext and the remaining values of xext are set to 0.
void Extend (const MultiGridCL& mg, const VecDescCL& x, VecDescCL& xext);

/// \brief Given a P1-vecdesc on the interface x and a P1-vecdesc xext on mg, the values of xext
///        are copied to x, where x has an unknown.
void Restrict (const MultiGridCL& mg, const VecDescCL& xext, VecDescCL& x);

/// \brief Helper for the accumulators: inserts the local matrix coup into M.
void update_global_matrix_P1 (MatrixBuilderCL& M, const double coup[4][4], const IdxT numr[4], const IdxT numc[4]);

/// \brief Writes a normal to each triangle of p to n. The normal is the normalized cross-product of two of the edges. Absdet is its norm.
template <class ResultContainerT>
  void
  resize_and_evaluate_piecewise_normal (const SurfacePatchCL& p, const TetraCL& t, ResultContainerT& n, std::valarray<double>* absdet= 0);

/// \brief Writes a normal to each facet of p to n. The normal has unit length. Absdet is the scaling of the volume from the reference tetra.
template <class ResultContainerT>
  void
  resize_and_evaluate_piecewise_normal (const SPatchCL<4>& p, const TetraPrismCL& prism, ResultContainerT& n, std::valarray<double>* absdet= 0);

/// \todo This should be a generic function somewhere in num or misc.
void P1Init (instat_scalar_fun_ptr icf, VecDescCL& ic, const MultiGridCL& mg, double t);


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

  public:
    const PrincipalLatticeCL& lat;
    LocalP1CL<> p1[4];

    const VecDescCL*   ls;      // P2-level-set
    const BndDataCL<>* lsetbnd; // boundary data for the level set function

    LocalP2CL<> locp2_ls;
    std::valarray<double> ls_loc;
    SurfacePatchCL surf;

    const InterfaceCommonDataP1CL& get_clone () const {
        const int tid= omp_get_thread_num();
        return tid == 0 ? *this : the_clones[tid][0];
    }

    bool empty () const { return surf.empty(); }

    InterfaceCommonDataP1CL (const VecDescCL& ls_arg, const BndDataCL<>& lsetbnd_arg)
        : lat( PrincipalLatticeCL::instance( 2)), ls( &ls_arg), lsetbnd( &lsetbnd_arg), ls_loc( lat.vertex_size())
    { p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; } // P1-Basis-Functions
    InterfaceCommonDataP1CL ()
        : lat( PrincipalLatticeCL::instance( 2)), ls( 0), lsetbnd( 0)
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

template <class LocalMatrixT>
class InterfaceMatrixAccuP1CL : public TetraAccumulatorCL
{
  private:
    const InterfaceCommonDataP1CL& cdata_;
    std::string name_;

    MatDescCL* mat_; // the matrix
    MatrixBuilderCL* M;

    LocalMatrixT local_mat;

    Uint lvl;
    IdxT numr[4],
         numc[4];

  public:
    InterfaceMatrixAccuP1CL (MatDescCL* Mmat, const LocalMatrixT& loc_mat, const InterfaceCommonDataP1CL& cdata,
                             std::string name= std::string())
        : cdata_( cdata), name_( name), mat_( Mmat), M( 0), local_mat( loc_mat) {}
    virtual ~InterfaceMatrixAccuP1CL () {}

    void set_name (const std::string& n) { name_= n; }

    virtual void begin_accumulation () {
        const IdxT num_rows= mat_->RowIdx->NumUnknowns();
        const IdxT num_cols= mat_->ColIdx->NumUnknowns();
        std::cout << "InterfaceMatrixAccuP1CL::begin_accumulation";
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
        std::cout << "InterfaceMatrixAccuP1CL::finalize_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\"";
        std::cout << ": " << mat_->Data.num_nonzeros() << " nonzeros." << std::endl;
    }

    virtual void visit (const TetraCL& t) {
        const InterfaceCommonDataP1CL& cdata= cdata_.get_clone();
        if (cdata.empty())
            return;
        local_mat.setup( t, cdata);
        GetLocalNumbP1NoBnd( numr, t, *mat_->RowIdx);
        GetLocalNumbP1NoBnd( numc, t, *mat_->ColIdx);
        update_global_matrix_P1( *M, local_mat.coup, numr, numc);
    }

    virtual InterfaceMatrixAccuP1CL* clone (int /*clone_id*/) { return new InterfaceMatrixAccuP1CL( *this); }
};

/// \brief Accumulate an interface-vector.
template <class LocalVectorT>
class InterfaceVectorAccuP1CL : public TetraAccumulatorCL
{
  private:
    const InterfaceCommonDataP1CL& cdata_;
    std::string name_;

    VecDescCL* y_;
    LocalVectorT local_vec;
    IdxT numry[4];

  public:
    InterfaceVectorAccuP1CL (VecDescCL* y, const LocalVectorT& loc_vec, const InterfaceCommonDataP1CL& cdata,
                             std::string name= std::string())
        : cdata_( cdata), name_( name), y_( y), local_vec( loc_vec) {}
    virtual ~InterfaceVectorAccuP1CL () {}

    void set_name (const std::string& n) { name_= n; }

    virtual void begin_accumulation () {
        std::cout << "InterfaceVectorAccuP1CL::begin_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\"";
        std::cout  << ": " << y_->RowIdx->NumUnknowns() << " rows.\n";
    }

    virtual void finalize_accumulation() {}

    virtual void visit (const TetraCL& t) {
        const InterfaceCommonDataP1CL& cdata= cdata_.get_clone();
        if (cdata.empty())
            return;
        GetLocalNumbP1NoBnd( numry, t, *y_->RowIdx);
        local_vec.setup( t, cdata, numry);
        for (int i= 0; i < 4; ++i)
            if (numry[i] != NoIdx)
                y_->Data[numry[i]]+= local_vec.vec[i];
    }

    virtual InterfaceVectorAccuP1CL* clone (int /*clone_id*/) { return new InterfaceVectorAccuP1CL( *this); }
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
  inline InterfaceVectorAccuP1CL<LocalMatVecP1CL< LocalMatrixT<DiscVelSolT> > >*
  make_wind_dependent_vectorP1_accu (VecDescCL* y, const VecDescCL* x, const InterfaceCommonDataP1CL& cdata, const DiscVelSolT& wind, std::string name= std::string())
{
    return new InterfaceVectorAccuP1CL<LocalMatVecP1CL< LocalMatrixT<DiscVelSolT> > >( y,
        LocalMatVecP1CL< LocalMatrixT<DiscVelSolT> >( LocalMatrixT<DiscVelSolT>( wind), x), cdata, name);
}

/// \brief Convenience-function to reduce the number of explicit template-parameters for the massdiv- and the convection-matrix.
template <template <class> class LocalMatrixT, class DiscVelSolT>
  inline InterfaceMatrixAccuP1CL< LocalMatrixT<DiscVelSolT> >*
  make_wind_dependent_matrixP1_accu (MatDescCL* mat, const InterfaceCommonDataP1CL& cdata, const DiscVelSolT& wind, std::string name= std::string())
{
    return new InterfaceMatrixAccuP1CL< LocalMatrixT<DiscVelSolT> >( mat,
        LocalMatrixT<DiscVelSolT>( wind), cdata, name);
}


class LocalInterfaceMassP1CL
{
  private:
    std::valarray<double> q[4];
    QuadDomain2DCL qdom;

  public:
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
    double coup[4][4];

    void setup (const TetraCL& t, const InterfaceCommonDataP1CL& cdata) {
        resize_and_evaluate_piecewise_normal( cdata.surf, t, n, &absdet);
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
    const DiscVelSolT& w_; // wind

    QuadDomain2DCL qdom;
    Point3DCL grad[4];
    double dummy;
    LocalP2CL<Point3DCL> w_loc;
    std::valarray<double> q[4];
    GridFunctionCL<Point3DCL> qw;

  public:
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
    const DiscVelSolT& w_;

    QuadDomain2DCL qdom;
    LocalP2CL<Point3DCL> w_loc;
    std::valarray<double> q[4];
    double dummy;
    SMatrixCL<3,3> T;
    GridFunctionCL<Point3DCL> n_tri,
                              n,
                              qgradp2i;
    std::valarray<double> qdivgamma_w;
    LocalP1CL<Point3DCL> gradrefp2[10],
                         gradp2[10];


  public:
    double coup[4][4];

    void setup (const TetraCL& t, const InterfaceCommonDataP1CL& cdata);

    LocalInterfaceMassDivP1CL (const DiscVelSolT& w)
        : w_( w) { P2DiscCL::GetGradientsOnRef( gradrefp2); }
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
    double omit_bound_;

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

///\brief Quadratic in space and linear in time space-time FE-function on a single TetraPrismCL.
/// The class holds references to two LocalP2CL-objects.
template <typename T= double>
class LocalSTP2P1ProxyCL
{
  public:
    typedef T value_type;
    typedef T (*instat_fun_ptr) (const Point3DCL&, double);
    typedef LocalP2CL<T> spatial_localfe_type;

  protected:
    typedef LocalSTP2P1ProxyCL<T> self_;

  private:
    const LocalP2CL<T>& v0_;
    const LocalP2CL<T>& v1_;

  public:
    LocalSTP2P1ProxyCL(const spatial_localfe_type& v0, const spatial_localfe_type& v1)
        : v0_( v0), v1_( v1) {}

    // pointwise evaluation in STCoordCL-coordinates
    inline value_type operator() (const STCoordCL& p) const {
        return (1. - p.t_ref)*v0_( p.x_bary) + p.t_ref*v1_( p.x_bary);
    }

    // The LocalP1-function at t0
    const spatial_localfe_type& at_t0 () const { return v0_; }
    // The LocalP1-function at t1
    const spatial_localfe_type& at_t1 () const { return v1_; }
};


class STP1P1DiscCL
{
  public:
    static LocalSTP1P1CL<>                 ref_val[8];
    static LocalSTP1P1CL<Point4DCL>        ref_grad[8];
    static LocalSTP1P1CL<Point3DCL>        ref_gradx[8];
    // static LocalSTP1P1CL< SMatrixCL<4,4> > ref_Hess[8];

    static void StaticInit();
    static void StaticDestruct() {}

    static void GetGradients (const TetraPrismCL& prism, LocalSTP1P1CL<Point4DCL> grad[8]) {
        SMatrixCL<3,3> M( Uninitialized);
        double det= 0.; // dummy
        GetTrafoTr( M, det, prism.t);
        const double dt_inv= 1./(prism.t1 - prism.t0);
        for (Uint i= 0; i < 4; ++i) {
            const Point3DCL p1grad( M*FE_P1CL::DHRef( i));
            for (Uint j= 0; j < 4; ++j) {
                grad[i    ].at_t0()[j]= MakePoint4D( p1grad[0], p1grad[1], p1grad[2], j == i ? -dt_inv : 0.);
                grad[i + 4].at_t1()[j]= MakePoint4D( p1grad[0], p1grad[1], p1grad[2], j == i ?  dt_inv : 0.);
                // gradx[i    ].at_t0()[j]= p1grad;
                // gradx[i + 4].at_t1()[j]= p1grad;
            }
        }
    }
};


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


/// \brief Collect indices of unknowns
///
/// This is convenient for discretisation of operators in the Setup-routines.
class LocalNumbSTP1P1CL
{
  private:
    size_t NumIniUnknowns_,
           NumFiniUnknowns_;

  public:
    /// \brief Field of unknown-indices; NoIdx, iff the degree of freedom lies
    /// on a boundary without unknowns.
    IdxT num[8];

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
    /// \briref True, iff the shape function exists on the new interface.
    bool IsFini (IdxT i) const { return num[i] >= NumIniUnknowns_ && num[i] < NumIniUnknowns_ + NumFiniUnknowns_; }
};


/// \brief P1-discretization and solution of the transport equation on the interface
class SurfactantcGdGP1CL : public SurfactantP1BaseCL
{
  public:
    MatrixCL A,    ///< ST-diffusion matrix
             Mder, ///< ST-material-derivative matrix
             Mdiv, ///< ST-mass-matrix with interface-divergence of velocity
             Mold; ///< mass matrix on old spatial interface.

  private:
    STP1P1IdxDescCL st_idx_;
    VectorCL st_oldic_, ///< the old solution represented in the space-time-FE-basis.
             st_ic_;    ///< the new solution represented in the space-time-FE-basis.
    MatrixCL L_; ///< sum of matrices

  public:
    SurfactantcGdGP1CL (MultiGridCL& mg,
        double theta, double D, VecDescCL* v, const VelBndDataT& Bnd_v, VecDescCL& lset_vd, const BndDataCL<>& lsetbnd,
        int iter= 1000, double tol= 1e-7, double omit_bound= -1.)
    : SurfactantP1BaseCL( mg, theta, D, v, Bnd_v, lset_vd, lsetbnd, iter, tol, omit_bound)
    {}

    /// save a copy of the old level-set and velocity; moves ic to oldic; must be called before DoStep.
    virtual void InitTimeStep () {
        SurfactantP1BaseCL::InitTimeStep();
    }

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

} // end of namespace DROPS

#include "surfactant/ifacetransp.tpp"

#endif

