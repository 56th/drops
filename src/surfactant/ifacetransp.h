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
#include "levelset/mgobserve.h"
#include "out/ensightOut.h"
#include "out/vtkOut.h"
#include "num/fe_repair.h"

#ifndef DROPS_IFACETRANSP_H
#define DROPS_IFACETRANSP_H

namespace DROPS {

/// \brief Given a P1-vecdesc on the interface x and a P1-vecdesc xext on mg, the values of x
///        are copied to xext and the remaining values of xext are set to 0.
void Extend (const MultiGridCL& mg, const VecDescCL& x, VecDescCL& xext);

/// \brief Given a P1-vecdesc on the interface x and a P1-vecdesc xext on mg, the values of xext
///        are copied to x, where x has an unknown.
void Restrict (const MultiGridCL& mg, const VecDescCL& xext, VecDescCL& x);

/// \brief The routine sets up the mass-matrix in matM on the interface defined by ls.
///        It belongs to the FE induced by standard P1-elements.
void SetupInterfaceMassP1 (const MultiGridCL& MG, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd, const double alpha= 1.);

/// \brief The routine sets up the mass-matrix in matM on the interface defined by ls.
///        The FE spaces associated with rows and columns can be mixed among P1 and P1X.
///        The jumping coefficient \a alpha is defined w.r.t. the positive(0) and negative(1) part of the domain. If omitted, \a alpha is assumed to be 1 on the whole domain.
void SetupInterfaceMassP1X (const MultiGridCL& MG, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd, const double alpha[2]= 0);

/// \brief The routine equips the matrices with the right FE spaces \a mass_idx and \a surf_idx and
///        sets up the ad/de-sorption terms for coupled mass/surfactant transport.
///        The jumping coefficients \a k_a, \a k_d are defined w.r.t. the positive(0) and negative(1) part of the domain.
void SetupInterfaceSorptionP1X (const MultiGridCL& MG, const VecDescCL& ls, const BndDataCL<>& lsetbnd,
        MatDescCL* R, MatDescCL* C, MatDescCL* R_i, MatDescCL* C_i, const IdxDescCL* mass_idx, const IdxDescCL* surf_idx, const double k_a[2], const double k_d[2]);

/// \brief The routine sets up the Laplace-Beltrami-matrix in mat on the interface defined by ls.
///        It belongs to the FE induced by standard P1-elements.
///
/// D is the diffusion-coefficient
void SetupLBP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsbnd, double D);

/// \brief The routine sets up the convection-matrix in mat on the interface defined by ls.
///        It belongs to the FE induced by standard P1-elements.
///
/// The template-parameter is only used to circumvent the exact type of the discrete
/// velocity solution in the Stokes classes.
template <class DiscVelSolT>
void SetupConvectionP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsbnd, const DiscVelSolT& v);

/// \brief Helper of SetupConvectionP1
void SetupConvectionP1OnTriangle (const BaryCoordCL triangle[3], double det,
    const LocalP1CL<> p1[4], Quad5_2DCL<> qp1[4],
    const LocalP2CL<Point3DCL>& u, Point3DCL grad[4], double coup[4][4]);

/// \brief The routine sets up the mass-matrix scaled with \f$ div_\Gamma v \f$ in mat on the interface
///        defined by ls. It belongs to the FE induced by standard P1-elements.
///
/// The template-parameter is only used to circumvent the exact type of the discrete
/// velocity solution in the Stokes classes.
template <class DiscVelSolT>
void SetupMassDivP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsbnd, const DiscVelSolT& v);

/// \brief Helper of SetupMassDivP1
void SetupMassDivP1OnTriangle (const BaryCoordCL triangle[3], double det,
    const LocalP1CL<> p1[4], Quad5_2DCL<> qp1[4],
    const LocalP2CL<Point3DCL>& u, LocalP1CL<Point3DCL> gradp2[10], const Point3DCL& n, double coup[4][4]);

/// \brief The routine sets up the mixed mass-matrix on the interface given by ls: The rows belong
///        to the new timestep, the columns to the old timestep. It belongs to the FE induced by
///        standard P1-elements.
void SetupMixedMassP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsbnd);

/// \brief The routine sets up the load-vector in v on the interface defined by ls.
///        It belongs to the FE induced by standard P1-elements.
void SetupInterfaceRhsP1 (const MultiGridCL& mg, VecDescCL* v,
    const VecDescCL& ls, const BndDataCL<>& lsbnd, instat_scalar_fun_ptr f);

/// \brief Short-hand for simple loops over the interface.
/// \param t  - Reference to a tetra
/// \param ls - Levelset-reference: Something that can be handed to InterfacePatchCL::Init as 2nd argument.
/// \param bnd  - ???
/// \param p  - The InterfacePatchCL that should be used.
/// \param n  - Name of the integer to reference the interface-triangles
#define DROPS_FOR_TETRA_INTERFACE_BEGIN( t, ls, bnd, p, n) \
    (p).Init( (t), (ls), (bnd)); \
    if ((p).Intersects()) { /*We are at the phase boundary.*/ \
        for (int ch__= 0; ch__ < 8; ++ch__) { \
            (p).ComputeForChild( ch__); \
            for (int n= 0; n < (p).GetNumTriangles(); ++n) \

#define DROPS_FOR_TETRA_INTERFACE_END }}

#define DROPS_FOR_TETRA_INTERFACE_COARSE_BEGIN( t, ls, bnd, p, n) \
    (p).Init( (t), (ls), (bnd)); \
    if ((p).IntersectsChild( 8)) { /*We are at the phase boundary.*/ \
        (p).ComputeForChild( 8); \
        for (int n= 0; n < (p).GetNumTriangles(); ++n) \

#define DROPS_FOR_TETRA_INTERFACE_COARSE_END }


/// \brief Short-hand for integral on the interface.
template <typename DiscP1FunT>
double Integral_Gamma (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const DROPS::BndDataCL<>& bnd,
    const DiscP1FunT& discsol)
{
    double d( 0.);
    const DROPS::Uint lvl = ls.GetLevel();
    DROPS::InterfaceTriangleCL triangle;
    DROPS::Quad5_2DCL<> qdiscsol;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        DROPS_FOR_TETRA_INTERFACE_BEGIN( *it, ls, bnd, triangle, tri) {
            qdiscsol.assign(  *it, &triangle.GetBary( tri), discsol);
            d+= qdiscsol.quad( triangle.GetAbsDet( tri));
        }
        DROPS_FOR_TETRA_INTERFACE_END
    }
    return d;
}

/// \brief Short-hand for integral on the interface, computed on the coarse interface.
template <typename DiscP1FunT>
double Integral_Gamma_Coarse (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const DROPS::BndDataCL<>& bnd,
    const DiscP1FunT& discsol)
{
    double d( 0.);
    const DROPS::Uint lvl = ls.GetLevel();
    DROPS::InterfaceTriangleCL triangle;
    DROPS::Quad5_2DCL<> qdiscsol;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        DROPS_FOR_TETRA_INTERFACE_COARSE_BEGIN( *it, ls, bnd, triangle, tri) {
            qdiscsol.assign(  *it, &triangle.GetBary( tri), discsol);
            d+= qdiscsol.quad( triangle.GetAbsDet( tri));
        }
        DROPS_FOR_TETRA_INTERFACE_COARSE_END
    }
    return d;
}

/// \brief Short-hand for integral on the interface, h^3-version through interpolation.
template <typename DiscP1FunT>
double Integral_Gamma_Extrapolate (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const DROPS::BndDataCL<>& bnd,
    const DiscP1FunT& discsol)
{
    return (4.*Integral_Gamma( mg, ls, bnd,discsol) - Integral_Gamma_Coarse( mg, ls, bnd,discsol))/3.;
}


/// \brief P1-discretization and solution of the transport equation on the interface
class SurfactantcGP1CL
{
  public:
    typedef BndDataCL<Point3DCL>                              VelBndDataT;
    typedef NoBndDataCL<>                                     BndDataT;
    typedef P1EvalCL<double, const BndDataT, VecDescCL>       DiscSolCL;
    typedef P1EvalCL<double, const BndDataT, const VecDescCL> const_DiscSolCL;

    IdxDescCL idx;
    VecDescCL ic; ///< concentration on the interface
    MatDescCL A,  ///< diffusion matrix
              M,  ///< mass matrix
              C,  ///< convection matrix
              Md, ///< mass matrix with interface-divergence of velocity
              M2; ///< mass matrix: new trial- and test- functions on old interface

  private:
    MatrixCL      L_;              ///< sum of matrices
    MultiGridCL&  MG_;
    double        D_,              ///< diffusion coefficient
                  theta_, dt_;     ///< time scheme parameter and time step

    BndDataT            Bnd_;
    const VelBndDataT&  Bnd_v_;  ///< Boundary condition for the velocity
    VecDescCL*          v_;      ///< velocity at current time step
    VecDescCL&          lset_vd_;///< levelset at current time step

    const BndDataCL<>&  lsetbnd_; ///< level set boundary

    IdxDescCL           oldidx_; ///< idx that corresponds to old time (and oldls_)
    VectorCL            oldic_;  ///< interface concentration at old time
    VecDescCL           oldls_;  ///< levelset at old time
    VecDescCL           oldv_;   ///< velocity at old time
    double              oldt_;   ///< old time

    GSPcCL                  pc_;
    GMResSolverCL<GSPcCL>   gm_;
    double omit_bound_; ///< not used atm

  public:
    SurfactantcGP1CL (MultiGridCL& mg, const VelBndDataT& Bnd_v,
        double theta, double D, VecDescCL* v, VecDescCL& lset_vd, const BndDataCL<>& lsetbnd,
        double dt, int iter= 1000, double tol= 1e-7, double omit_bound= -1.)
    : idx( P1IF_FE), MG_( mg), D_( D), theta_( theta), dt_( dt),
        Bnd_v_( Bnd_v), v_( v), lset_vd_( lset_vd), lsetbnd_( lsetbnd), oldidx_( P1IF_FE), gm_( pc_, 100, iter, tol, true),
        omit_bound_( omit_bound)
    { idx.GetXidx().SetBound( omit_bound); }

    const MultiGridCL& GetMG() const { return MG_; }
    GMResSolverCL<GSPcCL>& GetSolver() { return gm_; }

     /// initialize the interface concentration
    void Init (instat_scalar_fun_ptr);

    /// \remarks call SetupSystem \em before calling SetTimeStep!
    void SetTimeStep( double dt, double theta=-1);

    /// save a copy of the old level-set; must be called before DoStep.
    void InitOld ();

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
    VectorCL InitStep ();
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

} // end of namespace DROPS

#include "surfactant/ifacetransp.tpp"

#endif

