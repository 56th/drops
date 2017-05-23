/// \file poisson.h
/// \brief classes that constitute the poisson-problem
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Marcus Soemers, Volker Reichelt, Liang Zhang; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_POISSON_H
#define DROPS_POISSON_H

#include "misc/problem.h"
#include "misc/params.h"
#include "num/fe.h"
#include "num/discretize.h"
#include "num/bndData.h"
#include <deque>
#include <iostream>
#include <numeric>

#include "poisson/ale.h"
#include "geom/deformation.h"

namespace DROPS
{

typedef BndSegDataCL<> PoissonBndSegDataCL;
typedef BndDataCL<> PoissonBndDataCL;

class StripTimeCL
// converts time dependent function to one, that is time independent.
// i.e.:  scalar_instat_fun_ptr  -->  scalar_fun_ptr
// useful where one wants to use all the stuff written for the stationary problems
{
  private:
    static instat_scalar_fun_ptr func_;
    static double                t_;
  public:
    StripTimeCL( instat_scalar_fun_ptr func, double t)
      { func_= func; t_= t; }
    void SetTime( double t) { t_= t; }

    static double GetFunc( const Point3DCL& x)
      { return func_( x, t_); }
};

//Streamline diffusion stabilization class which can compute difference stabilization coefficient according to the grids type.
class SUPGCL
{
  private:
    double magnitude_;
    int    grids_;      // decide how to compute characteristic length to approximate the longest length in flow direction
    double longedge_;   // the longest edge for regular grids;
    bool   SUPG_;

  public:
    SUPGCL()
        : magnitude_(0), grids_(1), longedge_(0), SUPG_(false) {}
    SUPGCL(ParamCL para) { init(para); }
    void init(ParamCL para)
    {
        const double lx_= norm(para.get<Point3DCL>("Mesh.E1")),
                     ly_= norm(para.get<Point3DCL>("Mesh.E2")),
                     lz_= norm(para.get<Point3DCL>("Mesh.E3"));
        const int    nx_= para.get<int>("Mesh.N1"),
                     ny_= para.get<int>("Mesh.N2"),
                     nz_= para.get<int>("Mesh.N3");
        int Ref_=para.get<int>("Mesh.AdaptRef.FinestLevel");
        magnitude_ =para.get<double>("Stabilization.Magnitude");
        grids_     =para.get<double>("Stabilization.Grids");
        //pick up the longest edge
        if(grids_==1)
        {
            double dx_= lx_/(nx_*std::pow(2, Ref_));
            double dy_= ly_/(ny_*std::pow(2, Ref_));
            double dz_= lz_/(nz_*std::pow(2, Ref_));
            double m;
            if(dx_>=dy_)
                m = dx_;
            else
                m = dy_;
            if( m>=dz_)
                longedge_=m;
            else
                longedge_=dz_;
        }
        else
            longedge_ = 0;
        SUPG_ = para.get<int>("Stabilization.SUPG");
    }

    bool GetSUPG() const { return SUPG_; }

    double GetCharaLength(int grids) const
    {
        double h=0.;
        if(grids==1)
            h=longedge_;
        else
            std::cout<<"WARNING: The geometry type has not been implemented!\n";
        return h;
    }

    double Sta_Coeff(const DROPS::Point3DCL& Vel, double alpha) const
    { //Stabilization coefficient
        double Pec=0.;
        double h  =GetCharaLength(grids_);
        Pec=Vel.norm()*h/(2.*alpha);  //compute mesh Peclet number
        if (Pec<=1)
            return 0.0;
        else
            return magnitude_*h/(2.*Vel.norm())*(1.-1./Pec);
    }
};

/*class PoissonDeformationHelperCL
{
    private:
      bool ALE_;
      MultiGridCL& mg_;
    public:
      PoissonDeformationHelperCL(MultiGridCL& mg, bool ALE): mg_(mg), ALE_(ALE) {}
      Point3DCL GetVertexCoord(const VertexCL& vert) 
      {   
          if(ALE_)
            return mg_.GetMeshDeformation().GetTransformedVertexCoord(vert);
          else
            return vert.GetCoord();
      }
      Point3DCL GetEdgeBaryCenter(const EdgeCL& edge) 
      {   
          if(ALE_)
            return mg_.GetMeshDeformation().GetTransformedEdgeBaryCenter(edge);
          else
            return GetBaryCenter(edge);
      }    
}*/

template <class Coeff>
class PoissonP1CL : public ProblemCL<Coeff, PoissonBndDataCL>
{
  private:
    bool    adjoint_;
    SUPGCL& supg_;
    bool    ALE_;           // ALE method
    MeshDeformationCL& md_;

  public:
    typedef ProblemCL<Coeff, PoissonBndDataCL> base_;
    typedef typename base_::BndDataCL          BndDataCL;
    typedef typename base_::CoeffCL            CoeffCL;
    using                                      base_::MG_;
    using                                      base_::Coeff_;
    using                                      base_::BndData_;
    using                                      base_::GetBndData;
    using                                      base_::GetMG;

    typedef P1EvalCL<double, const BndDataCL, VecDescCL>       DiscSolCL;
    typedef P1EvalCL<double, const BndDataCL, const VecDescCL> const_DiscSolCL;
    typedef double (*est_fun)(const TetraCL&, const VecDescCL&, const BndDataCL&);

    MLIdxDescCL idx;
    VecDescCL   x;
    VecDescCL   b;
    VecDescCL   vU;
    MLMatDescCL A;
    MLMatDescCL M;
    MLMatDescCL U;


    PoissonP1CL(const MGBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata, SUPGCL& supg, bool ALE, bool adj=false)
        : base_( mgb, coeff, bdata), adjoint_( adj), supg_(supg), ALE_(ALE), md_(MeshDeformationCL::getInstance()), idx( P1_FE) {}

    PoissonP1CL(MultiGridCL& mg, const CoeffCL& coeff, const BndDataCL& bdata, SUPGCL& supg, bool ALE, bool adj=false)
        : base_( mg, coeff, bdata), adjoint_( adj), supg_(supg), ALE_(ALE), md_(MeshDeformationCL::getInstance()), idx( P1_FE) {}
    // numbering of unknowns
    void CreateNumbering( Uint level, MLIdxDescCL* idx)
        { idx->CreateNumbering( level, MG_, BndData_); }
    void DeleteNumbering( MLIdxDescCL* idx)
        { idx->DeleteNumbering( MG_); }
    void SetNumLvl( size_t n);

    bool usesALE() const { return ALE_; }

    // set up matrices and rhs
    void SetupSystem         (MLMatDescCL&, VecDescCL&) const;
    ///  \brief set up matrices (M is time independent)
    void SetupInstatSystem( MLMatDescCL& A, MLMatDescCL& M, double t) const;
    /// \brief set up matrix and couplings with bnd unknowns for convection term
    void SetupConvection( MLMatDescCL& U, VecDescCL& vU, double t) const;

    /// \brief Setup time dependent parts
    ///
    /// couplings with bnd unknowns, coefficient f(t)
    /// If the function is called with the same vector for some arguments,
    /// the vector will contain the sum of the results after the call
    void SetupInstatRhs( VecDescCL& vA, VecDescCL& vM, double t, VecDescCL& vf, double tSUPG= 0.) const;
    /// \brief Setup special source term including the gradient of a given P1 function
    void SetupGradSrc( VecDescCL& src, instat_scalar_fun_ptr T, instat_scalar_fun_ptr dalpha, double t= 0.) const;

    /// \brief Set initial value
    void Init( VecDescCL&, instat_scalar_fun_ptr, double t0= 0.) const;

    /// \brief check computed solution etc.
    double CheckSolution( const VecDescCL&, instat_scalar_fun_ptr, double t=0.) const;
    double CheckSolution( instat_scalar_fun_ptr Lsg, double t=0.) const { return CheckSolution(x, Lsg, t); }
    double CheckSolution( const VecDescCL&, scalar_tetra_function, double t=0.) const;
    
    void GetDiscError   ( const MLMatDescCL&, instat_scalar_fun_ptr, double =0.) const;
    void GetDiscError   ( instat_scalar_fun_ptr Lsg, double t=0.) const { GetDiscError(A, Lsg, t); }

    bool          EstimateError         ( const VecDescCL&, const double, double&, est_fun);
    static double ResidualErrEstimator  ( const TetraCL&, const VecDescCL&, const BndDataCL&);
    static double ResidualErrEstimatorL2( const TetraCL&, const VecDescCL&, const BndDataCL&);

    DiscSolCL GetSolution()
        { return DiscSolCL(&x, &GetBndData(), &GetMG()); }
    const_DiscSolCL GetSolution() const
        { return const_DiscSolCL(&x, &GetBndData(), &GetMG()); }
};

template <class Coeff>
class PoissonP2CL : public ProblemCL<Coeff, PoissonBndDataCL>
{
  private:
    bool ALE_;
    MeshDeformationCL& md_;

  public:
    typedef ProblemCL<Coeff, PoissonBndDataCL> base_;
    typedef typename base_::BndDataCL               BndDataCL;
    typedef typename base_::CoeffCL                 CoeffCL;
    using                                           base_::MG_;
    using                                           base_::Coeff_;
    using                                           base_::BndData_;
    using                                           base_::GetBndData;
    using                                           base_::GetMG;

    typedef P2EvalCL<double, const BndDataCL, VecDescCL>       DiscSolCL;
    typedef P2EvalCL<double, const BndDataCL, const VecDescCL> const_DiscSolCL;
    typedef double (*est_fun)(const TetraCL&, const VecDescCL&, const BndDataCL&);

    // new fields for the matrix A, the rhs b and the solution x
    MLIdxDescCL idx;
    VecDescCL   x;
    VecDescCL   b;
    MLMatDescCL A;

    MLMatDescCL U;   //Convection matrix
    MLMatDescCL M;   //Mass matrix
    VecDescCL   vU;  //Coupling with convection matrix

    //create an element of the class
    PoissonP2CL(const MGBuilderCL& mgb, const CoeffCL& coeff,
                const BndDataCL& bdata, bool ALE = false) : base_(mgb, coeff, bdata), ALE_(ALE), md_(MeshDeformationCL::getInstance()), idx(P2_FE) {}
    PoissonP2CL(MultiGridCL& mg, const CoeffCL& coeff, const BndDataCL& bdata, bool ALE = false)
        : base_( mg, coeff, bdata), ALE_(ALE), md_(MeshDeformationCL::getInstance()), idx( P2_FE) {}
    // numbering of unknowns
    void CreateNumbering( Uint level, MLIdxDescCL* idx)
        { idx->CreateNumbering( level, MG_, BndData_); }
    void DeleteNumbering( MLIdxDescCL* idx)
        { idx->DeleteNumbering( MG_); }
    void SetNumLvl( size_t n);

    bool usesALE() const { return ALE_; }

    // set up matrices and rhs
    void SetupSystem         ( MLMatDescCL&, VecDescCL&) const;

    ///  \brief set up matrices for instatProblem
    void SetupInstatSystem( MLMatDescCL& A, MLMatDescCL& M, double t) const;

    void SetupInstatRhs( VecDescCL& vA, VecDescCL& vM, double tA, VecDescCL& vf, double tf) const;

    //Set up convection
    void SetupConvection( MLMatDescCL&, VecDescCL&, double) const;

    //Set up initial value
    void Init( VecDescCL&, instat_scalar_fun_ptr, double t0= 0.) const;

    // check computed solution, etc.
    double CheckSolution( const VecDescCL&, instat_scalar_fun_ptr, double t=0.) const;
    double CheckSolution( instat_scalar_fun_ptr Lsg, double t=0.) const { return CheckSolution(x, Lsg, t); }
    double CheckSolution( const VecDescCL&, scalar_tetra_function, double t=0.) const;

    DiscSolCL GetSolution()
        { return DiscSolCL(&x, &GetBndData(), &GetMG()); }
    const_DiscSolCL GetSolution() const
        { return const_DiscSolCL(&x, &GetBndData(), &GetMG()); }
};


double SimpleGradEstimator ( const TetraCL& t, const VecDescCL& lsg, const PoissonBndDataCL&);

//==============================================================
//                  Marker classes
//==============================================================

template <class TetraEstT, class ProblemT>
class PoissonErrEstCL
{
  private:
    double    InitGlobErr_;
    double    RelReduction_;
    double    ConvExponent_;
    double    Meas_;
    double    ActGlobErr_;
    TetraEstT Estimator_;
    ProblemT& Problem_;
    Uint      NumLastMarkedForRef_;
    Uint      NumLastMarkedForDel_;
    bool      DoMark_;
    std::ostream*  outp_;

  public:
      PoissonErrEstCL(double RelReduction, double ConvExp, double Meas, bool DoMark, TetraEstT est,
                      ProblemT& problem, std::ostream* osp= &std::cout)
        : InitGlobErr_(0), RelReduction_(RelReduction), ConvExponent_(ConvExp), Meas_(Meas), ActGlobErr_(-1), Estimator_(est),
          Problem_(problem), NumLastMarkedForRef_(0), NumLastMarkedForDel_(0), DoMark_(DoMark), outp_(osp)
        {
#ifdef _PAR
            throw DROPSErrCL("This class has not been parallelized yet, sorry");
#endif
        }
    // default assignment-op, copy-ctor, dtor

    template <class BndData_, class VD_>
    void Init(const P1EvalCL<double, BndData_, VD_>&);

    double GetRelRed() const { return RelReduction_; }
    void SetRelRed(double newred) { RelReduction_= newred; }
    bool DoesMark() const { return DoMark_; }
    void SwitchMark() { DoMark_= DoMark_ ? false : true; }
    template <class BndData_, class VD_>
    bool Estimate(const P1EvalCL<double, BndData_, const VD_>&);
};


template <class TetraEstT, class ProblemT>
class DoerflerMarkCL
{
  private:
    double        InitGlobErr_;
    double        RelReduction_;
    double        min_tetra_ratio_;
    double        Threshold_;
    double        Meas_;
    double        ActGlobErr_;
    TetraEstT     Estimator_;
    ProblemT&     Problem_;
    Uint          NumLastMarkedForRef_;
    Uint          NumLastMarkedForDel_;
    bool          DoMark_;
    std::ostream* outp_;

  public:
  // the tetras are sorted: T_1 with biggest error, last T_n with smallest
  // a tetra T_i is marked for refinement, iff (a) or (b) holds:
  // (a) it is among "min_ratio" % of the tetras with biggest errors,
  // (b) the sum err_1+..+err_i accounts for less than "Threshold" % of the global error

      DoerflerMarkCL(double RelReduction, double min_ratio, double Threshold, double Meas, bool DoMark, TetraEstT est,
                     ProblemT& problem, std::ostream* osp= &std::cout)
        : InitGlobErr_(0), RelReduction_(RelReduction), min_tetra_ratio_(min_ratio), Threshold_(Threshold), Meas_(Meas), ActGlobErr_(-1), Estimator_(est),
          Problem_(problem), NumLastMarkedForRef_(0), NumLastMarkedForDel_(0), DoMark_(DoMark), outp_(osp)
        {
#ifdef _PAR
            throw DROPSErrCL("This class has not been parallelized yet, sorry");
#endif
        }
    // default assignment-op, copy-ctor, dtor

    template <class BndData_, class VD_>
    void Init(const P1EvalCL<double, BndData_, VD_>&);

    double GetRelRed() const { return RelReduction_; }
    void   SetRelRed(double newred) { RelReduction_= newred; }
    double GetThreshold() const { return Threshold_; }
    void   SetThreshold(double newThreshold) { Threshold_= newThreshold; }
    bool   DoesMark() const { return DoMark_; }
    void   SwitchMark() { DoMark_= DoMark_ ? false : true; }
    template <class BndData_, class VD_>
    bool Estimate(const P1EvalCL<double, BndData_, const VD_>&);
};


} // end of namespace DROPS

#include "poisson/poisson.tpp"

#endif
