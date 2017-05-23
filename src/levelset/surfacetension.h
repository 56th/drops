/// \file surfacetension.h
/// \brief compute the interfacial tension
/// \author LNM RWTH Aachen: Hieu Nguyen, Yuanjun Zhang; SC RWTH Aachen:

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

#ifndef DROPS_SURFACETENSION_H
#define DROPS_SURFACETENSION_H

#include "num/discretize.h"
#include "num/bndData.h"
#include "num/lattice-eval.h"

namespace DROPS
{

enum InputMethodT
/// different types of surface forces
{
    Sigma_X=0,  ///< surface tension as a function of position
    Sigma_C=1,  ///< surface tension as a function of mass concentration
    Sigma_S=2   ///< surface tension as a function of surfactant concentration
};

enum SurfaceModelT
/// different types of constitutive laws for surfactant
{
    LANGMUIR=0, ///< Langmuir model
    LINEAR  =1  ///< simple linear model  
};

class SurfaceTensionCL
{
  private:
    BndDataCL<> cBnd_;
    double C_[5], cp_;                  ///< coefficients for computing surface tension with Sigma_C
    instat_scalar_fun_ptr sigma_;       ///< variable surface tension with Sigma_X
    InputMethodT input_;
    SurfaceModelT surfmodel_;
    VecDescCL * c_;                     ///< mass concentration
    double smax_, T_;                   ///< smax_: maximum surfactant concentration (needed in Langmuir model); T_: temperature
    VecDescCL * s_;                     ///< surfactant concentration (on the interface)

    VecDescCL * sigma_vtk_;             ///< vtk output of the variable surface tension coefficient

    double sigma_c(double c) const;     ///< variable surface tension with Sigma_C
    double dsigma_dc(double c) const;   ///< variable surface tension gradient with Sigma_C
    double sigma_s(double s) const;     ///< variable surface tension with Sigma_S

  public:
    SurfaceTensionCL ( instat_scalar_fun_ptr sigma)
        :  cBnd_(0),  cp_(0.), sigma_(sigma), input_(Sigma_X), surfmodel_(LANGMUIR), c_(0), smax_(0.), T_(0.), s_(0), sigma_vtk_(0)
    { std::memset( C_, 0, 5*sizeof( double));}

    SurfaceTensionCL ( instat_scalar_fun_ptr sigma, BndDataCL<> cBnd)
    :  cBnd_(cBnd),  cp_(0.), sigma_(sigma), input_(Sigma_X), surfmodel_(LANGMUIR), c_(0), smax_(0.), T_(0.), s_(0), sigma_vtk_(0)
    { std::memset( C_, 0, 5*sizeof( double));}

    SurfaceTensionCL ( double C[5], double cp=0., BndDataCL<> cBnd = BndDataCL<>( 0),  VecDescCL* c =0)
    :  cBnd_(cBnd),  cp_(cp), sigma_(0), input_(Sigma_C), surfmodel_(LANGMUIR), c_(c), smax_(0.), T_(0.), s_(0), sigma_vtk_(0)
    { std::memcpy( C_, C, 5*sizeof( double));}

    SurfaceTensionCL ( double smax, instat_scalar_fun_ptr sigma, VecDescCL* s =0, double T =300, BndDataCL<> cBnd = BndDataCL<>( 0))
    :  cBnd_(cBnd),  cp_(0.), sigma_(sigma), input_(Sigma_S), surfmodel_(LANGMUIR), c_(0), smax_(smax), T_(T), s_(s), sigma_vtk_(0)
    { std::memset( C_, 0, 5*sizeof( double));}

    void SetInputMethod(InputMethodT input) {input_=input;}
    InputMethodT GetInputMethod() const {return input_;}

    void SetSurfaceModel(SurfaceModelT surfmodel) {surfmodel_=surfmodel;}
    SurfaceModelT GetSurfaceModel() const {return surfmodel_;}

    void SetConcentration(VecDescCL * c) {c_=c;}
    void SetBoundary (BndDataCL<> cBnd) {cBnd_ =  cBnd;}
    void SetTime(double time) { if (c_) c_->t=time; if (s_) s_->t=time; }
    void SetVtkOutput(VecDescCL * sigma_vtk) {sigma_vtk_=sigma_vtk;}

    instat_scalar_fun_ptr GetSigma() const {return sigma_;}
    void SetCoeff (double C[5], double cp) {std::memcpy( C_, C, 5*sizeof( double)); cp_= cp;}
    void SetIFCoeff (double smax, double T) {smax_=smax; T_=T;}

    /// \brief evaluate surface tension on triangle given by the first and second arguments
    void ComputeSF( const TetraCL&, const BaryCoordCL * const, Quad5_2DCL<>&) const;

    /// \brief Evaluate the surface tension on every vertex of q.
    template <class DomainT, class ResultIterT>
      inline ResultIterT
      evaluate_on_vertexes (const TetraCL& tet, const DomainT& dom, double t, ResultIterT result_iterator) const {
        if (input_ ==  Sigma_X) {
            return DROPS::evaluate_on_vertexes( sigma_, tet, dom, t, result_iterator);
        }
        else
            throw DROPSErrCL( "SurfaceTensionCL::evaluate_on_vertexes: Input method is not yet implemented.\n");

        return result_iterator;
    }
};


}// end of name space DROPS

#endif
