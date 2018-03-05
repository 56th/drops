/// \file auto_diff.h
/// \brief Simple wrapper for arithmetic types to perform automatic differentiation.
/// \Todo Provide more functions, add class for reverse mode, etc.
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

#ifndef DROPS_AUTO_DIFF_H
#define DROPS_AUTO_DIFF_H

#include <cmath>

#include "../misc/container.h"

namespace DROPS
{

namespace AutoDiff {

///\brief Forward-mode automatic differentiation with D independent variables.
/// Providing one argument to the constructor creates a constant.
/// Calling seed(x, i) turns an object into an instance of the ith independent variable with value x.
template <typename T, int D>
class ADFwdCL
{
  public:
    enum class FreeVariableE { free_variable };

    using value_type= T;
    using derivative_type= SVectorCL<D>;

  private:
    value_type       s;
     derivative_type ds;

  public:
    ADFwdCL ()
        : ADFwdCL (value_type(), derivative_type()) {}
    ADFwdCL (value_type sarg)
        : ADFwdCL (sarg, derivative_type()) {}
    ADFwdCL (value_type sarg, const derivative_type& dsarg)
        : s( sarg), ds( dsarg) {}

    ADFwdCL (value_type sarg, int i, FreeVariableE)
        : s (sarg), ds (derivative_type()) {
        ds[i]= 1;
    }

    ADFwdCL (const ADFwdCL& ad)= default;
    ADFwdCL (ADFwdCL&& ad)= default;
    ADFwdCL& operator= (const ADFwdCL& ad)= default;

    const value_type& value () const { return s; }
    value_type&       value ()       { return s; }
    const derivative_type& derivative () const { return ds; }
    derivative_type&       derivative ()       { return ds; }

    void seed (T sarg, int i) {
        s= sarg;
        std::fill_n (ds.begin(), D, value_type());
        ds[i]= 1;
    }
};

template <typename T, int D>
ADFwdCL<T, D> operator+ (ADFwdCL<T, D> f, ADFwdCL<T, D> g)
{
    return ADFwdCL<T, D>( f.value()      + g.value(),
                          f.derivative() + g.derivative());
}

template <typename T, int D>
ADFwdCL<T, D> operator+ (double f, ADFwdCL<T, D> g)
{
    return ADFwdCL<T, D>( f + g.value(),
                          g.derivative());
}

template <typename T, int D>
ADFwdCL<T, D> operator+ (ADFwdCL<T, D> f, double g)
{
    return ADFwdCL<T, D>( f.value() + g,
                          f.derivative());
}

template <typename T, int D>
ADFwdCL<T, D> operator- (ADFwdCL<T, D> f)
{
    return ADFwdCL<T, D>( -f.value(),
                          -f.derivative());
}

template <typename T, int D>
ADFwdCL<T, D> operator- (ADFwdCL<T, D> f, ADFwdCL<T, D> g)
{
    return ADFwdCL<T, D>( f.value()      - g.value(),
                          f.derivative() - g.derivative());
}

template <typename T, int D>
ADFwdCL<T, D> operator- (double f, ADFwdCL<T, D> g)
{
    return ADFwdCL<T, D>( f - g.value(),
                          -g.derivative());
}

template <typename T, int D>
ADFwdCL<T, D> operator- (ADFwdCL<T, D> f, double g)
{
    return ADFwdCL<T, D>( f.value() - g,
                          f.derivative());
}

template <typename T, int D>
ADFwdCL<T, D> operator* (ADFwdCL<T, D> f, ADFwdCL<T, D> g)
{
    return ADFwdCL<T, D>( f.value()*g.value(),
                          f.derivative()*g.value() + f.value()*g.derivative());
}

template <typename T, int D>
ADFwdCL<T, D> operator* (double f, ADFwdCL<T, D> g)
{
    return ADFwdCL<T, D>( f*g.value(),
                          f*g.derivative());
}

template <typename T, int D>
ADFwdCL<T, D> operator* (ADFwdCL<T, D> f, double g)
{
    return ADFwdCL<T, D>( f.value()*g,
                          f.derivative()*g);
}

template <typename T, int D>
ADFwdCL<T, D> operator/ (ADFwdCL<T, D> f, ADFwdCL<T, D> g)
{
    return ADFwdCL<T, D>( f.value()/g.value(),
                          (f.derivative()*g.value() - f.value()*g.derivative())/std::pow(g.value(), 2));
}

template <typename T, int D>
ADFwdCL<T, D> operator/ (double f, ADFwdCL<T, D> g)
{
    return ADFwdCL<T, D>( f/g.value(),
                          (-f)*g.derivative()/std::pow(g.value(), 2));
}

template <typename T, int D>
ADFwdCL<T, D> operator/ (ADFwdCL<T, D> f, double g)
{
    return ADFwdCL<T, D>( f.value()/g,
                          f.derivative()/g);
}

template <typename T, int D>
ADFwdCL<T, D> pow (ADFwdCL<T, D> f, double i)
{
    return ADFwdCL<T, D>( std::pow( f.value(), i),
                          i*std::pow( f.value(), i - 1)*f.derivative());
}
template <typename T, int D>
ADFwdCL<T, D> sqrt (ADFwdCL<T, D> f)
{
    return ADFwdCL<T, D>( std::pow( f.value(), 0.5),
                          0.5*std::pow( f.value(), -0.5)*f.derivative());
}
template <typename T, int D>
ADFwdCL<T, D> cos (ADFwdCL<T, D> f)
{
    return ADFwdCL<T, D>( std::cos( f.value()),
                         -std::sin( f.value())*f.derivative());
}

template <typename T, int D>
ADFwdCL<T, D> sin (ADFwdCL<T, D> f)
{
    return ADFwdCL<T, D>( std::sin( f.value()),
                          std::cos( f.value())*f.derivative());
}
template <typename T, int D>
ADFwdCL<T, D> exp (ADFwdCL<T, D> f)
{
    return ADFwdCL<T, D>( std::exp( f.value()),
                          std::exp( f.value())*f.derivative());
}




} // end of namespace DROPS::AutoDiff
} // end of namespace DROPS

#endif
