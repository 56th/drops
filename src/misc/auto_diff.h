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

namespace DROPS
{

///\brief Forward-mode automatic differentiation.
/// Providing one argument to the constructor creates a constant.
/// Calling seed() turns an object into an instance of the independent variable.
/// The parameter of seed (and the second constructor parameter) is useful for directional derivatives other than the partial derivatives.
template <typename T>
class ADFwdCL
{
  private:
    T s,
      ds;

  public:
    ADFwdCL (T sarg= T(), T dsarg= T())
        : s( sarg), ds( dsarg) {}

    T value () const { return s; }
    T derivative () const { return ds; }

    void seed (T dsarg= T( 1)) { ds= dsarg; }
};

template <typename T>
ADFwdCL<T> operator+ (ADFwdCL<T> f, ADFwdCL<T> g)
{
    return ADFwdCL<T>( f.value()      + g.value(),
                       f.derivative() + g.derivative());
}

template <typename T>
ADFwdCL<T> operator- (ADFwdCL<T> f)
{
    return ADFwdCL<T>( -f.value(),
                       -f.derivative());
}

template <typename T>
ADFwdCL<T> operator- (ADFwdCL<T> f, ADFwdCL<T> g)
{
    return ADFwdCL<T>( f.value()      - g.value(),
                       f.derivative() - g.derivative());
}

template <typename T>
ADFwdCL<T> operator* (ADFwdCL<T> f, ADFwdCL<T> g)
{
    return ADFwdCL<T>( f.value()*g.value(),
                       f.derivative()*g.value() + f.value()*g.derivative());
}

template <typename T>
ADFwdCL<T> operator/ (ADFwdCL<T> f, ADFwdCL<T> g)
{
    return ADFwdCL<T>( f.value()/g.value(),
                       (f.derivative()*g.value() - f.value()*g.derivative())/std::pow(g.value(), 2));
}

template <typename T>
ADFwdCL<T> pow (ADFwdCL<T> f, int i)
{
    return ADFwdCL<T>( std::pow( f.value(), i),
                       i*std::pow( f.value(), i - 1));
}

template <typename T>
ADFwdCL<T> cos (ADFwdCL<T> f)
{
    return ADFwdCL<T>( std::cos( f.value()),
                       -std::sin( f.value())*f.derivative());
}

template <typename T>
ADFwdCL<T> sin (ADFwdCL<T> f)
{
    return ADFwdCL<T>( std::sin( f.value()),
                       std::cos( f.value())*f.derivative());
}

} // end of namespace DROPS

#endif
