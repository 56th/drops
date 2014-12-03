/// \file omp_variable.h
/// \brief OpenMP-safe accumulation of scalars, vectors...
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

#ifndef DROPS_OMP_VARABLE_H
#define DROPS_OMP_VARABLE_H

#include <algorithm>

namespace DROPS {


///\brief Variable without update race in AccumulatorCL.
/// Usage: Define an instance f of a subclass of OpenMPVarBaseCL in the AccumulatorCL,
///        call scatter() in begin_accumulation,
///        call reduce() in finalize_accumulation, and
///        call make_reference_to( f) in clone() on the copy of f in the clone.
///        In visit(), refer to f.value( tid) in thread tid for the local version.
///        Outside the OpenMP-parallel region, refer to f.value() for the accumulated value.
/// Custom scatter and reduce can be provided by subclassing and writing
/// scatter_impl and reduce_impl (CRTP method).
template <typename T, typename Derived>
class OpenMPVarBaseCL
{
  public:
    typedef T value_type;

  protected:
    size_t n;   ///< maximum number of threads.
    bool owner; ///< If true, the object must delete xx.
    T x;        ///< reduced value (before scatter and after reduce)
    T* xx;      ///< scattered values for each thread-id (after scatter and before reduce).

  public:
    ///\brief Swap also swaps the ownership of xx.
    void swap (OpenMPVarBaseCL& r) {
        std::swap( n, r.n);
        std::swap( owner, r.owner);
        std::swap( x, r.x);
        std::swap( xx, r.xx);
    }

    ///\brief Constructs an object which owns xx.
    OpenMPVarBaseCL()
        : n( 0), owner( true), x(), xx( 0) {}
    ~OpenMPVarBaseCL() {
        if (owner && xx)
            delete[] xx;
    }
    ///\brief Deep copy.
    OpenMPVarBaseCL (const OpenMPVarBaseCL& r)
        : n( r.n), owner( r.owner), x( r.x), xx( !owner ? r.xx : (r.xx ? new T[n] : 0)) {
        if (owner && r.xx)
            std::copy( r.xx, r.xx + r.n, xx);
    }

    ///\brief Deep copy-assignment.
    /// Note that the argument is passed by value (i.e. copy constructed from the original).
    OpenMPVarBaseCL& operator= (OpenMPVarBaseCL r) {
        this->swap( r);
        return *this;
    }

    ///\brief Let xx point to r.xx;
    /// call on the clones after scatter has been called on the original (this fits to the model of AccumulatorCL).
    void make_reference_to (OpenMPVarBaseCL& r) {
        if (owner && xx) {
            delete[] xx;
            owner= false;
        }
        xx= r.xx;
        x= r.x;
    }

    const T& value () const {
        return x;
    }
    T& value () {
        return x;
    };
#ifdef _OPENMP
    const T& value (int tid) const {
        return xx[tid];
    };
    T& value (int tid) {
        return xx[tid];
    };
#else
    const T& value (int) const {
        return x;
    };
    T& value (int) {
        return x;
    };
#endif

    void scatter () {
#ifdef _OPENMP
        if (owner) {
            n= omp_get_max_threads();
            xx= new T[n];
        }
#endif
        static_cast<Derived*>( this)->scatter_impl();
    }
    void reduce  () {
        static_cast<Derived*>( this)->reduce_impl();
#ifdef _OPENMP
        if (owner) {
            delete[] xx;
            xx= 0;
            n= 0;
        }
#endif
    }
};

///\brief Scatter zero-initializes xx, reduce puts the sum of xx into x.
template <typename T>
class OpenMPVar_ZeroInit_Sum_CL : public OpenMPVarBaseCL<T, OpenMPVar_ZeroInit_Sum_CL<T> >
{
  public:
    typedef OpenMPVarBaseCL<T, OpenMPVar_ZeroInit_Sum_CL<T> > base_type;

    void scatter_impl () {
#ifdef _OPENMP
        std::fill( base_type::xx, base_type::xx + base_type::n, T());
#else
         base_type::x= T();
#endif
    }
    void reduce_impl () {
#ifdef _OPENMP
        base_type::x= std::accumulate( base_type::xx, base_type::xx + base_type::n, T());
#endif
    }

    using typename base_type::value_type;
    using base_type::swap;
    using base_type::make_reference_to;
    using base_type::value;
    using base_type::scatter;
    using base_type::reduce;
};

///\brief Scatter initializes xx with -std::numeric_limits<T>::max(), reduce puts the max of xx into x.
template <typename T>
class OpenMPVar_MinInit_Max_CL : public OpenMPVarBaseCL<T, OpenMPVar_MinInit_Max_CL<T> >
{
  public:
    typedef OpenMPVarBaseCL<T, OpenMPVar_MinInit_Max_CL<T> > base_type;

    void scatter_impl () {
#ifdef _OPENMP
        std::fill( base_type::xx, base_type::xx + base_type::n, -std::numeric_limits<T>::max());
#else
         base_type::x= -std::numeric_limits<T>::max();
#endif
    }
    void reduce_impl () {
#ifdef _OPENMP
        base_type::x= *std::max_element( base_type::xx, base_type::xx + base_type::n);
#endif
    }

    using typename base_type::value_type;
    using base_type::swap;
    using base_type::make_reference_to;
    using base_type::value;
    using base_type::scatter;
    using base_type::reduce;
};

///\brief Scatter initializes xx with std::numeric_limits<T>::max(), reduce puts the min of xx into x.
template <typename T>
class OpenMPVar_MaxInit_Min_CL : public OpenMPVarBaseCL<T, OpenMPVar_MaxInit_Min_CL<T> >
{
  public:
    typedef OpenMPVarBaseCL<T, OpenMPVar_MaxInit_Min_CL<T> > base_type;

    void scatter_impl () {
#ifdef _OPENMP
        std::fill( base_type::xx, base_type::xx + base_type::n, std::numeric_limits<T>::max());
#else
         base_type::x= std::numeric_limits<T>::max();
#endif
    }
    void reduce_impl () {
#ifdef _OPENMP
        base_type::x= *std::min_element( base_type::xx, base_type::xx + base_type::n);
#endif
    }

    using typename base_type::value_type;
    using base_type::swap;
    using base_type::make_reference_to;
    using base_type::value;
    using base_type::scatter;
    using base_type::reduce;
};

} // end of namespace DROPS

#endif
