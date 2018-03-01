/// \file container.h
/// \brief simple array in STL style, level-based array
/// \author LNM RWTH Aachen: Joerg Grande, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_CONTAINER_H
#define DROPS_CONTAINER_H

#include <vector>
#include <list>
#include <cmath>
#include <iostream>
#include <valarray>
#include "misc/utils.h"

namespace DROPS
{

//**************************************************************************
// Class:    DMatrixCL                                                     *
// Purpose:  dynamical storage for 2D data                                 *
// Remarks:  just the bare minimum                                         *
//**************************************************************************

template <typename T>
class DMatrixCL
{
  private:
    size_t Rows_, Cols_;
    T* Array_;

  public:
    DMatrixCL(size_t row, size_t col) : Rows_(row), Cols_(col), Array_(new T[row*col]) {}
    ~DMatrixCL() { delete[] Array_; }

    T& operator() (size_t row, size_t col)       {
        Assert(row<Rows_ && col<Cols_, DROPSErrCL("DMatrixCL::operator(): Invalide index"), DebugNumericC);
        return Array_[col*Rows_+row];
    }
    T  operator() (size_t row, size_t col) const {
        Assert(row<Rows_ && col<Cols_, DROPSErrCL("DMatrixCL::operator() const: Invalide index"), DebugNumericC);
        return Array_[col*Rows_+row];
    }
    T* GetCol       (size_t col)                 {
        Assert(col<Cols_, DROPSErrCL("DMatrixCL::GetCol: Invalide index"), DebugNumericC);
        return Array_+col*Rows_;
    }
    const T* GetCol (size_t col) const           {
        Assert(col<Cols_, DROPSErrCL("DMatrixCL::GetCol: Invalide index"), DebugNumericC);
        return Array_+col*Rows_;
    }
    Uint num_rows() const { return Rows_; }  // Zeilenzahl
    Uint num_cols() const { return Cols_; }  // Spaltenzahl
};

inline double
frobenius_norm_sq (const DMatrixCL<double>& a)
{
    double ret = 0;
    const double* d= a.GetCol( 0);
    for (Uint i= 0; i < a.num_rows()*a.num_cols(); ++i)
        ret += d[i]*d[i];
    return ret;
}


template <class T, Uint _Size>
  class SArrayCL;
template <class T, Uint _Size>
  class SBufferCL;
template <Uint _Size>
  class SVectorCL;
template <Uint _Rows, Uint _Cols>
  class SMatrixCL;
template <class T, Uint _Size>
  inline bool
  operator==(const SArrayCL<T, _Size>&, const SArrayCL<T, _Size>&);
template <class T, Uint _Size>
  inline bool
  operator<(const SArrayCL<T, _Size>&, const SArrayCL<T, _Size>&);
template <class T, Uint _Size>
  inline bool
  operator==(const SBufferCL<T, _Size>&, const SBufferCL<T, _Size>&);

/// Stores 2D coordinates
typedef SVectorCL<2> Point2DCL;
/// Stores 3D coordinates
typedef SVectorCL<3> Point3DCL;
/// Stores 4D coordinates
typedef SVectorCL<4> Point4DCL;
/// Stores barycentric coordinates
typedef SVectorCL<4> BaryCoordCL;
/// Stores barycentric coordinates of a space-time simplex (pentatope)
typedef SVectorCL<5> STBaryCoordCL;

/// Stores 4D coordinates
typedef SVectorCL<4> Point4DCL;


enum InitStateT { Uninitialized, Initialized };


//**************************************************************************
// Class:    SArrayCL                                                      *
// Purpose:  an array that remembers its size                              *
// Remarks:  All functions are inline, should be as fast as a "bare" array *
//**************************************************************************

template <class T, Uint _Size>
class SArrayCL
{
  private:
    T Array[_Size];

  public:
    typedef       T*       iterator;
    typedef const T* const_iterator;
    typedef       T&       reference;
    typedef const T& const_reference;
    typedef       T  value_type;

//    SArrayCL() {}
    explicit           SArrayCL(T val= T())        { std::fill_n(Array+0, _Size, val); }
    /*uninitialized memory; mainly for faster SVectorCL-math*/
    explicit SArrayCL(InitStateT) {}
    template<class In> explicit SArrayCL(In start) { std::copy(start, start+_Size, Array+0); }
    template<class In> SArrayCL(In start, In end)  { std::copy(start, end, Array+0); }
    // Default copy-ctor, assignment operator, dtor

    iterator       begin     ()             { return static_cast<T*>(Array); }
    const_iterator begin     ()       const { return static_cast<const T*>(Array); }
    iterator       end       ()             { return Array+_Size; }
    const_iterator end       ()       const { return Array+_Size; }
    reference      operator[](Uint i)       { Assert(i<_Size, DROPSErrCL("SArrayCL::operator[]: wrong index"), DebugContainerC);
                                              return Array[i]; }
    value_type     operator[](Uint i) const { Assert(i<_Size, DROPSErrCL("SArrayCL::operator[]: wrong index"), DebugContainerC);
                                              return Array[i]; }
    Uint           size      ()       const { return _Size; }

    friend bool operator==<>(const SArrayCL&, const SArrayCL&); // Component-wise equality
    friend bool operator< <>(const SArrayCL&, const SArrayCL&); // lexicographic ordering
};

template <class T, Uint _Size>
  inline bool
  operator==(const SArrayCL<T, _Size>& a0, const SArrayCL<T, _Size>& a1)
{
    for (Uint i=0; i<_Size; ++i)
        if (a0[i] != a1[i]) return false;
    return true;
}

template <class T, Uint _Size>
  inline bool
  operator<(const SArrayCL<T, _Size>& a0, const SArrayCL<T, _Size>& a1)
{
    for (Uint i=0; i<_Size; ++i)
        if (a0[i] < a1[i]) return true;
        else if ( a0[i] > a1[i]) return false;
    return false;
}

template <class T, Uint _Size>
  inline const T*
  Addr(const SArrayCL<T, _Size>& a)
{
    return a.begin();
}

template <class T, Uint _Size>
  inline T*
  Addr(SArrayCL<T, _Size>& a)
{
    return a.begin();
}

//**************************************************************************
// Class:    SBufferCL                                                     *
// Purpose:  A buffer or ring-buffer of fixed size with wrap-around indices*
//**************************************************************************
template <class T, Uint _Size>
class SBufferCL
{
private:
    T Array[_Size];
    int Front;

public:
    typedef       T&       reference;
    typedef const T& const_reference;
    typedef       T  value_type;

//    SBufferCL() {}
    explicit           SBufferCL(T val= T())        { std::fill_n(Array+0, _Size, val); Front= 0; }
    template<class In> explicit SBufferCL(In start) { std::copy(start, start+_Size, Array+0); Front= 0; }
    template<class In> SBufferCL(In start, In end)  { std::copy(start, end, Array+0); Front= 0; }
    SBufferCL(const SBufferCL& b) {
        std::copy( b.Array+0, b.Array+_Size, Array+0);
        Front= b.Front; }
    SBufferCL& operator=(const SBufferCL& b) {
        if (&b!=this) {
            std::copy( b.Array+0, b.Array+_Size, Array+0);
            Front= b.Front; }
        return *this; }
    // Default dtor

    reference  operator[](int i) {
        if (i<0) i+= _Size;
        else if (i >= static_cast<int>( _Size)) i-=_Size;
        Assert( (i>=0 && i<static_cast<int>( _Size)), DROPSErrCL("SBufferCL::operator[]: wrong index"), DebugContainerC);
        return Array[(i+Front < static_cast<int>( _Size)) ? i+Front : i+Front-_Size]; }
    value_type operator[](int i) const {
        if (i<0) i+= _Size;
        else if (i >= static_cast<int>( _Size)) i-=_Size;
        Assert( (i>=0 && i<_Size), DROPSErrCL("SBufferCL::operator[]: wrong index"), DebugContainerC);
        return Array[(i+Front < static_cast<int>( _Size)) ? i+Front : i+Front-_Size]; }
    Uint size() const { return _Size; }

    // first entry is overwritten.
    void push_back(const T& t) {
        Array[Front]= t;
        if (++Front == static_cast<int>( _Size)) Front-= _Size; }
    // for positive i, front element moves to the end.
    void rotate (int i= 1) {
        Assert(i<static_cast<int>( _Size), DROPSErrCL("SBufferCL::rotate: wrong index"), DebugContainerC);
        Front+= i;
        if (Front >= static_cast<int>( _Size)) Front-= _Size;
        else if (Front < 0) Front+= _Size; }


    friend bool operator==<>(const SBufferCL&, const SBufferCL&);
};

template <class T, Uint _Size>
  inline bool
  operator==(const SBufferCL<T, _Size>& a0, const SBufferCL<T, _Size>& a1)
{
    if (a0.Front != a1.Front) return false;
    for (Uint i=0; i<_Size; ++i)
        if (a0.Array[i] != a1.Array[i]) return false;
    return true;
}

//**************************************************************************
// Class:    SVectorCL                                                     *
// Purpose:  primitive vector class with templates for short vectors       *
// Remarks:  Many optimizations are possible.                              *
//**************************************************************************

template <Uint _Size>
class SVectorCL : public SArrayCL<double,_Size>
{
  public:
    typedef SArrayCL<double,_Size> base_type;

    SVectorCL()                                                             {}
    explicit           SVectorCL(InitStateT i)      : base_type( i)         {}
    explicit           SVectorCL(double val)        : base_type( val)       {}
    // constructors for Point3D/Point4D/BaryCoordCL 
    // (constructor throws Exception if number of arg. does not coincide with _Size)
    explicit SVectorCL(double arg1, double arg2, double arg3);
    explicit SVectorCL(double arg1, double arg2, double arg3, double arg4);
    template<Uint _Size2>
    explicit SVectorCL(SVectorCL<_Size2> arg1, double arg2);

    template<class In> explicit SVectorCL(In start) : base_type( start)     {}
    template<class In> SVectorCL(In start, In end)  : base_type( start,end) {}

    SVectorCL& operator+=(const SVectorCL&);
    SVectorCL& operator-=(const SVectorCL&);
    SVectorCL& operator*=(double);
    SVectorCL& operator/=(double);

    // komponentenweise Operatoren
    SVectorCL& operator*=(const SVectorCL&);
    SVectorCL& operator/=(const SVectorCL&);

    double norm_sq() const;
    double norm()    const { return std::sqrt(norm_sq()); }
};


template <Uint _Size>
template <Uint _Size2>
inline SVectorCL<_Size>::SVectorCL(SVectorCL<_Size2> arg1, double arg2) : base_type( Uninitialized )
{
    if (_Size != _Size2+1) // <- this is optimized away
    {
        throw DROPSErrCL("SVectorCL-classes have wrong size for constructor");
    }
    for (Uint i = 0; i < _Size2; i++)
        (*this)[i] = arg1[i];
    (*this)[_Size2] = arg2;
}

template <Uint _Size>
inline SVectorCL<_Size>::SVectorCL(double arg1, double arg2, double arg3) : base_type( Uninitialized )
{
    if (_Size != 3) // <- this is optimized away
        throw DROPSErrCL("SVectorCL-class has wrong size for constructor, number of args : 3");
    (*this)[0] = arg1;
    (*this)[1] = arg2;
    (*this)[2] = arg3;
}

template <Uint _Size>
inline SVectorCL<_Size>::SVectorCL(double arg1, double arg2, double arg3, double arg4) : base_type( Uninitialized )
{
    if (_Size != 4) // <- this is optimized away
        throw DROPSErrCL("SVectorCL-class has wrong size for constructor, number of args : 4");
    (*this)[0] = arg1;
    (*this)[1] = arg2;
    (*this)[2] = arg3;
    (*this)[3] = arg4;
}



template <Uint _Size>
SVectorCL<_Size>&
SVectorCL<_Size>::operator+=(const SVectorCL& v)
{
    for (Uint i=0; i!=_Size; ++i) (*this)[i]+= v[i];
    return *this;
}

template <Uint _Size>
SVectorCL<_Size>&
SVectorCL<_Size>::operator-=(const SVectorCL<_Size>& v)
{
    for (Uint i=0; i!=_Size; ++i) (*this)[i]-= v[i];
    return *this;
}

template <Uint _Size>
SVectorCL<_Size>&
SVectorCL<_Size>::operator*=(const SVectorCL<_Size>& v)
{
    for (Uint i=0; i!=_Size; ++i) (*this)[i]*= v[i];
    return *this;
}

template <Uint _Size>
SVectorCL<_Size>&
SVectorCL<_Size>::operator/=(const SVectorCL<_Size>& v)
{
    for (Uint i=0; i!=_Size; ++i) (*this)[i]/= v[i];
    return *this;
}

template <Uint _Size>
SVectorCL<_Size>&
SVectorCL<_Size>::operator*=(double s)
{
    for (Uint i=0; i!=_Size; ++i) (*this)[i]*= s;
    return *this;
}

template <Uint _Size>
SVectorCL<_Size>&
SVectorCL<_Size>::operator/=(double s)
{
    for (Uint i=0; i!=_Size; ++i) (*this)[i]/= s;
    return *this;
}

template <Uint _Size>
double SVectorCL<_Size>::norm_sq() const
{
    double x = 0.0;
    for (Uint i=0; i<_Size; ++i) x += (*this)[i]*(*this)[i];
    return x;
}

template <Uint _Size>
SVectorCL<_Size> BaryCenter(const SVectorCL<_Size>& v1, const SVectorCL<_Size>& v2)
{
    SVectorCL<_Size> tempv( Uninitialized);
    for (Uint i=0; i<_Size; ++i) tempv[i] = .5 * (v1[i] + v2[i]);
    return tempv;
}

template <Uint _Size>
SVectorCL<_Size> ConvexComb (double a,
                             const SVectorCL<_Size>& v1,
                             const SVectorCL<_Size>& v2)
{
    SVectorCL<_Size> tempv( Uninitialized);
    for (Uint i=0; i<_Size; ++i) tempv[i] = (1.0-a)*v1[i] + a*v2[i];
    return tempv;
}

template <Uint _Size>
SVectorCL<_Size> operator+(const SVectorCL<_Size>& v1,
                           const SVectorCL<_Size>& v2)
{
    SVectorCL<_Size> tempv( Uninitialized);
    for (Uint i=0; i<_Size; ++i) tempv[i]= v1[i] + v2[i];
    return tempv;
}

template <Uint _Size>
SVectorCL<_Size> operator-(const SVectorCL<_Size>& v1,
                           const SVectorCL<_Size>& v2)
{
    SVectorCL<_Size> tempv( Uninitialized);
    for (Uint i=0; i<_Size; ++i) tempv[i]= v1[i] - v2[i];
    return tempv;
}

template <Uint _Size>
SVectorCL<_Size> operator-(const SVectorCL<_Size>& v1)
{
    SVectorCL<_Size> tempv( Uninitialized);
    for (Uint i=0; i<_Size; ++i) tempv[i]= -v1[i];
    return tempv;
}

template <Uint _Size>
SVectorCL<_Size> operator*(double d, const SVectorCL<_Size>& v)
{
    SVectorCL<_Size> tempv( Uninitialized);
    for (Uint i=0; i<_Size; ++i) tempv[i]= d * v[i];
    return tempv;
}

template <Uint _Size>
SVectorCL<_Size> operator*(const SVectorCL<_Size>& v, double d)
{
    return d*v;
}

template <Uint _Size>
double inner_prod(const SVectorCL<_Size>& v1, const SVectorCL<_Size>& v2)
{
    double ret= 0.0;
    for (Uint i= 0; i <_Size; ++i) ret+= v1[i]*v2[i];
    return ret;
}

inline double
inner_prod(const SVectorCL<3u>& v1, const SVectorCL<3u>& v2)
{
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

template <Uint _Size>
SVectorCL<_Size> operator/(const SVectorCL<_Size>& v, double d)
{
    SVectorCL<_Size> tempv( Uninitialized);
    for (Uint i=0; i<_Size; ++i) tempv[i]= v[i]/d;
    return tempv;
}

template <Uint _Size>
SVectorCL<_Size> operator*(const SVectorCL<_Size>& v1,
                           const SVectorCL<_Size>& v2)
{
    SVectorCL<_Size> tempv( Uninitialized);
    for (Uint i=0; i<_Size; ++i) tempv[i]= v1[i] * v2[i];
    return tempv;
}

template <Uint _Size>
SVectorCL<_Size> operator/(const SVectorCL<_Size>& v1,
                           const SVectorCL<_Size>& v2)
{
    SVectorCL<_Size> tempv( Uninitialized);
    for (Uint i=0; i<_Size; ++i) tempv[i]= v1[i] / v2[i];
    return tempv;
}

template <Uint _Size>
SVectorCL<_Size> sqrt(const SVectorCL<_Size>& v)
{
    SVectorCL<_Size> tempv( Uninitialized);
    for (Uint i=0; i<_Size; ++i) tempv[i]= std::sqrt(v[i]);
    return tempv;
}

template <Uint _Size>
SVectorCL<_Size> fabs(const SVectorCL<_Size>& v)
{
    SVectorCL<_Size> tempv( Uninitialized);
    for (Uint i=0; i<_Size; ++i) tempv[i]= std::fabs(v[i]);
    return tempv;
}


template <Uint _Size>
std::ostream& operator<<(std::ostream& os, const SVectorCL<_Size>& v)
{
//    os << v.size() << "    ";
    for (Uint i=0; i<v.size(); ++i)
        os << v[i] << ' ';
    return os;
}


inline void
cross_product(Point3DCL& res, const Point3DCL& v0, const Point3DCL& v1)
// res= v0 x v1
{
    res[0]= v0[1]*v1[2] - v0[2]*v1[1];
    res[1]= v0[2]*v1[0] - v0[0]*v1[2];
    res[2]= v0[0]*v1[1] - v0[1]*v1[0];
}


//taken from thesis hollasch
inline void
cross_product(Point4DCL& res, const Point4DCL& u, const Point4DCL& v, const Point4DCL& w)
// res= X( v1, v2, v3); 
{
    const double a = v[0]*w[1]-v[1]*w[0];
    const double b = v[0]*w[2]-v[2]*w[0];
    const double c = v[0]*w[3]-v[3]*w[0];
    const double d = v[1]*w[2]-v[2]*w[1];
    const double e = v[1]*w[3]-v[3]*w[1];
    const double f = v[2]*w[3]-v[3]*w[2];
    res[0]=   u[1]*f - u[2]*e + u[3]*d;
    res[1]= - u[0]*f + u[2]*c - u[3]*b;
    res[2]=   u[0]*e - u[1]*c + u[3]*a;
    res[3]= - u[0]*d + u[1]*b - u[2]*a;
}

//taken from thesis hollasch
inline Point4DCL
cross_product(const Point4DCL& u, const Point4DCL& v, const Point4DCL& w)
// res= X( v1, v2, v3); 
{
    Point4DCL res;
    cross_product(res,u,v,w);
    return res;
}


inline double
l2norm_of_3D_vector(Point4DCL & in)
{
    return std::sqrt(in[0]*in[0]+in[1]*in[1]+in[2]*in[2]);
}

template <Uint _Size>
inline double
normalize_l2norm(SVectorCL<_Size> & in)
{
    double ret = in.norm();
    in /= ret;
    return ret;
}


// std_basis<n>(0)==Null, std_basis<n>(i)[j]==Delta_i-1_j
template <Uint _Size>
inline SVectorCL<_Size> std_basis(Uint i)
{
    SVectorCL<_Size> ret(0.);
    if (i>0) ret[i-1]= 1.;
    return ret;
}

/// By convention x_bary are barycentric coordinates with respect to a tetra; t_bary is a time in the (reference-) interval (0, 1);
struct STCoordCL
{
    BaryCoordCL x_bary;
    double      t_ref;

    STCoordCL () : x_bary(), t_ref( 0.) {} ///< Default init to 0.
    STCoordCL (const BaryCoordCL& b, double t)
        : x_bary( b), t_ref( t) {}

    STCoordCL& operator+= (const STCoordCL& v);
};

inline STCoordCL& STCoordCL::operator+= (const STCoordCL& v)
{
    x_bary+= v.x_bary;
    t_ref += v.t_ref;
    return *this;
}

inline STCoordCL operator+(const STCoordCL& v1,
                           const STCoordCL& v2)
{
    return STCoordCL( v1.x_bary + v2.x_bary, v1.t_ref + v2.t_ref);
}

inline STCoordCL operator*(double d, const STCoordCL& v)
{
    return STCoordCL( d*v.x_bary, d*v.t_ref);
}

inline STCoordCL operator*(const STCoordCL& v, double d)
{
    return STCoordCL( d*v.x_bary, d*v.t_ref);
}

inline STCoordCL ConvexComb (double a,
                               const STCoordCL& v1,
                               const STCoordCL& v2)
{
    return STCoordCL( ConvexComb( a, v1.x_bary, v2.x_bary), (1.0-a)*v1.t_ref + a*v2.t_ref);
}

inline STCoordCL MakeSTCoord (const BaryCoordCL& b, double t)
{
    return STCoordCL( b, t);
}

inline std::ostream& operator<< (std::ostream& os, const STCoordCL& st)
{
    return os << st.x_bary << st.t_ref;
}

inline STBaryCoordCL
MakeSTBaryCoord(double a, double b, double c, double d, double e)
{
    STBaryCoordCL ret( Uninitialized);
    ret[0]= a; ret[1]= b; ret[2]= c; ret[3]= d; ret[4]= e;
    return ret;
}

inline STBaryCoordCL
MakeSTBaryCoord(const Point4DCL& p)
{
    STBaryCoordCL ret( Uninitialized);
    ret[0]= 1-p[0]-p[1]-p[2]-p[3];
    ret[1] = p[0];
    ret[2] = p[1];
    ret[3] = p[2];
    ret[4] = p[3];
    return ret;
}


inline BaryCoordCL
MakeBaryCoord(double a, double b, double c, double d)
{
    BaryCoordCL ret( Uninitialized);
    ret[0]= a; ret[1]= b; ret[2]= c; ret[3]= d;
    return ret;
}

inline BaryCoordCL
MakeBaryCoord(double a, double b, double c)
{
    BaryCoordCL ret( Uninitialized);
    ret[0]= 1-a-b-c;
    ret[1] = a;
    ret[2] = b;
    ret[3] = c;
    return ret;
}


inline Point3DCL
BaryToRefCoord(const BaryCoordCL & p)
{
    Point3DCL ret( Uninitialized);
    ret[0]= p[1]; ret[1]= p[2]; ret[2]= p[3];
    return ret;
}

inline Point3DCL
RestrictToPoint3D(const Point4DCL & p)
{
    Point3DCL ret( Uninitialized);
    ret[0]= p[0]; ret[1]= p[1]; ret[2]= p[2];
    return ret;
}

inline BaryCoordCL
MakeBaryCoord(const Point3DCL& p)
{
    BaryCoordCL ret( Uninitialized);
    ret[0]= 1-p[0]-p[1]-p[2];
    ret[1] = p[0];
    ret[2] = p[1];
    ret[3] = p[2];
    return ret;
}


inline Point4DCL
MakePoint4D(Point3DCL a, double d) /*deprecated*/
{
    Point4DCL ret( Uninitialized);
    ret[0]= a[0]; ret[1]= a[1]; ret[2]= a[2]; ret[3]= d;
    return ret;
}

inline Point4DCL
MakePoint4D(double a, double b, double c, double d) /*deprecated*/
{
    Point4DCL ret( Uninitialized);
    ret[0]= a; ret[1]= b; ret[2]= c; ret[3]= d;
    return ret;
}

inline Point3DCL
MakePoint3D(double a, double b, double c)
{
    Point3DCL ret( Uninitialized);
    ret[0]= a; ret[1]= b; ret[2]= c;
    return ret;
}

inline Point2DCL
MakePoint2D(double a, double b)
{
    Point2DCL ret( Uninitialized);
    ret[0]= a; ret[1]= b;
    return ret;
}

template<class T>
SArrayCL<T, 2>
MakeSArray(T a, T b)
{
    SArrayCL<T, 2> ret( Uninitialized);
    ret[0]= a; ret[1]= b;
    return ret;
}

template<class T>
SArrayCL<T, 3>
MakeSArray(T a, T b, T c)
{
    SArrayCL<T, 3> ret( Uninitialized);
    ret[0]= a; ret[1]= b; ret[2]= c;
    return ret;
}

template<class T>
SArrayCL<T, 4>
MakeSArray(T a, T b, T c, T d)
{
    SArrayCL<T, 4> ret( Uninitialized);
    ret[0]= a; ret[1]= b; ret[2]= c; ret[3]= d;
    return ret;
}

template<class T>
SArrayCL<T, 5>
MakeSArray(T a, T b, T c, T d, T e)
{
    SArrayCL<T, 5> ret( Uninitialized);
    ret[0]= a; ret[1]= b; ret[2]= c; ret[3]= d; ret[4]= e;
    return ret;
}

template <class T, Uint _Size>
std::ostream& operator<<(std::ostream& os, const SArrayCL<T,_Size>& a)
{
//    os << v.size() << "    ";
    for (Uint i=0; i<a.size(); ++i)
        os << a[i] << ' ';
    return os;
}

template <Uint _Rows, Uint _Cols>
class SMatrixCL : public SVectorCL<_Rows*_Cols>
{
  public:
    typedef SVectorCL<_Rows*_Cols> _vec_base;

    SMatrixCL()                                                             {}
    explicit           SMatrixCL(InitStateT i)      : _vec_base( i)         {}
    explicit           SMatrixCL(double val)        : _vec_base( val)       {}
    template<class In> explicit SMatrixCL(In start) : _vec_base( start)     {}
    template<class In> SMatrixCL(In start, In end)  : _vec_base( start,end) {}

// Schreib- & Lesezugriff
    double& operator() (int row, int col)       { return (*this)[row*_Cols+col]; }// Matrix(i,j)
    double  operator() (int row, int col) const { return (*this)[row*_Cols+col]; }

    SVectorCL<_Rows> col( int) const;
    void             col( int, const SVectorCL<_Rows>&);

// Zuweisung & Co.
    SMatrixCL& operator+=(const SMatrixCL&);                // Matrix=Matrix+Matrix'
    SMatrixCL& operator-=(const SMatrixCL&);                // Matrix=Matrix-Matrix'
    SMatrixCL& operator*=(double s);                               // Matrix = c * Matrix
    SMatrixCL& operator/=(double s);                               // Matrix = Matrix'/c

// Dimensionen feststellen
    Uint num_rows() const { return _Rows; }                        // Zeilenzahl
    Uint num_cols() const { return _Cols; }                        // Spaltenzahl
};

template<Uint _Rows, Uint _Cols>
SVectorCL<_Rows>
SMatrixCL<_Rows, _Cols>::col (int c) const
{
    SVectorCL<_Rows> ret( Uninitialized);
    for (Uint i= 0; i != _Rows; ++i, c+= _Cols)
        ret[i]= (*this)[c];
    return ret;
}

template<Uint _Rows, Uint _Cols>
void
SMatrixCL<_Rows,_Cols>::col (int c, const SVectorCL<_Rows>& v)
{
    for (Uint i= 0; i != _Rows; ++i)
        (*this)( i, c)= v[i];
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>&
SMatrixCL<_Rows, _Cols>::operator+=(const SMatrixCL<_Rows, _Cols>& m)
{
    *static_cast<_vec_base*>(this)+= *static_cast<const _vec_base*>(&m);
    return *this;
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>&
SMatrixCL<_Rows, _Cols>::operator-=(const SMatrixCL<_Rows, _Cols>& m)
{
    *static_cast<_vec_base*>(this)-= *static_cast<const _vec_base*>(&m);
    return *this;
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>&
SMatrixCL<_Rows, _Cols>::operator*=(double d)
{
    *static_cast<_vec_base*>(this)*= d;
    return *this;
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>&
SMatrixCL<_Rows, _Cols>::operator/=(double d)
{
    *static_cast<_vec_base*>(this)/= d;
    return *this;
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>
operator+(const SMatrixCL<_Rows, _Cols>& m1, const SMatrixCL<_Rows, _Cols>& m2)
{
    SMatrixCL<_Rows, _Cols> ret( Uninitialized);
    *static_cast<typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&ret)
        = *static_cast<const typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&m1)
         +*static_cast<const typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&m2);
    return ret;
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>
operator-(const SMatrixCL<_Rows, _Cols>& m1, const SMatrixCL<_Rows, _Cols>& m2)
{
    SMatrixCL<_Rows, _Cols> ret( Uninitialized);
    *static_cast<typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&ret)
        = *static_cast<const typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&m1)
         -*static_cast<const typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&m2);
    return ret;
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>
operator-(const SMatrixCL<_Rows, _Cols>& m)
{
    SMatrixCL<_Rows, _Cols> ret( Uninitialized);
    *static_cast<typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&ret)
        = -*static_cast<const typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&m);
    return ret;
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>
operator*(double d, const SMatrixCL<_Rows, _Cols>& m)
{
    SMatrixCL<_Rows, _Cols> ret( Uninitialized);
    *static_cast<typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&ret)
        = d**static_cast<const typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&m);
    return ret;
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>
operator*(const SMatrixCL<_Rows, _Cols>& m, double d)
{
    SMatrixCL<_Rows, _Cols> ret( Uninitialized);
    *static_cast<typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&ret)
        = *static_cast<const typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&m)*d;
    return ret;
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>
operator/(const SMatrixCL<_Rows, _Cols>& m, double d)
{
    SMatrixCL<_Rows, _Cols> ret( Uninitialized);
    *static_cast<typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&ret)
        = *static_cast<const typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&m)/d;
    return ret;
}

template<Uint _RowsL, Uint _ColsR, Uint _Dim>
SMatrixCL<_RowsL, _ColsR>
operator*(const SMatrixCL<_RowsL, _Dim>& m1, const SMatrixCL<_Dim, _ColsR>& m2)
{
    SMatrixCL<_RowsL, _ColsR> ret(0.0);
    for (Uint row=0; row!=_RowsL; ++row)
        for (Uint col=0; col!=_ColsR; ++col)
            for (Uint i=0; i!=_Dim; ++i)
                ret(row, col)+= m1(row, i)*m2(i, col);
    return ret;
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Cols, _Cols>
GramMatrix(const SMatrixCL<_Rows, _Cols>& m)
/// Computes m^T*m
{
    SMatrixCL<_Cols, _Cols> ret( 0.0);
    for (Uint row= 0; row != _Cols; ++row) {
        for (Uint col= 0; col < row; ++col) {
            for (Uint i= 0; i != _Rows; ++i)
                ret( row, col)+= m( i, row)*m( i, col);
            ret( col, row)= ret( row, col);
        }
        for (Uint i= 0; i != _Rows; ++i)
            ret( row, row)+= std::pow( m( i, row), 2);
    }
    return ret;
}

template<Uint _Rows, Uint _Cols>
SVectorCL<_Cols>
transp_mul(const SMatrixCL<_Rows, _Cols>& m, const SVectorCL<_Rows>& v)
{
    SVectorCL<_Cols> ret(0.0);
    for (Uint col=0; col!=_Cols; ++col)
        for (Uint i=0; i!=_Rows; ++i)
                ret[col]+= m( i, col)*v[i];
    return ret;
}

template<Uint _Rows, Uint _Cols>
SVectorCL<_Rows>
operator*(const SMatrixCL<_Rows, _Cols>& m, const SVectorCL<_Cols>& v)
{
    SVectorCL<_Rows> ret(0.0);
    for (Uint row=0; row!=_Rows; ++row)
        for (Uint i=0; i!=_Cols; ++i)
                ret[row]+= m(row, i)*v[i];
    return ret;
}

inline SVectorCL<3>
operator*(const SMatrixCL<3, 3>& m, const SVectorCL<3>& v)
{
    SVectorCL<3> ret( Uninitialized);
    const double* const a= m.begin();
    ret[0]= a[0]*v[0] + a[1]*v[1] + a[2]*v[2];
    ret[1]= a[3]*v[0] + a[4]*v[1] + a[5]*v[2];
    ret[2]= a[6]*v[0] + a[7]*v[1] + a[8]*v[2];
    return ret;
}

inline SVectorCL<4>
operator*(const SMatrixCL<4, 4>& m, const SVectorCL<4>& v)
{
    SVectorCL<4> ret( Uninitialized);
    const double* const a= m.begin();
    ret[0]=  a[0]*v[0] +  a[1]*v[1] +  a[2]*v[2] +  a[3]*v[3];
    ret[1]=  a[4]*v[0] +  a[5]*v[1] +  a[6]*v[2] +  a[7]*v[3];
    ret[2]=  a[8]*v[0] +  a[9]*v[1] + a[10]*v[2] + a[11]*v[3];
    ret[3]= a[12]*v[0] + a[13]*v[1] + a[14]*v[2] + a[15]*v[3];
    return ret;
}


template <Uint Rows, Uint Cols>
double inner_prod(const SVectorCL<Rows>& v1, const SMatrixCL<Rows, Cols>& m, const SVectorCL<Cols>& v2)
{
    double ret= 0.0;
    for (Uint i= 0; i < Rows; ++i)
        for (Uint j= 0; j < Cols; ++j)
            ret+= v1[i]*m( i, j)*v2[j];
    return ret;
}

template <Uint _Rows, Uint _Cols>
std::ostream& operator << (std::ostream& os, const SMatrixCL<_Rows, _Cols>& m)
{
    const Uint M = m.num_rows();
    const Uint N = m.num_cols();

    os << M << ' ' << N << '\n' ;

    for (Uint row=0; row<M; ++row)
    {
        os << "  ";
        for (Uint col=0; col<N; ++col)
            os << m(row, col) << ' ';
        os << '\n';
    }
    return os;
}

template <Uint _Rows>
inline SMatrixCL<_Rows,_Rows>
outer_product (const SVectorCL<_Rows>& a, const SVectorCL<_Rows>& b)
{
    SMatrixCL<_Rows,_Rows> ret( Uninitialized);
    for (Uint i= 0; i < _Rows; ++i)
        for (Uint j= 0; j < _Rows; ++j)
            ret(i, j)= a[i]*b[j];
    return ret;
}

template <Uint _Rows>
inline double
frobenius_norm_sq (const SMatrixCL<_Rows, _Rows>& a)
{
    double ret = 0;
    for (Uint i= 0; i < _Rows*_Rows; ++i)
        ret += a[i]*a[i];
    return ret;
}

template <Uint _Rows>
inline double
trace (const SMatrixCL<_Rows, _Rows>& a)
{
    double ret= 0.;
    for (Uint i= 0; i < _Rows; ++i)
            ret+= a( i, i);
    return ret;
}

template <Uint _Rows>
inline SMatrixCL<_Rows,_Rows>&
assign_transpose (SMatrixCL<_Rows, _Rows>& out, const SMatrixCL<_Rows,_Rows>& in)
{
    for (Uint i= 0; i < _Rows; ++i) {
        for (Uint j= 0; j < i; ++j) {
            out( i, j)= in( j, i);
            out( j, i)= in( i, j);
        }
        out( i, i)= in( i, i);
    }
    return out;
}

template <Uint _Rows, Uint _Cols>
inline SMatrixCL<_Rows, _Cols> eye ()
{
    SMatrixCL<_Rows, _Cols> ret;

    const Uint m= std::min( _Rows, _Cols);
    for (Uint i= 0; i < m; ++i)
        ret( i, i)= 1.;
    return ret;
}

/// \brief \f$full_local+= (scalar_local^T) \operatorname{kroneckerproduct} Id_{3\times 3}\f$
///
/// This is the operation that distributes a scalar-valued operator over the block-diagonal of a vector-valued operator.
inline void
add_transpose_kronecker_id (SMatrixCL<3,3> full_local[10][10], const double scalar_local[10][10])
{
    for(int i= 0; i < 10; ++i)
        for(int j= 0; j < 10; ++j)
            for (int k= 0; k < 3; ++k)
                full_local[i][j]( k, k)+= scalar_local[j][i];
}

///\brief A small diagonal matrix. It is needed as distinct type for SparseMatBuilderCL for block diagonal sparse matrices.
template <Uint _Rows>
class SDiagMatrixCL : public SVectorCL<_Rows>
{
  public:
    typedef SVectorCL<_Rows> _vec_base;

    SDiagMatrixCL()                                                             {}
    explicit           SDiagMatrixCL(InitStateT i)      : _vec_base( i)         {}
    explicit           SDiagMatrixCL(double val)        : _vec_base( val)       {}
    template<class In> explicit SDiagMatrixCL(In start) : _vec_base( start)     {}
    template<class In> SDiagMatrixCL(In start, In end)  : _vec_base( start,end) {}

// Schreib- & Lesezugriff
    double& operator() (int row)       { return (*this)[row]; }// Matrix(i,i)
    double  operator() (int row) const { return (*this)[row]; }

// Zuweisung & Co.
    SDiagMatrixCL& operator+=(const SDiagMatrixCL&);                // Matrix=Matrix+Matrix'
    SDiagMatrixCL& operator-=(const SDiagMatrixCL&);                // Matrix=Matrix-Matrix'
    SDiagMatrixCL& operator*=(double s);                            // Matrix = c * Matrix
    SDiagMatrixCL& operator/=(double s);                            // Matrix = Matrix'/c

// Dimensionen feststellen
    Uint num_rows() const { return _Rows; }                        // Zeilenzahl
    Uint num_cols() const { return _Rows; }                        // Spaltenzahl
};

template<Uint _Rows>
SDiagMatrixCL<_Rows>&
SDiagMatrixCL<_Rows>::operator+=(const SDiagMatrixCL<_Rows>& m)
{
    *static_cast<_vec_base*>(this)+= *static_cast<const _vec_base*>(&m);
    return *this;
}

template<Uint _Rows>
SDiagMatrixCL<_Rows>&
SDiagMatrixCL<_Rows>::operator-=(const SDiagMatrixCL<_Rows>& m)
{
    *static_cast<_vec_base*>(this)-= *static_cast<const _vec_base*>(&m);
    return *this;
}

template<Uint _Rows>
SDiagMatrixCL<_Rows>&
SDiagMatrixCL<_Rows>::operator*=(double d)
{
    *static_cast<_vec_base*>(this)*= d;
    return *this;
}

template<Uint _Rows>
SDiagMatrixCL<_Rows>&
SDiagMatrixCL<_Rows>::operator/=(double d)
{
    *static_cast<_vec_base*>(this)/= d;
    return *this;
}

/// \brief A QR-factored, rectangular matrix, A=QR.
///
/// This allows for fast application of A^{-1} and A.
/// A can have more rows than columns. In this case, the least-squares-solution is computed.
template <Uint Rows_, Uint Cols_= Rows_>
class QRDecompCL
{
  private:
    SMatrixCL<Rows_,Cols_> a_;
    double d_[Cols_]; ///< The diagonal of R
    double beta_[Cols_]; ///< The reflections are R_j= I + beta_j*a_[j:Rows_-1][j]

  public:
    QRDecompCL () : a_( Uninitialized) {}
    template <class MatT>
      QRDecompCL (MatT m, bool assume_full_rank= true)
        : a_( m) { prepare_solve ( assume_full_rank); }

    SMatrixCL<Rows_,Cols_>&       GetMatrix ()       { return a_; }
    const SMatrixCL<Rows_,Cols_>& GetMatrix () const { return a_; }

    ///\brief Computes the factorization. By default, Throws if the matrix is numerically ranc-deficient.
    ///\return rank A < Rows_.
    bool prepare_solve (bool assume_full_rank= true);
    ///@{ Call only after prepare_solve; solves are inplace. For least-squares, the first Cols_ entries are the least squares solution, the remaining components of b are the residual-vector.
    void Solve (SVectorCL<Rows_>& b) const;
    void Solve (size_t n, SVectorCL<Rows_>* b) const;
    template <template<class> class SVecCont>
    void Solve (SVecCont<SVectorCL<Rows_> >& b) const;
    template <Uint Size>
    void Solve (SArrayCL<SVectorCL<Rows_>, Size>& b) const;
    template <Uint Size>
    void Solve (SMatrixCL<Rows_, Size>& b) const;

    ///@{ Apply (parts of) Q or Q^T
    void apply_reflection  (Uint j, SVectorCL<Rows_>& b) const;
    void apply_Q_transpose (SVectorCL<Rows_>& b) const;
    void apply_Q           (SVectorCL<Rows_>& b) const;
    SMatrixCL<Rows_, Rows_> get_Q () const;
    ///@}
    ///\brief Computes the determinant of R (stable). For Rows_ > Cols_, the determinant of the upper Cols_ x Cols_ block of R is returned.
    double Determinant_R () const;
    ///@}

    ///@{ Serialize and Deserialize a QR decomposition
    void Serialize(double*) const;
    void Deserialize(const double*);
    ///@}
};

template <Uint Rows_, Uint Cols_>
  inline void
    QRDecompCL<Rows_, Cols_>::apply_reflection (Uint j, SVectorCL<Rows_>& b) const
{
    double sp= 0.;
    for(Uint i= j; i < Rows_; ++i)
        sp+= a_(i, j) * b[i];
    sp*= beta_[j];
    for(Uint i= j; i < Rows_; ++i)
    b[i]+= a_(i, j)*sp;
}

template <Uint Rows_, Uint Cols_>
  inline void
    QRDecompCL<Rows_, Cols_>::apply_Q_transpose (SVectorCL<Rows_>& b) const
{
    for (Uint j= 0; j < Cols_; ++j)
        apply_reflection( j, b);
}

template <Uint Rows_, Uint Cols_>
  inline void
    QRDecompCL<Rows_, Cols_>::apply_Q (SVectorCL<Rows_>& b) const
{
    // The Q_j are symmetric, so to transpose Q^T one just has to reverse the order of application.
    for (Uint j= Cols_ - 1; j < Cols_; --j)
        apply_reflection( j, b);
}

template <Uint Rows_, Uint Cols_>
  SMatrixCL<Rows_, Rows_>
    QRDecompCL<Rows_, Cols_>::get_Q () const
{
    SMatrixCL<Rows_, Rows_> Q (eye<Rows_, Rows_> ());
    SVectorCL<Rows_> tmp;
    for (Uint j= 0; j < Rows_; ++j) {
        tmp= Q.col( j);
        apply_Q( tmp);
        Q.col( j, tmp);
    }
    return Q;
}

template <Uint Rows_, Uint Cols_>
  double
  QRDecompCL<Rows_, Cols_>::Determinant_R () const
{
    double tmp= 1.0;
    for(Uint i= 0; i < Cols_; ++i)
        tmp *= d_[i];
    return tmp;
}

template <Uint Rows_, Uint Cols_>
  bool
  QRDecompCL<Rows_, Cols_>::prepare_solve (bool assume_full_rank)
{
    // inplace Householder
    double sigma, sp;
    for (Uint j= 0; j < Cols_; ++j) {
        sigma = 0.;
        for(Uint i= j; i < Rows_; ++i)
            sigma+= std::pow( a_(i, j), 2);
        if(sigma == 0.) {
            if (assume_full_rank)
                throw DROPSErrCL( "QRDecompCL::prepare_solve: rank-deficient matrix\n");
            else
                return true;
        }
        d_[j]= (a_(j, j) < 0 ? 1. : -1.) * std::sqrt( sigma);
        beta_[j]= 1./(d_[j]*a_(j, j) - sigma);
        a_(j, j)-= d_[j];
        for(Uint k= j + 1; k < Cols_; ++k) { // Apply reflection in column j
            sp= 0.;
            for(Uint i= j; i < Rows_; ++i)
                sp+= a_(i, j) * a_(i, k);
            sp*= beta_[j];
            for(Uint i= j; i < Rows_; ++i)
                a_(i, k)+= a_(i, j)*sp;
        }
    }
    return false;
}

template <Uint Rows_, Uint Cols_>
 void
 QRDecompCL<Rows_, Cols_>::Solve (size_t n, SVectorCL<Rows_>* b) const
{
    for (Uint i= 0; i < n; ++i)
        Solve( b[i]);
}

template <Uint Rows_, Uint Cols_>
  void
  QRDecompCL<Rows_, Cols_>::Solve (SVectorCL<Rows_>& b) const
{
    apply_Q_transpose( b);
    for (Uint i= Cols_ - 1; i < Cols_; --i) { // backsolve
        for (Uint j= i + 1; j < Cols_; ++j)
            b[i]-= a_(i, j)*b[j];
        b[i]/= d_[i];
    }
}

template <Uint Rows_, Uint Cols_>
  template <template<class> class SVecCont>
  void
  QRDecompCL<Rows_, Cols_>::Solve (SVecCont<SVectorCL<Rows_> >& b) const
{
    for (Uint i= 0; i < b.size(); ++i)
        Solve( b[i]);
}

template <Uint Rows_, Uint Cols_>
  template <Uint Size>
    void
    QRDecompCL<Rows_, Cols_>::Solve (SArrayCL<SVectorCL<Rows_>, Size>& b) const
{
    for (Uint i= 0; i < Size; ++i)
        Solve( b[i]);
}

template <Uint Rows_, Uint Cols_>
  template <Uint Size>
    void
    QRDecompCL<Rows_, Cols_>::Solve (SMatrixCL<Rows_, Size>& b) const
{
    SVectorCL<Rows_> c( Uninitialized);
    for (Uint i= 0; i < Size; ++i) {
        c= b.col (i);
        Solve( c);
        b.col( i, c);
    }
}

/** Put the values of a_, d_ and beta_ in buffer. Note that buffer must be of size
    (Rows_ + 2)*Cols_
 */
template <Uint Rows_, Uint Cols_>
  void QRDecompCL<Rows_, Cols_>::Serialize(double* buffer) const
{
    std::copy( a_.begin(), a_.end(), buffer);
    std::copy( d_, d_ + Cols_, buffer + a_.size());
    std::copy( beta_, beta_ + Cols_, buffer + a_.size() + Cols_);
}

template <Uint Rows_, Uint Cols_>
  void QRDecompCL<Rows_, Cols_>::Deserialize( const double* buffer)
{
    std::copy( buffer, buffer + a_.size(), a_.begin());
    std::copy( buffer + a_.size(), buffer + a_.size() + Cols_, d_);
    std::copy( buffer + a_.size() + Cols_, buffer + a_.size() + 2*Cols_, beta_);
}

/// \brief Represents the Given-/Jacobi-rotation J=(c & s \\ -s & c) with respect to row/column-pairs p and q.
class GivensRotationCL
{
  private:
    size_t p_, q_; ///< Row/column-pairs i=j=p and i=j=q.
    double c_, s_; ///< Cosine and sine.

  public:

    /// \brief Let A be a symmetric matrix, B= J^TAJ, then, B_{pq} = B_{qp} = 0.
    template <typename Mat>
    void inline compute_symmetric_rotation (const Mat& A, size_t p, size_t q);

    ///\brief Compute A= J^TA.
    template <typename Mat>
    void apply_transpose_from_left (Mat& A) const;

    ///\brief Compute A= AJ.
    template <typename Mat>
    void apply_from_right (Mat& A) const;

    ///\brief Compute A= J^TAJ.
    template <typename Mat>
    void inline apply_conjugation (Mat& A) const;
};

template <typename Mat>
void inline GivensRotationCL::compute_symmetric_rotation (const Mat& A, size_t p, size_t q)
{
    c_= 1.;
    s_= 0.;
    p_= p;
    q_= q;
    if (A( p, q) == 0.)
        return;

    const double tau= (A( q, q) - A( p, p))/(2.*A( p, q)),
                 t= tau >= 0. ?  1./( tau + std::sqrt( 1. + tau*tau))
                              : -1./(-tau + std::sqrt( 1. + tau*tau));
    c_= 1./std::sqrt( 1. + t*t);
    s_= t*c_;
}

template <typename Mat>
void GivensRotationCL::apply_transpose_from_left (Mat& A) const
{
    double Apj, Aqj;
    const size_t n= A.num_cols();
    for (size_t j= 0; j < n; ++j) {
        Apj= A( p_, j);
        Aqj= A( q_, j);
        A( p_, j)= c_*Apj - s_*Aqj;
        A( q_, j)= s_*Apj + c_*Aqj;
    }
}

template <typename Mat>
void GivensRotationCL::apply_from_right (Mat& A) const
{
    double Aip, Aiq;
    const size_t n= A.num_rows();
    for (size_t i= 0; i < n; ++i) {
        Aip= A( i, p_);
        Aiq= A( i, q_);
        A( i, p_)= c_*Aip - s_*Aiq;
        A( i, q_)= s_*Aip + c_*Aiq;
    }
}

template <typename Mat>
void GivensRotationCL::apply_conjugation (Mat& A) const
{
    double Aip, Aiq, Api, Aqi;
    const size_t n= A.num_rows();
    for (size_t i= 0; i < n; ++i) {
        if (i == p_ || i == q_)
            continue;
        Aip= A( i, p_);
        Aiq= A( i, q_);
        Api= A( p_, i);
        Aqi= A( q_, i);
        A( i, p_)= c_*Aip - s_*Aiq;
        A( i, q_)= s_*Aip + c_*Aiq;
        A( p_, i)= c_*Api - s_*Aqi;
        A( q_, i)= s_*Api + c_*Aqi;
    }
    const double App= A( p_, p_), Apq= A( p_, q_), Aqp= A( q_, p_), Aqq= A( q_, q_),
                 N= c_*s_*(Apq + Aqp), D= c_*s_*(App - Aqq),
                 csq= c_*c_, ssq= s_*s_;
    A( p_, p_)= csq*App + ssq*Aqq - N;
    A( p_, q_)= csq*Apq - ssq*Aqp + D;
    A( q_, p_)= csq*Aqp - ssq*Apq + D;
    A( q_, q_)= csq*Aqq + ssq*App + N;
}

template <typename Mat>
double off_diagonal_sq (const Mat& A)
{
    double off= 0.;
    const size_t m= A.num_rows(),
                 n= A.num_cols();
    for (size_t i= 0; i < m; ++i)
        for (size_t j= 0; j < n; ++j)
            if (j != i)
                off+= std::pow( A( i, j), 2);
    return off;
}

///\brief Computes all eigenvalues of a symmetric matrix A up to tolerance tol; overwrites A.
///
/// If Q is given, then AQ=Q\Lambda, where \Lambda is the diagonal matrix returned in A and A is the original A.
/// Q is an orthogonal matrix.
template <typename Mat>
void cyclic_jacobi (Mat& A, double tol, Mat* Q= 0)
{
    const size_t n= A.num_rows();
    if (Q) {
        std::fill_n( &(*Q)(0, 0), n*n, 0.);
        for (size_t i= 0; i < n; ++i)
            (*Q)( i, i)= 1.;
    }
    GivensRotationCL J;
    const double eps= tol*std::sqrt( frobenius_norm_sq( A));
    while (std::sqrt( off_diagonal_sq( A)) > eps) {
        for (size_t p= 0; p < n - 1; ++p)
            for (size_t q= p + 1; q < n; ++q) {
                J.compute_symmetric_rotation( A, p, q);
                J.apply_conjugation( A);
                if (Q)
                    J.apply_from_right( *Q);
            }
    }
}

//**************************************************************************
// Class:   GlobalListCL                                                   *
// Purpose: A list that is subdivided in levels. For modifications, it can *
//          efficiently be split into std::lists per level and then merged *
//          after modifications.                                           *
// Remarks: Negative level-indices count backwards from end().             *
//**************************************************************************
template <class T>
class GlobalListCL
{
  public:
    typedef std::list<T>                  Cont;
    typedef std::list<T>                  LevelCont;
    typedef typename Cont::iterator       iterator;
    typedef typename Cont::const_iterator const_iterator;
    typedef typename LevelCont::iterator       LevelIterator;
    typedef typename LevelCont::const_iterator const_LevelIterator;

  private:
    Cont                        Data_;
    std::vector<iterator>       LevelStarts_;
    std::vector<const_iterator> const_LevelStarts_;
    std::vector<LevelCont*>     LevelViews_;
    bool                        modifiable_;

    int
    StdLevel (int lvl) const { return lvl >= 0 ? lvl : lvl + GetNumLevel(); }

  public:
    GlobalListCL (bool modifiable= true) : modifiable_( modifiable) {}
    // standard dtor

    Uint GetNumLevel () const {
        return modifiable_ ? LevelViews_.size()
            : (LevelStarts_.size() > 0 ? LevelStarts_.size()-1 : 0);
    }
    bool IsEmpty () const
        { return modifiable_ ? LevelViews_.empty() : LevelStarts_.empty(); }
    bool IsLevelEmpty (Uint Level) const {
        return modifiable_ ? LevelViews_[Level]->empty()
            : LevelStarts_[Level] == LevelStarts_[1+Level];
    }

    // Only useful, if modifiable_ == false, otherwise 0.
    size_t size () const { return Data_.size(); }
    // If modifiable_==true, Data_ is empty, thus these accessors are useless.
          iterator begin ()       { return Data_.begin(); }
          iterator end   ()       { return Data_.end(); }
    const_iterator begin () const { return Data_.begin(); }
    const_iterator end   () const { return Data_.end(); }

    iterator level_begin (int lvl)
        { return !modifiable_ ? LevelStarts_[StdLevel( lvl)] : LevelViews_[StdLevel( lvl)]->begin(); }
    iterator level_end   (int lvl)
        { return !modifiable_ ? LevelStarts_[StdLevel( lvl) + 1] : LevelViews_[StdLevel( lvl)]->end(); }
    const_iterator level_begin (int lvl) const
        { return !modifiable_ ? const_LevelStarts_[StdLevel( lvl)] : LevelViews_[StdLevel( lvl)]->begin(); }
    const_iterator level_end   (int lvl) const
        { return !modifiable_ ? const_LevelStarts_[StdLevel( lvl) + 1] : LevelViews_[StdLevel( lvl)]->end(); }

    // Split Data_ into level-wise lists in LevelViews_ or merge LevelViews_ into Data_.
    void PrepareModify();
    void FinalizeModify();

    void AppendLevel();
    void RemoveLastLevel();

    LevelCont& // Only usable, if modifiable_ == true
    operator[] (int lvl) {
        Assert( modifiable_, DROPSErrCL("GlobalListCL::operator[]: "
            "Data not modifiable."), DebugContainerC );
        return *LevelViews_[StdLevel( lvl)];
    }
};

template <class T>
  void // Split Data_ into level-wise lists in LevelViews_.
  GlobalListCL<T>::PrepareModify()
{
    Assert( !modifiable_, DROPSErrCL("GlobalListCL::PrepareModify:"
        "Data is already modifiable."), DebugContainerC );
    Assert( LevelViews_.empty(), DROPSErrCL("GlobalListCL::PrepareModify:"
        "Inconsistent LevelViews_."), DebugContainerC );
    LevelViews_.resize( GetNumLevel());
    for (Uint lvl= 0, numlvl= GetNumLevel(); lvl < numlvl; ++lvl) {
        LevelViews_[lvl]= new LevelCont();
        LevelViews_[lvl]->splice( LevelViews_[lvl]->end(), Data_,
            level_begin( lvl), level_begin( lvl+1));
    }
    LevelStarts_.clear();
    const_LevelStarts_.clear();
    Assert( Data_.empty(), DROPSErrCL("GlobalListCL::PrepareModify: "
        "Did not move all Data."), DebugContainerC );
    modifiable_= true;
}

template <class T>
  void // Merge LevelViews_ into Data_.
  GlobalListCL<T>::FinalizeModify()
{
    Assert( modifiable_,
        DROPSErrCL("GlobalListCL::FinalizeModify: Data is not modifiable."),
        DebugContainerC );
    Assert( LevelStarts_.empty(),
        DROPSErrCL("GlobalListCL::FinalizeModify: Inconsistent LevelStarts_."),
        DebugContainerC );
    Assert( Data_.empty(), DROPSErrCL("GlobalListCL::FinalizeModify:"
        "Inconsistent Data_."), DebugContainerC );
    modifiable_= false;
    if (LevelViews_.empty()) return;
    LevelStarts_.resize( LevelViews_.size() + 1);
    LevelStarts_[LevelViews_.size()]= Data_.end();
    const_LevelStarts_.resize( LevelViews_.size() + 1);
    const_LevelStarts_[LevelViews_.size()]= Data_.end();
    for (Uint lvl= LevelViews_.size(); lvl > 0; --lvl) {
        Data_.splice( Data_.begin(), *LevelViews_[lvl-1]);
        LevelStarts_[lvl-1]= Data_.begin();
        const_LevelStarts_[lvl-1]= Data_.begin();
        Assert( LevelViews_[lvl-1]->empty(),
            DROPSErrCL("GlobalListCL::FinalizeModify: Did not move all Data."),
            DebugContainerC );
        delete LevelViews_[lvl-1];
    }
    LevelViews_.clear();
}

template <class T>
  void
  GlobalListCL<T>::AppendLevel()
{
    Assert( modifiable_, DROPSErrCL("GlobalListCL::AppendLevel: "
        "Data not modifiable."), DebugContainerC );
    LevelViews_.push_back( new LevelCont());
}

template <class T>
  void
  GlobalListCL<T>::RemoveLastLevel()
{
    Assert( modifiable_, DROPSErrCL("GlobalListCL::RemoveLast: "
        "Data not modifiable."), DebugContainerC );
    Assert( LevelViews_.size() > 0, DROPSErrCL("GlobalListCL: RemoveLastLevel: "
        "There are no levels to be removed."), DebugContainerC);
    Assert( LevelViews_.back()->empty(), DROPSErrCL("GlobalListCL: RemoveLastLevel: "
        "Last level not empty"), DebugContainerC);
    delete LevelViews_.back();
    LevelViews_.pop_back();
}

//**************************************************************************
// Class:   MLDataCL                                                       *
//**************************************************************************
template <class T>
class MLDataCL : public std::list<T>
{
  public:
    explicit MLDataCL ()
        : std::list<T>() {}
    explicit MLDataCL (size_t n, const T& val= T())
        : std::list<T>( n, val) {}

    T&       GetFinest()      { return this->back();  }
    T&       GetCoarsest()    { return this->front(); }
    T*       GetFinestPtr()   { return &this->back(); }
    T*       GetCoarsestPtr() { return &this->front(); }
    const T& GetFinest()      const { return this->back();  }
    const T& GetCoarsest()    const { return this->front(); }
    const T* GetFinestPtr()   const { return &this->back(); }
    const T* GetCoarsestPtr() const { return &this->front(); }
    typename MLDataCL::iterator       GetFinestIter()         { return --this->end(); }
    typename MLDataCL::const_iterator GetFinestIter()   const { return --this->end(); }
    typename MLDataCL::iterator       GetCoarsestIter()       { return this->begin(); }
    typename MLDataCL::const_iterator GetCoarsestIter() const { return this->begin(); }
};

///\brief Designates the part of the domain, usually on tetras at the interface, one is interested in.
enum TetraSignEnum { AllTetraC, NegTetraC, PosTetraC };

} // end of namespace DROPS

#endif
