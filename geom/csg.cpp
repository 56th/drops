/// \file csg.cpp
/// \brief classes for simple constructive solid geometry (CSG) for level sets
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
 * Copyright 2012 LNM/SC RWTH Aachen, Germany
*/

#include "geom/csg.h"
#include "misc/funcmap.h"

namespace DROPS {
namespace CSG {

/// \brief Stack-type for the stack-machine.
typedef std::stack<const BodyCL*>  BodyStackT;


class ComplementBodyCL : public BodyCL
{
  private:
    const BodyCL* b_;

  public:
    ComplementBodyCL (BodyStackT& s) {
        if (s.empty())
            throw DROPSErrCL("ComplementBodyCL: Nothing to complement.\n");
        b_= s.top();
        s.pop();
    }

    double operator() (const Point3DCL& x, double t) const
        { return -(*b_)( x, t); }
};

class IntersectionBodyCL : public BodyCL
{
  private:
    const BodyCL* b0_;
    const BodyCL* b1_;

  public:
    IntersectionBodyCL (BodyStackT& s) {
        if (s.size() < 2)
            throw DROPSErrCL("IntersectionBodyCL: Less than 2 bodies on the stack.\n");
         b0_= s.top(); s.pop();
         b1_= s.top(); s.pop();
    }

    double operator() (const Point3DCL& x, double t) const {
        const double v0= (*b0_)( x, t), v1= (*b1_)( x, t);
        return std::max( v0, v1);
    }
};

class UnionBodyCL : public BodyCL
{
  private:
    const BodyCL* b0_;
    const BodyCL* b1_;

  public:
    UnionBodyCL (BodyStackT& s) {
        if (s.size() < 2)
            throw DROPSErrCL("UnionBodyCL: Less than 2 bodies on the stack.\n");
         b0_= s.top(); s.pop();
         b1_= s.top(); s.pop();
    }

    double operator() (const Point3DCL& x, double t) const {
        const double v0= (*b0_)( x, t), v1= (*b1_)( x, t);
        return std::min( v0, v1);
    }
};

/// n_^Tx + a_
/// Normalizes n.
class HalfspaceBodyCL : public BodyCL
{
  private:
    Point3DCL n_;
    const double a_;

  public:
    HalfspaceBodyCL (BodyStackT&, const ParamCL& p)
        : n_( p.get<Point3DCL>( "Normal")), a_( p.get<double>( "ValueAtOrigin")) {
        bool normalize= false;
        try {
            normalize= p.get<bool>( "MakeUnitNormal");
        } catch (DROPSErrCL) {}
        if (normalize == true)
            n_/= n_.norm();
    }

    double operator() (const Point3DCL& x, double) const {
        return inner_prod( n_, x) + a_;
    }
};

/// \|x -c_\| - r_
class SphereBodyCL : public BodyCL
{
  private:
    const Point3DCL c_;
    const double r_;

  public:
    SphereBodyCL (BodyStackT&, const ParamCL& p)
        : c_( p.get<Point3DCL>( "Center")), r_( p.get<double>( "Radius")) {}

    double operator() (const Point3DCL& x, double) const {
        return (x - c_).norm() - r_;
    }
};

/// Levelset
class LevelsetBodyCL : public BodyCL
{
  private:
    double (*ls_)(const Point3DCL&, double);

  public:
    LevelsetBodyCL (BodyStackT& s, const ParamCL& p)
        : ls_( InScaMap::getInstance()[p.get<std::string>( "Function")]) {}

    double operator() (const Point3DCL& x, double t) const {
        return (*ls_)( x, t);
    }
};


/// \brief Free-standing builder for a BodyCL-specialization.
/// Pops its arguments from s and pushes its product on s.
typedef void (*SimpleBodyBuilderT)(BodyStackT& s);
typedef SingletonMapCL<SimpleBodyBuilderT> SimpleBuilderMap;
typedef MapRegisterCL<SimpleBodyBuilderT>  RegisterSimpleBuilder;

/// \brief Free-standing builder for a BodyCL-specialization.
/// Pops its arguments from s and pushes its product on s. Additional arguments are in ParamCL p.
typedef void (*BodyBuilderT)(BodyStackT& s, const ParamCL& p);
typedef SingletonMapCL<BodyBuilderT> BuilderMap;
typedef MapRegisterCL<BodyBuilderT>  RegisterBuilder;

///\brief Templates for SimpleBodyBuilderT and BodyBuilderT
///@{
template <class T>
inline void
simple_make_new (BodyStackT& s)
{ s.push( new T( s)); }

template <class T>
inline void
make_new (BodyStackT& s, const ParamCL& p)
{ s.push( new T( s, p)); }
///@}

RegisterSimpleBuilder sreg00("Complement",   &simple_make_new<ComplementBodyCL>);
RegisterSimpleBuilder sreg01("Intersection", &simple_make_new<IntersectionBodyCL>);
RegisterSimpleBuilder sreg02("Union",        &simple_make_new<UnionBodyCL>);

RegisterBuilder reg00("Halfspace", &make_new<HalfspaceBodyCL>);
RegisterBuilder reg01("Sphere",    &make_new<SphereBodyCL>);
RegisterBuilder reg02("Levelset",  &make_new<LevelsetBodyCL>);


const BodyCL* body_builder (const ParamCL& p) // :-)
{
    BodyStackT s;
    BodyVectorT refs;

    if (p.get<std::string>( "Type") != std::string( "CSG-body (v0)"))
        throw DROPSErrCL( "body_builder: Unknown CSG-specification.\n");

    using boost::property_tree::ptree;
    for (ptree::const_iterator i= p.get_child( "Instructions").begin(); i != p.get_child( "Instructions").end(); ++i) {
        std::string op= i->second.get_value<std::string>(); // Operations without parameters (other than S) are represented by bare strings.
        if (!op.empty())
            SimpleBuilderMap::getInstance()[op]( s);
        else {
            op= i->second.get<std::string>( "Type");
            BuilderMap::getInstance()[op]( s, i->second);
        }
        refs.push_back( s.top());
    }

    refs.pop_back(); // If ret owns a pointer to itself, it will delete itself again...
    const BodyCL* ret= s.top();
    ret->make_owner_of( refs);
    if (s.size() != 1) {
        delete ret;
        throw DROPSErrCL( "body_builder: Stack not empty after popping the final product.\n");
    }
    return ret;
}

} // end of namespace DROPS::CSG
} // end of namespace DROPS