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

class ComplementBodyCL : public BodyCL
{
  private:
    const BodyCL* b_;

  public:
    ComplementBodyCL (BodyStackT& s, const std::string& /*instruction*/)
        { b_= pop_and_take_ownership( s); }

    double operator() (const Point3DCL& x, double t) const
        { return -(*b_)( x, t); }
};

class IntersectionBodyCL : public BodyCL
{
  private:
    const BodyCL* b0_;
    const BodyCL* b1_;
    bool  use_smooth_op_;

  public:
    IntersectionBodyCL (BodyStackT& s, const std::string& instruction) {
         use_smooth_op_= instruction == std::string( "SmoothIntersection");
         b0_= pop_and_take_ownership( s);
         b1_= pop_and_take_ownership( s);
    }

    double operator() (const Point3DCL& x, double t) const {
        const double v0= (*b0_)( x, t), v1= (*b1_)( x, t);
        return use_smooth_op_ ? v0 + v1 + hypot(v0, v1) : std::max( v0, v1);
    }
};

class UnionBodyCL : public BodyCL
{
  private:
    const BodyCL* b0_;
    const BodyCL* b1_;
    bool  use_smooth_op_;

  public:
    UnionBodyCL (BodyStackT& s, const std::string& instruction) {
         use_smooth_op_= instruction == std::string( "SmoothUnion");
         b0_= pop_and_take_ownership( s);
         b1_= pop_and_take_ownership( s);
    }

    double operator() (const Point3DCL& x, double t) const {
        const double v0= (*b0_)( x, t), v1= (*b1_)( x, t);
        return use_smooth_op_ ? v0 + v1 - hypot(v0, v1) : std::min( v0, v1);
    }
};

/// n_^Tx + a_
class HalfspaceBodyCL : public BodyCL
{
  private:
    Point3DCL n_;
    double a_;

  public:
    HalfspaceBodyCL (BodyStackT&, const ParamCL& p)
        : n_( p.get<Point3DCL>( "Normal")), a_( 0.) {
        try {
            a_= p.get<double>( "ValueAtOrigin");
        } catch (DROPSErrCL) {}
        bool normalize= false;
        try {
            normalize= p.get<bool>( "MakeUnitNormal");
        } catch (DROPSErrCL) {}
        if (normalize == true)
            n_/= n_.norm();
    }

    const VecOfStringT& known_keys () const {
        static std::string keys[] = {
            std::string( "Type"),
            std::string( "Normal"),
            std::string( "ValueAtOrigin"),
            std::string( "MakeUnitNormal") };
        static const VecOfStringT k( keys, keys + 4);
        return k;
    }

    double operator() (const Point3DCL& x, double) const
        { return inner_prod( n_, x) + a_; }
};

/// \|x -c_\| - r_
class SphereBodyCL : public BodyCL
{
  private:
    Point3DCL c_;
    double r_;

  public:
    SphereBodyCL (BodyStackT&, const ParamCL& p)
        : r_( 0.) {
        try {
            r_= p.get<double>( "Radius");
        } catch (DROPSErrCL) {}
        try {
            c_= p.get<Point3DCL>( "Center");
        } catch (DROPSErrCL) {}
    }

    const VecOfStringT& known_keys () const {
        static std::string keys[] = {
            std::string( "Type"),
            std::string( "Radius"),
            std::string( "Center")};
        static const VecOfStringT k( keys, keys + 3);
        return k;
    }

    double operator() (const Point3DCL& x, double) const
        { return (x - c_).norm() - r_; }
};

class InfiniteConeBodyCL : public BodyCL
{
  private:
    Point3DCL apex_,
              axis_;
    double c_; // cos(aperture/2)

  public:
    InfiniteConeBodyCL (BodyStackT&, const ParamCL& p)
        : axis_( std_basis<3>( 1)) {
        double a= M_PI/4.;
        try {
            a= p.get<double>( "SemiAperture");
        } catch (DROPSErrCL) {}
        c_= std::cos( a);
        try {
            apex_= p.get<Point3DCL>( "Apex");
        } catch (DROPSErrCL) {}
        try {
            axis_= p.get<Point3DCL>( "Axis");
        } catch (DROPSErrCL) {}
        axis_/=axis_.norm();
    }
 
   const VecOfStringT& known_keys () const {
        static std::string keys[] = {
            std::string( "Type"),
            std::string( "SemiAperture"),
            std::string( "Apex"),
            std::string( "Axis") };
        static const VecOfStringT k( keys, keys + 4);
        return k;
    }

    double operator() (const Point3DCL& x, double) const
        { return  c_*(x - apex_).norm() - inner_prod(x - apex_, axis_); }
};

class InfiniteCylinderBodyCL : public BodyCL
{
  private:
    Point3DCL origin_,
              axis_;
    double r_;

  public:
    InfiniteCylinderBodyCL (BodyStackT&, const ParamCL& p)
        : axis_( std_basis<3>( 1)), r_( 1.0) {
        try {
            r_= p.get<double>( "Radius");
        } catch (DROPSErrCL) {}
        try {
            origin_= p.get<Point3DCL>( "Origin");
        } catch (DROPSErrCL) {}
        try {
            axis_= p.get<Point3DCL>( "Axis");
        } catch (DROPSErrCL) {}
        axis_/=axis_.norm();
    }

    const VecOfStringT& known_keys () const {
        static std::string keys[] = {
            std::string( "Type"),
            std::string( "Radius"),
            std::string( "Origin"),
            std::string( "Axis") };
        static const VecOfStringT k( keys, keys + 4);
        return k;
    }

    double operator() (const Point3DCL& x, double) const
        { return  (x - origin_ - inner_prod(x - origin_, axis_)*axis_).norm() - r_; }
};


class TorusBodyCL : public BodyCL
{
  private:
    Point3DCL origin_,
              axis_;
    double R_,
           r_;

  public:
    TorusBodyCL (BodyStackT&, const ParamCL& p)
        : origin_( std_basis<3>( 0)), axis_( std_basis<3>( 1)), R_( 0.75), r_( 0.25) {
        try {
            R_= p.get<double>( "BigRadius");
        } catch (DROPSErrCL) {}
        try {
            r_= p.get<double>( "SmallRadius");
        } catch (DROPSErrCL) {}
        try {
            origin_= p.get<Point3DCL>( "Origin");
        } catch (DROPSErrCL) {}
        try {
            axis_= p.get<Point3DCL>( "Axis");
        } catch (DROPSErrCL) {}
        axis_/=axis_.norm();
    }

    const VecOfStringT& known_keys () const {
        static std::string keys[] = {
            std::string( "Type"),
            std::string( "BigRadius"),
            std::string( "SmallRadius"),
            std::string( "Origin"),
            std::string( "Axis") };
        static const VecOfStringT k( keys, keys + 5);
        return k;
    }

    double operator() (const Point3DCL& x, double) const {
        const double  z= inner_prod( axis_, x - origin_),
                     xy= std::sqrt( (x - origin_).norm_sq() - z*z);
        return  std::sqrt( z*z + std::pow( xy - R_, 2)) - r_; }
};

/// Levelset
class LevelsetBodyCL : public BodyCL
{
  private:
    double (*ls_)(const Point3DCL&, double);

  public:
    LevelsetBodyCL (BodyStackT&, const ParamCL& p)
        : ls_( InScaMap::getInstance()[p.get<std::string>( "Function")]) {}

    const VecOfStringT& known_keys () const {
        static std::string keys[] = {
            std::string( "Type"),
            std::string( "Function")};
        static const VecOfStringT k( keys, keys + 2);
        return k;
    }

    double operator() (const Point3DCL& x, double t) const
        { return (*ls_)( x, t); }
};

///\brief Load a body from another json-file.
/// It is constructed by a call to body_builder.
class ModuleBodyCL : public BodyCL
{
  private:
    const BodyCL* b_;

  public:
    ModuleBodyCL (BodyStackT& s, const ParamCL& p) {
        std::ifstream is( p.get<std::string>( "Path").c_str());
        ParamCL pm;
        is >> pm;
        std::string name;
        try {
            name= p.get<std::string>( "Name");
        } catch (DROPSErrCL) {}
        s.push( name == std::string( "") ? body_builder( pm)
                                         : body_builder( pm.get_child( name)));
        b_= pop_and_take_ownership( s);
    }

    const VecOfStringT& known_keys () const {
        static std::string keys[] = {
            std::string( "Type"),
            std::string( "Path"),
            std::string( "Name")};
        static const VecOfStringT k( keys, keys + 3);
        return k;
    }

    double operator() (const Point3DCL& x, double t) const
        { return (*b_)( x, t); }
};

class SimilarityTransformBodyCL : public BodyCL
{
  private:
    const BodyCL* b_;

    double s_;
    bool do_rotate_;
    SMatrixCL<3,3> r_;
    Point3DCL t_;

  public:
    SimilarityTransformBodyCL (BodyStackT& s, const ParamCL& p) 
        : s_(1.), do_rotate_( false), t_( 0.) {
        try {
            s_= p.get<double>( "Scaling");
        } catch (DROPSErrCL) {}
        bool haveangle= false;
        double a;
        try {
            a= p.get<double>( "RotationAngle");
            haveangle= true;
        } catch (DROPSErrCL) {}
        bool haveaxis= false;
        Point3DCL axis;
        try {
            axis= p.get<Point3DCL>( "RotationAxis");
            haveaxis= true;
        } catch (DROPSErrCL) {}
        if (haveaxis != haveangle)
            throw DROPSErrCL( "SimilarityTransformBodyCL: A rotation must have a \"RotationAxis\" and a \"RotationAngle\".\n");
        if (haveangle) { // Precompute the rotation-matrix.
            axis/=axis.norm();
            SMatrixCL<3,3> o;
            o(1,0)= -(o(0,1)= -axis[2]); // Matrix-rep of cross-product with axis.
            o(2,0)= -(o(0,2)=  axis[1]);
            o(2,1)= -(o(1,2)= -axis[0]);
            r_(0,0)= r_(1,1)= r_(2,2)= 1.;
            r_+= std::sin( a)*o + (1. - std::cos(a))*(o*o);
            do_rotate_= true;
        }
        try {
            t_= p.get<Point3DCL>( "Translation");
        } catch (DROPSErrCL) {}
        b_= pop_and_take_ownership( s);
    }

    const VecOfStringT& known_keys () const {
        static std::string keys[] = {
            std::string( "Type"),
            std::string( "Scaling"),
            std::string( "RotationAngle"),
            std::string( "RotationAxis"),
            std::string( "Translation")  };
        static const VecOfStringT k( keys, keys + 5);
        return k;
    }

    double operator() (const Point3DCL& x, double t) const {
        Point3DCL y= s_*x;
        if (do_rotate_)
            y= r_*y;
        y+= t_;
        return (*b_)( y, t)/s_; }
};

///\brief Refer to an already existing body. Therefore, do not take ownership of it.
/// This is an extension of a stack-machine.
/// By allowing references, one can construct a DAG instead of a tree. This
/// allows for smaller representations of complex geometries, as copies of
/// existing branches can be avoided.
class ReferenceBodyCL : public BodyCL
{
  private:
    const BodyCL* b_;

  public:
    ReferenceBodyCL (const BodyCL* b) : b_( b) {}

    double operator() (const Point3DCL& x, double t) const
        { return (*b_)( x, t); }
};

/// \brief Free-standing builder for a BodyCL-specialization.
/// Pops its arguments from s and pushes its product on s.
typedef void (*SimpleBodyBuilderT)(BodyStackT& s, const std::string& instruction);
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
make_new (BodyStackT& s, const std::string& instruction)
{ s.push( new T( s, instruction)); }

template <class T>
inline void
make_new (BodyStackT& s, const ParamCL& p)
{ s.push( new T( s, p)); }
///@}

// Note how C++ selects the right overload of make_new by the type of the MapRegisterCL<T>.
RegisterSimpleBuilder sreg00("Complement",         &make_new<ComplementBodyCL>);
RegisterSimpleBuilder sreg01("Intersection",       &make_new<IntersectionBodyCL>);
RegisterSimpleBuilder sreg02("Union",              &make_new<UnionBodyCL>);
RegisterSimpleBuilder sreg03("SmoothIntersection", &make_new<IntersectionBodyCL>);
RegisterSimpleBuilder sreg04("SmoothUnion",        &make_new<UnionBodyCL>);

RegisterBuilder reg00("Halfspace",               &make_new<HalfspaceBodyCL>);
RegisterBuilder reg01("Sphere",                  &make_new<SphereBodyCL>);
RegisterBuilder reg02("InfiniteCone",            &make_new<InfiniteConeBodyCL>);
RegisterBuilder reg03("InfiniteCylinder",        &make_new<InfiniteCylinderBodyCL>);
RegisterBuilder reg04("Torus",                   &make_new<TorusBodyCL>);
RegisterBuilder reg05("Levelset",                &make_new<LevelsetBodyCL>);
RegisterBuilder reg06("LoadFromModule",          &make_new<ModuleBodyCL>);
RegisterBuilder reg07("ApplySimilarityToDomain", &make_new<SimilarityTransformBodyCL>);


std::string unknown_keys (const ParamCL& p, const VecOfStringT& known_keys)
{
    std::string uk;
    using boost::property_tree::ptree;
    for (ptree::const_iterator i= p.begin(); i != p.end(); ++i)
        if (std::find( known_keys.begin(), known_keys.end(), i->first) == known_keys.end())
            uk+= (uk != std::string() ? std::string( " ") : std::string()) + i->first;
    return uk;
}

const BodyCL* body_builder (const ParamCL& p) // :-)
{
    BodyStackT s;
    std::map<std::string, const BodyCL*> named_refs;

    if (p.get<std::string>( "Type") != std::string( "CSG-body (v0)"))
        throw DROPSErrCL( "body_builder: Unknown CSG-specification.\n");

    size_t ic= 0; // Instruction counter for debugging.
    using boost::property_tree::ptree;
    try {
        for (ptree::const_iterator i= p.get_child( "Instructions").begin(); i != p.get_child( "Instructions").end(); ++i, ++ic) {
            std::string op= i->second.get_value<std::string>(); // Operations without parameters (other than s) are represented by bare strings.
            if (!op.empty())
                SimpleBuilderMap::getInstance()[op]( s, op);
            else {
                op= i->second.get<std::string>( "Type");
                if      (op == std::string( "CreateReference"))
                    named_refs[i->second.get<std::string>( "Name")]= s.top();
                else if (op == std::string( "PushReference"))
                    s.push( new ReferenceBodyCL( named_refs[i->second.get<std::string>( "Name")]));
                else {
                    BuilderMap::getInstance()[op]( s, i->second);
                    std::string uk= unknown_keys( i->second, s.top()->known_keys());
                    if (uk != std::string())
                        std::cerr << "body_builder: Warning: Ignoring unknown key(s) \"" << uk
                                  << "\" while proccessing instruction " << ic << ".\n";
                }
            }
        }
    } catch (...) {
        std::cerr << "body_builder: Error while proccessing instruction " << ic << ".\n";
        throw;
    }

    const BodyCL* ret= s.top();
    if (s.size() != 1) {
        delete ret;
        throw DROPSErrCL( "body_builder: Stack not empty after popping the final product.\n");
    }
    return ret;
}

} // end of namespace DROPS::CSG
} // end of namespace DROPS
