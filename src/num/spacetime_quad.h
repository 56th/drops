/// \file spacetime_quad.h
/// \brief helper functions to construct integration rules on (cutted) 4d-geometries (consisting of pentatopes)
/// \author LNM RWTH Aachen: Christoph Lehrenfeld

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

#ifndef DROPS_SPACETIME_QUAD_H
#define DROPS_SPACETIME_QUAD_H

#include "num/spacetime_geom.h"
#include "num/spacetime_map.h"
#include "num/discretize.h"
#include "misc/scopetimer.h"

namespace DROPS
{


typedef double (*instat_scalar_fun_ptr)(const DROPS::Point3DCL&, double);

template <class IntRuleData,class Geom>
void Gather4DIntegrationPoints(const std::vector<Geom> &, GridFunctionCL<Point4DCL> &);

template <class IntRuleData,class Geom>
void Gather4DIntegrationWeights(const std::vector<Geom> &, GridFunctionCL<double> &);

template <class IntRuleData>
void Gather4DNormals(const std::vector<Tetra4DCL> &, GridFunctionCL<Point4DCL> &);

template <class IntRuleData>
void Gather4DNu(const std::vector<Tetra4DCL> &, GridFunctionCL<double> &);

template <class IntRuleData,class Geom>
SArrayCL<double,IntRuleData::NumNodesC> Transform4DIntegrationWeights(const Geom &);

template <class IntRuleData,class Geom>
SArrayCL<Point4DCL,IntRuleData::NumNodesC> Transform4DIntegrationPoints(const Geom &);




// class decomposing given prism, evaluating level set, decomposing into
// decomposition of uncut geometries and afterwards gathering all the 
// data for integration on this geometry, including quadrature rules
// template argument is a 4D quadrature rule which is applied on each 
// uncut pentatope after the decomposition phase.
// Most work (decomposition, generating quad. rule) is done in Initialize
template <class IntRuleData4D>
class CompositeSTQuadCL
{
protected:
    Point4DContainer pcont;
    GridFunctionCL<Point4DCL> ips_neg;//int. point in neg. region
    GridFunctionCL<double> ipw_neg;   //int. point weights in neg. region
    GridFunctionCL<Point4DCL> ips_pos;//int. point in pos. region
    GridFunctionCL<double> ipw_pos;   //int. point weights in pos. region
    GridFunctionCL<Point4DCL> ips_if; //int. point at interface
    GridFunctionCL<double> ipw_if;    //int. point weights at interface
    GridFunctionCL<Point4DCL> ipn_if; //normals
    GridFunctionCL<Point3DCL> ipsn_if; //spacenormals
    GridFunctionCL<double> ipnu_if;   //nu values
    std::vector<PentatopeCL> negpentas;
    std::vector<PentatopeCL> pospentas;
    std::vector<Tetra4DCL> iftetras;
    bool  hasinterface;

    GridFunctionCL<double> Eval ( instat_scalar_fun_ptr f, const GridFunctionCL<Point4DCL> & points,
                                  const SpaceTimeMapping * map = &SpaceTimeIdentity::getInstance()) const;
    GridFunctionCL<Point3DCL> Eval ( instat_vector_fun_ptr f, const GridFunctionCL<Point4DCL> & points,
                                  const SpaceTimeMapping * map = &SpaceTimeIdentity::getInstance()) const;
    template <class T>
    GridFunctionCL<T> EvalLinear ( const LocalP2CL<T> & fold, const LocalP2CL<T> & fnew, 
                                        const GridFunctionCL<Point4DCL> & points) const;
    template <class T> 
    T Quad ( const GridFunctionCL<T> & f, const GridFunctionCL<double> & weights) const;

    std::valarray<bool> p1p1signs; //signs for P1(space)-P1(time)-ext. finite elements 
    // [TODO: find a place to put that into]

public:
    void Initialize(const TetraCL* tet, const TimeInterval* ti, 
                    const GeneralizedPrism4CL* refprism4,
                    const LocalP2CL<double>* lsetold, const LocalP2CL<double>* lsetnew, 
                    instat_scalar_fun_ptr f, Uint ints_per_space_edge = 2, Uint subtimeintervals = 1);

    CompositeSTQuadCL(const TetraCL& tet, const TimeInterval& ti, instat_scalar_fun_ptr f, 
                      Uint ints_per_space_edge = 2, Uint subtimeintervals = 1);

    CompositeSTQuadCL(const TetraCL& tet, const TimeInterval& ti, 
                      const LocalP2CL<double>& lsetold, const LocalP2CL<double>& lsetnew,
                      Uint ints_per_space_edge = 2, Uint subtimeintervals = 1);

    CompositeSTQuadCL(const GeneralizedPrism4CL& refprism4, instat_scalar_fun_ptr f, 
                      Uint ints_per_space_edge = 2, Uint subtimeintervals = 1);

    CompositeSTQuadCL(const GeneralizedPrism4CL& refprism4,
                      const LocalP2CL<double>& lsetold, const LocalP2CL<double>& lsetnew,
                      Uint ints_per_space_edge = 2, Uint subtimeintervals = 1);

    virtual ~CompositeSTQuadCL(){ };

    const GridFunctionCL<Point4DCL> & GetVolumeIntegrationPoints ( bool posPart) const;
    const GridFunctionCL<double> & GetVolumeIntegrationWeights ( bool posPart) const;
    Uint NumberOfIntegrationPoints ( bool posPart) const;
    const GridFunctionCL<Point4DCL> & GetInterfaceIntegrationPoints ( ) const;
    const GridFunctionCL<double> & GetInterfaceIntegrationWeights ( ) const;
    Uint NumberOfInterfaceIntegrationPoints ( ) const;

    GridFunctionCL<double> EvalOnPart ( instat_scalar_fun_ptr f, bool posPart,
                                        const SpaceTimeMapping * map = &SpaceTimeIdentity::getInstance()) const;
    GridFunctionCL<Point3DCL> EvalOnPart ( instat_vector_fun_ptr f, bool posPart,
                                        const SpaceTimeMapping * map = &SpaceTimeIdentity::getInstance()) const;
    GridFunctionCL<double> EvalOnInterface ( instat_scalar_fun_ptr f,
                                             const SpaceTimeMapping * map = &SpaceTimeIdentity::getInstance()) const;
    template <class T>
    GridFunctionCL<T> EvalLinearOnPart ( const LocalP2CL<T>& fold, 
                                         const LocalP2CL<T>& fnew, bool posPart) const;

    template <class T>
    GridFunctionCL<T> EvalLinearOnInterface ( const LocalP2CL<T>& fold, 
                                              const LocalP2CL<T>& fnew) const;
    const GridFunctionCL<Point4DCL> & GetNormalsOnInterface ( ) const;
    const GridFunctionCL<Point3DCL> & GetSpaceNormalsOnInterface ( ) const;
    const GridFunctionCL<double> & GetNuOnInterface ( ) const;
    bool HasInterface ( ) const;
    const std::valarray<bool> & GetVertexSigns () const; //only needed for P1P1-XFEM so far
                                                         //(could be removed if p1p1signs is removed)
    bool GetVertexSign (Uint i, bool newtime) const; //only needed for P1P1-XFEM so far
                                                     //(could be removed if p1p1signs is removed)

    template <class T> 
    T QuadOnPart ( const GridFunctionCL<T> & f, bool posPart) const;

    template <class T> 
    T QuadOnInterface ( const GridFunctionCL<T> & f) const;

    void Report(std::ostream & out, std::string linehead="", std::string linetail="") const;

};


template <class VolRule, class SurfRule>
struct Vol_Surf_DataCL{
    typedef VolRule Volume;
    typedef SurfRule Surface;
};

typedef Vol_Surf_DataCL<Quad3_4DDataCL,Quad3DataCL> Quad3_4D_and_3D_DataCL;
typedef Vol_Surf_DataCL<Quad5_4DDataCL,Quad5DataCL> Quad5_4D_and_3D_DataCL;

typedef Quad3_4D_and_3D_DataCL QuadRule;

// helper functions for pointwise operations on GridFunctions:

// e.g. for debugging
inline Point4DCL transform_points_in_space( const Point3DCL & offset, const SMatrixCL<3,3> & T, const Point4DCL & p)
{
    Point4DCL ret(p);
    for (Uint i = 0; i < 3; ++i)
    {
        ret[i] = offset[i];
        for (Uint j = 0; j < 3; ++j)
            ret[i] += T(i,j) * p[j];
    }
    return ret;
}

inline Point4DCL transform_points_in_spacetime( const Point3DCL & offset, const SMatrixCL<3,3> & T, const TimeInterval& ti,
                                                const Point4DCL & p)
{
    Point4DCL ret(p);
    for (Uint i = 0; i < 3; ++i)
    {
        ret[i] = offset[i];
        for (Uint j = 0; j < 3; ++j)
            ret[i] += T(i,j) * p[j];
    }
    ret[3] = ti.first + (ti.second-ti.first)*ret[3];
    return ret;
}

// e.g. for debugging
inline Point4DCL transform_points_in_time( const double t0, const double t1, const Point4DCL & p)
{
    Point4DCL ret(p);
    ret[3] = t0 + (t1-t0) * p[3];
    return ret;
}


// transformation of normals: n = |F| * F^-T * n_ref with F the jacobi of the mapping from ref. to phys. and n_ref
// the space time normal on the reference domain. Note that ||n|| is the measure ratio due to the transformation
// and especially it (||n||) is not 1.0!
inline Point4DCL transform_normals_in_space( const SMatrixCL<3,3> & T, const double & absdet, const Point4DCL & p)
{
    Point4DCL ret(0.0);
    for (Uint i = 0; i < 3; ++i)
    {
        for (Uint j = 0; j < 3; ++j)
            ret[i] += T(i,j) * p[j] * absdet;
    }
    ret[3] = p[3];
    return ret;
}

inline Point4DCL transform_normals_in_spacetime( const SMatrixCL<3,3> & T, const double & absdet, const double & dt, const Point4DCL & p)
{
    Point4DCL ret(0.0);
    for (Uint i = 0; i < 3; ++i)
    {
        for (Uint j = 0; j < 3; ++j)
            ret[i] += T(i,j) * p[j] * absdet * dt;
    }
    ret[3] = p[3] * absdet;
    return ret;
}


inline double calc_norm( const Point4DCL & p)
{
    return p.norm();
}

inline double calc_norm_sqr( const Point3DCL & p)
{
    const double a = p.norm();
    return a*a;
}

inline double calc_norm_and_normalize( Point4DCL & p)
{
    const double n = p.norm();
    p /= n;
    return n;
}

inline void normalize( Point3DCL & p)
{
    const double n = p.norm();
    p /= n;
}


inline double calc_nu_from_stnormal( const Point4DCL & p)
{
    const double nu = std::sqrt( p[0]*p[0]+p[1]*p[1]+p[2]*p[2] );
    return nu;
}

/*
  combination of calc_norm_and_normalize and calc_nu_from_normal
 */
inline double calc_nu_measure_and_normalize( Point4DCL & p)
{
    const double n = calc_norm_and_normalize(p);
    const double nu = calc_nu_from_stnormal(p);
    return n * nu;
}

inline Point3DCL calc_spacenormal_from_stnormal(const Point4DCL & p)
{
    Point3DCL ret(RestrictToPoint3D(p));
    const double nu = ret.norm();
    if (nu >= 1e-12)
        ret /= nu;
    else
    {
        throw DROPSErrCL("nu==0.0");
        ret = MakePoint3D(0.0, 0.0, 0.0);
    }
    return ret;
}


inline double square(const double& a)
{
    return a * a;
}





} // end of namespace DROPS

#endif
