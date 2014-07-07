/// \file spacetime_geom.h
/// \brief classes that constitute space time geometries - 
/// Note that in contrast to VertexCL, ..., SimplexCL, ... these geometries only appear element-local 
/// and their purpose is the decomposition of a cut space-time primitive into several uncut ones
/// main components:
///  * Point4DContainer (gathered collection of unique Points)
///  * geometric primitive classes:
///    * Tetra4DCL
///    * HyperTrigCL
///    * GeneralizedPrism4CL
///    * PentatopeCL
/// PentatopeCL basically does all the work. With the help of levelset values a cut Pentatope
/// can be decomposed into uncut Pentatopes, HyperTrigs and GeneralizedPrisms. Each class can
/// decompose itself into Pentatopes. The Pentatope class also decomposes the interface into 
/// Tetra4Ds. 
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

#ifndef DROPS_SPACETIME_GEOM_H
#define DROPS_SPACETIME_GEOM_H

#include <set>
#include <vector>
#include <iostream>
#include "misc/utils.h"
#include "geom/boundary.h"
#include "geom/topo.h"
#include "geom/simplex.h"
#include "geom/principallattice.h"

namespace DROPS
{

typedef double    (*instat_scalar_fun_ptr)(const DROPS::Point3DCL&, double);
typedef std::pair<double,double> TimeInterval;

/// struct which defines the relation a < b for Point4DCL (in order to use std::set-features)
struct Point4Dless {
  bool operator() (const Point4DCL& a, const Point4DCL& b) const {
    // if you want to merge point which are "the same up 14 digits..."
    const double EPS = 0.0; // 1e-14*abs(b[i])
    for (Uint i=0; i<4; i++)
    {
        if (a[i] < b[i] - EPS)
            return true;
        if (a[i] > b[i] + EPS)
            return false;
    }
    return false;
  }
  typedef Point4DCL first_argument_type;
  typedef Point4DCL second_argument_type;
  typedef bool result_type;
};


/// Container set constitutes a collection of Point4Ds 
/// main feature: the operator()(const Point4DCL & p)
/// The points in the container are owned and later 
/// released by Point4DContainer
class Point4DContainer
{
protected:
    std::set<Point4DCL,Point4Dless> pset;
#ifdef DEBUG
    size_t k;
#endif
public: 
    Point4DContainer()
    {
#ifdef DEBUG
        k=0;
#endif
        pset.clear();
    };
    
    /// Access operator to points
    /// Either point is already in the Container, 
    ///   then return pointer to that point
    /// or point is not in the Container yet,
    ///   then add point to container and return pointer to new Point4D
    /// The return value (pointer ) points to a Point4D which is owned
    /// and later released by Point4DContainer
    const Point4DCL* operator()(const Point4DCL & p)
    {
        std::set<Point4DCL,Point4Dless>::iterator it;
        it = pset.find(p);
        if (it == pset.end())
        {
            pset.insert(p);
            return &(*pset.find(p));
        }
        else
        {
#ifdef DEBUG
            k++;
#endif
            return &(*it);
        }
    }

    void Report(std::ostream & out) const
    {
        out << " Point4DContainer stored " << pset.size() << " points.\n";
#ifdef DEBUG
        out << " Point4DContainer rejected " << k << " points.\n";
#endif
    }

    ~Point4DContainer(){};
};

// Helper class constituting a four dimensional tetrahedra T (with meas_4(T) = 0)
class Tetra4DCL
{
protected:
    SArrayCL<const Point4DCL*,4> p;

    /// The calculation of absdet or the normal is only done once
    /// therefore the bool .._initialized checks if that has been done already or not
    /// the calculation is done in calc_normal() or calc_abs_determinant()
    mutable double absdet_of_trafo;
    mutable bool absdet_of_trafo_initialized;
    mutable Point4DCL normal;
    mutable bool normal_initialized;

    // helper members (e.g. to unify the orientation of the normal)
    mutable const Point4DCL * helppoint;
    mutable bool helpsign;
    mutable bool has_helppoint;

    Point4DContainer& pcont;

    /// calculate absolute value of determinant, i.e. sqrt(det(A^T A)) with A the trafo matrix
    inline void calc_abs_determinant() const;
    /// calculate absolute value of determinant, i.e. sqrt(det(A^T A)) with A the trafo matrix
    inline void calc_normal() const;
public:

    const SArrayCL<const Point4DCL*,4> & GetPoints() const {return p;}

    // reference tetra4D
    Tetra4DCL(Point4DContainer& pcont_in): absdet_of_trafo_initialized(false),normal_initialized(false),has_helppoint(false),pcont(pcont_in){
        p[0] = pcont(Point4DCL(0.0,0.0,0.0,0.0));
        p[1] = pcont(Point4DCL(1.0,0.0,0.0,0.0));
        p[2] = pcont(Point4DCL(0.0,1.0,0.0,0.0));
        p[3] = pcont(Point4DCL(0.0,0.0,1.0,0.0));
    }
    
    Tetra4DCL(Point4DContainer& pcont_in, const Point4DCL & points1, const Point4DCL & points2, const Point4DCL & points3, const Point4DCL & points4, bool points_already_in_pcont = false): absdet_of_trafo_initialized(false),normal_initialized(false),has_helppoint(false),pcont(pcont_in){
        p[0] = points_already_in_pcont ? &points1 : pcont(points1);
        p[1] = points_already_in_pcont ? &points2 : pcont(points2);
        p[2] = points_already_in_pcont ? &points3 : pcont(points3);
        p[3] = points_already_in_pcont ? &points4 : pcont(points4);
    };

    Tetra4DCL(const Tetra4DCL& tet4d):
        p( tet4d.p),
        absdet_of_trafo( tet4d.absdet_of_trafo),
        absdet_of_trafo_initialized( tet4d.absdet_of_trafo_initialized),
        normal( tet4d.normal),
        normal_initialized( tet4d.normal_initialized),
        helppoint( tet4d.helppoint),
        helpsign( tet4d.helpsign),
        has_helppoint( tet4d.has_helppoint),
        pcont( tet4d.pcont)
    {
    };

    void SetHelpPoint(const Point4DCL & p, bool sign, bool already_in_pcont = false) const
    {
        has_helppoint = true;
        helppoint = already_in_pcont ? &p : pcont(p);
        helpsign = sign;
    }
    
    Tetra4DCL& operator= (const Tetra4DCL &)
    {
        throw DROPSErrCL("Tetra4DCL& operator= called");
        return *this;
    }

    bool HasHelpPoint() { return has_helppoint; };

    double Measure() const;
    double AbsDet() const;

    // return space time normal
    Point4DCL Normal() const;

    // returns 1.0/(sqrt(1+w_n^2)) = ||n_s|| where n_s is the spatial part of the space time normal 
    double Nu() const;

    friend  std::ostream& operator<<(std::ostream& os, const Tetra4DCL &);

};

std::ostream& operator<<(std::ostream& os, const Tetra4DCL & penta);

/*
  Class constituting a pentatope (a.k.a. simplex-4, pentachoron, .. ) 
  as the convex hull of five Points(4D)
 */
class PentatopeCL
{
protected:
    Point4DContainer& pcont;

    /// vertices of the Pentatope:
    SArrayCL<const Point4DCL*,5> p;
    mutable double absdet_of_trafo;
    mutable bool absdet_of_trafo_initialized;

    /// calculate absolute value of determinant
    inline void calc_abs_determinant() const;
public:
    
    const SArrayCL<const Point4DCL*,5> & GetPoints() const {return p;}

    PentatopeCL(Point4DContainer& pcont_in, 
            const Point4DCL & points1, 
            const Point4DCL & points2, 
            const Point4DCL & points3, 
            const Point4DCL & points4, 
            const Point4DCL & points5, 
            bool already_in_pcont = false)
        : pcont(pcont_in),
        absdet_of_trafo_initialized(false)
        
    {
        p[0] = already_in_pcont ? &points1 : pcont(points1);
        p[1] = already_in_pcont ? &points2 : pcont(points2);
        p[2] = already_in_pcont ? &points3 : pcont(points3);
        p[3] = already_in_pcont ? &points4 : pcont(points4);
        p[4] = already_in_pcont ? &points5 : pcont(points5);
    };

    PentatopeCL(const PentatopeCL& penta): pcont( penta.pcont), p( penta.p), absdet_of_trafo_initialized(false)
    {
    };

    PentatopeCL(Point4DContainer& pcont_in, SArrayCL<Point4DCL,5> points): pcont(pcont_in), absdet_of_trafo_initialized(false){
        p[0] = pcont(points[0]);
        p[1] = pcont(points[1]);
        p[2] = pcont(points[2]);
        p[3] = pcont(points[3]);
        p[4] = pcont(points[4]);
    };
    // reference penta
    PentatopeCL(Point4DContainer& pcont_in): pcont(pcont_in), absdet_of_trafo_initialized(false){
        p[0] = pcont(Point4DCL(0.0,0.0,0.0,0.0));
        p[1] = pcont(Point4DCL(1.0,0.0,0.0,0.0));
        p[2] = pcont(Point4DCL(0.0,1.0,0.0,0.0));
        p[3] = pcont(Point4DCL(0.0,0.0,1.0,0.0));
        p[4] = pcont(Point4DCL(0.0,0.0,0.0,1.0));
    }

    double Measure() const;
    double AbsDet() const;

    PentatopeCL& operator= (const PentatopeCL & )
    {
        throw DROPSErrCL("PentatopeCL& operator= called");
        return *this;
    }

    // generate lsetvalues for all penta vertices
    template <class EvalObj, class T>
    void eval_timeinterpol_func_at_verts(EvalObj fx_at_t0, EvalObj fx_at_t1, 
                                                    SArrayCL<T,5> & ret)
    {
        BaryCoordCL spacebarycoord;
        double timecoord;
        T value0, value1;
        for (Uint i = 0; i < 5; ++i)
        {
            spacebarycoord = MakeBaryCoord(RestrictToPoint3D((*p[i])));
            timecoord = (*p[i])[3];
            value0 = fx_at_t0(spacebarycoord);
            value1 = fx_at_t1(spacebarycoord);
            ret[i] = (1-timecoord) * value0 + timecoord * value1;
        }
    }

    // generate lsetvalues for all penta vertices
    template <class EvalObj, class T>
    void eval_func_at_verts(EvalObj fxt, SArrayCL<T,5> & ret)
    {
        for (Uint i = 0; i < 5; ++i){
            ret[i] = fxt(MakePoint3D( (*p[i])[0], (*p[i])[1], (*p[i])[2]),(*p[i])[3]);
        }
    }

    // generate lsetvalues for all penta vertices
    // assumption: this pentatopeCL-object is part of the reference prism4 or reference pentatope
    template <class EvalObj, class T>
    void eval_func_at_verts(EvalObj fxt, const TetraCL & tet, const TimeInterval & tatb, SArrayCL<T,5> & ret)
    {
        for (Uint i = 0; i < 5; ++i){
            const double time = tatb.first * (1-(*p[i])[3]) + tatb.second * (*p[i])[3];
            Point3DCL x(0.0);
            BaryCoordCL b = MakeBaryCoord(RestrictToPoint3D((*p[i])));
            for (Uint j = 0; j < 4; ++j)
                x += b[j] * tet.GetVertex( j)->GetCoord();
            ret[i] = fxt(x,time);
        }
    }


    // decomposition into positive and negative pentatopes. The level set values already 
    // define a P1 representation. They could have been obtained via interpolation 
    // or L2-Projection or whatever... 
    void decompose_add_signed_pentas_4dtets(SArrayCL<double,5> & lsetvals, 
                                            std::vector<PentatopeCL> & negpentas, 
                                            std::vector<PentatopeCL> & pospentas, 
                                            std::vector<Tetra4DCL> & iftets);

    void decompose_signed_pentas(SArrayCL<double,5> & lsetvals, 
                                 std::vector<PentatopeCL> & negpentas, 
                                 std::vector<PentatopeCL> & pospentas, 
                                 std::vector<Tetra4DCL> & iftets)
    {
        negpentas.clear();
        pospentas.clear();
        iftets.clear();
        decompose_add_signed_pentas_4dtets(lsetvals, negpentas, pospentas, iftets);

    }

    virtual ~PentatopeCL(){};

    friend  std::ostream& operator<<(std::ostream& os, PentatopeCL);

};

std::ostream& operator<<(std::ostream& os, PentatopeCL penta);

/*
  helper class constituting a Generalized prism-4
  generalized prism-4: 4D geometry which defines the convex hull of 8 points. 
  Assumption: the eight points can be divided into two groups with 4 points
  each, where both groups define a valid (with pos. measure) tetrahedra.
*/

class GeneralizedPrism4CL
{
protected:
    SArrayCL<const Point4DCL*,4> x;
    SArrayCL<const Point4DCL*,4> y;
    Point4DContainer& pcont;
public:
    GeneralizedPrism4CL(Point4DContainer& pcont_in, SArrayCL<const Point4DCL*,4> px, SArrayCL<const Point4DCL*,4> py) : pcont(pcont_in)
    {
        for (Uint i = 0; i < 4; ++i)
        {
            x[i] = pcont(*px[i]);
            y[i] = pcont(*py[i]);
        }
    }

    GeneralizedPrism4CL(Point4DContainer& pcont_in, SArrayCL<const Point4DCL*,8> p) : pcont(pcont_in)
    {
        for (Uint i = 0; i < 4; ++i)
        {
            x[i] = pcont(*p[i]);
            y[i] = pcont(*p[i+4]);
        }
    }

    GeneralizedPrism4CL(Point4DContainer& pcont_in, SArrayCL<const Point4DCL*,4> px, Point4DCL py) : pcont(pcont_in)
    {
        for (Uint i = 0; i < 4; ++i)
        {
            x[i] = pcont(*px[i]);
            y[i] = pcont(*px[i]+py);
        }
    }


    // T x [0,1]
    GeneralizedPrism4CL(Point4DContainer& pcont_in, const TetraCL & tet, const TimeInterval & ti) : pcont(pcont_in)
    {
        for (Uint i = 0; i < 4; ++i)
        {
            Point4DCL xx(0.), yy(0.);
            for (Uint j = 0; j < 3; ++j){
                xx[j] = tet.GetVertex(i)->GetCoord()[j];
                yy[j] = tet.GetVertex(i)->GetCoord()[j];
            }
            xx[3] = ti.first;
            yy[3] = ti.second;
            
            x[i] = pcont(xx);
            y[i] = pcont(yy);
        }
    }

    // ref-T x [0,1]
    GeneralizedPrism4CL(Point4DContainer& pcont_in) : pcont(pcont_in)
    {
        Point4DCL zero(0.0);
        Point4DCL zerot(0.0);
        zerot[3] = 1.0;
        for (Uint i = 0; i < 4; ++i)
        {
            Point4DCL xx(0.), yy(0.);
            xx = zero;
            yy = zerot;
            if (i>0)
            {
                xx[i-1] = 1.0;
                yy[i-1] = 1.0;
            }
            x[i] = pcont(xx);
            y[i] = pcont(yy);
        }
    }

    // ref-T x [0,dt]
    GeneralizedPrism4CL(Point4DContainer& pcont_in, double dt) : pcont(pcont_in)
    {
        Point4DCL zero(0.0);
        Point4DCL zerot(0.0);
        zerot[3] = dt;
        for (Uint i = 0; i < 4; ++i)
        {
            Point4DCL xx(0.), yy(0.);
            xx = zero;
            yy = zerot;
            if (i>0)
            {
                xx[i-1] = 1.0;
                yy[i-1] = 1.0;
            }
            x[i] = pcont(xx);
            y[i] = pcont(yy);
        }
    }

    GeneralizedPrism4CL& operator= (const GeneralizedPrism4CL &)
    {
        throw DROPSErrCL("GeneralizedPrism4CL& operator= called");
        return *this;
    }

    void decompose_add_to_pentas (std::vector<PentatopeCL> & pentas, 
                                  int plattice_1d_els = 1, int time_1d_els = 1) const;
    
    void decompose_to_pentas (std::vector<PentatopeCL> & pentas, 
                              int plattice_1d_els = 1, int time_1d_els = 1) const
    {
        pentas.clear();
        decompose_add_to_pentas(pentas, plattice_1d_els, time_1d_els);
    }
 
    friend  std::ostream& operator<<(std::ostream& os, GeneralizedPrism4CL);

    virtual ~GeneralizedPrism4CL(){};

}; // end of class GeneralizedPrism4CL

std::ostream& operator<<(std::ostream& os, GeneralizedPrism4CL genprism4);


/**
   Helper class constituting a hyper triangle
 */
class HyperTrigCL
{
protected:
    SArrayCL<const Point4DCL*,3> u;
    SArrayCL<const Point4DCL*,3> v;
    SArrayCL<const Point4DCL*,3> w;
    Point4DContainer pcont;

public:
    HyperTrigCL(Point4DContainer& pcont_in,
                SArrayCL<const Point4DCL*,3> & iu, 
                SArrayCL<const Point4DCL*,3> & iv, 
                SArrayCL<const Point4DCL*,3> & iw)
        :u(iu),v(iv),w(iw), pcont(pcont_in)
    {
    };

    HyperTrigCL(Point4DContainer& pcont_in) : pcont(pcont_in){
        Point2DCL chi[3];
        SArrayCL<const Point4DCL*,3> * trigs[3];
        trigs[0] = &u; trigs[1] = &v; trigs[2] = &w;
        chi[0][0] = 0.0;   chi[0][1] = 0.0; 
        chi[1][0] = 1.0;   chi[1][1] = 0.0; 
        chi[2][0] = 0.0;   chi[2][1] = 1.0; 

        for (Uint i = 0; i < 3; ++i)
        {
            for (Uint j = 0; j < 3; ++j)
            {
                Point4DCL ps(0.);
                for (Uint k = 0; k < 2; ++k)
                    ps[k] = chi[j][k];
                for (Uint k = 0; k < 2; ++k)
                    ps[k+2] = chi[i][k];
                (*trigs[i])[j] = pcont(ps);
            }
        }
    }

    HyperTrigCL& operator= (const HyperTrigCL &)
    {
        throw DROPSErrCL("HyperTrigCL& operator= called");
        return *this;
    }

    void decompose_add_to_pentas (std::vector<PentatopeCL> & pentas);

    void decompose_to_pentas (std::vector<PentatopeCL> & pentas)
    {
        pentas.clear();
        decompose_add_to_pentas(pentas);
    }

    virtual ~HyperTrigCL(){};

    friend  std::ostream& operator<<(std::ostream& os, HyperTrigCL);

}; // end of class HyperTrig

std::ostream& operator<<(std::ostream& os, HyperTrigCL hypert);



} // end of namespace DROPS

#endif
