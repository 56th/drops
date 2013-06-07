/// \file boundary.h
/// \brief Boundary description
/// \author LNM RWTH Aachen:Joerg Peters, Volker Reichelt; SC RWTH Aachen:

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

#ifndef DROPS_BOUNDARY_H
#define DROPS_BOUNDARY_H


#include "misc/container.h"


namespace DROPS
{


class EdgeCL;
class BndPointCL;

typedef std::pair<Point2DCL, Point3DCL> BndPairCL;
typedef Usint BndIdxT;
const BndIdxT NoBndC= static_cast<BndIdxT>(-1);

//**************************************************************************
// Class:    BndSegCL                                                      *
// Purpose:  base class for the boundary descriptions used by DROPS        *
// Remarks:  To define your own boundary segment descriptions inherit from *
//           this class.                                                   *
//**************************************************************************

class BndSegCL
{
  private:
    bool           _IsNonPlanar;

  public:
    BndSegCL (bool IsNonPl)
        :_IsNonPlanar(IsNonPl) {}
    virtual ~BndSegCL () {}
    // Default copy-ctor, assignment-op.

    bool                  IsNonPlanar () const  { return _IsNonPlanar; }

    // true, iff p2D is in the domain of the boundary segment
    virtual bool      IsInBounds (const Point2DCL&) const = 0;
    // If IsInBounds(p2D)==true, the result is the image of p2D, else undefined
    virtual Point3DCL Map        (const Point2DCL&) const = 0;
    // Project p3D into to parameter space. It is not checked, whether p2D is in
    // the actual domain of the boundary segment. The last argument can be a hint,
    // where p2D lies - it is implementation dependant, whether it is used.
    virtual Point2DCL ProjectRaw (const Point3DCL&, const Point2DCL* = NULL) const = 0;
    // Map p3D into the domain
    virtual Point2DCL Project    (const Point3DCL&, const Point2DCL* = NULL) const = 0;
    // Calculate the midvertex coordinates of an edge, the endpoints of which have boundary-descriptions
    // bp1 & bp2. The midvertex shall lie on the boundary.
    virtual BndPairCL MidProject (const BndPointCL&, const BndPointCL&) const = 0;
};


//**************************************************************************
// Class:    BndPointCL                                                    *
// Purpose:  stores the boundary-description of a vertex                   *
// Remarks:  _Coord2D are the parametric coordinates of the vertex in the  *
//           boundary segment *_BndIdx.                                    *
//**************************************************************************

class BndPointCL
{
  private:
    BndIdxT   _BndIdx;
    Point2DCL _Coord2D;

  public:
    BndPointCL (BndIdxT BndIdx, const Point2DCL& p2D)
        : _BndIdx(BndIdx), _Coord2D(p2D) {}
    // Default copy-ctor, dtor, assignnment-op.

    // Returns the BndSegCL-object, to which these 2D-coordinates belong
    BndIdxT GetBndIdx() const { return _BndIdx; }
    // Returns the 2D-coordinates of the vertex in BndSegCL *GetBndIdx()
    const Point2DCL& GetCoord2D () const { return _Coord2D; }

#ifdef _PAR
    /// \brief Boundary-Points have to be transferred, so ParMultiGridCL needs access
    friend class ParMultiGridCL;
#endif
};


//**************************************************************************
// Class:    BndPointSegLessCL                                             *
// Purpose:  orders BndPointCL-objects by their BndSegCL-object-Id         *
// Remarks:  The boundary descriptions of a vertex are stored in this order*
//           to allow for the use of set_intersection in BuildMidVertex.   *
//**************************************************************************

class BndPointSegLessCL : public std::binary_function<BndPointCL,BndPointCL,bool>
{
  public:
    bool operator () (const BndPointCL& bp0, const BndPointCL& bp1) const
        { return bp0.GetBndIdx() < bp1.GetBndIdx(); }
};


//**************************************************************************
// Class:    BndPointSegEqCL                                               *
// Purpose:  Compare BndPointCL-objects by their BndSegCL-object-Id        *
// Remarks:                                                                *
//**************************************************************************

class BndPointSegEqCL : public std::unary_function<BndPointCL,bool>
{
  private:
    BndIdxT _bidx;

  public:
    BndPointSegEqCL(BndIdxT bidx) :_bidx(bidx) {}

    bool operator () (const BndPointCL& bp) const
        { return bp.GetBndIdx() == _bidx; }
};


//**************************************************************************
// Class:    AffineSquareCL                                                *
// Base:     BndSegCL                                                      *
// Purpose:  affine image of the 2D-unit-square (0,0), (1,0), (1,1), (0,1) *
// Remarks:                                                                *
//**************************************************************************

class AffineSquareCL : public BndSegCL
{
  private:
    Point3DCL _Orig, _d0, _d1;
    double    _d0d1, _d0sq, _d1sq, _Det;

  public:
    // Images of          (0,0),            (1,0),            (0,1).
    AffineSquareCL (const Point3DCL&, const Point3DCL&, const Point3DCL&);
    // Default copy-ctor, assignment-op.

    virtual bool      IsInBounds (const Point2DCL&) const;
    virtual Point3DCL Map        (const Point2DCL&) const;
    virtual Point2DCL ProjectRaw (const Point3DCL&, const Point2DCL* = NULL) const;
    virtual Point2DCL Project    (const Point3DCL&, const Point2DCL* = NULL) const;
    virtual BndPairCL MidProject (const BndPointCL&, const BndPointCL&) const;
};


//**************************************************************************
// Class:    AffineTriangleCL                                              *
// Base:     BndSegCL                                                      *
// Purpose:  affine image of the reference-triangle (0,0), (1,0), (0,1)    *
//**************************************************************************

class AffineTriangleCL : public BndSegCL
{
  private:
    Point3DCL _Orig, _d0, _d1;
    double    _d0d1, _d0sq, _d1sq, _Det;

  public:
    // Images of            (0,0),            (1,0),            (0,1).
    AffineTriangleCL (const Point3DCL&, const Point3DCL&, const Point3DCL&);
    // Default copy-ctor, assignment-op.

    virtual bool      IsInBounds (const Point2DCL&) const;
    virtual Point3DCL Map        (const Point2DCL&) const;
    virtual Point2DCL ProjectRaw (const Point3DCL&, const Point2DCL* = NULL) const;
    virtual Point2DCL Project    (const Point3DCL&, const Point2DCL* = NULL) const;
    virtual BndPairCL MidProject (const BndPointCL&, const BndPointCL&) const;
};


//**************************************************************************
// Class:    MeshBoundaryCL                                                *
// Base:     BndSegCL                                                      *
// Purpose:  Represents a topological boundary-segment read from a         *
//           mesh-file. Most operations are trivial.                       *
//**************************************************************************
// TODO: This will explode on refinement. We need 3D coordinates or better a vertex (2 vertices).
class MeshBoundaryCL : public BndSegCL
{
  private:
    Uint bc_; // Boundary condition as per mesh file. TODO: Use an enum.

  public:
    MeshBoundaryCL(Uint bc)
        :BndSegCL( false), bc_( bc) {}
    // Default copy-ctor, assignment-op.

    virtual bool      IsInBounds (const Point2DCL&) const {
        return true; }
    virtual Point3DCL Map        (const Point2DCL&) const {
        return Point3DCL(); }
    virtual Point2DCL ProjectRaw (const Point3DCL&, const Point2DCL* = 0) const {
        return Point2DCL(); }
    virtual Point2DCL Project    (const Point3DCL&, const Point2DCL* = 0) const {
        return Point2DCL(); }
    virtual BndPairCL MidProject (const BndPointCL&, const BndPointCL&) const {
        return std::make_pair( Point2DCL(), Point3DCL()); }
    Uint GetBC() const { return bc_; }
};


} // end of namespace DROPS

#endif
