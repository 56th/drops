/// \file geomid.h
/// \brief The DiST (short for distributed simplex type) module is responsible for the distributed geometric data structure on a parallel machine.
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Yuanjun Zhang; SC RWTH Aachen: Oliver Fortmeier, Daniel Medina Cardona.

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

#ifndef GEOMID_H_
#define GEOMID_H_

#include "misc/container.h"

namespace DROPS
{

class VertexCL;
class EdgeCL;
class FaceCL;
class TetraCL;

namespace DiST{

// \name Ask for the dimension of a simplex
//@{
template <typename SimplexT> inline Usint GetDim();
template <> inline Usint GetDim<VertexCL>() { return 0; }
template <> inline Usint GetDim<EdgeCL>()   { return 1; }
template <> inline Usint GetDim<FaceCL>()   { return 2; }
template <> inline Usint GetDim<TetraCL>()  { return 3; }
template <> inline Usint GetDim<const VertexCL>() { return 0; }
template <> inline Usint GetDim<const EdgeCL>()   { return 1; }
template <> inline Usint GetDim<const FaceCL>()   { return 2; }
template <> inline Usint GetDim<const TetraCL>()  { return 3; }
//@}

/// \brief Assign each simplex an unique geometric id
struct GeomIdCL
{
    Uint        level;      ///< level, the simplex occurs first
    Point3DCL   bary;       ///< barycenter of the simplex
    Usint       dim;        ///< Dimension of the simplex, i.e., 0 - vertex, 1 - edge, 2 - face, 3 - tetrahedron, 4 - uninitialized

    GeomIdCL() : level((Uint)(-1)), bary(), dim(4) {}
    GeomIdCL(Uint lvl, const Point3DCL& p, Usint dimension) : level(lvl), bary(p), dim(dimension) {}
    GeomIdCL( const GeomIdCL& h) : level(h.level), bary(h.bary), dim(h.dim) {}
    template <typename SimplexT>
    GeomIdCL(Uint lvl, const SimplexT& s) : level(lvl), bary( ComputeBaryCenter(s)), dim(GetDim<SimplexT>()) {}
    bool operator== (const GeomIdCL& h) const { return h.level == level && h.bary == bary;}
    bool operator!= (const GeomIdCL& h) const { return !(h==*this); }
    bool operator < (const GeomIdCL& h) const { return level < h.level && dim < h.dim && bary[0] < h.bary[0] && bary[1] < h.bary[1] && bary[2] < h.bary[2];}
};

inline std::ostream& operator << ( std::ostream& os, const GeomIdCL& h)
{
    static char scode[]= "VEFT?"; // simplex code for each dimension
    os << scode[h.dim] << h.level << " (" << h.bary << ')';
    return os;
}
const GeomIdCL NoGID= GeomIdCL( (Uint)(-1), Point3DCL(), 4);  // Dummy id, if simplex does not exist

#if __GNUC__ >= 4 || DROPS_WIN || __xlC__
struct Hashing : std::unary_function<GeomIdCL, size_t>
{
    size_t operator()(const GeomIdCL& h) const
    {
        size_t seed = 0;
        DROPS_STD_HASH<int> inthasher;
        DROPS_STD_HASH<double> doublehasher;
        // see boost::hash_combine
        seed ^= inthasher(h.level) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= doublehasher(h.bary[0]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= doublehasher(h.bary[1]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= doublehasher(h.bary[2]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};
#endif
} // end of namespace DiST
} // end of namespace DROPS


#endif /* GEOMID_H_ */
