/// \file deformation.h
/// \brief classes for describing global mesh deformations (for isoparam. els.)
/// \author LNM RWTH Aachen: Christoph Lehrenfeld, Liang Zhang

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

#ifndef DROPS_DEFORMATION_H
#define DROPS_DEFORMATION_H

//#include "geom/simplex.h"
//#include "misc/container.h"
//#include "misc/problem.h"
//#include "num/bndData.h"
//#include "geom/multigrid.h"
#include "num/discretize.h"

namespace DROPS
{

class MeshDeformationCL
{
private:

    MultiGridCL * mg_;
    MLIdxDescCL mlidx_; // perhaps MultiLevel Idx just for future use...
    VecDescCL pointsol_;
    BndDataCL<Point3DCL> * bnd_;
    size_t mgVersion_;
    std::map<const TetraCL*, bool> tet_is_curved;
    
    MeshDeformationCL() : mg_(0), mlidx_( vecP2_FE), pointsol_(&mlidx_), bnd_(0), mgVersion_(0)
    {}
    virtual ~MeshDeformationCL()
    { if (bnd_) delete bnd_; };
    void MaybeUpdateNumbering(); ///< depending on the version of the multigrid, the numbering of mlidx_ is updated

    static Point3DCL instat_Identity3D (const Point3DCL & a, const double)
    {  return a;  }

public:
    static MeshDeformationCL& getInstance();

    void SetMeshIdentity();
    void SetMeshTransformation(instat_vector_fun_ptr f, const double t, bool only_bnd_edges_curved = false, bool P2 = false);

    void CheckForCurved();
    void Initialize( MultiGridCL* mg);
    bool IsUsed() const { return mg_!=nullptr;}
    bool IsTetraCurved(const TetraCL& tet);
    LocalP2CL<Point3DCL> GetLocalP2Deformation( const TetraCL&) const;
    LocalP1CL<Point3DCL> GetLocalP1Deformation( const TetraCL&) const;
    void SetEdgeDeformation(const EdgeCL& edge, const Point3DCL & p);
    void SetInnerEdgesPlanar();
    Point3DCL GetTransformedVertexCoord( const VertexCL &) const;
    Point3DCL GetTransformedEdgeBaryCenter( const EdgeCL &) const;
    Point3DCL GetTransformedTetraBaryCenter( const TetraCL &) const;

};



/** calculates the transpose of the transformation  Tetra -> RefTetra
 *  if \f$ \Phi \f$ denotes the trafo from \f$ T_ref \f$ to \f$ T \f$
 *  then the result is \f$ T = (\nabla \Phi)^{-T} \f$ and 
 *  det \f$ = det(\nabla \Phi) \f$
 *  Inputs are p the point on the reference triangle and pt, the four points 
 *  of the transformed tetraeder.  */
inline void GetTrafoTr( SMatrixCL<3,3>& T, double& det, const LocalP1CL<Point3DCL>& pt)
{
    double M[3][3];
    const Point3DCL& pt0= pt[0];
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            M[j][i]= pt[i+1][j] - pt0[j];
    det=   M[0][0] * (M[1][1]*M[2][2] - M[1][2]*M[2][1])
         - M[0][1] * (M[1][0]*M[2][2] - M[1][2]*M[2][0])
         + M[0][2] * (M[1][0]*M[2][1] - M[1][1]*M[2][0]);

    T(0,0)= (M[1][1]*M[2][2] - M[1][2]*M[2][1])/det;
    T(0,1)= (M[2][0]*M[1][2] - M[1][0]*M[2][2])/det;
    T(0,2)= (M[1][0]*M[2][1] - M[2][0]*M[1][1])/det;
    T(1,0)= (M[2][1]*M[0][2] - M[0][1]*M[2][2])/det;
    T(1,1)= (M[0][0]*M[2][2] - M[2][0]*M[0][2])/det;
    T(1,2)= (M[2][0]*M[0][1] - M[0][0]*M[2][1])/det;
    T(2,0)= (M[0][1]*M[1][2] - M[1][1]*M[0][2])/det;
    T(2,1)= (M[1][0]*M[0][2] - M[0][0]*M[1][2])/det;
    T(2,2)= (M[0][0]*M[1][1] - M[1][0]*M[0][1])/det;
}

inline void GetTrafoTrDet(double& det, const LocalP1CL<Point3DCL>& pt)
{
    double M[3][3];
    const Point3DCL& pt0= pt[0];
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            M[j][i]= pt[i+1][j] - pt0[j];
    det=   M[0][0] * (M[1][1]*M[2][2] - M[1][2]*M[2][1])
         - M[0][1] * (M[1][0]*M[2][2] - M[1][2]*M[2][0])
         + M[0][2] * (M[1][0]*M[2][1] - M[1][1]*M[2][0]);
}


} // end of namespace DROPS

#endif
