/// \file deformation.cpp
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

// #include "geom/simplex.h"
// #include "geom/isoparamP2.h"
// #include "num/fe.h"I
#include "geom/deformation.h"

namespace DROPS
{

MeshDeformationCL& MeshDeformationCL::getInstance()
{
    static MeshDeformationCL instance;
    return instance;
}

void MeshDeformationCL::MaybeUpdateNumbering()
{
    if (mg_==nullptr)
        throw DROPSErrCL("MeshDeformationCL::MaybeUpdateNumbering: use Initialize( MultiGridCL*) first!");
    if (mg_->GetVersion() == mgVersion_)
        return;

    const Uint lvl= mg_->GetLastLevel();
    mlidx_.DeleteNumbering( *mg_);
    mlidx_.CreateNumbering( lvl , *mg_, *bnd_);
    pointsol_.SetIdx( &mlidx_);
    mgVersion_= mg_->GetVersion();
}

/// Apply a transformation to the reference mesh and fill the (P2)deformation vector 
/// accordingly. The flag "only_bnd_edges_curved" allows to basically reduce the number
/// of curved edges (and thus tetrahedra) by removing curvature information at inner edges
/// and replacing them with the mean of the vertex values.
void MeshDeformationCL::SetMeshTransformation(instat_vector_fun_ptr f, const double t, 
                                              bool only_bnd_edges_curved, bool P2){
    MaybeUpdateNumbering();
    const Uint pidx = mlidx_.GetIdx();
    DROPS_FOR_TRIANG_VERTEX( (*mg_), mg_->GetLastLevel(), it) {
        if (!it->Unknowns.Exist( pidx)) continue;
        const Point3DCL val(f(it->GetCoord(),t));
        for (Uint k = 0; k < 3; ++k)
            pointsol_.Data[it->Unknowns(pidx)+k] = val[k];
    }

    DROPS_FOR_TRIANG_EDGE( (*mg_), mg_->GetLastLevel(), it) {
        if (!it->Unknowns.Exist( pidx)) continue;
        if (P2 && (!only_bnd_edges_curved || it->IsOnBoundary())) // potentially curved elements
        {
            const Point3DCL val(f(GetBaryCenter(*it),t));
            for (Uint k = 0; k < 3; ++k)
                pointsol_.Data[it->Unknowns(pidx)+k] = val[k];
        }
        else // planar values, i.e. average of vertex values
        {
            const VertexCL& vt1 (*it->GetVertex(0));
            const VertexCL& vt2 (*it->GetVertex(1));
            if (!vt1.Unknowns.Exist( pidx)) continue;
            if (!vt2.Unknowns.Exist( pidx)) continue;
            if (!it->Unknowns.Exist( pidx)) continue;
            Point3DCL a,b;
            for (Uint k = 0; k < 3; ++k)
            {
                pointsol_.Data[it->Unknowns(pidx)+k] =
                    0.5 * pointsol_.Data[vt1.Unknowns(pidx)+k]
                  + 0.5 * pointsol_.Data[vt2.Unknowns(pidx)+k];
            }
        }
    }
    CheckForCurved();
}


// Fill the Deformation vector with the coordinates of the reference mesh
void MeshDeformationCL::SetMeshIdentity(){
    SetMeshTransformation(&instat_Identity3D,0.0);
}


// For inner edges curved elements are typically not necessary. For ease for implementation
// it might still be that an applied transformation lead to curved edges inside the domain.
// This function replaces the curved values on the edges with the average of the vertex values
// leading to a non-curved situation for all inner edges and thus for all inner elements
void MeshDeformationCL::SetInnerEdgesPlanar(){
    MaybeUpdateNumbering();
    DROPS_FOR_TRIANG_EDGE( (*mg_), mg_->GetLastLevel(), it) {
        if (it->IsOnBoundary()) continue; // the remainder is only done for inner edges!

        const Uint pidx = mlidx_.GetIdx();
        const VertexCL& vt1 (*it->GetVertex(0));
        const VertexCL& vt2 (*it->GetVertex(1));
        if (!vt1.Unknowns.Exist( pidx)) continue;
        if (!vt2.Unknowns.Exist( pidx)) continue;
        if (!it->Unknowns.Exist( pidx)) continue;
        
        Point3DCL a,b;
        for (Uint k = 0; k < 3; ++k)
        {
            pointsol_.Data[it->Unknowns(pidx)+k] =
                0.5 * pointsol_.Data[vt1.Unknowns(pidx)+k]
              + 0.5 * pointsol_.Data[vt2.Unknowns(pidx)+k];
        }
    }
}

/// Runs over all tetrahedra, checks if any aligned edge is curved and stores this in the mapping
/// tet_is_curved. Note that an edge is considered as curved if the angle is larger than a given
/// threshold, here 1e-8.
void MeshDeformationCL::CheckForCurved(){
    if (mg_==nullptr)
        throw DROPSErrCL("MeshDeformationCL::CheckForCurved: No MultiGridCL* given!");
    Ulint curvedels = 0;
    Ulint els = 0;
    DROPS_FOR_TRIANG_TETRA( (*mg_), mg_->GetLastLevel(), it) {
        const Uint pidx = mlidx_.GetIdx();
        bool curved = false;
        for (Uint i = 0; i < 6; ++i)
        {
            const EdgeCL& edge (*it->GetEdge(i));
            const VertexCL& vt1 (*edge.GetVertex(0));
            const VertexCL& vt2 (*edge.GetVertex(1));
            if (!vt1.Unknowns.Exist( pidx)) continue;
            if (!vt2.Unknowns.Exist( pidx)) continue;
            Point3DCL a,b;
            for (Uint k = 0; k < 3; ++k)
            {
                a[k] = pointsol_.Data[vt1.Unknowns(pidx)+k];
                b[k] = pointsol_.Data[vt2.Unknowns(pidx)+k];
            }
            const Point3DCL c = 0.5 * a + 0.5 * b; // edge midpoint (uncurved)
            const Point3DCL d = b - a; // difference of vert coords
            const double edgelength = d.norm();
            if (!edge.Unknowns.Exist( pidx)) continue;
            Point3DCL dc;
            for (Uint k = 0; k < 3; ++k)
                dc[k] = pointsol_.Data[edge.Unknowns(pidx)+k] - c[k];
            if (dc.norm() > 2e-8 * edgelength) // angle is larger than approx. 1e-8
            {
                curved = true;
                break;
            }
        }
        tet_is_curved[&(*it)] = curved;
        els ++;
        if (curved)
            curvedels ++;
    }
    std::cout << curvedels << " out of " << els << " tetrahedra are curved." << std::endl;
}

// Manipulate the deformation of a single edge midpoint
void MeshDeformationCL::SetEdgeDeformation(const EdgeCL& edge, const Point3DCL & p){
    MaybeUpdateNumbering();
    const Uint pidx = mlidx_.GetIdx();
    if (!edge.Unknowns.Exist( pidx)) return;
    for (Uint k = 0; k < 3; ++k)
        pointsol_.Data[edge.Unknowns(pidx)+k] = p[k];
}

bool MeshDeformationCL::IsTetraCurved(const TetraCL& tet)
{
     return tet_is_curved[&tet];
}

void MeshDeformationCL::Initialize( MultiGridCL* mg){
    if (mg==nullptr)
        throw DROPSErrCL("MeshDeformationCL::Initialize: No MultiGrid given!");
    std::cout << " Initializing MeshDeformationCL " << std::endl;

    mg_ = mg;
    const Usint nb = mg->GetBnd().GetNumBndSeg();
    if (bnd_) delete bnd_;
    bnd_ = new BndDataCL<Point3DCL>(nb);

    tet_is_curved.clear();
    DROPS_FOR_TRIANG_TETRA( (*mg_), mg_->GetLastLevel(), it) {
        tet_is_curved[&(*it)] = false;
    }
    mgVersion_= 0; // force update of numbering by SetMeshIdentity()
    SetMeshIdentity();
}

LocalP2CL<Point3DCL> MeshDeformationCL::GetLocalP2Deformation( const TetraCL& tet) const {
    return LocalP2CL<Point3DCL>(tet, pointsol_, *bnd_);
}

LocalP1CL<Point3DCL> MeshDeformationCL::GetLocalP1Deformation( const TetraCL& tet) const {
    return LocalP1CL<Point3DCL>(tet, pointsol_, *bnd_);
}

Point3DCL MeshDeformationCL::GetTransformedVertexCoord( const VertexCL &v) const
{
    Point3DCL ret;
    const Uint pidx( mlidx_.GetIdx());
    for (Uint k = 0; k < 3; ++k)
        ret[k] = pointsol_.Data[v.Unknowns(pidx)+k];
    return ret;
}

Point3DCL MeshDeformationCL::GetTransformedEdgeBaryCenter( const EdgeCL &v) const
{
    Point3DCL ret;
    const Uint pidx( mlidx_.GetIdx());
    for (Uint k = 0; k < 3; ++k)
        ret[k] =  pointsol_.Data[v.Unknowns(pidx)+k];
    return ret;
}

Point3DCL MeshDeformationCL::GetTransformedTetraBaryCenter( const TetraCL & v) const
{
    Point3DCL ret;
    const VertexCL& vt1 (*v.GetVertex(0));
    const VertexCL& vt2 (*v.GetVertex(1));
    const VertexCL& vt3 (*v.GetVertex(2));
    const VertexCL& vt4 (*v.GetVertex(3));
    const Uint pidx( mlidx_.GetIdx());
    for (Uint k = 0; k < 3; ++k)
        ret[k] =   0.25 * pointsol_.Data[vt1.Unknowns(pidx)+k]
                 + 0.25 * pointsol_.Data[vt2.Unknowns(pidx)+k]
                 + 0.25 * pointsol_.Data[vt3.Unknowns(pidx)+k]
                 + 0.25 * pointsol_.Data[vt4.Unknowns(pidx)+k];
    return ret;
}

} // end of namespace DROPS
