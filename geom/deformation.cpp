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
// #include "num/fe.h"
#include "geom/deformation.h"

namespace DROPS
{

MeshDeformationCL& MeshDeformationCL::getInstance()
{
    static MeshDeformationCL instance;
    return instance;
}

void MeshDeformationCL::SetMeshTransformation(instat_vector_fun_ptr f, const double t){
    if (mg_==NULL)
        throw DROPSErrCL("MeshDeformationCL::SetMeshIdentity: No MultiGridCL* given!");

    Uint pidx( mlidx_->GetIdx());
    DROPS_FOR_TRIANG_VERTEX( (*mg_), mg_->GetLastLevel(), it) {
        if (!it->Unknowns.Exist( pidx)) continue;
        const Point3DCL val(f(it->GetCoord(),t));
        for (Uint k = 0; k < 3; ++k)
            (*pointsol_).Data[it->Unknowns(pidx)+k] = val[k];
    }

    DROPS_FOR_TRIANG_EDGE( (*mg_), mg_->GetLastLevel(), it) {
        if (!it->Unknowns.Exist( pidx)) continue;
        const Point3DCL val(f(GetBaryCenter(*it),t));
        for (Uint k = 0; k < 3; ++k)
            (*pointsol_).Data[it->Unknowns(pidx)+k] = val[k];
    }
    std::cout << " (*pointsol_).Data = \n " << (*pointsol_).Data << std::endl;
}

void MeshDeformationCL::SetMeshIdentity(){
    SetMeshTransformation(&instat_Identity3D,0.0);
}

void MeshDeformationCL::CheckForCurved(){
    if (mg_==NULL)
        throw DROPSErrCL("MeshDeformationCL::CheckForCurved: No MultiGridCL* given!");

    std::cout << " MeshDeformationCL::CheckForCurved : not implemented yet" << std::endl;
    DROPS_FOR_TRIANG_TETRA( (*mg_), mg_->GetLastLevel(), it) {
        // check if Tetra has a curved Edge...
        tet_is_curved[&(*it)] = false;
    }
}


bool MeshDeformationCL::IsTetraCurved(const TetraCL& tet){
    return tet_is_curved[&tet];
}

void MeshDeformationCL::Initialize( MultiGridCL* mg){
    std::cout << " Initializing MeshDeformationCL " << std::endl;

    mg_ = mg;
    mlidx_ = new MLIdxDescCL( vecP2_FE);
    const Uint lvl = mg_->GetLastLevel();
    const Usint nb = mg->GetBnd().GetNumBndSeg();
    bnd_ = new BndDataCL<Point3DCL>(nb); 
    mlidx_->CreateNumbering( lvl , *mg_, *bnd_);
    pointsol_ = new VecDescCL( mlidx_);

    tet_is_curved.clear();
    DROPS_FOR_TRIANG_TETRA( (*mg_), mg_->GetLastLevel(), it) {
        tet_is_curved[&(*it)] = false;
    }

    SetMeshIdentity();
}

LocalP2CL<Point3DCL> MeshDeformationCL::GetLocalP2Deformation( const TetraCL& tet){
    return LocalP2CL<Point3DCL>(tet, *pointsol_, *bnd_);
}

LocalP1CL<Point3DCL> MeshDeformationCL::GetLocalP1Deformation( const TetraCL& tet){
    LocalP1CL<Point3DCL> ret;
    const Uint pidx( mlidx_->GetIdx());
    for (Uint i = 0; i < 4; ++i)
    {
        const VertexCL& v(*tet.GetVertex(i));
        if (!v.Unknowns.Exist(pidx)) throw DROPSErrCL(" MeshDeformationCL::GetLocalP1Deformation: no vertex value!");
        for (Uint k = 0; k < 3; ++k)
            ret[i][k] = (*pointsol_).Data[v.Unknowns(pidx)+k];
    }
    return ret;
}

Point3DCL MeshDeformationCL::GetTransformedVertexCoord( const VertexCL &v)
{
    Point3DCL ret;
    const Uint pidx( mlidx_->GetIdx());
    for (Uint k = 0; k < 3; ++k)
        ret[k] = (*pointsol_).Data[v.Unknowns(pidx)+k];
    return ret;
}

Point3DCL MeshDeformationCL::GetTransformedEdgeBaryCenter( const EdgeCL &v)
{
    Point3DCL ret;
    const Uint pidx( mlidx_->GetIdx());
    for (Uint k = 0; k < 3; ++k)
        ret[k] = (*pointsol_).Data[v.Unknowns(pidx)+k];
    return ret;
}

} // end of namespace DROPS
