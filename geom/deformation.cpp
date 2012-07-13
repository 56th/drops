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

void MeshDeformationCL::Initialize( MultiGridCL* mg){
    mg_ = mg;
    mlidx_ = new MLIdxDescCL( P2_FE);
    
    const Usint nb = mg->GetBnd().GetNumBndSeg();
    bnd_ = new BndDataCL<Point3DCL>(nb); 
    mlidx_->CreateNumbering( mg_->GetLastLevel() , *mg_, *bnd_);
    pointsol_ = new VecDescCL( mlidx_);
}

LocalP2CL<Point3DCL> MeshDeformationCL::GetLocalP2Deformation( const TetraCL& tet){
    return LocalP2CL<Point3DCL>(tet, *pointsol_, *bnd_);
}

// LocalP1CL<Point3DCL> MeshDeformationCL::GetLocalP1Deformation( const TetraCL& tet){
//     return LocalP1CL<Point3DCL>(LocalP2CL<Point3DCL>(tet, *pointsol_, *bnd_));
// }


} // end of namespace DROPS
