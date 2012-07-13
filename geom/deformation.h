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
#include "misc/problem.h"
#include "num/bndData.h"
#include "num/discretize.h"

namespace DROPS
{

class MeshDeformationCL
{
private:
    MultiGridCL * mg_;
    MLIdxDescCL * mlidx_; // perhaps MultiLevel Idx just for future use...
    VecDescCL * pointsol_;
    BndDataCL<Point3DCL> * bnd_;
    MeshDeformationCL() : mg_(0),mlidx_(0),pointsol_(0),bnd_(0){};
    virtual ~MeshDeformationCL(){
        if (mg_) delete mg_;
        if (mlidx_) delete mlidx_;
        if (pointsol_) delete pointsol_;
        if (bnd_) delete bnd_;
    };
public:
    static MeshDeformationCL& getInstance();

    void Initialize( MultiGridCL* mg);

    LocalP2CL<Point3DCL> GetLocalP2Deformation( const TetraCL&);
    // LocalP1CL<Point3DCL> GetLocalP1Deformation( const TetraCL&);
};

} // end of namespace DROPS

#endif
