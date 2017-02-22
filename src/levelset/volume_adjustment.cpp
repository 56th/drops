/// \file volume_adjustment.cpp
/// \brief Maintain and adjust the volume of the subdomains with positive/negative level set value.
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
 * Copyright 2017 LNM RWTH Aachen, Germany
*/

#include "levelset/volume_adjustment.h"
#include "levelset/levelset.h"
#include "misc/params.h"

namespace DROPS
{

std::unique_ptr<VolumeAdjustmentCL> VolumeAdjustmentCL::Create (LevelsetP2CL* lset, const ParamCL&)
{
    return std::unique_ptr<VolumeAdjustmentCL> (new GlobalVolumeAdjustmentCL (lset));
}

void VolumeAdjustmentCL::DebugOutput (std::ostream& os) const
{
    os << "VolumeAdjustmentCL: No volume adjustment.\n" << std::endl;
}

void GlobalVolumeAdjustmentCL::AdjustVolume ()
{
    if (global_reference_volume_ == -1.0)
        throw DROPSErrCL ("GlobalVolumeAdjustmentCL::AdjustVolume: The global_reference_volume_ has not been set, but it is required for this method (maybe you forgot SetGlobalReferenceVolume somewhere).\n");

    const double tol= tol_*global_reference_volume_;

    double v0=lset_->GetVolume(0., num_subdivision_) - global_reference_volume_;
    if (std::abs(v0)<=tol)
        return;

    // hint: surf(Kugel) = [3/4/pi*vol(Kugel)]^(2/3) * 4pi
    double d0=0,
           d1=v0*0.23/std::pow (global_reference_volume_, 2./3.);
    double v1=lset_->GetVolume(d1, num_subdivision_) - global_reference_volume_;
    if (std::abs(v1)<=tol) {
        dphi_= d1;
        goto converged;
    }

    // secant method for initial value
    while (v1*v0 > 0) // same sign
    {
        const double d2=d1-1.2*v1*(d1-d0)/(v1-v0);
        d0=d1; d1=d2; v0=v1; v1=lset_->GetVolume(d1, num_subdivision_) - global_reference_volume_;
        if (std::abs(v1)<=tol) {
            dphi_= d1;
            goto converged;
        }
    }

    // Anderson-Bjoerk for exact value
    while (true)
    {
        const double d2=(v1*d0-v0*d1)/(v1-v0),
                     v2=lset_->GetVolume(d2, num_subdivision_) - global_reference_volume_;
        if (std::abs(v2)<=tol) {
            dphi_= d2;
            goto converged;
        }

        if (v2*v1 < 0) // different signs
          { d0=d1; d1=d2; v0=v1; v1=v2; }
        else
          { const double c=1.0-v2/v1; d1=d2; v1=v2; v0*= c>0 ? c : 0.5; }
    }

    converged:
    lset_->Phi.Data+= dphi_;
}

void GlobalVolumeAdjustmentCL::DebugOutput (std::ostream& os) const
{
    os << "GlobalVolumeAdjustmentCL: shift for level set function dphi_: " << dphi_ << ", "
       "new rel. volume: " << lset_->GetVolume()/global_reference_volume_ << "." << std::endl;
}

} // end of namespace DROPS
