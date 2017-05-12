/// \file volume_adjustment.h
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

#ifndef DROPS_VOLUME_ADJUSTMENT_H
#define DROPS_VOLUME_ADJUSTMENT_H

#include <iosfwd>
#include <memory>
#include <valarray>
#include <vector>
#include "misc/container.h"
#include "num/spmat.h"
#include "out/vtkOut.h"

namespace DROPS
{

// forward declarations
class LevelsetP2CL;
class ParamCL;

class VolumeAdjustmentCL
{
  protected:
    LevelsetP2CL* lset_;

    double tol_= 1e-9;
    double global_reference_volume_= -1.0;
    int    num_subdivision_= 2;

    virtual void InitVolume_impl () {}

  public:
    static std::unique_ptr<VolumeAdjustmentCL> Create (LevelsetP2CL* lset, const ParamCL& P);

    VolumeAdjustmentCL (LevelsetP2CL* lset) : lset_ (lset) {}
    virtual ~VolumeAdjustmentCL() {}

    VolumeAdjustmentCL& SetTol (double tol) {
        tol_= tol;
        return *this;
    }
    VolumeAdjustmentCL& InitVolume (double refvol=-1.0) {
        global_reference_volume_= refvol;
        InitVolume_impl();
        return *this;
    }
    VolumeAdjustmentCL& SetNumSubdivision (int numsubdiv) {
        num_subdivision_= numsubdiv;
        return *this;
    }
    virtual void Repair () {}  // called in LevelsetRepairCL::post_refine_sequence
    virtual void AdjustVolume () {}
    virtual void DebugOutput (std::ostream&) const;
    virtual VTKVariableCL& make_VTKComponentMap(std::string var_name);
};

class GlobalVolumeAdjustmentCL : public VolumeAdjustmentCL
{
  private:
    double dphi_= 0.;

  public:
    GlobalVolumeAdjustmentCL (LevelsetP2CL* lset) : VolumeAdjustmentCL (lset) {}
    void AdjustVolume () override;
    void DebugOutput (std::ostream& os) const override;
};


/// \brief Individual volume correction for each connected component.
/// Automatically handles changes in topology.
///
/// Does not correct volumes of connected components which are too close together (distance has to be at least 4 P2 DOFs).
class ComponentBasedVolumeAdjustmentCL : public VolumeAdjustmentCL 
{
  public:
    using component_vector= std::vector<size_t>;

  private:
    ///@{ Data for each dof
    component_vector       component_of_dof_, ///< component_of_dof_[i] is the number of the component of dof i.
                           component_of_dof_backup_;
    std::vector<Point3DCL> coord_of_dof_,        ///< The points where the dofs live.
                           coord_of_dof_backup_; ///< The points where dofs_backup live (only different from coord_of_dof_ in Repair).
    ///@}

    ///@{ Data for each component
    std::vector<bool>                  doCorrection_;
    std::vector<int>                   sign_of_component_;
    std::vector<int>                   sign_of_component_backup_;
    std::vector<std::valarray<double>> indicator_functions_;
    std::vector<double>                Volumes; ///< volume per component
    std::vector<double>                targetVolumes; ///< target volume of each component
    std::vector<Point3DCL>             ReferencePoints; ///< markers for each connected component
    std::vector<Point3DCL>             ReferencePoints_backup; ///< old markers for each connected component
    ///@}

    /// \brief Initialize coord_of_dof_.
    void init_coord_of_dof ();

    /// \brief Helper of compute_indicator_functions. Extend all components (except 0) by one level (except where they would overlap).
    component_vector ExtendOneStep (const SparseMatBaseCL<unsigned char>& A, const component_vector& cp, std::vector<bool>& doCorrection) const;
    /// \brief Compute the indicator_functions_ of the extension of each connected component. Also sets doCorrection_ to false if the extensions would overlap.
    void compute_indicator_functions (const SparseMatBaseCL<unsigned char>&);
    /// \brief Compute the connected components of the level sets of lset_->Phi. This sets component_of_dof_, ReferencePoints, and calls compute_indicator_functions. It also computes Volumes.
    void FindComponents ();

     /// \brief Compute volume of component c; for c == 0, shift must be 0.
    double CalculateVolume(Uint c, double shift) const;
    /// \brief For each component, find a point in coord_of_dof_ with largest absolute value of the level set function. Sets ReferencePoints.
    void ComputeReferencePoints();
    /// \brief For each point in refpts, return the component_of_dof from the closest point in coord_of_dof with the same sign as sign_of_component_.
    component_vector component_of_point (const std::vector<Point3DCL>& refpts, const std::vector<int>& sign_of_component, const component_vector& component_of_dof, const std::vector<Point3DCL>& coord_of_dof, const std::vector<int>& sign_of_component2) const;
    /// \brief Recomputes targetVolumes_ based on a matching of old and new components.
    void MatchComponents();

    /// \brief Copy component_of_dof_ and ReferencePoints to their _backup siblings.
    void make_backup();

    Uint num_components() const { return Volumes.size(); }

    // initializes Split, Volumes and ReferencePoints plus their backups
    void InitVolume_impl() override;

  public:
    ComponentBasedVolumeAdjustmentCL(LevelsetP2CL* lset);
    ~ComponentBasedVolumeAdjustmentCL();

    // initializes Split and Volumes, matches with the old components and recomputes ReferencePoints ... called after grid adaption
    void Repair() override;

    void AdjustVolume () override;
    void DebugOutput (std::ostream& os) const override;
    VTKVariableCL& make_VTKComponentMap(std::string var_name);
};

} // end of namespace DROPS

#endif
