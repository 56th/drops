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

class ComponentBasedVolumeAdjustmentCL : public VolumeAdjustmentCL 
{
  private:
    std::vector<size_t>    component_of_dof_, ///< component_of_dof_[i] is the number of the component of dof i.
                           component_of_dof_backup_;
    std::vector<Point3DCL> coord_of_dof_; // The points where the dofs live.

    std::vector<bool>                  doCorrection_;
    std::vector<std::valarray<double>> indicator_functions_;
    std::vector<double>                Volumes; // volume per component
    std::vector<double>                targetVolumes; // target volume of each component
    std::vector<Point3DCL>             ReferencePoints; // markers for each connected component
    std::vector<Point3DCL>             ReferencePoints_backup; // old markers for each connected component

    /// \brief Initialize coord_of_dof_.
    void init_coord_of_dof ();

    /// \brief Renumbers the components in such a way, that the outer phase is always component 0
    void renumber_components();
    /// \brief Helper of compute_indicator_functions. Extend all components (except 0) by one level (except where they would overlap).
    std::vector<size_t> ExtendOneStep (const MatrixCL& A, const std::vector<size_t>& cp, std::vector<bool>& doCorrection) const;
    /// \brief Compute the indicator_functions_ of the extension of each connected component. Also sets doCorrection_ to false if the extensions would overlap.
    void compute_indicator_functions (const MatrixCL&);
    /// \brief Compute the connected components of the level sets of lset_->Phi. This sets component_of_dof_ and calls compute_indicator_functions.
    void FindComponents ();

    double CalculateVolume(Uint c, double shift) const; // Compute volume of component c; for c == 0, shift must be 0.
    void FindReferencePoints();
    void MatchComponents(); // uses the reference points to ensure a coherent numbering of the connected components between consecutive steps
    void make_backup();

    // Changes in Topology
    bool Handle_topo_change();


    double GetVolumeOfComponent(int i) { return Volumes[i]; }
    Uint num_components() const { return Volumes.size(); }
    Point3DCL GetReferencePoint(Uint i) {return ReferencePoints[i];}  

    // initializes Split, Volumes and ReferencePoints plus their backups
    void InitVolume_impl() override;

  public:
    ComponentBasedVolumeAdjustmentCL(LevelsetP2CL* lset);
    ~ComponentBasedVolumeAdjustmentCL();

    // initializes Split and Volumes, matches with the old components and recomputes ReferencePoints ... called after grid adaption
    void Repair() override;

    void AdjustVolume () override;
    void DebugOutput (std::ostream& os) const override;
};

} // end of namespace DROPS

#endif
