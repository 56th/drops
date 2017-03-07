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



// forward declarations
class GraphComponentsCL;



class ComponentBasedVolumeAdjustmentCL : public VolumeAdjustmentCL 
{
  private:
      GraphComponentsCL& Split; // determines connected components and numbers them
      std::valarray<double> Volumes; // volume per component
      std::valarray<double> Volumes_backup; // old volumes per component
      std::vector<Point3DCL> ReferencePoints; // markers for each connected component
      std::vector<Point3DCL> ReferencePoints_backup; // old markers for each connected component
      
      double ComputeComponentAdjustment (int compnumber);
      void CalculateInitialVolumes();
      void CalculateVolumes(std::valarray<double>& volumes, std::valarray<double>* epsilons) const;
      void FindReferencePoints();
      void MatchComponents(); // uses the reference points to ensure a coherent numbering of the connected components between consecutive steps
      void make_backup(bool complete=false);
      
      // Changes in Topology
      bool Handle_topo_change();
      
      
      double GetVolumeOfComponent(int i) {return Volumes[i];}
      int GetNumberOfComponents() const;
      Point3DCL GetReferencePoint(Uint i) {return ReferencePoints[i];}  
      GraphComponentsCL& GetSplit();
      
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
