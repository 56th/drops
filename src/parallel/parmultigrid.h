/// \file parmultigrid.h
/// \brief handling of a parallel multigrid
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#ifndef DROPS_PARMULTIGRID_H
#define DROPS_PARMULTIGRID_H



#include "parallel/loadbal.h"
#include "DiST/DiST.h"

#include "misc/utils.h"
#include "misc/problem.h"

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "geom/boundary.h"
#include "geom/topo.h"
#include "num/unknowns.h"

namespace DROPS
{

class ParMultiGridInitCL; //forward declaration

/****************************************************************************
* P A R  M U L T I  G R I D  C L A S S                                      *
****************************************************************************/
/// \brief Constitute a multigrid on several procs
///
/// This is the major class to constitute a multigrid on several procs.
/// It handles the transfers and the accumulation of numerical and geometric
/// values. This class makes intensive use of the DiST-Library.
///
/****************************************************************************
* P A R  M U L T I  G R I D  C L A S S                                      *
****************************************************************************/
class ParMultiGridCL
{
    friend class ParMultiGridInitCL;
    friend class MultiGridCL;

  private:    // member variables
    MultiGridCL*     mg_;                ///< Pointer to the multigrid
    DiST::ModifyCL*  modify_;            ///< Pointer to the DiST::ModifyCL, used for Modify environment
    DiST::TransferCL*transfer_;          ///< Pointer to the DiST::TransferCL, used for Transfer environment
    int              level_;             ///< triangulation level considered for transfer
    DiST::PrioListCL priosId_;           ///< priorities used for IdentifyVertex/Edge/Face

  private:
    /// \brief Create a ParMultiGridCL, due to the Singleton-Pattern, the constructor is private
    ParMultiGridCL();
    ParMultiGridCL( const ParMultiGridCL&); // copy ctr not implemented

    /// \name Interface communication handler
    class SanityCheckCL;                    ///< checking interprocess dependencies
    class HandlerAccMFRCL;                  ///< Make AccMFR consistent at process boundaries
    class HandlerRefMarkCL;                 ///< Communicate marks for refinement from ghost to master tetrahedra
    class RescueGhostVertsCL;               ///< Rescue the ghost vertices
    /// \name Determine priority of sub-simplices
    //@{
    class TreatHasGhostsCL;
    class RescueGhostsCL;
    class RescueMasterCL;
    //@}
    class AdaptMidVertexCL;           ///< Delete or set mid-vertex of edges
    //@}

    void AttachTo(MultiGridCL&);                        // attach the MultiGrid to the ParMultiGridCL
    void AdjustLevel();                                 // Apply all the same number of levels to all procs
    void MarkSimplicesForUnknowns();                    // Set flag "MayHaveUnknowns in the UnknownsHandleCL on all simplices that are able to store unknowns

    void AccumulateMFR(int Level=-1);                                    // accumulate mark for refinements on Edges on Level (-1==all levels!)
    void CommunicateRefMarks( Uint Level);                               // Tell everybody, if an tetra is marked for refinement
    void TreatGhosts (int Level= -1);                                    // Rescue all ghosts-subsimplices
    void RescueGhostVerts(Uint Level);                                   // Vertices are special
    void TreatHasGhosts (int Level= -1);                                 // Rescue all subsimplices of tetra that has ghosts
    void AdaptPrioOnSubs();                                              // Set prios of all subsimplices right
    void RescueSubs(TetraCL&);                                           // All subsimplices of tetra are rescued and get prio PrioMaster
    void AdaptMidVertex ();                                              // Delete or set mid-vertex for all edges

  public:
    /// \brief Get a reference to the ParMultiGridCL (Singleton-Pattern)
    static ParMultiGridCL& Instance()    { static ParMultiGridCL instance; return instance; }
    /// \brief Get a pointer to the ParMultiGridCL (Singleton-Pattern)
    static ParMultiGridCL* InstancePtr() { return &Instance(); }

    /// \name Functions concerning the MultiGrid
    // @{
    MultiGridCL& GetMG();                               // Get a reference on the MultiGrid
    const MultiGridCL& GetMG() const { return *mg_; }
    void Refine();                                      // Refine the MultiGrid
    void MarkAll();                                     // All Tetras of last level are marked
    // @}

    /// \name Modify environment
    // @{
    inline void ModifyBegin();
    inline void ModifyEnd();
    template<class SimplexT>
    inline void PrioChange(SimplexT* const, Priority Prio);              // Change prio of a simplex
    template<class SimplexT>
    inline void Delete(SimplexT* const);                                 // Delete simplex from distributed multigrid.
    template<class SimplexT>
    inline void Keep(SimplexT* const);                                   // Keep simplex in distributed multigrid (even though all adjacent tetras call "Delete()")
    // @}

    /// \name Transfer environment
    // @{
    void TransferBegin(int Level=-1);                                    // Call this everytime before using a transfer (or modify) command!
    void TransferEnd();                                                  // Call this everytime after all tetras are marked for transfer
    void Transfer(TetraCL&, int proc, Priority, bool del);               // Transfer a tetra to another proc
    // @}

    /// \name Identify environment
    // @{
    void IdentifyBegin() {}
    void IdentifyEnd()   {}
    void IdentifyVertex( const EdgeCL* Parent);                         // Identify an vertex by parent edge
    void IdentifyEdge( EdgeCL* Me, const EdgeCL* Parent);               // for subedge Me of parent edge Parent
    void IdentifyEdge( EdgeCL* Me, const FaceCL* Parent);               // for subedge Me in parent face Parent
    void IdentifyFace( FaceCL* Me, const FaceCL* Parent);
    // @}



    /// \name Checking and debug functions
    // @{
    bool IsSane(std::ostream&, int Level= -1) const;                 // Check if distributed edges and faces have the same subsimplices
    bool CheckMFR( int Level, std::ostream& os) const;               // Check local MFR on edges.
    bool ConsCheck(std::ostream&) const;                             // DiST consistency check

    void DebugInfo(std::ostream&) const;                             // writes useful infos onto output stream
    void Show(const DiST::GeomIdCL& gid, char *mesg, int proc= -1) const; // Show the simplex with a given GID
    double GetBalance() const;                                       // Calculate Imbalance of triangulation over procs on last triangulation-level --> Communication
    // @}
};

/// \brief Adaptor to provide access to ParMultiGridCL::Delete. Needed by MultiGridCL.
template<class SimplexT>
struct Delete_fun {
    void operator() (SimplexT& s) { ParMultiGridCL::Instance().Delete( &s); }
};

enum { MIG=0, REF=1 };

/// \brief Write Debug information into the file output/(proc_id)_MG_(REF|MIG)_(step).mg
void PrintMG(const ParMultiGridCL&, int type=MIG);

/// \brief Check parallel multigrid and write output into the file output/sane_(proc_id).chk
bool CheckParMultiGrid();

} // end of namespace DROPS

#include "parallel/parmultigrid.tpp"     // implementation of the template or/and inline functions
#endif
