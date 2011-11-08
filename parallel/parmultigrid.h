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

//******************************************************************************
// remark:                                                                     *
//  + no unknown on faces can be transfered. See doxygen!                      *
//******************************************************************************

#ifndef DROPS_PARMULTIGRID_H
#define DROPS_PARMULTIGRID_H



#include "parallel/loadbal.h"
#include "parallel/addeddata.h"
#include "parallel/interface.h"
#include "parallel/DiST.h"

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
/// This is the major class to constitute a multigrid on several procs
/// It handles the transfers and the accumulation of numerical and geometric
/// values. This class makes intensive use of the DDD-Library, which is
/// developed for UG.
/// See [],[],[] for more information and how to use this class within DROPS
///
/// This class cannot handle unknowns on faces. If you want to use them,
/// you have to change the MultiGridCL (adding an UnknownHandleCL to the FaceCL)
/// and this class by writing an HandlerFGather and HandlerFScatter.
/// \todo(of) AttachTo: fuer VecDesc und Bnd in eine Funktion packen. Geht das?
/****************************************************************************
* P A R  M U L T I  G R I D  C L A S S                                      *
****************************************************************************/
class ParMultiGridCL
{
    friend class ParMultiGridInitCL;

  public: // typedefs
      /// \brief Vector of pointers to  Vector-Describer-Classes.
      ///
      /// Used for storing the received numerical data and creating new numbering on the vertices, edges and tetras
      /// Each type of unknown has one VecDescCL* (like pressure or velocity)
    typedef std::vector<VecDescCL*>             VecDescPCT;
     /// \brief Buffer for double values
     ///
     /// Buffer type for storing double values. Used within refinement and migration
    typedef std::vector< double >               BufferCT;
     /// \brief Vector of TetraCL pointer for remembering special tetras
     ///
     /// Used for remembering Tetras, that have changed prio from ghost to master, to make MFR
     /// on edges consistent
    typedef std::vector<TetraCL*>               TetraPCT;
    /// \brief Vector of scalar boundary conditions
    typedef std::vector< const BndDataCL<double>* >    ScalBndCT;
    /// \brief Vector of vector boundary conditions
    typedef std::vector< const BndDataCL<Point3DCL>* > VecBndCT;

  private:    // variables
    MultiGridCL*     mg_;                ///< Pointer to the multigrid
    DiST::ModifyCL*  modify_;            ///< Pointer to the DiST::ModifyCL, used for Modify environment
    DiST::TransferCL*transfer_;          ///< Pointer to the DiST::TransferCL, used for Transfer environment
    int              level_;             ///< triangulation level considered for transfer
    VecDescPCT       _VecDesc;           ///< Vector of Pointers to the vector describers, where the unknowns are stored
    BufferCT         _RecvBuf;           ///< Buffer for the received numerical data (used for refinement and migration)!
    ScalBndCT        _ScalBnd;           ///< Store scalar boundary conditions
    VecBndCT         _VecBnd;            ///< Store vector boundary conditions
    VecDescCL*       _actualVec;         ///< actual index within HandleNewIdx (used within DDD-Gather and DDD-Scatter operation)
    bool             _UnkOnSimplex[3];   ///< Remember if there are unknowns on (Vertex,Edge,Tetra)
    IdxT             _RecvBufPos;        ///< last position in receive buffer that is not used
    DiST::PrioListCL priosId_;           ///< priorities used for IdentifyVertex/Edge/Face
    /// \todo _UnkOnSimplex wird nicht mehr unbedingt gebraucht. Besser direkt ueber _VecDesc abfragen.

/*
    static IFT _EdgeIF;                  // Interface for accumulate the "Marked for Refinement"
    static IFT _TetraIF;                 // Interface for communicateRefMarks() and TreatGhosts()

    static IFT           _FaceIF;           // |
    static std::ostream* _os;               //  >  for IsSane()
    static bool          _sane;             // |

    static TypeT _BndPtT,                // Boundary-Point-Type for Xfer_AddData                 (no header!)
                 _ChildPtrT;             // Childer-Pointer-Type for Xfer Children of a Tetra    (no header!)

    static IFT       NotMasterIF_;      // Vertices, edges and tetras that are not masters (for destroying unknowns after migration)

    static IFT       _GhToMaTetra_IF;   // Ghost tetras to master-tetras

    static MG_VertexContT*  _VertCont;      // Simplices stored on this proc
    static MG_EdgeContT*    _EdgeCont;
    static MG_FaceContT*    _FaceCont;
    static MG_TetraContT*   _TetraCont;



    static bool TransferMode,               // If an active DDD-Xfer-Mode is activ
                PrioChangeMode;             // if we are in a DDD-PrioChange enviroment

    static int _level;                      // all level to _level in XferStart() and XferEnd()


*/
    static ParMultiGridCL* instance_;       // only one instance of ParMultiGridCL may exist (Singleton-Pattern)

  private:
    /// \brief Create a ParMultiGridCL, due to the Singleton-Pattern, the constructor is private
    ParMultiGridCL();

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


  public:
    /// \brief Get a pointer to the ParMultiGridCL (Singleton-Pattern)
    static ParMultiGridCL* InstancePtr() { return instance_ ? instance_ : (instance_= new ParMultiGridCL); }
    /// \brief Get a reference to the ParMultiGridCL (Singleton-Pattern)
    static ParMultiGridCL& Instance()    { return *InstancePtr(); }
    ~ParMultiGridCL();


    /// \name Functions concerning the MultiGrid
    // @{
    void AttachTo(MultiGridCL&);                        // attach the MultiGrid to the ParMultiGridCL
    MultiGridCL& GetMG();                               // Get a reference on the MultiGrid
    const MultiGridCL& GetMG() const { return *mg_; }
    void Refine();                                      // Refine the MultiGrid
    void AdjustLevel();                                 // Apply all the same number of levels to all procs
    void MarkAll();                                     // All Tetras of last level are marked
    void MarkSimplicesForUnknowns();                    // Set flag "MayHaveUnknowns in the UnknownsHandleCL on all simplices that are able to store unknowns
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

    /// \name Functions concerning the handling of unknowns
    // @{
    template<typename BndT>
    void AttachTo(VecDescCL*, const BndT*);                                 // attach VecDescCL and boundary conditions
    void DeleteVecDesc();
    inline bool UnknownsOnSimplices();                               // are there unknowns on simplices
    void HandleUnknownsAfterRefine();                                       // Handle Unknowns after a refine operation
    void HandleNewIdx(IdxDescCL* old, VecDescCL* newVec);                   // Handle unknowns after migration and refinement
    /// \brief Handle unknowns after migration and refinement with a MLIdxDescCL. This is just a wrapper for the IdxDescCL
    void HandleNewIdx(MLIdxDescCL* old, VecDescCL* newVec) {
        HandleNewIdx( old->GetFinestPtr(), newVec);
    }
    void CompleteRepair(VecDescCL* newVec);                                 // After RepairAfterRefine[P1|P2] call this function
    void DelAllUnkRecv();                                                   // Clear all RecieveUnkowns-Marks
    void DeleteRecvBuffer();                                                // Delete recieve Buffer
    void DeleteUnksOnGhosts(int Level=-1);                           // Delete Unknowns on ghosts and vertical ghosts
    // @}


    /// \name Access to unknowns that are stored in a recieve buffer or an old index
    // @{
    template<typename SimplexT, typename DofT>
    struct GetDof{
      inline DofT operator() (const SimplexT&, Uint idx, int pos=-1);       // Get dof on a simplex (not implemented, see specializations down)
    };
    template<typename SimplexT>
    struct GetDof<SimplexT, Point3DCL>{
      inline Point3DCL operator() (const SimplexT&, Uint idx, int pos=-1);  // Spezialisation for vectorial data
    };
    template<typename SimplexT>
    struct GetDof<SimplexT, double>{
      inline double operator() (const SimplexT&, Uint idx, int pos=-1);     // Spezialisation for scalar data
    };
    // @}


    /// \name Checking and debug functions
    // @{
    bool IsSane(std::ostream&, int Level= -1) const;                 // Check if distributed edges and faces have the same subsimplices
    bool CheckMFR( int Level, std::ostream& os) const;                    // Check local MFR on edges.

    void DebugInfo(std::ostream&) const;                             // writes usefull infos onto outputstream
    void Show(const DiST::Helper::GeomIdCL& gid, char *mesg, int proc= -1);                // Show the simplex with a given GID
    void ConsCheck();                                                // DDD-Consisty-Check
    double GetBalance();                                             // Calculate Imbalance of triangulation over procs on last triangulation-level --> Communication
    void ShowVecDesc();
    size_t GetRecvBufferSize();                                      // Get the size of _RecvBuf
    // @}


    void AccumulateMFR(int Level=-1);                                    // accumulate mark for refinements on Edges on Level (-1==all levels!)
    void CommunicateRefMarks( Uint Level);                               // Tell everybody, if an tetra is marked for refinement
    void TreatGhosts (int Level= -1);                                    // Rescue all ghosts-subsimplices
    void RescueGhostVerts(Uint Level);                                   // Vertices are special
    void TreatHasGhosts (int Level= -1);                                 // Rescue all subsimplices of tetra that has ghosts
    void AdaptPrioOnSubs();                                              // Set prios of all subsimplices right
    void RescueSubs(TetraCL&);                                           // All subsimplices of tetra are rescued and get prio PrioMaster
    void AdaptMidVertex ();                                              // Delete or set mid-vertex for all edges

  private:
    // functions concerning the interfaces and handlers
    // ------------------------------------------------

    void VXfer(VertexCL&, PROCT, PrioT, bool del);     // transfer a vertex, this can damage the MultiGrid structure
    void EXfer(EdgeCL&, PROCT, PrioT, bool del);       // transfer an edge, this can damage the MultiGrid structure
    void FXfer(FaceCL&, PROCT, PrioT, bool del);       // not implemented

    // declare and definetypes from this class
    void DeclareBndPtT();                                    // Declare BoundaryPointer-Type
    void DeclareChildPtrT();                                 // Declare ChildPointer-Type
    void DefineBndPtT();                                     // Define BoundaryPointer-Type
    void DefineChildPtrT();                                  // Define ChildPointer-Type


    // functions concerning the internal handling of unknowns
    // ------------------------------------------------------

    // checking and size-estimating functions
    //@{
    inline bool VecDescRecv();                               // VecDesc revieved?
    Uint NumberOfUnknownsOnTetra();                          // Get overall number of unknowns on a tetrahedron
    Uint NumberOfUnknownsOnVertex();                         // Get number of unknowns on a vertex
    Uint NumberOfUnknownsOnEdge();                           // Get number of unknowns on an edge
    //@}


    // Copy values and doing linear interpolation after refine/migrate
    //@{
    inline Uint GetStorePos(const IdxDescCL*);                       // Get position where the IdxDesc is internally stored
    template<typename BndT>
    void AttachTo(const IdxDescCL*, const BndT*);                           // attach boundary conditions
    template<typename BndT>
    inline bool LinearInterpolation(const EdgeCL&, Uint, const BndT*, const VectorCL&, typename BndT::bnd_type& new_dof);
    void RescueUnknownsOnEdges();                                    // Rescue unknowns on edges, that are deleted
    void PutData(MultiGridCL::const_VertexIterator&,
                 const VectorCL* const, VectorCL*, const Uint,
                 const Uint, const IdxDescCL*);                             // Copy unknowns on a vertex into a new datafield
    template<typename BndT>
    void PutData(MultiGridCL::const_EdgeIterator&,
                 const VectorCL* const, VectorCL*, const Uint,
                 const Uint, const IdxDescCL*, const BndT*);                // Copy unknowns on an edge into a new datafield
    template<typename BndT>
    inline const BndT* GetBndCond(const IdxDescCL*);                 // Get boundary condition to store VecDesCL
    //@}

    // Send and receive unknowns
    //@{
    template<class SimplexT>
    inline void SendUnknowns(SimplexT*, TypeT, void*, int);       // Send Unknwons within the Handler<SimplexT>Gather
    template<class SimplexT>
    inline void RecvUnknowns(SimplexT*, TypeT, void*, int);       // Recieve Unknwons within the Handler<SimplexT>Gather
    void EnlargeReceiveBuffer();                                     // Enlarge the _RecvBuf
    //@}


    // set and get functions for vectors
    //@{
    template<typename SimplexT>
    void SetDof(const SimplexT&, Uint, VectorCL&, const Point3DCL&); // Put data into a given vector
    template<typename SimplexT>
    void SetDof(const SimplexT&, Uint, VectorCL&, const double&);    // Put data into a given vector

    template<typename SimplexT>
    void PutDofIntoRecvBuffer(SimplexT&, Uint, const Point3DCL&);    // Put a value of unknown into the Recieve Buffer
    template<typename SimplexT>
    void PutDofIntoRecvBuffer(SimplexT&, Uint, const double&);       // Put values of unknown into the Recieve Buffer

    template<typename SimplexT, typename ContainerT, typename DofT>
    struct GetDofOutOfVector{
      inline DofT operator() (const SimplexT&, Uint, const ContainerT&);    // Get values of unknown out of a given vector
    };
    template<typename SimplexT, typename ContainerT>
    struct GetDofOutOfVector<SimplexT, ContainerT, Point3DCL>{
      inline Point3DCL operator() (const SimplexT&, Uint, const ContainerT&);// Get values of unknown out of a given vector
    };
    template<typename SimplexT, typename ContainerT>
    struct GetDofOutOfVector<SimplexT, ContainerT, double>{
      inline double operator() (const SimplexT&, Uint, const ContainerT&);  // Get values of unknown out of a given vector
    };
    //@}


  public:
    /// \name DDD-Handler (should be private)
    //@{
    void DeleteObj(void *buffer, size_t size, int ddd_typ);          // see cpp-File for further information
    template<class SimplexT> static void HandlerDelete( OBJT);           // DDD-Handler for DynamicDataInterfaceExtraCL::XferDeleteObj

    OBJT HandlerVConstructor( size_t, PrioT, ATTRT level);  // how the reciever can construct a VertexCL
    void HandlerVXfer( OBJT , PROCT, PrioT);                // how the sender can send a VertexCL
    void HandlerVGather( OBJT, int, TypeT, void*);             // how to send the boundary-information and numerical Data
    void HandlerVScatter( OBJT, int, TypeT, void*, int);       // how to recieve the boundary-information and numerical Data

    OBJT HandlerEConstructor( size_t, PrioT, ATTRT level);  // How to construct an edge,
    void HandlerEXfer(OBJT, PROCT, PrioT);                  // to send one,
    void HandlerEGather( OBJT, int, TypeT, void*);             // to pack numerical data to message,
    void HandlerEScatter( OBJT, int, TypeT, void*, int);       // and to unpack this data

    OBJT HandlerFConstructor( size_t, PrioT, ATTRT level);  // How to construct a face
    void HandlerFXfer(OBJT, PROCT, PrioT);                  // not implemented

    OBJT HandlerTConstructor( size_t, PrioT, ATTRT level);  // How to construct a tetraeder,
    void HandlerTXfer(OBJT, PROCT, PrioT);                  // to send one
    void HandlerTGather( OBJT, int, TypeT, void*);             // to pack numerical data to message,
    void HandlerTScatter( OBJT, int, TypeT, void*, int);       // to unpack this data
    void HandlerTUpdate(OBJT);                                    // to link tetra to faces
    void HandlerTObjMkCons( OBJT, int);                           // to set right MFR values onto the edges
    void HandlerTSetPrio(OBJT, PrioT);                         // to set the priority

    int GatherEdgeMFR( OBJT, void*);                              // send MFR from Edge
    int ScatterEdgeMFR( OBJT, void*);                             // collect MFR on Edge

    int GatherTetraRestrictMarks( OBJT, void*);                   // send if a tetra is marked for refinement
    int ScatterTetraRestrictMarks( OBJT, void*);                  // recieve, if a tetra is marked for refinement

    int ExecGhostRescue( OBJT);                                   // _TetraIF calls this functions with TreatGhosts()
    int ExecGhVertRescue( OBJT);                                  // Vertices has to be treaten special!
    int ExecHasGhost( OBJT);                                      // _TetraIF calls this functions with TreatHasGhosts()
    int ExecAdaptVGhostMidVertex( OBJT);                          // _EdgeIF calls this function with AdaptMidVertex()

    int GatherEdgeSane( OBJT, void*, PROCT, ATTRT);         // Check if Edges are sane
    int ScatterEdgeSane( OBJT, void*, PROCT, ATTRT);
    int GatherFaceSane( OBJT, void*, PROCT, ATTRT);         // Check if Faces are sane
    int ScatterFaceSane( OBJT, void*, PROCT, ATTRT);

    template<class SimplexT> int DestroyUnksOnSimplex(OBJT);      // Destroy Unks on a simplex
    int GatherUnknownsRef (OBJT, void*);                          // Gather  all unknowns on tetra with subsimplices for RepairAfterRefine
    int ScatterUnknownsRef(OBJT, void*);                          // Scatter all unknowns on tetra with subsimplices for RepairAfterRefine

    int GatherUnknownsMigV (OBJT, void*);                         // Gather  all unknowns on a vertex for sending unknowns to a vertex that has changed prio from vghos/ghost to PrioHasUnk
    int ScatterUnknownsMigV(OBJT, void*);                         // Scatter all unknowns on a vertex for sending unknowns on a vertex that has changed prio from vghos/ghost to PrioHasUnk
    int GatherUnknownsMigE (OBJT, void*);                         // Gather  all unknowns on a edge for sending unknowns to a edge that has changed prio from vghos/ghost to PrioHasUnk
    int ScatterUnknownsMigE(OBJT, void*);                         // Scatter all unknowns on a egde for sending unknowns on a edge that has changed prio from vghos/ghost to PrioHasUnk

    template<typename SimplexT>
    int GatherInterpolValues (OBJT, void*);                       // Gather  unknowns of an interpolated simplex
    template<typename SimplexT>
    int ScatterInterpolValues(OBJT, void*);                       // Scatter unknowns of an interpolated simplex
    //@}
};

/// \brief Manage construction and destruction of ParMultiGridCL
class ParMultiGridInitCL
{
  public:
    ParMultiGridInitCL()  { ParMultiGridCL::Instance(); }
    ~ParMultiGridInitCL() { if (ParMultiGridCL::InstancePtr()) delete ParMultiGridCL::InstancePtr(); }
};

/// \brief Adaptor to provide access to ParMultiGridCL::Delete. Needed by MultiGridCL.
template<class SimplexT>
struct Delete_fun {
    void operator() (SimplexT& s) { ParMultiGridCL::Instance().Delete( &s); }
};

// Declaration of specialized template functions
//----------------------------------------------
template<>
void ParMultiGridCL::AttachTo<BndDataCL<double> >(const IdxDescCL*, const  BndDataCL<double>*);
template<>
void ParMultiGridCL::AttachTo<BndDataCL<Point3DCL> >(const IdxDescCL*, const  BndDataCL<Point3DCL>*);
template<>
  const BndDataCL<double>* ParMultiGridCL::GetBndCond<BndDataCL<double> >( const IdxDescCL*);
template<>
  const BndDataCL<Point3DCL>* ParMultiGridCL::GetBndCond<BndDataCL<Point3DCL> >( const IdxDescCL*);


// Declaration of DDD-Handlers used within the template functions
//---------------------------------------------------------------
extern "C" int GatherInterpolValuesVC (OBJT o, void* b);
extern "C" int ScatterInterpolValuesVC(OBJT o, void* b);
extern "C" int GatherUnknownsMigVC (OBJT o, void* b);
extern "C" int ScatterUnknownsMigVC(OBJT o, void* b);
extern "C" int GatherUnknownsMigEC (OBJT o, void* b);
extern "C" int ScatterUnknownsMigEC(OBJT o, void* b);


//  Helper classes and functions
//------------------------------

/// \brief For transfering unknowns of a killed ghost a marker has to be sent as well
///        to distinguish if this is a valid unknown
struct TransferUnkT
{
    double val;         ///< value
    Uint   idx;         ///< index
    bool   mark;        ///< marker
};

const int MIG=0, REF=1;         // constants for choosing the the filename in PrintMG

/// \brief Write Debug information into the file output/(proc_id)_MG_(REF|MIG)_(step).mg
void PrintMG(const ParMultiGridCL&, int type=MIG);

/// \brief Check parallel multigrid and write output into the file output/sane_(proc_id).chk
bool CheckParMultiGrid(const DROPS::ParMultiGridCL&);

} // end of namespace DROPS

#include "parallel/parmultigrid.tpp"     // implementation of the template or/and inline functions
#endif
