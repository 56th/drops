/// \file interfacePatch.h
/// \brief Computes 2D patches and 3D cuts of tetrahedra and interface
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Martin Horsky, Eva Loch; SC RWTH Aachen:

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

#ifndef DROPS_INTERFACEPATCH_H
#define DROPS_INTERFACEPATCH_H

#include "num/discretize.h"

namespace DROPS
{

class InterfacePatchCL
/// Computes approximation of interface.
/** Computes the planar interface patches, which are the intersection of a child T' of
 *  a tetrahedron T and the zero level of I(phi), where I(phi) is the linear interpolation
 *  of the level set function phi on T'. With LinearEdgeIntersection==false the planar interface
 *  patch on T' is determined by computing the roots of the quadratic level set function phi
 *  on each edge of T' where phi changes its sign.
 */
{
  public:
    typedef SArrayCL<BaryCoordCL,4> SubTetraT;

  protected:
    const RefRuleCL RegRef_;
    int intersec_;                    ///< number of zeros found
    int ch_;                          ///< number of the child of T (as given to GetChildData(): 0-7: regular children,
                                      ///< 8: the tetra itself, -1: class uninitialized, -2: use the subtetra st_ (defined below).
    static BaryCoordCL BaryDoF_[10];  ///< barycentric coordinates of the (P2-) levelset-dofs.
    int num_sign_[3];                 ///< 0/1/2 = -/0/+
    int sign_[10];                    ///< sign of the levelset-function in the dofs of T
    BaryCoordCL Bary_[4];             ///< barycentric coordinates of the zeros with respect to T
    int Edge_[4];                     ///< number of the cut edges (as in VertOfEdgeAr)
    Point3DCL PQRS_[4];
    int innersec_;                    ///< number of edge-intersections
    bool barysubtetra_;               ///< If true, children of the subtetra st_ of T will be considered.
                                      ///< The columns of st_ are the barycentric coordinates of the vertices of this subtetra.
                                      ///< All barycentric coordinates computed by this class are with respect to T
    SubTetraT st_;        ///< Consider the subtetra st_ of T instead of T itself.
    LocalP2CL<> PhiLoc_;  ///< levelset-function on T
    int numtriangles_;    ///< number of triangles in the intersection with the interface (0, 1, 2);
    Point3DCL Coord_[10]; ///< coordinates of the vertices and edge-barycenters of T
    bool cut_point_on_face[4][4];       ///< if (i,j)-th value is true, the i-th point of the cut is on j-th face

  private:
    static BaryCoordCL AllEdgeBaryCenter_[10][10]; ///< barycenters of all edges of the reference-tetra.
                                                    ///< \todo These fields belong into geom/topo.h or a new geom/geom.h.
    static const double approxZero_;              ///< smaller absolute values of the levelset-function are assumed to be zero.
    static const bool LinearEdgeIntersection;     ///< if true, compute the zeros of the levelset-function by linear interpolation on the edges.
                                                    ///< If false, compute the zero of the quadratic levelset-function.
  public:
      InterfacePatchCL();

      static int Sign( double phi) { return std::abs(phi)<approxZero_ ? 0 : (phi>0 ? 1 : -1); } ///< returns -1/0/1
      inline static double EdgeIntersection (Uint v0, Uint v1, const LocalP2CL<>& philoc);            ///< Compute the root of the LS-Function restricted to the edge (v0,v1) as barycentric coordinate on this edge.
      void Init( const TetraCL& t, const VecDescCL& ls, const BndDataCL<>& lsetbnd, double translation= 0.);
      void Init( const TetraCL& t, const LocalP2CL<double>& ls, double translation= 0.);
      void Init( const TetraCL& t, const SubTetraT& st, const LocalP2CL<double>& ls, double translation);
      ///< Wird nur von masstransport P1X verwendet      
      void Init( const SubTetraT& st, const LocalP2CL<double>& ls, double translation);

      /// \name Use after Init
      /// \remarks The following functions are only valid, if Init(...) was called before! They refer to T. If st_ was given to Init, they refer to the transformation of T.
      ///@{
      int GetSign( Uint DoF) const { return sign_[DoF]; }              ///< returns -1/0/1
      double GetPhi( Uint DoF) const { return PhiLoc_[DoF]; }          ///< returns value of level set function
      const Point3DCL& GetPoint( Uint i) const { return PQRS_[i]; }
      bool Intersects() const                                          ///< returns whether a patch exists on one of the 8 regular children (i.e. interface intersects tetra)
        { for(int i=1; i<10; ++i) if (sign_[0]!=sign_[i]) return true; return false; }
      bool IntersectsInterior() const                                  ///<  returns whether patch exists on one of the 8 regular children, which is not subset of a face
        { for(int i=0; i<9; ++i) for (int j=i+1; j<10; ++j) if (sign_[i]*sign_[j]==-1) return true; return false; }
      bool IntersectsChild( Uint ch) const;                            ///< returns whether a patch exists on an arbitrary child ch.
      int GetNumIntersectedSubTetras () const; ///< returns the number of intersected regular child tetrahedra

      bool ComputeVerticesOfCut( Uint ch, bool compute_PQRS= false);   ///< returns true, iff a patch exists for the child ch. If st_ was given to Init, the child of the transformed tetra T will be considered.
          ///@}

      /// \name Use after ComputeVerticesOfCut
      /// \remarks The following functions are only valid, if ComputeVerticesOfCut( ch, ...) was called before! They return information on the patch in the child ch.
      ///@{
      int                GetNumTriangles() const { return numtriangles_; } ///< returns, how many triangles form the intersection of the child and the interface.
      bool               IsQuadrilateral() const { return intersec_==4; }
      bool               EqualToFace() const { return num_sign_[1]>=3; }   ///< returns true, if patch is shared by two tetras
      Uint               GetNumPoints() const { return intersec_; }
      const BaryCoordCL& GetBary ( Uint i) const { return Bary_[i]; }      ///< The first three points are the vertices of the triangular patch;
                                                                           ///< if the patch is quadrilateral, the last three points are the vertices of the second triangle.
                                                                           ///< The barycentric coordinates are with respect to T, even if a transformation of T was applied.
      int                GetNumSign ( int sign)const { return num_sign_[sign+1]; } ///< returns number of child points with given sign, where sign is in {-1, 0, 1}

      void               WriteGeom( std::ostream&) const;                          ///< Geomview output for debugging
      void               DebugInfo( std::ostream&, bool InfoOnChild= false) const;
    ///@}
};

class InterfaceTetraCL : public InterfacePatchCL
{
  private:
    std::vector<SubTetraT> posTetras, negTetras;
    std::vector<Uint>      posChildIdx, negChildIdx;

    SubTetraT TransformToSubTetra(const SubTetraT& tetra); ///< Transforms tetra to the subtetra st_ given in InterfacePatchCL: the columns of tetra are interpreted as barycentric coordinates w.r.t. st_, the return values are the barycentric coordinates of these points w. r. t. T.
    void InsertSubTetra(const SubTetraT& BaryCoords, bool pos, Uint child); ///< modifies its first argument, if barysubtetra_ == true, because the transformation to the subtetra st_ is performed here

  public:
    bool   ComputeCutForChild( Uint ch); ///< returns true, if a patch exists for this child
    void   ComputeSubTets( Uint ch, bool clearTetras= true); ///< Computes a tetrahedralization of \f$\{\varphi<0\}\cap (\mbox{child ch of }T)\f$ and \f$\{\varphi>0\}\cap (\mbox{child ch of} T)\f$; For ch == 8 this gives the tetrahedralization of T itself. The call with clearTetras == false appends the computed tetras to posTetras and negTetras.
    void   ComputeSubTets(bool subdivide_first = true); ///< Computes a tetrahedralization of \f$\{\varphi<0\}\cap T\f$ and \f$\{\varphi>0\}\cap T\f$; the regular children of T are triangulated.
    ///@}

    /// \name Use after ComputeSubTets
    /// \remarks The following functions are only valid, if ComputeSubTets() was called before!
    ///@{
    const SubTetraT& GetTetra (Uint i)  const { return i < negTetras.size() ? negTetras[i] : posTetras[i-negTetras.size()];} ///< returns sub tetra with index \a i
    Uint  GetChildIdx         (Uint i)  const { return i < negChildIdx.size() ? negChildIdx[i] : posChildIdx[i-negChildIdx.size()];} ///< returns index of child containing sub tetra \a i
    Uint  GetNumTetra()         const {return negTetras.size() + posTetras.size();} ///< returns number of sub tetras
    Uint  GetNumNegTetra()      const {return negTetras.size();}                    ///< returns number of tetras with level set function < 0
    Uint  GetNumPosTetra()      const {return posTetras.size();}                    ///< returns number of tetras with level set function > 0
    ///@}

    /// \name Use after ComputeCutForChild
    /// \remarks The following functions are only valid, if ComputeCutForChild(...) was called before!
    ///@{
    template<class ValueT>
    ValueT quad( const LocalP2CL<ValueT>&, double absdet, bool posPart= true); ///< integrate on pos./neg. part
    template<class ValueT>
    void   quadBothParts( ValueT& int_pos, ValueT& int_neg, const LocalP2CL<ValueT>&, double absdet); ///< integrate on pos. and neg. part
    ///@}
};


class InterfaceTriangleCL : public InterfacePatchCL
{
  private:
    double          DetA_;   // determinant = 2 * area of triangle
    Point3DCL       B_[3];
    Point2DCL       ab_;

    BaryCoordCL TransformToSubTetra (const BaryCoordCL& b); ///< compute st_*b \todo remove this by introducing a column-oriented small matrix class

  public:
    bool ComputeForChild( Uint ch);                            ///< returns true, if a patch exists for this child
    double GetAbsDet( Uint tri= 0) const { return DetA_*(tri==0 ? 1.0 : GetAreaFrac()); } ///< Returns the Determinant for surface integration on the triangle \p tri.
    double GetAreaFrac()   const { return intersec_==4 ? ab_[0]+ab_[1]-1 : 0; }                   ///< Quotient of the areas of the first and the second triangle.
    template<class ValueT>
    ValueT quad2D( const LocalP2CL<ValueT>&, Uint tri= 0) const;  ///< integrate on triangle \p tri, quadrature exact up to degree 2
    const Point3DCL& GetGradId( Uint i) const { return B_[i]; }   ///< Returns the projection of the i-th standard-basis-vector of \f$R^3\f$ on the patch.    
    Point3DCL GetNormal () const;                         ///< Returns the unit normal to the linear approximation of \f$\Gamma\f$, that points from \f$\{\varphi<0\}\f$ to \f$\{\varphi<0\}\f$.
    Quad5_2DCL<Point3DCL> GetImprovedNormal(Uint) const;  ///< Returns the improved unit normal
    Point3DCL ApplyProj( const Point3DCL& grad) const { return grad[0]*B_[0] + grad[1]*B_[1] + grad[2]*B_[2]; }
};

//***********************************************************************************************************
// InterfaceLineCL is used to handle the cases in which the interface intersects with slip/symmetry boundary.
// The number of contact line segments are computed and the information of the end nodes of the segments are
// computed, too. In addition, the normal vectors at the contact line are computed.
//***********************************************************************************************************
class InterfaceLineCL : public InterfacePatchCL
{
  private:
    Uint     numMCL_;	                   ///< number of moving contact lines (MCL)
    Uint     IdxMCL_[4][2];                ///< the edge index for each contact line
    BndCondT  BC_Face_[4];                 ///< boundary condition type for all four faces
    BndCondT  BC_Edge_[6];                 ///< boundary condition type for all six edges
    instat_vector_fun_ptr outnormal_;      ///< the outer normal of the (slip) boundary
    bool SymmType[4];                      ///< store if a contact line segment is symmetric

  public:
    bool ComputeMCLForChild(Uint ch);                          ///< returns true, if a moving contact line exists for this child
    Uint GetNumMCL();                                          ///< returns number of MCL segments
    void SetBndCondT(const TetraCL& tet, const BndDataCL<Point3DCL>& BndData);     ///< set the boundary condition type of the tetra
    void SetBndOutNormal(instat_vector_fun_ptr outnormal);     ///< set the outer normal of the slip boundary
    bool IsSymmType(Uint i) {return SymmType[i];}              ///<return if a contact line segment is on symmetry boundary
    double GetInfoMCL(Uint v, BaryCoordCL& bary0, BaryCoordCL& bary1, Point3DCL& pt0, Point3DCL& pt1); ///< return the length of the MCL
                                                                                                       /// set the BaryCoord and Point3D of two end nodes
    Quad9_1DCL<Point3DCL> GetImprovedNormalAtMCL(Uint v) const;  ///< Returns the improved unit normal of the interface at the moving contact line
    Quad9_1DCL<Point3DCL> GetImprovedMCLNormalOnSlipBnd(const TetraCL& t, Uint v) const;    ///< Returns the unit outer normal at the contact line which lies on the slip boundary                                                           ///computed in an improved way using levelset function
                                                                            /// v denotes the v-th contact line. 
                                                                            /// It must be called after SetBndOutNormal() 
    Quad9_1DCL<double> GetDynamicCtAngle(const TetraCL& t, Uint v) const;   ///< Returns the contact angle

};


LocalP2CL<double> ProjectIsoP2ChildToParentP1 (LocalP2CL<double> lpin, Uint child);


} // end of namespace DROPS

#include "num/interfacePatch.tpp"

#endif
