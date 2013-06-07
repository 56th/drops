/// \file fe.tpp
/// \brief Description of various finite-element functions
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt ; SC RWTH Aachen:

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

/// content of the file (longer version)

/// \file
/// \brief Description of various finite-element functions

#ifdef _PAR
#  include "parallel/parallel.h"
#  include "parallel/parmultigrid.h"
#endif

namespace DROPS
{


//**************************************************************************
// Class:   FE_P1CL                                                        *
//**************************************************************************
inline double
FE_P1CL::H(Uint dof, double v1, double v2, double v3)
{
    switch(dof) {
      case 0: return H0( v1, v2, v3);
      case 1: return H1( v1, v2, v3);
      case 2: return H2( v1, v2, v3);
      case 3: return H3( v1, v2, v3);
      default: throw DROPSErrCL("FE_P1CL::H: Invalid shape function.");
    };
}


//**************************************************************************
// Class:   FE_P1DCL                                                       *
//**************************************************************************
inline double
FE_P1DCL::H(Uint dof, double v1, double v2, double v3)
{
    switch(dof) {
      case 0: return H0( v1, v2, v3);
      case 1: return H1( v1, v2, v3);
      case 2: return H2( v1, v2, v3);
      case 3: return H3( v1, v2, v3);
      default: throw DROPSErrCL("FE_P1DCL::H: Invalid shape function.");
    };
}

inline double
FE_P1DCL::H(Uint dof, const BaryCoordCL& p)
{
    switch(dof) {
      case 0: return H0( p);
      case 1: return H1( p);
      case 2: return H2( p);
      case 3: return H3( p);
      default: throw DROPSErrCL("FE_P1DCL::H: Invalid shape function.");
    };
}


//**************************************************************************
// Class:   FE_P2CL                                                        *
//**************************************************************************
inline double
FE_P2CL::H(Uint dof, double v1, double v2, double v3)
{
    switch(dof) {
      case 0: return H0( v1, v2, v3);
      case 1: return H1( v1, v2, v3);
      case 2: return H2( v1, v2, v3);
      case 3: return H3( v1, v2, v3);
      case 4: return H4( v1, v2, v3);
      case 5: return H5( v1, v2, v3);
      case 6: return H6( v1, v2, v3);
      case 7: return H7( v1, v2, v3);
      case 8: return H8( v1, v2, v3);
      case 9: return H9( v1, v2, v3);
      default: throw DROPSErrCL("FE_P2CL::H: Invalid shape function.");
    };
}

inline double
FE_P2CL::H(Uint dof, const BaryCoordCL& p)
{
    switch(dof) {
      case 0: return H0( p);
      case 1: return H1( p);
      case 2: return H2( p);
      case 3: return H3( p);
      case 4: return H4( p);
      case 5: return H5( p);
      case 6: return H6( p);
      case 7: return H7( p);
      case 8: return H8( p);
      case 9: return H9( p);
      default: throw DROPSErrCL("FE_P2CL::H: Invalid shape function.");
    };
}

inline SVectorCL<3>
FE_P2CL::DH0Ref(double v1, double v2, double v3)
{
    SVectorCL<3> ret(4.*(v1 + v2 + v3) -3.);
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DH1Ref(double v1, double, double)
{
    SVectorCL<3> ret;
    ret[0]= 4.*v1 -1.; ret[1]= ret[2]= 0.;
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DH2Ref(double, double v2, double)
{
    SVectorCL<3> ret;
    ret[0]= ret[2]= 0.; ret[1]= 4*v2 -1.;
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DH3Ref(double, double, double v3)
{
    SVectorCL<3> ret;
    ret[0]= 0.; ret[1]= 0.; ret[2]= 4.*v3 -1.;
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DH4Ref(double v1, double v2, double v3)
{
    SVectorCL<3> ret;
    ret[0]= 4.*( 1. - (2.*v1 + v2 + v3) ); ret[1]= ret[2]= -4.*v1;
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DH5Ref(double v1, double v2, double v3)
{
    SVectorCL<3> ret;
    ret[0]= ret[2]= -4.*v2; ret[1]= 4.*( 1. -(2.*v2 + v1 + v3) );
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DH6Ref(double v1, double v2, double)
{
    SVectorCL<3> ret;
    ret[0]= 4.*v2; ret[1]= 4.*v1; ret[2]= 0.;
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DH7Ref(double v1, double v2, double v3)
{
    SVectorCL<3> ret;
    ret[0]= ret[1]= -4.*v3; ret[2]= 4.*( 1. -(2.*v3 + v1 + v2) );
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DH8Ref(double v1, double, double v3)
{
    SVectorCL<3> ret;
    ret[0]= 4.*v3; ret[1]= 0.; ret[2]= 4.*v1;
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DH9Ref(double, double v2, double v3)
{
    SVectorCL<3> ret;
    ret[0]= 0.; ret[1]= 4.*v3; ret[2]= 4.*v2;
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DHRef(Uint dof, double v1, double v2, double v3)
{
    switch (dof)
    {
      case 0: return DH0Ref(v1, v2, v3);
      case 1: return DH1Ref(v1, v2, v3);
      case 2: return DH2Ref(v1, v2, v3);
      case 3: return DH3Ref(v1, v2, v3);
      case 4: return DH4Ref(v1, v2, v3);
      case 5: return DH5Ref(v1, v2, v3);
      case 6: return DH6Ref(v1, v2, v3);
      case 7: return DH7Ref(v1, v2, v3);
      case 8: return DH8Ref(v1, v2, v3);
      case 9: return DH9Ref(v1, v2, v3);
      default: throw DROPSErrCL("FE_P2CL::DHRef: Invalid shape function.");
    };
}

inline double
FE_P2CL::Laplace(Uint dof, const SMatrixCL<3,3>& M)
{
    double ret= 0.;
    for (Uint i=0; i<3; ++i)
      for(Uint j=0; j<3; ++j)
        for (Uint k=0; k<3; ++k)
        {
            ret+= M(i,j)*M(i,k)*D2HRef(dof, j, k);
        }
    return ret;
}


//**************************************************************************
// Class:   P1EvalCL                                                       *
//**************************************************************************
template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P1EvalCL<Data, _BndData, _VD>::GetDoF(const VertexCL& s, _Cont& c) const
{
    c[0]= GetDoF( s);
}

template<class Data, class _BndData, class _VD>
  inline Data
  P1EvalCL<Data, _BndData, _VD>::val(const VertexCL& s) const
{
    return GetDoF( s);
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline Data
P1EvalCL<Data, _BndData, _VD>::val(const _Cont& c) const
{
    return  c[0];
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P1EvalCL<Data, _BndData, _VD>::GetDoF(const EdgeCL& s, _Cont& c) const
{
    c[0]= GetDoF( *s.GetVertex( 0));
    c[1]= GetDoF( *s.GetVertex( 1));
}

template<class Data, class _BndData, class _VD>
  inline Data
  P1EvalCL<Data, _BndData, _VD>::val(const EdgeCL& s, double v1) const
{
    return  GetDoF( *s.GetVertex( 0))*FE_P1CL::H0( v1)
          + GetDoF( *s.GetVertex( 1))*FE_P1CL::H1( v1);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P1EvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1) const
{
    return c[0]*FE_P1CL::H0( v1) + c[1]*FE_P1CL::H1( v1);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P1EvalCL<Data, _BndData, _VD>::GetDoF(const TetraCL& s, _Cont& c) const
{
    for (Uint i= 0; i < NumVertsC; ++i)
        c[i]= GetDoF( *s.GetVertex( i));
}

template<class Data, class _BndData, class _VD>
  inline Data
  P1EvalCL<Data, _BndData, _VD>::val(const TetraCL& s, double v1, double v2, double v3) const
{
    return  GetDoF( *s.GetVertex( 0))*FE_P1CL::H0( v1, v2, v3)
           +GetDoF( *s.GetVertex( 1))*FE_P1CL::H1( v1, v2, v3)
           +GetDoF( *s.GetVertex( 2))*FE_P1CL::H2( v1, v2, v3)
           +GetDoF( *s.GetVertex( 3))*FE_P1CL::H3( v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
    inline Data
    P1EvalCL<Data, _BndData, _VD>::val(const TetraCL& s, const BaryCoordCL& p) const
{
    return  GetDoF( *s.GetVertex( 0))*FE_P1CL::H0( p)
           +GetDoF( *s.GetVertex( 1))*FE_P1CL::H1( p)
           +GetDoF( *s.GetVertex( 2))*FE_P1CL::H2( p)
           +GetDoF( *s.GetVertex( 3))*FE_P1CL::H3( p);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P1EvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1, double v2, double v3) const
{
    return c[0] * FE_P1CL::H0(v1, v2, v3) + c[1] * FE_P1CL::H1(v1, v2, v3)
         + c[2] * FE_P1CL::H2(v1, v2, v3) + c[3] * FE_P1CL::H3(v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
  inline void
  P1EvalCL<Data, _BndData, _VD>::SetDoF(const VertexCL& s, const Data& d)
{
    Assert( !_bnd->IsOnDirBnd(s),
        DROPSErrCL("P1EvalBaseCL::SetDoF: Trying to assign to"
        "Dirichlet-boundary-vertex."), DebugNumericC);
    DoFHelperCL<Data, typename VecDescT::DataType>::set(
        _sol->Data, s.Unknowns(_sol->RowIdx->GetIdx()), d);
}


//**************************************************************************
// Class:   P1DEvalCL                                                      *
//**************************************************************************
template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P1DEvalCL<Data, _BndData, _VD>::GetDoF(const TetraCL& s, _Cont& c) const
{
    for (Uint i= 0; i < NumFacesC; ++i)
        c[i]= GetDoF( *s.GetFace( i));
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P1DEvalCL<Data, _BndData, _VD>::GetDoFP1(const TetraCL& s, _Cont& c) const
{
        c[0]= val( s, 0., 0., 0.);
        c[1]= val( s, 1., 0., 0.);
        c[2]= val( s, 0., 1., 0.);
        c[3]= val( s, 0., 0., 1.);
}

template<class Data, class _BndData, class _VD>
  inline Data
  P1DEvalCL<Data, _BndData, _VD>::val(const TetraCL& s, double v1, double v2, double v3) const
{
    return  GetDoF( *s.GetFace( 0))*FE_P1DCL::H0( v1, v2, v3)
           +GetDoF( *s.GetFace( 1))*FE_P1DCL::H1( v1, v2, v3)
           +GetDoF( *s.GetFace( 2))*FE_P1DCL::H2( v1, v2, v3)
           +GetDoF( *s.GetFace( 3))*FE_P1DCL::H3( v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P1DEvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1, double v2, double v3) const
{
    return c[0] * FE_P1DCL::H0(v1, v2, v3) + c[1] * FE_P1DCL::H1(v1, v2, v3)
         + c[2] * FE_P1DCL::H2(v1, v2, v3) + c[3] * FE_P1DCL::H3(v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
  inline void
  P1DEvalCL<Data, _BndData, _VD>::SetDoF(const FaceCL& s, const Data& d)
{
    Assert( !_bnd->IsOnDirBnd(s),
        DROPSErrCL("P1EvalBaseCL::SetDoF: Trying to assign to"
        "Dirichlet-boundary-vertex."), DebugNumericC);
    DoFHelperCL<Data, typename VecDescT::DataType>::set(
        _sol->Data, s.Unknowns(_sol->RowIdx->GetIdx()), d);
}


//**************************************************************************
// Class:   P2EvalCL                                                       *
//**************************************************************************
template<class Data, class _BndData, class _VD>
  inline void
  P2EvalCL<Data, _BndData, _VD>::SetDoF(const VertexCL& s, const Data& d)
{
    Assert( !_bnd->IsOnDirBnd(s),
            DROPSErrCL( "P2EvalBaseCL::SetDoF: Trying to assign to Dirichlet-boundary-vertex."),
            DebugNumericC);
    DoFHelperCL<Data, typename VecDescT::DataType>::set(
        _sol->Data, s.Unknowns( _sol->RowIdx->GetIdx()), d);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P2EvalCL<Data, _BndData, _VD>::GetDoF(const VertexCL& s, _Cont& c) const
{
    c[0]= GetDoF( s);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P2EvalCL<Data, _BndData, _VD>::val(const _Cont& c) const
{
    return c[0];
}

template<class Data, class _BndData, class _VD>
  inline Data
  P2EvalCL<Data, _BndData, _VD>::val(const VertexCL& s) const
{
    return GetDoF( s);
}


template<class Data, class _BndData, class _VD>
  inline void
  P2EvalCL<Data, _BndData, _VD>::SetDoF(const EdgeCL& s, const Data& d)
{
    Assert( !_bnd->IsOnDirBnd(s),
            DROPSErrCL( "P2EvalBaseCL::SetDoF: Trying to assign to Dirichlet-boundary-edge."),
            DebugNumericC);
    DoFHelperCL<Data, typename VecDescT::DataType>::set(
        _sol->Data, s.Unknowns( _sol->RowIdx->GetIdx()), d);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P2EvalCL<Data, _BndData, _VD>::GetDoF(const EdgeCL& s, _Cont& c) const
{
    c[0]= GetDoF( *s.GetVertex( 0));
    c[1]= GetDoF( *s.GetVertex( 1));
    c[2]= GetDoF( s);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P2EvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1) const
{
    return c[0] * FE_P2CL::H0( v1) + c[1] * FE_P2CL::H1( v1) + c[2] * FE_P2CL::H2( v1);
}

template<class Data, class _BndData, class _VD>
  inline Data
  P2EvalCL<Data, _BndData, _VD>::val(const EdgeCL& s, double v1) const
{
    return  GetDoF( *s.GetVertex(0))*FE_P2CL::H0( v1)
          + GetDoF( *s.GetVertex(1))*FE_P2CL::H1( v1)
          + GetDoF( s)*FE_P2CL::H2( v1);
}

template<class Data, class _BndData, class _VD>
  inline Data
  P2EvalCL<Data, _BndData, _VD>::val(const EdgeCL& s) const
{
    return GetDoF( s);
}


template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P2EvalCL<Data, _BndData, _VD>::GetDoF(const TetraCL& s, Uint face, _Cont& c) const
{
    for(Uint i= 0; i < 3; ++i) {
        c[i]= GetDoF( *s.GetVertex( VertOfFace( face, i)));
        c[i+3]= GetDoF( *s.GetEdge( EdgeOfFace( face, i)));
    }
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P2EvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1, double v2) const
{
    return c[0] * FE_P2CL::H0(v1, v2) + c[1] * FE_P2CL::H1(v1, v2)
         + c[2] * FE_P2CL::H2(v1, v2) + c[3] * FE_P2CL::H3(v1, v2)
         + c[4] * FE_P2CL::H4(v1, v2) + c[5] * FE_P2CL::H5(v1, v2);
}


template<class Data, class _BndData, class _VD>
  inline Data
  P2EvalCL<Data, _BndData, _VD>::val(const TetraCL& s, Uint face, double v1, double v2) const
{
    return  GetDoF( *s.GetVertex( VertOfFace( face, 0)))*FE_P2CL::H0( v1, v2)
           +GetDoF( *s.GetVertex( VertOfFace( face, 1)))*FE_P2CL::H1( v1, v2)
           +GetDoF( *s.GetVertex( VertOfFace( face, 2)))*FE_P2CL::H2( v1, v2)
           +GetDoF( *s.GetEdge( EdgeOfFace( face, 0)))*FE_P2CL::H3( v1, v2)
           +GetDoF( *s.GetEdge( EdgeOfFace( face, 1)))*FE_P2CL::H4( v1, v2)
           +GetDoF( *s.GetEdge( EdgeOfFace( face, 2)))*FE_P2CL::H5( v1, v2);
}


template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P2EvalCL<Data, _BndData, _VD>::GetDoF(const TetraCL& s, _Cont& c) const
{
    for (Uint i= 0; i < NumVertsC; ++i)
        c[i]= GetDoF( *s.GetVertex( i));
    for (Uint i= 0; i < NumEdgesC; ++i)
        c[i+NumVertsC]= GetDoF( *s.GetEdge( i));
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P2EvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1, double v2, double v3) const
{
    return  c[0] * FE_P2CL::H0( v1, v2, v3) + c[1] * FE_P2CL::H1( v1, v2, v3)
          + c[2] * FE_P2CL::H2( v1, v2, v3) + c[3] * FE_P2CL::H3 (v1, v2, v3)
          + c[4] * FE_P2CL::H4( v1, v2, v3) + c[5] * FE_P2CL::H5( v1, v2, v3)
          + c[6] * FE_P2CL::H6( v1, v2, v3) + c[7] * FE_P2CL::H7( v1, v2, v3)
          + c[8] * FE_P2CL::H8( v1, v2, v3) + c[9] * FE_P2CL::H9( v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
  inline Data
  P2EvalCL<Data, _BndData, _VD>::val(const TetraCL& s, double v1, double v2, double v3) const
{
    return  GetDoF( *s.GetVertex( 0))*FE_P2CL::H0( v1, v2, v3)
           +GetDoF( *s.GetVertex( 1))*FE_P2CL::H1( v1, v2, v3)
           +GetDoF( *s.GetVertex( 2))*FE_P2CL::H2( v1, v2, v3)
           +GetDoF( *s.GetVertex( 3))*FE_P2CL::H3( v1, v2, v3)
           +GetDoF( *s.GetEdge( 0))*FE_P2CL::H4( v1, v2, v3)
           +GetDoF( *s.GetEdge( 1))*FE_P2CL::H5( v1, v2, v3)
           +GetDoF( *s.GetEdge( 2))*FE_P2CL::H6( v1, v2, v3)
           +GetDoF( *s.GetEdge( 3))*FE_P2CL::H7( v1, v2, v3)
           +GetDoF( *s.GetEdge( 4))*FE_P2CL::H8( v1, v2, v3)
           +GetDoF( *s.GetEdge( 5))*FE_P2CL::H9( v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
  inline Data
  P2EvalCL<Data, _BndData, _VD>::val(const TetraCL& s, const BaryCoordCL& p) const
{
    return  GetDoF( *s.GetVertex( 0))*FE_P2CL::H0( p)
           +GetDoF( *s.GetVertex( 1))*FE_P2CL::H1( p)
           +GetDoF( *s.GetVertex( 2))*FE_P2CL::H2( p)
           +GetDoF( *s.GetVertex( 3))*FE_P2CL::H3( p)
           +GetDoF( *s.GetEdge( 0))*FE_P2CL::H4( p)
           +GetDoF( *s.GetEdge( 1))*FE_P2CL::H5( p)
           +GetDoF( *s.GetEdge( 2))*FE_P2CL::H6( p)
           +GetDoF( *s.GetEdge( 3))*FE_P2CL::H7( p)
           +GetDoF( *s.GetEdge( 4))*FE_P2CL::H8( p)
           +GetDoF( *s.GetEdge( 5))*FE_P2CL::H9( p);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P2EvalCL<Data, _BndData, _VD>::val(const _Cont& c, const BaryCoordCL& p)
{
    return  c[0] * FE_P2CL::H0( p) + c[1] * FE_P2CL::H1( p)
          + c[2] * FE_P2CL::H2( p) + c[3] * FE_P2CL::H3( p)
          + c[4] * FE_P2CL::H4( p) + c[5] * FE_P2CL::H5( p)
          + c[6] * FE_P2CL::H6( p) + c[7] * FE_P2CL::H7( p)
          + c[8] * FE_P2CL::H8( p) + c[9] * FE_P2CL::H9( p);
}

template<class Data, class _BndData, class _VD>
inline bool P2EvalCL<Data, _BndData, _VD>::UnknownsMissing(const TetraCL& t) const
{
    const Uint idx= _sol->RowIdx->GetIdx();
    for (TetraCL::const_VertexPIterator it= t.GetVertBegin(), end= t.GetVertEnd();
         it!=end; ++it)
        if ( !IsDefinedOn( **it)) return true;
    for (TetraCL::const_EdgePIterator it= t.GetEdgesBegin(), end= t.GetEdgesEnd();
         it!=end; ++it)
        if ( !(_bnd->IsOnDirBnd( **it) || (*it)->Unknowns.Exist(idx) )) return true;
    return false;
}


template<class Data, class _BndData, class _VD>
inline bool P2EvalCL<Data, _BndData, _VD>::IsDefinedOn(const VertexCL& v) const
{
    return _bnd->IsOnDirBnd( v)
           || (v.Unknowns.Exist() && v.Unknowns.Exist( _sol->RowIdx->GetIdx()));
}

template<class Data, class _BndData, class _VD>
inline bool P2EvalCL<Data, _BndData, _VD>::IsDefinedOn(const EdgeCL& e) const
{
    return IsDefinedOn( *e.GetVertex( 0)) && IsDefinedOn( *e.GetVertex( 1))
           && (_bnd->IsOnDirBnd( e)
               || (e.Unknowns.Exist() && e.Unknowns.Exist( _sol->RowIdx->GetIdx())));
}

template<class Data, class _BndData, class _VD>
inline bool P2EvalCL<Data, _BndData, _VD>::IsDefinedOn(
    const TetraCL& t, const Uint face) const
{
    const Uint idx= _sol->RowIdx->GetIdx();
    for (Uint i=0; i<3; ++i) {
        const VertexCL& v= *t.GetVertex( VertOfFace( face, i));
        if (!IsDefinedOn( v)) return false;
        const EdgeCL* const ep= t.GetEdge( EdgeOfFace( face, i));
        if (!(_bnd->IsOnDirBnd( *ep)
              || (ep->Unknowns.Exist() && ep->Unknowns.Exist( idx))))
            return false;
    }
    return true;
}


template<class Data, class _BndData, class _VD>
inline bool P2EvalCL<Data, _BndData, _VD>::IsDefinedOn(const TetraCL& t) const
{
    for (Uint i=0; i<NumVertsC; ++i)
        if (!IsDefinedOn( *t.GetVertex( i))) return false;
    const Uint idx= _sol->RowIdx->GetIdx();
    for (Uint i=0; i<NumEdgesC; ++i) {
        const EdgeCL* const ep= t.GetEdge( i);
        if (!(_bnd->IsOnDirBnd( *ep)
              || (ep->Unknowns.Exist() && ep->Unknowns.Exist( idx))))
            return false;
    }
    return true;
}


//**************************************************************************
// Class:   P1BubbleEvalCL                                                 *
//**************************************************************************
template<class Data, class _BndData, class _VD>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::GetDoF(const VertexCL& s) const
{
    return _bnd->IsOnDirBnd(s) ? _bnd->GetDirBndValue(s)
        : DoFHelperCL<Data,typename VecDescT::DataType>::get(
        _sol->Data, s.Unknowns(_sol->RowIdx->GetIdx()));
}


template<class Data, class _BndData, class _VD> template<class _Cont>
inline void
P1BubbleEvalCL<Data, _BndData, _VD>::GetDoF(const VertexCL& s, _Cont& c) const
{
    c[0]= GetDoF(s);
}


template<class Data, class _BndData, class _VD>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::GetDoF(const TetraCL& s) const
{
    return DoFHelperCL<Data,typename VecDescT::DataType>::get(
        _sol->Data, s.Unknowns(_sol->RowIdx->GetIdx()));
}


template<class Data, class _BndData, class _VD>
inline void
P1BubbleEvalCL<Data, _BndData, _VD>::SetDoF(const VertexCL& s, const Data& d)
{
    Assert( !_bnd->IsOnDirBnd(s),
        DROPSErrCL("P1EvalCL::SetDoF: Trying to assign to Dirichlet-boundary-vertex."),
        DebugNumericC);
    DoFHelperCL<Data, typename VecDescT::DataType>::set(
        _sol->Data, s.Unknowns(_sol->RowIdx->GetIdx()), d);
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::val(const _Cont& c) const
{
    return  c[0];
}

template<class Data, class _BndData, class _VD>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::val(const VertexCL& s) const
{
    return GetDoF(s);
}


template<class Data, class _BndData, class _VD> template<class _Cont>
inline void
P1BubbleEvalCL<Data, _BndData, _VD>::GetDoF(const EdgeCL& s, _Cont& c) const
{
    c[0]= GetDoF(*s.GetVertex(0));
    c[1]= GetDoF(*s.GetVertex(1));
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1) const
{
    return c[0] * FE_P1BubbleCL::H0(v1) + c[1] * FE_P1BubbleCL::H1(v1);
}

template<class Data, class _BndData, class _VD>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::val(const EdgeCL& s, double v1) const
{
    return  GetDoF(*s.GetVertex(0))*FE_P1BubbleCL::H0(v1)
          + GetDoF(*s.GetVertex(1))*FE_P1BubbleCL::H1(v1);
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline void
P1BubbleEvalCL<Data, _BndData, _VD>::GetDoF(const TetraCL& s, _Cont& c) const
{
    for(Uint i= 0; i < NumVertsC; ++i)
        c[i]= GetDoF(*s.GetVertex(0));
    c[4]= GetDoF(s);
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1, double v2, double v3) const
{
    return c[0] * FE_P1BubbleCL::H0(v1, v2, v3) + c[1] * FE_P1BubbleCL::H1(v1, v2, v3)
         + c[2] * FE_P1BubbleCL::H2(v1, v2, v3) + c[3] * FE_P1BubbleCL::H3(v1, v2, v3)
         + c[4] * FE_P1BubbleCL::H4(v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::val(const TetraCL& s, double v1, double v2, double v3) const
{
    return  GetDoF(*s.GetVertex(0))*FE_P1BubbleCL::H0(v1, v2, v3)
           +GetDoF(*s.GetVertex(1))*FE_P1BubbleCL::H1(v1, v2, v3)
           +GetDoF(*s.GetVertex(2))*FE_P1BubbleCL::H2(v1, v2, v3)
           +GetDoF(*s.GetVertex(3))*FE_P1BubbleCL::H3(v1, v2, v3)
           +GetDoF(s)*FE_P1BubbleCL::H4(v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::lin_val(const TetraCL& s, double v1, double v2, double v3) const
{
    return  GetDoF(*s.GetVertex(0))*FE_P1BubbleCL::H0(v1, v2, v3)
           +GetDoF(*s.GetVertex(1))*FE_P1BubbleCL::H1(v1, v2, v3)
           +GetDoF(*s.GetVertex(2))*FE_P1BubbleCL::H2(v1, v2, v3)
           +GetDoF(*s.GetVertex(3))*FE_P1BubbleCL::H3(v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::bubble_val(const TetraCL& s, double v1, double v2, double v3) const
{
    return  GetDoF(s)*FE_P1BubbleCL::H4(v1, v2, v3);
}


template<class Data, class _BndData, class _VD>
void
Interpolate(P1EvalCL<Data, _BndData, _VD>& sol, const P1EvalCL<Data, _BndData, const _VD>& old_sol)
// This only works, if Interpolate is called after every refinement of the multigrid.
// Take care, that x and old_x are on successive triangulations.
{
#ifdef _PAR
    throw DROPSErrCL("InterpolateP1: The functionality of interpolating is done in RepairP1 in parallel!");
#endif
    typedef typename P1EvalCL<Data, _BndData, _VD>::BndDataCL BndCL;
    const BndCL* const _bnd= old_sol.GetBndData();

    const MultiGridCL& _MG= old_sol.GetMG();
    const Uint old_level= old_sol.GetLevel();
    //const Uint level= sol.GetLevel();
    const Uint old_idx= old_sol.GetSolution()->RowIdx->GetIdx();
    Uint counter1= 0, counter2= 0;

//    Assert( level==old_level || old_level==level-1, DROPSErrCL("Interpolate: successive triangs are expected\n"));
    // Iterate over all edges, interpolate values on new mid vertices
    for (MultiGridCL::const_EdgeIterator sit= _MG.GetAllEdgeBegin(old_level), theend= _MG.GetAllEdgeEnd(old_level);
         sit!=theend; ++sit)
        if ( sit->IsRefined() && !_bnd->IsOnDirBnd(*sit->GetMidVertex())  ) // only new non-boundary vertices are interpolated
        {
            sol.SetDoF(*sit->GetMidVertex(), (old_sol.val(*sit->GetVertex(0)) + old_sol.val(*sit->GetVertex(1)))/2.);
            ++counter2;
        }
    // Iterate over the vertices of the old triangulation and copy values
    for (MultiGridCL::const_TriangVertexIteratorCL sit= _MG.GetTriangVertexBegin(old_level), theend= _MG.GetTriangVertexEnd(old_level);
         sit!=theend; ++sit)
        if ( !_bnd->IsOnDirBnd(*sit) )
        {
// @@@ Handle missing unknowns: Do it right here?
            if (sit->Unknowns.Exist(old_idx))
                sol.SetDoF(*sit, old_sol.val(*sit));
//else they were set in the for-loop before !
            ++counter1;
        }

    std::cout << "Interpolate: " << counter1 << " vertex-DoF of "
              << old_sol.GetSolution()->Data.size() << " copied, "
              << counter2 << " new Mid-vertex-DoF interpolated." << std::endl;
}

template<class Data, class _BndData, class _VD>
void InterpolateChildren( const TetraCL& t, P2EvalCL<Data, _BndData, _VD>& sol, const P2EvalCL<Data, _BndData, const _VD>& old_sol)
{
    typedef typename P2EvalCL<Data, _BndData, _VD>::BndDataCL BndCL;
    const BndCL* const _bnd= old_sol.GetBndData();

    const double edgebary[3][2]=
        { {0.25, 0.25},
          {0.5 , 0.25},
          {0.25, 0.5}
        };

    // Hole des Tetraeders RefRule; gehe ueber alle Kinder ;-): durchlaufe die edges in
    // der childrule des kindes: falls IsSubEdge(edge): finde mit ParentEdge & NumOfSubEdge heraus,
    // von welcher Kante des tetras, (falls der wert bei edge != 0) und interpoliere ueber kante.
    // sonst, falls IsSubInParFace(subedge): WhichEdgeInFace(subedge, face, pos) und interpoliere ueber face;
    // sonst, behandele raumdiagonale;
    const RefRuleCL& refrule= t.GetRefData();
    TetraCL::const_ChildPIterator child= t.GetChildBegin();
    const TetraCL::const_ChildPIterator childend= t.GetChildEnd();
    for (Uint childnum=0; child!=childend; ++childnum, ++child)
    {
        const ChildDataCL& childdata= GetChildData(refrule.Children[childnum]);
        for (Uint chedge=0; chedge<NumEdgesC; ++chedge)
        {
            const EdgeCL* const edgep= (*child)->GetEdge(chedge);
            const Uint chedgeinparent= childdata.Edges[chedge];
            if (!_bnd->IsOnDirBnd(*edgep))
            {
                if ( IsSubEdge(chedgeinparent) )
                {
                    const Uint paredge= ParentEdge(chedgeinparent);
                    const Uint num= NumOfSubEdge(chedgeinparent);
                    sol.SetDoF( *edgep, old_sol.val(*t.GetEdge(paredge), 0.25+0.5*num) );
                }
                else if ( IsSubInParFace(chedgeinparent) )
                {
                    Uint parface;
                    Uint pos;
                    WhichEdgeInFace(chedgeinparent, parface, pos);
                    sol.SetDoF( *edgep, old_sol.val(t, parface, edgebary[pos][0], edgebary[pos][1]) );
                }
                else
                {
                    sol.SetDoF( *edgep, old_sol.val(t, 0.25, 0.25, 0.25) );
                }
            }
        }
    }
}


template<class Data, class _BndData, class _VD>
inline void CopyValues( const TetraCL& t, P2EvalCL<Data, _BndData, _VD>& sol, const P2EvalCL<Data, _BndData, const _VD>& old_sol)
{
    typedef typename P2EvalCL<Data, _BndData, _VD>::BndDataCL BndCL;
    const BndCL* const _bnd= old_sol.GetBndData();

    for (TetraCL::const_VertexPIterator it= t.GetVertBegin(), end= t.GetVertEnd();
        it!=end; ++it)
        if (!_bnd->IsOnDirBnd( **it) )
            sol.SetDoF( **it, old_sol.val( **it) );
    for (TetraCL::const_EdgePIterator it= t.GetEdgesBegin(), end= t.GetEdgesEnd();
        it!=end; ++it)
        if (!_bnd->IsOnDirBnd( **it) )
            sol.SetDoF( **it, old_sol.val( **it, 0.5) );
}

template<class Data, class _BndData, class _VD>
void Adapt( P2EvalCL<Data, _BndData, _VD>& sol, const P2EvalCL<Data, _BndData, const _VD>& old_sol)
{
// Adapt a solution on a triangulation of a certain level, that has changed during the refinement.
// Notation: T = old triang, T' = new triang. Both T and T' are of the same level.
// This change can be classified in several cases (not complete...): Fot tetra t in T:
//    a) t is missing not only in T' but in all new triang levels
//       -> information is lost
//    b) t is missing in T', but is member of the new MultiGrid
//       -> t was refined, children of t are members of T', for these information is interpolated from t!
//    c) t and its brotherhood were replaced by other children
//       -> change of the parents' refinement rule, restore information in parent and interpolate!
// TODO: missing: handling of unrefined Tetras

// Adapt should be very robust in all occuring situations!!!

    const MultiGridCL& _MG= old_sol.GetMG();
    const Uint level= sol.GetLevel();
    const Uint old_idx= old_sol.GetSolution()->RowIdx->GetIdx();

    Assert( level==old_sol.GetLevel(),
        DROPSErrCL("Adapt: Same triang levels are expected\n"), ~0);

    // 1. Iterate tetras of triangulation: Interpolate missing unknowns
    // 2. Iterate tetras of triangulation: Copy known values
    // => known values override interpolated values

    for (MultiGridCL::const_TriangTetraIteratorCL sit= _MG.GetTriangTetraBegin(level), theend= _MG.GetTriangTetraEnd(level);
         sit!=theend; ++sit)
    {
         if ( old_sol.UnknownsMissing( *sit) )
         // Tetra war in T vorhanden -> neu in T'! Fall b)c)
             if ( sol.UnknownsMissing( *sit) )
             {
                 // Ergaenze evtl. fehlenden Idx des Vaters auf Edges! (nur fuer c) noetig)
                 for (TetraCL::const_EdgePIterator it= sit->GetParent()->GetEdgesBegin(), end= sit->GetParent()->GetEdgesEnd();
                     it!=end; ++it)
                     if ((*it)->IsRefined() && (*it)->GetMidVertex()->Unknowns.Exist(old_idx) )
                     {
                         // evtl. UnknownIdx fuer old_idx anlegen
                         (*it)->Unknowns.Prepare( old_idx);
                         // Indexwerte von MidVertex kopieren
                         //for (Uint i=0; i<NumUnknowns; ++i)
                             (*it)->Unknowns(old_idx)= (*it)->GetMidVertex()->Unknowns(old_idx);
                     }
                 // Interpoliere Werte vom Vater
                 InterpolateChildren( *sit->GetParent(), sol, old_sol);
             }
    }

    for (MultiGridCL::const_TriangTetraIteratorCL sit= _MG.GetTriangTetraBegin(level), theend= _MG.GetTriangTetraEnd(level);
         sit!=theend; ++sit)
    {
        if ( !old_sol.UnknownsMissing(*sit) )
            CopyValues( *sit, sol, old_sol);
    }
}




template<class Data, class _BndData, class _VD>
void Interpolate(P2EvalCL<Data, _BndData, _VD>& sol, const P2EvalCL<Data, _BndData, const _VD>& old_sol)
// This only works, if Interpolate is called after every refinement of the multigrid.
// Take care, that x and old_x are on successive triangulations.
{
#ifdef _PAR
    throw DROPSErrCL("Interpolate of P2-Elements not yet implemented for parallel use");
#endif
    typedef typename P2EvalCL<Data, _BndData, _VD>::BndDataCL BndCL;
    const BndCL* const _bnd= old_sol.GetBndData();
    const Uint old_idx= old_sol.GetSolution()->RowIdx->GetIdx();

    // All velocity-components use the same row-index and the same trianglevel
    const MultiGridCL& _MG= old_sol.GetMG();
    const Uint old_level= old_sol.GetLevel();
    const Uint level= sol.GetLevel();
    Uint num_vert_copy= 0, num_edge_copy= 0, num_child_edge= 0;

//    Assert( level==old_level || old_level==level-1, DROPSErrCL("Interpolate: successive triangs are expected\n"));

    // Iterate over all tetras in old_level, interpolate edge-DoF on diagonal, if there is one;
    // interpolate edge-DoF's of the children's non-copied edges
    for (MultiGridCL::const_TriangTetraIteratorCL sit= _MG.GetTriangTetraBegin(old_level), theend= _MG.GetTriangTetraEnd(old_level);
         sit!=theend; ++sit)
    {
        // If *sit is unrefined, all interpolation will be done via copying DoF on edges and vertices later
        if ( !sit->IsUnrefined() && (*sit->GetChildBegin())->IsInTriang(level) )
            InterpolateChildren( *sit, sol, old_sol);
    }
    // Iterate over all edges, interpolate values on new mid-vertices and edge-copies (plain copying)
    for (MultiGridCL::const_EdgeIterator sit= _MG.GetAllEdgeBegin(level), theend= _MG.GetAllEdgeEnd(level);
         sit!=theend; ++sit)
    {
        if ( sit->IsRefined() && sit->IsInTriang(old_level)
             && sit->GetMidVertex()->IsInTriang(level) && !_bnd->IsOnDirBnd(*sit->GetMidVertex()) ) // only new non-boundary vertices are interpolated
        {
            if (!sit->Unknowns.Exist(old_idx)) continue;
                sol.SetDoF( *sit->GetMidVertex(), old_sol.val(*sit, 0.5) );
                ++num_edge_copy;
        }
        else if ( sit->IsInTriang(old_level) && sit->IsInTriang(level) && !_bnd->IsOnDirBnd(*sit) )
            {
                if (!sit->Unknowns.Exist(old_idx)) continue;
                    sol.SetDoF( *sit, old_sol.val(*sit, 0.5) );
                    ++num_edge_copy;
            }
    }
    // Iterate over the vertices of the old triangulation and copy values
    for (MultiGridCL::const_TriangVertexIteratorCL sit= _MG.GetTriangVertexBegin(old_level), theend= _MG.GetTriangVertexEnd(old_level);
         sit!=theend; ++sit)
        if ( !_bnd->IsOnDirBnd(*sit) )
        {
            if (!sit->Unknowns.Exist(old_idx)) continue;
                sol.SetDoF(*sit, old_sol.val(*sit) );
                ++num_vert_copy;
        }

        std::cout << "Interpolate: " << num_vert_copy << " vertex-DoF copied, "
                                     << num_edge_copy << " edge-DoF copied, "
                                     << num_child_edge << " edge-DoF interpolated."
                                     << std::endl;
}

template <class VecDescT, class BndDataT, class Cont>
void RestrictP2(const TetraCL& s, const VecDescT& vd, const BndDataT& bnd, Cont& c)
{
    const Uint slvl= s.GetLevel();
    const Uint flvl= vd.GetLevel();
    Assert( slvl<=flvl, DROPSErrCL("RestrictP2: Tetra is on a finer level"
            "than the function."), ~0);

    typedef typename VecDescT::DataType VecT;
    typedef DoFHelperCL< typename BndDataT::bnd_type, VecT> DoFT;
    const VecT& v= vd.Data;
    const Uint idx= vd.RowIdx->GetIdx();
    for (Uint i= 0; i<NumVertsC; ++i)
        c[i]= !bnd.IsOnDirBnd( *s.GetVertex( i))
                ? DoFT::get( v, s.GetVertex( i)->Unknowns( idx))
                : bnd.GetDirBndValue( *s.GetVertex( i), vd.t);
    for (Uint i= 0; i<NumEdgesC; ++i) {
        const EdgeCL& e= *s.GetEdge( i);
        c[i+NumVertsC]= (slvl < flvl && e.IsRefined())
            ? ( !bnd.IsOnDirBnd( *e.GetMidVertex())
                ? DoFT::get( v, e.GetMidVertex()->Unknowns( idx))
                : bnd.GetDirBndValue( *e.GetMidVertex(), vd.t))
            : ( !bnd.IsOnDirBnd( e)
                ? DoFT::get( v, e.Unknowns( idx))
                : bnd.GetDirBndValue( e, vd.t));
    }
}


} // end of namespace DROPS
