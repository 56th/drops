/// \file surfacestokes_utils.h
/// \brief Error and initialization routines for surfacestokes.cpp
/// \author LNM RWTH Aachen: Thomas Jankuhn

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


#ifndef DROPS_SURFACESTOKES_UTILS_H
#define DROPS_SURFACESTOKES_UTILS_H

#include "surfactant/ifacetransp.h"
#include "surfactant/surfacestokes_funcs.h"

using namespace DROPS;

//////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////Error and Initialization/////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


template<class DiscFunType>
double L2_error (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
    const DiscFunType& discsol, DROPS::instat_scalar_fun_ptr extsol, double t= 0.)
{
    double d( 0.);
    const DROPS::Uint lvl = ls.GetLevel();
    DROPS::InterfaceTriangleCL triangle;
    DROPS::Quad5_2DCL<> qsol, qdiscsol, qpow;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
    	triangle.Init( *it, ls, lsbnd);
        if (triangle.Intersects()) { // We are at the phase boundary.
            for (int ch= 0; ch < 8; ++ch) {
            	triangle.ComputeForChild( ch);
                for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri) {
                    qsol.assign( *it, &triangle.GetBary( tri), extsol, t);
                    qdiscsol.assign(  *it, &triangle.GetBary( tri), discsol);
                    // d+= DROPS::Quad5_2DCL<>( std::pow( qdiscsol - qsol, 2)).quad( triangle.GetAbsDet( tri));
                    qpow = (qdiscsol-qsol)*(qdiscsol-qsol);
                    d+= qpow.quad( triangle.GetAbsDet( tri));
                }
            }
        }
    }
    return std::sqrt( d);
}

template<class DiscFunType>
double L2_Vector_error (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
    const DiscFunType& discsol, DROPS::instat_vector_fun_ptr extsol, double t= 0.)
{
    double d( 0.);
    const DROPS::Uint lvl = ls.GetLevel();
    DROPS::InterfaceTriangleCL triangle;
    DROPS::Quad5_2DCL<DROPS::Point3DCL> qsol, qdiscsol;
    DROPS::Quad5_2DCL<double> qpow;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        triangle.Init( *it, ls, lsbnd);
        if (triangle.Intersects()) { // We are at the phase boundary.
            for (int ch= 0; ch < 8; ++ch) {
                triangle.ComputeForChild( ch);
                for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri) {
                    qsol.assign( *it, &triangle.GetBary( tri), extsol, t);
                    qdiscsol.assign(  *it, &triangle.GetBary( tri), discsol);
                    // d+= DROPS::Quad5_2DCL<>( std::pow( qdiscsol - qsol, 2)).quad( triangle.GetAbsDet( tri));
                    qdiscsol-= qsol;
                    qpow = dot(qdiscsol, qdiscsol);
//                    qpow = dot(qdiscsol, qdiscsol) - dot(2*qdiscsol, qsol) + dot(qsol, qsol);
                    d+= qpow.quad( triangle.GetAbsDet( tri));
                }
            }
        }
    }
    return std::sqrt( d);
}

double L2_norm (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
    DROPS::instat_scalar_fun_ptr extsol, double t= 0.)
{
    double d( 0.);
    const DROPS::Uint lvl = ls.GetLevel();
    DROPS::InterfaceTriangleCL triangle;
    DROPS::Quad5_2DCL<> qsol;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
    	triangle.Init( *it, ls, lsbnd);
        if (triangle.Intersects()) { // We are at the phase boundary.
            for (int ch= 0; ch < 8; ++ch) {
            	triangle.ComputeForChild( ch);
                for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri) {
                    qsol.assign( *it, &triangle.GetBary( tri), extsol, t);
                    d+= DROPS::Quad5_2DCL<>( qsol*qsol).quad( triangle.GetAbsDet( tri));
                }
            }
        }
    }
    return std::sqrt( d);
}

void LinearLSInit (const DROPS::MultiGridCL& mg, DROPS::VecDescCL& ls, DROPS::instat_scalar_fun_ptr d, double t=0.)
{
    const DROPS::Uint lvl= ls.GetLevel(),
                      idx= ls.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= d( it->GetCoord(), t);

    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= ls.Data[it->Unknowns( idx)]=
            0.5*(ls.Data[it->GetVertex( 0)->Unknowns( idx)] + ls.Data[it->GetVertex( 1)->Unknowns( idx)]);
}

void InitScalar (const MultiGridCL& mg, VecDescCL& ic, instat_scalar_fun_ptr icf, double t= 0)
{
    const Uint lvl= ic.GetLevel(),
               idx= ic.RowIdx->GetIdx();
    ic.t= t;

    if (ic.RowIdx->NumUnknownsVertex())
        DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
            if (it->Unknowns.Exist( idx)) {
                ic.Data[it->Unknowns( idx)]= icf( it->GetCoord(), t);
            }
        }
    if (ic.RowIdx->NumUnknownsEdge())
        DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it) {
            if (it->Unknowns.Exist( idx))
                ic.Data[it->Unknowns( idx)]= icf( GetBaryCenter(*it) , t);
        }

}

// Does not work for Dirichlet Boundaries
void InterpolateP2Vec (const MultiGridCL& mg, const VecDescCL& vP1, VecDescCL& vP2)
{
    const DROPS::Uint lvl= vP2.GetLevel(),
                      idxP2= vP2.RowIdx->GetIdx(),
                      idxP1= vP1.RowIdx->GetIdx();


    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
        if (it->Unknowns.Exist( idxP2)) {
            for (int k=0; k<3; ++k)
                vP2.Data[it->Unknowns( idxP2)+k] = vP1.Data[it->Unknowns( idxP1)+k];
        }
    }

    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it) {
        if (it->Unknowns.Exist( idxP2)) {
            for (int k=0; k<3; ++k)
                vP2.Data[it->Unknowns( idxP2)+k] = 0.5*(vP1.Data[it->GetVertex( 0)->Unknowns( idxP1)+k] + vP1.Data[it->GetVertex( 1)->Unknowns( idxP1)+k]);
        }
    }
}

void InitVector(const MultiGridCL& MG, VecDescCL& vec, instat_vector_fun_ptr LsgVel, double t0 = 0)
{
    VectorCL& lsgvel= vec.Data;
    vec.t = t0;
    const Uint lvl  = vec.GetLevel(),
               vidx = vec.RowIdx->GetIdx();

    if( vec.RowIdx->NumUnknownsVertex()) {
        for (MultiGridCL::const_TriangVertexIteratorCL sit= MG.GetTriangVertexBegin( lvl),
            send= MG.GetTriangVertexEnd( lvl);
            sit != send; ++sit) {
            if ( sit->Unknowns.Exist( vidx))
                DoFHelperCL<Point3DCL, VectorCL>::set( lsgvel, sit->Unknowns( vidx),
                    LsgVel(sit->GetCoord(), t0));
        }
    }
    if( vec.RowIdx->NumUnknownsEdge()) {
        for (MultiGridCL::const_TriangEdgeIteratorCL sit= MG.GetTriangEdgeBegin( lvl),
             send= MG.GetTriangEdgeEnd( lvl);
            sit != send; ++sit) {
            if ( sit->Unknowns.Exist( vidx))
                DoFHelperCL<Point3DCL, VectorCL>::set( lsgvel, sit->Unknowns( vidx),
                    LsgVel( GetBaryCenter(*sit), t0));
        }
    }
}

// Scalar H1-error
double H1_error (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
                 const DROPS::VecDescCL& discsol, DROPS::instat_vector_fun_ptr extsol_grad, bool surfgradmark, double t=0.)
{
  double d( 0.), d_surf( 0.);
  const DROPS::Uint lvl = ls.GetLevel();

  LocalP1CL<Point3DCL> GradRefLP1[10], GradLP1[10];
  SMatrixCL<3,3> T;
  double det;

  IdxT num[10];

  Quad5_2DCL<Point3DCL> qsol_grad, qsol_surfgrad, qfullgrad, qsurfgrad;
  Quad5_2DCL<> qH1_full, qH1_surfgrad;

  P2DiscCL::GetGradientsOnRef( GradRefLP1);

  DROPS::InterfaceTriangleCL triangle;

  DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
    triangle.Init( *it, ls, lsbnd);
    if (triangle.Intersects()) { // We are at the phase boundary.
      GetLocalNumbP2NoBnd( num, *it, *(discsol.RowIdx));
      GetTrafoTr( T, det, *it);
      P2DiscCL::GetGradients( GradLP1, GradRefLP1, T);

      for (int ch= 0; ch < 8; ++ch) {
        triangle.ComputeForChild( ch);
        for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri) {
          Quad5_2DCL<Point3DCL> surfGrad[10];
          Quad5_2DCL<Point3DCL> Grad[10];
          
          for (int i=0; i < 10; ++i) {
            Grad[i].assign( GradLP1[i], &triangle.GetBary( tri));
            surfGrad[i].assign( GradLP1[i], &triangle.GetBary( tri));
            surfGrad[i].apply( triangle, &InterfaceTriangleCL::ApplyProj);
          }
          
          qsol_grad.assign( *it, &triangle.GetBary( tri), extsol_grad, t);
          qsol_surfgrad.assign( *it, &triangle.GetBary( tri), extsol_grad, t);
          qsol_surfgrad.apply( triangle, &InterfaceTriangleCL::ApplyProj);

          qfullgrad = qsol_grad;
          qsurfgrad = qsol_surfgrad;
          for (int i=0; i<10; i++) {
            qfullgrad -= discsol.Data[num[i]]*Grad[i];
            qsurfgrad -= discsol.Data[num[i]]*surfGrad[i];
          }

          qH1_full = dot(qfullgrad,qfullgrad);
          qH1_surfgrad = dot(qsurfgrad,qsurfgrad);
          d += qH1_full.quad( triangle.GetAbsDet( tri));
          d_surf += qH1_surfgrad.quad( triangle.GetAbsDet( tri));
        }
      }
    }
  }
  if(!surfgradmark)
    return std::sqrt( d);
  else
    return std::sqrt( d_surf);
}


void H1_Vector_error_P2 (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
                 const DROPS::VecDescCL& discsol, const DROPS::BndDataCL<Point3DCL>& bnddata, DROPS::instat_vector_fun_ptr extsol, DROPS::instat_vector_fun_ptr extsol_grad1, DROPS::instat_vector_fun_ptr extsol_grad2, DROPS::instat_vector_fun_ptr extsol_grad3, double eta, double& H1, double& simple_surfH1, double& advanced_surfH1, double& normal_velocity, double t=0.)
{
    double d( 0.), d_surf( 0.), d_advanced_surf( 0.), d_normal_velocity( 0.);
    const DROPS::Uint lvl = ls.GetLevel();

    LocalP2CL<> P2Hat[10];
    LocalP2CL<Point3DCL> loc_uh;
    LocalP1CL<Point3DCL> GradRefP2[10], GradP2[10];
    SMatrixCL<3,3> T;
    double det;

    LocalNumbP2CL n;

    Quad5_2DCL<Point3DCL> qsol_grad[3], qsol_surfgrad[3], qfullgrad[3], qsurfgrad[3], qerror, q_uh;
    Quad5_2DCL<> qH1_full, qH1_surfgrad, qH1_advanced_surfgrad, qL2_normal_velocity;

    P2DiscCL::GetP2Basis( P2Hat);
    P2DiscCL::GetGradientsOnRef( GradRefP2);

    DROPS::InterfaceTriangleCL triangle;
    Quad5_2DCL<> Hat[10];
    Quad5_2DCL<Point3DCL> surfGrad[10];
    Quad5_2DCL<Point3DCL> Grad[10];
    Quad5_2DCL<Point3DCL> normal;
    Quad5_2DCL<> qsurfgrad_comp[3][3], qerror_comp[3], qsol_comp[3];
    GridFunctionCL<> normal_comp[3];

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        triangle.Init( *it, ls, lsbnd);
        if (triangle.Intersects()) { // We are at the phase boundary.
            n.assign( *it, *(discsol.RowIdx), discsol.RowIdx->GetBndInfo());
            loc_uh.assign( *it, discsol, bnddata);
            GetTrafoTr( T, det, *it);
            P2DiscCL::GetGradients( GradP2, GradRefP2, T);

            for (int ch= 0; ch < 8; ++ch) {
                triangle.ComputeForChild( ch);
                for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri) {
                    normal = triangle.GetImprovedNormal(tri);
                    for (int i=0; i<3; ++i) {
                        normal_comp[i].resize(normal.size());
                        ExtractComponent(normal, normal_comp[i], i);
                    }

                    for (int i=0; i < 10; ++i) {
                        Hat[i].assign( P2Hat[i], &triangle.GetBary( tri));
                        Grad[i].assign( GradP2[i], &triangle.GetBary( tri));
                        surfGrad[i]= Grad[i];
                        surfGrad[i]-= dot( surfGrad[i], normal)*normal;
                    }
                    q_uh.assign( loc_uh, &triangle.GetBary( tri));
                    qerror.assign( *it, &triangle.GetBary( tri), extsol, t);
                    qsol_grad[0].assign( *it, &triangle.GetBary( tri), extsol_grad1, t);
                    qsol_grad[1].assign( *it, &triangle.GetBary( tri), extsol_grad2, t);
                    qsol_grad[2].assign( *it, &triangle.GetBary( tri), extsol_grad3, t);
                    for(int i=0; i<3; ++i) {
                        qsol_surfgrad[i] = qsol_grad[i];
                        qsol_surfgrad[i]-= dot( qsol_surfgrad[i], normal)*normal;
                    }
                    //qsol_surfgrad.apply( triangle, &InterfaceTriangleCL::ApplyProj);

                    for (int i=0; i < 3; ++i) {
                        qfullgrad[i] = qsol_grad[i];
                        qsurfgrad[i] = qsol_surfgrad[i];
                        ExtractComponent( qerror, qerror_comp[i], i);
                    }
                    for (int i=0; i<10; i++) {
                        for (int j=0; j<3; ++j) {
                            if(n.WithUnknowns(i) || !(n.num[i]==NoIdx)) {
                                qsol_comp[j] += discsol.Data[n.num[i]+j]*Hat[i];
                                qerror_comp[j] -= discsol.Data[n.num[i]+j]*Hat[i];
                                qfullgrad[j] -= discsol.Data[n.num[i]+j]*Grad[i];
                                qsurfgrad[j] -= discsol.Data[n.num[i]+j]*surfGrad[i];
                            }
                        }
                    }

                    //if(surfgradmark) {
                        for (int i=0; i<3; ++i) {
                            for (int j=0; j<3; ++j) {
                                ExtractComponent( qsurfgrad[j], qsurfgrad_comp[j][i], i);
                            }
                        }
                        for (int i=0; i<3; ++i) {
                            for(int j=0; j<3; ++j) {
                                for(int k=0; k<3; ++k) {
                                    qsurfgrad_comp[j][i] -= qsurfgrad_comp[k][i]*normal_comp[j]*normal_comp[k];
                                }
                            }
                        }
                    //}
                    qH1_full = 0;
                    qH1_surfgrad = 0;
                    qH1_advanced_surfgrad = 0;
                    qL2_normal_velocity = dot( normal, q_uh);
                    qL2_normal_velocity*= qL2_normal_velocity;
                    for (int k=0; k<3; ++k) {
                        qH1_full += dot( qfullgrad[k], qfullgrad[k]);
                        qH1_advanced_surfgrad += qerror_comp[k]*qerror_comp[k] + eta*qerror_comp[k]*normal_comp[k]*qerror_comp[k]*normal_comp[k];
                        for(int i=0; i<3; ++i) {
//                            qL2_normal_velocity += qsol_comp[k]*normal_comp[k]*qsol_comp[i]*normal_comp[i];
                            qH1_surfgrad += qsurfgrad_comp[k][i]*qsurfgrad_comp[k][i];
                            qH1_advanced_surfgrad += qsurfgrad_comp[k][i]*qsurfgrad_comp[k][i] + qsurfgrad_comp[k][i]*qsurfgrad_comp[i][k] + qsurfgrad_comp[k][i]*qsurfgrad_comp[i][k] + qsurfgrad_comp[i][k]*qsurfgrad_comp[i][k];
                        }
                    }
                    d += qH1_full.quad( triangle.GetAbsDet( tri));
                    d_surf += qH1_surfgrad.quad( triangle.GetAbsDet( tri));
                    d_advanced_surf += qH1_advanced_surfgrad.quad( triangle.GetAbsDet( tri));
                    d_normal_velocity += qL2_normal_velocity.quad(triangle.GetAbsDet( tri));
                }
            }
        }
    }
    H1 = std::sqrt( d);
    simple_surfH1 = std::sqrt( d_surf);
    advanced_surfH1 = std::sqrt( d_advanced_surf);
    normal_velocity = std::sqrt( d_normal_velocity);
}


void H1_Vector_error_P1 (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
                 const DROPS::VecDescCL& discsol, const DROPS::BndDataCL<Point3DCL>& bnddata, DROPS::instat_vector_fun_ptr extsol, DROPS::instat_vector_fun_ptr extsol_grad1, DROPS::instat_vector_fun_ptr extsol_grad2, DROPS::instat_vector_fun_ptr extsol_grad3, double eta, double& H1, double& simple_surfH1, double& advanced_surfH1, double& normal_velocity, double t=0.)
{
    double d( 0.), d_surf( 0.), d_advanced_surf( 0.), d_normal_velocity( 0.);
    const DROPS::Uint lvl = ls.GetLevel();

    LocalP1CL<> P1Hat[4];
    LocalP1CL<Point3DCL> loc_uh;
    Point3DCL GradP1[4];
    SMatrixCL<3,3> T;
    double det;

    //IdxT num[4];
    LocalNumbP1CL n;

    Quad5_2DCL<Point3DCL> qsol_grad[3], qsol_surfgrad[3], qfullgrad[3], qsurfgrad[3], qerror, q_uh;
    Quad5_2DCL<> qH1_full, qH1_surfgrad, qH1_advanced_surfgrad, qL2_normal_velocity;

    DROPS::InterfaceTriangleCL triangle;
    Quad5_2DCL<> Hat[4];
    Quad5_2DCL<Point3DCL> surfGrad[4];
    Quad5_2DCL<Point3DCL> Grad[4];
    Quad5_2DCL<Point3DCL> normal;
    Quad5_2DCL<> qsurfgrad_comp[3][3], qerror_comp[3], qsol_comp[3];
    GridFunctionCL<> normal_comp[3];

    P1DiscCL::GetP1Basis( P1Hat);

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        triangle.Init( *it, ls, lsbnd);
        if (triangle.Intersects()) { // We are at the phase boundary.
            //GetLocalNumbP1NoBnd( num, *it, *(discsol.RowIdx));
            n.assign( *it, *(discsol.RowIdx), discsol.RowIdx->GetBndInfo());
            loc_uh.assign( *it, discsol, bnddata);
            GetTrafoTr( T, det, *it);
            P1DiscCL::GetGradients( GradP1, T);

            for (int ch= 0; ch < 8; ++ch) {
                triangle.ComputeForChild( ch);
                for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri) {
                    normal = triangle.GetImprovedNormal(tri);
                    for (int i=0; i<3; ++i) {
                        normal_comp[i].resize(normal.size());
                        ExtractComponent(normal, normal_comp[i], i);
                    }

                    for (int i=0; i < 4; ++i) {
                        Hat[i].assign( P1Hat[i], &triangle.GetBary( tri));
                        Grad[i] = GradP1[i];
                        surfGrad[i]= Grad[i];
                        surfGrad[i]-= dot( surfGrad[i], normal)*normal;
                        //surfGrad[i].apply( triangle, &InterfaceTriangleCL::ApplyProj);
                    }
                    q_uh.assign( loc_uh, &triangle.GetBary( tri));
                    qerror.assign( *it, &triangle.GetBary( tri), extsol, t);
                    qsol_grad[0].assign( *it, &triangle.GetBary( tri), extsol_grad1, t);
                    qsol_grad[1].assign( *it, &triangle.GetBary( tri), extsol_grad2, t);
                    qsol_grad[2].assign( *it, &triangle.GetBary( tri), extsol_grad3, t);
                    for(int i=0; i<3; ++i) {
                        qsol_surfgrad[i] = qsol_grad[i];
                        qsol_surfgrad[i]-= dot( qsol_surfgrad[i], normal)*normal;
                    }
                    //qsol_surfgrad.apply( triangle, &InterfaceTriangleCL::ApplyProj);

                    for (int i=0; i < 3; ++i) {
                        qfullgrad[i] = qsol_grad[i];
                        qsurfgrad[i] = qsol_surfgrad[i];
                        ExtractComponent( qerror, qerror_comp[i], i);
                    }
                    for (int i=0; i<4; i++) {
                        for (int j=0; j<3; ++j) {
                            if(n.WithUnknowns(i) || !(n.num[i]==NoIdx)) {
                                qsol_comp[j] += discsol.Data[n.num[i]+j]*Hat[i];
                                qerror_comp[j] -= discsol.Data[n.num[i]+j]*Hat[i];
                                qfullgrad[j] -= discsol.Data[n.num[i]+j]*Grad[i];
                                qsurfgrad[j] -= discsol.Data[n.num[i]+j]*surfGrad[i];
                            }
                        }
                    }
                    //if(surfgradmark) {
                        for (int i=0; i<3; ++i) {
                            for (int j=0; j<3; ++j) {
                                ExtractComponent( qsurfgrad[j], qsurfgrad_comp[j][i], i);
                            }
                        }
                        for (int i=0; i<3; ++i) {
                            for(int j=0; j<3; ++j) {
                                for(int k=0; k<3; ++k) {
                                    qsurfgrad_comp[j][i] -= qsurfgrad_comp[k][i]*normal_comp[j]*normal_comp[k];
                                }
                            }
                        }
                    //}
                    qH1_full = 0;
                    qH1_surfgrad = 0;
                    qH1_advanced_surfgrad = 0;
                    qL2_normal_velocity = dot( normal, q_uh);
                    qL2_normal_velocity*= qL2_normal_velocity;
                    for (int k=0; k<3; ++k) {
                        qH1_full += dot( qfullgrad[k], qfullgrad[k]);
                        qH1_advanced_surfgrad += qerror_comp[k]*qerror_comp[k] + eta*qerror_comp[k]*normal_comp[k]*qerror_comp[k]*normal_comp[k];
                        for(int i=0; i<3; ++i) {
//                            qL2_normal_velocity += qsol_comp[k]*normal_comp[k]*qsol_comp[i]*normal_comp[i];
                            qH1_surfgrad += qsurfgrad_comp[k][i]*qsurfgrad_comp[k][i];
                            qH1_advanced_surfgrad += qsurfgrad_comp[k][i]*qsurfgrad_comp[k][i] + qsurfgrad_comp[k][i]*qsurfgrad_comp[i][k] + qsurfgrad_comp[k][i]*qsurfgrad_comp[i][k] + qsurfgrad_comp[i][k]*qsurfgrad_comp[i][k];
                        }
                    }
                    d += qH1_full.quad( triangle.GetAbsDet( tri));
                    d_surf += qH1_surfgrad.quad( triangle.GetAbsDet( tri));
                    d_advanced_surf += qH1_advanced_surfgrad.quad( triangle.GetAbsDet( tri));
                    d_normal_velocity += qL2_normal_velocity.quad( triangle.GetAbsDet( tri));
                    if(std::isnan( d_normal_velocity)) {
                        std::cout << "NaN" << std::endl;
                    }
                }
            }
        }
    }
    H1 = std::sqrt( d);
    simple_surfH1 = std::sqrt( d_surf);
    advanced_surfH1 = std::sqrt( d_advanced_surf);
    normal_velocity = std::sqrt( d_normal_velocity);
}

// parse iterations from stream
void ComputeAverageIterations( std::stringstream& stream, double& average)
{
    int max = 0;
    double counter = 0.;
    while( !stream.eof()) {
        std::string row;
        std::getline(stream, row);
        if( row.empty()) break;
        int posofr = row.find('r') - 5;
        std::string stringvalue = row.substr(5, posofr);
        int value = std::stoi(stringvalue);
        if( ((value == 1) && (max == 1)) || ((value == 1) && (max > 1))) {
            average += max;
            //std::cout << max << std::endl;
            max = 1;
            counter++;
        }
        if( value > max) {
            max = value;
        }
    }
    if( !(counter == 0.))
        average = average/counter;
}

void ComputeVariationFromAverageIterations( std::stringstream& stream, double average, double& variation)
{
    int max = 0;
    double counter = 0.;
    while( !stream.eof()) {
        std::string row;
        std::getline(stream, row);
        if( row.empty()) break;
        int posofr = row.find('r') - 5;
        std::string stringvalue = row.substr(5, posofr);
        int value = std::stoi(stringvalue);
        if( ((value == 1) && (max == 1)) || ((value == 1) && (max > 1))) {
            if(std::abs(average - max) > variation) {
                variation = std::abs(average - max);
            }
            max = 1;
            counter++;
        }
        if( value > max) {
            max = value;
        }
    }
}

#endif

