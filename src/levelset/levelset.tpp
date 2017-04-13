/// \file levelset.tpp
/// \brief levelset equation for two phase flow problems
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#include <fstream>
#include "num/accumulator.h"
#include "misc/progressaccu.h"
#include "misc/scopetimer.h"

namespace DROPS
{

template<class DiscVelSolT>
void LevelsetP2CL::GetInfo( double& maxGradPhi, double& Volume, Point3DCL& bary, Point3DCL& vel, const DiscVelSolT& velsol, Point3DCL& minCoord, Point3DCL& maxCoord, double& surfArea) const
/**
 * - \p maxGradPhi is the maximal 2-norm of the gradient of the level set function. This can be used as an indicator to decide
 *   whether a reparametrization should be applied.
 * - \p Volume is the volume inside the approximate interface consisting of planar segments.
 * - \p bary is the barycenter of the droplet.
 * - \p vel is the velocity of the barycenter of the droplet.
 * - The entries of \p minCoord store the minimal x, y and z coordinates of the approximative interface, respectively.
 * - The entries of \p maxCoord store the maximal x, y and z coordinates of the approximative interface, respectively.
 * - \p surfArea is the surface area of the approximative interface
 */
{
    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    SMatrixCL<3,3> T;
    double det, absdet;
    InterfaceTetraCL tetra;
    InterfaceTriangleCL triangle;

    P2DiscCL::GetGradientsOnRef( GradRef);
    maxGradPhi= -1.;
    Volume= 0.;
    surfArea= 0.;
    bary[0]= bary[1]= bary[2]= 0;
    vel[0]= vel[1]= vel[2]= 0;
    minCoord[0]= minCoord[1]= minCoord[2]= 1e99;
    maxCoord[0]= maxCoord[1]= maxCoord[2]= -1e99;
    LocalP2CL<double> ones( 1.);
    LocalP2CL<Point3DCL> Coord, Vel;

    for (MultiGridCL::const_TriangTetraIteratorCL it=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(idx.TriangLevel()), end=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(idx.TriangLevel());
        it!=end; ++it)
    {
        GetTrafoTr( T, det, *it);
        absdet= std::abs( det);
        P2DiscCL::GetGradients( Grad, GradRef, T); // Gradienten auf aktuellem Tetraeder

        tetra.Init( *it, Phi, BndData_);
        triangle.Init( *it, Phi, BndData_);

        // compute maximal norm of grad Phi
        Quad2CL<Point3DCL> gradPhi;
        for (int v=0; v<10; ++v) // init gradPhi, Coord
        {
            gradPhi+= tetra.GetPhi(v)*Grad[v];
            Coord[v]= v<4 ? it->GetVertex(v)->GetCoord() : GetBaryCenter( *it->GetEdge(v-4));
        }
        Vel.assign( *it, velsol);
        VectorCL normGrad( 5);
        for (int v=0; v<5; ++v) // init normGrad
            normGrad[v]= norm( gradPhi[v]);
        const double maxNorm= normGrad.max();
        if (maxNorm > maxGradPhi) maxGradPhi= maxNorm;

        for (int ch=0; ch<8; ++ch)
        {
            // compute volume, barycenter and velocity
            tetra.ComputeCutForChild(ch);
            Volume+= tetra.quad( ones, absdet, false);
            bary+= tetra.quad( Coord, absdet, false);
            vel+= tetra.quad( Vel, absdet, false);

            // find minimal/maximal coordinates of interface
            if (!triangle.ComputeForChild(ch)) // no patch for this child
                continue;
            for (int tri=0; tri<triangle.GetNumTriangles(); ++tri)
                surfArea+= triangle.GetAbsDet(tri);
            for (Uint i=0; i<triangle.GetNumPoints(); ++i)
            {
                const Point3DCL p= triangle.GetPoint(i);
                for (int j=0; j<3; ++j)
                {
                    if (p[j] < minCoord[j]) minCoord[j]= p[j];
                    if (p[j] > maxCoord[j]) maxCoord[j]= p[j];
                }
            }
        }
    }
#ifdef _PAR
    // Globalization of  data
    // -----
    const Point3DCL local_bary(bary), local_vel(vel), local_minCoord(minCoord), local_maxCoord(maxCoord);
    maxGradPhi= ProcCL::GlobalMax(maxGradPhi);
    Volume    = ProcCL::GlobalSum(Volume);
    surfArea  = ProcCL::GlobalSum(surfArea);
    ProcCL::GlobalSum(Addr(local_bary), Addr(bary), 3);
    ProcCL::GlobalSum(Addr(local_vel), Addr(vel), 3);
    ProcCL::GlobalMin(Addr(local_minCoord), Addr(minCoord), 3);
    ProcCL::GlobalMax(Addr(local_maxCoord), Addr(maxCoord), 3);
#endif

    bary/= Volume;
    vel/= Volume;
    surfArea*= 0.5;
}

template<class DiscVelSolT>
void LevelsetP2CL::GetFilmInfo( double& maxGradPhi, double& Volume, Point3DCL& vel, const DiscVelSolT& velsol,  double& x, double& z, double& h) const
/*
 */
{
    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    SMatrixCL<3,3> T;
    double det, absdet;
    InterfaceTetraCL tetra;
    InterfaceTriangleCL triangle;
    double lamda[3];
    bool   triangle_found = false;

    P2DiscCL::GetGradientsOnRef( GradRef);
    maxGradPhi= -1.;
    Volume= 0.;
    vel[0]= vel[1]= vel[2]= 0;
    LocalP2CL<double> ones( 1.);
    LocalP2CL<Point3DCL> Coord, Vel;

    for (MultiGridCL::const_TriangTetraIteratorCL it=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(idx.TriangLevel()), end=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(idx.TriangLevel());
        it!=end; ++it)
    {
        GetTrafoTr( T, det, *it);
        absdet= std::abs( det);
        P2DiscCL::GetGradients( Grad, GradRef, T); // Gradienten auf aktuellem Tetraeder

        tetra.Init( *it, Phi, BndData_);
        triangle.Init( *it, Phi, BndData_);

        // compute maximal norm of grad Phi
        Quad2CL<Point3DCL> gradPhi;
        for (int v=0; v<10; ++v) // init gradPhi, Coord
        {
            gradPhi+= tetra.GetPhi(v)*Grad[v];
            Coord[v]= v<4 ? it->GetVertex(v)->GetCoord() : GetBaryCenter( *it->GetEdge(v-4));
        }
        Vel.assign( *it, velsol);
        VectorCL normGrad( 5);
        for (int v=0; v<5; ++v) // init normGrad
            normGrad[v]= norm( gradPhi[v]);
        const double maxNorm= normGrad.max();
        if (maxNorm > maxGradPhi) maxGradPhi= maxNorm;

        for (int ch=0; ch<8; ++ch)
        {
            // compute volume, barycenter and velocity
            tetra.ComputeCutForChild(ch);
            Volume+= tetra.quad( ones, absdet, false);
            vel+= tetra.quad( Vel, absdet, false);

            // find minimal/maximal coordinates of interface
            if (!triangle.ComputeForChild(ch)) // no patch for this child
                continue;
            if(!triangle_found)
            {
                Point3DCL p1= triangle.GetPoint(0);
                Point3DCL p2= triangle.GetPoint(1);
                Point3DCL p3= triangle.GetPoint(2);
                double det = (p1[0] - p3[0]) * (p2[2] - p3[2]) - (p2[0] - p3[0]) * (p1[2] - p3[2]);
                 //compute the barycentric coordinates for x, z in the projection of the triangle p1p2p3 to xz plane
                lamda[0] = (p2[2] - p3[2]) * (x - p3[0]) + (p3[0] - p2[0]) * (z - p3[2]);
                lamda[0]/= det;
                lamda[1] = (p3[2] - p1[2]) * (x - p3[0]) + (p1[0] - p3[0]) * (z - p3[2]);
                lamda[1]/= det;
                lamda[2]= 1. - lamda[0] - lamda[1];
                //To see if point ( x, z ) in the projected triangle
                if( !(lamda[0]< 0) && !(lamda[1]< 0) &&!(lamda[2]< 0))
                    triangle_found = true;
                if( (!triangle_found)&&triangle.GetNumPoints() == 4)
                {
                    p1= triangle.GetPoint(1);
                    p2= triangle.GetPoint(2);
                    p3= triangle.GetPoint(3);
                    double det = (p1[0] - p3[0]) * (p2[2] - p3[2]) - (p2[0] - p3[0]) * (p1[2] - p3[2]);
                    lamda[0] = (p2[2] - p3[2]) * (x - p3[0]) + (p3[0] - p2[0]) * (z - p3[2]);
                    lamda[0]/= det;
                    lamda[1] = (p3[2] - p1[2]) * (x - p3[0]) + (p1[0] - p3[0]) * (z - p3[2]);
                    lamda[1]/= det;
                    lamda[2]= 1. - lamda[0] - lamda[1];
                    if( !(lamda[0]< 0) && !(lamda[1]< 0) &&!(lamda[2]< 0))
                        triangle_found = true;
                }
                //if the triangle founded, compute the y;
                if(triangle_found)
                {
                   //std::cout<<"Here"<<std::endl;
                   Point3DCL n(0.);
                   Point3DCL a = p2 -p1;
                   Point3DCL b = p3 -p1;
                   n[0] = a[1] * b[2] - a[2] * b[1];
                   n[1] = a[2] * b[0] - a[0] * b[2];
                   n[2] = a[0] * b[1] - a[1] * b[0];
                   h = p1[0] * n[0] + p1[1] * n[1] + p1[2] * n[2] - n[0] * x - n[2] * z;
                   h /= n[1];
                }
            }
        }
    }
#ifdef _PAR
    // Globalization of  data
    // -----
    const Point3DCL local_vel(vel);
    maxGradPhi= ProcCL::GlobalMax(maxGradPhi);
    Volume    = ProcCL::GlobalSum(Volume);
    ProcCL::GlobalSum(Addr(local_vel), Addr(vel), 3);
#endif
    vel/= Volume;
}

/// \brief Accumulator to set up the matrices E and H for the level set equation.
template<class DiscVelSolT>
class LevelsetAccumulator_P2CL : public TetraAccumulatorCL
{
    LevelsetP2CL& ls_;
    const DiscVelSolT& vel_;
    const double SD_;
    SparseMatBuilderCL<double> *bE_, *bH_;

    Quad5CL<Point3DCL> Grad[10], GradRef[10], u_loc;
    Quad5CL<double> u_Grad[10]; // fuer u grad v_i
    SMatrixCL<3,3> T;
    LocalNumbP2CL n;

  public:
    LevelsetAccumulator_P2CL( LevelsetP2CL& ls, const DiscVelSolT& vel, double SD, __UNUSED__ double dt)
      : ls_(ls), vel_(vel), SD_(SD)
    { P2DiscCL::GetGradientsOnRef( GradRef); }

    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();
    ///\brief Do setup of E, H on given tetra
    void visit (const TetraCL&);

    TetraAccumulatorCL* clone (int /*tid*/) { return new LevelsetAccumulator_P2CL ( *this); };
};

template<class DiscVelSolT>
void LevelsetAccumulator_P2CL<DiscVelSolT>::begin_accumulation ()
{
    const IdxT num_unks= ls_.Phi.RowIdx->NumUnknowns();
    bE_= new SparseMatBuilderCL<double>(&ls_.E, num_unks, num_unks);
    bH_= new SparseMatBuilderCL<double>(&ls_.H, num_unks, num_unks);

#ifndef _PAR
    __UNUSED__ const IdxT allnum_unks= num_unks;
#else
    __UNUSED__ const IdxT allnum_unks= DROPS::ProcCL::GlobalSum(num_unks);
#endif
    Comment("entering LevelsetAccumulator_P2CL: " << allnum_unks << " levelset unknowns.\n", DebugDiscretizeC);
}

template<class DiscVelSolT>
void LevelsetAccumulator_P2CL<DiscVelSolT>::finalize_accumulation ()
{
    bE_->Build();
    delete bE_;
    bH_->Build();
    delete bH_;
}

template<class DiscVelSolT>
void LevelsetAccumulator_P2CL<DiscVelSolT>::visit (const TetraCL& t)
/**Sets up the stiffness matrices: <br>
   E is of mass matrix type:    E_ij = ( v_j       , v_i + SD * u grad v_i ) <br>
   H describes the convection:  H_ij = ( u grad v_j, v_i + SD * u grad v_i ) <br>
   where v_i, v_j denote the ansatz functions.
   \todo: implementation of other boundary conditions
*/
{
    double det;
    GetTrafoTr( T, det, t);
    P2DiscCL::GetGradients( Grad, GradRef, T);
    const double absdet= std::fabs( det),
            h_T= std::pow( absdet, 1./3.);

    // save information about the edges and verts of the tetra in Numb
    n.assign( t, *ls_.Phi.RowIdx, ls_.GetBndData());

    // save velocities inside tetra for quadrature in u_loc
    u_loc.assign( t, vel_);

    for(int i=0; i<10; ++i)
        u_Grad[i]= dot( u_loc, Grad[i]);

    double maxV = 0.; // scaling of SD parameter (cf. master thesis of Rodolphe Prignitz)
    //const double limit= (h_T*h_T)/dt;
    const double limit= 1e-3;
    for(int i=0; i<Quad5DataCL::NumNodesC; ++i)
        maxV = std::max( maxV, u_loc[i].norm());
    maxV= std::max( maxV, limit/h_T); // no scaling for extremely small velocities
    //double maxV= 1; // no scaling

    SparseMatBuilderCL<double> &bE= *bE_, &bH= *bH_;
    for(int i=0; i<10; ++i)    // assemble row Numb[i]
        for(int j=0; j<10; ++j)
        {
            // E is of mass matrix type:    E_ij = ( v_j       , v_i + SD * u grad v_i )
            bE( n.num[i], n.num[j])+= P2DiscCL::GetMass(i,j) * absdet
                                 + u_Grad[i].quadP2(j, absdet)*SD_/maxV*h_T;

            // H describes the convection:  H_ij = ( u grad v_j, v_i + SD * u grad v_i )
            bH( n.num[i], n.num[j])+= u_Grad[j].quadP2(i, absdet)
                                 + Quad5CL<>(u_Grad[i]*u_Grad[j]).quad( absdet) * SD_/maxV*h_T;
        }
}


/// \brief Accumulator to set up the matrices E and H for the discontinous P2 level set equation.
/// takes care of volume integarls only
/// creates SparseMatBuilder for E and H, but creates E only
///creating H is accomplished in the FaceAccumulator
template<class DiscVelSolT>
class LevelsetTetraAccumulator_P2DCL : public TetraAccumulatorCL
{
    LevelsetP2DiscontCL& ls_;
    const DiscVelSolT& vel_;
    SparseMatBuilderCL<double> *bE_, *bH_;

    Quad5CL<Point3DCL> Grad[10], GradRef[10], u_loc;
    Quad5CL<double> u_Grad[10]; // fuer u grad v_i
    SMatrixCL<3,3> T;
    LocalNumbP2CL n;

  public:
    LevelsetTetraAccumulator_P2DCL( LevelsetP2CL& ls, const DiscVelSolT& vel)
      : ls_(dynamic_cast<LevelsetP2DiscontCL&>(ls)), vel_(vel)
    { P2DiscCL::GetGradientsOnRef( GradRef); }

    SparseMatBuilderCL<double> * GetHMatrix() {return bH_;};
    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();
    ///\brief Do setup of E, H on given tetra
    void visit (const TetraCL&);

    TetraAccumulatorCL* clone (int /*tid*/) { return new LevelsetTetraAccumulator_P2DCL ( *this); };
};


template<class DiscVelSolT>
void LevelsetTetraAccumulator_P2DCL<DiscVelSolT>::begin_accumulation ()
{
    const IdxT num_unks= ls_.Phi.RowIdx->NumUnknowns();
    bE_= new SparseMatBuilderCL<double>(&ls_.E, num_unks, num_unks);
    bH_= new SparseMatBuilderCL<double>(&ls_.H, num_unks, num_unks);

#ifndef _PAR
    __UNUSED__ const IdxT allnum_unks= num_unks;
#else
    __UNUSED__ const IdxT allnum_unks= DROPS::ProcCL::GlobalSum(num_unks);
#endif
    Comment("entering LevelsetTetraAccumulator_P2DCL: " << allnum_unks << " levelset unknowns.\n", DebugDiscretizeC);
}

template<class DiscVelSolT>
void LevelsetTetraAccumulator_P2DCL<DiscVelSolT>::finalize_accumulation ()
{
    bE_->Build();
    delete bE_;
    // bH_->Build(); // FaceAccumulator takes care of bH
    // delete bH_;
#ifndef _PAR
    Comment(ls_.E.num_nonzeros() << " nonzeros in E, "<< ls_.H.num_nonzeros() << " nonzeros in H! " << std::endl, DebugDiscretizeC);
#endif
}

template<class DiscVelSolT>
void LevelsetTetraAccumulator_P2DCL<DiscVelSolT>::visit (const TetraCL& t)
/**Sets up the stiffness matrices: <br>
   E is of mass matrix type:    E_ij =   (  v_j ,        v_i ) <br>
   H describes the convection:  H_ij = - (  v_j , u grad v_i)  <br>
   where v_i, v_j denote the ansatz functions.
   \todo: implementation of other boundary conditions
*/
{
    double det;
    GetTrafoTr( T, det, t);
    P2DiscCL::GetGradients( Grad, GradRef, T);
    const double absdet= std::fabs( det);

    // save information about the edges and verts of the tetra in Numb
    n.assign( t, *ls_.Phi.RowIdx, ls_.GetBndData()); //BndData is dummy
    // save velocities inside tetra for quadrature in u_loc
    u_loc.assign( t, vel_);

    for(int i=0; i<10; ++i)
        u_Grad[i]= dot( u_loc, Grad[i]);

    SparseMatBuilderCL<double> &bE= *bE_, &bH= *bH_;
    for(int i=0; i<10; ++i)    // assemble row Numb[i]
        for(int j=0; j<10; ++j)
        {
            // E is of mass matrix type:    E_ij = ( v_j       , v_i )
            bE( n.num[i], n.num[j])+= P2DiscCL::GetMass(i,j) * absdet;
            // H describes the convection:  H_ij = - (  v_j , u grad v_i)
            bH( n.num[i], n.num[j])-= u_Grad[i].quadP2(j, absdet);
        }
}



/// \brief Accumulator to set up the matrix H for the level set equation.
/// matrix builder is assumend to be created before
/// matric H is build in finalize_accumulation()
/// face integrals are added to volume integrals
template<class DiscVelSolT>
class LevelsetFaceAccumulator_P2DCL : public FaceAccumulatorCL
{
    LevelsetP2DiscontCL& ls_;
    const DiscVelSolT& vel_;

    SparseMatBuilderCL<double> *bH_;
    VecDescCL & rhs_;

    Uint lvl;
    LocalNumbP2CL n[2];

  public:
    LevelsetFaceAccumulator_P2DCL( LevelsetP2CL& ls, const DiscVelSolT& vel, SparseMatBuilderCL<double> * bH, VecDescCL & rhs)
        : ls_(dynamic_cast<LevelsetP2DiscontCL&>(ls)), vel_(vel), bH_(bH), rhs_(rhs), lvl(ls.idx.TriangLevel())
        { ; }

    ///\brief Initializes load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrix H
    void finalize_accumulation();
    ///\brief Do setup of H on given face
    void visit (const FaceCL&);

    FaceAccumulatorCL* clone (int ) { return new LevelsetFaceAccumulator_P2DCL ( *this); };
};

template<class DiscVelSolT>
void LevelsetFaceAccumulator_P2DCL<DiscVelSolT>::begin_accumulation ()
{
    // TetraAccumulator takes care of MatrixBuilder construction
    Comment("entering LevelsetFaceAccumulator_P2DCL: " << ls_.Phi.RowIdx->NumUnknowns() << " levelset unknowns.\n", DebugDiscretizeC);
    rhs_.SetIdx(ls_.Phi.RowIdx);
    rhs_.Data = 0. ;
}

template<class DiscVelSolT>
void LevelsetFaceAccumulator_P2DCL<DiscVelSolT>::finalize_accumulation ()
{
    bH_->Build();
    delete bH_;
    std::cout << ls_.E.num_nonzeros() << " nonzeros in E, "<< ls_.H.num_nonzeros() << " nonzeros in H! " << std::endl;
#ifndef _PAR
    Comment(ls_.E.num_nonzeros() << " nonzeros in E, "<< ls_.H.num_nonzeros() << " nonzeros in H! " << std::endl, DebugDiscretizeC);
#endif
}

template<class DiscVelSolT>
void LevelsetFaceAccumulator_P2DCL<DiscVelSolT>::visit (const FaceCL& face)
/**Sets up the stiffness matrices: <br>
   H describes the convection:  H_ij(face) = (un v_i,h(v_j)), h upwind flux operator <br>
   where v_i, v_j denote the ansatz functions.
   \todo: implementation of other boundary conditions
*/
{

// search for aligned tets of face
// if only one tet is found we are on the boundary
   /* static int mycnt = 0;
    mycnt++;
    std::cout << "mycnt = " << mycnt << std::endl;
    std::cout << "&face= " << &face << std::endl;*/

    const TetraCL* tets[2];
    Uint facnum[2];
    Uint cnt=0;
    for (int i = 0; i < 4; ++i)
    {
        const TetraCL* tet ( face.GetNeighbor(i));
        if (tet==NULL || cnt==2) continue;

        if ( tet->IsInTriang ( lvl ) )
        {
            tets[cnt] = tet;
            facnum[cnt++] = face.GetFaceNumInTetra(tet);
        }

    }
    if (cnt == 2){
    // description of face vertex as a barycoord of each neighbor tet
    // takes care of diffent oriantation of the tets
        BaryCoordCL BaryC_tet[2][3]; // description of face vertex as a barycoord of each neighbor tet
        //Point3DCL m(0.75), p3zero(0.), q[3],qm[3];

        for (int i = 0; i < 3; ++i) //vertices of face
            for (int k = 0; k < 2; ++k) // neighbor tets
                for (int j = 0; j < 4; ++j) //vertices of (both) tets
                    if (tets[k]->GetVertex(j)==face.GetVertex(i))
                        BaryC_tet[k][i][j] = 1.0;

        Quad5_2DCL<Point3DCL> velocity;//, velocity_extern;
        velocity.assign( *tets[0], BaryC_tet[0], vel_);
        Point3DCL normal;
        tets[0]->GetOuterNormal(facnum[0],normal);
        Quad5_2DCL<double> vn (dot(normal,velocity));
        LocalP2CL<> lp2;
        LocalP2CL<> mlp2;
        // shape_as_q5: shape function at Quad5 quadrature points on face
        Quad5_2DCL<> shape_as_q5[2][10];
        //upwshape_as_q5: upwind operator applied to shape function at Quad5 quadrature points on face
        Quad5_2DCL<> upwshape_as_q5[2][10];
        for (int j = 0; j < 10; ++j)
        {
            if (j>0){
                lp2[j-1] = 0.0;
                mlp2[j-1] = 0.0;
            }
            lp2[j]  =  1.0;
            mlp2[j] = -1.0;
            for (int i = 0; i < 2; ++i)
                upwshape_as_q5[i][j].assign(lp2, BaryC_tet[i]);
            shape_as_q5[0][j].assign(lp2, BaryC_tet[0]);
            shape_as_q5[1][j].assign(mlp2, BaryC_tet[1]);
        }

      //  Point3DCL normal;
      //  tets[0]->GetOuterNormal(facnum[0],normal);

       // Quad5_2DCL<Point3DCL> velocity;
       // velocity.assign( *tets[0], BaryC_tet[0], vel_);
        // flow direction relative to face at Quad5 quadrature points
        //Quad5_2DCL<double> vn (dot(normal,velocity));

        //Quad5_2DCL<double> v2 (dot(velocity,velocity));
        //const double v_max = std::max(v2.max(),-v2.min());
        // applying upwind operator
        for (int j = 0; j < 10; ++j)
            for (int k = 0; k < Quad5_2DDataCL::NumNodesC; ++k)
                if (vn[k] > 0)
                    upwshape_as_q5[1][j][k] = 0.0;
                else
                    upwshape_as_q5[0][j][k] = 0.0;

        for (int i = 0; i < 2; ++i)
        {
            n[i].assign( *tets[i], *ls_.Phi.RowIdx, ls_.GetBndData()); //BndData is dummy
        }

        const VertexCL* v[3];
        for (Uint i= 0; i < 3; ++i)
            v[i]= face.GetVertex( i);

        const double absdet_face=FuncDet2D(v[1]->GetCoord()-v[0]->GetCoord(),
                                           v[2]->GetCoord()-v[0]->GetCoord());


        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                for (int k = 0; k < 10; ++k)
                    for (int l = 0; l < 10; ++l)
                    {
                        const double integral = Quad5_2DCL<>(vn*upwshape_as_q5[j][l]*shape_as_q5[i][k]).quad(absdet_face);
                        // test integarl for beeing close to zero to reduce matrix entries
                        //if (std::abs(integral) > 1e-17*v_max)
                            (*bH_)( n[i].num[k], n[j].num[l]) += integral;
                    }
    }
    else // weakly imposed boundary condition case
    {
        if (!tets[0]->IsBndSeg(facnum[0]))
            throw DROPSErrCL("not a boundary --- should not happen");
        BaryCoordCL BaryC_tet[3]; // description of face vertex as a barycoord of each neighbor tet

        for (int i = 0; i < 3; ++i) //vertices of face
            for (int j = 0; j < 4; ++j) //vertices of (both) tets
                if (tets[0]->GetVertex(j)==face.GetVertex(i))
                    BaryC_tet[i][j] = 1.0;

        LocalP2CL<> lp2;

        Point3DCL normal;
        tets[0]->GetOuterNormal(facnum[0],normal);

        Quad5_2DCL<Point3DCL> velocity;
        velocity.assign( *tets[0], BaryC_tet, vel_);
        //velocity.assign( *tets[0], BaryC_tet, fixed_velocity);
        Quad5_2DCL<double> vn (dot(normal,velocity));

        Quad5_2DCL<> shape_as_q5[10];
        Quad5_2DCL<> upwshape_as_q5[10];

        for (int j = 0; j < 10; ++j)
        {
            if (j>0){
                lp2[j-1] = 0.0;
            }
            lp2[j]  =  1.0;
            upwshape_as_q5[j].assign(lp2, BaryC_tet);
            shape_as_q5[j].assign(lp2, BaryC_tet);
        }

        for (int j = 0; j < 10; ++j)
            for (int k = 0; k < Quad5_2DDataCL::NumNodesC; ++k)
                if (vn[k] <= 0)
                    upwshape_as_q5[j][k] = 0.0;

        n[0].assign( *tets[0], *ls_.Phi.RowIdx, ls_.GetBndData()); //BndData is dummy

        const VertexCL* v[3];
        for (Uint i= 0; i < 3; ++i)
            v[i]= face.GetVertex( i);

        const double absdet_face=FuncDet2D(v[1]->GetCoord()-v[0]->GetCoord(),
                                           v[2]->GetCoord()-v[0]->GetCoord());


        for (int k = 0; k < 10; ++k)
            for (int l = 0; l < 10; ++l)
                (*bH_)( n[0].num[k], n[0].num[l])
                    += Quad5_2DCL<>(vn*upwshape_as_q5[l]*shape_as_q5[k]).quad(absdet_face);


        Quad5_2DCL<> bndf;
        bndf.assign( *tets[0], BaryC_tet, ls_.GetBndData().GetBndFun(tets[0]->GetBndIdx(facnum[0])));

        for (int k = 0; k < Quad5_2DDataCL::NumNodesC; ++k)
            if (vn[k] > 0)
                bndf[k] = 0.0;
        // create rhs due to boundary conditions
        for (int k = 0; k < 10; ++k)
                rhs_.Data[ n[0].num[k] ] -= Quad5_2DCL<>(vn*bndf*shape_as_q5[k]).quad(absdet_face);

    }

    // const double absdet_face=FuncDet2D(v[1]->GetCoord()-v[0]->GetCoord(),
    //                                    v[2]->GetCoord()-v[0]->GetCoord());

    // Quad5_2DCL(const TetraCL&, const BaryCoordCL* const, const PFunT&);
}

/// because virtual functions can not be templatized, cases are catched individually
/// and implemented with hard casts depending on FE type
template<class DiscVelSolT>
void LevelsetP2CL::SetupSystem( const DiscVelSolT& vel, __UNUSED__ const double dt)
/// Setup level set matrices E, H
{
    ScopeTimerCL scope("Levelset SetupSystem");
    if (!IsDG)
    {
      LevelsetP2ContCL* asp2cont = dynamic_cast<LevelsetP2ContCL*>(this);
      if (asp2cont)
        asp2cont->SetupSystem(vel,dt);
      else
         throw DROPSErrCL("SetupSystem: cast failed");
    }
    else if (IsDG)
    {
      LevelsetP2DiscontCL* asp2discont = dynamic_cast<LevelsetP2DiscontCL*>(this);
      if (asp2discont)
        asp2discont->SetupSystem(vel,dt);
      else
         throw DROPSErrCL("SetupSystem: cast failed");
    }
    else
      throw DROPSErrCL("SetupSystem: unknown FETYPE");
}

template<class DiscVelSolT>
void LevelsetP2ContCL::SetupSystem( const DiscVelSolT& vel, const double dt)
/// Setup level set matrices E, H
{
    LevelsetAccumulator_P2CL<DiscVelSolT> accu( *this, vel, SD_, dt);
    TetraAccumulatorTupleCL accus;
    MaybeAddProgressBar(MG_, "Levelset Setup", accus, Phi.RowIdx->TriangLevel());
    accus.push_back( &accu);
    accumulate( accus, MG_, Phi.RowIdx->TriangLevel(), Phi.RowIdx->GetBndInfo());
    rhs.SetIdx( &idx.GetFinest());
    // rhs.Data.resize(Phi.RowIdx->NumUnknowns());
}

template<class DiscVelSolT>
void LevelsetP2DiscontCL::SetupSystem( const DiscVelSolT& vel, __UNUSED__ const double dt)
/// Setup level set matrices E, H
/// contains face and tetra accumulators
{
//    throw DROPSErrCL("LevelsetP2DiscontCL::SetupSystem : not yet");

    LevelsetTetraAccumulator_P2DCL<DiscVelSolT> accu_tet( *this, vel);
    TetraAccumulatorTupleCL accus_tet;
    accus_tet.push_back( &accu_tet);
    accumulate( accus_tet, MG_, Phi.RowIdx->TriangLevel(), Phi.RowIdx->GetBndInfo());

    LevelsetFaceAccumulator_P2DCL<DiscVelSolT> accu_face( *this, vel, accu_tet.GetHMatrix(), rhs);
    FaceAccumulatorTupleCL accus_face;
    accus_face.push_back( &accu_face);
    accumulate_faces( accus_face, MG_, Phi.RowIdx->TriangLevel(), Phi.RowIdx->GetBndInfo());
}

template <class DiscVelSolT>
PermutationT LevelsetP2CL::downwind_numbering (const DiscVelSolT& vel, IteratedDownwindCL dw)
{
    ScopeTimerCL scope("Downwind Numbering");
    std::cout << "LevelsetP2CL::downwind_numbering:\n";
    std::cout << "...accumulating convection matrix...\n";
    MatrixCL C;
    DownwindAccu_P2CL accu( *vel.GetBndData(), *vel.GetSolution(), idx.GetFinest(), C);
    TetraAccumulatorTupleCL accus;
    MaybeAddProgressBar(MG_, "Downwind Numbering", accus, Phi.RowIdx->TriangLevel());
    accus.push_back( &accu);
    accumulate( accus, this->GetMG(), idx.TriangLevel(), idx.GetBndInfo());

    const PermutationT& p= dw.downwind_numbering( C);
    permute_fe_basis( GetMG(), idx.GetFinest(), p);
    permute_Vector( Phi.Data, p);
    std::cout << "...downwind numbering finished.\n";

    return p;
}

} // end of namespace DROPS

