/// \file spacetime_setup.h
/// \brief setup routines for space time xfem
/// \author LNM RWTH Aachen: Christoph Lehrenfeld

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
 * Copyright 2012 LNM/SC RWTH Aachen, Germany
*/

#ifndef DROPS_SPACETIME_SETUP_H
#define DROPS_SPACETIME_SETUP_H

#include "misc/problem.h"
#include "misc/params.h"
#include "num/spacetime_geom.h"
#include "num/spacetime_quad.h"
#include "spacetimetransp/stxfem.h"
#include "num/interfacePatch.h"
#include "num/discretize.h"
#include "num/accumulator.h"
#include "num/quadrature.h"
#include "levelset/levelset.h"

#include <iostream>
#include <fstream>

namespace DROPS
{

const Uint quad1d_max_n_points = 4;

const Uint quad1d_p3_n_points = 4;
const double quad1d_p3_points[4] =
    { 0.0, 1.0/3.0,  2.0/3.0, 1.0 };
const double quad1d_p3_modal_psik_times_psi[4][2] =
    { {13.0/120.0 ,  1.0/ 60.0},
      { 3.0/ 10.0 ,  3.0/ 40.0},
      { 3.0/ 40.0 ,  3.0/ 10.0},
      { 1.0/ 60.0 , 13.0/120.0}};     // int u*v, u in P1, v in P2
const double quad1d_p3_modal_psik_times_psi_i_psi_j[4][2][2] =
    { { { 1./ 10.0, 1./120.0} , {1./120.0, 1./120.0} },
      { { 9./ 40.0, 3./ 40.0} , {3./ 40.0,       0.} },
      { {      0.0, 3./ 40.0} , {3./ 40.0, 9./ 40.0} },
      { { 1./120.0, 1./120.0} , {1./120.0, 1./ 10.0} } }; // int u*v*w, u in P1, v in P1, w in P2
//naiv backup:
const double quad1d_p3_intweights[4] =
    { 1.0/8.0, 3.0/8.0, 3.0/8.0, 1.0/8.0};

const Uint quad1d_p2_n_points = 3;
const double quad1d_p2_points[3] =
    { 0.0, 0.5, 1.0 };
const double quad1d_p2_modal_psik_times_psi[3][2] =
    { {1.0/6.0,     0.0},
      {1.0/3.0, 1.0/3.0},
      {    0.0, 1.0/6.0} };     // int u*v, u in P1, v in P2
const double quad1d_p2_modal_psik_times_psi_i_psi_j[3][2][2] =
    { { { 9./60.0, 1./60.0} , {1./60.0, -1./60.0} },
      { {12./60.0, 8./60.0} , {8./60.0, 12./60.0} },
      { {-1./60.0, 1./60.0} , {1./60.0,  9./60.0} } }; // int u*v*w, u in P1, v in P1, w in P2
const double quad1d_p2_intweights[3] =
    { 1.0/6.0, 2.0/3.0, 1.0/6.0 };

const Uint quad1d_p1_n_points = 2;
const double quad1d_p1_points[2] =
    { 0.0, 1.0};
const double quad1d_p1_modal_psik_times_psi[2][2] =
    { {1.0/3.0, 1.0/6.0} , {1.0/6.0, 1.0/3.0}};     //     u * v
const double quad1d_p1_modal_psik_times_psi_i_psi_j[2][2][2] =
    { { { 1.0/4.0, 1.0/12.0} , {1.0/12.0, 1.0/12.0} },      // u * v * w
      { {1.0/12.0, 1.0/12.0} , {1.0/12.0, 1.0/4.00} } };
const double quad1d_p1_intweights[2] =
    { 0.5, 0.5};

const Uint quad1d_p0_n_points = 1;
const double quad1d_p0_points[1] =
    { 0.5};
const double quad1d_p0_modal_psik_times_psi[1][2] =
    { {1.0/2.0, 1.0/2.0} };     //     u * v
const double quad1d_p0_modal_psik_times_psi_i_psi_j[1][2][2] =
    { { { 1.0/3.0, 1.0/6.0} , {1.0/6.0, 1.0/3.0} } };
const double quad1d_p0_intweights[1] =
    { 1.0};

class MomentsQuad1D
{
public:
    
    MomentsQuad1D(){;}

    static int GetNPoints(int order){
        return order+1;
    }

    static double GetPoint(int order, int i){
        switch(order)
        {
        case 0:
            return quad1d_p0_points[i];
            break;
        case 2:
            return quad1d_p2_points[i];
            break;
        case 3:
            return quad1d_p3_points[i];
            break;
        case 1:
        default:
            return quad1d_p1_points[i];
        }
    }

    static double GetIntegralPkP1P1(int order, int l, int i, int j) // i,j P1 - l Pk
    {
        switch(order)
        {
        case 0:
            return quad1d_p0_modal_psik_times_psi_i_psi_j[l][i][j];
            break;
        case 1:
            return quad1d_p1_modal_psik_times_psi_i_psi_j[l][i][j];
            break;
        case 2:
            return quad1d_p2_modal_psik_times_psi_i_psi_j[l][i][j];
            break;
        case 3:
            return quad1d_p3_modal_psik_times_psi_i_psi_j[l][i][j];
            break;        
        default:
            throw DROPSErrCL("MomentsQuad1D::GetIntegralPkP1P1: Unknown integration order!");
            return 0;
        }
    }

    static double GetIntegralPkP1(int order, int l, int i) // i P1 - l Pk
    {
        switch(order)
        {
        case 0:
            return quad1d_p0_modal_psik_times_psi[l][i];
            break;
        case 1:
            return quad1d_p1_modal_psik_times_psi[l][i];
            break;
        case 2:
            return quad1d_p2_modal_psik_times_psi[l][i];
            break;
        case 3:
            return quad1d_p3_modal_psik_times_psi[l][i];
            break;
        default:
            throw DROPSErrCL("MomentsQuad1D::GetIntegralPkP1: Unknown integration order!");
            return 0;
        }
    }

    static double GetIntegrationWeight(int order, int i)
    {
        switch(order)
        {
        case 0:
            return quad1d_p0_intweights[i];
            break;
        case 2:
            return quad1d_p2_intweights[i];
            break;
        case 3:
            return quad1d_p3_intweights[i];
            break;
        case 1:
        default:
            return quad1d_p1_intweights[i];
        }
    }
};


const double mass_1D_p1p1[2][2] = { {1.0/3.0, 1.0/6.0} , {1.0/6.0, 1.0/3.0}};     //     u * v
const double dudtv_1D_p1p1[2][2] = { {-1.0/2.0, 1.0/2.0} , {-1.0/2.0, 1.0/2.0}};  // du/dt * v
const double dudtdvdt_1D_p1p1[2][2] = { {1.0, -1.0} , {-1.0, 1.0}};  // du/dt * dv/dt
const double m_udvdt_1D_p1p1[2][2] = { {1.0/2.0, 1.0/2.0} , {-1.0/2.0, -1.0/2.0}};//    -u * dv/dt
const double dudtv_0_5_m_udvdt_0_5_1D_p1p1[2][2] = 
    { { 0.0, 0.5}, 
      {-0.5, 0.0} };   // 0.5*u * dv/dt - 0.5 * v * du/dt
const double udvdt_1D_p1p1[2][2] = 
    { {-1.0/2.0, -1.0/2.0}, 
      { 1.0/2.0,  1.0/2.0} };  //     u * dv/dt
const double mass_1D_p1p1p1[2][2][2] =
    { { { 1.0/4.0, 1.0/12.0} , {1.0/12.0, 1.0/12.0} },      // u * v * w
      { {1.0/12.0, 1.0/12.0} , {1.0/12.0, 1.0/4.00} } };
const double mass_1D_p2p1[3][2] = 
    { {1.0/6.0,     0.0}, 
      {1.0/3.0, 1.0/3.0},
      {    0.0, 1.0/6.0} };     //     u * v, u in P1, v in P2

const double mass_1D_p2p1p1[3][2][2] = 
    { { { 9./60.0, 1./60.0} , {1./60.0, -1./60.0} },
      { {12./60.0, 8./60.0} , {8./60.0, 12./60.0} },
      { {-1./60.0, 1./60.0} , {1./60.0,  9./60.0} } }; // u*v*w, u in P1, v in P1, w in P2


class STVelocityContainer
{
  typedef BndDataCL<Point3DCL>                              VelBndDataT;
  typedef P2EvalCL<SVectorCL<3>, const VelBndDataT, const VecDescCL> const_DiscVelSolCL;  
  private:
    VecDescCL *vold_;   
    VecDescCL *vnew_;   
    const VelBndDataT*    Bnd_v_;     
    MultiGridCL*   MG_;  
    instat_vector_fun_ptr vfptr_; 
    const_DiscVelSolCL * asp2_old;
    const_DiscVelSolCL * asp2_new;
  public:
    STVelocityContainer(VecDescCL & v_old, VecDescCL & v_new,
                      const VelBndDataT& Bnd_v,MultiGridCL& MG)
        :vold_(&v_old),vnew_(&v_new),Bnd_v_(&Bnd_v), MG_(&MG),vfptr_(0)
    {
      asp2_old = new const_DiscVelSolCL( vold_, Bnd_v_, MG_);
      asp2_new = new const_DiscVelSolCL( vnew_, Bnd_v_, MG_);
    };
    
    STVelocityContainer(instat_vector_fun_ptr v):vold_(0),vnew_(0),Bnd_v_(0),MG_(0),vfptr_(v),asp2_old(0),asp2_new(0){};
    
    ~STVelocityContainer()
    {
      if (asp2_old) delete asp2_old;
      if (asp2_new) delete asp2_new;
    }
    
    const_DiscVelSolCL & GetVelocityAsP2_Old() const
        { 
          if (!(vold_ && Bnd_v_))
            throw DROPSErrCL("velocity not prescribed as a const_DiscVelSolCL");
          return *asp2_old; 
        }
        
    const_DiscVelSolCL & GetVelocityAsP2_New() const
        { 
          if (!(vnew_ && Bnd_v_))
            throw DROPSErrCL("velocity not prescribed as a const_DiscVelSolCL");
          return *asp2_new; 
        }
        
    const_DiscVelSolCL GetVelocityAsP2_Old(const VecDescCL& vel) const
        { 
          if (!(vold_ && Bnd_v_))
            throw DROPSErrCL("velocity not prescribed as a const_DiscVelSolCL");
          return const_DiscVelSolCL( &vel, Bnd_v_, MG_); 
        }

    instat_vector_fun_ptr GetVelocityAsFunctionPointer() const
    {
      if (!vfptr_)
        throw DROPSErrCL("velocity not prescribed as a function(pointer)");
      return vfptr_;
    }
    
    bool VelocityAsP2() const {
      return (vold_ && vnew_ && Bnd_v_);
    }
     
    bool VelocityAsFunctionPointer() const {
      return (vfptr_);
    }
 
};



template <int D>
inline void composite_elcontrib_to_xelcontrib(double coup_s_s_neg[D][D], double coup_s_s_pos[D][D],
                                      double elvec_neg[D], double elvec_pos[D],
                                      bool *sign,
                                      double intf_coup_s_s[D][D], double intf_coup_s_x[D][D],
                                      double intf_coup_x_s[D][D], double intf_coup_x_x[D][D],
                                      double coup_s_s[D][D], double coup_s_x[D][D],
                                      double coup_x_s[D][D], double coup_x_x[D][D],
                                      double elvec[D], double elvec_x[D])
{
    for(Uint i=0; i<D; ++i)
    {
        if (elvec_neg != NULL || elvec_pos != NULL)
        {
            elvec[i] = elvec_neg[i] + elvec_pos[i];
            elvec_x[i] = sign[i] ? -elvec_neg[i] : elvec_pos[i];
        }
         for(Uint j=0; j<D; ++j)
        {
            coup_s_s[i][j] = coup_s_s_neg[i][j] + coup_s_s_pos[i][j];
            coup_x_s[i][j] = sign[i] ? -coup_s_s_neg[i][j] : coup_s_s_pos[i][j];
            coup_s_x[i][j] = sign[j] ? -coup_s_s_neg[i][j] : coup_s_s_pos[i][j];
            if (sign[i] == sign[j])
                coup_x_x[i][j] = sign[j] /*( == sign[i] )*/ ? coup_s_s_neg[i][j] : coup_s_s_pos[i][j];
            else
                coup_x_x[i][j] = 0.0;
        }
    }

    for(Uint i=0; i<D; ++i)
    {
        for(Uint j=0; j<D; ++j)
        {
            coup_s_s[i][j] += intf_coup_s_s[i][j];
            coup_s_x[i][j] += intf_coup_s_x[i][j];
            coup_x_s[i][j] += intf_coup_x_s[i][j];
            coup_x_x[i][j] += intf_coup_x_x[i][j];
        }
    }
}


// Accumulator for space time integrals - testing surface area, space-time surface area etc... 
class STGeomApproxTestAccumulatorCL : public TetraAccumulatorCL
{
    protected:   
    const MultiGridCL& MG_;

    const LevelsetP2CL * lsetp2old;
    const LevelsetP2CL * lsetp2new;
    instat_scalar_fun_ptr lset_fpt;

    double det;
    double absdet;
    SMatrixCL<3,3> T;
    bool sign[8]; //sign of levelset at vertices (4(space) x 2(time))
    const Uint ints_per_space_edge;
    const Uint subtimeintervals;
    const double told;
    const double tnew;
    double iscut;
    double * surfarea;
    double * stsurfarea;
    double * negvol;
    double * posvol;
    Point3DCL * intnorm;
    Point4DCL * intnorm4d;

public: 
    STGeomApproxTestAccumulatorCL (const MultiGridCL& MG, const LevelsetP2CL * lsetp2old_in,
                                   const LevelsetP2CL * lsetp2new_in, 
                                   instat_scalar_fun_ptr lset_fpt, const double t1, const double t2,
                                   const ParamCL::ptree_type & P);
    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();
    virtual void visit (const TetraCL&); 
    virtual TetraAccumulatorCL* clone (int /*tid*/) { return new STGeomApproxTestAccumulatorCL( *this); }
};



// Accumulator for space time integrals - measure error in jump [ beta  u ], different versions..
class InterfaceJumpAccumulatorCL : public TetraAccumulatorCL
{
    protected:
    const MultiGridCL& MG_;

    const LevelsetP2CL * lsetp2old;
    const LevelsetP2CL * lsetp2new;
    instat_scalar_fun_ptr lset_fpt;

    double det;
    double absdet;
    SMatrixCL<3,3> T;
    bool sign[8]; //sign of levelset at vertices (4(space) x 2(time))

    const Uint ints_per_space_edge;
    const Uint subtimeintervals;

    const double told;
    const double tnew;

    double iscut;
    double * ifjump;
    double * ifjump_past;
    double * st_ifjump;
    double * st_ifjump_nu;

    const VecDescCL & old_sol_neg;
    const VecDescCL & old_sol_pos;
    const VecDescCL & new_sol_neg;
    const VecDescCL & new_sol_pos;
    const BndDataCL<> & Bnd_neg;
    const BndDataCL<> & Bnd_pos;

    const double beta_neg;
    const double beta_pos;
    const double lambda_stab;

public:
    InterfaceJumpAccumulatorCL (const MultiGridCL& MG, const LevelsetP2CL * lsetp2old_in,
                                const LevelsetP2CL * lsetp2new_in,
                                instat_scalar_fun_ptr lset_fpt, const double t1, const double t2,
                                const VecDescCL & oldsol_neg_in, const VecDescCL & oldsol_pos_in,
                                const VecDescCL & newsol_neg_in, const VecDescCL & newsol_pos_in,
                                const BndDataCL<> & Bnd_neg_in,
                                const BndDataCL<> & Bnd_pos_in,
                                const ParamCL::ptree_type & P);
    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();
    virtual void visit (const TetraCL&);
    virtual TetraAccumulatorCL* clone (int /*tid*/) { return new InterfaceJumpAccumulatorCL( *this); }
};




// Accumulator for space time integrals
class STVolumeAccumulator_P1SP1TXCL : public TetraAccumulatorCL
{
    protected:   
    const MultiGridCL& MG_;
    const BndDataCL<> * BndData_neg; 
    const BndDataCL<> * BndData_pos; 
    MatrixCL* Amat_; 
    VecDescCL* b_; 
    const IdxDescCL& RowIdx_; 
    const IdxDescCL& ColIdx_; 

    const LevelsetP2CL * lsetp2old;
    const LevelsetP2CL * lsetp2new;
    
    MatrixBuilderCL * A_;    
    
    const ExtIdxDescCL& Xidx;

    //local informations
    
    // - sharable (for future changes)
    Point3DCL G[4];
    double det;
    double absdet;
    SMatrixCL<3,3> T;
    IdxT UnknownIdx[4];  // the std. unknowns on one time level (the other follows by +1 (stride.. ))

    LocalP2CL<double> phi[4];
    Quad5CL<double> q_phi[4];
    LocalP2CL<double> phi_i_phi_j[4][4];

    // - not sharable (for future changes)

    double coup_s_s_neg[8][8];
    double coup_s_s_pos[8][8];

    double coup_s_s_space[4][4];

    bool sign[8]; //sign of levelset at vertices (4(space) x 2(time))

    double coup_s_s[8][8];
    double coup_x_s[8][8];
    double coup_s_x[8][8];
    double coup_x_x[8][8];

    double intf_coup_s_s[8][8];
    double intf_coup_x_s[8][8];
    double intf_coup_s_x[8][8];
    double intf_coup_x_x[8][8];

    double elvec_neg[8];
    double elvec_pos[8];
    double elvec[8];
    double elvec_x[8];

    /* QuadCL<> U_Grad[4]; */
    
    const Uint lvl;
    const Uint idx;
    
    const Uint ints_per_space_edge;
    const Uint subtimeintervals;

    const double told;
    const double tnew;

    const double beta_neg;
    const double beta_pos;

    double iscut;

    void transform_comp_elmats_to_xelmats();
    void update_global_matrix_without_x();
    void update_global_matrix_with_x();
    void update_coupling(const TetraCL& sit, bool with_x);
    
public: 
    STVolumeAccumulator_P1SP1TXCL (const MultiGridCL& MG, const BndDataCL<> * BndData_neg_in, 
                                   const BndDataCL<> * BndData_pos_in, 
                                   const LevelsetP2CL * lsetp2old_in,
                                   const LevelsetP2CL * lsetp2new_in,
                                   MatrixCL* Amat, VecDescCL* b, 
                                   const IdxDescCL& RowIdx, const IdxDescCL& ColIdx, 
                                   const double t1, const double t2,
                                   const ParamCL::ptree_type & P);

    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    virtual void local_setup (const TetraCL& sit) = 0;

    virtual void visit (const TetraCL&); 
    void output_elmats (); // just for debugging... Can be removed later..

    virtual TetraAccumulatorCL* clone (int /*tid*/) = 0; //{ return new STVolumeAccumulator_P1SP1TXCL ( *this); }
    
};

// lets-see-matrix 
class MassTestAccumulator_P1SP1TXCL : public STVolumeAccumulator_P1SP1TXCL
{
    protected:
    typedef STVolumeAccumulator_P1SP1TXCL base_;
    using                           base_::MG_;
    using                           base_::BndData_neg; 
    using                           base_::BndData_pos; 
    using                           base_::Amat_; 
    using                           base_::b_; 
    using                           base_::RowIdx_; 
    using                           base_::ColIdx_; 
    using                           base_::lsetp2old; 
    using                           base_::lsetp2new; 
    using                           base_::A_;              //MassTests matrix
    using                           base_::G;
    using                           base_::det;
    using                           base_::absdet;
    using                           base_::T;
    using                           base_::UnknownIdx;
    using                           base_::phi;
    using                           base_::q_phi;
    using                           base_::phi_i_phi_j;
    using                           base_::coup_s_s_neg;
    using                           base_::coup_s_s_pos;
    using                           base_::coup_s_s_space;
    using                           base_::coup_s_s;
    using                           base_::coup_s_x;
    using                           base_::coup_x_s;
    using                           base_::coup_x_x;
    using                           base_::intf_coup_s_s;
    using                           base_::intf_coup_s_x;
    using                           base_::intf_coup_x_s;
    using                           base_::intf_coup_x_x;
    using                           base_::sign;
    using                           base_::lvl;
    using                           base_::idx;
    using                           base_::ints_per_space_edge;
    using                           base_::subtimeintervals;
    using                           base_::told;
    using                           base_::tnew;
    using                           base_::iscut;
    instat_scalar_fun_ptr rhs_neg;
    instat_scalar_fun_ptr rhs_pos;

    public:
    MassTestAccumulator_P1SP1TXCL (const MultiGridCL& MG, const BndDataCL<> * BndData_neg_in, 
                                   const BndDataCL<> * BndData_pos_in, 
                                   const LevelsetP2CL * lsetp2old_in,
                                   const LevelsetP2CL * lsetp2new_in,
                                   instat_scalar_fun_ptr rhs_neg_in,
                                   instat_scalar_fun_ptr rhs_pos_in,
                                   MatrixCL* Amat, VecDescCL* b, 
                                   const IdxDescCL& RowIdx, const IdxDescCL& ColIdx, 
                                   const double t1, const double t2,
                                   const ParamCL::ptree_type & P)
        : base_(MG,BndData_neg_in,BndData_pos_in,lsetp2old_in,lsetp2new_in,Amat,b,RowIdx,ColIdx,t1,t2,P),
        rhs_neg(rhs_neg_in),rhs_pos(rhs_pos_in){}
    virtual void local_setup (const TetraCL& sit);
    virtual TetraAccumulatorCL* clone (int /*tid*/) { return new MassTestAccumulator_P1SP1TXCL ( *this); }
    
};


// diffusion term
class SpatialLaplaceAccumulator_P1SP1TXCL : public STVolumeAccumulator_P1SP1TXCL
{
    protected:
    typedef STVolumeAccumulator_P1SP1TXCL base_;
    using                           base_::MG_;
    using                           base_::BndData_neg; 
    using                           base_::BndData_pos; 
    using                           base_::Amat_; 
    using                           base_::b_; 
    using                           base_::RowIdx_; 
    using                           base_::ColIdx_; 
    using                           base_::lsetp2old; 
    using                           base_::lsetp2new; 
    using                           base_::A_;              //SpatialLaplaces matrix
    using                           base_::G;
    using                           base_::det;
    using                           base_::absdet;
    using                           base_::T;
    using                           base_::UnknownIdx;
    using                           base_::phi;
    using                           base_::q_phi;
    using                           base_::phi_i_phi_j;
    using                           base_::coup_s_s_neg;
    using                           base_::coup_s_s_pos;
    using                           base_::coup_s_s_space;
    using                           base_::coup_s_s;
    using                           base_::coup_s_x;
    using                           base_::coup_x_s;
    using                           base_::coup_x_x;
    using                           base_::intf_coup_s_s;
    using                           base_::intf_coup_s_x;
    using                           base_::intf_coup_x_s;
    using                           base_::intf_coup_x_x;
    using                           base_::sign;
    using                           base_::lvl;
    using                           base_::idx;
    using                           base_::ints_per_space_edge;
    using                           base_::subtimeintervals;
    using                           base_::told;
    using                           base_::tnew;
    using                           base_::iscut;
    instat_scalar_fun_ptr rhs_neg;
    instat_scalar_fun_ptr rhs_pos;

    /* QuadCL<> U_Grad[4]; */

    public:
    SpatialLaplaceAccumulator_P1SP1TXCL (const MultiGridCL& MG, const BndDataCL<> * BndData_neg_in, 
                                         const BndDataCL<> * BndData_pos_in, 
                                         const LevelsetP2CL * lsetp2old_in,
                                         const LevelsetP2CL * lsetp2new_in,
                                         instat_scalar_fun_ptr rhs_neg_in,
                                         instat_scalar_fun_ptr rhs_pos_in,
                                         MatrixCL* Amat, VecDescCL* b, 
                                         const IdxDescCL& RowIdx, const IdxDescCL& ColIdx, 
                                         const double t1, const double t2,
                                         const ParamCL::ptree_type & P)
        : base_(MG,BndData_neg_in,BndData_pos_in,lsetp2old_in,lsetp2new_in,Amat,b,RowIdx,ColIdx,t1,t2,P),
        rhs_neg(rhs_neg_in),rhs_pos(rhs_pos_in){}
    virtual void local_setup (const TetraCL& sit);

    virtual TetraAccumulatorCL* clone (int /*tid*/ ) { return new SpatialLaplaceAccumulator_P1SP1TXCL ( *this); }
    
};

// STTransport
class STTransportVolumeAccumulator_P1SP1TXCL : public STVolumeAccumulator_P1SP1TXCL
{
    protected:
    typedef STVolumeAccumulator_P1SP1TXCL base_;
    using                           base_::MG_;
    using                           base_::BndData_neg; 
    using                           base_::BndData_pos; 
    using                           base_::Amat_; 
    using                           base_::b_; 
    using                           base_::RowIdx_; 
    using                           base_::ColIdx_; 
    using                           base_::lsetp2old; 
    using                           base_::lsetp2new; 
    using                           base_::A_;              //SpatialLaplaces matrix
    using                           base_::G;
    using                           base_::det;
    using                           base_::absdet;
    using                           base_::T;
    using                           base_::UnknownIdx;
    using                           base_::phi;
    using                           base_::q_phi;
    using                           base_::phi_i_phi_j;
    using                           base_::coup_s_s_neg;
    using                           base_::coup_s_s_pos;
    using                           base_::coup_s_s_space;
    using                           base_::coup_s_s;
    using                           base_::coup_s_x;
    using                           base_::coup_x_s;
    using                           base_::coup_x_x;
    using                           base_::intf_coup_s_s;
    using                           base_::intf_coup_s_x;
    using                           base_::intf_coup_x_s;
    using                           base_::intf_coup_x_x;
    using                           base_::sign;
    using                           base_::lvl;
    using                           base_::idx;
    using                           base_::ints_per_space_edge;
    using                           base_::subtimeintervals;
    using                           base_::told;
    using                           base_::tnew;
    using                           base_::beta_neg;
    using                           base_::beta_pos;
    using                           base_::iscut;
    instat_scalar_fun_ptr lset_fpt;
    instat_scalar_fun_ptr rhs_neg;
    instat_scalar_fun_ptr rhs_pos;
    
    /* instat_vector_fun_ptr convection; */
    STVelocityContainer & convection;

    const double alpha_neg;
    const double alpha_pos;
    /* const double beta_neg; */
    /* const double beta_pos; */
    const double delta_beta;
    const double lambda_stab;

    const double vmax;
    /* QuadCL<> U_Grad[4]; */
    const VecDescCL & oldsolution_neg;
    const VecDescCL & oldsolution_pos;
    bool use_stkappa_rule;
    bool use_space_kappa_rule;
    bool antisym_conv;
    double SD_stab;

    const int quad1d_rhs_order;
    const int quad1d_conv_order;

    bool adjoint;

    public:
    STTransportVolumeAccumulator_P1SP1TXCL (const MultiGridCL& MG, const BndDataCL<> * BndData_neg_in, 
                                            const BndDataCL<> * BndData_pos_in, 
                                            const LevelsetP2CL * lsetp2old_in,
                                            const LevelsetP2CL * lsetp2new_in,
                                            instat_scalar_fun_ptr lset_fpt_in,
                                            instat_scalar_fun_ptr rhs_neg_in,
                                            instat_scalar_fun_ptr rhs_pos_in,
                                            STVelocityContainer & convection_in,
                                            MatrixCL* Amat, VecDescCL* b, 
                                            const IdxDescCL& RowIdx, const IdxDescCL& ColIdx, 
                                            const double t1, const double t2,
                                            const VecDescCL & oldsol_neg, const VecDescCL & oldsol_pos,
                                            const ParamCL::ptree_type & P);
    virtual void local_setup (const TetraCL& sit);

    virtual TetraAccumulatorCL* clone (int /*tid*/ ) { return new STTransportVolumeAccumulator_P1SP1TXCL ( *this); }
    
};


}//end of namespace DROPS

#endif
