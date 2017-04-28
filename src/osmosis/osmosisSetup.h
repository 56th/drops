/// \file osmosisSetup.h
/// \brief Classes that constitute a osmosis problem.
/// \author Christoph Lehrenfeld IGPM, Thorolf Schulte

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

#ifndef DROPS_TRANSPORTNITSCHE_H
#define DROPS_TRANSPORTNITSCHE_H

#include "levelset/levelset.h"
#include "levelset/mgobserve.h"
#include "misc/problem.h"
#include "num/fe.h"
#include "num/discretize.h"
#include "num/krylovsolver.h"
#include "num/precond.h"
#include "num/bndData.h"
#include "stokes/instatstokes2phase.h"
#include "misc/xfem.h"
#include "misc/params.h"
#include <iostream>
#include <numeric>
#include <cstring>

extern DROPS::ParamCL P;
typedef double (*instat_fun_ptr)(const DROPS::Point3DCL&, double);


namespace DROPS
{
  
class VelocityContainer
{
  typedef BndDataCL<Point3DCL>                              VelBndDataT;
  typedef P2EvalCL<SVectorCL<3>, const VelBndDataT, const VecDescCL> const_DiscVelSolCL;  
  private:
    VecDescCL *v_;   
    const VelBndDataT*    Bnd_v_;     
    MultiGridCL*   MG_;  
    instat_vector_fun_ptr vfptr_; 
    const_DiscVelSolCL * asp2;
  public:
    VelocityContainer(VecDescCL & v,const VelBndDataT& Bnd_v,MultiGridCL& MG):v_(&v),Bnd_v_(&Bnd_v), MG_(&MG),vfptr_(0)
    {
      asp2 = new const_DiscVelSolCL( v_, Bnd_v_, MG_);
    };
    
    VelocityContainer(instat_vector_fun_ptr v):v_(0),Bnd_v_(0),MG_(0),vfptr_(v),asp2(0){};
    
    ~VelocityContainer()
    {
      if (asp2) delete asp2;
    }
    
    const_DiscVelSolCL & GetVelocityAsP2() const
        { 
          if (!(v_ && Bnd_v_))
            throw DROPSErrCL("velocity not prescribed as a const_DiscVelSolCL");
          return *asp2; 
        }
        
    const_DiscVelSolCL GetVelocityAsP2(const VecDescCL& vel) const
        { 
          if (!(v_ && Bnd_v_))
            throw DROPSErrCL("velocity not prescribed as a const_DiscVelSolCL");
          return const_DiscVelSolCL( &vel, Bnd_v_, MG_); 
        }

    instat_vector_fun_ptr GetVelocityAsFunctionPointer() const
    {
      if (!vfptr_)
        throw DROPSErrCL("velocity not prescribed as a function(pointer)");
      return vfptr_;
    }
    
    bool hasVelocityAsP2() const {
      return (v_ && Bnd_v_);
    }
     
    bool hasVelocityAsFunctionPointer() const {
      return (vfptr_);
    }
 
};

class OsmosisP1CL
{
  public:
    typedef BndDataCL<>                                       BndDataT;
    typedef BndDataCL<Point3DCL>                              VelBndDataT;
    typedef P1EvalCL<double, const BndDataT, VecDescCL>       DiscSolCL;
    typedef P1EvalCL<double, const BndDataT, const VecDescCL> const_DiscSolCL;
    typedef P2EvalCL<SVectorCL<3>, const VelBndDataT, const VecDescCL> const_DiscVelSolCL;
    typedef P1EvalCL<SVectorCL<3>, const VelBndDataT, const VecDescCL> const_DiscVelP1SolCL;

    double oldt_, t_;
    MLIdxDescCL idx, oldidx;
    IdxDescCL Velidx, oldVelidx;
    VecDescCL conc, oldconc;   ///<concentration
    VecDescCL Vn_, oldVn_;    ///< velocity field around interface
    MLMatDescCL A,         ///< diffusion matrix
                M,         ///< mass matrix
                C;        ///< convection matrix
    VecDescCL	b;
    
    MatDescCL MVn, ///< Matrix for interface velocity field problem
    		  MV; ///< Matrix for globalised velocity field problem
    VecDescCL bVn, ///< RHS vector for interface velocity field problem
    		  bV; ///< RHS vector for globalised velocity field problem


  private:
    MLMatrixCL     L_;              ///< sum of matrices
    MatrixCL	   LsetL_;
    MultiGridCL&   MG_;
    BndDataT&      Bnd_;     ///< Boundary condition for the concentration
    const VelBndDataT&   Bnd_v_;
    LsetBndDataCL& Bnd_ls_;          ///< Boundary condition for the level set function
    VelocityContainer v_;
    double         theta_,          ///< time scheme parameter
                   dt_,             ///< time step 
                   D_;           ///< diffusion constant
    LevelsetP2CL &lset_,  &oldlset_;  ///< levelset at current time step and previous step
    SSORPcCL				ssorpc_;
    //GMResSolverCL<GSPcCL>   gm_;
    //DummyPcCL               pc_;
    GMResSolverCL<GSDiag0PcCL>   gm_;
    /* GMResSolverCL<GSDiag0PcCL>   gm_; */
    PCGSolverCL<SSORPcCL> 	cg_;
    GSDiag0PcCL                  pc_;
    instat_scalar_fun_ptr f_;         ///<source term
    instat_scalar_fun_ptr c_;        ///<mass/reaction term
    const double omit_bound_;

    void SetupInstatSystem(MatrixCL&, MatrixCL&, MatrixCL&, VecDescCL*, IdxDescCL&, const double) const;
    void SetupInstatMixedMassMatrix( MatrixCL&, VecDescCL*, IdxDescCL&, IdxDescCL&, const double) const;
    void SetupInstatMixedSystem(MatrixCL&, VecDescCL*, MatrixCL&, VecDescCL*, MatrixCL&, VecDescCL*, VecDescCL*,
        IdxDescCL&, IdxDescCL&, const double) const;

  public:
    OsmosisP1CL( MultiGridCL& mg, BndDataT& Bnd, VelBndDataT& BndVel, VelocityContainer& v, LsetBndDataCL& Bnd_ls,
        LevelsetP2CL& lset, LevelsetP2CL& oldlset,
        DROPS::ParamCL& P, double initialtime=0, instat_scalar_fun_ptr reac=0, instat_scalar_fun_ptr rhs=0)
        : oldt_(initialtime), t_( initialtime), 
        idx( P1_FE, 1, Bnd), oldidx( P1_FE, 1, Bnd),
        Velidx( vecP1_FE, BndVel),
        MG_( mg), Bnd_( Bnd), Bnd_v_(BndVel), Bnd_ls_(Bnd_ls), v_ (v),
        theta_( P.get<double>("Time.Theta")), dt_( P.get<int>("Time.NumSteps")!=0 ? P.get<double>("Time.FinalTime")/P.get<int>("Time.NumSteps") : 0),
        D_( P.get<double>("Osmosis.Diffusivity")),
        lset_( lset), oldlset_(oldlset),
        gm_( pc_, 20, P.get<int>("Solver.Iter"), P.get<double>("Solver.Tol"), false, false, RightPreconditioning),
        cg_( ssorpc_, P.get<int>("Solver.Iter"), P.get<double>("Solver.Tol"), false),
        f_(rhs), c_(reac),
        omit_bound_(P.get<double>("Osmosis.XFEMReduced"))
    {
        if (theta_ >= 1.0 || theta_ <= 0.0)
        {
        	std::cerr<< "OsmosisP1CL::OsmosisP1CL: Sorry. Cannot deal with theta != (0,1). Overwriting: Using Theta = 0.5 now! \n";
        	theta_ = 0.5;
        }

    }

    const VelocityContainer & GetVelocity() const {return v_;}

    MultiGridCL& GetMG()               { return MG_; }
    const MultiGridCL& GetMG() const   { return MG_; }
    //GMResSolverCL<GSPcCL>& GetSolver() { return gm_; }
    GMResSolverCL<GSDiag0PcCL>& GetSolver() { return gm_; }
    const BndDataT& GetBndData() const { return Bnd_; }

    void CreateNumbering( Uint level, MLIdxDescCL* idx1, MLIdxDescCL *oldidx1, IdxDescCL *velidx, IdxDescCL *oldvelidx, VecDescCL& lset1, VecDescCL& oldlset1) {
         idx1->CreateNumbering( level, MG_, &lset1, &Bnd_ls_);
         oldidx1->CreateNumbering( level, MG_, &oldlset1, &Bnd_ls_);
         velidx->CreateNumbering( MG_.GetLastLevel(), MG_, &lset1, &Bnd_ls_);
         oldvelidx->CreateNumbering( MG_.GetLastLevel(), MG_, &lset1, &Bnd_ls_);
    }
    void DeleteNumbering( MLIdxDescCL* idx1) { idx1->DeleteNumbering( MG_); }
    /// initialize transformed concentration function
    void Init( instat_scalar_fun_ptr, double t=0.0);
    // scale...
    void InitWithScaling( instat_scalar_fun_ptr, double totalconcentration, double t=0.0);
    /// Set one P1X Function as the other with according scaling parameters in the positive and negative domain parts

    void SetTimeStep( double dt, double theta=-1);
    void SetTwoStepIdx();
    void SetNewIdx();
    void SetupInstatSystem( MLMatDescCL&, MLMatDescCL&, MLMatDescCL&, VecDescCL&, const double) const;
    void SetupInstatMixedMassMatrix( MLMatDescCL&, VecDescCL&, const double) const;
    void SetupInstatMixedSystem( MLMatDescCL&, VecDescCL&, MLMatDescCL&, VecDescCL&, MLMatDescCL&, VecDescCL&, VecDescCL&, const double) const;
    double GetXFEMOmitBound(){ return omit_bound_; }

    //Quad3CL<Point3DCL> GetImprovedTetraNormal(TetraCL& t) const;

    //Velocity field at interface
    void SetupInterfaceVelSystem(MatDescCL& , VecDescCL& , IdxDescCL&, double );
    //Stabilisation for velocity field
    void SetupInterfaceStabilisation(TetraCL& sit, double detTet, MatrixBuilderCL& matM, LocalNumbP1CL& locIdx, Point3DCL Grad[4]);

    /// perform one time step
    void DoStep( double new_t);
    void DoInterfaceStep(double& t);
    void DoLevelsetStep(double& t);
    void SetRHS( instat_scalar_fun_ptr rhs) {f_= rhs;}
    const_DiscSolCL GetSolution() const
        { return const_DiscSolCL( &conc, &Bnd_, &MG_); }
    const_DiscSolCL GetSolution( const VecDescCL& Myc) const
        { return const_DiscSolCL( &Myc, &Bnd_, &MG_);}
    const_DiscVelSolCL GetVelSolution( const VecDescCL& Myc) const
		{ return const_DiscVelSolCL( &Myc, &Bnd_v_, &MG_);}
    const_DiscVelP1SolCL GetVelP1Solution( const VecDescCL& Myc) const
		{ return const_DiscVelP1SolCL( &Myc, &Bnd_v_, &MG_);}
                         
 
    void CheckDivergence();

   
    /// \name For internal use only
    /// The following member functions are added to enable an easier implementation
    /// of the locling navstokes-levelset. They should not be called by a common user.
    /// Use LevelsetP2CL::DoStep() instead.
    ///@{
    void InitStep (VectorCL&);
    void DoStep (const VectorCL&);
    void CommitStep ();
    ///@}
    /// Get FE type for pressure space
    VecDescCL& GetLevelset ()  {return lset_.Phi;}
    const VecDescCL& GetLevelset ()  const {return lset_.Phi;}
    LsetBndDataCL& GetLevelsetBnd ()  {return Bnd_ls_;}
    const LsetBndDataCL& GetLevelsetBnd ()  const {return Bnd_ls_;}
    VecDescCL& GetOldLevelset ()  {return oldlset_.Phi;}
    const VecDescCL& GetOldLevelset ()  const {return oldlset_.Phi;}
    //@}
    void SetupInstatRhs( VecDescCL & b,  const double time) const;
    double MeanDropConcentration(bool absolute=false);
    void VolumeAndMass(double& vol, double& con);
    void SolutionErrorCase1(double& , double& , double );
    void SolutionErrorCase2(double& , double& , double );
    void SolutionErrorCase3(double& , double& , double );
    void SolutionErrorCase4(double& , double& , double );
    double CheckSolution(instat_scalar_fun_ptr Lsgn, instat_scalar_fun_ptr Lsgp, double time);

    //special osmosis functions
    instat_fun_ptr InterfaceVel(instat_fun_ptr kappa, instat_fun_ptr ct) const;
};

class GNUPlotCL
{
	public:
		double relErr;
		double absErr;

		GNUPlotCL(const char* s, OsmosisP1CL* osm, void (OsmosisP1CL::*f)(double&, double&, double))
		{
			gnup.open(s, std::ios::trunc  | std::ios::out);
			gnup << "# Time \t " << "# Out1 \t " << "# Out2 \t " << std::endl;
			dataFun = f;
			osmosis = osm;
		}
		void Write(double t)
		{
			(osmosis->*dataFun)(relErr, absErr, t);
			gnup << t << "\t " << relErr << "\t " << absErr << "\t " << std::endl;
		}

		void Close()
		{
			gnup.close();
		}

	private:
		std::fstream gnup;
		OsmosisP1CL* osmosis;
		void (OsmosisP1CL::*dataFun)(double&, double&, double);
};

/// \brief Observes the MultiGridCL-changes by AdapTriangCL to repair the function c.ct.
///
/// The actual work is done in post_refine().

class OsmosisRepairCL : public MGObserverCL
{
  private:
    OsmosisP1CL& c_;
    std::unique_ptr<P1XRepairCL> oldp1xrepair_;
    Uint mylvl;
  public:
    OsmosisRepairCL (OsmosisP1CL& c, Uint amylvl)
        : c_( c), mylvl(amylvl) {}

    void pre_refine  () {}
    void post_refine ();

    void pre_refine_sequence  (); 
    void post_refine_sequence (); 
    const IdxDescCL* GetIdxDesc() const;
};

class VelOsmosisRepairCL : public MGObserverCL
{
  public:
    typedef BndDataCL<Point3DCL>                                          VelBndDataT;
    typedef P2EvalCL<SVectorCL<3>, const VelBndDataT, const VecDescCL>    const_DiscVelSolCL;
  private:
    MultiGridCL& mg_;
    VecDescCL& v_;
    const VelBndDataT& Bnd_v_;
    IdxDescCL& vidx_;
    double time_;

  public:
    VelOsmosisRepairCL (VecDescCL& v, MultiGridCL& mg, const VelBndDataT& Bnd_v, IdxDescCL& vidx, double t )
        :  mg_(mg), v_(v), Bnd_v_(Bnd_v), vidx_(vidx) , time_(t){}

    void pre_refine  () {}
    void post_refine ();

    void pre_refine_sequence  () {}
    void post_refine_sequence () {}
    void SetTime(double t) {time_=t;}
    const IdxDescCL* GetIdxDesc() const;
};


/// \todo: think about where to put this in a more systematical way

class P1FEGridfunctions{
  private:
    LocalP1CL<> p1[4];
    LocalP2CL<> p2[4];
    Quad3CL<> q3_p[4];
    Quad5CL<> q5_p[4];
    LocalP2CL<> pipj[4][4];
  public:
    static void SetupPiPj(LocalP2CL<>pipj[4][4])
    {    
        for(int i= 0; i < 4; ++i) {
            for(int j= 0; j < i; ++j) {
                pipj[j][i][EdgeByVert( i, j) + 4]= 0.25;
                pipj[i][j][EdgeByVert( j, i) + 4]= 0.25;
            }
            pipj[i][i][i]= 1.;
            for (int vert= 0; vert < 3; ++vert)
                pipj[i][i][EdgeByVert( i, VertOfFace( i, vert)) + 4]= 0.25;
        }
    }

    P1FEGridfunctions(){
      for(int i= 0; i < 4; ++i) {
          p1[i][i]=1.;
          q3_p[i].assign(p1[i]);
          q5_p[i].assign(p1[i]);
          p2[i].assign(p1[i]);
       }
      SetupPiPj(pipj);  
    }
    
    LocalP1CL<>& GetShapeAsLocalP1(int i) {return p1[i];}
    LocalP2CL<>& GetShapeAsLocalP2(int i) {return p2[i];}
    LocalP2CL<>& GetProductShapeAsLocalP2(int i, int j) {return pipj[i][j];}
    Quad3CL<>& GetShapeAsQuad3(int i) {return q3_p[i];}
    Quad5CL<>& GetShapeAsQuad5(int i) {return q5_p[i];}
};

class LocalConvDiffReacCoefficients; //forward declaration

class TransformedP1FiniteElement{
  private:
    SMatrixCL<3,4> G;
    SMatrixCL<4,4> GTG;
    double det, absdet;
    double vol;
  protected:
    TetraCL* tet;
    bool has_trafo_base;
    bool has_Gram;
    bool oninterface;
    P1FEGridfunctions& p1fegfs;
    BaryCoordCL* nodes;
    Quad3CL<> q3_baseshape[4];
  public:
    TransformedP1FiniteElement(P1FEGridfunctions& ap1fegfs, TetraCL* atet = nullptr):tet(atet), p1fegfs(ap1fegfs){
      has_trafo_base=false;    
      has_Gram=false;    
      oninterface = false; 
      nodes = nullptr;
    }
    
    virtual ~TransformedP1FiniteElement(){
      if (nodes) delete[] nodes;
    }
    
    P1FEGridfunctions& GetGridfunctions(){
      return p1fegfs;
    }
    
    void SetTetra(TetraCL& atet){
      tet = &atet;
      has_trafo_base=false;    
      has_Gram=false;    
      oninterface = false;    
    }
    
    virtual void SetLocal(TetraCL& atet, LocalConvDiffReacCoefficients& , bool ){    
      SetTetra(atet);
    }

    void SetSubTetra(const SArrayCL<BaryCoordCL,4>& cutT){
      vol = VolFrac(cutT) * absdet * 1.0/6.0;
      if (nodes) delete[] nodes;
      nodes = Quad3CL<>::TransformNodes(cutT);
      oninterface = true;
      
      for (int i = 0; i < 4; i ++){
        q3_baseshape[i] = Quad3CL<>(p1fegfs.GetShapeAsLocalP2(i), nodes);      
      }      
    }

    TetraCL& GetTetra() const{
      if(tet == nullptr) throw DROPSErrCL("TransformedP1FiniteElement::GetTetra - No TetraCL object given!");
      return *tet;
    }
    
    void CalcTrafoBase(){
      P1DiscCL::GetGradients(G, det, GetTetra());
      absdet = std::fabs( det);
      vol = absdet * 1.0/6.0;
      has_trafo_base = true;
    }
    
    BaryCoordCL* GetNodes() const{
      if (!oninterface)
        throw DROPSErrCL("GetNodes should only be called if a tetra is subdivided!");
      else
        return nodes;
    }
    
    double GetDeterminant(){
      if (!has_trafo_base) CalcTrafoBase();
      return det;
    }

    double GetAbsDeterminant(){
      if (!has_trafo_base) CalcTrafoBase();
      return absdet;
    }

    double GetVolume(){
      if (!has_trafo_base) CalcTrafoBase();
      return vol;
    }

    SMatrixCL<3,4> & GetDShape(){
      if (!has_trafo_base) CalcTrafoBase();
      return G;
    }

    SMatrixCL<4,4> & GetGramShape(){
      if (!has_Gram){
        if (!has_trafo_base) CalcTrafoBase();
        GTG = GramMatrix(G);
        has_Gram = true;
      } 
      return GTG;
    }

    Quad3CL<> & GetBaseShapeAsQuad3CL(int i) {
      if (!oninterface)
        return p1fegfs.GetShapeAsQuad3(i);
      else
        return q3_baseshape[i];
    }

};

class GlobalConvDiffReacCoefficients{
  friend class LocalConvDiffReacCoefficients;
  private:
    typedef BndDataCL<Point3DCL> VelBndDataT;
    typedef P2EvalCL<SVectorCL<3>, const VelBndDataT, const VecDescCL> const_DiscVelSolCL;  
    double D_;
    double time;    
    const VelocityContainer& vel;
    instat_scalar_fun_ptr source;
    instat_scalar_fun_ptr mass;
  public:

    GlobalConvDiffReacCoefficients(const double D, const VelocityContainer& u, instat_scalar_fun_ptr c, instat_scalar_fun_ptr f, double atime)
      :  D_( D), time(atime), vel(u), source(f), mass(c){}

    double GetDiffusionCoef(){
    	return D_;
    }
	
	instat_scalar_fun_ptr GetSourceAsFuntionPointer(){return source;}
	instat_scalar_fun_ptr GetMassAsFuntionPointer(){return mass;}
};

class LocalConvDiffReacCoefficients{
  private:
    GlobalConvDiffReacCoefficients& gcdcoefs;  
    Quad5CL<> q5_source;
    Quad3CL<> q3_source;
    Quad5CL<> q5_mass;
    Quad3CL<> q3_mass;
    Quad3CL<Point3DCL> *q3_velocity;
    LocalP2CL<Point3DCL> *lp2_velocity;  
  public:
    LocalConvDiffReacCoefficients(GlobalConvDiffReacCoefficients& agcdcoefs, TetraCL& tet):gcdcoefs(agcdcoefs), 
        q5_source(tet,gcdcoefs.source,gcdcoefs.time), q3_source(tet,gcdcoefs.source,gcdcoefs.time),
        q5_mass(tet,gcdcoefs.mass,gcdcoefs.time), q3_mass(tet,gcdcoefs.mass,gcdcoefs.time)
    {
      if (gcdcoefs.vel.hasVelocityAsP2()){
        q3_velocity = new Quad3CL<Point3DCL>(tet, gcdcoefs.vel.GetVelocityAsP2());
        lp2_velocity = new LocalP2CL<Point3DCL>(tet, gcdcoefs.vel.GetVelocityAsP2());
      }
      else
      {
        q3_velocity = new Quad3CL<Point3DCL>(tet, gcdcoefs.vel.GetVelocityAsFunctionPointer(),gcdcoefs.time);
        lp2_velocity = new LocalP2CL<Point3DCL>(tet, gcdcoefs.vel.GetVelocityAsFunctionPointer(),gcdcoefs.time);
      }
    }
    
    ~LocalConvDiffReacCoefficients(){
      delete q3_velocity;
      delete lp2_velocity;
    }

    double GetDiffusionCoef(){
      return gcdcoefs.GetDiffusionCoef();
    }    
	
	instat_scalar_fun_ptr GetSourceAsFuntionPointer(){
	  return gcdcoefs.GetSourceAsFuntionPointer();
	}
	double GetTime(){return gcdcoefs.time;}
	Quad5CL<>& GetSourceAsQuad5(){return q5_source;}
	Quad3CL<>& GetSourceAsQuad3(){return q3_source;}
	Quad5CL<>& GetReactionAsQuad5(){return q5_mass;}
	Quad3CL<>& GetReactionAsQuad3(){return q3_mass;}
	Quad3CL<Point3DCL>& GetVelocityAsQuad3(){return *q3_velocity;}
	LocalP2CL<Point3DCL>& GetVelocityAsLocalP2(){return *lp2_velocity;}

};





typedef double Elmat4x4[4][4];
typedef double Elvec4[4];

class ConvDiffElementMatrices{
  public:
    Elmat4x4 A;
    Elmat4x4 M;
    Elmat4x4 Mr;
    Elmat4x4 C;

    void ResetUnsigned(){
      std::memset( A,0, 4*4*sizeof(double));
      std::memset( M,0, 4*4*sizeof(double));
      std::memset( Mr,0, 4*4*sizeof(double));
      std::memset( C,0, 4*4*sizeof(double));
    }

    void ResetAll(){
      ResetUnsigned();
    }

    ConvDiffElementMatrices(){
      ResetAll();
    }
};



} // end of namespace DROPS

#endif
