/// \file integrTime.h
/// \brief classes that perform time-integration steps
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Marcus Soemers, Volker Reichelt, Liang Zhang; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_POI_INTEGRTIME_H
#define DROPS_POI_INTEGRTIME_H

#include "num/poissonsolverfactory.h"
#include "misc/problem.h"
#include "poisson/poisson.h"

/// \todo FracStepScheme fuer instat. Poisson


namespace DROPS
{


/*****************************************************************************
*   for solving the instat. Poisson equation of type PoissonT with a
*   1-step-theta-scheme: theta=1   -> impl. Euler (BDF1, backward Euler)
*                        theta=1/2 -> Crank-Nicholson (Trapezregel)
*
*   Inner stat. Poisson-type problems are solved with a SolverT-solver.
*   The matrices A, M and the rhs b of the Poisson class have
*   to be set properly! After construction, SetTimeStep has to be called once.
*   Then every DoStep performs one step in time. Changing time steps require
*   further calls to SetTimeStep.
******************************************************************************/

template <class PoissonT, class SolverT>
class InstatPoissonThetaSchemeCL
{
  private:
    PoissonT&   _Poisson;
    SolverT&    _solver;
    const ParamCL& _param;
    
    VecDescCL *_b, *_old_b;             // rhs
    VecDescCL *_cplA, *_old_cplA;       // couplings with poisson matrix A
    VecDescCL *_cplM, *_old_cplM;       // couplings with mass matrix M
    VecDescCL *_cplU;                   // couplings with convection matrix U
    VectorCL  _rhs;
    VectorCL  _Zeta;                    // one help sequence for generalized theta scheme
    MLMatrixCL  _Lmat;                  // M + theta*dt*nu*A  = linear part

    double _theta, _dt;
    bool   _Convection;
    bool   _supg;
    bool   _ale;
    bool   _firstStep;                  


  public:
    InstatPoissonThetaSchemeCL( PoissonT& Poisson, SolverT& solver, const ParamCL& param)
    : _Poisson( Poisson), _solver( solver), _param(param),
      _b( &Poisson.b), _old_b( new VecDescCL),
      _cplA( new VecDescCL), _old_cplA( new VecDescCL),
      _cplM( new VecDescCL), _old_cplM( new VecDescCL),
      _cplU( new VecDescCL),
      _rhs( Poisson.b.RowIdx->NumUnknowns()), 
      _Zeta( Poisson.b.RowIdx->NumUnknowns()),
      _theta( _param.get<double>("Time.Theta")),
      _Convection(_param.get<int>("Poisson.Coeff.withConvection")), _supg(_param.get<int>("Stabilization.SUPG")),
      _ale(_param.get<int>("ALE.wavy")),_firstStep(true)
    {
      _old_b->SetIdx( _b->RowIdx);
      _cplA->SetIdx( _b->RowIdx); _old_cplA->SetIdx( _b->RowIdx);
      _cplM->SetIdx( _b->RowIdx); _old_cplM->SetIdx( _b->RowIdx);
      _Poisson.SetupInstatRhs( *_old_cplA, *_old_cplM, _Poisson.x.t, *_old_b, _Poisson.x.t);
      if (_Convection)
      {
        _cplU->SetIdx( _b->RowIdx);
        _Poisson.SetupConvection( _Poisson.U, *_cplU, _Poisson.x.t);
      }
    }

    ~InstatPoissonThetaSchemeCL()
    {
      if (_old_b == &_Poisson.b)
        delete _b;
      else
        delete _old_b;
      delete _cplA; delete _old_cplA;
      delete _cplM; delete _old_cplM;
      delete _cplU;
    }

    double GetTheta()    const { return _theta; }
    double GetTimeStep() const { return _dt; }

    void SetTimeStep( double dt)
    {
      _dt= dt;
      if (!_Convection)
        _Lmat.LinComb( 1., _Poisson.M.Data, _theta*_dt, _Poisson.A.Data);
    }

    void DoStep( VecDescCL& v);
    void StadScheme( VecDescCL& v);        //standard theta scheme with mass matrix not dependant on time
    void GeneralScheme   ( VecDescCL& v);  //generalized scheme with mass matirx depends on time
};




//=================================
//     template definitions
//=================================


template <class PoissonT, class SolverT>
void InstatPoissonThetaSchemeCL<PoissonT,SolverT>::DoStep( VecDescCL& v)
{
    if(_supg||_ale)
    {   
        GeneralScheme(v);       
    }
    else 
        StadScheme(v);
}

template <class PoissonT, class SolverT>
void InstatPoissonThetaSchemeCL<PoissonT,SolverT>::StadScheme( VecDescCL& v)
{
  _Poisson.x.t+= _dt;

  _Poisson.SetupInstatRhs( *_cplA, *_cplM, _Poisson.x.t, *_b, _Poisson.x.t);
  
  _rhs = _Poisson.A.Data * v.Data;
  _rhs*= -_dt*(1.0-_theta); 

  _rhs+=  _dt*(1.0-_theta)*(_old_b->Data)
         +_dt*(1.0-_theta)* _old_cplA->Data;
         
  _rhs +=_Poisson.M.Data*v.Data
          -_old_cplM->Data
          +_dt*_theta*_b->Data
          + _dt*_theta*_cplA->Data
          + _cplM->Data;
  

  if (_Convection)
  {
      _rhs+= (_dt*(1.0-_theta)) * (_cplU->Data - _Poisson.U.Data * v.Data );
      _Poisson.SetupConvection( _Poisson.U, *_cplU, _Poisson.x.t);
      _rhs+= (_dt*_theta) * _cplU->Data;
      MLMatrixCL AU;
      AU.LinComb( 1, _Poisson.A.Data, 1, _Poisson.U.Data);
      _Lmat.LinComb( 1, _Poisson.M.Data, _dt*_theta, AU);
  }
  _solver.Solve( _Lmat, v.Data, _rhs, v.RowIdx->GetEx());

  std::swap( _b, _old_b);
  std::swap( _cplA, _old_cplA);
  std::swap( _cplM, _old_cplM);    
    
    
}

template <class PoissonT, class SolverT>
///implemented correctly only in the case Dirichlet boundary conditions are constant in time;
void InstatPoissonThetaSchemeCL<PoissonT,SolverT>::GeneralScheme( VecDescCL& v)
{
  _Poisson.x.t+= _dt;
  if(_firstStep)
  {
      VectorCL  _innerb( _Poisson.b.RowIdx->NumUnknowns());
      _innerb = _old_b->Data - _Poisson.A.Data * v.Data + _old_cplA->Data;
      if (_Convection)
      {
                _innerb+= _cplU->Data - _Poisson.U.Data * v.Data;
      } 
      SSORPcCL pc(1.0);
      GMResSolverCL<SSORPcCL> solver(pc, 50, 1000, 1.0e-12);
      solver.Solve( _Poisson.M.Data , _Zeta, _innerb, _Poisson.idx.GetEx());
      std::cout << " o Solved system with:\n"
                << "   - iterations    " << solver.GetIter()  << '\n'
                << "   - residuum      " << solver.GetResid() << '\n';
        
     _firstStep = false;  
  }
 
  VectorCL _old_sol( _Poisson.b.RowIdx->NumUnknowns());
  _old_sol = v.Data;
  //update mass matrix, stiffness matrix for SUPG and ALE, since the test functions change in different steps
  _Poisson.SetupInstatSystem( _Poisson.A, _Poisson.M, _Poisson.x.t);
  
  _Poisson.SetupInstatRhs( *_cplA, *_cplM, _Poisson.x.t, *_b, _Poisson.x.t);
  

  //_Poisson.SetupInstatRhs( *_old_cplA, *_old_cplM, _Poisson.x.t-_dt, *_old_b, _Poisson.x.t);
  

  _rhs = (1. - _theta) * _dt * (_Poisson.M.Data * _Zeta);
  _rhs += _Poisson.M.Data * v.Data
          //-_old_cplM->Data
          +_dt*_theta*_b->Data
          + _dt*_theta*_cplA->Data;
          //+ _cplM->Data;
  

  if (_Convection)
  {
      _Poisson.SetupConvection( _Poisson.U, *_cplU, _Poisson.x.t);
      _rhs+= (_dt*_theta) * _cplU->Data;
      MLMatrixCL AU;
      AU.LinComb( 1, _Poisson.A.Data, 1, _Poisson.U.Data);
      _Lmat.LinComb( 1, _Poisson.M.Data, _dt*_theta, AU);
  }
  _solver.Solve( _Lmat, v.Data, _rhs, v.RowIdx->GetEx());
  //Update zeta sequence
  VectorCL _old_Zeta( _Poisson.b.RowIdx->NumUnknowns());
  _old_Zeta = _Zeta;
  _Zeta = 1./(_dt * _theta) * (_Poisson.x.Data - _old_sol) - (1. - _theta)/_theta * _old_Zeta;
  std::swap( _b, _old_b);
  std::swap( _cplA, _old_cplA);
  std::swap( _cplM, _old_cplM);    
    
}

}    // end of namespace DROPS

#endif
