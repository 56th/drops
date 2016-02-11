/// \file integrTime2phase.h
/// \brief classes that perform time-integration steps for two phase flows with given levelset function
/// \author LNM RWTH Aachen: Thomas Ludescher

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

#ifndef DROPS_STO_INTEGRTIME2PHASE_H
#define DROPS_STO_INTEGRTIME2PHASE_H

#include "stokes/stokes.h"
#include "stokes/instatstokes2phase.h"
//#include "stokes/integrTime.h"
#ifdef _PAR
# include "misc/problem.h"
#endif
#include <typeinfo>

//extern DROPS::ParamCL P;
typedef std::numeric_limits<double> dbllimit;

namespace DROPS
{

template< class StokesT, class SolverT, class LsetT>
class TimeDiscStokes2phaseCL
{
  protected:
    StokesT& _Stokes;
    SolverT& _solver;
    LsetT& _lset;

    VelVecDescCL *_b, *_old_b;        // rhs + couplings with poisson matrix A
    VelVecDescCL *_cplM, *_old_cplM;  // couplings with mass matrix M
    VectorCL      _rhs;
    MLMatrixCL    _mat;               // M + theta*dt*A

    double _theta, _dt;

  public:
    TimeDiscStokes2phaseCL( StokesT& Stokes, SolverT& solver, LsetT& lset, double theta= 0.5)
        : _Stokes( Stokes), _solver( solver), _lset(lset), _b( &Stokes.b), _old_b( new VelVecDescCL),
        _cplM( new VelVecDescCL), _old_cplM( new VelVecDescCL), _rhs( Stokes.b.RowIdx->NumUnknowns()),
        _theta( theta)
        {
            _old_b->SetIdx( _b->RowIdx); _cplM->SetIdx( _b->RowIdx); _old_cplM->SetIdx( _b->RowIdx);
            //_Stokes.SetupInstatRhs( _old_b, &_Stokes.c, _old_cplM, _Stokes.v.t, _old_b, _Stokes.v.t);

            // only _old_cplM is required here, if theta!=1 _old_b is needed as well

            //_Stokes.SetupSystem1( &_Stokes.A, &_Stokes.M, _old_b, _old_b, _old_cplM, _lset, _Stokes.v.t);

        };

    virtual ~TimeDiscStokes2phaseCL()
    {
        if (_old_b == &_Stokes.b)
            delete _b;
        else
            delete _old_b;
        delete _cplM; delete _old_cplM;
    }

    double GetTheta()    const { return _theta; }
    double GetTime()     const { return _Stokes.v.t; }
    double GetTimeStep() const { return _dt; }

    /// \brief Get constant reference on solver
    const SolverT& GetSolver() const { return _solver; }
    /// \brief Get reference on solver
    SolverT& GetSolver()       { return _solver; }

    virtual void SetTimeStep( double dt)= 0;

    virtual void SetTimeStep( double dt, double theta)= 0;

    virtual void DoStep( VectorCL& v, VectorCL& p)= 0;

    MLMatrixCL* GetUpperLeftBlock() { return &_mat; }
    const MLMatrixCL* GetUpperLeftBlock() const { return &_mat; }
};

template < class StokesT, class SolverT, class LsetT>
class InstatStokes2phaseThetaSchemeCL : public TimeDiscStokes2phaseCL< StokesT, SolverT, LsetT>
//**************************************************************************
//  for solving the instationary Stokes equation of type StokesT with a
//  1-step-theta-scheme: theta=1   -> impl. Euler (BDF1, backward Euler)
//                       theta=1/2 -> Crank-Nicholson (Trapezregel)
//
//  Inner stationary Stokes-type problems are solved with a SolverT-solver.
//  The matrices A, B, M and the rhs b, c of the Stokes class have to be set
//  properly! After construction, SetTimeStep has to be called once. Then
//  every DoStep performs one step in time. Changing time steps require
//  further calls to SetTimeStep.
//**************************************************************************
{
  private:
    typedef TimeDiscStokes2phaseCL< StokesT, SolverT, LsetT> base_;
    using base_:: _Stokes;
    using base_:: _solver;
    using base_::  _lset;

    using base_:: _b;        using base_::_old_b;        // rhs + couplings with poisson matrix A
    using base_:: _cplM;     using base_::_old_cplM;  // couplings with mass matrix M
    using base_:: _rhs;
    using base_:: _mat;               // M + theta*dt*A

    using base_:: _theta;    using base_:: _dt;


  public:
    InstatStokes2phaseThetaSchemeCL( StokesT& Stokes, SolverT& solver, LsetT& lset, double theta= 0.5 )
        :base_(Stokes, solver, lset, theta) {}

    ~InstatStokes2phaseThetaSchemeCL(){}

    void SetTimeStep( double dt)
    {
        _dt= dt;        
        _mat.LinComb( 1./_dt, _Stokes.M.Data, _theta, _Stokes.A.Data);        
    }

    void SetTimeStep( double dt, double theta)
    {        
        _theta= theta;
        SetTimeStep( dt );
    }

    void DoStep(VectorCL &v, VectorCL &p);

    void CommitStep( const VecDescCL &un );

};

template <class StokesT, class SolverT, class LsetT>
void InstatStokes2phaseThetaSchemeCL<StokesT,SolverT,LsetT>::DoStep( VectorCL& v, VectorCL& p)
{
    double oldTime = _Stokes.v.t;
    _Stokes.v.t+= _dt;

    //update XFEM numberings ... only needed for moving interface
    if ( _Stokes.UsesXFEM() )
    {
        _Stokes.UpdateXNumbering( &_Stokes.pr_idx, _lset);
        _Stokes.UpdatePressure( &_Stokes.p);
    }

    // The MatrixBuilderCL's method of determining when to reuse the pattern
    // is not safe for P1X-elements.
    _Stokes.ClearMat();
    _Stokes.SetIdx();
    _mat.clear();


    //MLIdxDescCL *vidx= &_Stokes.vel_idx;
    IdxDescCL *vidx =  &(_Stokes.vel_idx.GetFinest() );


    _rhs.resize( vidx->NumUnknowns() );
    _b->SetIdx( vidx );    
    _cplM->SetIdx( vidx );
    _old_cplM->SetIdx( vidx );        

    //get _old_cplM    
    _Stokes.SetupCplM( _old_cplM, _lset, oldTime );

    /*
     * for testing purposes only
    // very inefficient routine .... A not needed at all ...
    VelVecDescCL *tmpb = new VelVecDescCL( vidx );
    //tmpb->SetIdx( vidx );
    _Stokes.SetupSystem1( &_Stokes.A, &_Stokes.M, tmpb, tmpb, _cpl_tst, _lset, oldTime);
    delete tmpb;


    VectorCL v1 =  _cpl_tst->Data ;
    VectorCL v2 =  _old_cplM->Data ;
    VectorCL res( v1.size() );
    res = v1-v2;    

    double scalprod = norm(res);

    std::cout << "\n\ncplM test: " << scalprod << std::endl;
    std::cout << norm(v1) << std::endl;
    std::cout << norm(v2) << std::endl << std::endl;
    */


    // if system matrices do not change (i.e. stationary interface) quite inefficient
    // but cplM is needed if boundary conditions are time dependent; cplA is inside b
    _Stokes.SetupSystem1( &_Stokes.A, &_Stokes.M, _b, _b, _cplM, _lset, _Stokes.v.t);
    _Stokes.SetupSystem2( &_Stokes.B, &_Stokes.C, &_Stokes.c, _lset, _Stokes.v.t );

    // if moving interface: uncomment the following lines        
    _Stokes.SetupPrMass ( &_Stokes.prM, _lset );
    _Stokes.SetupPrStiff( &_Stokes.prA, _lset, P.get<double>("Stokes.lambda",1.0) );


    // create matrix mat
    SetTimeStep(_dt);

    // compute surface tension integral for rhs    
    VelVecDescCL surften( &_Stokes.vel_idx );
    _lset.AccumulateBndIntegral( surften );
    _b->Data += surften.Data;

    _rhs = (1./_dt) * (_Stokes.M.Data*v + _cplM->Data - _old_cplM->Data ) + _theta* _b->Data
            + ( 1. - _theta ) * _old_b->Data;


    /*
    _rhs=  _Stokes.A.Data * v;
    _rhs*= (_theta-1.);
    _rhs+= (1./_dt)*(_Stokes.M.Data*v + _cplM->Data - _old_cplM->Data)
    //_rhs = (1./_dt)*(_Stokes.M.Data*v + _cplM->Data - _old_cplM->Data)
         +  _theta*_b->Data + (1.-_theta)*_old_b->Data;
    //_rhs = (1./_dt) * ( _Stokes.M.Data * v ) + _b->Data;
    //_rhs = _Stokes.b.Data + (1./_dt) * ( _Stokes.M.Data * v );
    */

    // solve stabilized version for two-phase flow
    std::cout << "dostep: solver.solve" << std:: endl;
    _solver.Solve( _mat, _Stokes.B.Data, _Stokes.C.Data, v, p, _rhs, _Stokes.c.Data, _Stokes.vel_idx.GetEx(), _Stokes.pr_idx.GetEx());
    //_solver.Solve( _mat, _Stokes.B.Data, v, p, _rhs, _Stokes.c.Data, _Stokes.vel_idx.GetEx(), _Stokes.pr_idx.GetEx());

    //std::swap( _b, _old_b);
    //std::swap( _cplM, _old_cplM);
}

// grid is now adapted to level set. old coupling terms need to be retrieved
template <class StokesT, class SolverT, class LsetT>
void InstatStokes2phaseThetaSchemeCL<StokesT,SolverT,LsetT>::CommitStep( const VecDescCL &un  )
{
    if( _theta == 1.0 ) return;

    MLIdxDescCL *vidx= &_Stokes.vel_idx;

    _Stokes.ClearMat();
    _Stokes.SetIdx();

    _old_b->SetIdx( vidx );

    _Stokes.SetupRhs1( _old_b, _lset, _Stokes.v.t );

    VelVecDescCL newCplA( vidx );
    _Stokes.SetupAdotU( &newCplA, un, _lset, _Stokes.v.t );

    _old_b->Data -= newCplA.Data;


    // get old surface tension force
    VelVecDescCL surften( vidx );
    _lset.AccumulateBndIntegral( surften );
    _old_b->Data += surften.Data;

    return;


    /*
     * for testing only
     *
    VelVecDescCL cplA( vidx );
    VectorCL oldRes( vidx->NumUnknowns() );
    VectorCL diff( oldRes.size() );

    VelVecDescCL *tmpb = new VelVecDescCL( vidx );
    _Stokes.SetupSystem1( &_Stokes.A, &_Stokes.M, tmpb, &cplA, tmpb, _lset, _Stokes.v.t);
    delete tmpb;

    oldRes = _Stokes.A.Data * un.Data - cplA.Data;
    VectorCL newRes = newCplA.Data;
    diff = newRes - oldRes;
    std::cout << "scalar prod" << std::endl;
    double scalprod = norm(diff);

    std::cout << "\n\ncplA test: " << scalprod << std::endl;

    std::cout.precision(dbllimit::digits10);
    std::cout << std::scientific;
    std::cout << norm(oldRes) << std::endl;
    std::cout << norm(newRes) << std::endl << std::endl;
    */


}

}    // end of namespace DROPS

#endif
