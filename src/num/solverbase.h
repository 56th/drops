/// \file solverbase.h
/// \brief iterative solvers
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Helmut Jarausch, Volker Reichelt; SC RWTH Aachen:

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
 * Copyright 2012--2017 LNM/SC RWTH Aachen, Germany
*/

#ifndef DROPS_SOLVERBASE_H
#define DROPS_SOLVERBASE_H

namespace DROPS
{

/// What every iterative solver should have
class SolverBaseCL
{
  protected:
    int            maxiter_;
    mutable int    iter_= -1;
    double         tol_;
    mutable double res_= -1.;
    bool           rel_;

    mutable std::ostream* output_;

    SolverBaseCL (int maxiter, double tol, bool rel= false, std::ostream* output= 0)
        : maxiter_( maxiter), tol_( tol), rel_( rel), output_( output)  {}
    virtual ~SolverBaseCL() {}

  public:
    virtual void   SetTol     (double tol) { tol_= tol; }
    virtual void   SetMaxIter (int iter)   { maxiter_= iter; }
    virtual void   SetRelError(bool rel)   { rel_= rel; }

    virtual double GetTol     () const { return tol_; }
    virtual int    GetMaxIter () const { return maxiter_; }
    virtual double GetResid   () const { return res_; }
    virtual int    GetIter    () const { return iter_; }
    virtual bool   GetRelError() const { return rel_; }

    virtual void   SetOutput( std::ostream* os) { output_=os; }
};

} // end of namespace DROPS

#endif
