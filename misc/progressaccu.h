/// \file progressaccu.h
/// \brief accumulator which displays progress of accumulation 
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

#ifndef DROPS_PROGRESS_ACCU_H
#define DROPS_PROGRESS_ACCU_H

#include "num/accumulator.h"
#include "num/discretize.h"
#include "misc/problem.h"

namespace DROPS
{


class ProgressBarTetraAccumulatorCL : public TetraAccumulatorCL
{
    protected:   
    Uint * tetprog;
    int * prog;
    Uint ntet;
    const std::string name;
    bool isterm;
public: 

    static bool active;

    static void Activate(){ active = true; }
    static void Deactivate(){ active = false; }
    
    ProgressBarTetraAccumulatorCL (const MultiGridCL& MG, 
                                   const std::string aname = std::string("accumulator"), int lvl = -1);

    ProgressBarTetraAccumulatorCL (int tets, 
                                   const std::string aname = std::string("accumulator"));

    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    virtual void visit (const TetraCL&); 

    virtual TetraAccumulatorCL* clone (int /*tid*/) { return new ProgressBarTetraAccumulatorCL ( *this); }
    
};

void MaybeAddProgressBar ( const MultiGridCL & MG, const std::string name, TetraAccumulatorTupleCL & accutup, int lvl);
void MaybeAddMLProgressbar( const MultiGridCL & MG, const std::string name, MLTetraAccumulatorTupleCL & accutup, int lastlevel);

}

#endif
