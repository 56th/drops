/// \file progressaccu.cpp
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

#include "misc/progressaccu.h"
#include "misc/scopetimer.h"

#ifndef DROPS_WIN
  #include <unistd.h>
#else
  #include <io.h>
#endif

namespace DROPS
{


void ProgressBarTetraAccumulatorCL::begin_accumulation ()
{
  if (active)
  {
    // std::cout  //<< "\n"
    //            // << std::setw(20) 
    //            << name << ": ";
    tetprog = new Uint(0);
    prog = new int(-1);
  }
}

void ProgressBarTetraAccumulatorCL::finalize_accumulation ()
{
  if (active)
  {
    if (isterm)
    {
      std::cout                     
        << "\r"
        // << std::setw(20) 
        << name << ": " << std::setw(3) << 100 << "% [ "
        << std::setw(8) << ntet << " / " << std::setw(8) << ntet << " ]"
        << " - [ finished ]" << std::endl;
    }
    else
      std::cout << "] - [ finished ]" << std::endl;
    delete tetprog;
    delete prog;
  }
}

void ProgressBarTetraAccumulatorCL::visit (const TetraCL& )
{
  ScopeTimerCL scope("ProgressBarTetraAccumulator");
  if (active)
  {
    bool doout = false; //default: keine ausgabe
    int newprog = 0; //local copy of prog
    Uint mytetprog = 0; //local copy of tetprog

#pragma omp critical(tetprog)
    mytetprog = (*tetprog)++;

    newprog = (Uint) ( (100.0 * mytetprog) / (double)ntet);

#pragma omp critical(prog)
    if (newprog > *prog){
      *prog=newprog;

      if (isterm || newprog%5==0)
        doout = true;
    }

    if (doout) //output of (possibly outdated [who cares?!] ) local copies...
    {
#pragma omp critical(doout)
      {
        if (isterm){
          std::cout // << "\r                          "
            << "\r"
            // << std::setw(20) 
            << name << ": " << std::setw(3) << newprog << "% [ "
            << std::setw(8) << mytetprog << " / " << std::setw(8) << ntet << " ]" << " - [";
          for (int i = 0; i < newprog/10; ++i)
            std::cout << "#";
          for (int i = newprog/10; i < 10; ++i)
            std::cout << " ";
          std::cout << "]" << std::flush;
        }
        else //poor mans output (compatible to output files)
        {
          if (newprog==0) 
            std::cout << name << ": ["<< std::flush;
          else
            std::cout << "#" << std::flush;
        }
      }
    }
  }
}

ProgressBarTetraAccumulatorCL::ProgressBarTetraAccumulatorCL(const MultiGridCL& MG, const std::string aname, int lvl)
  :name(aname)
{
#ifndef DROPS_WIN
    isterm = isatty(fileno(stdout));
#else
    isterm = _isatty( _fileno(stdout));
#endif
    if (active)
    {
      ntet=0;
      DROPS_FOR_TRIANG_CONST_TETRA( MG, lvl, It) ++ntet;
    }
}

ProgressBarTetraAccumulatorCL::ProgressBarTetraAccumulatorCL(int tets, const std::string aname)
  :ntet(tets), name(aname) {  }


void MaybeAddProgressBar ( const MultiGridCL & MG, const std::string name, TetraAccumulatorTupleCL & accutup, int lvl){
  if (ProgressBarTetraAccumulatorCL::active)
    accutup.push_back_acquire( new ProgressBarTetraAccumulatorCL(MG, name, lvl));
}

void MaybeAddMLProgressbar( const MultiGridCL & MG, const std::string name, MLTetraAccumulatorTupleCL & mlaccutup, int lastlevel){
  if (ProgressBarTetraAccumulatorCL::active)
  {
    int dist = std::distance(mlaccutup.begin(),mlaccutup.end());
    MLTetraAccumulatorTupleCL::iterator accutup = mlaccutup.begin();
    for (int lvl=lastlevel+1-dist; accutup != mlaccutup.end(); ++accutup, lvl++)
    {
      accutup->push_back_acquire( new ProgressBarTetraAccumulatorCL(MG, name, lvl));
    }
  }
}

bool ProgressBarTetraAccumulatorCL::active = false;

}
