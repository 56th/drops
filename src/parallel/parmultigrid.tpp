/// \file parmultigrid.tpp
/// \brief handling of a parallel multigrid
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

namespace DROPS{

template<class IterT>
void ShowSimplex( IterT begin, IterT end, const DiST::GeomIdCL& gid, char *mesg, int proc= -1)
{
    if (proc!=ProcCL::MyRank() && proc!=-1)
        return;
    for (IterT it= begin; it!=end; ++it)
        if (it->GetGID()==gid)
        {
            std::cerr << "Show " << mesg << ":\n";
            it->DebugInfo( std::cerr);
        }
}

void ParMultiGridCL::ModifyBegin()
{
    if (transfer_) // modify module inside transfer module is just ignored
        return;
    Assert( !modify_, DROPSErrCL("ParMultiGridCL::ModifyBegin: already called ModifyBegin"), DebugParallelC);
    modify_= new DiST::ModifyCL(*mg_, false); modify_->Init();
}

void ParMultiGridCL::ModifyEnd()
{
    if (transfer_)  // modify module inside transfer module is just ignored
        return;
    Assert( modify_ , DROPSErrCL("ParMultiGridCL::ModifyEnd: not in Modify mode"), DebugParallelC);
    modify_->Finalize();
    delete modify_; modify_=0;
}

/// \brief Change the priority of a simplex
template<>
  void ParMultiGridCL::PrioChange<TetraCL>(TetraCL* const Tp, Priority Prio);


/// \brief Change the priority of a simplex
template<class SimplexT>
  void ParMultiGridCL::PrioChange(SimplexT* const Tp, Priority Prio)
{
    Assert( modify_, DROPSErrCL("ParMultiGridCL::PrioChange: There must be an active Transfer or Modify module to run this procedure"), DebugParallelC);
    modify_->ChangePrio( *Tp, Prio);
}

template<>
  void ParMultiGridCL::Keep<TetraCL>(TetraCL* const Tp); // not defined, as Keep() should not be called for a tetra


template<class SimplexT>
  void ParMultiGridCL::Keep(SimplexT* const Tp)
{
    Assert( modify_, DROPSErrCL("ParMultiGridCL::PrioChange: There must be an active Transfer or Modify module to run this procedure"), DebugParallelC);
    modify_->Keep( *Tp);
}

/// \brief Delete simplex from distributed multigrid
template<class SimplexT>
  void ParMultiGridCL::Delete(SimplexT* const Tp)
{
    Assert( modify_, DROPSErrCL("ParMultiGridCL::Delete: There must be an active Transfer or Modify module to run this procedure"), DebugParallelC);
    modify_->Delete( *Tp);
}

} // end of namespace DROPS
