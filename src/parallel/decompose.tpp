/// \file decompose.tpp
/// \brief determining a partitioning of a distributed hierarchy of triangulations.
/// \author LNM RWTH Aachen: Patrick Esser, Sven Gross; SC RWTH Aachen: Oliver Fortmeier

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
*/

namespace DROPS{

bool LbIteratorCL::IsInLbSet(const TetraCL& t) const
/// Test whether a tetrahedra represents a multi-node of the LoadBalanceSet
{
    if (t.HasGhost() )
        return false;
    if (t.IsUnrefined() )
    {
        if (t.GetLevel()==0)
            return true;
        else
            return false;   // there is a parent, that represents this tetra
    }

    // check for children in this triangulation level, if there is a child, than this tetra is in the LoadBalanceSet
    for (TetraCL::const_ChildPIterator ch(t.GetChildBegin()), chend(t.GetChildEnd()); ch!=chend; ++ch)
        if ((*ch)->IsInTriang(TriLevel_) )
            return true;

    // If non of the above cases happens, this is not in the LoadBalanceSet
    return false;
}


LbIteratorCL& LbIteratorCL::operator ++ ()
/** Increase the position, until we reach an element of the LoadBalanceSet or the
    end of all tetraeders in the triangulation level _TriLevel
*/
{
    do
    {
        ++pos_;                                         // increase pos

        while (pos_ == mg_->GetTetrasEnd( Level_) )     // if we reaches the end of a level
        {
            if (Level_ < TriLevel_)
                pos_= mg_->GetTetrasBegin( ++Level_);   // goto the next level, if in the right triangulation
            else
                return *this;                           // return end
        }
    } while (!IsInLbSet() );        // until reached an element in the LoadBalanceSet
    return *this;                   // return the found element
}


LbIteratorCL LbIteratorCL::operator ++ (int)
{
    LbIteratorCL tmp(*this);
    ++(*this);
    return tmp;
}
}
