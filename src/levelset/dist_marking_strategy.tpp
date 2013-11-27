/// \file dist_marking_strategy.tpp
/// \brief 
/// \author LNM RWTH Aachen: Matthias Kirchhart

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
 * Copyright 2011 LNM/SC RWTH Aachen, Germany
*/

namespace DROPS
{

template <class DistFctT>
void DistMarkingStrategyCL<DistFctT>::visit( const TetraCL& t )
{
    double d = std::numeric_limits<double>::max();
    int num_pos = 0;
    for ( Uint j = 0; j < 4; ++j )
    {
        const double dist = GetValue( *t.GetVertex(j) );
        if ( dist >= 0 ) ++num_pos;
        d = std::min( d, std::abs(dist) );
    }
    
    for ( Uint j = 0; j < 6; ++j )
    {
        const double dist = GetValue( *t.GetEdge(j) );
        if ( dist >= 0 ) ++num_pos;
        d = std::min( d, std::abs(dist) );
    }

    d = std::min( d, std::abs(GetValue(t)) );

    const bool vzw = num_pos!=0 && num_pos!=10; // change of sign
    const Uint l = t.GetLevel();

    // In the shell:      level should be f_level_.
    // Outside the shell: level should be c_level_.
    const Uint soll_level = ( d <= width_ || vzw ) ? f_level_ : c_level_;

    if ( l !=  soll_level || (l == soll_level && ! t.IsRegular()) )
    {
        // tetra will be marked for refinement/removement
        if ( l <= soll_level ) t.SetRegRefMark();
        else t.SetRemoveMark();
        modified_ = true;
    }
}

void DistMarkingStrategyCL<instat_scalar_fun_ptr>::visit( const TetraCL& t )
{
    double d = std::numeric_limits<double>::max();
    int num_pos = 0;
    for ( Uint j = 0; j < 4; ++j )
    {
        const double dist = GetValue( *t.GetVertex(j) );
        if ( dist >= 0 ) ++num_pos;
        d = std::min( d, std::abs(dist) );
    }
    
    for ( Uint j = 0; j < 6; ++j )
    {
        const double dist = GetValue( *t.GetEdge(j) );
        if ( dist >= 0 ) ++num_pos;
        d = std::min( d, std::abs(dist) );
    }

    d = std::min( d, std::abs(GetValue(t)) );

    const bool vzw = num_pos!=0 && num_pos!=10; // change of sign
    const Uint l = t.GetLevel();

    // In the shell:      level should be f_level_.
    // Outside the shell: level should be c_level_.
    const Uint soll_level = ( d <= width_ || vzw ) ? f_level_ : c_level_;

    if ( l !=  soll_level || (l == soll_level && ! t.IsRegular()) )
    {
        // tetra will be marked for refinement/removement
        if ( l <= soll_level ) t.SetRegRefMark();
        else t.SetRemoveMark();
        modified_ = true;
    }
}

}

