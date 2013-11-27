/// \file marking_strategy.cpp
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
#include "levelset/marking_strategy.h"

namespace DROPS
{

void StrategyCombinerCL::visit( const TetraCL& t )
{
    bool remove = false;

    for ( Uint i = 0; i < strategies_.size(); ++i )
    {
        // Every input strategy expects Refmark to be on NoRef
        // on input. By default, this is true for any leaf.
	t.SetNoRefMark();
	
        MarkingStrategyCL *s = strategies_[ i ];
        s->visit( t );

        Uint NewMark = t.GetRefMark();
        switch ( NewMark )
        {
        case NoRefMarkC:
             // I don't care case. Don't do anything.
             break;
        case RegRefMarkC:
             // Regular refinement always wins. Stop here.
             modified_ = true;
             return;
        case RemoveMarkC:
             // Continue loop. Some other strategy might still want to refine.
             remove = true;
             break;
        }
    }

    if ( remove )
    {
        t.SetRemoveMark();
        modified_ = true;
    }
}

}

