/// \file bndData.cpp
/// \brief Classes for storing and handling boundary data.
/// \author LNM RWTH Aachen: Sven Gross, Joerg Grande

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
 * Copyright 2013 LNM, Germany
*/


#include "num/bndData.h"

namespace DROPS
{

void BndCondInfo (BndCondT bc, std::ostream& os)
{
    switch(bc)
    {
      case Dir0BC: /* WallBC has the same number */
                         os << "hom. Dirichlet BC / wall\n"; break;
      case DirBC:        os << "inhom. Dirichlet BC / inflow\n"; break;
      case Per1BC:       os << "periodic BC\n"; break;
      case Per2BC:       os << "periodic BC, correspondent\n"; break;
      case Nat0BC: /* OutflowBC has the same number */
                         os << "hom. Natural BC / outflow\n"; break;
      case NatBC:        os << "inhom. Natural BC\n"; break;
      case NoBC:         os << "no boundary\n"; break;
      case UndefinedBC_: os << "WARNING! unknown BC from ReadMeshBuilderCL\n"; break;
      default:           os << "WARNING! unknown BC\n";
    }
}

} //end of namespace DROPS
