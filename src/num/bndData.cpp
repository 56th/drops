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
#include "misc/params.h"
#include "misc/singletonmap.h"
#include "geom/multigrid.h"
#include "num/discretize.h"

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
      case Slip0BC:      os << "hom. slip BC\n"; break;
      case SlipBC:       os << "inhom. slip BC\n"; break;
      case NoBC:         os << "no boundary\n"; break;
      case UndefinedBC_: os << "WARNING! unknown BC from ReadMeshBuilderCL\n"; break;
      default:           os << "WARNING! unknown BC\n";
    }
}

BndCondT string_to_BndCondT (std::string s)
{
    const char* names[]=    { "Dir0BC", "DirBC", "Per1BC", "Per2BC", "Nat0BC", "NatBC", "Slip0BC", "SlipBC", "SymmBC", "OutflowBC", "WallBC", "NoBC", "UndefinedBC_", "MaxBC_" };
    const BndCondT types[]= {  Dir0BC,   DirBC,   Per1BC,   Per2BC,   Nat0BC,   NatBC,   Slip0BC,   SlipBC,   SymmBC,   OutflowBC,   WallBC,   NoBC,   UndefinedBC_,   MaxBC_  };
    const size_t num_types= sizeof( names)/sizeof( const char*);

    Uint i;
    for (i= 0; i < num_types; ++i)
        if (s == std::string( names[i]))
            break;
    if (i == num_types)
        throw DROPSErrCL( "string_to_enum: The string '" + s + "' does not matrch a boundary condition.\n");
    return types[i];
}


BndCondCL::BndCondCL (BndIdxT numbndseg, const BndCondT* bc)
{
    BndCond_.resize( numbndseg);
    for (Uint i=0; i<numbndseg; ++i)
        BndCond_[i]= bc ? bc[i] : Nat0BC;
}

void assignZeroFunc( instat_scalar_fun_ptr& f)
{
    f= SingletonMapCL<instat_scalar_fun_ptr>::getInstance()["Zero"];
}

void assignZeroFunc( instat_vector_fun_ptr& f)
{
    f= SingletonMapCL<instat_vector_fun_ptr>::getInstance()["ZeroVel"];
}


/// \brief Read an json-array with 1 or 2 members. The first is the BndCondT, the second possibly is a function-name.
/// Implementation detail of read_BndData.
template <class BndValFunT>
  void
  read_BndData1 (const ParamCL& P, BndCondT& type, BndValFunT& fun)
{
    ParamCL::ptree_const_iterator_type it= P.begin();
    type= string_to_BndCondT( it->second.get_value<std::string>());
    if (type != DirBC && type != NatBC) {
        assignZeroFunc( fun);
        return;
    }
    if (++it == P.end())
        throw DROPSErrCL( "read_BndData1: DirBC and NatBC must specify a function.");
    fun= SingletonMapCL<BndValFunT>::getInstance()[it->second.get_value<std::string>()];
}

template <class T>
  void
  read_BndData (BndDataCL<T>& bnddata, const MultiGridCL& mg, const ParamCL& P)
{
    typedef typename BndDataCL<T>::bnd_val_fun BndValFunT;
    const BoundaryCL& bnd= mg.GetBnd();
    const BndIdxT num_bnd= bnd.GetNumBndSeg();

    // Try to read and set the default.
    BndCondT   default_type= UndefinedBC_;
    BndValFunT default_fun=  0;
    const ParamCL::ptree_type* child= 0;
    try {
        child= &P.get_child( "Default");
    } catch (DROPSParamErrCL e) {}
    if (child != 0)
        try {
            read_BndData1 ( *child, default_type, default_fun);
        } catch (DROPSErrCL e) {
            std:: cerr << "read_BndData: While processing 'Default'...\n";
            throw e;
        }
    std::vector<BndCondT>   bnd_type( num_bnd, default_type);
    std::vector<BndValFunT> bnd_fun(  num_bnd, default_fun);


    // Read data for the boundary segments.
    for (ParamCL::ptree_const_iterator_type it= P.begin(), end= P.end(); it != end; ++it) {
        const std::string key= it->first;
        if (key == std::string( "Default") || key == std::string( "PeriodicMatching"))
            continue;

        BndIdxT i;
        std::istringstream iss( key);
        iss >> i;
        if (!iss) // As BndIdxT is unsigned, the 2nd test is redundant.
            throw DROPSErrCL( "read_BndData: Invalid boundary segment '" + key + "'.\n");

        BndCondT type=  UndefinedBC_;
        BndValFunT fun= 0;
        try {
            read_BndData1 ( it->second, type,fun);
        } catch (DROPSParamErrCL e) {
            std:: cerr << "read_BndData: While processing key '" << key << "'...\n";
            throw e;
        }
        bnd_type[i]= type;
        bnd_fun[i]=  fun;
    }

    // No UndefinedBC_ are allowed.
    std::vector<BndCondT>::const_iterator undef= std::find( bnd_type.begin(), bnd_type.end(), UndefinedBC_);
    if (undef != bnd_type.end()) {
        std::ostringstream os;
        os <<  "read_BndData: UndefinedBC_ set for boundary segment " << undef - bnd_type.begin() << ".\n",
        throw DROPSErrCL( os.str());
    }

    // warn if PeriodicMatching is set (will be ignored), this should be done in Mesh.PeriodicBnd
    child= 0;
    try {
        child= &P.get_child( "PeriodicMatching");
    } catch (DROPSParamErrCL e) {}
    if (child != 0)
        std:: cerr << "read_BndData: Warning: 'PeriodicMatching' is ignored and should be specified in the param section 'Mesh.PeriodicBnd'!\n";

    // Enter the data to bnddata
    std::vector<BndCondInfoCL> bnd_info( num_bnd);
    for (Uint i= 0; i < num_bnd; ++i)
        bnd_info[i]= bnd_type[i];
    bnddata.Init( bnd_info, bnd_fun);
}

/// Explicit instantiations for T = double and T = Point3DCL.
template void read_BndData<double>    (BndDataCL<double>& bnddata,    const MultiGridCL& mg, const ParamCL& P);
template void read_BndData<Point3DCL> (BndDataCL<Point3DCL>& bnddata, const MultiGridCL& mg, const ParamCL& P);

} //end of namespace DROPS
