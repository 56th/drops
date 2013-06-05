/// \file staticinit.h
/// \brief Nifty counter idiom. Reference-counted initialization and destruction of static data (e.g. RefTetraPartitionCL).
/// \author LNM RWTH Aachen: Joerg Grande; SC RWTH Aachen:

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
 * Copyright 2012 LNM, RWTH Aachen, Germany
*/

#ifndef DROPS_STATICINIT_H
#define DROPS_STATICINIT_H

namespace DROPS {

/// \brief Initialization and destruction of static data for a class T (e.g. RefTetraPartitionCL).
/// Usage: Define a variable StaticInitializerCL<T> in an unnamed namespace in the .cpp-file for the class behind all members that need initialization.
template<typename T>
class StaticInitializerCL
{
  public:
    StaticInitializerCL () {
            T::StaticInit();
    }
    ~StaticInitializerCL () {
            T::StaticDestruct();
    }
};

/// \brief Reference-counted initialization and destruction of static data for a class T.
/// This is for uses, where the order accross translation units matters.
/// Usage: Define a variable RefCountedStaticInitializerCL<T> in an unnamed namespace in the header,
/// where T is defined. Define the member RefCountedStaticInitializerCL<T>::counter_= 0 in the corresponding .cpp-file.
///
/// Caveat: Static class-type objects, that have a constructor, defeat this scheme, as the constructor will be run after the call to StaticInit.
///
/// The static member T::StaticInit() is called once at programm startup as early as neccessary. The static member T::StaticDestruct() is called once at programm exit as late as neccessary.
template<typename T>
class RefCountedStaticInitializerCL
{
  private:
     static size_t counter_;

  public:
    RefCountedStaticInitializerCL () {
        if (counter_++ == 0)
            T::StaticInit();
    }
    ~RefCountedStaticInitializerCL () {
        if (--counter_ == 0)
            T::StaticDestruct();
    }
};

} // end of namespace DROPS
#endif
