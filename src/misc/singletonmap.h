/// \file singletonmap.h
/// \brief Template class for a SingletonMap string-->T.
/// \author LNM RWTH Aachen: Jens Berger
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
 * Copyright 2012, 2013 LNM RWTH Aachen, Germany
*/

#ifndef SINGLETONMAP_H
#define	SINGLETONMAP_H

#include "misc/utils.h"

namespace DROPS
{

template<class T>
class SingletonMapCL : public std::map<std::string, T>
{
  private:
    SingletonMapCL() {}                         // von aussen keine Instanzen erzeugbar
    SingletonMapCL(const SingletonMapCL&) : std::map<std::string, T>() { }  // nicht kopierbar
    ~SingletonMapCL() {}
  public:
    static SingletonMapCL& getInstance();
    T& operator[](std::string s);
    void PrintAll() {
        std::cout << " map contains : \n";
        for ( typename SingletonMapCL<T>::const_iterator iter = this->begin();
              iter != this->end(); ++iter )
            std::cout << iter->first << std::endl;

    }
};


template <class T>
class MapRegisterCL
{
  public:
    MapRegisterCL(std::string name, T t) {
        SingletonMapCL<T>::getInstance().insert(std::make_pair(name, t));
    }
};


template<class T>
SingletonMapCL<T>& SingletonMapCL<T>::getInstance()
{
    static SingletonMapCL instance;
    return instance;
}

template<class T>
T& SingletonMapCL<T>::operator[](std::string s)
{
    if (this->find(s) == this->end())
        throw DROPSErrCL(std::string("SingletonMapCL::operator[]: function with the name \"") + s + std::string("\" not found in container.\n"));
    return this->find(s)->second;
}

}
#endif	/* SINGLETONMAP_H */

