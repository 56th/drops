/* 
 * File:   singletonmap.h
 * Author: berger
 *
 * Created on July 25, 2012, 10:15 AM
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
};

template <class T>
class MapRegisterCL
{
  public:
    MapRegisterCL(std::string name, T t){
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
T& SingletonMapCL<T>::operator[](std::string s){
    if (this->find(s) == this->end())
        throw DROPSErrCL(std::string("SingletonMapCL::operator[]: function with the name \"") + s + std::string("\" not found in container!"));
    return this->find(s)->second;
}





}
#endif	/* SINGLETONMAP_H */

