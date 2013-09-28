/// \file dynamicload.cpp
/// \brief dynamic load functionality
/// \author LNM RWTH Aachen: Niklas Fischer, Christoph Lehrenfeld; SC RWTH Aachen:
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
 * Copyright 2013 LNM/SC RWTH Aachen, Germany
*/

#include "dynamicload.h"

#ifdef DROPS_WIN
#define SHAREDLIB_SUFFIX ".dll"
#else
#define SHAREDLIB_SUFFIX ".so"
#endif


namespace DROPS{

void read_file(const std::string& path);

#ifndef DROPS_WIN
void read_directory(const std::string &path);
#endif

void dynamicLoad(const std::string& prefix, const std::vector<std::string> & paths){
    for(std::vector<std::string>::const_iterator it = paths.begin(); it != paths.end(); it++)
        read_file(prefix + *it + SHAREDLIB_SUFFIX);
}
    


#ifndef DROPS_WIN
void read_directory(const std::string &path){
    /* not needed in the current implementation */
    DIR* dirfd= opendir(path.c_str());
    
    dirent* enumerator;
    while((enumerator = readdir(dirfd))){
        std::string filename(enumerator->d_name);
        if (filename.length() > 3 && filename.substr(filename.length()-3) == ".so") 
            read_file(path + "/" + filename);
    }
    if(errno) //relates to readdir()
        throw DROPSErrCL(strerror(0));
}
#endif

void read_file(const std::string& path){
#ifndef DROPS_WIN
    void* handle = dlopen(path.c_str(), RTLD_LAZY);
    
    if(!handle) throw DROPSErrCL((std::string)dlerror());
    
    dlerror(); //clear errors
    
    /* in case we ever want to load function symbols, this is how to do it
    
    fnc_t load = (fnc_t) dlsym(handle, "load");
    const char *dlsym_error = dlerror();
    
    if(dlsym_error)
        throw DROPSErrCL(dlsym_error);
    
    load(handle);
    
    */
    
    //no dlclose in our case so far (means dyn lib. gets cleaned when program stops..)
    //dlclose(handle);
#else
    HINSTANCE handle = LoadLibrary (path.c_str());
    if(!handle) throw DROPSErrCL("windows: couldn't load dll...");
#endif
}

}
