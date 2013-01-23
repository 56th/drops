#pragma once

#include <errno.h> //errno
#include <string.h> //various string functions (strcat…)
#include <vector>
#include <string>
#include <algorithm>
#include <misc/utils.h>

#ifndef DROPS_WIN

#include <dirent.h> //for readdir, DIR
#include <dlfcn.h> //dlopen

#else

#include <windows.h> 

#endif

namespace DROPS{
    void dynamicLoad(const std::string&, const std::vector<std::string>&);
}
