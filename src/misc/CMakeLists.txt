#CMakeLists.txt for misc package


#LIBRARIES
add_library(misc-compensight SHARED compensight)
add_library(misc-funcmap SHARED funcmap)
add_library(misc-params SHARED params)
add_library(misc-problem SHARED problem)
add_library(misc-scalarFunctions SHARED scalarFunctions)
add_library(misc-scopetimer SHARED scopetimer)
add_library(misc-utils SHARED utils)
add_library(misc-vectorFunctions SHARED vectorFunctions)
add_library(misc-xfem SHARED xfem)

target_link_libraries(misc-problem num-interfacePatch)

#EXECUTABLES
if(NOT MPI)

        add_executable(compensight.exe compensight)
        target_link_libraries(compensight.exe misc-compensight misc-utils)
        
endif(NOT MPI)

