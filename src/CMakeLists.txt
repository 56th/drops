cmake_minimum_required(VERSION 2.8)

project(Drops)

#config types definieren
#if(CMAKE_CONFIGURATION_TYPES)
#	set(CMAKE_CONFIGURATION_TYPES Debug Release)
#	set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING
#		"Reset the configurations to what we need"
#		FORCE)
#endif()

cmake_policy(SET CMP0015 NEW) #allows us to give relative paths in ${PARMETIS_HOME} and ${METIS_HOME} (when given relative paths in librarys these will be interpreted as relative to the root source directory)
include(CMakeSettings.txt)

#tell cmake to always include the home dir
include_directories(${CMAKE_SOURCE_DIR})

#use openMP
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

#use MPI
if(MPI)
	message(STATUS "###  MPI active ###")
	find_package(MPI REQUIRED)
	
	set(METIS_INCLUDE ${METIS_HOME}/include)
	set(METIS_LIBRARY ${METIS_HOME}/build/Linux-x86_64/libmetis/libmetis.a)
	set(PARMETIS_INCLUDE ${PARMETIS_HOME}/include)
	set(PARMETIS_LIBRARY ${PARMETIS_HOME}/build/Linux-x86_64/libparmetis/libparmetis.a)

	include_directories(${MPI_INCLUDE_PATH} ${PARMETIS_INCLUDE} ${METIS_INCLUDE})
	
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS} -DMPICH_IGNORE_CXX_SEEK")
	set(PARMETIS_LIBRARIES ${MPI_CXX_LIBRARIES} ${METIS_LIBRARY} ${PARMETIS_LIBRARY} m)
	add_definitions("-D_PAR=1")
else(MPI)
	message(STATUS "### MPI disabled ###")
endif(MPI)

#compiler specific options
if(CMAKE_CXX_COMPILER MATCHES icpc)
	message(STATUS "### setting flags for INTEL compiler ###")
	set(CMAKE_CXX_FLAGS_DEBUG "-DVALARRAY_BUG -wd1572")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ip -xhost -fp-model fast=2")
elseif(CMAKE_COMPILER_IS_GNUCXX)
	message(STATUS "### setting flags for GCC compiler ###")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed")
	set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-as-needed")
	set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
	set(CMAKE_CXX_FLAGS_RELEASE "-O2 -finline-limit=600")
endif(CMAKE_CXX_COMPILER MATCHES icpc)

###############################################################################
# libs* functions are specializations of add_library, which allow us to add
# multiple files in a single line. The ${HOME} directory is automatically 
# prepended.
# THESE FUNCTIONS REQUIRE ${HOME}

function(libs)
	foreach(lib ${ARGN})
		add_library(${HOME}-${lib} SHARED ${lib})
	endforeach(lib)
endfunction(libs)

function(libs_ser)
	if(NOT MPI)
		foreach(lib ${ARGN})
			add_library(${HOME}-${lib} SHARED ${lib})
		endforeach(lib)
	endif(NOT MPI)
endfunction(libs_ser)

function(libs_par)
	if(MPI)
		foreach(lib ${ARGN})
			add_library(${HOME}-${lib} SHARED ${lib})
		endforeach(lib)
	endif(MPI)
endfunction(libs_par)

###############################################################################
# target_link_libraries_par appends ${PAR_OBJ} ${MPI_CXX_LIBRARIES} to a call
# of target_link_libraries

function(target_link_libraries_par)
	target_link_libraries(${ARGN} ${PAR_OBJ} ${MPI_CXX_LIBRARIES} ${PARMETIS_LIBRARIES})
endfunction(target_link_libraries_par)

###############################################################################
# the exec command can be used to conveniently add an executable. 
# Arguments: filename <link_objects...>
# First a library and an executable target with suitable names are created.
# Afterwards we add <link_objects...> as a requirement of the executable 
# THIS FUNCTION REQUIRES ${HOME} and ${PAR_OBJ}

function(exec f)
	libs(${f})
	add_executable(${f} ${f})
	target_link_libraries_par(${f} ${HOME}-${f} ${ARGN})
endfunction(exec)

function(exec_ser f)
	if(NOT MPI)
		libs_ser(${f})
		add_executable(${f} ${f})
		target_link_libraries_par(${f} ${HOME}-${f} ${ARGN})
	endif(NOT MPI)
endfunction(exec_ser)

function(exec_par f)
	if(MPI)
		libs_par(${f})
		add_executable(${f} ${f})
		target_link_libraries_par(${f} ${HOME}-${f} ${ARGN})
	endif(MPI)
endfunction(exec_par)

set(PACKAGES parallel levelset poisson DiST geom num out misc stokes navstokes surfactant transport partests tests)

option(TESTS "compile tests" OFF)
if(TESTS)
	list(APPEND PACKAGES tests)
	enable_testing()
endif(TESTS)

foreach(package ${PACKAGES})
	add_subdirectory(${package})
	file(GLOB jsonfiles ${package}/*.json)
	foreach(file ${jsonfiles})
		  configure_file(${file} ${package}/)
		  #message(STATUS "copying ${file} ${package}/")
	endforeach(file)
endforeach(${package})

message(STATUS "### GLOBAL COMPILE FLAGS: ${CMAKE_CXX_FLAGS} ###")
message(STATUS "### GLOBAL LINKER FLAGS: ${CMAKE_EXE_LINKER_FLAGS} ###")	
