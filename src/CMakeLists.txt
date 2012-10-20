cmake_minimum_required(VERSION 2.8)

project(Drops)

#config types definieren
#if(CMAKE_CONFIGURATION_TYPES)
#	set(CMAKE_CONFIGURATION_TYPES Debug Release)
#	set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING
#		"Reset the configurations to what we need"
#		FORCE)
#endif()

include(CMakeSettings.txt)

#tell cmake to always include the home dir
include_directories(${CMAKE_SOURCE_DIR})

#use openMP
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

#use MPI
option(MPI "use MPI" OFF)
if(MPI)
	message(STATUS "###  MPI active ###")
	find_package(MPI REQUIRED)
	include_directories(${MPI_INCLUDE_PATH})
	
	execute_process(COMMAND mpicxx --showme:compile OUTPUT_VARIABLE MPI_CXX OUTPUT_STRIP_TRAILING_WHITESPACE)
	execute_process(COMMAND mpicxx --showme:link OUTPUT_VARIABLE MPI_LD OUTPUT_STRIP_TRAILING_WHITESPACE)
	
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX} -I${PARMETIS_HOME}/include -I${METIS_HOME}/include -I${ZOLTAN_HOME} -I${SCOTCH_HOME}/include -I${HYPRE_HOME}/src/hypre/include -DMPICH_IGNORE_CXX_SEEK")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LD} -L${PARMETIS_HOME}/build/Linux-x86_64/libparmetis -L${METIS_HOME}/build/Linux-x86_64/libmetis -lm -lparmetis -lmetis")
	add_definitions("-D_PAR=1")
else(MPI)
	message(STATUS "### MPI disabled ###")
endif(MPI)


#compiler specific options
if(CMAKE_CXX_COMPILER MATCHES icpc)
	message(STATUS "### setting flags for INTEL compiler ###")
	set(CMAKE_CXX_FLAGS_DEBUG "-DVALARRAY_BUG -wd1572")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ip -axSSE4.2,SSSE3,SSE2 -fp-model fast=2 -vec-report0")
elseif(CMAKE_COMPILER_IS_GNUCXX)
	message(STATUS "### setting flags for GCC compiler ###")
	set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
	set(CMAKE_CXX_FLAGS_RELEASE "-O2 -finline-limit=600")
endif(CMAKE_CXX_COMPILER MATCHES icpc)

function(libs)
	foreach(lib ${ARGN})
		add_library(${HOME}-${lib} SHARED ${lib})
	endforeach(lib)
endfunction(libs)

function(libs_ser)
	if(NOT ${MPI})
		foreach(lib ${ARGN})
			add_library(${HOME}-${lib} SHARED ${lib})
		endforeach(lib)
	endif(NOT ${MPI})
endfunction(libs_ser)

function(libs_par)
	if(${MPI})
		foreach(lib ${ARGN})
			add_library(${HOME}-${lib} SHARED ${lib})
		endforeach(lib)
	endif(${MPI})
endfunction(libs_par)

function(target_link_libraries_par)
	target_link_libraries(${ARGN} ${PAR_OBJ} ${MPI_CXX_LIBRARIES})
endfunction(target_link_libraries_par)

set(PACKAGES parallel levelset poisson DiST geom num out misc stokes navstokes surfactant transport partests)

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
