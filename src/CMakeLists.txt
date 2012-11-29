

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
	set(DROPS_CXX_FLAGS_DEBUG "${IPC_CXX_FLAGS_DEBUG}")
	set(DROPS_CXX_FLAGS_RELEASE "${IPC_CXX_FLAGS_RELEASE}")
elseif(CMAKE_COMPILER_IS_GNUCXX)
	message(STATUS "### setting flags for GCC compiler ###")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed")
	set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-as-needed")
	set(DROPS_CXX_FLAGS_DEBUG "${GCC_CXX_FLAGS_DEBUG}")
	set(DROPS_CXX_FLAGS_RELEASE "${GCC_CXX_FLAGS_RELEASE}")
endif(CMAKE_CXX_COMPILER MATCHES icpc)

#apply build type specific options
if(DROPS_BUILD_TYPE STREQUAL "RELEASE")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DROPS_CXX_FLAGS_RELEASE}")
elseif(DROPS_BUILD_TYPE STREQUAL "DEBUG")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DROPS_CXX_FLAGS_DEBUG}")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DROPS_CXX_FLAGS}")

message(STATUS "### GLOBAL COMPILER FLAGS:       ${CMAKE_CXX_FLAGS} ###")
message(STATUS "### GLOBAL SHARED LINKER FLAGS: ${CMAKE_SHARED_LINKER_FLAGS} ###")	
message(STATUS "### GLOBAL EXE LINKER FLAGS:    ${CMAKE_EXE_LINKER_FLAGS} ###")	

###############################################################################
# libs* functions are specializations of add_library, which allow us to add
# multiple files in a single line. The ${HOME} directory is automatically 
# prepended.
# THESE FUNCTIONS REQUIRE ${HOME}

function(libs_ when)
	if((when STREQUAL "ALL") OR (${when} STREQUAL "PAR" AND MPI) OR (${when} STREQUAL "SER" AND NOT MPI))
		foreach(lib ${ARGN})
			add_library(${HOME}-${lib} SHARED ${lib})
		endforeach(lib)
	endif()
endfunction(libs_)

function(libs)
	libs_("ALL" ${ARGN})
endfunction(libs)

function(libs_ser)
	libs_("SER" ${ARGN})
endfunction(libs_ser)

function(libs_par)
	libs_("PAR" ${ARGN})
endfunction(libs_par)

###############################################################################
# target_link_libraries_par appends ${PAR_OBJ} ${MPI_CXX_LIBRARIES} to a call
# of target_link_libraries

function(target_link_libraries_ when)
	if((${when} STREQUAL "ALL" OR ${when} STREQUAL "PAR") AND MPI)
		target_link_libraries(${ARGN} ${PAR_OBJ} ${MPI_CXX_LIBRARIES} ${PARMETIS_LIBRARIES})
	elseif((${when} STREQUAL "ALL" OR ${when} STREQUAL "SER") AND NOT MPI)
		target_link_libraries(${ARGN})
	endif()
endfunction(target_link_libraries_)

function(target_link_libraries_par)
	target_link_libraries_("PAR" ${ARGN})
endfunction(target_link_libraries_par)

function(target_link_libraries_ser)
	target_link_libraries_("SER" ${ARGN})
endfunction(target_link_libraries_ser)

function(set_ when)
    if((${when} STREQUAL "ALL") OR (${when} STREQUAL "PAR" AND MPI) OR (${when} STREQUAL "SER" AND NOT MPI))
	    message("SET ${HOME} ${when} ${ARGN}")
	    set(${ARGN})
    endif()
endfunction(set_)

function(set_par)
	set_("PAR" ${ARGN})
endfunction(set_par)

function(set_ser)
	set_("SER" ${ARGN})
endfunction(set_ser)

###############################################################################
# the exec command can be used to conveniently add an executable. 
# Arguments: filename <link_objects...>
# First a library and an executable target with suitable names are created.
# Afterwards we add <link_objects...> as a requirement of the executable 
# THIS FUNCTION REQUIRES ${HOME} and ${PAR_OBJ}

function(exec_ when f)
	if((${when} STREQUAL "ALL") OR (${when} STREQUAL "PAR" AND MPI) OR (${when} STREQUAL "SER" AND NOT MPI))
		libs_(${when} ${f})
		add_executable(${f} ${f})
		target_link_libraries_(${when} ${f} ${HOME}-${f} ${ARGN})
	endif()
endfunction(exec_)

function(exec f)
	exec_("ALL" ${f} ${ARGN})
endfunction(exec)

function(exec_ser f)
	exec_("SER" ${f} ${ARGN})
endfunction(exec_ser)

function(exec_par f)
	exec_("PAR" ${f} ${ARGN})
endfunction(exec_par)

set(PACKAGES parallel levelset poisson DiST geom num out misc stokes navstokes surfactant transport partests tests)
string(REPLACE ";" " " PACKAGES_STRING "${PACKAGES}")

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

add_custom_target(stat COMMAND ls */*.cpp */*.tpp */*.h | xargs wc -l WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

#find_package(Doxygen)
#if(DOXYGEN_FOUND)
#	add_custom_target(doc COMMAND ${DOXYGEN_EXECUTABLE} dox.cfg WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMENT "Generating API documentation with Doxygen" VERBATIM)
#endif(DOXYGEN_FOUND)

find_package(Doxygen)
if(DOXYGEN_FOUND)
	FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/dox.cfg ${CMAKE_CURRENT_BINARY_DIR}/doxygen/dox.cfg)
	add_custom_target(doc
		${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxygen/dox.cfg
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
		COMMENT "Generating API documentation with Doxygen" VERBATIM
	)
endif()
