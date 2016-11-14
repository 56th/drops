cmake_minimum_required(VERSION 2.8.1)

project(drops)

cmake_policy(SET CMP0015 NEW) #allows us to give relative paths in ${PARMETIS_HOME} and ${METIS_HOME} (when given relative paths in librarys these will be interpreted as relative to the root source directory)
include(CMakeSettings.txt)

#tell cmake to always include the home dir
include_directories(${CMAKE_SOURCE_DIR})

#use openMP
message(STATUS "### GENERATOR is ${CMAKE_GENERATOR}")
if(CMAKE_GENERATOR MATCHES "Visual Studio 10")
  set(OpenMP_CXX_FLAGS "/openmp")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
else()
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    find_package(OpenMP REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

#use MPI
if(MPI)
    message(STATUS "###  MPI active ###")
    find_package(MPI REQUIRED)

    set(METIS_INCLUDE ${METIS_HOME}/include)
    set(METIS_LIBRARY ${METIS_HOME}/build/Linux-x86_64/libmetis/libmetis.a)
    set(PARMETIS_INCLUDE ${PARMETIS_HOME}/include)
    set(PARMETIS_LIBRARY ${PARMETIS_HOME}/build/Linux-x86_64/libparmetis/libparmetis.a)

    include_directories(${MPI_INCLUDE_PATH} ${PARMETIS_INCLUDE} ${METIS_INCLUDE})

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS} -DMPICH_IGNORE_CXX_SEEK -DOMPI_SKIP_MPICXX")
    set(PARMETIS_LIBRARIES ${MPI_CXX_LIBRARIES} ${METIS_LIBRARY} ${PARMETIS_LIBRARY} m)
    add_definitions("-D_PAR=1")
else(MPI)
    message(STATUS "### MPI disabled ###")
endif(MPI)

if(WIN32)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DDROPS_WIN  /I${BOOST_HOME}")
else(WIN32)
  if (BOOST_HOME STREQUAL "")
  else (BOOST_HOME STREQUAL "")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${BOOST_HOME}")
  endif (BOOST_HOME STREQUAL "")

  if (MKL_HOME STREQUAL "")
  else (MKL_HOME STREQUAL "")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDROPS_PARDISO -DMKL_ILP64 -m64 -I${MKL_HOME}/include -Wno-long-long")
  endif (MKL_HOME STREQUAL "")

endif(WIN32)

#compiler specific options
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    message(STATUS "### setting flags for INTEL compiler ###")
    set(DROPS_CXX_FLAGS_DEBUG "${IPC_CXX_FLAGS_DEBUG}")
    set(DROPS_CXX_FLAGS_RELEASE "${IPC_CXX_FLAGS_RELEASE}")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "### setting flags for GCC compiler ###")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--as-needed")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--as-needed")
    set(DROPS_CXX_FLAGS_DEBUG "${GCC_CXX_FLAGS_DEBUG}")
    set(DROPS_CXX_FLAGS_RELEASE "${GCC_CXX_FLAGS_RELEASE}")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "XL")
    message(STATUS "### setting flags for IBM XL compiler ###")
    set(DROPS_CXX_FLAGS_DEBUG "${XL_CXX_FLAGS_DEBUG}")
    set(DROPS_CXX_FLAGS_RELEASE "${XL_CXX_FLAGS_RELEASE}")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(STATUS "### setting flags for CLANG compiler ###")
    set(DROPS_CXX_FLAGS_DEBUG "${CLANG_CXX_FLAGS_DEBUG}")
    set(DROPS_CXX_FLAGS_RELEASE "${CLANG_CXX_FLAGS_RELEASE}")
else()
    message(STATUS "### could not determine compiler ${CMAKE_CXX_COMPILER_ID}, it seems to be neither gcc nor intel")
endif()

if(NOT DROPS_BUILD_TYPE)
    set(DROPS_BUILD_TYPE "DEBUG")
    message(STATUS "### Setting DEFAULT build type: ${DROPS_BUILD_TYPE} ###")
else(NOT DROPS_BUILD_TYPE)
    message(STATUS "### BUILD TYPE:                 ${DROPS_BUILD_TYPE} ###")
endif(NOT DROPS_BUILD_TYPE)

#apply build type specific options
if(DROPS_BUILD_TYPE STREQUAL "RELEASE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DROPS_CXX_FLAGS_RELEASE}")
elseif(DROPS_BUILD_TYPE STREQUAL "DEBUG")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DROPS_CXX_FLAGS_DEBUG}")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DROPS_CXX_FLAGS} -fPIC")

message(STATUS "### CXX_COMPILER:               ${CMAKE_CXX_COMPILER}")
message(STATUS "### GLOBAL COMPILER FLAGS:      ${CMAKE_CXX_FLAGS} ###")
message(STATUS "### GLOBAL SHARED LINKER FLAGS: ${CMAKE_SHARED_LINKER_FLAGS} ###")
message(STATUS "###    GLOBAL EXE LINKER FLAGS: ${CMAKE_EXE_LINKER_FLAGS} ###")

###############################################################################
# libs* functions are specializations of add_library, which allow us to add
# multiple files in a single line. The ${HOME} directory is automatically
# prepended.
# THESE FUNCTIONS REQUIRE ${HOME}

function(libs)
    foreach(lib ${ARGN})
        add_library(${HOME}-${lib} STATIC ${lib})
    endforeach(lib)
endfunction(libs)

function(libs_ser)
    if (NOT MPI)
        libs(${ARGN})
    endif()
endfunction(libs_ser)

function(libs_par)
    if (MPI)
        libs(${ARGN})
    endif()
endfunction(libs_par)

###############################################################################
# target_link_libraries_par appends ${PAR_OBJ} ${MPI_CXX_LIBRARIES}
# ${PARMETIS_LIBRARIES} to a call of target_link_libraries

function(target_link_libraries_ser)
    if (NOT MPI)
        target_link_libraries(${ARGN})
    endif()
endfunction(target_link_libraries_ser)

function(target_link_libraries_par)
    if (MPI)
        target_link_libraries(${ARGN} ${PAR_OBJ} ${MPI_CXX_LIBRARIES} ${PARMETIS_LIBRARIES})
    endif()
endfunction(target_link_libraries_par)

###############################################################################
# Set a variable depending on MPI

function(set_par)
    if (MPI)
        set(${ARGN} PARENT_SCOPE)
    endif()
endfunction(set_par)

if (WIN32)
#  set(BEGIN_STATIC_LIBS -OPT:NOREF)
#  set(END_STATIC_LIBS  -OPT:REF)
else(WIN32)
    set(BEGIN_STATIC_LIBS -Wl,--whole-archive)
    set(END_STATIC_LIBS  -Wl,--no-whole-archive)
endif(WIN32)

###############################################################################
# the exec command can be used to conveniently add an executable.
# Arguments: filename <link_objects...>
# First a library and an executable target with suitable names are created.
# Afterwards we add <link_objects...> as a requirement of the executable
# THIS FUNCTION REQUIRES ${HOME} and ${PAR_OBJ}

function(exec_ f)
#    libs(${f})
    add_executable(${f} ${f})
endfunction(exec_)

function(exec_ser f)
    if(NOT MPI)
        exec_(${f})
        target_link_libraries(${f} ${ARGN} ${${f}-staticlibs})
    endif()
endfunction(exec_ser)

function(exec_par f)
    if(MPI)
        exec_(${f})
        target_link_libraries_par(${f} ${ARGN} ${${f}-staticlibs})
    endif()
endfunction(exec_par)

function(exec f)
    exec_ser(${f} ${ARGN})
    exec_par(${f} ${ARGN})
endfunction(exec)

function(add_my_custom_targets f)
  add_custom_target(${f}-param-files 
    echo \"copy $f json files...\" \;
    cp -f ${CMAKE_CURRENT_SOURCE_DIR}/*.json ${CMAKE_CURRENT_BINARY_DIR}/ 2> /dev/null \; 
    COMMENT "copy param files from ${CMAKE_CURRENT_SOURCE_DIR}/ to ${CMAKE_CURRENT_BINARY_DIR}/"
    )
  add_custom_target(${f}-stat COMMAND ls ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/*.tpp ${CMAKE_CURRENT_SOURCE_DIR}/*.h | xargs wc -l WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_custom_target(${f}-coeffs DEPENDS misc-scalarFunctions misc-vectorFunctions levelset-twophaseCoeff levelset-filmCoeff poisson-poissonCoeff stokes-stokesCoeff)
endfunction(add_my_custom_targets f)

set(PACKAGES parallel levelset poisson DiST geom num out misc osmosis stokes navstokes surfactant spacetimetransp transport partests tests)
string(REPLACE ";" " " PACKAGES_STRING "${PACKAGES}")

option(TESTS "compile tests" OFF)
if(TESTS)
    list(APPEND PACKAGES tests)
    enable_testing()
endif(TESTS)

foreach(package ${PACKAGES})
    add_subdirectory(${package})
endforeach(${package})

add_custom_target(param-files 
    COMMAND for PACK in ${PACKAGES} \; do
              echo \"copy \$$PACK json files...\" \;
              cp -f ${CMAKE_CURRENT_SOURCE_DIR}/\$$PACK/*.json ${CMAKE_CURRENT_BINARY_DIR}/\$$PACK/ 2> /dev/null \; 
            done
    COMMENT "copy param files from ${CMAKE_CURRENT_SOURCE_DIR} to ${CMAKE_CURRENT_BINARY_DIR}"
)

add_custom_target(stat COMMAND ls */*.cpp */*.tpp */*.h | xargs wc -l WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_custom_target(distclean COMMAND rm -rf ${CMAKE_CURRENT_BINARY_DIR}/*)

add_custom_target(ctagsdb COMMAND ctags -R --tag-relative=yes --c++-types=+px --excmd=pattern --exclude=Makefile --exclude= --langmap=c++:+.tpp -f ${CMAKE_CURRENT_BINARY_DIR}/ctagsdb WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

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

add_custom_target(coeffs DEPENDS misc-scalarFunctions misc-vectorFunctions levelset-twophaseCoeff levelset-filmCoeff poisson-poissonCoeff stokes-stokesCoeff)
