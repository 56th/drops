# DROPS

DROPS is a research-driven software package for the simulation of 
two-phase flows, mass and surfactant transport in two-phase systems developed 
at the Chair of Numerical Mathematics at RWTH Aachen University. 
The discretization is based on special finite element methods (XFEM/CutFEM, traceFEM).

## How to build DROPS

### Prerequisites

You will need
* a C++11 compiler (e.g., gcc-5.2)
* CMake
* the [Boost library](http://www.boost.org/users/history/), preferably version 
  1.58. From version 1.59 on, using C++-style comments in JSON files (which we use 
  as parameter files) is not possible anymore. You can strip comments from JSON
  files by using `g++ -E -x c++ myparam.json > myparam-stripped.json` and removing
  the first lines starting with `#`.
  
### UH branch

You will need the following third party libraries (TPLs): Trilinos, VTK, and Matlab. To build Trilinos, you may use the following script (modify paths for your local machine):

```
cmake \
-DTPL_ENABLE_MPI=ON \
-DMPI_BASE_DIR=/usr/lib/x86_64-linux-gnu/openmpi \
-DCMAKE_INSTALL_PREFIX=~/trilinos/build \
-DBUILD_SHARED_LIBS=ON \
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION=ON \
-DTrilinos_ENABLE_Epetra=ON \
-DTrilinos_ENABLE_AztecOO=ON \
-DTrilinos_ENABLE_Ifpack=ON \
-DTrilinos_ENABLE_Ifpack2=ON \
-DTrilinos_ENABLE_Belos=ON \
-DTrilinos_ENABLE_Kokkos=ON \
~/trilinos/source

make -j2 install
```
Use Paraview version 5.7.0-RC1 or newer.

For UH compute nodes guide, check [Ilya's presentation](https://www.math.uh.edu/~ilya/gs/talk_comput.pdf) and [wiki](https://sites.google.com/view/josiclabwiki/home).

### Generating Makefiles

CMake is used to generate the build system based on Makefiles.

First, generate a build directory. Within this build directory, call 
`cmake /path/to/drops/src`. For example, the build directory could be a folder 
`bin` inside the `drops` folder:

```bash
mkdir bin
cd bin
cmake ../src
```

This will populate a directory structure inside `bin` with associated Makefiles.
These can be used to compile the different executables.

#### MPI version

To build the MPI version of DROPS (to be used on a parallel computer),
the graph partitioner library ParMETIS is used. Download ParMETIS 4.0.2
and install it on your machine. After that, modify the path to ParMETIS and METIS
(which is part of ParMETIS) in the `CMakeSettings.txt` in the `src` directory.

Of course, an MPI library is needed, which is often chosen by using a corresponding
`module load` command. The MPI library will then be automatically selected by 
CMake in the following. As DROPS is able to use OpenMP within a parallel node 
(i.e., hybrid parallelization), a thread-safe MPI library is the best choice.

As next step, you have to add the option `-DMPI=1` to the `cmake` command like this:
```bash
mkdir bin-par
cd bin-par
cmake -DMPI=1 ../src
```

### Compiling Executables

Inside the build directory, go to the executable's directory and use `make`
to build the executable. E.g., to compile `twophasedrops`:

```bash
cd levelset
make twophasedrops
```

Calling `make all` in the build directory will build all DROPS executables.

### Running Executable

Usually, you have to provide a parameter file in JSON format as a command line
option when calling the executable. Some example parameter files can be found
in the directory `param` (on the same level as `src`) and has the same structure
as the `src` directory. If the code for some executable `xyz` is located in 
`src/abc`, then the corresponding parameter files are stored in `param/abc/xyz/`.

For example, if the executable `twophasedrops`
resides in the directory `drops/bin/levelset`, from within that directory call

```bash
mkdir vtk  # for VTK output
./twophasedrops ../../param/levelset/twophasedrops/risingbutanoldroplet.json
```
