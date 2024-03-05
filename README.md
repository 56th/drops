# DROPS

DROPS is a research-driven software package for the simulation of 
two-phase flows, mass and surfactant transport in two-phase systems developed 
at the Chair of Numerical Mathematics at RWTH Aachen University. 
The discretization is based on special finite element methods (XFEM/CutFEM, traceFEM).

Here you can find the models which were added by a research group at the University of Houston. Also there is a list of related published papers in the end of this page. 

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

## UH branch instructions

This is a set of instructions for reproducing the results from papers listed below.

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

For UH compute nodes guide, check [wiki](https://sites.google.com/view/josiclabwiki/home).


### Publications
1. V.  Yushutin, A. Quaini, S. Majd, M. Olshanskii, A computational study of lateral phase separation in biological membranes, International Journal for Numerical Methods in Biomedical Engineering, V. 35 (2019), e3182;  doi: 10.1002/cnm.3181;
2. M. Olshanskii, A.  Reusken, A.  Zhiliakov, Inf-sup stability of the trace P2-P1 Taylor-Hood elements for surface PDEs, Mathematics of Computation V.  90 (2021), p. 1527-1555; doi: 10.1090/mcom/3551;
3. A.  Zhiliakov, Y.  Wang, A. Quaini, M. Olshanskii, S. Majd, Experimental validation of a phase-field model to predict coarsening dynamics of lipid domains in multicomponent membranes, BBA - Biomembranes, V.  1863 (2021), Article 183446;   doi: 10.1016/j.bbamem.2020.183446;
4. Y.  Palzhanov, A.  Zhiliakov, A.  Quaini, M.  Olshanskii, A decoupled, stable, and linear FEM for a phase-field model of variable density two-phase incompressible surface flow, Computer Methods in Applied Mechanics and Engineering, V.  387 (2021), Article 114167;   doi: 10.1016/j.cma.2021.114167;
5. M. Olshanskii, Y.  Palzhanov, Q.  Sun, A comparison of Cahn-Hilliard and Navier-Stokes-Cahn-Hilliard models on manifolds, Vietnam Journal of Mathematics,   doi: 10.1007/s10013-022-00564-5;
6. M. Olshanskii, A.  Zhiliakov, Recycling augmented Lagrangian preconditioner in an incompressible fluid solver, Numerical Linear Algebra with Applications , V.  1864 (2022), Article 183898; doi: 10.1002/nla.2415;
7. Y.  Wang, Y.  Palzhanov, A. Quaini, M. Olshanskii, S. Majd, Lipid domain coarsening and fluidity in multicomponent lipid vesicles: A continuum based model and its experimental validation,   BBA - Biomembranes, V.  29 (2022), Article e2415;   doi: 10.1016/j.bbamem.2022.183898;
8. M. Olshanskii, A. Reusken, A. Zhilyakov, Tangential Navier-Stokes equations on evolving surfaces: Analysis and simulations, Mathematical Models and Methods in Applied Sciences V.  14 (2022), 2817-2852; doi: 10.1142/S0218202522500658;
9. M. Olshanskii, Y. Palzhanov, A. Quaini, A scalar auxiliary variable unfitted FEM for the surface Cahn-Hilliard equation, Journal of Scientific Computing V.  97 (2023), Article 57; doi: 10.1007/s10915-023-02370-8;
10. Y. Wang, Y. Palzhanov, D Dang, A. Quaini, M. Olshanskii, S. Majd, On fusogenicity of positively charged phased-separated lipid vesicles: experiments and computational simulations, Biomolecules V.  13 (2023), Article 1473; doi: 10.3390/biom13101473; 
