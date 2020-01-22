#!/bin/bash

mkdir bin
cd bin

cmake \
	-D VTK_DIR:FILEPATH="/usr/local/vtk" \
	-D TRILINOS_PATH:FILEPATH="/home/alexander/trilinos/build" \
  	-D Matlab_ROOT_DIR:FILEPATH="/usr/local/MATLAB/R2019a" \
	-D DROPS_BUILD_TYPE:STRING="RELEASE" \
	-D MPI:BOOL=OFF \
	../src/
