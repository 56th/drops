#!/bin/bash

rm -r bin
mkdir bin
cd bin

cmake \
	-D VTK_DIR:FILEPATH="/home/yerbol/Packages/VTK-build" \
	-D TRILINOS_PATH:FILEPATH="/home/yerbol/Packages/Trilinos-build" \
  	-D Matlab_ROOT_DIR:FILEPATH="/home/yerbol/MATLAB" \
	-D DROPS_BUILD_TYPE:STRING="RELEASE" \
	-D MPI:BOOL=OFF \
	../src/
