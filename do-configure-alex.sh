#!/bin/bash

mkdir bin
cd bin

cmake \
	-D VTK_DIR:FILEPATH="/usr/lib/cmake/vtk-7.1" \
	-D TRILINOS_PATH:FILEPATH="/home/alexander/Documents/Amanzi/amanzi/install/tpls/trilinos-12-12-1" \
	-D Matlab_ROOT_DIR:FILEPATH="/usr/local/MATLAB/R2019a" \
	-D DROPS_BUILD_TYPE:STRING="RELEASE" \
	../src/
