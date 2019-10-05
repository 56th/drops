#!/bin/bash

mkdir bin
cd bin

cmake \
	-D TRILINOS_PATH:FILEPATH="/home/alexander/Documents/Amanzi/amanzi/install/tpls/trilinos-12-12-1" \
	-D Matlab_ROOT_DIR:FILEPATH="/usr/local/MATLAB/R2019a" \
	../src/
