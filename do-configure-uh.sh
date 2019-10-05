#!/bin/bash

mkdir bin
cd bin

cmake \
	-D TRILINOS_PATH:FILEPATH="/home/alexander/Documents/Amanzi/amanzi/install/tpls/trilinos-12-12-1" \
	-D Matlab_ROOT_DIR:FILEPATH="/opt/matlab" \
	../src/

cd surfnavierstokes
cat > runml << EOF
	#/bin/sh
	{ time -p nice -n 15 ./surfnavierstokes ; } |& mail -s "Process done!" alex@math.uh.edu &
EOF
chmod +x runml
