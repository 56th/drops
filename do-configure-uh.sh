#!/bin/bash

rm -r bin
mkdir bin
cd bin

cmake \
	-D TRILINOS_PATH:FILEPATH="/shared/surfpde-nb/alex/trilinos/build" \
	-D VTK_DIR:FILEPATH="/usr/lib64/cmake/vtk" \
	-D Matlab_ROOT_DIR:FILEPATH="/opt/matlab" \
        -D DROPS_BUILD_TYPE:STRING="RELEASE" \
	-D MPI:BOOL=OFF \
	../src/

cd surfnavierstokes
cat > runml << EOF
	#/bin/sh
	{ time -p nice -n 15 ./surfnavierstokes ; } |& mail -s "Process done!" alex@math.uh.edu &
EOF
chmod +x runml
