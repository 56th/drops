#CMakeLists.txt for DiST package


#LIBRARIES
if(MPI)
	add_library(DiST-DiST SHARED DiST)
	add_library(DiST-mpistream SHARED mpistream)
	add_library(DiST-remotedata SHARED remotedata)
endif(MPI)
#EXECUTABLES
