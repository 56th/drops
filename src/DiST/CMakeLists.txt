set(HOME DiST)
libs_par(DiST mpistream remotedata)

if(${MPI})
	target_link_libraries(DiST-remotedata DiST-DiST)
	target_link_libraries(DiST-DiST geom-simplex num-unknowns)
endif(${MPI})
