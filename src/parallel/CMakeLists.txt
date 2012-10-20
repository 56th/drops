set(HOME parallel)

libs_par(decompose exchange loadbal migrateunknowns parallel parmultigrid partime)

if(${MPI})
	target_link_libraries(parallel-parmultigrid geom-multigrid)
endif(${MPI})
