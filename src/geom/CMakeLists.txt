set(HOME geom)
libs(boundary builder geomselect multigrid principallattice reftetracut simplex subtriangulation topo)

if(${MPI})
	target_link_libraries(geom-simplex parallel-parmultigrid geom-topo)
	target_link_libraries(geom-multigrid parallel-partime)
endif(${MPI})
