set(HOME partests)

if(${MPI})
	SET(PAR_OBJ DiST-DiST DiST-mpistream DiST-remotedata
		parallel-parallel parallel-partime
		parallel-exchange parallel-parmultigrid
		parallel-loadbal parallel-migrateunknowns
		parallel-decompose)
		
	libs_par(dist_streamtest)

	add_executable(dist_streamtest dist_streamtest)
	target_link_libraries(dist_streamtest partests-dist_streamtest geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer ${PAR_OBJ} ${MPI_CXX_LIBRARIES})
endif(${MPI})
