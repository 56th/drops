set(HOME partests)

if(${MPI})
	SET(PAR_OBJ DiST-DiST DiST-mpistream DiST-remotedata
		parallel-parallel parallel-partime
		parallel-exchange parallel-parmultigrid
		parallel-loadbal parallel-migrateunknowns
		parallel-decompose)
		
	libs_par(dist_streamtest dist_remotedatatest dist_sendtetraofbrick)

	add_executable(dist_streamtest dist_streamtest)
	target_link_libraries_par(dist_streamtest partests-dist_streamtest geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer)
	
	add_executable(dist_remotedatatest dist_remotedatatest)
	target_link_libraries_par(dist_remotedatatest partests-dist_remotedatatest geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer)
	
	add_executable(dist_sendtetraofbrick dist_sendtetraofbrick)
	target_link_libraries_par(dist_sendtetraofbrick partests-dist_sendtetraofbrick geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer)
endif(${MPI})
