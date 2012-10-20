set(HOME partests)

if(${MPI})
	SET(PAR_OBJ DiST-DiST DiST-mpistream DiST-remotedata
		parallel-parallel parallel-partime
		parallel-exchange parallel-parmultigrid
		parallel-loadbal parallel-migrateunknowns
		parallel-decompose)
		
	exec(dist_streamtest geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer)
	
	exec(dist_remotedatatest geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer)
	
	exec(dist_sendtetraofbrick geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer)
endif(${MPI})
