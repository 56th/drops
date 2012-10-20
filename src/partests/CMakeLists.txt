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
	
	exec(dist_modify geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer)
	
  exec(dist_simplexsend geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer)
    
  exec(dist_interface geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer)

	exec(dist_ref geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize out-vtkOut parallel-parmultigrid num-fe num-interfacePatch levelset-levelset levelset-fastmarch levelset-surfacetension geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scopetimer)

	exec(dist_migration geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize out-vtkOut parallel-parmultigrid num-fe num-interfacePatch levelset-levelset levelset-fastmarch levelset-surfacetension geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-funcmap misc-vectorFunctions misc-scalarFunctions stokes-instatstokes2phase misc-scopetimer)

	exec(dist_TestRefPar geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize out-vtkOut parallel-parmultigrid num-fe num-interfacePatch levelset-levelset levelset-fastmarch levelset-surfacetension misc-params out-output geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scopetimer)

	exec(dist_ref_rising_drop geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize out-vtkOut parallel-parmultigrid num-fe num-interfacePatch levelset-levelset levelset-fastmarch levelset-surfacetension misc-params out-output geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scopetimer)

	exec(dist_exchange geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer)

	exec(dist_sdropsP2 geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem misc-utils num-discretize num-interfacePatch num-fe geom-reftetracut geom-principallattice misc-scopetimer)

	exec(dist_interp2 partests-dist_interp2 geom-simplex geom-multigrid geom-topo out-output num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils num-discretize)
	
	exec(dist_interp1 geom-simplex geom-multigrid geom-topo out-output num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils num-discretize)
endif(${MPI})
