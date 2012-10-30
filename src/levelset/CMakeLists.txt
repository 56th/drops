set(HOME levelset)

libs(coupling levelset fastmarch surfacetension twophaseCoeff twophaseutils)
libs_ser(filmCoeff)

target_link_libraries(levelset-levelset levelset-fastmarch)

if(${MPI})
	set(PAR_OBJ DiST-DiST DiST-mpistream DiST-remotedata
           parallel-parallel parallel-partime
           parallel-exchange parallel-parmultigrid
           parallel-loadbal parallel-migrateunknowns
           parallel-decompose)
else(${MPI})

	exec(reparam levelset-fastmarch levelset-levelset geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe out-ensightOut num-interfacePatch levelset-surfacetension geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)
	
	exec(lsshear geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize levelset-fastmarch num-fe out-ensightOut stokes-instatstokes2phase navstokes-instatnavstokes2phase num-MGsolver num-interfacePatch levelset-surfacetension levelset-coupling misc-funcmap geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scalarFunctions num-renumber)
	
	exec(surfTens geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize misc-params levelset-fastmarch stokes-instatstokes2phase navstokes-instatnavstokes2phase num-fe out-ensightOut num-oseenprecond num-MGsolver misc-funcmap num-interfacePatch levelset-surfacetension levelset-coupling misc-vectorFunctions geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scalarFunctions num-renumber)

	exec(prJump geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize misc-params levelset-fastmarch stokes-instatstokes2phase num-MGsolver num-fe  out-ensightOut num-oseenprecond num-interfacePatch levelset-surfacetension levelset-coupling misc-vectorFunctions misc-funcmap geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scalarFunctions num-renumber)
	
	exec(film geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize navstokes-instatnavstokes2phase misc-params levelset-fastmarch stokes-instatstokes2phase num-MGsolver num-fe out-ensightOut out-vtkOut num-oseenprecond num-interfacePatch levelset-surfacetension levelset-coupling misc-vectorFunctions misc-funcmap levelset-filmCoeff geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scalarFunctions num-renumber)

	exec(brick_transp geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize misc-params levelset-fastmarch num-fe poisson-transport2phase out-ensightOut out-vtkOut num-interfacePatch levelset-surfacetension stokes-instatstokes2phase misc-vectorFunctions misc-funcmap levelset-twophaseCoeff geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scalarFunctions)
endif(${MPI})

exec(twophasedrops geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize navstokes-instatnavstokes2phase misc-params levelset-fastmarch stokes-instatstokes2phase num-MGsolver num-fe out-ensightOut num-oseenprecond poisson-transport2phase num-interfacePatch out-vtkOut surfactant-ifacetransp levelset-surfacetension geom-geomselect  levelset-twophaseutils num-hypre levelset-coupling misc-vectorFunctions misc-scalarFunctions misc-funcmap levelset-twophaseCoeff geom-principallattice geom-reftetracut geom-subtriangulation num-renumber num-quadrature)
