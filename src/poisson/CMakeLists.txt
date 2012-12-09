set(HOME poisson)

if(MPI)
set(PAR_OBJ DiST-DiST DiST-mpistream DiST-remotedata parallel-parallel parallel-partime parallel-exchange parallel-parmultigrid parallel-loadbal parallel-migrateunknowns parallel-decompose)
endif(MPI)

libs(ale poissonCoeff poisson transport2phase)

exec(cdrdrops poisson-poissonCoeff poisson-ale poisson-cdrdrops geom-boundary geom-builder out-ensightOut geom-simplex geom-multigrid num-unknowns geom-topo poisson-poisson misc-problem misc-utils out-output num-fe num-discretize num-interfacePatch geom-geomselect misc-params out-vtkOut misc-funcmap misc-scalarFunctions misc-vectorFunctions geom-deformation misc-scopetimer)
