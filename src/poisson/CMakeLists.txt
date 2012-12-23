set(HOME poisson)

set_par(PAR_OBJ DiST-DiST DiST-mpistream DiST-remotedata parallel-parallel parallel-partime parallel-exchange parallel-parmultigrid parallel-loadbal parallel-migrateunknowns parallel-decompose)


libs(ale poissonCoeff poisson transport2phase)

set(cdrdrops-staticlibs ${BEGIN_STATIC_LIBS} misc-scalarFunctions misc-vectorFunctions poisson-poissonCoeff ${END_STATIC_LIBS})
exec(cdrdrops poisson-ale poisson-cdrdrops geom-boundary geom-builder out-ensightOut geom-simplex geom-multigrid num-unknowns geom-topo poisson-poisson misc-problem misc-utils out-output num-fe num-discretize num-interfacePatch geom-geomselect misc-params out-vtkOut misc-funcmap geom-deformation misc-scopetimer)
