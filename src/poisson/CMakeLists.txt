set(HOME poisson)

set_par(PAR_OBJ DiST-DiST DiST-mpistream DiST-remotedata parallel-parallel parallel-partime parallel-exchange parallel-parmultigrid parallel-loadbal parallel-migrateunknowns parallel-decompose)


libs(ale poisson transport2phase poissonParam)

add_library(poisson-poissonCoeff SHARED poissonCoeff)
#set_property(TARGET poisson-poissonCoeff PROPERTY SUFFIX "")


if(NOT WIN32)
    target_link_libraries(poisson-poissonCoeff)
endif(NOT WIN32)


exec(cdrdrops poisson-ale geom-boundary geom-builder out-ensightOut geom-simplex geom-multigrid num-unknowns geom-topo poisson-poisson misc-problem misc-progressaccu misc-utils out-output num-fe num-discretize num-interfacePatch geom-geomselect misc-params out-vtkOut misc-funcmap geom-deformation misc-scopetimer poisson-poissonParam misc-dynamicload)
add_dependencies(cdrdrops misc-scalarFunctions misc-vectorFunctions poisson-poissonCoeff)

add_my_custom_targets(poisson)
