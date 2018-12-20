set(HOME surfnavierstokes)

add_library(${HOME}-ifacetransp STATIC ../surfactant/ifacetransp)

target_link_libraries(${HOME}-ifacetransp levelset-levelsetmapper)

exec_ser(surfnavierstokes geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset levelset-marking_strategy levelset-adaptriang misc-scopetimer misc-progressaccu misc-utils out-output num-discretize misc-params num-interfacePatch levelset-fastmarch surfactant-ifacetransp num-fe out-ensightOut levelset-surfacetension out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload misc-funcmap misc-scopetimer)

if(NOT MPI)
add_dependencies(surfactant misc-scalarFunctions misc-vectorFunctions)
add_dependencies(surfnavierstokes misc-scalarFunctions misc-vectorFunctions levelset-twophaseCoeff)
endif(NOT MPI)

add_my_custom_targets(surfnavierstokes)
