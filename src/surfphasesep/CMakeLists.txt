set(HOME surfphasesep)

add_library(${HOME}-ifacetransp STATIC ../surfactant/ifacetransp)

target_link_libraries(${HOME}-ifacetransp levelset-levelsetmapper)

exec_ser(surfphasesep geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset levelset-marking_strategy levelset-adaptriang misc-scopetimer misc-progressaccu misc-utils out-output num-discretize misc-params num-interfacePatch levelset-fastmarch surfactant-ifacetransp num-fe out-ensightOut levelset-surfacetension out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload misc-funcmap misc-scopetimer)

exec_ser(separation geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset levelset-marking_strategy levelset-adaptriang misc-scopetimer misc-progressaccu misc-utils out-output num-discretize misc-params num-interfacePatch levelset-fastmarch surfactant-ifacetransp num-fe out-ensightOut levelset-surfacetension out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload misc-funcmap misc-scopetimer num-gradient_recovery levelset-levelsetmapper)

if(NOT MPI)
add_dependencies(surfactant misc-scalarFunctions misc-vectorFunctions)
add_dependencies(surfphasesep misc-scalarFunctions misc-vectorFunctions levelset-twophaseCoeff)
add_dependencies(separation misc-scalarFunctions misc-vectorFunctions levelset-twophaseCoeff)
endif(NOT MPI)

add_my_custom_targets(surfphasesep)
add_my_custom_targets(separation)

add_definitions(-Wfatal-errors)
