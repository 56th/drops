set(HOME transport)
libs_ser(transportNitsche)

add_library(transport-transportCoeff SHARED transportCoeff)

#exec_ser(ns_transp transport-transportCoeff geom-boundary geom-builder geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize navstokes-instatnavstokes2phase misc-scalarFunctions misc-vectorFunctions misc-funcmap misc-params levelset-fastmarch stokes-instatstokes2phase num-MGsolver num-fe out-ensightOut transport-transportNitsche num-interfacePatch levelset-fastmarch num-fe out-vtkOut levelset-surfacetension geom-simplex surfactant-ifacetransp geom-geomselect  levelset-twophaseutils num-hypre levelset-coupling geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature num-renumber num-oseenprecond)

exec_ser(ns_transp
    transport-transportCoeff
    geom-boundary geom-builder geom-multigrid
    geom-deformation
    num-unknowns geom-topo num-fe misc-problem levelset-levelset
    levelset-adaptriang levelset-marking_strategy
    misc-scopetimer misc-progressaccu
    misc-utils out-output num-discretize
    navstokes-instatnavstokes2phase misc-scalarFunctions
    misc-vectorFunctions misc-funcmap misc-params
    levelset-fastmarch stokes-instatstokes2phase num-MGsolver
    num-fe out-ensightOut transport-transportNitsche
    num-interfacePatch levelset-fastmarch num-fe
    out-vtkOut levelset-surfacetension geom-simplex surfactant-ifacetransp
    geom-geomselect  levelset-twophaseutils num-hypre levelset-coupling
    geom-principallattice geom-reftetracut geom-subtriangulation
    num-quadrature num-renumber num-oseenprecond
    misc-dynamicload stokes-slipBndOnePhase)
if(NOT MPI)
add_dependencies(ns_transp transport-transportCoeff misc-scalarFunctions misc-vectorFunctions spacetimetransp-sttranspCoeff)
endif(NOT MPI)

add_my_custom_targets(transport)
