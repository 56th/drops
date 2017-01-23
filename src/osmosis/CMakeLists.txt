set(HOME osmosis)
libs_ser( osmosisSetup)

add_library(osmosis-osmosisCoeff SHARED osmosisCoeff)

if(NOT MPI)
    if(NOT WIN32)
        target_link_libraries(osmosis-osmosisCoeff)
    endif(NOT WIN32)
endif(NOT MPI)

exec_ser(osmosis misc-scopetimer geom-boundary geom-builder geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset levelset-marking_strategy levelset-adaptriang misc-utils out-output num-discretize navstokes-instatnavstokes2phase misc-funcmap misc-funcmap misc-params levelset-fastmarch stokes-instatstokes2phase num-MGsolver num-fe out-ensightOut stokes-integrTime osmosis-osmosisSetup num-interfacePatch levelset-fastmarch num-fe out-vtkOut levelset-surfacetension geom-simplex surfactant-ifacetransp geom-geomselect levelset-twophaseutils num-hypre levelset-coupling geom-principallattice geom-reftetracut geom-subtriangulation num-renumber num-quadrature misc-progressaccu misc-dynamicload stokes-slipBndOnePhase dl)

if(NOT MPI)
    add_dependencies(osmosis misc-scalarFunctions misc-vectorFunctions misc-csgFunctions osmosis-osmosisCoeff)
endif(NOT MPI)

add_my_custom_targets(osmosis)
