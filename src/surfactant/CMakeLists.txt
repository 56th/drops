set(HOME surfactant)

libs(ifacetransp)

exec_ser(surfactant geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-scopetimer misc-progressaccu misc-utils out-output num-discretize misc-params num-interfacePatch levelset-fastmarch surfactant-ifacetransp num-fe out-ensightOut levelset-surfacetension out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload misc-funcmap misc-scopetimer num-gradient_recovery levelset-levelsetmapper)

add_dependencies(surfactant misc-scalarFunctions misc-vectorFunctions)

exec_ser(surfacenorms geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-scopetimer misc-progressaccu misc-utils out-output num-discretize misc-params num-interfacePatch levelset-fastmarch surfactant-ifacetransp num-fe out-ensightOut levelset-surfacetension out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload misc-scopetimer)

add_dependencies(surfacenorms misc-scalarFunctions misc-vectorFunctions)