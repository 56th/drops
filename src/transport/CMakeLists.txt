set(HOME transport)
libs_ser(transportCoeff transportNitsche)

exec_ser(ns_transp transport-transportCoeff transport-ns_transp geom-boundary geom-builder geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize navstokes-instatnavstokes2phase misc-scalarFunctions misc-vectorFunctions misc-funcmap misc-params levelset-fastmarch stokes-instatstokes2phase num-MGsolver num-fe out-ensightOut transport-transportNitsche num-interfacePatch levelset-fastmarch num-fe out-vtkOut levelset-surfacetension geom-simplex surfactant-ifacetransp geom-geomselect  levelset-twophaseutils num-hypre levelset-coupling geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature num-renumber num-oseenprecond)
