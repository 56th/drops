set(HOME surfactant)

libs(ifacetransp)

exec_ser(surfactant surfactant-surfactant geom-boundary geom-builder geom-simplex geom-multigrid 
    num-unknowns geom-topo num-fe misc-problem levelset-levelset 
    misc-utils out-output num-discretize misc-params num-interfacePatch 
    levelset-fastmarch surfactant-ifacetransp 
    num-fe out-ensightOut levelset-surfacetension out-vtkOut 
    geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)
    
exec_ser(surfacenorms
    surfactant-surfacenorms geom-boundary geom-builder geom-simplex geom-multigrid 
    num-unknowns geom-topo num-fe misc-problem levelset-levelset 
    misc-utils out-output num-discretize misc-params num-interfacePatch 
    levelset-fastmarch surfactant-ifacetransp 
    num-fe out-ensightOut levelset-surfacetension out-vtkOut 
    geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)
