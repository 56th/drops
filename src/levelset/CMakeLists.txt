set(HOME levelset)

libs(coupling levelset fastmarch surfacetension twophaseCoeff twophaseutils)
libs_ser(filmCoeff)

target_link_libraries(levelset-levelset levelset-fastmarch misc-progressaccu misc-scopetimer)

set_par(PAR_OBJ DiST-DiST DiST-mpistream DiST-remotedata parallel-parallel parallel-partime parallel-exchange parallel-parmultigrid parallel-loadbal parallel-migrateunknowns parallel-decompose)

exec_ser(reparam levelset-fastmarch levelset-levelset misc-progressaccu misc-scopetimer geom-deformation geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe out-ensightOut num-interfacePatch levelset-surfacetension geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

exec_ser(lsshear geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-progressaccu misc-scopetimer misc-utils out-output num-discretize levelset-fastmarch num-fe out-ensightOut stokes-instatstokes2phase navstokes-instatnavstokes2phase num-MGsolver num-interfacePatch levelset-surfacetension levelset-coupling misc-funcmap geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature num-renumber)

set(prJump-staticlibs ${BEGIN_STATIC_LIBS} misc-scalarFunctions misc-vectorFunctions ${END_STATIC_LIBS})
exec_ser(prJump geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-progressaccu misc-scopetimer misc-utils out-output num-discretize misc-params levelset-fastmarch stokes-instatstokes2phase num-MGsolver num-fe  out-ensightOut num-oseenprecond num-interfacePatch levelset-surfacetension levelset-coupling misc-funcmap geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature num-renumber)

set(film-staticlibs ${BEGIN_STATIC_LIBS} misc-scalarFunctions misc-vectorFunctions levelset-filmCoeff ${END_STATIC_LIBS})
exec_ser(film geom-boundary geom-builder geom-deformation geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-progressaccu misc-scopetimer misc-utils out-output num-discretize navstokes-instatnavstokes2phase misc-params levelset-fastmarch stokes-instatstokes2phase num-MGsolver num-fe out-ensightOut out-vtkOut num-oseenprecond num-interfacePatch levelset-surfacetension levelset-coupling misc-funcmap geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature num-renumber)

set(brick_transp-staticlibs ${BEGIN_STATIC_LIBS} misc-scalarFunctions misc-vectorFunctions levelset-twophaseCoeff ${END_STATIC_LIBS})
exec_ser(brick_transp geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-progressaccu misc-scopetimer misc-utils out-output num-discretize misc-params levelset-fastmarch num-fe poisson-transport2phase out-ensightOut out-vtkOut num-interfacePatch levelset-surfacetension stokes-instatstokes2phase misc-funcmap geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature )

set(twophasedrops-staticlibs ${BEGIN_STATIC_LIBS} misc-scalarFunctions misc-vectorFunctions levelset-twophaseCoeff ${END_STATIC_LIBS})
exec(twophasedrops geom-boundary geom-builder geom-deformation geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-progressaccu misc-scopetimer misc-utils out-output num-discretize navstokes-instatnavstokes2phase misc-params levelset-fastmarch stokes-instatstokes2phase num-MGsolver num-fe out-ensightOut num-oseenprecond poisson-transport2phase num-interfacePatch out-vtkOut surfactant-ifacetransp levelset-surfacetension geom-geomselect levelset-twophaseutils num-hypre levelset-coupling misc-funcmap geom-principallattice geom-reftetracut geom-subtriangulation num-renumber num-quadrature)

