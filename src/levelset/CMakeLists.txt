set(HOME levelset)

libs(coupling levelset fastmarch surfacetension twophaseutils marking_strategy adaptriang volume_adjustment)

target_link_libraries(levelset-levelset levelset-volume_adjustment levelset-fastmarch misc-progressaccu misc-scopetimer)
target_link_libraries(levelset-volume_adjustment out-vtkOut)

add_library(levelset-twophaseCoeff SHARED twophaseCoeff)
add_library(levelset-filmCoeff SHARED filmCoeff)
#set_property(TARGET levelset-filmCoeff PROPERTY SUFFIX "")

if(NOT WIN32)
    target_link_libraries(levelset-twophaseCoeff)
    target_link_libraries(levelset-filmCoeff)
else(NOT WIN32)
    target_link_libraries(levelset-twophaseCoeff misc-funcmap misc-params)
endif(NOT WIN32)

if(MPI)
    set(PAR_OBJ DiST-DiST DiST-mpistream DiST-remotedata parallel-parallel parallel-partime parallel-exchange parallel-parmultigrid parallel-loadbal parallel-migrateunknowns parallel-decompose)
endif(MPI)

exec_ser(reparam levelset-fastmarch levelset-levelset misc-progressaccu misc-scopetimer geom-deformation geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe out-ensightOut num-interfacePatch levelset-surfacetension geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

exec_ser(lsshear geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-progressaccu misc-scopetimer misc-utils out-output num-discretize levelset-fastmarch num-fe out-ensightOut stokes-instatstokes2phase navstokes-instatnavstokes2phase num-MGsolver num-interfacePatch levelset-surfacetension levelset-coupling misc-funcmap geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature num-renumber stokes-slipBndOnePhase)

#set(prJump-staticlibs ${BEGIN_STATIC_LIBS} misc-scalarFunctions misc-vectorFunctions ${END_STATIC_LIBS})
exec_ser(prJump geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-progressaccu misc-scopetimer misc-utils out-output num-discretize misc-params levelset-fastmarch levelset-adaptriang levelset-marking_strategy stokes-instatstokes2phase num-MGsolver num-fe out-ensightOut out-vtkOut num-oseenprecond num-interfacePatch levelset-surfacetension levelset-coupling misc-funcmap geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature num-renumber misc-dynamicload stokes-slipBndOnePhase)
# convention: make sure that those coefficient-defining shared objects are built which ared used in the default configuration of the executable..
if(NOT MPI)
add_dependencies(prJump levelset-twophaseCoeff misc-scalarFunctions misc-vectorFunctions)
endif(NOT MPI)

#exec_ser(film geom-boundary geom-builder geom-deformation geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset levelset-adaptriang levelset-marking_strategy misc-progressaccu misc-scopetimer misc-utils out-output num-discretize navstokes-instatnavstokes2phase misc-params levelset-fastmarch stokes-instatstokes2phase num-MGsolver num-fe out-ensightOut out-vtkOut num-oseenprecond num-interfacePatch levelset-surfacetension levelset-coupling misc-funcmap geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature num-renumber misc-dynamicload stokes-slipBndOnePhase)
# convention: make sure that those coefficient-defining shared objects are built which ared used in the default configuration of the executable..
#if(NOT MPI)
#add_dependencies(film levelset-filmCoeff misc-scalarFunctions misc-vectorFunctions)
#endif(NOT MPI)

exec_ser(brick_transp geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-progressaccu misc-scopetimer misc-utils out-output num-discretize misc-params levelset-fastmarch num-fe poisson-transport2phase out-ensightOut out-vtkOut num-interfacePatch levelset-surfacetension stokes-instatstokes2phase misc-funcmap geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload stokes-slipBndOnePhase)
# convention: make sure that those coefficient-defining shared objects are built which ared used in the default configuration of the executable..
if(NOT MPI)
add_dependencies(brick_transp misc-scalarFunctions misc-vectorFunctions levelset-twophaseCoeff)
endif(NOT MPI)

exec_par(twophasedrops geom-boundary geom-builder geom-deformation geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset levelset-adaptriang levelset-marking_strategy misc-progressaccu misc-scopetimer misc-utils out-output num-discretize navstokes-instatnavstokes2phase misc-params levelset-fastmarch stokes-instatstokes2phase stokes-slipBndOnePhase num-MGsolver num-fe out-ensightOut num-oseenprecond poisson-transport2phase num-interfacePatch surfactant-ifacetransp levelset-surfacetension geom-geomselect levelset-twophaseutils num-hypre levelset-coupling misc-funcmap geom-principallattice geom-reftetracut geom-subtriangulation num-renumber num-quadrature misc-dynamicload out-vtkOut)

exec_ser(twophasedrops geom-boundary geom-builder geom-deformation geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset levelset-adaptriang levelset-marking_strategy misc-progressaccu misc-scopetimer misc-utils out-output num-discretize navstokes-instatnavstokes2phase misc-params levelset-fastmarch stokes-instatstokes2phase stokes-slipBndOnePhase num-MGsolver num-fe out-ensightOut num-oseenprecond poisson-transport2phase num-interfacePatch surfactant-ifacetransp levelset-surfacetension geom-geomselect levelset-twophaseutils num-hypre levelset-coupling misc-funcmap geom-principallattice geom-reftetracut geom-subtriangulation num-renumber num-quadrature misc-dynamicload out-vtkOut num-stokespardiso)


# convention: make sure that those coefficient-defining shared objects are built which ared used in the default configuration of the executable..
add_dependencies(twophasedrops levelset-twophaseCoeff levelset-filmCoeff misc-scalarFunctions misc-vectorFunctions misc-csgFunctions)

add_my_custom_targets(levelset)
