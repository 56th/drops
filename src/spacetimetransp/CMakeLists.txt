set(HOME spacetimetransp)

libs(spacetime_error spacetime_setup spacetime_sol stxfem indicator)

add_library(spacetimetransp-sttranspCoeff SHARED sttranspCoeff)

if(NOT MPI)
  if(NOT WIN32)
      target_link_libraries(spacetimetransp-spacetime_setup num-spacetime_quad num-spacetime_map)
      target_link_libraries(spacetimetransp-spacetime_error num-spacetime_quad num-spacetime_map)
      target_link_libraries(spacetimetransp-spacetime_sol spacetimetransp-spacetime_setup spacetimetransp-stxfem)
  endif(NOT WIN32)
endif(NOT MPI)

target_link_libraries_par(spacetimetransp-spacetime_setup num-spacetime_quad num-spacetime_map)
target_link_libraries_par(spacetimetransp-spacetime_error num-spacetime_quad num-spacetime_map)
target_link_libraries_par(spacetimetransp-spacetime_sol spacetimetransp-spacetime_setup spacetimetransp-stxfem)


exec_ser(st_transp geom-boundary geom-builder geom-deformation geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-progressaccu misc-scopetimer misc-utils out-output num-discretize navstokes-instatnavstokes2phase misc-params levelset-adaptriang levelset-marking_strategy levelset-fastmarch stokes-instatstokes2phase num-MGsolver num-fe out-ensightOut num-oseenprecond poisson-transport2phase num-interfacePatch out-vtkOut surfactant-ifacetransp levelset-surfacetension geom-geomselect levelset-twophaseutils num-hypre levelset-coupling misc-funcmap geom-principallattice geom-reftetracut geom-subtriangulation num-renumber num-quadrature num-spacetime_geom num-spacetime_quad num-spacetime_map spacetimetransp-spacetime_error spacetimetransp-spacetime_setup spacetimetransp-spacetime_sol misc-dynamicload spacetimetransp-indicator num-stokespardiso stokes-slipBndOnePhase)

if(MPI)
    set(PAR_OBJ DiST-DiST DiST-mpistream DiST-remotedata parallel-parallel parallel-partime parallel-exchange parallel-parmultigrid parallel-loadbal parallel-migrateunknowns parallel-decompose)
endif(MPI)

exec_par(st_transp geom-boundary geom-builder geom-deformation geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-progressaccu misc-scopetimer misc-utils out-output num-discretize navstokes-instatnavstokes2phase misc-params levelset-adaptriang levelset-marking_strategy levelset-fastmarch stokes-instatstokes2phase num-MGsolver num-fe out-ensightOut num-oseenprecond poisson-transport2phase num-interfacePatch out-vtkOut surfactant-ifacetransp levelset-surfacetension geom-geomselect levelset-twophaseutils num-hypre levelset-coupling misc-funcmap geom-principallattice geom-reftetracut geom-subtriangulation num-renumber num-quadrature num-spacetime_geom num-spacetime_quad num-spacetime_map spacetimetransp-spacetime_error spacetimetransp-spacetime_setup spacetimetransp-spacetime_sol misc-dynamicload spacetimetransp-indicator stokes-slipBndOnePhase)

add_dependencies(st_transp misc-scalarFunctions misc-vectorFunctions misc-csgFunctions spacetimetransp-sttranspCoeff)

if(NOT MPI)
# st_transp needs Pardiso from MKL
	if (MKL_HOME STREQUAL "")
	else(MKL_HOME STREQUAL "")
		target_link_libraries_ser(st_transp -fopenmp -Wl,--start-group ${MKL_HOME}/lib/intel64/libmkl_intel_ilp64.a ${MKL_HOME}/lib/intel64/libmkl_gnu_thread.a ${MKL_HOME}/lib/intel64/libmkl_core.a -Wl,--end-group -ldl -lpthread -lm) 
	endif(MKL_HOME STREQUAL "")
endif(NOT MPI)

add_my_custom_targets(spacetimetransp)
