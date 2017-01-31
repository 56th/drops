set(HOME partests)

SET(PAR_OBJ DiST-DiST DiST-mpistream DiST-remotedata parallel-parallel parallel-partime parallel-exchange parallel-parmultigrid parallel-loadbal parallel-migrateunknowns parallel-decompose)

exec_par(dist_streamtest geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer geom-deformation)

exec_par(dist_remotedatatest geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize geom-deformation num-fe num-interfacePatch misc-scopetimer)

exec_par(dist_sendtetraofbrick geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize geom-deformation num-fe num-interfacePatch misc-scopetimer)

exec_par(dist_modify geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize geom-deformation num-fe num-interfacePatch misc-scopetimer)

exec_par(dist_simplexsend geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize geom-deformation num-fe num-interfacePatch misc-scopetimer)

exec_par(dist_interface geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer geom-deformation)

exec_par(dist_ref geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize out-vtkOut parallel-parmultigrid num-fe num-interfacePatch levelset-levelset levelset-fastmarch levelset-surfacetension geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scopetimer geom-deformation)

exec_par(dist_migration geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize out-vtkOut parallel-parmultigrid num-fe num-interfacePatch levelset-levelset levelset-fastmarch levelset-surfacetension geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-funcmap stokes-instatstokes2phase stokes-slipBndOnePhase misc-scopetimer geom-deformation misc-dynamicload)
if (MPI)
  add_dependencies(dist_migration misc-scalarFunctions misc-vectorFunctions)
endif (MPI)

exec_par(dist_TestRefPar geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize out-vtkOut parallel-parmultigrid num-fe num-interfacePatch levelset-levelset levelset-fastmarch levelset-surfacetension misc-params out-output geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scopetimer geom-deformation misc-progressaccu)

exec_par(dist_ref_rising_drop geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize out-vtkOut parallel-parmultigrid num-fe num-interfacePatch levelset-levelset levelset-fastmarch levelset-surfacetension misc-params out-output geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scopetimer geom-deformation)

exec_par(dist_exchange geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer geom-deformation)

exec_par(dist_sdropsP2 geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem misc-utils num-discretize num-interfacePatch num-fe geom-reftetracut geom-principallattice misc-scopetimer geom-deformation stokes-slipBndOnePhase)

exec_par(dist_interp2 geom-simplex geom-multigrid geom-topo out-output num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils num-discretize geom-deformation)

exec_par(dist_interp1 geom-simplex geom-multigrid geom-topo out-output num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils num-discretize geom-deformation)

add_my_custom_targets(partests)
