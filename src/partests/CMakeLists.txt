#CMakeLists.txt for partests package

if(MPI)

#LIBRARIES
add_library(partests-dist_exchange SHARED dist_exchange)
add_library(partests-dist_interface SHARED dist_interface)
add_library(partests-dist_interp1 SHARED dist_interp1)
add_library(partests-dist_interp2 SHARED dist_interp2)
add_library(partests-dist_migration SHARED dist_migration)
add_library(partests-dist_modify SHARED dist_modify)
add_library(partests-dist_ref SHARED dist_ref)
add_library(partests-dist_ref_rising_drop SHARED dist_ref_rising_drop)
add_library(partests-dist_remotedatatest SHARED dist_remotedatatest)
add_library(partests-dist_sdropsP2 SHARED dist_sdropsP2)
add_library(partests-dist_sendtetraofbrick SHARED dist_sendtetraofbrick)
add_library(partests-dist_simplexsend SHARED dist_simplexsend)
add_library(partests-dist_streamtest SHARED dist_streamtest)
add_library(partests-dist_TestRefPar SHARED dist_TestRefPar)
add_library(partests-TestPoissonPar SHARED TestPoissonPar)

add_library(parallel ../DiST/DiST ../DiST/mpistream ../DiST/remotedata
           ../parallel/parallel ../parallel/partime
           ../parallel/exchange ../parallel/parmultigrid
           ../parallel/loadbal ../parallel/migrateunknowns
           ../parallel/decompose)

#EXECUTABLES
add_executable(dist_streamtest.exe dist_streamtest)
target_link_libraries(dist_streamtest.exe partests-dist_streamtest geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer)

add_executable(dist_remotedatatest.exe dist_remotedatatest)
target_link_libraries(dist_remotedatatest.exe partests-dist_remotedatatest geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer)

add_executable(dist_sendtetraofbrick.exe dist_sendtetraofbrick)
target_link_libraries(dist_sendtetraofbrick.exe partests-dist_sendtetraofbrick geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer)

add_executable(dist_modify.exe dist_modify)
target_link_libraries(dist_modify.exe partests-dist_modify geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer)

add_executable(dist_simplexsend.exe dist_simplexsend)
target_link_libraries(dist_simplexsend.exe partests-dist_simplexsend geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer)

add_executable(dist_interface.exe dist_interface)
target_link_libraries(dist_interface.exe partests-dist_interface geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer)

add_executable(dist_ref.exe dist_ref)
target_link_libraries(dist_ref.exe partests-dist_ref geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize out-vtkOut parallel-parmultigrid num-fe num-interfacePatch levelset-levelset levelset-fastmarch levelset-surfacetension geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scopetimer)

add_executable(dist_migration.exe dist_migration)
target_link_libraries(dist_migration.exe partests-dist_migration geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize out-vtkOut parallel-parmultigrid num-fe num-interfacePatch levelset-levelset levelset-fastmarch levelset-surfacetension geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-funcmap misc-vectorFunctions misc-scalarFunctions stokes-instatstokes2phase misc-scopetimer)

add_executable(dist_TestRefPar.exe dist_TestRefPar)
target_link_libraries(dist_TestRefPar.exe partests-dist_TestRefPar geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize out-vtkOut parallel-parmultigrid num-fe num-interfacePatch levelset-levelset levelset-fastmarch levelset-surfacetension misc-params out-output geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scopetimer)

add_executable(dist_ref_rising_drop.exe dist_ref_rising_drop)
target_link_libraries(dist_ref_rising_drop.exe partests-dist_ref_rising_drop geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize out-vtkOut parallel-parmultigrid num-fe num-interfacePatch levelset-levelset levelset-fastmarch levelset-surfacetension misc-params out-output geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scopetimer)

add_executable(dist_exchange.exe dist_exchange)
target_link_libraries(dist_exchange.exe partests-dist_exchange geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch misc-scopetimer)

add_executable(dist_sdropsP2.exe dist_sdropsP2)
target_link_libraries(dist_sdropsP2.exe partests-dist_sdropsP2 geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem misc-utils num-discretize num-interfacePatch num-fe geom-reftetracut geom-principallattice misc-scopetimer)

add_executable(dist_interp2.exe dist_interp2)
target_link_libraries(dist_interp2.exe partests-dist_interp2 geom-simplex geom-multigrid geom-topo out-output num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils num-discretize)

add_executable(dist_interp1.exe dist_interp1)
target_link_libraries(dist_interp1.exe partests-dist_interp1 geom-simplex geom-multigrid geom-topo out-output num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils num-discretize)

endif(MPI)