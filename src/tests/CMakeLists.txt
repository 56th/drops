set(HOME tests)

exec_ser(reftest tests-reftest geom-boundary geom-builder geom-deformation geom-simplex geom-multigrid num-unknowns misc-utils geom-topo out-output misc-problem num-interfacePatch num-fe)

exec_ser(interp1 geom-simplex geom-deformation geom-multigrid geom-topo out-output num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils)

#exec_ser(interfacepatch_doublecut misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils)

exec_ser(interfacepatch_doublecut
    misc-utils geom-builder geom-simplex
    geom-multigrid geom-boundary geom-topo num-unknowns misc-problem
    num-fe num-discretize num-interfacePatch geom-deformation)

exec_ser(interp2 geom-simplex geom-multigrid geom-deformation geom-multigrid geom-topo out-output num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils num-discretize)

exec_ser(prolongationp2test tests-prolongationp2test geom-deformation geom-simplex geom-multigrid geom-topo num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils)

exec_ser(tetrabuildertest geom-deformation geom-simplex geom-multigrid geom-topo num-unknowns geom-builder geom-boundary misc-utils misc-problem num-interfacePatch num-fe)

exec_ser(mattest misc-utils)

exec_ser(testfe tests-testfe misc-utils geom-topo geom-boundary geom-deformation misc-problem num-interfacePatch num-fe geom-simplex geom-multigrid num-unknowns geom-builder)

exec_ser(sbuffer tests-sbuffer misc-utils)

exec_ser(minres tests-minres misc-utils)

exec_ser(meshreader geom-deformation geom-simplex geom-multigrid geom-topo num-unknowns geom-builder geom-boundary misc-utils out-output misc-problem num-interfacePatch num-fe)

exec_ser(restrictp2 geom-deformation geom-simplex geom-multigrid geom-topo out-output num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils)

exec_ser(vectest tests-vectest misc-utils)

exec_ser(p2local misc-utils geom-deformation geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns out-output misc-problem geom-builder num-interfacePatch num-discretize num-fe)

exec_ser(quadbase tests-quadbase misc-utils geom-builder geom-deformation geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem num-fe num-discretize num-interfacePatch) 

exec_ser(globallist tests-globallist misc-utils)

exec_ser(serialization misc-utils misc-utils geom-builder geom-simplex geom-multigrid geom-deformation geom-boundary geom-topo num-unknowns out-output num-fe misc-problem num-interfacePatch)

exec_ser(triang tests-triang misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns geom-deformation misc-problem num-interfacePatch num-fe)

exec_ser(quadCut tests-quadCut misc-utils geom-builder geom-deformation geom-simplex geom-multigrid misc-scopetimer misc-progressaccu geom-boundary geom-topo num-unknowns misc-problem num-interfacePatch levelset-levelset levelset-fastmarch num-discretize num-fe levelset-surfacetension geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

exec_ser(bicgstab tests-bicgstab misc-utils)

exec_ser(gcr tests-gcr misc-utils)

exec_ser(blockmat tests-blockmat misc-utils)

exec_ser(mass tests-mass misc-utils)

exec_ser(quad5 misc-utils geom-builder geom-deformation geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem num-fe num-discretize num-interfacePatch geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

exec_ser(downwind misc-utils geom-multigrid geom-simplex geom-topo num-unknowns num-fe geom-deformation num-interfacePatch misc-problem)

exec_ser(quad5_2D tests-quad5_2D misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem geom-deformation num-fe num-discretize num-interfacePatch)

exec_ser(interfaceP1FE misc-utils geom-builder geom-simplex geom-multigrid geom-deformation misc-scopetimer misc-progressaccu geom-boundary geom-topo num-unknowns misc-problem num-interfacePatch levelset-levelset levelset-fastmarch num-discretize num-fe levelset-surfacetension geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

set(xfem-staticlibs ${BEGIN_STATIC_LIBS} misc-vectorFunctions ${END_STATIC_LIBS})
exec_ser(xfem misc-utils geom-builder geom-simplex geom-multigrid geom-deformation misc-scopetimer misc-progressaccu geom-boundary geom-topo num-unknowns misc-problem num-interfacePatch num-fe num-discretize levelset-levelset levelset-fastmarch stokes-instatstokes2phase out-ensightOut levelset-surfacetension misc-funcmap misc-vectorFunctions geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

exec_ser(f_Gamma geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation misc-scopetimer misc-progressaccu num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize num-interfacePatch misc-params levelset-fastmarch stokes-instatstokes2phase levelset-surfacetension misc-funcmap misc-vectorFunctions geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload)

exec_ser(neq misc-utils) 

exec_ser(extendP1onChild tests-extendP1onChild num-discretize misc-utils geom-topo num-fe misc-problem geom-deformation geom-simplex geom-multigrid num-unknowns num-interfacePatch)

exec_ser(splitboundary geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns misc-utils geom-topo out-output num-fe misc-problem num-interfacePatch)

exec_ser(reparam_init geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation misc-scopetimer misc-progressaccu num-unknowns misc-utils geom-topo out-output misc-problem levelset-levelset num-fe misc-utils num-discretize num-interfacePatch levelset-fastmarch levelset-surfacetension out-ensightOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

#reparam renamed to reparamtest
exec_ser(reparamtest levelset-fastmarch levelset-levelset geom-deformation misc-scopetimer misc-progressaccu geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe out-ensightOut num-interfacePatch levelset-surfacetension misc-params out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload) ##${PAR_OBJ}

exec_ser(principallattice tests-principallattice misc-utils geom-principallattice num-discretize geom-deformation geom-topo num-fe num-interfacePatch misc-problem num-unknowns geom-simplex geom-multigrid geom-builder geom-boundary num-unknowns geom-reftetracut geom-subtriangulation num-quadrature)

#exec_ser(quad_extra tests-quad_extra misc-utils geom-builder geom-simplex geom-multigrid geom-deformation geom-boundary geom-topo num-unknowns misc-problem num-fe num-discretize num-interfacePatch geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

exec_ser(csgtest levelset-levelset levelset-fastmarch geom-deformation misc-funcmap misc-scopetimer misc-progressaccu geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch levelset-surfacetension misc-params out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature geom-csg)

exec(dynload geom-simplex misc-funcmap misc-dynamicload)

add_dependencies(dynload misc-scalarFunctions misc-vectorFunctions)

if(EXISTS ${SUITSPARSE_HOME})
	lib_ser(directsolver)
	target_link_libraries_ser(tests-directsolver ${SUITESPARSE_HOME}/UMFPACK/Lib/libumfpack.a ${SUITESPARSE_HOME}/CHOLMOD/Lib/libcholmod.a ${SUITESPARSE_HOME}/AMD/Lib/libamd.a ${SUITESPARSE_HOME}/COLAMD/Lib/libcolamd.a ${SUITESPARSE_HOME}/CCOLAMD/Lib/libccolamd.a ${SUITESPARSE_HOME}/CAMD/Lib/libcamd.a ${SUITESPARSE_HOME}/metis-4.0/libmetis.a blas lapack)
	set_target_properties(tests-directsolver PROPERTIES COMPILE_FLAGS -I${SUITESPARSE_HOME}/CHOLMOD/Include -I${SUITESPARSE_HOME}/UFconfig -I${SUITESPARSE_HOME}/UMFPACK/Include -I${SUITESPARSE_HOME}/AMD/Include)
endif()

############### TEST FAILS TO COMPILE ####################### 
#exec_ser(quad_extra misc-utils geom-builder geom-simplex geom-multigrid 
# geom-boundary geom-topo num-unknowns misc-problem 
# num-fe num-discretize num-interfacePatch
# geom-principallattice geom-reftetracut geom-subtriangulation num-#quadrature)

#### BEGIN HACKY VERSION OF TEST SUITE ###

add_subdirectory(scripts)
add_subdirectory(specifications)

#### END HACKY VERSION OF TEST SUITE ###
