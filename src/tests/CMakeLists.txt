set(HOME tests)

exec_ser(reftest geom-boundary geom-builder geom-deformation geom-simplex geom-multigrid num-unknowns misc-utils geom-topo out-output misc-problem num-interfacePatch num-fe)

exec_ser(interp1 geom-simplex geom-deformation geom-multigrid geom-topo out-output num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils)

#exec_ser(interfacepatch_doublecut misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils)

exec_ser(interfacepatch_doublecut
    misc-utils geom-builder geom-simplex
    geom-multigrid geom-boundary geom-topo num-unknowns misc-problem
    num-fe num-discretize num-interfacePatch geom-deformation)

exec_ser(interp2 geom-simplex geom-multigrid geom-deformation geom-multigrid geom-topo out-output num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils num-discretize)

exec_ser(prolongationp2test geom-deformation geom-simplex geom-multigrid geom-topo num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils)

exec_ser(tetrabuildertest geom-deformation geom-simplex geom-multigrid geom-topo num-unknowns geom-builder geom-boundary misc-utils misc-problem num-interfacePatch num-fe)

exec_ser(mattest misc-utils)

exec_ser(testfe misc-utils geom-topo geom-boundary geom-deformation misc-problem num-interfacePatch num-fe geom-simplex geom-multigrid num-unknowns geom-builder)

exec_ser(sbuffer misc-utils)

exec_ser(minres misc-utils)

exec_ser(meshreader geom-deformation geom-simplex geom-multigrid geom-topo num-unknowns geom-builder geom-boundary misc-utils out-output misc-problem num-interfacePatch num-fe misc-params)

exec_ser(restrictp2 geom-deformation geom-simplex geom-multigrid geom-topo out-output num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils)

exec_ser(vectest misc-utils)

exec_ser(spacetime_decomp misc-utils geom-deformation geom-simplex geom-boundary geom-topo num-spacetime_geom num-spacetime_map num-spacetime_quad misc-scopetimer geom-multigrid geom-builder geom-principallattice geom-topo)

exec_ser(spacetime_surf geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-spacetime_geom num-spacetime_map num-spacetime_quad misc-scopetimer misc-progressaccu num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize num-interfacePatch misc-params levelset-fastmarch stokes-instatstokes2phase levelset-surfacetension misc-funcmap misc-vectorFunctions geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload)

exec_ser(p2local misc-utils geom-deformation geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns out-output misc-problem geom-builder num-interfacePatch num-discretize num-fe)

exec_ser(quadbase misc-utils geom-builder geom-deformation geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem num-fe num-discretize num-interfacePatch)

exec_ser(globallist misc-utils)

exec_ser(serialization misc-utils misc-params geom-builder geom-simplex geom-multigrid geom-deformation geom-boundary geom-topo num-unknowns out-output num-fe misc-problem num-interfacePatch)

exec_ser(triang misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns geom-deformation misc-problem num-interfacePatch num-fe)

exec_ser(combiner misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo levelset-adaptriang levelset-marking_strategy out-output out-vtkOut)

exec_ser(quadCut misc-utils geom-builder geom-deformation geom-simplex geom-multigrid misc-scopetimer misc-progressaccu geom-boundary geom-topo num-unknowns misc-problem num-interfacePatch levelset-levelset levelset-fastmarch num-discretize num-fe levelset-surfacetension geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

exec_ser(bicgstab misc-utils)

if (MKL_HOME STREQUAL "")
else (MKL_HOME STREQUAL "")
  exec_ser(pardiso misc-utils)
  target_link_libraries_ser(pardiso -fopenmp -Wl,--start-group ${MKL_HOME}/lib/intel64/libmkl_intel_ilp64.a ${MKL_HOME}/lib/intel64/libmkl_gnu_thread.a ${MKL_HOME}/lib/intel64/libmkl_core.a -Wl,--end-group -ldl -lpthread -lm) 
endif (MKL_HOME STREQUAL "")


exec_ser(gcr misc-utils)

exec_ser(blockmat misc-utils)

exec_ser(quad5 misc-utils geom-builder geom-deformation geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem num-fe num-discretize num-interfacePatch geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

exec_ser(downwind misc-utils geom-multigrid geom-simplex geom-topo num-unknowns num-fe geom-deformation num-interfacePatch misc-problem)

exec_ser(quad5_2D misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem geom-deformation num-fe num-discretize num-interfacePatch)

exec_ser(interfaceP1FE misc-utils geom-builder geom-simplex geom-multigrid geom-deformation misc-scopetimer misc-progressaccu geom-boundary geom-topo num-unknowns misc-problem num-interfacePatch levelset-levelset levelset-fastmarch num-discretize num-fe levelset-surfacetension geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

#set(xfem-staticlibs ${BEGIN_STATIC_LIBS} misc-vectorFunctions ${END_STATIC_LIBS})
exec_ser(xfem misc-utils geom-builder geom-simplex geom-multigrid geom-deformation misc-scopetimer misc-progressaccu geom-boundary geom-topo num-unknowns misc-problem num-interfacePatch num-fe num-discretize levelset-levelset levelset-fastmarch stokes-instatstokes2phase out-output out-ensightOut levelset-surfacetension misc-funcmap misc-vectorFunctions geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload stokes-slipBndOnePhase)

exec_ser(f_Gamma geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation misc-scopetimer misc-progressaccu num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize num-interfacePatch misc-params levelset-fastmarch stokes-instatstokes2phase levelset-surfacetension misc-funcmap misc-vectorFunctions geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload stokes-slipBndOnePhase)
if(NOT MPI)
    add_dependencies(f_Gamma levelset-twophaseCoeff misc-scalarFunctions misc-vectorFunctions misc-csgFunctions)
    add_dependencies(xfem levelset-twophaseCoeff misc-scalarFunctions misc-vectorFunctions)
endif(NOT MPI)

exec_ser(neq misc-utils)

exec_ser(extendP1onChild num-discretize misc-utils geom-topo num-fe misc-problem geom-deformation geom-simplex geom-multigrid num-unknowns num-interfacePatch)

exec_ser(splitboundary geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns misc-utils geom-topo out-output num-fe misc-problem num-interfacePatch)

exec_ser(reparam_init geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation misc-scopetimer misc-progressaccu num-unknowns misc-utils geom-topo out-output misc-problem levelset-levelset levelset-adaptriang levelset-marking_strategy num-fe misc-utils misc-params num-discretize num-interfacePatch levelset-fastmarch levelset-surfacetension out-ensightOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

#reparam renamed to reparamtest
exec_ser(reparamtest levelset-fastmarch levelset-levelset levelset-adaptriang levelset-marking_strategy geom-deformation misc-scopetimer misc-progressaccu geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe out-ensightOut num-interfacePatch levelset-surfacetension misc-params out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload) ##${PAR_OBJ}
if(NOT MPI)
    add_dependencies(reparamtest levelset-twophaseCoeff misc-scalarFunctions misc-vectorFunctions)
endif(NOT MPI)

exec_ser(principallattice misc-utils geom-principallattice num-discretize geom-deformation geom-topo num-fe num-interfacePatch misc-problem num-unknowns geom-simplex geom-multigrid geom-builder geom-boundary num-unknowns geom-reftetracut geom-subtriangulation num-quadrature)

exec_ser(bndTrianglePartition misc-utils geom-principallattice num-discretize geom-deformation geom-topo num-fe num-interfacePatch misc-problem num-unknowns geom-simplex geom-multigrid geom-builder geom-boundary num-unknowns geom-reftetracut geom-subtriangulation num-quadrature geom-maketopo_helper)

exec_ser(quad_extra misc-utils geom-builder geom-simplex geom-multigrid geom-deformation geom-boundary geom-topo num-unknowns misc-problem num-fe num-discretize num-interfacePatch geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

exec_ser(csgtest levelset-levelset levelset-fastmarch geom-deformation misc-funcmap misc-scopetimer misc-progressaccu geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe num-interfacePatch levelset-surfacetension misc-params out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature geom-csg)

exec_ser(dynload geom-simplex misc-funcmap misc-dynamicload)

exec_ser(show_domain geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation misc-scopetimer num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize num-interfacePatch misc-params levelset-fastmarch stokes-instatstokes2phase levelset-surfacetension misc-funcmap misc-vectorFunctions geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload out-vtkOut)
if(NOT MPI)
    add_dependencies(show_domain levelset-twophaseCoeff misc-scalarFunctions misc-vectorFunctions misc-csgFunctions)
endif(NOT MPI)

if(NOT MPI)
    add_dependencies(dynload misc-scalarFunctions misc-vectorFunctions)

    if(EXISTS ${SUITESPARSE_HOME}) 
        exec_ser(directsolver misc-utils)
        #${SUITESPARSE_HOME}/metis-4.0/libmetis.a not used for now
        find_package(BLAS)
        if (BLAS_FOUND)
            target_link_libraries_ser(directsolver ${SUITESPARSE_HOME}/UMFPACK/Lib/libumfpack.a ${SUITESPARSE_HOME}/CHOLMOD/Lib/libcholmod.a ${SUITESPARSE_HOME}/AMD/Lib/libamd.a ${SUITESPARSE_HOME}/COLAMD/Lib/libcolamd.a ${SUITESPARSE_HOME}/CCOLAMD/Lib/libccolamd.a ${SUITESPARSE_HOME}/CAMD/Lib/libcamd.a ${SUITESPARSE_HOME}/SuiteSparse_config/libsuitesparseconfig.a ${BLAS_LIBRARIES} lapack)
            set_target_properties(directsolver PROPERTIES COMPILE_FLAGS "-I${SUITESPARSE_HOME}/CHOLMOD/Include -I${SUITESPARSE_HOME}/SuiteSparse_config -I${SUITESPARSE_HOME}/UMFPACK/Include -I${SUITESPARSE_HOME}/AMD/Include")
        else(BLAS_FOUND)
            add_custom_target(GiveBlasWarning COMMENT "Warning: BLAS library not found, directsolver will not link properly.")
            target_link_libraries(directsolver GiveBlasWarning)
        endif(BLAS_FOUND)    
    else()
        add_custom_target(directsolver COMMENT "Warning: SuiteSparse library not found, skipped compilation of directsolver.")
    endif()
endif(NOT MPI)


#### BEGIN HACKY VERSION OF TEST SUITE ###

add_subdirectory(scripts)
add_subdirectory(specifications)

#### END HACKY VERSION OF TEST SUITE ###

add_my_custom_targets(tests)
