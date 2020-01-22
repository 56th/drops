set(HOME surfphasesep)

add_library(${HOME}-ifacetransp STATIC ../surfactant/ifacetransp)

target_link_libraries(${HOME}-ifacetransp levelset-levelsetmapper)

exec_ser(surfphasesep geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset levelset-marking_strategy levelset-adaptriang misc-scopetimer misc-progressaccu misc-utils out-output num-discretize misc-params num-interfacePatch levelset-fastmarch surfactant-ifacetransp num-fe out-ensightOut levelset-surfacetension out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload misc-funcmap misc-scopetimer)

exec_ser(separation geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset levelset-marking_strategy levelset-adaptriang misc-scopetimer misc-progressaccu misc-utils out-output num-discretize misc-params num-interfacePatch levelset-fastmarch surfactant-ifacetransp num-fe out-ensightOut levelset-surfacetension out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload misc-funcmap misc-scopetimer num-gradient_recovery levelset-levelsetmapper)

# logger
add_library(${HOME}-logger STATIC ${LOGGER_PATH}/SingletonLogger.cpp)
target_link_libraries(${HOME} ${HOME}-logger)

# vtk
if (VTK_FOUND)
    target_link_libraries(${HOME} ${VTK_LIBRARIES})
    #vtk_module_autoinit(
    #    TARGETS ${HOME}
    #    MODULES ${VTK_LIBRARIES}
    #)
endif()

# trilinos
if (Trilinos_FOUND)
    target_link_libraries(${HOME}
            ${Trilinos_LIBRARIES}
            ${Trilinos_TPL_LIBRARIES}
            )
endif()

# matlab for sparse mtx .mat binary output
if(Matlab_FOUND)
    add_library(${HOME}-matlab1 SHARED IMPORTED)
    set_target_properties(${HOME}-matlab1 PROPERTIES
            IMPORTED_LOCATION "${Matlab_ROOT_DIR}/extern/bin/glnxa64/libMatlabDataArray.so"
            INTERFACE_INCLUDE_DIRECTORIES "${Matlab_ROOT_DIR}/extern/include"
            )
    add_library(${HOME}-matlab2 SHARED IMPORTED)
    set_target_properties(${HOME}-matlab2 PROPERTIES
            IMPORTED_LOCATION "${Matlab_ROOT_DIR}/extern/bin/glnxa64/libMatlabEngine.so"
            INTERFACE_INCLUDE_DIRECTORIES "${Matlab_ROOT_DIR}/extern/include"
            )
    target_link_libraries(${HOME} -lpthread ${HOME}-matlab1 ${HOME}-matlab2)
endif()

if(NOT MPI)
add_dependencies(surfactant misc-scalarFunctions misc-vectorFunctions)
add_dependencies(surfphasesep misc-scalarFunctions misc-vectorFunctions levelset-twophaseCoeff)
add_dependencies(separation misc-scalarFunctions misc-vectorFunctions levelset-twophaseCoeff)
endif(NOT MPI)

add_my_custom_targets(surfphasesep)
add_my_custom_targets(separation)

add_definitions(-Wfatal-errors)
