set(HOME1 surfphasesep)
set(HOME2 separation)

#add_library(${HOME1}-ifacetransp STATIC ../surfactant/ifacetransp)
#target_link_libraries(${HOME1}-ifacetransp levelset-levelsetmapper)

exec_ser(${HOME1} geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset levelset-marking_strategy levelset-adaptriang misc-scopetimer misc-progressaccu misc-utils out-output num-discretize misc-params num-interfacePatch levelset-fastmarch surfactant-ifacetransp num-fe out-ensightOut levelset-surfacetension out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload misc-funcmap misc-scopetimer)
exec_ser(${HOME2} geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset levelset-marking_strategy levelset-adaptriang misc-scopetimer misc-progressaccu misc-utils out-output num-discretize misc-params num-interfacePatch levelset-fastmarch surfactant-ifacetransp num-fe out-ensightOut levelset-surfacetension out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload misc-funcmap misc-scopetimer num-gradient_recovery levelset-levelsetmapper)

# logger
add_library(${HOME1}-logger STATIC ${LOGGER_PATH}/SingletonLogger.cpp)
target_link_libraries(${HOME1} ${HOME1}-logger)
target_link_libraries(${HOME2} ${HOME1}-logger)

# vtk
if (VTK_FOUND)
    target_link_libraries(${HOME1} ${VTK_LIBRARIES})
    target_link_libraries(${HOME2} ${VTK_LIBRARIES})
    #vtk_module_autoinit(
    #    TARGETS ${HOME1}
    #    MODULES ${VTK_LIBRARIES}
    #)
endif()

# trilinos
if (Trilinos_FOUND)
    target_link_libraries(${HOME1}
            ${Trilinos_LIBRARIES}
            ${Trilinos_TPL_LIBRARIES}
    )
    target_link_libraries(${HOME2}
            ${Trilinos_LIBRARIES}
            ${Trilinos_TPL_LIBRARIES}
    )
endif()

# matlab for sparse mtx .mat binary output
if(Matlab_FOUND)
    add_library(${HOME1}-matlab1 SHARED IMPORTED)
    set_target_properties(${HOME1}-matlab1 PROPERTIES
            IMPORTED_LOCATION "${Matlab_ROOT_DIR}/extern/bin/glnxa64/libMatlabDataArray.so"
            INTERFACE_INCLUDE_DIRECTORIES "${Matlab_ROOT_DIR}/extern/include"
            )
    add_library(${HOME1}-matlab2 SHARED IMPORTED)
    set_target_properties(${HOME1}-matlab2 PROPERTIES
            IMPORTED_LOCATION "${Matlab_ROOT_DIR}/extern/bin/glnxa64/libMatlabEngine.so"
            INTERFACE_INCLUDE_DIRECTORIES "${Matlab_ROOT_DIR}/extern/include"
            )
    target_link_libraries(${HOME1} -lpthread ${HOME1}-matlab1 ${HOME1}-matlab2)
    target_link_libraries(${HOME2} -lpthread ${HOME1}-matlab1 ${HOME1}-matlab2)
endif()

if(NOT MPI)
add_dependencies(surfactant misc-scalarFunctions misc-vectorFunctions)
add_dependencies(${HOME1} misc-scalarFunctions misc-vectorFunctions levelset-twophaseCoeff)
add_dependencies(${HOME2} misc-scalarFunctions misc-vectorFunctions levelset-twophaseCoeff)
endif(NOT MPI)

add_my_custom_targets(${HOME1})
add_my_custom_targets(${HOME2})

add_definitions(-Wfatal-errors)
