set(HOME geom)
libs(boundary builder geomselect multigrid principallattice reftetracut simplex subtriangulation maketopo_helper deformation csg)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/topo.cpp) #otherwise it is impossible to build topo
    libs(topo)
else()
    message(STATUS "topo.cpp does not exist. Please build it using `make topo` on a serial build. Building anything else before may fail.")
endif()

target_link_libraries_par(geom-simplex parallel-parmultigrid geom-topo)
target_link_libraries_par(geom-multigrid parallel-partime)

target_link_libraries(geom-builder geom-boundary misc-utils misc-params)
target_link_libraries(geom-simplex geom-deformation num-unknowns)
target_link_libraries(geom-geomselect geom-builder)
target_link_libraries(geom-multigrid geom-simplex num-bndData misc-params)
target_link_libraries(geom-deformation misc-problem)

if(NOT MPI)
    exec(maketopo geom-maketopo_helper misc-utils)
    set_target_properties(maketopo PROPERTIES EXCLUDE_FROM_ALL 1)

    add_custom_target(topo COMMAND ${CMAKE_CURRENT_BINARY_DIR}/maketopo ${CMAKE_CURRENT_SOURCE_DIR}/maketopo_helper.cpp ${CMAKE_CURRENT_SOURCE_DIR}/topo.cpp DEPENDS maketopo)
endif()

add_my_custom_targets(geom)
