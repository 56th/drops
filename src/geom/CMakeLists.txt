set(HOME geom)
libs(boundary builder geomselect multigrid principallattice reftetracut simplex subtriangulation maketopo_helper)

if(MPI)
	target_link_libraries(geom-simplex parallel-parmultigrid geom-topo)
	target_link_libraries(geom-multigrid parallel-partime)
endif(MPI)

exec(maketopo geom-maketopo_helper misc-utils)
add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/topo.cpp COMMAND ${CMAKE_CURRENT_BINARY_DIR}/maketopo ${CMAKE_CURRENT_SOURCE_DIR}/maketopo_helper.cpp ${CMAKE_CURRENT_SOURCE_DIR}/topo.cpp DEPENDS maketopo)
add_library(geom-topo SHARED topo.cpp)

set_directory_properties(PROPERTY ADDITIONAL_MAKE_CLEAN_FILES topo.cpp)
