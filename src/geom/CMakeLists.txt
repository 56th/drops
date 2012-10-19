#CMakeLists.txt for geom package


#LIBRARIES
add_library(geom-boundary SHARED boundary)
add_library(geom-builder SHARED builder)
add_library(geom-geomselect SHARED geomselect)
add_library(geom-maketopo SHARED maketopo)
add_library(geom-multigrid SHARED multigrid)
add_library(geom-principallattice SHARED principallattice)
add_library(geom-reftetracut SHARED reftetracut)
add_library(geom-simplex SHARED simplex)
add_library(geom-subtriangulation SHARED subtriangulation)
add_library(geom-topo SHARED topo)
add_library(geom-maketopo_helper SHARED maketopo_helper)

target_link_libraries(geom-builder geom-boundary)
target_link_libraries(geom-geomselect geom-builder)
target_link_libraries(geom-simplex num-unknowns)
target_link_libraries(geom-multigrid geom-simplex)

#EXECUTABLES
if(NOT MPI)

        add_executable(maketopo.exe maketopo)
        target_link_libraries(maketopo.exe misc-utils geom-maketopo_helper)
    
endif(NOT MPI)
