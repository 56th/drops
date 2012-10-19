#CMakeLists.txt for navstokes package


#LIBRARIES
add_library(navstokes-instatnavstokes2phase SHARED instatnavstokes2phase)
target_link_libraries(navstokes-instatnavstokes2phase num-renumber)

#EXECUTABLES
if(NOT MPI)
	add_library(navstokes-nsdrops SHARED nsdrops)
	add_executable(nsdrops.exe nsdrops)
	target_link_libraries(nsdrops.exe geom-boundary geom-builder navstokes-nsdrops geom-simplex geom-multigrid num-unknowns geom-topo num-fe num-interfacePatch misc-problem misc-utils out-output num-discretize geom-principallattice geom-reftetracut)
	
	add_library(navstokes-nsdrops_begehung SHARED nsdrops_begehung)
	add_executable(nsdrops_begehung.exe nsdrops_begehung)
	target_link_libraries(nsdrops_begehung.exe geom-boundary geom-builder navstokes-nsdrops_begehung geom-simplex geom-multigrid num-unknowns geom-topo num-fe num-interfacePatch misc-problem misc-utils out-output num-discretize geom-principallattice geom-reftetracut)

	add_library(navstokes-insdrops SHARED insdrops)
	add_executable(insdrops.exe insdrops)
	target_link_libraries(insdrops.exe geom-boundary geom-builder navstokes-insdrops geom-simplex geom-multigrid num-unknowns geom-topo num-fe num-interfacePatch misc-problem misc-utils out-output num-discretize out-ensightOut geom-principallattice geom-reftetracut)

	add_library(navstokes-insadrops SHARED insadrops)
	add_executable(insadrops.exe insadrops)
	target_link_libraries(insadrops.exe geom-boundary geom-builder navstokes-insadrops geom-simplex geom-multigrid num-unknowns geom-topo num-fe num-interfacePatch misc-problem misc-utils out-output num-discretize out-ensightOut geom-principallattice geom-reftetracut)

endif(NOT MPI)
