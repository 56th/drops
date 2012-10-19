#CMakeLists.txt for surfactant package


#LIBRARIES
add_library(surfactant-ifacetransp SHARED ifacetransp)

#EXECUTABLES
if(NOT MPI)
	add_library(surfactant-surfactant SHARED surfactant)
	add_executable(surfactant.exe surfactant)
	target_link_libraries(surfactant.exe surfactant-surfactant geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize misc-params num-interfacePatch levelset-fastmarch surfactant-ifacetransp num-fe out-ensightOut levelset-surfacetension out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

	add_library(surfactant-surfacenorms SHARED surfacenorms)
	add_executable(surfacenorms.exe surfacenorms)
	target_link_libraries(surfacenorms.exe surfactant-surfacenorms geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize misc-params num-interfacePatch levelset-fastmarch surfactant-ifacetransp num-fe out-ensightOut levelset-surfacetension out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

endif(NOT MPI)
