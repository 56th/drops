#CMakeLists.txt for levelset package


#LIBRARIES
add_library(levelset-levelset SHARED levelset)
add_library(levelset-fastmarch SHARED fastmarch)
add_library(levelset-surfacetension SHARED surfacetension)
add_library(levelset-twophaseutils SHARED twophaseutils)
add_library(levelset-coupling SHARED coupling)
add_library(levelset-twophaseCoeff SHARED twophaseCoeff)
#~ add_library(levelset-twophasedrops SHARED twophasedrops)

target_link_libraries(levelset-levelset geom-reftetracut geom-subtriangulation num-quadrature levelset-fastmarch levelset-surfacetension geom-principallattice)

#~ target_link_libraries(levelset-twophasedrops levelset-twophasedrops geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize navstokes-instatnavstokes2phase misc-params levelset-fastmarch stokes-instatstokes2phase num-MGsolver num-fe out-ensightOut num-oseenprecond poisson-transport2phase num-interfacePatch out-vtkOut surfactant-ifacetransp levelset-surfacetension geom-geomselect levelset-twophaseutils num-hypre levelset-coupling misc-vectorFunctions misc-scalarFunctions misc-funcmap levelset-twophaseCoeff geom-principallattice geom-reftetracut geom-subtriangulation num-renumber num-quadrature)

#~ if(NOT MPI)
	#~ add_library(levelset-film SHARED film)
	add_library(levelset-filmCoeff SHARED filmCoeff)
	#~ add_library(levelset-brick_transp SHARED brick_transp)
	#~ add_library(levelset-lsshear SHARED lsshear)
	#~ add_library(levelset-prJump SHARED prJump)
	#~ add_library(levelset-reparam SHARED reparam)
	#~ add_library(levelset-surfTens SHARED surfTens)
		
	#~ target_link_libraries(levelset-film levelset-film geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize navstokes-instatnavstokes2phase misc-params levelset-fastmarch stokes-instatstokes2phase num-MGsolver num-fe out-ensightOut out-vtkOut num-oseenprecond num-interfacePatch levelset-surfacetension levelset-coupling misc-vectorFunctions misc-funcmap levelset-filmCoeff geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scalarFunctions num-renumber)
	#~ target_link_libraries(levelset-brick_transp levelset-brick_transp geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize misc-params levelset-fastmarch num-fe poisson-transport2phase out-ensightOut out-vtkOut num-interfacePatch levelset-surfacetension stokes-instatstokes2phase misc-vectorFunctions misc-funcmap levelset-twophaseCoeff geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scalarFunctions)
	#~ target_link_libraries(levelset-lsshear levelset-lsshear geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize levelset-fastmarch num-fe out-ensightOut stokes-instatstokes2phase navstokes-instatnavstokes2phase num-MGsolver num-interfacePatch levelset-surfacetension levelset-coupling misc-funcmap geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scalarFunctions num-renumber)
	#~ target_link_libraries(levelset-prJump levelset-prJump geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize misc-params levelset-fastmarch stokes-instatstokes2phase num-MGsolver num-fe out-ensightOut num-oseenprecond num-interfacePatch levelset-surfacetension levelset-coupling misc-vectorFunctions misc-funcmap geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scalarFunctions num-renumber)
	#~ target_link_libraries(levelset-reparam levelset-reparam levelset-fastmarch levelset-levelset geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe out-ensightOut num-interfacePatch levelset-surfacetension geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)
	#~ target_link_libraries(levelset-surfTens levelset-surfTens geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize misc-params levelset-fastmarch stokes-instatstokes2phase navstokes-instatnavstokes2phase num-fe out-ensightOut num-oseenprecond num-MGsolver misc-funcmap num-interfacePatch levelset-surfacetension levelset-coupling misc-vectorFunctions geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scalarFunctions num-renumber)
#~ endif(NOT MPI)

#EXECUTABLES
add_executable(twophasedrops.exe twophasedrops)
target_link_libraries(twophasedrops.exe geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize navstokes-instatnavstokes2phase misc-params levelset-fastmarch stokes-instatstokes2phase num-MGsolver num-fe out-ensightOut num-oseenprecond poisson-transport2phase num-interfacePatch out-vtkOut surfactant-ifacetransp levelset-surfacetension geom-geomselect levelset-twophaseutils num-hypre levelset-coupling misc-vectorFunctions misc-scalarFunctions misc-funcmap levelset-twophaseCoeff geom-principallattice geom-reftetracut geom-subtriangulation num-renumber num-quadrature)
if(MPI)
    target_link_libraries(twophasedrops.exe DiST-DiST DiST-mpistream DiST-remotedata parallel-parallel parallel-partime parallel-exchange parallel-parmultigrid parallel-loadbal parallel-migrateunknowns parallel-decompose)
endif(MPI)
	
if(NOT MPI)
	add_executable(reparam.exe reparam)
	add_executable(lsshear.exe lsshear)
	add_executable(surfTens.exe surfTens)
	add_executable(prJump.exe prJump)
	add_executable(film.exe film)
	add_executable(brick_transp.exe brick_transp)

	target_link_libraries(reparam.exe levelset-fastmarch levelset-levelset geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe out-ensightOut num-interfacePatch levelset-surfacetension geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)
	target_link_libraries(lsshear.exe geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize levelset-fastmarch num-fe out-ensightOut stokes-instatstokes2phase navstokes-instatnavstokes2phase num-MGsolver num-interfacePatch levelset-surfacetension levelset-coupling misc-funcmap geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scalarFunctions num-renumber)
	target_link_libraries(surfTens.exe geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize misc-params levelset-fastmarch stokes-instatstokes2phase navstokes-instatnavstokes2phase num-fe out-ensightOut num-oseenprecond num-MGsolver misc-funcmap num-interfacePatch levelset-surfacetension levelset-coupling misc-vectorFunctions geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scalarFunctions num-renumber)
	target_link_libraries(prJump.exe geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize misc-params levelset-fastmarch stokes-instatstokes2phase num-MGsolver num-fe out-ensightOut num-oseenprecond num-interfacePatch levelset-surfacetension levelset-coupling misc-vectorFunctions misc-funcmap geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scalarFunctions num-renumber)
	target_link_libraries(film.exe geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize navstokes-instatnavstokes2phase misc-params levelset-fastmarch stokes-instatstokes2phase num-MGsolver num-fe out-ensightOut out-vtkOut num-oseenprecond num-interfacePatch levelset-surfacetension levelset-coupling misc-vectorFunctions misc-funcmap levelset-filmCoeff geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scalarFunctions num-renumber)
	target_link_libraries(brick_transp.exe geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize misc-params levelset-fastmarch num-fe poisson-transport2phase out-ensightOut out-vtkOut num-interfacePatch levelset-surfacetension stokes-instatstokes2phase misc-vectorFunctions misc-funcmap levelset-twophaseCoeff geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-scalarFunctions)
endif(NOT MPI)
