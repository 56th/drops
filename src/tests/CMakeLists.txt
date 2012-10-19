#CMakeLists.txt for tests package


#LIBRARIES
if(NOT MPI)
	add_library(tests-bicgstab SHARED bicgstab)
	add_library(tests-blockmat SHARED blockmat)
#	add_library(tests-directsolver SHARED directsolver)
	add_library(tests-downwind SHARED downwind)
	add_library(tests-extendP1onChild SHARED extendP1onChild)
	add_library(tests-f_Gamma SHARED f_Gamma)
	add_library(tests-gcr SHARED gcr)
	add_library(tests-globallist SHARED globallist)
	add_library(tests-interfaceP1FE SHARED interfaceP1FE)
	add_library(tests-interfacepatch_doublecut SHARED interfacepatch_doublecut)
	add_library(tests-interp1 SHARED interp1)
	add_library(tests-interp2 SHARED interp2)
	add_library(tests-mass SHARED mass)
	add_library(tests-mattest SHARED mattest)
	add_library(tests-meshreader SHARED meshreader)
	add_library(tests-minres SHARED minres)
	add_library(tests-neq SHARED neq)
	add_library(tests-p2local SHARED p2local)
	add_library(tests-principallattice SHARED principallattice)
	add_library(tests-prolongationp2test SHARED prolongationp2test)
	add_library(tests-quad5_2D SHARED quad5_2D)
	add_library(tests-quad5 SHARED quad5)
	add_library(tests-quadbase SHARED quadbase)
	add_library(tests-quadCut SHARED quadCut)
	add_library(tests-quad_extra SHARED quad_extra)
	add_library(tests-reftest SHARED reftest)
	add_library(tests-reparam SHARED reparam)
	add_library(tests-reparam_init SHARED reparam_init)
	add_library(tests-restrictp2 SHARED restrictp2)
	add_library(tests-sbuffer SHARED sbuffer)
	add_library(tests-serialization SHARED serialization)
	add_library(tests-splitboundary SHARED splitboundary)
	add_library(tests-testfe SHARED testfe)
	add_library(tests-tetrabuildertest SHARED tetrabuildertest)
	add_library(tests-triang SHARED triang)
	add_library(tests-vectest SHARED vectest)
	add_library(tests-xfem SHARED xfem)
endif(NOT MPI)

#EXECUTABLES
if(NOT MPI)

        add_executable(reftest.exe reftest)
        target_link_libraries(reftest.exe tests-reftest geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns misc-utils geom-topo out-output misc-problem num-interfacePatch num-fe )
        add_test(reftest ${CMAKE_CURRENT_BINARY_DIR}/reftest.exe)

if(FAIL)
        add_executable(ip1test.exe interp1)
        target_link_libraries(ip1test.exe tests-interp1 geom-simplex geom-multigrid geom-topo out-output num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils )
        add_test(interp1 ${CMAKE_CURRENT_BINARY_DIR}/interp1.exe)
endif(FAIL)

#        add_executable(doublecut.exe doublecut)
#        target_link_libraries(doublecut.exe misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem num-fe num-discretize num-interfacePatch )
#		 add_test(doublecut ${CMAKE_CURRENT_BINARY_DIR}/doublecut.exe)
if(FAIL)
        add_executable(ip2test.exe interp2)
        target_link_libraries(ip2test.exe tests-interp2 geom-simplex geom-multigrid geom-topo out-output num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils num-discretize )
        add_test(interp2 ${CMAKE_CURRENT_BINARY_DIR}/interp2.exe)
endif(FAIL)
        add_executable(prolongationp2test.exe prolongationp2test)
        target_link_libraries(prolongationp2test.exe tests-prolongationp2test geom-simplex geom-multigrid geom-topo num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils )
        add_test(prolongationp2test ${CMAKE_CURRENT_BINARY_DIR}/prolongationp2test.exe)

        add_executable(tetrabuildertest.exe tetrabuildertest)
        target_link_libraries(tetrabuildertest.exe tests-tetrabuildertest geom-simplex geom-multigrid geom-topo num-unknowns geom-builder geom-boundary misc-utils )
		add_test(tetrabuildertest ${CMAKE_CURRENT_BINARY_DIR}/tetrabuildertest.exe)

        add_executable(mattest.exe mattest)
        target_link_libraries(mattest.exe tests-mattest misc-utils )

        add_executable(testfe.exe testfe)
        target_link_libraries(testfe.exe tests-testfe misc-utils geom-topo geom-boundary geom-simplex geom-multigrid num-unknowns geom-builder )
    	add_test(testfe ${CMAKE_CURRENT_BINARY_DIR}/testfe.exe)

        add_executable(sbuffer.exe sbuffer)
        target_link_libraries(sbuffer.exe tests-sbuffer misc-utils )
    	add_test(sbuffer ${CMAKE_CURRENT_BINARY_DIR}/sbuffer.exe)
    	
        add_executable(minres.exe minres)
        target_link_libraries(minres.exe tests-minres misc-utils )
       	add_test(minres ${CMAKE_CURRENT_BINARY_DIR}/minres.exe)
       	
if(FAIL)
        add_executable(meshreader.exe meshreader)
        target_link_libraries(meshreader.exe tests-meshreader geom-simplex geom-multigrid geom-topo num-unknowns geom-builder geom-boundary misc-utils out-output misc-problem num-interfacePatch num-fe )
        add_test(meshreader ${CMAKE_CURRENT_BINARY_DIR}/meshreader.exe)
endif(FAIL)

        add_executable(restrictp2.exe restrictp2)
        target_link_libraries(restrictp2.exe tests-restrictp2 geom-simplex geom-multigrid geom-topo out-output num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils )
		add_test(restrictp ${CMAKE_CURRENT_BINARY_DIR}/restrictp.exe)
		
        add_executable(vectest.exe vectest)
        target_link_libraries(vectest.exe tests-vectest misc-utils )
		add_test(vectest ${CMAKE_CURRENT_BINARY_DIR}/vectest.exe)
if(FAIL)
        add_executable(p2local.exe p2local)
        target_link_libraries(p2local.exe geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns out-output misc-problem geom-builder num-interfacePatch num-discretize num-fe )
		add_test(p2local ${CMAKE_CURRENT_BINARY_DIR}/p2local.exe)
endif(FAIL)
        add_executable(quadbase.exe quadbase)
        target_link_libraries(quadbase.exe tests-quadbase tests-quadbase misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem num-fe num-discretize num-interfacePatch )
		add_test(quadbase ${CMAKE_CURRENT_BINARY_DIR}/quadbase.exe)
		
        add_executable(globallist.exe globallist)
        target_link_libraries(globallist.exe tests-globallist misc-utils )
		add_test(globallist ${CMAKE_CURRENT_BINARY_DIR}/globallist.exe)
		
if(FAIL)
        add_executable(serialization.exe serialization)
        target_link_libraries(serialization.exe tests-serialization misc-utils misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns out-output num-fe misc-problem num-interfacePatch )
   		add_test(serialization ${CMAKE_CURRENT_BINARY_DIR}/serialization.exe)
endif(FAIL)
        add_executable(triang.exe triang)
        target_link_libraries(triang.exe tests-triang misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns )
		add_test(triang ${CMAKE_CURRENT_BINARY_DIR}/triang.exe)
		
        add_executable(quadCut.exe quadCut)
        target_link_libraries(quadCut.exe tests-quadCut misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem num-interfacePatch levelset-levelset levelset-fastmarch num-discretize num-fe levelset-surfacetension geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature )
		add_test(quadCut ${CMAKE_CURRENT_BINARY_DIR}/quadCut.exe)
		
        add_executable(bicgstab.exe bicgstab)
        target_link_libraries(bicgstab.exe tests-bicgstab misc-utils )
		add_test(bicgstab ${CMAKE_CURRENT_BINARY_DIR}/bicgstab.exe)
		
        add_executable(gcr.exe gcr)
        target_link_libraries(gcr.exe tests-gcr misc-utils )
		add_test(gcr ${CMAKE_CURRENT_BINARY_DIR}/gcr.exe)
		
        add_executable(blockmat.exe blockmat)
        target_link_libraries(blockmat.exe tests-blockmat misc-utils )
		add_test(blockmat ${CMAKE_CURRENT_BINARY_DIR}/blockmat.exe)
		
        add_executable(mass.exe mass)
        target_link_libraries(mass.exe tests-mass misc-utils )
		add_test(mass ${CMAKE_CURRENT_BINARY_DIR}/mass.exe)
		
        add_executable(quad5.exe quad5)
        target_link_libraries(quad5.exe tests-quad5 misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem num-fe num-discretize num-interfacePatch geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature )
		add_test(quad5 ${CMAKE_CURRENT_BINARY_DIR}/quad5.exe)
        
        add_executable(downwind.exe downwind)
        target_link_libraries(downwind.exe tests-downwind misc-utils geom-multigrid geom-simplex geom-topo num-unknowns )
		add_test(downwind ${CMAKE_CURRENT_BINARY_DIR}/downwind.exe)
		
#        add_executable(2D.exe 2D)
#       target_link_libraries(2D.exe misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem num-fe num-discretize num-interfacePatch )
#		add_test(minres ${CMAKE_CURRENT_BINARY_DIR}/minres.exe)

        add_executable(interfaceP1FE.exe interfaceP1FE)
        target_link_libraries(interfaceP1FE.exe tests-interfaceP1FE misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem num-interfacePatch levelset-levelset levelset-fastmarch num-discretize num-fe levelset-surfacetension geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature )
		add_test(interfaceP1FE ${CMAKE_CURRENT_BINARY_DIR}/interfaceP1FE.exe)
		
        add_executable(xfem.exe xfem)
        target_link_libraries(xfem.exe tests-xfem misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem num-interfacePatch num-fe num-discretize levelset-levelset levelset-fastmarch stokes-instatstokes2phase out-ensightOut levelset-surfacetension misc-funcmap misc-vectorFunctions geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature )
		add_test(xfem ${CMAKE_CURRENT_BINARY_DIR}/xfem.exe)
		
#        add_executable(Gamma.exe Gamma)
#        target_link_libraries(Gamma.exe geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize num-interfacePatch misc-params levelset-fastmarch stokes-instatstokes2phase levelset-surfacetension misc-funcmap misc-vectorFunctions geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature )
#		add_test(Gamma ${CMAKE_CURRENT_BINARY_DIR}/Gamma.exe)

if(FAIL)
        add_executable(neq.exe neq)
        target_link_libraries(tests-neq misc-utils )
		add_test(neq ${CMAKE_CURRENT_BINARY_DIR}/neq.exe)
endif(FAIL)

        add_executable(extendP1onChild.exe extendP1onChild)
        target_link_libraries(extendP1onChild.exe tests-extendP1onChild num-discretize misc-utils geom-topo num-fe misc-problem geom-simplex geom-multigrid num-unknowns num-interfacePatch )
		add_test(extendP1onChild ${CMAKE_CURRENT_BINARY_DIR}/extendP1onChild.exe)

        add_executable(splitboundary.exe splitboundary)
        target_link_libraries(splitboundary.exe tests-splitboundary geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns misc-utils geom-topo out-output num-fe misc-problem num-interfacePatch )
		add_test(splitboundary ${CMAKE_CURRENT_BINARY_DIR}/splitboundary.exe)
		
#        add_executable(init.exe init)
#        target_link_libraries(init.exe geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns misc-utils geom-topo out-output misc-problem levelset-levelset num-fe misc-utils num-discretize num-interfacePatch levelset-fastmarch levelset-surfacetension out-ensightOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature )
#		add_test(minres ${CMAKE_CURRENT_BINARY_DIR}/minres.exe)

        add_executable(test-reparam.exe reparam)
        target_link_libraries(test-reparam.exe tests-reparam levelset-fastmarch levelset-levelset geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe out-ensightOut num-interfacePatch levelset-surfacetension misc-params out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature )
		add_test(test-reparam ${CMAKE_CURRENT_BINARY_DIR}/test-reparam.exe)
		
        add_executable(principallattice.exe principallattice)
        target_link_libraries(principallattice.exe tests-principallattice misc-utils geom-principallattice num-discretize geom-topo num-fe num-interfacePatch misc-problem num-unknowns geom-simplex geom-multigrid geom-builder geom-boundary num-unknowns geom-reftetracut geom-subtriangulation num-quadrature )
		add_test(principallattice ${CMAKE_CURRENT_BINARY_DIR}/principallattice.exe)
		
#       add_executable(extra.exe extra)
#       target_link_libraries(extra.exe misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem num-fe num-discretize num-interfacePatch geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature )
#		add_test(extra ${CMAKE_CURRENT_BINARY_DIR}/extra.exe)

#DIRECTSOLVER REQUIRES CHOLMOD/UMFPACK CHECK THAT
#       add_executable(directsolver.exe directsolver)
#       target_link_libraries(directsolver.exe tests-directsolver misc-utils )
#		add_test(directsolver ${CMAKE_CURRENT_BINARY_DIR}/directsolver.exe)

endif(NOT MPI)
