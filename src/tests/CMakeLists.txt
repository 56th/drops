set(HOME tests)

exec_ser(reftest tests-reftest geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns misc-utils geom-topo out-output misc-problem num-interfacePatch num-fe)

exec_ser(interp1 geom-simplex geom-multigrid geom-topo out-output num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils)

exec_ser(interfacepatch_doublecut tests-interfacepatch_doublecut misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem num-fe num-discretize num-interfacePatch)

exec_ser(interp2 geom-simplex geom-multigrid geom-topo out-output num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils num-discretize)

exec_ser(prolongationp2test tests-prolongationp2test geom-simplex geom-multigrid geom-topo num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils)

exec_ser(tetrabuildertest tests-tetrabuildertest geom-simplex geom-multigrid geom-topo num-unknowns geom-builder geom-boundary misc-utils)

exec_ser(mattest misc-utils)

exec_ser(testfe tests-testfe misc-utils geom-topo geom-boundary geom-simplex geom-multigrid num-unknowns geom-builder)

exec_ser(sbuffer tests-sbuffer misc-utils)

exec_ser(minres tests-minres misc-utils)

exec_ser(meshreader tests-meshreader geom-simplex geom-multigrid geom-topo num-unknowns geom-builder geom-boundary misc-utils out-output misc-problem num-interfacePatch num-fe)

exec_ser(restrictp2 tests-restrictp2 geom-simplex geom-multigrid geom-topo out-output num-unknowns geom-builder misc-problem num-interfacePatch num-fe geom-boundary misc-utils)

exec_ser(vectest tests-vectest misc-utils)

exec_ser(p2local tests-p2local misc-utils geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns out-output misc-problem geom-builder num-interfacePatch num-discretize num-fe)

exec_ser(quadbase tests-quadbase misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem num-fe num-discretize num-interfacePatch) 

exec_ser(globallist tests-globallist misc-utils)

exec_ser(serialization tests-serialization misc-utils misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns out-output num-fe misc-problem num-interfacePatch)

exec_ser(triang tests-triang misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns)

exec_ser(quadCut tests-quadCut misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem num-interfacePatch levelset-levelset levelset-fastmarch num-discretize num-fe levelset-surfacetension geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

exec_ser(bicgstab tests-bicgstab misc-utils)

exec_ser(gcr tests-gcr misc-utils)

exec_ser(blockmat tests-blockmat misc-utils)

exec_ser(mass tests-mass misc-utils)

exec_ser(quad5 tests-quad5 misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem num-fe num-discretize num-interfacePatch geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

exec_ser(downwind tests-downwind misc-utils geom-multigrid geom-simplex geom-topo num-unknowns)

exec_ser(quad5_2D tests-quad5_2D misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem num-fe num-discretize num-interfacePatch)

exec_ser(interfaceP1FE tests-interfaceP1FE misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem num-interfacePatch levelset-levelset levelset-fastmarch num-discretize num-fe levelset-surfacetension geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

exec_ser(xfem tests-xfem misc-utils geom-builder geom-simplex geom-multigrid geom-boundary geom-topo num-unknowns misc-problem num-interfacePatch num-fe num-discretize levelset-levelset levelset-fastmarch stokes-instatstokes2phase out-ensightOut levelset-surfacetension misc-funcmap misc-vectorFunctions geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

exec_ser(f_Gamma tests-f_Gamma geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem levelset-levelset misc-utils out-output num-discretize num-interfacePatch misc-params levelset-fastmarch stokes-instatstokes2phase levelset-surfacetension misc-funcmap misc-vectorFunctions geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature )

exec_ser(neq tests-neq misc-utils) 

exec_ser(extendP1onChild tests-extendP1onChild num-discretize misc-utils geom-topo num-fe misc-problem geom-simplex geom-multigrid num-unknowns num-interfacePatch)

exec_ser(splitboundary tests-splitboundary geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns misc-utils geom-topo out-output num-fe misc-problem num-interfacePatch)

exec_ser(reparam_init tests-reparam_init geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns misc-utils geom-topo out-output misc-problem levelset-levelset num-fe misc-utils num-discretize num-interfacePatch levelset-fastmarch levelset-surfacetension out-ensightOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature)

#reparam renamed to reparamtest
exec_ser(reparamtest levelset-fastmarch levelset-levelset geom-simplex geom-multigrid geom-builder geom-topo geom-boundary num-unknowns misc-utils misc-problem num-discretize num-fe out-ensightOut num-interfacePatch levelset-surfacetension misc-params out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature) ##${PAR_OBJ}

exec_ser(principallattice tests-principallattice misc-utils geom-principallattice num-discretize geom-topo num-fe num-interfacePatch misc-problem num-unknowns geom-simplex geom-multigrid geom-builder geom-boundary num-unknowns geom-reftetracut geom-subtriangulation num-quadrature)

if(EXISTS ${UMFPACK_HOME})
    
endif()

######################## MISSING ############################
# directsolver

############### TEST FAILS TO COMPILE ####################### 
#exec_ser(quad_extra misc-utils geom-builder geom-simplex geom-multigrid 
# geom-boundary geom-topo num-unknowns misc-problem 
# num-fe num-discretize num-interfacePatch
# geom-principallattice geom-reftetracut geom-subtriangulation num-#quadrature)
