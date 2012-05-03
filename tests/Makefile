#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = tests

EXEC =  ip1test ip2test mattest prolongationp2test reftest testfe \
        tetrabuildertest sbuffer minres meshreader restrictp2 vectest \
        p2local quadbase globallist triang quadCut bicgstab gcr blockmat \
        mass quad5 downwind quad5_2D interfaceP1FE serialization xfem \
        directsolver f_Gamma neq splitboundary reparam_init reparam \
        extendP1onChild principallattice quad_extra dist_streamtest dist_remotedatatest \
        dist_sendtetraofbrick dist_simplexsend dist_interface dist_exchange \
        dist_poisson dist_sdropsP2 dist_twophasedrops dist_ref dist_modify \
        dist_TestRefPar dist_ref_rising_drop dist_migration dist_interp2 dist_interp1

DELETE = $(EXEC) *.out *.diff *.off *.mg *.dat

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)

# include common rules:
include ../Make.common

PAR_OBJ_ = ../parallel/parallel.o ../parallel/partime.o  \
           ../parallel/exchange.o  ../parallel/parmultigrid.o \
           ../parallel/logger.o ../parallel/DiST.o ../parallel/mpistream.o \
           ../parallel/loadbal.o ../parallel/migrateunknowns.o \
           ../parallel/decompose.o

PAR_OBJ = $(if $(PAR_BUILD),$(PAR_OBJ_),)

TESTNAME = drivencavity-adaptiv-schur drivencavity-adaptiv-uzawa \
    navierstokes-stationaer-uniform \
    poisson-instationaer-uniform poisson-stationaer-adaptiv \
    poisson-stationaer-multigrid poisson-stationaer-neumann \
    poisson-stationaer-uniform stokes-instationaer-uniform \
    stokes-stationaer-uniform

# rules

default: ClrScr

dist: dist_streamtest dist_remotedatatest dist_sendtetraofbrick dist_modify dist_simplexsend dist_interface dist_exchange dist_poisson dist_sdropsP2 dist_twophasedrops \
    dist_ref dist_TestRefPar dist_migration dist_interp2 dist_interp1

check: $(TESTNAME:%=check_%) ip1check
	@echo "--> All testcases have been executed!"

check_%:
	cd ../$$(./testinfo "makedir" $*) && $(MAKE) $$(../${DIR}/testinfo "maketarget" $*)
	./runtest $*

reftest: \
    ../tests/reftest.o ../geom/boundary.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../num/unknowns.o ../misc/utils.o ../geom/topo.o ../out/output.o ../misc/problem.o \
    ../num/interfacePatch.o ../num/fe.o 
	$(CXX) -o $@ $^ $(LFLAGS)

ip1test: \
    ../tests/interp1.o ../geom/simplex.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o ../num/interfacePatch.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

interfacepatch_doublecut: \
    ../tests/interfacepatch_doublecut.o ../misc/utils.o ../geom/builder.o ../geom/simplex.o \
    ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o ../num/unknowns.o ../misc/problem.o \
    ../num/fe.o ../num/discretize.o ../num/interfacePatch.o
	$(CXX) -o $@ $^ $(LFLAGS)

ip1check: ip1test
	@ ./ip1test > /dev/null 2>&1 || echo "WARNING: P1-repair failed."

ip2test: \
    ../tests/interp2.o ../geom/simplex.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o ../num/interfacePatch.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o ../num/discretize.o
	$(CXX) -o $@ $^ $(LFLAGS)

prolongationp2test: \
    ../tests/prolongationp2test.o ../geom/simplex.o ../geom/multigrid.o ../geom/topo.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o ../num/interfacePatch.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

tetrabuildertest: \
    ../tests/tetrabuildertest.o ../geom/simplex.o ../geom/multigrid.o ../geom/topo.o \
    ../num/unknowns.o ../geom/builder.o ../geom/boundary.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

mattest: \
    ../tests/mattest.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

testfe: \
    ../tests/testfe.o ../misc/utils.o ../geom/topo.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o
	$(CXX) -o $@ $^ $(LFLAGS)

sbuffer: \
    ../tests/sbuffer.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

minres: \
    ../tests/minres.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

meshreader: \
    ../tests/meshreader.o ../geom/simplex.o ../geom/multigrid.o ../geom/topo.o \
    ../num/unknowns.o ../geom/builder.o ../geom/boundary.o ../misc/utils.o \
    ../out/output.o ../misc/problem.o ../num/interfacePatch.o ../num/fe.o 
	$(CXX) -o $@ $^ $(LFLAGS)

restrictp2: \
    ../tests/restrictp2.o ../geom/simplex.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o ../num/interfacePatch.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

vectest: \
    ../tests/vectest.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

p2local: ../tests/p2local.o ../misc/utils.o \
  ../geom/simplex.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o ../num/unknowns.o \
  ../out/output.o  ../misc/problem.o ../geom/builder.o ../num/interfacePatch.o \
  ../num/discretize.o ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

quadbase: \
    ../tests/quadbase.o ../misc/utils.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o ../misc/problem.o \
    ../num/fe.o ../num/discretize.o ../num/interfacePatch.o
	$(CXX) -o $@ $^ $(LFLAGS)

globallist: \
    ../tests/globallist.o  ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

serialization: \
    ../tests/serialization.o  ../misc/utils.o  ../misc/utils.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o ../out/output.o ../num/fe.o ../misc/problem.o ../num/interfacePatch.o
	$(CXX) -o $@ $^ $(LFLAGS)

triang: \
    ../tests/triang.o  ../misc/utils.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o
	$(CXX) -o $@ $^ $(LFLAGS)

quadCut: \
    ../tests/quadCut.o  ../misc/utils.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o ../misc/problem.o ../num/interfacePatch.o \
    ../levelset/levelset.o ../levelset/fastmarch.o ../num/discretize.o ../num/fe.o ../levelset/surfacetension.o \
    ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o ../num/quadrature.o
	$(CXX) -o $@ $^ $(LFLAGS)

bicgstab: \
    ../tests/bicgstab.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

gcr: \
    ../tests/gcr.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)
blockmat: \
    ../tests/blockmat.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

mass: \
    ../tests/mass.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

quad5: \
    ../tests/quad5.o ../misc/utils.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o ../misc/problem.o \
    ../num/fe.o ../num/discretize.o ../num/interfacePatch.o\
    ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o ../num/quadrature.o
	$(CXX) -o $@ $^ $(LFLAGS)

downwind: \
    ../tests/downwind.o ../misc/utils.o ../geom/multigrid.o ../geom/simplex.o ../geom/topo.o ../num/unknowns.o
	$(CXX) -o $@ $^ $(LFLAGS)

quad5_2D: \
    ../tests/quad5_2D.o ../misc/utils.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o ../misc/problem.o \
    ../num/fe.o ../num/discretize.o ../num/interfacePatch.o
	$(CXX) -o $@ $^ $(LFLAGS)

interfaceP1FE: \
    ../tests/interfaceP1FE.o  ../misc/utils.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o ../misc/problem.o ../num/interfacePatch.o \
    ../levelset/levelset.o ../levelset/fastmarch.o ../num/discretize.o ../num/fe.o ../levelset/surfacetension.o \
    ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o ../num/quadrature.o
	$(CXX) -o $@ $^ $(LFLAGS)

xfem: \
    ../tests/xfem.o ../misc/utils.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o ../misc/problem.o ../num/interfacePatch.o \
    ../num/fe.o ../num/discretize.o ../levelset/levelset.o ../levelset/fastmarch.o \
    ../stokes/instatstokes2phase.o ../out/ensightOut.o ../levelset/surfacetension.o \
    ../misc/bndmap.o ../geom/bndVelFunctions.o \
    ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o ../num/quadrature.o 
	$(CXX) -o $@ $^ $(LFLAGS)

f_Gamma: \
    ../tests/f_Gamma.o ../geom/boundary.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../num/interfacePatch.o \
    ../misc/params.o ../levelset/fastmarch.o ../stokes/instatstokes2phase.o \
    ../levelset/surfacetension.o ../misc/bndmap.o ../geom/bndVelFunctions.o \
    ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o ../num/quadrature.o 
	$(CXX) -o $@ $^ $(LFLAGS)

neq: \
    ../tests/neq.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

extendP1onChild: \
	../tests/extendP1onChild.o ../num/discretize.o ../misc/utils.o ../geom/topo.o ../num/fe.o ../misc/problem.o \
	../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../num/interfacePatch.o
	$(CXX) -o $@ $^ $(LFLAGS)

splitboundary: \
    ../tests/splitboundary.o ../geom/boundary.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../num/unknowns.o ../misc/utils.o ../geom/topo.o ../out/output.o ../num/fe.o ../misc/problem.o ../num/interfacePatch.o
	$(CXX) -o $@ $^ $(LFLAGS)

reparam_init: \
    ../tests/reparam_init.o ../geom/boundary.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../num/unknowns.o ../misc/utils.o ../geom/topo.o ../out/output.o \
    ../misc/problem.o ../levelset/levelset.o ../num/fe.o \
    ../misc/utils.o ../num/discretize.o ../num/interfacePatch.o \
    ../levelset/fastmarch.o ../levelset/surfacetension.o ../out/ensightOut.o \
    ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o ../num/quadrature.o
	$(CXX) -o $@ $^ $(LFLAGS)

reparam: \
    ../tests/reparam.o ../levelset/fastmarch.o ../levelset/levelset.o \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o \
    ../num/fe.o ../out/ensightOut.o ../num/interfacePatch.o ../levelset/surfacetension.o \
    ../misc/params.o ../out/vtkOut.o ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o \
    ../num/quadrature.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

principallattice: \
    ../tests/principallattice.o ../misc/utils.o ../geom/principallattice.o ../num/discretize.o ../geom/topo.o \
    ../num/fe.o ../num/interfacePatch.o ../misc/problem.o ../num/unknowns.o ../geom/simplex.o \
    ../geom/multigrid.o ../geom/builder.o ../geom/boundary.o ../num/unknowns.o ../geom/reftetracut.o \
    ../geom/subtriangulation.o ../num/quadrature.o
	$(CXX) -o $@ $^ $(LFLAGS)

quad_extra: \
    ../tests/quad_extra.o ../misc/utils.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o ../misc/problem.o \
    ../num/fe.o ../num/discretize.o ../num/interfacePatch.o\
    ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o ../num/quadrature.o
	$(CXX) -o $@ $^ $(LFLAGS)

dist_streamtest: \
    ../tests/dist_streamtest.o \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o \
    ../num/fe.o ../num/interfacePatch.o ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_remotedatatest: \
    ../tests/dist_remotedatatest.o \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o \
    ../num/fe.o ../num/interfacePatch.o ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_sendtetraofbrick: \
    ../tests/dist_sendtetraofbrick.o \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o \
    ../num/fe.o ../num/interfacePatch.o ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_modify: \
    ../tests/dist_modify.o \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o \
    ../num/fe.o ../num/interfacePatch.o ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_simplexsend: \
    ../tests/dist_simplexsend.o  \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o \
    ../num/fe.o ../num/interfacePatch.o ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_interface: \
    ../tests/dist_interface.o  \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o \
    ../num/fe.o ../num/interfacePatch.o ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_ref: \
    ../tests/dist_ref.o \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o ../out/vtkOut.o ../parallel/parmultigrid.o \
    ../num/fe.o ../num/interfacePatch.o ../levelset/levelset.o ../levelset/fastmarch.o ../levelset/surfacetension.o \
    ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o ../num/quadrature.o \
    ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_migration: \
    ../tests/dist_migration.o \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o ../out/vtkOut.o ../parallel/parmultigrid.o \
    ../num/fe.o ../num/interfacePatch.o ../levelset/levelset.o ../levelset/fastmarch.o ../levelset/surfacetension.o \
    ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o ../num/quadrature.o $(PAR_OBJ) \
    ../misc/bndmap.o ../geom/bndVelFunctions.o ../geom/bndScalarFunctions.o \
    ../stokes/instatstokes2phase.o ../misc/scopetimer.o
	$(CXX) -o $@ $^ $(LFLAGS)

dist_TestRefPar: \
    ../tests/dist_TestRefPar.o \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o ../out/vtkOut.o ../parallel/parmultigrid.o \
    ../num/fe.o ../num/interfacePatch.o ../levelset/levelset.o ../levelset/fastmarch.o ../levelset/surfacetension.o \
    ../misc/params.o ../out/output.o ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o \
    ../num/quadrature.o ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_ref_rising_drop: \
    ../tests/dist_ref_rising_drop.o \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o ../out/vtkOut.o ../parallel/parmultigrid.o \
    ../num/fe.o ../num/interfacePatch.o ../levelset/levelset.o ../levelset/fastmarch.o ../levelset/surfacetension.o \
    ../misc/params.o ../out/output.o ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o ../num/quadrature.o \
    ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)


dist_exchange: \
    ../tests/dist_exchange.o  \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o \
    ../num/fe.o ../num/interfacePatch.o ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_poisson: \
    ../tests/dist_poisson.o ../geom/boundary.o ../geom/builder.o\
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o \
    ../poisson/poisson.o ../misc/problem.o ../misc/utils.o ../out/output.o \
    ../num/fe.o ../num/discretize.o ../num/interfacePatch.o ../geom/geomselect.o\
    ../misc/params.o ../out/vtkOut.o ../misc/bndmap.o ../geom/bndScalarFunctions.o\
    ../geom/bndVelFunctions.o ../poisson/poissonCoeff.o \
    ../num/hypre.o ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_sdropsP2: \
    ../tests/dist_sdropsP2.o ../geom/boundary.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o \
    ../misc/utils.o ../num/discretize.o ../num/interfacePatch.o ../num/stokessolver.o \
    ../num/fe.o ../stokes/integrTime.o ../geom/reftetracut.o ../geom/principallattice.o \
    ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_twophasedrops: \
    ../tests/dist_twophasedrops.o ../geom/boundary.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../navstokes/instatnavstokes2phase.o \
    ../misc/params.o ../levelset/fastmarch.o ../stokes/instatstokes2phase.o ../num/MGsolver.o \
    ../num/fe.o ../stokes/integrTime.o  ../num/interfacePatch.o ../out/vtkOut.o ../levelset/surfacetension.o \
    ../geom/geomselect.o  ../levelset/twophaseutils.o ../num/hypre.o ../levelset/coupling.o ../geom/bndVelFunctions.o \
    ../geom/bndScalarFunctions.o ../misc/bndmap.o ../levelset/twophaseCoeff.o  \
    ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o ../num/quadrature.o \
    ../num/renumber.o ../surfactant/ifacetransp.o \
    ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)
	
dist_interp2: \
    ../tests/dist_interp2.o ../geom/simplex.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o ../num/interfacePatch.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o ../num/discretize.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)
	
dist_interp1: \
    ../tests/dist_interp1.o ../geom/simplex.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o ../num/interfacePatch.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o ../num/discretize.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

../tests/directsolver.o:
ifneq ($(DSOLVER_CXXFLAGS),)
	$(CXX) $(CXXFLAGS) $(DSOLVER_CXXFLAGS) -o $*.o -c $*.cpp
else
	@echo "compilation skipped due to missing CHOLMOD/UMFPACK installation"
endif

ifneq ($(DSOLVER_LFLAGS),)
directsolver: \
    ../tests/directsolver.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS) $(DSOLVER_LFLAGS)
else
directsolver:
	@echo "compilation skipped due to missing CHOLMOD/UMFPACK installation"
endif

