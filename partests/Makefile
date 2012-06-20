#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = parpartests

DELETE = $(EXEC) *.out *.diff *.off *.mg *.dat
CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)

# include common rules:
include ../Make.common

EXEC =  dist_streamtest dist_remotedatatest dist_sendtetraofbrick dist_simplexsend \
        dist_interface dist_exchange dist_sdropsP2 dist_ref dist_modify \
        dist_TestRefPar dist_ref_rising_drop dist_migration dist_interp2 dist_interp1
            
PAR_OBJ_ = ../DiST/DiST.o ../DiST/mpistream.o ../DiST/remotedata.o \
           ../parallel/parallel.o ../parallel/partime.o \
           ../parallel/exchange.o ../parallel/parmultigrid.o \
           ../parallel/loadbal.o ../parallel/migrateunknowns.o \
           ../parallel/decompose.o

PAR_OBJ = $(if $(PAR_BUILD),$(PAR_OBJ_),)

# rules

all: $(EXEC)

default: ClrScr

dist: dist_streamtest dist_remotedatatest dist_sendtetraofbrick dist_modify dist_simplexsend dist_interface dist_exchange dist_poisson dist_sdropsP2 dist_twophasedrops \
    dist_ref dist_TestRefPar dist_migration dist_interp2 dist_interp1 dist_ref_rising_drop


dist_streamtest: \
    ../partests/dist_streamtest.o \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o \
    ../num/fe.o ../num/interfacePatch.o ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_remotedatatest: \
    ../partests/dist_remotedatatest.o \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o \
    ../num/fe.o ../num/interfacePatch.o ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_sendtetraofbrick: \
    ../partests/dist_sendtetraofbrick.o \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o \
    ../num/fe.o ../num/interfacePatch.o ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_modify: \
    ../partests/dist_modify.o \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o \
    ../num/fe.o ../num/interfacePatch.o ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_simplexsend: \
    ../partests/dist_simplexsend.o  \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o \
    ../num/fe.o ../num/interfacePatch.o ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_interface: \
    ../partests/dist_interface.o  \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o \
    ../num/fe.o ../num/interfacePatch.o ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_ref: \
    ../partests/dist_ref.o \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o ../out/vtkOut.o ../parallel/parmultigrid.o \
    ../num/fe.o ../num/interfacePatch.o ../levelset/levelset.o ../levelset/fastmarch.o ../levelset/surfacetension.o \
    ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o ../num/quadrature.o \
    ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_migration: \
    ../partests/dist_migration.o \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o ../out/vtkOut.o ../parallel/parmultigrid.o \
    ../num/fe.o ../num/interfacePatch.o ../levelset/levelset.o ../levelset/fastmarch.o ../levelset/surfacetension.o \
    ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o ../num/quadrature.o $(PAR_OBJ) \
    ../misc/bndmap.o ../geom/bndVelFunctions.o ../geom/bndScalarFunctions.o \
    ../stokes/instatstokes2phase.o ../misc/scopetimer.o
	$(CXX) -o $@ $^ $(LFLAGS)

dist_TestRefPar: \
    ../partests/dist_TestRefPar.o \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o ../out/vtkOut.o ../parallel/parmultigrid.o \
    ../num/fe.o ../num/interfacePatch.o ../levelset/levelset.o ../levelset/fastmarch.o ../levelset/surfacetension.o \
    ../misc/params.o ../out/output.o ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o \
    ../num/quadrature.o ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_ref_rising_drop: \
    ../partests/dist_ref_rising_drop.o \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o ../out/vtkOut.o ../parallel/parmultigrid.o \
    ../num/fe.o ../num/interfacePatch.o ../levelset/levelset.o ../levelset/fastmarch.o ../levelset/surfacetension.o \
    ../misc/params.o ../out/output.o ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o ../num/quadrature.o \
    ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)


dist_exchange: \
    ../partests/dist_exchange.o  \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o \
    ../num/fe.o ../num/interfacePatch.o ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_sdropsP2: \
    ../partests/dist_sdropsP2.o ../geom/boundary.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o \
    ../misc/utils.o ../num/discretize.o ../num/interfacePatch.o \
    ../num/fe.o ../geom/reftetracut.o ../geom/principallattice.o \
    ../misc/scopetimer.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

dist_interp2: \
    ../partests/dist_interp2.o ../geom/simplex.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o ../num/interfacePatch.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o ../num/discretize.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)
	
dist_interp1: \
    ../partests/dist_interp1.o ../geom/simplex.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o ../num/interfacePatch.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o ../num/discretize.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)
