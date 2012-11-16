#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = navstokes

EXEC = nsdrops nsdrops_begehung insdrops insadrops

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)


# rules

default: ClrScr

nsdrops: \
    ../geom/boundary.o ../geom/builder.o ../geom/deformation.o ../navstokes/nsdrops.o \
	../misc/scopetimer.o ../misc/progressaccu.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o ../num/fe.o ../num/interfacePatch.o \
    ../misc/problem.o ../misc/utils.o ../out/output.o ../num/discretize.o ../geom/principallattice.o \
    ../geom/reftetracut.o
	$(CXX) -o $@ $^ $(LFLAGS)

nsdrops_begehung: \
    ../geom/boundary.o ../geom/builder.o ../geom/deformation.o ../navstokes/nsdrops_begehung.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o ../num/fe.o ../num/interfacePatch.o \
    ../misc/problem.o ../misc/utils.o ../out/output.o ../num/discretize.o ../geom/principallattice.o \
    ../geom/reftetracut.o
	$(CXX) -o $@ $^ $(LFLAGS)

insdrops: \
    ../geom/boundary.o ../geom/builder.o  ../geom/deformation.o ../navstokes/insdrops.o \
	../misc/scopetimer.o ../misc/progressaccu.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o ../num/fe.o ../num/interfacePatch.o \
    ../misc/problem.o ../misc/utils.o ../out/output.o ../num/discretize.o \
    ../out/ensightOut.o ../geom/principallattice.o ../geom/reftetracut.o
	$(CXX) -o $@ $^ $(LFLAGS)

insadrops: \
    ../geom/boundary.o ../geom/builder.o  ../geom/deformation.o ../navstokes/insadrops.o \
	../misc/scopetimer.o ../misc/progressaccu.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o ../num/fe.o ../num/interfacePatch.o \
    ../misc/problem.o ../misc/utils.o ../out/output.o ../num/discretize.o \
    ../out/ensightOut.o ../geom/principallattice.o ../geom/reftetracut.o
	$(CXX) -o $@ $^ $(LFLAGS)


# include common rules:
include ../Make.common
