#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = poisson

EXEC = cdrdrops

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)

# include common rules:
include ../Make.common

PAR_OBJ_ = ../DiST/DiST.o ../DiST/mpistream.o ../DiST/remotedata.o \
           ../parallel/parallel.o ../parallel/partime.o \
           ../parallel/exchange.o ../parallel/parmultigrid.o \
           ../parallel/loadbal.o ../parallel/migrateunknowns.o \
           ../parallel/decompose.o

PAR_OBJ = $(if $(PAR_BUILD),$(PAR_OBJ_),)

# rules

default: ClrScr

ipdropsAD: \
    ../poisson/poissonCoeff.o \
    ../geom/boundary.o ../geom/builder.o ../poisson/ipdropsAD.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o ../num/interfacePatch.o \
    ../poisson/poisson.o ../misc/problem.o ../misc/utils.o ../misc/funcmap.o \
    ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

cdrdrops:\
    ../poisson/poissonCoeff.o \
    ../poisson/ale.o \
    ../poisson/cdrdrops.o ../geom/boundary.o ../geom/builder.o ../out/ensightOut.o\
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o \
    ../poisson/poisson.o ../misc/problem.o ../misc/utils.o ../out/output.o \
    ../num/fe.o ../num/discretize.o ../num/interfacePatch.o ../geom/geomselect.o\
    ../misc/params.o \
    ../out/vtkOut.o ../misc/funcmap.o ../misc/scalarFunctions.o\
    ../misc/vectorFunctions.o\
    ../geom/deformation.o \
    ../misc/scopetimer.o \
    $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

matlab/ipdrops:\
    ../poisson/poissonCoeff.o \
    ../geom/boundary.o ../geom/builder.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o ../num/interfacePatch.o \
    ../poisson/poisson.o ../misc/problem.o ../out/output.o ../misc/utils.o
	$(CXX) $(CXXFLAGS) matlab/ipdrops.cpp -o $@ $^ $(LFLAGS)

matlab/ipfilm:\
    ../poisson/poissonCoeff.o \
    ../geom/boundary.o ../geom/builder.o ../num/interfacePatch.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o \
    ../poisson/poisson.o ../misc/problem.o ../out/output.o ../misc/utils.o
	$(CXX) $(CXXFLAGS) matlab/ipfilm.cpp -o $@ $^ $(LFLAGS)


clean: matlab/clean

matlab/clean:
	cd matlab && $(MAKE) clean
