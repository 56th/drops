#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = parallel
EXEC =
DELETE = $(EXEC)
CPP = $(wildcard *.cpp)
OBJ = $(CPP:%.cpp=%.o)

# rules:

default:

all: $(EXEC)

PAR_OBJ = parallel.o parmultigrid.o partime.o loadbal.o exchange.o logger.o

# include common rules:
include ../Make.common
