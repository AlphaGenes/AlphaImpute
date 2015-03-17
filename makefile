# General variables
NAME:=AlphaImpute
VERSION:=1
SUBVERSION:=0
PROGRAM:=${NAME}${VERSION}.${SUBVERSION}

# Compiler
FC:=ifort
#FC:=gfortran

# Options
FFLAGS:=-O3 -m64 -openmp -heap-arrays
#FFLAGS:=-O3 -m64 -fopenmp -ffree-line-length-0

all: executable

debug: FFLAGS += -DDEBUG -g -O0 -openmp -traceback -heap-arrays
#debug: FFLAGS =  -DDEBUG -g -ffree-line-length-0 -O0 -fopenmp
debug: executable

OBJS:=global.o hmm.o

%.o:%.f90
	${FC} ${FFLAGS} -c $<

executable: ${OBJS}
	${FC} AlphaImpute.f90 ${OBJS} ${FFLAGS} -o ${PROGRAM}

clean:
	rm -f *.o *.mod *~ 

.PHONY: make clean

