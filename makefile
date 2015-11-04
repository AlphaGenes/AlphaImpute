# General variables
NAME:=AlphaImpute
VERSION:=1
SUBVERSION:=2.0
PROGRAM:=${NAME}${VERSION}.${SUBVERSION}
GP:=GeneProbForAlphaImpute

# Compiler
FC:=ifort
#FC:=gfortran

# Set precompilation options by default
CLUSTER?=0
DEBUG?=0

# Options
FFLAGS:=-O3 -m64 -openmp -fpp -DCLUSTER=$(CLUSTER) -openmp-link=static -static-intel
#FFLAGS:=-O3 -m64 -fopenmp -ffree-line-length-0

GPFFLAGS:=-O3 -m64 -static-intel

all: executable

debug: FFLAGS = -DDEBUG=${DEBUG} -g -O0 -openmp -check bounds -fpp -static-intel -openmp-link=static -traceback
#debug: FFLAGS =  -DDEBUG=${DEBUG} -g -ffree-line-length-0 -O0 -fopenmp
debug: executable

OBJS:=global.o par_zig_mod.o random.o hmm.o Output.o

%.o:%.f90
	${FC} ${FFLAGS} -c $<

executable: ${OBJS}
	export OMP_STACKSIZE=" 128 M"
	export OMP_NUM_THREADS=4
	${FC} AlphaImpute.f90 ${OBJS} ${FFLAGS} -o ${PROGRAM}
	${FC} GeneProbForAlphaImpute.f90 ${GPFFLAGS} -o ${GP}

clean:
	rm -f *.o *.mod *~ 

.PHONY: make clean

