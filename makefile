# General variables
NAME:=AlphaImpute
VERSION:=1
SUBVERSION:=0
PROGRAM:=${NAME}${VERSION}.${SUBVERSION}

# Compiler
FC:=ifort
#FC:=gfortran

# Options
FFLAGS:=-O3 -m64 -openmp -fpp
#FFLAGS:=-O3 -m64 -fopenmp -ffree-line-length-0

all: executable

debug: FFLAGS = -DDEBUG -g -O0 -openmp -check bounds -fpp
#debug: FFLAGS =  -DDEBUG -g -ffree-line-length-0 -O0 -fopenmp
debug: executable

OBJS:=global.o par_zig_mod.o random.o hmm.o

%.o:%.f90
	${FC} ${FFLAGS} -c $<

executable: ${OBJS}
	export OMP_STACKSIZE=" 128 M"
	export OMP_NUM_THREADS=4
	${FC} AlphaImpute.f90 ${OBJS} ${FFLAGS} -o ${PROGRAM}

clean:
	rm -f *.o *.mod *~ 

.PHONY: make clean

