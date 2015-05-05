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

ifeq ($(DEBUG),) # If DEBUG option is not set
    DEBUG:=0
endif

ifeq ($(CLUSTER),) # If CLUSTER option is not set
    CLUSTER:=0
endif

debug: FFLAGS = -DDEBUG=${DEBUG} -g -O0 -openmp -check bounds -fpp
#debug: FFLAGS =  -DDEBUG=${DEBUG} -g -ffree-line-length-0 -O0 -fopenmp
debug: executable

OBJS:=global.o par_zig_mod.o random.o hmm.o

%.o:%.f90
	${FC} ${FFLAGS} -c $<

executable: ${OBJS}
	$(eval FFLAGS += -DCLUSTER=$(CLUSTER))
	export OMP_STACKSIZE=" 128 M"
	export OMP_NUM_THREADS=4
	${FC} AlphaImpute.f90 ${OBJS} ${FFLAGS} -o ${PROGRAM}


clean:
	rm -f *.o *.mod *~ 

.PHONY: make clean

