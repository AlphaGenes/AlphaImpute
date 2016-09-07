# General variables
NAME:=AlphaImpute
VERSION:=$(shell git rev-parse --short HEAD)
MASTERVERSION:=$(shell git describe --tag)
# SUBVERSION:=2.0
# PROGRAM:=${NAME}${VERSION}.${SUBVERSION}
PROGRAM:=$(NAME)$(MASTERVERSION)
GP:=GeneProbForAlphaImpute

# Compiler
FC:=ifort
#FC:=gfortran

# Options
FFLAGS:=-O3 -m64 -DVERS=""commit-$(VERSION)""

# Set precompilation options by default
CLUSTER?=0
DEBUG?=0

# MS Windows
ifeq ($(OS), Windows_NT)
	OSFLAG := "OS_WIN"
	FFLAGS := $(FFLAGS) /static /fpp  /D $(OSFLAG) -D CLUSTER=0 /Qopenmp /libs:static

	obj:= .obj
	exe:= .exe

	DEL:= del
else
	OSFLAG := "OS_UNIX"
	FFLAGS:= $(FFLAGS) -qopenmp -static-intel -fpp -qopenmp-link=static  -D $(OSFLAG) -D CLUSTER=$(CLUSTER)

	obj:= .o
	exe:= 

	DEL:= rm -rf
endif

# Options
# FFLAGS:=-O3 -m64 -openmp -fpp -DCLUSTER=$(CLUSTER) -openmp-link=static -static-intel
# FFLAGS:=-O3 -m64 -fopenmp -ffree-line-length-0

GPFFLAGS:=-O3 -m64 -static-intel

all: executable

debug: FFLAGS = -DDEBUG=${DEBUG} -g -O0 -qopenmp -check bounds -fpp -static-intel -qopenmp-link=static -traceback -D $(OSFLAG)
#debug: FFLAGS =  -DDEBUG=${DEBUG} -g -ffree-line-length-0 -O0 -fopenmp
debug: PROGRAM:=$(NAME)_debug
debug: executable

OBJS:=global$(obj) par_zig_mod$(obj) random$(obj) hmmvariables$(obj) Output$(obj) hmmHaplotype$(obj) utils$(obj) recombination$(obj) hmm$(obj) PhaseRounds$(obj) HaplotypeBits$(obj) Imputation$(obj)

ifeq ($(OS), Windows_NT)
	 OBJS:= $(OBJS) manageWindows.obj
endif

%$(obj):%.f90
	$(FC) $(FFLAGS) -c $<

executable: $(OBJS)
	$(FC) AlphaImpute.f90 $(OBJS) $(FFLAGS) -o $(PROGRAM)$(exe)
	${FC} GeneProbForAlphaImpute.f90 ${GPFFLAGS} -o ${GP}

clean:
	$(DEL) *$(obj) *.mod *~

veryclean:
	$(DEL) *$(obj)*.mod *~ AlphaImpute$(exe)

.PHONY: make clean veryclean all

