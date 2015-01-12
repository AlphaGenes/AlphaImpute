# General variables
NAME:=AlphaImpute
VERSION:=1
SUBVERSION:=0
PROGRAM:=${NAME}${VERSION}.${SUBVERSION}

# Compiler
FC:=ifort

# Options
FFLAGS:=-O3 -m64 -ffree-line-length-0

SRCS=$(wildcard *.f90)

OBJS:=$(SRCS:.o=.f90)

all: ${OBJS}
	${FC} AlphaImpute.f90 ${FFLAGS} -o ${PROGRAM}

clean:
	rm -f *.o *.mod *~ 

.PHONY: make clean

