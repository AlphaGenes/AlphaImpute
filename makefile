# General variables
NAME:=AlphaImpute
VERSION:=1
SUBVERSION:=0
PROGRAM:=${NAME}${VERSION}.${SUBVERSION}

# Compiler
FC:=ifort

# Options
FFLAGS:=-O3 -m64

all:executable

debug: FFLAGS += -DDEBUG -g -O0
debug: executable

SRCS=$(wildcard *.f90)

OBJS:=$(SRCS:.o=.f90)

executable: ${OBJS}
	${FC} AlphaImpute.f90 ${FFLAGS} -o ${PROGRAM}

clean:
	rm -f *.o *.mod *~ 

.PHONY: make clean

