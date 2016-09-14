.PHONY: git make clean veryclean all
# General variables
NAME:=AlphaImpute
VERSION:= $(shell git rev-parse --short HEAD)
SUBVERSION:=0
MASTERVERSION:=$(shell git describe --tag)
PROGRAM:=$(NAME)$(MASTERVERSION)
GP:=GeneProbForAlphaImpute

# Set the default compiler to iFort
FC:=ifort
FFLAGS:=-O3 -DCOMMIT=" $(VERSION)"

# Set precompilation options by default
CLUSTER?=0
DEBUG?=0

ifeq ($(OS), Windows_NT)
	# MS Windows
	SRCDIR      := src/
	BUILDDIR    :=
	TARGETDIR   :=
	PFUNIT :=
	TESTDIR :=
	OSFLAG := "OS_WIN"
	FFLAGS := $(FFLAGS) /static /i8 /fpp /D $(OSFLAG) -D CLUSTER=0 /Qopenmp /libs:static
	GPFFLAGS := -O3 /static /i8 /libs:static

	obj := .obj

	MAKEDIR :=
	exe := .exe
	CC := cl
	CFLAGS := /EHsc

	DEL := del
else
	# Linux or Mac OSX
	SRCDIR      :=src/
	BUILDDIR    := objs/
	TARGETDIR   := bin/
	PFUNIT:=/usr/local/pFUnit_serial
	TESTDIR:=tests/
	obj := .o
	OSFLAG := "OS_UNIX"
	ABOPT := -qopenmp -static-intel -qopenmp-link=static
	exe :=
	FFLAGS:= $(FFLAGS) -m64 -fpp -module $(BUILDDIR) -D $(OSFLAG) $(ABOPT) -D CLUSTER=$(CLUSTER)
	GPFFLAGS := $(FFLAGS) -m64 $(ABOPT) -module $(BUILDDIR)
	uname := $(shell uname)
	MAKEDIR := @mkdir -p
	DEL := rm -rf

	DEBUG_FLAGS:= -traceback -g -qopenmp -check bounds -check format -check output_conversion -check pointers -check uninit -DDEBUG=${DEBUG}
endif

# THIS_DIR :=$(shell pwd)/

TES:= $(THIS_DIR)$(TESTDIR)

# FILES:=$(wildcard $(THIS_DIR)$(SRCDIR)*.f90)
FILES:= $(THIS_DIR)$(SRCDIR)global.f90 \
		$(THIS_DIR)$(SRCDIR)par_zig_mod.f90 \
		$(THIS_DIR)$(SRCDIR)random.f90 \
		$(THIS_DIR)$(SRCDIR)hmmvariables.f90 \
		$(THIS_DIR)$(SRCDIR)Output.f90 \
		$(THIS_DIR)$(SRCDIR)hmmHaplotype.f90 \
		$(THIS_DIR)$(SRCDIR)utils.f90 \
		$(THIS_DIR)$(SRCDIR)recombination.f90 \
		$(THIS_DIR)$(SRCDIR)hmm.f90 \
		$(THIS_DIR)$(SRCDIR)PhaseRounds.f90 \
		$(THIS_DIR)$(SRCDIR)HaplotypeBits.f90 \
		$(THIS_DIR)$(SRCDIR)Imputation.f90 \
		$(THIS_DIR)$(SRCDIR)AlphaImpute.f90

PFTESTS:=$(wildcard $(TES)*.pf)

F90TESTS:=$(PFTESTS:.pf=.F90)

BUILDDATE=$(shell date +%Y%m%d-%H:%M:%S)

OBJS:= global$(obj) par_zig_mod$(obj) random$(obj) hmmvariables$(obj) Output$(obj) hmmHaplotype$(obj) utils$(obj) recombination$(obj) hmm$(obj) PhaseRounds$(obj) HaplotypeBits$(obj) Imputation$(obj)

ifeq ($(OS), Windows_NT)
	 OBJS:= $(OBJS) manageWindows.obj
endif

# all: files preprocess_test clean

# all_git: files git preprocess_test clean

# debug:  files_debug preprocess_test clean

# debug_git:  files_debug git preprocess_test clean

# preprocess_test: $(F90TESTS)
# 	cd $(TES)
# 	$(FC) -o tests.x $^ -I$(THIS_DIR)$(BUILDDIR) $(PFUNIT)/include/driver.f90 -I$(PFUNIT)/mod -I$(PFUNIT)/include -L$(PFUNIT)/lib -lpfunit -module $(BUILDDIR)
# 	./tests.x

all: files clean

all_git: files clean

debug:  files_debug clean

debug_git:  files_debug git clean

clean:
	$(RM) $(THIS_DIR)$(TESTDIR)*.F90
	$(RM) $(BUILDDIR)*.mod
	$(RM) $(BUILDDIR)tests.x
	$(RM) $(BUILDDIR)testSuites.inc

files: $(FILES)
	$(FC) -o $(TARGETDIR)$(PROGRAM) $^ $(FFLAGS)
	$(FC) -o $(TARGETDIR)$(GP) $(GPFFLAGS) $(THIS_DIR)$(SRCDIR)GeneProbForAlphaImpute.f90

files_debug: $(FILES)
	$(FC) -o $(TARGETDIR)$(PROGRAM) $^ $(FFLAGS) $(DEBUG_FLAGS)

git: $(FILES) $(PFTESTS)
	git add $^
	git commit --allow-empty -m "Successful compilation $(BUILDDATE)"

%.F90: %.pf
	$(PFUNIT)/bin/pFUnitParser.py $< $@

%$(obj):%.f90
	$(FC) $(FFLAGS) -c -o $(BUILDDIR)$@ $<






