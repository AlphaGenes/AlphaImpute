CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT(AlphaImpute)
enable_language (Fortran)
FIND_PACKAGE(Doxygen)

OPTION(RUN_TESTS "Run and compile tests"
       OFF)
OPTION(USE_HDF5 "Use hdf"
       OFF)

set(CLUSTER 0 CACHE STRING "cluster size")
# Add our local modlues to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/alphacmakefiles/cmake/modules")
set(CMAKE_MACOSX_RPATH 0)

INCLUDE(${CMAKE_MODULE_PATH}/setVersion.cmake)



add_definitions(-DCLUSTER=${CLUSTER})
# Define the executable name
SET(ALPHAIMPUTEEXE AlphaImpute)
SET(GENEPROBEXE GeneProbForAlphaImpute)
MESSAGE("Prog name ${ALPHAIMPUTEEXE}")
# Define the library name
SET(AHLIB alphahouse)
SET(GPLIB geneprob)
SET(APLIB alphaphase)
SET(AMLPLIB alphamlp)
SET(AHMMLIB alphahmm)
SET(FERLIB ferdosi)
# Define some directories
SET(SRC ${CMAKE_SOURCE_DIR}/src)
SET(OBJ ${CMAKE_SOURCE_DIR}/objs)
SET(BIN ${CMAKE_SOURCE_DIR}/bin)
SET(TESTS ${CMAKE_SOURCE_DIR}/tests)
SET(SRCALPHAIMPUTE ${SRC}/)


set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN})

INCLUDE(${CMAKE_MODULE_PATH}/findDoxygen.cmake)
# Defined libary source

# Uncomment if it is required that Fortran 90 is supported
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
   MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Set some options the user may choose
# Uncomment the below if you want the user to choose a parallelization library
OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
OPTION(USE_OPENMP "Use OpenMP for parallelization" ON)


if (${RUN_TESTS})
SET(testDeps ${SRC}/AlphaImputeInMod.f90
${SRC}/global.f90
)
endif (${RUN_TESTS})

# set(_files "test.pf")


INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

INCLUDE(${CMAKE_MODULE_PATH}/findMKL.cmake)
# Setup the LAPACK libraries.  This also takes care of peculiarities, such as
# the fact the searching for MKL requires a C compiler, and that the results
# are not stored in the cache.
#INCLUDE(${CMAKE_MODULE_PATH}/SetUpLAPACK.cmake)
if (${USE_HDF5})
	INCLUDE(${CMAKE_MODULE_PATH}/findHDF5.cmake)
endif(${USE_HDF5})

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING.  You should  review this file and make sure the flags
# are to your liking.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)
SET(SRCAH ./alphahouse/src/)
SET(SRCGP ./geneprob/src/)
SET(SRCAP ./alphaphase/src/)
SET(SRCAMLP ./alphamlp/src/)
SET(SRCAHMM ./alphahmm/src/)
SET(SRCFER ./ferdosi/src/)
# Have the .mod files placed in the OBJ folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${OBJ})


ADD_SUBDIRECTORY(${SRCAH} ${OBJ}/AH)
ADD_SUBDIRECTORY(${SRCGP} ${OBJ}/GP)
ADD_SUBDIRECTORY(${SRCAP} ${OBJ}/AP)
ADD_SUBDIRECTORY(${SRCAHMM} ${OBJ}/AHMM)
ADD_SUBDIRECTORY(${SRCFER} ${OBJ}/FER)
#ADD_SUBDIRECTORY(${SRCAMLP} ${OBJ}/AMLP)

# TODO add tests to correct location
ADD_SUBDIRECTORY(${SRCALPHAIMPUTE} ${BIN})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_MODULE_PATH}/distclean.cmake
)

if (${RUN_TESTS})
	INCLUDE(${CMAKE_MODULE_PATH}/findPFUnit.cmake)
endif(${RUN_TESTS})
