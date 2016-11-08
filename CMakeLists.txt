# CMak project file for FOO

##################################################
# Define the project and the depencies that it has
##################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT(AlphaImpute)
enable_language (Fortran)
SET(VERSION 1.6)

OPTION(RUN_TESTS "Run and compile tests"
       OFF)
OPTION(USE_HDF5 "Use hdf"
       OFF)

set(CLUSTER 0 CACHE STRING "cluster size")
# Add our local modlues to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
set(CMAKE_MACOSX_RPATH 0)


add_definitions(-DCLUSTER=${CLUSTER})

############################################################
# Define the actual files and folders that make up the build
############################################################

# Define the executable name
SET(ALPHAIMPUTEEXE AlphaImpute${VERSION})
SET(GENEPROBEXE GeneProbForAlphaImpute)
# Define the library name
SET(AHLIB alphahouse)

# Define some directories
SET(SRC ${CMAKE_SOURCE_DIR}/src)
SET(LIB ${CMAKE_SOURCE_DIR}/objs)
SET(BIN ${CMAKE_SOURCE_DIR}/bin)
SET(TESTS ${CMAKE_SOURCE_DIR}/tests)
SET(SRCALPHAIMPUTE ${SRC}/)

# Defined libary source

# Uncomment if it is required that Fortran 90 is supported
#IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
#    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
#ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Set some options the user may choose
# Uncomment the below if you want the user to choose a parallelization library
OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
OPTION(USE_OPENMP "Use OpenMP for parallelization" ON)

 

INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)
# Locate and set parallelization libraries.  There are some CMake peculiarities
# taken care of here, such as the fact that the FindOpenMP routine doesn't know
# about Fortran.
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
SET(SRCAH ../alphahouse/src/)

# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${LIB})

# The source for the BAR library and have it placed in the lib folder
ADD_SUBDIRECTORY(${SRCAH} ${LIB}/AH)


# TODO add tests to correct location
ADD_SUBDIRECTORY(${SRCALPHAIMPUTE} ${BIN})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)

if (${RUN_TESTS})
	INCLUDE(${CMAKE_MODULE_PATH}/findPFUnit.cmake) 
endif(${RUN_TESTS})

# ADD_CUSTOM_TARGET(tests 
#     COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/t.cmake
# )