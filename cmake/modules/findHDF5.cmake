# INCLUDE(${CMAKE_MODULE_PATH}/SetCompileFlags.cmake)



set(HDF5_USE_STATIC_LIBRARIES ON)
# find_package(HDF5 REQUIRED)
SET (HDF5PATH "/usr/local/hdf5_2015/")

# if (HDF5PATH)
# 	SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS}"
#                 	Fortran "-L${HDF5PATH}lib -I${HDF5PATH}include -lhdf5_fortran -lhdf5"
#                 	)
# endif(HDF5PATH) 
