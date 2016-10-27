
SET (PFUNIT_DIR "/usr/local/pFUnit_2015")


enable_testing()

# SET(_files ${TESTS}/TestASIndividualModule.pf)


if(DEFINED PFUNIT_DIR)
    # message(STATUS "Manual setup of variable PFUNIT_DIR : ${PFUNIT_DIR}")
    set(PFUNIT_DIR ${PFUNIT_DIR})
else(DEFINED PFUNIT_DIR)
    include(ExternalProject)

    set(ExternalProjectCMakeArgs
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external/pfunit
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        )
    ExternalProject_Add(pfunit
        DOWNLOAD_COMMAND git submodule update
        DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/pfunit
        BINARY_DIR ${PROJECT_BINARY_DIR}/external/pfunit-build
        STAMP_DIR ${PROJECT_BINARY_DIR}/external/pfunit-stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/external/pfunit-tmp
        INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        CMAKE_ARGS ${ExternalProjectCMakeArgs}
        )
    include_directories(${PROJECT_BINARY_DIR}/external/pfunit/mod)
    set(PFUNIT_DIR ${PROJECT_BINARY_DIR}/external/pfunit)
endif(DEFINED PFUNIT_DIR)

file(MAKE_DIRECTORY ${LIB}/generated)
file(WRITE ${LIB}/generated/testSuites.inc "")

include_directories(${PFUNIT_DIR}/source)

include_directories(
    # ${TESTS}
    ${SRC}
    # ${LIB}
    ${PFUNIT_DIR}/mod
    ${LIB}/generated
    )
message(STATUS "Manual setup of variable PFUNIT_DIR : ${SRC}/*.f90")
set(_test_sources)
if (${_files})
    set(files ${_files})
else()
    file(GLOB files "${TESTS}/*.pf")
endif()
    foreach(_file ${files})

    get_filename_component (_test ${_file} NAME_WE)
    set(test_dependency ${TESTS}/${_test}.pf ${SRC}/*.f90)
    add_custom_command(
        OUTPUT ${LIB}/generated/${_test}.F90
        COMMAND python ${PFUNIT_DIR}/bin/pFUnitParser.py ${TESTS}/${_test}.pf ${LIB}/generated/${_test}.F90
        DEPENDS ${test_dependency}
        )
    set(_test_sources ${_test_sources} ${LIB}/generated/${_test}.F90)
    file(APPEND ${LIB}/generated/testSuites.inc "ADD_TEST_SUITE(${_test}_suite)\n")
endforeach()

set_source_files_properties(${PFUNIT_DIR}/include/driver.F90 PROPERTIES GENERATED 1)

SET(testDeps ${SRCALPHASIM}/AlphaSimInputModule.f90
    ${SRCALPHASIM}/GenomeModule.f90   
    ${SRCALPHASIM}/HaplotypesModule.f90   
        ${SRCALPHASIM}/AlphaSimIndividualModule.f90      
)
# file(GLOB testDeps "${SRC}/*.f90")
add_executable(
    pftest_alltests
    ${testDeps}
    ${PFUNIT_DIR}/include/driver.F90
    ${_test_sources}
    
    )

target_link_libraries(
    pftest_alltests
    ${PFUNIT_DIR}/lib/libpfunit.a
    )

TARGET_LINK_LIBRARIES(pftest_alltests ${AHLIB})
TARGET_LINK_LIBRARIES(pftest_alltests ${MACSLIB})

add_test(pftest_alltests ${BIN}/pftest_alltests)
INSTALL(TARGETS pftest_alltests RUNTIME DESTINATION bin)