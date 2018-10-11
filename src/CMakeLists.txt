SET(ALPHAIMPUTE_src
            AlphaImputeSpecFileModule.f90
            global.f90
            AlphaImputeInputOutput.f90
            Imputation.f90
            HeuristicGeneprob.f90
            AlphaImpute.f90
            HeuristicGeneprob.f90
            main.f90

)



if (ALPHAIMPUTEEXE)      
    add_definitions(-D COMMIT=${VERSION})
    ADD_EXECUTABLE(${ALPHAIMPUTEEXE} ${ALPHAIMPUTE_src})

    # Link to external libraries
    TARGET_LINK_LIBRARIES(${ALPHAIMPUTEEXE} ${AHLIB})

    TARGET_LINK_LIBRARIES(${ALPHAIMPUTEEXE} ${APLIB})
    TARGET_LINK_LIBRARIES(${ALPHAIMPUTEEXE} ${AHMMLIB})
    TARGET_LINK_LIBRARIES(${ALPHAIMPUTEEXE} ${FERLIB})
endif()
