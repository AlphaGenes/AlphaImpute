SET(ALPHAIMPUTE_src
            AlphaImputeSpecFileModule.f90
            global.f90
            HaplotypeBits.f90
            hmm.f90
            AlphaImputeInputOutput.f90
            hmmvariables.f90
            hmmHaplotype.f90      
            informationModule.f90    
            Imputation.f90
            utils.f90
            AlphaImpute.f90
            main.f90

)
set (HMM_src
    global.f90
    par_zig_mod.f90
    random.f90   
    AlphaImputeSpecFileModule.f90
    hmm.f90
    utils.f90
    hmmvariables.f90
     AlphaImputeInputOutput.f90
            hmmHaplotype.f90  )


if (ALPHAIMPUTEEXE)      
    add_definitions(-D COMMIT=${VERSION})
    ADD_EXECUTABLE(${ALPHAIMPUTEEXE} ${ALPHAIMPUTE_src})
    # ADD_EXECUTABLE(${GENEPROBEXE} ${SRCALPHAIMPUTE}GeneProbForAlphaImpute.f90)
    # add external links here
    TARGET_LINK_LIBRARIES(${ALPHAIMPUTEEXE} ${AHLIB})
    TARGET_LINK_LIBRARIES(${ALPHAIMPUTEEXE} ${GPLIB})
    #TARGET_LINK_LIBRARIES(${ALPHAIMPUTEEXE} ${AMLPLIB})
    TARGET_LINK_LIBRARIES(${ALPHAIMPUTEEXE} ${APLIB})
    
endif()

# if (HMMLIB)
#     ADD_LIBRARY(${HMMLIB} STATIC ${HMM_src})
#     TARGET_LINK_LIBRARIES(${HMMLIB} ${AHLIB})
# endif()
#INSTALL(TARGETS ${ALPHAIMPUTEEXE} RUNTIME DESTINATION bin)
