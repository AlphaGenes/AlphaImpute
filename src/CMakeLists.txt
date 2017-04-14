SET(ALPHAIMPUTE_src
            ${SRCALPHAIMPUTE}AlphaImputeInMod.f90
            ${SRCALPHAIMPUTE}global.f90
            ${SRCALPHAIMPUTE}HaplotypeBits.f90
            ${SRCALPHAIMPUTE}hmm.f90
            ${SRCALPHAIMPUTE}Output.f90
            ${SRCALPHAIMPUTE}hmmvariables.f90
            ${SRCALPHAIMPUTE}hmmHaplotype.f90      
            ${SRCALPHAIMPUTE}informationModule.f90    
            ${SRCALPHAIMPUTE}Imputation.f90
            ${SRCALPHAIMPUTE}InputMod.f90
            ${SRCALPHAIMPUTE}utils.f90
            ${SRCALPHAIMPUTE}AlphaImpute.f90

)
set (HMM_src
    global.f90
    par_zig_mod.f90
    random.f90   
    AlphaImputeInMod.f90
    hmm.f90
    utils.f90
    hmmvariables.f90
     InputMod.f90
            hmmHaplotype.f90  )


if (ALPHAIMPUTEEXE)      
    add_definitions(-D COMMIT=${VERSION})
    ADD_EXECUTABLE(${ALPHAIMPUTEEXE} ${ALPHAIMPUTE_src})
    # ADD_EXECUTABLE(${GENEPROBEXE} ${SRCALPHAIMPUTE}GeneProbForAlphaImpute.f90)
    # add external links here
    TARGET_LINK_LIBRARIES(${ALPHAIMPUTEEXE} ${AHLIB})
    TARGET_LINK_LIBRARIES(${ALPHAIMPUTEEXE} ${GPLIB})
    TARGET_LINK_LIBRARIES(${ALPHAIMPUTEEXE} ${APLIB})
    
endif()

# if (HMMLIB)
#     ADD_LIBRARY(${HMMLIB} STATIC ${HMM_src})
#     TARGET_LINK_LIBRARIES(${HMMLIB} ${AHLIB})
# endif()
#INSTALL(TARGETS ${ALPHAIMPUTEEXE} RUNTIME DESTINATION bin)
