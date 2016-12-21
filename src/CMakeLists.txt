SET(ALPHAIMPUTE_src
            ${SRCALPHAIMPUTE}AlphaImputeInMod.f90
            ${SRCALPHAIMPUTE}global.f90
            ${SRCALPHAIMPUTE}HaplotypeBits.f90
            ${SRCALPHAIMPUTE}hmm.f90
            ${SRCALPHAIMPUTE}Output.f90
            ${SRCALPHAIMPUTE}hmmvariables.f90
            ${SRCALPHAIMPUTE}hmmHaplotype.f90          
            ${SRCALPHAIMPUTE}Imputation.f90
            ${SRCALPHAIMPUTE}informationModule.f90
            ${SRCALPHAIMPUTE}InputMod.f90
            ${SRCALPHAIMPUTE}manageWindows.f90
            ${SRCALPHAIMPUTE}par_zig_mod.f90
            ${SRCALPHAIMPUTE}PhaseRounds.f90
            ${SRCALPHAIMPUTE}random.f90
            ${SRCALPHAIMPUTE}recombination.f90
            ${SRCALPHAIMPUTE}utils.f90
            ${SRCALPHAIMPUTE}AlphaImpute.f90

)

add_definitions(-D COMMIT=${VERSION})
ADD_EXECUTABLE(${ALPHAIMPUTEEXE} ${ALPHAIMPUTE_src})

SET(GENEPROB_src ${SRCALPHAIMPUTE}GeneProbForAlphaImpute.f90)
ADD_EXECUTABLE(${GENEPROBEXE} ${GENEPROB_src})

set(TESTPROG testProg)


# add external links here
TARGET_LINK_LIBRARIES(${ALPHAIMPUTEEXE} ${AHLIB})
TARGET_LINK_LIBRARIES(${ALPHAIMPUTEEXE} ${RESULT})

# TARGET_LINK_LIBRARIES(${TESTPROG} ${RESULT})
INSTALL(TARGETS ${ALPHAIMPUTEEXE} RUNTIME DESTINATION bin)
INSTALL(TARGETS ${GENEPROBEXE} RUNTIME DESTINATION bin)
