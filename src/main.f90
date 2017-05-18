!######################################################################

program AlphaImpute
    ! The program is intended to be run numerous times
    ! The first time AlphaImpute is run in order to calculate Genotype probabilities (GenoProb)
    ! The second time, it is run to phase the individuals
    ! Finally, AlphaImpute should be run to impute genotypes.
    ! RestartOptions handle this situation, so:
    !   * inputParams%restartOption=0 => Passes through the whole process: GenoProb, Phasing and Imputing
    !   * inputParams%restartOption=1 => Makes only GenoProbReadInParameterFile
    !   * inputParams%restartOption=2 => Makes GenoProb and Phasing
    !   * inputParams%restartOption=3 => Makes only Imputation. This implies AlphaImpute has to be run
    !                           already in order to get GenoProb done or Genotype Probabilities
    !                           Genotype Probabilities have to be edited by hand
    !   * inputParams%restartOption=4 =>
    use Global
    use AlphaImputeModule
    use informationModule
    use GlobalVariablesHmmMaCH
    use Output
    use AlphaImputeInMod
    use Imputation
    use InputMod
    use GeneProbModule
    use AlphaPhaseResultsDefinition
    implicit none

    character(len=4096) :: cmd, SpecFile

    inputParams => defaultInput
    if (Command_Argument_Count() > 0) then
        call get_command_argument(1,cmd)
        if (cmd(1:2) .eq. "-v") then
            call PrintVersion
            call exit(0)
        end if
    end if

    if (Command_Argument_Count() > 0) then
        call Get_Command_Argument(1,SpecFile)
    else
        specfile="AlphaImputeSpec.txt"
    end if


    call Titles

    ! use default input, TODO this can be changed
    allocate(defaultInput)
    call defaultInput%ReadInParameterFile(SpecFile)

    inputParams => defaultInput

    if (inputParams%hmmoption /= RUN_HMM_NGS) then
        if (inputParams%restartOption<OPT_RESTART_PHASING) call MakeDirectories(RUN_HMM_NULL)

        call ReadInData
        !call cpu_time(start)
        call SnpCallRate
        call CheckParentage

        if (inputParams%MultiHD/=0) call ClassifyAnimByChips

        call FillInSnp

        call FillInBasedOnOffspring
        call InternalEdit


        if (inputparams%cluster) then
            call MakeFiles
        else if (inputParams%PreProcess==.true.) then
            print*, "  "
            print*, "  ","ERROR: PREPROCESSING OPTION IS NO LONGER AVAILABLE WITH CLUSTER MODE DISABLED"
            stop
        endif
    else

        call MakeDirectories(RUN_HMM_NGS)
        call ReadInData
        call SnpCallRate
        allocate(Reads(ped%nGenotyped,inputParams%nsnp))
        allocate(ImputeGenos(0:ped%pedigreeSize,inputParams%nsnpRaw))
        allocate(ImputePhase(0:ped%pedigreeSize,inputParams%nsnpRaw,2))
        allocate(SnpIncluded(inputParams%nsnp))
        call CheckParentage
        call ReadSeq(inputParams%GenotypeFileUnit)
    endif

    if (inputParams%hmmoption == RUN_HMM_NGS) then

        call MaCHController(inputParams%hmmoption)
        call FromHMM2ImputePhase
        call WriteOutResults

    else if (inputParams%hmmoption==RUN_HMM_ONLY) then

        print*, ""
        print*, "Bypass calculation of probabilities and phasing"

    else ! if hmm option is not ngs or only
        write(6,*) " "
        write(6,*) " ","Data editing completed"

        if (inputParams%SexOpt==0) then
            select case (inputParams%bypassgeneprob)

                case (0)

                    if (inputParams%restartOption== OPT_RESTART_ALL .or. inputParams%restartOption== OPT_RESTART_GENEPROB) Then
                        ! call ped%addGenotypeInformation(imputeGenos)
                        ! WriteOutResults is a piece of shit and makes life hard
                        ! call ped%addGenotypeInformation(Genos)
                        if (inputParams%cluster) then
#if CLUSTER==1
                           call GeneProbManagement
#else

                            write(6,*) ""
                            write(6,*) "Restart option 0 or 1 stops program before Geneprobs jobs have been submitted with cluster option set"
                            stop
#endif                           
                        endif
                        print *, "Calling geneprob"
                        call runGeneProbAlphaImpute(1, inputParams%nsnp, ped, GenosProbs, MAF)
                        print *, "writing probabilities"
                        call WriteProbabilitiesFull("./GeneProb/GenotypeProbabilities.txt", GenosProbs, ped,ped%pedigreeSize-ped%nDummys)
                        call WriteProbabilities("./Results/GenotypeProbabilities.txt", GenosProbs, ped,ped%pedigreeSize-ped%nDummys, inputParams%nsnp)
                    endif

                    ! deallocate(GenosProbs)

                    if (inputParams%restartOption==OPT_RESTART_GENEPROB) then
                        call ped%writeOutGenotypes("./GeneProb/individualGenotypes.txt")
                        call WriteProbabilitiesFull("./GeneProb/GenotypeProbabilities.txt", GenosProbs, ped,ped%pedigreeSize-ped%nDummys)
                        write(6,*) "Restart option 1 stops program after Geneprobs jobs have finished"
                        stop
                    endif
                    write(6,*) " "
                    write(6,*) " ","Genotype probabilities calculated"
                    

                case (1)
                    write(6,*) " "
                    write(6,*) " ","Genotype probabilities bypassed."
                    if (inputParams%restartOption==OPT_RESTART_GENEPROB) then
                        write(6,*) "Restart option 1 stops program after Geneprobs jobs have finished"
                        write(6,*) "Warning - BYPASSGENEPROB has been given yes. Thus these options are incompatible."
                        stop
                    endif
                    ! TODO check this 
                case (2)

                    write(6,*) "Restart option 1 stops program after genotype probabilities have been outputted"
                    stop
                case default
                    print *, "ERROR: BYPASS GENEPROB SET INCORRECTLY"
                    stop 1
            end select
        endif



    if (inputParams%managephaseon1off0==1) then

        if (inputParams%restartOption<OPT_RESTART_IMPUTATION) Then
            call PhasingManagementNew(APResults)

        endif

        if (inputParams%restartOption==OPT_RESTART_PHASING) then
            block
                use OutputParametersDefinition
                use InputOutput
                integer :: i
                type(OutputParameters) :: oParams
                oParams = newOutputParametersImpute()
                do i=1, apResults%nResults
                    write(oParams%outputDirectory,'("./Phasing/Phase"i0)') i
                    call writeAlphaPhaseResults(APResults%results(i), ped, oParams)

                enddo
                write(6,*) "Restart option 2 stops program after Phasing has been managed"
                stop
            end block
        endif
    endif


endif

if (inputParams%hmmoption/=RUN_HMM_NGS) then
        if (inputParams%restartOption> OPT_RESTART_PHASING) Then
            print *,"Reading in Phasing information"
            ! Read back in geneprob data

            if (inputParams%BypassGeneProb == 0) then
                if (inputParams%cluster) then
                    call readProbabilitiesFullCluster(GenosProbs,ped%pedigreeSize-ped%nDummys, inputParams%nsnp,inputparams,GpIndex)
                else
                    call readProbabilitiesFull("./GeneProb/GenotypeProbabilities.txt",GenosProbs,ped%pedigreeSize-ped%nDummys, inputParams%nsnp)
                endif
            endif
            block 

                use OutputParametersDefinition
                use InputOutput
                use AlphaPhaseResultsDefinition
                integer :: i
                type(OutputParameters) :: oParams
                oParams = newOutputParametersImpute()
                ApResults%nResults = size(inputParams%CoreLengths)
                allocate(ApResults%results(ApResults%nResults))
                do i=1, ApResults%nResults
                    write(oParams%outputDirectory,'("./Phasing/Phase"i0)') i
                    call readAlphaPhaseResults(ApResults%results(i), oParams, ped)
                enddo
            end block
        endif
        print *, "Phasing Completed"

    ! If we only want to phase data, then skip all the imputation steps
    if (inputParams%PhaseTheDataOnly==0) Then
        call ImputationManagement

        call WriteOutResults

        if (inputparams%restartOption == OPT_RESTART_IMPUTATION) then
            write(*,*) "Restart option 3 stops program after Iterate Geneprob jobs have been finished"
            stop
        endif

#ifdef DEBUG
        write(0,*) 'DEBUG: Model Recombination'
#endif
        ! WARNING: Skip the modelling the recombination because it interferes with HMM propabilites
        ! TODO:
        if (inputParams%hmmoption==RUN_HMM_NO) call ModelRecomb

#ifdef DEBUG

        write(0,*) 'DEBUG: Final Checker'
#endif

        if (inputParams%TrueGenos1None0==1) then
            call FinalChecker

        endif
    endif
endif
call ped%destroyPedigree()
call PrintTimerTitles


end program AlphaImpute

