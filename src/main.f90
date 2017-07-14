#ifdef _WIN32

#define STRINGIFY(x)#x
#define TOSTRING(x) STRINGIFY(x)

#DEFINE DASH "\"
#DEFINE COPY "copy"
#DEFINE MD "md"
#DEFINE RMDIR "RMDIR /S /Q"
#DEFINE RM "del"
#DEFINE RENAME "MOVE /Y"
#DEFINE SH "BAT"
#DEFINE EXE ".exe"
#DEFINE NULL " >NUL"


#else

#define STRINGIFY(x)#x
#define TOSTRING(x) STRINGIFY(x)

#DEFINE DASH "/"
#DEFINE COPY "cp"
#DEFINE MD "mkdir"
#DEFINE RMDIR "rm -r"
#DEFINE RM "rm"
#DEFINE RENAME "mv"
#DEFINE SH "sh"
#DEFINE EXE ""
#DEFINE NULL ""


#endif
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
    use AlphaImputeInputOutputModule
    use AlphaImputeSpecFileModule
    use Imputation
    use ModuleRunFerdosi

    use AlphaPhaseResultsModule    
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
        ! call InitialiseArrays
        !call cpu_time(start)
        call SnpCallRate
        call CheckParentage
        if (inputParams%MultiHD/=0) then 
            call ClassifyAnimByChips
        endif
        
        call FillInSnp

        call FillInBasedOnOffspring
        call InternalEdit

        if (inputParams%PreProcess==.true.) then
            print*, "Data preprocessed"
            stop
        endif
        allocate(GlobalWorkPhase(0:ped%pedigreeSize,inputParams%nsnpraw,2))
        call InitialiseArrays

    else

        call MakeDirectories(RUN_HMM_NGS)
        call ReadInData
        call SnpCallRate
        allocate(SnpIncluded(inputParams%nsnp))
        call CheckParentage
        call ped%addSequenceFromFile(inputparams%GenotypeFile, inputParams%nsnpRaw, MAX_READS_COUNT)
    endif

        

    if (inputParams%hmmoption == RUN_HMM_NGS) then


        block
            use AlphaHmmInMod
            use ExternalHMMWrappers
            type (AlphaHMMinput) :: inputParamsHMM
            integer(kind=1) ,dimension(:,:), allocatable :: res

            inputParamsHMM%nsnp = inputParams%nsnp
            inputParamsHMM%nHapInSubH = inputParams%nHapInSubH
            inputParamsHMM%HmmBurnInRound = inputParams%HmmBurnInRound
            inputParamsHMM%nRoundsHmm = inputParams%nRoundsHmm
            inputParamsHMM%useProcs = inputParams%useProcs
            inputParamsHMM%imputedThreshold = inputParams%imputedThreshold
            inputParamsHMM%phasedThreshold = inputParams%phasedThreshold
            inputParamsHMM%HapList = inputParams%HapList

            res = ped%getGenotypesAsArray()
            call AlphaImputeHMMRunner(inputParamsHMM, ped, ProbImputeGenosHmm, ProbImputePhaseHmm, GenosCounts, FullH)


        end block
        call FromHMM2ImputePhase
        call WriteOutResults

    else if (inputParams%hmmoption==RUN_HMM_ONLY) then

        print*, ""
        print*, "Bypass calculation of probabilities and phasing"

    else ! if hmm option is not ngs or only
        write(6,*) " "
        write(6,*) " ","Data editing completed"

        if (inputParams%useFerdosi) then

            call doFerdosi(ped)
        endif



    if (inputParams%managephaseon1off0==1) then


        if (inputParams%restartOption<OPT_RESTART_IMPUTATION) Then
            call PhasingManagementNew(APResults)

        endif

        if (inputParams%restartOption==OPT_RESTART_PHASING) then
            block
                use OutputParametersModule
                use InputOutput
                integer :: i
                type(OutputParameters) :: oParams
                oParams = newOutputParametersImpute()
                do i=1, apResults%nResults
                    write(oParams%outputDirectory,'("."a"Phasing",a,"Phase"i0)') DASH,DASH, i
                    call writeAlphaPhaseResults(APResults%results(i), ped, oParams)

                enddo
                write(6,*) "Restart option 1 stops program after Phasing has been managed"
                stop
            end block
        endif
    endif


endif



if (inputParams%hmmoption/=RUN_HMM_NGS) then
        if (inputParams%restartOption> OPT_RESTART_PHASING) Then
            print *,"Reading in Phasing information"
            ! Read back in geneprob data
            block 

                use OutputParametersModule
                use InputOutput
                use AlphaPhaseResultsModule
                integer :: i
                type(OutputParameters) :: oParams
                oParams = newOutputParametersImpute()
                ApResults%nResults = size(inputParams%CoreLengths)*2
                allocate(ApResults%results(ApResults%nResults))
                do i=1, ApResults%nResults
                     write(oParams%outputDirectory,'("."a"Phasing",a,"Phase"i0)') DASH,DASH, i
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
            write(*,*) "Restart option 2 stops program after Imputation has finished"
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
            ! call FinalChecker
            ! TODO remimplement final checker

        endif
    endif
endif
call ped%destroyPedigree()
call PrintTimerTitles


end program AlphaImpute

