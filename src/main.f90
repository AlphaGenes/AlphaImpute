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
    use informationModule
    use GlobalVariablesHmmMaCH
    use Output
    use AlphaImputeInMod
    use Imputation
    use InputMod
    implicit none

    integer :: markers
    double precision, allocatable :: GenosProbs(:,:,:)
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

        !call cpu_time(start)
        call CountInData
        !call cpu_time(finish)
        !print '("Time ReadInData= ",f6.3," seconds.")',finish-start

        !call cpu_time(start)
        call ReadInData
        !call cpu_time(finish)
        !print '("Time ReadInData= ",f6.3," seconds.")',finish-start

        !call cpu_time(start)
        call SnpCallRate
        !call cpu_time(finish)
        !print '("Time SnpCallRate= ",f6.3," seconds.")',finish-start

        !call cpu_time(start)
        call CheckParentage
        !call cpu_time(finish)
        !print '("Time CheckParentage= ",f6.3," seconds.")',finish-start

        if (inputParams%MultiHD/=0) call ClassifyAnimByChips

        !call cpu_time(start)
        call FillInSnp
        !call cpu_time(finish)
        !print '("Time FillInSnp= ",f6.3," seconds.")',finish-start

        !call cpu_time(start)
        call FillInBasedOnOffspring
        !call cpu_time(finish)
        !print '("Time FillInBasedOnOffspring= ",f6.3," seconds.")',finish-start

        !call cpu_time(start)
        call InternalEdit
        !call cpu_time(finish)
        !print '("Time InternalEdit= ",f6.3," seconds.")',finish-start

        !call cpu_time(start)
        call MakeFiles
        !call cpu_time(finish)
        !print '("Time MakeFiles= ",f6.3," seconds.")',finish-start

    else

        call MakeDirectories(RUN_HMM_NGS)
        call CountInData
        call ReadInData
        call SnpCallRate
        allocate(Reads(nAnisG,inputParams%nsnp))
        allocate(ImputeGenos(0:nAnisG,inputParams%nsnp))
        allocate(ImputePhase(0:nAnisG,inputParams%nsnp,2))
        allocate(SnpIncluded(inputParams%nsnp))
        call CheckParentage
        call ReadSeq(inputParams%GenotypeFileUnit)
    endif

    if (inputParams%hmmoption == RUN_HMM_NGS) then

#ifdef DEBUG
        write(0,*) 'DEBUG: HMM NGS'
#endif
        call MaCHController(inputParams%hmmoption)
        call FromHMM2ImputePhase
        call WriteOutResults

    else if (inputParams%hmmoption==RUN_HMM_ONLY) then

#ifdef DEBUG
        write(0,*) 'DEBUG: HMM only'
#endif

        print*, ""
        print*, "Bypass calculation of probabilities and phasing"

    else
        write(6,*) " "
        write(6,*) " ","Data editing completed"

        if (inputParams%SexOpt==0) then
            select case (inputParams%bypassgeneprob)
            !        if (inputParams%bypassgeneprob==0) then
        case (0)

#ifdef DEBUG
            write(0,*) 'DEBUG: Calculate Genotype Probabilites'
#endif

            if (inputParams%restartOption== OPT_RESTART_ALL .or. inputParams%restartOption== OPT_RESTART_GENEPROB) Then

#ifndef _WIN32
#if CLUSTER==2
                write(6,*) ""
                write(6,*) "Restart option 1 stops program before Geneprobs jobs have been submitted"
                stop
#else
                call GeneProbManagement
#endif
#else
                !call cpu_time(start)
                call GeneProbManagementWindows

                !call cpu_time(finish)
                !print '("Time GeneProbManagementWindows= ",f6.3," seconds.")',finish-start

#endif

            endif

            markers = inputParams%nsnp
            if (inputParams%outopt==1) then
                markers = inputParams%nSnpRaw
            end if
            allocate(GenosProbs(ped%nDummys + nAnisP, markers, 2))
            call ReReadIterateGeneProbs(GenosProbs, .FALSE., nAnisP)
            call WriteProbabilities("./Results/GenotypeProbabilities.txt", GenosProbs, ped,nAnisP, inputParams%nsnp)
            deallocate(GenosProbs)

            if (inputParams%restartOption==OPT_RESTART_GENEPROB) then
#if CLUSTER==1
                write(6,*) "Restart option 1 stops program before Geneprobs jobs have finished"
#elif CLUSTER==0
                write(6,*) "Restart option 1 stops program after Geneprobs jobs have finished"
#endif
                stop
            endif
            write(6,*) " "
            write(6,*) " ","Genotype probabilities calculated"
            !        endif
        case (2)
            markers = inputParams%nsnp
            if (inputParams%outopt==1) then
                markers = inputParams%nSnpRaw
            end if
            allocate(GenosProbs(ped%nDummys+ nAnisP, markers, 2))
            call ReReadIterateGeneProbs(GenosProbs, .FALSE., nAnisP)
            call WriteProbabilities("./Results/GenotypeProbabilities.txt", GenosProbs, ped,nAnisP, inputParams%nsnp)
            deallocate(GenosProbs)
            write(6,*) "Restart option 1 stops program after genotype probabilities have been outputted"
            stop
        end select
    endif



    if (inputParams%managephaseon1off0==1) then

#ifdef DEBUG
        write(0,*) 'DEBUG: Phase haplotypes with AlphaPhase'
#endif

        if (inputParams%restartOption<OPT_RESTART_IMPUTATION) Then
#ifndef _WIN32
#if CLUSTER==2
            write(6,*) ""
            write(6,*) "Restart option 1 stops program before Phasing has been managed"
            stop
#else
            call PhasingManagement
#endif
#else
            !call cpu_time(start)
            call PhasingManagementWindows
            !call cpu_time(finish)
            !print '("Time PhasingManagementWindows= ",f6.3," seconds.")',finish-start
#endif
        endif

        if (inputParams%restartOption==OPT_RESTART_PHASING) then
#if CLUSTER==1
            write(6,*) "Restart option 2 stops program before Phasing has finished"
#elif CLUSTER==0
            write(6,*) "Restart option 2 stops program after Phasing has been managed"
#endif
            stop
        endif
    endif

    ! print*, " "
    ! print*, " ","Phasing completed"

    ! This is not necessary, already output in subroutine PhasingManagement
    if ((inputParams%restartOption/=OPT_RESTART_ALL).and.(inputParams%restartOption<OPT_RESTART_IMPUTATION)) then
        write(6,*) "Restart option 2 stops program after Phasing has been managed"
        stop
    endif
endif

if (inputParams%hmmoption/=RUN_HMM_NGS) then
    ! If we only want to phase data, then skip all the imputation steps
    if (inputParams%PhaseTheDataOnly==0) Then
        !call cpu_time(start)
        call ImputationManagement
        !call cpu_time(finish)
        !print '("Time ImputationManagement= ",f6.3," seconds.")',finish-start

#ifdef DEBUG
        write(0,*) 'DEBUG: Write results'
#endif

        call WriteOutResults

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
            !call cpu_time(start)
            call FinalChecker
            !call cpu_time(finish)
            !print '("Time FinalChecker= ",f6.3," seconds.")',finish-start
        endif
        ! call Cleaner
    endif
endif
!call cpu_time(start)
call PrintTimerTitles
!call cpu_time(finish)
!print '("Time call PrintTimerTitles= ",f6.3," seconds.")',finish-start

if (inputParams%restartOption > OPT_RESTART_IMPUTATION) then
    call system(RM // " Tmp2345678.txt")
end if

end program AlphaImpute
