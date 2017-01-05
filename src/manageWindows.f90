!#############################################################################################################################################################################################################################
subroutine GeneProbManagementWindows
    use Global
    use alphaimputeModule
    use AlphaImputeInMod
    implicit none

    integer :: i
    character(len=300) :: filout
    type(AlphaImputeInput), pointer :: inputParams

    inputParams => defaultInput
    open (unit=109,file="TempGeneProb.BAT",status="unknown")

    print*, " "
    print*, " ","       Calculating genotype probabilities"

    ! Create bash script for run GeneProb subprocesses
    do i=1,inputParams%nProcessors
        write (filout,'("cd GeneProb\GeneProb"i0)')i
        write (109,*) trim(filout)
        ! Call the external package GeneProbForAlphaImpute
        if (GeneProbPresent==0) write (109,*) "start /b GeneProbForAlphaImpute.exe > out 2>&1"
        if (GeneProbPresent==1) write (109,*) "start /b .\GeneProbForAlphaImpute.exe > out 2>&1"
        write (109,*) "cd ../.."
    enddo

    close (109)
    call system("start ""GeneProbs"" .\TempGeneProb.BAT >NUL")

    ! Check that every process has finished before going on
    call CheckGeneProbFinished(inputParams%nProcessors)

end subroutine GeneProbManagementWindows

!#############################################################################################################################################################################################################################

subroutine PhasingManagementWindows
    use Global
    use AlphaImputeInMod
    implicit none

    type(AlphaImputeInput), pointer :: inputParams
    integer :: i,StartJob,Tmp,ProcUsed
    integer, dimension(:), allocatable :: JobsDone, JobsStarted
    character(len=300) :: filout,infile
    logical :: FileExists

    inputParams => defaultInput
    print*, " "
    print*, " ","       Performing the phasing of the data"

    allocate(JobsDone(inputParams%nPhaseInternal))
    allocate(JobsStarted(inputParams%nPhaseInternal))
    if (inputParams%nProcessors<inputParams%nPhaseInternal) then
        print*, "ERROR - To use this Restart option you need as many processors as phasing internal jobs"
        stop
    endif

#ifdef OS_UNIX
#else
    do i=1,inputParams%nProcessors
        write (filout,'("Phasing\Phase"i0,"\PhasingResults\Timer.txt")')i
        inquire(file=trim(filout),exist=FileExists)
        if (FileExists .eqv. .true.) then
            write (filout,'("Phasing\Phase"i0"\PhasingResults")')i
            call system("rmdir /s /q " // filout // " >NUL")
        endif
    enddo
#endif


    open (unit=107,file="TempPhase1.BAT",status="unknown")
    JobsStarted=0
    ProcUsed=0
    do i=1,inputParams%nProcessors
        ProcUsed=ProcUsed+1
        write (infile,'("cd Phasing\Phase"i0)')i
        write (107,*) trim(infile)
        if (AlphaPhasePresent==0) write (107,*) "start /b AlphaPhase.exe > out 2>&1"
        if (AlphaPhasePresent==1) write (107,*) "start /b .\AlphaPhase.exe > out 2>&1"
        write (107,*) "cd ..\.."
        JobsStarted(i)=1
        if (ProcUsed==inputParams%nPhaseInternal) exit
    enddo
    StartJob=ProcUsed
    close (107)
    call system("start ""Phasing"" .\TempPhase1.BAT >NUL")
    Tmp=inputParams%nProcessors
    JobsDone(:)=0

    ! Check that every process has finished before go on
    do
        do i=1,inputParams%nPhaseInternal
            write (filout,'("Phasing\Phase"i0,"\PhasingResults\Timer.txt")')i

            inquire(file=trim(filout),exist=FileExists)
            if ((FileExists .eqv. .true.).and.(JobsDone(i)==0)) then
                print *," ","      AlphaPhase job ",i," done"
                JobsDone(i)=1
                if ((sum(JobsStarted(:))<inputParams%nPhaseInternal).and.(sum(JobsDone(:))<inputParams%nPhaseInternal)) then
                    Tmp=Tmp+1
                    JobsStarted(Tmp)=1
                    write (filout,'("TempPhase"i0,".BAT")')Tmp
                    open (unit=107,file=trim(filout),status="unknown")
                    write (infile,'("cd Phasing\Phase"i0)')Tmp
                    write (107,*) trim(infile)
                    if (AlphaPhasePresent==0) write (107,*) "start /b AlphaPhase.exe > out 2>&1"
                    if (AlphaPhasePresent==1) write (107,*) "start /b .\AlphaPhase.exe > out 2>&1"
                    close(107)
                    call system("start """" .\\" // filout)
                endif
            endif
        enddo
        call sleep(SleepParameter)
        if (sum(JobsDone(:))==inputParams%nPhaseInternal) exit
    enddo
    call system("del TempPhase*.BAT")

    deallocate(JobsDone)
    deallocate(JobsStarted)
end subroutine PhasingManagementWindows

!#############################################################################################################################################################################################################################
