!#############################################################################################################################################################################################################################
subroutine GeneProbManagementWindows
use Global
implicit none

integer :: i,JobsDone(nProcessors)
character(len=300) :: filout,f
logical :: FileExists

open (unit=109,file="TempGeneProb.BAT",status="unknown")

print*, " "
print*, " ","       Calculating genotype probabilities"

! Create bash script for run GeneProb subprocesses
do i=1,nProcessors
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
call CheckGeneProbFinished(nProcessors)

end subroutine GeneProbManagementWindows

!#############################################################################################################################################################################################################################

subroutine PhasingManagementWindows
use Global
implicit none

integer :: i,j,JobsDone(nPhaseInternal),StartJob,Tmp,StartNewJob,ProcUsed,JobsStarted(nPhaseInternal)
character(len=300) :: filout,infile,f
logical :: FileExists

print*, " "
print*, " ","       Performing the phasing of the data"

if (nProcessors<nPhaseInternal) then
    print*, "ERROR - To use this Restart option you need as many processors as phasing internal jobs"
    stop
endif

#ifdef OS_UNIX
#else
do i=1,nProcessors
    write (f,'("Phasing\Phase"i0,"\PhasingResults\Timer.txt")')i
    inquire(file=trim(f),exist=FileExists)
    if (FileExists .eqv. .true.) then
        write (f,'("Phasing\Phase"i0"\PhasingResults")')i
        call system("rmdir /s /q " // f // " >NUL")
    endif
enddo
#endif


open (unit=107,file="TempPhase1.BAT",status="unknown")
JobsStarted=0
ProcUsed=0
do i=1,nProcessors
    ProcUsed=ProcUsed+1
    write (infile,'("cd Phasing\Phase"i0)')i
    write (107,*) trim(infile)
    if (AlphaPhasePresent==0) write (107,*) "start /b AlphaPhase1.1.exe > out 2>&1"
    if (AlphaPhasePresent==1) write (107,*) "start /b .\AlphaPhase1.1.exe > out 2>&1"
    write (107,*) "cd ..\.."
    JobsStarted(i)=1
    if (ProcUsed==nPhaseInternal) exit
enddo
StartJob=ProcUsed
close (107)
call system("start ""Phasing"" .\TempPhase1.BAT >NUL")
Tmp=nProcessors
JobsDone(:)=0

! Check that every process has finished before go on
do
    do i=1,nPhaseInternal
        write (filout,'("Phasing\Phase"i0,"\PhasingResults\Timer.txt")')i

        inquire(file=trim(filout),exist=FileExists)
        if ((FileExists .eqv. .true.).and.(JobsDone(i)==0)) then
            print *," ","      AlphaPhase job ",i," done"
            JobsDone(i)=1
            if ((sum(JobsStarted(:))<nPhaseInternal).and.(sum(JobsDone(:))<nPhaseInternal)) then
                Tmp=Tmp+1
                JobsStarted(Tmp)=1
                write (filout,'("TempPhase"i0,".BAT")')Tmp
                open (unit=107,file=trim(filout),status="unknown")
                write (infile,'("cd Phasing\Phase"i0)')Tmp
                write (107,*) trim(infile)
                if (AlphaPhasePresent==0) write (107,*) "start /b AlphaPhase1.1.exe > out 2>&1"
                if (AlphaPhasePresent==1) write (107,*) "start /b .\AlphaPhase1.1.exe > out 2>&1"
                close(107)
                call system("start """" .\" // filout)
            endif
        endif
    enddo
    call sleep(SleepParameter)
    if (sum(JobsDone(:))==nPhaseInternal) exit
enddo
call system("del TempPhase*.BAT")

end subroutine PhasingManagementWindows 

!#############################################################################################################################################################################################################################
