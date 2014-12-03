!#############################################################################################################################################################################################################################

module Global
implicit none


integer,parameter :: WindowsLinux=0     !If 1 then compile for Windows / If 0 then compile for Linux

integer,parameter :: TestVersion=0  !If 1 then this is a development version with intermediate checking, if 0 it is not

integer,parameter :: PicVersion=0   !If 1 then this is a PIC version with suitability for their system, if 0 it is not

integer,parameter :: SleepParameter=1!00

integer,parameter :: lengan=20,MissingGenotypeCode=9,OffspringFillMin=10
integer,parameter :: ImputeFromHDLibraryCountThresh=1,ImputeFromHDPhaseThresh=1
integer,parameter :: ImputeFromParentCountThresh=1,ImputeFromGrandParentCountThresh=1
real,parameter :: DisagreeThreshold=0.05,GeneProbThresh=0.99

character (len=300) :: PedigreeFile,PhasePath,TrueGenosFile,GenotypeFile,GenderFile

integer :: nAnisG,nAnisRawPedigree,nSnp,nAnisP,IntEditStat,nPhaseInternal,nPhaseExternal,OutOpt,SexOpt,HetGameticStatus,HomGameticStatus
integer :: nProcessors,nProcessGeneProb,nProcessAlphaPhase,ManagePhaseOn1Off0,CountRawGenos,InternalIterations,nAnisInGenderFile
integer :: nAgreeImputeHDLib,nAgreeParentPhaseElim,nAgreeInternalHapLibElim,MaxLeftRightSwitch,MinSpan,ConservativeHapLibImputation
integer :: TrueGenos1None0,nSnpRaw,nObsDataRaw,nAgreePhaseElim,nAgreeGrandParentPhaseElim,PreProcess,UseGP,BypassGeneProb,HMMOption
integer :: nSnpIterate,NoPhasing,AlphaPhasePresent,GeneProbPresent,PrePhased,UserDefinedHD,PedFreePhasing,PhaseTheDataOnly,RestartOption

real :: PercGenoForHD,PercSnpMiss,SecondPercGenoForHD,GenotypeErrorPhase,WellPhasedThresh

integer(kind=1),allocatable,dimension (:) :: SnpIncluded,RecIdHDIndex,GenderRaw,RecGender,IndivIsGenotyped
integer(kind=1),allocatable,dimension (:,:) :: Genos,TempGenos,TmpGenos,MSTermInfo
integer(kind=1),allocatable,dimension (:,:) :: ImputeGenos,SireDam
integer(kind=1),allocatable,dimension (:,:,:) :: ImputePhase,TmpPhase,GlobalWorkPhase
integer,allocatable :: RecPed(:,:),Setter(:),CoreAndTailLengths(:),CoreLengths(:),GpIndex(:,:),BaseAnimals(:),GlobalTmpCountInf(:,:)
integer,allocatable :: GlobalHmmID(:)
real,allocatable,dimension (:) :: Maf
real,allocatable,dimension (:,:) :: ProbImputeGenos
real,allocatable,dimension (:,:,:) :: ProbImputePhase
character*(lengan),allocatable :: GenotypeId(:),GenderId(:)

end module Global

!#############################################################################################################################################################################################################################


module GlobalPedigree
use Global
implicit none

real(kind=4),allocatable :: xnumrelmatHold(:)
integer :: NRMmem, shell, shellmax, shellWarning
integer,allocatable:: seqid(:),seqsire(:),seqdam(:),RecodeGenotypeId(:),passedorder(:)
character*(lengan),allocatable :: ped(:,:),Id(:),sire(:),dam(:)
integer :: GlobalExtraAnimals       !Change John Hickey

end module GlobalPedigree

!########################################################################################################################################################################

module GlobalVariablesHmmMaCH
implicit none

character(len=300) :: GenotypeFileName,CheckPhaseFileName,CheckGenoFileName
integer :: nIndHmmMaCH,GlobalRoundHmm,nSnpHmm
integer :: nHapInSubH,idum,nRoundsHmm,HmmBurnInRound
double precision :: Theta
integer(kind=1),allocatable,dimension(:,:) :: GenosHmmMaCH,SubH
integer(kind=1),allocatable,dimension(:,:,:) :: FullH
integer,allocatable,dimension(:) :: ErrorUncertainty,ErrorMatches,ErrorMismatches,Crossovers,GlobalHmmHDInd
double precision,allocatable,dimension(:) :: Thetas,Epsilon
double precision,allocatable,dimension(:,:) :: ForwardProbs
double precision,allocatable,dimension(:,:,:) :: Penetrance
real,allocatable,dimension (:,:) :: ProbImputeGenosHmm

end module GlobalVariablesHmmMaCH

!########################################################################################################################################################################
 
program AlphaImpute
! The program is intended to be run numerous times
! The first time AlphaImpute is run in order to calculate Genotype probabilities (GenoProb)
! The second time, it is run to phase the individuals
! Finally, AlphaImpute should be run to impute genotypes.
! RestartOptions handle this situation, so:
!   * RestartOption=0 => Passes through the whole process: GenoProb, Phasing and Imputing
!   * RestartOption=1 => Makes only GenoProb
!   * RestartOption=2 => Makes GenoProb and Phasing
!   * RestartOption=3 => Makes only Imputation. This implies AlphaImpute has to be run 
!                           already in order to get GenoProb done or Genotype Probabilities 
!                           Genotype Probabilities have to be edited by hand
!   * RestartOption=4 =>
use Global
implicit none

call Titles
call ReadInParameterFile
if (RestartOption<2) call MakeDirectories
call CountInData
call ReadInData
call CheckParentage
call FillInSnp
call FillInBasedOnOffspring
call InternalEdit
call MakeFiles


print*, " "
print*, " ","Data editing completed"
if (SexOpt==0) then
    if (BypassGeneProb==0) then
        if (RestartOption<2) call GeneProbManagement
        print*, " "
        print*, " ","Genotype probabilities calculated"
    endif   
endif

if (ManagePhaseOn1Off0==1) then
    if (RestartOption<3) call PhasingManagement
endif   
print*, " "
print*, " ","Phasing completed"

! This is not necessary, already output in subroutine PhasingManagement
if ((RestartOption/=0).and.(RestartOption<3)) then
    print*, "Restart option 2 stops program after Phasing has been managed"
    stop
endif

if (PhaseTheDataOnly==0) then
    call ImputationManagement
    call WriteOutResults
    call ModelRecomb
    if (TrueGenos1None0==1) call FinalChecker
    call Cleaner
endif

call PrintTimerTitles

end program AlphaImpute

!#############################################################################################################################################################################################################################

subroutine ImputationManagement
use Global
use GlobalPedigree
implicit none

integer :: i,j,loop,dum

allocate(SireDam(0:nAnisP,2))
SireDam=0
do i=1,nAnisP
    do j=1,2
        SireDam(RecPed(i,j+1),j)=1 ! SireDam tells if individuals are Sires or Dams
    enddo
enddo

! WARNING: Need to discuss this part of code with John. Nonsense going on here!

if (RestartOption==4) then
    allocate(ImputeGenos(0:nAnisP,nSnp))
    allocate(ImputePhase(0:nAnisP,nSnp,2))
else
    if (SexOpt==0) then
        ! Impute initial genotypes from calculated genotype probabilities
        if (BypassGeneProb==0) then
            allocate(ImputeGenos(0:nAnisP,nSnp))
            allocate(ImputePhase(0:nAnisP,nSnp,2))
            allocate(GlobalWorkPhase(0:nAnisP,nSnp,2))
            call ReReadGeneProbs
        else
            ! Phase in the homozygous case for the SEX CHROMOSOME
            ! WARNING: NOTHING IS DONE!!
            call InsteadOfReReadGeneProb
        endif   

        ! Get Genotype information
        call InitialiseArrays       ! This is similar to InsteadOfReReadGeneProb subroutine but allocating ImputePhase
        call GeneProbPhase          ! Recover and store information about which and how many alleles/SNPs have been genotyped/phased 
    else
        allocate(MSTermInfo(nAnisP,2))
        MSTermInfo=0
    endif   

    if (NoPhasing==1) then
        call BaseAnimalFillIn
        if (PrePhased==1) call ReadInPrePhasedData
        if (SexOpt==1) call EnsureHetGametic

        ! General imputation procedures
        call GeneralFillIn

        if (HMMOption==3) then
            call MaCHController

        else
            print*, " "
            print*, " ","Imputation of base animals completed"
            do loop=1,InternalIterations
                print*, " "
                print*, "Performing imputation loop",loop 

                call PhaseElimination
                if (SexOpt==1) call EnsureHetGametic
                call GeneralFillIn
                print*, " "
                print*, " ","Parent of origin assigmnent of high density haplotypes completed"  

                call ParentPhaseElimination
                if (SexOpt==1) call EnsureHetGametic
                call GeneralFillIn
                call RestrictedWorkLeftRight
                if (SexOpt==1) call EnsureHetGametic
                call GeneralFillIn
                print*, " "
                print*, " ","Imputation from high-density parents completed"

                call ImputeFromHDLibrary
                if (SexOpt==1) call EnsureHetGametic
                call GeneralFillIn
                call RestrictedWorkLeftRight
                if (SexOpt==1) call EnsureHetGametic
                call GeneralFillIn
                print*, " "
                print*, " ","Haplotype library imputation completed"

                call InternalParentPhaseElim
                if (SexOpt==1) call EnsureHetGametic
                call GeneralFillIn
                call RestrictedWorkLeftRight
                if (SexOpt==1) call EnsureHetGametic
                call GeneralFillIn
                print*, " "
                print*, " ","Internal imputation from parents haplotype completed"

                call InternalHapLibImputation
                if (SexOpt==1) call EnsureHetGametic
                call GeneralFillIn
                if (SexOpt==1) call EnsureHetGametic
                call RestrictedWorkLeftRight
                call GeneralFillIn
                print*, " "
                print*, " ","Internal haplotype library imputation completed"
            enddo
            call ManageWorkLeftRight
        endif
    endif

    if (SexOpt==1) call EnsureHetGametic
    call GeneralFillIn

    deallocate(GlobalWorkPhase)
endif


end subroutine ImputationManagement

!#############################################################################################################################################################################################################################

subroutine ReadInParameterFile
use Global
use GlobalPedigree
use GlobalVariablesHmmMaCH

implicit none
integer :: k,i,resid,Changer,nLines
character (len=300) :: dumC,IntEdit,PhaseDone,OutputOptions,PreProcessOptions,TempOpt,TempHetGameticStatus
character (len=300) :: UserDefinedHDAnimalsFile,PrePhasedAnimalFile,PedigreeFreePhasing,PhasingOnlyOptions
character (len=300) :: UseGeneProb,ConservHapLibImp,CharBypassGeneProb,TmpHmmOption

open (unit=1,file="AlphaImputeSpec.txt",status="old")

nLines=0
do
    read (1,*,iostat=k) dumC
    nLines=nLines+1
    if (k/=0) then
        nLines=nLines-1
        exit
    endif
enddo
rewind(1)

if (nLines/=24) then
    print*, "   ","There are some lines missing from AlphaImputeSpec.txt"
    print*, "   ","HINT - maybe you are using the Spec file from the beta version which is out of date"
    stop
endif

read (1,*) dumC,PedigreeFile
read (1,*) dumC,GenotypeFile
read (1,*) dumC,TempOpt
SexOpt=9
HetGameticStatus=9
HomGameticStatus=9
if (trim(TempOpt)=="Yes") then
    backspace(1)
    read (1,*) dumC,TempOpt,GenderFile,TempHetGameticStatus
    HetGameticStatus=9
    if (trim(TempHetGameticStatus)=="Male") then        ! Species  with heterogametic males
        HetGameticStatus=1                              ! My father is heterogametic
        HomGameticStatus=2                              ! My mother is homogametic
    endif   
    if (trim(TempHetGameticStatus)=="Female") then      ! Species with heterogametic females
        HetGameticStatus=2                              ! My mother is heterogametic
        HomGameticStatus=1                              ! My father is homogametic
    endif   
    if (HetGameticStatus==9) then
        print*, "Warning - heterogametic status is misspecified"
        stop
    endif
    SexOpt=1
endif   

!Not sex chrom
if (trim(TempOpt)=="No") then
    SexOpt=0
endif
if (SexOpt==9) then
    print*, "Warning - Sex chromosome status is misspecified"
    stop
endif


read (1,*) dumC,nSnp
if (nSnp>240000) then
    print*, "Contact John Hickey if you want to do more than 240,000 SNP"
    stop
endif

read (1,*) dumC,IntEdit
if (trim(IntEdit)=='Yes') IntEditStat=1
if (trim(IntEdit)=='No') IntEditStat=0
if ((trim(IntEdit)/='Yes').and.(trim(IntEdit)/='No')) then
    write (*,*) "Specify editing status properly"               ! write (*,*) "You fucking bastard! X-)
endif

OutOpt=9
if (IntEditStat==1) then
    read (1,*) dumC,PercGenoForHD,PercSnpMiss,SecondPercGenoForHD,OutputOptions
    if (trim(OutputOptions)=="AllSnpOut") OutOpt=1
    if (trim(OutputOptions)=="EditedSnpOut") OutOpt=0
    if (OutOpt==9) then
        print*, "Output options incorrectly specified"
        print*, "Beware!!!!! AlphaImpute is case sensitive"
        stop
    endif
else
    PercGenoForHD=90.0
    read (1,*) dumC
    OutOpt=1
endif
PercGenoForHD=PercGenoForHD/100
PercSnpMiss=PercSnpMiss/100
SecondPercGenoForHD=SecondPercGenoForHD/100

read (1,*) dumC,PhaseDone
NoPhasing=1
if (trim(PhaseDone)=="PhaseDone") then
    ManagePhaseOn1Off0=0
    rewind (1)
    do i=1,6
        read (1,*) dumC
    enddo    
    read (1,*) dumC,PhaseDone,PhasePath,nPhaseInternal
    NoPhasing=1
elseif (trim(PhaseDone)=="NoPhase") then
    NoPhasing=0
    ManagePhaseOn1Off0=0
    rewind (1)
    do i=1,6
        read (1,*) dumC
    enddo    
    read (1,*) dumC
else
    ManagePhaseOn1Off0=1
    rewind (1)
    do i=1,6
        read (1,*) dumC
    enddo
    read (1,*) dumC,nPhaseExternal
endif

if (trim(PhaseDone)/="PhaseDone" .and. trim(PhaseDone)/="NoPhase") then
    if (nPhaseExternal>40) then
            print*, "Too many phasing runs required, at most you can do 40"
            stop
    endif
    if (nPhaseExternal<2) then
            print*, "Not enough phasing runs required, you must do 2 at least, 10 is better"
            stop
    endif

    nPhaseInternal=2*nPhaseExternal
    allocate(CoreAndTailLengths(nPhaseExternal))
    allocate(CoreLengths(nPhaseExternal))
    read (1,*) dumC,CoreAndTailLengths(:)
    read (1,*) dumC,CoreLengths(:)
    read (1,*) dumC,PedigreeFreePhasing
    if (trim(PedigreeFreePhasing)=="No") then
        PedFreePhasing=0
    else
        if (trim(PedigreeFreePhasing)=="Yes") then
            PedFreePhasing=1
        else
            print*, "Stop - Pedigree free phasing option incorrectly specified"
            stop
        endif
    endif
    read (1,*) dumC,GenotypeErrorPhase
else
    read (1,*) dumC
    read (1,*) dumC
    read (1,*) dumC
    read (1,*) dumC
endif

!read (1,*) dumC,UseGeneProb
!if (trim(UseGeneProb)=='Yes') UseGP=1
!if (trim(UseGeneProb)=='No') UseGP=0
!if ((trim(UseGeneProb)/='Yes').and.(trim(UseGeneProb)/='No')) then
!   write (*,*) "Specify editing status properly" 
!endif

    
read (1,*) dumC,nProcessors
read (1,*) dumC,InternalIterations
read (1,*) dumC,PreProcessOptions
if (PreProcessOptions=="No") then
    PreProcess=0
else
    if (PreProcessOptions=="Yes") then
        PreProcess=1
    else
        print*, "Stop - Preprocess of data option incorrectly specified"
        stop
    endif
endif

read (1,*) dumC,PhasingOnlyOptions
if (PhasingOnlyOptions=="No") then
    PhaseTheDataOnly=0
else
    if (PhasingOnlyOptions=="Yes") then
        PhaseTheDataOnly=1
    else
        print*, "Stop - Phasing only option incorrectly specified"
        stop
    endif
endif

ConservativeHapLibImputation=-1
read (1,*) dumC,ConservHapLibImp
if (trim(ConservHapLibImp)=="No") then
    ConservativeHapLibImputation=0
endif
if (trim(ConservHapLibImp)=="Yes") then
        ConservativeHapLibImputation=1
endif
if (ConservativeHapLibImputation==-1) then
    print*, "ConservativeHaplotypeLibraryUse not correctly specified"
    stop
endif

read (1,*) dumC,WellPhasedThresh

read (1,*) dumC,UserDefinedHDAnimalsFile
if (UserDefinedHDAnimalsFile=="None") then
    UserDefinedHD=0
else
    UserDefinedHD=1
    open (unit=46,file=trim(UserDefinedHDAnimalsFile),status="old")
endif

read (1,*) dumC,PrePhasedAnimalFile
if (PrePhasedAnimalFile=="None") then
    PrePhased=0
else
    PrePhased=1
    open (unit=47,file=trim(PrePhasedAnimalFile),status="old")
endif

BypassGeneProb=-1
read (1,*) dumC,CharBypassGeneProb
if (trim(CharBypassGeneProb)=="No") then
    BypassGeneProb=0
endif
if (trim(CharBypassGeneProb)=="Yes") then
    BypassGeneProb=1
endif
if (BypassGeneProb==-1) then
    print*, "BypassGeneProb not correctly specified"
    stop
endif

read (1,*) dumC,RestartOption

read (1,*) dumC,TmpHmmOption
HMMOption=0
if (trim(TmpHmmOption)=='No') HMMOption=1
if (trim(TmpHmmOption)=='Yes') HMMOption=2
if (trim(TmpHmmOption)=='Only') HMMOption=3
if (HMMOption==0) then
    print*, "HMMOption not correctly specified"
    stop
endif
read (1,*) dumC,nHapInSubH,HmmBurnInRound,nRoundsHMM,idum

read (1,*) dumC,TrueGenosFile
if (TrueGenosFile=="None") then
    TrueGenos1None0=0
else
    TrueGenos1None0=1
endif


open (unit=2,file=trim(PedigreeFile),status="old")
open (unit=3,file=trim(GenotypeFile),status="old")
if (SexOpt==1) open (unit=4,file=trim(GenderFile),status="old")

nProcessAlphaPhase=nProcessors-nProcessGeneProb ! Never used!

if (nPhaseInternal==2) then
    nAgreeImputeHDLib=1
    nAgreeParentPhaseElim=1
    nAgreePhaseElim=1
    nAgreeInternalHapLibElim=1
endif
if (nPhaseInternal==4) then
    nAgreeImputeHDLib=2
    nAgreeParentPhaseElim=2
    nAgreePhaseElim=2
    nAgreeInternalHapLibElim=2
endif
if (nPhaseInternal==6) then
    nAgreeImputeHDLib=3
    nAgreeParentPhaseElim=3
    nAgreePhaseElim=3
    nAgreeInternalHapLibElim=3
endif
if (nPhaseInternal>6) then
    nAgreeImputeHDLib=4
    nAgreeParentPhaseElim=4
    nAgreeGrandParentPhaseElim=4
    nAgreePhaseElim=4
    nAgreeInternalHapLibElim=4
endif

GlobalExtraAnimals=0

nSnpRaw=nSnp

!$  CALL OMP_SET_NUM_THREADS(nProcessors)

end subroutine ReadInParameterFile

!#############################################################################################################################################################################################################################

subroutine MaCHController
use Global
use GlobalPedigree
use GlobalVariablesHmmMaCH

implicit none

integer :: i,j


call ParseMaCHData
call SetUpEquations

open (unit=6,form='formatted',CARRIAGECONTROL='FORTRAN') 
do GlobalRoundHmm=1,nRoundsHmm
    write(6, 100) "   HMM Round   ",GlobalRoundHmm
    100 format ('+', a17,i10)    
    
    do j=1,nIndHmmMaCH
        call MaCHForInd(j)
    enddo
enddo
close (6)

ProbImputeGenosHmm=ProbImputeGenosHmm/(nRoundsHmm-HmmBurnInRound)

end subroutine MaCHController

!#############################################################################################################################################################################################################################

subroutine ParseMaCHData
use Global
use GlobalPedigree
use GlobalVariablesHmmMaCH

implicit none
integer :: i,j,k

nSnpHmm=nSnp
nIndHmmMaCH=nAnisG
allocate(GenosHmmMaCH(nIndHmmMaCH,nSnp))
allocate(GlobalHmmID(nIndHmmMaCH))
allocate(GlobalHmmHDInd(nIndHmmMaCH))
allocate(ProbImputeGenosHmm(nIndHmmMaCH,nSnp))
ProbImputeGenosHmm=0.0

GlobalHmmHDInd=0
k=0
do i=1,nAnisP
    if (IndivIsGenotyped(i)==1) then
        k=k+1
        GenosHmmMaCH(k,:)=ImputeGenos(i,:)
        GlobalHmmID(k)=i
        if ((float(count(GenosHmmMaCH(k,:)==9))/nSnp)<0.10) then
            GlobalHmmHDInd(k)=1
        endif
        do j=1,nSnp
            if ((GenosHmmMaCH(k,j)<0).or.(GenosHmmMaCH(k,j)>2)) GenosHmmMaCH(k,j)=3
        enddo   
    endif   
enddo
if (k/=nAnisG) then
    print*, "Error in ParseMaCHData"
    stop
endif

if (nHapInSubH>2*sum(GlobalHmmHDInd(:))) then
    print*, "Data set is too small for the number of Haplotypes in Sub H specified"
    stop
endif

end subroutine ParseMaCHData

!########################################################################################################################################################################
 
subroutine MaCHForInd(CurrentInd)
use GlobalVariablesHmmMaCH
implicit none

integer :: HapCount,CurrentInd,ShuffleInd1,ShuffleInd2
integer,allocatable,dimension(:) :: Shuffle1,Shuffle2

!Extract SubH
allocate(Shuffle1(nIndHmmMaCH))
allocate(Shuffle2(nIndHmmMaCH))
call RandomOrder(Shuffle1,nIndHmmMaCH,idum)
call RandomOrder(Shuffle2,nIndHmmMaCH,idum)

HapCount=0
ShuffleInd1=0
ShuffleInd2=0
do while (HapCount<nHapInSubH)
    if (mod(HapCount,2)==0) then
        ShuffleInd1=ShuffleInd1+1
        if ((Shuffle1(ShuffleInd1)/=CurrentInd).and.(GlobalHmmHDInd(ShuffleInd1)==1)) then
            HapCount=HapCount+1
            SubH(HapCount,:)=FullH(Shuffle1(ShuffleInd1),:,1)
        endif   
    else
        ShuffleInd2=ShuffleInd2+1   
        if ((Shuffle2(ShuffleInd2)/=CurrentInd).and.(GlobalHmmHDInd(ShuffleInd2)==1)) then
            HapCount=HapCount+1
            SubH(HapCount,:)=FullH(Shuffle2(ShuffleInd2),:,2)       
        endif
    endif
enddo

allocate(ForwardProbs(nHapInSubH*nHapInSubH,nSnpHmm))

call ForwardAlgorithm(CurrentInd)
call SampleChromosomes(CurrentInd)

if (GlobalRoundHmm>HmmBurnInRound) ProbImputeGenosHmm(CurrentInd,:)=ProbImputeGenosHmm(CurrentInd,:)+FullH(CurrentInd,:,1)+FullH(CurrentInd,:,2)

deallocate(ForwardProbs)


end subroutine MaCHForInd

!########################################################################################################################################################################

subroutine SampleChromosomes(CurrentInd)
use GlobalVariablesHmmMaCH
implicit none

integer :: i,j,k,l,SuperJ,CurrentInd,Index,OffOn,State1,State2,TmpJ,TopBot,FirstState,SecondState,Tmp
double precision :: Summer,ran1,Choice,Sum00,Sum01,Sum10,Sum11

Summer=0.0
Index=0
do i=1,nHapInSubH
    do j=1,i
        Index=Index+1
        Summer=Summer+ForwardProbs(Index,nSnpHmm)
    enddo
enddo   

Choice=ran1(idum)*Summer    

Summer=0.0
Index=0
OffOn=0
do i=1,nHapInSubH
    do j=1,i
        Index=Index+1
        Summer=Summer+ForwardProbs(Index,nSnpHmm)
        if (Summer>Choice) then
            State1=i
            State2=j
            OffOn=1
            exit
        endif
    enddo
    if (OffOn==1) exit
enddo   

SuperJ=nSnpHmm
do while (SuperJ>1)
    SuperJ=SuperJ-1
    call ImputeAlleles(CurrentInd,SuperJ+1,State1,State2)
    TmpJ=SuperJ
    Theta=Thetas(SuperJ)
    do while ((GenosHmmMaCH(CurrentInd,SuperJ)==3).and.SuperJ>1)
        SuperJ=SuperJ-1
        Theta=Theta+Thetas(SuperJ)-Theta*Thetas(SuperJ)
    enddo
        
    !When examining the previous location we consider three alternatives:
    !states that could be reached when both haplotypes recombine (11),
    !states that can be reached when the first (10) or second (01) haplotype recombines,
    !and the states that can be reached without recombination.
    Sum00=0.0
    Sum01=0.0
    Sum10=0.0
    Sum11=0.0   

    Index=0
    do k=1,nHapInSubH
        do l=1,k
            Index=Index+1
            Sum11=Sum11+ForwardProbs(Index,SuperJ)
            if ((State1==k).or.(State1==l)) Sum01=Sum01+ForwardProbs(Index,SuperJ)  
            if ((State2==k).or.(State2==l)) Sum10=Sum10+ForwardProbs(Index,SuperJ)  
            if (((State1==k).and.(State2==l)).or.((State1==l).and.(State2==k))) Sum00=Sum00+ForwardProbs(Index,SuperJ)  
        enddo
    enddo   

    Summer=Sum11*Theta*Theta/(nHapInSubH*nHapInSubH)+(Sum10+Sum01)*Theta*(1.0-Theta)/nHapInSubH+Sum00*(1.0-Theta)*(1.0-Theta)
    if (SuperJ==1) exit
    

    !Sample number and decide how many state changes occurred between the
    !two positions
    Choice=ran1(idum)*Summer

    !The most likely outcome is that no changes occur ...
    Choice=Choice-(Sum00*(1.0-Theta)*(1.0-Theta))

    if (Choice<=0.0) then
        !Record outcomes for intermediate, uninformative, positions
        TopBot=1
        call FillPath(CurrentInd,SuperJ,TmpJ+1,State1,TopBot)
        TopBot=2
        call FillPath(CurrentInd,SuperJ,TmpJ+1,State2,TopBot)
        cycle
    endif

    !But perhaps the first or second haplotype recombined
    Choice=Choice-(Sum10*Theta*(1.0-Theta)/nHapInSubH)

    if (Choice<=0.0) then
        !The first haplotype changed ...
        Choice=Choice*nHapInSubH/(Theta*(1.0-Theta))
        !Record the original state
        FirstState=State1                       
        do while (State1<nHapInSubH)
            State1=State1+1                                                         !Check with Carl
            if (State1>=State2) then
                Choice=Choice+ForwardProbs(State1*(State1-1)/2+State2,SuperJ)
            else
                Choice=Choice+ForwardProbs(State2*(State2-1)/2+State1,SuperJ)           
            endif
            if (Choice>=0.0) exit
        enddo

        !Record outcomes for intermediate, uninformative, positions
        TopBot=1  
        call SamplePath(CurrentInd,SuperJ,TmpJ+1,State1,FirstState,TopBot) 
        TopBot=2
        call FillPath(CurrentInd,SuperJ,TmpJ+1,State2,TopBot)      
        cycle
    endif

    Choice=Choice-(Sum01*Theta*(1.0-Theta)/nHapInSubH)

    if (Choice<=0.0) then
        !The second haplotype changed ...
        Choice=Choice*nHapInSubH/(Theta*(1.0-Theta))
        !Save the original state
        SecondState=State2
        do while (State2<nHapInSubH)
            State2=State2+1                                                         !Check with Carl
            if (State1>=State2) then
                Choice=Choice+ForwardProbs(State1*(State1-1)/2+State2,SuperJ)       
            else
               Choice=Choice+ForwardProbs(State2*(State2-1)/2+State1,SuperJ)
            endif    
            if (Choice>=0.0) exit
        enddo                                                                       !Check with Carl                            

        !Record outcomes for intermediate, uninformative, positions
        TopBot=1    
        call FillPath(CurrentInd,SuperJ,TmpJ+1,State1,TopBot)                   !Check with Carl
        TopBot=2            
        call SamplePath(CurrentInd,SuperJ,TmpJ+1,State2,SecondState,TopBot)     !Check with Carl    
        cycle

    endif

    !Try to select any other state
    Choice=Choice*nHapInSubH*nHapInSubH/(Theta*Theta)

    !Save the original states
    FirstState=State1
    SecondState=State2

    Summer=0.0
    Index=0
    OffOn=0
    do i=1,nHapInSubH
        do j=1,i
            Index=Index+1
            Summer=Summer+ForwardProbs(Index,nSnpHmm)
            if (Summer>Choice) then
                State1=i
                State2=j
                OffOn=1
                exit
            endif
        enddo
        if (OffOn==1) exit
    enddo   
    
    if (ran1(idum)>0.5) then
        Tmp=State1
        State2=State1
        State2=Tmp
    endif
    !Record outcomes for intermediate, uninformative, positions
    TopBot=1    
    call SamplePath(CurrentInd,SuperJ,TmpJ+1,State1,FirstState,TopBot) 
    TopBot=2
    call SamplePath(CurrentInd,SuperJ,TmpJ+1,State2,SecondState,TopBot)     

enddo

call ImputeAlleles(CurrentInd,1,State1,State2)

end subroutine SampleChromosomes

!########################################################################################################################################################################
 
subroutine FillPath(CurrentInd,FromMarker,ToMarker,State,TopBot)
use GlobalVariablesHmmMaCH
implicit none

integer :: j,CurrentInd,FromMarker,ToMarker,State,TopBot

do j=FromMarker+1,ToMarker-1
    call ImputeAllele(CurrentInd,j,State,TopBot)
enddo

end subroutine FillPath

!########################################################################################################################################################################
 
subroutine SamplePath(CurrentInd,FromMarker,ToMarker,FromState,ToState,TopBot)
use GlobalVariablesHmmMaCH
implicit none

integer :: i,j,CurrentInd,TopBot,FromMarker,ToMarker,FromState,ToState
double precision :: R,Theta1,ran1

Theta=0.0
do j=FromMarker,ToMarker-1
    !Calculate overall recombination fraction for the interval
    Theta=Thetas(j)+Theta-Theta*Thetas(j)
enddo

!Impute a path between the two end markers, assuming no genotypes
!are observed -- the only constraint is that we must start at
!fromState and end at toState with at least one intervening recombinant

do while (FromMarker<ToMarker-1)
    R=ran1(idum)*Theta
    Theta1=Thetas(FromMarker)
    if (Theta < 0.9) then
        !Fast closed formula
        Theta=(Theta-Theta1)/(1.0-Theta1)
    else
        Theta = 0.0
        !More accurate, iterative formula
        do i=FromMarker+1,ToMarker-1
            Theta=Thetas(i)+Theta-Theta*Thetas(i)
        enddo
    endif

    if (R>theta1) then
        !No recombinant in the in first interval
        FromMarker=FromMarker+1
        call ImputeAllele(CurrentInd,FromMarker,ToState,TopBot)                 
        cycle
    endif
    Crossovers(FromMarker)=Crossovers(FromMarker)+1

    if (R<Theta1*(1.0-Theta)) then
        !No recombinant in the second interval
        call FillPath(CurrentInd,FromMarker,ToMarker,ToState,TopBot);
        return
    else
    
        FromMarker=FromMarker+1
        ToState=int(ran1(idum)*nHapInSubH)+1
        call ImputeAllele(CurrentInd,FromMarker,ToState,TopBot)                 
    endif

enddo

!If we get here, record obligate recombinant between two consecutive markers
Crossovers(FromMarker)=Crossovers(FromMarker)+1


end subroutine SamplePath

!########################################################################################################################################################################
 
subroutine ImputeAlleles(CurrentInd,CurrentMarker,State1,State2)
use GlobalVariablesHmmMaCH
implicit none

integer :: CurrentInd,CurrentMarker,State1,State2,Imputed1,Imputed2,Genotype,Differences
double precision :: ran1


Imputed1=SubH(State1,CurrentMarker)
Imputed2=SubH(State2,CurrentMarker)

Genotype=GenosHmmMaCH(CurrentInd,CurrentMarker)

if ((Genotype/=0).and.(Genotype/=2)) then
    FullH(CurrentInd,CurrentMarker,1)=Imputed1
    FullH(CurrentInd,CurrentMarker,2)=Imputed2
endif

if (Genotype==3) return

Differences=abs(Genotype - (Imputed1+Imputed2))

if ((Genotype==1).and.(Differences==0)) then
    ErrorUncertainty(CurrentMarker)=ErrorUncertainty(CurrentMarker)+1
else
    ErrorMatches(CurrentMarker)=ErrorMatches(CurrentMarker)+(2-Differences)
    ErrorMismatches(CurrentMarker)=ErrorMismatches(CurrentMarker)+Differences
endif 

if (Genotype/=1) return

if (Imputed1==Imputed2) then
    if (ran1(idum)>=0.5) then
        FullH(CurrentInd,CurrentMarker,1)=abs(Imputed1-1)       
    else    
        FullH(CurrentInd,CurrentMarker,1)=abs(Imputed2-1)                   
    endif
endif

end subroutine ImputeAlleles

!########################################################################################################################################################################
 
subroutine ImputeAllele(CurrentInd,CurrentMarker,State,TopBot)
use GlobalVariablesHmmMaCH
implicit none

integer :: CurrentInd,CurrentMarker,State,TopBot

FullH(CurrentInd,CurrentMarker,TopBot)=SubH(State,CurrentMarker)

end subroutine ImputeAllele

!########################################################################################################################################################################
 
subroutine ForwardAlgorithm(CurrentInd)
use GlobalVariablesHmmMaCH
implicit none

integer :: i,j,CurrentInd,PrecedingMarker

call SetUpPrior

j=1
call ConditionOnData(CurrentInd,j)

Theta=0.0
PrecedingMarker=1
do j=2,nSnpHmm
    Theta=Theta+Thetas(j-1)-Theta*Thetas(j-1)
    if ((GenosHmmMaCH(CurrentInd,j)/=3).or.(j==nSnpHmm)) then
        call Transpose(j,PrecedingMarker)
        call ConditionOnData(CurrentInd,j)
        PrecedingMarker=j
        Theta=0.0
    endif
enddo

end subroutine ForwardAlgorithm

!########################################################################################################################################################################
 
subroutine Transpose(CurrentMarker,PrecedingMarker)
use GlobalVariablesHmmMaCH
implicit none

integer :: i,j,CurrentMarker,PrecedingMarker,Index
double precision :: Summer,Marginals(nHapInSubH),NoChange,OneChange,TwoChange

if (Theta==0.0) then
    ForwardProbs(:,CurrentMarker)=ForwardProbs(:,PrecedingMarker)
else
    Summer=0.0
    Index=0
    Marginals(:)=0.0
    do i=1,nHapInSubH
        do j=1,i-1      !Karl suggestion
            Index=Index+1
            Summer=Summer+ForwardProbs(Index,PrecedingMarker)
            Marginals(i)=Marginals(i)+ForwardProbs(Index,PrecedingMarker)
            Marginals(j)=Marginals(j)+ForwardProbs(Index,PrecedingMarker)           
        enddo
        Index=Index+1
        Summer=Summer+ForwardProbs(Index,PrecedingMarker)
        Marginals(i)=Marginals(i)+(ForwardProbs(Index,PrecedingMarker)*2.0)     
    enddo       
    
    NoChange=(1.0-Theta)*(1.0-Theta)
    OneChange=(1.0-Theta)*Theta/nHapInSubH
    TwoChange=Summer*Theta*Theta/(nHapInSubH*nHapInSubH)
    
    !Automatically rescale likelihoods when they get too small
    if (Summer < 1e-15) then
        NoChange=NoChange*1e30
        OneChange=OneChange*1e30
        TwoChange=TwoChange*1e30
    endif
    
    !This final loop actually transposes the probabilities for each state
    Index=0
    do i=1,nHapInSubH
        do j=1,i-1      !Karl suggestion
            Index=Index+1
            ForwardProbs(Index,CurrentMarker)=(ForwardProbs(Index,PrecedingMarker)*NoChange)+(Marginals(i)*OneChange)+(Marginals(j)*OneChange)+(2*TwoChange)
        enddo
        Index=Index+1
        ForwardProbs(Index,CurrentMarker)=(ForwardProbs(Index,PrecedingMarker)*NoChange)+(Marginals(i)*OneChange)+(2*TwoChange)
    enddo       
endif

end subroutine Transpose

!########################################################################################################################################################################

 
subroutine ConditionOnData(CurrentInd,Marker)
use GlobalVariablesHmmMaCH
implicit none

integer :: i,j,CurrentInd,Marker,Index
double precision :: Factors(0:1)

if (GenosHmmMaCH(CurrentInd,Marker)==3) then
    return
else
    Index=0
    do i=1,nHapInSubH
        Factors(0)=Penetrance(Marker,SubH(i,Marker),GenosHmmMaCH(CurrentInd,Marker))
        Factors(1)=Penetrance(Marker,SubH(i,Marker)+1,GenosHmmMaCH(CurrentInd,Marker))
        do j=1,i
            Index=Index+1
            ForwardProbs(Index,Marker)=ForwardProbs(Index,Marker)*Factors(SubH(j,Marker))
        enddo       
    enddo
endif       

end subroutine ConditionOnData

!########################################################################################################################################################################
 
subroutine CalcPenetrance
use GlobalVariablesHmmMaCH
implicit none

integer :: j

allocate(Penetrance(nSnpHmm,0:2,0:2))

do j=1,nSnpHmm
    Penetrance(j,0,0)=(1.0-Epsilon(j))**2
    Penetrance(j,0,1)=2.0*(1.0-Epsilon(j))*Epsilon(j)   
    Penetrance(j,0,2)=(Epsilon(j)**2)       
    Penetrance(j,1,0)=(1.0-Epsilon(j))*Epsilon(j)   
    Penetrance(j,1,1)=((1.0-Epsilon(j))**2)+(Epsilon(j)**2)
    Penetrance(j,1,2)=(1.0-Epsilon(j))*Epsilon(j)   
    Penetrance(j,2,0)=(Epsilon(j)**2)   
    Penetrance(j,2,1)=2.0*(1.0-Epsilon(j))*Epsilon(j)   
    Penetrance(j,2,2)=(1.0-Epsilon(j))**2   
enddo

end subroutine CalcPenetrance

!########################################################################################################################################################################
 
subroutine SetUpPrior
use GlobalVariablesHmmMaCH
implicit none

ForwardProbs(:,1)=1.0/(nHapInSubH*nHapInSubH)

end subroutine SetUpPrior

!########################################################################################################################################################################

subroutine SetUpEquations
use GlobalVariablesHmmMaCH
implicit none

integer :: i,j
double precision :: ran1

allocate(FullH(nIndHmmMaCH,nSnpHmm,2))
allocate(SubH(nHapInSubH,nSnpHmm))
allocate(Epsilon(nSnpHmm))
allocate(Thetas(nSnpHmm-1))
allocate(ErrorUncertainty(nSnpHmm))
allocate(ErrorMatches(nSnpHmm))
allocate(ErrorMismatches(nSnpHmm))
allocate(Crossovers(nSnpHmm))

Epsilon=0.00000001
Thetas=0.01

!Initialise FullH
do i=1,nIndHmmMaCH
    do j=1,nSnpHmm

        if (GenosHmmMaCH(i,j)==0) then
            FullH(i,j,:)=0
        endif
        
        if (GenosHmmMaCH(i,j)==2) then
            FullH(i,j,:)=1
        endif   
        
        if (GenosHmmMaCH(i,j)==1) then
            if (ran1(idum)>=0.5) then           
                FullH(i,j,1)=0
                FullH(i,j,2)=1
            else
                FullH(i,j,1)=1
                FullH(i,j,2)=0
            endif
        endif   
        if (GenosHmmMaCH(i,j)==3) then
            if (ran1(idum)>=0.5) then           
                FullH(i,j,1)=0
            else    
                FullH(i,j,1)=1
            endif   
            if (ran1(idum)>=0.5) then           
                FullH(i,j,2)=0
            else    
                FullH(i,j,2)=1
            endif   
        endif   

    enddo        
enddo

call CalcPenetrance

ErrorUncertainty(:)=0
ErrorMatches(:)=0
ErrorMismatches(:)=0
Crossovers(:)=0

end subroutine SetUpEquations

!########################################################################################################################################################################

subroutine ReadInPrePhasedData
! Impute phase information from pre-phased file. Count the number of pre-phased individuals
use Global
use GlobalPedigree

integer :: h,i,j,k,nAnisPrePhased,WorkPhase(nSnpRaw,2),CountPrePhased
character(len=300) :: dumC

! Count animals prephased
nAnisPrePhased=0
do
    read (47,*,iostat=k) dumC
    nAnisPrePhased=nAnisPrePhased+1
    if (k/=0) then
        nAnisPrePhased=nAnisPrePhased-1
        exit
    endif
enddo
rewind(47)
nAnisPrePhased=nAnisPrePhased/2         ! Two haplotypes per animal have been read

CountPrePhased=0
do k=1,nAnisPrePhased
    read (47,*) dumC,WorkPhase(:,1)     ! Paternal haplotype
    read (47,*) dumC,WorkPhase(:,2)     ! Maternal haplotype
    do i=1,nAnisP
        if (trim(dumC)==trim(Id(i))) then   ! Check if any animal in the file agrees with the animals in the pedigree
            h=0
            do j=1,nSnpRaw
                if (SnpIncluded(j)==1) then ! Check if this SNP has to be considered (may be it has been removed during the edition step)
                    h=h+1
                    ! Impute phase only if this locus is phased (in the PrePhased file)
                    if ((WorkPhase(j,1)==0).or.(WorkPhase(j,1)==1)) ImputePhase(i,h,1)=WorkPhase(j,1)
                    if ((WorkPhase(j,2)==0).or.(WorkPhase(j,2)==1)) ImputePhase(i,h,2)=WorkPhase(j,2)
                endif
            enddo
            CountPrePhased=CountPrePhased+1
            exit
        endif
    enddo   
enddo

print*, " "
print*, " ",CountPrePhased," valid pre-phased indiviudals in the user specified pre-phased file"

end subroutine ReadInPrePhasedData

!#############################################################################################################################################################################################################################

subroutine ManageWorkLeftRight
use Global

implicit none

MaxLeftRightSwitch=4; MinSpan=200
call WorkLeftRight
if (SexOpt==1) call EnsureHetGametic
call GeneralFillIn

MaxLeftRightSwitch=3; MinSpan=200
call WorkLeftRight
if (SexOpt==1) call EnsureHetGametic
call GeneralFillIn

MaxLeftRightSwitch=4; MinSpan=100
call WorkLeftRight
if (SexOpt==1) call EnsureHetGametic
call GeneralFillIn

MaxLeftRightSwitch=4; MinSpan=50
call WorkLeftRight
if (SexOpt==1) call EnsureHetGametic
call GeneralFillIn

end subroutine ManageWorkLeftRight

!#############################################################################################################################################################################################################################

subroutine GeneProbManagement
use Global
implicit none

integer :: i,JobsDone(nProcessors),JDone(nProcessors)
character(len=300) :: filout,f
logical :: FileExists

open (unit=109,file="TempGeneProb.sh",status="unknown")

if (PicVersion==0) then         ! GeneProbForAlphaImpute algorithm to calculate probabilities of genotype
    ! Create a bash script
    do i=1,nProcessors
        write (filout,'("cd GeneProb/GeneProb"i0)')i
        write (109,*) trim(filout)
        ! Call the external package GeneProbForAlphaImpute
        if (GeneProbPresent==0) write (109,*) "nohup sh -c ""GeneProbForAlphaImpute > out 2>&1"" >/dev/null &"
        if (GeneProbPresent==1) write (109,*) "nohup sh -c ""./GeneProbForAlphaImpute > out 2>&1"" >/dev/null &"    
        write (109,*) "cd ../.."
    enddo

    print*, " "
    print*, " ","       Calculating genotype probabilities"
    close (109)
    call system("chmod +x TempGeneProb.sh")
    call system("./TempGeneProb.sh")
    JobsDone(:)=0
    JDone(:)=0
    if (RestartOption/=1) then
        do  ! This should be a while do loop
            do i=1,nProcessors
                write (filout,'("./GeneProb/GeneProb"i0,"/GpDone.txt")')i
                inquire(file=trim(filout),exist=FileExists)
                if (FileExists .eqv. .true.) then
                    if (JDone(i)==0) print*, " ","      GeneProb job ",i," done"
                    JobsDone(i)=1
                    JDone(i)=1
                endif
            enddo
            call sleep(SleepParameter)
            if (sum(JobsDone(:))==nProcessors) exit
        enddo
    endif
else            ! PIC company algorithm to calculate probabilities of genotype
    write (filout,'("cd GeneProb/")')
    write(f,'(i0)') nProcessors
    write (109,*) trim(filout)
    write(109,*) "cp ../../../SharedFiles/AlphaImpute_scripts/runSubmitGeneProb.sh ."
    write(109,*) "cp ../../../SharedFiles/AlphaImpute_scripts/submitGeneProb.sh ."
    write(109,'(a,a)') "./runSubmitGeneProb.sh ",trim(f)
    print*, " "
    print*, " ","       Calculating genotype probabilities"
    close (109)
    call system("chmod +x TempGeneProb.sh")
    call system("./TempGeneProb.sh")
    call system("rm TempGeneProb.sh")

    JobsDone(:)=0
    if (RestartOption/=1) then
        do
            do i=1,nProcessors
                write (filout,'("./GeneProb/GeneProb"i0,"/GpDone.txt")')i
                inquire(file=trim(filout),exist=FileExists)
                if (FileExists .eqv. .true.) then
                    if (JDone(i)==0) print*, " ","      GeneProb job ",i," done"
                    JobsDone(i)=1
                endif
            enddo
            call sleep(SleepParameter)
            if (sum(JobsDone(:))==nProcessors) exit
        enddo
    endif   
endif

if (RestartOption==1) then
    print*, "Restart option 1 stops program after Geneprobs jobs have been submitted"
    stop
endif

end subroutine GeneProbManagement

!#############################################################################################################################################################################################################################

subroutine PhasingManagement
use Global
implicit none

integer :: i,j,JobsDone(nPhaseInternal),StartJob,Tmp,StartNewJob,ProcUsed,JobsStarted(nPhaseInternal)
character(len=300) :: filout,infile,f
logical :: FileExists

print*, " "
print*, " ","       Performing the phasing of the data"
if (PicVersion==0) then
    open (unit=107,file="TempPhase1.sh",status="unknown")
    JobsStarted=0
    ProcUsed=0
    do i=1,nProcessors
        ProcUsed=ProcUsed+1
        write (infile,'("cd Phasing/Phase"i0)')i
        write (107,*) trim(infile)
        if (AlphaPhasePresent==0) write (107,*) "nohup sh -c ""AlphaPhase1.1 > out 2>&1"" >/dev/null &"
        if (AlphaPhasePresent==1) write (107,*) "nohup sh -c ""./AlphaPhase1.1 > out 2>&1"" >/dev/null &"
        write (107,*) "cd ../.."
        JobsStarted(i)=1
        if (ProcUsed==nPhaseInternal) exit
    enddo
    StartJob=ProcUsed
    close (107)
    call system("chmod +x TempPhase*.sh")
    call system("./TempPhase1.sh")
    Tmp=nProcessors
    JobsDone(:)=0
    if (RestartOption==2) then
        if (nProcessors<nPhaseInternal) then
            print*, "ERROR - To use this Restart option you need as many processors as phasing internal jobs"
            stop
        endif
    else    
        do ! This should be a while do loop
            do i=1,nPhaseInternal
                write (filout,'("./Phasing/Phase"i0,"/PhasingResults/Timer.txt")')i
        
                inquire(file=trim(filout),exist=FileExists)     
                if ((FileExists .eqv. .true.).and.(JobsDone(i)==0)) then
                    print*, " ","AlphaPhase job ",i," done"
                    JobsDone(i)=1
                    if ((sum(JobsStarted(:))<nPhaseInternal).and.(sum(JobsDone(:))<nPhaseInternal)) then
                        Tmp=Tmp+1
                        JobsStarted(Tmp)=1
                        write (filout,'("TempPhase"i0,".sh")')Tmp
                        open (unit=107,file=trim(filout),status="unknown")
                        write (infile,'("cd Phasing/Phase"i0)')Tmp
                        write (107,*) trim(infile) 
                        if (AlphaPhasePresent==0) write (107,*) "nohup sh -c ""AlphaPhase1.1 > out 2>&1"" >/dev/null &"
                        if (AlphaPhasePresent==1) write (107,*) "nohup sh -c ""./AlphaPhase1.1 > out 2>&1"" >/dev/null &"
                        close(107)
                        call system("chmod +x TempPhase*.sh")
                        call system("./" // filout)
                    endif   
                endif
            enddo
            call sleep(SleepParameter)
            if (sum(JobsDone(:))==nPhaseInternal) exit 
        enddo
    endif       
    call system("rm TempPhase*.sh")
else
    open (unit=107,file="TempPhase1.sh",status="unknown")
    write (filout,'("cd Phasing/")')
    write(f,'(i0)') nPhaseInternal
    write (107,*) trim(filout)
    write(107,*) "cp ../../../SharedFiles/AlphaImpute_scripts/runSubmitPhasing.sh ."
    write(107,*) "cp ../../../SharedFiles/AlphaImpute_scripts/submitPhasing.sh ."
    write(107,'(a,a)') "./runSubmitPhasing.sh ",trim(f)
    print*, " "
    close (107)
    call system("chmod +x TempPhase1.sh")
    call system("./TempPhase1.sh")
    call system("rm TempPhase1.sh")

    if (RestartOption/=2) then
        JobsDone(:)=0
        do
            do i=1,nPhaseInternal
                write (filout,'("./Phasing/Phase"i0,"/PhasingResults/Timer.txt")')i
                inquire(file=trim(filout),exist=FileExists)
                if ((FileExists .eqv. .true.).and.(JobsDone(i)==0)) then
                    print*, " ","AlphaPhase job ",i," done"
                    JobsDone(i)=1
                endif
            enddo
            call sleep(SleepParameter)
            if (sum(JobsDone(:))==nPhaseInternal) exit
        enddo
    endif   
endif

if (RestartOption==2) then
    print*, "Restart option 2 stops program after AlphaPhase jobs have been submitted"
    stop
endif

end subroutine PhasingManagement 

!#############################################################################################################################################################################################################################

subroutine IterateGeneProbs
use Global

implicit none

integer :: h,i,j,k,dum,tmp,JobsDone(nProcessors),JDone(nProcessors),StSnp,EnSnp,OnOff
real,allocatable :: PatAlleleProb(:,:),MatAlleleProb(:,:),HetProb(:),GeneProbWork(:,:)
character(len=300) :: filout,f
logical :: FileExists

if (OutOpt==0) nSnpIterate=nSnp
if (OutOpt==1) nSnpIterate=nSnpRaw

allocate(PatAlleleProb(nSnpIterate,2))
allocate(MatAlleleProb(nSnpIterate,2))
allocate(HetProb(nSnpIterate))
allocate(GeneProbWork(nSnpIterate,4))
allocate(ProbImputeGenos(0:nAnisP,nSnpIterate))
allocate(ProbImputePhase(0:nAnisP,nSnpIterate,2))
deallocate(GpIndex)
allocate(GpIndex(nProcessors,2))

ProbImputeGenos(0,:)=0.0
ProbImputePhase(0,:,:)=0.0
ProbImputeGenos(1:nAnisP,:)=-9.0
ProbImputePhase(1:nAnisP,:,:)=-9.0

Tmp=int(float(nSnpIterate)/nProcessors)
GpIndex(1,1)=1
GpIndex(1,2)=Tmp
if (nProcessors>1) then
    do i=2,nProcessors
        GpIndex(i,1)=GpIndex(i-1,1)+Tmp     
        GpIndex(i,2)=GpIndex(i-1,2)+Tmp
    enddo
endif
GpIndex(nProcessors,2)=nSnpIterate


do i=1,nProcessors
    write (filout,'("./IterateGeneProb/GeneProb"i0,"/GeneProbSpec.txt")')i
    open (unit=108,file=trim(filout),status='unknown')
    write (108,*) "nAnis        ,",nAnisP 
    write (108,*) "nSnp     ,",nSnpIterate
    write (108,*) "InputFilePath    ,",'"../IterateGeneProbInput.txt"'
    write (108,*) "OutputFilePath   ,",'"GeneProbs.txt"'
    write (108,*) "StartSnp     ,",GpIndex(i,1) 
    write (108,*) "EndSnp       ,",GpIndex(i,2) 
    close(108)
    write (filout,'("GeneProb"i0)')i
    if(PicVersion==0) call system ("cp GeneProbForAlphaImpute IterateGeneProb/" // filout)
enddo   

open (unit=109,file="TempIterateGeneProb.sh",status="unknown")

if (RestartOption/=4) then
    if (PicVersion==0) then
        do i=1,nProcessors
            write (filout,'("cd IterateGeneProb/GeneProb"i0)')i
            write (109,*) trim(filout) 
            if (GeneProbPresent==0) write (109,*) "nohup sh -c ""GeneProbForAlphaImpute > out 2>&1"" >/dev/null &"
            if (GeneProbPresent==1) write (109,*) "nohup sh -c ""./GeneProbForAlphaImpute > out 2>&1"" >/dev/null &"
            write (109,*) "cd ../.."
            call flush(109)
        enddo
        close(109)
        call system("chmod +x TempIterateGeneProb.sh")
        call system("./TempIterateGeneProb.sh")
    
        print*, " "
        print*, " ","       Calculating genotype probabilities"
        JobsDone(:)=0
        JDone(:)=0
        if (RestartOption/=3) then
            do
                do i=1,nProcessors
                    write (filout,'("./IterateGeneProb/GeneProb"i0,"/GpDone.txt")')i
                    inquire(file=trim(filout),exist=FileExists)
                    if (FileExists .eqv. .true.) then
                        if (JDone(i)==0) print*, " ","      GeneProb job ",i," done"
                        JobsDone(i)=1
                        JDone(i)=1
                    endif
                enddo
                if (sum(JobsDone(:))==nProcessors) exit
            enddo
        endif   
    else
        write (filout,'("cd IterateGeneProb/")')
        write(f,'(i0)') nProcessors
        write (109,*) trim(filout)
        write(109,*) "cp ../../../SharedFiles/AlphaImpute_scripts/runSubmitGeneProb.sh ."
        write(109,*) "cp ../../../SharedFiles/AlphaImpute_scripts/submitGeneProb.sh ."
        write(109,'(a,a)') "./runSubmitGeneProb.sh ",trim(f)
        print*, " "
        close (109)
        call system("chmod +x TempIterateGeneProb.sh")
        call system("./TempIterateGeneProb.sh")
        call system("rm TempIterateGeneProb.sh")
    
        JobsDone(:)=0
        if (RestartOption/=3) then  
            do
                do i=1,nProcessors
                    write (filout,'("./IterateGeneProb/GeneProb"i0,"/GpDone.txt")')i
                    inquire(file=trim(filout),exist=FileExists)
                    if ((FileExists .eqv. .true.).and.(JobsDone(i)==0)) then
                        print*, " ","IterateGeneProb job ",i," done"
                        JobsDone(i)=1
                    endif
                enddo
                call sleep(SleepParameter)
                if (sum(JobsDone(:))==nProcessors) exit
            enddo
        endif   
    endif
    
    close (109)
    
    if (RestartOption==3) then
        open (unit=109,file="Tmp2345678.txt",status="unknown")
        do i=1,nAnisP
            write (109,'(i10,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ImputePhase(i,:,1)
            write (109,'(i10,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ImputePhase(i,:,2)          
            write (109,'(i10,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ImputeGenos(i,:)
        enddo
        close (109) 
        print*, "Restart option 3 stops program after Iterate Geneprob jobs have been submitted"
        stop
    endif
endif

call IterateGeneProbPhase
call IterateParentHomoFill
call PhaseComplement
call IterateMakeGenotype

do i=1,nAnisP
    do j=1,nSnpIterate
        do k=1,2
            if (ImputePhase(i,j,k)/=9) ProbImputePhase(i,j,k)=float(ImputePhase(i,j,k))
        enddo

        if (ImputeGenos(i,j)/=9) then
            ProbImputeGenos(i,j)=float(ImputeGenos(i,j))
        else
            ProbImputeGenos(i,j)=sum(ProbImputePhase(i,j,:))
        endif
    enddo
enddo

ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

end subroutine IterateGeneProbs

!#############################################################################################################################################################################################################################

subroutine IterateInsteadOfGeneProbs
use Global
implicit none

integer :: e,i,j,k,Counter,ParId
real,allocatable,dimension(:) :: TempAlleleFreq

if (SexOpt==1) then
    if (OutOpt==0) nSnpIterate=nSnp
    if (OutOpt==1) nSnpIterate=nSnpRaw
    
    allocate(ProbImputeGenos(0:nAnisP,nSnpIterate))
    allocate(ProbImputePhase(0:nAnisP,nSnpIterate,2))
    allocate(TempAlleleFreq(nSnpIterate))
    
    TempAlleleFreq=0.0
    do j=1,nSnpIterate
        Counter=0
        do i=1,nAnisP
            if (ImputeGenos(i,j)/=9) then
                TempAlleleFreq(j)=TempAlleleFreq(j)+ImputeGenos(i,j)
                Counter=Counter+2
            endif
        enddo
        if (Counter/=0) then
            TempAlleleFreq(j)=TempAlleleFreq(j)/Counter
        else
            TempAlleleFreq(j)=9.0
        endif
    enddo   
    
    ProbImputePhase(0,:,1)=TempAlleleFreq(:)
    ProbImputePhase(0,:,2)=TempAlleleFreq(:)
    ProbImputeGenos(0,:)=2*TempAlleleFreq(:)
    ProbImputeGenos(1:nAnisP,:)=-9.0
    ProbImputePhase(1:nAnisP,:,:)=-9.0
    
    call IterateParentHomoFill
    call PhaseComplement
    call IterateMakeGenotype
    
    do i=1,nAnisP
        do j=1,nSnpIterate
            do e=1,2
                if (ImputePhase(i,j,e)/=9) ProbImputePhase(i,j,e)=float(ImputePhase(i,j,e))
            enddo
        enddo   
    enddo   
    
    do i=1,nAnisP
        do e=1,2
            ParId=RecPed(i,e+1)
            if (ParId==0) then
                do j=1,nSnpIterate
                    if (ImputePhase(i,j,e)==9) ProbImputePhase(i,j,e)=TempAlleleFreq(j)
                enddo
            endif
            if (RecGender(i)==HomGameticStatus) then
                do j=1,nSnpIterate
                    if (ImputePhase(i,j,e)==9) then
                        ProbImputePhase(i,j,e)=(sum(ProbImputePhase(ParId,j,:))/2)
                    endif
                enddo
            endif
        enddo
        if (RecGender(i)==HetGameticStatus) then 
            ParId=RecPed(i,HomGameticStatus+1)
            do j=1,nSnpIterate
                if (ImputePhase(i,j,1)==9) then
                    ProbImputePhase(i,j,:)=(sum(ProbImputePhase(ParId,j,:))/2)
                endif
            enddo
        endif
    enddo
    
    do i=1,nAnisP
        do j=1,nSnpIterate
            do k=1,2
                if (ImputePhase(i,j,k)/=9) ProbImputePhase(i,j,k)=float(ImputePhase(i,j,k))
            enddo
            if (ImputeGenos(i,j)/=9) then
                ProbImputeGenos(i,j)=float(ImputeGenos(i,j))
            else
                ProbImputeGenos(i,j)=sum(ProbImputePhase(i,j,:))
            endif
        enddo
    enddo
    
    if (SexOpt==1) then
        do i=1,nAnisP
            if (RecGender(i)==HetGameticStatus) then
                do j=1,nSnpIterate
                    if ((ImputePhase(i,j,1)==9).and.(ImputePhase(i,j,2)/=9)) ImputePhase(i,j,1)=ImputePhase(i,j,2)
                    if ((ImputePhase(i,j,2)==9).and.(ImputePhase(i,j,1)/=9)) ImputePhase(i,j,2)=ImputePhase(i,j,1)      
                    if ((ProbImputePhase(i,j,1)==-9.0).and.(ProbImputePhase(i,j,2)/=-9.0)) ProbImputePhase(i,j,1)=ProbImputePhase(i,j,2)
                    if ((ProbImputePhase(i,j,2)==-9.0).and.(ProbImputePhase(i,j,1)/=-9.0)) ProbImputePhase(i,j,2)=ProbImputePhase(i,j,1)        
                enddo
            endif
        enddo
    endif
    
    do i=1,nAnisP
        do j=1,nSnpIterate
            if (ProbImputeGenos(i,j)==-9.0) ProbImputeGenos(i,j)=sum(ProbImputePhase(i,j,:))
            if (ProbImputeGenos(i,j)>1.999) ImputeGenos(i,j)=2
            if (ProbImputeGenos(i,j)<0.0001) ImputeGenos(i,j)=0
            if ((ProbImputeGenos(i,j)>0.999).and.(ProbImputeGenos(i,j)<1.00001)) ImputeGenos(i,j)=1
        enddo   
    enddo   
    
    
    do j=1,nSnpIterate
        if (TempAlleleFreq(j)==9.0) then
            ProbImputeGenos(:,j)=9.0
            ProbImputePhase(:,j,:)=9.0
        endif
    enddo
    
    ImputePhase(0,:,:)=9
    ImputeGenos(0,:)=9

else

    if (OutOpt==0) nSnpIterate=nSnp
    if (OutOpt==1) nSnpIterate=nSnpRaw
    
    allocate(ProbImputeGenos(0:nAnisP,nSnpIterate))
    allocate(ProbImputePhase(0:nAnisP,nSnpIterate,2))
    allocate(TempAlleleFreq(nSnpIterate))
    
    TempAlleleFreq=0.0
    do j=1,nSnpIterate
        Counter=0
        do i=1,nAnisP
            if (ImputeGenos(i,j)/=9) then
                TempAlleleFreq(j)=TempAlleleFreq(j)+ImputeGenos(i,j)
                Counter=Counter+2
            endif
        enddo
        if (Counter/=0) then
            TempAlleleFreq(j)=TempAlleleFreq(j)/Counter
        else
            TempAlleleFreq(j)=9.0
        endif
    enddo   
    
    ProbImputePhase(0,:,1)=TempAlleleFreq(:)
    ProbImputePhase(0,:,2)=TempAlleleFreq(:)
    ProbImputeGenos(0,:)=2*TempAlleleFreq(:)
    ProbImputeGenos(1:nAnisP,:)=-9.0
    ProbImputePhase(1:nAnisP,:,:)=-9.0
    
    call IterateParentHomoFill
    call PhaseComplement
    call IterateMakeGenotype
    
    do i=1,nAnisP
        do j=1,nSnpIterate
            do e=1,2
                if (ImputePhase(i,j,e)/=9) ProbImputePhase(i,j,e)=float(ImputePhase(i,j,e))
            enddo
        enddo   
    enddo   
    
    do i=1,nAnisP
        do e=1,2
            ParId=RecPed(i,e+1)
            if (ParId==0) then
                do j=1,nSnpIterate
                    if (ImputePhase(i,j,e)==9) ProbImputePhase(i,j,e)=TempAlleleFreq(j)
                enddo
            endif
            do j=1,nSnpIterate
                if (ImputePhase(i,j,e)==9) then
                    ProbImputePhase(i,j,e)=(sum(ProbImputePhase(ParId,j,:))/2)
                endif
            enddo
        enddo
    enddo
    
    do i=1,nAnisP
        do j=1,nSnpIterate
            do k=1,2
                if (ImputePhase(i,j,k)/=9) ProbImputePhase(i,j,k)=float(ImputePhase(i,j,k))
            enddo
            if (ImputeGenos(i,j)/=9) then
                ProbImputeGenos(i,j)=float(ImputeGenos(i,j))
            else
                ProbImputeGenos(i,j)=sum(ProbImputePhase(i,j,:))
            endif
        enddo
    enddo
    
    do i=1,nAnisP
        do j=1,nSnpIterate
            if (ProbImputeGenos(i,j)==-9.0) ProbImputeGenos(i,j)=sum(ProbImputePhase(i,j,:))
            if (ProbImputeGenos(i,j)>1.999) ImputeGenos(i,j)=2
            if (ProbImputeGenos(i,j)<0.0001) ImputeGenos(i,j)=0
            if ((ProbImputeGenos(i,j)>0.999).and.(ProbImputeGenos(i,j)<1.00001)) ImputeGenos(i,j)=1
        enddo   
    enddo   
    
    do j=1,nSnpIterate
        if (TempAlleleFreq(j)==9.0) then
            ProbImputeGenos(:,j)=9.0
            ProbImputePhase(:,j,:)=9.0
        endif
    enddo
    
    ImputePhase(0,:,:)=9
    ImputeGenos(0,:)=9

endif   

end subroutine IterateInsteadOfGeneProbs

!#############################################################################################################################################################################################################################

subroutine WriteOutResults
use Global
use GlobalPedigree
use GlobalVariablesHmmMaCH
implicit none

character(len=7) :: cm !use for formatting output - allows for up to 1 million SNPs
integer :: i,j,k,l,WorkTmp(nSnpRaw)
double precision :: ImputationQuality(nAnisP,6)
character(len=300) :: TmpId 

write(cm,'(I7)') nSnpRaw !for formatting
cm = adjustl(cm)

open (unit=33,file="./Results/ImputePhase.txt",status="unknown")
open (unit=34,file="./Results/ImputeGenotypes.txt",status="unknown")
open (unit=40,file="./Results/ImputePhaseProbabilities.txt",status="unknown")
open (unit=41,file="./Results/ImputeGenotypeProbabilities.txt",status="unknown")
open (unit=50,file="./Results/ImputationQualityIndividual.txt",status="unknown")
open (unit=51,file="./Results/ImputationQualitySnp.txt",status="unknown")
open (unit=52,file="./Results/WellPhasedIndividuals.txt",status="unknown")

if (OutOpt==0) then

    if (SexOpt==0) then
    
        do i=1,nAnisP
            do j=1,nSnp
                if (ImputePhase(i,j,1)<0) ImputePhase(i,j,1)=9 
                if (ImputePhase(i,j,1)>1) ImputePhase(i,j,1)=9      
                if (ImputePhase(i,j,2)<0) ImputePhase(i,j,2)=9 
                if (ImputePhase(i,j,2)>1) ImputePhase(i,j,2)=9      
                if (ImputeGenos(i,j)>2) ImputeGenos(i,j)=9  
                if (ImputeGenos(i,j)<0) ImputeGenos(i,j)=9                      
            enddo
        enddo

        open (unit=39,file="IterateGeneProb/IterateGeneProbInput.txt")
        do i=1,nAnisP
            write (39,'(3i10,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') RecPed(i,:),ImputeGenos(i,:)
        enddo
        call flush(39)
        
        if (BypassGeneProb==0) then 
            call IterateGeneProbs
        else
            call IterateInsteadOfGeneProbs
        endif   
    else

        do i=1,nAnisP
            do j=1,nSnp
                if (ImputePhase(i,j,1)<0) ImputePhase(i,j,1)=9 
                if (ImputePhase(i,j,1)>1) ImputePhase(i,j,1)=9      
                if (ImputePhase(i,j,2)<0) ImputePhase(i,j,2)=9 
                if (ImputePhase(i,j,2)>1) ImputePhase(i,j,2)=9      
                if (ImputeGenos(i,j)>2) ImputeGenos(i,j)=9  
                if (ImputeGenos(i,j)<0) ImputeGenos(i,j)=9                      
            enddo
        enddo
    
    
        call IterateInsteadOfGeneProbs
    endif

    do i=GlobalExtraAnimals+1,nAnisP
         !write (33,'(a,'//cm//'(1x,i1))') Id(i),ImputePhase(i,:,1)
         !write (33,'(a,'//cm//'(1x,i1))') Id(i),ImputePhase(i,:,2)
         !write (34,'(a,'//cm//'(1x,i1))') Id(i),ImputeGenos(i,:)
         write (33,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') Id(i),ImputePhase(i,:,1)
         write (33,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') Id(i),ImputePhase(i,:,2)
         write (34,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') Id(i),ImputeGenos(i,:)
         
    enddo

    do i=GlobalExtraAnimals+1,nAnisP
         !write (40,'(a,'//cm//'(1x,f5.2))') Id(i),ProbImputePhase(i,:,1)
         !write (40,'(a,'//cm//'(1x,f5.2))') Id(i),ProbImputePhase(i,:,2)
         !write (41,'(a,'//cm//'(1x,f5.2))') Id(i),ProbImputeGenos(i,:)
         write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') Id(i),ProbImputePhase(i,:,1)
         write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') Id(i),ProbImputePhase(i,:,2)
         write (41,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') Id(i),ProbImputeGenos(i,:)

    enddo

    if (SexOpt==1) then
        allocate(Maf(nSnp))
        do j=1,nSnp
            Maf(j)=sum(ProbImputeGenos(:,j))/(2*nAnisP)
        enddo
        if (WindowsLinux==1) then
            open (unit=111,file=".\Miscellaneous\MinorAlleleFrequency.txt",status="unknown")
        else
            open (unit=111,file="./Miscellaneous/MinorAlleleFrequency.txt",status="unknown")
        endif

        do j=1,nSnpRaw
            write (111,*) j,Maf(j)
        enddo
        close(111)      
    endif
    ImputationQuality(:,1)=sum(2*Maf(:))/nSnp

    ImputationQuality(:,2)=0.0
    do i=GlobalExtraAnimals+1,nAnisP
        do j=1,nSnp
            !ImputationQuality(i,2)=ImputationQuality(i,2)+abs(ProbImputeGenos(i,j)-(2*Maf(j)))
            ImputationQuality(i,2)=ImputationQuality(i,2)+abs(ProbImputeGenos(i,j)-((Maf(j)**4)+(4*(Maf(j)**2)*((1.0-Maf(j))**2))+(1.0-Maf(j)**4)))
        enddo   
        ImputationQuality(i,2)=ImputationQuality(i,2)/nSnp
        ImputationQuality(i,3)=float(nSnp-count(ImputePhase(i,:,1)==9))/nSnp
        ImputationQuality(i,4)=float(nSnp-count(ImputePhase(i,:,2)==9))/nSnp    
        ImputationQuality(i,5)=(ImputationQuality(i,3)+ImputationQuality(i,4))/2                    
        ImputationQuality(i,6)=float(nSnp-count(ImputeGenos(i,:)==9))/nSnp  
        write (50,'(a20,20000f7.2)') Id(i),ImputationQuality(i,:)
    enddo

    do j=1,nSnp
        write (51,'(i10,20000f7.2)') j,float(((nAnisP-(GlobalExtraAnimals+1))+1)-count(ImputeGenos(GlobalExtraAnimals+1:nAnisP,j)==9))/((nAnisP-(GlobalExtraAnimals+1))+1)
    enddo

    WellPhasedThresh=WellPhasedThresh/100   
    do i=GlobalExtraAnimals+1,nAnisP
        if (ImputationQuality(i,5)>=WellPhasedThresh) then
            write (52,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') Id(i),ImputePhase(i,:,1)
            write (52,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') Id(i),ImputePhase(i,:,2)
        endif   
    enddo

else

    do i=1,nAnisP
        do j=1,nSnp
            if (ImputePhase(i,j,1)<0) ImputePhase(i,j,1)=9 
            if (ImputePhase(i,j,1)>1) ImputePhase(i,j,1)=9      
            if (ImputePhase(i,j,2)<0) ImputePhase(i,j,2)=9 
            if (ImputePhase(i,j,2)>1) ImputePhase(i,j,2)=9      
            if (ImputeGenos(i,j)>2) ImputeGenos(i,j)=9  
            if (ImputeGenos(i,j)<0) ImputeGenos(i,j)=9                      
        enddo
    enddo


    open (unit=42,file=trim(GenotypeFile),status='old')
    allocate(TmpGenos(0:nAnisP,nSnpRaw))
    allocate(TmpPhase(0:nAnisP,nSnpRaw,2))
    TmpGenos=9  
    TmpPhase=9

    l=0 
    do j=1,nSnpRaw
        if (SnpIncluded(j)==1) then
            l=l+1
            TmpGenos(:,j)=ImputeGenos(:,l)
            TmpPhase(:,j,1)=ImputePhase(:,l,1)
            TmpPhase(:,j,2)=ImputePhase(:,l,2)
        endif
    enddo

    do i=1,nAnisG
        read (42,*) TmpId,WorkTmp(:)
        do k=1,nAnisP
            if (trim(Id(k))==trim(TmpId)) then
                do j=1,nSnpRaw
                    if (SnpIncluded(j)==0) then
                        if (WorkTmp(j)==1) TmpGenos(k,j)=1
                        if (WorkTmp(j)==0) then
                            TmpPhase(k,j,:)=0
                            TmpGenos(k,j)=0
                        endif
                        if (WorkTmp(j)==2) then 
                            TmpPhase(k,j,:)=1
                            TmpGenos(k,j)=2
                        endif
                    endif
                enddo
                exit
            endif
        enddo
    enddo
    close(42)

    do i=GlobalExtraAnimals+1,nAnisP
         write (33,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') Id(i),TmpPhase(i,:,1)
         write (33,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') Id(i),TmpPhase(i,:,2)
         write (34,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') Id(i),TmpGenos(i,:)
    enddo
    if (SexOpt==0) then
        open (unit=39,file="IterateGeneProb/IterateGeneProbInput.txt")
        do i=1,nAnisP
            write (39,'(3i10,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') RecPed(i,:),TmpGenos(i,:)
        enddo
        call flush(39)
    endif

    deallocate(ImputePhase)
    deallocate(ImputeGenos)
    allocate(ImputeGenos(0:nAnisP,nSnpRaw))
    allocate(ImputePhase(0:nAnisP,nSnpRaw,2))
    ImputeGenos=TmpGenos
    ImputePhase=TmpPhase
!REMOVE THIS WHEN HMM IS FINALISED
    if (HMMOption/=3) then
        if (SexOpt==0) then
            if (BypassGeneProb==0) then 
                call IterateGeneProbs
            else
                call IterateInsteadOfGeneProbs
            endif   
        endif
        if (SexOpt==1) call IterateInsteadOfGeneProbs
    endif
!REMOVE THIS

    if (HMMOption==3) then
        l=0 
        do j=1,nSnpRaw
            if (SnpIncluded(j)==1) then
                l=l+1
                do i=1,nAnisG
                    ProbImputeGenos(GlobalHmmID(i),l)=ProbImputeGenosHmm(i,l)
                enddo
            endif
        enddo
    endif

    do i=GlobalExtraAnimals+1,nAnisP
         write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') Id(i),ProbImputePhase(i,:,1)
         write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') Id(i),ProbImputePhase(i,:,2)
         write (41,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') Id(i),ProbImputeGenos(i,:)
    enddo
    if ((SexOpt==1).or.(BypassGeneProb==1)) then
        allocate(Maf(nSnpRaw))
        do j=1,nSnpRaw
            Maf(j)=sum(ProbImputeGenos(:,j))/(2*nAnisP)
        enddo
        if (WindowsLinux==1) then
            open (unit=111,file=".\Miscellaneous\MinorAlleleFrequency.txt",status="unknown")
        else
            open (unit=111,file="./Miscellaneous/MinorAlleleFrequency.txt",status="unknown")
        endif

        do j=1,nSnpRaw
            write (111,*) j,Maf(j)
        enddo
        close(111)
    endif

    ImputationQuality(:,1)=sum(2*Maf(:))/nSnpRaw
    ImputationQuality(:,2)=0.0
    do i=GlobalExtraAnimals+1,nAnisP
        do j=1,nSnpRaw
            !ImputationQuality(i,2)=ImputationQuality(i,2)+abs(ProbImputeGenos(i,j)-(2*Maf(j)))
            ImputationQuality(i,2)=ImputationQuality(i,2)+abs(ProbImputeGenos(i,j)-((Maf(j)**4)+(4*(Maf(j)**2)*((1.0-Maf(j))**2))+(1.0-Maf(j)**4)))         
        enddo   
        ImputationQuality(i,2)=ImputationQuality(i,2)/nSnpRaw
        ImputationQuality(i,3)=float(nSnpRaw-count(ImputePhase(i,:,1)==9))/nSnpRaw
        ImputationQuality(i,4)=float(nSnpRaw-count(ImputePhase(i,:,2)==9))/nSnpRaw
        ImputationQuality(i,5)=(ImputationQuality(i,3)+ImputationQuality(i,4))/2
        ImputationQuality(i,6)=float(nSnpRaw-count(ImputeGenos(i,:)==9))/nSnpRaw                    
        write (50,'(a20,20000f7.2)') Id(i),ImputationQuality(i,:)
    enddo

    do j=1,nSnpRaw
        write (51,'(i10,20000f7.2)') j,float(((nAnisP-(GlobalExtraAnimals+1))+1)-count(ImputeGenos(GlobalExtraAnimals+1:nAnisP,j)==9))/((nAnisP-(GlobalExtraAnimals+1))+1)
    enddo

    WellPhasedThresh=WellPhasedThresh/100   
    do i=GlobalExtraAnimals+1,nAnisP
        if (ImputationQuality(i,5)>=WellPhasedThresh) then
            write (52,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') Id(i),TmpPhase(i,:,1)
            write (52,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') Id(i),TmpPhase(i,:,2)
        endif   
    enddo
endif

close (33)
close (34)
close (40)
close (41)
close (50)
close (51)
close (52)

end subroutine WriteOutResults

!#############################################################################################################################################################################################################################

subroutine ModelRecomb
use Global
use GlobalPedigree
implicit none

integer :: e,i,j,h,k,l,SuperJ,StartDisFound,EndDisFound,HetEnd,HetStart,RSide,LSide,PatMat,SireDamRL,nSnpFinal,StSnp,EnSnp,Counter
integer :: StartDisPrev,EndDisPrev,RecombOnOff
integer :: GamA,GamB,Tmp,StartDisOld,StartDisTmp
integer :: CountRightSwitch,CountLeftSwitch,StartPt,EndPt,PedId,StartDis,EndDis,dum,StartJ,nRec
integer(kind=1),allocatable,dimension(:,:,:) :: WorkPhase,TempWork
integer,allocatable,dimension(:) :: WorkLeft,WorkRight,TempVec,StR,EnR,StRNarrow,EnRNarrow
real,allocatable,dimension(:) :: LengthVec
real,allocatable,dimension(:,:) :: PatAlleleProb,MatAlleleProb,GeneProbWork
real :: Val
character(len=300) :: filout
character(len=7) :: cm 
write(cm,'(I7)') nSnpRaw !for formatting
cm = adjustl(cm)

open (unit=42,file="Results/RecombinationInformation.txt",status="unknown")
open (unit=43,file="Results/RecombinationInformationNarrow.txt",status="unknown")
open (unit=44,file="Results/NumberRecombinations.txt",status="unknown")
open (unit=45,file="Results/RecombinationInformationR.txt",status="unknown")
open (unit=46,file="Results/RecombinationInformationNarrowR.txt",status="unknown")


if (OutOpt==0) then
    nSnpFinal=nSnp
else
    nSnpFinal=nSnpRaw
endif

Tmp=int(float(nSnp)/nProcessors)
GpIndex(1,1)=1
GpIndex(1,2)=Tmp
if (nProcessors>1) then
    do i=2,nProcessors
        GpIndex(i,1)=GpIndex(i-1,1)+Tmp     
        GpIndex(i,2)=GpIndex(i-1,2)+Tmp
    enddo
endif
GpIndex(nProcessors,2)=nSnp

allocate(GlobalWorkPhase(0:nAnisP,nSnp,2))
allocate(WorkPhase(0:nAnisP,nSnpFinal,2))
allocate(TempVec(nSnpFinal))
allocate(LengthVec(nSnpFinal))
allocate(WorkLeft(nSnpFinal))
allocate(WorkRight(nSnpFinal))
allocate(TempWork(0:nAnisP,nSnpFinal,2))
allocate(PatAlleleProb(nSnpFinal,2))
allocate(MatAlleleProb(nSnpFinal,2))
allocate(GeneProbWork(nSnpFinal,4))
allocate(StR(nSnpFinal))
allocate(EnR(nSnpFinal))
allocate(StRNarrow(nSnpFinal))
allocate(EnRNarrow(nSnpFinal))

WorkPhase=9
if (SexOpt==0) then
    if (BypassGeneProb==0) then
        call ReReadGeneProbs
    else
        call InsteadOfReReadGeneProb
    endif
endif           
if (SexOpt==1) call InsteadOfReReadGeneProb

l=0 
do j=1,nSnpFinal

    if (SnpIncluded(j)==1) then
        l=l+1
        WorkPhase(:,j,1)=GlobalWorkPhase(:,l,1)
        WorkPhase(:,j,2)=GlobalWorkPhase(:,l,2)
    endif
enddo

do i=1,nAnisP
    HetEnd=-1
    HetStart=-1
    WorkRight(:)=9
    WorkLeft(:)=9
    do e=1,2
        nRec=0  
        StR=-9
        EnR=-9
        StRNarrow=-9
        EnRNarrow=-9

        PatMat=e
        SireDamRL=e+1
        CountLeftSwitch=0
        CountRightSwitch=0
        PedId=RecPed(i,SireDamRL)
        if ((SexOpt==1).and.(RecGender(PedId)==HetGameticStatus)) cycle
        if ((PedId>0).and.((float(count(ImputePhase(PedId,:,:)==9))/(2*nSnpFinal))<0.30)) then          !(RecIdHDIndex(PedId)==1)   
            WorkRight=9
            RSide=9
            do j=1,nSnpFinal
                if ((ImputePhase(PedId,j,1)/=ImputePhase(PedId,j,2)).and.(ImputePhase(PedId,j,1)/=9).and.(ImputePhase(PedId,j,2)/=9))  then
                    HetStart=j
                    if (ImputePhase(i,HetStart,PatMat)==ImputePhase(PedId,HetStart,1)) then
                        WorkRight(HetStart)=1 
                        RSide=1
                        exit
                    endif
                    if (ImputePhase(i,HetStart,PatMat)==ImputePhase(PedId,HetStart,2)) then
                        WorkRight(HetStart)=2
                        RSide=2 
                        exit
                    endif
                endif  
            enddo
            if (RSide/=9) then
                do j=HetStart+1,nSnpFinal
                    if ((ImputePhase(i,j,PatMat)/=ImputePhase(PedId,j,RSide)).and.(ImputePhase(PedId,j,RSide)/=9).and.(ImputePhase(i,j,PatMat)/=9)) then
                        RSide=abs((RSide-1)-1)+1
                        CountRightSwitch=CountRightSwitch+1     
                    endif
                    WorkRight(j)=RSide
                enddo
            endif

            LSide=9
            do j=nSnpFinal,1,-1
                if ((ImputePhase(PedId,j,1)/=ImputePhase(PedId,j,2)).and.(ImputePhase(PedId,j,1)/=9).and.(ImputePhase(PedId,j,2)/=9))  then
                    HetEnd=j
                    if (ImputePhase(i,HetEnd,PatMat)==ImputePhase(PedId,HetEnd,1)) then
                        WorkLeft(HetEnd)=1  !£$$$$
                        LSide=1
                        exit
                    endif
                    if (ImputePhase(i,HetEnd,PatMat)==ImputePhase(PedId,HetEnd,2)) then
                        WorkLeft(HetEnd)=2  !£$$$$
                        LSide=2 
                        exit
                    endif
                endif  
            enddo
            if (LSide/=9) then
                do j=HetEnd-1,1,-1
                    if ((ImputePhase(i,j,PatMat)/=ImputePhase(PedId,j,LSide)).and.(ImputePhase(PedId,j,LSide)/=9).and.(ImputePhase(i,j,PatMat)/=9)) then
                        LSide=abs((LSide-1)-1)+1
                        CountLeftSwitch=CountLeftSwitch+1
                    endif
                    WorkLeft(j)=LSide
                enddo
            endif
            StartDis=-9
            EndDis=-9
            TempVec=9
            LengthVec=0.0
            StartJ=1

            !prototype start

            StartDis=1
            StartDisOld=1
            EndDis=nSnpFinal
            nRec=0
            SuperJ=0
            TempVec=9
            LengthVec=0.0

            StartDisPrev=StartDis                   
            EndDisPrev=-9                                           

            do while (SuperJ<nSnp) 
                SuperJ=SuperJ+1             
                
                !Finding StartDis and Moving it left and EndDis and Movie it Right
                
                !Find StartDis
                StartDisFound=0
                if (abs(WorkLeft(SuperJ)-WorkRight(SuperJ))==1) then
                    StartDisOld=StartDis
                    StartDis=SuperJ
                    StartDisFound=1
                endif

                if (StartDisFound==1) then
                    nRec=nRec+1
                                        
                    !Find EndDis
                    EndDisFound=0
                    do k=StartDis+1,nSnpFinal
                        if (WorkLeft(k)==WorkRight(k)) then
                            EndDis=k
                            SuperJ=k
                            EndDisFound=1
                            exit
                        endif
                    enddo   
                    if (EndDisFound==0) then
                        EndDis=nSnpFinal
                        SuperJ=nSnpFinal
                    endif   
                
                    !Move StartDis to the left informative marker
                    StartDisTmp=StartDis            
                    do k=StartDis,1,-1
                        if ((WorkPhase(PedId,k,1)+WorkPhase(PedId,k,2))==1) then
                            if (WorkPhase(i,k,e)/=9) then
                                StartDis=k
                                exit
                            endif
                        endif
                    enddo
                    if (StartDis<=StartDisOld) StartDis=StartDisTmp
                    if (StartDis<1) StartDis=1
                    if (StartDis<EndDisPrev) StartDis=EndDisPrev
                    
                    !Move EndDis to the right informative marker            
                    do k=EndDis,nSnpFinal
                        if ((WorkPhase(PedId,k,1)+WorkPhase(PedId,k,2))==1) then
                            if (WorkPhase(i,k,e)/=9) then
                                EndDis=k
                                SuperJ=k            
                                exit
                            endif
                        endif                           
                    enddo
                    
                    StR(nRec)=StartDis                  
                    EnR(nRec)=EndDis
                    LengthVec(StartDis:EndDis)=1.0/(EndDis-StartDis)                
                    StartDisPrev=StartDis                   
                    EndDisPrev=EndDis   
                    
                endif
            enddo
            !prototype end (temp)
            StR(nRec+1)=-9
            EnR(nRec+1)=-9

            GamA=9
            GamB=9
            RecombOnOff=0
            k=1
            do j=1,nSnpFinal

                if (j==StR(k)) then
                    k=k+1
                    RecombOnOff=1
                    Counter=0
                    GamA=WorkRight(j)
                
                    If (GamA==1) GamB=2
                    If (GamA==2) GamB=1                 
                endif
                if (j==EnR(k-1)) then
                    RecombOnOff=0       
                endif
                if (GamA==9) cycle
                if (RecombOnOff==1) then
                    Counter=Counter+1
                    if (Counter/=1) then
                        if (ImputePhase(PedId,j,1)/=ImputePhase(PedId,j,2)) then
                            if ((ImputePhase(PedId,j,1)/=9).and.(ImputePhase(PedId,j,2)/=9)) then
                                ImputePhase(i,j,e)=9
                                ImputeGenos(i,j)=9
                                ProbImputePhase(i,j,e)=((1.0-(LengthVec(j)*Counter))*ImputePhase(PedId,j,GamA))+(LengthVec(j)*Counter*ImputePhase(PedId,j,GamB))
                                ProbImputeGenos(i,j)=ProbImputePhase(i,j,1)+ProbImputePhase(i,j,2)      

                            endif
                        endif       
                    endif
                endif   
            enddo           
        endif
        if (i>GlobalExtraAnimals) then
            write (44,'(a20,20000i20)') Id(i),nRec
            write (42,'(a20,20000i20)') Id(i),nRec,StR(1:nRec)
            write (42,'(a20,20000i20)') Id(i),nRec,EnR(1:nRec)      
!
!           write (43,'(a20,20000i20)') Id(i),nRec,StRNarrow(1:nRec)
!           write (43,'(a20,20000i20)') Id(i),nRec,EnRNarrow(1:nRec)
            do j=1,nRec
                write (45,'(a20,20000i20)') Id(i),e,StR(j),EnR(j)   
!               write (46,'(a20,20000i20)') Id(i),e,StRNarrow(j),EnRNarrow(j)   
            enddo       

        endif   
    enddo
    WorkPhase(i,:,:)=ImputePhase(i,:,:)

enddo

open (unit=33,file="Results/ImputePhase.txt",status="unknown")
open (unit=34,file="Results/ImputeGenotypes.txt",status="unknown")
open (unit=40,file="Results/ImputePhaseProbabilities.txt",status="unknown")
open (unit=41,file="Results/ImputeGenotypeProbabilities.txt",status="unknown")

do i=GlobalExtraAnimals+1,nAnisP
     write (33,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') Id(i),ImputePhase(i,:,1)
     write (33,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') Id(i),ImputePhase(i,:,2)
     write (34,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') Id(i),ImputeGenos(i,:)
enddo

do i=GlobalExtraAnimals+1,nAnisP
     write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') Id(i),ProbImputePhase(i,:,1)
     write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') Id(i),ProbImputePhase(i,:,2)
     write (41,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') Id(i),ProbImputeGenos(i,:)
enddo

print*, " "
print*, " ","Imputation by detection of recombination events completed"

end subroutine ModelRecomb

!#############################################################################################################################################################################################################################

subroutine IterateGeneProbPhase
use Global

implicit none

integer :: h,i,j,dum,StSnp,EnSnp,counter
real :: PatAlleleProb(nSnpIterate,2),MatAlleleProb(nSnpIterate,2),HetProb(nSnpIterate),GeneProbWork(nSnpIterate,4)
character(len=300) :: filout

allocate(Maf(nSnpIterate))

if (RestartOption==4) then
    open (unit=109,file="Tmp2345678.txt",status="unknown")
    do i=1,nAnisP
        read (109,*) ImputePhase(i,:,1)
        read (109,*) ImputePhase(i,:,2)         
        read (109,*) ImputeGenos(i,:)
    enddo   
    call system("rm Tmp2345678.txt")
endif

counter=0
do h=1,nProcessors
    write (filout,'("./IterateGeneProb/GeneProb"i0,"/GeneProbs.txt")')h         !here
    open (unit=110,file=trim(filout),status="unknown")
    write (filout,'("./IterateGeneProb/GeneProb"i0,"/MinorAlleleFrequency.txt")')h          !here
    open (unit=111,file=trim(filout),status="unknown")

    StSnp=GpIndex(h,1)
    EnSnp=GpIndex(h,2)
    do i=1,nAnisP
        do j=1,4
            read (110,*) dum,GeneProbWork(StSnp:EnSnp,j)
        enddo
        PatAlleleProb(StSnp:EnSnp,1)=GeneProbWork(StSnp:EnSnp,1)&
                +GeneProbWork(StSnp:EnSnp,2)
        PatAlleleProb(StSnp:EnSnp,2)=GeneProbWork(StSnp:EnSnp,3)&
                +GeneProbWork(StSnp:EnSnp,4)
        MatAlleleProb(StSnp:EnSnp,1)=GeneProbWork(StSnp:EnSnp,1)&
                +GeneProbWork(StSnp:EnSnp,3)
        MatAlleleProb(StSnp:EnSnp,2)=GeneProbWork(StSnp:EnSnp,2)&
                +GeneProbWork(StSnp:EnSnp,4)

        HetProb(StSnp:EnSnp)=GeneProbWork(StSnp:EnSnp,2)+GeneProbWork(StSnp:EnSnp,3) 

        do j=StSnp,EnSnp
            if ((PatAlleleProb(j,1)>=GeneProbThresh).and.(ImputePhase(i,j,1)==9)) ImputePhase(i,j,1)=0
            if ((PatAlleleProb(j,2)>=GeneProbThresh).and.(ImputePhase(i,j,1)==9)) ImputePhase(i,j,1)=1
            if ((MatAlleleProb(j,1)>=GeneProbThresh).and.(ImputePhase(i,j,2)==9)) ImputePhase(i,j,2)=0
            if ((MatAlleleProb(j,2)>=GeneProbThresh).and.(ImputePhase(i,j,2)==9)) ImputePhase(i,j,2)=1
            if (ImputePhase(i,j,1)==9) then
                ProbImputePhase(i,j,1)=GeneProbWork(j,3)+GeneProbWork(j,4)
            else
                ProbImputePhase(i,j,1)=float(ImputePhase(i,j,1))
            endif
            if (ImputePhase(i,j,2)==9) then
                ProbImputePhase(i,j,2)=GeneProbWork(j,2)+GeneProbWork(j,4)
            else
                ProbImputePhase(i,j,2)=float(ImputePhase(i,j,2))
            endif
            if (HetProb(j)>=GeneProbThresh) ProbImputeGenos(i,j)=1
        enddo
        ProbImputeGenos(i,StSnp:EnSnp)=ProbImputePhase(i,StSnp:EnSnp,1)+ProbImputePhase(i,StSnp:EnSnp,2)
    enddo

    do j=StSnp,EnSnp
        read (111,*) Maf(j)
    enddo

    close(110)
    close(111)
enddo

if (WindowsLinux==1) then
        open (unit=111,file=".\Miscellaneous\MinorAlleleFrequency.txt",status="unknown")
else
        open (unit=111,file="./Miscellaneous/MinorAlleleFrequency.txt",status="unknown")
endif

do j=1,nSnpIterate
    write (111,*) j,Maf(j)
enddo
close(111)

call IterateMakeGenotype

ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

end subroutine IterateGeneProbPhase
!#############################################################################################################################################################################################################################

subroutine IterateMakeGenotype
use Global
implicit none

integer :: i,j

do i=1,nAnisP
    do j=1,nSnpIterate
        if (ImputeGenos(i,j)==9) then
            if ((ImputePhase(i,j,1)/=9).and.(ImputePhase(i,j,2)/=9)) ImputeGenos(i,j)=sum(ImputePhase(i,j,:))
        endif
    enddo
enddo

ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

end subroutine IterateMakeGenotype

!#############################################################################################################################################################################################################################

subroutine IteratePhaseComplement
use Global
implicit none

integer :: i,j

do i=1,nAnisP
    do j=1,nSnpIterate
        if (ImputeGenos(i,j)/=9) then
            if ((ImputePhase(i,j,1)/=9).and.(ImputePhase(i,j,2)==9)) ImputePhase(i,j,2)=ImputeGenos(i,j)-ImputePhase(i,j,1)
            if ((ImputePhase(i,j,2)/=9).and.(ImputePhase(i,j,1)==9)) ImputePhase(i,j,1)=ImputeGenos(i,j)-ImputePhase(i,j,2)
        endif
    enddo
enddo

ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

end subroutine IteratePhaseComplement

!#############################################################################################################################################################################################################################

subroutine IterateParentHomoFill
use Global
implicit none

integer :: e,i,j,PedLoc,PatMat

if (SexOpt==0) then
    do i=1,nAnisP
        do e=1,2
            PedLoc=e+1
            do j=1,nSnpIterate
                if (ImputePhase(i,j,e)==9) then
                    if ((ImputePhase(RecPed(i,PedLoc),j,1)==ImputePhase(RecPed(i,PedLoc),j,2)).and.&
                        (ImputePhase(RecPed(i,PedLoc),j,1)/=9)) ImputePhase(i,j,e)=ImputePhase(RecPed(i,PedLoc),j,1)
                endif
            enddo
        enddo
    enddo
else
    do i=1,nAnisP
        if (RecGender(i)==HomGameticStatus) then  
            do e=1,2
                PedLoc=RecPed(i,e+1)
                do j=1,nSnpIterate
                    if (ImputePhase(i,j,e)==9) then
                        if ((ImputePhase(PedLoc,j,1)==ImputePhase(PedLoc,j,2)).and.(ImputePhase(PedLoc,j,1)/=9)) ImputePhase(i,j,e)=ImputePhase(PedLoc,j,1)
                    endif
                enddo
            enddo
        else
            PedLoc=RecPed(i,HomGameticStatus+1)
            do j=1,nSnpIterate
                if (ImputePhase(i,j,1)==9) then     !Comment From John Hickey I changed what was indexed e to 1 I think it is ok same thing in analogous routine
                                                    !There is no do loop for e here
                    if ((ImputePhase(PedLoc,j,1)==ImputePhase(PedLoc,j,2)).and.(ImputePhase(PedLoc,j,1)/=9)) ImputePhase(i,j,:)=ImputePhase(PedLoc,j,1)
                endif
            enddo
        endif
    enddo
endif   
ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

end subroutine IterateParentHomoFill

!######################################################################################################################################################################################

subroutine ReReadGeneProbs
! Read genotype probabilities from files and phase allele based in these probabilities. 
! This files should have been already created during previous calls to AlphaImpute (RestartOption<3)
! Phasing information is store in the variable GlobalWorkPhase
use Global

implicit none

integer :: h,i,j,dum,StSnp,EnSnp
real :: PatAlleleProb(nSnp,2),MatAlleleProb(nSnp,2),HetProb(nSnp),GeneProbWork(nSnp,4)
character(len=300) :: filout

GlobalWorkPhase=9
do h=1,nProcessors
    write (filout,'("GeneProb/GeneProb"i0,"/GeneProbs.txt")')h          !here
    open (unit=110,file=trim(filout),status="unknown")
    StSnp=GpIndex(h,1)          ! Where SNPs start
    EnSnp=GpIndex(h,2)          ! Where SNPs end
    do i=1,nAnisP                                           ! The number of lines of GeneProbs.txt files is = nAnisP x 4
        do j=1,4                                            ! where 4 stands for the two paternal and the two maternal haplotypes
            read (110,*) dum,GeneProbWork(StSnp:EnSnp,j)
        enddo

        ! GeneProbWork(:,1) == Probability 0-0 = Prob00
        ! GeneProbWork(:,2) == Probability 0-1 = Prob01
        ! GeneProbWork(:,3) == Probability 1-0 = Prob10 
        ! GeneProbWork(:,4) == Probability 1-1 = Prob11 
        PatAlleleProb(StSnp:EnSnp,1)=GeneProbWork(StSnp:EnSnp,1)+GeneProbWork(StSnp:EnSnp,2)    ! PatAlleleProb(:,1) == Probability Paternal allele is 0 = Prob00 + Prob01
        PatAlleleProb(StSnp:EnSnp,2)=GeneProbWork(StSnp:EnSnp,3)+GeneProbWork(StSnp:EnSnp,4)    ! PatAlleleProb(:,2) == Probability Paternal allele is 1 = Prob10 + Prob11
        MatAlleleProb(StSnp:EnSnp,1)=GeneProbWork(StSnp:EnSnp,1)+GeneProbWork(StSnp:EnSnp,3)    ! PatAlleleProb(:,3) == Probability Maternal allele is 0 = Prob00 + Prob10
        MatAlleleProb(StSnp:EnSnp,2)=GeneProbWork(StSnp:EnSnp,2)+GeneProbWork(StSnp:EnSnp,4)    ! PatAlleleProb(:,4) == Probability Maternal allele is 1 = Prob01 + Prob11

        do j=StSnp,EnSnp
            if (PatAlleleProb(j,1)>=GeneProbThresh) GlobalWorkPhase(i,j,1)=0
            if (PatAlleleProb(j,2)>=GeneProbThresh) GlobalWorkPhase(i,j,1)=1
            if (MatAlleleProb(j,1)>=GeneProbThresh) GlobalWorkPhase(i,j,2)=0
            if (MatAlleleProb(j,2)>=GeneProbThresh) GlobalWorkPhase(i,j,2)=1
        enddo
    enddo
    close(110)
enddo
GlobalWorkPhase(0,:,:)=9

end subroutine ReReadGeneProbs

!######################################################################################################################################################################################

subroutine InsteadOfReReadGeneProb
! Phase alleles in the SEX CHROMOSOME whenever it is possible (homozygous case). 
! Phasing information is store in the variable GlobalWorkPhase
use Global
implicit none

integer :: e,i,j,ParId

if (SexOpt==1) then                                         ! Sex chromosome
    deallocate(GlobalWorkPhase)
    allocate(GlobalWorkPhase(0:nAnisP,nSnp,2))
    Genos(0,:)=9                                                        ! Erase any possible information about the phantom parents
    GlobalWorkPhase=9
    do i=1,nAnisP
        do j=1,nSnp                                                     ! Phase alleles in the homozygous case
            if (Genos(i,j)==0) GlobalWorkPhase(i,j,:)=0
            if (Genos(i,j)==2) GlobalWorkPhase(i,j,:)=1
        enddo
        if (RecGender(i)/=HetGameticStatus) then
            do e=1,2                                                    ! Phase alleles for homogametic individuals in the homozygous case
                ParId=RecPed(i,e+1)
                do j=1,nSnp
                    if (Genos(ParId,j)==0) GlobalWorkPhase(i,j,e)=0
                    if (Genos(ParId,j)==2) GlobalWorkPhase(i,j,e)=1
                enddo
            enddo
        else
            ParId=RecPed(i,HomGameticStatus+1)                          ! Phase alleles for heterogametic individuals in the homozygous case 
            do j=1,nSnp
                if (Genos(ParId,j)==0) GlobalWorkPhase(i,j,:)=0
                if (Genos(ParId,j)==2) GlobalWorkPhase(i,j,:)=1
            enddo   
        endif
    enddo
    GlobalWorkPhase(0,:,:)=9
else                                ! Nothing is done in other chromosomes
    !! WARNING: This should be some copied, pasted and erased stuff
endif   

end subroutine InsteadOfReReadGeneProb

!######################################################################################################################################################################################

subroutine InsteadOfGeneProb
! Phase haplotypes whenever there is enough information from the parents (homozygous case)
use Global
implicit none

integer :: e,i,j,ParId

if (SexOpt==1) then                                                     ! Sex chromosome
    allocate(GlobalWorkPhase(0:nAnisP,nSnp,2))
    allocate(ImputeGenos(0:nAnisP,nSnp))
    allocate(ImputePhase(0:nAnisP,nSnp,2))

    Genos(0,:)=9
    GlobalWorkPhase=9
    do i=1,nAnisP
        do j=1,nSnp                                                     ! Phase in the homozygous case
            if (Genos(i,j)==0) GlobalWorkPhase(i,j,:)=0
            if (Genos(i,j)==2) GlobalWorkPhase(i,j,:)=1
        enddo
        if (RecGender(i)/=HetGameticStatus) then                        ! Am I homogametic?
            do e=1,2                                                    ! Phase a single haplotype whenever my parents are homozygous
                ParId=RecPed(i,e+1)
                do j=1,nSnp
                    if (Genos(ParId,j)==0) GlobalWorkPhase(i,j,e)=0
                    if (Genos(ParId,j)==2) GlobalWorkPhase(i,j,e)=1
                enddo
            enddo
        else                                                            ! Am I heterogametic?
            ParId=RecPed(i,HomGameticStatus+1)  
            do j=1,nSnp                                                 ! Phase the two haplotypes whenever my homogametic parent is homozygous
                if (Genos(ParId,j)==0) GlobalWorkPhase(i,j,:)=0
                if (Genos(ParId,j)==2) GlobalWorkPhase(i,j,:)=1
            enddo   
        endif
    enddo
    
    GlobalWorkPhase(0,:,:)=9
    ImputeGenos=9
    ImputePhase=9
    ImputeGenos(:,:)=Genos(:,:)
    ImputePhase(:,:,:)=GlobalWorkPhase(:,:,:)
    
    allocate(GlobalTmpCountInf(nAnisP,6))
    GlobalTmpCountInf(:,:)=0

else                                                                    ! Other chromosome
    allocate(GlobalWorkPhase(0:nAnisP,nSnp,2))
    allocate(ImputeGenos(0:nAnisP,nSnp))
    allocate(ImputePhase(0:nAnisP,nSnp,2))

    Genos(0,:)=9
    GlobalWorkPhase=9
    do i=1,nAnisP
        do j=1,nSnp                                                     ! Phase in the homozygous case
            if (Genos(i,j)==0) GlobalWorkPhase(i,j,:)=0
            if (Genos(i,j)==2) GlobalWorkPhase(i,j,:)=1
        enddo
        ParId=RecPed(i,2)   
        do j=1,nSnp                                                     ! Phase if my father is homozygous
            if (Genos(ParId,j)==0) GlobalWorkPhase(i,j,1)=0
            if (Genos(ParId,j)==2) GlobalWorkPhase(i,j,1)=1
        enddo   
        ParId=RecPed(i,3)   
        do j=1,nSnp                                                     ! Phase if my mother is homozygous
            if (Genos(ParId,j)==0) GlobalWorkPhase(i,j,2)=0
            if (Genos(ParId,j)==2) GlobalWorkPhase(i,j,2)=1
        enddo   

    enddo
    
    GlobalWorkPhase(0,:,:)=9
    ImputeGenos=9
    ImputePhase=9
    ImputeGenos(:,:)=Genos(:,:)
    ImputePhase(:,:,:)=GlobalWorkPhase(:,:,:)

endif   

end subroutine InsteadOfGeneProb

!######################################################################################################################################################################################

subroutine RestrictedWorkLeftRight
use Global

implicit none

integer :: e,i,j,HetEnd,HetStart,WorkRight(nSnp),WorkLeft(nSnp),RSide,LSide,PatMat,SireDamRL
integer :: CountRightSwitch,CountLeftSwitch,StartPt,EndPt,PedId

integer :: StartDis,EndDis,dum,StartJ,k
integer,allocatable,dimension(:) :: TempVec
real,allocatable,dimension(:) :: LengthVec

allocate(TempVec(nSnp))
allocate(LengthVec(nSnp))


ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

MaxLeftRightSwitch=4; MinSpan=200

do i=1,nAnisP
    HetEnd=-1
    HetStart=-1
    do e=1,2
        PatMat=e
        SireDamRL=e+1
        CountLeftSwitch=0
        CountRightSwitch=0
        PedId=RecPed(i,SireDamRL)
        if ((SexOpt==1).and.(RecGender(i)==HetGameticStatus).and.(RecGender(PedId)==HetGameticStatus)) cycle
        if ((PedId>0).and.((float(count(ImputePhase(PedId,:,:)==9))/(2*nSnp))<0.07)) then           !(RecIdHDIndex(PedId)==1)   
            WorkRight=9
            RSide=9
            do j=1,nSnp
                if ((ImputePhase(PedId,j,1)/=ImputePhase(PedId,j,2)).and.&
                        (ImputePhase(PedId,j,1)/=9).and.(ImputePhase(PedId,j,2)/=9))  then
                    HetStart=j
                    if (ImputePhase(i,HetStart,PatMat)==ImputePhase(PedId,HetStart,1)) then
                        WorkRight(HetStart)=1 
                        RSide=1
                        exit
                    endif
                    if (ImputePhase(i,HetStart,PatMat)==ImputePhase(PedId,HetStart,2)) then
                        WorkRight(HetStart)=2
                        RSide=2 
                        exit
                    endif
                endif  
            enddo
                if (RSide/=9) then
                do j=HetStart+1,nSnp
                    if ((ImputePhase(i,j,PatMat)/=ImputePhase(PedId,j,RSide)).and.&
                            (ImputePhase(PedId,j,RSide)/=9).and.(ImputePhase(i,j,PatMat)/=9)) then
                        RSide=abs((RSide-1)-1)+1
                        CountRightSwitch=CountRightSwitch+1
                    endif
                    WorkRight(j)=RSide
                enddo
            endif
            WorkLeft=9
            LSide=9
            do j=nSnp,1,-1
                if ((ImputePhase(PedId,j,1)/=ImputePhase(PedId,j,2)).and.&
                        (ImputePhase(PedId,j,1)/=9).and.(ImputePhase(PedId,j,2)/=9))  then
                    HetEnd=j
                    if (ImputePhase(i,HetEnd,PatMat)==ImputePhase(PedId,HetEnd,1)) then
                        WorkRight(HetEnd)=1 
                        LSide=1
                        exit
                    endif
                    if (ImputePhase(i,HetEnd,PatMat)==ImputePhase(PedId,HetEnd,2)) then
                        WorkRight(HetEnd)=2
                        LSide=2 
                        exit
                    endif
                endif  
            enddo
            if (LSide/=9) then
                do j=HetEnd-1,1,-1
                    if ((ImputePhase(i,j,PatMat)/=ImputePhase(PedId,j,LSide)).and.&
                        (ImputePhase(PedId,j,LSide)/=9).and.(ImputePhase(i,j,PatMat)/=9) ) then
                        LSide=abs((LSide-1)-1)+1
                        CountLeftSwitch=CountLeftSwitch+1
                    endif
                    WorkLeft(j)=LSide
                enddo
            endif
            
!$$$$$$$$$$$$$$$$$$$



            StartDis=-9
            EndDis=-9
            TempVec=9
            LengthVec=0.0
            StartJ=1
            do j=StartJ,nSnp
                if (StartDis==-9) then
                    if (abs(WorkLeft(j)-WorkRight(j))==1) then
                        StartDis=j
                        TempVec(StartDis)=1
                    endif
                endif
                if ((WorkLeft(j)==WorkRight(j)).and.(WorkLeft(j)/=9).and.(StartDis/=-9).and.(EndDis==-9)) then
                    EndDis=j-1
                    TempVec(EndDis)=2
                endif   

                if (StartDis/=-9) then
                    do k=StartDis,1,-1
                        if ((GlobalWorkPhase(PedId,k,1)+GlobalWorkPhase(PedId,k,2))==1) then
                            if (GlobalWorkPhase(i,k,e)/=9) then
                                exit
                            endif
                        endif
                    enddo
                    StartDis=k
                    if (StartDis<1) StartDis=1
                    TempVec(StartDis)=1 
                endif
                if (EndDis/=-9) then
                    do k=EndDis,nSnp
                        if ((GlobalWorkPhase(PedId,k,1)+GlobalWorkPhase(PedId,k,2))==1) then
                            if (GlobalWorkPhase(i,k,e)/=9) then
                                exit
                            endif
                        endif
                    enddo
                    EndDis=k
                    if (EndDis>nSnp) EndDis=nSnp
                    TempVec(EndDis)=2
                endif

                if ((StartDis/=-9).and.(EndDis/=-9)) then
                    TempVec(StartDis+1:EndDis-1)=3
                    LengthVec(StartDis+1:EndDis-1)=1.0/(((EndDis-1)-(StartDis+1))+1)
                    StartJ=EndDis+1
                    StartDis=-9
                    EndDis=-9
                endif   
            enddo
            
            do j=1,nSnp
                if (TempVec(j)==3) then
                    if (ImputePhase(PedId,j,1)/=ImputePhase(PedId,j,2)) then
                        if ((ImputePhase(PedId,j,1)/=9).and.(ImputePhase(PedId,j,2)/=9)) then
                            ImputePhase(i,j,e)=9
                            ImputeGenos(i,j)=9
                        endif
                    endif
                endif   
            enddo
            GlobalWorkPhase(i,:,:)=ImputePhase(i,:,:)

!$$$$$$$$$$$$$$$$$$$            
            
            
            if ((CountLeftSwitch+CountRightSwitch)<(2*MaxLeftRightSwitch)) then
                do j=1,nSnp
                    if (ImputePhase(i,j,PatMat)==9) then
                        if ((WorkLeft(j)==9).and.(WorkRight(j)/=9)) &
                            ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkRight(j))
                        if ((WorkLeft(j)/=9).and.(WorkRight(j)==9)) &
                            ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkLeft(j))
                        if ((WorkLeft(j)/=9).and.(WorkRight(j)==WorkLeft(j))) &
                            ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkLeft(j))
                    endif
                enddo
            else
                EndPt=0
                StartPt=0
                do while ((StartPt<(nSnp-MinSpan)).and.(EndPt<(nSnp-MinSpan)))
                    do j=EndPt+1,nSnp
                        if ((WorkLeft(j)/=9).and.(WorkRight(j)==WorkLeft(j)))  then
                            StartPt=j
                            exit
                        endif
                        if (j==nSnp) StartPt=j
                    enddo
                    do j=StartPt,nSnp
                        if ((WorkLeft(j)==9).or.(WorkRight(j)/=WorkLeft(j)))  then
                            EndPt=j
                            exit
                        endif
                        if (j==nSnp) EndPt=j
                    enddo
                    if (((EndPt-StartPt)+1)>MinSpan) then
                        do j=StartPt,EndPt
                            if (ImputePhase(PedId,j,WorkRight(j))/=9)&
                                ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkRight(j))
                        enddo
                    endif               
                enddo
            endif
        endif
        
    enddo
enddo

do i=1,nAnisP
    if ((SexOpt==1).and.(RecGender(i)==HetGameticStatus)) then
        ImputePhase(i,:,HetGameticStatus)=ImputePhase(i,:,HomGameticStatus)     !JohnHickey changed the j to :
        GlobalWorkPhase(i,:,:)=ImputePhase(i,:,:)
    endif
enddo

ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

end subroutine RestrictedWorkLeftRight

!#############################################################################################################################################################################################################################

subroutine WorkLeftRight
use Global

implicit none

integer :: e,i,j,HetEnd,HetStart,WorkRight(nSnp),WorkLeft(nSnp),RSide,LSide,PatMat,SireDamRL
integer :: CountRightSwitch,CountLeftSwitch,StartPt,EndPt,PedId

integer :: StartDis,EndDis,dum,StartJ,k
integer,allocatable,dimension(:) :: TempVec
real,allocatable,dimension(:) :: LengthVec

allocate(TempVec(nSnp))
allocate(LengthVec(nSnp))

ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

do i=1,nAnisP
    HetEnd=-1
    HetStart=-1
    do e=1,2
        PatMat=e
        SireDamRL=e+1
        CountLeftSwitch=0
        CountRightSwitch=0
        PedId=RecPed(i,SireDamRL)
        if ((SexOpt==1).and.(RecGender(i)==HetGameticStatus).and.(RecGender(PedId)==HetGameticStatus)) cycle
        if (PedId>0) then
            WorkRight=9
            RSide=9
            do j=1,nSnp
                if ((ImputePhase(PedId,j,1)/=ImputePhase(PedId,j,2)).and.&
                        (ImputePhase(PedId,j,1)/=9).and.(ImputePhase(PedId,j,2)/=9))  then
                    HetStart=j
                    if (ImputePhase(i,HetStart,PatMat)==ImputePhase(PedId,HetStart,1)) then
                        WorkRight(HetStart)=1 
                        RSide=1
                        exit
                    endif
                    if (ImputePhase(i,HetStart,PatMat)==ImputePhase(PedId,HetStart,2)) then
                        WorkRight(HetStart)=2
                        RSide=2 
                        exit
                    endif
                endif  
            enddo
                if (RSide/=9) then
                do j=HetStart+1,nSnp
                    if ((ImputePhase(i,j,PatMat)/=ImputePhase(PedId,j,RSide)).and.&
                            (ImputePhase(PedId,j,RSide)/=9).and.(ImputePhase(i,j,PatMat)/=9)) then
                        RSide=abs((RSide-1)-1)+1
                        CountRightSwitch=CountRightSwitch+1
                    endif
                    WorkRight(j)=RSide
                enddo
            endif
            WorkLeft=9
            LSide=9
            do j=nSnp,1,-1
                if ((ImputePhase(PedId,j,1)/=ImputePhase(PedId,j,2)).and.&
                        (ImputePhase(PedId,j,1)/=9).and.(ImputePhase(PedId,j,2)/=9))  then
                    HetEnd=j
                    if (ImputePhase(i,HetEnd,PatMat)==ImputePhase(PedId,HetEnd,1)) then
                        WorkRight(HetEnd)=1 
                        LSide=1
                        exit
                    endif
                    if (ImputePhase(i,HetEnd,PatMat)==ImputePhase(PedId,HetEnd,2)) then
                        WorkRight(HetEnd)=2
                        LSide=2 
                        exit
                    endif
                endif  
            enddo
            if (LSide/=9) then
                do j=HetEnd-1,1,-1
                    if ((ImputePhase(i,j,PatMat)/=ImputePhase(PedId,j,LSide)).and.&
                        (ImputePhase(PedId,j,LSide)/=9).and.(ImputePhase(i,j,PatMat)/=9) ) then
                        LSide=abs((LSide-1)-1)+1
                        CountLeftSwitch=CountLeftSwitch+1
                    endif
                    WorkLeft(j)=LSide
                enddo
            endif
            
            !$$$$$$$$$$$$$$$$$$$



            StartDis=-9
            EndDis=-9
            TempVec=9
            LengthVec=0.0
            StartJ=1
            do j=StartJ,nSnp
                if (StartDis==-9) then
                    if (abs(WorkLeft(j)-WorkRight(j))==1) then
                        StartDis=j
                        TempVec(StartDis)=1
                    endif
                endif
                if ((WorkLeft(j)==WorkRight(j)).and.(WorkLeft(j)/=9).and.(StartDis/=-9).and.(EndDis==-9)) then
                    EndDis=j-1
                    TempVec(EndDis)=2
                endif   

                if (StartDis/=-9) then
                    do k=StartDis,1,-1
                        if ((GlobalWorkPhase(PedId,k,1)+GlobalWorkPhase(PedId,k,2))==1) then
                            if (GlobalWorkPhase(i,k,e)/=9) then
                                exit
                            endif
                        endif
                    enddo
                    StartDis=k
                    if (StartDis<1) StartDis=1
                    TempVec(StartDis)=1 
                endif
                if (EndDis/=-9) then
                    do k=EndDis,nSnp
                        if ((GlobalWorkPhase(PedId,k,1)+GlobalWorkPhase(PedId,k,2))==1) then
                            if (GlobalWorkPhase(i,k,e)/=9) then
                                exit
                            endif
                        endif
                    enddo
                    EndDis=k
                    if (EndDis>nSnp) EndDis=nSnp
                    TempVec(EndDis)=2
                endif

                if ((StartDis/=-9).and.(EndDis/=-9)) then
                    TempVec(StartDis+1:EndDis-1)=3
                    LengthVec(StartDis+1:EndDis-1)=1.0/(((EndDis-1)-(StartDis+1))+1)
                    StartJ=EndDis+1
                    StartDis=-9
                    EndDis=-9
                endif   
            enddo
            
            do j=1,nSnp
                if (TempVec(j)==3) then
                    if (ImputePhase(PedId,j,1)/=ImputePhase(PedId,j,2)) then
                        if ((ImputePhase(PedId,j,1)/=9).and.(ImputePhase(PedId,j,2)/=9)) then
                            ImputePhase(i,j,e)=9
                            ImputeGenos(i,j)=9
                        endif
                    endif
                endif   
            enddo
            GlobalWorkPhase(i,:,:)=ImputePhase(i,:,:)

!$$$$$$$$$$$$$$$$$$$            

            if ((CountLeftSwitch+CountRightSwitch)<(2*MaxLeftRightSwitch)) then
                do j=1,nSnp
                    if (ImputePhase(i,j,PatMat)==9) then
                        if ((WorkLeft(j)==9).and.(WorkRight(j)/=9)) &
                            ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkRight(j))
                        if ((WorkLeft(j)/=9).and.(WorkRight(j)==9)) &
                            ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkLeft(j))
                        if ((WorkLeft(j)/=9).and.(WorkRight(j)==WorkLeft(j))) &
                            ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkLeft(j))
                    endif
                enddo
            else
                EndPt=0
                StartPt=0
                do while ((StartPt<(nSnp-MinSpan)).and.(EndPt<(nSnp-MinSpan)))
                    do j=EndPt+1,nSnp
                        if ((WorkLeft(j)/=9).and.(WorkRight(j)==WorkLeft(j)))  then
                            StartPt=j
                            exit
                        endif
                        if (j==nSnp) StartPt=j
                    enddo
                    do j=StartPt,nSnp
                        if ((WorkLeft(j)==9).or.(WorkRight(j)/=WorkLeft(j)))  then
                            EndPt=j
                            exit
                        endif
                        if (j==nSnp) EndPt=j
                    enddo
                    if (((EndPt-StartPt)+1)>MinSpan) then
                        do j=StartPt,EndPt
                            if (ImputePhase(PedId,j,WorkRight(j))/=9)&
                                ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkRight(j))
                        enddo
                    endif               
                enddo
            endif
        endif
    enddo
enddo

do i=1,nAnisP
    if ((SexOpt==1).and.(RecGender(i)==HetGameticStatus)) then
        !ImputePhase(i,j,HetGameticStatus)=ImputePhase(i,j,HomGameticStatus)
        ImputePhase(i,:,HetGameticStatus)=ImputePhase(i,:,HomGameticStatus)     !JohnHickey changed the j to :
        GlobalWorkPhase(i,:,:)=ImputePhase(i,:,:)
    endif
enddo

ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

end subroutine WorkLeftRight

!#############################################################################################################################################################################################################################

subroutine InternalParentPhaseElim
use Global
implicit none

integer :: m,f,e,h,g,i,j,k,l,nCore,nHap,nGlobalLoop,CoreLength,CoreStart,CoreEnd,InLib,NotHere,CompPhase,Count0
integer :: LoopStart,Offset,AnimalOn(nAnisP,2)
integer :: Count1,Work(nSnp,2),Ban,CompPhasePar,GamA,GamB

integer,allocatable,dimension (:,:) :: CoreIndex,HapLib,LoopIndex,HapElim
integer(kind=1),allocatable,dimension (:,:,:,:) :: Temp

nGlobalLoop=25
allocate(LoopIndex(nGlobalLoop,2))

LoopIndex(1,1)=400
LoopIndex(2,1)=300
LoopIndex(3,1)=250
LoopIndex(4,1)=200
LoopIndex(5,1)=175
LoopIndex(6,1)=150
LoopIndex(7,1)=125
LoopIndex(8,1)=115
LoopIndex(9,1)=100
LoopIndex(10,1)=80
LoopIndex(11,1)=70
LoopIndex(12,1)=60
LoopIndex(13,1)=50
LoopIndex(14,1)=40
LoopIndex(15,1)=35
LoopIndex(16,1)=30
LoopIndex(17,1)=25
LoopIndex(18,1)=20
LoopIndex(19,1)=15
LoopIndex(20,1)=12
LoopIndex(21,1)=10
LoopIndex(22,1)=8
LoopIndex(23,1)=6
LoopIndex(24,1)=5
LoopIndex(25,1)=4

if(nSnp<=50) return
if(nSnp>50) LoopStart=24
if(nSnp>100) LoopStart=21
if(nSnp>200) LoopStart=20
if(nSnp>400) LoopStart=17
if(nSnp>800) LoopStart=17
if(nSnp>1000) LoopStart=15
if(nSnp>1500) LoopStart=13
if(nSnp>2000) LoopStart=10
if(nSnp>3000) LoopStart=6
if(nSnp>4000) LoopStart=3
if(nSnp>5000) LoopStart=1

do i=1,nGlobalLoop
    LoopIndex(i,2)=int(float(nSnp)/LoopIndex(i,1))  
enddo

allocate(Temp(nAnisP,nSnp,2,2))
Temp=0
AnimalOn=0
do m=1,2
        do l=LoopStart,nGlobalLoop
            if (m==1) then
            nCore=nSnp/LoopIndex(l,2) 
                CoreStart=1
                    CoreEnd=LoopIndex(l,2)
            else
                OffSet=int(float(LoopIndex(l,1))/2)
                nCore=(nSnp-(2*OffSet))/LoopIndex(l,2)
                CoreStart=1+Offset
                CoreEnd=LoopIndex(l,2)+Offset
            endif    
            do g=1,nCore
                if ((m==1).and.(g==nCore)) CoreEnd=nSnp
                if ((m==2).and.(g==nCore)) CoreEnd=nSnp-OffSet
            CoreLength=(CoreEnd-CoreStart)+1
            if (CoreLength<10) exit

            !# PARALLEL DO SHARED (nAnisP,RecPed,ImputePhase,CoreStart,CoreEnd,AnimalOn,Temp) private(i,e,CompPhase,GamA,j,GamB)
            do i=1,nAnisP
                do e=1,2
                    if ((SexOpt==1).and.(RecGender(i)==HetGameticStatus).and.(RecGender(RecPed(i,e+1))==HetGameticStatus)) cycle
                    if (RecPed(i,e+1)>0) then
                        CompPhase=1
                        if (count(ImputePhase(i,CoreStart:CoreEnd,e)==9)>0) CompPhase=0
                        if (CompPhase==0) then
                            GamA=1
                            do j=CoreStart,CoreEnd
                                if ((ImputePhase(i,j,e)/=9).and.&
                                    (ImputePhase(RecPed(i,e+1),j,1)&
                                        /=ImputePhase(i,j,e)).and.&
                                    (ImputePhase(RecPed(i,e+1),j,1)/=9)) then       
                                    GamA=0
                                    exit
                                endif
                            enddo
                            GamB=1
                            do j=CoreStart,CoreEnd
                                if ((ImputePhase(i,j,e)/=9).and.&
                                    (ImputePhase(RecPed(i,e+1),j,2)&
                                        /=ImputePhase(i,j,e)).and.&
                                    (ImputePhase(RecPed(i,e+1),j,2)/=9)) then       
                                    GamB=0
                                    exit
                                endif
                            enddo
                            if ((GamA==1).and.(GamB==0)) then
                                AnimalOn(i,e)=1
                                do j=CoreStart,CoreEnd
                                    if (ImputePhase(i,j,e)==9) then
                                        if (ImputePhase(RecPed(i,e+1),j,1)==0)&
                                            Temp(i,j,e,1)=Temp(i,j,e,1)+1
                                        if (ImputePhase(RecPed(i,e+1),j,1)==1)&
                                            Temp(i,j,e,2)=Temp(i,j,e,2)+1
                                    endif
                                enddo
                            endif
                            if ((GamA==0).and.(GamB==1)) then
                                AnimalOn(i,e)=1
                                do j=CoreStart,CoreEnd
                                    if (ImputePhase(i,j,e)==9) then
                                        if (ImputePhase(RecPed(i,e+1),j,2)==0)&
                                            Temp(i,j,e,1)=Temp(i,j,e,1)+1
                                        if (ImputePhase(RecPed(i,e+1),j,2)==1)&
                                            Temp(i,j,e,2)=Temp(i,j,e,2)+1
                                    endif
                                enddo
                            endif
                        endif
                    endif
                enddo
            enddo
            !# END PARALLEL DO 
            
            CoreStart=CoreStart+LoopIndex(l,2)
            CoreEnd=CoreEnd+LoopIndex(l,2)
            if ((m==2).and.(g==nCore)) exit
        enddo
    enddo
enddo

do i=1,nAnisP
    do e=1,2
        if (AnimalOn(i,e)==1) then
            if ((SexOpt==0).or.(RecGender(i)==HomGameticStatus)) then
                do j=1,nSnp
                    if (ImputePhase(i,j,e)==9) then
                        if ((Temp(i,j,e,1)>nAgreeInternalHapLibElim).and.(Temp(i,j,e,2)==0)) ImputePhase(i,j,e)=0
                        if ((Temp(i,j,e,1)==0).and.(Temp(i,j,e,2)>nAgreeInternalHapLibElim)) ImputePhase(i,j,e)=1
                    endif
                enddo
            endif
        endif
    enddo
    if ((SexOpt==1).and.(RecGender(i)==HetGameticStatus)) then
        if (AnimalOn(i,HomGameticStatus)==1) then
            do j=1,nSnp
                if (ImputePhase(i,j,HomGameticStatus)==9) then
                    if ((Temp(i,j,HomGameticStatus,1)>nAgreeInternalHapLibElim).and.(Temp(i,j,HomGameticStatus,2)==0)) ImputePhase(i,j,:)=0
                    if ((Temp(i,j,HomGameticStatus,1)==0).and.(Temp(i,j,HomGameticStatus,2)>nAgreeInternalHapLibElim)) ImputePhase(i,j,:)=1
                endif
            enddo
        endif
    endif   
enddo
deallocate(Temp)

ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

end subroutine InternalParentPhaseElim

!#############################################################################################################################################################################################################################

subroutine InternalHapLibImputation
use Global
implicit none

integer :: f,e,h,g,i,j,k,l,nCore,nHap,nGlobalLoop,CoreLength,CoreStart,CoreEnd,InLib,NotHere,CompPhase,Count0,Count1,Work(nSnp,2)
integer :: Counter,BanBoth(2),Ban(2),AnimalOn(nAnisP,2)
integer :: LoopStart,OffSet

integer,allocatable,dimension (:,:) :: CoreIndex,HapLib,LoopIndex,HapElim
integer(kind=1),allocatable,dimension (:,:,:,:) :: Temp

nGlobalLoop=25
allocate(LoopIndex(nGlobalLoop,2))

LoopIndex(1,1)=400
LoopIndex(2,1)=300
LoopIndex(3,1)=250
LoopIndex(4,1)=200
LoopIndex(5,1)=175
LoopIndex(6,1)=150
LoopIndex(7,1)=125
LoopIndex(8,1)=115
LoopIndex(9,1)=100
LoopIndex(10,1)=80
LoopIndex(11,1)=70
LoopIndex(12,1)=60
LoopIndex(13,1)=50
LoopIndex(14,1)=40
LoopIndex(15,1)=35
LoopIndex(16,1)=30
LoopIndex(17,1)=25
LoopIndex(18,1)=20
LoopIndex(19,1)=15
LoopIndex(20,1)=12
LoopIndex(21,1)=10
LoopIndex(22,1)=8
LoopIndex(23,1)=6
LoopIndex(24,1)=5
LoopIndex(25,1)=4

if(nSnp<=50) return
if(nSnp>50) LoopStart=24
if(nSnp>100) LoopStart=21
if(nSnp>200) LoopStart=20
if(nSnp>400) LoopStart=17
if(nSnp>800) LoopStart=17
if(nSnp>1000) LoopStart=15
if(nSnp>1500) LoopStart=13
if(nSnp>2000) LoopStart=10
if(nSnp>3000) LoopStart=6
if(nSnp>4000) LoopStart=3
if(nSnp>5000) LoopStart=1

do i=LoopStart,nGlobalLoop
    LoopIndex(i,2)=int(float(nSnp)/LoopIndex(i,1))  
enddo

allocate(Temp(nAnisP,nSnp,2,2))
Temp=0
AnimalOn=0
do f=1,2
    allocate(HapLib(nAnisP*2,nSnp))
        do l=LoopStart,nGlobalLoop
            if (f==1) then
            nCore=nSnp/LoopIndex(l,2) 
                CoreStart=1
                    CoreEnd=LoopIndex(l,2)
            else
                OffSet=int(float(LoopIndex(l,1))/2)
                nCore=(nSnp-(2*OffSet))/LoopIndex(l,2)
                CoreStart=1+Offset
                CoreEnd=LoopIndex(l,2)+Offset
            endif    
            do g=1,nCore
                if ((f==1).and.(g==nCore)) CoreEnd=nSnp
                if ((f==2).and.(g==nCore)) CoreEnd=nSnp-OffSet
            CoreLength=(CoreEnd-CoreStart)+1
            if (CoreLength<10) exit
            nHap=0

            HapLib=9
            
            !# PARALLEL DO SHARED (nAnisP,ImputePhase,CoreStart,CoreEnd,nHap,HapLib) private(i,e,CompPhase,InLib,h,NotHere,j)
            do i=1,nAnisP
                do e=1,2
                    if ((ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle 
                    CompPhase=1
                    if (count(ImputePhase(i,CoreStart:CoreEnd,e)==9)>0) CompPhase=0
                    if (CompPhase==1) then
                        if (nHap==0) then
                            HapLib(1,CoreStart:CoreEnd)=ImputePhase(i,CoreStart:CoreEnd,e)
                            nHap=1
                        else
                            InLib=0
                            do h=1,nHap
                                NotHere=1
                                do j=CoreStart,CoreEnd
                                    if(ImputePhase(i,j,e)/=HapLib(h,j)) then
                                        NotHere=0   
                                        exit
                                    endif
                                enddo
                                if (NotHere==1) then
                                    InLib=1
                                    exit                
                                endif
                            enddo
                            if (InLib==0) then
                                nHap=nHap+1             
                                HapLib(nHap,CoreStart:CoreEnd)=ImputePhase(i,CoreStart:CoreEnd,e)           
                            endif
                        endif
                    endif
                enddo
            enddo
            !# END PARALLEL DO 
            
            !# PARALLEL DO SHARED (nAnisP,ImputePhase,CoreStart,CoreEnd,CoreLength,nHap,HapLib,ImputeGenos,AnimalOn,Temp) private(i,HapElim,Work,BanBoth,e,h,j,Counter,Count0,Count1,Ban)
            do i=1,nAnisP
                allocate(HapElim(nAnisP*2,2))
                HapElim=1
                Work=9
                BanBoth=0
                do e=1,2
                    if ((ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle 
                    if ((count(ImputePhase(i,CoreStart:CoreEnd,e)==9)/=CoreLength)&
                        .and.(count(ImputePhase(i,CoreStart:CoreEnd,e)/=9)/=CoreLength)) then
                        do h=1,nHap
                            do j=CoreStart,CoreEnd
                                if ((ImputePhase(i,j,e)/=9)&
                                        .and.(ImputePhase(i,j,e)/=HapLib(h,j))) then
                                    HapElim(h,e)=0
                                    exit    
                                endif
                            enddo
                        enddo
                        Counter=count(HapElim(1:nHap,e)==1)
                        if (float(Counter)<(float(nHap)*0.25)) then
                            BanBoth(e)=1
                            do j=CoreStart,CoreEnd
                                Count0=0
                                Count1=0
                                do h=1,nHap
                                    if (HapElim(h,e)==1) then
                                        if (HapLib(h,j)==0) Count0=Count0+1
                                        if (HapLib(h,j)==1) Count1=Count1+1
                                        if ((Count0>0).and.(Count1>0)) exit
                                    endif                       
                                enddo
                                if (Count0>0) then
                                    if (Count1==0) Work(j,e)=0
                                else
                                    if (Count1>0) Work(j,e)=1
                                endif
                            enddo
                        endif
                    endif
                enddo
                Ban=0
                if (BanBoth(1)==1) Ban(1)=1
                if (BanBoth(2)==1) Ban(2)=1
                if (sum(BanBoth(:))==2) then
                    do j=CoreStart,CoreEnd
                        if (ImputeGenos(i,j)/=9) then
                            if ((Work(j,1)/=9).and.(Work(j,2)/=9)) then
                                if (ImputeGenos(i,j)/=(Work(j,1)+Work(j,2))) then
                                    Ban=0
                                    exit
                                endif
                            endif
                        endif
                    enddo 
                endif   
                do e=1,2
                    if ((ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle 
                    if (Ban(e)==1) then
                        AnimalOn(i,e)=1
                        do j=CoreStart,CoreEnd
                            if (Work(j,e)==0) Temp(i,j,e,1)=Temp(i,j,e,1)+1
                            if (Work(j,e)==1) Temp(i,j,e,2)=Temp(i,j,e,2)+1
                        enddo
                    endif
                enddo   
                deallocate(HapElim)
            enddo
            !# END PARALLEL DO
            
            CoreStart=CoreStart+LoopIndex(l,2)
            CoreEnd=CoreEnd+LoopIndex(l,2)
            if ((f==2).and.(g==nCore)) exit
        enddo
    enddo
    deallocate(HapLib)
enddo

do i=1,nAnisP
    do e=1,2
        if ((ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle 
        if (AnimalOn(i,e)==1) then
            do j=1,nSnp
                if (ImputePhase(i,j,e)==9) then
                    if ((Temp(i,j,e,1)>nAgreeInternalHapLibElim).and.(Temp(i,j,e,2)==0)) ImputePhase(i,j,e)=0
                    if ((Temp(i,j,e,1)==0).and.(Temp(i,j,e,2)>nAgreeInternalHapLibElim)) ImputePhase(i,j,e)=1
                endif
            enddo
        endif   
    enddo
enddo
deallocate(Temp)

ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

end subroutine InternalHapLibImputation

!#############################################################################################################################################################################################################################

subroutine PhaseElimination
use Global
use GlobalPedigree
implicit none

integer :: e,g,h,i,j,nCore,CoreLength,dum,GamA,GamB,nAnisHD,PosHDInd
integer :: StartSnp,EndSnp,Gam1,Gam2,TempCount,AnimalOn(nAnisP,2)
integer,allocatable,dimension (:) :: PosHD
integer,allocatable,dimension (:,:) :: CoreIndex
integer,allocatable,dimension (:,:,:) :: PhaseHD
integer(kind=1),allocatable,dimension (:,:,:,:) :: Temp

character(len=1000) :: FileName,dumC

nAnisHD=(count(Setter(:)==1))

allocate(Temp(nAnisP,nSnp,2,2))
allocate(PhaseHD(nAnisHD,nSnp,2))
allocate(PosHD(nAnisP))
PosHD=0
Temp=0
AnimalOn=0
do h=1,nPhaseInternal
    if (ManagePhaseOn1Off0==0) then 
        write (FileName,'(a,"Phase",i0,"/PhasingResults/CoreIndex.txt")') trim(PhasePath),h
    else
        write (FileName,'("./Phasing/Phase",i0,"/PhasingResults/CoreIndex.txt")')h  
    endif   
    call CountLines(FileName,nCore)
    allocate(CoreIndex(nCore,2))
    open (unit=2001,file=trim(FileName),status="old")
    do g=1,nCore
        read (2001,*) dum,CoreIndex(g,:)
    enddo
    close(2001)
    if (ManagePhaseOn1Off0==0) then 
        write (FileName,'(a,"Phase",i0,"/PhasingResults/FinalPhase.txt")') trim(PhasePath),h
    else
        write (FileName,'("./Phasing/Phase",i0,"/PhasingResults/FinalPhase.txt")')h 
    endif   
    open (unit=2001,file=trim(FileName),status="old")
    do i=1,nAnisHD
        read (2001,*) dumC,PhaseHD(i,:,1)
        read (2001,*) dumC,PhaseHD(i,:,2)
        do j=1,nAnisP
            if (trim(dumC)==trim(Id(j))) then
                PosHD(j)=i
                exit
            endif
        enddo
    enddo
    close(2001) 
    do g=1,nCore
        StartSnp=CoreIndex(g,1)
        EndSnp=CoreIndex(g,2)
        !# PARALLEL DO SHARED (nAnisP,PosHD,RecPed,ImputePhase,StartSnp,EndSnp,PhaseHD,AnimalOn,Temp) private(i,PosHDInd,Gam1,Gam2,e,GamA,GamB,TempCount,j)
        do i=1,nAnisP
            PosHDInd=PosHD(RecPed(i,1))
            if ((count(ImputePhase(i,StartSnp:EndSnp,:)==9)/=0).and.(PosHDInd>0)) then
                if (count(ImputePhase(i,StartSnp:EndSnp,1)/=ImputePhase(i,StartSnp:EndSnp,2))>0) then
                    Gam1=0
                    Gam2=0
                    do e=1,2
                        GamA=1
                        GamB=1
                        TempCount=0
                        do j=StartSnp,EndSnp                
                            if (ImputePhase(i,j,e)/=9) then
                                if (ImputePhase(i,j,e)/=PhaseHD(PosHDInd,j,1)) then 
                                    TempCount=TempCount+1
                                    if (ImputeFromHDPhaseThresh==TempCount) then
                                        GamA=0
                                        exit
                                    endif
                                endif
                            endif
                        enddo
                        TempCount=0
                        do j=StartSnp,EndSnp
                            if (ImputePhase(i,j,e)/=9) then
                                if (ImputePhase(i,j,e)/=PhaseHD(PosHDInd,j,2)) then
                                    TempCount=TempCount+1
                                    if (ImputeFromHDPhaseThresh==TempCount) then
                                        GamB=0
                                        exit
                                    endif
                                endif
                            endif
                        enddo
                        if ((e==1).and.(GamA==1).and.(GamB==0)) Gam1=1
                        if ((e==1).and.(GamA==0).and.(GamB==1)) Gam1=2
                        if ((e==2).and.(GamA==1).and.(GamB==0)) Gam2=1
                        if ((e==2).and.(GamA==0).and.(GamB==1)) Gam2=2
                    enddo
                    if (Gam1/=Gam2) then
                        AnimalOn(i,:)=1
                        if (Gam1/=0) then
                            do j=StartSnp,EndSnp
                                if (ImputePhase(i,j,1)==9) then
                                    if(PhaseHD(PosHDInd,j,Gam1)==0)&
                                                 Temp(i,j,1,1)=Temp(i,j,1,1)+1
                                    if(PhaseHD(PosHDInd,j,Gam1)==1)&
                                                 Temp(i,j,1,2)=Temp(i,j,1,2)+1  
                                endif
                            enddo
                        endif
                        if (Gam2/=0) then
                            do j=StartSnp,EndSnp
                                if (ImputePhase(i,j,2)==9) then
                                    if(PhaseHD(PosHDInd,j,Gam2)==0)&
                                                 Temp(i,j,2,1)=Temp(i,j,2,1)+1
                                    if(PhaseHD(PosHDInd,j,Gam2)==1)&
                                                 Temp(i,j,2,2)=Temp(i,j,2,2)+1  
                                endif
                            enddo
                        endif
                    endif
                endif
            endif
        enddo
        !# END PARALLEL DO
    enddo
    deallocate(CoreIndex)
enddo

do i=1,nAnisP
    do e=1,2
        if (AnimalOn(i,e)==1) then
            do j=1,nSnp
                if (ImputePhase(i,j,e)==9) then
                    if ((Temp(i,j,e,1)>nAgreeInternalHapLibElim).and.(Temp(i,j,e,2)==0)) ImputePhase(i,j,e)=0
                    if ((Temp(i,j,e,1)==0).and.(Temp(i,j,e,2)>nAgreeInternalHapLibElim)) ImputePhase(i,j,e)=1
                endif
            enddo
        endif
    enddo
enddo
deallocate(Temp)    

ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

end subroutine PhaseElimination

!#############################################################################################################################################################################################################################

subroutine ParentPhaseElimination
use Global
use GlobalPedigree
implicit none

integer :: e,g,h,i,j,nCore,CoreLength,dum,PedId,GamA,GamB,nAnisHD,PosHDInd
integer :: StartSnp,EndSnp,TempCount,AnimalOn(nAnisP,2)
integer,allocatable,dimension (:) :: PosHD
integer,allocatable,dimension (:,:) :: CoreIndex
integer,allocatable,dimension (:,:,:) :: PhaseHD
integer(kind=1),allocatable,dimension (:,:,:,:) :: Temp

character(len=1000) :: FileName,dumC

nAnisHD=(count(Setter(:)==1))

allocate(Temp(nAnisP,nSnp,2,2))
allocate(PhaseHD(nAnisHD,nSnp,2))
allocate(PosHD(nAnisP))
PosHD=0
Temp=0
AnimalOn=0
do h=1,nPhaseInternal
    if (ManagePhaseOn1Off0==0) then 
        write (FileName,'(a,"Phase",i0,"/PhasingResults/CoreIndex.txt")') trim(PhasePath),h
    else
        write (FileName,'("./Phasing/Phase",i0,"/PhasingResults/CoreIndex.txt")')h  
    endif   
    call CountLines(FileName,nCore)
    allocate(CoreIndex(nCore,2))
    open (unit=2001,file=trim(FileName),status="old")
    do g=1,nCore
        read (2001,*) dum,CoreIndex(g,:)
    enddo
    close(2001)
    if (ManagePhaseOn1Off0==0) then 
        write (FileName,'(a,"Phase",i0,"/PhasingResults/FinalPhase.txt")') trim(PhasePath),h
    else
        write (FileName,'("./Phasing/Phase",i0,"/PhasingResults/FinalPhase.txt")')h 
    endif   
    open (unit=2001,file=trim(FileName),status="old")
    do i=1,nAnisHD
        read (2001,*) dumC,PhaseHD(i,:,1)
        read (2001,*) dumC,PhaseHD(i,:,2)
        do j=1,nAnisP
            if (trim(dumC)==trim(Id(j))) then
                PosHD(j)=i
                exit
            endif
        enddo
    enddo
    close(2001) 
    do g=1,nCore
        StartSnp=CoreIndex(g,1)
        EndSnp=CoreIndex(g,2)
        
        !# PARALLEL DO SHARED (RecPed,PosHD,ImputePhase,StartSnp,EndSnp,PhaseHD,AnimalOn,Temp) private(i,e,PedId,PosHDInd,GamA,GamB,TempCount,j)
        do i=1,nAnisP
            do e=1,2
                PedId=e+1
                if ((SexOpt==1).and.(RecGender(i)==HetGameticStatus).and.(RecGender(RecPed(i,PedId))==HetGameticStatus)) cycle
                if (RecPed(i,PedId)>0) then
                    PosHDInd=PosHD(RecPed(i,PedId))
                    if ((count(ImputePhase(i,StartSnp:EndSnp,e)==9)/=0).and.(PosHDInd>0)) then
                        GamA=1
                        GamB=1
                        TempCount=0
                        do j=StartSnp,EndSnp
                            if (ImputePhase(i,j,e)/=9) then
                                if (ImputePhase(i,j,e)/=PhaseHD(PosHDInd,j,1)) then
                                    TempCount=TempCount+1
                                    if (ImputeFromParentCountThresh==TempCount) then
                                        GamA=0  
                                        exit
                                    endif
                                endif
                            endif
                        enddo
                        TempCount=0
                        do j=StartSnp,EndSnp
                            if (ImputePhase(i,j,e)/=9) then
                                if (ImputePhase(i,j,e)/=PhaseHD(PosHDInd,j,2)) then
                                    TempCount=TempCount+1
                                    if (ImputeFromParentCountThresh==TempCount) then
                                        GamB=0  
                                        exit
                                    endif
                                endif
                            endif
                        enddo
                        if ((GamA==1).and.(GamB==0)) then
                            AnimalOn(i,e)=1
                            do j=StartSnp,EndSnp
                                if (ImputePhase(i,j,e)==9) then
                                    if(PhaseHD(PosHDInd,j,1)==0) Temp(i,j,e,1)=Temp(i,j,e,1)+1
                                    if(PhaseHD(PosHDInd,j,1)==1) Temp(i,j,e,2)=Temp(i,j,e,2)+1  
                                endif
                            enddo
                        endif
                        if ((GamA==0).and.(GamB==1)) then
                            AnimalOn(i,e)=1
                            do j=StartSnp,EndSnp
                                if (ImputePhase(i,j,e)==9) then
                                    if(PhaseHD(PosHDInd,j,2)==0) Temp(i,j,e,1)=Temp(i,j,e,1)+1
                                    if(PhaseHD(PosHDInd,j,2)==1) Temp(i,j,e,2)=Temp(i,j,e,2)+1  
                                endif
                            enddo
                        endif
                    endif
                endif
            enddo
        enddo
        !# END PARALLEL DO 
        
    enddo
    deallocate(CoreIndex)
enddo

do i=1,nAnisP
    do e=1,2
        if (AnimalOn(i,e)==1) then
            if ((SexOpt==0).or.(RecGender(i)==HomGameticStatus)) then
                do j=1,nSnp
                    if (ImputePhase(i,j,e)==9) then
                        if ((Temp(i,j,e,1)>nAgreeInternalHapLibElim).and.(Temp(i,j,e,2)==0))&
                             ImputePhase(i,j,e)=0
                        if ((Temp(i,j,e,1)==0).and.(Temp(i,j,e,2)>nAgreeInternalHapLibElim))&
                             ImputePhase(i,j,e)=1
                    endif
                enddo
            endif
        endif
    enddo
    if ((SexOpt==1).and.(RecGender(i)==HetGameticStatus)) then
        if (AnimalOn(i,HomGameticStatus)==1) then
            do j=1,nSnp
                if (ImputePhase(i,j,HomGameticStatus)==9) then
                    if ((Temp(i,j,HomGameticStatus,1)>nAgreeInternalHapLibElim).and.(Temp(i,j,HomGameticStatus,2)==0)) ImputePhase(i,j,:)=0
                    if ((Temp(i,j,HomGameticStatus,1)==0).and.(Temp(i,j,HomGameticStatus,2)>nAgreeInternalHapLibElim)) ImputePhase(i,j,:)=1
                endif
            enddo
        endif
    endif   
enddo
deallocate(Temp)

ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

end subroutine ParentPhaseElimination

!#############################################################################################################################################################################################################################

subroutine ImputeFromHDLibrary
use Global

implicit none

integer :: l,i,j,k,h,e,f,g,nCore,dum,CoreLength,nHap,CountAB(nSnp,0:1),Work(nSnp,2),TempCount
integer :: StartSnp,EndSnp,PatMatDone(2),Counter,BanBoth(2),Ban(2),AnimalOn(nAnisP,2)
integer,allocatable,dimension (:,:,:,:) :: Temp
integer,allocatable,dimension (:,:) :: CoreIndex
integer(kind=1),allocatable,dimension (:,:) :: HapLib,HapCand

character(len=1000) :: FileName,dumC

allocate(Temp(0:nAnisP,nSnp,2,2))
Temp=0
AnimalOn=0
do h=1,nPhaseInternal
    if (ManagePhaseOn1Off0==0) then 
        write (FileName,'(a,"Phase",i0,"/PhasingResults/CoreIndex.txt")') trim(PhasePath),h
    else
        write (FileName,'("./Phasing/Phase",i0,"/PhasingResults/CoreIndex.txt")')h  
    endif   
    call CountLines(FileName,nCore)
    allocate(CoreIndex(nCore,2))
    open (unit=2001,file=trim(FileName),status="old")
    do g=1,nCore
        read (2001,*) dum,CoreIndex(g,:)
    enddo
    close(2001)
    do g=1,nCore
        StartSnp=CoreIndex(g,1)
        EndSnp=CoreIndex(g,2)
        CoreLength=(EndSnp-StartSnp)+1
        if (ManagePhaseOn1Off0==0) then 
            write (FileName,'(a,"Phase",i0,"/PhasingResults/HaplotypeLibrary/HapLib",i0,".bin")') trim(PhasePath),h,g
        else
            write (FileName,'("./Phasing/Phase",i0,"/PhasingResults/HaplotypeLibrary/HapLib",i0,".bin")')h,g
        endif   
        open (unit=2001,file=trim(FileName),status="old",form="unformatted")
        read(2001) nHap,CoreLength
        if(nHap/=0) then    
            allocate(HapLib(nHap,CoreLength))
            do l=1,nHap
                    read(2001) HapLib(l,:)
            enddo
            close (2001)

            !$OMP PARALLEL DO SHARED (nAnisP,ImputePhase,StartSnp,EndSnp,CoreLength,nHap,HapLib) private(i,e,f,j,k,TempCount,CountAB,Counter,PatMatDone,Work,BanBoth,Ban,HapCand)
            do i=1,nAnisP
                allocate(HapCand(nHap,2))   
                PatMatDone=0
                if (count(ImputePhase(i,StartSnp:EndSnp,1)==9)/=CoreLength) PatMatDone(1)=1
                if (count(ImputePhase(i,StartSnp:EndSnp,2)==9)/=CoreLength) PatMatDone(2)=1
                HapCand=1
                Work=9
                BanBoth=0
                do e=1,2
                    if ((ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle 
                    if ((PatMatDone(e)==1).and.(count(ImputePhase(i,StartSnp:EndSnp,e)/=9)/=CoreLength)) then
                        do f=1,nHap
                            k=0
                            TempCount=0
                            do j=StartSnp,EndSnp
                                k=k+1
                                if ((ImputePhase(i,j,e)/=9).and.(ImputePhase(i,j,e)/=HapLib(f,k))) then
                                    TempCount=TempCount+1
                                    if (ImputeFromHDLibraryCountThresh==TempCount) then
                                        HapCand(f,e)=0
                                        exit
                                    endif   
                                endif
                            enddo
                        enddo
                        CountAB=0
                        Counter=count(HapCand(:,e)==1)
                        if (float(Counter)<(float(nHap)*0.25)) then
                            BanBoth(e)=1
                            do f=1,nHap
                                if (HapCand(f,e)==1) then
                                    k=0
                                    do j=StartSnp,EndSnp
                                        k=k+1
                                        CountAB(j,HapLib(f,k))=CountAB(j,HapLib(f,k))+1     
                                    enddo
                                endif
                            enddo
                            do j=StartSnp,EndSnp
                                if (CountAB(j,0)>0) then
                                    if (CountAB(j,1)==0) Work(j,e)=0
                                else
                                    if (CountAB(j,1)>0) Work(j,e)=1
                                endif
                            enddo
                        endif
                    endif
                enddo
                if (sum(PatMatDone(:))>0) then
                    Ban=0   
                    if (BanBoth(1)==1) Ban(1)=1
                    if (BanBoth(2)==1) Ban(2)=1
                    if (sum(BanBoth(:))==2) then
                        TempCount=0
                        do j=StartSnp,EndSnp
                            if (ImputeGenos(i,j)/=9) then
                                if (ImputeGenos(i,j)/=(Work(j,1)+Work(j,2))) then
                                    TempCount=TempCount+1
                                    if (ImputeFromHDLibraryCountThresh==TempCount) then
                                        Ban=0
                                        exit
                                    endif
                                endif
                            endif
                        enddo 
                    endif   
                    do e=1,2
                        if ((ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle 
                        if (Ban(e)==1) then
                            AnimalOn(i,e)=1
                            do j=StartSnp,EndSnp
                                if (Work(j,e)==0) Temp(i,j,e,1)=Temp(i,j,e,1)+1
                                if (Work(j,e)==1) Temp(i,j,e,2)=Temp(i,j,e,2)+1
                            enddo
                        endif   
                    enddo   
                endif
                deallocate(HapCand)
            enddo
            !$OMP END PARALLEL DO 

            deallocate(HapLib)
        endif
    enddo
    deallocate(CoreIndex)
enddo

do i=1,nAnisP
    do e=1,2
        if ((ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle 
        if (AnimalOn(i,e)==1) then
            do j=1,nSnp
                if (ImputePhase(i,j,e)==9) then
                    if ((Temp(i,j,e,1)>nAgreeInternalHapLibElim).and.(Temp(i,j,e,2)==0))&
                         ImputePhase(i,j,e)=0
                    if ((Temp(i,j,e,1)==0).and.(Temp(i,j,e,2)>nAgreeInternalHapLibElim))&
                         ImputePhase(i,j,e)=1
                endif
            enddo
        endif
    enddo
enddo
deallocate(Temp)

ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

end subroutine ImputeFromHDLibrary

!#############################################################################################################################################################################################################################

subroutine BaseAnimalFillIn
! Impute phase for Base Animals (without any pedigree information).
! The internal phasing has been calculated by means of AlphaPhase using both 'with' and 'without' shift parameter. 
! The function select the middle phase run without shifting the cores and the middle core within this phase run, and then 
! selects the complementary phase run with the shift cores. Shifted and non shifted cores allows overlapping and so the
! function goes through SNPs from left to right and from right to left from the MiddleCore jumping to shifted and non shifted
! cores.
! The function ends when all the SNPs have been phased. If the haplotypes have not been fully phased means that there is some
! recombination in that haplotype.

use Global
use GlobalPedigree

implicit none

integer :: e,h,i,g,l,j,nCoreA,nCoreB,MiddlePhaseRun,dum,MiddleCoreA,MiddleCoreB,Ban,nHap,CoreLength,nAnisHD,PosHDInd,CountDisagree
integer :: CompPhaseRun,CompJump,StartSnp,EndSnp,UptoRightSnp,UptoLeftSnp,UpToCoreA,UpToCoreB,C1,C2,C3,C4,Recmb,CompLength,RL
integer :: UpToSnp,StPt,EndPt,FillInSt,FillInEnd
integer,allocatable,dimension (:) :: PosHD
integer,allocatable,dimension (:,:) :: CoreIndexA,CoreIndexB,AnimRecomb
integer,allocatable,dimension (:,:,:,:) :: PhaseHD
character(len=1000) :: FileName,dumC

nAnisHD=(count(Setter(:)==1))

allocate(PhaseHD(nAnisHD,nSnp,2,2))     ! HIGH DENSITY PHASING: PhaseHD = (Animals, SNPs, Haplotypes, Nonshifted and Shifted phasing)
allocate(AnimRecomb(nAnisP,2))
allocate(PosHD(nAnisP))
AnimRecomb=0
PosHD=0

! Select the phasing in the middle    
MiddlePhaseRun=(int(nPhaseInternal)/4)
if (MiddlePhaseRun==0) MiddlePhaseRun=1
CompJump=int(nPhaseInternal)/2
if (CompJump==0) CompJump=1             ! This is only possible if no phasing info is available
CompPhaseRun=MiddlePhaseRun+CompJump

! The different "PhasingResults/CoreIndex.txt" files are created 
! by AlphaPhase when it is called during PhasingManagement

! Get cores corresponding to this phasing step
if (ManagePhaseOn1Off0==0) then
    ! Phasing has been done independently to AlphaImpute
    write (FileName,'(a,"Phase",i0,"/PhasingResults/CoreIndex.txt")') trim(PhasePath),MiddlePhaseRun
else
    ! Phasing has been managed by previous calls to AlphaImpute
    write (FileName,'("./Phasing/Phase",i0,"/PhasingResults/CoreIndex.txt")') MiddlePhaseRun 
endif
call CountLines(FileName,nCoreA)
allocate(CoreIndexA(nCoreA,2))
! Get the start and end of cores of this phasing step
open (unit=2001,file=trim(FileName),status="old")
do g=1,nCoreA
    read (2001,*) dum,CoreIndexA(g,:)
enddo
close(2001)
! Select the core in the middle
MiddleCoreA=int(nCoreA)/2
if (MiddleCoreA==0) MiddleCoreA=1

! Get the core corresponding to the complementary phasing step
if (ManagePhaseOn1Off0==0) then 
    write (FileName,'(a,"Phase",i0,"/PhasingResults/CoreIndex.txt")') trim(PhasePath),CompPhaseRun
else
    write (FileName,'("./Phasing/Phase",i0,"/PhasingResults/CoreIndex.txt")')CompPhaseRun   
endif   
call CountLines(FileName,nCoreB)
allocate(CoreIndexB(nCoreB,2))
open (unit=2001,file=trim(FileName),status="old")
do g=1,nCoreB
    read (2001,*) dum,CoreIndexB(g,:)
enddo
close(2001)
MiddleCoreB=int(nCoreB)/2
if (MiddleCoreB==0) MiddleCoreB=1

! Get HIGH DENSITY phase information of this phasing step
! WARNING: If I only want to phase base animals, why do I need to read the whole file?
if (ManagePhaseOn1Off0==0) then 
    write (FileName,'(a,"Phase",i0,"/PhasingResults/FinalPhase.txt")') trim(PhasePath),MiddlePhaseRun
else
    write (FileName,'("./Phasing/Phase",i0,"/PhasingResults/FinalPhase.txt")')MiddlePhaseRun    
endif   
open (unit=2001,file=trim(FileName),status="old")
do i=1,nAnisHD
    read (2001,*) dumC,PhaseHD(i,:,1,1)
    read (2001,*) dumC,PhaseHD(i,:,2,1)
    do j=1,nAnisP
        ! Match individuals with high density phase information in the FinalPhase output file
        if (trim(dumC)==trim(Id(j))) then
            PosHD(j)=i
            exit
        endif
    enddo
enddo
close(2001) 

if (ManagePhaseOn1Off0==0) then 
    write (FileName,'(a,"Phase",i0,"/PhasingResults/FinalPhase.txt")') trim(PhasePath),CompPhaseRun
else
    write (FileName,'("./Phasing/Phase",i0,"/PhasingResults/FinalPhase.txt")')CompPhaseRun  
endif   
open (unit=2001,file=trim(FileName),status="old")
do i=1,nAnisHD
    read (2001,*) dumC,PhaseHD(i,:,1,2)
    read (2001,*) dumC,PhaseHD(i,:,2,2)
    ! It is assumed that there is no difference in both FinalPhase files in terms of animals 
enddo
close(2001) 

! Impute HD phase of the middle core of the middle phasing step
! WARNING: Why to impute phase information only for this case?
do g=MiddleCoreA,MiddleCoreA        ! Why is it necessary a loop here?
    StartSnp=CoreIndexA(g,1)
    EndSnp=CoreIndexA(g,2)
    CoreLength=(EndSnp-StartSnp)+1
    do i=1,nAnisP
        ! If I have no parents and if I am somebody
        if ((BaseAnimals(i)==1).and.(PosHD(i)/=0)) then
            CountDisagree=0
            ! Check if the two haplotypes are equal
            do j=StartSnp,EndSnp 
                if (ImputePhase(i,j,1)/=ImputePhase(i,j,2)) then
                    CountDisagree=CountDisagree+1
                    if (CountDisagree>1) exit
                endif
            enddo
            ! If haplotypes are equal
            ! Impute High Density phase
            if (CountDisagree==0) then
                ImputePhase(i,StartSnp:EndSnp,1)=PhaseHD(PosHD(i),StartSnp:EndSnp,1,1)
                ImputePhase(i,StartSnp:EndSnp,2)=PhaseHD(PosHD(i),StartSnp:EndSnp,2,1)
            endif
        endif
    enddo
enddo

UpToRightSnp=EndSnp
UpToLeftSnp=StartSnp
UpToCoreA=MiddleCoreA

! Go through SNPs from left to right (1) and from right to left (2) from the MiddleCore
! The internal phasing are calculated both with and without shift what allows to have 
! two different cores overlapping.
! e variable controls in which direction the phasing is being performed:
!   - e=1 Left -> Right
!   - e=2 Right -> Left
! h variable swaps between the two different phasing steps with and without shift, so the 
! the imputation can go on through all the SNPs
do e=1,2
    if (e==1) then              ! Going forward (from left to right)
        UpToSnp=UpToRightSnp
        RL=2                    ! Set UpToSnp at the end of the core in the next step

    else
        UpToSnp=UpToLeftSnp     ! Going backward (from right to left)
        RL=1                    ! Set UpToSnp at the start of the core in the next step
    endif

    h=0
    do ! Repeat till all SNPs have been covered 
        if ((nCoreA==1).and.(nCoreB==1)) exit   ! If the number of cores is 1, EXIT
                                                ! This will force the subroutine to finish 
                                                ! since it will be exit from both DO statements
        h=h+1
        if (mod(h,2)/=0) then                   ! If ODD
            do g=1,nCoreB
                if ((CoreIndexB(g,1)<UptoSnp).and.(CoreIndexB(g,2)>UptoSnp)) then
                    UpToCoreB=g
                    exit
                endif
            enddo
            if (e==1) then
                StartSnp=CoreIndexB(UpToCoreB,1)
                EndSnp=CoreIndexB(UpToCoreB,2)
                StPt=StartSnp
                EndPt=UpToSnp
                FillInSt=StartSnp
                FillInEnd=EndSnp
            else
                EndSnp=CoreIndexB(UpToCoreB,1)
                StartSnp=CoreIndexB(UpToCoreB,2)
                StPt=UpToSnp
                EndPt=StartSnp
                FillInSt=EndSnp
                FillInEnd=StartSnp
            endif
            
            CompLength=abs(UpToSnp-StartSnp)+1

            do i=1,nAnisP
                ! if ((BaseAnimal(i)==1).and.(PosHD(i)/=0).and.((i,RL)==0)) then
                if ((RecPed(i,2)==0).and.(RecPed(i,3)==0).and.(PosHD(i)/=0).and.(AnimRecomb(i,RL)==0)) then
                    C1=0
                    C2=0
                    C3=0
                    C4=0
                    Recmb=1
                    do j=StPt,EndPt
                        ! NOTE: ImputePhase array is a HD phase data because the SNPs are within the MiddleCoreA core
                        if (PhaseHD(PosHD(i),j,1,2)==ImputePhase(i,j,1)) C1=C1+1
                        if (PhaseHD(PosHD(i),j,1,2)==ImputePhase(i,j,2)) C2=C2+1
                        if (PhaseHD(PosHD(i),j,2,2)==ImputePhase(i,j,1)) C3=C3+1
                        if (PhaseHD(PosHD(i),j,2,2)==ImputePhase(i,j,2)) C4=C4+1        
                    enddo
                    if ((CompLength==C1).and.(CompLength/=C3)) then     ! If one haplotype is the same as the paternal, impute
                        ImputePhase(i,FillInSt:FillInEnd,1)=PhaseHD(PosHD(i),FillInSt:FillInEnd,1,2)
                        Recmb=0
                    endif
                    if ((CompLength/=C1).and.(CompLength==C3)) then     ! If one haplotype is the same as the paternal, impute
                        ImputePhase(i,FillInSt:FillInEnd,1)=PhaseHD(PosHD(i),FillInSt:FillInEnd,2,2)
                        Recmb=0
                    endif
                    if ((CompLength==C2).and.(CompLength/=C4)) then     ! If one haplotype is the same as the maternal, impute
                        ImputePhase(i,FillInSt:FillInEnd,2)=PhaseHD(PosHD(i),FillInSt:FillInEnd,1,2)
                        Recmb=0
                    endif
                    if ((CompLength/=C2).and.(CompLength==C4)) then     ! If one haplotype is the same as the maternal, impute
                        ImputePhase(i,FillInSt:FillInEnd,2)=PhaseHD(PosHD(i),FillInSt:FillInEnd,2,2)
                        Recmb=0
                    endif
                    if (Recmb==1) AnimRecomb(i,RL)=1    ! There is bridge in the phase in the RL direction. Nothing more will be done
                endif
            enddo
            UpToSnp=CoreIndexB(UpToCoreB,RL)
        else                                    ! if EVEN
            do g=1,nCoreA
                if ((CoreIndexA(g,1)<UptoSnp).and.(CoreIndexA(g,2)>UptoSnp)) then
                    UpToCoreA=g
                    exit
                endif
            enddo
            if (e==1) then
                StartSnp=CoreIndexA(UpToCoreA,1)
                EndSnp=CoreIndexA(UpToCoreA,2)
                StPt=StartSnp
                EndPt=UpToSnp
                FillInSt=StartSnp
                FillInEnd=EndSnp
            else
                EndSnp=CoreIndexA(UpToCoreA,1)
                StartSnp=CoreIndexA(UpToCoreA,2)
                StPt=UpToSnp
                EndPt=StartSnp
                FillInSt=EndSnp
                FillInEnd=StartSnp
            endif
            CompLength=abs(UpToSnp-StartSnp)+1
            do i=1,nAnisP
                if ((RecPed(i,2)==0).and.(RecPed(i,3)==0).and.(PosHD(i)/=0).and.(AnimRecomb(i,RL)==0)) then
                    C1=0
                    C2=0
                    C3=0
                    C4=0
                    Recmb=1
                    do j=StPt,EndPt
                        if (PhaseHD(PosHD(i),j,1,1)==ImputePhase(i,j,1)) C1=C1+1
                        if (PhaseHD(PosHD(i),j,1,1)==ImputePhase(i,j,2)) C2=C2+1
                        if (PhaseHD(PosHD(i),j,2,1)==ImputePhase(i,j,1)) C3=C3+1
                        if (PhaseHD(PosHD(i),j,2,1)==ImputePhase(i,j,2)) C4=C4+1
                    enddo
                    if ((CompLength==C1).and.(CompLength/=C3)) then
                        ImputePhase(i,FillInSt:FillInEnd,1)=PhaseHD(PosHD(i),FillInSt:FillInEnd,1,1)
                        Recmb=0
                    endif
                    if ((CompLength/=C1).and.(CompLength==C3)) then
                        ImputePhase(i,FillInSt:FillInEnd,1)=PhaseHD(PosHD(i),FillInSt:FillInEnd,2,1)
                        Recmb=0
                    endif
                    if ((CompLength==C2).and.(CompLength/=C4)) then
                        ImputePhase(i,FillInSt:FillInEnd,2)=PhaseHD(PosHD(i),FillInSt:FillInEnd,1,1)
                        Recmb=0
                    endif
                    if ((CompLength/=C2).and.(CompLength==C4)) then
                        ImputePhase(i,FillInSt:FillInEnd,2)=PhaseHD(PosHD(i),FillInSt:FillInEnd,2,1)
                        Recmb=0
                    endif
                    if (Recmb==1) (i,RL)=1
                endif
            enddo
            UpToSnp=CoreIndexA(UpToCoreA,RL)
        endif

        ! Exit condition
        if (e==1) then
            if (UpToSnp>=nSnp) exit     ! Exit if we've reached the last SNP
        else
            if (UpToSnp<=1) exit        ! Exit if we've reached the first SNP
        endif
    enddo
enddo

ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

end subroutine BaseAnimalFillIn

!#############################################################################################################################################################################################################################

subroutine InitialiseArrays
! Impute phase information for homozygous cases
use Global
implicit none

integer :: i,j,dum

ImputeGenos=9
ImputePhase=9

! Get information from RecodedGeneProbInput.txt which has been created in Makefiles subroutine
! WARNING: Why don't read information from Geno(:,:) that has been used to feed RecodedGeneProbInput.txt instead??
!          Read from file is always slower!
open (unit=43,file='./InputFiles/RecodedGeneProbInput.txt',status='old')
do i=1,nAnisP
    read (43,*) dum,dum,dum,ImputeGenos(i,:)
    do j=1,nSnp
        if (ImputeGenos(i,j)==0) ImputePhase(i,j,:)=0
        if (ImputeGenos(i,j)==2) ImputePhase(i,j,:)=1
    enddo
enddo

! WARNING: Close statement is missing?

end subroutine InitialiseArrays

!#############################################################################################################################################################################################################################

subroutine GeneralFillIn
! This function implements the four Minor sub-steps explained in Hickey et al. (2012; Appendix A)
use Global
implicit none

call ParentHomoFill                     ! Minor sub-step 1. Parent Homozygous fill in
call PhaseComplement                    ! Minor sub-step 2. Phase Complement
call ImputeParentByProgenyComplement    ! Minor sub-step 3. Impute Parents from Progeny Complement
call MakeGenotype                       ! Minor sub-step 4. Make Genotype
if (TestVersion==1) call CurrentYield 
if (TestVersion==1) call Checker

ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

end subroutine GeneralFillIn

!#############################################################################################################################################################################################################################

subroutine CurrentYield
use Global
implicit none

integer :: i,CountPatAl,CountMatAl,CountGeno
real :: PropPatAl,PropMatAl,PropGeno,NotKnownStart,NotKnownEnd


CountPatAl=count(ImputePhase(:,:,1)==9)
CountMatAl=count(ImputePhase(:,:,2)==9)
CountGeno=count(ImputeGenos(:,:)/=9)

PropPatAl=100*(float(CountPatAl)/(nAnisP*nSnp))
PropMatAl=100*(float(CountMatAl)/(nAnisP*nSnp))

NotKnownStart=(nAnisP*nSnp)-CountRawGenos
NotKnownEnd=(nAnisP*nSnp)-CountGeno
PropGeno=100*((NotKnownStart-NotKnownEnd)/NotKnownStart)

print*, " "
print*, "           ","Proportion not imputed:"
write(*,'(a10,1x,a15,1x,a15,1x,a33)') "           ","Paternal allele","Maternal allele","Proportion missing now genotyped"
write (*,'(a10,1x,2f15.2,1x,f33.2)') "           ",PropPatAl,PropMatAl,PropGeno

end subroutine CurrentYield

!#############################################################################################################################################################################################################################

subroutine EnsureHetGametic
! Impute phase to Y chromosome from X chromosome for heterogametic individuals 
use Global
implicit none

integer :: i,j

do i=1,nAnisP
    if (RecGender(i)==HetGameticStatus) then
        do j=1,nSnp
            if ((ImputePhase(i,j,1)==9).and.(ImputePhase(i,j,2)/=9)) ImputePhase(i,j,1)=ImputePhase(i,j,2)
            if ((ImputePhase(i,j,2)==9).and.(ImputePhase(i,j,1)/=9)) ImputePhase(i,j,2)=ImputePhase(i,j,1)      
        enddo
    endif
enddo

end subroutine EnsureHetGametic

!#############################################################################################################################################################################################################################

subroutine MakeGenotype
! Any individual that has a missing genotype information but has both alleles
! known, has its genotype filled in as the sum of the two alleles
use Global
implicit none

integer :: i,j

do i=1,nAnisP
    do j=1,nSnp
        if (ImputeGenos(i,j)==9) then
            if ((ImputePhase(i,j,1)/=9).and.(ImputePhase(i,j,2)/=9)) ImputeGenos(i,j)=sum(ImputePhase(i,j,:))
        endif
    enddo
enddo

ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

end subroutine MakeGenotype

!#############################################################################################################################################################################################################################

subroutine PhaseComplement
! If the genotype at a locus for an individual is known and one of its alleles has been determined
! then impute the missing allele as the complement of the genotype and the known phased allele
use Global
implicit none

integer :: i,j

do i=1,nAnisP
    do j=1,nSnp
        if (ImputeGenos(i,j)/=9) then
            if ((ImputePhase(i,j,1)/=9).and.(ImputePhase(i,j,2)==9)) ImputePhase(i,j,2)=ImputeGenos(i,j)-ImputePhase(i,j,1)
            if ((ImputePhase(i,j,2)/=9).and.(ImputePhase(i,j,1)==9)) ImputePhase(i,j,1)=ImputeGenos(i,j)-ImputePhase(i,j,2)
        endif
    enddo
enddo

ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

end subroutine PhaseComplement

!#############################################################################################################################################################################################################################

subroutine ParentHomoFill
! Fill in the allele of an offspring of a parent that has both its 
! alleles filled in and has a resulting genotype that is homozygous
use Global
implicit none

integer :: e,i,j,PatMat,ParId

do i=1,nAnisP
    ! WARNING: (SexOpt==0).or.((SexOpt==1) is always TRUE
    if ((SexOpt==0).or.((SexOpt==1).and.(RecGender(i)/=HetGameticStatus))) then     ! If individual is homogametic
        do e=1,2
            ParId=RecPed(i,e+1)
            do j=1,nSnp
                if (ImputePhase(i,j,e)==9) then                                     ! Always that the SNP is not genotyped
                    if ((ImputePhase(ParId,j,1)==ImputePhase(ParId,j,2)).and.(ImputePhase(ParId,j,1)/=9)) &
                        ImputePhase(i,j,e)=ImputePhase(ParId,j,1)                   ! Imput phase if parent is homozygous
                endif
            enddo
        enddo
    else
        ParId=RecPed(i,HomGameticStatus+1)      ! The homogametic parent
        do j=1,nSnp
            if (ImputePhase(i,j,1)==9) then     !Comment from John Hickey see analogous iterate subroutine
                if ((ImputePhase(ParId,j,1)==ImputePhase(ParId,j,2)).and.(ImputePhase(ParId,j,1)/=9)) &
                    ImputePhase(i,j,:)=ImputePhase(ParId,j,1)                      ! Imput phase to the two haplotypes if parent is homozygous 
            endif
        enddo
    endif
enddo
ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

end subroutine ParentHomoFill

!#############################################################################################################################################################################################################################

subroutine ImputeParentByProgenyComplement
! If one of the parental allele is known and the other missing, the fill
! in the missing allele in the parent if, at least, one of its offspring
! is known to carry an allele that does not match the known allele in the parent    
use Global

implicit none

integer :: i,j,k,l,Count1,Count0

do i=1,nAnisP
    do j=1,2
        if (SireDam(i,j)==1) then       ! We are only interested in sires because they have more progeny

            ! Sex chromosome
            if (SexOpt==1) then
                do k=1,nSnp

                    ! Mat gamete missing -> fill if offspring suggest heterozygous
                    ! WARNING: This was comment the other way around in the original version of the code
                    if ((ImputePhase(i,k,1)/=9).and.(ImputePhase(i,k,2)==9)) then
                        Count1=0
                        Count0=0
                        if (ImputePhase(i,k,1)==1) Count1=1
                        if (ImputePhase(i,k,1)==0) Count0=1

                        ! Look for the individual progeny and count their phase
                        do l=1,nAnisP

                            ! WARNING: This is the only difference with the SexOpt=0 code below. Duplicating
                            !           the code can be avoided by including the IF statement here instead than
                            !           outside the SNPs loop.
                            if ((RecGender(i)==HetGameticStatus).and.(RecGender(l)==HetGameticStatus)) cycle

                            if (RecPed(l,j+1)==i) then
                                if (ImputePhase(l,k,j)==0) Count0=Count0+1
                                if (ImputePhase(l,k,j)==1) Count1=Count1+1
                            endif   
                            if ((Count0>0).and.(Count1>0)) then                 !Consider increasing the number of offspring required to all ow for genotyping error
                                if (ImputePhase(i,k,1)==0) ImputePhase(i,k,2)=1
                                if (ImputePhase(i,k,1)==1) ImputePhase(i,k,2)=0
                                exit
                            endif
                        enddo
                    endif
    
                    !Pat gamete missing -> fill if offspring suggest heterozygous
                    ! WARNING: This comment was the other way around in the original version of the code
                    if ((ImputePhase(i,k,2)/=9).and.(ImputePhase(i,k,1)==9)) then  
                        Count1=0
                        Count0=0
                        if (ImputePhase(i,k,2)==1) Count1=1
                        if (ImputePhase(i,k,2)==0) Count0=1
                        do l=1,nAnisP

                            ! WARNING: This is the only difference with the SexOpt=0 code below. Duplicating
                            !           the code can be avoided by including the IF statement here instead than
                            !           outside the SNPs loop.
                            if ((RecGender(i)==HetGameticStatus).and.(RecGender(l)==HetGameticStatus)) cycle

                            if (RecPed(l,j+1)==i) then
                                if (ImputePhase(l,k,j)==0) Count0=Count0+1
                                if (ImputePhase(l,k,j)==1) Count1=Count1+1
                            endif   
                            if ((Count0>0).and.(Count1>0)) then                 !Consider increasing the number of offspring required to all ow for genotyping error
                                if (ImputePhase(i,k,2)==0) ImputePhase(i,k,1)=1
                                if (ImputePhase(i,k,2)==1) ImputePhase(i,k,1)=0
                                exit
                            endif
                        enddo
                    endif
                enddo

            ! Generic chromosome
            else
                do k=1,nSnp
                    if ((ImputePhase(i,k,1)/=9).and.(ImputePhase(i,k,2)==9)) then               !Pat gamete missing fill if offspring suggest heterozygous
                        Count1=0
                        Count0=0
                        if (ImputePhase(i,k,1)==1) Count1=1
                        if (ImputePhase(i,k,1)==0) Count0=1
                        do l=1,nAnisP
                            if (RecPed(l,j+1)==i) then
                                if (ImputePhase(l,k,j)==0) Count0=Count0+1
                                if (ImputePhase(l,k,j)==1) Count1=Count1+1
                            endif   
                            if ((Count0>0).and.(Count1>0)) then                 !Consider increasing the number of offspring required to all ow for genotyping error
                                if (ImputePhase(i,k,1)==0) ImputePhase(i,k,2)=1
                                if (ImputePhase(i,k,1)==1) ImputePhase(i,k,2)=0
                                exit
                            endif
                        enddo
                    endif
    
                    if ((ImputePhase(i,k,2)/=9).and.(ImputePhase(i,k,1)==9)) then               !Mat gamete missing fill if offspring suggest heterozygous  
                        Count1=0
                        Count0=0
                        if (ImputePhase(i,k,2)==1) Count1=1
                        if (ImputePhase(i,k,2)==0) Count0=1
                        do l=1,nAnisP
                            if (RecPed(l,j+1)==i) then
                                if (ImputePhase(l,k,j)==0) Count0=Count0+1
                                if (ImputePhase(l,k,j)==1) Count1=Count1+1
                            endif   
                            if ((Count0>0).and.(Count1>0)) then                 !Consider increasing the number of offspring required to all ow for genotyping error
                                if (ImputePhase(i,k,2)==0) ImputePhase(i,k,1)=1
                                if (ImputePhase(i,k,2)==1) ImputePhase(i,k,1)=0
                                exit
                            endif
                        enddo
                    endif
                enddo
            endif   
        endif   
    enddo
enddo

ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

end subroutine ImputeParentByProgenyComplement

!#############################################################################################################################################################################################################################

subroutine GeneProbPhase
!! Provides information about:
!!   * if the paternal and maternal haplotypes have been phased
!!   * how many and which alleles of my parents have been phased
!!   * how many and which SNPs of my grandparents are heterozygous
!! and stores all this information in two files:
!!   * IndividualSnpInformativeness.txt
!!   * IndividualMendelianInformativeness.txt

use Global
use GlobalPedigree
implicit none

integer :: h,i,j,k,m,dum,StSnp,EnSnp
real :: PatAlleleProb(nSnp,2),MatAlleleProb(nSnp,2),HetProb(nSnp),GeneProbWork(nSnp,4)
integer :: Informativeness(nSnp,6),TmpInfor(nSnp,6),GrandPar
character(len=300) :: filout

if (BypassGeneProb==0) then
    ! Get information from GeneProb
    do h=1,nProcessors
        write (filout,'("./GeneProb/GeneProb"i0,"/GeneProbs.txt")')h            !here
        open (unit=110,file=trim(filout),status="unknown")
        StSnp=GpIndex(h,1)
        EnSnp=GpIndex(h,2)
        do i=1,nAnisP
            do j=1,4
                read (110,*) dum,GeneProbWork(StSnp:EnSnp,j)
            enddo
            PatAlleleProb(StSnp:EnSnp,1)=GeneProbWork(StSnp:EnSnp,1)+GeneProbWork(StSnp:EnSnp,2)
            PatAlleleProb(StSnp:EnSnp,2)=GeneProbWork(StSnp:EnSnp,3)+GeneProbWork(StSnp:EnSnp,4)
            MatAlleleProb(StSnp:EnSnp,1)=GeneProbWork(StSnp:EnSnp,1)+GeneProbWork(StSnp:EnSnp,3)
            MatAlleleProb(StSnp:EnSnp,2)=GeneProbWork(StSnp:EnSnp,2)+GeneProbWork(StSnp:EnSnp,4)
    
            ! Probability of heterozygosity
            HetProb(StSnp:EnSnp)=GeneProbWork(StSnp:EnSnp,2)+GeneProbWork(StSnp:EnSnp,3)
    
            do j=StSnp,EnSnp
                if (PatAlleleProb(j,1)>=GeneProbThresh) ImputePhase(i,j,1)=0
                if (PatAlleleProb(j,2)>=GeneProbThresh) ImputePhase(i,j,1)=1
                if (MatAlleleProb(j,1)>=GeneProbThresh) ImputePhase(i,j,2)=0
                if (MatAlleleProb(j,2)>=GeneProbThresh) ImputePhase(i,j,2)=1
                if (HetProb(j)>=GeneProbThresh) ImputeGenos(i,j)=1
            enddo
        enddo
        close(110)
    enddo
endif
    
ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

if (WindowsLinux==1) then
        open (unit=102,file=".\Miscellaneous\IndividualSnpInformativeness.txt",status="unknown")
else
        open (unit=102,file="./Miscellaneous/IndividualSnpInformativeness.txt",status="unknown")
endif

if (WindowsLinux==1) then
        open (unit=103,file=".\Miscellaneous\IndividualMendelianInformativeness.txt",status="unknown")
else
        open (unit=103,file="./Miscellaneous/IndividualMendelianInformativeness.txt",status="unknown")
endif

allocate(GlobalTmpCountInf(nAnisP,8))
allocate(MSTermInfo(nAnisP,2))

MSTermInfo=0
do i=1,nAnisP
    if (IndivIsGenotyped(i)==1) MSTermInfo(i,:)=1
    TmpInfor(:,:)=-99
    GlobalTmpCountInf(i,:)=0
    Informativeness(:,:)=9 ! What the hell is this variable for??
    j=0
    ! Check whether my parents and grandparents are heterozygous
    do m=1,nSnpRaw
        if (SnpIncluded(m)==1) then                     ! Whether to consider this SNP
            j=j+1                                       ! Number of SNPs included so far
            if (ImputeGenos(i,j)==1) then               ! If heterozygous

                ! My father is heterozygous
                if (ImputeGenos(RecPed(i,2),j)==1) then
                    ! And have my father haplotype phased
                    if ((ImputePhase(i,j,1)==0).or.(ImputePhase(i,j,1)==1)) then
                        Informativeness(j,1)=1
                        GlobalTmpCountInf(i,1)=GlobalTmpCountInf(i,1)+1     ! Count the number of SNPs phased of my father haplotype
                        TmpInfor(GlobalTmpCountInf(i,1),1)=j                ! The SNP no. GlobalTmpCountInf(i,1) is my SNP no. j
                    endif               
                endif

                ! My mother is heterozygous
                if (ImputeGenos(RecPed(i,3),j)==1) then
                    ! And have my mother haplotype phased
                    if ((ImputePhase(i,j,2)==0).or.(ImputePhase(i,j,2)==1)) then
                        Informativeness(j,2)=1
                        GlobalTmpCountInf(i,2)=GlobalTmpCountInf(i,2)+1
                        TmpInfor(GlobalTmpCountInf(i,2),2)=j                
                    endif
                endif   
    
                ! My father haplotype is phased
                if ((ImputePhase(i,j,1)==0).or.(ImputePhase(i,j,1)==1)) then
                    ! If my paternal GranSire is heterozygous
                    GrandPar=RecPed(RecPed(i,2),2)
                    if (ImputeGenos(GrandPar,j)==1) then
                        Informativeness(j,3)=1
                        GlobalTmpCountInf(i,3)=GlobalTmpCountInf(i,3)+1
                        TmpInfor(GlobalTmpCountInf(i,3),3)=j                
                    endif
                    ! If my maternal GranDam is heterozygous
                    GrandPar=RecPed(RecPed(i,2),3) 
                    if (ImputeGenos(GrandPar,j)==1) then
                        Informativeness(j,4)=1
                        GlobalTmpCountInf(i,4)=GlobalTmpCountInf(i,4)+1
                        TmpInfor(GlobalTmpCountInf(i,4),4)=j                
                    endif
                endif
    
                ! My mother haplotype is phased
                if ((ImputePhase(i,j,2)==0).or.(ImputePhase(i,j,2)==1)) then
                    ! If my maternal GranSire is heterozygous
                    GrandPar=RecPed(RecPed(i,3),2) 
                    if (ImputeGenos(GrandPar,j)==1) then
                        Informativeness(j,5)=1
                        GlobalTmpCountInf(i,5)=GlobalTmpCountInf(i,5)+1
                        TmpInfor(GlobalTmpCountInf(i,5),5)=j                
                    endif
                    ! If my maternal GranDam is heterozygous
                    GrandPar=RecPed(RecPed(i,3),3) 
                    if (ImputeGenos(GrandPar,j)==1) then
                        Informativeness(j,6)=1
                        GlobalTmpCountInf(i,6)=GlobalTmpCountInf(i,6)+1
                        TmpInfor(GlobalTmpCountInf(i,6),6)=j                
                    endif
                endif
            endif   
        endif   
    enddo

    GlobalTmpCountInf(i,7)=count(ImputePhase(i,:,1)/=9)         ! Count the number of genotyped allele in the paternal haplotype
    GlobalTmpCountInf(i,8)=count(ImputePhase(i,:,2)/=9)         ! Count the number of genotyped allele in the maternal haplotype
    if (GlobalTmpCountInf(i,1)>0) MSTermInfo(i,1)=1             ! If the paternal haplotype is phased
    if (GlobalTmpCountInf(i,2)>0) MSTermInfo(i,2)=1             ! If the maternal haplotype is phased
    do k=1,6
        write(102,'(a20,i2,i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7)') Id(i), Setter(i), GlobalTmpCountInf(i,k), TmpInfor(1:GlobalTmpCountInf(i,k),k)
    enddo   
    write(102,'(a20,i2,i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7)') Id(i), Setter(i), GlobalTmpCountInf(i,7)
    write(102,'(a20,i2,i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7)') Id(i), Setter(i), GlobalTmpCountInf(i,8)    
    write(103,'(a20,2i10)') Id(i), MSTermInfo(i,:)
enddo
close(102)
close(103)

end subroutine GeneProbPhase

!######################################################################################################################################################################################

subroutine MakeFiles
use Global
use GlobalPedigree

implicit none

integer :: i,j,TempCore(nPhaseInternal),TempCplusT(nPhaseInternal)
integer :: Tmp
character(len=7) :: cm 
character(len=300) :: filout,FileCheck
logical :: FileExists

allocate(GpIndex(nProcessors,2))

! WARNING: This code is not necessary
write(cm,'(I7)') nSnpRaw !for formatting   
cm = adjustl(cm)                          
! end WARNING

do i=1,nPhaseExternal
    TempCore(i)=CoreLengths(i)
    TempCore(i+nPhaseExternal)=CoreLengths(i)
    TempCplusT(i)=CoreAndTailLengths(i)
    TempCplusT(i+nPhaseExternal)=CoreAndTailLengths(i)
    if (TempCore(i)>nSnp) TempCore(i)=nSnp 
    if (TempCore(i+nPhaseExternal)>nSnp) TempCore(i+nPhaseExternal)=nSnp
    if (TempCplusT(i)>nSnp) TempCplusT(i)=nSnp
    if (TempCplusT(i+nPhaseExternal)>nSnp) TempCplusT(i+nPhaseExternal)=nSnp
enddo

if (WindowsLinux==1) then
     open (unit=103,file=".\InputFiles\AlphaPhaseInputPedigree.txt",status="unknown")
     open (unit=104,file=".\InputFiles\RecodedGeneProbInput.txt",status="unknown")
     open (unit=105,file=".\InputFiles\AlphaPhaseInputGenotypes.txt",status="unknown")
else
     open (unit=103,file="./InputFiles/AlphaPhaseInputPedigree.txt",status="unknown") 
     open (unit=104,file="./InputFiles/RecodedGeneProbInput.txt",status="unknown")
     open (unit=105,file="./InputFiles/AlphaPhaseInputGenotypes.txt",status="unknown")
endif

do i=1,nAnisP
     write (104,'(i16,1x,i16,1x,i16,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') RecPed(i,:),Genos(i,:)        
     if (Setter(i)==1) write (105,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') Id(i),Genos(i,:)           
enddo
call flush(104)
call flush(105)
close(104)
close(105)

CountRawGenos=count(Genos(:,:)/=9)

do i=1,nObsDataRaw
    write (103,'(4a20)') Ped(i,:)
enddo
call flush(103)
close(103)


!ANDREAS AND JOHN CHANGE ON FRIDAY - POSSIBLY REMOVE
if (SexOpt==1) then                 ! Sex chromosome
    call InsteadOfGeneProb          ! Calculate Genotype Probabilities
else                                ! Not sex chromosome
    if (BypassGeneProb==1) then     ! Do I have to bypass Genotype Probabilities?
        call InsteadOfGeneProb
    else
        deallocate(Genos)
    endif
endif

! Check whether AlphaPhase is present
write (FileCheck,'("AlphaPhase1.1")')
inquire(file=trim(FileCheck),exist=FileExists)      
if (FileExists .eqv. .true.) then
    AlphaPhasePresent=1
else
    AlphaPhasePresent=0
    print*, " "
    print*, " ","AlphaPhase1.1 not present and copied, version in the bin directory used"
    ! WARNING: What if there is any version of this software in the bin either?
endif

! Check whether GeneProbForAlphaImpute is present
write (FileCheck,'("GeneProbForAlphaImpute")')
inquire(file=trim(FileCheck),exist=FileExists)      
if (FileExists .eqv. .true.) then
    GeneProbPresent=1
else
    GeneProbPresent=0
    print*, " "
    print*, " ","GeneProb not present and copied, version in the bin directory used"
    ! WARNING: What if there is any version of this software in the bin either?
endif

! Create AlphaPhaseSpec file
do i=1,nPhaseInternal           ! Phasing is done in parallel
    if (WindowsLinux==1) then
        ! WARNING: Apparently, AlphaImpute does not work for Windows systems
    else
        write (filout,'("./Phasing/Phase"i0,"/AlphaPhaseSpec.txt")')i
        open (unit=106,file=trim(filout),status='unknown')
        if (PedFreePhasing==0) then
            if (SexOpt==0) write (106,*) 'PedigreeFile              ,"../../InputFiles/AlphaPhaseInputPedigree.txt"'
        else
            if (SexOpt==0) write (106,*) 'PedigreeFile                      ,"NoPedigree"'
        endif
        if (SexOpt==1) write (106,*) 'PedigreeFile                      ,"NoPedigree"'

        write (106,'(a97)') &
                'GenotypeFile                       ,"../../InputFiles/AlphaPhaseInputGenotypes.txt",GenotypeFormat'
        write (106,*) 'NumberOfSnp                      ,',nSnp
        write (106,*) 'GeneralCoreAndTailLength     ,',TempCplusT(i)
        if(i<=nPhaseInternal/2) then  
            write (106,*) 'GeneralCoreLength        ,',TempCore(i),',Offset'
        else
            write (106,*) 'GeneralCoreLength        ,',TempCore(i),',NotOffset'
        endif
        write (106,*) 'UseThisNumberOfSurrogates    ,',10
        write (106,*) 'PercentageSurrDisagree       ,',10.00
        write (106,*) 'PercentageGenoHaploDisagree  ,',GenotypeErrorPhase
        write (106,*) 'GenotypeMissingErrorPercentage   ,',0.00
        write (106,*) 'NrmThresh                        ,',0.00
        write (106,*) 'FullOutput                   ,0'
        write (106,*) 'Graphics             ,0'
        write (106,*) 'Simulation           ,0'
        write (106,*) 'TruePhaseFile            ,None'
        call flush(106)
        close(106)
        write (filout,'("Phase"i0)')i
        if (AlphaPhasePresent==1) call system ("cp AlphaPhase1.1 Phasing/" // filout) 
    endif
enddo

Tmp=int(float(nSnp)/nProcessors)
GpIndex(1,1)=1
GpIndex(1,2)=Tmp
if (nProcessors>1) then
    do i=2,nProcessors
        GpIndex(i,1)=GpIndex(i-1,1)+Tmp
        GpIndex(i,2)=GpIndex(i-1,2)+Tmp
    enddo
endif
GpIndex(nProcessors,2)=nSnp

! Create GeneProbSpec file
do i=1,nProcessors
    write (filout,'("./GeneProb/GeneProb"i0,"/GeneProbSpec.txt")')i

    open (unit=108,file=trim(filout),status='unknown')
    write (108,*) "nAnis        ,",nAnisP 
    write (108,*) "nSnp     ,",nSnp
    write (108,*) "InputFilePath    ,",'"../../InputFiles/RecodedGeneProbInput.txt"'
    write (108,*) "OutputFilePath   ,",'"GeneProbs.txt"'
    write (108,*) "StartSnp     ,",GpIndex(i,1) 
    write (108,*) "EndSnp       ,",GpIndex(i,2) 
    call flush(108)
    close(108)

    write (filout,'("GeneProb"i0)')i
    if (GeneProbPresent==1) call system ("cp GeneProbForAlphaImpute GeneProb/" // filout)
enddo

if (PreProcess==1) then
    print*, "  "
    print*, "  ","The program has preprocessed the data and now it stops"
    stop    
endif

end subroutine MakeFiles

!#############################################################################################################################################################################################################################

subroutine InternalEdit
use Global
use GlobalPedigree

implicit none

integer :: i,j,k,CountMiss,CountHD,nSnpR,dum,Counter(nSnp)
real :: SnpSummary(nSnp),TempFreq(nSnp)
character (len=300) :: dumC,FileName

allocate(SnpIncluded(nSnp))
allocate(Setter(0:nAnisP))

SnpIncluded(:)=0
if ((ManagePhaseOn1Off0==0).and.(NoPhasing==1)) then
    write (FileName,'(a,"/EditingSnpSummary.txt")') trim(PhasePath)
    open(unit=111,file=trim(FileName),status="old")
    do i=1,nSnp
        read (111,*) dum,SnpSummary(i),SnpIncluded(i)
    enddo
endif

if (NoPhasing==0) SnpIncluded(:)=1

if (UserDefinedHD==0) then
    Setter(0)=0
    Setter(1:nAnisP)=1
    RecIdHDIndex(0)=0
    RecIdHDIndex(1:nAnisP)=1
    do i=1,nAnisP
        CountMiss=count(TempGenos(i,:)==9)
        if ((float(CountMiss)/nSnp)>(1.0-PercGenoForHD)) then
            Setter(i)=0
            RecIdHDIndex(i)=0
        endif
    enddo
    CountHD=count(Setter(:)==1)
else
    Setter(0)=0
    Setter(1:nAnisP)=0
    RecIdHDIndex(0)=0
    RecIdHDIndex(1:nAnisP)=0

    CountHD=0
    do
        read (46,*,iostat=k) dumC
        CountHD=CountHD+1
        if (k/=0) then
            CountHD=CountHD-1
            exit
        endif
    enddo
    rewind(46)
    do k=1,CountHD
        read (46,*) dumC
        do i=1,nAnisP
            if (trim(dumC)==trim(Id(i))) then
                Setter(i)=1
                RecIdHDIndex(i)=1
                exit
            endif
        enddo   
    enddo
    CountHD=count(Setter(:)==1)
    print*, " "
    print*, " ",CountHD," valid indiviudals in the user specified AlphaPhase1.1 file"
endif


if (WindowsLinux==1) then
        open (unit=102,file=".\Miscellaneous\EditingSnpSummary.txt",status="unknown")
else
        open (unit=102,file="./Miscellaneous/EditingSnpSummary.txt",status="unknown")
endif

if (ManagePhaseOn1Off0==1) then
    TempFreq(:)=0.0
    Counter(:)=0
    do i=1,nAnisP
        do j=1,nSnp
            if (TempGenos(i,j)/=9) then
                TempFreq(j)=TempFreq(j)+float(TempGenos(i,j))
                Counter(j)=Counter(j)+2
            endif   
        enddo
    enddo
    do j=1,nSnp
        if (Counter(j)>0.000000) TempFreq(j)=TempFreq(j)/Counter(j)
    enddo
endif

if (ManagePhaseOn1Off0==1) then
    SnpSummary=0.0
    do i=1,nAnisP
        if (Setter(i)==1) then
            do j=1,nSnp
                if (TempGenos(i,j)==9) SnpSummary(j)=SnpSummary(j)+1.0
            enddo
        endif
    enddo   
endif

if (IntEditStat==0) then
    nSnpR=nSnp
    allocate(Genos(0:nAnisP,nSnp))
    Genos=TempGenos
    deallocate(TempGenos)
    if (ManagePhaseOn1Off0==1) SnpIncluded(:)=1
else
    if (ManagePhaseOn1Off0==1) then
        SnpSummary(:)=SnpSummary(:)/CountHD
        nSnpR=0
        do j=1,nSnp
            if ((SnpSummary(j)<PercSnpMiss).and.((TempFreq(j)>0.00000001).and.(TempFreq(j)<0.9999999))) nSnpR=nSnpR+1
        enddo       
    else
        nSnpR=count(SnpIncluded(:)==1)
    endif
    if (nSnpR==nSnp) then
        allocate(Genos(0:nAnisP,nSnp))
        Genos=TempGenos
        deallocate(TempGenos)
        SnpIncluded(:)=1
    else
        allocate(Genos(0:nAnisP,nSnpR))
        Genos(0,:)=9
        if (ManagePhaseOn1Off0==1) then
            k=0
            do j=1,nSnp
                if ((SnpSummary(j)<PercSnpMiss).and.((TempFreq(j)>0.00000001).and.(TempFreq(j)<0.9999999))) then
                    k=k+1
                    Genos(:,k)=TempGenos(:,j)
                    SnpIncluded(j)=1
                endif   
            enddo
            deallocate(TempGenos)
            nSnp=nSnpR  
        else
            k=0
            do j=1,nSnp
                if (SnpIncluded(j)==1) then
                    k=k+1
                    Genos(:,k)=TempGenos(:,j)
                endif   
            enddo
            deallocate(TempGenos)
            nSnp=nSnpR  
        endif
    endif
    if (UserDefinedHD==0) then
        Setter(1:nAnisP)=1
        RecIdHDIndex(1:nAnisP)=1
        do i=1,nAnisP
            CountMiss=count(Genos(i,:)==9)
            if ((float(CountMiss)/nSnp)>(1.0-SecondPercGenoForHD)) then
                Setter(i)=0
                RecIdHDIndex(i)=0
            endif
        enddo
        CountHD=count(Setter(:)==1)
    else
        do i=1,nAnisP
            if (Setter(i)==1) then
                CountMiss=count(Genos(i,:)==9)
                if ((float(CountMiss)/nSnp)>(1.0-SecondPercGenoForHD)) then
                    Setter(i)=0
                    RecIdHDIndex(i)=0
                endif
            endif
        enddo
        CountHD=count(Setter(:)==1)
    endif
endif

do j=1,nSnpRaw
    write (102,*) j,SnpSummary(j),SnpIncluded(j)        !'(i,1x,f5.3,1x,i)'
enddo
close(102)

open (unit=112,file="./Phasing/EditingSnpSummary.txt",status="unknown")
do j=1,nSnpRaw
    write (112,*) j,SnpSummary(j),SnpIncluded(j)        !'(i,1x,f5.3,1x,i)'
enddo
close(112)

print*, " "
print*, " "
print*, " ",CountHD," indiviudals passed to AlphaPhase1.1"
print*, " ",nSnp," snp remain after editing"

end subroutine InternalEdit

!#############################################################################################################################################################################################################################

subroutine FillInBasedOnOffspring
! Genotype SNPs based on the genetic information of my offsprings
use Global

implicit none

integer :: i,j,k,Count0(nSnp),Count1(nSnp),Count2(nSnp)

do i=1,nAnisP ! These are parents
    ! This three variables will count the different number of genotypes of the offsprings
    Count0=0
    Count1=0
    Count2=0
    do j=1,nAnisP ! These are offsprings
        if ((RecPed(j,2)==i).or.(RecPed(j,3)==i)) then
            if ((SexOpt==1).and.(RecGender(i)==HetGameticStatus).and.(RecGender(j)==HetGameticStatus)) cycle
            do k=1,nSnp
                if (TempGenos(i,k)==9) then ! If my parent is not genotyped
                    if (TempGenos(j,k)==0) Count0(k)=Count0(k)+1    ! Number of offspring genotype as 0
                    if (TempGenos(j,k)==1) Count1(k)=Count1(k)+1    ! Number of offspring genotype as 1
                    if (TempGenos(j,k)==2) Count2(k)=Count2(k)+1    ! Number of offspring genotype as 2
                endif
            enddo
        endif
    enddo

    do k=1,nSnp
        if ((Count0(k)+Count1(k)+Count2(k))>OffspringFillMin) then
            if (Count0(k)==0) TempGenos(i,k)=2                       ! This is the most likely thing, but it might be not true
            if (Count2(k)==0) TempGenos(i,k)=0                       ! This is the most likely thing, but it might be not true
            if ((SexOpt==1).and.(RecGender(i)==HetGameticStatus)) cycle
            if ((Count0(k)>2).and.(Count2(k)>2)) TempGenos(i,k)=1    ! This is the most likely thing, but it might be not true
        endif
    enddo
enddo

end subroutine FillInBasedOnOffspring

!#############################################################################################################################################################################################################################

subroutine FillInSnp
! Genotype SNPs based on the pedigree information
use Global

implicit none

integer :: i,j,k,TurnOn

do i=1,nAnisP
    do k=2,3
        TurnOn=1
        ! if the proband is heterogametic, and 
        ! considering the heterogametic parent, then avoid!!
        if ((SexOpt==1).and.(RecGender(i)==HetGameticStatus).and.((k-1)==HetGameticStatus)) TurnOn=0

        ! Homogametic individuals and the homogametic parent of a heterogametic individual
        if (TurnOn==1) then
            do j=1,nSnp
                if ((TempGenos(i,j)==0).and.(TempGenos(RecPed(i,k),j)==2)) then 
                    TempGenos(i,j)=9
                    TempGenos(RecPed(i,k),j)=9
                endif
                if ((TempGenos(i,j)==2).and.(TempGenos(RecPed(i,k),j)==0)) then
                    TempGenos(i,j)=9
                    TempGenos(RecPed(i,k),j)=9
                endif
            enddo
        endif
    enddo
enddo 

! WARNING: This can be refactored
do i=1,nAnisP
    do j=1,nSnp
        if (TempGenos(i,j)==9) then 
            if ((TempGenos(RecPed(i,2),j)==0).and.(TempGenos(RecPed(i,3),j)==0)) TempGenos(i,j)=0
            if ((TempGenos(RecPed(i,2),j)==2).and.(TempGenos(RecPed(i,3),j)==2)) TempGenos(i,j)=2
            if (SexOpt==1) then
                if (RecGender(i)/=HetGameticStatus) then
                    if ((TempGenos(RecPed(i,2),j)==0).and.(TempGenos(RecPed(i,3),j)==2)) TempGenos(i,j)=1  
                    if ((TempGenos(RecPed(i,2),j)==2).and.(TempGenos(RecPed(i,3),j)==0)) TempGenos(i,j)=1
                else
                    ! HomGameticSatus(1 or 2) +1 = sire (2) or dam (3)
                    if (TempGenos(RecPed(i,(HomGameticStatus+1)),j)==0) TempGenos(i,j)=0
                    if (TempGenos(RecPed(i,(HomGameticStatus+1)),j)==2) TempGenos(i,j)=2
                endif
            else
                if ((TempGenos(RecPed(i,2),j)==0).and.(TempGenos(RecPed(i,3),j)==2)) TempGenos(i,j)=1  
                if ((TempGenos(RecPed(i,2),j)==2).and.(TempGenos(RecPed(i,3),j)==0)) TempGenos(i,j)=1
            endif
        endif
    enddo   
enddo

end subroutine FillInSnp

!#############################################################################################################################################################################################################################

subroutine CheckParentage
! Check which is the parentage relation between animals. 
! Set the vector baseline that specify whether an animal 
! has no parents or its parents have been pruned
use GlobalPedigree

use Global

implicit none

integer :: e,i,j,k,CountBothGeno,CountDisagree,CountChanges,GenoYesNo(nAnisRawPedigree,3),IndId,ParId,flag,SumPruned,ParPos
integer :: TurnOn
integer,allocatable,dimension (:) :: Genotyped,Pruned
logical,allocatable,dimension (:) :: IsParent

if (WindowsLinux==1) then
    open (unit=101,file=".\Miscellaneous\PedigreeMistakes.txt",status="unknown")
else
    open (unit=101,file="./Miscellaneous/PedigreeMistakes.txt",status="unknown")
endif

GenoYesNo=0         ! Matrix (nAnisRawPedigree x 3). 
                    ! It basically says which is the proband's genotype (GenoYesNo(:,1)) but also
                    ! the genotype of proband's sire (GenoYesNo(:,2)) and dam (GenoYesNo(:,3))
do i=1,nAnisRawPedigree
    do j=1,nAnisG
        if (trim(GenotypeId(j))==trim(Ped(i,1))) then
            GenoYesNo(i,1)=j
            exit
        endif       
    enddo
    do j=1,nAnisG
        if (trim(GenotypeId(j))==trim(Ped(i,2))) then
            GenoYesNo(i,2)=j
            exit
        endif       
    enddo
    do j=1,nAnisG
        if (trim(GenotypeId(j))==trim(Ped(i,3))) then
            GenoYesNo(i,3)=j
            exit
        endif       
    enddo
enddo

CountChanges=0
do e=1,2                    ! Do whatever this does, first on males and then on females
    ParPos=e+1              ! Index in the Genotype and Pedigree matrices for sires and dams
    do i=1,nAnisRawPedigree
        IndId=GenoYesNo(i,1)            ! My Id
        ParId=GenoYesNo(i,ParPos)       ! Paternal Id,
        TurnOn=1
        ! GenderRaw: Says whether the proband is a male or a female
        ! If there are sex cromosome information, and
        ! if the proband is heterogametic, and 
        ! I am considering the heterogametic parent, then avoid!!
        ! That is, avoid males and their sires (hetero=1), or females and their dams (hetero=2) 
        if ((SexOpt==1).and.(GenderRaw(i)==HetGameticStatus).and.((ParPos-1)==HetGameticStatus)) TurnOn=0

        ! Considere the Homogametic probands and the heterogametic proband with homogametic parent
        if ((IndId/=0).and.(ParId/=0).and.(TurnOn==1)) then
            CountBothGeno=0                         
            CountDisagree=0

            ! Look for mendelenian errors
            do j=1,nSnp
                if ((Genos(IndId,j)/=9).and.(Genos(ParId,j)/=9)) then
                    CountBothGeno=CountBothGeno+1
                    if ((Genos(IndId,j)==0).and.(Genos(ParId,j)==2)) CountDisagree=CountDisagree+1
                    if ((Genos(IndId,j)==2).and.(Genos(ParId,j)==0)) CountDisagree=CountDisagree+1
                endif
            enddo
            if ((float(CountDisagree)/CountBothGeno)>DisagreeThreshold) then ! Mendelenian error
                Ped(i,ParPos)='0'
                write (101,*) Ped(i,1)
                CountChanges=CountChanges+1
            else
                ! Remove genotype of proband and parent
                do j=1,nSnp
                    if ((Genos(IndId,j)/=9).and.(Genos(ParId,j)/=9)) then
                        if ((Genos(IndId,j)==0).and.(Genos(ParId,j)==2)) then
                            Genos(IndId,j)=9
                            Genos(ParId,j)=9    
                        endif
                        if ((Genos(IndId,j)==2).and.(Genos(ParId,j)==0)) then
                            Genos(IndId,j)=9
                            Genos(ParId,j)=9    
                        endif
                    endif
                enddo
            endif
        endif 
    enddo
enddo
write (101,*) CountChanges," changes were made to the pedigree"
close (101)
print*, " ",CountChanges," errors in the pedigree due to Mendelian inconsistencies"

! Sort sires and dams, and look for mistakes (bisexuality,...).
call PVseq(nAnisRawPedigree,nAnisP)

allocate(RecPed(0:nAnisP,3))
allocate(RecIdHDIndex(0:nAnisP))
allocate(RecGender(0:nAnisP))
allocate(IndivIsGenotyped(nAnisP))

RecIdHDIndex=0
RecPed(0,:)=0
do i=1,nAnisP
    RecPed(i,1)=i
enddo
RecPed(1:nAnisP,2)=seqsire(1:nAnisP)
RecPed(1:nAnisP,3)=seqdam(1:nAnisP)

if (WindowsLinux==1) then
        open (unit=101,file=".\Miscellaneous\InternalDataRecoding.txt",status="unknown")
else
        open (unit=101,file="./Miscellaneous/InternalDataRecoding.txt",status="unknown")
endif

do i=1,nAnisP
    write (101,'(3i20,a20)') RecPed(i,:),trim(Id(i))
enddo
close (101)

RecGender=9     ! It basically says the gender of the individual
if (SexOpt==1) then
    do i=1,nAnisP
        TurnOn=1
        do j=1,nAnisInGenderFile
            if (trim(Id(i))==trim(GenderId(j))) then
                RecGender(i)=GenderRaw(j)
                TurnOn=0
                exit
            endif
        enddo   
        if (TurnOn==1) then
            do j=1,nAnisP
                if (i==RecPed(j,2)) then
                    RecGender(i)=1
                    exit
                endif
                if (i==RecPed(j,3)) then
                    RecGender(i)=2
                    exit
                endif
            enddo
        endif
    enddo
endif

deallocate(seqid)
deallocate(seqsire)
deallocate(seqdam)
allocate(Genotyped(nAnisP))
allocate(Pruned(0:nAnisP))
allocate(IsParent(nAnisP))
allocate(BaseAnimals(nAnisP))
allocate(TempGenos(0:nAnisP,nSnp))

TempGenos=9
Genotyped=0
Pruned=0
Pruned(0)=1
do i=1,nAnisG
    do j=1,nAnisP
        if (trim(GenotypeId(i))==trim(Id(j))) then
            TempGenos(j,:)=Genos(i,:)
            if (count(TempGenos(j,:)/=9)>0) Genotyped(j)=1
            exit
        endif
    enddo
enddo
deallocate(Genos)
IndivIsGenotyped(:)=Genotyped(:)

Pruned=0
do i=1,nAnisP
    If(RecPed(i,2)==0 .and. RecPed(i,3)==0 .and. Genotyped(RecPed(i,1))==0 ) Pruned(RecPed(i,1)) = 1
enddo
SumPruned=Sum(Pruned)
flag=1
k=0
do while (flag==1)
    flag=0
    k=k+1
    do i=1,nAnisP
        If(RecPed(i,2)==0 .and. Pruned(RecPed(i,3))==1 .and. Genotyped(RecPed(i,1))==0 ) Pruned(RecPed(i,1)) = 1
        If(RecPed(i,3)==0 .and. Pruned(RecPed(i,2))==1 .and. Genotyped(RecPed(i,1))==0 ) Pruned(RecPed(i,1)) = 1
    enddo
    if(sum(Pruned)>SumPruned) flag=1
    SumPruned=Sum(Pruned)
enddo

IsParent=.false.
do i=1,nAnisP
    If(RecPed(i,2)>0) IsParent(RecPed(i,2))=.true.
    If(RecPed(i,3)>0) IsParent(RecPed(i,3))=.true.
enddo

flag=1
k=0
do while (flag==1)
    flag=0
    k=k+1
    do i=1,nAnisP
        If(.not. IsParent(RecPed(i,1)) .and. Genotyped(RecPed(i,1))==0 ) Pruned(RecPed(i,1)) = 1
    enddo
    if(sum(Pruned)>SumPruned) then
        flag=1
        SumPruned=Sum(Pruned)
        IsParent=.false.
        do i=1,nAnisP
            if(Pruned(RecPed(i,1))==0) then
                If(RecPed(i,2)>0) IsParent(RecPed(i,2))=.true.
                If(RecPed(i,3)>0) IsParent(RecPed(i,3))=.true.
            endif
        enddo
    endif
enddo

BaseAnimals=0

! Some cases are considered twice, i.e.:
!       if (RecPed(i,2)==0) then
!         if (RecPed(i,3)==0) BaseAnimals(i)=1
! and:
!       if (RecPed(i,3)==0) then
!         if (RecPed(i,2)==0) BaseAnimals(i)=1
! can be transformed into:
! if (RecPed(i,3)==0.and.RecPed(i,2)==0) BaseAnimals(i)=1

! Determine whether an animal is a base animals (==without parents)
do i=1,nAnisP
    if (Genotyped(i)==1) then                               ! If genotyped
        ! If no sire
        if (RecPed(i,2)==0) then
            ! If no dam
            if (RecPed(i,3)==0) BaseAnimals(i)=1
            ! If dam has been pruned
            if (Pruned(RecPed(i,3))==1) BaseAnimals(i)=1
        endif
        ! If sire has been pruned
        if (Pruned(RecPed(i,2))==1) then
            ! If no dam
            if (RecPed(i,3)==0) BaseAnimals(i)=1
            ! If dam has been pruned
            if (Pruned(RecPed(i,3))==1) BaseAnimals(i)=1
        endif
        ! If no dam
        if (RecPed(i,3)==0) then
            ! If no sire
            if (RecPed(i,2)==0) BaseAnimals(i)=1
            ! If sire has been pruned
            if (Pruned(RecPed(i,2))==1) BaseAnimals(i)=1
        endif
        ! If dam has been pruned
        if (Pruned(RecPed(i,3))==1) then
            ! If no sire
            if (RecPed(i,2)==0) BaseAnimals(i)=1
            ! If sire has been pruned
            if (Pruned(RecPed(i,2))==1) BaseAnimals(i)=1
        endif
    endif
enddo

end subroutine CheckParentage

!#############################################################################################################################################################################################################################

subroutine CountInData
! Count the number of individuals genotyped and the number of individuals in the pedigree.
! If no pedigree information is available, then these to counts are equal to the number of 
! individuals genotyped.

use Global
implicit none

integer :: k
character (len=300) :: dumC

do
    read (2,*,iostat=k) dumC
    nAnisRawPedigree=nAnisRawPedigree+1
    if (k/=0) then
        nAnisRawPedigree=nAnisRawPedigree-1
        exit            ! This forces to exit if an error is found 
    endif
enddo
rewind(2)

print*, " ",nAnisRawPedigree," individuals in the pedigree file"
nObsDataRaw=nAnisRawPedigree

do
    read (3,*,iostat=k) dumC
    nAnisG=nAnisG+1
    if (k/=0) then
        nAnisG=nAnisG-1
        exit
    endif
enddo

rewind(3)

print*, " ",nAnisG," individuals in the genotype file"

! This is incoherent with functions ReadInData and ReadInParameterFile
if (trim(PedigreeFile)=="NoPedigree") nAnisRawPedigree=nAnisG

end subroutine CountInData

!#############################################################################################################################################################################################################################

subroutine ReadInData
use GlobalPedigree
use Global
implicit none

integer :: i,j,k,Temp(nSnp),CountLinesGender,GenCode,AnimalPresent
character(len=300) :: dumC

allocate(GenotypeId(nAnisG))
allocate(Ped(nAnisRawPedigree,3))
allocate(Genos(0:nAnisG,nSnp))
allocate(GenderId(nAnisRawPedigree))
allocate(GenderRaw(nAnisRawPedigree))

Genos(0,:)=9

! Read the pedigree information 
do i=1,nAnisRawPedigree
    read(2,*) ped(i,:)
enddo

do i=1,nAnisG
    read (3,*) GenotypeId(i),Temp(:)
    do j=1,nSnp
        if ((Temp(j)<0).or.(Temp(j)>2)) Temp(j)=9
    enddo
    Genos(i,:)=Temp(:)
enddo
close(2)
close(3)

GenderRaw=9
if (SexOpt==1) then
    CountLinesGender=0
    do
        read (4,*,iostat=k) dumC
        CountLinesGender=CountLinesGender+1
        if (k/=0) then
            CountLinesGender=CountLinesGender-1
            exit
        endif
    enddo
    rewind(4)
    nAnisInGenderFile=CountLinesGender
    if (CountLinesGender/=nAnisRawPedigree) then
        print*, "Warning - number of lines in Gender file not the same as in pedigree file"
        stop
    endif
    do j=1,nAnisRawPedigree                 ! For each individual in the file
        read (4,*) dumC,GenCode 
        if ((GenCode/=1).and.(GenCode/=2)) then
            print*, "Warning - Gender code incorrect for at least one animal"
            stop    
        endif
        AnimalPresent=0
        do i=1,nAnisRawPedigree             ! For each individual in the pedigree
            if (trim(dumC)==trim(ped(i,1))) then
                GenderId(i)=dumC
                GenderRaw(i)=GenCode
                AnimalPresent=1
                exit
            endif
        enddo
        if (AnimalPresent==0) then
            print*, "Warning - Animal missing in gender file"
            stop
        endif
    enddo   
endif

end subroutine ReadInData

!#############################################################################################################################################################################################################################

subroutine MakeDirectories
use global
implicit none

integer :: i
character(len=300) :: FolderName

print*, ""

call rmdir("Miscellaneous")
call rmdir("Phasing")
call rmdir("Results")
call rmdir("InputFiles")
call rmdir("GeneProb")      !here
call rmdir("IterateGeneProb")   !here

call system("mkdir Phasing")
call system("mkdir Miscellaneous")

call system("mkdir

 Results")
call system("mkdir InputFiles")
call system("mkdir GeneProb")       !here
call system("mkdir IterateGeneProb")    !here

if (WindowsLinux==1) then


else

    do i=1,nProcessors
        write (FolderName,'("GeneProb"i0)')i
        call system ("mkdir GeneProb/" // FolderName)       !here
    enddo
    do i=1,nPhaseInternal
        write (FolderName,'("Phase"i0)')i
        call system ("mkdir Phasing/" // FolderName)
    enddo
    do i=1,nProcessors
        write (FolderName,'("GeneProb"i0)')i            !here
        call system ("mkdir IterateGeneProb/" // FolderName)    !here
    enddo
endif


end subroutine MakeDirectories


!#############################################################################################################################################################################################################################

subroutine CountLines(FileName,nLines)
! Return the number of lines (nLines) in FileName
implicit none

character(len=*) :: FileName
character(len=300) :: dumC
integer :: nLines,f

nLines=0
open (unit=101,file=trim(FileName),status="old")
do
    read (101,*,iostat=f) dumC
    nLines=nLines+1
    if (f/=0) then
        nLines=nLines-1
        exit
    endif
enddo
close(101)

end subroutine CountLines

!#############################################################################################################################################################################################################################

subroutine PVseq(nObs,nAnisPedigree)

USE GlobalPedigree
implicit none
character (LEN=lengan), ALLOCATABLE :: holdsireid(:), holddamid(:)
character (LEN=lengan), ALLOCATABLE :: holdid(:), SortedId(:), SortedSire(:), SortedDam(:)
character (LEN=lengan)              :: IDhold
integer, ALLOCATABLE                :: SortedIdIndex(:), SortedSireIndex(:), SortedDamIndex(:)
integer, ALLOCATABLE                :: OldN(:), NewN(:), holdsire(:), holddam(:)
INTEGER :: mode    ! mode=1 to generate dummy ids where one parent known.  Geneprob->1  Matesel->0
INTEGER :: i, j, k, kk, newid, itth, itho, ihun, iten, iunit
integer :: nsires, ndams, newsires, newdams, nbisexuals, flag
INTEGER :: ns, nd, iextra, oldnobs, kn, kb, oldkn, ks, kd
INTEGER :: Noffset, Limit, Switch, ihold, ipoint
integer :: nObs,nAnisPedigree,verbose
character (LEN=lengan) :: path

mode=1 
allocate(id(0:nobs),sire(nobs),dam(nobs),seqid(nobs),seqsire(nobs),seqdam(nobs))

do i=1,nobs
        id(i)=ped(i,1)
        sire(i)=ped(i,2)
        dam(i)=ped(i,3)
end do

nAnisPedigree=nObs
path=".\"
Verbose=1

! Initialize and standarize
do j = 1, nobs
  If (dam(j) == ''.or. dam(j) == '0'.or. dam(j) == '#'.or. dam(j) == '*' .or. dam(j) == '.') Then
    dam(j) = '0'
    seqdam(j)=0
  endif
  If (sire(j) == ''.or.sire(j) == '0'.or.sire(j) == '#'.or.sire(j) == '*'.or.sire(j) == '.') Then
    sire(j) = '0'
    seqsire(j)=0
  endif
enddo

! Insert dummy IDs (mode=1)
if(mode.eq.1) then
  !PRINT*,  ' Inserting dummy IDs ... '
  newid=0
  do j = 1, nobs
    ! Count individuals with a single parent
    if(((sire(j)=='0').and.(dam(j).ne.'0')).or.((sire(j).ne.'0').and.(dam(j)=='0'))) then
       newid=newid+1
       
       if(newid.gt.99999) then
       !         PRINT*, newid, ' ...'
           stop 'too many dummy single parent IDs'
       endif

       ! Give dummy Id to missing parents
       itth=int(newid/10000)
       itho=int(newid/1000)-10*itth
       ihun=int(newid/100)-10*itho-100*itth
       iten=int(newid/10)-10*ihun-100*itho-1000*itth
       iunit=newid-10*iten-100*ihun-1000*itho-10000*itth
       if(sire(j)=='0') sire(j)='dum'//achar(48+itth)//achar(48+itho)//achar(48+ihun)//achar(48+iten)//achar(48+iunit)
       if( dam(j)=='0')  dam(j)='dum'//achar(48+itth)//achar(48+itho)//achar(48+ihun)//achar(48+iten)//achar(48+iunit)
    endif
  enddo
endif

!PRINT*,  ' Sorting Sires ... '
ALLOCATE  (SortedId(nobs), SortedIdIndex(nobs))
SortedId(1:nobs) = Sire(1:nobs)
Noffset = INT(nobs/2)
DO WHILE (Noffset>0)
    Limit = nobs - Noffset
    switch=1
    DO WHILE (Switch.ne.0)
       Switch = 0
       do i = 1, Limit
          IF (SortedId(i).gt.SortedId(i + Noffset)) THEN
               IDhold=SortedId(i)
               SortedId(i)=SortedId(i + Noffset)
               SortedId(i + Noffset)=IDhold
               Switch = i
          endif
       enddo
       Limit = Switch - Noffset
    enddo
    Noffset = INT(Noffset/2)
enddo

! Count the number of sires
nsires=0
IF(SortedId(1) /= '0') nsires=1
do i=2,nobs
    IF(SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') nsires=nsires+1
end do

ALLOCATE (SortedSire(0:nsires), SortedSireIndex(nsires))

! Sort sires
SortedSire(0)='0'
nsires=0
IF(SortedId(1) /= '0') THEN
    nsires=1
    SortedSire(1) = SortedId(1)
ENDIF
do i=2,nobs
  IF(SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') then
   nsires=nsires+1
   SortedSire(nsires) = SortedId(i)
  ENDIF
end do

!PRINT*,  ' Sorting Dams ... '
SortedId(1:nobs) = Dam(1:nobs)
  Noffset = INT(nobs/2)
  DO WHILE (Noffset>0)
      Limit = nobs - Noffset
      switch=1
    DO WHILE (Switch.ne.0)
       Switch = 0
       do i = 1, Limit
          IF (SortedId(i).gt.SortedId(i + Noffset)) THEN
               IDhold=SortedId(i)
               SortedId(i)=SortedId(i + Noffset)
               SortedId(i + Noffset)=IDhold
               Switch = i
          endif
       enddo
       Limit = Switch - Noffset
    enddo
    Noffset = INT(Noffset/2)
  enddo

! Count the number of dams
nDams=0
IF(SortedId(1) /= '0') nDams=1
do i=2,nobs
  IF(SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') nDams=nDams+1
end do

ALLOCATE (SortedDam(0:nDams), SortedDamIndex(ndams))

! Sort sires
SortedDam(0)='0'
nDams=0

IF(SortedId(1) /= '0') THEN
 nDams=1
 SortedDam(1) = SortedId(1)
ENDI
do i=2,nobs
  IF(SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') then
   nDams=nDams+1
   SortedDam(nDams) = SortedId(i)
  ENDIF
end do
! NOTE: Have to go through the number of Observations several times. 
!       We must find a way to avoid this and do it only a few times

!PRINT*,  ' Sorting IDs ... '
SortedId(1:nobs) = ID(1:nobs)
do i=1,nobs
 SortedIdIndex(i) = i
end do

Noffset = INT(nobs/2)
DO WHILE (Noffset>0)
    Limit=nobs-Noffset
    switch=1
    DO WHILE (Switch.ne.0)
        Switch=0
        do i=1,Limit
            IF (SortedId(i).gt.SortedId(i + Noffset)) THEN
               IDhold=SortedId(i)
               SortedId(i)=SortedId(i + Noffset)
               SortedId(i + Noffset)=IDhold
               ihold=SortedIdIndex(i)
               SortedIdIndex(i)=SortedIdIndex(i + Noffset)
               SortedIdIndex(i + Noffset)=ihold
               Switch = i
            endif
        enddo
        Limit = Switch - Noffset
    enddo
    Noffset = INT(Noffset/2)
enddo

!PRINT*,  ' Check for duplicate IDs ... '
flag = -1
Do i = 2, nobs
  If (SortedID(i) == SortedID(i - 1)) Then
   If (flag == -1) Then
     open (1,FILE='ID_err.txt',STATUS = 'unknown')
     WRITE(1,*) 'Duplicated IDs ...'
     flag = 0
   End If
   WRITE(1,*) SortedID(i)
   flag = flag + 1
  End If
enddo
If (flag > -1) Then
  Close (1)
  ! PRINT*, flag,' case(s) of duplicate ID. See ID_ERR.TXT        <------------ WARNING !!!'
End If


!PRINT*,  ' Males ... '
!PRINT*,  '  Find or set sire indices ... '
newsires = 0
do j=1,nsires
   ! check if already listed as an individual (Dichotomic search!)
   ipoint=INT(nobs/2)
   Noffset = INT(ipoint/2)
   do while (Noffset>1)
     IF (SortedSire(j).lt.SortedId(ipoint)) THEN
       ipoint = ipoint-Noffset
       Noffset = INT(Noffset/2)
     else
       ipoint = ipoint+Noffset
       Noffset = INT(Noffset/2)
     endif
   enddo

   kn=0
   if (SortedSire(j)==SortedId(ipoint)) kn=1 ! We've found j is listed as individual

   ! This is nosense! The Dichotomic search should be implemented in such a way that this search
   ! should not be necessary
   ! Looking forwards
   do while (ipoint<nobs.and.(kn==0).and.SortedSire(j)>SortedId(ipoint))
     ipoint=ipoint+1
   enddo
   
   if (SortedSire(j)==SortedId(ipoint)) kn=1    ! Found?
   
   ! Looking backwards
   do while ((ipoint>1).and.(kn==0).and.SortedSire(j)<SortedId(ipoint))
     ipoint=ipoint-1
   enddo
   
   if (SortedSire(j)==SortedId(ipoint)) kn=1    ! Found?
   
   ! If found, set sire index
   IF(kn==1) then
     SortedSireIndex(j) = SortedIdIndex(ipoint)
   else    ! otherwise, sire is unlisted base sire: Found a new sire
     newsires = newsires + 1
     SortedSireIndex(j) = nobs + newsires ! for now. 
   endif
end do !j

ALLOCATE  (holdsireid(newsires))
kn=0
  do j=1,nsires
    if (SortedSireIndex(j) > nobs) then
      kn=kn+1
      holdsireid(SortedSireIndex(j)-nobs) = SortedSire(j)
    end if
  enddo
IF(kn /= newsires) stop'newsires error'

!PRINT*,  '  Find seqsire ... '
do j = 1, nobs
  If (sire(j) == '0') Then
    seqsire(j)=0
  else
    ! Find sire (dichotomic search)
    ipoint=INT(nsires/2)
    Noffset = INT(ipoint/2)
    do while (Noffset>1)
      IF (Sire(j).lt.SortedSire(ipoint)) THEN
        ipoint = ipoint - Noffset
        Noffset = INT(Noffset/2)
      else
        ipoint = ipoint + Noffset
        Noffset = INT(Noffset/2)
      endif
    enddo
    kn=0
    if (Sire(j)==SortedSire(ipoint)) kn=1
      do while (ipoint<nsires .and. kn==0 .and. Sire(j) > SortedSire(ipoint))
        ipoint=ipoint+1
      enddo
      if (Sire(j)==SortedSire(ipoint)) kn=1
      do while (ipoint>1 .and. kn==0 .and. Sire(j) < SortedSire(ipoint))
        ipoint=ipoint-1
      enddo
      if (Sire(j)==SortedSire(ipoint)) kn=1
      IF(kn==1) then
        seqsire(j) = SortedSireIndex(ipoint)
      else
      !PRINT*, ' Error: Sire missing: ', Sire(j)
      stop
    endif
  endif
ENDDO !j

!PRINT*,  '  Sires: ',newsires,' unlisted, ',nsires,' in total'
!PRINT*,  ' Females ... '
!PRINT*,  '  Find or set dam indices ... '

newdams = 0
nbisexuals = 0
do j=1,ndams
! check if already listed as an individual
   ipoint=INT(nobs/2)
   Noffset = INT(ipoint/2)
   do while (Noffset>1)
    IF (Sorteddam(j).lt.SortedId(ipoint)) THEN
     ipoint = ipoint - Noffset
     Noffset = INT(Noffset/2)
    else
     ipoint = ipoint + Noffset
     Noffset = INT(Noffset/2)
    endif
   enddo
    kn=0
    if (Sorteddam(j)==SortedId(ipoint)) kn=ipoint  ! store ipoint here as ipoint can change with bisexuals
    do while (ipoint<nobs .and. kn==0 .and. Sorteddam(j) > SortedId(ipoint))
     ipoint=ipoint+1
    enddo
    if (Sorteddam(j)==SortedId(ipoint)) kn=ipoint
    do while (ipoint>1 .and. kn==0 .and. Sorteddam(j) < SortedId(ipoint))
     ipoint=ipoint-1
    enddo
    if (Sorteddam(j)==SortedId(ipoint)) kn=ipoint
! check if already listed as a sire (and therefore bisexual)
   ipoint=INT(nsires/2)
   Noffset = INT(ipoint/2)
   do while (Noffset>1)
    IF (SortedDam(j).lt.SortedSire(ipoint)) THEN
     ipoint = ipoint - Noffset
     Noffset = INT(Noffset/2)
    else
     ipoint = ipoint + Noffset
     Noffset = INT(Noffset/2)
    endif
   enddo
    kb=0
    if (SortedDam(j)==SortedSire(ipoint)) kb=1
    do while (ipoint<nsires .and. kb==0 .and. SortedDam(j) > SortedSire(ipoint))
     ipoint=ipoint+1
    enddo
    if (SortedDam(j)==SortedSire(ipoint)) kb=1
    do while (ipoint>1 .and. kb==0 .and. SortedDam(j) < SortedSire(ipoint))
     ipoint=ipoint-1
    enddo
    if (SortedDam(j)==SortedSire(ipoint)) kb=1

    IF(kb==1) then
      nbisexuals = nbisexuals + 1
      open (1,FILE='bisex.txt',position = 'append')
       WRITE(1,*) SortedDam(j)
      close(1)
    endif
    if (kb==1) then
     SorteddamIndex(j) = SortedSireIndex(ipoint)
    elseif (kn>=1) then
     SorteddamIndex(j) = SortedIdIndex(kn)
    else    ! dam is unlisted base dam
     newdams = newdams + 1
     SorteddamIndex(j) = nobs + newsires + newdams ! for now
    endif
end do !j

If (nbisexuals > 0)  PRINT*, nbisexuals,' bisexual parent(s) found. See file bisex.txt.  <------------ WARNING !!!'
 ALLOCATE  (holddamid(newdams))
 kn=0
 do j=1,ndams
  if (SortedDamIndex(j) > nobs+newsires) then
   kn=kn+1
   holddamid(SortedDamIndex(j)-nobs-newsires) = SortedDam(j)
  end if
 enddo
 IF(kn /= newdams) stop'newdams error'

!PRINT*,  '  Find seqdam ... '
do j = 1, nobs
  If (dam(j) == '0') Then
    seqdam(j)=0
  else
    ipoint=INT(ndams/2)
    Noffset = INT(ipoint/2)
   do while (Noffset>1)
    IF (dam(j).lt.Sorteddam(ipoint)) THEN
     ipoint = ipoint - Noffset
     Noffset = INT(Noffset/2)
    else
     ipoint = ipoint + Noffset
     Noffset = INT(Noffset/2)
    endif
   enddo
    kn=0
    if (dam(j)==Sorteddam(ipoint)) kn=1
    do while (ipoint<ndams .and. kn==0 .and. dam(j) > Sorteddam(ipoint))
     ipoint=ipoint+1
    enddo
    if (dam(j)==Sorteddam(ipoint)) kn=1
    do while (ipoint>1 .and. kn==0 .and. dam(j) < Sorteddam(ipoint))
     ipoint=ipoint-1
    enddo
    if (dam(j)==Sorteddam(ipoint)) kn=1
    IF(kn==1) then
     seqdam(j) = SorteddamIndex(ipoint)
    else
    ! PRINT*, ' Error: dam missing: ', dam(j)
     stop
    endif
  endif
ENDDO !j

!PRINT*,  '  Dams: ',newdams,' unlisted, ',ndams,' in total'
!PRINT*,  ' Arranging unlisted base parents ... '
iextra = newsires + newdams
If (iextra > 0) then
     ! PRINT*, ' ', iextra, ' unlisted base parents found.'
 ! SortedId and SortedIdIndex just used as a holder while redimensioning
 SortedId(1:nobs)=id(1:nobs)
 deallocate (id)
 ALLOCATE(id(nobs+iextra))
 id(1+iextra:nobs+iextra)=SortedId(1:nobs)

 SortedId(1:nobs)=sire(1:nobs)
 deallocate (sire)
 ALLOCATE(sire(nobs+iextra))
 sire(1+iextra:nobs+iextra)=SortedId(1:nobs)

 SortedId(1:nobs)=dam(1:nobs)
 deallocate (dam)
 ALLOCATE(dam(nobs+iextra))
 dam(1+iextra:nobs+iextra)=SortedId(1:nobs)

 SortedIdIndex(1:nobs)=seqsire(1:nobs)
 deallocate (seqsire)
 ALLOCATE(seqsire(nobs+iextra))
 seqsire(1+iextra:nobs+iextra)=SortedIdIndex(1:nobs)

 SortedIdIndex(1:nobs)=seqdam(1:nobs)
 deallocate (seqdam)
 ALLOCATE(seqdam(nobs+iextra))
 seqdam(1+iextra:nobs+iextra)=SortedIdIndex(1:nobs)
endif

!PRINT*, ' Inserting unlisted base parents ...'

oldnobs = nobs
nobs = nobs + iextra
!PRINT*, ' Total number of animals = ',nobs

ALLOCATE (passedorder(nobs))
passedorder=0

do i = 1+iextra, nobs
 passedorder(i)= i-iextra

 If (sire(i) == '0')then
   seqsire(i) = 0
 Else
   seqsire(i) = iextra + seqsire(i)
   If (seqsire(i) > nobs)  seqsire(i) = seqsire(i) - nobs  ! for unlisted sires
 End If

 If (dam(i) == '0') Then
   seqdam(i) = 0
  Else
   seqdam(i) = iextra + seqdam(i)
   If (seqdam(i) > nobs)  seqdam(i) = seqdam(i) - nobs
  End If
ENDDO !i

do i = 1, newsires
 ID(i) = holdsireid(i)
 passedorder(i)=0
 seqsire(i) = 0
 seqdam(i) = 0
ENDDO !i
do i = newsires + 1, newsires + newdams
 ID(i) = holddamid(i - newsires)
 passedorder(i)=0
 seqsire(i) = 0
 seqdam(i) = 0
ENDDO !i

DEALLOCATE(holdsireid, holddamid, SortedIdIndex, SortedId)

flag = 0
Do i = 1, nobs
If (i <= seqsire(i) .Or. i <= seqdam(i) ) flag = 1
enddo !i
!If (flag == 0) !PRINT*, 'not needed'!return

!PRINT*, ' Re-Ordering pedigree ...'


Allocate ( OldN(0:nobs), NewN(0:nobs) )
ALLOCATE ( holdid(0:nobs), holdsire(nobs), holddam(nobs) )

OldN(0) = 0
NewN=0
!seqsire(0) = 0 !not needed !
!seqdam(0) = 0

holdid(1:nobs) = ID(1:nobs)
holdsire = seqsire
holddam = seqdam

!Find base ancestors ...
kn = 0
do i = 1, nobs
 If (seqsire(i) == 0 .And. seqdam(i) == 0) Then
      kn = kn + 1
      NewN(i) = kn
      OldN(kn) = i
 End If
ENDDO !i

!Re-order pedigree ...
NewN(0) = nobs + 1
flag = 0
Do While (kn < nobs)
 oldkn = kn
 do i = 1, nobs
  If (NewN(i) == 0) Then !And ID(i) <> 'UniqueNULL' Then
    Ks = seqsire(i)
    Kd = seqdam(i)
    If (NewN(Ks) > 0 .And. NewN(Kd) > 0) Then
      kn = kn + 1
      NewN(i) = kn
      OldN(kn) = i
    End If
  End If
 enddo !i

 ! to avoid hang on unexpected problem ...
 If (kn == oldkn) Then
  flag = flag + 1
 Else
  flag = 0
 endif

 If (flag > 10) Then
   open(1,file='ped_err.txt',status='unknown')
   write(1,*) 'Pedigree errors found involving two or more of the following relationships ...'
   write(1,*)
   write(1,*) '       Index numbers are followed by names.'
   write(1,*) '       Index number 0 means unknown, whence name is blank.'
   write(1,*)
   do i = 1, nobs
    If (NewN(i) == 0) Then
     write(1,*) 'Individual:',          i, ':  ', ID(i)
     write(1,*) '    Father:', seqsire(i), ':  ', ID(seqsire(i))
     write(1,*) '    Mother:',  seqdam(i), ':  ', ID(seqdam(i))
     write(1,*)
    End If
   ENDDO !i
   Close (1)
   PRINT*,  'Logical error when re-ordering pedigree - see details in file PED_ERR.TXT'
   stop
 End If
ENDDO !while

NewN(0) = 0

do i = 1, nobs
 ID(i) = holdid(OldN(i))
enddo

do i = 1, nobs
seqsire(i) = NewN(holdsire(OldN(i)))
seqdam(i) = NewN(holddam(OldN(i)))
If (i <= NewN(holdsire(OldN(i))) .Or. i <= NewN(holddam(OldN(i)))) then
   !PRINT*,  'out of order'
   stop
endif
ENDDO !i

DO i = 1, nobs
  holdsire(i) = passedorder(i)  ! holdsire just because it is free
enddo

DO i = 1, nobs
  passedorder(i) = holdsire(OldN(i))
enddo

deallocate ( OldN, NewN, holdid, holdsire, holddam) ! holdrec)

!do i = 1, nobs
! PRINT'(3i5,2x,3a4,i5)', i, seqsire(i), seqdam(i), id(i), sire(i), dam(i), passedorder(i)
!enddo
nAnisPedigree=nObs
GlobalExtraAnimals=iextra   !Change John Hickey

end subroutine PVseq
!#############################################################################################################################################################################################################################

subroutine rmdir(tmpdir)
 
character(len=*) :: tmpdir

open(unit=1000,file=".rmdirsh",status="unknown")
write(1000,*) "if [ -d "// trim(tmpdir) //" ]"
write(1000,*) "then rm -r " // trim(tmpdir)
write(1000,*) "fi"
close(1000)

call system("chmod a+x .rmdirsh")
call system("./.rmdirsh")
call system("rm .rmdirsh")
 
end subroutine rmdir
!#############################################################################################################################################################################################################################

subroutine Checker
use Global
use GlobalPedigree
implicit none

integer :: h,i,j,k,l,f,nAnisTest,Work(nSnpRaw),WorkTmp(nSnpRaw),dum,GenoStratIndex(nAnisP),CountCatTest(6),CountCorr
integer :: SummaryStats(3,6),Div,CountLen,Counter
real :: SummaryProps(3,6),SumPat(6),SumMat(6)
character(len=300) :: Names(6),FileName,dumC
integer,allocatable,dimension(:) :: RecTestId,FinalSetter
integer,allocatable,dimension(:,:) :: TrueGenos,RawGenos,TestMat
real,allocatable,dimension(:) :: Correlations
real,allocatable,dimension(:,:) :: AnisSummary,WorkVec,RealTestGenos
character(len=lengan),allocatable,dimension(:) :: TrueGenosId

FileName=trim(TrueGenosFile)
call CountLines(FileName,nAnisTest)

if (WindowsLinux==1) then

        call system("rmdir /s /q TempTestAlphaImpute")

else

    call rmdir("TempTestAlphaImpute")
endif
call system("mkdir TempTestAlphaImpute")

open (unit=35,file=trim(TrueGenosFile),status="old")
open (unit=36,file=trim(GenotypeFile),status="unknown")
open (unit=37,file="./TempTestAlphaImpute/IndividualAnimalAccuracy.txt",status="unknown")
open (unit=38,file="./TempTestAlphaImpute/SummaryAnimalAccuracy.txt",status="unknown")
open (unit=44,file="./TempTestAlphaImpute/IndividualSummaryAccuracy.txt",status="unknown")
open (unit=45,file="./TempTestAlphaImpute/IndividualSummaryYield.txt",status="unknown")

Names(1)="Both Parents Genotyped"
Names(2)="Sire and Maternal GrandSire Genotyped"
Names(3)="Dam and Paternal Grandsire Genotyped"
Names(4)="Sire Genotyped"
Names(5)="Dam Genotyped"
Names(6)="Other Relatives Genotyped"

allocate(FinalSetter(0:nAnisP))
FinalSetter=0
do i=1,nAnisG
    read (36,*) dumC,WorkTmp(:)
    Counter=0
    do j=1,nSnpRaw
        if ((WorkTmp(j)>=0).and.(WorkTmp(j)<=2)) Counter=Counter+1
    enddo
    if (float(Counter)>(float(nSnpRaw)/2)) then
        do k=1,nAnisP
            if (trim(Id(k))==dumC) then
                FinalSetter(k)=1
                exit
            endif
        enddo
    endif
enddo
rewind(36)

if (OutOpt==0) then

    allocate(TrueGenos(nAnisTest,nSnp))
    allocate(TrueGenosId(nAnisTest))
    allocate(RawGenos(nAnisTest,nSnp))
    allocate(TestMat(nAnisTest,nSnp))
    allocate(RecTestId(nAnisTest))
    allocate(AnisSummary(nAnisTest,5))
    allocate(Correlations(6))
    allocate(RealTestGenos(nAnisTest,nSnp))
    do i=1,nAnisTest
        read (35,*) TrueGenosId(i),Work(:)
        k=0
        do j=1,nSnpRaw
            if (SnpIncluded(j)/=0) then
                k=k+1
                TrueGenos(i,k)=Work(j)
            endif
        enddo
    enddo
    do i=1,nAnisTest
        do j=1,nAnisP
            if (trim(TrueGenosId(i))==trim((Id(j)))) then
                RecTestId(i)=j
                exit
            endif
        enddo
    enddo
    do i=1,nAnisG
        read (36,*) dumC,WorkTmp(:)
        do j=1,nAnisTest
            if (trim(TrueGenosId(j))==dumC) then
                k=0
                do l=1,nSnpRaw
                    if (SnpIncluded(l)==1) then
                        k=k+1
                        RawGenos(j,k)=WorkTmp(l)
                    endif
                enddo
                exit
            endif
        enddo
    enddo
    GenoStratIndex(:)=0
    do i=1,nAnisP
        if (FinalSetter(i)/=1) then
            GenoStratIndex(i)=6
            if (FinalSetter(RecPed(i,3))==1) then
                GenoStratIndex(i)=5
                if (FinalSetter(RecPed(RecPed(i,2),2))==1) then
                    GenoStratIndex(i)=3
                endif
            endif
            if (FinalSetter(RecPed(i,2))==1) then
                GenoStratIndex(i)=4
                if (FinalSetter(RecPed(RecPed(i,3),2))==1) then
                    GenoStratIndex(i)=2
                endif
            endif
            if ((FinalSetter(RecPed(i,2))==1).and.(FinalSetter(RecPed(i,3))==1)) then
                GenoStratIndex(i)=1
            endif
        endif
    enddo
    TestMat=4
    CountCatTest=0
    AnisSummary=0.0
    do i=1,nAnisTest
        do j=1,nSnp
            if ((TrueGenos(i,j)<0).or.(TrueGenos(i,j)>2)) then
                TestMat(i,j)=5 
            else
                if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
                    if (TrueGenos(i,j)==ImputeGenos(RecTestId(i),j)) TestMat(i,j)=1
                    if (TrueGenos(i,j)/=ImputeGenos(RecTestId(i),j)) TestMat(i,j)=2
                    if (ImputeGenos(RecTestId(i),j)==9) TestMat(i,j)=3
                endif   
            endif
        enddo
        write (37,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),TestMat(i,:)
        CountCatTest(GenoStratIndex(RecTestId(i)))=CountCatTest(GenoStratIndex(RecTestId(i)))+1
        Div=(count(TestMat(i,:)==1))+(count(TestMat(i,:)==2))+(count(TestMat(i,:)==3))
        AnisSummary(i,1)=100*(float(count(TestMat(i,:)==1))/Div)
        AnisSummary(i,2)=100*(float(count(TestMat(i,:)==2))/Div)
        AnisSummary(i,3)=100*(float(count(TestMat(i,:)==3))/Div)
        AnisSummary(i,4)=100*(float(count(ImputePhase(RecTestId(i),:,1)/=9))/nSnp)
        AnisSummary(i,5)=100*(float(count(ImputePhase(RecTestId(i),:,2)/=9))/nSnp)
        write (44,'(a20,i3,5f7.2)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),AnisSummary(i,:)
    enddo
    SummaryStats=0
    SumPat=0.0
    SumMat=0.0
    do i=1,nAnisTest
        do j=1,6
            if (GenoStratIndex(RecTestId(i))==j) then
                SummaryStats(1,j)=SummaryStats(1,j)+(count(TestMat(i,:)==1))
                SummaryStats(2,j)=SummaryStats(2,j)+(count(TestMat(i,:)==2))
                SummaryStats(3,j)=SummaryStats(3,j)+(count(TestMat(i,:)==3))
                SumPat(j)=SumPat(j)+AnisSummary(i,4)
                SumMat(j)=SumMat(j)+AnisSummary(i,5)
            endif
        enddo
    enddo

    SummaryProps=0.0
    do i=1,3
        do j=1,6
            if (CountCatTest(j)==0) then
                SummaryProps(i,j)=0.0
            else
                SummaryProps(i,j)=float(SummaryStats(i,j))/sum(SummaryStats(:,j))
                SummaryProps(i,j)=SummaryProps(i,j)*100
            endif
        enddo
    enddo

    do h=1,6
        CountLen=0
        do i=1,nAnisTest
            if(GenoStratIndex(RecTestId(i))==h) then
                do j=1,nSnp
                    if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
                        if ((TrueGenos(i,j)>=0).and.(TrueGenos(i,j)<=2)) then
                            CountLen=CountLen+1
                        endif
                    endif
                enddo
            endif
        enddo
    enddo

    do i=1,6
        SumPat(i)=SumPat(i)/CountCatTest(i)
        SumMat(i)=SumMat(i)/CountCatTest(i)
        write (38,'(5f7.2,i7,a40)') SummaryProps(:,i),SumPat(i),SumMat(i),CountCatTest(i),trim(Names(i))
    enddo

    print*, " "
    do i=1,6
        if (CountCatTest(i)>0) write (*,'(3f7.2,a3,a3,2f7.2,i7,a40)') SummaryProps(:,i),"   "&
                    ,"   ",SumPat(i),SumMat(i),CountCatTest(i),trim(Names(i))
    enddo

    do i=GlobalExtraAnimals+1,nAnisP
        write (45,'(a25,i3,2f7.2)') Id(i),FinalSetter(i),float(count(ImputePhase(i,:,1)/=9))/nSnp &
                            ,float(count(ImputePhase(i,:,2)/=9))/nSnp
    enddo
else
    allocate(TrueGenos(nAnisTest,nSnpRaw))
    allocate(TrueGenosId(nAnisTest))
    allocate(RawGenos(nAnisTest,nSnpRaw))
    allocate(TestMat(nAnisTest,nSnpRaw))
    allocate(RecTestId(nAnisTest))
    allocate(AnisSummary(nAnisTest,5))

    do i=1,nAnisTest
        read (35,*) TrueGenosId(i),TrueGenos(i,:)
    enddo

    do i=1,nAnisTest
        do j=1,nAnisP
            if (trim(TrueGenosId(i))==trim((Id(j)))) then
                RecTestId(i)=j
                exit
            endif
        enddo
    enddo

    do i=1,nAnisG
        read (36,*) dumC,WorkTmp(:)
        do j=1,nAnisTest
            if (trim(TrueGenosId(j))==dumC) then
                RawGenos(j,:)=WorkTmp(:)
                exit
            endif
        enddo
    enddo
    GenoStratIndex(:)=0
    do i=1,nAnisP
        if (FinalSetter(i)/=1) then
            GenoStratIndex(i)=6
            if (FinalSetter(RecPed(i,3))==1) then
                GenoStratIndex(i)=5
                if (FinalSetter(RecPed(RecPed(i,2),2))==1) then
                    GenoStratIndex(i)=3
                endif
            endif
            if (FinalSetter(RecPed(i,2))==1) then
                GenoStratIndex(i)=4
                if (FinalSetter(RecPed(RecPed(i,3),2))==1) then
                    GenoStratIndex(i)=2
                endif
            endif
            if ((FinalSetter(RecPed(i,2))==1).and.(FinalSetter(RecPed(i,3))==1)) then
                GenoStratIndex(i)=1
            endif
        endif
    enddo
    TestMat=4
    CountCatTest=0
    AnisSummary=0.0
    do i=1,nAnisTest
        do j=1,nSnpRaw
            if ((TrueGenos(i,j)<0).or.(TrueGenos(i,j)>2)) then
                TestMat(i,j)=5 
            else    
                if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
                    if (TrueGenos(i,j)==ImputeGenos(RecTestId(i),j)) TestMat(i,j)=1
                    if (TrueGenos(i,j)/=ImputeGenos(RecTestId(i),j)) TestMat(i,j)=2
                    if (ImputeGenos(RecTestId(i),j)==9) TestMat(i,j)=3
                endif
            endif   
        enddo
        write (37,'(a20,i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),TestMat(i,:)
        CountCatTest(GenoStratIndex(RecTestId(i)))=CountCatTest(GenoStratIndex(RecTestId(i)))+1
        Div=(count(TestMat(i,:)==1))+(count(TestMat(i,:)==2))+(count(TestMat(i,:)==3))
        AnisSummary(i,1)=100*(float(count(TestMat(i,:)==1))/Div)
        AnisSummary(i,2)=100*(float(count(TestMat(i,:)==2))/Div)
        AnisSummary(i,3)=100*(float(count(TestMat(i,:)==3))/Div)
        AnisSummary(i,4)=100*(float(count(ImputePhase(RecTestId(i),:,1)/=9))/nSnpRaw)
        AnisSummary(i,5)=100*(float(count(ImputePhase(RecTestId(i),:,2)/=9))/nSnpRaw)
        write (44,'(a20,i3,5f7.2)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),AnisSummary(i,:)
    enddo

    SummaryStats=0
    SumPat=0.0
    SumMat=0.0
    do i=1,nAnisTest
        do j=1,6
            if (GenoStratIndex(RecTestId(i))==j) then
                SummaryStats(1,j)=SummaryStats(1,j)+(count(TestMat(i,:)==1))
                SummaryStats(2,j)=SummaryStats(2,j)+(count(TestMat(i,:)==2))
                SummaryStats(3,j)=SummaryStats(3,j)+(count(TestMat(i,:)==3))
                SumPat(j)=SumPat(j)+AnisSummary(i,4)
                SumMat(j)=SumMat(j)+AnisSummary(i,5)
            endif
        enddo
    enddo

    SummaryProps=0.0
    do i=1,3
        do j=1,6
            if (CountCatTest(j)==0) then
                SummaryProps(i,j)=0.0
            else
                SummaryProps(i,j)=float(SummaryStats(i,j))/sum(SummaryStats(:,j))
                SummaryProps(i,j)=SummaryProps(i,j)*100
            endif
        enddo
    enddo

    do h=1,6
        CountLen=0
        do i=1,nAnisTest
            if(GenoStratIndex(RecTestId(i))==h) then
                do j=1,nSnpRaw
                    if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
                        if ((TrueGenos(i,j)>=0).and.(TrueGenos(i,j)<=2)) then
                            CountLen=CountLen+1
                        endif
                    endif
                enddo
            endif
        enddo
    enddo

    do i=1,6
        SumPat(i)=SumPat(i)/CountCatTest(i)
        SumMat(i)=SumMat(i)/CountCatTest(i)
        write (38,'(5f7.2,i7,a40)') SummaryProps(:,i),SumPat(i),SumMat(i),CountCatTest(i),trim(Names(i))
    enddo

    print*, " "
    do i=1,6
        if (CountCatTest(i)>0) write (*,'(3f7.2,a3,a3,2f7.2,i7,a40)') SummaryProps(:,i),"   "&
                    ,"   ",SumPat(i),SumMat(i),CountCatTest(i),trim(Names(i))
    enddo

    do i=GlobalExtraAnimals+1,nAnisP
        write (45,'(a25,i3,2f7.2)') Id(i),FinalSetter(i),float(count(ImputePhase(i,:,1)/=9))/nSnpRaw &
                            ,float(count(ImputePhase(i,:,2)/=9))/nSnpRaw
    enddo

endif
close(35)
close(36)
close(37)
close(38)
close(44)
close(45)


end subroutine Checker

!#############################################################################################################################################################################################################################

subroutine FinalChecker
use Global
use GlobalPedigree
implicit none

integer :: h,i,j,k,l,f,nAnisTest,Work(nSnpRaw),WorkTmp(nSnpRaw),dum,GenoStratIndex(nAnisP),CountCatTest(6),CountCorr
integer :: SummaryStats(3,6),Div,CountLen,Counter,Top1,Top2,Top3,Top4,Bot,ContSnpCor,CountValAnim(6)
double precision :: SummaryProps(3,6),SumPat(6),SumMat(6),MeanCorPerInd(6),StdDevPerGrp(6),AveCategoryInformativeness(6,6)
double precision :: Tmpave,Tmpadev,Tmpvar,Tmpskew,Tmpcurt
character(len=300) :: Names(6),FileName,dumC
integer,allocatable,dimension(:) :: RecTestId,FinalSetter
integer,allocatable,dimension(:,:) :: TrueGenos,RawGenos,TestMat,TestAnimInformativeness
double precision,allocatable,dimension(:) :: Correlations,CorrelationPerAnimal,TmpVarPerGrp
double precision,allocatable,dimension(:,:) :: AnisSummary,WorkVec,RealTestGenos,CalcCorPerAnimal
character(len=lengan),allocatable,dimension(:) :: TrueGenosId

FileName=trim(TrueGenosFile)
call CountLines(FileName,nAnisTest)

if (WindowsLinux==1) then

     call system("rmdir /s /q TestAlphaImpute")

else

    call rmdir("TestAlphaImpute")
endif
call system("mkdir TestAlphaImpute")

open (unit=35,file=trim(TrueGenosFile),status="old")
open (unit=36,file=trim(GenotypeFile),status="unknown")
open (unit=37,file="./TestAlphaImpute/IndividualAnimalAccuracy.txt",status="unknown")
open (unit=38,file="./TestAlphaImpute/SummaryAnimalAccuracy.txt",status="unknown")
open (unit=44,file="./TestAlphaImpute/IndividualSummaryAccuracy.txt",status="unknown")
open (unit=45,file="./TestAlphaImpute/IndividualSummaryYield.txt",status="unknown")
open (unit=48,file="./TestAlphaImpute/IndividualSnpAccuracy.txt",status="unknown")

Names(1)="Both Parents Genotyped"
Names(2)="Sire and Maternal GrandSire Genotyped"
Names(3)="Dam and Paternal Grandsire Genotyped"
Names(4)="Sire Genotyped"
Names(5)="Dam Genotyped"
Names(6)="Other Relatives Genotyped"

allocate(FinalSetter(0:nAnisP))
FinalSetter=0
do i=1,nAnisG
    read (36,*) dumC,WorkTmp(:)
    Counter=0
    do j=1,nSnpRaw
        if ((WorkTmp(j)>=0).and.(WorkTmp(j)<=2)) Counter=Counter+1
    enddo
    if (float(Counter)>(float(nSnpRaw)/2)) then
        do k=1,nAnisP
            if (trim(Id(k))==dumC) then
                FinalSetter(k)=1
                exit
            endif
        enddo
    endif
enddo
rewind(36)

allocate(TestAnimInformativeness(nAnisTest,6))

if (OutOpt==0) then

    allocate(TrueGenos(nAnisTest,nSnp))
    allocate(TrueGenosId(nAnisTest))
    allocate(RawGenos(nAnisTest,nSnp))
    allocate(TestMat(nAnisTest,nSnp))
    allocate(RecTestId(nAnisTest))
    allocate(AnisSummary(nAnisTest,5))
    allocate(Correlations(6))
    allocate(RealTestGenos(nAnisTest,nSnp))
    allocate(CalcCorPerAnimal(nSnp,2))
    allocate(CorrelationPerAnimal(nAnisTest))
    allocate(TmpVarPerGrp(nAnisTest))

    do i=1,nAnisTest
        read (35,*) TrueGenosId(i),Work(:)
        k=0
        do j=1,nSnpRaw
            if (SnpIncluded(j)/=0) then
                k=k+1
                TrueGenos(i,k)=Work(j)
            endif
        enddo
    enddo

    RecTestId(:)=-99
    do i=1,nAnisTest
        do j=1,nAnisP
            if (trim(TrueGenosId(i))==trim((Id(j)))) then
                RecTestId(i)=j
                TestAnimInformativeness(i,:)=GlobalTmpCountInf(j,1:6)
                exit
            endif
        enddo
    enddo
    if (count(RecTestId(:)==-99)>0) print*, "Error - There seems to be unidentifiablecount ",count(RecTestId(:)==-99)," individuals in the test file"

    do i=1,nAnisG
        read (36,*) dumC,WorkTmp(:)
        do j=1,nAnisTest
            if (trim(TrueGenosId(j))==dumC) then
                k=0
                do l=1,nSnpRaw
                    if (SnpIncluded(l)==1) then
                        k=k+1
                        RawGenos(j,k)=WorkTmp(l)
                    endif
                enddo
                exit
            endif
        enddo
    enddo

    GenoStratIndex(:)=0
    do i=1,nAnisP
        if (FinalSetter(i)/=1) then
            GenoStratIndex(i)=6
            if (FinalSetter(RecPed(i,3))==1) then
                GenoStratIndex(i)=5
                if (FinalSetter(RecPed(RecPed(i,2),2))==1) then
                    GenoStratIndex(i)=3
                endif
            endif
            if (FinalSetter(RecPed(i,2))==1) then
                GenoStratIndex(i)=4
                if (FinalSetter(RecPed(RecPed(i,3),2))==1) then
                    GenoStratIndex(i)=2
                endif
            endif
            if ((FinalSetter(RecPed(i,2))==1).and.(FinalSetter(RecPed(i,3))==1)) then
                GenoStratIndex(i)=1
            endif
        endif
    enddo

    TestMat=4
    CountCatTest=0
    AnisSummary=0.0

    do i=1,nAnisTest
        do j=1,nSnp
            if ((TrueGenos(i,j)<0).or.(TrueGenos(i,j)>2)) then
                TestMat(i,j)=5 
            else
                if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
                    if (TrueGenos(i,j)==ImputeGenos(RecTestId(i),j)) TestMat(i,j)=1
                    if (TrueGenos(i,j)/=ImputeGenos(RecTestId(i),j)) TestMat(i,j)=2
                    if (ImputeGenos(RecTestId(i),j)==9) TestMat(i,j)=3
                endif   
            endif
        enddo
        RealTestGenos(i,:)=ProbImputeGenos(RecTestId(i),:)
        write (37,'(a20,i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),TestMat(i,:)
        CountCatTest(GenoStratIndex(RecTestId(i)))=CountCatTest(GenoStratIndex(RecTestId(i)))+1
        Div=(count(TestMat(i,:)==1))+(count(TestMat(i,:)==2))+(count(TestMat(i,:)==3))
        AnisSummary(i,1)=100*(float(count(TestMat(i,:)==1))/Div)
        AnisSummary(i,2)=100*(float(count(TestMat(i,:)==2))/Div)
        AnisSummary(i,3)=100*(float(count(TestMat(i,:)==3))/Div)
        AnisSummary(i,4)=100*(float(count(ImputePhase(RecTestId(i),:,1)/=9))/nSnp)
        AnisSummary(i,5)=100*(float(count(ImputePhase(RecTestId(i),:,2)/=9))/nSnp)
    enddo

    SummaryStats=0
    SumPat=0.0
    SumMat=0.0
    do i=1,nAnisTest
        do j=1,6
            if (GenoStratIndex(RecTestId(i))==j) then
                SummaryStats(1,j)=SummaryStats(1,j)+(count(TestMat(i,:)==1))
                SummaryStats(2,j)=SummaryStats(2,j)+(count(TestMat(i,:)==2))
                SummaryStats(3,j)=SummaryStats(3,j)+(count(TestMat(i,:)==3))
                SumPat(j)=SumPat(j)+AnisSummary(i,4)
                SumMat(j)=SumMat(j)+AnisSummary(i,5)
            endif
        enddo
    enddo

    SummaryProps=0.0
    do i=1,3
        do j=1,6
            if (CountCatTest(j)==0) then
                SummaryProps(i,j)=0.0
            else
                SummaryProps(i,j)=float(SummaryStats(i,j))/sum(SummaryStats(:,j))
                SummaryProps(i,j)=SummaryProps(i,j)*100
            endif
        enddo
    enddo

    MeanCorPerInd(:)=0.0
    CountValAnim(:)=0
    AveCategoryInformativeness(:,:)=0.0 
    do h=1,6
        CountLen=0
        do i=1,nAnisTest
            if(GenoStratIndex(RecTestId(i))==h) then
                AveCategoryInformativeness(h,:)=AveCategoryInformativeness(h,:)+TestAnimInformativeness(i,:)
                do j=1,nSnp
                    if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
                        if ((TrueGenos(i,j)>=0).and.(TrueGenos(i,j)<=2)) then
                            CountLen=CountLen+1
                        endif
                    endif
                enddo
            endif
        enddo
        allocate(WorkVec(CountLen,2))
        CountLen=0
        do i=1,nAnisTest
            if(GenoStratIndex(RecTestId(i))==h) then
                ContSnpCor=0
                do j=1,nSnp
                    if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
                        if ((TrueGenos(i,j)>=0).and.(TrueGenos(i,j)<=2)) then
                            ContSnpCor=ContSnpCor+1                                         !#HereToday
                            CountLen=CountLen+1
                            WorkVec(CountLen,1)=float(TrueGenos(i,j))-(2*maf(j))
                            WorkVec(CountLen,2)=RealTestGenos(i,j)-(2*maf(j))
                            CalcCorPerAnimal(ContSnpCor,1)=float(TrueGenos(i,j))-(2*maf(j))
                            CalcCorPerAnimal(ContSnpCor,2)=RealTestGenos(i,j)-(2*maf(j))                            
                        endif
                    endif
                enddo
                if (ContSnpCor>5) then
                    call Pearsn (CalcCorPerAnimal(1:ContSnpCor,1),CalcCorPerAnimal(1:ContSnpCor,2),ContSnpCor,CorrelationPerAnimal(i))
                    MeanCorPerInd(h)=MeanCorPerInd(h)+CorrelationPerAnimal(i)
                    CountValAnim(h)=CountValAnim(h)+1
                    TmpVarPerGrp(CountValAnim(h))=CorrelationPerAnimal(i)

                else
                    CorrelationPerAnimal(i)=-99.0
                endif

            endif
        enddo
        if (CountLen>5) then
            call Pearsn (WorkVec(:,1),WorkVec(:,2),CountLen,Correlations(h))
        else
            Correlations(h)=-99.0
        endif
        deallocate(WorkVec)
        if (CountValAnim(h)>5) then
            MeanCorPerInd(h)=MeanCorPerInd(h)/CountValAnim(h)
            call moment(TmpVarPerGrp(1:CountValAnim(h)),CountValAnim(h),Tmpave,Tmpadev,StdDevPerGrp(h),Tmpvar,Tmpskew,Tmpcurt) 
        else
            MeanCorPerInd(h)=-99.0      
        endif   
    enddo

    do i=1,6
        SumPat(i)=SumPat(i)/CountCatTest(i)
        SumMat(i)=SumMat(i)/CountCatTest(i)
        AveCategoryInformativeness(i,:)=AveCategoryInformativeness(i,:)/CountCatTest(i)
        write (38,'(14f7.2,i7,a40)') SummaryProps(:,i),Correlations(i),MeanCorPerInd(i),StdDevPerGrp(i),SumPat(i),SumMat(i),AveCategoryInformativeness(i,:),CountCatTest(i),trim(Names(i))
    enddo

    print*, " "
    do i=1,6
        if (CountCatTest(i)>0) write (*,'(3f7.2,a3,3f7.2,a3,8f7.2,i7,a40)') SummaryProps(:,i),"   ",Correlations(i),MeanCorPerInd(i),StdDevPerGrp(i)&
                    ,"   ",SumPat(i),SumMat(i),AveCategoryInformativeness(i,:),CountCatTest(i),trim(Names(i))
    enddo

    do i=GlobalExtraAnimals+1,nAnisP
        write (45,'(a25,i3,2f7.2)') Id(i),FinalSetter(i),float(count(ImputePhase(i,:,1)/=9))/nSnp &
                            ,float(count(ImputePhase(i,:,2)/=9))/nSnp
    enddo
    
    do j=1,nSnp
        Bot=count(TestMat(:,j)/=4)      
        Top1=count(TestMat(:,j)==1) 
        Top2=count(TestMat(:,j)==2) 
        Top3=count(TestMat(:,j)==3) 
        Top4=count(TestMat(:,j)==5)     
        write (48,'(2i10,4f9.2)') j,Bot,100*(float(Top1)/Bot),100*(float(Top2)/Bot),100*(float(Top3)/Bot),100*(float(Top4)/Bot)
    enddo
    
    do i=1,nAnisTest
        write (44,'(a20,i3,6f7.2,6i10)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),AnisSummary(i,:),CorrelationPerAnimal(i),TestAnimInformativeness(i,:)
    enddo
    
    
else
    allocate(TrueGenos(nAnisTest,nSnpRaw))
    allocate(TrueGenosId(nAnisTest))
    allocate(RawGenos(nAnisTest,nSnpRaw))
    allocate(TestMat(nAnisTest,nSnpRaw))
    allocate(RecTestId(nAnisTest))
    allocate(AnisSummary(nAnisTest,5))
    allocate(Correlations(6))
    allocate(RealTestGenos(nAnisTest,nSnpRaw))
    allocate(CalcCorPerAnimal(nSnpRaw,2))   
    allocate(CorrelationPerAnimal(nAnisTest))
    allocate(TmpVarPerGrp(nAnisTest))

    do i=1,nAnisTest
        read (35,*) TrueGenosId(i),TrueGenos(i,:)
    enddo

    RecTestId(:)=-99
    do i=1,nAnisTest
        do j=1,nAnisP
            if (trim(TrueGenosId(i))==trim((Id(j)))) then
                RecTestId(i)=j
                TestAnimInformativeness(i,:)=GlobalTmpCountInf(j,1:6)
                exit
            endif
        enddo
    enddo
    if (count(RecTestId(:)==-99)>0) print*, "Error - There seems to be unidentifiablecount ",count(RecTestId(:)==-99)," individuals in the test file"

    do i=1,nAnisG
        read (36,*) dumC,WorkTmp(:)
        do j=1,nAnisTest
            if (trim(TrueGenosId(j))==dumC) then
                RawGenos(j,:)=WorkTmp(:)
                exit
            endif
        enddo
    enddo
    GenoStratIndex(:)=0
    do i=1,nAnisP
        if (FinalSetter(i)/=1) then
            GenoStratIndex(i)=6
            if (FinalSetter(RecPed(i,3))==1) then
                GenoStratIndex(i)=5
                if (FinalSetter(RecPed(RecPed(i,2),2))==1) then
                    GenoStratIndex(i)=3
                endif
            endif
            if (FinalSetter(RecPed(i,2))==1) then
                GenoStratIndex(i)=4
                if (FinalSetter(RecPed(RecPed(i,3),2))==1) then
                    GenoStratIndex(i)=2
                endif
            endif
            if ((FinalSetter(RecPed(i,2))==1).and.(FinalSetter(RecPed(i,3))==1)) then
                GenoStratIndex(i)=1
            endif
        endif
    enddo
    TestMat=4
    CountCatTest=0
    AnisSummary=0.0
    do i=1,nAnisTest
        do j=1,nSnpRaw
            if ((TrueGenos(i,j)<0).or.(TrueGenos(i,j)>2)) then
                TestMat(i,j)=5 
            else    
                if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
                    if (TrueGenos(i,j)==ImputeGenos(RecTestId(i),j)) TestMat(i,j)=1
                    if (TrueGenos(i,j)/=ImputeGenos(RecTestId(i),j)) TestMat(i,j)=2
                    if (ImputeGenos(RecTestId(i),j)==9) TestMat(i,j)=3
                endif
            endif   
        enddo
        RealTestGenos(i,:)=ProbImputeGenos(RecTestId(i),:)
        write (37,'(a20,i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),TestMat(i,:)
        CountCatTest(GenoStratIndex(RecTestId(i)))=CountCatTest(GenoStratIndex(RecTestId(i)))+1
        Div=(count(TestMat(i,:)==1))+(count(TestMat(i,:)==2))+(count(TestMat(i,:)==3))
        AnisSummary(i,1)=100*(float(count(TestMat(i,:)==1))/Div)
        AnisSummary(i,2)=100*(float(count(TestMat(i,:)==2))/Div)
        AnisSummary(i,3)=100*(float(count(TestMat(i,:)==3))/Div)
        AnisSummary(i,4)=100*(float(count(ImputePhase(RecTestId(i),:,1)/=9))/nSnpRaw)
        AnisSummary(i,5)=100*(float(count(ImputePhase(RecTestId(i),:,2)/=9))/nSnpRaw)
    enddo

    SummaryStats=0
    SumPat=0.0
    SumMat=0.0
    do i=1,nAnisTest
        do j=1,6
            if (GenoStratIndex(RecTestId(i))==j) then
                SummaryStats(1,j)=SummaryStats(1,j)+(count(TestMat(i,:)==1))
                SummaryStats(2,j)=SummaryStats(2,j)+(count(TestMat(i,:)==2))
                SummaryStats(3,j)=SummaryStats(3,j)+(count(TestMat(i,:)==3))
                SumPat(j)=SumPat(j)+AnisSummary(i,4)
                SumMat(j)=SumMat(j)+AnisSummary(i,5)
            endif
        enddo
    enddo

    SummaryProps=0.0
    do i=1,3
        do j=1,6
            if (CountCatTest(j)==0) then
                SummaryProps(i,j)=0.0
            else
                SummaryProps(i,j)=float(SummaryStats(i,j))/sum(SummaryStats(:,j))
                SummaryProps(i,j)=SummaryProps(i,j)*100
            endif
        enddo
    enddo

    MeanCorPerInd(:)=0.0
    CountValAnim(:)=0
    AveCategoryInformativeness(:,:)=0.0
    do h=1,6
        CountLen=0
        do i=1,nAnisTest
            if(GenoStratIndex(RecTestId(i))==h) then
                AveCategoryInformativeness(h,:)=AveCategoryInformativeness(h,:)+TestAnimInformativeness(i,:)
                do j=1,nSnpRaw
                    if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
                        if ((TrueGenos(i,j)>=0).and.(TrueGenos(i,j)<=2)) then
                            CountLen=CountLen+1
                        endif
                    endif
                enddo
            endif
        enddo
        allocate(WorkVec(CountLen,2))
        CountLen=0
        do i=1,nAnisTest
            if(GenoStratIndex(RecTestId(i))==h) then
                ContSnpCor=0
                do j=1,nSnpRaw
                    if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
                        if ((TrueGenos(i,j)>=0).and.(TrueGenos(i,j)<=2)) then
                            ContSnpCor=ContSnpCor+1
                            CountLen=CountLen+1
                            WorkVec(CountLen,1)=float(TrueGenos(i,j))-(2*maf(j))
                            WorkVec(CountLen,2)=RealTestGenos(i,j)-(2*maf(j))
                            CalcCorPerAnimal(ContSnpCor,1)=float(TrueGenos(i,j))-(2*maf(j))
                            CalcCorPerAnimal(ContSnpCor,2)=RealTestGenos(i,j)-(2*maf(j))                            
                        endif
                    endif
                enddo
                if (ContSnpCor>5) then
                    call Pearsn (CalcCorPerAnimal(1:ContSnpCor,1),CalcCorPerAnimal(1:ContSnpCor,2),ContSnpCor,CorrelationPerAnimal(i))
                    MeanCorPerInd(h)=MeanCorPerInd(h)+CorrelationPerAnimal(i)
                    CountValAnim(h)=CountValAnim(h)+1
                    TmpVarPerGrp(CountValAnim(h))=CorrelationPerAnimal(i)
                else
                    CorrelationPerAnimal(i)=-99.0
                endif               
            endif
        enddo               
        
        if (CountLen>5) then
            call Pearsn(WorkVec(:,1),WorkVec(:,2),CountLen,Correlations(h))
        else
            Correlations(h)=-99.0
        endif

        deallocate(WorkVec)

        if (CountValAnim(h)>5) then
            MeanCorPerInd(h)=MeanCorPerInd(h)/CountValAnim(h)
            call moment(TmpVarPerGrp(1:CountValAnim(h)),CountValAnim(h),Tmpave,Tmpadev,StdDevPerGrp(h),Tmpvar,Tmpskew,Tmpcurt)          
        else
            MeanCorPerInd(h)=-99.0      
        endif   
    enddo

    do i=1,6
        SumPat(i)=SumPat(i)/CountCatTest(i)
        SumMat(i)=SumMat(i)/CountCatTest(i)
        AveCategoryInformativeness(i,:)=AveCategoryInformativeness(i,:)/CountCatTest(i)
        write (38,'(14f7.2,i7,a40)') SummaryProps(:,i),Correlations(i),MeanCorPerInd(i),StdDevPerGrp(i),SumPat(i),SumMat(i),AveCategoryInformativeness(i,:),CountCatTest(i),trim(Names(i))
    enddo

    print*, " "
    do i=1,6
        if (CountCatTest(i)>0) write (*,'(3f7.2,a3,3f7.2,a3,8f7.2,i7,a40)') SummaryProps(:,i),"   ",Correlations(i),MeanCorPerInd(i),StdDevPerGrp(i)&
                    ,"   ",SumPat(i),SumMat(i),AveCategoryInformativeness(i,:),CountCatTest(i),trim(Names(i))
    enddo

    do i=GlobalExtraAnimals+1,nAnisP
        write (45,'(a25,i3,2f7.2)') Id(i),FinalSetter(i),float(count(ImputePhase(i,:,1)/=9))/nSnpRaw &
                            ,float(count(ImputePhase(i,:,2)/=9))/nSnpRaw
    enddo

    do j=1,nSnp
        Bot=count(TestMat(:,j)/=4)      
        Top1=count(TestMat(:,j)==1) 
        Top2=count(TestMat(:,j)==2) 
        Top3=count(TestMat(:,j)==3) 
        Top4=count(TestMat(:,j)==5)     
        write (48,'(2i10,4f9.2)') j,Bot,100*(float(Top1)/Bot),100*(float(Top2)/Bot),100*(float(Top3)/Bot),100*(float(Top4)/Bot)                 
    enddo
    
    do i=1,nAnisTest
        write (44,'(a20,i3,6f7.2,6i10)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),AnisSummary(i,:),CorrelationPerAnimal(i),TestAnimInformativeness(i,:)  
    enddo
    
endif

end subroutine FinalChecker

!#############################################################################################################################################################################################################################

subroutine Pearsn (x,y,n,r)

implicit none
integer n
double precision :: prob,r,z,x(n),y(n),TINY
parameter (tiny=1.e-20)
integer j
double precision :: ax,ay,df,sxx,sxy,syy,t,xt,betai,yt

ax=0.0
ay=0.0
DO j=1,n
        ax=ax+x(j)
        ay=ay+y(j)
END DO
ax=ax/n
ay=ay/n
sxx=0.
syy=0.
sxy=0.
DO j=1,n
        xt=x(j)-ax
        yt=y(j)-ay
        sxx=sxx+xt**2
        syy=syy+yt**2
        sxy=sxy+xt*yt
END DO
r=sxy/(SQRT(sxx*syy)+TINY)
z=0.5*LOG(((1.+r)+TINY)/((1.-r)+TINY))
df=n-2
t=r*SQRT(df/(((1.-r)+TINY)*((1.+r)+TINY)))
!prob=betai(0.5*df,0.5,df/(df+t**2))
!prob=erfcc(ABS(z*SQRT(n-1.))/1.4142136)
prob=0
return

end subroutine Pearsn

!###########################################################################################################################################################

SUBROUTINE moment(DATA,n,ave,adev,sdev,var,skew,curt)
IMPLICIT NONE
INTEGER n
DOUBLE PRECISION adev,ave,curt,sdev,skew,var,DATA(n)
INTEGER j
DOUBLE PRECISION p,s,ep
IF (n.le.1) PAUSE 'n must be at least 2 in moment'
s=0
DO j= 1,n
        s=s+DATA(j)
END DO

ave=s/n
adev=0
var=0
skew=0
curt=0
ep=0

DO j=1,n
        s=DATA(j)-ave
        ep=ep+s
        adev=adev+ABS(s)
        p=s*s
        var=var+p
        p=p*s
        skew=skew+p
        p=p*s
        curt=curt+p
END DO

adev=adev/n
var=(var-ep**2/n)/(n-1)
sdev=SQRT(var)
IF(var.ne.0)then
        skew=skew/(n*sdev**3)
        curt=curt/(n*var**2)-3
ELSE
        !PRINT*, 'no skew or kurtosis when zero variance in moment'
        !PAUSE 'no skew or kurtosis when zero variance in moment'
END IF
RETURN
END SUBROUTINE moment

!#############################################################################################################################################################################################################################

subroutine Cleaner
use Global
use GlobalPedigree
implicit none

call rmdir("GeneProb")      
call rmdir("IterateGeneProb")   

if (SexOpt==0) call system(" rm TempGeneProb.sh")
if (SexOpt==0) call system("rm TempIterateGeneProb.sh") 

end subroutine Cleaner

!#############################################################################################################################################################################################################################

subroutine Titles

print*, ""
print*, "                              ***********************                         "
print*, "                              *                     *                         "
print*, "                              *     AlphaImpute     *                         "
print*, "                              *      Beta 1.21      *                         "
print*, "                              *                     *                         "
print*, "                              ***********************                         "
print*, "                                                                              "
print*, "                    Software For Phasing and Imputing Genotypes               "
print*, ""
print*, "  Written by John Hickey, Matt Cleveland, Andreas Kranis, and Brian Kinghorn  "
print*, ""
print*, ""
print*, ""
print*, ""

end subroutine Titles

!#############################################################################################################################################################################################################################

subroutine PrintTimerTitles
use Global
implicit none

real :: etime          ! Declare the type of etime()
real :: elapsed(2)     ! For receiving user and system time
real :: total,Minutes,Hours,Seconds

print*, ""
print*, ""
print*, ""
print*, ""
print*, "                              ***********************                         "
print*, "                              *                     *                         "
print*, "                              *     AlphaImpute     *                         "
print*, "                              *                     *                         "
print*, "                              ***********************                         "
print*, "                                                                              "
print*, "                    Software For Phasing and Imputing Genotypes               "
print*, ""
print*, "  Written by John Hickey, Matt Cleveland, Andreas Kranis, and Brian Kinghorn  "
print*, ""
print*, "                                  No Liability"
print*, "                          Bugs to John.Hickey@une.edu.au"
print*, ""
print*, "                Analysis Finished                         "

total=etime(elapsed)
Minutes=total/60
Seconds=Total-(INT(Minutes)*60)
Hours=Minutes/60
Minutes=INT(Minutes)-(INT(Hours)*60)
print '(A107,A7,I3,A9,I3,A9,F6.2)', "Time Elapsed","Hours", INT(Hours),"Minutes",INT(Minutes),"Seconds",Seconds

if (WindowsLinux==1) then
        open (unit=32,file=".\Miscellaneous\Timer.txt",status="unknown")
else
        open (unit=32,file="./Miscellaneous/Timer.txt",status="unknown")
endif

write(32,'(A27,A7,I3,A9,I3,A9,F6.2)') "Time Elapsed","Hours", INT(Hours),"Minutes",INT(Minutes),"Seconds",Seconds

end subroutine PrintTimerTitles

!#############################################################################################################################################################################################################################

subroutine READINJUNK

use Global

integer :: i,dum


open (3001,file="fort.2008",status="old")
do i=1,nAnisP
    read (3001,*) dum,ImputePhase(i,:,1)
    read (3001,*) dum,ImputePhase(i,:,2)
enddo

do i=1,nAnisP
    do j=1,nSnp
        if ((ImputePhase(i,j,1)/=9).and.(ImputePhase(i,j,2)/=9)) then
            if (ImputeGenos(i,j)==9) ImputeGenos(i,j)=ImputePhase(i,j,1)+ImputePhase(i,j,2)
        endif
    enddo
enddo

end subroutine READINJUNK

!#############################################################################################################################################################################################################################

subroutine RandomOrder(order,n,idum)
 implicit none

!     Generate a random ordering of the integers 1 ... n.

integer, INTENT(IN)  :: n
integer, INTENT(OUT) :: order(n)
integer :: idum
double precision ran1

!     Local variables

integer :: i, j, k
double precision    :: wk

do i = 1, n
  order(i) = i
end do

!     Starting at the end, swap the current last indicator with one
!     randomly chosen from those preceeding it.

do i = n, 2, -1
  wk=ran1(idum)
  j = 1 + i * wk
  if (j < i) then
    k = order(i)
    order(i) = order(j)
    order(j) = k
  end if
end do

RETURN
end subroutine RandomOrder

!#############################################################################################################################################################################################################################

! This Function returns a uniform random deviate between 0.0 and 1.0.
! Set IDUM to any negative value to initialize or reinitialize the sequence.
!MODIFIED FOR REAL

FUNCTION ran1(idum)
IMPLICIT NONE
 INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
 DOUBLE PRECISION ran1,AM,EPS,RNMX
 PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
 INTEGER j,k,iv(NTAB),iy
 SAVE iv,iy
 DATA iv /NTAB*0/, iy /0/
  IF (idum.le.0.or.iy.eq.0) then
      idum=max(-idum,1)
  DO 11 j=NTAB+8,1,-1
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
  IF (idum.lt.0) idum=idum+IM
  IF (j.le.NTAB) iv(j)=idum

11 CONTINUE
     iy=iv(1)
  END IF
     k=idum/IQ
     idum=IA*(idum-k*IQ)-IR*k
  IF (idum.lt.0) idum=idum+IM
     j=1+iy/NDIV
     iy=iv(j)
     iv(j)=idum
     ran1=min(AM*iy,RNMX)
  RETURN
END

!  (C) Copr. 1986-92 Numerical Recipes Software 6

!#############################################################################################################################################################################################################################
