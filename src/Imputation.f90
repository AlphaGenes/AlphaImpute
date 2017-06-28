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


MODULE Imputation
    use Global

    use AlphaImputeSpecFileModule
    use AlphaImputeInputOutputModule
    use ConstantModule
    implicit none

    type(AlphaImputeInput), pointer :: inputParams
CONTAINS

    SUBROUTINE ImputationManagement
        use omp_lib
        use AlphaImputeInputOutputModule, only : ReadInPrePhasedData, ReReadGeneProbs,readgeneprobscluster
        use AlphaPhaseResultsModule

        integer :: loop
        character(len=150) :: timeOut

        inputParams => defaultInput



        call ped%writeOutGenotypes("Results/" //"startimpMan")    



        ! WARNING: Need to disuss this part of code with John. Nonsense going on here!

        if (inputParams%HMMOption==RUN_HMM_ONLY) then ! Avoid any adulteration of genotypes with imputation subroutines

#ifdef DEBUG
write(0,*) 'DEBUG: Allocate memory for genotypes and haplotypes'
#endif     

            allocate(GlobalTmpCountInf(ped%pedigreeSize,8))

            allocate(MSTermInfo(ped%pedigreeSize,2))

#ifdef DEBUG
write(0,*) 'DEBUG: Call Mach'
#endif

            block
                use AlphaHmmInMod
                use ExternalHMMWrappers
                type (AlphaHMMinput) :: inputParamsHMM
                integer(kind=1) ,dimension(:,:), allocatable :: genos
                integer(kind=1) ,dimension(:,:,:), allocatable :: phase
                
                inputParamsHMM%nsnp = inputParams%nsnp
                inputParamsHMM%nHapInSubH = inputParams%nHapInSubH
                inputParamsHMM%HmmBurnInRound = inputParams%HmmBurnInRound
                inputParamsHMM%nRoundsHmm = inputParams%nRoundsHmm
                inputParamsHMM%useProcs = inputParams%useProcs
                inputParamsHMM%imputedThreshold = inputParams%imputedThreshold
                inputParamsHMM%phasedThreshold = inputParams%phasedThreshold
                inputParamsHMM%HapList = inputParams%HapList


                call AlphaImputeHMMRunner(inputParamsHMM, ped, ProbImputeGenosHmm, ProbImputePhaseHmm, GenosCounts, FullH)


            end block
            call FromHMM2ImputePhase

#ifdef DEBUG
write(0,*) 'DEBUG: Mach Finished'
#endif

        else

            if (inputParams%sexopt==0) then    
                if (inputParams%BypassGeneProb==0) then


                    if (inputParams%cluster) then
                        call readGeneProbsCluster(GlobalWorkPhase,ped,GpIndex, inputParams,GeneProbThresh)
                    else   
                        call ReReadGeneProbs(globalworkphase, ped,"./Results/GenotypeProbabilities.txt", inputParams%nsnp,GeneProbThresh)
                    endif
                else
                    ! Phase in the homozygous case for the SEX CHROMOSOME
                    ! WARNING: NOTHING IS DONE!!
                    call InsteadOfReReadGeneProb
                endif


                ! Get Genotype information
                call GeneProbPhase          ! Recover and store information about which and how many alleles/SNPs have been genotyped/phased
            else
                allocate(MSTermInfo(ped%pedigreeSize,2))
                MSTermInfo=0
            endif

            if (inputParams%NoPhasing==1) then
                ! Major sub-step 2 as explained in Hickey et al. (2012; Appendix A)
                call BaseAnimalFillIn
                 call ped%writeOutGenotypes("Results/" //"afterbase1")  
                    ! Impute phase whenever a pre-phase file exists
                    if (inputParams%PrePhased==1) call ReadInPrePhasedData

                    ! Impute phase in the sex chromosome
                    if (inputParams%sexopt==1) call EnsureHetGametic

                    ! General imputation procedures
                    call GeneralFillInInit

                    call ped%writeOutGenotypes("Results/" //"aftergeneral1")  



                        if (inputParams%HMMOption==RUN_HMM_PREPHASE) Then
                            block
                                use AlphaHmmInMod
                                use ExternalHMMWrappers
                                type (AlphaHMMinput) :: inputParamsHMM

                                inputParamsHMM%nsnp = inputParams%nsnp
                                inputParamsHMM%nHapInSubH = inputParams%nHapInSubH
                                inputParamsHMM%HmmBurnInRound = inputParams%HmmBurnInRound
                                inputParamsHMM%nRoundsHmm = inputParams%nRoundsHmm
                                inputParamsHMM%useProcs = inputParams%useProcs
                                inputParamsHMM%imputedThreshold = inputParams%imputedThreshold
                                inputParamsHMM%phasedThreshold = inputParams%phasedThreshold
                                inputParamsHMM%HapList = inputParams%HapList



                                call AlphaImputeHMMRunner(inputParamsHMM, ped, ProbImputeGenosHmm, ProbImputePhaseHmm, GenosCounts, FullH)


                            end block
                        else
                            print*, " "
                            print*, " ","Imputation of base animals completed"


                            call ped%WriteoutPhase("Results/" //"beforeLoops")
                            call ped%writeOutGenotypes("Results/" //"beforeLoopsgeno")    

                            do loop=1,inputParams%InternalIterations
                                print*, " "
                                print*, "Performing imputation loop",loop

                                call PhaseElimination                   ! Major Sub-Step 5 (Hickey et al., 2012; Appendix A)
                                

                                    call ped%WriteoutPhase("Results/" //"afterPhase")

                                if (inputParams%sexopt==1) then
                                    call EnsureHetGametic
                                end if
                                call GeneralFillIn

                                    call ped%WriteoutPhase("Results/" //"afterGeneral")

                                print*, " "
                                print*, " ","Parent of origin assigmnent of high density haplotypes completed"

                                
                                call ParentPhaseElimination             ! Major Sub-Step 4 (Hickey et al., 2012; Appendix A)
                                if (inputParams%sexopt==1) then
                                    call EnsureHetGametic
                                end if

                                 call ped%WriteoutPhase("Results/" //"after4")


                                call GeneralFillIn
                                call RestrictedWorkLeftRight            ! Major Sub-Step 8 (Hickey et al., 2012; Appendix A)
                                if (inputParams%sexopt==1) then
                                    call EnsureHetGametic
                                end if
                                call GeneralFillIn

                                 call ped%WriteoutPhase("Results/" //"after81")

                                print*, " "
                                CALL DATE_AND_TIME(time=timeOut)
                                print*, " ","Imputation from high-density parents completed at: ",trim(timeOut)

                                call ImputeFromHDLibrary                ! Major Sub-Step 3 (Hickey et al., 2012; Appendix A)

                                if (inputParams%sexopt==1) then
                                    call EnsureHetGametic
                                end if

                                 call ped%WriteoutPhase("Results/" //"after3")

                                call GeneralFillIn
                                call RestrictedWorkLeftRight            ! Major Sub-Step 8 (Hickey et al., 2012; Appendix A)
                                call ped%WriteoutPhase("Results/" //"after82")
                                if (inputParams%sexopt==1) then
                                    call EnsureHetGametic
                                end if
                                call GeneralFillIn
                                CALL DATE_AND_TIME(time=timeOut)

                                 

                                print*, " "
                                print*, " ","Haplotype library imputation completed at: ",trim(timeOut)

                                call InternalParentPhaseElim            ! Major Sub-Step 7 (Hickey et al., 2012; Appendix A)

                                if (inputParams%sexopt==1) then
                                    call EnsureHetGametic
                                end if
                                call GeneralFillIn

                                 call ped%WriteoutPhase("Results/" //"after7")


                                call RestrictedWorkLeftRight            ! Major Sub-Step 8 (Hickey et al., 2012; Appendix A)
                                                                 call ped%WriteoutPhase("Results/" //"after83")
                                if (inputParams%sexopt==1) then
                                    call EnsureHetGametic
                                end if
                                call GeneralFillIn
                                print*, " "
                                CALL DATE_AND_TIME(time=timeOut)


                                print*, " ","Internal imputation from parents haplotype completed at: ",timeOut

                                call InternalHapLibImputationOLD           ! Major Sub-Step 6 (Hickey et al., 2012; Appendix A)
                                if (inputParams%sexopt==1) then
                                    call EnsureHetGametic
                                end if
                                call GeneralFillIn
                                if (inputParams%sexopt==1) then
                                    call EnsureHetGametic
                                end if


                                 call ped%WriteoutPhase("Results/" //"after6")

                                call RestrictedWorkLeftRight            ! Major Sub-Step 8 (Hickey et al., 2012; Appendix A)
                                call GeneralFillIn
                                print*, " "
                                CALL DATE_AND_TIME(time=timeOut)
                                
                                 call ped%WriteoutPhase("Results/" //"after84")

                                print*, " ","Internal haplotype library imputation completed at: ", timeOut

                                call ped%WriteoutPhase("Results/" //"endLoop")    
                                 call ped%writeOutGenotypes("Results/" //"endLoopGeno")      
                             enddo

                            call ManageWorkLeftRight

                        endif


                        if (inputParams%sexopt==1) then
                            call EnsureHetGametic
                        end if
                        call GeneralFillIn

                        if (inputParams%HMMOption==RUN_HMM_YES) Then

                            block
                                use AlphaHmmInMod
                                use ExternalHMMWrappers
                                type (AlphaHMMinput) :: inputParamsHMM

                                inputParamsHMM%nsnp = inputParams%nsnp
                                inputParamsHMM%nHapInSubH = inputParams%nHapInSubH
                                inputParamsHMM%HmmBurnInRound = inputParams%HmmBurnInRound
                                inputParamsHMM%nRoundsHmm = inputParams%nRoundsHmm
                                inputParamsHMM%useProcs = inputParams%useProcs
                                inputParamsHMM%imputedThreshold = inputParams%imputedThreshold
                                inputParamsHMM%phasedThreshold = inputParams%phasedThreshold
                                inputParamsHMM%HapList = inputParams%HapList
                                call AlphaImputeHMMRunner(inputParamsHMM, ped, ProbImputeGenosHmm, ProbImputePhaseHmm, GenosCounts, FullH)


                            end block
                            call FromHMM2ImputePhase
                        endif
                        deallocate(GlobalWorkPhase)
                    endif
                endif

                call ped%writeOutGenotypes("Results/" //"endimp")    
            END SUBROUTINE ImputationManagement

            subroutine InternalParentPhaseElim
                ! Internal imputation of single alleles based on parental phase.
                ! Several different core lengths are used to define the length of the haplotypes and to ensure to
                ! prevent the use of phasing errors that originate from LRPHLI. For each core of each round of the
                ! LRPHLI algorithm, all haplotypes that have been found and stored in the haplotype library.
                ! Candidate haplotypes are restricted to the two haplotypes that have been identified for each of
                ! the individual parents with high-density genotype information by the LRPHLI algorithm for the true
                ! haplotype that an individual carries on its gametes. Within the core, all alleles that are known
                ! are compared to corresponding alleles in each of the individual haplotypes. The candidate
                ! haplotypes are checked for locations that have unanimous agreement about a particular allele. For
                ! alleles with complete agreement, a count of the suggested allele is incremented. Alleles are
                ! imputed if, at the end of passing across each core and each round of the LRPHLI algorithm, the
                ! count of whether the alleles are 0 or 1 is above a threshold in one direction and below a
                ! threshold in the other. This helps to prevent the use of phasing errors that originate from
                ! LRPHLI.
                ! This subroutine corresponds to Major sub-step 7 from Hickey et al., 2012 (Appendix A)
                use, intrinsic :: ISO_Fortran_Env
                use individualModule

                implicit none

                integer :: m,e,g,i,j,nCore,nGlobalLoop,CoreLength,CoreStart,CoreEnd
                integer :: LoopStart,Offset,AnimalOn(ped%pedigreeSize,2)
                integer ::GamA,GamB

                integer,allocatable,dimension (:,:) :: LoopIndex
                integer :: l
                integer(kind=1),allocatable,dimension (:,:,:,:) :: Temp
                type(individual), pointer :: parent
                type(haplotype) :: subset 
                inputParams => defaultInput
                ! WARNING: This should go in a function since it is the same code as InternalParentPhaseElim subroutine
                nGlobalLoop=25

                ! LoopeIndex is a matrix with two columns that will serve to define:
                !   * 1.- nCores
                !   * 2.- Core lengths
                allocate(LoopIndex(nGlobalLoop,2))
                call setLoopIndex(LoopIndex)
                ! WARNING: This can be better arrange with a ELSEIF statement and should go in a function since it
                !          is the same code as InternalParentPhaseElim subroutine
                ! LoopStart indicates which is the first loop the algorithm should treat. The bigger the number of
                ! SNPs, the more the loops to be considered

                LoopStart = getLoopStart(inputParams%nsnp)
                if (LoopStart == 0) return

                ! Assumed that LoopIndex(:,1) are the numbers of cores for each phase step, LoopIndex):,2) are the core lengths
                do i=1,nGlobalLoop
                    LoopIndex(i,2)=int(float(inputParams%nsnp)/LoopIndex(i,1))
                enddo

                allocate(Temp(ped%pedigreeSize,inputParams%nsnp,2,2))
                Temp=0
                AnimalOn=0

                ! SIMULATE PHASING
                ! m is a variable to simulate shift or no-shift phasing
                do m=1,2
                    do l=LoopStart,nGlobalLoop

                        ! Simulate phase without shift
                        if (m==1) then
                            nCore=inputParams%nsnp/LoopIndex(l,2)
                            CoreStart=1
                            CoreEnd=LoopIndex(l,2)
                        else ! Simulate phase with shift
                            OffSet=int(float(LoopIndex(l,2))/2)
                            nCore=(inputParams%nsnp-(2*OffSet))/LoopIndex(l,2)
                            CoreStart=1+Offset
                            CoreEnd=LoopIndex(l,2)+Offset
                        endif

                           do g=1,nCore
                            ! Make sure that cores ends correctly
                            if ((m==1).and.(g==nCore)) CoreEnd=inputParams%nsnp
                            if ((m==2).and.(g==nCore)) CoreEnd=inputParams%nsnp-OffSet

                            ! Exit if the corelength is too small
                            CoreLength=(CoreEnd-CoreStart)+1
                            if (CoreLength<10) exit

                            ! PARALLELIZATION BEGINS
                            !# PARALLEL DO SHARED (nAnisP,RecPed,CoreStart,CoreEnd,AnimalOn,Temp) private(i,e,CompPhase,GamA,j,GamB)
                            do i=1,ped%pedigreeSize-ped%ndummys
                                do e=1,2
                                    ! Skip if, in the case of sex chromosome, me and my parent are heterogametic
                                    if (ped%pedigree(i)%isDummyBasedOnIndex(e+1)) cycle
                                    parent => ped%pedigree(i)%getSireDamObjectByIndex(e+1)
                                    if ((inputParams%SexOpt==1).and.(ped%pedigree(i)%gender==inputParams%HetGameticStatus).and.(parent%gender==inputParams%HetGameticStatus)) cycle

                                    ! If not a Base Animal
                                    if (.not. ped%pedigree(i)%founder) then
                                        ! If the haplotype for this core is not completely phased
                                        subset = ped%pedigree(i)%individualPhase(e)%subset(coreStart,coreEnd)
                                        if (.not. subset%fullyPhased())  then

                                            ! Check is this haplotype is the very same that the paternal haplotype
                                            ! of the parent of the individual
                                            GamA=1


                                            if (.not. subset%compatible(ped%pedigree(parent%id)%individualPhase(1)%subset(coreStart,CoreEnd),0,1)) then
                                                GamA=0
                                            endif
                                              ! Check is this haplotype is the very same that the maternal haplotype
                                            ! of the parent of the individual
                                            GamB=1

                                            if (.not. subset%compatible(ped%pedigree(parent%id)%individualPhase(2)%subset(coreStart,CoreEnd),0,1)) then
                                                GamB=0
                                            endif

                                            ! This haplotype is the paternal haplotype of the individual's parent
                                            ! Then count the number of occurrences a particular phase is impute in a
                                            ! a particular allele across the cores and across the internal phasing steps
                                            ! WARNING: This chunk of code and the next chunk can be colapse in a
                                            !          DO statement. Look in InternalHapLibImputation for an example
                                            if ((GamA==1).and.(GamB==0)) then
                                                AnimalOn(i,e)=1
                                                do j=CoreStart,CoreEnd
                                                    if (ped%pedigree(i)%individualPhase(e)%isMissing(j)) then
                                                    
                                                        block

                                                        integer(kind=1) :: parentPhase

                                                        parentPhase = parent%individualPhase(1)%getPhase(j)
                                                        
                                                        if (parentPhase==0) then
                                                            !$OMP ATOMIC
                                                            Temp(i,j,e,1)=Temp(i,j,e,1)+1
                                                        else if (parentPhase==1) then
                                                            !$OMP ATOMIC
                                                            Temp(i,j,e,2)=Temp(i,j,e,2)+1
                                                        endif
                                                        end block
                                                    endif
                                                enddo
                                            endif

                                            ! This haplotype is the maternal haplotype of the individual's parent
                                            ! Then count the number of occurrences a particular phase is impute in a
                                            ! a particular allele across the cores and across the internal phasing steps
                                            ! WARNING: This chunk of code and the previous chunk can be colapse in a
                                            !          DO statement. Look in InternalHapLibImputation for an example
                                            if ((GamA==0).and.(GamB==1)) then
                                                AnimalOn(i,e)=1
                                                do j=CoreStart,CoreEnd
                                                    if (ped%pedigree(i)%individualPhase(e)%isMissing(j)) then
                                                        block

                                                        integer(kind=1) :: parentPhase

                                                        parentPhase = parent%individualPhase(1)%getPhase(j)

                                                        if (parentPhase==0) then
                                                            Temp(i,j,e,1)=Temp(i,j,e,1)+1
                                                        else if (parentPhase==1) then
                                                            Temp(i,j,e,2)=Temp(i,j,e,2)+1
                                                        endif
                                                        end block
                                                    endif
                                                enddo
                                            endif
                                        endif
                                    endif
                                enddo
                            enddo
                            !# END PARALLEL DO

                            ! Prepare the core for the next cycle
                            CoreStart=CoreStart+LoopIndex(l,2)
                            CoreEnd=CoreEnd+LoopIndex(l,2)
                            if ((m==2).and.(g==nCore)) exit
                        enddo
                    enddo
                enddo

                ! If all alleles across the cores and across the internal phasing steps have been phased the same way, impute
                ! The individual has two haplotypes
            do i=1,ped%pedigreeSize- ped%nDummys
                do e=1,2
                
                    if (AnimalOn(i,e)==1) then
                        if ((inputParams%sexopt==0).or.(ped%pedigree(i)%gender==inputParams%HomGameticStatus)) then
                            do j=1,inputParams%nsnp
                                    
                                if (ped%pedigree(i)%individualPhase(e)%isMissing(j)) then
                                    if ((Temp(i,j,e,1)>inputParams%nAgreeInternalHapLibElim).and.(Temp(i,j,e,2)==0)) then
                                        call ped%pedigree(i)%individualPhase(e)%setPhase(j,0)
                                    else if ((Temp(i,j,e,1)==0).and.(Temp(i,j,e,2)>inputParams%nAgreeInternalHapLibElim)) then
                                        call ped%pedigree(i)%individualPhase(e)%setPhase(j,1)
                                    endif
                                endif
                                enddo
                        end if
                    endif
            

                ! The individual is has one haplotype: In Sex Chromosome, the heterogametic case
                
                    if ((inputParams%sexopt==1).and.(ped%pedigree(i)%gender==inputParams%HetGameticStatus)) then
                        if (AnimalOn(i,inputParams%HomGameticStatus)==1) then
                            do j=1,inputParams%nsnp
                                if (ped%pedigree(i)%individualPhase(inputParams%HomGameticStatus)%isMissing(j)) then
                                    if ((Temp(i,j,inputParams%HomGameticStatus,1)>inputParams%nAgreeInternalHapLibElim).and.(Temp(i,j,inputParams%HomGameticStatus,2)==0))&
                                        call ped%pedigree(i)%individualPhase(1)%setPhase(j,0)
                                        call ped%pedigree(i)%individualPhase(2)%setPhase(j,0)
                                    if ((Temp(i,j,inputParams%HomGameticStatus,1)==0).and.(Temp(i,j,inputParams%HomGameticStatus,2)>inputParams%nAgreeInternalHapLibElim))&
                                        call ped%pedigree(i)%individualPhase(1)%setPhase(j,1)
                                        call ped%pedigree(i)%individualPhase(2)%setPhase(j,1)
                                endif
                            enddo
                        end if
                    end if
                enddo
            enddo


                deallocate(Temp)


            end subroutine InternalParentPhaseElim

            !#############################################################################################################################################################################################################################


subroutine InternalHapLibImputationOld
! Internal candidate haplotype library imputation of alleles.q
! Haplotype libraries are internally built using the information that has been previously imputed.
! Several different core lengths are used to define the length of the haplotypes and to ensure to
! prevent the use of phasing errors that originate from LRPHLI. For each core, all haplotypes that
! have been found and stored in the haplotype library are initially considered to be candidates for
! the true haplotype that an individual carries on its gametes. Within the core, all alleles that
! are known are compared to corresponding alleles in each of the haplotypes in the library.
! Haplotypes that have a number of disagreements greater than a small error threshold have their
! candidacy rejected. At the end of this loop, the surviving candidate haplotypes are checked for
! locations that have unanimous agreement about a particular allele. For alleles with complete
! agreement, a count of the suggested allele is incremented. Alleles are imputed if, at the end of
! passing across each core and each round of the LRPHLI algorithm, the count of whether the alleles
! are 0 or 1 is above a threshold in one direction and below a threshold in the other.
! This subroutine corresponds to Major sub-step 6 from Hickey et al., 2012 (Appendix A)

use Global
use HaplotypeModule
use HaplotypeLibraryModule
implicit none

integer :: f,e,g,i,j,l,nCore,nGlobalLoop,CoreLength,CoreStart,CoreEnd,CompPhase
integer :: BanBoth(2),Ban(2),AnimalOn(ped%pedigreesize-ped%ndummys,2)
integer :: LoopStart,OffSet

integer,allocatable,dimension (:,:) :: LoopIndex
integer(kind=1),allocatable,dimension (:,:,:,:) :: Temp

type(HaplotypeLibrary) :: hapLib

type(Haplotype), dimension(2) :: workHap
type(Haplotype) :: tmpHap
type(Genotype) :: workGeno
integer, dimension(:),allocatable :: matches
integer :: id
integer :: phase

inputParams => defaultInput

! WARNING: This should go in a function since it is the same code as InternalParentPhaseElim subroutine
nGlobalLoop=25


! LoopeIndex is a matrix with two columns that will serve to define:
!   * 1.- nCores
!   * 2.- Core lengths
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

! WARNING: This can be better arrange with a ELSEIF statement and should go in a function since it
!          is the same code as InternalParentPhaseElim subroutine
! LoopStart indicates which is the first loop the algorithm should treat. The bigger the number of
! SNPs, the more the loops to be considered
if(inputParams%nsnp<=50) return
if(inputParams%nsnp>50) LoopStart=24
if(inputParams%nsnp>100) LoopStart=21
if(inputParams%nsnp>200) LoopStart=20
if(inputParams%nsnp>400) LoopStart=17
if(inputParams%nsnp>800) LoopStart=17
if(inputParams%nsnp>1000) LoopStart=15
if(inputParams%nsnp>1500) LoopStart=13
if(inputParams%nsnp>2000) LoopStart=10
if(inputParams%nsnp>3000) LoopStart=6
if(inputParams%nsnp>4000) LoopStart=3
if(inputParams%nsnp>5000) LoopStart=1

! Assumed that LoopIndex(:,1) are the numbers of cores for each phase step, LoopIndex):,2) are the core lengths
do i=LoopStart,nGlobalLoop
    LoopIndex(i,2)=int(float(inputParams%nsnp)/LoopIndex(i,1))
enddo

allocate(Temp(ped%pedigreeSize,inputParams%nsnp,2,2))
Temp=0
AnimalOn=0



! SIMULATE PHASING
! f is a variable to simulate shift or no-shift phasing
do f=1,2



    do l=LoopStart,nGlobalLoop



        ! Simulate phase without shift
        if (f==1) then
            nCore=inputParams%nsnp/LoopIndex(l,2)
            CoreStart=1
            CoreEnd=LoopIndex(l,2)

        ! Simulate phase with shift
        else
            OffSet=int(float(LoopIndex(l,2))/2)
            nCore=(inputParams%nsnp-(2*OffSet))/LoopIndex(l,2)
            CoreStart=1+Offset
            CoreEnd=LoopIndex(l,2)+Offset
        endif

        do g=1,nCore

            ! Make sure that cores ends correctly
            if ((f==1).and.(g==nCore)) CoreEnd=inputParams%nsnp
            if ((f==2).and.(g==nCore)) CoreEnd=inputParams%nsnp-OffSet

            ! Exit if the corelength is too small
            CoreLength=(CoreEnd-CoreStart)+1
            if (CoreLength<10) exit

            hapLib = newHaplotypeLibrary(nsnps=CoreLength, storeSize=500, stepSize=500)


            !$!OMP PARALLEL DO &
            !$!OMP DEFAULT(SHARED) &
            !$!OMP PRIVATE(i,e,tmphap) &
            !$!OMP FIRSTPRIVATE(compPhase)
            do i=1,ped%pedigreesize-ped%ndummys
                do e=1,2
                    ! WARNING: If GeneProbPhase has been executed, that is, if not considering the Sex Chromosome, then MSTermInfo={0,1}.
                    !          Else, if Sex Chromosome, then MSTermInfo is 0 always
                    !          So, if a Conservative imputation of haplotypes is selected, this DO statement will do nothing
                    if ((inputParams%ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle

                    ! Check if the haplotype for this core is completely phased
                    CompPhase=1
                    
                    tmpHap = ped%pedigree(i)%individualPhase(e)%subset(coreStart,coreEnd)
                    if (.not. tmpHap%fullyPhased()) CompPhase=0

                    ! If haplotype is completely phased, populate HapLib
                    ! NOTE: Since there is code in order to populate the Haplotype Library in
                    !       in AlphaPhase, it can be convenient to create a share procedure in
                    !       AlphaHouse
                    !$!OMP CRITICAL
                    if (CompPhase==1) then
                            id = hapLib%hasHap(tmphap)
                            if (id == 0) then 
                             ! If haplotype is not in the library, then
                             ! a new haplotype has been found, then populate the library
                                
                                id = hapLib%addHap(tmphap)
                            endif
                    endif
                      !$!OMP END CRITICAL
                enddo
            enddo
            !$!OMP END PARALLEL DO

            ! WARNING: This code does not match the corresponding code of the subroutine ImputeFromHDLibrary
            !          In ImputeFromHDLibrary, there are two steps, counting agreements and impute
            !          across candidate haplotypes, and counting agreements and impute across cores
            !          and phasing steps.

            !$OMP PARALLEL DO &
            !$OMP DEFAULT(SHARED) &
            !$OMP PRIVATE(i,j,e,BanBoth,matches,Ban,workGeno,phase,tmpHap) &
            !$OMP FIRSTPRIVATE(workHap)
            do i=1,ped%pedigreesize-ped%ndummys            
                BanBoth=0
                do e=1,2
                    tmpHap = ped%pedigree(i)%individualPhase(e)%subset(coreStart,coreEnd)

                    ! WARNING: If GeneProbPhase has been executed, that is, if not considering the Sex Chromosome, then MSTermInfo={0,1}.
                    !          Else, if Sex Chromosome, then MSTermInfo is 0 always
                    !          So, if a Conservative imputation of haplotypes is selected, this DO statement will do nothing
                    if ((inputParams%ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle

                    ! If haplotype is partially phased  

                    if (.not. tmpHap%allMissingOrError() .and. .not. tmpHap%fullyPhased()) then
                        ! Identify matching haplotypes in library]
                        ! THIS SHOULD PROBABLY USE A NEW FUNCTION IN HAPLIB (match) THAT DOESN'T TEST ERRORS BUT FOR NOW...
                        matches = hapLib%matchWithError(tmpHap,0)


                        ! If the number of candidate haplotypes is less than the 25% of the Library,
                        ! then impute if all alleles have been phased the same way
                        if (size(matches)<(float(hapLib%size)*0.25)) then
                            ! Ban this haplotype will be phased here and nowhere else
                            BanBoth(e)=1

                            workHap(e) = haplib%getConsensusHap(matches)
                        endif
                    endif
                enddo

                ! Has any haplotype been previously banned/phased?
                Ban=0
                if (BanBoth(1)==1) Ban(1)=1
                if (BanBoth(2)==1) Ban(2)=1
                ! If both gametes have been previously banned/phased, and
                ! if allele phases disagree with the genotype, then unban haplotypes
                if (sum(BanBoth(:))==2) then

                    workGeno = ped%pedigree(i)%individualGenotype%subset(corestart,coreEnd)
                    if (.not. workGeno%compatibleHaplotypes(workHap(1),workHap(2), 0)) then
                        Ban=0
                    endif
                endif

                ! Count the number of occurrences a particular phase is impute in a particular
                ! allele across the cores and across the internal phasing steps
                ! This implies occurrences across all the haplotype libraries of the different internal phasing steps
                do e=1,2
                    if ((inputParams%ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle
                    if (Ban(e)==1) then
                        AnimalOn(i,e)=1
                        do j=CoreStart,CoreEnd                                
                            phase = workHap(e)%getPhase((j- corestart)+1)
                            if (phase == 0) then
                                !$OMP atomic
                                Temp(i,j,e,1)=Temp(i,j,e,1)+1
                            else if (phase == 1) then
                                !$OMP atomic
                                Temp(i,j,e,2)=Temp(i,j,e,2)+1
                            endif
                        enddo
                    endif
                enddo
            enddo  
            !$OMP END PARALLEL DO

            ! Prepare the core for the next cycle
            CoreStart=CoreStart+LoopIndex(l,2)
            CoreEnd=CoreEnd+LoopIndex(l,2)
            call hapLib%destroyHaplotypeLibrary()
            if ((f==2).and.(g==nCore)) exit
        enddo
    
        
    enddo
enddo

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(i,j,e)
do i=1,ped%pedigreesize-ped%ndummys
    do e=1,2
        ! WARNING: If GeneProbPhase has been executed, that is, if not considering the Sex Chromosome, then MSTermInfo={0,1}.
        !          Else, if Sex Chromosome, then MSTermInfo is 0 always
        !          So, if a Conservative imputation of haplotypes is selected, this DO statement will do nothing
        if ((inputParams%ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle

        ! If all alleles across the cores and across the internal phasing steps have been phased the same way, impute
        if (AnimalOn(i,e)==1) then
            do j=1,inputParams%nsnp
                if (ped%pedigree(i)%individualPhase(e)%isMissing(j)) then
                    if ((Temp(i,j,e,1)>inputParams%nAgreeInternalHapLibElim).and.(Temp(i,j,e,2)==0))&
                        call ped%pedigree(i)%individualPhase(e)%setPhase(j,0)
                    if ((Temp(i,j,e,1)==0).and.(Temp(i,j,e,2)>inputParams%nAgreeInternalHapLibElim))&
                        call ped%pedigree(i)%individualPhase(e)%setPhase(j,1)
                endif
            enddo
        endif
    enddo
enddo
!$OMP END PARALLEL DO

! call hapLib%destroyHaplotypeLibrary()
deallocate(Temp)


end subroutine InternalHapLibImputationOld

            ! #############################################################################################################################################################################################################################


            SUBROUTINE PhaseElimination
                ! Candidate haplotype library imputation of alleles.
                ! For each core of each round of the LRPHLI algorithm, all haplotypes that have been found and
                ! stored in the haplotype library. Candidate haplotypes are restricted to the two haplotypes that
                ! have been identified for the individual by the LRPHLI algorithm for the true haplotype that an
                ! individual carries on its gametes. Within the core, all alleles that are known are compared to
                ! corresponding alleles in each of the individual haplotypes. The candidate haplotypes are checked
                ! for locations that have unanimous agreement about a particular allele. For alleles with complete
                ! agreement, a count of the suggested allele is incremented. Alleles are imputed if, at the end of
                ! passing across each core and each round of the LRPHLI algorithm, the count of whether the
                ! alleles are 0 or 1 is above a threshold in one direction and below a threshold in the other.
                ! This helps to prevent the use of phasing errors that originate from LRPHLI.
                ! This subroutine corresponds to Major sub-step 5 from Hickey et al., 2012 (Appendix A)

                use ISO_Fortran_Env
                use Global


                use AlphaImputeSpecFileModule
                use AlphaPhaseResultsModule
                implicit none

                integer :: e,g,i,j,GamA,GamB,PosHDInd
                integer :: StartSnp,EndSnp,Gam1,Gam2,AnimalOn(ped%pedigreeSize,2)
                integer(kind=1),allocatable,dimension (:,:,:,:) :: Temp
                integer :: unknownFreeIterator
                integer(kind=1) :: tmpHDPhase

                type(haplotype) :: tmpPhase(2)
                inputParams => defaultInput
                ! Number of animals that have been HD phased

                allocate(Temp(ped%pedigreeSize,inputParams%nsnpraw,2,2))

                Temp=0
                AnimalOn=0
                ! FOR EACH CORE OF EACH ROUND OF THE LRPHLI
                do unknownFreeIterator=1,apresults%nResults

                    startSnp = 1
                    EndSnp = 0

                    do g=1,size(apResults%results(unknownFreeIterator)%cores)
                        ! Initialize Start and End snps of the cores
                        StartSnp=apResults%results(unknownFreeIterator)%startIndexes(g)
                        EndSnp=apResults%results(unknownFreeIterator)%endIndexes(g)


                        !$!OMP PARALLEL DO &
                        !$!OMP DEFAULT(SHARED) &
                        !$!OMP PRIVATE(i,j,e,Gam1,Gam2,GamA,GamB,tmpPhase,tmpHDPhase,PosHDInd)
                        do i=1,ped%nHd
                            ! Look for possible gametes through the Haplotype
                            ! Library constructed during the phasing step
                            PosHDInd=ped%hdMap(i)         ! Index of the individual in the HD phase information
                            
                            tmpPhase(1) = ped%pedigree(posHDInd)%individualPhase(1)%subset(startSnp,endSnp)
                            tmpPhase(2) = ped%pedigree(posHDInd)%individualPhase(2)%subset(startSnp,endSnp)
                            ! If there is one allele phased at least
                               if (.not. (tmpPhase(1)%fullyPhased() .and. tmpPhase(2)%fullyPhased()))  then
                                    ! If at least one locus is heterozygous
                                    if (.not. tmpPhase(1)%equalHap(tmpPhase(2))) then
                                        Gam1=0
                                        Gam2=0
                                        do e=1,2
                                            GamA=1
                                            GamB=1

                                            if (tmpPhase(e)%mismatches(apResults%results(unknownFreeIterator)%cores(g)%phase(i,1)) /= 0) then
                                                gamA = 0
                                            endif
                                            if (tmpPhase(e)%mismatches(apResults%results(unknownFreeIterator)%cores(g)%phase(i,2)) /= 0) then
                                                gamB = 0
                                            endif


                                        


                                            ! Paternal haplotype (gamete) is strictly the same as my paternal haplotype from the Hap Library
                                            if ((e==1).and.(GamA==1).and.(GamB==0)) Gam1=1
                                            ! Paternal haplotype (gamete) is strictly the same as my maternal haplotype from the Hap Library
                                            if ((e==1).and.(GamA==0).and.(GamB==1)) Gam1=2
                                            ! Maternal haplotype (gamete) is strictly the same as my paternal haplotype from the Hap Library
                                            if ((e==2).and.(GamA==1).and.(GamB==0)) Gam2=1
                                            ! Maternal haplotype (gamete) is strictly the same as my maternal haplotype from the Hap Library
                                            if ((e==2).and.(GamA==0).and.(GamB==1)) Gam2=2

                                            ! Basically the important thing is that haplotype e is present in the Haplotype library. It is not
                                            ! important which haplotype it is, whether the paternal or the maternal.
                                        enddo

                                        ! If the paternal and maternal gametes are different
                                        if (Gam1/=Gam2) then
                                            AnimalOn(i,:)=1             ! Consider this animal in further steps

                                            ! Paternal gamete is in the Hap Library
                                            if (Gam1/=0) then
                                                do j=StartSnp,EndSnp
                                                    ! Count the number of alleles coded with 0 and 1
                                                    if (ped%pedigree(posHdind)%individualPhase(1)%isMissing(j)) then
                                                        tmpHDPhase = apResults%results(unknownFreeIterator)%cores(g)%phase(i,gam1)%getPhase(j-StartSnp+1)
                                                        if(tmpHdPhase==0) then
                                                            !$OMP ATOMIC
                                                            Temp(PosHDInd,j,1,1)=Temp(PosHDInd,j,1,1)+1
                                                        else if(tmpHdPhase==1) then
                                                            !$OMP ATOMIC
                                                            Temp(PosHDInd,j,1,2)=Temp(PosHDInd,j,1,2)+1
                                                        endif
                                                    endif
                                                enddo
                                            endif
                                            ! Maternal gamete is in the Hap Library
                                            if (Gam2/=0) then
                                                do j=StartSnp,EndSnp
                                                    ! Count the number of alleles coded with 0 and 1
                                                        
                                                    if (ped%pedigree(posHdind)%individualPhase(2)%isMissing(j)) then
                                                        tmpHDPhase = apResults%results(unknownFreeIterator)%cores(g)%phase(i,gam2)%getPhase(j-StartSnp+1)
                                                        if(tmpHdPhase==0) then
                                                            !$OMP ATOMIC
                                                            Temp(PosHDInd,j,2,1)=Temp(PosHDInd,j,2,1)+1
                                                        else if(tmpHdPhase==1) then
                                                            !$OMP ATOMIC
                                                            Temp(PosHDInd,j,2,2)=Temp(PosHDInd,j,2,2)+1
                                                        endif
                                                    endif
                                                enddo
                                            endif
                                        endif
                                    endif
                                endif

                        enddo
                        !$!OMP END PARALLEL DO

                    enddo
                enddo

                do e=1,2
                    do j=1,inputParams%nsnp
                        do i=1,ped%pedigreeSize- ped%nDummys
                            if (AnimalOn(i,e)==1) then
                                if (ped%pedigree(i)%individualPhase(e)%isMissing(j)) then
                                    ! Impute phase allele with the most significant code for that allele across haplotypes
                                    ! only if the other codification never happens
                                    if ((Temp(i,j,e,1)>inputParams%nAgreeInternalHapLibElim).and.(Temp(i,j,e,2)==0)) then
                                        call ped%pedigree(i)%individualPhase(e)%setPhase(j,0)
                                    end if
                                    if ((Temp(i,j,e,1)==0).and.(Temp(i,j,e,2)>inputParams%nAgreeInternalHapLibElim)) then
                                        call ped%pedigree(i)%individualPhase(e)%setPhase(j,1)
                                    end if
                                endif
                            endif
                        enddo
                    enddo
                enddo
                deallocate(Temp)

            END SUBROUTINE PhaseElimination

            !#############################################################################################################################################################################################################################

            SUBROUTINE ParentPhaseElimination
                ! Imputation of single alleles based on parental phase.
                ! For each core of each round of the LRPHLI algorithm, all haplotypes that have been found and
                ! stored in the haplotype library. Candidate haplotypes are restricted to the two haplotypes that
                ! have been identified for each of the individual parents with high-density genotype information by
                ! the LRPHLI algorithm for the true haplotype that an individual carries on its gametes. Within the
                ! core, all alleles that are known are compared to corresponding alleles in each of the individual
                ! haplotypes. The candidate haplotypes are checked for locations that have unanimous agreement about
                ! a particular allele. For alleles with complete agreement, a count of the suggested allele is
                ! incremented. Alleles are imputed if, at the end of passing across each core and each round of the
                ! LRPHLI algorithm, the count of whether the alleles are 0 or 1 is above a threshold in one
                ! direction and below a threshold in the other. This helps to prevent the use of phasing errors that
                ! originate from LRPHLI.
                ! This subroutine corresponds to Major sub-step 4 from Hickey et al., 2012 (Appendix A)

                use Global

                use AlphaImputeSpecFileModule
                implicit none

                integer :: e,g,h,i,j,PedId,GamA,GamB,PosHDInd
                integer :: StartSnp,EndSnp,AnimalOn(ped%pedigreeSize,2)
                integer(kind=1),allocatable,dimension (:,:,:,:) :: Temp

                integer :: tempCount, tmpPhase


                inputParams => defaultInput
                ped%nHd=(count(Setter(:)==1))

                allocate(Temp(ped%pedigreeSize,inputParams%nsnpraw,2,2))
                Temp=0
                AnimalOn=0

                do h=1,apResults%nResults



                    ! Get phase information

                    do g=1,size(apResults%results(h)%cores)
                        ! Initialize Start and End snps of the cores
                        StartSnp= apResults%results(h)%startIndexes(g)
                        EndSnp=apResults%results(h)%endIndexes(g)                        
                        
                        block
                            use individualModule
                            type(individual) ,pointer :: parent
                            type(Haplotype) :: tmpHap

                            !$OMP PARALLEL DO &
                            !$OMP DEFAULT(SHARED) &
                            !$OMP PRIVATE(i,j,e,PedId,PosHDInd,GamA,GamB,parent,tmpHap,TempCount,tmpPhase)
                            do i=1,ped%pedigreeSize- ped%nDummys

                                do e=1,2
                                    PedId=e+1
                                    if (ped%pedigree(i)%isDummyBasedOnIndex(pedId)) cycle !checked this one makes sense
                                    parent => ped%pedigree(i)%getSireDamObjectByIndex(pedId)

                                    ! parent => ped%pedigree((ped%pedigree(i)%getSireDamNewIDByIndex(pedID)))

                                    ! Skip if, in the case of sex chromosome, me and my parent are heterogametic
                                   if ((inputParams%SexOpt==1 .and.ped%pedigree(i)%gender==inputParams%HetGameticStatus .and.(parent%gender==inputParams%HetGameticStatus))) cycle

                                    ! We look for gamete through those individuals that have parents with HD genotype information
                                    if (associated(parent)) then
                                        ! We look for possible gametes within the haplotypes identified to each of the individual's parents constructed during the phasing step
                                        posHdInd = ped%hdDictionary%getValue(parent%originalId)
                                        tmpHap = ped%pedigree(i)%individualPhase(e)%subset(startsnp,endSnp)
                                        ! If there is one allele phased at least
                                        if ((.not. tmpHap%fullyPhased()).and.(PosHDInd>0)) then
                                            GamA=1
                                            GamB=1
                                            TempCount=0

                                            if (apResults%results(h)%cores(g)%phase(PosHDInd,1)%mismatches(tmpHap) >= ImputeFromParentCountThresh) then
                                                gamA=0
                                            endif
                                            if (apResults%results(h)%cores(g)%phase(PosHDInd,2)%mismatches(tmpHap) >= ImputeFromParentCountThresh) then
                                                GamB=0
                                            endif
                                            ! NOTE: [..."and the candidate haplotypes for each individual's gametes are restricted
                                            !       to the two haplotypes that have been identified for each of its parents..."]
                                            !   If GamA==0, then the haplotype is different from individual's paternal haplotype
                                            !   If GamB==0, then the haplotype is different from individual's maternal haplotype

                                            ! e is the individual's paternal haplotype
                                            if ((GamA==1).and.(GamB==0)) then
                                                ! Consider the haplotype of this animal in further steps
                                                AnimalOn(i,e)=1
                                                do j=StartSnp,EndSnp
                                                    ! Count the number of alleles coded with 0 and 1
                                                   
                                                    if ( ped%pedigree(i)%individualPhase(e)%ismissing(j)) then
                                                        tmpPhase = apResults%results(h)%cores(g)%phase(PosHDInd,1)%getPhase(j-startsnp+1)
                                                        if(tmpPhase==0) then
                                                            !$OMP ATOMIC
                                                            Temp(i,j,e,1)=Temp(i,j,e,1)+1
                                                        else if(tmpPhase==1) then
                                                            !$OMP ATOMIC
                                                            Temp(i,j,e,2)=Temp(i,j,e,2)+1
                                                        endif
                                                    endif
                                                enddo
                                            endif

                                            ! e is the individual's maternal haplotype
                                            if ((GamA==0).and.(GamB==1)) then
                                                ! Consider the haplotype of this animal in further steps
                                                AnimalOn(i,e)=1
                                                do j=StartSnp,EndSnp
                                                    ! Count the number of alleles coded with 0 and 1
                                                    if ( ped%pedigree(i)%individualPhase(e)%ismissing(j)) then
                                                        tmpPhase = apResults%results(h)%cores(g)%phase(PosHDInd,1)%getPhase(j-startSnp+1)
                                                        if(tmpPhase==0) then
                                                            !$OMP ATOMIC
                                                            Temp(i,j,e,1)=Temp(i,j,e,1)+1
                                                        else if(tmpPhase==1) then
                                                            !$OMP ATOMIC
                                                            Temp(i,j,e,2)=Temp(i,j,e,2)+1
                                                        end if
                                                    endif
                                                enddo
                                            endif
                                        endif
                                    endif
                                enddo
                            enddo
                            !$OMP END PARALLEL DO
                        end block
                    enddo
                enddo


                do i=1,ped%pedigreeSize- ped%nDummys
                    do e=1,2
                    
                        
                            if (AnimalOn(i,e)==1) then
                                if ((inputParams%sexopt==0).or.(ped%pedigree(i)%gender==inputParams%HomGameticStatus)) then
                                    do j=1,inputParams%nsnp
                                        if (ped%pedigree(i)%individualPhase(e)%isMissing(j)) then
                                            ! Impute phase allele with the most significant code for that allele across haplotypes
                                            ! only if the other codification never happens
                                            if ((Temp(i,j,e,1)>inputParams%nAgreeInternalHapLibElim).and.(Temp(i,j,e,2)==0)) then
                                                call ped%pedigree(i)%individualPhase(e)%setPhase(j,0)
                                            end if
                                            if ((Temp(i,j,e,1)==0).and.(Temp(i,j,e,2)>inputParams%nAgreeInternalHapLibElim)) then
                                                call ped%pedigree(i)%individualPhase(e)%setPhase(j,1)
                                            end if
                                        end if
                                    end do
                                end if
                            end if
                        
                    end do
                
                        if ((inputParams%sexopt==1).and.(ped%pedigree(i)%gender==inputParams%HetGameticStatus)) then
                            if (AnimalOn(i,inputParams%HomGameticStatus)==1) then
                                do j=1,inputParams%nsnp
                                    if (ped%pedigree(i)%individualPhase(inputParams%HomGameticStatus)%isMissing(j)) then
                                        ! Impute phase allele with the most significant code for that allele across haplotypes
                                        ! only if the other codification never happens
                                        if ((Temp(i,j,inputParams%HomGameticStatus,1)>inputParams%nAgreeInternalHapLibElim).and.&
                                            (Temp(i,j,inputParams%HomGameticStatus,2)==0)) then
                                            call ped%pedigree(i)%individualPhase(1)%setPhase(j,0)
                                            call ped%pedigree(i)%individualPhase(2)%setPhase(j,0)
                                        end if
                                        if ((Temp(i,j,inputParams%HomGameticStatus,1)==0).and.&
                                            (Temp(i,j,inputParams%HomGameticStatus,2)>inputParams%nAgreeInternalHapLibElim)) then
                                            call ped%pedigree(i)%individualPhase(1)%setPhase(j,1)
                                            call ped%pedigree(i)%individualPhase(2)%setPhase(j,1)
                                        end if
                                    endif
                                end do
                            end if
                        end if
                    
                enddo
                deallocate(Temp)

            END SUBROUTINE ParentPhaseElimination

            !#############################################################################################################################################################################################################################

            SUBROUTINE ImputeFromHDLibrary
                ! Candidate haplotype library imputation of alleles.
                ! For each core of each round of the LRPHLI algorithm, all haplotypes that have been found and
                ! stored in the haplotype library are initially considered to be candidates for the true haplotype
                ! that an individual carries on its gametes. Within the core, all alleles that are known are
                ! compared to corresponding alleles in each of the haplotypes in the library. Haplotypes that have a
                ! number of disagreements greater than a small error threshold have their candidacy rejected. At the
                ! end of this loop, the surviving candidate haplotypes are checked for locations that have unanimous
                ! agreement about a particular allele. For alleles with complete agreement, a count of the suggested
                ! allele is incremented. Alleles are imputed if, at the end of passing across each core and each
                ! round of the LRPHLI algorithm, the count of whether the alleles are 0 or 1 is above a threshold in
                ! one direction and below a threshold in the other. This helps to prevent the use of phasing errors
                ! that originate from LRPHLI.
                ! This subroutine corresponds to Major sub-step 3 from Hickey et al., 2012 (Appendix A)

                use Global

                use HaplotypeModule
                use AlphaImputeSpecFileModule
                use AlphaPhaseResultsModule
                implicit none


                integer :: i,j,k,h,e,f,g,CoreLength,nHap,TempCount
                integer :: StartSnp,EndSnp,Counter,BanBoth(2),Ban(2),AnimalOn(ped%pedigreeSize,2)
                logical :: PatMatDone(2)
                integer,allocatable,dimension (:,:,:,:) :: Temp

                integer :: phase

                type(Haplotype), dimension(:),allocatable :: workHap
                type(Haplotype) :: tmpHap
                type(Genotype) :: workGeno
                integer, dimension(:),allocatable :: matches

                inputParams => defaultInput
                ! Temp(ped%pedigreeSize, inputParams%nsnps, PatHap, Phase)
                allocate(Temp(0:ped%pedigreeSize,inputParams%nsnpraw,2,2))
                Temp=0

                AnimalOn=0

                do h=1,apResults%nResults
                    ! Get HIGH DENSITY phase information of this phasing step and information
                    ! of core indexes

                    ! Get core information of number of cores and allocate start and end cores information
                    do g=1,size(apResults%results(h)%cores)
                        ! Initialize Start and End snps of the cores
                        StartSnp=apResults%results(h)%startIndexes(g)
                        EndSnp=apResults%results(h)%endIndexes(g)

                        CoreLength=(EndSnp-StartSnp)+1


                        nHap = apResults%results(h)%libraries(g)%size
                        coreLength = apResults%results(h)%libraries(g)%nsnps



                        ! Binary haplib is apResults%results(h)%libraries(g)%newStore(z)

                        !$OMP PARALLEL DO &
                        !$OMP DEFAULT(SHARED) &
                        !$OMP PRIVATE(i,e,f,j,k,TempCount,Counter,PatMatDone,BanBoth,Ban,phase,tmpHap, matches,workGeno, workHap)
                        do i=1,ped%pedigreeSize- ped%nDummys
                            ! The number of candidate haplotypes is the total number of haps in the library times
                            ! 2 (Paternal and maternal candidates)
                            allocate(workHap(2))

                            PatMatDone=.false.

                            BanBoth=0

                            ! For each haplotype
                            do e=1,2
                                ! If GeneProbPhase has been executed, that is, if not considering the Sex Chromosome, then MSTermInfo={0,1}.
                                ! Else, if Sex Chromosome, then MSTermInfo is 0 always
                                ! So, if a Conservative imputation of haplotypes is selected, this DO statement will do nothing
                                if ((inputParams%ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle
                                    tmpHap = ped%pedigree(i)%individualPhase(e)%subset(startSnp,endSnp)
                                    !  If haplotype is partially phased
                                     if (.not. tmpHap%allMissingOrError()) then
                                        ! TODO - this can probably be removed - and this should be made a logical
                                        PatMatDone(e) = .true.
                                        if (.not. tmpHap%fullyPhased()) then
                                                                               
                                            matches = apResults%results(h)%libraries(g)%matchWithError(tmpHap,ImputeFromHDLibraryCountThresh)

                        
                                        ! CountAB(inputParams%nsnps,0:1) matrix indicating how many haplotypes has been phased as 0
                                            ! or 1 in a particular allele

                                            ! If the number of candidate haplotypes is less than the 25% of the Library,
                                            ! then impute if all alleles have been phased the same way
                                            if (float(size(matches))<(float(nHap)*0.25)) then
                                                ! Ban this haplotype will be phased here and nowhere else
                                                BanBoth(e)=1

                                                ! Count the occurrences in phasing of alleles across candidate haplotypes
                                                WorkHap(e) = apResults%results(h)%libraries(g)%getConsensusHap(matches)
                                            endif
                                        endif
                                    endif
                            enddo

                            ! If one of the haplotypes is partially phased
                            if (any(PatMatDone == .true.)) then
                                Ban=0
                                ! Any haplotype has been previously banned/phased?
                                if (BanBoth(1)==1) Ban(1)=1
                                if (BanBoth(2)==1) Ban(2)=1

                                ! If both gametes have been previously banned/phased,
                                ! check whether the phase given agrees with genotype
                                if (sum(BanBoth(:))==2) then
                                    workGeno = ped%pedigree(i)%individualGenotype%subset(StartSnp,endSnp)
                                    if (.not. workGeno%compatibleHaplotypes(workHap(1),workHap(2), 0)) then
                                        Ban=0
                                    endif
                                endif

                                ! Count the number of occurrences a phase is impute in a particular
                                ! allele across the cores across the internal phasing steps
                                ! This implies occurrences across all the haplotype libraries of the different
                                ! internal phasi#ng steps
                                do e=1,2
                                    ! If GeneProbPhase has been executed, that is, if not considering the Sex Chromosome, then MSTermInfo={0,1}.
                                    ! Else, if Sex Chromosome, then MSTermInfo is 0 always
                                    ! So, if a Conservative imputation of haplotypes is selected, this DO statement will do nothing
                                    if ((inputParams%ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle
                                    if (Ban(e)==1) then
                                        AnimalOn(i,e)=1
                                        do j=StartSnp,EndSnp
                                            phase = workHap(e)%getPhase((j- startSnp)+1)
                                            if (phase==0) Temp(i,j,e,1)=Temp(i,j,e,1)+1
                                            if (phase==1) Temp(i,j,e,2)=Temp(i,j,e,2)+1
                                        enddo
                                    endif
                                enddo
                            endif
                            deallocate(workHap)
                        enddo
                        !$OMP END PARALLEL DO
                        
                    enddo

                enddo


                do j=1,inputParams%nsnp
                    ! If GeneProbPhase has been executed, that is, if not considering the Sex Chromosome, then MSTermInfo={0,1}.
                    ! Else, if Sex Chromosome, then MSTermInfo is 0 always
                    ! So, if a Conservative imputation of haplotypes is selected, this DO statement will do nothing

                    do i=1,ped%pedigreeSize- ped%nDummys
                        do e=1,2
                            if (ped%pedigree(i)%isDummyBasedOnIndex(e)) cycle
                            if (AnimalOn(i,e)==1) then
                                if ((inputParams%ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle
                                ! If all alleles across the cores across the internal phasing steps have been phased the
                                ! same way, impute

                                if ( ped%pedigree(i)%individualPhase(e)%ismissing(j)) then
                                    if ((Temp(i,j,e,1)>inputParams%nAgreeInternalHapLibElim).and.(Temp(i,j,e,2)==0)) then
                                         call ped%pedigree(i)%individualPhase(e)%setPhase(j,0)
                                    end if
                                    if ((Temp(i,j,e,1)==0).and.(Temp(i,j,e,2)>inputParams%nAgreeInternalHapLibElim)) then
                                        call ped%pedigree(i)%individualPhase(e)%setPhase(j,1)
                                    end if
                                endif
                            end if
                        end do
                    enddo
                enddo
                deallocate(Temp)


            END SUBROUTINE ImputeFromHDLibrary

            !#############################################################################################################################################################################################################################

            SUBROUTINE BaseAnimalFillIn
                ! Impute phase for Base Animals (without any pedigree information).
                ! The internal phasing has been calculated by means of AlphaPhase using both 'with' and 'without' shift parameter.
                ! The function select the middle phase run without shifting the cores and the middle core within this phase run, and then
                ! selects the complementary phase run with the shift cores. Shifted and non shifted cores allows overlapping and so the
                ! function goes through SNPs from left to right and from right to left from the MiddleCore jumping to shifted and non shifted
                ! cores.
                ! The function ends when all the SNPs have been phased. If the haplotypes have not been fully phased means that there is some
                ! recombination in that haplotype.

                use Global
                use ISO_Fortran_Env

                use AlphaImputeSpecFileModule
                use AlphaPhaseResultsModule

                implicit none

                integer :: e,h,i,g,MiddleResult,MiddleResultNoShift,CoreLength
                integer :: CompJump,StartSnp,EndSnp,UptoRightSnp,UptoLeftSnp,UpToCoreA,UpToCoreB,Recmb,CompLength,RL
                integer :: UpToSnp,StPt,EndPt,FillInSt,FillInEnd,hdAnimid
                ! integer,allocatable,dimension (:,:) :: CoreIndexA,CoreIndexB,AnimRecomb
                integer,allocatable,dimension (:,:) :: AnimRecomb
                integer :: middleCoreIndex,middleCoreIndexShift
                type(individual), pointer :: tmpAnim
                logical :: C1,C2,C3,C4
                inputParams => defaultInput

                ! TODO This should be removed
                ! ped%nHd=(count(Setter(:)==1))

                allocate(AnimRecomb(ped%pedigreeSize,2))
                AnimRecomb=0

                ! Get core information of number of cores and allocate start and end cores information

                ! Select the core in the middle
                MiddleResult=apresults%nResults/4
                if (MiddleResult==0) MiddleResult=1

                CompJump = apresults%nResults/2
                if (CompJump==0) CompJump=1
                MiddleResultNoShift = MiddleResult + CompJump

                ! Get HIGH DENSITY phase information of this phasing step
                ! WARNING: If I only want to phase base animals, why do I need to read the whole file?
               

                ! Impute HD phase of the middle core of the middle phasing step
                ! WARNING: Why to impute phase information only for this case?
                middleCoreIndex = apresults%results(MiddleResult)%nCores/2
                if (middleCoreIndex == 0) then
                    middleCoreIndex = 1
                endif

                middleCoreIndexShift = apresults%results(MiddleResultNoShift)%nCores/2
                if (middleCoreIndexShift == 0) then
                    middleCoreIndexShift = 1
                endif

                StartSnp=apresults%results(MiddleResult)%startIndexes(middleCoreIndex)
                EndSnp=apresults%results(MiddleResult)%endIndexes(middleCoreIndex)
                CoreLength=(EndSnp-StartSnp)+1
                ! do i=1,ped%pedigreeSize- ped%nDummys
                do i=1,ped%nHd
                    ! If I have no parents and if I am somebody

                    tmpAnim => ped%pedigree(ped%hdMap(i))
                    if (tmpAnim%founder) then
                        block
                        type(Haplotype) :: tmpHap
                        tmpHap = tmpAnim%individualPhase(1)%subset(startSnp,endSnp)
                        ! Check if the two haplotypes are equal
                        if (tmpHap%equalHap(tmpAnim%individualPhase(2)%subset(startSnp,EndSnp))) then
                            call tmpAnim%individualPhase(1)%setSubset(apResults%results(MiddleResult)%cores(middleCoreIndex)%phase(i,1),startSnp)
                            call tmpAnim%individualPhase(2)%setSubset(apResults%results(MiddleResult)%cores(middleCoreIndex)%phase(i,2),startSnp)
                        end if
                        end block
                        
                    endif
                enddo

                UpToRightSnp=EndSnp
                UpToLeftSnp=StartSnp
                UpToCoreA=middleCoreIndex

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
                        if ((apResults%results(MiddleResult)%nCores==1)&
                            .AND.(apResults%results(MiddleResultNoShift)%nCores==1)) exit ! If the number of cores is 1, EXIT
                        ! This will force the subroutine to finish
                        ! since it will be exit from both DO statements
                        h=h+1
                        if (mod(h,2)/=0) then                   ! If ODD
                            do g=1,apresults%results(MiddleResultNoShift)%nCores
                                if ((apresults%results(MiddleResultNoShift)%startIndexes(g)<UptoSnp)&
                                    .AND.(apresults%results(MiddleResultNoShift)%endIndexes(g)>UptoSnp)) then

                                    UpToCoreB=g
                                    exit
                                endif
                            enddo
                            if (e==1) then
                                StartSnp=apresults%results(MiddleResultNoShift)%startIndexes(UpToCoreB)
                                EndSnp=apresults%results(MiddleResultNoShift)%endIndexes(UpToCoreB)

                                StPt=StartSnp
                                EndPt=UpToSnp
                                FillInSt=StartSnp
                                FillInEnd=EndSnp
                            else
                                StartSnp=apresults%results(MiddleResultNoShift)%endIndexes(UpToCoreB)
                                EndSnp=apresults%results(MiddleResultNoShift)%startIndexes(UpToCoreB)
                                StPt=UpToSnp
                                EndPt=StartSnp
                                FillInSt=EndSnp
                                FillInEnd=StartSnp
                            endif

                            CompLength=abs(UpToSnp-StartSnp)+1

                            do i=1,ped%pedigreeSize- ped%nDummys
                                hdAnimid = ped%hdDIctionary%getValue(ped%pedigree(i)%originalId)
                                if (ped%pedigree(i)%founder .and.(hdAnimid/=DICT_NULL).and.(AnimRecomb(i,RL)==0)) then

                                    block 

                                    type(haplotype) :: phase1,phase2, phase1Hd,phase2HD,tmp1,tmp2

                                    phase1= ped%pedigree(i)%individualPhase(1)%subset(stpt,endpt)
                                    phase2= ped%pedigree(i)%individualPhase(2)%subset(stpt,endpt)
                                    tmp1 = apResults%results(MiddleResultNoShift)%getfullphasesinglehap(hdAnimid,1)
                                    tmp2 = apResults%results(MiddleResultNoShift)%getfullphasesinglehap(hdAnimid,2)
                                    phase1HD = tmp1%subset(stpt,endpt)
                                    phase2HD = tmp2%subset(stpt,endpt)
                                    
                                    C1 = phase1HD%equalHap(phase1)
                                    C2 = phase1HD%equalHap(phase2)
                                    C3 = phase2HD%equalHap(phase1)
                                    C4 = phase2HD%equalHap(phase2)

                                    Recmb=1
                                    ! If one haplotype is the same as the paternal, impute
                                        if (C1 .and. .not. C3) then                                        
                                            call ped%pedigree(i)%individualPhase(1)%setSubset(phase1HD,FillInSt)
                                            Recmb=0

                                        ! If one haplotype is the same as the paternal, impute
                                        else if (.not. c1 .and. c3) then
                                            call ped%pedigree(i)%individualPhase(1)%setSubset(phase2HD,FillInSt)
                                            Recmb=0
                                        endif

                                        ! If one haplotype is the same as the maternal, impute
                                        if (c2 .and. .not. c4) then
                                            call ped%pedigree(i)%individualPhase(2)%setSubset(phase1HD,FillInSt)
                                            Recmb=0

                                        ! If one haplotype is the same as the maternal, impute
                                        else if (.not. c2 .and. c4) then
                                            call ped%pedigree(i)%individualPhase(2)%setSubset(phase2HD,FillInSt)
                                            Recmb=0
                                        endif

                                        ! There is bridge in the phase in the RL direction. Nothing more will be done
                                        if (Recmb==1) then
                                            AnimRecomb(i,RL)=1
                                        end if
                                    end block
                                endif
                                
                            enddo
                            if (RL == 1) then
                                UpToSnp=apresults%results(MiddleResultNoShift)%startIndexes(UpToCoreB)
                            else
                                UpToSnp=apresults%results(MiddleResultNoShift)%endIndexes(UpToCoreB)
                            end if
                        else                                    ! if EVEN
                            do g=1,apresults%results(MiddleResult)%nCores
                                if ((apresults%results(MiddleResult)%startIndexes(g)<UptoSnp)&
                                    .AND.(apresults%results(MiddleResult)%endIndexes(g))>UptoSnp) then
                                    UpToCoreA=g
                                    exit
                                endif
                            enddo

                            if (e==1) then
                                startSnp = apresults%results(MiddleResult)%startIndexes(UpToCoreA)
                                endSnp = apresults%results(MiddleResult)%endIndexes(UpToCoreA)
                                StPt=StartSnp
                                EndPt=UpToSnp
                                FillInSt=StartSnp
                                FillInEnd=EndSnp
                            else
                                startSnp = apresults%results(MiddleResult)%endIndexes(UpToCoreA)
                                endSnp = apresults%results(MiddleResult)%startIndexes(UpToCoreA)
                                StPt=UpToSnp
                                EndPt=StartSnp
                                FillInSt=EndSnp
                                FillInEnd=StartSnp
                            endif
                            CompLength=abs(UpToSnp-StartSnp)+1
                            do i=1,ped%pedigreeSize- ped%nDummys

                                hdAnimid = ped%hdDIctionary%getValue(ped%pedigree(i)%originalID)
                                if ( ped%pedigree(i)%founder .and.(hdAnimid/= DICT_NULL).and.(AnimRecomb(i,RL)==0)) then
                                   
                                    block 

                                    type(haplotype) :: phase1,phase2, phase1Hd,phase2HD,tmp1,tmp2
                                    phase1= ped%pedigree(i)%individualPhase(1)%subset(stpt,endpt)
                                    phase2= ped%pedigree(i)%individualPhase(2)%subset(stpt,endpt)

                                      tmp1 = apResults%results(MiddleResult)%getfullphasesinglehap(hdAnimid,1)
                                    tmp2 = apResults%results(MiddleResult)%getfullphasesinglehap(hdAnimid,2)
                                    phase1HD = tmp1%subset(stpt,endpt)
                                    phase2HD = tmp2%subset(stpt,endpt)

                                    C1 = phase1HD%equalHap(phase1)
                                    C2 = phase1HD%equalHap(phase2)
                                    C3 = phase2HD%equalHap(phase1)
                                    C4 = phase2HD%equalHap(phase2)

                                    Recmb=1
                                    ! If one haplotype is the same as the paternal, impute
                                        if (C1 .and. .not. C3) then                                        
                                            call ped%pedigree(i)%individualPhase(1)%setSubset(phase1HD,FillInSt)
                                            Recmb=0

                                        ! If one haplotype is the same as the paternal, impute
                                        else if (.not. c1 .and. c3) then
                                            call ped%pedigree(i)%individualPhase(1)%setSubset(phase2HD,FillInSt)
                                            Recmb=0
                                        endif

                                        ! If one haplotype is the same as the maternal, impute
                                        if (c2 .and. .not. c4) then
                                            call ped%pedigree(i)%individualPhase(2)%setSubset(phase1HD,FillInSt)
                                            Recmb=0

                                        ! If one haplotype is the same as the maternal, impute
                                        else if (.not. c2 .and. c4) then
                                            call ped%pedigree(i)%individualPhase(2)%setSubset(phase1HD,FillInSt)
                                            Recmb=0
                                        endif

                                        ! There is bridge in the phase in the RL direction. Nothing more will be done
                                        if (Recmb==1) then
                                            AnimRecomb(i,RL)=1
                                        end if
                                    end block
                                    if (Recmb==1) then
                                        AnimRecomb(i,RL)=1
                                    end if
                                endif
                            enddo

                            if (RL == 1) then
                                UpToSnp= apresults%results(MiddleResult)%startIndexes(UpToCoreA)

                            else
                                UpToSnp= apresults%results(MiddleResult)%endIndexes(UpToCoreA)
                            end if
                        endif

                        ! Exit condition
                        if (e==1) then
                            if (UpToSnp>=inputParams%nsnp) exit     ! Exit if we've reached the last SNP
                        else
                            if (UpToSnp<=1) exit        ! Exit if we've reached the first SNP
                        endif
                    enddo
                enddo

            END SUBROUTINE BaseAnimalFillIn

            !#############################################################################################################################################################################################################################

            subroutine InitialiseArrays
                ! Impute phase information for homozygous cases
                use AlphaImputeSpecFileModule
                use global, only :ped
                implicit none

                integer :: i
                inputParams => defaultInput

                do i=1, ped%pedigreeSize
                    call ped%pedigree(i)%individualGenotype%setHaplotypeFromGenotype(ped%pedigree(i)%individualPhase(1))
                    call ped%pedigree(i)%individualGenotype%setHaplotypeFromGenotype(ped%pedigree(i)%individualPhase(2))
                enddo

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
                ! if (TestVersion==1) call CurrentYield
                ! if (TestVersion==1) call Checker


            end subroutine GeneralFillIn




            subroutine GeneralFillInInit
                ! This function implements the four Minor sub-steps explained in Hickey et al. (2012; Appendix A)
                use Global

                implicit none

                call ParentHomoFill                     ! Minor sub-step 1. Parent Homozygous fill in
                

                    call ped%writeOutGenotypes("Results/" //"afters1")
                call PhaseComplement                    ! Minor sub-step 2. Phase Complement
                    call ped%writeOutGenotypes("Results/" //"afters2")
                call ImputeParentByProgenyComplement    ! Minor sub-step 3. Impute Parents from Progeny Complement
                    call ped%writeOutGenotypes("Results/" //"afters3")
                call MakeGenotype                       ! Minor sub-step 4. Make Genotype
                    call ped%writeOutGenotypes("Results/" //"afters4")
                ! if (TestVersion==1) call CurrentYield
                ! if (TestVersion==1) call Checker


            end subroutine GeneralFillInInit


            !#############################################################################################################################################################################################################################

            SUBROUTINE EnsureHetGametic
                ! Impute phase to Y chromosome from X chromosome for heterogametic individuals
                use Global
                use AlphaImputeSpecFileModule
                implicit none

                integer :: i

                inputParams => defaultInput

                    do i=1,ped%pedigreeSize- ped%nDummys
                        if (ped%pedigree(i)%gender==inputParams%HetGameticStatus) then
                            call ped%pedigree(i)%individualPhase(1)%setFromOtherIfMissing(ped%pedigree(i)%individualPhase(2))
                            call ped%pedigree(i)%individualPhase(2)%setFromOtherIfMissing(ped%pedigree(i)%individualPhase(1))
                        endif
                    enddo
                

            END SUBROUTINE EnsureHetGametic

            !#############################################################################################################################################################################################################################

            SUBROUTINE MakeGenotype
                ! Any individual that has a missing genotype information but has both alleles
                ! known, has its genotype filled in as the sum of the two alleles
                use Global
                use AlphaImputeSpecFileModule

                implicit none
                integer :: i, j

                integer :: phase1,phase2

                do i=1,ped%pedigreeSize-ped%nDummys

                    do j=1, inputParams%nsnpRaw
                    ! call ped%pedigree(i)%IndividualGenotype%setFromHaplotypesIfMissing(ped%pedigree(i)%individualPhase(1),ped%pedigree(i)%individualPhase(2))

                        if (ped%pedigree(i)%individualGenotype%isMissing(j)) then
                            phase1 = ped%pedigree(i)%individualPhase(1)%getPhase(j)
                            phase2 = ped%pedigree(i)%individualPhase(2)%getPhase(j)
                            if (phase1 /= 9 .and. phase2 /= 9) then
                                call ped%pedigree(i)%individualGenotype%setGenotype(j,(phase1 + phase2))
                            endif
                        endif
                    enddo
                enddo 


            END SUBROUTINE MakeGenotype

            !#############################################################################################################################################################################################################################

            subroutine PhaseComplement
                ! If the genotype at a locus for an individual is known and one of its alleles has been determined
                ! then impute the missing allele as the complement of the genotype and the known phased allele
                use Global
                implicit none

                integer :: i
                type(haplotype) :: comp1, comp2

                do i=1,ped%pedigreeSize- ped%nDummys  
                    comp2 = ped%pedigree(i)%individualGenotype%complement(ped%pedigree(i)%individualPhase(1))
                    comp1 = ped%pedigree(i)%individualGenotype%complement(ped%pedigree(i)%individualPhase(2))

                    call ped%pedigree(i)%individualPhase(1)%setErrorToMissing()
                    call ped%pedigree(i)%individualPhase(2)%setErrorToMissing()
                enddo

            end subroutine PhaseComplement

            !#############################################################################################################################################################################################################################

            subroutine ParentHomoFill
                ! Fill in the allele of an offspring of a parent that has both its
                ! alleles filled in and has a resulting genotype that is homozygous
                use Global
                use AlphaImputeSpecFileModule
                implicit none

                integer :: e,i,ParId
                type(Genotype) :: tmpGeno

                inputParams => defaultInput
                do i=1,ped%pedigreeSize- ped%nDummys
                    if (inputParams%sexopt==0 .or. (inputParams%sexopt==1 .and. ped%pedigree(i)%gender/=inputParams%HetGameticStatus) ) then     ! If individual is homogametic
                        do e=1,2
                            ParId=ped%pedigree(i)%getSireDamNewIDByIndex(e+1)
                            if (parId == 0) cycle
                            tmpGeno = newGenotypeHap(ped%pedigree(parId)%IndividualPhase(1),ped%pedigree(parId)%individualPhase(2))
                                call tmpGeno%setHaplotypeFromGenotypeIfMissing(ped%pedigree(i)%individualPhase(e))                                
                        enddo
                    else
                        ParId= ped%pedigree(i)%getSireDamNewIDByIndex(inputParams%HomGameticStatus+1) !the homogametic parent
                        if (parId == 0) cycle
                        tmpGeno = newGenotypeHap(ped%pedigree(parId)%IndividualPhase(1),ped%pedigree(parId)%individualPhase(2))
                        call tmpGeno%setHaplotypeFromGenotypeIfMissing(ped%pedigree(i)%individualPhase(1)) 
                        call tmpGeno%setHaplotypeFromGenotypeIfMissing(ped%pedigree(i)%individualPhase(2)) 
                    endif
                enddo

            end subroutine ParentHomoFill

            !#############################################################################################################################################################################################################################

            subroutine ImputeParentByProgenyComplement
                ! If one of the parental allele is known and the other missing, the fill
                ! in the missing allele in the parent if, at least, one of its offspring
                ! is known to carry an allele that does not match the known allele in the parent
                use Global
                use AlphaImputeSpecFileModule
                use individualModule
                implicit none

                integer :: i,k,l,Count1,Count0
                integer :: sireDam
                type(individual), pointer :: tmpChild
                integer(kind=1) :: phase1,phase2, childphase
                inputParams => defaultInput
                ! TODO maybe change this to pedigreeSIze? 
                do i=1,ped%pedigreeSize- ped%nDummys
                    if (ped%pedigree(i)%nOffs /= 0) then       ! check that animal i,j is a sire or a dam

                        ! Sex chromosome
                        if (inputParams%sexopt==1) then
                            do k=1,inputParams%nsnp
                                phase1 = ped%pedigree(i)%individualPhase(1)%getPhase(k)
                                phase2 = ped%pedigree(i)%individualPhase(2)%getPhase(k)
                                ! Mat gamete missing -> fill if offspring suggest heterozygous
                                ! WARNING: This was comment the other way around in the original version of the code
                                if ((phase1/=9).and.(phase2==9)) then
                                    Count1=0
                                    Count0=0
                                    if (phase1==1) Count1=1
                                    if (phase1==0) Count0=1

                                    ! Look for the individual progeny and count their phase
                                    do l=1,ped%pedigree(i)%nOffs

                                        tmpChild => ped%pedigree(i)%offsprings(l)%p
                                        childPhase = ped%pedigree(tmpChild%id)%individualPhase(sireDam)%getPhase(k)
                                        ! This is the only difference with the inputParams%sexopt=0 code below. Duplicating
                                        ! the code can be avoided by including the IF statement here instead than
                                        ! outside the SNPs loop.
                                        if ((ped%pedigree(i)%gender ==inputParams%HetGameticStatus).and.(tmpChild%gender==inputParams%HetGameticStatus)) cycle
                                        if (tmpChild%sirePointer == ped%pedigree(i)) then
                                                sireDam = 1
                                            else
                                                sireDam = 2
                                            endif
                                        if (childPhase==0) Count0=Count0+1
                                        if (childPhase==1) Count1=Count1+1

                                        if ((Count0>0).and.(Count1>0)) then                 !Consider increasing the number of offspring required to all ow for genotyping error
                                            if (phase1==0) then
                                                call ped%pedigree(i)%individualPhase(2)%setPhase(k,1)
                                            else if (phase1==1) then
                                                call ped%pedigree(i)%individualPhase(2)%setPhase(k,0)
                                            endif
                                        endif
                                    enddo
                                endif

                                !Pat gamete missing -> fill if offspring suggest heterozygous
                                ! WARNING: This comment was the other way around in the original version of the code
                                if ((phase2/=9).and.(phase2==9)) then
                                    Count1=0
                                    Count0=0
                                    if (phase2==1) Count1=1
                                    if (phase2==0) Count0=1

                                    do l=1,ped%pedigree(i)%nOffs

                                        tmpChild => ped%pedigree(i)%offsprings(l)%p
                                        childPhase = tmpChild%individualPhase(sireDam)%getPhase(k)
                                        ! This is the only difference with the inputParams%sexopt=0 code below. Duplicating
                                        ! the code can be avoided by including the IF statement here instead than
                                        ! outside the SNPs loop.
                                        if ((ped%pedigree(i)%gender ==inputParams%HetGameticStatus).and.(tmpChild%gender==inputParams%HetGameticStatus)) cycle
                                        
                                          if (tmpChild%sirePointer == ped%pedigree(i)) then
                                                sireDam = 1
                                            else
                                                sireDam = 2
                                            endif
                                        
                                        if (childPhase==0) Count0=Count0+1
                                        if (childPhase==1) Count1=Count1+1

                                        if ((Count0>0).and.(Count1>0)) then                 !Consider increasing the number of offspring required to all ow for genotyping error
                                            if (phase2==0) then
                                                call ped%pedigree(i)%individualPhase(1)%setPhase(k,1)
                                            else if (phase2==1) then
                                                call ped%pedigree(i)%individualPhase(1)%setPhase(k,0)
                                            endif
                                        endif
                                    enddo
                                endif
                            enddo

                            ! Generic chromosome
                        else
                            
                            do k=1,inputParams%nsnp
                                phase1 = ped%pedigree(i)%individualPhase(1)%getPhase(k)
                                phase2 = ped%pedigree(i)%individualPhase(2)%getPhase(k)
                                if ((phase1/=9).and.(phase2==9)) then               !Pat gamete missing fill if offspring suggest heterozygous
                                    Count1=0
                                    Count0=0
                                    

                                        if (phase1==1) then
                                            Count1=1
                                        else if (phase1==0) then 
                                            count0=1
                                        endif
                                    

                                    do l=1,ped%pedigree(i)%nOffs
                                        tmpChild => ped%pedigree(i)%offsprings(l)%p
                                        
                                        if (tmpChild%sirePointer == ped%pedigree(i)) then
                                            sireDam = 1
                                        else
                                            sireDam = 2
                                        endif
                                        childPhase = tmpChild%individualPhase(sireDam)%getPhase(k)
                                        if (childPhase==0) Count0=Count0+1
                                        if (childPhase==1) Count1=Count1+1
                                        if ((Count0>0).and.(Count1>0)) then                 !Consider increasing the number of offspring required to all ow for genotyping error
                                            if (phase1==0) then
                                                call ped%pedigree(i)%individualPhase(2)%setPhase(k,1)
                                            else if (phase2==1) then
                                                call ped%pedigree(i)%individualPhase(2)%setPhase(k,0)
                                            endif
                                        endif
                                    enddo
                                
                                else if ((phase2/=9).and.(phase1==9)) then               !Mat gamete missing fill if offspring suggest heterozygous
                                    Count1=0
                                    Count0=0
                                    if (phase2==1) Count1=1
                                    if (phase2==0) Count0=1
                                    do l=1,ped%pedigree(i)%nOffs
                                        tmpChild => ped%pedigree(i)%offsprings(l)%p
                                        
                                        if (tmpChild%sirePointer == ped%pedigree(i)) then
                                            sireDam = 1
                                        else
                                            sireDam = 2
                                        endif
                                        childPhase = ped%pedigree(tmpChild%id)%individualPhase(sireDam)%getPhase(k)

                                        if (childPhase==0) Count0=Count0+1
                                        if (childPhase==1) Count1=Count1+1

                                        if ((Count0>0).and.(Count1>0)) then                 !Consider increasing the number of offspring required to all ow for genotyping error
                                            if (phase2==0) then
                                                call ped%pedigree(i)%individualPhase(1)%setPhase(k,1)
                                            else if (phase2==1) then
                                                call ped%pedigree(i)%individualPhase(1)%setPhase(k,0)
                                            endif
                                        endif
                                    enddo
                                endif
                            enddo
                        endif
                    endif
                enddo

    

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
                ! TODOphase this all needs redone
                use AlphaImputeSpecFileModule
                use AlphaImputeInputOutputModule, only: readProbabilitiesFull, readProbabilitiesFullCluster

                implicit none

                integer :: i,j,k,m,StSnp,EnSnp
                real :: PatAlleleProb(inputParams%nsnp,2),MatAlleleProb(inputParams%nsnp,2),HetProb(inputParams%nsnp)
                integer :: Informativeness(inputParams%nsnp,6),TmpInfor(inputParams%nsnp,6),GrandPar
                integer :: phase1, phase2 
                inputParams => defaultInput


                if (inputParams%BypassGeneProb==0) then
                    ! Get information from GeneProb


                    if (inputParams%restartOption /= OPT_RESTART_ALL) then
                        if (inputParams%cluster) then
                            call readProbabilitiesFullCluster(GenosProbs,ped%pedigreeSize-ped%nDummys, inputParams%nsnp, inputParams, GpIndex)
                        else
                            call readProbabilitiesFull("./GeneProb/GenotypeProbabilities.txt",GenosProbs,ped%pedigreeSize-ped%nDummys, inputParams%nsnp)
                        endif
                    endif 
                    StSnp=1
                    EnSnp=inputParams%nSnp

                    do i=1,ped%pedigreeSize- ped%nDummys

                        PatAlleleProb(StSnp:EnSnp,2)=GenosProbs(i,StSnp:EnSnp,3)+GenosProbs(i,StSnp:EnSnp,4)
                        MatAlleleProb(StSnp:EnSnp,1)=GenosProbs(i,StSnp:EnSnp,1)+GenosProbs(i,StSnp:EnSnp,3)
                        MatAlleleProb(StSnp:EnSnp,2)=GenosProbs(i,StSnp:EnSnp,2)+GenosProbs(i,StSnp:EnSnp,4)

                        ! Probability of heterozygosity
                        HetProb(StSnp:EnSnp)=GenosProbs(i,StSnp:EnSnp,2)+GenosProbs(i,StSnp:EnSnp,3)

                        do j=StSnp,EnSnp
                            if (PatAlleleProb(j,1)>=GeneProbThresh) then
                                call ped%pedigree(i)%individualPhase(1)%setPhase(j,0)
                            else if (PatAlleleProb(j,2)>=GeneProbThresh) then
                                call ped%pedigree(i)%individualPhase(1)%setPhase(j,1)
                            endif

                            if (MatAlleleProb(j,1)>=GeneProbThresh) then
                                call ped%pedigree(i)%individualPhase(2)%setPhase(j,0)
                            else if (MatAlleleProb(j,2)>=GeneProbThresh) then 
                                call ped%pedigree(i)%individualPhase(2)%setPhase(j,1)
                            endif

                            if (HetProb(j)>=GeneProbThresh) call ped%pedigree(i)%individualGenotype%setGenotype(j,1)
                        enddo
                    enddo
                endif


                open(unit=102,file="." // DASH // "Miscellaneous" // "IndividualSnpInformativeness.txt", status="unknown")
                open(unit=103,file="." // DASH // "Miscellaneous" // "IndividualMendelianInformativeness.txt", status="unknown")

                allocate(GlobalTmpCountInf(ped%pedigreeSize,8))
                allocate(MSTermInfo(ped%pedigreeSize,2))

                MSTermInfo=0
                do i=1,ped%pedigreeSize- ped%nDummys
                    if (ped%pedigree(i)%Genotyped) MSTermInfo(i,:)=1
                    TmpInfor(:,:)=-99
                    GlobalTmpCountInf(i,:)=0
                    Informativeness(:,:)=9 ! What the hell is this variable for??
                    j=0

                    ! Check whether my parents and grandparents are heterozygous
                    do m=1,inputParams%nsnpRaw
                        if (SnpIncluded(m)==1) then                     ! Whether to consider this SNP
                            j=j+1                                       ! Number of SNPs included so far

                            if (ped%pedigree(i)%individualGenotype%getgenotype(j)==1) then               ! If heterozygous
                                phase1 =  ped%pedigree(i)%individualPhase(1)%getPhase(j)
                                phase2 =  ped%pedigree(i)%individualPhase(2)%getPhase(j)
                                ! TODO - now looking at dummy animals as parents - why not? 
                                if (ped%pedigree(i)%getSireDamNewIDByIndex(2) /= 0) then
                                    ! My father is heterozygous

                                    
                                    if (ped%pedigree(ped%pedigree(i)%getSireDamNewIDByIndex(2))%individualGenotype%getGenotype(j)==1) then
                                        ! And have my father haplotype phased

                                        
                                        if ((phase1==0).or.(phase1==1)) then
                                            Informativeness(j,1)=1
                                            GlobalTmpCountInf(i,1)=GlobalTmpCountInf(i,1)+1     ! Count the number of SNPs phased of my father haplotype
                                            TmpInfor(GlobalTmpCountInf(i,1),1)=j                ! The SNP no. GlobalTmpCountInf(i,1) is my SNP no. j
                                        endif
                                    endif
                                endif

                                ! My mother is heterozygous
                                if (ped%pedigree(i)%getSireDamNewIDByIndex(3) /= 0) then
                                    if (ped%pedigree(ped%pedigree(i)%getSireDamNewIDByIndex(3))%individualGenotype%getGenotype(j)==1) then
                                        ! And have my mother haplotype phased
                                        if ((phase2==0).or.(phase2==1)) then
                                            Informativeness(j,2)=1
                                            GlobalTmpCountInf(i,2)=GlobalTmpCountInf(i,2)+1
                                            TmpInfor(GlobalTmpCountInf(i,2),2)=j
                                        endif
                                    endif
                                endif


                                ! My father haplotype is phased
                                if ((phase1==0).or.(phase1==1)) then
                                    ! If my paternal GranSire is heterozygous
                                    GrandPar=ped%pedigree(i)%getPaternalGrandSireRecodedIndexNoDummy()
                                    if (ped%pedigree(grandPar)%individualGenotype%getgenotype(j)==1) then
                                        Informativeness(j,3)=1
                                        GlobalTmpCountInf(i,3)=GlobalTmpCountInf(i,3)+1
                                        TmpInfor(GlobalTmpCountInf(i,3),3)=j
                                    endif
                                    ! If my maternal GranDam is heterozygous
                                    GrandPar=ped%pedigree(i)%getPaternalGrandSireRecodedIndexNoDummy()
                                    if (ped%pedigree(grandPar)%individualGenotype%getgenotype(j)==1) then
                                        Informativeness(j,4)=1
                                        GlobalTmpCountInf(i,4)=GlobalTmpCountInf(i,4)+1
                                        TmpInfor(GlobalTmpCountInf(i,4),4)=j
                                    endif
                                endif

                                ! My mother haplotype is phased
                                if ((phase2==0).or.(phase2==1)) then
                                    ! If my maternal GranSire is heterozygous
                                    ! if (ped%pedigree(i)%hasDummyParentsOrGranparents()) cycle
                                    ! TODO make this return 0
                                    GrandPar= ped%pedigree(i)%getPaternalGrandSireRecodedIndexNoDummy()

                                    if (ped%pedigree(grandPar)%individualGenotype%getgenotype(j)==1) then
                                        Informativeness(j,5)=1
                                        GlobalTmpCountInf(i,5)=GlobalTmpCountInf(i,5)+1
                                        TmpInfor(GlobalTmpCountInf(i,5),5)=j
                                    endif
                                    ! If my maternal GranDam is heterozygous
                                    GrandPar=ped%pedigree(i)%getmaternalGrandDamRecodedIndexNoDummy()
                                    if (ped%pedigree(grandPar)%individualGenotype%getgenotype(j)==1) then
                                        Informativeness(j,6)=1
                                        GlobalTmpCountInf(i,6)=GlobalTmpCountInf(i,6)+1
                                        TmpInfor(GlobalTmpCountInf(i,6),6)=j
                                    endif
                                endif
                            endif
                        endif
                    enddo

                    GlobalTmpCountInf(i,7)=  ped%pedigree(i)%individualPhase(1)%numberNotMissingOrError()         ! Count the number of genotyped allele in the paternal haplotype
                    GlobalTmpCountInf(i,8)=  ped%pedigree(i)%individualPhase(2)%numberNotMissingOrError()         ! Count the number of genotyped allele in the maternal haplotype
                    if (GlobalTmpCountInf(i,1)>0) MSTermInfo(i,1)=1             ! If the paternal haplotype is phased
                    if (GlobalTmpCountInf(i,2)>0) MSTermInfo(i,2)=1             ! If the maternal haplotype is phased
                    do k=1,6
                        write(102,'(a20,i2,i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7)') ped%pedigree(i)%originalID, Setter(i), GlobalTmpCountInf(i,k), TmpInfor(1:GlobalTmpCountInf(i,k),k)
                    enddo
                    write(102,'(a20,i2,i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7)') ped%pedigree(i)%originalID, Setter(i), GlobalTmpCountInf(i,7)
                    write(102,'(a20,i2,i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7)') ped%pedigree(i)%originalID, Setter(i), GlobalTmpCountInf(i,8)
                    write(103,'(a20,2i10)') ped%pedigree(i)%id, MSTermInfo(i,:)
                enddo
                close(102)
            end subroutine GeneProbPhase

            !#############################################################################################################################################################################################################################

            subroutine ManageWorkLeftRight
                ! Major sub-step 8 is iterated a number of times with increasingly relaxed
                ! restrictions. After each iteration, the minor sub-steps are also carried out.
                use Global
                use AlphaImputeSpecFileModule
                implicit none

                inputParams => defaultInput
                MaxLeftRightSwitch=4; MinSpan=200
                call WorkLeftRight


                if (inputParams%sexopt==1) then
                    call EnsureHetGametic
                end if
                call GeneralFillIn

                MaxLeftRightSwitch=3; MinSpan=200
                call WorkLeftRight
                if (inputParams%sexopt==1) then
                    call EnsureHetGametic
                end if
                call GeneralFillIn

                MaxLeftRightSwitch=4; MinSpan=100
                call WorkLeftRight
                if (inputParams%sexopt==1) then
                    call EnsureHetGametic
                end if
                call GeneralFillIn

                MaxLeftRightSwitch=4; MinSpan=50
                call WorkLeftRight
                if (inputParams%sexopt==1) then
                    call EnsureHetGametic
                end if
                call GeneralFillIn

            end subroutine ManageWorkLeftRight

            function getLoopStart(snps) Result(LoopStart)
                integer, intent(in) :: snps
                integer :: LoopStart

                select case (snps)
                case (:50)
                    LoopStart = 0
                case (51:100)
                    LoopStart = 24
                case (101:200)
                    LoopStart = 21
                case (201:400)
                    LoopStart = 20
                case (401:800)
                    LoopStart = 17
                case (801:1000)
                    LoopStart = 17
                case (1001:1500)
                    LoopStart = 15
                case (1501:2000)
                    LoopStart = 13
                case (2001:3000)
                    LoopStart = 10
                case (3001:4000)
                    LoopStart = 6
                case (4001:5000)
                    LoopStart = 3
                case (5001:)
                    LoopStart = 1
                end select
            end function getLoopStart

            subroutine setLoopIndex(LoopIndex)
                integer, dimension(:,:), intent(out) :: LoopIndex

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
            end subroutine setLoopIndex

            !#############################################################################################################################################################################################################################

            subroutine RestrictedWorkLeftRight
                ! Imputation based on identifying where recombination occurs during inheritance from parent to offspring.
                ! Each gamete of an individual is examined from the beginning to the end and from the end to the
                ! beginning of the chromosome. In each direction, at loci where both the individual and its parent are
                ! heterozygous and have phase information resolved, this information is used to determine which of the
                ! parental gametes the individual received. Loci for which this cannot be determined but which are
                ! between two loci that (a) can be determined and (b) come from the same parental gamete, are assumed
                ! to come from this gamete (i.e. no double recombination event in between). Alleles are imputed in the
                ! individual when analysis in both directions of the chromosome has identified the same inherited gamete
                ! and when the parent is phased for this locus in the suggested gamete, subject to the restrictions that
                ! he number of recombination events for the individuals is less than a threshold and that the region in
                ! which two recombination events occurred exceeds a threshold length.
                ! This subroutine corresponds to WorkLeftRight (Major sub-step 8 from Hickey et al., 2012 (Appendix A))
                ! but with two main restrictions:
                !   * The number of recombinations is fixed to MaxLeftRightSwitch=4; and
                !     the threshold lenght for recombination is fixed to MinSpan=200
                !   * The number of unphased alleles has to be lower than a threshold ((2*inputParams%nsnp)*0.07))

                use Global
                use AlphaImputeSpecFileModule
                implicit none

                integer :: e,i,j,HetEnd,HetStart,RSide,LSide,PatMat,SireDamRL
                integer :: CountRightSwitch,CountLeftSwitch,StartPt,EndPt,PedId

                integer, dimension(:), allocatable :: WorkRight, WorkLeft

                integer :: StartDis,EndDis,StartJ,k
                integer,allocatable,dimension(:) :: TempVec
                real,allocatable,dimension(:) :: LengthVec
                type(AlphaImputeInput), pointer :: inputParams
                integer(kind=1) :: phase(2),tmpPhase


                inputParams => defaultInput

                allocate(WorkRight(inputParams%nsnp))
                allocate(WorkLeft(inputParams%nsnp))
                allocate(TempVec(inputParams%nsnp))
                allocate(LengthVec(inputParams%nsnp))

                MaxLeftRightSwitch=4; MinSpan=200

                do i=1,ped%pedigreeSize- ped%nDummys
                    HetEnd=-1
                    HetStart=-1
                    ! For each gamete
                    do e=1,2
                        PatMat=e
                        SireDamRL=e+1
                        CountLeftSwitch=0
                        CountRightSwitch=0
                        block
                            integer :: tmpGender
                            PedId=ped%pedigree(i)%getSireDamNewIDByIndex(SireDamRL)



                            if (PedId /= 0) then
                                tmpGender = ped%pedigree(PedId)%gender
                                ! if (ped%pedigree(pedId)%isDummy) then
                                !     cycle
                                ! endif
                            else
                                cycle
                            endif
                            ! Skip if, in the case of sex chromosome, me and my parent are heterogametic
                            if ((inputParams%SexOpt==1).and.(ped%pedigree(i)%gender==inputParams%hetGameticStatus).and.(tmpGender==inputParams%hetGameticStatus)) cycle
                        endblock
                        !! SCAN HAPLOTYPE IN TWO DIRECTIONS: L->R AND R->L
                        ! If not a base animal and the number of unphased alleles is lower than a threshold
                        ! WARNING: WHAT IS THIS THRESHOLD?
                        
                        if (((float(ped%pedigree(PedId)%individualPhase(1)%numberMissing()+ped%pedigree(PedId)%individualPhase(2)%numberMissing())/(2*inputParams%nsnp))<0.07)) then           !(RecIdHDIndex(PedId)==1)
                            WorkRight=9
                            RSide=9

                            ! Go throught haplotype from Left to Right
                            ! finding the first heterozygous allele of this parent, and...
                            do j=1,inputParams%nsnp
                                phase(1) = ped%pedigree(PedId)%individualPhase(1)%getPhase(j)
                                phase(2) = ped%pedigree(pedId)%individualPhase(2)%getPhase(j)
                                if ((phase(1)/=phase(2)).and.&
                                    (phase(1)/=9).and.(phase(2)/=9))  then
                                    HetStart=j
                                    tmpPhase = ped%pedigree(i)%individualPhase(PatMat)%getPhase(HetStart)
                                    ! Check if this allele corresponds to my parent paternal haplotype
                                    if (tmpPhase==phase(1)) then
                                        WorkRight(HetStart)=1   ! HetStart allele corresponds to Pat Haplotype in the LR direction
                                        RSide=1                 ! We are actually in the Paternal haplotype for the LR direction
                                        exit
                                    endif
                                    ! Check if this allele corresponds to my parent maternal haplotype
                                    if (tmpPhase==phase(2)) then
                                        WorkRight(HetStart)=2   ! HetStart allele corresponds to Mat Haplotype in the LR direction
                                        RSide=2                 ! We are actually in the Maternal haplotype for the LR direction
                                        exit
                                    endif
                                endif
                            enddo

                            ! ... Identifying recombinations
                            if (RSide/=9) then
                                do j=HetStart+1,inputParams%nsnp
                                    ! If this allele has different phased as the current haplotype, then
                                    ! Change haplotype and increase the number of recombinations of this haplotype
                                    phase(1) = ped%pedigree(pedId)%individualPhase(RSide)%getPhase(j)
                                    tmpPhase = ped%pedigree(i)%individualPhase(patMat)%getPhase(j)
                                    if ((tmpPhase/=phase(1) ).and.&
                                        (phase(1) /=9).and.(tmpPhase/=9)) then
                                        RSide=abs((RSide-1)-1)+1
                                        CountRightSwitch=CountRightSwitch+1
                                    endif
                                    ! Wich paternal gamete the individual has received
                                    WorkRight(j)=RSide
                                enddo
                            endif

                            WorkLeft=9
                            LSide=9

                            ! Go through haplotype from Right to Left
                            ! finding the first heterozygous allele of this parent, and...
                            do j=inputParams%nsnp,1,-1
                                phase(1) = ped%pedigree(PedId)%individualPhase(1)%getPhase(j)
                                phase(2) = ped%pedigree(pedId)%individualPhase(2)%getPhase(j)
                                if ((phase(1) /=phase(2)) .and.&
                                    (phase(1) /=9).and.(phase(2) /=9))  then
                                    HetEnd=j
                                    tmpPhase = ped%pedigree(i)%individualPhase(patMat)%getPhase(HetEnd)
                                    if (tmpPhase==phase(1)) then
                                        WorkRight(HetEnd)=1
                                        LSide=1
                                        exit
                                    endif
                                    if (tmpPhase==phase(2)) then
                                        WorkRight(HetEnd)=2
                                        LSide=2
                                        exit
                                    endif
                                endif
                            enddo

                            ! ... Identifying recombinations
                            if (LSide/=9) then
                                do j=HetEnd-1,1,-1

                                    phase(1) = ped%pedigree(i)%individualPhase(patMat)%getPhase(j)
                                    phase(2) = ped%pedigree(pedId)%individualPhase(LSide)%getPhase(j)
                                    if ((phase(1) /=phase(2)).and.&
                                        (phase(2)/=9).and.(phase(1) /=9) ) then
                                        LSide=abs((LSide-1)-1)+1
                                        CountLeftSwitch=CountLeftSwitch+1
                                    endif
                                    ! Wich paternal gamete the individual has received
                                    WorkLeft(j)=LSide
                                enddo
                            endif

                            !$$$$$$$$$$$$$$$$$$$

                            ! UNPHASE THOSE ALLELES WITH SOME AMBIGUITY DUE TO RECOMBINATION
                            ! Let's be (StartDis:EndDis) the SNPs of the two direction disagree
                            StartDis=-9
                            EndDis=-9
                            TempVec=9
                            LengthVec=0.0
                            StartJ=1
                            do j=StartJ,inputParams%nsnp
                                ! Initalize variables StartDis and EndDis
                                ! StartDis is the first allele where different directions differ
                                if (StartDis==-9) then
                                    if (abs(WorkLeft(j)-WorkRight(j))==1) then
                                        StartDis=j
                                        TempVec(StartDis)=1
                                    endif
                                endif
                                ! EndDis is the last allele where different directions agree.
                                !   (StartDis/=-9) guarantees that EndDis > StartDis
                                !   (EndDis==-9) guarantees that EndDis is not updated
                                if ((WorkLeft(j)==WorkRight(j)).and.(WorkLeft(j)/=9).and.(StartDis/=-9).and.(EndDis==-9)) then
                                    EndDis=j-1
                                    TempVec(EndDis)=2
                                endif

                                ! Move StartDis to the first phased allele (from left) that comes from a heterozygous case
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

                                ! Move EndDis to the last phased allele (from left) that comes from a heterozygous case
                                if (EndDis/=-9) then
                                    do k=EndDis,inputParams%nsnp
                                        if ((GlobalWorkPhase(PedId,k,1)+GlobalWorkPhase(PedId,k,2))==1) then
                                            if (GlobalWorkPhase(i,k,e)/=9) then
                                                exit
                                            endif
                                        endif
                                    enddo
                                    EndDis=k
                                    if (EndDis>inputParams%nsnp) EndDis=inputParams%nsnp
                                    TempVec(EndDis)=2
                                endif

                                ! If StartDis==9 means that there haplotype is the same from Left to Right than from Right to Left
                                !   * In this case, EndDis==9 too, then DO NOTHING
                                ! If EndDis==9 means that
                                !   * StartDis==9 or
                                !   * WorkLeft(j)==9, which means that
                                !       - The whole haplotype is homozygous
                                !       - parent is not completely phased for that allele
                                if ((StartDis/=-9).and.(EndDis/=-9)) then
                                    ! WARNING: 3 is the only value that is used for TempVec
                                    TempVec(StartDis+1:EndDis-1)=3
                                    LengthVec(StartDis+1:EndDis-1)=1.0/(((EndDis-1)-(StartDis+1))+1)
                                    StartJ=EndDis+1
                                    StartDis=-9
                                    EndDis=-9
                                endif
                            enddo

                            ! Remove phase and genotype for those alleles with no explanation due to heterozygosity and recombination
                            do j=1,inputParams%nSnp
                                if ( TempVec(j)==3) then
                                    phase(1) = ped%pedigree(pedId)%individualPhase(1)%getPhase(j)
                                    phase(2)  = ped%pedigree(pedId)%individualPhase(2)%getPhase(j)
                                    if (phase(1) /=phase(2) .AND. &
                                    (phase(1) /=9) .and. phase(2) /=9) then

                                        call ped%pedigree(i)%individualPhase(e)%setPhase(j,9)
                                    if (Setter(i)/=1) then ! Skip HD individuals
                                        call ped%pedigree(i)%individualGenotype%setGenotype(j,9)
                                    end if
                                endif

                                endif
                            enddo
                            !$$$$$$$$$$$$$$$$$$$

                            !! IMPUTE PHASE WHETHER IT IS POSSIBLE
                            ! WARNING: From Hickey et al. 2012 (Appendix A):
                            !          ["Alleles are imputed ... subject to the restriction that the number of recombinations events for the individuals is less than a threshold, AND
                            !          that the region in which two recombination events occurred exceeds a threshold lenght."]
                            !          What it is coded is ["... than a threshold, OR that the region..."]
                            ! The number of recombinations in total (LR + RL) is less than a threshold
                            if ((CountLeftSwitch+CountRightSwitch)<(2*MaxLeftRightSwitch)) then
                                do j=1,inputParams%nsnp
                                    if (ped%pedigree(i)%individualPhase(Patmat)%isMissing(j)) then

                                        ! WARNING: This can be coded in a conciser way
                                        ! Phase if the allele in one of the two directions is missing
                                        if ((WorkLeft(j)==9).and.(WorkRight(j)/=9)) then
                                            call ped%pedigree(i)%individualPhase(patMat)%setPhase(j, ped%pedigree(pedId)%individualPhase(WorkRight(j))%getPhase(j))
                                        else if ((WorkLeft(j)/=9).and.(WorkRight(j)==9)) then

                                            call ped%pedigree(i)%individualPhase(patMat)%setPhase(j, ped%pedigree(pedId)%individualPhase(WorkLeft(j))%getPhase(j))

                                        ! Phase if alleles is the two directions agree
                                        else if ((WorkLeft(j)/=9).and.(WorkRight(j)==WorkLeft(j))) then
                                            call ped%pedigree(i)%individualPhase(patMat)%setPhase(j, ped%pedigree(pedId)%individualPhase(WorkLeft(j))%getPhase(j))
                                        endif
                                    endif
                                enddo
                            else
                                ! Let's be (StartPt:EndPt) the SNPs in the two direction agree
                                EndPt=0
                                StartPt=0
                                ! TODO add openmp
                                do while ((StartPt<(inputParams%nsnp-MinSpan)).and.(EndPt<(inputParams%nsnp-MinSpan)))      ! If EndPt >(inputParams%nsnp-MinSpan), then recombination events does not exceed the threshold MinSpan
                                    do j=EndPt+1,inputParams%nsnp
                                        if ((WorkLeft(j)/=9).and.(WorkRight(j)==WorkLeft(j)))  then
                                            StartPt=j
                                            exit
                                        endif
                                        if (j==inputParams%nsnp) StartPt=j
                                    enddo
                                    do j=StartPt,inputParams%nsnp
                                        if ((WorkLeft(j)==9).or.(WorkRight(j)/=WorkLeft(j)))  then
                                            EndPt=j
                                            exit
                                        endif
                                        if (j==inputParams%nsnp) EndPt=j
                                    enddo
                                    ! The region in which two recombination events occurred exceeds a threshold lenght
                                    if (((EndPt-StartPt)+1)>MinSpan) then
                                        do j=StartPt,EndPt
                                            ! WARNING: This condition is supposed to be meet always since SNPs in (StartPt:EndPt)
                                            !          meet the condition (WorkLeft(j)/=9).and.(WorkRight(j)==WorkLeft(j))
                                            tmpPhase = ped%pedigree(pedId)%individualPhase(workRight(j))%getPhase(j)
                                            if (tmpPhase/=9) then
                                                call ped%pedigree(i)%individualPhase(patMat)%setPhase(j,tmpPhase)
                                            endif
                                        enddo
                                    endif
                                enddo
                            endif
                        endif
                    enddo
                enddo

                ! Impute phase for the Heterogametic chromosome from the Homogametic one, which has been already phased
                do i=1,ped%pedigreeSize- ped%nDummys
                    if ((inputParams%SexOpt==1).and.(ped%pedigree(i)%gender==inputParams%hetGameticStatus)) then
                         ped%pedigree(i)%individualPhase(inputParams%hetGameticStatus) = ped%pedigree(i)%individualPhase(inputParams%HomGameticStatus)
                        
                    endif
                enddo

                GlobalWorkPhase=ped%getPhaseAsArray()

                deallocate(WorkRight)
                deallocate(WorkLeft)
                deallocate(TempVec)
                deallocate(LengthVec)

            end subroutine RestrictedWorkLeftRight


            !#############################################################################################################################################################################################################################

            subroutine WorkLeftRight
                ! Imputation based on identifying where recombination occurs during inheritance from parent to offspring.
                ! Each gamete of an individual is examined from the beginning to the end and from the end to the
                ! beginning of the chromosome. In each direction, at loci where both the individual and its parent
                ! are heterozygous and have phase information resolved, this information is used to determine which
                ! of the parental gametes the individual received. Loci for which this cannot be determined but
                ! which are between two loci that (a) can be determined and (b) come from the same parental gamete,
                ! are assumed to come from this gamete (i.e. no double recombination event in between). Alleles are
                ! imputed in the individual when analysis in both directions of the chromosome has identified the
                ! same inherited gamete and when the parent is phased for this locus in the suggested gamete,
                ! subject to the restrictions that he number of recombination events for the individuals is less
                ! than a threshold and that the region in which two recombination events occurred exceeds a
                ! threshold length. Major sub-step 8 is iterated a number of times with increasingly relaxed
                ! restrictions. After each iteration, the minor sub-steps are also carried out.
                ! This subroutine corresponds to Major sub-step 8 from Hickey et al., 2012 (Appendix A)

                use Global
                use AlphaImputeSpecFileModule

                implicit none

                integer :: e,i,j,HetEnd,HetStart,RSide,LSide,PatMat,SireDamRL
                integer,dimension(:), allocatable :: WorkRight,WorkLeft
                integer :: CountRightSwitch,CountLeftSwitch,StartPt,EndPt,PedId

                integer :: StartDis,EndDis,StartJ,k
                integer,allocatable,dimension(:) :: TempVec
                real,allocatable,dimension(:) :: LengthVec
                type(AlphaImputeInput), pointer :: inputParams
                integer(kind=1) :: phase(2), tmpPhase

                inputParams => defaultInput

                allocate(WorkRight(inputParams%nsnp))
                allocate(WorkLeft(inputParams%nsnp))
                allocate(TempVec(inputParams%nsnp))
                allocate(LengthVec(inputParams%nsnp))


                do i=1,ped%pedigreeSize- ped%nDummys
                    HetEnd=-1
                    HetStart=-1

                    ! For each gamete
                    do e=1,2
                        PatMat=e
                        SireDamRL=e+1
                        CountLeftSwitch=0
                        CountRightSwitch=0
                        pedID=ped%pedigree(i)%getSireDamNewIDByIndex(SireDamRL)
                        if (ped%isDummy(pedID)) cycle
                        ! TODO can  probably skip if value is 0 too
                        ! Skip if, in the case of sex chromosome, me and my parent are heterogametic
                        if ((inputParams%SexOpt==1).and.(ped%pedigree(i)%gender==inputParams%hetGameticStatus).and.(ped%pedigree(i)%getParentGenderBasedOnIndex(SireDamRL)==inputParams%hetGameticStatus)) cycle

                        !! SCAN HAPLOTYPE IN TWO DIRECTIONS: L->R AND R->L
                        ! If not a base animal
                        if (PedId>0) then
                            WorkRight=9
                            RSide=9

                            ! Go through haplotype from Left to Right
                            ! finding the first heterozygous allele of this parent, and...
                            do j=1,inputParams%nsnp
                                phase(1) = ped%pedigree(pedId)%individualPhase(1)%getPhase(j)
                                phase(2) = ped%pedigree(pedId)%individualPhase(2)%getPhase(j)
                                if ((phase(1)/=phase(2)).and.&
                                    (phase(1)/=9).and.(phase(2)/=9))  then
                                    HetStart=j

                                    tmpPhase = ped%pedigree(i)%individualPhase(patmat)%getPhase(HetStart)
                                    ! Check if this allele corresponds to my parent paternal haplotype
                                    if (tmpPhase==phase(1)) then
                                        WorkRight(HetStart)=1   ! HetStart allele corresponds to Pat Haplotype in the LR direction
                                        RSide=1                 ! We are actually in the Paternal haplotype for the LR direction
                                        exit
                                    endif
                                    ! Check if this allele corresponds to my parent maternal haplotype
                                    if (tmpPhase==phase(2)) then
                                        WorkRight(HetStart)=2   ! HetStart allele corresponds to Mat Haplotype in the LR direction
                                        RSide=2                 ! We are actually in the Maternal haplotype for the LR direction
                                        exit
                                    endif
                                endif
                            enddo

                            ! ... Identifying recombinations
                            if (RSide/=9) then
                                do j=HetStart+1,inputParams%nsnp
                                    phase(1) = ped%pedigree(pedId)%individualPhase(RSide)%getPhase(j)
                                    tmpPhase = ped%pedigree(i)%individualPhase(patMat)%getPhase(j)
                                    ! If this allele has different phased as the current haplotype, then
                                    ! Change haplotype and increase the number of recombinations of this haplotype
                                    if ((tmpPhase/=phase(1)).and.&
                                        ((phase(1)==9).and.(tmpPhase/=9))) then
                                        RSide=abs((RSide-1)-1)+1
                                        CountRightSwitch=CountRightSwitch+1
                                    endif
                                    ! Wich paternal gamete the individual has received
                                    WorkRight(j)=RSide
                                enddo
                            endif

                            WorkLeft=9
                            LSide=9

                            ! Go through haplotype from Right to Left
                            ! finding the first heterozygous allele of this parent, and...
                            do j=inputParams%nsnp,1,-1
                                phase(1) = ped%pedigree(pedId)%individualPhase(1)%getPhase(j)
                                phase(2)  = ped%pedigree(pedId)%individualPhase(2)%getPhase(j)
                                if ((phase(1)/=phase(2)).and.&
                                    (phase(1)/=9).and.(phase(2)/=9))  then
                                    HetEnd=j
                                    tmpPhase = ped%pedigree(i)%individualPhase(PatMat)%getPhase(hetEnd)
                                    if (tmpPhase==phase(1)) then
                                        WorkRight(HetEnd)=1
                                        LSide=1
                                        exit
                                    endif
                                    if (tmpPhase==phase(2)) then
                                        WorkRight(HetEnd)=2
                                        LSide=2
                                        exit
                                    endif
                                endif
                            enddo

                            ! ... Identifying recombinations
                            if (LSide/=9) then
                                do j=HetEnd-1,1,-1
                                    phase(1) = ped%pedigree(i)%individualPhase(patMat)%getPhase(j)
                                    phase(2) = ped%pedigree(pedId)%individualPhase(LSide)%getPhase(j)
                                    if ((phase(1)/=phase(2) ).and.&
                                        (phase(2) /=9).and.(phase(1)/=9) ) then
                                        LSide=abs((LSide-1)-1)+1
                                        CountLeftSwitch=CountLeftSwitch+1
                                    endif
                                    ! Wich paternal gamete the individual has received
                                    WorkLeft(j)=LSide
                                enddo
                            endif

                            !$$$$$$$$$$$$$$$$$$$

                            ! UNPHASE THOSE ALLELES WITH SOME AMBIGUITY DUE TO RECOMBINATION
                            ! Let's be (StartDis:EndDis) the SNPs of the two direction disagree
                            StartDis=-9
                            EndDis=-9
                            TempVec=9
                            LengthVec=0.0
                            StartJ=1


                            do j=StartJ,inputParams%nsnp
                                ! Initalize variables StartDis and EndDis
                                ! StartDis is the first allele where different directions differ
                                if (StartDis==-9) then
                                    if (abs(WorkLeft(j)-WorkRight(j))==1) then
                                        StartDis=j
                                        TempVec(StartDis)=1
                                    endif
                                endif
                                ! EndDis is the last allele where different directions agree.
                                !   (StartDis/=-9) guarantees that EndDis > StartDis
                                !   (EndDis==-9) guarantees that EndDis is not updated
                                if ((WorkLeft(j)==WorkRight(j)).and.(WorkLeft(j)/=9).and.(StartDis/=-9).and.(EndDis==-9)) then
                                    EndDis=j-1
                                    TempVec(EndDis)=2
                                endif

                                ! Move StartDis to the first phased allele (from left) that comes from a heterozygous case
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

                                ! Move EndDis to the last phased allele (from left) that comes from a heterozygous case
                                if (EndDis/=-9) then
                                    do k=EndDis,inputParams%nsnp
                                        if ((GlobalWorkPhase(PedId,k,1)+GlobalWorkPhase(PedId,k,2))==1) then
                                            if (GlobalWorkPhase(i,k,e)/=9) then
                                                exit
                                            endif
                                        endif
                                    enddo
                                    EndDis=k
                                    if (EndDis>inputParams%nsnp) EndDis=inputParams%nsnp
                                    TempVec(EndDis)=2
                                endif

                                ! If StartDis==9 means that there haplotype is the same from Left to Right than from Right to Left
                                !   * In this case, EndDis==9 too, then DO NOTHING
                                ! If EndDis==9 means that
                                !   * StartDis==9 or
                                !   * WorkLeft(j)==9, which means that
                                !       - The whole haplotype is homozygous
                                !       - parent is not completely phased for that allele
                                if ((StartDis/=-9).and.(EndDis/=-9)) then
                                    ! WARNING: 3 is the only value that is used for TempVec
                                    TempVec(StartDis+1:EndDis-1)=3
                                    LengthVec(StartDis+1:EndDis-1)=1.0/(((EndDis-1)-(StartDis+1))+1)
                                    StartJ=EndDis+1
                                    StartDis=-9
                                    EndDis=-9
                                endif
                            enddo

                            ! Remove phase and genotype for those alleles with no explanation due to heterozygosity and recombination
                            block 
                                integer(kind=1) :: phase1,phase2
                                do j=1,inputParams%nSnp
                                    if ( TempVec(j)==3) then
                                        phase1 = ped%pedigree(pedId)%individualPhase(1)%getPhase(j)
                                        phase2 = ped%pedigree(pedId)%individualPhase(2)%getPhase(j)
                                        if ((phase1/=phase2) .AND. &
                                        (phase1/=9) .and. (phase2/=9) ) then

                                            call ped%pedigree(i)%IndividualPhase(e)%setPhase(j,9)
                                            if (Setter(i)/=1) then ! Skip HD individuals
                                                call ped%pedigree(i)%IndividualGenotype%setGenotype(j,9)
                                            end if
                                        endif
                                    endif
                                enddo
                        
                            end block



                            !$$$$$$$$$$$$$$$$$$$

                            !! IMPUTE PHASE WHETHER IT IS POSSIBLE
                            ! WARNING: From Hickey et al. 2012 (Appendix A):
                            !          ["Alleles are imputed ... subject to the restriction that the number of recombinations events for the individuals is less than a threshold, AND
                            !          that the region in which two recombination events occurred exceeds a threshold lenght."]
                            !          What it is coded is ["... than a threshold, OR that the region..."]
                            ! The number of recombinations in total (LR + RL) is less than a threshold
                            if ((CountLeftSwitch+CountRightSwitch)<(2*MaxLeftRightSwitch)) then
                                do j=1,inputParams%nsnp
                                
                                    if (ped%pedigree(i)%individualPhase(PatMat)%isMissing(j)) then

                                        ! WARNING: This can be coded in a conciser way
                                        ! Phase if the allele in one of the two directions is missing
                                        if ((WorkLeft(j)==9).and.(WorkRight(j)/=9)) then
                                            call ped%pedigree(i)%individualPhase(patMat)%setPhase(j,ped%pedigree(pedId)%individualPhase(workRight(j))%getPhase(j))
                                        else if ((WorkLeft(j)/=9).and.(WorkRight(j)==9)) then
                                            call ped%pedigree(i)%individualPhase(patMat)%setPhase(j,ped%pedigree(pedId)%individualPhase(WorkLeft(j))%getPhase(j))
                                        ! Phase if alleles is the two directions agree
                                        else if ((WorkLeft(j)/=9).and.(WorkRight(j)==WorkLeft(j))) then
                                            call ped%pedigree(i)%individualPhase(patMat)%setPhase(j,ped%pedigree(pedId)%individualPhase(WorkLeft(j))%getPhase(j))
                                        endif
                                    endif
                                enddo
                            else
                                ! Let's be (StartPt:EndPt) the SNPs in the two direction agree
                                EndPt=0
                                StartPt=0
                                do while ((StartPt<(inputParams%nsnp-MinSpan)).and.(EndPt<(inputParams%nsnp-MinSpan)))      ! If EndPt >(inputParams%nsnp-MinSpan), then recombination events does not exceed the threshold MinSpan
                                    do j=EndPt+1,inputParams%nsnp
                                        if ((WorkLeft(j)/=9).and.(WorkRight(j)==WorkLeft(j)))  then
                                            StartPt=j
                                            exit
                                        endif
                                        if (j==inputParams%nsnp) StartPt=j
                                    enddo
                                    do j=StartPt,inputParams%nsnp
                                        if ((WorkLeft(j)==9).or.(WorkRight(j)/=WorkLeft(j)))  then
                                            EndPt=j
                                            exit
                                        endif
                                        if (j==inputParams%nsnp) EndPt=j
                                    enddo
                                    ! The region in which two recombination events occurred exceeds a threshold lenght
                                    if (((EndPt-StartPt)+1)>MinSpan) then
                                        do j=StartPt,EndPt
                                            ! WARNING: This condition is supposed to be meet always since SNPs in (StartPt:EndPt)
                                            !          meet the condition (WorkLeft(j)/=9).and.(WorkRight(j)==WorkLeft(j))
                                            
                                            if (.not. ped%pedigree(pedId)%individualPhase(WorkRight(j))%isMissing(j)) then
                                                call ped%pedigree(i)%individualPhase(patMat)%setPhase(j,ped%pedigree(pedId)%individualPhase(WorkRight(j))%getPhase(j))
                                            endif
                                        enddo
                                    endif
                                enddo
                            endif
                        endif
                    enddo
                enddo

                ! Impute phase for the Heterogametic chromosome from the Homogametic one, which has been already phased
                do i=1,ped%pedigreeSize- ped%nDummys
                    if ((inputParams%SexOpt==1).and.(ped%pedigree(i)%gender==inputParams%hetGameticStatus)) then
                        ped%pedigree(i)%individualPhase(inputParams%hetGameticStatus) = ped%pedigree(i)%individualPhase(inputParams%HomGameticStatus)
                    
                    endif
                enddo
                
                GlobalWorkPhase =ped%getPhaseAsArray()

            end subroutine WorkLeftRight



    subroutine FromHMM2ImputePhase
        ! Impute alleles from HMM dosage probabilities
        use Global

        use AlphaImputeSpecFileModule

        implicit none

        integer :: i,j,k
        type(AlphaImputeInput), pointer :: inputParams

        inputParams => defaultInput

        do i=1,ped%nGenotyped
            do j=1,inputParams%nsnp
                do k=1,2
                    if (FullH(i,j,k)<0.001.and.FullH(i,j,k)>=0.0) Then
                        call ped%pedigree(ped%genotypeMap(i))%individualPhase(k)%setPhase(j,0)
                    elseif (FullH(i,j,k)>0.999.and.FullH(i,j,k)<=1.0) then
                        call ped%pedigree(ped%genotypeMap(i))%individualPhase(k)%setPhase(j,1)
                    else
                        call ped%pedigree(ped%genotypeMap(i))%individualPhase(k)%setPhase(j,9)
                    endif
                enddo
            enddo
        enddo

    end subroutine FromHMM2ImputePhase


        subroutine InsteadOfReReadGeneProb
        ! Phase alleles in the SEX CHROMOSOME whenever it is possible (homozygous case).
        ! Phasing information is store in the variable GlobalWorkPhase
        use Global
        use AlphaImputeSpecFileModule
        implicit none

        type(AlphaImputeInput), pointer :: inputParams
        integer :: e,i,j,ParId
        integer, dimension(:,:) , allocatable :: Genos !  Temp variable


        genos = ped%getGenotypesAsArray()

        inputParams => defaultInput
        if (defaultInput%SexOpt==1) then                                         ! Sex chromosome
            deallocate(GlobalWorkPhase)
            allocate(GlobalWorkPhase(0:ped%pedigreeSize-ped%nDummys,inputParams%nsnp,2))
            GlobalWorkPhase=9
            do i=1,ped%pedigreeSize
                do j=1,inputParams%nsnp                                                     ! Phase alleles in the homozygous case
                    if (Genos(i,j)==0) GlobalWorkPhase(i,j,:)=0
                    if (Genos(i,j)==2) GlobalWorkPhase(i,j,:)=1
                enddo
                if (ped%pedigree(i)%gender/=inputParams%hetGameticStatus) then
                    do e=1,2                                                    ! Phase alleles for homogametic individuals in the homozygous case
                        ParId=ped%pedigree(i)%getSireDamNewIDByIndex(e+1)
                        do j=1,inputParams%nsnp
                            if (Genos(ParId,j)==0) GlobalWorkPhase(i,j,e)=0
                            if (Genos(ParId,j)==2) GlobalWorkPhase(i,j,e)=1
                        enddo
                    enddo
                else
                    ParId=ped%pedigree(i)%getSireDamNewIDByIndex(inputParams%HomGameticStatus+1)                          ! Phase alleles for heterogametic individuals in the homozygous case
                    do j=1,inputParams%nsnp
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



        END MODULE Imputation
