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
        use informationModule
        use AlphaImputeInputOutputModule, only : ReadInPrePhasedData, ReReadGeneProbs,readgeneprobscluster
        use AlphaPhaseResultsModule

        integer :: loop
        character(len=150) :: timeOut

        inputParams => defaultInput



        print *," in imputation"
        ! WARNING: Need to discuss this part of code with John. Nonsense going on here!

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

                inputParamsHMM%nsnp = inputParams%nsnp
                inputParamsHMM%nHapInSubH = inputParams%nHapInSubH
                inputParamsHMM%HmmBurnInRound = inputParams%HmmBurnInRound
                inputParamsHMM%nRoundsHmm = inputParams%nRoundsHmm
                inputParamsHMM%useProcs = inputParams%useProcs
                inputParamsHMM%imputedThreshold = inputParams%imputedThreshold
                inputParamsHMM%phasedThreshold = inputParams%phasedThreshold
                inputParamsHMM%HapList = inputParams%HapList
                call AlphaImputeHMMRunner(inputParamsHMM, ImputeGenos, ImputePhase, ped, ProbImputeGenosHmm, ProbImputePhaseHmm, GenosCounts, FullH)


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

                    ! Impute phase whenever a pre-phase file exists
                    if (inputParams%PrePhased==1) call ReadInPrePhasedData

                    ! Impute phase in the sex chromosome
                    if (inputParams%sexopt==1) call EnsureHetGametic

                    ! General imputation procedures
                    call GeneralFillInInit
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

                                call AlphaImputeHMMRunner(inputParamsHMM, ImputeGenos, ImputePhase, ped, ProbImputeGenosHmm, ProbImputePhaseHmm, GenosCounts, FullH)


                            end block
                        else
                            print*, " "
                            print*, " ","Imputation of base animals completed"

                            do loop=1,inputParams%InternalIterations
                                print*, " "
                                print*, "Performing imputation loop",loop

                                call PhaseElimination                   ! Major Sub-Step 5 (Hickey et al., 2012; Appendix A)
                                if (inputParams%sexopt==1) then
                                    call EnsureHetGametic
                                end if
                                call GeneralFillIn
                                print*, " "
                                print*, " ","Parent of origin assigmnent of high density haplotypes completed"

                                
                                call ParentPhaseElimination             ! Major Sub-Step 4 (Hickey et al., 2012; Appendix A)
                                
                                if (inputParams%sexopt==1) then
                                    call EnsureHetGametic
                                end if


                                call GeneralFillIn
                                call RestrictedWorkLeftRight            ! Major Sub-Step 8 (Hickey et al., 2012; Appendix A)
                                if (inputParams%sexopt==1) then
                                    call EnsureHetGametic
                                end if
                                call GeneralFillIn
                                print*, " "
                                CALL DATE_AND_TIME(time=timeOut)
                                print*, " ","Imputation from high-density parents completed at: ",trim(timeOut)

                                call ImputeFromHDLibrary                ! Major Sub-Step 3 (Hickey et al., 2012; Appendix A)

                                if (inputParams%sexopt==1) then
                                    call EnsureHetGametic
                                end if
                                call GeneralFillIn
                                call RestrictedWorkLeftRight            ! Major Sub-Step 8 (Hickey et al., 2012; Appendix A)
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
                                call RestrictedWorkLeftRight            ! Major Sub-Step 8 (Hickey et al., 2012; Appendix A)
                                if (inputParams%sexopt==1) then
                                    call EnsureHetGametic
                                end if
                                call GeneralFillIn
                                print*, " "
                                CALL DATE_AND_TIME(time=timeOut)
                                print*, " ","Internal imputation from parents haplotype completed at: ",timeOut

                                call InternalHapLibImputationOld           ! Major Sub-Step 6 (Hickey et al., 2012; Appendix A)
                                if (inputParams%sexopt==1) then
                                    call EnsureHetGametic
                                end if
                                call GeneralFillIn
                                if (inputParams%sexopt==1) then
                                    call EnsureHetGametic
                                end if
                                call RestrictedWorkLeftRight            ! Major Sub-Step 8 (Hickey et al., 2012; Appendix A)
                                call GeneralFillIn
                                print*, " "
                                CALL DATE_AND_TIME(time=timeOut)
                                print*, " ","Internal haplotype library imputation completed at: ", timeOut
                                                      
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
                                call AlphaImputeHMMRunner(inputParamsHMM, ImputeGenos, ImputePhase, ped, ProbImputeGenosHmm, ProbImputePhaseHmm, GenosCounts, FullH)


                            end block
                            call FromHMM2ImputePhase
                        endif
                        deallocate(GlobalWorkPhase)
                    endif
                endif


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
                use HaplotypeBits

                implicit none

                integer :: m,e,g,i,j,nCore,nGlobalLoop,CoreLength,CoreStart,CoreEnd,CompPhase
                integer :: LoopStart,Offset,AnimalOn(ped%pedigreeSize,2)
                integer ::GamA,GamB

                integer,allocatable,dimension (:,:) :: LoopIndex
                integer(kind=1),allocatable,dimension (:,:,:,:) :: Temp

                integer(kind=int64), allocatable, dimension(:,:,:) :: BitImputePhase, MissImputePhase

                type(BitSection) :: Section
                integer :: numSections, curSection, curPos,l


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

                            Section = BitSection((CoreEnd - CoreStart + 1), 64)
                            numSections = Section%numSections

                            if (allocated(BitImputePhase)) then
                                deallocate(BitImputePhase)
                                deallocate(MissImputePhase)
                            end if

                            allocate(BitImputePhase(numSections,2,0:ped%pedigreeSize))
                            allocate(MissImputePhase(numSections,2,0:ped%pedigreeSize))

                            BitImputePhase = 0
                            MissImputePhase = 0

                            do i = 1, ped%pedigreeSize - ped%nDummys
                                do e = 1, 2
                                    curSection = 1
                                    curPos = 1
                                    do j = CoreStart, CoreEnd
                                        select case (ImputePhase(i, j, e))
                                        case (1)
                                            BitImputePhase(curSection, e, i) = ibset(BitImputePhase(curSection, e, i), curPos)
                                        case (9)
                                            MissImputePhase(curSection, e, i) = ibset(MissImputePhase(curSection, e, i), curPos)
                                        end select

                                        curPos = curPos + 1
                                        if (curPos == 65) then
                                            curPos = 1
                                            curSection = curSection + 1
                                        end if
                                    end do
                                end do
                            end do
                            block
                                use individualModule
                                type(individual), pointer :: parent

                                !$OMP PARALLEL DO &
                                !$OMP DEFAULT(SHARED) &
                                !$OMP PRIVATE(i,j,e,CompPhase,GamA,GamB,curPos,curSection, parent)
                                do i=1,ped%pedigreeSize-ped%ndummys

                                    do e=1,2

                                        if (ped%pedigree(i)%isDummyBasedOnIndex(e+1)) cycle
                                        parent => ped%pedigree(i)%getSireDamObjectByIndex(e+1)
                                        ! parent => ped%pedigree(ped%pedigree(i)%getSireDamNewIDByIndex(e+1))
                                        ! Skip if, in the case of sex chromosome, me and my parent are heterogametic
                                        if ((inputParams%sexopt==1).and.(ped%pedigree(i)%gender==inputParams%HetGameticStatus).and.(parent%gender==inputParams%HetGameticStatus)) cycle
                                        ! If not a Base Animal
                                        if (associated(parent)) then
                                            CompPhase=1
                                            ! If the haplotype for this core is not completely phased
                                            if (Section%BitCompletePhased(MissImputePhase(:,e,i)) == .FALSE.) then
                                                ! Check is this haplotype is the very same that the paternal haplotype
                                                ! of the parent of the individual
                                                GamA=1
                                                if (.NOT. Section%compareHaplotype(BitImputePhase(:,1,parent%id), BitImputePhase(:,e,i), &
                                                    MissImputePhase(:,1,parent%id), MissImputePhase(:,e,i)) ) then
                                                    GamA = 0
                                                end if

                                                ! Check is this haplotype is the very same that the maternal haplotype
                                                ! of the parent of the individual
                                                GamB=1
                                                if (.NOT. Section%compareHaplotype(BitImputePhase(:,2,parent%id), BitImputePhase(:,e,i), &
                                                    MissImputePhase(:,2,parent%id), MissImputePhase(:,e,i))) then
                                                    GamB = 0
                                                end if

                                                ! This haplotype is the paternal haplotype of the individual's parent
                                                ! Then count the number of occurrences a particular phase is impute in a
                                                ! a particular allele across the cores and across the internal phasing steps
                                                if ((GamA==1).and.(GamB==0)) then
                                                    AnimalOn(i,e)=1
                                                    curPos = 1
                                                    curSection = 1

                                                    do j=CoreStart,CoreEnd
                                                        if ( BTEST(MissImputePhase(curSection,e,i), curPos) == .TRUE. ) then
                                                            if ( BTEST(BitImputePhase(curSection,1,parent%id), curPos) == .FALSE. .AND. &
                                                                BTEST(MissImputePhase(curSection,1,parent%id), curPos)  == .FALSE. ) then
                                                                !$OMP ATOMIC
                                                                Temp(i,j,e,1)=Temp(i,j,e,1)+1
                                                            end if

                                                            if ( BTEST(BitImputePhase(curSection,1,parent%id), curPos) == .TRUE. ) then
                                                                !$OMP ATOMIC
                                                                Temp(i,j,e,2)=Temp(i,j,e,2)+1
                                                            end if
                                                        end if

                                                        curPos = curPos + 1
                                                        if (curPos == 65) then
                                                            curPos = 1
                                                            curSection = curSection + 1
                                                        end if
                                                    enddo

                                                endif

                                                ! This haplotype is the maternal haplotype of the individual's parent
                                                ! Then count the number of occurrences a particular phase is impute in a
                                                ! a particular allele across the cores and across the internal phasing steps
                                                if ((GamA==0).and.(GamB==1)) then
                                                    AnimalOn(i,e)=1
                                                    curPos = 1
                                                    curSection = 1
                                                    do j=CoreStart,CoreEnd

                                                        if ( BTEST(MissImputePhase(curSection,e,i), curPos) == .TRUE. ) then

                                                        ! if phase is 0
                                                            if ( BTEST(BitImputePhase(curSection,2,parent%id), curPos) == .FALSE. .AND. &
                                                                BTEST(MissImputePhase(curSection,2,parent%id), curPos) == .FALSE. ) then
                                                                !$OMP ATOMIC
                                                                Temp(i,j,e,1)=Temp(i,j,e,1)+1
                                                            end if

                                                            ! if  phase is 1
                                                            if ( BTEST(BitImputePhase(curSection,2, parent%id), curPos) == .TRUE. ) then
                                                                !$OMP ATOMIC
                                                                Temp(i,j,e,2)=Temp(i,j,e,2)+1
                                                            end if
                                                        end if

                                                        curPos = curPos + 1
                                                        if (curPos == 65) then
                                                            curPos = 1
                                                            curSection = curSection + 1
                                                        end if
                                                    enddo
                                                endif
                                            endif
                                        endif

                                    enddo
                                enddo
                                !$OMP END PARALLEL DO
                            end block
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
                                if (ImputePhase(i,j,e)==9) then
                                    if ((Temp(i,j,e,1)>inputParams%nAgreeInternalHapLibElim).and.(Temp(i,j,e,2)==0)) ImputePhase(i,j,e)=0
                                    if ((Temp(i,j,e,1)==0).and.(Temp(i,j,e,2)>inputParams%nAgreeInternalHapLibElim)) ImputePhase(i,j,e)=1
                                endif
                                enddo
                        end if
                    endif
            

                ! The individual is has one haplotype: In Sex Chromosome, the heterogametic case
                
                    if ((inputParams%sexopt==1).and.(ped%pedigree(i)%gender==inputParams%HetGameticStatus)) then
                        if (AnimalOn(i,inputParams%HomGameticStatus)==1) then
                            do j=1,inputParams%nsnp
                                if (ImputePhase(i,j,inputParams%HomGameticStatus)==9) then
                                    if ((Temp(i,j,inputParams%HomGameticStatus,1)>inputParams%nAgreeInternalHapLibElim).and.(Temp(i,j,inputParams%HomGameticStatus,2)==0))&
                                        ImputePhase(i,j,:)=0
                                    if ((Temp(i,j,inputParams%HomGameticStatus,1)==0).and.(Temp(i,j,inputParams%HomGameticStatus,2)>inputParams%nAgreeInternalHapLibElim))&
                                        ImputePhase(i,j,:)=1
                                endif
                            enddo
                        end if
                    end if
                enddo
            enddo


                deallocate(Temp)

                ImputePhase(0,:,:)=9
                ImputeGenos(0,:)=9

            end subroutine InternalParentPhaseElim

            !#############################################################################################################################################################################################################################

            subroutine InternalHapLibImputation
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
                use HaplotypeBits
                use AlphaImputeSpecFileModule
                implicit none

                integer :: f,e,h,g,i,j,l,nCore,nHap,nGlobalLoop,CoreLength,CoreStart,CoreEnd,InLib,Count0,Count1
                integer :: Counter,BanBoth(2),Ban(2),AnimalOn(ped%pedigreeSize,2)
                integer :: LoopStart,OffSet

                integer,allocatable,dimension (:,:) :: HapLib,LoopIndex,HapElim
                integer(kind=1),allocatable,dimension (:,:,:,:) :: Temp

                integer(kind=8), allocatable, dimension(:,:) :: BitHapLib, MissHapLib, BitWork, MissWork
                integer(kind=8), allocatable, dimension(:,:,:) :: BitImputePhase, MissImputePhase

                type(BitSection) :: Section
                integer :: numSections, curSection, curPos
                integer :: BitGeno

                inputParams => defaultInput
                ! WARNING: This should go in a function since it is the same code as InternalParentPhaseElim subroutine
                nGlobalLoop=25

                ! LoopeIndex is a matrix with two columns that will serve to define:
                !   * 1.- nCores
                !   * 2.- Core lengths
                allocate(LoopIndex(nGlobalLoop,2))
                call setLoopIndex(LoopIndex)

                ! LoopStart indicates which is the first loop the algorithm should treat. The bigger the number of
                ! SNPs, the more the loops to be considered
                LoopStart = getLoopStart(inputParams%nsnp)
                if (LoopStart == 0) return

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
                    ! Allocate the Internal Haplotype Library
                    allocate(HapLib(ped%pedigreeSize*2,inputParams%nsnp))
                    do l=LoopStart,nGlobalLoop

                        ! Simulate phase without shift
                        if (f==1) then
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
                            if ((f==1).and.(g==nCore)) CoreEnd=inputParams%nsnp
                            if ((f==2).and.(g==nCore)) CoreEnd=inputParams%nsnp-OffSet

                            ! Exit if the corelength is too small
                            CoreLength=(CoreEnd-CoreStart)+1
                            if (CoreLength<10) exit

                            nHap=0
                            HapLib=9

                            Section = BitSection((CoreEnd - CoreStart + 1), 64)
                            numSections = Section%numSections

                            if (allocated(BitHapLib)) then
                                deallocate(BitHapLib)
                                deallocate(MissHapLib)
                                deallocate(BitImputePhase)
                                deallocate(MissImputePhase)
                            end if

                            allocate(BitHapLib(ped%pedigreeSize*2,numSections))
                            allocate(MissHapLib(ped%pedigreeSize*2,numSections))
                            allocate(BitImputePhase(0:ped%pedigreeSize,numSections,2))
                            allocate(MissImputePhase(0:ped%pedigreeSize,numSections,2))

                            BitHapLib = 0
                            MissHapLib = 0
                            BitImputePhase = 0
                            MissImputePhase = 0

                            curSection = 1
                            curPos = 1
                            do j = CoreStart, CoreEnd
                                do e = 1, 2
                                    do i = 1, ped%pedigreeSize- ped%nDummys
                                        select case (ImputePhase(i, j, e))
                                        case (1)
                                            BitImputePhase(i, curSection, e) = ibset(BitImputePhase(i, curSection, e), curPos)
                                        case (9)
                                            MissImputePhase(i, curSection, e) = ibset(MissImputePhase(i, curSection, e), curPos)
                                        end select
                                    end do
                                end do

                                do i = 1, (ped%pedigreeSize- ped%nDummys)*2
                                    MissHapLib(i, curSection) = IBSET(MissHapLib(i, curSection), curPos)
                                end do

                                curPos = curPos + 1
                                if (curPos == 65) then
                                    curPos = 1
                                    curSection = curSection + 1
                                end if
                            end do

                            !$!OMP PARALLEL DO ORDERED&
                            !$!OMP DEFAULT(SHARED) &
                            !$!OMP PRIVATE(i,j,e,h,InLib)
                            do i=1,ped%pedigreeSize- ped%nDummys
                                do e=1,2
                                    ! If GeneProbPhase has been executed, that is, if not considering the Sex Chromosome, then MSTermInfo={0,1}.
                                    ! Else, if Sex Chromosome, then MSTermInfo is 0 always
                                    ! So, if a Conservative imputation of haplotypes is selected, this DO statement will do nothing
                                    if ((inputParams%ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle

                                    ! Check if the haplotype for this core is completely phased and populate HapLib
                                    if (Section%BitCompletePhased(MissImputePhase(i,:,e)) ) then
                                        !$OMP ORDERED
                                        if (nHap==0) then       ! The first haplotype in the library
                                            HapLib(1,CoreStart:CoreEnd)=ImputePhase(i,CoreStart:CoreEnd,e)
                                            BitHapLib(1,:) = BitImputePhase(i,:,e)
                                            MissHapLib(1,:) = MissImputePhase(i,:,e)
                                            nHap=1
                                        else
                                            InLib=0
                                            do h=1,nHap
                                                if (Section%compareHaplotypeAllowMissing(BitHapLib(h,:), BitImputePhase(i,:,e), &
                                                    MissHapLib(h,:), MissImputePhase(i,:,e))) then
                                                    InLib = 1
                                                    exit
                                                end if
                                            enddo

                                            ! If haplotype is not in the library,
                                            ! a new haplotype has been found, then populate the library
                                            if (InLib==0) then
                                                !$!OMP ATOMIC
                                                nHap=nHap+1
                                                HapLib(nHap,CoreStart:CoreEnd)=ImputePhase(i,CoreStart:CoreEnd,e)
                                                BitHapLib(nHap,:) = BitImputePhase(i,:,e)
                                                MissHapLib(nHap,:) = MissImputePhase(i,:,e)
                                            endif
                                        endif
                                        !$OMP END ORDERED
                                    endif
                                enddo
                            enddo
                            !$!OMP END PARALLEL DO ORDERED

                            ! WARNING: This code does not match the corresponding code of the subroutine ImputeFromHDLibrary
                            !          In ImputeFromHDLibrary, there are two steps, counting agreements and impute
                            !          across candidate haplotypes, and counting agreements and impute across cores
                            !          and phasing steps.
                            allocate(BitWork(numSections,2))
                            allocate(MissWork(numSections,2))
                            allocate(HapElim(ped%pedigreeSize*2,2))
                            HapElim = 0 !  changed from HapElim=1 to improve speed
                            !$!OMP PARALLEL DO &
                            !$!OMP DEFAULT(SHARED) &
                            !$!OMP PRIVATE(i,j,e,h,HapElim,BanBoth,counter,Count0,Count1,Ban,curPos,curSection,BitWork,MissWork,BitGeno)
                            do i=1,ped%pedigreeSize- ped%nDummys


                                BanBoth=0
                                BitWork = 0
                                MissWork = 0

                                do e=1,2

                                    counter = nHap*2
                                    ! If GeneProbPhase has been executed, that is, if not considering the Sex Chromosome, then MSTermInfo={0,1}.
                                    ! Else, if Sex Chromosome, then MSTermInfo is 0 always
                                    ! So, if a Conservative imputation of haplotypes is selected, this DO statement will do nothing
                                    if ((inputParams%ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle

                                    ! If haplotype is partially phased
                                    if ( (Section%BitCompletePhased(MissImputePhase(i,:,e)) == .FALSE.) .AND. &
                                        (Section%BitCompleteMissing(MissImputePhase(i,:,e)) == .FALSE.) ) then

                                        ! Identify and reject the candidate haplotypes if it does not explain the whole haplotype
                                        do h=1,nHap
                                            if ( .NOT. Section%compareHaplotypeAllowMissing(BitHapLib(h,:), BitImputePhase(i,:,e), &
                                                MissHapLib(h,:), MissImputePhase(i,:,e))) then
                                                HapElim(h,e)=i
                                                counter = counter - 1 !same as previous but without having to acutally use count
                                            end if
                                        enddo

                                        ! If the number of candidate haplotypes is less than the 25% of the Library,
                                        ! then impute if all alleles have been phased the same way
                                        ! Counter=count(HapElim(1:nHap,e)==1)
                                        if (float(Counter)<(float(nHap)*0.25)) then
                                            ! Ban this haplotype will be phased here and nowhere else
                                            BanBoth(e)=1
                                            curPos = 1
                                            curSection = 1

                                            do j=CoreStart,CoreEnd
                                                ! How many haplotypes has been phased as 0 or 1 in this particular allele?
                                                Count0=0
                                                Count1=0

                                                ! Count the occurrences in phasing of alleles across candidate haplotypes
                                                do h=1,nHap
                                                    if (HapElim(h,e)/=i) then
                                                        if (BTEST(BitHapLib(h,curSection), curPos) == .FALSE. .AND. &
                                                            BTEST(MissHapLib(h,curSection), curPos) == .FALSE.) then
                                                            Count0=Count0+1
                                                        end if
                                                        if (BTEST(BitHapLib(h,curSection), curPos) == .TRUE.) then
                                                            Count1=Count1+1
                                                        end if
                                                        if ((Count0>0).and.(Count1>0)) exit
                                                    endif
                                                enddo

                                                ! If all alleles across the candidate haplotypes have been phased the same way, impute
                                                if (Count0 == 0 .AND. Count1 > 0) then
                                                    BitWork(curSection, e) = IBSET(BitWork(curSection, e),curPos)
                                                else if (Count0 > 0 .AND. Count1 == 0) then
                                                    BitWork(curSection, e) = IBCLR(BitWork(curSection, e),curPos)
                                                else
                                                    MissWork(curSection, e) = IBSET(MissWork(curSection, e),curPos)
                                                end if

                                                curPos = curPos + 1
                                                if (curPos == 65) then
                                                    curPos = 1
                                                    curSection = curSection + 1
                                                end if
                                            enddo
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
                                    curPos = 1
                                    curSection = 1
                                    do j=CoreStart,CoreEnd
                                        ! If Working genotype is the same as the imputed one
                                        if (ImputeGenos(i,j)/=9) then
                                            if ( BTEST(MissWork(curSection,1), curPos) == .FALSE. .AND. &
                                                BTEST(MissWork(curSection,2), curPos) == .FALSE. ) then
                                                BitGeno = 0
                                                if (BTEST(BitWork(curSection,1), curPos)) BitGeno = BitGeno + 1
                                                if (BTEST(BitWork(curSection,2), curPos)) BitGeno = BitGeno + 1
                                                if (ImputeGenos(i,j)/=BitGeno) then
                                                    Ban=0
                                                    exit
                                                end if
                                            endif
                                        end if
                                        curPos = curPos + 1
                                        if (curPos == 65) then
                                            curPos = 1
                                            curSection = curSection + 1
                                        end if
                                    enddo
                                endif

                                ! Count the number of occurrences a particular phase is impute in a particular
                                ! allele across the cores and across the internal phasing steps
                                ! This implies occurrences across all the haplotype libraries of the different internal phasing steps
                                do e=1,2
                                    if ((inputParams%ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle
                                    if (Ban(e)==1) then
                                        AnimalOn(i,e)=1
                                        curPos = 1
                                        curSection = 1
                                        do j=CoreStart,CoreEnd
                                            if (BTEST(BitWork(curSection,e), curPos) == .FALSE. .AND. &
                                                BTEST(MissWork(curSection,e), curPos) == .FALSE.) then
                                                !$OMP ATOMIC
                                                Temp(i,j,e,1)=Temp(i,j,e,1)+1
                                            end if
                                            if (BTEST(BitWork(curSection,e), curPos) == .TRUE.) then
                                                !$OMP ATOMIC
                                                Temp(i,j,e,2)=Temp(i,j,e,2)+1
                                            end if

                                            curPos = curPos + 1
                                            if (curPos == 65) then
                                                curPos = 1
                                                curSection = curSection + 1
                                            end if
                                        enddo
                                    endif
                                enddo

                            enddo
                            !$!OMP END PARALLEL DO


                            deallocate(HapElim)
                            deallocate(BitWork)
                            deallocate(MissWork)

                            ! Prepare the core for the next cycle
                            CoreStart=CoreStart+LoopIndex(l,2)
                            CoreEnd=CoreEnd+LoopIndex(l,2)
                            if ((f==2).and.(g==nCore)) exit
                        enddo
                    enddo
                    deallocate(HapLib)
                enddo

                !$!OMP PARALLEL DEFAULT(SHARED)
                do e=1,2
                    !$!OMP DO PRIVATE(i,j)
                    do j=1,inputParams%nsnp
                        do i=1,ped%pedigreeSize- ped%nDummys
                            ! If GeneProbPhase has been executed, that is, if not considering the Sex Chromosome, then MSTermInfo={0,1}.
                            ! Else, if Sex Chromosome, then MSTermInfo is 0 always
                            ! So, if a Conservative imputation of haplotypes is selected, this DO statement will do nothing
                            if ((inputParams%ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle

                            ! If all alleles across the cores and across the internal phasing steps have been phased the same way, impute
                            if (AnimalOn(i,e)==1) then
                                if (ImputePhase(i,j,e)==9) then
                                    if ((Temp(i,j,e,1)>inputParams%nAgreeInternalHapLibElim).and.(Temp(i,j,e,2)==0)) then
                                        ImputePhase(i,j,e)=0
                                    end if
                                    if ((Temp(i,j,e,1)==0).and.(Temp(i,j,e,2)>inputParams%nAgreeInternalHapLibElim)) then
                                        ImputePhase(i,j,e)=1
                                    end if
                                endif
                            end if
                        enddo
                    enddo
                    !$!OMP END DO
                enddo
                !$!OMP END PARALLEL

                deallocate(Temp)

                ImputePhase(0,:,:)=9
                ImputeGenos(0,:)=9

            END SUBROUTINE InternalHapLibImputation



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
implicit none

integer :: f,e,h,g,i,j,k,l,nCore,nHap,nGlobalLoop,CoreLength,CoreStart,CoreEnd,InLib,NotHere,CompPhase,Count0,Count1
integer, dimension(:,:), allocatable :: Work
integer :: Counter,BanBoth(2),Ban(2),AnimalOn(ped%pedigreesize-ped%ndummys,2)
integer :: LoopStart,OffSet

integer,allocatable,dimension (:,:) :: CoreIndex,HapLib,LoopIndex,HapElim
integer(kind=1),allocatable,dimension (:,:,:,:) :: Temp


inputParams => defaultInput

! WARNING: This should go in a function since it is the same code as InternalParentPhaseElim subroutine
nGlobalLoop=25


allocate(work(inputParams%nsnp,2))
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
    ! Allocate the Internal Haplotype Library
    allocate(HapLib(ped%pedigreeSize*2,inputParams%nsnp))
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

            nHap=0
            HapLib=9

            ! THE FIRST PARALLELIZATION BEGINS: POPULATE THE INTERNAL HAPLOTYPE LIBRARY
            !# PARALLEL DO SHARED (ped%pedigreesize-ped%ndummys,ImputePhase,CoreStart,CoreEnd,nHap,HapLib) private(i,e,CompPhase,InLib,h,NotHere,j)
            do i=1,ped%pedigreesize-ped%ndummys
                do e=1,2
                    ! WARNING: If GeneProbPhase has been executed, that is, if not considering the Sex Chromosome, then MSTermInfo={0,1}.
                    !          Else, if Sex Chromosome, then MSTermInfo is 0 always
                    !          So, if a Conservative imputation of haplotypes is selected, this DO statement will do nothing
                    if ((inputParams%ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle

                    ! Check if the haplotype for this core is completely phased
                    CompPhase=1
                    if (count(ImputePhase(i,CoreStart:CoreEnd,e)==9)>0) CompPhase=0

                    ! If haplotype is completely phased, populate HapLib
                    ! NOTE: Since there is code in order to populate the Haplotype Library in
                    !       in AlphaPhase, it can be convenient to create a share procedure in
                    !       AlphaHouse
                    if (CompPhase==1) then
                        if (nHap==0) then       ! The first haplotype in the library
                            HapLib(1,CoreStart:CoreEnd)=ImputePhase(i,CoreStart:CoreEnd,e)
                            nHap=1
                        else
                            InLib=0
                            do h=1,nHap
                                NotHere=1
                                ! Check if haplotype is in the Library
                                do j=CoreStart,CoreEnd
                                    if(ImputePhase(i,j,e)/=HapLib(h,j)) then
                                        NotHere=0
                                        exit
                                    endif
                                enddo
                                ! If haplotype in the library, do nothing
                                if (NotHere==1) then
                                    InLib=1
                                    exit
                                endif
                            enddo
                            ! If haplotype is not in the library, then
                            ! a new haplotype has been found, then populate the library
                            if (InLib==0) then
                                nHap=nHap+1
                                HapLib(nHap,CoreStart:CoreEnd)=ImputePhase(i,CoreStart:CoreEnd,e)
                            endif
                        endif
                    endif
                enddo
            enddo
            !# END PARALLEL DO

            ! THE SECOND PARALLELIZATION BEGINS
            !# PARALLEL DO SHARED (ped%pedigreesize-ped%ndummys,ImputePhase,CoreStart,CoreEnd,CoreLength,nHap,HapLib,ImputeGenos,AnimalOn,Temp) private(i,HapElim,Work,BanBoth,e,h,j,Counter,Count0,Count1,Ban)

            ! WARNING: This code does not match the corresponding code of the subroutine ImputeFromHDLibrary
            !          In ImputeFromHDLibrary, there are two steps, counting agreements and impute
            !          across candidate haplotypes, and counting agreements and impute across cores
            !          and phasing steps.
            do i=1,ped%pedigreesize-ped%ndummys
                allocate(HapElim(ped%pedigreesize-ped%ndummys*2,2))
                HapElim=1
                Work=9
                BanBoth=0
                do e=1,2
                    ! WARNING: If GeneProbPhase has been executed, that is, if not considering the Sex Chromosome, then MSTermInfo={0,1}.
                    !          Else, if Sex Chromosome, then MSTermInfo is 0 always
                    !          So, if a Conservative imputation of haplotypes is selected, this DO statement will do nothing
                    if ((inputParams%ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle

                    ! If haplotype is partially phased
                    if ((count(ImputePhase(i,CoreStart:CoreEnd,e)==9)/=CoreLength)&
                            .and.(count(ImputePhase(i,CoreStart:CoreEnd,e)/=9)/=CoreLength)) then

                        ! Identify and reject the candidate haplotypes if it does not explain the whole haplotype
                        do h=1,nHap
                            do j=CoreStart,CoreEnd
                                if ((ImputePhase(i,j,e)/=9)&
                                        .and.(ImputePhase(i,j,e)/=HapLib(h,j))) then
                                    HapElim(h,e)=0
                                    exit
                                endif
                            enddo
                        enddo

                        ! If the number of candidate haplotypes is less than the 25% of the Library,
                        ! then impute if all alleles have been phased the same way
                        Counter=count(HapElim(1:nHap,e)==1)
                        if (float(Counter)<(float(nHap)*0.25)) then
                            ! Ban this haplotype will be phased here and nowhere else
                            BanBoth(e)=1
                            do j=CoreStart,CoreEnd
                                ! How many haplotypes has been phased as 0 or 1 in this particular allele?
                                Count0=0
                                Count1=0

                                ! Count the occurrences in phasing of alleles across candidate haplotypes
                                do h=1,nHap
                                    if (HapElim(h,e)==1) then
                                        if (HapLib(h,j)==0) Count0=Count0+1
                                        if (HapLib(h,j)==1) Count1=Count1+1
                                        if ((Count0>0).and.(Count1>0)) exit
                                    endif
                                enddo

                                ! If all alleles across the candidate haplotypes have been phased the same way, impute
                                if (Count0>0) then
                                    if (Count1==0) Work(j,e)=0
                                else
                                    if (Count1>0) Work(j,e)=1
                                endif
                            enddo
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

                ! Count the number of occurrences a particular phase is impute in a particular
                ! allele across the cores and across the internal phasing steps
                ! This implies occurrences across all the haplotype libraries of the different internal phasing steps
                do e=1,2
                    if ((inputParams%ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle
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

            ! Prepare the core for the next cycle
            CoreStart=CoreStart+LoopIndex(l,2)
            CoreEnd=CoreEnd+LoopIndex(l,2)
            if ((f==2).and.(g==nCore)) exit
        enddo
    enddo
    deallocate(HapLib)
enddo


do i=1,ped%pedigreesize-ped%ndummys
    do e=1,2
        ! WARNING: If GeneProbPhase has been executed, that is, if not considering the Sex Chromosome, then MSTermInfo={0,1}.
        !          Else, if Sex Chromosome, then MSTermInfo is 0 always
        !          So, if a Conservative imputation of haplotypes is selected, this DO statement will do nothing
        if ((inputParams%ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle

        ! If all alleles across the cores and across the internal phasing steps have been phased the same way, impute
        if (AnimalOn(i,e)==1) then
            do j=1,inputParams%nsnp
                if (ImputePhase(i,j,e)==9) then
                    if ((Temp(i,j,e,1)>inputParams%nAgreeInternalHapLibElim).and.(Temp(i,j,e,2)==0))&
                        ImputePhase(i,j,e)=0
                    if ((Temp(i,j,e,1)==0).and.(Temp(i,j,e,2)>inputParams%nAgreeInternalHapLibElim))&
                        ImputePhase(i,j,e)=1
                endif
            enddo
        endif
    enddo
enddo

deallocate(Temp)

ImputePhase(0,:,:)=9
ImputeGenos(0,:)=9

end subroutine InternalHapLibImputationOld

            !#############################################################################################################################################################################################################################


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


                use HaplotypeBits
                use AlphaImputeSpecFileModule
                use AlphaPhaseResultsModule
                implicit none

                integer :: e,g,i,j,GamA,GamB,PosHDInd
                integer :: StartSnp,EndSnp,Gam1,Gam2,AnimalOn(ped%pedigreeSize,2)
                integer,allocatable,dimension (:,:,:) :: PhaseHD
                integer(kind=8), allocatable, dimension(:,:,:) :: BitPhaseHD, BitImputePhase, MissPhaseHD, MissImputePhase
                integer(kind=1),allocatable,dimension (:,:,:,:) :: Temp
                integer :: unknownFreeIterator

                type(BitSection) :: Section



                integer :: numSections, curSection, curPos


                inputParams => defaultInput
                ! Number of animals that have been HD phased

                allocate(Temp(ped%pedigreeSize,inputParams%nsnpraw,2,2))
                allocate(PhaseHD(ped%nHd,inputParams%nsnpraw,2))

                Temp=0
                AnimalOn=0
                ! FOR EACH CORE OF EACH ROUND OF THE LRPHLI
                do unknownFreeIterator=1,apresults%nResults
                    !    TODOPhase read in phase info here if not red
                    ! Get phase information
                    ! call ReadPhased(ped%nHd, FileNamePhase, Ped, PhaseHD, PosHD)

                    startSnp = 1
                    EndSnp = 0

                    phaseHd = apResults%results(unknownFreeIterator)%getFullPhaseIntArray()
                    do g=1,size(apResults%results(unknownFreeIterator)%cores)
                        ! Initialize Start and End snps of the cores
                        StartSnp=apResults%results(unknownFreeIterator)%startIndexes(g)
                        EndSnp=apResults%results(unknownFreeIterator)%endIndexes(g)

                        Section = BitSection((EndSnp - StartSnp + 1), 64)
                        numSections = Section%numSections

                        allocate(BitPhaseHD(ped%nHd,numSections,2))
                        allocate(BitImputePhase(0:ped%pedigreeSize,numSections,2))
                        allocate(MissPhaseHD(ped%nHd,numSections,2))
                        allocate(MissImputePhase(0:ped%pedigreeSize,numSections,2))

                        BitPhaseHD = 0
                        MissPhaseHD = 0
                        BitImputePhase = 0
                        MissImputePhase = 0

                        do e=1,2
                            curSection = 1
                            curPos = 1

                            do j = StartSnp, EndSnp

                                do i = 1, ped%nHd
                                    select case (PhaseHD(i, j, e))
                                    case (1)
                                        ! set that phase information exists
                                        BitPhaseHD(i, curSection, e) = ibset(BitPhaseHD(i, curSection, e), curPos)
                                    case (9)
                                        ! set that missing iformation does not
                                        MissPhaseHD(i, curSection, e) = ibset(MissPhaseHD(i, curSection, e), curPos)
                                    end select
                                end do

                                do i = 1, ped%pedigreeSize- ped%nDummys
                                    select case (ImputePhase(i, j, e))
                                    case (1)
                                        BitImputePhase(i, curSection, e) = ibset(BitImputePhase(i, curSection, e), curPos)
                                    case (9)
                                        MissImputePhase(i, curSection, e) = ibset(MissImputePhase(i, curSection, e), curPos)
                                    end select
                                end do

                                curPos = curPos + 1
                                if (curPos == 65) then
                                    curPos = 1
                                    curSection = curSection + 1
                                end if
                            end do
                        end do

                        !$OMP PARALLEL DO &
                        !$OMP DEFAULT(SHARED) &
                        !$OMP FIRSTPRIVATE(i,j,e,Gam1,Gam2,GamA,GamB) &
                        !$OMP PRIVATE(PosHDInd)
                        do i=1,ped%nHd
                            ! Look for possible gametes through the Haplotype
                            ! Library constructed during the phasing step
                            PosHDInd=ped%hdMap(i)         ! Index of the individual in the HD phase information

                            ! If there is one allele phased at least
                            if ((Section%BitCountAllelesImputed(MissImputePhase(PosHDInd,:,1)) + &
                                Section%BitCountAllelesImputed(MissImputePhase(PosHDInd,:,2))) > 0) then
                                ! If at least one locus is heterozygous
                                if (.NOT. Section%compareHaplotype(BitImputePhase(PosHDInd,:,1), BitImputePhase(PosHDInd,:,2), &
                                    MissImputePhase(PosHDInd,:,1), MissImputePhase(PosHDInd,:,2))) then
                                    Gam1=0
                                    Gam2=0
                                    do e=1,2
                                        GamA=18
                                        GamB=1

                                        if (.NOT. Section%compareHaplotypeAllowMissing(BitPhaseHD(i,:,1), BitImputePhase(PosHDInd,:,e), &
                                            MissPhaseHD(i,:,1), MissImputePhase(PosHDInd,:,e), ImputeFromHDPhaseThresh)) then
                                            GamA = 0
                                        end if

                                        if (.NOT. Section%compareHaplotypeAllowMissing(BitPhaseHD(i,:,2), BitImputePhase(PosHDInd,:,e), &
                                            MissPhaseHD(i,:,2), MissImputePhase(PosHDInd,:,e), ImputeFromHDPhaseThresh)) then
                                            GamB = 0
                                        end if


                                        if (e == 1) Then
                                            ! Paternal haplotype is strictly my paternal haplotype from the Hap Library
                                            if (GamA==1 .and. GamB == 0) Then
                                                gam1 = 1
                                                ! Paternal haplotype is strictly my maternal haplotype from the Hap Library
                                            else if (GamA==0 .and. GamB == 1) Then
                                                gam2 = 0
                                            endif

                                        else if (e==2) Then
                                        ! Maternal haplotype is strictly my paternal haplotype from the Hap Library
                                            if (GamA==1 .and. GamB == 0) Then
                                                gam2 = 1
                                        ! Maternal haplotype is strictly my maternal haplotype from the Hap Library
                                            else if (GamA == 0 .and. GamB == 1) Then
                                                gam2 = 2
                                            endif
                                        endif
                                        ! Basically the important thing is that haplotype e is present in the Haplotype
                                        ! library. It is not important which haplotype it is, whether the paternal or the
                                        ! maternal.
                                    enddo

                                    ! If the paternal and maternal gametes are different
                                    if (Gam1/=Gam2) then
                                        AnimalOn(PosHDInd,:)=1             ! Consider this animal in further steps

                                        ! Paternal gamete is in the Hap Library
                                        if (Gam1/=0) then
                                            do j=StartSnp,EndSnp
                                                ! Count the number of alleles coded with 0 and 1
                                                if (ImputePhase(PosHDInd,j,1)==9) then
                                                    if(PhaseHD(i,j,Gam1)==0) then
                                                        !$OMP ATOMIC
                                                        Temp(PosHDInd,j,1,1)=Temp(PosHDInd,j,1,1)+1
                                                    end if
                                                    if(PhaseHD(i,j,Gam1)==1) then
                                                        !$OMP ATOMIC
                                                        Temp(PosHDInd,j,1,2)=Temp(PosHDInd,j,1,2)+1
                                                    end if
                                                endif
                                            enddo
                                        endif

                                        ! Maternal gamete is in the Hap Library
                                        if (Gam2/=0) then
                                            do j=StartSnp,EndSnp
                                                ! Count the number of alleles coded with 0 and 1
                                                if (ImputePhase(PosHDInd,j,2)==9) then
                                                    if(PhaseHD(i,j,Gam2)==0) then
                                                        !$OMP ATOMIC
                                                        Temp(PosHDInd,j,2,1)=Temp(PosHDInd,j,2,1)+1
                                                    end if
                                                    if(PhaseHD(i,j,Gam2)==1) then
                                                        !$OMP ATOMIC
                                                        Temp(PosHDInd,j,2,2)=Temp(PosHDInd,j,2,2)+1
                                                    end if
                                                endif
                                            enddo
                                        endif

                                    endif

                                endif
                            endif

                        enddo
                        !$OMP END PARALLEL DO

                        deallocate(BitPhaseHD)
                        deallocate(BitImputePhase)
                        deallocate(MissPhaseHD)
                        deallocate(MissImputePhase)
                    enddo
                enddo

                do e=1,2
                    do j=1,inputParams%nsnp
                        do i=1,ped%pedigreeSize- ped%nDummys
                            if (AnimalOn(i,e)==1) then
                                if (ImputePhase(i,j,e)==9) then
                                    ! Impute phase allele with the most significant code for that allele across haplotypes
                                    ! only if the other codification never happens
                                    if ((Temp(i,j,e,1)>inputParams%nAgreeInternalHapLibElim).and.(Temp(i,j,e,2)==0)) then
                                        ImputePhase(i,j,e)=0
                                    end if
                                    if ((Temp(i,j,e,1)==0).and.(Temp(i,j,e,2)>inputParams%nAgreeInternalHapLibElim)) then
                                        ImputePhase(i,j,e)=1
                                    end if
                                endif
                            endif
                        enddo
                    enddo
                enddo
                deallocate(Temp)

                ImputePhase(0,:,:)=9
                ImputeGenos(0,:)=9

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

                use HaplotypeBits

                use AlphaImputeSpecFileModule
                implicit none

                integer :: e,g,h,i,j,PedId,GamA,GamB,PosHDInd
                integer :: StartSnp,EndSnp,AnimalOn(ped%pedigreeSize,2)
                integer,allocatable,dimension (:,:,:) :: PhaseHD
                integer(kind=1),allocatable,dimension (:,:,:,:) :: Temp
                integer(kind=8), allocatable, dimension(:,:,:) :: BitPhaseHD, BitImputePhase, MissPhaseHD, MissImputePhase

                type(BitSection) :: Section
                integer :: numSections, curSection, curPos




                inputParams => defaultInput
                ped%nHd=(count(Setter(:)==1))

                allocate(Temp(ped%pedigreeSize,inputParams%nsnpraw,2,2))
                allocate(PhaseHD(ped%nHd,inputParams%nsnpraw,2))
                Temp=0
                AnimalOn=0


                do h=1,apResults%nResults



                    ! Get phase information

                    do g=1,size(apResults%results(h)%cores)
                        ! Initialize Start and End snps of the cores
                        StartSnp= apResults%results(h)%startIndexes(g)
                        EndSnp=apResults%results(h)%endIndexes(g)


                        Section = BitSection((EndSnp - StartSnp + 1), 64)
                        numSections = Section%numSections

                        allocate(BitPhaseHD(ped%nHd,numSections,2))
                        allocate(BitImputePhase(0:ped%pedigreeSize,numSections,2))
                        allocate(MissPhaseHD(ped%nHd,numSections,2))
                        allocate(MissImputePhase(0:ped%pedigreeSize,numSections,2))

                        BitPhaseHD = 0
                        MissPhaseHD = 0
                        BitImputePhase = 0
                        MissImputePhase = 0

                        do e=1,2
                            curSection = 1
                            curPos = 1

                            do j = StartSnp, EndSnp

                                do i = 1, ped%nHd
                                    select case (PhaseHD(i, j, e))
                                    case (1)
                                        BitPhaseHD(i, curSection, e) = ibset(BitPhaseHD(i, curSection, e), curPos)
                                    case (9)
                                        MissPhaseHD(i, curSection, e) = ibset(MissPhaseHD(i, curSection, e), curPos)
                                    end select
                                end do

                                do i = 1, ped%pedigreeSize- ped%nDummys
                                    select case (ImputePhase(i, j, e))
                                    case (1)
                                        BitImputePhase(i, curSection, e) = ibset(BitImputePhase(i, curSection, e), curPos)
                                    case (9)
                                        MissImputePhase(i, curSection, e) = ibset(MissImputePhase(i, curSection, e), curPos)
                                    end select
                                end do

                                curPos = curPos + 1
                                if (curPos == 65) then
                                    curPos = 1
                                    curSection = curSection + 1
                                end if
                            end do
                        end do

                        block
                            use individualModule
                            type(individual) ,pointer :: parent

                            !$OMP PARALLEL DO &
                            !$OMP DEFAULT(SHARED) &
                            !$OMP PRIVATE(i,j,e,PedId,PosHDInd,GamA,GamB,parent)
                            do i=1,ped%pedigreeSize- ped%nDummys

                                do e=1,2
                                    PedId=e+1
                                    if (ped%pedigree(i)%isDummyBasedOnIndex(pedId)) cycle !checked this one makes sense
                                    parent => ped%pedigree(i)%getSireDamObjectByIndex(pedId)

                                    ! parent => ped%pedigree((ped%pedigree(i)%getSireDamNewIDByIndex(pedID)))

                                    ! Skip if, in the case of sex chromosome, me and my parent are heterogametic
                                    if ((inputParams%sexopt==1).and.(ped%pedigree(i)%gender==inputParams%HetGameticStatus).and.&
                                        (parent%gender == inputParams%HetGameticStatus)) then
                                        cycle
                                    end if

                                    ! We look for gamete through those individuals that have parents with HD genotype information
                                    if (associated(parent)) then
                                        ! We look for possible gametes within the haplotypes identified to each of the
                                        ! individual's parents constructed during the phasing step
                                     ! Index of the parent in the HD phase information
                                        posHDInd = ped%hdDIctionary%getValue(parent%originalId)

                                        ! If there is one allele phased at least
                                        if ((Section%BitCountAllelesImputed(MissImputePhase(i,:,1)) + &
                                            Section%BitCountAllelesImputed(MissImputePhase(i,:,2))) > 0 .AND. PosHDInd/= DICT_NULL) then

                                            GamA=1
                                            GamB=1

                                            if (.NOT. Section%compareHaplotypeAllowMissing(BitPhaseHD(PosHDInd,:,1), BitImputePhase(i,:,e), &
                                                MissPhaseHD(PosHDInd,:,1), MissImputePhase(i,:,e), ImputeFromParentCountThresh)) then
                                                GamA = 0
                                            end if

                                            if (.NOT. Section%compareHaplotypeAllowMissing(BitPhaseHD(PosHDInd,:,2), BitImputePhase(i,:,e), &
                                                MissPhaseHD(PosHDInd,:,2), MissImputePhase(i,:,e), ImputeFromHDPhaseThresh)) then
                                                GamB = 0
                                            end if


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
                                                    if (ImputePhase(i,j,e)==9) then
                                                        if(PhaseHD(PosHDInd,j,1)==0) then
                                                            !$OMP ATOMIC
                                                            Temp(i,j,e,1)=Temp(i,j,e,1)+1
                                                        end if
                                                        if(PhaseHD(PosHDInd,j,1)==1) then
                                                            !$OMP ATOMIC
                                                            Temp(i,j,e,2)=Temp(i,j,e,2)+1
                                                        end if
                                                    endif
                                                enddo
                                            endif

                                            ! e is the individual's maternal haplotype
                                            if ((GamA==0).and.(GamB==1)) then
                                                ! Consider the haplotype of this animal in further steps
                                                AnimalOn(i,e)=1
                                                do j=StartSnp,EndSnp
                                                    ! Count the number of alleles coded with 0 and 1
                                                    if (ImputePhase(i,j,e)==9) then
                                                        if(PhaseHD(PosHDInd,j,2)==0) then
                                                            !$OMP ATOMIC
                                                            Temp(i,j,e,1)=Temp(i,j,e,1)+1
                                                        end if
                                                        if(PhaseHD(PosHDInd,j,2)==1) then
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
                        deallocate(BitPhaseHD)
                        deallocate(BitImputePhase)
                        deallocate(MissPhaseHD)
                        deallocate(MissImputePhase)

                    enddo
                enddo

                do e=1,2
                    do j=1,inputParams%nsnp
                        do i=1,ped%pedigreeSize- ped%nDummys
                            if (AnimalOn(i,e)==1) then
                                if ((inputParams%sexopt==0).or.(ped%pedigree(i)%gender==inputParams%HomGameticStatus)) then
                                    if (ImputePhase(i,j,e)==9) then
                                        ! Impute phase allele with the most significant code for that allele across haplotypes
                                        ! only if the other codification never happens
                                        if ((Temp(i,j,e,1)>inputParams%nAgreeInternalHapLibElim).and.(Temp(i,j,e,2)==0)) then
                                            ImputePhase(i,j,e)=0
                                        end if
                                        if ((Temp(i,j,e,1)==0).and.(Temp(i,j,e,2)>inputParams%nAgreeInternalHapLibElim)) then
                                            ImputePhase(i,j,e)=1
                                        end if
                                    end if
                                end if
                            end if
                        end do
                    end do
                end do

                do j=1,inputParams%nsnp
                    do i = 1, ped%pedigreeSize- ped%nDummys
                        if ((inputParams%sexopt==1).and.(ped%pedigree(i)%gender==inputParams%HetGameticStatus)) then
                            if (AnimalOn(i,inputParams%HomGameticStatus)==1) then

                                if (ImputePhase(i,j,inputParams%HomGameticStatus)==9) then
                                    ! Impute phase allele with the most significant code for that allele across haplotypes
                                    ! only if the other codification never happens
                                    if ((Temp(i,j,inputParams%HomGameticStatus,1)>inputParams%nAgreeInternalHapLibElim).and.&
                                        (Temp(i,j,inputParams%HomGameticStatus,2)==0)) then
                                        ImputePhase(i,j,:)=0
                                    end if
                                    if ((Temp(i,j,inputParams%HomGameticStatus,1)==0).and.&
                                        (Temp(i,j,inputParams%HomGameticStatus,2)>inputParams%nAgreeInternalHapLibElim)) then
                                        ImputePhase(i,j,:)=1
                                    end if
                                endif
                            end if
                        end if
                    enddo
                enddo
                deallocate(Temp)

                ImputePhase(0,:,:)=9
                ImputeGenos(0,:)=9

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

                use HaplotypeBits
                use AlphaImputeSpecFileModule
                use AlphaPhaseResultsModule
                implicit none


                integer :: i,j,k,h,e,f,g,CoreLength,nHap,CountAB(inputParams%nsnpraw,0:1),Work(inputParams%nsnpraw,2),TempCount
                integer :: StartSnp,EndSnp,PatMatDone(2),Counter,BanBoth(2),Ban(2),AnimalOn(ped%pedigreeSize,2)
                integer,allocatable,dimension (:,:,:,:) :: Temp
                integer(kind=1),allocatable,dimension (:,:) :: HapLib,HapCand

                integer(kind=8), allocatable, dimension(:,:,:) :: BitImputePhase, MissImputePhase
                integer(kind=8), allocatable, dimension(:,:) :: BitHapLib, MissHapLib


                type(BitSection) :: Section
                integer :: numSections, curSection, curPos
                integer :: z

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
                        Section = BitSection(CoreLength, 64)
                        numSections = Section%numSections

                        allocate(BitImputePhase(0:ped%pedigreeSize,numSections,2))
                        allocate(MissImputePhase(0:ped%pedigreeSize,numSections,2))

                        BitImputePhase = 0
                        MissImputePhase = 0

                        do e=1,2
                            curSection = 1
                            curPos = 1
                            do j = StartSnp, EndSnp
                                do i = 1, ped%pedigreeSize- ped%nDummys
                                    select case (ImputePhase(i, j, e))
                                    case (1)
                                        BitImputePhase(i, curSection, e) = ibset(BitImputePhase(i, curSection, e), curPos)
                                    case (9)
                                        MissImputePhase(i, curSection, e) = ibset(MissImputePhase(i, curSection, e), curPos)
                                    end select
                                end do

                                curPos = curPos + 1
                                if (curPos == 65) then
                                    curPos = 1
                                    curSection = curSection + 1
                                end if
                            end do
                        end do

                        nHap = apResults%results(h)%libraries(g)%size
                        coreLength = apResults%results(h)%libraries(g)%nsnps
                        allocate(HapLib(nHap, CoreLength))
                        do z=1, apResults%results(h)%libraries(g)%size
                            hapLib(z,:) = apResults%results(h)%libraries(g)%newStore(z)%toIntegerArray()
                        enddo


                        allocate(BitHapLib(nHap,numSections))
                        allocate(MissHapLib(nHap,numSections))

                        BitHapLib = 0
                        MissHapLib = 0

                        curSection = 1
                        curPos = 1
                        do j = 1, CoreLength
                            do i = 1, nHap
                                select case (HapLib(i, j))
                                case (1)
                                    BitHapLib(i, curSection) = ibset(BitHapLib(i, curSection), curPos)
                                case (9)
                                    MissHapLib(i, curSection) = ibset(MissHapLib(i, curSection), curPos)
                                end select
                            end do

                            curPos = curPos + 1
                            if (curPos == 65) then
                                curPos = 1
                                curSection = curSection + 1
                            end if
                        end do

                        !$OMP PARALLEL DO &
                        !$OMP DEFAULT(SHARED) &
                        !$OMP PRIVATE(i,e,f,j,k,TempCount,CountAB,Counter,PatMatDone,Work,BanBoth,Ban,HapCand)
                        do i=1,ped%pedigreeSize- ped%nDummys
                            ! The number of candidate haplotypes is the total number of haps in the library times
                            ! 2 (Paternal and maternal candidates)
                            allocate(HapCand(nHap,2))

                            PatMatDone=0
                            if (.not. Section%BitCompleteMissing(MissImputePhase(i,:,1))) then
                                PatMatDone(1) = 1
                            end if
                            if (.not. Section%BitCompleteMissing(MissImputePhase(i,:,2))) then
                                PatMatDone(2) = 1
                            end if

                            HapCand=1
                            Work=9
                            BanBoth=0

                            ! For each haplotype
                            do e=1,2
                                ! If GeneProbPhase has been executed, that is, if not considering the Sex Chromosome, then MSTermInfo={0,1}.
                                ! Else, if Sex Chromosome, then MSTermInfo is 0 always
                                ! So, if a Conservative imputation of haplotypes is selected, this DO statement will do nothing
                                if ((inputParams%ConservativeHapLibImputation==1).and.(MSTermInfo(i,e)==0)) cycle

                                ! If haplotype is partially phased
                                if ((PatMatDone(e)==1).and.&
                                    Section%BitCompletePhased(MissImputePhase(i,:,e)) == .FALSE.) then

                                    do f=1,nHap
                                        if (.not. Section%compareHaplotypeAllowMissing(BitHapLib(f,:), BitImputePhase(i,:,e), &
                                            MissHapLib(f,:), MissImputePhase(i,:,e))) then
                                            HapCand(f,e)=0
                                        end if
                                    enddo

                                    ! CountAB(inputParams%nsnps,0:1) matrix indicating how many haplotypes has been phased as 0
                                    ! or 1 in a particular allele
                                    CountAB=0

                                    ! If the number of candidate haplotypes is less than the 25% of the Library,
                                    ! then impute if all alleles have been phased the same way
                                    Counter=count(HapCand(:,e)==1)
                                    if (float(Counter)<(float(nHap)*0.25)) then
                                        ! Ban this haplotype will be phased here and nowhere else
                                        BanBoth(e)=1

                                        ! Count the occurrences in phasing of alleles across candidate haplotypes
                                        do f=1,nHap
                                            if (HapCand(f,e)==1) then
                                                k=0
                                                do j=StartSnp,EndSnp
                                                    k=k+1
                                                    ! Count occurrence of phase code HapLib(f,k)={0,1}
                                                    CountAB(j,HapLib(f,k))=CountAB(j,HapLib(f,k))+1
                                                enddo
                                            endif
                                        enddo

                                        ! If all alleles across the candidate haplotypes have been phased the same way, impute
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

                            ! If one of the haplotypes is partially phased
                            if (sum(PatMatDone(:))>0) then
                                Ban=0
                                ! Any haplotype has been previously banned/phased?
                                if (BanBoth(1)==1) Ban(1)=1
                                if (BanBoth(2)==1) Ban(2)=1

                                ! If both gametes have been previously banned/phased,
                                ! check whether the phase given agrees with genotype
                                if (sum(BanBoth(:))==2) then
                                    TempCount=0
                                    do j=StartSnp,EndSnp
                                        if (ImputeGenos(i,j)/=9) then
                                            ! If disagreement is greater than a threshold, unban haplotypes
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
                        deallocate(BitHapLib)
                        deallocate(MissHapLib)

                        deallocate(BitImputePhase)
                        deallocate(MissImputePhase)
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

                                if (ImputePhase(i,j,e)==9) then
                                    if ((Temp(i,j,e,1)>inputParams%nAgreeInternalHapLibElim).and.(Temp(i,j,e,2)==0)) then
                                        ImputePhase(i,j,e)=0
                                    end if
                                    if ((Temp(i,j,e,1)==0).and.(Temp(i,j,e,2)>inputParams%nAgreeInternalHapLibElim)) then
                                        ImputePhase(i,j,e)=1
                                    end if
                                endif
                            end if
                        end do
                    enddo
                enddo
                deallocate(Temp)

                ImputePhase(0,:,:)=9
                ImputeGenos(0,:)=9

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

                integer :: e,h,i,g,j,MiddleResult,MiddleResultNoShift,CoreLength,CountDisagree
                integer :: CompJump,StartSnp,EndSnp,UptoRightSnp,UptoLeftSnp,UpToCoreA,UpToCoreB,C1,C2,C3,C4,Recmb,CompLength,RL
                integer :: UpToSnp,StPt,EndPt,FillInSt,FillInEnd,hdAnimid
                ! integer,allocatable,dimension (:,:) :: CoreIndexA,CoreIndexB,AnimRecomb
                integer,allocatable,dimension (:,:) :: AnimRecomb
                integer,allocatable,dimension (:,:,:,:) :: PhaseHD
                integer :: middleCoreIndex,middleCoreIndexShift


                inputParams => defaultInput

                ! TODO This should be removed
                ! ped%nHd=(count(Setter(:)==1))

                allocate(PhaseHD(ped%nHd,inputParams%nsnpraw,2,2))     ! HIGH DENSITY PHASING: PhaseHD = (Animals, SNPs, Haplotypes, Nonshifted and Shifted phasing)
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

                phaseHD(:,:,:,1) = apResults%results(MiddleResult)%getFullPhaseIntArray()
                phaseHD(:,:,:,2) = apResults%results(MiddleResultNoShift)%getFullPhaseIntArray()

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
                    if (ped%pedigree(ped%hdMap(i))%founder) then
                        CountDisagree=0
                        ! Check if the two haplotypes are equal
                        do j=StartSnp,EndSnp
                            if (ImputePhase(ped%hdMap(i),j,1)/=ImputePhase(ped%hdMap(i),j,2)) then
                                CountDisagree=CountDisagree+1
                                if (CountDisagree>1) exit
                            endif
                        enddo
                        ! If haplotypes are equal
                        ! Impute High Density phase
                        if (CountDisagree==0) then
                            ImputePhase(ped%hdMap(i),StartSnp:EndSnp,1)=PhaseHD(i,StartSnp:EndSnp,1,1)
                            ImputePhase(ped%hdMap(i),StartSnp:EndSnp,2)=PhaseHD(i,StartSnp:EndSnp,2,1)
                        endif
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
                                    C1=0
                                    C2=0
                                    C3=0
                                    C4=0
                                    Recmb=1
                                    do j=StPt,EndPt
                                        ! NOTE: ImputePhase array is a HD phase data because the SNPs are within the MiddleResult core

                                        
                                        if (PhaseHD(hdAnimid,j,1,2)==ImputePhase(i,j,1)) C1=C1+1
                                        if (PhaseHD(hdAnimid,j,1,2)==ImputePhase(i,j,2)) C2=C2+1
                                        if (PhaseHD(hdAnimid,j,2,2)==ImputePhase(i,j,1)) C3=C3+1
                                        if (PhaseHD(hdAnimid,j,2,2)==ImputePhase(i,j,2)) C4=C4+1
                                    enddo

                                    ! If one haplotype is the same as the paternal, impute
                                    if ((CompLength==C1).and.(CompLength/=C3)) then
                                        ImputePhase(i,FillInSt:FillInEnd,1)=PhaseHD(hdAnimid,FillInSt:FillInEnd,1,2)
                                        Recmb=0
                                    endif

                                    ! If one haplotype is the same as the paternal, impute
                                    if ((CompLength/=C1).and.(CompLength==C3)) then
                                        ImputePhase(i,FillInSt:FillInEnd,1)=PhaseHD(hdAnimid,FillInSt:FillInEnd,2,2)
                                        Recmb=0
                                    endif

                                    ! If one haplotype is the same as the maternal, impute
                                    if ((CompLength==C2).and.(CompLength/=C4)) then
                                        ImputePhase(i,FillInSt:FillInEnd,2)=PhaseHD(hdAnimid,FillInSt:FillInEnd,1,2)
                                        Recmb=0
                                    endif

                                    ! If one haplotype is the same as the maternal, impute
                                    if ((CompLength/=C2).and.(CompLength==C4)) then
                                        ImputePhase(i,FillInSt:FillInEnd,2)=PhaseHD(hdAnimid,FillInSt:FillInEnd,2,2)
                                        Recmb=0
                                    endif

                                    ! There is bridge in the phase in the RL direction. Nothing more will be done
                                    if (Recmb==1) then
                                        AnimRecomb(i,RL)=1
                                    end if
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
                                    C1=0
                                    C2=0
                                    C3=0
                                    C4=0
                                    Recmb=1
                                    do j=StPt,EndPt
                                        if (PhaseHD(hdAnimid,j,1,1)==ImputePhase(i,j,1)) C1=C1+1
                                        if (PhaseHD(hdAnimid,j,1,1)==ImputePhase(i,j,2)) C2=C2+1
                                        if (PhaseHD(hdAnimid,j,2,1)==ImputePhase(i,j,1)) C3=C3+1
                                        if (PhaseHD(hdAnimid,j,2,1)==ImputePhase(i,j,2)) C4=C4+1
                                    enddo
                                    if ((CompLength==C1).and.(CompLength/=C3)) then
                                        ImputePhase(i,FillInSt:FillInEnd,1)=PhaseHD(hdAnimid,FillInSt:FillInEnd,1,1)
                                        Recmb=0
                                    endif
                                    if ((CompLength/=C1).and.(CompLength==C3)) then
                                        ImputePhase(i,FillInSt:FillInEnd,1)=PhaseHD(hdAnimid,FillInSt:FillInEnd,2,1)
                                        Recmb=0
                                    endif
                                    if ((CompLength==C2).and.(CompLength/=C4)) then
                                        ImputePhase(i,FillInSt:FillInEnd,2)=PhaseHD(hdAnimid,FillInSt:FillInEnd,1,1)
                                        Recmb=0
                                    endif
                                    if ((CompLength/=C2).and.(CompLength==C4)) then
                                        ImputePhase(i,FillInSt:FillInEnd,2)=PhaseHD(hdAnimid,FillInSt:FillInEnd,2,1)
                                        Recmb=0
                                    endif
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

                ImputePhase(0,:,:)=9
                ImputeGenos(0,:)=9

            END SUBROUTINE BaseAnimalFillIn

            !#############################################################################################################################################################################################################################

            subroutine InitialiseArrays
                ! Impute phase information for homozygous cases
                use AlphaImputeSpecFileModule
                use global, only :ped
                implicit none

                integer :: i,j
                ! ImputeGenos=9
                ImputePhase=9
                inputParams => defaultInput

                ! imputeGenos(:,:) = ped%getGenotypesAsArray()
                imputeGenos(:,:) = ped%getGenotypesAsArrayWitHMissing()
                do i=1, ped%pedigreeSize
                    do j=1,inputParams%nsnp
                        if (ImputeGenos(i,j)==0) ImputePhase(i,j,:)=0
                        if (ImputeGenos(i,j)==2) ImputePhase(i,j,:)=1
                    enddo
                enddo



            end subroutine InitialiseArrays

            !#############################################################################################################################################################################################################################

            subroutine GeneralFillIn
                ! This function implements the four Minor sub-steps explained in Hickey et al. (2012; Appendix A)
                use Global
                use informationModule

                implicit none

                call ParentHomoFill                     ! Minor sub-step 1. Parent Homozygous fill in
                call PhaseComplement                    ! Minor sub-step 2. Phase Complement
                call ImputeParentByProgenyComplement    ! Minor sub-step 3. Impute Parents from Progeny Complement
                call MakeGenotype                       ! Minor sub-step 4. Make Genotype
                ! if (TestVersion==1) call CurrentYield
                ! if (TestVersion==1) call Checker

                ImputePhase(0,:,:)=9
                ImputeGenos(0,:)=9

            end subroutine GeneralFillIn



            subroutine GeneralFillInInit
                ! This function implements the four Minor sub-steps explained in Hickey et al. (2012; Appendix A)
                use Global
                use informationModule

                implicit none

                call ParentHomoFill                     ! Minor sub-step 1. Parent Homozygous fill in
                call PhaseComplement                    ! Minor sub-step 2. Phase Complement
                call ImputeParentByProgenyComplement    ! Minor sub-step 3. Impute Parents from Progeny Complement
                call MakeGenotype                       ! Minor sub-step 4. Make Genotype
                ! if (TestVersion==1) call CurrentYield
                ! if (TestVersion==1) call Checker

                ImputePhase(0,:,:)=9
                ImputeGenos(0,:)=9

            end subroutine GeneralFillInInit
            !#############################################################################################################################################################################################################################

            SUBROUTINE EnsureHetGametic
                ! Impute phase to Y chromosome from X chromosome for heterogametic individuals
                use Global
                use AlphaImputeSpecFileModule
                implicit none

                integer :: i,j

                inputParams => defaultInput

                do j=1,inputParams%nsnp
                    do i=1,ped%pedigreeSize- ped%nDummys
                        if (ped%pedigree(i)%gender==inputParams%HetGameticStatus) then
                            if ((ImputePhase(i,j,1)==9).and.(ImputePhase(i,j,2)/=9)) then
                                ImputePhase(i,j,1)=ImputePhase(i,j,2)
                            end if
                            if ((ImputePhase(i,j,2)==9).and.(ImputePhase(i,j,1)/=9)) then
                                ImputePhase(i,j,2)=ImputePhase(i,j,1)
                            end if
                        end if
                    enddo
                enddo

            END SUBROUTINE EnsureHetGametic

            !#############################################################################################################################################################################################################################

            SUBROUTINE MakeGenotype
                ! Any individual that has a missing genotype information but has both alleles
                ! known, has its genotype filled in as the sum of the two alleles
                use Global
                use AlphaImputeSpecFileModule

                implicit none

                integer :: i,j

                inputParams => defaultInput
                do i=1,ped%pedigreeSize- ped%nDummys
                    do j=1,inputParams%nsnp         
                        if (ImputeGenos(i,j)==9) then
                            if ((ImputePhase(i,j,1)/=9).and.(ImputePhase(i,j,2)/=9)) then
                                ImputeGenos(i,j)=sum(ImputePhase(i,j,:))
                            end if
                        endif
                    enddo
                enddo

                ImputePhase(0,:,:)=9
                ImputeGenos(0,:)=9

            END SUBROUTINE MakeGenotype

            !#############################################################################################################################################################################################################################

            subroutine PhaseComplement
                ! If the genotype at a locus for an individual is known and one of its alleles has been determined
                ! then impute the missing allele as the complement of the genotype and the known phased allele
                use Global
                implicit none

                integer :: i,j

                do i=1,ped%pedigreeSize- ped%nDummys
                    do j=1,inputParams%nsnp 
                        if (ImputeGenos(i,j)/=9) then
                            ! TODO check if this is  correct
                            ! if (imputeGenos(i,j) /= 0) then
                                if ((ImputePhase(i,j,1)/=9).and.(ImputePhase(i,j,2)==9)) then
                                    ImputePhase(i,j,2)=ImputeGenos(i,j)-ImputePhase(i,j,1)
                                    if (ImputePhase(i,j,2) < 0) Then
                                        ImputePhase(i,j,2) = 0
                                    endif
                                end if
                                if ((ImputePhase(i,j,2)/=9).and.(ImputePhase(i,j,1)==9)) then
                                    ImputePhase(i,j,1)=ImputeGenos(i,j)-ImputePhase(i,j,2)
                                    if (ImputePhase(i,j,2) < 0) Then
                                        ImputePhase(i,j,2) = 0
                                    endif
                                end if
                            ! endif
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
                use AlphaImputeSpecFileModule
                implicit none

                integer :: e,i,j,ParId

                inputParams => defaultInput
                do i=1,ped%pedigreeSize- ped%nDummys
                    if (inputParams%sexopt==0 .or. (inputParams%sexopt==1 .and. ped%pedigree(i)%gender/=inputParams%HetGameticStatus) ) then     ! If individual is homogametic
                        do e=1,2
                            if (ped%pedigree(i)%isDummyBasedOnIndex(e+1)) cycle
                            ParId=ped%pedigree(i)%getSireDamNewIDByIndex(e+1)
                            do j=1,inputParams%nsnp
                                if (ImputePhase(i,j,e)==9) then                 ! Always that the SNP is not genotyped
                                    if ((ImputePhase(ParId,j,1)==ImputePhase(ParId,j,2)).and. &
                                        (ImputePhase(ParId,j,1)/=9)) then
                                        ! Imput phase if parent is homozygous
                                        ImputePhase(i,j,e)=ImputePhase(ParId,j,1)
                                    end if
                                endif
                            enddo
                        enddo
                    else
                        ParId= ped%pedigree(i)%getSireDamNewIDByIndex(inputParams%HomGameticStatus+1) !the homogametic parent
                        do j=1,inputParams%nsnp
                            if (ImputePhase(i,j,1)==9) then !Comment from John Hickey see analogous iterate subroutine
                                if ((ImputePhase(ParId,j,1)==ImputePhase(ParId,j,2)).and. &
                                    (ImputePhase(ParId,j,1)/=9)) then
                                    ! Imput phase to the two haplotypes if parent is homozygous
                                    ImputePhase(i,j,:)=ImputePhase(ParId,j,1)
                                end if
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
                use AlphaImputeSpecFileModule
                use individualModule
                implicit none

                integer :: i,k,l,Count1,Count0
                integer :: sireDam
                type(individual), pointer :: tmpChild
                inputParams => defaultInput
                ! TODO maybe change this to pedigreeSIze? 
                do i=1,ped%pedigreeSize- ped%nDummys
                    if (ped%pedigree(i)%nOffs /= 0) then       ! check that animal i,j is a sire or a dam

                        ! Sex chromosome
                        if (inputParams%sexopt==1) then
                            do k=1,inputParams%nsnp

                                ! Mat gamete missing -> fill if offspring suggest heterozygous
                                ! WARNING: This was comment the other way around in the original version of the code
                                if ((ImputePhase(i,k,1)/=9).and.(ImputePhase(i,k,2)==9)) then
                                    Count1=0
                                    Count0=0
                                    if (ImputePhase(i,k,1)==1) Count1=1
                                    if (ImputePhase(i,k,1)==0) Count0=1

                                    ! Look for the individual progeny and count their phase
                                    do l=1,ped%pedigree(i)%nOffs


                                        
                                        tmpChild => ped%pedigree(i)%offsprings(l)%p

                                        ! This is the only difference with the inputParams%sexopt=0 code below. Duplicating
                                        ! the code can be avoided by including the IF statement here instead than
                                        ! outside the SNPs loop.
                                        if ((ped%pedigree(i)%gender ==inputParams%HetGameticStatus).and.(tmpChild%gender==inputParams%HetGameticStatus)) cycle
                                        if (tmpChild%sirePointer == ped%pedigree(i)) then
                                                sireDam = 1
                                            else
                                                sireDam = 2
                                            endif
                                        if (ImputePhase(tmpChild%id,k,sireDam)==0) Count0=Count0+1
                                        if (ImputePhase(tmpChild%id,k,sireDam)==1) Count1=Count1+1

                                        if ((Count0>0).and.(Count1>0)) then                 !Consider increasing the number of offspring required to all ow for genotyping error
                                            if (ImputePhase(i,k,1)==0) ImputePhase(i,k,2)=1
                                            if (ImputePhase(i,k,1)==1) ImputePhase(i,k,2)=0
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

                                    do l=1,ped%pedigree(i)%nOffs

                                        tmpChild => ped%pedigree(i)%offsprings(l)%p
                                        ! This is the only difference with the inputParams%sexopt=0 code below. Duplicating
                                        ! the code can be avoided by including the IF statement here instead than
                                        ! outside the SNPs loop.
                                        if ((ped%pedigree(i)%gender ==inputParams%HetGameticStatus).and.(tmpChild%gender==inputParams%HetGameticStatus)) cycle
                                        
                                          if (tmpChild%sirePointer == ped%pedigree(i)) then
                                                sireDam = 1
                                            else
                                                sireDam = 2
                                            endif
                                        
                                        if (ImputePhase(tmpChild%id,k,sireDam)==0) Count0=Count0+1
                                        if (ImputePhase(tmpChild%id,k,sireDam)==1) Count1=Count1+1

                                        if ((Count0>0).and.(Count1>0)) then                 !Consider increasing the number of offspring required to all ow for genotyping error
                                            if (ImputePhase(i,k,2)==0) ImputePhase(i,k,1)=1
                                            if (ImputePhase(i,k,2)==1) ImputePhase(i,k,1)=0
                                        endif
                                    enddo
                                endif
                            enddo

                            ! Generic chromosome
                        else
                            do k=1,inputParams%nsnp
                                if ((ImputePhase(i,k,1)/=9).and.(ImputePhase(i,k,2)==9)) then               !Pat gamete missing fill if offspring suggest heterozygous
                                    Count1=0
                                    Count0=0
                                    if (ImputePhase(i,k,1)==1) Count1=1
                                    if (ImputePhase(i,k,1)==0) Count0=1
                                    

                                    do l=1,ped%pedigree(i)%nOffs
                                        tmpChild => ped%pedigree(i)%offsprings(l)%p

                                          if (tmpChild%sirePointer == ped%pedigree(i)) then
                                                sireDam = 1
                                            else
                                                sireDam = 2
                                            endif

                                        if (ImputePhase(tmpChild%id,k,sireDam)==0) Count0=Count0+1
                                        if (ImputePhase(tmpChild%id,k,sireDam)==1) Count1=Count1+1
                                        if ((Count0>0).and.(Count1>0)) then                 !Consider increasing the number of offspring required to all ow for genotyping error
                                            if (ImputePhase(i,k,1)==0) ImputePhase(i,k,2)=1
                                            if (ImputePhase(i,k,1)==1) ImputePhase(i,k,2)=0
                                        endif
                                    enddo
                                
                                ! TODO change to else 
                                else if ((ImputePhase(i,k,2)/=9).and.(ImputePhase(i,k,1)==9)) then               !Mat gamete missing fill if offspring suggest heterozygous
                                    Count1=0
                                    Count0=0
                                    if (ImputePhase(i,k,2)==1) Count1=1
                                    if (ImputePhase(i,k,2)==0) Count0=1
                                    do l=1,ped%pedigree(i)%nOffs
                                        tmpChild => ped%pedigree(i)%offsprings(l)%p
                                        if (tmpChild%sirePointer == ped%pedigree(i)) then
                                            sireDam = 1
                                        else
                                            sireDam = 2
                                        endif
                                        if (ImputePhase(Tmpchild%id,k,sireDam)==0) Count0=Count0+1
                                        if (ImputePhase(Tmpchild%id,k,sireDam)==1) Count1=Count1+1

                                        if ((Count0>0).and.(Count1>0)) then                 !Consider increasing the number of offspring required to all ow for genotyping error
                                            if (ImputePhase(i,k,2)==0) ImputePhase(i,k,1)=1
                                            if (ImputePhase(i,k,2)==1) ImputePhase(i,k,1)=0
                                        endif
                                    enddo
                                endif
                            enddo
                        endif
                    endif
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
                ! TODOphase this all needs redone
                use AlphaImputeSpecFileModule
                use AlphaImputeInputOutputModule, only: readProbabilitiesFull, readProbabilitiesFullCluster

                implicit none

                integer :: i,j,k,m,StSnp,EnSnp
                real :: PatAlleleProb(inputParams%nsnp,2),MatAlleleProb(inputParams%nsnp,2),HetProb(inputParams%nsnp)
                integer :: Informativeness(inputParams%nsnp,6),TmpInfor(inputParams%nsnp,6),GrandPar

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
                            if (PatAlleleProb(j,1)>=GeneProbThresh) ImputePhase(i,j,1)=0
                            if (PatAlleleProb(j,2)>=GeneProbThresh) ImputePhase(i,j,1)=1
                            if (MatAlleleProb(j,1)>=GeneProbThresh) ImputePhase(i,j,2)=0
                            if (MatAlleleProb(j,2)>=GeneProbThresh) ImputePhase(i,j,2)=1
                            if (HetProb(j)>=GeneProbThresh) ImputeGenos(i,j)=1
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

                            if (ImputeGenos(i,j)==1) then               ! If heterozygous

                                if (.not. ped%pedigree(i)%isDummyBasedOnIndex(2)) then
                                    ! My father is heterozygous
                                    if (ImputeGenos(ped%pedigree(i)%getSireDamNewIDByIndex(2),j)==1) then
                                        ! And have my father haplotype phased
                                        if ((ImputePhase(i,j,1)==0).or.(ImputePhase(i,j,1)==1)) then
                                            Informativeness(j,1)=1
                                            GlobalTmpCountInf(i,1)=GlobalTmpCountInf(i,1)+1     ! Count the number of SNPs phased of my father haplotype
                                            TmpInfor(GlobalTmpCountInf(i,1),1)=j                ! The SNP no. GlobalTmpCountInf(i,1) is my SNP no. j
                                        endif
                                    endif
                                endif

                                ! My mother is heterozygous
                                if (.not. ped%pedigree(i)%isDummyBasedOnIndex(3)) then
                                    if (ImputeGenos(ped%pedigree(i)%getSireDamNewIDByIndex(3),j)==1) then
                                        ! And have my mother haplotype phased
                                        if ((ImputePhase(i,j,2)==0).or.(ImputePhase(i,j,2)==1)) then
                                            Informativeness(j,2)=1
                                            GlobalTmpCountInf(i,2)=GlobalTmpCountInf(i,2)+1
                                            TmpInfor(GlobalTmpCountInf(i,2),2)=j
                                        endif
                                    endif
                                endif


                                ! My father haplotype is phased
                                if ((ImputePhase(i,j,1)==0).or.(ImputePhase(i,j,1)==1)) then
                                    ! If my paternal GranSire is heterozygous
                                    GrandPar=ped%pedigree(i)%getPaternalGrandSireRecodedIndexNoDummy()
                                    if (ImputeGenos(GrandPar,j)==1) then
                                        Informativeness(j,3)=1
                                        GlobalTmpCountInf(i,3)=GlobalTmpCountInf(i,3)+1
                                        TmpInfor(GlobalTmpCountInf(i,3),3)=j
                                    endif
                                    ! If my maternal GranDam is heterozygous
                                    GrandPar=ped%pedigree(i)%getPaternalGrandSireRecodedIndexNoDummy()
                                    if (ImputeGenos(GrandPar,j)==1) then
                                        Informativeness(j,4)=1
                                        GlobalTmpCountInf(i,4)=GlobalTmpCountInf(i,4)+1
                                        TmpInfor(GlobalTmpCountInf(i,4),4)=j
                                    endif
                                endif

                                ! My mother haplotype is phased
                                if ((ImputePhase(i,j,2)==0).or.(ImputePhase(i,j,2)==1)) then
                                    ! If my maternal GranSire is heterozygous
                                    ! if (ped%pedigree(i)%hasDummyParentsOrGranparents()) cycle
                                    ! TODO make this return 0
                                    GrandPar= ped%pedigree(i)%getPaternalGrandSireRecodedIndexNoDummy()

                                    if (ImputeGenos(GrandPar,j)==1) then
                                        Informativeness(j,5)=1
                                        GlobalTmpCountInf(i,5)=GlobalTmpCountInf(i,5)+1
                                        TmpInfor(GlobalTmpCountInf(i,5),5)=j
                                    endif
                                    ! If my maternal GranDam is heterozygous
                                    GrandPar=ped%pedigree(i)%getmaternalGrandDamRecodedIndexNoDummy()
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
                        write(102,'(a20,i2,i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7)') ped%pedigree(i)%originalID, Setter(i), GlobalTmpCountInf(i,k), TmpInfor(1:GlobalTmpCountInf(i,k),k)
                    enddo
                    write(102,'(a20,i2,i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7)') ped%pedigree(i)%originalID, Setter(i), GlobalTmpCountInf(i,7)
                    write(102,'(a20,i2,i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7,20000i7)') ped%pedigree(i)%originalID, Setter(i), GlobalTmpCountInf(i,8)
                    write(103,'(a20,2i10)') ped%pedigree(i)%id, MSTermInfo(i,:)
                enddo
                close(102)
                close(103)

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


                inputParams => defaultInput

                allocate(WorkRight(inputParams%nsnp))
                allocate(WorkLeft(inputParams%nsnp))
                allocate(TempVec(inputParams%nsnp))
                allocate(LengthVec(inputParams%nsnp))


                ImputePhase(0,:,:)=9
                ImputeGenos(0,:)=9

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
                                if (ped%pedigree(pedId)%isDummy) then
                                    cycle
                                endif
                            else
                                tmpGender = 0
                            endif
                            ! Skip if, in the case of sex chromosome, me and my parent are heterogametic
                            if ((inputParams%SexOpt==1).and.(ped%pedigree(i)%gender==inputParams%hetGameticStatus).and.(tmpGender==inputParams%hetGameticStatus)) cycle
                        endblock
                        !! SCAN HAPLOTYPE IN TWO DIRECTIONS: L->R AND R->L
                        ! If not a base animal and the number of unphased alleles is lower than a threshold
                        ! WARNING: WHAT IS THIS THRESHOLD?
                        if ((PedId>0).and.((float(count(ImputePhase(PedId,:,:)==9))/(2*inputParams%nsnp))<0.07)) then           !(RecIdHDIndex(PedId)==1)
                            WorkRight=9
                            RSide=9

                            ! Go throught haplotype from Left to Right
                            ! finding the first heterozygous allele of this parent, and...
                            do j=1,inputParams%nsnp
                                if ((ImputePhase(PedId,j,1)/=ImputePhase(PedId,j,2)).and.&
                                    (ImputePhase(PedId,j,1)/=9).and.(ImputePhase(PedId,j,2)/=9))  then
                                    HetStart=j
                                    ! Check if this allele corresponds to my parent paternal haplotype
                                    if (ImputePhase(i,HetStart,PatMat)==ImputePhase(PedId,HetStart,1)) then
                                        WorkRight(HetStart)=1   ! HetStart allele corresponds to Pat Haplotype in the LR direction
                                        RSide=1                 ! We are actually in the Paternal haplotype for the LR direction
                                        exit
                                    endif
                                    ! Check if this allele corresponds to my parent maternal haplotype
                                    if (ImputePhase(i,HetStart,PatMat)==ImputePhase(PedId,HetStart,2)) then
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
                                    if ((ImputePhase(i,j,PatMat)/=ImputePhase(PedId,j,RSide)).and.&
                                        (ImputePhase(PedId,j,RSide)/=9).and.(ImputePhase(i,j,PatMat)/=9)) then
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

                            ! ... Identifying recombinations
                            if (LSide/=9) then
                                do j=HetEnd-1,1,-1
                                    if ((ImputePhase(i,j,PatMat)/=ImputePhase(PedId,j,LSide)).and.&
                                        (ImputePhase(PedId,j,LSide)/=9).and.(ImputePhase(i,j,PatMat)/=9) ) then
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
                                if ( TempVec(j)==3 .AND. &
                                    (ImputePhase(PedId,j,1)/=ImputePhase(PedId,j,2)) .AND. &
                                    (ImputePhase(PedId,j,1)/=9) .and. (ImputePhase(PedId,j,2)/=9) ) then

                                    ImputePhase(i,j,e)=9
                                    if (Setter(i)/=1) then ! Skip HD individuals
                                        ImputeGenos(i,j)=9
                                    end if

                                endif
                            enddo

                            GlobalWorkPhase(i,:,:)=ImputePhase(i,:,:)

                            !$$$$$$$$$$$$$$$$$$$

                            !! IMPUTE PHASE WHETHER IT IS POSSIBLE
                            ! WARNING: From Hickey et al. 2012 (Appendix A):
                            !          ["Alleles are imputed ... subject to the restriction that the number of recombinations events for the individuals is less than a threshold, AND
                            !          that the region in which two recombination events occurred exceeds a threshold lenght."]
                            !          What it is coded is ["... than a threshold, OR that the region..."]
                            ! The number of recombinations in total (LR + RL) is less than a threshold
                            if ((CountLeftSwitch+CountRightSwitch)<(2*MaxLeftRightSwitch)) then
                                do j=1,inputParams%nsnp
                                    if (ImputePhase(i,j,PatMat)==9) then

                                        ! WARNING: This can be coded in a conciser way
                                        ! if ( (WorkLeft(j)/=9) .and. ( (WorkRight(j)==WorkLeft(j)).or.(WorkRight(j)==9) )  ) ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkLeft(j))
                                        ! if ( (WorkRight(j)/=9) .and. ( (WorkRight(j)==WorkLeft(j)).or.(WorkLeft(j)==9) )  ) ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkRight(j))
                                        ! Phase if the allele in one of the two directions is missing
                                        if ((WorkLeft(j)==9).and.(WorkRight(j)/=9)) &
                                            ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkRight(j))
                                        if ((WorkLeft(j)/=9).and.(WorkRight(j)==9)) &
                                            ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkLeft(j))

                                        ! Phase if alleles is the two directions agree
                                        if ((WorkLeft(j)/=9).and.(WorkRight(j)==WorkLeft(j))) &
                                            ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkLeft(j))
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
                                            if (ImputePhase(PedId,j,WorkRight(j))/=9)&
                                                ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkRight(j))
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
                        ImputePhase(i,:,inputParams%hetGameticStatus)=ImputePhase(i,:,inputParams%HomGameticStatus)     !JohnHickey changed the j to :
                        GlobalWorkPhase(i,:,:)=ImputePhase(i,:,:)
                    endif
                enddo

                ImputePhase(0,:,:)=9
                ImputeGenos(0,:)=9

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


                inputParams => defaultInput

                allocate(WorkRight(inputParams%nsnp))
                allocate(WorkLeft(inputParams%nsnp))
                allocate(TempVec(inputParams%nsnp))
                allocate(LengthVec(inputParams%nsnp))


                ImputePhase(0,:,:)=9
                ImputeGenos(0,:)=9

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
                                if ((ImputePhase(PedId,j,1)/=ImputePhase(PedId,j,2)).and.&
                                    (ImputePhase(PedId,j,1)/=9).and.(ImputePhase(PedId,j,2)/=9))  then
                                    HetStart=j
                                    ! Check if this allele corresponds to my parent paternal haplotype
                                    if (ImputePhase(i,HetStart,PatMat)==ImputePhase(PedId,HetStart,1)) then
                                        WorkRight(HetStart)=1   ! HetStart allele corresponds to Pat Haplotype in the LR direction
                                        RSide=1                 ! We are actually in the Paternal haplotype for the LR direction
                                        exit
                                    endif
                                    ! Check if this allele corresponds to my parent maternal haplotype
                                    if (ImputePhase(i,HetStart,PatMat)==ImputePhase(PedId,HetStart,2)) then
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
                                    if ((ImputePhase(i,j,PatMat)/=ImputePhase(PedId,j,RSide)).and.&
                                        (ImputePhase(PedId,j,RSide)/=9).and.(ImputePhase(i,j,PatMat)/=9)) then
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

                            ! ... Identifying recombinations
                            if (LSide/=9) then
                                do j=HetEnd-1,1,-1
                                    if ((ImputePhase(i,j,PatMat)/=ImputePhase(PedId,j,LSide)).and.&
                                        (ImputePhase(PedId,j,LSide)/=9).and.(ImputePhase(i,j,PatMat)/=9) ) then
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
                                if ( TempVec(j)==3 .AND. &
                                    (ImputePhase(PedId,j,1)/=ImputePhase(PedId,j,2)) .AND. &
                                    (ImputePhase(PedId,j,1)/=9) .and. (ImputePhase(PedId,j,2)/=9) ) then

                                    ImputePhase(i,j,e)=9
                                    if (Setter(i)/=1) then ! Skip HD individuals
                                        ImputeGenos(i,j)=9
                                    end if

                                endif
                            enddo

                            GlobalWorkPhase(i,:,:)=ImputePhase(i,:,:)

                            !$$$$$$$$$$$$$$$$$$$

                            !! IMPUTE PHASE WHETHER IT IS POSSIBLE
                            ! WARNING: From Hickey et al. 2012 (Appendix A):
                            !          ["Alleles are imputed ... subject to the restriction that the number of recombinations events for the individuals is less than a threshold, AND
                            !          that the region in which two recombination events occurred exceeds a threshold lenght."]
                            !          What it is coded is ["... than a threshold, OR that the region..."]
                            ! The number of recombinations in total (LR + RL) is less than a threshold
                            if ((CountLeftSwitch+CountRightSwitch)<(2*MaxLeftRightSwitch)) then
                                do j=1,inputParams%nsnp
                                    if (ImputePhase(i,j,PatMat)==9) then

                                        ! WARNING: This can be coded in a conciser way
                                        ! if ( (WorkLeft(j)/=9) .and. ( (WorkRight(j)==WorkLeft(j)).or.(WorkRight(j)==9) )  ) ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkLeft(j))
                                        ! if ( (WorkRight(j)/=9) .and. ( (WorkRight(j)==WorkLeft(j)).or.(WorkLeft(j)==9) )  ) ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkRight(j))
                                        ! Phase if the allele in one of the two directions is missing
                                        if ((WorkLeft(j)==9).and.(WorkRight(j)/=9)) &
                                            ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkRight(j))
                                        if ((WorkLeft(j)/=9).and.(WorkRight(j)==9)) &
                                            ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkLeft(j))

                                        ! Phase if alleles is the two directions agree
                                        if ((WorkLeft(j)/=9).and.(WorkRight(j)==WorkLeft(j))) &
                                            ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkLeft(j))
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
                                            if (ImputePhase(PedId,j,WorkRight(j))/=9)&
                                                ImputePhase(i,j,PatMat)=ImputePhase(PedId,j,WorkRight(j))
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
                        !ImputePhase(i,j,inputParams%hetGameticStatus)=ImputePhase(i,j,inputParams%HomGameticStatus)
                        ImputePhase(i,:,inputParams%hetGameticStatus)=ImputePhase(i,:,inputParams%HomGameticStatus)     !JohnHickey changed the j to :
                        GlobalWorkPhase(i,:,:)=ImputePhase(i,:,:)
                    endif
                enddo

                ImputePhase(0,:,:)=9
                ImputeGenos(0,:)=9

            end subroutine WorkLeftRight


        END MODULE Imputation
