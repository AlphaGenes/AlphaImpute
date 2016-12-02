#ifdef OS_UNIX
#DEFINE DASH "/"
#else
#DEFINE DASH "\"
#endif

MODULE Imputation
  use Global
  
  use GlobalVariablesHmmMaCH
  use AlphaImputeInMod
  use InputMod
  implicit none

  type(AlphaImputeInput), pointer :: inputParams
CONTAINS

  SUBROUTINE ImputationManagement
    use omp_lib

    integer :: i,j,loop
    character(len=150) :: timeOut

    inputParams => defaultInput


    ! WARNING: Need to discuss this part of code with John. Nonsense going on here!

    if (inputParams%HMMOption==RUN_HMM_ONLY) then ! Avoid any adulteration of genotypes with imputation subroutines

#ifdef DEBUG
      write(0,*) 'DEBUG: Allocate memory for genotypes and haplotypes'
#endif

      ! TODO: This hack avoid mem allocation problems with ImputeGenos and ImputePhase
      !       up in the code: at InsteadOfGeneProb in MakeFiles subroutine.
      !       Something has to be done with InsteadOfGeneProb cos' it is causing lots
      !       of problems!!
      if (allocated(ImputeGenos)) Then
          deallocate(ImputeGenos)
      endif
      allocate(ImputeGenos(0:nAnisP,inputParams%nsnp))
      if (allocated(ImputePhase)) Then
          deallocate(ImputePhase)
      endif
      allocate(ImputePhase(0:nAnisP,inputParams%nsnp,2))
      ImputePhase=9

      allocate(GlobalTmpCountInf(nAnisP,8))

      allocate(MSTermInfo(nAnisP,2))

#ifdef DEBUG
      write(0,*) 'DEBUG: Read Genotypes'
#endif

      call ped%addGenotypeInformation(inputParams%GenotypeFile,inputParams%nsnp,NanisG)



#ifdef DEBUG
      write(0,*) 'DEBUG: Call Mach'
#endif

      call MaCHController(inputParams%HMMOption)
      call FromHMM2ImputePhase

#ifdef DEBUG
      write(0,*) 'DEBUG: Mach Finished'
#endif

    else

      if (inputParams%RestartOption==4) then
        allocate(ImputeGenos(0:nAnisP,inputParams%nsnp))
        allocate(ImputePhase(0:nAnisP,inputParams%nsnp,2))
      else
        if (inputParams%sexopt==0) then
          ! Impute initial genotypes from calculated genotype probabilities
          if (inputParams%BypassGeneProb==0) then
            allocate(ImputeGenos(0:nAnisP,inputParams%nsnp))
            allocate(ImputePhase(0:nAnisP,inputParams%nsnp,2))
            allocate(GlobalWorkPhase(0:nAnisP,inputParams%nsnp,2))
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

        if (inputParams%NoPhasing==1) then
          ! Major sub-step 2 as explained in Hickey et al. (2012; Appendix A)
          call BaseAnimalFillIn

          ! Impute phase whenever a pre-phase file exists
          if (inputParams%PrePhased==1) call ReadInPrePhasedData

          ! Impute phase in the sex chromosome
          if (inputParams%sexopt==1) call EnsureHetGametic

          ! General imputation procedures
          call GeneralFillIn

          if (inputParams%HMMOption==RUN_HMM_PREPHASE) Then
            call MaCHController(inputParams%HMMOption)
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

              call InternalHapLibImputation           ! Major Sub-Step 6 (Hickey et al., 2012; Appendix A)
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
        endif

        if (inputParams%sexopt==1) then
          call EnsureHetGametic
        end if
        call GeneralFillIn

        if (inputParams%HMMOption==RUN_HMM_YES) Then
          call MaCHController(inputParams%HMMOption)
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
    integer :: LoopStart,Offset,AnimalOn(nAnisP,2)
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

    allocate(Temp(nAnisP,inputParams%nsnp,2,2))
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

          allocate(BitImputePhase(numSections,2,0:nAnisP))
          allocate(MissImputePhase(numSections,2,0:nAnisP))

          BitImputePhase = 0
          MissImputePhase = 0

          do i = 1, nAnisP
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
            do i=1,nAnisP
              
              do e=1,2

                  if (ped%pedigree(i)%isDummyBasedOnIndex(e+1)) cycle
                  parent => ped%pedigree(i)%getSireDamObjectByIndex(e+1)
                ! Skip if, in the case of sex chromosome, me and my parent are heterogametic
                if ((inputParams%sexopt==1).and.(ped%pedigree(i)%gender==HetGameticStatus).and.(parent%gender==HetGameticStatus)) cycle
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
                            if ( BTEST(BitImputePhase(curSection,2,parent%id), curPos) == .FALSE. .AND. &
                                 BTEST(MissImputePhase(curSection,2,parent%id), curPos) == .FALSE. ) then
                              !$OMP ATOMIC
                              Temp(i,j,e,1)=Temp(i,j,e,1)+1
                            end if

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
    do e=1,2
      do j=1,inputParams%nsnp
        do i=1,nAnisP
          if (AnimalOn(i,e)==1) then
              if ((inputParams%sexopt==0).or.(ped%pedigree(i)%gender==HomGameticStatus)) then

                if (ImputePhase(i,j,e)==9) then
                    if ((Temp(i,j,e,1)>inputParams%nAgreeInternalHapLibElim).and.(Temp(i,j,e,2)==0)) ImputePhase(i,j,e)=0
                    if ((Temp(i,j,e,1)==0).and.(Temp(i,j,e,2)>inputParams%nAgreeInternalHapLibElim)) ImputePhase(i,j,e)=1
                endif
              end if
          endif
        end do
      enddo
    enddo

    ! The individual is has one haplotype: In Sex Chromosome, the heterogametic case
    do j=1,inputParams%nsnp
      do i =1, nAnisP
        if ((inputParams%sexopt==1).and.(ped%pedigree(i)%gender==HetGameticStatus)) then
          if (AnimalOn(i,HomGameticStatus)==1) then
            if (ImputePhase(i,j,HomGameticStatus)==9) then
                if ((Temp(i,j,HomGameticStatus,1)>inputParams%nAgreeInternalHapLibElim).and.(Temp(i,j,HomGameticStatus,2)==0))&
                    ImputePhase(i,j,:)=0
                if ((Temp(i,j,HomGameticStatus,1)==0).and.(Temp(i,j,HomGameticStatus,2)>inputParams%nAgreeInternalHapLibElim))&
                    ImputePhase(i,j,:)=1
            endif
          end if
        end if
      end do
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
    use alphaimputeinmod
    implicit none

    integer :: f,e,h,g,i,j,l,nCore,nHap,nGlobalLoop,CoreLength,CoreStart,CoreEnd,InLib,Count0,Count1
    integer :: Counter,BanBoth(2),Ban(2),AnimalOn(nAnisP,2)
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

    allocate(Temp(nAnisP,inputParams%nsnp,2,2))
    Temp=0
    AnimalOn=0

    ! SIMULATE PHASING
    ! f is a variable to simulate shift or no-shift phasing
    do f=1,2
      ! Allocate the Internal Haplotype Library
      allocate(HapLib(nAnisP*2,inputParams%nsnp))
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

          allocate(BitHapLib(nAnisP*2,numSections))
          allocate(MissHapLib(nAnisP*2,numSections))
          allocate(BitImputePhase(0:nAnisP,numSections,2))
          allocate(MissImputePhase(0:nAnisP,numSections,2))

          BitHapLib = 0
          MissHapLib = 0
          BitImputePhase = 0
          MissImputePhase = 0

          curSection = 1
          curPos = 1
          do j = CoreStart, CoreEnd
            do e = 1, 2
              do i = 1, nAnisP
                select case (ImputePhase(i, j, e))
                  case (1)
                    BitImputePhase(i, curSection, e) = ibset(BitImputePhase(i, curSection, e), curPos)
                  case (9)
                    MissImputePhase(i, curSection, e) = ibset(MissImputePhase(i, curSection, e), curPos)
                end select
             end do
            end do

            do i = 1, nAnisP*2
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
          do i=1,nAnisP
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

          !$OMP PARALLEL DO &
          !$OMP DEFAULT(SHARED) &
          !$OMP PRIVATE(i,j,e,h,HapElim,BanBoth,Counter,Count0,Count1,Ban,curPos,curSection,BitWork,MissWork,BitGeno)
          do i=1,nAnisP
            allocate(BitWork(numSections,2))
            allocate(MissWork(numSections,2))
            allocate(HapElim(nAnisP*2,2))
            HapElim=1
            BanBoth=0
            BitWork = 0
            MissWork = 0
            do e=1,2
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
                    HapElim(h,e)=0
                  end if
                enddo

                ! If the number of candidate haplotypes is less than the 25% of the Library,
                ! then impute if all alleles have been phased the same way
                Counter=count(HapElim(1:nHap,e)==1)
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
                      if (HapElim(h,e)==1) then
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
            deallocate(HapElim)
            deallocate(BitWork)
            deallocate(MissWork)
          enddo
          !$OMP END PARALLEL DO


          ! Prepare the core for the next cycle
          CoreStart=CoreStart+LoopIndex(l,2)
          CoreEnd=CoreEnd+LoopIndex(l,2)
          if ((f==2).and.(g==nCore)) exit
        enddo
      enddo
      deallocate(HapLib)
    enddo

    !$OMP PARALLEL DEFAULT(SHARED)
    do e=1,2
      !$OMP DO PRIVATE(i,j)
      do j=1,inputParams%nsnp
        do i=1,nAnisP
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
      !$OMP END DO
    enddo
    !$OMP END PARALLEL

    deallocate(Temp)

    ImputePhase(0,:,:)=9
    ImputeGenos(0,:)=9

  END SUBROUTINE InternalHapLibImputation

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
    
    use PhaseRounds
    use Utils
    use HaplotypeBits
    use alphaimputeinmod
    implicit none

    integer :: e,g,h,i,j,GamA,GamB,nAnisHD,PosHDInd
    integer :: StartSnp,EndSnp,Gam1,Gam2,AnimalOn(nAnisP,2)
    integer,allocatable,dimension (:) :: PosHD
    integer,allocatable,dimension (:,:,:) :: PhaseHD
    integer(kind=8), allocatable, dimension(:,:,:) :: BitPhaseHD, BitImputePhase, MissPhaseHD, MissImputePhase
    integer(kind=1),allocatable,dimension (:,:,:,:) :: Temp

    character(len=1000) :: FileName,FileNamePhase
    type(CoreIndex) :: CoreI
    type(BitSection) :: Section

    integer :: numSections, curSection, curPos

    inputParams => defaultInput
    ! Number of animals that have been HD phased
    nAnisHD=(count(Setter(:)==1))

    allocate(Temp(nAnisP,inputParams%nsnp,2,2))
    allocate(PhaseHD(nAnisHD,inputParams%nsnp,2))
    allocate(PosHD(nAnisP))
    PosHD=0
    Temp=0
    AnimalOn=0
    ! FOR EACH CORE OF EACH ROUND OF THE LRPHLI
    do h=1,inputParams%nPhaseInternal
      ! Get HIGH DENSITY phase information of this phasing step and information
      ! of core indexes
      FileName=""
      FileNamePhase=""
      if (inputParams%ManagePhaseOn1Off0==0) then
        FileName = getFileNameCoreIndex(trim(inputParams%phasePath),h)
        FileNamePhase = getFileNameFinalPhase(trim(inputParams%phasePath),h)
      else
        FileName = getFileNameCoreIndex(h)
        FileNamePhase = getFileNameFinalPhase(h)
      end if

      ! Get core information of number of cores and allocate start and end cores information
      CoreI = ReadCores(FileName)

      ! Get phase information
      call ReadPhased(nAnisHD, FileNamePhase, Ped, PhaseHD, PosHD)

      do g=1,CoreI%nCores
        ! Initialize Start and End snps of the cores
        StartSnp=CoreI%StartSnp(g)
        EndSnp=CoreI%EndSnp(g)

        Section = BitSection((EndSnp - StartSnp + 1), 64)
        numSections = Section%numSections

        allocate(BitPhaseHD(nAnisHD,numSections,2))
        allocate(BitImputePhase(0:nAnisP,numSections,2))
        allocate(MissPhaseHD(nAnisHD,numSections,2))
        allocate(MissImputePhase(0:nAnisP,numSections,2))

        BitPhaseHD = 0
        MissPhaseHD = 0
        BitImputePhase = 0
        MissImputePhase = 0

        do e=1,2
          curSection = 1
          curPos = 1

          do j = StartSnp, EndSnp

            do i = 1, nAnisHD
              select case (PhaseHD(i, j, e))
                case (1)
                  BitPhaseHD(i, curSection, e) = ibset(BitPhaseHD(i, curSection, e), curPos)
                case (9)
                  MissPhaseHD(i, curSection, e) = ibset(MissPhaseHD(i, curSection, e), curPos)
              end select
            end do

            do i = 1, nAnisP
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
        !$OMP PRIVATE(i,j,e,PosHDInd,Gam1,Gam2,GamA,GamB)
        do i=1,nAnisP
          ! Look for possible gametes through the Haplotype
          ! Library constructed during the phasing step
          PosHDInd=PosHD(i)         ! Index of the individual in the HD phase information

          ! If there is one allele phased at least
          if ((Section%BitCountAllelesImputed(MissImputePhase(i,:,1)) + &
               Section%BitCountAllelesImputed(MissImputePhase(i,:,2))) > 0 .AND. PosHDInd>0) then
            ! If at least one locus is heterozygous
              if (.NOT. Section%compareHaplotype(BitImputePhase(i,:,1), BitImputePhase(i,:,2), &
                  MissImputePhase(i,:,1), MissImputePhase(i,:,2))) then
              Gam1=0
              Gam2=0
              do e=1,2
                GamA=1
                GamB=1

                if (.NOT. Section%compareHaplotypeAllowMissing(BitPhaseHD(PosHDInd,:,1), BitImputePhase(i,:,e), &
                    MissPhaseHD(PosHDInd,:,1), MissImputePhase(i,:,e), ImputeFromHDPhaseThresh)) then
                  GamA = 0
                end if

                if (.NOT. Section%compareHaplotypeAllowMissing(BitPhaseHD(PosHDInd,:,2), BitImputePhase(i,:,e), &
                    MissPhaseHD(PosHDInd,:,2), MissImputePhase(i,:,e), ImputeFromHDPhaseThresh)) then
                  GamB = 0
                end if

                ! Paternal haplotype is strictly my paternal haplotype from the Hap Library
                if ((e==1).and.(GamA==1).and.(GamB==0)) Gam1=1
                ! Paternal haplotype is strictly my maternal haplotype from the Hap Library
                if ((e==1).and.(GamA==0).and.(GamB==1)) Gam1=2
                ! Maternal haplotype is strictly my paternal haplotype from the Hap Library
                if ((e==2).and.(GamA==1).and.(GamB==0)) Gam2=1
                ! Maternal haplotype is strictly my maternal haplotype from the Hap Library
                if ((e==2).and.(GamA==0).and.(GamB==1)) Gam2=2

                ! Basically the important thing is that haplotype e is present in the Haplotype
                ! library. It is not important which haplotype it is, whether the paternal or the
                ! maternal.
              enddo

              ! If the paternal and maternal gametes are different
              if (Gam1/=Gam2) then
                AnimalOn(i,:)=1             ! Consider this animal in further steps

                ! Paternal gamete is in the Hap Library
                if (Gam1/=0) then
                  do j=StartSnp,EndSnp
                    ! Count the number of alleles coded with 0 and 1
                    if (ImputePhase(i,j,1)==9) then
                      if(PhaseHD(PosHDInd,j,Gam1)==0) then
                        !$OMP ATOMIC
                        Temp(i,j,1,1)=Temp(i,j,1,1)+1
                      end if
                      if(PhaseHD(PosHDInd,j,Gam1)==1) then
                        !$OMP ATOMIC
                        Temp(i,j,1,2)=Temp(i,j,1,2)+1
                      end if
                    endif
                  enddo
                endif

                ! Maternal gamete is in the Hap Library
                if (Gam2/=0) then
                  do j=StartSnp,EndSnp
                    ! Count the number of alleles coded with 0 and 1
                    if (ImputePhase(i,j,2)==9) then
                      if(PhaseHD(PosHDInd,j,Gam2)==0) then
                        !$OMP ATOMIC
                        Temp(i,j,2,1)=Temp(i,j,2,1)+1
                      end if
                      if(PhaseHD(PosHDInd,j,Gam2)==1) then
                        !$OMP ATOMIC
                        Temp(i,j,2,2)=Temp(i,j,2,2)+1
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
        do i=1,nAnisP
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
    
    use PhaseRounds
    use HaplotypeBits
    use Utils
    use alphaimputeinmod
    implicit none

    integer :: e,g,h,i,j,PedId,GamA,GamB,nAnisHD,PosHDInd
    integer :: StartSnp,EndSnp,AnimalOn(nAnisP,2)
    integer,allocatable,dimension (:) :: PosHD
    integer,allocatable,dimension (:,:,:) :: PhaseHD
    integer(kind=1),allocatable,dimension (:,:,:,:) :: Temp
    integer(kind=8), allocatable, dimension(:,:,:) :: BitPhaseHD, BitImputePhase, MissPhaseHD, MissImputePhase

    character(len=1000) :: FileName,FileNamePhase
    type(CoreIndex) :: CoreI
    type(BitSection) :: Section
    integer :: numSections, curSection, curPos

    inputParams => defaultInput
    nAnisHD=(count(Setter(:)==1))

    allocate(Temp(nAnisP,inputParams%nsnp,2,2))
    allocate(PhaseHD(nAnisHD,inputParams%nsnp,2))
    allocate(PosHD(nAnisP))
    PosHD=0
    Temp=0
    AnimalOn=0

    do h=1,inputParams%nPhaseInternal
      ! Get HIGH DENSITY phase information of this phasing step and information
      ! of core indexes
      if (inputParams%ManagePhaseOn1Off0==0) then
        FileName = getFileNameCoreIndex(trim(inputParams%phasePath),h)
        FileNamePhase = getFileNameFinalPhase(trim(inputParams%phasePath),h)
      else
        FileName = getFileNameCoreIndex(h)
        FileNamePhase = getFileNameFinalPhase(h)
      end if

      ! Get core information of number of cores and allocate start and end cores information
      CoreI = ReadCores(FileName)

      ! Get phase information
      call ReadPhased(nAnisHD, FileNamePhase, ped, PhaseHD, PosHD)

      do g=1,CoreI%nCores
        ! Initialize Start and End snps of the cores
        StartSnp=CoreI%StartSnp(g)
        EndSnp=CoreI%EndSnp(g)

        Section = BitSection((EndSnp - StartSnp + 1), 64)
        numSections = Section%numSections

        allocate(BitPhaseHD(nAnisHD,numSections,2))
        allocate(BitImputePhase(0:nAnisP,numSections,2))
        allocate(MissPhaseHD(nAnisHD,numSections,2))
        allocate(MissImputePhase(0:nAnisP,numSections,2))

        BitPhaseHD = 0
        MissPhaseHD = 0
        BitImputePhase = 0
        MissImputePhase = 0

        do e=1,2
          curSection = 1
          curPos = 1

          do j = StartSnp, EndSnp

            do i = 1, nAnisHD
              select case (PhaseHD(i, j, e))
                case (1)
                  BitPhaseHD(i, curSection, e) = ibset(BitPhaseHD(i, curSection, e), curPos)
                case (9)
                  MissPhaseHD(i, curSection, e) = ibset(MissPhaseHD(i, curSection, e), curPos)
              end select
            end do

            do i = 1, nAnisP
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
          do i=1,nAnisP
            
            do e=1,2
              PedId=e+1
              if (ped%pedigree(i)%isDummyBasedOnIndex(pedId)) cycle !checked this one makes sense
              parent => ped%pedigree(i)%getSireDamObjectByIndex(pedId)
              ! Skip if, in the case of sex chromosome, me and my parent are heterogametic
              if ((inputParams%sexopt==1).and.(ped%pedigree(i)%gender==HetGameticStatus).and.&
                (parent%gender == HetGameticStatus)) then
                    cycle
              end if

              ! We look for gamete through those individuals that have parents with HD genotype information
              if (associated(parent)) then
                ! We look for possible gametes within the haplotypes identified to each of the
                ! individual's parents constructed during the phasing step
                PosHDInd=PosHD(parent%id)   ! Index of the parent in the HD phase information 
                !TODO check- SHOULD THIS BE INDIVIDUAL id (above) rather than parent ID?

                ! If there is one allele phased at least
                if ((Section%BitCountAllelesImputed(MissImputePhase(i,:,1)) + &
                     Section%BitCountAllelesImputed(MissImputePhase(i,:,2))) > 0 .AND. PosHDInd>0) then

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
        do i=1,nAnisP
          if (AnimalOn(i,e)==1) then
            if ((inputParams%sexopt==0).or.(ped%pedigree(i)%gender==HomGameticStatus)) then
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
      do i = 1, nAnisP
        if ((inputParams%sexopt==1).and.(ped%pedigree(i)%gender==HetGameticStatus)) then
          if (AnimalOn(i,HomGameticStatus)==1) then

            if (ImputePhase(i,j,HomGameticStatus)==9) then
              ! Impute phase allele with the most significant code for that allele across haplotypes
              ! only if the other codification never happens
              if ((Temp(i,j,HomGameticStatus,1)>inputParams%nAgreeInternalHapLibElim).and.&
                  (Temp(i,j,HomGameticStatus,2)==0)) then
                    ImputePhase(i,j,:)=0
              end if
              if ((Temp(i,j,HomGameticStatus,1)==0).and.&
                  (Temp(i,j,HomGameticStatus,2)>inputParams%nAgreeInternalHapLibElim)) then
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
    use PhaseRounds
    use Utils
    use HaplotypeBits
    use alphaimputeinmod
    implicit none


    integer :: l,i,j,k,h,e,f,g,CoreLength,nHap,CountAB(inputParams%nsnp,0:1),Work(inputParams%nsnp,2),TempCount
    integer :: StartSnp,EndSnp,PatMatDone(2),Counter,BanBoth(2),Ban(2),AnimalOn(nAnisP,2)
    integer,allocatable,dimension (:,:,:,:) :: Temp
    integer(kind=1),allocatable,dimension (:,:) :: HapLib,HapCand

    integer(kind=8), allocatable, dimension(:,:,:) :: BitImputePhase, MissImputePhase
    integer(kind=8), allocatable, dimension(:,:) :: BitHapLib, MissHapLib

    integer :: UHLib
    character(len=1000) :: FileName
    type(CoreIndex) :: CoreI

    type(BitSection) :: Section
    integer :: numSections, curSection, curPos

    inputParams => defaultInput
    ! Temp(nAnisP, inputParams%nsnps, PatHap, Phase)
    allocate(Temp(0:nAnisP,inputParams%nsnp,2,2))
    Temp=0

    AnimalOn=0

    do h=1,inputParams%nPhaseInternal
      ! Get HIGH DENSITY phase information of this phasing step and information
      ! of core indexes
      if (inputParams%ManagePhaseOn1Off0==0) then
        FileName = getFileNameCoreIndex(trim(inputParams%phasePath),h)
      else
        FileName = getFileNameCoreIndex(h)
      end if

      ! Get core information of number of cores and allocate start and end cores information
      CoreI = ReadCores(FileName)

      do g=1,CoreI%nCores
        ! Initialize Start and End snps of the cores
        StartSnp=CoreI%StartSnp(g)
        EndSnp=CoreI%EndSnp(g)

        CoreLength=(EndSnp-StartSnp)+1
        Section = BitSection(CoreLength, 64)
        numSections = Section%numSections

        allocate(BitImputePhase(0:nAnisP,numSections,2))
        allocate(MissImputePhase(0:nAnisP,numSections,2))

        BitImputePhase = 0
        MissImputePhase = 0

        do e=1,2
          curSection = 1
          curPos = 1
          do j = StartSnp, EndSnp
            do i = 1, nAnisP
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

        if (inputParams%ManagePhaseOn1Off0==0) then
          FileName = getFileNameHapLib(trim(inputParams%phasePath),h,g)
        else
          FileName = getFileNameHapLib(h,g)
        end if
        open (newunit=UHLib,file=trim(FileName),status="old",form="unformatted")

        ! Read the number of Hap in the library and how long they are
        read(UHLib) nHap,CoreLength

        if(nHap/=0) then
          ! Allocate the Haplotype Library and read it from file
          allocate(HapLib(nHap, CoreLength))
          do l=1,nHap
            read(UHLib) HapLib(l,:)
          enddo
          close (UHLib)

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
          do i=1,nAnisP
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
              ! internal phasing steps
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
        endif
        deallocate(BitImputePhase)
        deallocate(MissImputePhase)
      enddo
    enddo
    
    
      do j=1,inputParams%nsnp
               ! If GeneProbPhase has been executed, that is, if not considering the Sex Chromosome, then MSTermInfo={0,1}.
        ! Else, if Sex Chromosome, then MSTermInfo is 0 always
        ! So, if a Conservative imputation of haplotypes is selected, this DO statement will do nothing
  
          do i=1,nAnisP
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
    
    use PhaseRounds
    use Utils
    use alphaimputeinmod

    implicit none

    integer :: e,h,i,g,j,MiddlePhaseRun,MiddleCoreA,MiddleCoreB,CoreLength,nAnisHD,CountDisagree
    integer :: CompPhaseRun,CompJump,StartSnp,EndSnp,UptoRightSnp,UptoLeftSnp,UpToCoreA,UpToCoreB,C1,C2,C3,C4,Recmb,CompLength,RL
    integer :: UpToSnp,StPt,EndPt,FillInSt,FillInEnd
    integer,allocatable,dimension (:) :: PosHD
    ! integer,allocatable,dimension (:,:) :: CoreIndexA,CoreIndexB,AnimRecomb
    integer,allocatable,dimension (:,:) :: AnimRecomb
    integer,allocatable,dimension (:,:,:,:) :: PhaseHD
    character(len=1000) :: FileName,dumC
    type(CoreIndex) :: CoreIA, CoreIB

    inputParams => defaultInput
    nAnisHD=(count(Setter(:)==1))

    allocate(PhaseHD(nAnisHD,inputParams%nsnp,2,2))     ! HIGH DENSITY PHASING: PhaseHD = (Animals, SNPs, Haplotypes, Nonshifted and Shifted phasing)
    allocate(AnimRecomb(nAnisP,2))
    allocate(PosHD(nAnisP))
    AnimRecomb=0
    PosHD=0

    ! Select the phasing in the middle
    MiddlePhaseRun=(int(inputParams%nPhaseInternal)/4)
    if (MiddlePhaseRun==0) then
      MiddlePhaseRun=1
    end if
    CompJump=int(inputParams%nPhaseInternal)/2
    if (CompJump==0) then
      CompJump=1             ! This is only possible if no phasing info is available
    end if
    CompPhaseRun=MiddlePhaseRun+CompJump

    ! The different "PhasingResults/CoreIndex.txt" files are created
    ! by AlphaPhase when it is called during PhasingManagement

    ! Get HIGH DENSITY phase information of this phasing step and information
    ! of core indexes
    if (inputParams%ManagePhaseOn1Off0==0) then
      FileName = getFileNameCoreIndex(trim(inputParams%phasePath),MiddlePhaseRun)
    else
      FileName = getFileNameCoreIndex(MiddlePhaseRun)
    end if

    ! Get core information of number of cores and allocate start and end cores information
    CoreIA = ReadCores(FileName)

    ! Select the core in the middle
    MiddleCoreA=int(CoreIA%nCores)/2
    if (MiddleCoreA==0) MiddleCoreA=1

    ! Get HIGH DENSITY phase information of this phasing step and information
    ! of core indexes
    if (inputParams%ManagePhaseOn1Off0==0) then
      FileName = getFileNameCoreIndex(trim(inputParams%phasePath),CompPhaseRun)
    else
      FileName = getFileNameCoreIndex(CompPhaseRun)
    end if

    ! Get core information of number of cores and allocate start and end cores information
    CoreIB = ReadCores(FileName)

    ! Select the core in the middle
    MiddleCoreB=int(CoreIB%nCores)/2
    if (MiddleCoreB==0) MiddleCoreB=1

    ! Get HIGH DENSITY phase information of this phasing step
    ! WARNING: If I only want to phase base animals, why do I need to read the whole file?
#ifdef OS_UNIX
if (inputParams%ManagePhaseOn1Off0==0) then
    write (FileName,'(a,"Phase",i0,"/PhasingResults/FinalPhase.txt")') trim(inputParams%phasePath),MiddlePhaseRun
else
    write (FileName,'("./Phasing/Phase",i0,"/PhasingResults/FinalPhase.txt")')MiddlePhaseRun
endif
#else
if (inputParams%ManagePhaseOn1Off0==0) then
    write (FileName,'(a,"Phase",i0,"\PhasingResults\FinalPhase.txt")') trim(inputParams%phasePath),MiddlePhaseRun
else
    write (FileName,'(".\Phasing\Phase",i0,"\PhasingResults\FinalPhase.txt")')MiddlePhaseRun
endif
#endif
    block
      integer :: tmpIndex
      open (unit=2001,file=trim(FileName),status="old")
      do i=1,nAnisHD
        read (2001,*) dumC,PhaseHD(i,:,1,1)
        read (2001,*) dumC,PhaseHD(i,:,2,1)

        tmpIndex = ped%dictionary%getValue(dumC)
        if (tmpIndex /= DICT_NULL) then
          PosHD(tmpIndex)=i        

        else 
          write(error_unit,*) "WARNING - HD Phase information for animal not in the pedigree"
        endif 
      enddo
    endblock
    close(2001)

#ifdef OS_UNIX
if (inputParams%ManagePhaseOn1Off0==0) then
    write (FileName,'(a,"Phase",i0,"/PhasingResults/FinalPhase.txt")') trim(inputParams%phasePath),CompPhaseRun
else
    write (FileName,'("./Phasing/Phase",i0,"/PhasingResults/FinalPhase.txt")')CompPhaseRun
endif
#else
if (inputParams%ManagePhaseOn1Off0==0) then
    write (FileName,'(a,"Phase",i0,"\PhasingResults\FinalPhase.txt")') trim(inputParams%phasePath),CompPhaseRun
else
    write (FileName,'(".\Phasing\Phase",i0,"\PhasingResults\FinalPhase.txt")')CompPhaseRun
endif
#endif

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
      StartSnp=CoreIA%StartSnp(g)
      EndSnp=CoreIA%EndSnp(g)
      CoreLength=(EndSnp-StartSnp)+1
      do i=1,nAnisP
        ! If I have no parents and if I am somebody
        if (ped%pedigree(i)%founder .and.(PosHD(i)/=0)) then
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
        if ((CoreIA%nCores==1).and.(CoreIB%nCores==1)) exit   ! If the number of cores is 1, EXIT
                                                ! This will force the subroutine to finish
                                                ! since it will be exit from both DO statements
        h=h+1
        if (mod(h,2)/=0) then                   ! If ODD
          do g=1,CoreIB%nCores
            if ((CoreIB%StartSnp(g)<UptoSnp).and.(CoreIB%EndSnp(g)>UptoSnp)) then
              UpToCoreB=g
              exit
            endif
          enddo
          if (e==1) then
            StartSnp=CoreIB%StartSnp(UpToCoreB)
            EndSnp=CoreIB%EndSnp(UpToCoreB)

            StPt=StartSnp
            EndPt=UpToSnp
            FillInSt=StartSnp
            FillInEnd=EndSnp
          else
            StartSnp=CoreIB%EndSnp(UpToCoreB)
            EndSnp=CoreIB%StartSnp(UpToCoreB)
            StPt=UpToSnp
            EndPt=StartSnp
            FillInSt=EndSnp
            FillInEnd=StartSnp
          endif

          CompLength=abs(UpToSnp-StartSnp)+1

          do i=1,nAnisP
            if (ped%pedigree(i)%founder .and.(PosHD(i)/=0).and.(AnimRecomb(i,RL)==0)) then
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

              ! If one haplotype is the same as the paternal, impute
              if ((CompLength==C1).and.(CompLength/=C3)) then
                ImputePhase(i,FillInSt:FillInEnd,1)=PhaseHD(PosHD(i),FillInSt:FillInEnd,1,2)
                Recmb=0
              endif

              ! If one haplotype is the same as the paternal, impute
              if ((CompLength/=C1).and.(CompLength==C3)) then
                ImputePhase(i,FillInSt:FillInEnd,1)=PhaseHD(PosHD(i),FillInSt:FillInEnd,2,2)
                Recmb=0
              endif

              ! If one haplotype is the same as the maternal, impute
              if ((CompLength==C2).and.(CompLength/=C4)) then
                ImputePhase(i,FillInSt:FillInEnd,2)=PhaseHD(PosHD(i),FillInSt:FillInEnd,1,2)
                Recmb=0
              endif

              ! If one haplotype is the same as the maternal, impute
              if ((CompLength/=C2).and.(CompLength==C4)) then
                ImputePhase(i,FillInSt:FillInEnd,2)=PhaseHD(PosHD(i),FillInSt:FillInEnd,2,2)
                Recmb=0
              endif

              ! There is bridge in the phase in the RL direction. Nothing more will be done
              if (Recmb==1) then
                AnimRecomb(i,RL)=1
              end if
            endif
          enddo
          if (RL == 1) then
            UpToSnp=CoreIB%StartSnp(UpToCoreB)
          else
            UpToSnp=CoreIB%EndSnp(UpToCoreB)
          end if
        else                                    ! if EVEN
          do g=1,CoreIA%nCores
            if ((CoreIA%StartSnp(g)<UptoSnp).and.(CoreIA%EndSnp(g)>UptoSnp)) then
              UpToCoreA=g
              exit
            endif
          enddo
          if (e==1) then
            StartSnp=CoreIA%StartSnp(UpToCoreA)
            EndSnp=CoreIA%EndSnp(UpToCoreA)

            StPt=StartSnp
            EndPt=UpToSnp
            FillInSt=StartSnp
            FillInEnd=EndSnp
          else
            StartSnp=CoreIA%StartSnp(UpToCoreA)
            EndSnp=CoreIA%EndSnp(UpToCoreA)

            StPt=UpToSnp
            EndPt=StartSnp
            FillInSt=EndSnp
            FillInEnd=StartSnp
          endif
          CompLength=abs(UpToSnp-StartSnp)+1
          do i=1,nAnisP
            if ( ped%pedigree(i)%founder .and.(PosHD(i)/=0).and.(AnimRecomb(i,RL)==0)) then
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
              if (Recmb==1) then
                AnimRecomb(i,RL)=1
              end if
            endif
          enddo

          if (RL == 1) then
            UpToSnp=CoreIA%StartSnp(UpToCoreA)
          else
            UpToSnp=CoreIA%EndSnp(UpToCoreA)
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
    use alphaimputeinmod
    implicit none

    integer :: i,j,dum

    ImputeGenos=9
    ImputePhase=9
    inputParams => defaultInput
    ! Get information from RecodedGeneProbInput.txt which has been created in Makefiles subroutine
    ! WARNING: Why don't read information from Geno(:,:) that has been used to feed RecodedGeneProbInput.txt instead??
    !          Read from file is always slower!
    open (unit=43,file='.' // DASH // 'InputFiles' // DASH // 'RecodedGeneProbInput.txt',status='old')
    do i=1,nAnisP
      read (43,*) dum,dum,dum,ImputeGenos(i,:)
      do j=1,inputParams%nsnp
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

  SUBROUTINE EnsureHetGametic
    ! Impute phase to Y chromosome from X chromosome for heterogametic individuals
    use Global
    use alphaimputeinmod
    implicit none

    integer :: i,j

    inputParams => defaultInput

    do j=1,inputParams%nsnp
      do i=1,nAnisP
        if (ped%pedigree(i)%gender==HetGameticStatus) then
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
    use alphaimputeinmod

    implicit none

    integer :: i,j

    inputParams => defaultInput
    do j=1,inputParams%nsnp
      do i=1,nAnisP
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

    do j=1,inputParams%nsnp
      do i=1,nAnisP
        if (ImputeGenos(i,j)/=9) then
            if ((ImputePhase(i,j,1)/=9).and.(ImputePhase(i,j,2)==9)) then
              ImputePhase(i,j,2)=ImputeGenos(i,j)-ImputePhase(i,j,1)
            end if
            if ((ImputePhase(i,j,2)/=9).and.(ImputePhase(i,j,1)==9)) then
              ImputePhase(i,j,1)=ImputeGenos(i,j)-ImputePhase(i,j,2)
            end if
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
    use alphaimputeinmod
    implicit none

    integer :: e,i,j,ParId

    inputParams => defaultInput
    do i=1,nAnisP
      if (inputParams%sexopt==0 .or. (inputParams%sexopt==1 .and. ped%pedigree(i)%gender/=HetGameticStatus) ) then     ! If individual is homogametic
        do e=1,2
          if (ped%pedigree(i)%isDummyBasedOnIndex(e)) cycle
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
        ParId= ped%pedigree(i)%getSireDamNewIDByIndex(HomGameticStatus+1) !the homogametic parent
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
    use alphaimputeinmod
    use individualModule
    implicit none

    integer :: i,j,k,l,Count1,Count0
    type(individual), pointer :: tmpChild
    inputParams => defaultInput
    do i=1,nAnisP
      do j=1,2
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
                    if ((ped%pedigree(i)%gender ==HetGameticStatus).and.(tmpChild%gender==HetGameticStatus)) cycle 
                    if (ImputePhase(tmpChild%id,k,j)==0) Count0=Count0+1
                    if (ImputePhase(tmpChild%id,k,j)==1) Count1=Count1+1
            
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

                do l=1,ped%pedigree(i)%nOffs

                  tmpChild => ped%pedigree(i)%offsprings(l)%p

                  ! This is the only difference with the inputParams%sexopt=0 code below. Duplicating
                  ! the code can be avoided by including the IF statement here instead than
                  ! outside the SNPs loop.
                  if ((ped%pedigree(i)%gender ==HetGameticStatus).and.(tmpChild%gender==HetGameticStatus)) cycle 
                  if (ImputePhase(tmpChild%id,k,j)==0) Count0=Count0+1
                  if (ImputePhase(tmpChild%id,k,j)==1) Count1=Count1+1

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
            do k=1,inputParams%nsnp
              if ((ImputePhase(i,k,1)/=9).and.(ImputePhase(i,k,2)==9)) then               !Pat gamete missing fill if offspring suggest heterozygous
                Count1=0
                Count0=0
                if (ImputePhase(i,k,1)==1) Count1=1
                if (ImputePhase(i,k,1)==0) Count0=1

                do l=1,ped%pedigree(i)%nOffs
                    tmpChild => ped%pedigree(i)%offsprings(l)%p
                    if (ImputePhase(tmpChild%id,k,j)==0) Count0=Count0+1
                    if (ImputePhase(tmpChild%id,k,j)==1) Count1=Count1+1
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
                do l=1,ped%pedigree(i)%nOffs
                  tmpChild => ped%pedigree(i)%offsprings(l)%p
                  if (ImputePhase(Tmpchild%id,k,j)==0) Count0=Count0+1
                  if (ImputePhase(Tmpchild%id,k,j)==1) Count1=Count1+1
                
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
    
    use alphaimputeinmod

    implicit none

    integer :: h,i,j,k,m,dum,StSnp,EnSnp
    real :: PatAlleleProb(inputParams%nsnp,2),MatAlleleProb(inputParams%nsnp,2),HetProb(inputParams%nsnp),GeneProbWork(inputParams%nsnp,4)
    integer :: Informativeness(inputParams%nsnp,6),TmpInfor(inputParams%nsnp,6),GrandPar
    character(len=300) :: filout

    inputParams => defaultInput
    if (inputParams%BypassGeneProb==0) then
    ! Get information from GeneProb
      do h=1,inputParams%nProcessors
#ifdef OS_UNIX
        write (filout,'("./GeneProb/GeneProb"i0,"/GeneProbs.txt")')h            !here
#else
        write (filout,'(".\GeneProb\GeneProb"i0,"\GeneProbs.txt")')h            !here
#endif
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

    open(unit=102,file="." // DASH // "Miscellaneous" // "IndividualSnpInformativeness.txt", status="unknown")
    open(unit=103,file="." // DASH // "Miscellaneous" // "IndividualMendelianInformativeness.txt", status="unknown")

    allocate(GlobalTmpCountInf(nAnisP,8))
    allocate(MSTermInfo(nAnisP,2))

    MSTermInfo=0
    do i=1,nAnisP
      if (ped%pedigree(i)%Genotyped) MSTermInfo(i,:)=1
      TmpInfor(:,:)=-99
      GlobalTmpCountInf(i,:)=0
      Informativeness(:,:)=9 ! What the hell is this variable for??
      j=0
  

    ! if (ped%pedigree(i)%hasDummyParent()) cycle
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


            ! if (ped%pedigree(i)%hasDummyParentsOrGranparents()) cycle
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
    use alphaimputeinmod

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

END MODULE Imputation