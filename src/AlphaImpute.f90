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

!#############################################################################################################################################################################################################################

module AlphaImputeModule
    implicit none    
contains




    subroutine PhasingManagementNew(results)
        use AlphaImputeInMod
        use omp_lib
        use AlphaPhaseParametersDefinition
        use AlphaPhaseFunctions
        use AlphaPhaseResultsDefinition
        use Global, only : ped
        implicit none

        type(AlphaImputeInput), pointer :: inputParams
        type(AlphaPhaseParameters) :: params
        type(AlphaPhaseResultsContainer), intent(out) :: results
        integer :: nCoreLengths,i

        inputParams=> defaultInput
        nCoreLengths = size(inputParams%CoreAndTailLengths)
        results%nResults = nCoreLengths
                ! todo we need to handle shifted and non shifted outputs!
        allocate(results%results(nCoreLengths))
        params = newParameters()

    
        call omp_set_nested(.true.)




        print *, "Running AlphaPhase"
       !$OMP parallel DO schedule(dynamic) &
       !$OMP FIRSTPRIVATE(params)
        do i= 1, nCoreLengths

            params%CoreAndTailLength = inputParams%CoreAndTailLengths(i)
            params%itterateNumber = inputParams%PhaseNIterations
            params%jump = inputParams%CoreAndTailLengths(i)
            params%numsurrdisagree = 10
            params%useSurrsN = 10
            params%PercGenoHaploDisagree = inputParams%GenotypeErrorPhase
            results%results(i) = phaseAndCreateLibraries(ped, params, quiet=.true.)
        enddo 
        
        !$omp end parallel do
        print *, "Finished Running AlphaPhase"


    end subroutine PhasingManagementNew
    !#############################################################################################################################################################################################################################

    subroutine IterateGeneProbsNew(GenosProbs)

        use iso_fortran_env
        use global, only : Maf,ImputePhase,ImputeGenos, ped,ProbImputePhase, ProbImputeGenos
        use AlphaImputeInMod
        use GeneProbModule , only : runGeneProbAlphaImpute
        use imputation, only : PhaseComplement

        integer :: i,j,k

        integer :: startSnp, endsnp, nSnp
        type(AlphaImputeInput), pointer :: inputParams
        real(kind=real64), intent(out), allocatable :: GenosProbs(:,:,:)
        inputParams=> defaultInput

        if (inputParams%outopt==0) nSnp=inputParams%nsnp
        if (inputParams%outopt==1) nSnp=inputParams%nSnpRaw


        ! This will be MPI jobs in future, as right now parallelization is run using openMP inside Geneprob itself
        ! TODO add tihs in
        ! do i=1,inputParams%nprocessors
            ! startsnp = GpIndex(i,1)
            ! endsnp = GpIndex(i,2)

            startSnp = 1
            endsnp = inputParams%nSnp
            ! No writes should be made to pedigree, so this can be shared. GenosProbs and MAF need to be combined
            call runGeneProbAlphaImpute(StartSnp, endsnp, ped, GenosProbs, MAF)

        ! enddo


        call IterateGeneProbPhase
        call IterateParentHomoFill
        call PhaseComplement
        call IterateMakeGenotype

        do i=1,ped%pedigreeSize-ped%nDummys
            do j=1,nSnp
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

    end subroutine IterateGeneProbsNew


    !#############################################################################################################################################################################################################################



    !#############################################################################################################################################################################################################################

    subroutine IterateInsteadOfGeneProbs
        use Global
        use Imputation
        use alphaimputeinmod
        implicit none

        integer :: e,i,j,k,Counter,ParId
        real,allocatable,dimension(:) :: TempAlleleFreq

        inputParams => defaultInput
        if (inputParams%SexOpt==1) then
            if (inputParams%outopt==0) nSnpIterate=inputParams%nsnp
            if (inputParams%outopt==1) nSnpIterate=inputParams%nSnpRaw

            allocate(ProbImputeGenos(0:ped%pedigreeSize-ped%nDummys,nSnpIterate))
            allocate(ProbImputePhase(0:ped%pedigreeSize-ped%nDummys,nSnpIterate,2))
            allocate(TempAlleleFreq(nSnpIterate))

            TempAlleleFreq=0.0
            do j=1,nSnpIterate
                Counter=0
                do i=1,ped%pedigreeSize-ped%nDummys
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
            ProbImputeGenos(1:ped%pedigreeSize-ped%nDummys,:)=-9.0
            ProbImputePhase(1:ped%pedigreeSize-ped%nDummys,:,:)=-9.0

            call IterateParentHomoFill
            call PhaseComplement
            call IterateMakeGenotype

            do i=1,ped%pedigreeSize-ped%nDummys
                do j=1,nSnpIterate
                    do e=1,2
                        if (ImputePhase(i,j,e)/=9) ProbImputePhase(i,j,e)=float(ImputePhase(i,j,e))
                    enddo
                enddo
            enddo

            do i=1,ped%pedigreeSize-ped%nDummys
                do e=1,2
                    parId = ped%pedigree(i)%getSireDamNewIDByIndex(e+1)
                    if (ParId==0) then
                        do j=1,nSnpIterate
                            if (ImputePhase(i,j,e)==9) ProbImputePhase(i,j,e)=TempAlleleFreq(j)
                        enddo
                    endif
                    if (ped%pedigree(i)%gender==inputParams%HomGameticStatus) then
                        do j=1,nSnpIterate
                            if (ImputePhase(i,j,e)==9) then
                                ProbImputePhase(i,j,e)=(sum(ProbImputePhase(ParId,j,:))/2)
                            endif
                        enddo
                    endif
                enddo
                if (ped%pedigree(i)%gender==inputParams%hetGameticStatus) then
                    ParId=ped%pedigree(i)%getSireDamNewIDByIndex(inputParams%hetGameticStatus+1)
                    do j=1,nSnpIterate
                        if (ImputePhase(i,j,1)==9) then
                            ProbImputePhase(i,j,:)=(sum(ProbImputePhase(ParId,j,:))/2)
                        endif
                    enddo
                endif
            enddo

            do i=1,ped%pedigreeSize-ped%nDummys
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

            if (inputParams%SexOpt==1) then
                do i=1,ped%pedigreeSize-ped%nDummys
                    if (ped%pedigree(i)%gender==inputParams%hetGameticStatus) then
                        do j=1,nSnpIterate
                            if ((ImputePhase(i,j,1)==9).and.(ImputePhase(i,j,2)/=9)) ImputePhase(i,j,1)=ImputePhase(i,j,2)
                            if ((ImputePhase(i,j,2)==9).and.(ImputePhase(i,j,1)/=9)) ImputePhase(i,j,2)=ImputePhase(i,j,1)
                            if ((ProbImputePhase(i,j,1)==-9.0).and.(ProbImputePhase(i,j,2)/=-9.0)) ProbImputePhase(i,j,1)=ProbImputePhase(i,j,2)
                            if ((ProbImputePhase(i,j,2)==-9.0).and.(ProbImputePhase(i,j,1)/=-9.0)) ProbImputePhase(i,j,2)=ProbImputePhase(i,j,1)
                        enddo
                    endif
                enddo
            endif

            do i=1,ped%pedigreeSize-ped%nDummys
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

            if (inputParams%outopt==0) nSnpIterate=inputParams%nsnp
            if (inputParams%outopt==1) nSnpIterate=inputParams%nSnpRaw

            allocate(ProbImputeGenos(0:ped%pedigreeSize-ped%nDummys,nSnpIterate))
            allocate(ProbImputePhase(0:ped%pedigreeSize-ped%nDummys,nSnpIterate,2))
            allocate(TempAlleleFreq(nSnpIterate))

            TempAlleleFreq=0.0
            do j=1,nSnpIterate
                Counter=0
                do i=1,ped%pedigreeSize-ped%nDummys
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
            ProbImputeGenos(1:ped%pedigreeSize-ped%nDummys,:)=-9.0
            ProbImputePhase(1:ped%pedigreeSize-ped%nDummys,:,:)=-9.0

            call IterateParentHomoFill
            call PhaseComplement
            call IterateMakeGenotype

            do i=1,ped%pedigreeSize-ped%nDummys
                do j=1,nSnpIterate
                    do e=1,2
                        if (ImputePhase(i,j,e)/=9) ProbImputePhase(i,j,e)=float(ImputePhase(i,j,e))
                    enddo
                enddo
            enddo

            do i=1,ped%pedigreeSize-ped%nDummys
                do e=1,2
                    parID=ped%pedigree(i)%getSireDamNewIDByIndex(e+1)
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

            do i=1,ped%pedigreeSize-ped%nDummys
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

            do i=1,ped%pedigreeSize-ped%nDummys
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

        use output
        use GlobalVariablesHmmMaCH
        use AlphaImputeInMod

        use output
        implicit none

        character(len=7) :: cm !use for formatting output - allows for up to 1 million SNPs
        integer :: i,j,l
        integer,allocatable,dimension(:):: WorkTmp
        double precision :: ImputationQuality(ped%pedigreeSize-ped%nDummys,6)
        
        type(AlphaImputeInput), pointer :: inputParams
        character(len=300) :: TmpId
        integer :: n0, n1, n2



        inputParams => defaultInput

#ifdef DEBUG
        write(0,*) 'DEBUG: WriteOutResults'
#endif

        allocate(WorkTmp(inputParams%nSnpRaw))



        write(cm,'(I7)') inputParams%nSnpRaw !for formatting
        cm = adjustl(cm)

        open (unit=33,file="." // DASH// "Results" // DASH // "ImputePhase.txt",status="unknown")
        open (unit=34,file="." // DASH// "Results" // DASH // "ImputeGenotypes.txt",status="unknown")
        open (unit=40,file="." // DASH// "Results" // DASH // "ImputePhaseProbabilities.txt",status="unknown")
        open (unit=41,file="." // DASH// "Results" // DASH // "ImputeGenotypeProbabilities.txt",status="unknown")
        open (unit=50,file="." // DASH// "Results" // DASH // "ImputationQualityIndividual.txt",status="unknown")
        open (unit=51,file="." // DASH// "Results" // DASH // "ImputationQualitySnp.txt",status="unknown")
        open (unit=52,file="." // DASH// "Results" // DASH // "WellPhasedIndividuals.txt",status="unknown")

        open (unit=53,file="." // DASH// "Results" // DASH // "ImputePhaseHMM.txt",status="unknown")
        open (unit=54,file="." // DASH// "Results" // DASH // "ImputeGenotypesHMM.txt",status="unknown")


#ifdef DEBUG
        write(0,*) 'DEBUG: output=0 [WriteOutResults]'
#endif

        if (inputParams%outopt==0) then

            if (inputParams%SexOpt==0) then

                call CheckImputationInconsistencies(ImputeGenos, ImputePhase, ped%pedigreeSize-ped%nDummys, inputParams%nsnp)

                open (unit=39, file="IterateGeneProb" // DASH // "IterateGeneProbInput.txt")
                do i=1,ped%pedigreeSize-ped%nDummys
                    write (39,'(i16,1x,i16,1x,i16,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(i)%getIntegerVectorOfRecodedIdsNoDummy(),ImputeGenos(i,:)
                enddo
                call flush(39)
                ! close (39)

                if (inputParams%hmmoption == RUN_HMM_NO) then
                    if (inputParams%bypassgeneprob==0) then
                        call IterateGeneProbsNew(GenosProbs)
                    else
                        call IterateInsteadOfGeneProbs
                    endif
                end if
            else

                call CheckImputationInconsistencies(ImputeGenos, ImputePhase, ped%pedigreeSize-ped%nDummys, inputParams%nsnp)

                call IterateInsteadOfGeneProbs
            endif

            call CheckImputationInconsistencies(ImputeGenos, ImputePhase, ped%pedigreeSize-ped%nDummys, inputParams%nsnp)

            do i=1, ped%pedigreeSize-ped%nDummys
                if (ped%pedigree(i)%isDummy) then
                    exit
                endif
                write (33,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(i)%originalID,ImputePhase(i,:,1)
                write (33,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(i)%originalID,ImputePhase(i,:,2)
                write (34,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(i)%originalID,ImputeGenos(i,:)
            enddo


            if (inputParams%hmmoption /= RUN_HMM_NO) then

        ! #ifdef DEBUG
        !         write(0,*) 'DEBUG: Write HMM results [WriteOutResults]'
        ! #endif

                if (allocated(ImputeGenos)) then
                    deallocate(ImputeGenos)
                end if
                if (allocated(ImputePhase)) then
                    deallocate(ImputePhase)
                end if
                if (allocated(ProbImputeGenos)) then
                    deallocate(ProbImputeGenos)
                end if
                if (allocated(ProbImputePhase)) then
                    deallocate(ProbImputePhase)
                end if
                if (allocated(maf)) then
                    deallocate(maf)
                end if
                allocate(ImputeGenos(0:ped%nGenotyped,inputParams%nSnp))
                allocate(ImputePhase(0:ped%nGenotyped,inputParams%nSnp,2))
                allocate(ProbImputeGenos(0:ped%nGenotyped,inputParams%nSnp))
                allocate(ProbImputePhase(0:ped%nGenotyped,inputParams%nSnp,2))
                allocate(Maf(inputParams%nSnp))

                ImputeGenos = 9
                ImputePhase = 9
                ProbImputeGenos(1:ped%pedigreeSize-ped%nDummys,:) = 9.0
                ProbImputePhase(1:ped%pedigreeSize-ped%nDummys,:,:) = 9.0

                ! Feed Impute and Phase probabilites
                l=0
                do j=1,inputParams%nsnp
                        l=l+1
                        do i=1,ped%nGenotyped
                            ! if (GlobalHmmID(i) > ped%pedigreeSize-ped%nDummys ) then
                            !     GlobalHmmID(i) = 0
                            !     ! TODO this means animal is a dummy - need to deal with this
                            ! endif
                            ProbImputeGenos(i,j)   = ProbImputeGenosHmm(i,j)
                            ProbImputePhase(i,j,1) = ProbImputePhaseHmm(i,j,1)
                            ProbImputePhase(i,j,2) = ProbImputePhaseHmm(i,j,2)
                        enddo
                enddo

#ifdef DEBUG
                write(0,*) 'DEBUG: Impute alleles and genotypes based on HMM genotypes probabilities [WriteOutResults]'
#endif
                ! Impute the most likely genotypes. (Most frequent genotype)
                do i=1,ped%nGenotyped
                    do j=1,inputParams%nSnp
                        n2 = GenosCounts(i,j,2)                           ! Homozygous: 2 case
                        n1 = GenosCounts(i,j,1)                           ! Heterozygous
                        n0 = (inputParams%nRoundsHmm-inputParams%HmmBurnInRound) - n1 - n2        ! Homozygous: 0 case
                        if ((n0>n1).and.(n0>n2)) then
                            ImputeGenos(i,j)   = 0
                            ImputePhase(i,j,:) = 0
                        elseif (n1>n2) then
                            ImputeGenos(GlobalHmmID(i),j) = 1
                            if (ProbImputePhaseHmm(i,j,1) > ProbImputePhaseHmm(i,j,2) ) then
                                ImputePhase(i,j,1) = 1
                                ImputePhase(i,j,2) = 0
                            else
                                ImputePhase(i,j,1) = 0
                                ImputePhase(i,j,2) = 1
                            endif
                        else
                            ImputeGenos(i,j)   = 2
                            ImputePhase(i,j,:) = 1
                        endif
                    enddo
                enddo

                ! call CheckImputationInconsistencies(ImputeGenos, ImputePhase, ped%nGenotyped, inputParams%nSnp)

                BLOCK
                    integer :: hmmID
                    do i=1,ped%nGenotyped
                        ! TODO: Remove this variable with the isssue is really fixed
                        ! hmmID = GlobalHmmID(i)
                        hmmID = ped%genotypeMap(i)
                        write (53,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(hmmID)%originalID,ImputePhase(i,:,1)
                        write (53,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(hmmID)%originalID,ImputePhase(i,:,2)
                        write (54,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(hmmID)%originalID,ImputeGenos(i,:)
                        write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(hmmID)%originalID,ProbImputePhase(i,:,1)
                        write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(hmmID)%originalID,ProbImputePhase(i,:,2)
                        write (41,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(hmmID)%originalID,ProbImputeGenos(i,:)
                    enddo
                END BLOCK
            else
                do i=1, ped%pedigreeSize-ped%nDummys
                    if (ped%pedigree(i)%isDummy) then
                        exit
                    endif
                    write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(i)%originalID,ProbImputePhase(i,:,1)
                    write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(i)%originalID,ProbImputePhase(i,:,2)
                    write (41,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(i)%originalID,ProbImputeGenos(i,:)
                enddo
            endif

            if (allocated(Maf)) Then
                deallocate(Maf)
            endif
            allocate(Maf(inputParams%nSnp))

            if (inputParams%SexOpt==1) then
                allocate(Maf(inputParams%nsnp))
                do j=1,inputParams%nsnp
                    Maf(j)=sum(ProbImputeGenos(:,j))/(2*ped%pedigreeSize-ped%nDummys)
                enddo
                open(unit=111,file="." // DASH // "Miscellaneous" // DASH // "MinorAlleleFrequency.txt", status="unknown")

                do j=1,inputParams%nSnpRaw
                    write (111,*) j,Maf(j)
                enddo
                close(111)
            endif
            ImputationQuality(:,1)=sum(2*Maf(:))/inputParams%nsnp

            ImputationQuality(:,2)=0.0
            do i=1, ped%pedigreeSize-ped%nDummys
                if (ped%pedigree(i)%isDummy) then
                    exit
                endif
                do j=1,inputParams%nsnp
                    !ImputationQuality(i,2)=ImputationQuality(i,2)+abs(ProbImputeGenos(i,j)-(2*Maf(j)))
                    ImputationQuality(i,2)=ImputationQuality(i,2)+abs(ProbImputeGenos(i,j)-((Maf(j)**4)+(4*(Maf(j)**2)*((1.0-Maf(j))**2))+(1.0-Maf(j)**4)))
                enddo
                ImputationQuality(i,2)=ImputationQuality(i,2)/inputParams%nsnp
                ImputationQuality(i,3)=float(inputParams%nsnp-count(ImputePhase(i,:,1)==9))/inputParams%nsnp
                ImputationQuality(i,4)=float(inputParams%nsnp-count(ImputePhase(i,:,2)==9))/inputParams%nsnp
                ImputationQuality(i,5)=(ImputationQuality(i,3)+ImputationQuality(i,4))/2
                ImputationQuality(i,6)=float(inputParams%nsnp-count(ImputeGenos(i,:)==9))/inputParams%nsnp
                write (50,'(a20,20000f7.2)') ped%pedigree(i)%originalID,ImputationQuality(i,:)
            enddo

            do j=1,inputParams%nsnp
                write (51,'(i10,20000f7.2)') j,float(((ped%pedigreeSize-(ped%nDummys+1))+1)-count(ImputeGenos(ped%nDummys+1:ped%pedigreeSize,j)==9))/((ped%pedigreeSize-(ped%nDummys+1))+1)
            enddo

            call CheckImputationInconsistencies(ImputeGenos, ImputePhase, ped%pedigreeSize-ped%nDummys, inputParams%nsnp)
            inputParams%WellPhasedThresh=inputParams%WellPhasedThresh/100
            do i=1, ped%pedigreeSize-ped%nDummys
                if (ped%pedigree(i)%isDummy) then
                    exit
                endif
                if (ImputationQuality(i,5)>=inputParams%WellPhasedThresh) then
                    write (52,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(i)%originalID,ImputePhase(i,:,1)
                    write (52,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(i)%originalID,ImputePhase(i,:,2)
                endif
            enddo

        else

#ifdef DEBUG
            write(0,*) 'DEBUG: Unphase wrong alleles [WriteOutResults]'
#endif

            call CheckImputationInconsistencies(ImputeGenos, ImputePhase, ped%pedigreeSize-ped%nDummys, inputParams%nsnp)


! TODO this might be erronous
            if (inputParams%outopt == 1) then
                allocate(TmpGenos(0:ped%pedigreeSize-ped%nDummys,inputParams%nSnpRaw))
                allocate(TmpPhase(0:ped%pedigreeSize-ped%nDummys,inputParams%nSnpRaw,2))
            else
                allocate(TmpGenos(0:ped%pedigreeSize-ped%nDummys,inputParams%nSnp))
                allocate(TmpPhase(0:ped%pedigreeSize-ped%nDummys,inputParams%nSnp,2))
            endif
            TmpGenos=9
            TmpPhase=9

            if (inputParams%hmmoption==RUN_HMM_NGS) SnpIncluded=1

            l=0
            do j=1,inputParams%nSnpRaw
                if (SnpIncluded(j)==1) then
                    l=l+1
                    TmpGenos(:,j)=ImputeGenos(:,l)
                    TmpPhase(:,j,1)=ImputePhase(:,l,1)
                    TmpPhase(:,j,2)=ImputePhase(:,l,2)
                endif
            enddo

            block
                integer :: tmpIDInt
                if (inputParams%inteditstat == 1) then
                    open (unit=42,file=trim(inputParams%GenotypeFile),status='old')
                    do i=1,ped%nGenotyped
                        read (42,*) TmpId,WorkTmp(:)

                        tmpIDInt = ped%dictionary%getValue(trim(tmpID))

                        if (tmpIDInt /= DICT_NULL) then
                            do j=1,inputParams%nSnpRaw
                                if (SnpIncluded(j)==0) then
                                    if (WorkTmp(j)==1) TmpGenos(tmpIDInt,j)=1
                                    if (WorkTmp(j)==0) then
                                        TmpPhase(tmpIDInt,j,:)=0
                                        TmpGenos(tmpIDInt,j)=0
                                    endif
                                    if (WorkTmp(j)==2) then
                                        TmpPhase(tmpIDInt,j,:)=1
                                        TmpGenos(tmpIDInt,j)=2
                                    endif
                                endif
                            enddo
                            exit
                        endif

                    enddo

                    close(42)
                end if
            endblock
            call CheckImputationInconsistencies(TmpGenos, TmpPhase, ped%pedigreeSize-ped%nDummys, inputParams%nsnp)
            do i=1, ped%pedigreeSize-ped%nDummys
                if (ped%pedigree(i)%isDummy) then
                    exit
                endif
                write (33,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(i)%originalID,TmpPhase(i,:,1)
                write (33,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(i)%originalID,TmpPhase(i,:,2)
                write (34,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(i)%originalID,TmpGenos(i,:)
            enddo
            if (inputParams%SexOpt==0 .and. inputParams%hmmoption/=RUN_HMM_NGS) then
                open (unit=39, file="IterateGeneProb" // DASH // "IterateGeneProbInput.txt")

                do i=1,ped%pedigreeSize-ped%nDummys
                    write (39,'(i16,1x,i16,1x,i16,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(i)%getIntegerVectorOfRecodedIdsNoDummy(),TmpGenos(i,:)
                enddo
                call flush(39)
                close(39)
            endif

            !REMOVE THIS WHEN HMM IS FINALISED
            if (inputParams%hmmoption==RUN_HMM_NO) then
                deallocate(ImputePhase)
                deallocate(ImputeGenos)
                allocate(ImputeGenos(0:ped%pedigreeSize-ped%nDummys,inputParams%nSnpRaw))
                allocate(ImputePhase(0:ped%pedigreeSize-ped%nDummys,inputParams%nSnpRaw,2))
                ImputeGenos=TmpGenos
                ImputePhase=TmpPhase
                if (inputParams%SexOpt==0) then
                    if (inputParams%bypassgeneprob==0) then
                        call IterateGeneProbsNew(GenosProbs)
                    else
                        call IterateInsteadOfGeneProbs
                    endif
                endif
                if (inputParams%SexOpt==1) call IterateInsteadOfGeneProbs
            endif
            !REMOVE THIS

            call CheckImputationInconsistencies(ImputeGenos, ImputePhase, ped%pedigreeSize-ped%nDummys, inputParams%nsnp)

            !if (inputParams%hmmoption==RUN_HMM_ONLY.or.inputParams%hmmoption==RUN_HMM_PREPHASE) then
            if (inputParams%hmmoption/=RUN_HMM_NO) then

#ifdef DEBUG
                write(0,*) 'DEBUG: Write HMM results [WriteOutResults]'
#endif
                ! nSnpIterate=inputParams%nsnp
                if (allocated(ImputeGenos)) then
                    deallocate(ImputeGenos)
                end if
                allocate(ImputeGenos(0:ped%nGenotyped,inputParams%nSnpRaw))
                if (allocated(ImputePhase)) then
                    deallocate(ImputePhase)
                end if
                allocate(ImputePhase(0:ped%nGenotyped,inputParams%nSnpRaw,2))
                if (allocated(ProbImputeGenos)) then
                    deallocate(ProbImputeGenos)
                end if
                allocate(ProbImputeGenos(0:ped%pedigreeSize-ped%nDummys,inputParams%nSnpRaw))
                if (allocated(ProbImputePhase)) then
                    deallocate(ProbImputePhase)
                end if
                allocate(ProbImputePhase(0:ped%pedigreeSize-ped%nDummys,inputParams%nSnpRaw,2))
                if (allocated(Maf)) then
                    deallocate(Maf)
                end if
                allocate(Maf(inputParams%nSnpRaw))
                ProbImputeGenos(1:ped%pedigreeSize-ped%nDummys,:)= 9.0
                ProbImputePhase(1:ped%pedigreeSize-ped%nDummys,:,:)= 9.0
                ImputeGenos = 9
                ImputePhase = 9

                l=0
                do j=1,inputParams%nsnp
                        l=l+1
                        do i=1,ped%nGenotyped

                            ProbImputeGenos(i,j)   = ProbImputeGenosHmm(i,j)
                            ProbImputePhase(i,j,1) = ProbImputePhaseHmm(i,j,1)
                            ProbImputePhase(i,j,2) = ProbImputePhaseHmm(i,j,2)
                        enddo
                enddo

#ifdef DEBUG
                write(0,*) 'DEBUG: Impute alleles and genotypes based on HMM genotypes probabilities [WriteOutResults]'
#endif
                ! Impute the most likely genotypes. (Most frequent genotype)
                do i=1,ped%nGenotyped
                ! TODO check if this should be nsnp or nsnpraw
                    do j=1,inputParams%nSnp
                        n2 = GenosCounts(i,j,2)                           ! Homozygous: 2 case
                        n1 = GenosCounts(i,j,1)                           ! Heterozygous
                        n0 = (inputParams%nRoundsHmm-inputParams%HmmBurnInRound) - n1 - n2        ! Homozygous: 0 case
                        if ((n0>n1).and.(n0>n2)) then
                            ImputeGenos(i,j)   = 0
                            ImputePhase(i,j,:) = 0
                        elseif (n1>n2) then
                            ImputeGenos(i,j) = 1
                            if (ProbImputePhaseHmm(i,j,1) > ProbImputePhaseHmm(i,j,2) ) then
                                ImputePhase(i,j,1) = 1
                                ImputePhase(i,j,2) = 0
                            else
                                ImputePhase(i,j,1) = 0
                                ImputePhase(i,j,2) = 1
                            endif
                        else
                            ImputeGenos(i,j)   = 2
                            ImputePhase(i,j,:) = 1
                        endif
                    enddo
                enddo

#ifdef DEBUG
                write(0,*) 'DEBUG: Write phase, genotypes and probabilities into files [WriteOutResults]'
#endif
                ! call CheckImputationInconsistencies(ImputeGenos, ImputePhase, nAnisP, inputParams%nsnp)
                !  TODO remove
                BLOCK
                    integer :: hmmID
                    do i=1, ped%nGenotyped
                        ! TODO: Remove this variable when this issue is really fixed
                        ! hmmID = GlobalHmmID(i)
                        hmmID =i
                        write (53,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(hmmID)%originalID,ImputePhase(i,:,1)
                        write (53,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(hmmID)%originalID,ImputePhase(i,:,2)
                        write (54,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(hmmID)%originalID,ImputeGenos(i,:)

                        write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(hmmID)%originalID,ProbImputePhase(i,:,1)
                        write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(hmmID)%originalID,ProbImputePhase(i,:,2)
                        write (41,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(hmmID)%originalID,ProbImputeGenos(i,:)
                    enddo
                END BLOCK
            else
                do i=1, ped%pedigreeSize-ped%nDummys
                    if (ped%pedigree(i)%isDummy) then
                        exit
                    endif
                    write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(i)%originalID,ProbImputePhase(i,:,1)
                    write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(i)%originalID,ProbImputePhase(i,:,2)
                    write (41,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(i)%originalID,ProbImputeGenos(i,:)
                enddo
            endif

            if (inputParams%hmmoption/=RUN_HMM_NO) then
                ! call WriteProbabilities("./Results/GenotypeProbabilities.txt", GlobalHmmID, ID, ped%nGenotyped, inputParams%nsnp)
                call WriteProbabilities("./Results/GenotypeProbabilities.txt", GlobalHmmID, ped%nGenotyped, inputParams%nsnp)
            else
                if (inputParams%bypassgeneprob==0) then
                    allocate(GenosProbs(ped%pedigreeSize-ped%nDummys,nSnpIterate,2))
                    ! TODOGENEPROB geneprob should not have been bypassed here so check that it indeed is not
                    call WriteProbabilities("./Results/GenotypeProbabilities.txt", GenosProbs, ped,ped%pedigreeSize-ped%nDummys, inputParams%nsnp)
                endif
            endif

            if ((inputParams%SexOpt==1).or.(inputParams%bypassgeneprob==1)) then

#ifdef DEBUG
                write(0,*) 'DEBUG: Bypass genotype probabilities [WriteOutResults]'
#endif

                if (inputParams%hmmoption/=RUN_HMM_NO) Then
                    deallocate(Maf)
                endif
                allocate(Maf(inputParams%nSnpRaw))
                do j=1,inputParams%nSnpRaw
                    Maf(j)=sum(ProbImputeGenos(:,j))/(2*ped%pedigreeSize-ped%nDummys)
                enddo
                open(unit=111,file="." // DASH // "Miscellaneous" // DASH // "MinorAlleleFrequency.txt", status="unknown")


                do j=1,inputParams%nSnpRaw
                    write (111,*) j,Maf(j)
                enddo
                close(111)
            endif

#ifdef DEBUG
            write(0,*) 'DEBUG: Imputation Quality [WriteOutResults]'
#endif

            ImputationQuality(:,1)=sum(2*Maf(:))/inputParams%nSnpRaw
            ImputationQuality(:,2)=0.0
            do i=1, ped%pedigreeSize-ped%nDummys
                if (ped%pedigree(i)%isDummy) then
                    exit
                endif
                do j=1,inputParams%nSnp
                    !ImputationQuality(i,2)=ImputationQuality(i,2)+abs(ProbImputeGenos(i,j)-(2*Maf(j)))
                    ImputationQuality(i,2)=ImputationQuality(i,2)+abs(ProbImputeGenos(i,j)-((Maf(j)**4)+(4*(Maf(j)**2)*((1.0-Maf(j))**2))+(1.0-Maf(j)**4)))
                enddo
                ImputationQuality(i,2)=ImputationQuality(i,2)/inputParams%nSnpRaw
                ImputationQuality(i,3)=float(inputParams%nSnpRaw-count(ImputePhase(i,:,1)==9))/inputParams%nSnpRaw
                ImputationQuality(i,4)=float(inputParams%nSnpRaw-count(ImputePhase(i,:,2)==9))/inputParams%nSnpRaw
                ImputationQuality(i,5)=(ImputationQuality(i,3)+ImputationQuality(i,4))/2
                ImputationQuality(i,6)=float(inputParams%nSnpRaw-count(ImputeGenos(i,:)==9))/inputParams%nSnpRaw
                write (50,'(a20,20000f7.2)') ped%pedigree(i)%originalID,ImputationQuality(i,:)
            enddo

#ifdef DEBUG
            write(0,*) 'DEBUG: Write [WriteOutResults]'
#endif

            do j=1,inputParams%nSnpRaw
                write (51,'(i10,20000f7.2)') j,float(((ped%pedigreeSize-ped%nDummys-(ped%nDummys+1))+1)-count(ImputeGenos(ped%nDummys+1:ped%pedigreeSize-ped%nDummys,j)==9))/((ped%pedigreeSize-ped%nDummys-(ped%nDummys+1))+1)
            enddo

            call CheckImputationInconsistencies(TmpGenos, TmpPhase, ped%pedigreeSize-ped%nDummys, inputParams%nsnp)
            inputParams%WellPhasedThresh=inputParams%WellPhasedThresh/100
            do i=1, ped%pedigreeSize-ped%nDummys
                if (ped%pedigree(i)%isDummy) then
                    exit
                endif
                if (ImputationQuality(i,5)>=inputParams%WellPhasedThresh) then
                    write (52,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(i)%originalID,TmpPhase(i,:,1)
                    write (52,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(i)%originalID,TmpPhase(i,:,2)
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

        close (53)
        close (54)
        close (60)

    end subroutine WriteOutResults

    !#############################################################################################################################################################################################################################

    subroutine ModelRecomb
        use Global
        use informationModule, only : InsteadOfReReadGeneProb
        use alphaimputeinmod
        implicit none

        integer :: e,i,j,k,l,SuperJ,StartDisFound,EndDisFound,HetEnd,HetStart,RSide,LSide,PatMat,SireDamRL,nSnpFinal,Counter
        integer :: StartDisPrev,EndDisPrev,RecombOnOff
        integer :: GamA,GamB,Tmp,StartDisOld,StartDisTmp
        integer :: CountRightSwitch,CountLeftSwitch,PedId,StartDis,EndDis,StartJ,nRec
        integer(kind=1),allocatable,dimension(:,:,:) :: WorkPhase,TempWork
        integer,allocatable,dimension(:) :: WorkLeft,WorkRight,TempVec,StR,EnR,StRNarrow,EnRNarrow
        real,allocatable,dimension(:) :: LengthVec
        real,allocatable,dimension(:,:) :: PatAlleleProb,MatAlleleProb,GeneProbWork
        character(len=7) :: cm
        type(AlphaImputeInput), pointer :: inputParams


        inputParams => defaultInput


        write(cm,'(I7)') inputParams%nSnpRaw !for formatting
        cm = adjustl(cm)

        open (unit=42, file="Results" // DASH // "RecombinationInformation.txt")
        open (unit=43, file="Results" // DASH // "RecombinationInformationNarrow.txt")
        open (unit=44, file="Results" // DASH // "NumberRecombinations.txt")
        open (unit=45, file="Results" // DASH // "RecombinationInformationR.txt")
        open (unit=46, file="Results" // DASH // "RecombinationInformationNarrowR.txt")

        ! Check whether to consider all the raw snps or only the snps left after the edition procedure
        ! If EditedSnpOut in Spec file
        if (inputParams%outopt==0) then
            nSnpFinal=inputParams%nsnp
            ! If AllSnpOut in Spec file
        else
            nSnpFinal=inputParams%nSnpRaw
        endif

        ! Divide haplotypes into chunks of the same length.
        ! Each chunk will be treated separately in different processors
        Tmp=int(float(inputParams%nsnp)/inputParams%nProcessors)
        GpIndex(1,1)=1
        GpIndex(1,2)=Tmp
        if (inputParams%nProcessors>1) then
            do i=2,inputParams%nProcessors
                GpIndex(i,1)=GpIndex(i-1,1)+Tmp
                GpIndex(i,2)=GpIndex(i-1,2)+Tmp
            enddo
        endif
        GpIndex(inputParams%nProcessors,2)=inputParams%nsnp

        allocate(GlobalWorkPhase(0:ped%pedigreeSize-ped%nDummys,inputParams%nsnp,2))
        allocate(WorkPhase(0:ped%pedigreeSize-ped%nDummys,nSnpFinal,2))
        allocate(TempVec(nSnpFinal))
        allocate(LengthVec(nSnpFinal))
        allocate(WorkLeft(nSnpFinal))
        allocate(WorkRight(nSnpFinal))
        allocate(TempWork(0:ped%pedigreeSize-ped%nDummys,nSnpFinal,2))
        allocate(PatAlleleProb(nSnpFinal,2))
        allocate(MatAlleleProb(nSnpFinal,2))
        allocate(GeneProbWork(nSnpFinal,4))
        allocate(StR(nSnpFinal))
        allocate(EnR(nSnpFinal))
        allocate(StRNarrow(nSnpFinal))
        allocate(EnRNarrow(nSnpFinal))

        WorkPhase=9
        if (inputParams%SexOpt==0) then
            if (inputParams%bypassgeneprob==0) then
                ! TODO read in gene prob info here
            else
                call InsteadOfReReadGeneProb
            endif
        endif
        if (inputParams%SexOpt==1) call InsteadOfReReadGeneProb

        l=0
        do j=1,nSnpFinal

            if (SnpIncluded(j)==1) then
                l=l+1
                WorkPhase(:,j,1)=GlobalWorkPhase(:,l,1)
                WorkPhase(:,j,2)=GlobalWorkPhase(:,l,2)
            endif
        enddo

        do i=1,ped%pedigreeSize-ped%nDummys
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
                pedID=ped%pedigree(i)%getSireDamNewIDByIndex(SireDamRL)
                if ((inputParams%SexOpt==1).and.(ped%pedigree(PedId)%gender==inputParams%hetGameticStatus)) cycle
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

                    do while (SuperJ<inputParams%nsnp)
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
                        if (k>1) then
                            if (j==EnR(k-1)) then
                                RecombOnOff=0
                            endif
                        endif
                        if (GamA==9) cycle
                        if (RecombOnOff==1) then
                            Counter=Counter+1
                            if (Counter/=1) then
                                if (ImputePhase(PedId,j,1)/=ImputePhase(PedId,j,2)) then
                                    if ((ImputePhase(PedId,j,1)/=9).and.(ImputePhase(PedId,j,2)/=9)) then
                                        ImputePhase(i,j,e)=9
                                        ImputeGenos(i,j)=9
                                        ProbImputePhase(i,j,e)&
                                            =((1.0-(LengthVec(j)*Counter))*ImputePhase(PedId,j,GamA))&
                                            +(LengthVec(j)*Counter*ImputePhase(PedId,j,GamB))
                                        ProbImputeGenos(i,j)=ProbImputePhase(i,j,1)+ProbImputePhase(i,j,2)
                                    endif
                                endif
                            endif
                        endif
                    enddo
                endif
                write (44,'(a20,20000i20)') ped%pedigree(i)%originalID,nRec
                write (42,'(a20,20000i20)') ped%pedigree(i)%originalID,nRec,StR(1:nRec)
                write (42,'(a20,20000i20)') ped%pedigree(i)%originalID,nRec,EnR(1:nRec)
                !
                !           write (43,'(a20,20000i20)') ped%pedigree(i)%originalID,nRec,StRNarrow(1:nRec)
                !           write (43,'(a20,20000i20)') ped%pedigree(i)%originalID,nRec,EnRNarrow(1:nRec)
                do j=1,nRec
                    write (45,'(a20,20000i20)') ped%pedigree(i)%originalID,e,StR(j),EnR(j)
                    !               write (46,'(a20,20000i20)') ped%pedigree(i)%originalID,e,StRNarrow(j),EnRNarrow(j)
                enddo

            enddo
            WorkPhase(i,:,:)=ImputePhase(i,:,:)

        enddo

        open (unit=33,file="Results" // DASH // "ImputePhase.txt",status="unknown")
        open (unit=34,file="Results" // DASH // "ImputeGenotypes.txt",status="unknown")
        open (unit=40,file="Results" // DASH // "ImputePhaseProbabilities.txt",status="unknown")
        open (unit=41,file="Results" // DASH // "ImputeGenotypeProbabilities.txt",status="unknown")

        do i=1, ped%pedigreeSize-ped%nDummys
            if (ped%pedigree(i)%isDummy) then
                exit
            endif
            write (33,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(i)%originalID,ImputePhase(i,:,1)
            write (33,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(i)%originalID,ImputePhase(i,:,2)
            write (34,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(i)%originalID,ImputeGenos(i,:)
            write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(i)%originalID,ProbImputePhase(i,:,1)
            write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(i)%originalID,ProbImputePhase(i,:,2)
            write (41,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(i)%originalID,ProbImputeGenos(i,:)
        enddo

        print*, " "
        print*, " ","Imputation by detection of recombination events completed"

    end subroutine ModelRecomb

    !#############################################################################################################################################################################################################################

    subroutine IterateGeneProbPhase
        use Global
        use alphaimputeinmod
        implicit none

        integer :: i,j, fileUnit
        type(AlphaImputeInput), pointer :: inputParams

        inputParams => defaultInput

        if (allocated(Maf)) then
            deallocate(maf)
        endif
        allocate(Maf(nSnpIterate))


        if (inputParams%restartOption==4) then
            open (unit=fileUnit,file="Tmp2345678.txt",status="old")
            do i=1,ped%pedigreeSize-ped%nDummys
                read (fileUnit,*) ImputePhase(i,:,1)  
                read (fileUnit,*) ImputePhase(i,:,2)
                read (fileUnit,*) ImputeGenos(i,:)
            enddo
            close (fileUnit)
        endif

        ! StSnp=GpIndex(h,1)
        ! EnSnp=GpIndex(h,2)
        ! TODO once geneprob is MPI'd use this
        open(newunit=fileUnit,file="." // DASH // "Miscellaneous" // "MinorAlleleFrequency.txt", status="unknown")
        do j=1,nSnpIterate
            read (fileUnit,*) Maf(j)
        enddo

        close(fileUnit)


        call IterateMakeGenotype

        ImputePhase(0,:,:)=9
        ImputeGenos(0,:)=9

    end subroutine IterateGeneProbPhase
    !#############################################################################################################################################################################################################################

    subroutine IterateMakeGenotype
        use Global
        implicit none

        integer :: i,j

        do i=1,ped%pedigreeSize-ped%nDummys
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

        do i=1,ped%pedigreeSize-ped%nDummys
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
        use alphaimputeinmod
        implicit none

        integer :: e,i,j,PedLoc, id
        type(AlphaImputeInput), pointer :: inputParams

        inputParams => defaultInput

        if (inputParams%SexOpt==0) then
            do i=1,ped%pedigreeSize-ped%nDummys
                do e=1,2
                    PedLoc=e+1
                    do j=1,nSnpIterate
                        if (ImputePhase(i,j,e)==9) then
                            id = ped%pedigree(i)%getSireDamNewIDByIndex(pedLoc)
                            if ((ImputePhase(id,j,1)==ImputePhase(id,j,2)).and.&
                                (ImputePhase(id,j,1)/=9)) ImputePhase(i,j,e)=ImputePhase(id,j,1)
                        endif
                    enddo
                enddo
            enddo
        else
            do i=1,ped%pedigreeSize-ped%nDummys
                if (ped%pedigree(i)%gender==inputParams%HomGameticStatus) then
                    do e=1,2
                        id = ped%pedigree(i)%getSireDamNewIDByIndex(e+1)
                        do j=1,nSnpIterate
                            if (ImputePhase(i,j,e)==9) then
                                if ((ImputePhase(id,j,1)==ImputePhase(id,j,2)).and.(ImputePhase(id,j,1)/=9)) ImputePhase(i,j,e)=ImputePhase(id,j,1)
                            endif
                        enddo
                    enddo
                else
                    id = ped%pedigree(i)%getSireDamNewIDByIndex(inputParams%HomGameticStatus+1)
                    do j=1,nSnpIterate
                        if (ImputePhase(i,j,1)==9) then     !Comment From John Hickey I changed what was indexed e to 1 I think it is ok same thing in analogous routine
                            !There is no do loop for e here
                            if ((ImputePhase(id,j,1)==ImputePhase(id,j,2)).and.(ImputePhase(id,j,1)/=9)) ImputePhase(i,j,:)=ImputePhase(id,j,1)
                        endif
                    enddo
                endif
            enddo
        endif
        ImputePhase(0,:,:)=9
        ImputeGenos(0,:)=9

    end subroutine IterateParentHomoFill


    !######################################################################################################################################################################################

    subroutine InsteadOfGeneProb
        ! Phase haplotypes whenever there is enough information from the parents (homozygous case)
        use Global
        use AlphaImputeInMod
        implicit none

        integer :: e,i,j,ParId
        type(AlphaImputeInput), pointer :: inputParams

        inputParams => defaultInput
        if (inputParams%SexOpt==1) then                                                     ! Sex chromosome
            allocate(GlobalWorkPhase(0:ped%pedigreeSize-ped%nDummys,inputParams%nsnp,2))
            allocate(ImputeGenos(0:ped%pedigreeSize-ped%nDummys,inputParams%nsnp))
            allocate(ImputePhase(0:ped%pedigreeSize-ped%nDummys,inputParams%nsnp,2))

            GlobalWorkPhase=9
            do i=1,ped%pedigreeSize-ped%nDummys
                do j=1,inputParams%nsnp                                                     ! Phase in the homozygous case
                    
                    if (ped%pedigree(i)%individualGenotype%getGenotype(j)==0) then
                        GlobalWorkPhase(i,j,:)=0
                    else if (ped%pedigree(i)%individualGenotype%getGenotype(j)==2) then
                        GlobalWorkPhase(i,j,:)=1
                    endif
                enddo
                if (ped%pedigree(i)%gender/=inputParams%hetGameticStatus) then                        ! Am I homogametic?
                    do e=1,2                                                    ! Phase a single haplotype whenever my parents are homozygous
                        ParId=ped%pedigree(i)%getSireDamNewIDByIndex(e+1)
                        do j=1,inputParams%nsnp
                            if (ped%pedigree(parId)%individualGenotype%getGenotype(j)==0) then
                                GlobalWorkPhase(i,j,e)=0

                            else if (ped%pedigree(parId)%individualGenotype%getGenotype(j)==2) then
                                GlobalWorkPhase(i,j,e)=1
                            endif
                        enddo
                    enddo
                else                                                            ! Am I heterogametic?
                    ParId=ped%pedigree(i)%getSireDamNewIDByIndex(inputParams%HomGameticStatus+1)
                    do j=1,inputParams%nsnp                                                 ! Phase the two haplotypes whenever my homogametic parent is homozygous
                        if (ped%pedigree(parId)%individualGenotype%getGenotype(j)==0) then
                            GlobalWorkPhase(i,j,:)=0
                        else if (ped%pedigree(parId)%individualGenotype%getGenotype(j)==2) then
                            GlobalWorkPhase(i,j,:)=1
                        endif 
                    enddo
                endif
            enddo

            GlobalWorkPhase(0,:,:)=9
            
            ImputeGenos(:,:)=ped%getGenotypesAsArray()
            ImputePhase(:,:,:)=GlobalWorkPhase(:,:,:)

            allocate(GlobalTmpCountInf(ped%pedigreeSize-ped%nDummys,6))
            GlobalTmpCountInf(:,:)=0

        else                                                                    ! Other chromosome
            allocate(GlobalWorkPhase(0:ped%pedigreeSize-ped%nDummys,inputParams%nsnp,2))
            allocate(ImputeGenos(0:ped%pedigreeSize-ped%nDummys,inputParams%nsnp))
            allocate(ImputePhase(0:ped%pedigreeSize-ped%nDummys,inputParams%nsnp,2))

            GlobalWorkPhase=9
            do i=1,ped%pedigreeSize-ped%nDummys
                do j=1,inputParams%nsnp                                                     ! Phase in the homozygous case
                    if (ped%pedigree(i)%individualGenotype%getGenotype(j)==0) then
                        GlobalWorkPhase(i,j,:)=0
                    else if (ped%pedigree(i)%individualGenotype%getGenotype(j)==2) then
                        GlobalWorkPhase(i,j,:)=1
                    endif
                enddo
                if (associated(ped%pedigree(i)%sirePointer)) then
                    if (ped%pedigree(i)%sirePointer%isDummy) exit
                    ParId=ped%pedigree(i)%sirePointer%id

                    do j=1,inputParams%nsnp                                                     ! Phase if my father is homozygous
                        if (ped%pedigree(ParId)%individualGenotype%getGenotype(j)==0) then
                            GlobalWorkPhase(i,j,1)=0
                        else if (ped%pedigree(ParId)%individualGenotype%getGenotype(j)==2) then
                            GlobalWorkPhase(i,j,1)=1
                        endif 
                    enddo
                endif
                if (associated(ped%pedigree(i)%damPointer)) then
                    if (ped%pedigree(i)%damPointer%isDummy) exit
                    ParId=ped%pedigree(i)%damPointer%id
                    do j=1,inputParams%nsnp                                                     ! Phase if my mother is homozygous
                        if (ped%pedigree(ParId)%individualGenotype%getGenotype(j)==0) then
                            GlobalWorkPhase(i,j,2)=0
                        else if (ped%pedigree(ParId)%individualGenotype%getGenotype(j)==2) then
                            GlobalWorkPhase(i,j,2)=1
                        endif 
                    enddo
                endif

            enddo

            GlobalWorkPhase(0,:,:)=9
            ImputeGenos(:,:)=ped%getGenotypesAsArray()
            ImputePhase(:,:,:)=GlobalWorkPhase(:,:,:)

        endif

    end subroutine InsteadOfGeneProb

    !######################################################################################################################################################################################




    !#############################################################################################################################################################################################################################

    subroutine MakeFiles
        use Global

        use alphaimputeinmod
        implicit none

        integer :: i
        integer,dimension(:), allocatable :: TempCore,TempCplusT
        integer :: Tmp
        character(len=7) :: cm
        type(AlphaImputeInput), pointer :: inputParams


        inputParams => defaultInput

        allocate(TempCore(inputParams%nPhaseInternal))
        allocate(TempCplusT(inputParams%nPhaseInternal))

        allocate(GpIndex(inputParams%nprocessors,2))

        ! WARNING: This code is not necessary
        write(cm,'(I7)') inputParams%nSnpRaw !for formatting
        cm = adjustl(cm)
        ! end WARNING

        do i=1,inputParams%nPhaseExternal
            TempCore(i)=inputParams%CoreLengths(i)
            TempCore(i+inputParams%nPhaseExternal)=inputParams%CoreLengths(i)
            TempCplusT(i)=inputParams%CoreAndTailLengths(i)
            TempCplusT(i+inputParams%nPhaseExternal)=inputParams%CoreAndTailLengths(i)
            if (TempCore(i)>inputParams%nsnp) TempCore(i)=inputParams%nsnp
            if (TempCore(i+inputParams%nPhaseExternal)>inputParams%nsnp) TempCore(i+inputParams%nPhaseExternal)=inputParams%nsnp
            if (TempCplusT(i)>inputParams%nsnp) TempCplusT(i)=inputParams%nsnp
            if (TempCplusT(i+inputParams%nPhaseExternal)>inputParams%nsnp) TempCplusT(i+inputParams%nPhaseExternal)=inputParams%nsnp
        enddo

        open(unit=103,file="." // DASH // "InputFiles" // DASH // "AlphaPhaseInputPedigree.txt", status="unknown")
        open(unit=104,file="." // DASH // "InputFiles" // DASH // "RecodedGeneProbInput.txt", status="unknown")
        open(unit=105,file="." // DASH // "InputFiles" // DASH // "AlphaPhaseInputGenotypes.txt", status="unknown")


        do i=1,ped%pedigreeSize-ped%nDummys
            write (104,'(i16,1x,i16,1x,i16,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(i)%getIntegerVectorOfRecodedIdsNoDummy(), ped%pedigree(i)%individualGenotype%toIntegerArray()
            if (Setter(i)==1) write (105,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(i)%originalID,ped%pedigree(i)%individualGenotype%toIntegerArray()
        enddo
        call flush(104)
        call flush(105)
        close(104)
        close(105)

        CountRawGenos=ped%getNumGenotypesMissing()

        do i=1,nObsDataRaw
            write (103,'(4a20)') Ped%pedigree(i)%originalID,Ped%pedigree(i)%sireID,Ped%pedigree(i)%damID
        enddo
        call flush(103)
        close(103)

        if (inputParams%SexOpt==1) then                 ! Sex chromosome
            call InsteadOfGeneProb          ! Calculate Genotype Probabilities
        else                                ! Not sex chromosome
            if (inputParams%bypassgeneprob==1 .OR. inputParams%hmmoption == RUN_HMM_ONLY) then
                call InsteadOfGeneProb
            endif
        endif
        Tmp=int(float(inputParams%nsnp)/inputParams%nprocessors)
        GpIndex(1,1)=1
        GpIndex(1,2)=Tmp
        if (inputParams%nprocessors>1) then
            do i=2,inputParams%nprocessors
                GpIndex(i,1)=GpIndex(i-1,1)+Tmp
                GpIndex(i,2)=GpIndex(i-1,2)+Tmp
            enddo
        endif
        GpIndex(inputParams%nprocessors,2)=inputParams%nsnp


        deallocate(TempCore)
        deallocate(TempCplusT)
    end subroutine MakeFiles

    !#############################################################################################################################################################################################################################
    subroutine ClassifyAnimByChips
        ! Classify animals according to the HD chip information with a margin of missing markers
        ! The condition for an animal to be classify with a particular HD snp chip is:
        ! If the missing markers are not above a threshold and the number of markers is below the
        ! nominal number of markers for that chip
        ! LD animals or HD animals with missing markers above a threshold are not assign to any
        ! snp panel

        use Global

        use alphaimputeinmod
        use ISO_Fortran_Env
        implicit none

        integer :: i, j, CountMiss, UOutputs
        logical, allocatable :: printed(:)
        type(AlphaImputeInput), pointer :: inputParams

        inputParams => defaultInput
        open(newunit=UOutputs, file="." // DASH // "Miscellaneous" // DASH // "SnpCallRateByAnimalByChip.txt",status='unknown')
        allocate(animChip(ped%pedigreeSize-ped%nDummys))
        animChip(:)=0

        allocate(printed(ped%pedigreeSize-ped%nDummys))
        printed=.FALSE.

        do i=1,ped%pedigreeSize-ped%nDummys
            CountMiss=count(TempGenos(i,:)==9)
            do j=1,inputParams%MultiHD
                if ( (CountMiss-(inputParams%nsnp-inputParams%nSnpByChip(j))) < (1.0-inputParams%PercGenoForHD)*inputParams%nSnpByChip(j)&
                    .and. (inputParams%nsnp-CountMiss)<inputParams%nSnpByChip(j)&
                    .and. ped%pedigree(i)%genotyped) then
                    animChip(i)=j
                    write(UOutputs,'(a20,6f5.1)') ped%pedigree(i)%originalID, (inputParams%nsnp-CountMiss)*100/real(inputParams%nSnpByChip(j))
                    exit
                endif
                if ((CountMiss-(inputParams%nsnp-inputParams%nSnpByChip(j))) > (1.0-inputParams%PercGenoForHD)*inputParams%nSnpByChip(j)&
                    ! .and. animChip(i)/=0&
                .and. printed(i)==.false.&
                    .and. ped%pedigree(i)%genotyped) Then
                    write(UOutputs,'(a20,6f5.1)') ped%pedigree(i)%originalID, (inputParams%nsnp-CountMiss)*100/real(inputParams%nSnpByChip(j))
                    printed(i)=.true.
                end if
            enddo
        enddo
        close(UOutputs)
    end subroutine ClassifyAnimByChips

    !#############################################################################################################################################################################################################################
    subroutine InternalEdit
        use Global
        use alphaimputeinmod
        implicit none

        integer :: i,j,k,CountMiss,CountHD,nSnpR,dum
        real, allocatable, dimension(:) :: SnpSummary, TempFreq,Counter
        character (len=300) :: dumC,FileName
        type(AlphaImputeInput), pointer :: inputParams

        inputParams => defaultInput

        allocate(SnpSummary(inputParams%nsnp))
        allocate(TempFreq(inputParams%nsnp))
        allocate(Counter(inputParams%nsnp))
        allocate(SnpIncluded(inputParams%nsnp))

        allocate(Setter(0:ped%pedigreeSize-ped%nDummys))

        SnpIncluded(:)=0
        if ((inputParams%managephaseon1off0==0).and.(inputParams%NoPhasing==1)) then
            FileName = trim(inputParams%phasePath) // DASH // "EditingSnpSummary.txt"
            open(unit=111,file=trim(FileName),status="old")
            do i=1,inputParams%nsnp
                read (111,*) dum,SnpSummary(i),SnpIncluded(i)
            enddo
        endif

        if (inputParams%NoPhasing==0) SnpIncluded(:)=1

        ! I user do not specify any file with HD individuals
        if (inputParams%UserDefinedHD==0) then
            Setter(0)=0
            Setter(1:ped%pedigreeSize-ped%nDummys)=1
            RecIdHDIndex(0)=0
            RecIdHDIndex(1:ped%pedigreeSize-ped%nDummys)=1
            do i=1,ped%pedigreeSize-ped%nDummys
                CountMiss=count(TempGenos(i,:)==9)
                if (inputParams%MultiHD/=0) then
                    ! Disregard animals at LD or those HD animals with a number of markers missing
                    if (animChip(i)==0) then
                        Setter(i)=0
                        RecIdHDIndex(i)=0
                    endif
                else
                    if ((float(CountMiss)/inputParams%nsnp)>(1.0-inputParams%PercGenoForHD)) then
                        Setter(i)=0
                        RecIdHDIndex(i)=0
                    endif
                endif
            enddo
            CountHD=count(Setter(:)==1)
        else                                ! User has specified HD individuals
            Setter(0)=0
            Setter(1:ped%pedigreeSize-ped%nDummys)=0
            RecIdHDIndex(0)=0
            RecIdHDIndex(1:ped%pedigreeSize-ped%nDummys)=0

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

            block
                integer :: tmpID
                do k=1,CountHD
                    read (46,*) dumC

                    tmpID = ped%dictionary%getValue(dumC)
                    if (tmpID /= DICT_NULL) then
                        Setter(tmpID)=1
                        RecIdHDIndex(tmpID)=1
                        exit
                    endif
                enddo
            end block
            CountHD=count(Setter(:)==1)
            print*, " "
            print*, " ",CountHD," valid indiviudals in the user specified AlphaPhase file"
        endif

        open (unit=102,file="." // DASH // "Miscellaneous" // DASH // "EditingSnpSummary.txt",status="unknown")

        if (inputParams%managephaseon1off0==1) then
            TempFreq(:)=0.0
            Counter(:)=0
            SnpSummary=0.0
            do j=1,inputParams%nsnp
                do i=1, ped%pedigreeSize-ped%nDummys
                    if (TempGenos(i,j)/=9) then
                        TempFreq(j)=TempFreq(j)+float(TempGenos(i,j))
                        Counter(j)=Counter(j)+2
                    else
                        if (Setter(i) == 1) then
                            SnpSummary(j)=SnpSummary(j)+1.0
                        endif
                    end if
                end do
                if (Counter(j)>0.000000) TempFreq(j)=TempFreq(j)/Counter(j)
            end do
        endif

        if (inputParams%MultiHD/=0 .or. inputParams%IntEditStat==0) then
            nSnpR=inputParams%nsnp
            if (inputParams%managephaseon1off0==1) SnpIncluded(:)=1
        else
            if (inputParams%managephaseon1off0==1) then
                SnpSummary(:)=SnpSummary(:)/CountHD
                nSnpR=0
                do j=1,inputParams%nsnp
                    if ((SnpSummary(j)<inputParams%PercSnpMiss).and.((TempFreq(j)>0.00000001).and.(TempFreq(j)<0.9999999))) nSnpR=nSnpR+1
                enddo
            else
                nSnpR=count(SnpIncluded(:)==1)
            endif
            if (nSnpR==inputParams%nsnp) then
                SnpIncluded(:)=1
            else
                if (inputParams%managephaseon1off0==1) then
                    k=0
                    do j=1,inputParams%nsnp
                        if ((SnpSummary(j)<inputParams%PercSnpMiss).and.((TempFreq(j)>0.00000001).and.(TempFreq(j)<0.9999999))) then
                            k=k+1
                            SnpIncluded(j)=1
                        endif
                    enddo
                    inputParams%nsnp=nSnpR
                endif
            endif
            if (inputParams%UserDefinedHD==0) then
                Setter(1:ped%pedigreeSize-ped%nDummys)=1
                RecIdHDIndex(1:ped%pedigreeSize-ped%nDummys)=1
                do i=1,ped%nGenotyped
                    CountMiss=ped%pedigree(ped%genotypeMap(i))%individualGenotype%numMissing()
                    if ((float(CountMiss)/inputParams%nsnp)>(1.0-inputParams%SecondPercGenoForHD)) then
                        Setter(i)=0
                        RecIdHDIndex(i)=0
                    endif
                enddo
                CountHD=count(Setter(:)==1)
            else
                do i=1,ped%nGenotyped
                    if (Setter(i)==1) then
                        CountMiss=ped%pedigree(ped%genotypeMap(i))%individualGenotype%numMissing()
                        if ((float(CountMiss)/inputParams%nsnp)>(1.0-inputParams%SecondPercGenoForHD)) then
                            Setter(i)=0
                            RecIdHDIndex(i)=0
                        endif
                    endif
                enddo
                CountHD=count(Setter(:)==1)
            endif
        endif

        do j=1,inputParams%nSnpRaw
            write (102,*) j,SnpSummary(j),SnpIncluded(j)        !'(i,1x,f5.3,1x,i)'
        enddo
        close(102)

        if ((inputParams%ManagePhaseOn1Off0==0).and.(inputParams%NoPhasing==1)) then
            FileName = trim(inputParams%PhasePath) // DASH // "EditingSnpSummary.txt"
        else
            FileName = "." // DASH // "Phasing" // DASH // "EditingSnpSummary.txt"
        end if

        open (unit=112,file=FileName,status="unknown")
        do j=1,inputParams%nSnp
            write (112,*) j,SnpSummary(j),SnpIncluded(j)        !'(i,1x,f5.3,1x,i)'
        enddo
        close(112)

        print*, " "
        print*, " "
        print*, " ",CountHD," indiviudals passed to AlphaPhase"
        print*, " ",inputParams%nsnp," snp remain after editing"

        ! TODO clean this whole subroutine
        

        ! we add animals to hd list here
        do i=1, ped%pedigreeSize - ped%nDummys
            if (setter(i) == 1) then
                call ped%setAnimalAsHD(i)
            endif
        enddo 
        deallocate(SnpSummary)
        deallocate(TempFreq)
        deallocate(Counter)
    end subroutine InternalEdit

    !#############################################################################################################################################################################################################################

    subroutine FillInBasedOnOffspring
        ! Genotype SNPs based on the genetic information of my offsprings
        use Global
        use alphaimputeinmod
        implicit none

        integer :: i,j,k
        integer, allocatable, dimension(:) :: count0,count1,count2
        type(AlphaImputeInput), pointer :: inputParams
        type(individual) ,pointer :: tmpOff
        inputParams => defaultInput

        allocate(Count0(inputParams%nsnp))
        allocate(Count1(inputParams%nsnp))
        allocate(Count2(inputParams%nsnp))
        do i=1,ped%pedigreeSize-ped%nDummys ! These are parents
            ! This three variables will count the different number of genotypes of the offsprings
            Count0=0
            Count1=0
            Count2=0
            do j=1,ped%pedigree(i)%nOffs ! These are offsprings
                tmpOff => ped%pedigree(i)%offsprings(j)%p
                if ((inputParams%SexOpt==1).and.(ped%pedigree(i)%gender==inputParams%hetGameticStatus).and.(tmpOff%gender==inputParams%hetGameticStatus)) cycle
                do k=1,inputParams%nsnp
                    if (TempGenos(i,k)==9) then ! If my parent is not genotyped
                        if (TempGenos(tmpOff%id,k)==0) Count0(k)=Count0(k)+1    ! Number of offspring genotype as 0
                        if (TempGenos(tmpOff%id,k)==1) Count1(k)=Count1(k)+1    ! Number of offspring genotype as 1
                        if (TempGenos(tmpOff%id,k)==2) Count2(k)=Count2(k)+1    ! Number of offspring genotype as 2
                    endif
                enddo

            enddo

            do k=1,inputParams%nsnp
                if ((Count0(k)+Count1(k)+Count2(k))>OffspringFillMin) then
                    if (Count0(k)==0) TempGenos(i,k)=2                       ! This is the most likely thing, but it might be not true
                    if (Count2(k)==0) TempGenos(i,k)=0                       ! This is the most likely thing, but it might be not true
                    if ((inputParams%SexOpt==1).and.(ped%pedigree(i)%gender==inputParams%hetGameticStatus)) cycle
                    if ((Count0(k)>2).and.(Count2(k)>2)) TempGenos(i,k)=1    ! This is the most likely thing, but it might be not true
                endif
            enddo
        enddo

    end subroutine FillInBasedOnOffspring

    !#############################################################################################################################################################################################################################

    subroutine FillInSnp
        ! Genotype SNPs based on the pedigree information
        use Global
        use AlphaImputeInMod
        implicit none

        integer :: i,j,k,TurnOn
        type(AlphaImputeInput), pointer :: inputParams
        integer :: tmpParentId
        inputParams => defaultInput

        do i=1,ped%pedigreeSize-ped%nDummys
            do k=2,3
                TurnOn=1
                tmpParentId = ped%pedigree(i)%getSireDamNewIDByIndex(k)
                ! if the proband is heterogametic, and
                ! considering the heterogametic parent, then avoid!!
                if ((inputParams%SexOpt==1).and.(ped%pedigree(i)%gender==inputParams%hetGameticStatus).and. ((k-1)==inputParams%hetGameticStatus)) TurnOn=0
                if (tmpParentId /= 0) then
                    if (ped%pedigree(tmpParentId)%isDummy) cycle
                    ! if ((inputParams%SexOpt==1).and.(ped%pedigree(i)%gender==inputParams%hetGameticStatus).and. (pedigree%(i)%parent(k-1)%gender==inputParams%hetGameticStatus)) TurnOn=0
                    ! Homogametic individuals and the homogametic parent of a heterogametic individual
                    if (TurnOn==1) then
                        do j=1,inputParams%nsnp
                            if ((TempGenos(i,j)==0).and.(TempGenos(tmpParentId,j)==2)) then
                                TempGenos(i,j)=9
                                TempGenos(ped%pedigree(i)%getSireDamNewIDByIndex(k),j)=9
                            endif
                            if ((TempGenos(i,j)==2).and.(TempGenos(tmpParentId,j)==0)) then
                                TempGenos(i,j)=9
                                TempGenos(tmpParentId,j)=9
                            endif
                        enddo
                    endif
                endif
            enddo
        enddo

        ! WARNING: This can be refactored
        do i=1,ped%pedigreeSize-ped%nDummys
            do j=1,inputParams%nsnp
                if (TempGenos(i,j)==9 .and. .not. ped%pedigree(i)%hasDummyParent()) then
                    if ((TempGenos(ped%pedigree(i)%getSireDamNewIDByIndex(2),j)==0).and.(TempGenos(ped%pedigree(i)%getSireDamNewIDByIndex(3),j)==0)) then
                        TempGenos(i,j)=0
                    else if ((TempGenos(ped%pedigree(i)%getSireDamNewIDByIndex(2),j)==2).and.(TempGenos(ped%pedigree(i)%getSireDamNewIDByIndex(3),j)==2)) then
                        TempGenos(i,j)=2
                    endif
                    if (inputParams%SexOpt==1) then
                        if (ped%pedigree(i)%gender/=inputParams%hetGameticStatus) then
                            if ((TempGenos(ped%pedigree(i)%getSireDamNewIDByIndex(2),j)==0).and.(TempGenos(ped%pedigree(i)%getSireDamNewIDByIndex(3),j)==2)) TempGenos(i,j)=1
                            if ((TempGenos(ped%pedigree(i)%getSireDamNewIDByIndex(2),j)==2).and.(TempGenos(ped%pedigree(i)%getSireDamNewIDByIndex(3),j)==0)) TempGenos(i,j)=1
                        else
                            ! HomGameticSatus(1 or 2) +1 = sire (2) or dam (3)
                            if (TempGenos(ped%pedigree(i)%getSireDamNewIDByIndex(inputParams%HomGameticStatus+1),j)==0) TempGenos(i,j)=0
                            if (TempGenos(ped%pedigree(i)%getSireDamNewIDByIndex(inputParams%HomGameticStatus+1),j)==2) TempGenos(i,j)=2
                        endif
                    else
                        if ((TempGenos(ped%pedigree(i)%getSireDamNewIDByIndex(2),j)==0).and.(TempGenos(ped%pedigree(i)%getSireDamNewIDByIndex(3),j)==2)) TempGenos(i,j)=1
                        if ((TempGenos(ped%pedigree(i)%getSireDamNewIDByIndex(2),j)==2).and.(TempGenos(ped%pedigree(i)%getSireDamNewIDByIndex(3),j)==0)) TempGenos(i,j)=1
                    endif
                endif
            enddo
        enddo

    end subroutine FillInSnp

    !#############################################################################################################################################################################################################################
    !#############################################################################################################################################################################################################################

    subroutine CheckParentage
        ! Check which is the parentage relation between animals.
        ! Set the vector baseline that specify whether an animal
        ! has no parents or its parents have been pruned

        use alphaimputeinmod
        use Global

        implicit none

        integer :: e,i,j,CountBothGeno,CountDisagree,CountChanges,IndId,ParId,ParPos
        integer :: TurnOn
        type(AlphaImputeInput), pointer :: inputParams
        integer :: tmpID,dumId
        integer :: nHomoParent, nBothHomo

        inputParams => defaultInput
        open (unit=101,file="." // DASH // "Miscellaneous" // DASH // "PedigreeMistakes.txt",status="unknown")




        CountChanges=0
        nHomoParent = 0
        nBothHomo = 0
        do e=1,2                    ! Do whatever this does, first on males and then on females
            ParPos=e+1              ! Index in the Genotype and Pedigree matrices for sires and dams
            do i=1,ped%pedigreeSize
                IndId=ped%pedigree(i)%id           ! My Id
                ParId=ped%pedigree(i)%getSireDamNewIDByIndex(e+1)       ! Paternal Id,
                TurnOn=1
                ! GenderRaw: Says whether the proband is a male or a female
                ! If there are sex cromosome information, and
                ! if the proband is heterogametic, and
                ! I am considering the heterogametic parent, then avoid!!
                ! That is, avoid males and their sires (hetero=1), or females and their dams (hetero=2)
                if ((inputParams%SexOpt==1).and.(ped%pedigree(i)%gender==inputParams%hetGameticStatus).and.((ParPos-1)==inputParams%hetGameticStatus)) TurnOn=0

                ! Consider the Homogametic probands and the heterogametic proband with homogametic parent
                if (parId /= 0 ) then
                    if ((ped%pedigree(IndId)%genotyped.and. ped%pedigree(parId)%genotyped).and.(TurnOn==1)) then
                        CountBothGeno=0
                        CountDisagree=0
                        nHomoParent = 0
                        nBothHomo = 0

                        ! Look for mendelenian errors
                        do j=1,inputParams%nsnp
                            
                            ! TODO can probably do the following better with new genotype class
                            if ((ped%pedigree(IndId)%individualGenotype%getGenotype(j)/=9).and.(ped%pedigree(parId)%individualGenotype%getGenotype(j)/=9)) then
                                CountBothGeno=CountBothGeno+1
                                if (ped%pedigree(parId)%individualGenotype%getGenotype(j)/=1) then
                                    nHomoParent = nHomoParent+1
                                    if ( ped%pedigree(IndId)%individualGenotype%getGenotype(j)/=1) then
                                        nBothHomo =  nBothHomo + 1
                                    end if
                                end if
                                if ((ped%pedigree(IndId)%individualGenotype%getGenotype(j)==0).and.(ped%pedigree(parId)%individualGenotype%getGenotype(j)==2)) CountDisagree=CountDisagree+1
                                if ((ped%pedigree(IndId)%individualGenotype%getGenotype(j)==2).and.(ped%pedigree(parId)%individualGenotype%getGenotype(j)==0)) CountDisagree=CountDisagree+1
                            endif
                        enddo
                        if ((float(CountDisagree)/CountBothGeno)>DisagreeThreshold) then ! Mendelenian error
                            write (101,'(2a20,4I,3f5.3)') &
                                Ped%pedigree(i)%originalID, Ped%pedigree(i)%getSireDamByIndex(ParPos), CountDisagree, CountBothGeno, nHomoParent, nBothHomo, &
                                float(CountDisagree)/CountBothGeno, float(CountDisagree)/nHomoParent, float(CountDisagree)/nBothHomo
                            CountChanges=CountChanges+1
                            if (parPos == 2) then
                                call ped%pedigree(i)%sirePointer%removeOffspring(ped%pedigree(i))
                                call ped%createDummyAnimalAtEndOfPedigree(dumId, i)
                                ! ped%pedigree(i)%sireID = '0'
                            else if (parPos == 3) then
                                call ped%pedigree(i)%damPointer%removeOffspring(ped%pedigree(i))
                                call ped%createDummyAnimalAtEndOfPedigree(dumId, i)
                                ! ped%pedigree(i)%damPointer => null()
                                ! ped%pedigree(i)%damId = '0'
                            endif
                        else
                            ! Remove genotype of proband and parent
                            do j=1,inputParams%nsnp
                                if ((ped%pedigree(IndId)%individualGenotype%getGenotype(j)/=9).and.(ped%pedigree(parId)%individualGenotype%getGenotype(j)/=9)) then
                                    if ((ped%pedigree(IndId)%individualGenotype%getGenotype(j)==0).and.(ped%pedigree(parId)%individualGenotype%getGenotype(j)==2)) then
                                        call ped%pedigree(IndId)%individualGenotype%setGenotype(j,9)
                                        call ped%pedigree(parId)%individualGenotype%setGenotype(j,9)
                                    endif
                                    if ((ped%pedigree(IndId)%individualGenotype%getGenotype(j)==2).and.(ped%pedigree(parId)%individualGenotype%getGenotype(j)==0)) then
                                        call ped%pedigree(IndId)%individualGenotype%setGenotype(j,9)
                                        call ped%pedigree(parId)%individualGenotype%setGenotype(j,9)
                                    endif
                                endif
                            enddo
                        endif
                    endif
                endif
            enddo
        enddo
        write (101,*) CountChanges," changes were made to the pedigree"
        close (101)
        print*, " ",CountChanges," errors in the pedigree due to Mendelian inconsistencies"

        ! Sort sires and dams, and look for mistakes (bisexuality,...).

        allocate(RecIdHDIndex(0:ped%pedigreeSize-ped%nDummys))

        RecIdHDIndex=0


        call ped%outputSortedPedigreeInAlphaImputeFormat("." // DASH // "Miscellaneous" // DASH // "InternalDataRecoding.txt")


        call ped%sortPedigreeAndOverwrite()
        if (inputParams%SexOpt==1) then
            do j=1,nAnisInGenderFile
                tmpId = ped%dictionary%getValue(genderId(j))
                if (tmpID /= dict_null) then
                    ped%pedigree(tmpId)%gender=GenderRaw(j)
                    TurnOn=0
                endif
            enddo
        endif

        allocate(TempGenos(0:ped%pedigreeSize,inputParams%nsnp))
        tempGenos(0,:) = 9
        tempGenos(1:ped%nGenotyped,:) = ped%getGenotypesAsArray()

    end subroutine CheckParentage
    !#############################################################################################################################################################################################################################

    subroutine MakeDirectories(HMM)
        use global
        use GlobalVariablesHmmMaCH
        use alphaimputeinmod

        implicit none

        integer, intent(in) :: HMM
        integer :: i
        character(len=300) :: FolderName
        type(AlphaImputeInput), pointer :: inputParams

        inputParams => defaultInput

        if (HMM == RUN_HMM_NGS) then
            ! call rmdir("Results")
            ! call rmdir("Miscellaneous")
            call system(RMDIR // " Results")
            call system(RMDIR // " Miscellaneous")
            ! call system("mkdir Results")
            ! call system("mkdir Miscellaneous")
            call system(MD // " Results")
            call system(MD // " Miscellaneous")

        else
            print*, ""
            call system(RMDIR // " Miscellaneous")
            call system(RMDIR // " Phasing")
            call system(RMDIR // " Results")
            call system(RMDIR // " InputFiles")
            call system(RMDIR // " GeneProb")
            call system(RMDIR // " IterateGeneProb")

            call system(MD // " Phasing")
            call system(MD // " Miscellaneous")
            call system(MD // " Results")
            call system(MD // " InputFiles")
            call system(MD // " GeneProb")
            call system(MD // " IterateGeneProb")

            do i=1,inputParams%nprocessors
                write (FolderName,'("GeneProb"i0)')i
                ! call system ("mkdir GeneProb/" // FolderName)       !here
                call system(MD // " GeneProb" // DASH // FolderName)
            enddo
            do i=1,inputParams%nPhaseInternal
                write (FolderName,'("Phase"i0)')i
                ! call system ("mkdir Phasing/" // FolderName)
                call system(MD // " Phasing" // DASH // FolderName)
            enddo
            do i=1,inputParams%nprocessors
                write (FolderName,'("GeneProb"i0)')i            !here
                ! call system ("mkdir IterateGeneProb/" // FolderName)    !here
                call system(MD // " IterateGeneProb" // DASH // FolderName)
            enddo
        endif


    end subroutine MakeDirectories

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



    !#############################################################################################################################################################################################################################

    subroutine FinalChecker
        use Global

        use Utils
        use alphahouseMod, only : countLines
        use alphastatmod, only : pearsn, moment
        use alphaimputeinmod
        implicit none

        integer :: h,i,j,k,l,nAnisTest,CountCatTest(6),tmpIDInt
        integer :: SummaryStats(3,6),Div,CountLen,Counter,Top1,Top2,Top3,Top4,Bot,ContSnpCor,CountValAnim(6)
        double precision :: SummaryProps(3,6),SumPat(6),SumMat(6),MeanCorPerInd(6),StdDevPerGrp(6),AveCategoryInformativeness(6,6)
        double precision :: Tmpave,Tmpadev,Tmpvar,Tmpskew,Tmpcurt
        character(len=300) :: Names(6),FileName,dumC
        integer,allocatable,dimension(:) :: RecTestId,FinalSetter,WorkTmp,GenoStratIndex,Work
        integer,allocatable,dimension(:,:) :: TrueGenos,RawGenos,TestMat,TestAnimInformativeness
        double precision,allocatable,dimension(:) :: Correlations,CorrelationPerAnimal,TmpVarPerGrp
        double precision,allocatable,dimension(:,:) :: AnisSummary,WorkVec,RealTestGenos,CalcCorPerAnimal
        type(AlphaImputeInput), pointer :: inputParams

        character(len=lengan),allocatable,dimension(:) :: TrueGenosId

        inputParams =>defaultInput

        allocate(Work(inputParams%nSnpRaw))
        allocate(WorkTmp(inputParams%nSnpRaw))
        allocate(GenoStratIndex(ped%pedigreeSize-ped%nDummys))

        FileName=trim(inputParams%TrueGenotypeFile)
        ! call CountLines(FileName,nAnisTest)
        nAnisTest = CountLines(FileName)

        ! if (WindowsLinux==1) then
        !      call system("rmdir /s /q TestAlphaImpute")
        ! else
        !     call rmdir("TestAlphaImpute")
        ! endif
        ! call system("mkdir TestAlphaImpute")

        call system(RMDIR // " TestAlphaImpute")
        call system(MD // " TestAlphaImpute")

        open (unit=35,file=trim(inputParams%TrueGenotypeFile),status="old")
        open (unit=36,file=trim(inputParams%GenotypeFile),status="unknown")
        ! open (unit=37,file="./TestAlphaImpute/IndividualAnimalAccuracy.txt",status="unknown")
        ! open (unit=38,file="./TestAlphaImpute/SummaryAnimalAccuracy.txt",status="unknown")
        ! open (unit=44,file="./TestAlphaImpute/IndividualSummaryAccuracy.txt",status="unknown")
        ! open (unit=45,file="./TestAlphaImpute/IndividualSummaryYield.txt",status="unknown")
        ! open (unit=48,file="./TestAlphaImpute/IndividualSnpAccuracy.txt",status="unknown")
        open (unit=37, file="." // DASH // "TestAlphaImpute" // DASH // "IndividualAnimalAccuracy.txt", status="unknown")
        open (unit=38, file="." // DASH // "TestAlphaImpute" // DASH // "SummaryAnimalAccuracy.txt", status="unknown")
        open (unit=44, file="." // DASH // "TestAlphaImpute" // DASH // "IndividualSummaryAccuracy.txt", status="unknown")
        open (unit=45, file="." // DASH // "TestAlphaImpute" // DASH // "IndividualSummaryYield.txt", status="unknown")
        open (unit=48, file="." // DASH // "TestAlphaImpute" // DASH // "IndividualSnpAccuracy.txt", status="unknown")

        Names(1)="Both Parents Genotyped"
        Names(2)="Sire and Maternal GrandSire Genotyped"
        Names(3)="Dam and Paternal Grandsire Genotyped"
        Names(4)="Sire Genotyped"
        Names(5)="Dam Genotyped"
        Names(6)="Other Relatives Genotyped"

        if (allocated(GlobalTmpCountInf)==.FALSE.) then
            allocate(GlobalTmpCountInf(ped%pedigreeSize-ped%nDummys,6))
            GlobalTmpCountInf(:,:)=0
        endif

        allocate(FinalSetter(0:ped%pedigreeSize-ped%nDummys))
        FinalSetter=0
        do i=1,ped%nGenotyped
            read (36,*) dumC,WorkTmp(:)
            Counter=0
            do j=1,inputParams%nSnpRaw
                if ((WorkTmp(j)>=0).and.(WorkTmp(j)<=2)) Counter=Counter+1
            enddo
            if (float(Counter)>(float(inputParams%nSnpRaw)/2)) then
                tmpIDInt = ped%dictionary%getValue(dumC)
                if (tmpIDInt /= DICT_NULL) then
                    FinalSetter(tmpIDInt)=1
                endif

            endif
        enddo
        rewind(36)

        allocate(TestAnimInformativeness(nAnisTest,6))

        if (inputParams%outopt==0) then

            allocate(TrueGenos(nAnisTest,inputParams%nsnp))
            allocate(TrueGenosId(nAnisTest))
            allocate(RawGenos(nAnisTest,inputParams%nsnp))
            allocate(TestMat(nAnisTest,inputParams%nsnp))
            allocate(RecTestId(nAnisTest))
            allocate(AnisSummary(nAnisTest,5))
            allocate(Correlations(6))
            allocate(RealTestGenos(nAnisTest,inputParams%nsnp))
            allocate(CalcCorPerAnimal(inputParams%nsnp,2))
            allocate(CorrelationPerAnimal(nAnisTest))
            allocate(TmpVarPerGrp(nAnisTest))

            do i=1,nAnisTest
                read (35,*) TrueGenosId(i),Work(:)
                k=0
                do j=1,inputParams%nSnpRaw
                    if (SnpIncluded(j)/=0) then
                        k=k+1
                        TrueGenos(i,k)=Work(j)
                    endif
                enddo
            enddo

            RecTestId(:)=-99
            do i=1,nAnisTest
                do j=1,ped%pedigreeSize-ped%nDummys
                    if (trim(TrueGenosId(i))==trim((ped%pedigree(j)%originalID))) then
                        RecTestId(i)=j
                        TestAnimInformativeness(i,:)=GlobalTmpCountInf(j,1:6)
                        exit
                    endif
                enddo
            enddo
            if (count(RecTestId(:)==-99)>0) print*, "Error - There seems to be unidentifiablecount ",count(RecTestId(:)==-99)," individuals in the test file"

            do i=1,ped%nGenotyped
                read (36,*) dumC,WorkTmp(:)
                do j=1,nAnisTest
                    if (trim(TrueGenosId(j))==dumC) then
                        k=0
                        do l=1,inputParams%nSnpRaw
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
            do i=1,ped%pedigreeSize-ped%nDummys
                if (FinalSetter(i)/=1) then
                    GenoStratIndex(i)=6
                    if (FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(3))==1) then
                        GenoStratIndex(i)=5
                        if (FinalSetter(ped%pedigree(i)%getPaternalGrandSireRecodedIndex())==1) then
                            GenoStratIndex(i)=3
                        endif
                    endif
                    if (FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(2))==1) then
                        GenoStratIndex(i)=4
                        if (FinalSetter(ped%pedigree(i)%getMaternalGrandSireRecodedIndex())==1) then
                            GenoStratIndex(i)=2
                        endif
                    endif
                    if ((FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(2))==1).and.(FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(2))==1)) then
                        GenoStratIndex(i)=1
                    endif
                endif
            enddo

            TestMat=4
            CountCatTest=0
            AnisSummary=0.0

            do i=1,nAnisTest
                do j=1,inputParams%nsnp
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
                AnisSummary(i,4)=100*(float(count(ImputePhase(RecTestId(i),:,1)/=9))/inputParams%nsnp)
                AnisSummary(i,5)=100*(float(count(ImputePhase(RecTestId(i),:,2)/=9))/inputParams%nsnp)
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
                        do j=1,inputParams%nsnp
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
                        do j=1,inputParams%nsnp
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

            do i=1, ped%pedigreeSize-ped%nDummys
                if (ped%pedigree(i)%isDummy) then
                    exit
                endif
                write (45,'(a25,i3,2f7.2)') ped%pedigree(i)%originalID,FinalSetter(i),float(count(ImputePhase(i,:,1)/=9))/inputParams%nsnp &
                    ,float(count(ImputePhase(i,:,2)/=9))/inputParams%nsnp
            enddo

            do j=1,inputParams%nsnp
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
            allocate(TrueGenos(nAnisTest,inputParams%nSnpRaw))
            allocate(TrueGenosId(nAnisTest))
            allocate(RawGenos(nAnisTest,inputParams%nSnpRaw))
            allocate(TestMat(nAnisTest,inputParams%nSnpRaw))
            allocate(RecTestId(nAnisTest))
            allocate(AnisSummary(nAnisTest,5))
            allocate(Correlations(6))
            allocate(RealTestGenos(nAnisTest,inputParams%nSnpRaw))
            allocate(CalcCorPerAnimal(inputParams%nSnpRaw,2))
            allocate(CorrelationPerAnimal(nAnisTest))
            allocate(TmpVarPerGrp(nAnisTest))

            do i=1,nAnisTest
                read (35,*) TrueGenosId(i),TrueGenos(i,:)
            enddo

            RecTestId(:)=-99
            do i=1,nAnisTest
                do j=1,ped%pedigreeSize-ped%nDummys
                    ! print *, i,j, TrueGenosId(i), ped%pedigree(j)%originalID
                    if (trim(TrueGenosId(i))==trim((ped%pedigree(j)%originalID))) then
                        RecTestId(i)=j
                        TestAnimInformativeness(i,:)=GlobalTmpCountInf(j,1:6)
                        exit
                    endif
                enddo
            enddo
            if (count(RecTestId(:)==-99)>0) print*, "Error - There seems to be unidentifiablecount ",count(RecTestId(:)==-99)," individuals in the test file"

            do i=1,ped%nGenotyped
                read (36,*) dumC,WorkTmp(:)
                do j=1,nAnisTest
                    if (trim(TrueGenosId(j))==dumC) then
                        RawGenos(j,:)=WorkTmp(:)
                        exit
                    endif
                enddo
            enddo
            GenoStratIndex(:)=0
            do i=1,ped%pedigreeSize-ped%nDummys
                if (FinalSetter(i)/=1) then
                    GenoStratIndex(i)=6
                    if (FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(3))==1) then
                        GenoStratIndex(i)=5
                        if (FinalSetter(ped%pedigree(i)%getPaternalGrandSireRecodedIndex())==1) then
                            GenoStratIndex(i)=3
                        endif
                    endif
                    if (FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(2))==1) then
                        GenoStratIndex(i)=4
                        if (FinalSetter(ped%pedigree(i)%getMaternalGrandSireRecodedIndex())==1) then
                            GenoStratIndex(i)=2
                        endif
                    endif
                    if ((FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(2))==1).and.(FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(3))==1)) then
                        GenoStratIndex(i)=1
                    endif
                endif
            enddo
            TestMat=4
            CountCatTest=0
            AnisSummary=0.0
            do i=1,nAnisTest
                do j=1,inputParams%nSnpRaw
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
                AnisSummary(i,4)=100*(float(count(ImputePhase(RecTestId(i),:,1)/=9))/inputParams%nSnpRaw)
                AnisSummary(i,5)=100*(float(count(ImputePhase(RecTestId(i),:,2)/=9))/inputParams%nSnpRaw)
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
                        do j=1,inputParams%nSnpRaw
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
                        do j=1,inputParams%nSnpRaw
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

            do i=1, ped%pedigreeSize-ped%nDummys
                if (ped%pedigree(i)%isDummy) then
                    exit
                endif
                write (45,'(a25,i3,2f7.2)') ped%pedigree(i)%originalID,FinalSetter(i),float(count(ImputePhase(i,:,1)/=9))/inputParams%nSnpRaw &
                    ,float(count(ImputePhase(i,:,2)/=9))/inputParams%nSnpRaw
            enddo

            do j=1,inputParams%nsnp
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

        deallocate(Work)
        deallocate(WorkTmp)
        deallocate(GenoStratIndex)
    end subroutine FinalChecker


    !#############################################################################################################################################################################################################################

    subroutine Cleaner
        use Global

        use alphaimputeinmod
        implicit none

        type(AlphaImputeInput), pointer :: inputParams

        inputParams=>defaultInput
        ! call rmdir("GeneProb")
        ! call rmdir("IterateGeneProb")
        call system(RMDIR // " GeneProb")
        call system(RMDIR // " IterateGeneProb")

        ! if (inputParams%SexOpt==0) call system(" rm TempGeneProb.sh")
        ! if (inputParams%SexOpt==0) call system("rm TempIterateGeneProb.sh")
        if (inputParams%SexOpt==0) call system(RM // " TempGeneProb." // SH)
        if (inputParams%SexOpt==0) call system(RM // " TempIterateGeneProb." // SH)

    end subroutine Cleaner

    !#############################################################################################################################################################################################################################

    subroutine READINJUNK

        use Global
        use alphaimputeinmod

        type(AlphaImputeInput), pointer :: inputParams
        integer :: i,dum,j

        inputParams => defaultInput

        open (3001,file="fort.2008",status="old")
        do i=1,ped%pedigreeSize-ped%nDummys
            read (3001,*) dum,ImputePhase(i,:,1)
            read (3001,*) dum,ImputePhase(i,:,2)
        enddo

        do i=1,ped%pedigreeSize-ped%nDummys
            do j=1,inputParams%nsnp
                if ((ImputePhase(i,j,1)/=9).and.(ImputePhase(i,j,2)/=9)) then
                    if (ImputeGenos(i,j)==9) ImputeGenos(i,j)=ImputePhase(i,j,1)+ImputePhase(i,j,2)
                endif
            enddo
        enddo

    end subroutine READINJUNK

    !#############################################################################################################################################################################################################################

    ! SUBROUTINE ClusterIndivByChip(nClusters,ClusterMemberIndv,Centroid)
    SUBROUTINE ClusterIndivByChip(nClusters)
        use Global
        use alphaimputeinmod
        ! use GlobalClustering
        implicit none

        type(AlphaImputeInput), pointer :: inputParams
        integer, intent(IN) :: nClusters              ! Number of different SNP chips
        ! integer, intent(OUT) :: ClusterMemberIndv(:), Centroid(:)

        integer, allocatable :: res(:)             ! The output
        integer :: k                               ! The number of unique elements
        integer :: i, j

        ! integer :: nChips=3, SurrCounter, nClusters, dist, nIterations, NewCluster
        integer :: SurrCounter,dist,nIterations,NewCluster
        integer, allocatable :: ClusterMember(:), nSurrPerCluster(:)!, Centroid(:), ClusterMemberIndv(:)
        logical :: moved


        inputParams => defaultInput
        allocate(res(ped%nGenotyped))
        k = 1
        res(1) = nSnpsAnimal(1)

        outer: do i=2,size(nSnpsAnimal)
            do j=1,k
                if (res(j) == nSnpsAnimal(i)) then
                    ! Found a match so start looking again
                    cycle outer
                end if
            end do
            ! No match found so add it to the output
            k = k + 1
            res(k) = nSnpsAnimal(i)
        end do outer
        ! write(*,advance='no',fmt='(a,i0,a)') 'Unique list has ',k,' elements: '
        ! write(*,*) res(1:k)

        ! nClusters=nChips
        SurrCounter=k
        allocate(ClusterMember(SurrCounter))
        allocate(ClusterMemberIndv(ped%nGenotyped))
        allocate(centroid(nClusters))
        allocate(nSurrPerCluster(nClusters))

        ! First clusterization fixing centroids equally distant. Recalculation of centroids
        Centroid=0
        nSurrPerCluster=0

        do j=1,nClusters
            Centroid(j)=((j*inputParams%nsnp/nClusters)+((j-1)*inputParams%nsnp/nClusters))/2
        enddo

        nIterations=0
        moved=.TRUE.
        ClusterMember=0
        do while (moved==.TRUE. .and. nIterations<=100)
            nIterations=nIterations+1
            moved=.FALSE.

            do i=1,SurrCounter
                ! print *, ClusterMember(i)
                dist=inputParams%nsnp
                do j=1,nClusters
                    if(abs(res(i)-Centroid(j))<dist) then
                        dist=abs(res(i)-Centroid(j))
                        NewCluster=j
                    endif
                enddo
                if (NewCluster/=ClusterMember(i)) then
                    moved=.TRUE.
                    ClusterMember(i)=NewCluster
                endif
            enddo
            nSurrPerCluster=0
            Centroid=0
            do i=1,SurrCounter
                do j=1,nClusters
                    if (ClusterMember(i)==j) then
                        Centroid(j)=Centroid(j)+res(i)
                        nSurrPerCluster(j)=nSurrPerCluster(j)+1
                    endif
                enddo
            enddo
            do j=1,nClusters
                if (nSurrPerCluster(j)/=0) Centroid(j)=Centroid(j)/nSurrPerCluster(j)
            enddo
        enddo

        ! Assign individuals to clusters
        NewCluster=0
        dist=0
        print *, size(ClusterMemberIndv)
        ClusterMemberIndv=0

        ! open (unit=2222,file='nSnpsAnimalCluster.txt',status='unknown')
        do i=1,ped%nGenotyped
            dist=inputParams%nsnp
            do j=1,nClusters
                if(abs(nSnpsAnimal(i)-Centroid(j))<dist) then
                    dist=abs(nSnpsAnimal(i)-Centroid(j))
                    NewCluster=j
                endif
            enddo
            if (NewCluster/=ClusterMemberIndv(i)) then
                ClusterMemberIndv(i)=NewCluster
            endif
            ! write(2222,*) i, ClusterMemberIndv(i), Centroid(ClusterMemberIndv(i))
        enddo
        ! close(2222)

    end SUBROUTINE ClusterIndivByChip

    !#############################################################################################################################################################################################################################
    SUBROUTINE SnpCallRate()
        use Global

        use alphaimputeinmod
        use ISO_Fortran_Env

        implicit none

        integer :: i, CountMiss, UOutputs
        type(AlphaImputeInput), pointer :: inputParams

        inputParams => defaultInput


        open(newunit=UOutputs, file="." // DASH // "Miscellaneous" // DASH // "SnpCallRateByAnimal.txt",status='unknown')

        do i=1,ped%nGenotyped
            countMiss = ped%pedigree(ped%genotypeMap(i))%individualGenotype%numMissing()
            write(UOutputs,'(a20,6f5.1)') ped%pedigree(ped%genotypeMap(i))%originalID, (inputParams%nsnp-CountMiss)*100/real(inputParams%nsnp)
        end do
        close(UOutputs)

    END SUBROUTINE SnpCallRate

    !#############################################################################################################################################################################################################################
    SUBROUTINE CheckImputationInconsistencies(ImpGenos, ImpPhase, n, m)
        implicit none

        integer, intent(in) :: n, m
        integer(kind=1), dimension (:,:), intent(inout) :: ImpGenos
        integer(kind=1), dimension (:,:,:), intent(inout) :: ImpPhase

        integer :: i, j, k

        do j = 1, m
            do i = 1, n
                do k = 1, 2
                    if (ImpPhase(i, j, k) < 0) then
                        ImpPhase(i, j, k) = 9
                        ! ImpGenos(i, j) = 9
                    end if
                    if (ImpPhase(i, j, k) > 1) then
                        ImpPhase(i, j, k) = 9
                        ! ImpGenos(i, j) = 9
                    end if
                enddo
                if (ImpGenos(i, j) < 0) then
                    ImpGenos(i, j) = 9
                end if
                if (ImpGenos(i, j) > 2) then
                    ImpGenos(i, j) = 9
                end if
            enddo
        enddo

    END SUBROUTINE CheckImputationInconsistencies

    !#############################################################################################################################################################################################################################

    subroutine Titles

        call PrintVersion
        print *, ""
        print *, ""
        print *, ""

    end subroutine Titles

    !#############################################################################################################################################################################################################################

    subroutine Header

        print *, ""
        print *, "                              ***********************                         "
        print *, "                              *                     *                         "
        print *, "                              *     AlphaImpute     *                         "
        print *, "                              *                     *                         "
        print *, "                              ***********************                         "
        print *, "                                                                              "
        print *, "                    Software For Phasing and Imputing Genotypes               "

    end subroutine Header

    !#############################################################################################################################################################################################################################

    subroutine PrintVersion

        call Header
        print *, ""
        print *, "                              Commit:   "//TOSTRING(COMMIT),"                     "
        print *, "                              Compiled: "//__DATE__//", "//__TIME__
        print *, ""

    end subroutine PrintVersion

    !#############################################################################################################################################################################################################################

    subroutine PrintTimerTitles
        use Global
        use iso_fortran_env
        implicit none

        real(kind=real64) :: etime          ! Declare the type of etime()
        real(kind=real64) :: elapsed(2)     ! For receiving user and system time
        real(kind=real64) :: total,Minutes,Hours,Seconds

        print *, ""
        print *, ""
        call Header
        print*, ""
        print*, "                                  No Liability"
        print*, ""
        print*, "                Analysis Finished                         "

        total=etime(elapsed)
        Minutes=total/60
        Seconds=Total-(INT(Minutes)*60)
        Hours=Minutes/60
        Minutes=INT(Minutes)-(INT(Hours)*60)
        print '(A107,A7,I3,A9,I3,A9,F6.2)', "Time Elapsed","Hours", INT(Hours),"Minutes",INT(Minutes),"Seconds",Seconds

        open (unit=32,file="." // DASH // "Miscellaneous" // DASH // "Timer.txt",status="unknown")

        write(32,'(A27,A7,I3,A9,I3,A9,F6.2)') "Time Elapsed","Hours", INT(Hours),"Minutes",INT(Minutes),"Seconds",Seconds

    end subroutine PrintTimerTitles

end module AlphaImputeModule




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

    integer :: markers
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

        call MakeFiles

    else

        call MakeDirectories(RUN_HMM_NGS)
        call ReadInData
        call SnpCallRate
        allocate(Reads(ped%nGenotyped,inputParams%nsnp))
        allocate(ImputeGenos(0:ped%nGenotyped,inputParams%nsnp))
        allocate(ImputePhase(0:ped%nGenotyped,inputParams%nsnp,2))
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
                        
                        call runGeneProbAlphaImpute(1, inputParams%nsnp, ped, GenosProbs, MAF)
                        call WriteProbabilities("./Results/GenotypeProbabilities.txt", GenosProbs, ped,ped%pedigreeSize-ped%nDummys, inputParams%nsnp)
                    endif

                    ! deallocate(GenosProbs)

                    if (inputParams%restartOption==OPT_RESTART_GENEPROB) then
                        write(6,*) "Restart option 1 stops program after Geneprobs jobs have finished"
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
                    allocate(GenosProbs(ped%nDummys+ ped%pedigreeSize-ped%nDummys, markers, 2))
                    call readInGeneProbData(GenosProbs)
                    ! TODO may not need
                    call WriteProbabilities("./Results/GenotypeProbabilities.txt", GenosProbs, ped,ped%pedigreeSize-ped%nDummys, inputParams%nsnp)
                    deallocate(GenosProbs)
                    write(6,*) "Restart option 1 stops program after genotype probabilities have been outputted"
                    stop
                case default
                    print *, "ERROR: BYPAS GENEPROB SET INCORRECTLY"
                    stop 1
            end select
        endif



    if (inputParams%managephaseon1off0==1) then



        if (inputParams%restartOption<OPT_RESTART_IMPUTATION) Then
            call PhasingManagementNew(APResults)

        endif

        if (inputParams%restartOption==OPT_RESTART_PHASING) then
        ! TODO need to write out phasing results
            write(6,*) "Restart option 2 stops program after Phasing has been managed"
            stop
        endif
    endif


endif

if (inputParams%hmmoption/=RUN_HMM_NGS) then
    ! If we only want to phase data, then skip all the imputation steps
    if (inputParams%PhaseTheDataOnly==0) Then
        call ImputationManagement

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
            call FinalChecker

        endif
    endif
endif

call PrintTimerTitles
! if (inputParams%restartOption > OPT_RESTART_IMPUTATION) then
!     call system(RM // " Tmp2345678.txt")
! end if

end program AlphaImpute

