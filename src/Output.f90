module Output
    ! use global

    implicit none


    INTERFACE WriteProbabilities
        module procedure :: WriteProbabilitiesHMM
        module procedure :: WriteProbabilitiesGeneProb

    end INTERFACE WriteProbabilities
contains


    ! subroutine WriteProbabilitiesHMM(outFile, Indexes, Ids, nAnims, nSnps)
    subroutine WriteProbabilitiesHMM(outFile, Indexes, nAnims, nSnps)
        use global
        use GlobalVariablesHmmMaCH
        use alphaimputeinmod
        character(len=*), intent(IN) :: outFile
        integer, intent(IN) :: nAnims, nSnps
        integer, intent(IN) :: Indexes(nAnims)
        type(AlphaImputeInput), pointer :: inputParams
        ! character*(20), intent(IN) :: Ids(:)
        ! Local variables
        integer :: i,j, n0, n1, n2
        real :: d
        real, allocatable :: Probs0(:), Probs1(:)

        inputParams => defaultInput
        open (unit=55,file=outFile,status="unknown")

        allocate(Probs0(nSnps))
        allocate(Probs1(nSnps))

        n0=0
        n1=0
        n2=0
        d=0.0

        do i=1,nAnims
            do j=1,nSnps
                n1 = GenosCounts(i,j,1)                           ! Heterozygous
                n2 = GenosCounts(i,j,2)                           ! Homozygous: 2 case
                n0 = (GlobalRoundHmm-inputParams%HmmBurnInRound) - n1 - n2     ! Homozygous: 0 case
                d = n0 + n1 + n2
                Probs0(j)=n0/d
                Probs1(j)=n1/d
            enddo
            write (55,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%Pedigree(indexes(i))%originalID,Probs0(:)
            write (55,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%Pedigree(indexes(i))%originalID,Probs1(:)
            n0=0
            n1=0
            n2=0
            d=0.0
        enddo

        close(55)
        ! write(0,*) "Wrote out file", outFile, "with genotype probabilities"
    end subroutine WriteProbabilitiesHMM

    subroutine WriteProbabilitiesGeneProb(outFile, GenosProbs, ped,nAnims, nSnps)
        use PedigreeModule
        character(len=*), intent(IN) :: outFile
        integer, intent(IN) :: nSnps,nAnims
        type(pedigreeHolder), intent(IN) :: ped
        double precision, intent(IN) :: GenosProbs(ped%pedigreesize-ped%nDummys,nSnps,2)

        ! Local Variable
        integer :: i!,j,k, n0, n1, n2
        ! real, allocatable :: Probs0(:), Probs1(:)

        open (unit=55,file=outFile,status="unknown")

        ! allocate(Probs0(nSnps))
        ! allocate(Probs1(nSnps))

        do i=1,nAnims
            write (55,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%Pedigree(i)%originalID,GenosProbs(i,:,1)
            write (55,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%Pedigree(i)%originalID,GenosProbs(i,:,2)
            ! enddo
        enddo
    end subroutine WriteProbabilitiesGeneProb



    subroutine ReReadIterateGeneProbs(GenosProbs, IterGeneProb, nAnis)
        ! Read genotype probabilities from files and phase allele based in these probabilities.
        ! This files should have been already created during previous calls to AlphaImpute (RestartOption<3)
        ! The subroutine outputs the genotype probabilities of the homozygous genotype of the reference allele,
        ! G00, and the heterozygous genotype, Gh = G10 + G01. The homozygous genotype for the alternative allele can be inferred
        ! from the these two as G11 = 1 - G00 - Gh
        use Global
        use alphaimputeinmod
        implicit none

        logical, intent(IN) :: IterGeneProb
        integer, intent(IN) :: nAnis
        !double precision, dimension(:,:,:), intent(INOUT) :: GenosProbs(nAnis,markers,2)
        double precision, dimension(:,:,:), intent(INOUT) :: GenosProbs
        type(AlphaImputeInput), pointer :: inputParams
        ! Local variables
        integer :: h,i,j,dum,StSnp,EnSnp
        double precision, allocatable :: GeneProbWork(:,:)
        character(len=300) :: inFile

        inputParams => defaultInput

        do h=1,inputParams%nProcessors
            if (IterGeneProb) then
                write (inFile,'("IterateGeneProb/GeneProb"i0,"/GeneProbs.txt")')h          !here
            else
                write (inFile,'("GeneProb/GeneProb"i0,"/GeneProbs.txt")')h          !here
            end if

            open (unit=110,file=trim(inFile),status="unknown")
            StSnp=GpIndex(h,1)          ! Where SNPs start
            EnSnp=GpIndex(h,2)          ! Where SNPs end
            allocate(GeneProbWork(EnSnp-StSnp+1,4))
            GeneProbWork=9
            do i=1,nAnis                                           ! The number of lines of GeneProbs.txt files is = nAnisP x 4
                do j=1,4                                            ! where 4 stands for the two paternal and the two maternal haplotypes
                    read (110,*) dum,GeneProbWork(:,j)
                enddo

                GenosProbs(i,StSnp:EnSnp,1) = GeneProbWork(:,1)
                GenosProbs(i,StSnp:EnSnp,2) = GeneProbWork(:,2) + GeneProbWork(:,3)
            enddo
            deallocate(GeneProbWork)
            close(110)
        enddo

    end subroutine ReReadIterateGeneProbs



    subroutine ReadInPrePhasedData
        ! Impute phase information from pre-phased file. Count the number of pre-phased individuals
        use Global

        use AlphaImputeInMod

        integer :: h,j,k,nAnisPrePhased,CountPrePhased,tmpID
        integer, allocatable,dimension(:,:) :: WorkPhase
        character(len=300) :: dumC
        type(AlphaImputeInput), pointer :: inputParams

        inputParams => defaultInput
        ! Count animals prephased
        nAnisPrePhased=0
        allocate(WorkPhase(inputParams%nSnpRaw,2))
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

            tmpID = ped%dictionary%getValue(trim(dumC))
            if (tmpID /= DICT_NULL) then   ! Check if any animal in the file agrees with the animals in the pedigree
                h=0
                do j=1,inputParams%nSnpRaw
                    if (SnpIncluded(j)==1) then ! Check if this SNP has to be considered (may be it has been removed during the edition step)
                        h=h+1
                        ! Impute phase only if this locus is phased (in the PrePhased file)
                        if ((WorkPhase(j,1)==0).or.(WorkPhase(j,1)==1)) ImputePhase(tmpID,h,1)=WorkPhase(j,1)
                        if ((WorkPhase(j,2)==0).or.(WorkPhase(j,2)==1)) ImputePhase(tmpID,h,2)=WorkPhase(j,2)
                    endif
                enddo
                CountPrePhased=CountPrePhased+1
            else

            endif
        enddo

        print*, " "
        print*, " ",CountPrePhased," valid pre-phased indiviudals in the user specified pre-phased file"

    end subroutine ReadInPrePhasedData


    !######################################################################################################################################################################################

    subroutine ReReadGeneProbs
        ! Read genotype probabilities from files and phase allele based in these probabilities.
        ! This files should have been already created during previous calls to AlphaImpute (inputParams%restartOption<3)
        ! Phasing information is store in the variable GlobalWorkPhase
        use Global
        use AlphaImputeInMod
        implicit none


        integer :: h,i,j,dum,StSnp,EnSnp
        real, allocatable, dimension(:,:) :: PatAlleleProb,MatAlleleProb,GeneProbWork
        character(len=300) :: filout
        type(AlphaImputeInput), pointer :: inputParams



        inputParams => defaultInput

        allocate(PatAlleleProb(inputParams%nsnp,2))
        allocate(MatAlleleProb(inputParams%nsnp,2))
        allocate(GeneProbWork(inputParams%nsnp,4))
        GlobalWorkPhase=9
        do h=1,inputParams%nProcessors
#ifndef _WIN32
            write (filout,'("GeneProb/GeneProb"i0,"/GeneProbs.txt")')h          !here
#else
            write (filout,'("GeneProb\GeneProb"i0,"\GeneProbs.txt")')h          !here
#endif
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

        deallocate(PatAlleleProb)
        deallocate(MatAlleleProb)
        deallocate(GeneProbWork)
    end subroutine ReReadGeneProbs



end module Output
