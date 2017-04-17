#ifdef OS_UNIX

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

#else

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
#endif

module Output
    ! use global
    use iso_fortran_env
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
        integer, intent(IN) :: Indexes(:)
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
        double precision, intent(IN) :: GenosProbs(:,:,:)
        double precision, allocatable :: GenosProbsTmp(:,:,:)
        ! Local Variable
        integer :: i!,j,k, n0, n1, n2
        ! real, allocatable :: Probs0(:), Probs1(:)

        allocate(GenosProbsTmp(nAnims,nsnps,2))
        GenosProbsTmp(:,:,1) = GenosProbs(:,:, 1)
        GenosProbsTmp(:,:,2) =  GenosProbs(:,:,2) +  GenosProbs(:,:,3)
        open (unit=55,file=outFile,status="unknown")

        ! allocate(Probs0(nSnps))
        ! allocate(Probs1(nSnps))

        do i=1,nAnims
            write (55,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%Pedigree(i)%originalID,GenosProbsTmp(i,:,1)
            write (55,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%Pedigree(i)%originalID,GenosProbsTmp(i,:,2)
            ! enddo
        enddo
    end subroutine WriteProbabilitiesGeneProb


    subroutine WriteProbabilitiesFull(outFile, GenosProbs, ped,nAnims)
        use PedigreeModule
        character(len=*), intent(IN) :: outFile
        integer, intent(IN) :: nAnims
        type(pedigreeHolder), intent(IN) :: ped
        double precision, intent(IN) :: GenosProbs(:,:,:)
        ! Local Variable
        integer :: i!,j,k, n0, n1, n2
        ! real, allocatable :: Probs0(:), Probs1(:)
        open(unit=55,file=outFile,status="unknown")

        do i=1,nAnims
            write (55,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%Pedigree(i)%originalID,GenosProbs(i,:,1)
            write (55,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%Pedigree(i)%originalID,GenosProbs(i,:,2)
            write (55,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%Pedigree(i)%originalID,GenosProbs(i,:,3)
            write (55,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%Pedigree(i)%originalID,GenosProbs(i,:,4)
            ! enddo
        enddo
    end subroutine WriteProbabilitiesFull

    subroutine readProbabilitiesGeneProb(file, GenosProbs,nAnims)
        use constantModule
        character(len=*), intent(IN) :: file
        integer, intent(IN) :: nAnims
        double precision, intent(out) :: GenosProbs(:,:,:)
        integer :: fileUnit        
        integer :: i
        character(len=IDLENGTH) :: dum

        open (newunit=fileUnit,file=file,status="unknown")

        do i=1,nAnims
            read(fileUnit,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') dum,GenosProbs(i,:,1)
            read(fileUnit,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') dum,GenosProbs(i,:,2)
        enddo
    end subroutine readProbabilitiesGeneProb


        subroutine readProbabilitiesGeneProbCluster(GenosProbs,nAnims, inputParams, GpIndex)
        use constantModule
        use AlphaImputeInMod
        type(AlphaImputeInput), intent(in) :: inputParams
        integer, intent(IN) :: nAnims
        integer, intent(in) :: GpIndex(:,:)
        double precision, intent(out) :: GenosProbs(:,:,:)
        integer :: fileUnit        
        integer :: i,StSnp,EnSnp,h
        character(len=300) ::filout
        character(len=IDLENGTH) :: dum

        do h=1,inputParams%nProcessors
#ifndef _WIN32
            write (filout,'("./GeneProb/GeneProb"i0,"/GeneProbs.txt")')h            
#else
            write (filout,'(".\GeneProb\GeneProb"i0,"\GeneProbs.txt")')h            
#endif
            open (newUnit=fileUnit,file=trim(filout),status="unknown")
            StSnp=GpIndex(h,1)
            EnSnp=GpIndex(h,2)

            do i=1,nAnims
                read(fileUnit,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') dum,GenosProbs(i,StSnp:EnSnp,1)
                read(fileUnit,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') dum,GenosProbs(i,StSnp:EnSnp,2)
            enddo
        enddo
    end subroutine readProbabilitiesGeneProbCluster

    subroutine readProbabilitiesFull(file, GenosProbs,nAnims, nSnps)
        use constantModule
        character(len=*), intent(IN) :: file
        integer, intent(IN) :: nSnps,nAnims
        real(kind=real64),allocatable, intent(out) :: GenosProbs(:,:,:)
        integer :: fileUnit        
        integer :: i
        character(len=IDLENGTH) :: dum

        open (newunit=fileUnit,file=file,status="unknown")

        ! allocate(Probs0(nSnps))
        ! allocate(Probs1(nSnps))
        if (allocated(GenosProbs)) then
            deallocate(Genosprobs)
        endif
        allocate(GenosProbs(nAnims,nSnps, 4))

        do i=1,nAnims
            read(fileUnit,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') dum,GenosProbs(i,:,1)
            read(fileUnit,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') dum,GenosProbs(i,:,2)
            read(fileUnit,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') dum,GenosProbs(i,:,3)
            read(fileUnit,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') dum,GenosProbs(i,:,4)
            ! enddo
        enddo
    end subroutine readProbabilitiesFull

        subroutine readProbabilitiesFullCluster(GenosProbs,nAnims, nSnps, inputParams, GpIndex)
        use constantModule
        use AlphaImputeInMod
        type(AlphaImputeInput), intent(in) :: inputParams
        integer, intent(IN) :: nSnps,nAnims
        integer, intent(in) :: GpIndex(:,:)
        real(kind=real64),allocatable, intent(out) :: GenosProbs(:,:,:)
        integer :: fileUnit        
        integer :: i,h,StSnp, EnSnp
        character(len=IDLENGTH) :: dum
        character(len=300) :: filout
        do h=1,inputParams%nProcessors
#ifndef _WIN32
            write (filout,'("./GeneProb/GeneProb"i0,"/GeneProbs.txt")')h            
#else
            write (filout,'(".\GeneProb\GeneProb"i0,"\GeneProbs.txt")')h            
#endif
            open (newUnit=fileUnit,file=trim(filout),status="unknown")
            StSnp=GpIndex(h,1)
            EnSnp=GpIndex(h,2)



            ! allocate(Probs0(nSnps))
            ! allocate(Probs1(nSnps))
            if (allocated(GenosProbs)) then
                deallocate(Genosprobs)
            endif
            allocate(GenosProbs(nAnims,nSnps, 4))
            
            do i=1,nAnims
              read(fileUnit,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') dum,GenosProbs(i,StSnp:EnSnp,1)
                read(fileUnit,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') dum,GenosProbs(i,StSnp:EnSnp,2)
                read(fileUnit,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') dum,GenosProbs(i,StSnp:EnSnp,3)
                read(fileUnit,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') dum,GenosProbs(i,StSnp:EnSnp,4)
                ! enddo
            enddo
            close(fileUnit)
        enddo
    end subroutine readProbabilitiesFullCluster




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

    subroutine readGeneProbsCluster(GlobalWorkPhase,ped,GpIndex, inputParams,GeneProbThresh)
        use PedigreeModule
        use iso_fortran_env
        use AlphaImputeInMod
        use constantModule
        type(PedigreeHolder) :: ped
        integer(kind=1), intent(out) :: GlobalWorkPhase(:,:,:)
        integer, intent(in) :: GpIndex(:,:)
        real, intent(in) :: GeneProbThresh
        type(AlphaImputeInput), pointer,intent(in) :: inputParams
        character(len=300) :: filout
        real(real64), allocatable, dimension(:,:) :: PatAlleleProb,MatAlleleProb,GeneProbWork
        integer :: unit,h,StSnp,EnSnp,i,j
        character(len=IDLENGTH) :: dum
        

        allocate(GeneProbWork(inputParams%nsnp,4))
        allocate(PatAlleleProb(inputParams%nsnp,2))
        allocate(MatAlleleProb(inputParams%nsnp,2))
        do h=1,inputParams%nProcessors
#ifndef _WIN32
            write (filout,'("./GeneProb/GeneProb"i0,"/GeneProbs.txt")')h            
#else
            write (filout,'(".\GeneProb\GeneProb"i0,"\GeneProbs.txt")')h            
#endif
            open (newUnit=unit,file=trim(filout),status="unknown")
            StSnp=GpIndex(h,1)
            EnSnp=GpIndex(h,2)
            do i=1,ped%pedigreeSize-ped%nDummys
        
                do j=1,4
                    read (unit,*) dum,GeneProbWork(StSnp:EnSnp,j)
                enddo
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
            close(unit)
        enddo

        deallocate(PatAlleleProb)
        deallocate(MatAlleleProb)

        



    end subroutine readGeneProbsCluster

    subroutine ReReadGeneProbs(GlobalWorkPhase,ped,path, nsnp, GeneProbThresh) !TODO this needs rewritten
        ! Read genotype probabilities from files and phase allele based in these probabilities.
        ! This files should have been already created during previous calls to AlphaImpute (inputParams%restartOption<3)
        ! Phasing information is store in the variable GlobalWorkPhase
        use AlphaImputeInMod
        use PedigreeModule
        implicit none

        type(PedigreeHolder),intent(in) :: ped
        real, intent(in) :: GeneProbThresh
        integer(kind=1), intent(out) :: GlobalWorkPhase(:,:,:)
        integer, intent(in) :: nsnp
        integer :: i,j,dum,fileUnit
        real, allocatable, dimension(:,:) :: PatAlleleProb,MatAlleleProb,GeneProbWork
        character(len=*), intent(in) :: path
        type(AlphaImputeInput), pointer :: inputParams



        inputParams => defaultInput

        allocate(PatAlleleProb(inputParams%nsnp,2))
        allocate(MatAlleleProb(inputParams%nsnp,2))
        allocate(GeneProbWork(inputParams%nsnp,4))
        GlobalWorkPhase=9

            ! TODOgeneprob info read here 
        open (newunit=fileUnit,file=path,status="unknown")
        do i=1,ped%pedigreeSize-ped%nDummys                                           ! The number of lines of GeneProbs.txt files is = nAnisP x 4
            do j=1,2                                            ! where 4 stands for the two paternal and the two maternal haplotypes
                read (fileUnit,*) dum,GeneProbWork(1:nsnp,j)
            enddo

            ! GeneProbWork(:,1) == Probability 0-0 = Prob00
            ! GeneProbWork(:,2) == Probability 0-1 = Prob01
            ! GeneProbWork(:,3) == Probability 1-0 = Prob10
            ! GeneProbWork(:,4) == Probability 1-1 = Prob11
            PatAlleleProb(:,1)=GeneProbWork(:,1)+GeneProbWork(:,2)    ! PatAlleleProb(:,1) == Probability Paternal allele is 0 = Prob00 + Prob01
            PatAlleleProb(:,2)=GeneProbWork(:,3)+GeneProbWork(:,4)    ! PatAlleleProb(:,2) == Probability Paternal allele is 1 = Prob10 + Prob11
            MatAlleleProb(:,1)=GeneProbWork(:,1)+GeneProbWork(:,3)    ! PatAlleleProb(:,3) == Probability Maternal allele is 0 = Prob00 + Prob10
            MatAlleleProb(:,2)=GeneProbWork(:,2)+GeneProbWork(:,4)    ! PatAlleleProb(:,4) == Probability Maternal allele is 1 = Prob01 + Prob11

            do j=1,nsnp
                if (PatAlleleProb(j,1)>=GeneProbThresh) GlobalWorkPhase(i,j,1)=0
                if (PatAlleleProb(j,2)>=GeneProbThresh) GlobalWorkPhase(i,j,1)=1
                if (MatAlleleProb(j,1)>=GeneProbThresh) GlobalWorkPhase(i,j,2)=0
                if (MatAlleleProb(j,2)>=GeneProbThresh) GlobalWorkPhase(i,j,2)=1
            enddo
        enddo
        close(110)
        
        deallocate(PatAlleleProb)
        deallocate(MatAlleleProb)
        deallocate(GeneProbWork)
    end subroutine ReReadGeneProbs


end module Output
