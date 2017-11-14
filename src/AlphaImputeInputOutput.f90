#ifndef _WIN32

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

module AlphaImputeInputOutputModule
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
        use AlphaImputeSpecFileModule
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




    subroutine ReadInPrePhasedData
        ! Impute phase information from pre-phased file. Count the number of pre-phased individuals
        use Global

        use AlphaImputeSpecFileModule

        integer :: h,j,k,nAnisPrePhased,CountPrePhased,tmpID
        integer(kind=1), allocatable,dimension(:,:) :: WorkPhase
        character(len=300) :: dumC
        type(AlphaImputeInput), pointer :: inputParams

        inputParams => defaultInput
        ! Count animals prephased
        nAnisPrePhased=0
        allocate(WorkPhase(inputParams%nSnpRaw,2))
        do
            read (inputParams%prePhasedFileUnit,*,iostat=k) dumC
            nAnisPrePhased=nAnisPrePhased+1
            if (k/=0) then
                nAnisPrePhased=nAnisPrePhased-1
                exit
            endif
        enddo
        rewind(inputParams%prePhasedFileUnit)
        nAnisPrePhased=nAnisPrePhased/2         ! Two haplotypes per animal have been read

        CountPrePhased=0
        do k=1,nAnisPrePhased
            read (inputParams%prePhasedFileUnit,*) dumC,WorkPhase(:,1)     ! Paternal haplotype
            read (inputParams%prePhasedFileUnit,*) dumC,WorkPhase(:,2)     ! Maternal haplotype

            tmpID = ped%dictionary%getValue(trim(dumC))
            if (tmpID /= DICT_NULL) then   ! Check if any animal in the file agrees with the animals in the pedigree
                h=0
                do j=1,inputParams%nSnpRaw
                    if (SnpIncluded(j)/=0) then ! Check if this SNP has to be considered (may be it has been removed during the edition step)
                        h=h+1
                        ! Impute phase only if this locus is phased (in the PrePhased file)
                        if ((WorkPhase(j,1)==0).or.(WorkPhase(j,1)==1)) then
                            call ped%pedigree(tmpId)%individualPhase(1)%setPhase(h,WorkPhase(j,1))
                        else if ((WorkPhase(j,2)==0).or.(WorkPhase(j,2)==1)) then
                            call ped%pedigree(tmpId)%individualPhase(2)%setPhase(h,WorkPhase(j,2))
                        endif
                    endif
                enddo
                CountPrePhased=CountPrePhased+1
            else

            endif
        enddo
        rewind(inputParams%prePhasedFileUnit)

        print*, " "
        print*, " ",CountPrePhased," valid pre-phased indiviudals in the user specified pre-phased file"

    end subroutine ReadInPrePhasedData



        !#############################################################################################################################################################################################################################

    subroutine MakeDirectories(HMM)
        use global
        use AlphaImputeSpecFileModule

        implicit none

        integer, intent(in) :: HMM
        integer :: i
        character(len=300) :: FolderName
        type(AlphaImputeInput), pointer :: inputParams

        inputParams => defaultInput

        if (HMM == RUN_HMM_NGS) then
            ! call rmdir("Miscellaneous")
            call system(RMDIR // " " // trim(inputparams%resultFolderPath))
            call system(RMDIR // " Miscellaneous")
            ! call system("mkdir Miscellaneous")
            call system(MD // " " // trim(inputParams%resultFolderPath))
            call system(MD // " Miscellaneous")

        else
            print*, ""
            call system(RMDIR // " Miscellaneous")
            call system(RMDIR // " Phasing")
            call system(RMDIR // " " // trim(inputParams%resultFolderPath))

            call system(MD // " Phasing")
            call system(MD // " Miscellaneous")
            call system(MD // " " // trim(inputParams%resultFolderPath))

            do i=1,inputParams%nPhaseInternal
                write (FolderName,'("Phase"i0)')i
                ! call system ("mkdir Phasing/" // FolderName)
                call system(MD // " Phasing" // DASH // FolderName)
            enddo
        endif


    end subroutine MakeDirectories



    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Read files
    !
    !> @details    This subroutine reads the pedigree data and the gender file
    !>             in case of Sex Chromosome, as well as genotype information
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date       Nov 10, 2016
    !
    ! PARAMETERS:
    !> @param[in]
    !---------------------------------------------------------------------------
    subroutine ReadInData

        use Global
        use AlphaImputeSpecFileModule
        implicit none

        type(AlphaImputeInput), pointer :: inputParams
        inputParams=> defaultInput
        
        ! Read the pedigree information

        if (trim(inputParams%pedigreefile) /= "NoPedigree") then
            if (inputParams%SexOpt==1) then
                call initPedigree(ped,inputParams%pedigreefile,genderfile=inputParams%genderFile, nsnps=inputParams%nsnp)
            else 
                call initPedigree(ped,inputParams%pedigreefile, nsnps=inputParams%nsnp)
            endif

            if (inputParams%hmmoption /= RUN_HMM_NGS) then
                call ped%addGenotypeInformationFromFile(inputParams%GenotypeFile,inputParams%nsnp,initAll=1)
            endif
        else

           if (inputParams%hmmoption /= RUN_HMM_NGS) then
                ! init pedigree from genotype file
                call initPedigreeGenotypeFiles(ped,inputParams%GenotypeFile, nsnp=inputParams%nsnp,initAll=1)
            endif
        endif

    end subroutine ReadInData


end module AlphaImputeInputOutputModule
