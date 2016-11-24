module Output
! use global

implicit none

INTERFACE
  SUBROUTINE ReReadIterateGeneProbs(GenosProbs, IterGeneProb, nAnis)
    use Global
    logical, intent(IN) :: IterGeneProb
    integer, intent(IN) :: nAnis
!    double precision, intent(OUT) :: GenosProbs(nAnisP,markers,2)
    double precision, dimension(:,:,:), intent(INOUT) :: GenosProbs
  END SUBROUTINE ReReadIterateGeneProbs
END INTERFACE
end module Output



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

subroutine WriteProbabilitiesGeneProb(outFile, GenosProbs, ped, nExtraAnims, nAnisP, nSnps)
use PedigreeModule
character(len=*), intent(IN) :: outFile
integer, intent(IN) :: nExtraAnims, nAnisP, nSnps
double precision, intent(IN) :: GenosProbs(nAnisP,nSnps,2)
type(pedigreeHolder), intent(IN) :: ped

! Local Variable
integer :: i,k!,j,k, n0, n1, n2
! real, allocatable :: Probs0(:), Probs1(:)

open (unit=55,file=outFile,status="unknown")

! allocate(Probs0(nSnps))
! allocate(Probs1(nSnps))

do i=nExtraAnims+1,nAnisP
    k=i-nExtraAnims
    write (55,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%Pedigree(i)%originalID,GenosProbs(k,:,1)
    write (55,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%Pedigree(i)%originalID,GenosProbs(k,:,2)
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

