module Output
! use global

implicit none

! contains
end module Output

subroutine WriteProbabilitiesHMM(outFile, nExtraAnims, Ids, nAnisP, nSnps)
use GlobalVariablesHmmMaCH

character(len=*), intent(IN) :: outFile
integer, intent(IN) :: nExtraAnims, nAnisP, nSnps
character*(20), intent(IN) :: Ids(nAnisP)

! Local variables
integer :: i,j,k, n0, n1, n2
real :: d
real, allocatable :: Probs0(:), Probs1(:)

open (unit=55,file=outFile,status="unknown")

allocate(Probs0(nSnps))
allocate(Probs1(nSnps))

n0=0
n1=0
n2=0
d=0.0

do i=nExtraAnims+1,nAnisP
    k=i-nExtraAnims
    do j=1,nSnps
        ! call GetCounts(GenosCounts, n0, n1, n2, nAnisP, nSnps)
        n1 = GenosCounts(k,j,1)                           ! Heterozygous
        n2 = GenosCounts(k,j,2)                           ! Homozygous: 2 case
        n0 = (GlobalRoundHmm-HmmBurnInRound) - n1 - n2     ! Homozygous: 0 case
        d = n0 + n1 + n2
        Probs0(j)=n0/d
        Probs1(j)=n1/d
    enddo
    write (55,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') Ids(i),Probs0(:)
    write (55,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') Ids(i),Probs1(:)
    n0=0
    n1=0
    n2=0
    d=0.0
enddo

close(55)
! write(0,*) "Wrote out file", outFile, "with genotype probabilities"
end subroutine WriteProbabilitiesHMM

subroutine WriteProbabilitiesGeneProb(outFile, GenosProbs, Ids, nExtraAnims, nAnisP, nSnps)

character(len=*), intent(IN) :: outFile
integer, intent(IN) :: nExtraAnims, nAnisP, nSnps
double precision, intent(IN) :: GenosProbs(nAnisP,nSnps,4)
character*(20), intent(IN) :: Ids(nAnisP)

! Local Variable
integer :: i,k!,j,k, n0, n1, n2
! real, allocatable :: Probs0(:), Probs1(:)

open (unit=55,file=outFile,status="unknown")

! allocate(Probs0(nSnps))
! allocate(Probs1(nSnps))

do i=nExtraAnims+1,nAnisP
    k=i-nExtraAnims
    ! do j=1,nSnps
    ! Probs0=GenosProbs(k,:,1)
    ! Probs1=GenosProbs(k,:,2)
    write (55,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') Ids(i),GenosProbs(k,:,1)
    write (55,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') Ids(i),GenosProbs(k,:,2)
    ! enddo
enddo
end subroutine WriteProbabilitiesGeneProb



subroutine ReReadIterateGeneProbs(GenosProbs)
! Read genotype probabilities from files and phase allele based in these probabilities. 
! This files should have been already created during previous calls to AlphaImpute (RestartOption<3)
! Phasing information is store in the variable GlobalWorkPhase
use Global

implicit none

double precision, intent(OUT) :: GenosProbs(nAnisP,nSnp,2)

! Local variables
integer :: h,i,j,dum,StSnp,EnSnp
double precision :: PatAlleleProb(nSnp,2),MatAlleleProb(nSnp,2),HetProb(nSnp),GeneProbWork(nSnp,4)
character(len=300) :: inFile

GeneProbWork=9
do h=1,nProcessors
    write (inFile,'("IterateGeneProb/GeneProb"i0,"/GeneProbs.txt")')h          !here
    open (unit=110,file=trim(inFile),status="unknown")
    StSnp=GpIndex(h,1)          ! Where SNPs start
    EnSnp=GpIndex(h,2)          ! Where SNPs end
    do i=1,nAnisP                                           ! The number of lines of GeneProbs.txt files is = nAnisP x 4
        do j=1,4                                            ! where 4 stands for the two paternal and the two maternal haplotypes
            read (110,*) dum,GeneProbWork(StSnp:EnSnp,j)
        enddo

        GenosProbs(i,StSnp:EnSnp,1) = GeneProbWork(StSnp:EnSnp,1)
        GenosProbs(i,StSnp:EnSnp,2) = GeneProbWork(StSnp:EnSnp,2) + GeneProbWork(StSnp:EnSnp,3)

    enddo
    close(110)
enddo

end subroutine ReReadIterateGeneProbs

