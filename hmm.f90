!######################################################################

module GlobalVariablesHmmMaCH
implicit none

integer, parameter :: GENOTYPE_MISSING=3
integer, parameter :: ALLELE_MISSING=3
integer, parameter :: READ_MISSING=0
integer            :: MISSING=3

double precision, parameter :: SEQUENCING_ERROR=0.01
double precision, parameter :: EPSILON_ERROR=0.00000001


character(len=300) :: GenotypeFileName,CheckPhaseFileName,CheckGenoFileName
integer :: nIndHmmMaCH,GlobalRoundHmm,nSnpHmm,nPhased
integer :: nHapInSubH,useProcs,nRoundsHmm,HmmBurnInRound,phasedThreshold,idum
integer,allocatable,dimension(:,:) :: GenosHmmMaCH,SubH
integer(kind=1),allocatable,dimension(:,:,:) :: PhaseHmmMaCH,FullH
integer,allocatable,dimension(:) :: ErrorUncertainty,ErrorMatches,ErrorMismatches,Crossovers
integer,allocatable,dimension(:) :: GlobalHmmHDInd
logical,allocatable,dimension(:) :: GlobalHmmPhasedInd
double precision,allocatable,dimension(:) :: Thetas,Epsilon
double precision,allocatable,dimension(:,:) :: ForwardProbs
double precision,allocatable,dimension(:,:,:) :: Penetrance, ShotgunErrorMatrix
real,allocatable,dimension (:,:) :: ProbImputeGenosHmm
real,allocatable,dimension (:,:,:) :: ProbImputePhaseHmm
integer, allocatable :: GenosCounts(:,:,:)
!$omp threadprivate(ForwardProbs, SubH)

end module GlobalVariablesHmmMaCH

!######################################################################
subroutine MaCHController(HMM)
use Global
use GlobalVariablesHmmMaCH
use Par_Zig_mod
use omp_lib

implicit none
integer, intent(in) :: HMM
integer :: i, nprocs, nthreads, j
real(4) :: r
real(8) :: t1, t2, tT
double precision :: Theta
integer, allocatable :: seed(:)
integer :: grainsize, count, secs, seed0
integer :: n0, n1, n2

#ifdef DEBUG
    write(0,*) 'DEBUG: [MaCHController] Allocate memory'
#endif

! Number of SNPs and genotyped animals for the HMM algorithm
nSnpHmm=nSnp
nIndHmmMaCH=nAnisG

! ALLOCATE MEMORY
! Allocate a matrix to store the diploids of every Animal
! Template Diploids Library
! NOTE: GenosHmmMaCH can contain either genotype or reads information (if working with sequence data NGS)
allocate(GenosHmmMaCH(nIndHmmMaCH,nSnp))
allocate(PhaseHmmMaCH(nIndHmmMaCH,nSnp,2))
! Allocate memory to store Animals contributing to the Template
! Haplotype Library
allocate(GlobalHmmID(nIndHmmMaCH))

! Allocate memory to store Animals Highly Dense Genotyped
allocate(GlobalHmmHDInd(nIndHmmMaCH))

! Allocate memory to store Animals Highly Dense Genotyped
allocate(GlobalHmmPhasedInd(nIndHmmMaCH))
! No animal has been HD genotyped YET
! WARNING: If this variable only stores 1 and 0, then its type should
!          logical: GlobalHmmHDInd=.false.
GlobalHmmHDInd=0

! Allocate a matrix to store probabilities of genotypes and 
! alleles for each animal
allocate(ProbImputeGenosHmm(nIndHmmMaCH,nSnp))
allocate(ProbImputePhaseHmm(nIndHmmMaCH,nSnp,2))
! Initialise probabilities to 0
ProbImputeGenosHmm=0.0
ProbImputePhaseHmm=0.0

! The full Template Haplotype Library, H (Li et al. 2010, Appendix)
allocate(FullH(nIndHmmMaCH,nSnpHmm,2))

! Vector of Combination of genotyping error (Li et al. 2010, Appendix)
! Epsilon is related with the Penetrance Matrix of the HMM which gives
! the emision probabilities for each state/marker/snp.
allocate(Epsilon(nSnpHmm))
allocate(Penetrance(nSnpHmm,0:2,0:2))

! Vector of Combination of population recombination (Li et al. 2010, Appendix)
! Thetas is related with the transition Matrix of the HMM. Since there
! are nSnpHmm states, there are nSnpHmm transitions between states.
! WARNING: Is this correctly implemented throughout the hmm process??
allocate(Thetas(nSnpHmm-1))

allocate(ErrorUncertainty(nSnpHmm))
allocate(ErrorMatches(nSnpHmm))
allocate(ErrorMismatches(nSnpHmm))

! Crossover parameter in order to maximize the investigation of
! different mosaic configurations
! WARNING: crossovers are related with the transition matrix of the HMM.
!          If there are nSnpHmm states and nSnpHmm-1 transitions
!          between states, why the number of crossovers is nSnpHmm??
allocate(Crossovers(nSnpHmm-1))

if (HMM==RUN_HMM_NGS) then
    ! Set the value MISSING of reads
    MISSING = READ_MISSING

    allocate(ShotgunErrorMatrix(0:2,0:MAX_READS_COUNT,0:MAX_READS_COUNT))
endif

! Set up number of process to be used
nprocs = OMP_get_num_procs()
call OMP_set_num_threads(useProcs)
nthreads = OMP_get_num_threads()

allocate(seed(useProcs))

! Warm up the random seed generator
call system_clock(count)
secs = mod(count,int(1e4))
do i = 1,secs
    call random_number(r)
enddo

! Set up random process seeds for threads
! Feed seed as a function of the milliseconds of the system clock
call system_clock(count)
secs = mod(count,int(1e6))
seed0 = secs*1e5
idum=-abs(seed0)
do i = 1,useProcs
    call random_number(r)
    seed(i) = seed0*r
enddo

grainsize = 32
call par_zigset(useProcs, seed, grainsize)

#ifdef DEBUG
    write(0,*) 'DEBUG: [ParseMaCHData] ...'
#endif

call ParseMaCHData(HMM)

! do i=1,nIndHmmMaCH
!     do j=1,nSnp
!         if (GenosHmmMaCH(i,j)==0) then
!             PhaseHmmMaCH(i,j,:)=0
!         elseif (GenosHmmMaCH(i,j)==2) then
!             PhaseHmmMaCH(i,j,:)=1
!         else
!             PhaseHmmMaCH(i,j,:)=3
!         endif
!     enddo
! enddo

! Initialization of HMM parameters
Epsilon=EPSILON_ERROR
Thetas=0.01

do i=1,nSnpHmm
    call SetPenetrance(i, EPSILON_ERROR)
enddo

if (HMM==RUN_HMM_NGS) then
    call SetShotgunError(SEQUENCING_ERROR)
endif

#ifdef DEBUG
    write(0,*) 'DEBUG: [SetUpEquations] ...'
#endif

if (HMM==RUN_HMM_NGS) then
    call SetUpEquationsReads

else if (HMM==RUN_HMM_ONLY) then
    call SetUpEquationsGenotypesDiploid
else
    call SetUpEquationsGenotypesHaploid
endif

open (unit=6,form='formatted')

print*, ""
print*, " Impute genotypes by HMM"
print*, "    Using", useProcs, "processors of", nprocs

! Allocate and set up variables storing allele frequencies
allocate(GenosCounts(nIndHmmMaCH,nSnpHmm,2))
GenosCounts=0

tT=0.0
do GlobalRoundHmm=1,nRoundsHmm
    write(6, 100, advance="no") char(13),"   HMM Round   ",GlobalRoundHmm
    !write(6, 100) char(13),"   HMM Round   ",GlobalRoundHmm
    100 format (a1, a17, i10)
    call flush(6)

    call ResetCrossovers

    ! Parallelise the HMM process in animals
#if DEBUG.EQ.1
        t1 = omp_get_wtime()
        write(0,*) 'DEBUG: Begin paralellisation [MaCHController]'
#endif
    !$OMP PARALLEL DO DEFAULT(shared)
    !$!OMP DO 
    do i=1,nIndHmmMaCH
        call MaCHForInd(i)
    enddo
    !$!OMP END DO
    !$OMP END PARALLEL DO
#if DEBUG.EQ.1
        t2 = omp_get_wtime()
        tT = tT + (t2-t1)
        write(0,*) 'DEBUG: End paralellisation. Partial time:', t2-t1
#endif

    Theta = 0.01

    ! Update emission probabilities of the HMM process
    call UpdateThetas
    ! Update transition probabilities of the HMM process
    call UpdateErrorRate(Theta)
enddo

#ifdef DEBUG
    write(0,*) 'DEBUG: End paralellisation'
#endif


! Average genotype probability of the different hmm processes
ProbImputeGenosHmm=ProbImputeGenosHmm/(nRoundsHmm-HmmBurnInRound)
ProbImputePhaseHmm=ProbImputePhaseHmm/(nRoundsHmm-HmmBurnInRound)

! Most likely genotype is the genotype that has been sampled most frequently
!IndHmmMaCH
!    do j=1,nSnpHmm
!        n2 = GenosCounts(i,j,3)                           ! Homozygous: 2 case
!        n1 = GenosCounts(i,j,2)                           ! Heterozygous
!        n0 = (nRoundsHmm-HmmBurnInRound) - n1 - n2      ! Homozygous: 0 case
!        if ((n0>n1).and.(n0>n2)) then
!            ProbImputeGenosHmm(i,j)=0
!        elseif (n1>n2) then
!            ProbImputeGenosHmm(i,j)=1
!        else
!            ProbImputeGenosHmm(i,j)=2
!        endif
!    enddo
!enddo

!deallocate(GenosCounts)

end subroutine MaCHController
!######################################################################
subroutine GenosToImputeGenos
use Global
use GlobalVariablesHmmMaCH

implicit none
integer :: i

#ifdef DEBUG
    write(0,*) 'DEBUG: [GenosToImputeGenos]'
#endif


! do i=1,nAnisP
!     ImputeGenos(i,:) = 
! enddo

end subroutine GenosToImputeGenos
!######################################################################
subroutine ParseMaCHData(HMM)
use Global
use GlobalVariablesHmmMaCH

implicit none
integer, intent(in) :: HMM

if (HMM == RUN_HMM_NGS) then
    call ParseMaCHDataNGS
else
    call ParseMaCHDataGenos
endif

end subroutine ParseMaCHData

!######################################################################
subroutine ParseMaCHDataNGS
use Global
use GlobalVariablesHmmMaCH

implicit none
integer :: i

#ifdef DEBUG
    write(0,*) 'DEBUG: [ParseMaCHDataNGS]'
#endif

do i=1,nAnisG
    ! Add animal's diploid to the Diploids Library
    GenosHmmMaCH(i,:)=reads(i,:)
    GlobalHmmID(i)=i
    ! Find individuals sequenced with high coverage
    if ((float(count(GenosHmmMaCH(i,:)/=MISSING))/nSnp)>0.90) then
        ! WARNING: If this variable only stores 1 and 0, then its
        !          type should logical: GlobalHmmHDInd=.true.
        GlobalHmmHDInd(i)=1
    endif
enddo

! Check if the number of Haplotypes the user has considered in the
! Spec file, Sub H (MaCH paper: Li et al. 2010), is reached.
if (nHapInSubH>2*sum(GlobalHmmHDInd(:))) then
    print*, "WARNING! Number of individuals highly-covered is too small"
    print*, "         for the number of Haplotypes in Sub H specified."
    print*, "         Reference haplotypes will be taken from the whole population"
    GlobalHmmHDInd=1
    ! stop
endif

end subroutine ParseMaCHDataNGS

!######################################################################
subroutine ParseMaCHDataGenos
! subroutine ParseMaCHData
use Global
use GlobalVariablesHmmMaCH

implicit none
integer :: i,j,k
integer :: maxHaps      ! Maximum number of haplotypes possible

INTERFACE
  FUNCTION CountPhasedGametes RESULT( gametesPhased )
    use Global
    use GlobalVariablesHmmMaCH

    INTEGER :: gametesPhased, ind
    INTEGER :: CountPhasedAlleles
    INTEGER, ALLOCATABLE :: gamete(:)
  END FUNCTION CountPhasedGametes
END INTERFACE

#ifdef DEBUG
    write(0,*) 'DEBUG: [ParseMaCHDataGenos] ...'
#endif

GlobalHmmPhasedInd=.FALSE.
k=0

! Count the number of phased gametes
nGametesPhased=0
nGametesPhased = CountPhasedGametes()

print *,nGametesPhased

if (nGametesPhased/float(2*nAnisP) > phasedThreshold/100.0) then
! If the number of phased gametes is above a threshold, then
! read the phased information of AlphaImpute, store it in
! PhaseHmmMaCH, and keep track of which gamete is phased (GlobalHmmPhasedInd)
    ! print *, nGametesPhased/float(2*nAnisP), phasedThreshold/100.0
    do i=1,nAnisP
        ! Check if individual is in the genotype file
       if (IndivIsGenotyped(i)==1) then
            k=k+1
            GlobalHmmID(k)=i

            ! Take the phased information from AlphaImpute
            PhaseHmmMaCH(k,:,1)=ImputePhase(i,:,1)
            PhaseHmmMaCH(k,:,2)=ImputePhase(i,:,2)

            ! Clean the alleles from possible coding errors
            do j=1,nSnp
                if ((PhaseHmmMaCH(k,j,1)/=0) .AND. (PhaseHmmMaCH(k,j,1)/=1)) PhaseHmmMaCH(k,j,1)=3
                if ((PhaseHmmMaCH(k,j,2)/=0) .AND. (PhaseHmmMaCH(k,j,2)/=1)) PhaseHmmMaCH(k,j,2)=3
            enddo

            ! Check if this individual has its haplotypes phased
            if (float(count(PhaseHmmMaCH(k,:,1)/=3))/nSnpHmm >= (WellPhasedThresh/100.0)) Then
                GlobalHmmPhasedInd(k,1)=.TRUE.
            endif

            if (float(count(PhaseHmmMaCH(k,:,2)/=3))/nSnpHmm >= (WellPhasedThresh/100.0)) Then
                GlobalHmmPhasedInd(k,2)=.TRUE.
            endif

            ! Count the number of phased animals
            if ((GlobalHmmPhasedInd(k,1)==.TRUE.).AND.(GlobalHmmPhasedInd(k,2)==.TRUE.)) Then
                nAnimPhased=nAnimPhased+1
            endif
        endif

    enddo
    maxHaps = 2*nAnimPhased
else
! If the number of phased gametes is below a threshold, then
! read both the phased information of AlphaImpute and high-denisty genotypes,
! store it in PhaseHmmMaCH and GenosHmmMaCH, and keep track of  which
! gamete is phased (GlobalHmmPhasedInd) and the high-denisity
! genotyped animal (GlobalHmmHDInd)

    do i=1,nAnisP

        ! Check if individual is in the genotype file
        if (IndivIsGenotyped(i)==1) then
            k=k+1
            GlobalHmmID(k)=i

            ! Add animal's diploid to the Diploids Library
            GenosHmmMaCH(k,:)=ImputeGenos(i,:)

            ! Take the phased information from AlphaImpute
            PhaseHmmMaCH(k,:,1)=ImputePhase(i,:,1)
            PhaseHmmMaCH(k,:,2)=ImputePhase(i,:,2)

            ! Check if this animal is Highly Dense genotyped
            if ((float(count(GenosHmmMaCH(k,:)==9))/nSnp)<0.10) then
                GlobalHmmHDInd(k)=1
            endif

            ! Clean the genotypes and alleles from possible coding errors
            do j=1,nSnp
                if ((GenosHmmMaCH(k,j)<0).or.(GenosHmmMaCH(k,j)>2)) GenosHmmMaCH(k,j)=MISSING
                if ((PhaseHmmMaCH(k,j,1)/=0) .AND. (PhaseHmmMaCH(k,j,1)/=1)) PhaseHmmMaCH(k,j,1)=3
                if ((PhaseHmmMaCH(k,j,2)/=0) .AND. (PhaseHmmMaCH(k,j,2)/=1)) PhaseHmmMaCH(k,j,2)=3
            enddo

            ! Check if this individual has its haplotypes phased
            if (float(count(PhaseHmmMaCH(k,:,1)/=3))/nSnpHmm >= (WellPhasedThresh/100.0)) Then
                GlobalHmmPhasedInd(k,1)=.TRUE.
            endif
            if (float(count(PhaseHmmMaCH(k,:,2)/=3))/nSnpHmm >= (WellPhasedThresh/100.0)) Then
                GlobalHmmPhasedInd(k,2)=.TRUE.
            endif

            ! Count the number of phased animals
            if ((GlobalHmmPhasedInd(k,1)==.TRUE.).AND.(GlobalHmmPhasedInd(k,2)==.TRUE.)) Then
                nAnimPhased=nAnimPhased+1
            endif
        endif
    enddo

    maxHaps = 2*sum(GlobalHmmHDInd(:))
endif


! Check if the number of genotyped animals is correct
if (k/=nAnisG) then
    print*, "Error in ParseMaCHDataGenos"
    stop
endif

! Check if the number of Haplotypes the user has considered in the
! Spec file, Sub H (MaCH paper: Li et al. 2010), is reached.
print *, nHapInSubH, maxHaps
if (nHapInSubH>maxHaps) then
    print*, "Data set is too small for the number of Haplotypes in Sub H specified"
    stop
endif

! end subroutine ParseMaCHData
end subroutine ParseMaCHDataGenos

!######################################################################
subroutine MaCHForInd(CurrentInd)
! Create a Template Haplotype Library, H, and create HMM for each
! individual

use Global
use GlobalVariablesHmmMaCH
use random
use Par_Zig_mod
use omp_lib

implicit none

integer, intent(in) :: CurrentInd

! Local variables
!integer :: HapCount, ShuffleInd1, ShuffleInd2, states, thread
integer :: genotype, i, states, thread
integer :: Shuffle1(nIndHmmMaCH), Shuffle2(nIndHmmMaCH)

INTERFACE
  SUBROUTINE ForwardAlgorithmForHaplotype(CurrentInd, hap)
    use Global
    use GlobalVariablesHmmMaCH
    use Par_Zig_mod
    use omp_lib
    implicit none
    integer, intent(in) :: CurrentInd, hap
    ! double precision, intent(IN), allocatable :: Thetas(:)
  END SUBROUTINE ForwardAlgorithmForHaplotype

  SUBROUTINE SampleHaplotypeSource(CurrentInd, hap)
    use Global
    use GlobalVariablesHmmMaCH
    use Par_Zig_mod
    use omp_lib

    implicit none
    integer,intent(IN) :: CurrentInd, hap
  END SUBROUTINE SampleHaplotypeSource

  SUBROUTINE SampleSegmentHaplotypeSource(CurrentInd,hap,snpStart,snpStop)
    use Global
    use GlobalVariablesHmmMaCH
    use Par_Zig_mod
    use omp_lib

    implicit none
    integer,intent(IN) :: CurrentInd, hap, snpStart, snpStop
  END SUBROUTINE SampleSegmentHaplotypeSource
END INTERFACE



! The number of parameters of the HMM are:
!   nHapInSubH = Number of haplotypes in the template haplotype set, H
!   nHapInSubH*(nHapInSubH+1)/2 = Number of states, N = H^2 (Li et al, 2010)
! ForwardProbs are the accumulated probabilities, and
! ForwardPrbos(:,1) are the prior probabilities
! Allocate all possible state sequencies
states = nHapInSubH*(nHapInSubH+1)/2
allocate(ForwardProbs(states,nSnpHmm))
allocate(SubH(nHapInSubH,nSnpHmm))

! Create vectors of random indexes
! Serial
!#ifdef DEBUG
!    write(0,*) 'DEBUG: Shuffle Individuals [MaCHForInd]'
!#endif
!call RandomOrder(Shuffle1,nIndHmmMaCH,idum)
!call RandomOrder(Shuffle2,nIndHmmMaCH,idum)

#if DEBUG.EQ.1
    write(0,*) "DEBUG: RandomOrderPar [MaCHForInd]"
#endif
! Parallel
thread=omp_get_thread_num()
call RandomOrderPar(Shuffle1,nIndHmmMaCH,thread)
call RandomOrderPar(Shuffle2,nIndHmmMaCH,thread)

! Extract haps template (SubH) ...
if (nGametesPhased/float(2*nAnisP)>phasedThreshold/100.0) then
    ! If the number of phased gametes with AlphaImpute is above
    ! a threshold, then template is populated with the phased data
    call ExtractTemplateByHaps(CurrentInd,Shuffle1,Shuffle2)
else
    ! Otherwise, the template is populated with haplotypes at random
    ! from all the HD animals
    call ExtractTemplateHaps(CurrentInd,Shuffle1,Shuffle2)
endif

! ... selecting pairs of haplotypes at random
!call ExtractTemplateHapsByAnimals(CurrentInd,Shuffle1)


! WARNING: This code is something to change according to the Hybrid paper
! if (GlobalHmmPhasedInd(currentInd)==.FALSE.) then
! if (GlobalHmmHDInd(currentInd)==0) then
if (phasedThreshold>nPhased/nAnisP) then ! Training mode 1
    call ForwardAlgorithm(CurrentInd)
    call SampleChromosomes(CurrentInd)
else
    call ForwardAlgorithmForHaplotype(currentInd,1)
    call SampleHaplotypeSource(CurrentInd,1)
    call ForwardAlgorithmForHaplotype(currentInd,2)
    call SampleHaplotypeSource(CurrentInd,2)
    ! TODO: We have to do something here to training mode 3
    ! call ForwardAlgorithm(CurrentInd)
    ! call SampleChromosomes(CurrentInd)
endif

! WARNING: The idea of not to use the first HmmBurnInRound rounds suggests
!          the imputation in those rounds aren't accurate enough, which
!          also suggests that each round a better solution is found.
!          Better solutions are obtained by improving previous solutions
!          by means of update HMM parameters (recombinations  rates,
!          Thetas, and the genotyping errors) as implemented in MaCH
!          code with functions UpdateThetas, UpdateErrorRate and
!          TotalCrossovers.
!
!          However, each time the MaCHForInd subroutine is called is
!          independent from the previous call and so, HMM solutions
!          given by ForwardAlgorithm and SampleChromosomes are also
!          independent.

#if DEBUG.EQ.1
    write(0,*) 'DEBUG: Calculate genotype counts [MaCHForInd]'
#endif
if (GlobalRoundHmm>HmmBurnInRound) then
    do i=1,nSnpHmm
        genotype = FullH(CurrentInd,i,1)+FullH(CurrentInd,i,2)
        if (genotype==2) then
            GenosCounts(CurrentInd,i,2)=GenosCounts(CurrentInd,i,2)+1
        elseif (genotype==1) then
            GenosCounts(CurrentInd,i,1)=GenosCounts(CurrentInd,i,1)+1
        endif
        ! GenosCounts(CurrentInd,i,1)=(GlobalRoundHmm-HmmBurnInRound) - GenosCounts(CurrentInd,i,2) - GenosCounts(CurrentInd,i,3)
    enddo
endif

! Cumulative genotype probabilities through hmm processes
#if DEBUG.EQ.1
    write(0,*) 'DEBUG: Calculate genotype dosages [MaCHForInd]'
#endif
if (GlobalRoundHmm>HmmBurnInRound) then
    ProbImputeGenosHmm(CurrentInd,:)=ProbImputeGenosHmm(CurrentInd,:)&
        +FullH(CurrentInd,:,1)+FullH(CurrentInd,:,2)
    ProbImputePhaseHmm(CurrentInd,:,1)=ProbImputePhaseHmm(CurrentInd,:,1)&
        +FullH(CurrentInd,:,1)
    ProbImputePhaseHmm(CurrentInd,:,2)=ProbImputePhaseHmm(CurrentInd,:,2)&
        +FullH(CurrentInd,:,2)
endif

#if DEBUG.EQ.1
    write(0,*) 'DEBUG: Deallocate Forward variable and Haplotype Template [MaCHForInd]'
#endif
deallocate(ForwardProbs)
deallocate(SubH)

end subroutine MaCHForInd

!######################################################################
subroutine SampleChromosomes(CurrentInd)
use Global
use GlobalVariablesHmmMaCH
use omp_lib
use Par_Zig_mod
implicit none

integer,intent(in) :: CurrentInd

! Local variables
integer :: i,j,k,l,SuperJ,Index,OffOn,State1,State2,TmpJ,TopBot,FirstState,SecondState,Tmp,Thread
double precision :: Summer,ran1,Choice,Sum00,Sum01,Sum10,Sum11
double precision :: Probs(nHapInSubH*(nHapInSubH+1)/2)
double precision :: Theta

Summer=0.0
Index=0
Probs = ForwardProbs(:,nSnpHmm)
Thread = omp_get_thread_num()

! Calculate sum over all states. The sum corresponds to all the
! forward probabilities, that is, the probability of the sequence of
! observed genotypes of animal CurrentInd.
do i=1,nHapInSubH
    do j=1,i
        Index=Index+1
        Summer=Summer+Probs(Index)
    enddo
enddo

! Sample number at random and select state: (State1,State2)
!Choice=ran1(idum)*Summer
Choice = par_uni(Thread)*Summer
Summer=0.0
Index=0
OffOn=0

State1 = 1
State2 = 1

Probs = ForwardProbs(:,nSnpHmm)
do i=1,nHapInSubH
    do j=1,i
        Index=Index+1
        Summer=Summer+Probs(Index)
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
    if (HMMOption/=RUN_HMM_NGS) then
        call ImputeAlleles(CurrentInd,SuperJ+1,State1,State2)
    else
        call ImputeAllelesNGS(CurrentInd,SuperJ+1,State1,State2)
    endif
    TmpJ=SuperJ

    ! Cumulative recombination fraction allows us to skip over
    ! uninformative positions: Alleles with missing genotype are skipped
    ! but the recombination information (Thetas(SuperJ) is accumulated
    ! and used in the next location.
    Theta=Thetas(SuperJ)
    do while ((GenosHmmMaCH(CurrentInd,SuperJ)==MISSING).and.SuperJ>1)
        SuperJ=SuperJ-1
        Theta=Theta+Thetas(SuperJ)-Theta*Thetas(SuperJ)
    enddo

    ! When examining the previous location we consider three alternatives:
    !   * states that could be reached when both haplotypes recombine (11),
    !   * states that can be reached when the first (10) or second (01)
    !       haplotype recombines, and
    !   * the states that can be reached without recombination (00).
    Sum00=0.0
    Sum01=0.0
    Sum10=0.0
    Sum11=0.0

    Index=0
    Probs = ForwardProbs(:,SuperJ)
    do k=1,nHapInSubH
        do l=1,k
            Index=Index+1
            Sum11=Sum11+Probs(Index)
            if ((State1==k).or.(State1==l))&
                Sum01=Sum01+Probs(Index)
            if ((State2==k).or.(State2==l))&
                Sum10=Sum10+Probs(Index)
            if (((State1==k).and.(State2==l))&
                .or.((State1==l).and.(State2==k)))&
                    Sum00=Sum00+Probs(Index)
        enddo
    enddo

    Summer=Sum11*Theta*Theta/(nHapInSubH*nHapInSubH)&
          +(Sum10+Sum01)*Theta*(1.0-Theta)/nHapInSubH&
          +Sum00*(1.0-Theta)*(1.0-Theta)

    ! WARNING: Why is this assignment here?!?! In case it has to exit,
    !          shouldn't it exit before?
    if (SuperJ==1) exit


    !Sample number and decide how many state changes occurred between the
    !two positions
    !Choice=ran1(idum)*Summer
    Choice = par_uni(Thread)*Summer

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

    ! Did the first hap recombine?
    if (Choice<=0.0) then
        !The first haplotype changed ...
        Choice=Choice*nHapInSubH/(Theta*(1.0-Theta))
        !Record the original state
        FirstState=State1

        ! Sample number at random and decide haplotype
!        do while (State1<nHapInSubH)
!            State1=State1+1
        do State1=1,nHapInSubH
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

    ! Did the second hap recombine?
    if (Choice<=0.0) then
        !The second haplotype changed ...
        Choice=Choice*nHapInSubH/(Theta*(1.0-Theta))
        !Save the original state
        SecondState=State2

        ! Sample number at random and decide haplotype
        ! WARNING: State2 variable should be set to 0 before the loop!?!?
!        do while (State2<nHapInSubH)
!            State2=State2+1
        do State2=1,nHapInSubH
            if (State1>=State2) then
                Choice=Choice+ForwardProbs(State1*(State1-1)/2+State2,SuperJ)
            else
               Choice=Choice+ForwardProbs(State2*(State2-1)/2+State1,SuperJ)
            endif
            if (Choice>=0.0) exit
        enddo

        !Record outcomes for intermediate, uninformative, positions
        TopBot=1
        call FillPath(CurrentInd,SuperJ,TmpJ+1,State1,TopBot)
        TopBot=2
        call SamplePath(CurrentInd,SuperJ,TmpJ+1,State2,SecondState,TopBot)
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
    Probs=ForwardProbs(:,nSnpHmm)
    do i=1,nHapInSubH
        do j=1,i
            Index=Index+1
            Summer=Summer+Probs(Index)
            if (Summer>Choice) then
                State1=i
                State2=j
                OffOn=1
                exit
            endif
        enddo
        if (OffOn==1) exit
    enddo

    ! Shuffle haplotypes at random
    !if (ran1(idum)>0.5) then
    if (par_uni(Thread)>0.5) then
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

if (HMMOption/=RUN_HMM_NGS) then
    call ImputeAlleles(CurrentInd,1,State1,State2)
else
    call ImputeAllelesNGS(CurrentInd,1,State1,State2)
endif
 
! call ImputeAlleles(CurrentInd,1,State1,State2)

end subroutine SampleChromosomes

!######################################################################
subroutine FillPath(CurrentInd,FromMarker,ToMarker,State,TopBot)
! Impute alleles to a haplotype region. The region is from FromMarker
! to ToMarker, none of them included.

use GlobalVariablesHmmMaCH
implicit none

integer,intent(in) :: CurrentInd,FromMarker,ToMarker,State,TopBot

! Local variable
integer :: j

do j=FromMarker+1,ToMarker-1
    call ImputeAllele(CurrentInd,j,State,TopBot)
enddo

end subroutine FillPath

!######################################################################
subroutine SamplePath(CurrentInd,FromMarker,ToMarker,FromState,ToState,TopBot)
! Impute a path between the two end markers, assuming no genotypes
! are observed -- the only constraint is that we must start at
! fromState and end at toState with at least one intervening recombinant

use GlobalVariablesHmmMaCH
use omp_lib
use Par_Zig_mod
implicit none

integer,intent(in) :: CurrentInd,TopBot,FromMarker,ToMarker,FromState,ToState

! Local variables
integer :: i,State,FromMarkerLocal, Thread
double precision :: Recomb,Theta1,ran1
double precision :: Theta

Thread = omp_get_thread_num()

FromMarkerLocal=FromMarker
Theta=0.0
! Calculate overall recombination fraction for the interval,
! FromMarker excluded
do i=FromMarkerLocal,ToMarker-1
    Theta=Thetas(i)+Theta-Theta*Thetas(i)
enddo


! Impute a path between the two end markers
do while (FromMarkerLocal<ToMarker-1)
    ! Random recombination
    !Recomb=ran1(idum)*Theta
    Recomb = par_uni(Thread)*Theta

    ! Recombination fraction of the FromMarkerLocal
    Theta1=Thetas(FromMarkerLocal)

    ! Calculate overall recombination fraction for the interval,
    ! FromMarkerLocal included
    if (Theta < 0.9) then
        !Fast closed formula
        Theta=(Theta-Theta1)/(1.0-Theta1)
    else
        Theta = 0.0
        !More accurate, iterative formula
        do i=FromMarkerLocal+1,ToMarker-1
            Theta=Thetas(i)+Theta-Theta*Thetas(i)
        enddo
    endif

    if (Recomb>Theta1) then
        ! No recombinant in the first interval =>
        !    => Recombinant in second interval
        FromMarkerLocal=FromMarkerLocal+1
        ! WARNING: Here was a bug. It there isn't recombination,
        !          imputation has to be done in the FromState haplotype
        !call ImputeAllele(CurrentInd,FromMarkerLocal,ToState,TopBot)
        call ImputeAllele(CurrentInd,FromMarkerLocal,FromState,TopBot)
        cycle
    endif

    ! If there is no recombinant in the second interval, then
    ! there is recombinant in the first...
    !$OMP ATOMIC
    Crossovers(FromMarkerLocal)=Crossovers(FromMarkerLocal)+1

    if (Recomb<Theta1*(1.0-Theta)) then
        ! No recombinant in the second interval
        call FillPath(CurrentInd,FromMarkerLocal,ToMarker,ToState,TopBot);
        return
    else
        ! Recombinants in both intervals, so we must sample
        ! an intervening state
        FromMarkerLocal=FromMarkerLocal+1

        ! WARNING: Not really a bug, but not used a coherent notation
        !ToState=int(ran1(idum)*nHapInSubH)+1
        !call ImputeAllele(CurrentInd,FromMarkerLocal,ToState,TopBot)
        !State=int(ran1(idum)*nHapInSubH)+1
        State=int(par_uni(Thread)*nHapInSubH)+1

        call ImputeAllele(CurrentInd,FromMarkerLocal,State,TopBot)
    endif

enddo

!If we get here, record obligate recombinant between two consecutive markers
!$OMP ATOMIC
Crossovers(FromMarkerLocal)=Crossovers(FromMarkerLocal)+1

end subroutine SamplePath

!######################################################################
subroutine ImputeAlleles(CurrentInd,CurrentMarker,State1,State2)
! Impute alleles to haplotypes based on the HMM information. Count the
! number of uncertainties, matches and mismatches of the imputed
! alleles according to the genotype information that the individual
! carries.
use GlobalVariablesHmmMaCH
use omp_lib
use Par_Zig_mod
implicit none

integer,intent(in) :: CurrentInd,CurrentMarker,State1,State2

! Local variables
integer :: Imputed1,Imputed2,Genotype,Differences, Thread
double precision :: ran1

Thread = omp_get_thread_num()
! These will be the observed imputed alleles defined by the state:
! Sj=(State1, State2)
Imputed1=SubH(State1,CurrentMarker)
Imputed2=SubH(State2,CurrentMarker)

! This is the individual observed genotype
Genotype=GenosHmmMaCH(CurrentInd,CurrentMarker)

if ((Genotype/=0).and.(Genotype/=2)) then
    FullH(CurrentInd,CurrentMarker,1)=Imputed1
    FullH(CurrentInd,CurrentMarker,2)=Imputed2
endif

! If genotype is missing, skip
if (Genotype==3) return

! Difference between the observed genotype and the gentoype implied by
! the state S=(State1, State2)
Differences=abs(Genotype - (Imputed1+Imputed2))

! If allele is heterozygous, there is uncertainty
if ((Genotype==1).and.(Differences==0)) then
    !$OMP ATOMIC
    ErrorUncertainty(CurrentMarker)=ErrorUncertainty(CurrentMarker)+1

! If allele is homozygous or the genotype does not agree the observation
else
    ! count the number of alleles matching
    !$OMP ATOMIC
    ErrorMatches(CurrentMarker)=ErrorMatches(CurrentMarker)+(2-Differences)
    ! count the number of mismatching alleles
    !$OMP ATOMIC
    ErrorMismatches(CurrentMarker)=ErrorMismatches(CurrentMarker)+Differences
endif

! If gentoype is homozygous or missing, the skip
if (Genotype/=1) return

! If the observed allele is homozygous but the genotype is heterozygous
if (Imputed1==Imputed2) then
    ! Impute the allele to the paternal or maternal haplotype at random
    !if (ran1(idum)>=0.5) then
    if (par_uni(Thread)>=0.5) then
        FullH(CurrentInd,CurrentMarker,1)=abs(Imputed1-1)
    else
        ! WARNING: Bug fixed!
        !          Imputation was done to the paternal haplotype also here
        FullH(CurrentInd,CurrentMarker,2)=abs(Imputed2-1)
    endif
endif

end subroutine ImputeAlleles

!######################################################################
subroutine ImputeAllele(CurrentInd,CurrentMarker,State,TopBot)
use GlobalVariablesHmmMaCH
implicit none

integer :: CurrentInd,CurrentMarker,State,TopBot

FullH(CurrentInd,CurrentMarker,TopBot)=SubH(State,CurrentMarker)

end subroutine ImputeAllele

!######################################################################
subroutine ForwardAlgorithm(CurrentInd)
! Update the forward variable of the HMM model

use GlobalVariablesHmmMaCH
use omp_lib

implicit none

integer, intent(in) :: CurrentInd
double precision :: Theta

! Local variables
integer :: i,j,PrecedingMarker

! Setup the initial state distributions

call SetUpPrior

j=1
! CurrentInd is the individual being studied and it is necessary to
! obtain the genotype of this individual in ConditionaOnData subroutine
! For j=1, ConditionaOnData will initialize the variable
! ForwardProbs(:,1) with the Prior Probabilities
#if DEBUG.EQ.1
    write(0,*) 'DEBUG: Update forward variable with emission probabilities [ForwardAlgorithm]'
#endif
call ConditionOnData(CurrentInd,j)

! WARNING: This variable, Theta, should be considered as local as is
!          global through out the HMM code for different purposes.
!          Look at subroutines Transpose, SampleChromosomes and
!          SamplePath
Theta=0.0

PrecedingMarker=1
#if DEBUG.EQ.1
    write(0,*) 'DEBUG: Calculate Forward variables [ForwardAlgorithm]'
#endif

do j=2,nSnpHmm
    ! Cumulative recombination fraction allows us to skip uninformative positions
    Theta=Theta+Thetas(j-1)-Theta*Thetas(j-1)
    ! Skip over uninformative positions to save time
    if ((GenosHmmMaCH(CurrentInd,j)/=MISSING).or.(j==nSnpHmm)) then
        call Transpose(j,PrecedingMarker,Theta)
        call ConditionOnData(CurrentInd,j)
        PrecedingMarker=j
        Theta=0.0
    endif
enddo

end subroutine ForwardAlgorithm

!######################################################################
subroutine Transpose(CurrentMarker,PrecedingMarker,Theta)
! Calculates the probability of get a particular state at CurrentMarker
! from any other state at PrecedingMarker using the transition probabilities.
! It basically calculates the first term of the forward variable.
use GlobalVariablesHmmMaCH
implicit none

integer,intent(in) :: CurrentMarker,PrecedingMarker
double precision,intent(in) :: Theta

!Local variables
integer :: i,j,Index
double precision :: Summer,Marginals(nHapInSubH),NoChange,OneChange,TwoChange

#if DEBUG.EQ.1
    write(0,*) 'DEBUG: [Transpose]'
#endif

if (Theta==0.0) then
    ForwardProbs(:,CurrentMarker)=ForwardProbs(:,PrecedingMarker)
    return
endif

Summer=0.0
Index=0
Marginals(:)=0.0

! Calculate the sum of all forward variables and the marginal of a
! given state.
#if DEBUG.EQ.1
    write(0,*) 'DEBUG: Calculate the acumulate sum of forward variables and the marginals [Transpose]'
#endif
do i=1,nHapInSubH
    do j=1,i-1
        Index=Index+1
        Summer=Summer+ForwardProbs(Index,PrecedingMarker)

        ! Given two states sharing a haplotype, there is only one way to
        ! get one state from the other changing only one haplotype. This
        ! is valid for whatever the shared haplotype. Let's consider the
        ! the sum of all the forward variables corresponding
        ! to the states the state (k1,k2) can be reached from changing
        ! only the first haplotype: Sum{h=1,K}(h,k2). In a similar way,
        ! let's define this other Sum{h=1,k}(k1,h). Since the states
        ! are symmetric and we are working in the upper triangular matrix
        ! both marginal are the same, and each forward variable is summed
        ! twice. Each sum can be called marginal in the first or second
        ! haplotype. Again, since we consider the upper triangular
        ! matrix of the forward variables, we need to define the marginals
        ! in a different way: Sum{h=1,k1}(h,k2)+Sum{h=k1,K}(k1,k). In
        ! this way, the sum of the two marginals will give as a result
        ! that each forward variable is summed twice.
        Marginals(i)=Marginals(i)+ForwardProbs(Index,PrecedingMarker)
        Marginals(j)=Marginals(j)+ForwardProbs(Index,PrecedingMarker)
    enddo
    Index=Index+1
    Summer=Summer+ForwardProbs(Index,PrecedingMarker)

    ! The state (k,k) has to be summed twice
    Marginals(i)=Marginals(i)+(ForwardProbs(Index,PrecedingMarker)*2.0)
enddo

#if DEBUG.EQ.1
    write(0,*) 'DEBUG: Calculate Probabilities of number of changes [Transpose]'
#endif
! NOTE: to the transition probabilities:
!   (1-Theta) is the probability that a haplotype does not change =>
!       => Theta is the probability that a haplotype does change =>
!       => Theta/nHapInSubH is the probability that a haplotype change
!          to a particular haplotype
!
! Given those probabilities, then:
!   * None hapltoyped have changed
NoChange=(1.0-Theta)*(1.0-Theta)
!   * Only one haplotype has changed
OneChange=(1.0-Theta)*Theta/nHapInSubH
!   * The two haplotypes have changed
TwoChange=Summer*Theta*Theta/(nHapInSubH*nHapInSubH)

!Automatically rescale likelihoods when they get too small
#if DEBUG.EQ.1
    write(0,*) 'DEBUG: Rescale [Transpose]'
#endif
if (Summer < 1e-15) then
    NoChange=NoChange*1e30
    OneChange=OneChange*1e30
    TwoChange=TwoChange*1e30
endif

! This final loop actually transposes the probabilities for each state,
! that is, calculates the final probabilities of getting a particular state.
Index=0
#if DEBUG.EQ.1
    write(0,*) 'DEBUG: Calculate the final probabilities [Transpose]'
#endif
do i=1,nHapInSubH
    do j=1,i-1
        Index=Index+1
        ! The new forward probability will be the probability of no change,
        ! plus the probability of one change, plus two times the
        ! probability of two changes (prob of (k1,k2) and prob of (k2,k1))
        ForwardProbs(Index,CurrentMarker)&
            = (ForwardProbs(Index,PrecedingMarker)*NoChange)&
            + (Marginals(i)*OneChange)&
            + (Marginals(j)*OneChange)&
            + (2*TwoChange)
    enddo
    Index=Index+1
    ForwardProbs(Index,CurrentMarker)&
        = (ForwardProbs(Index,PrecedingMarker)*NoChange)&
        + (Marginals(i)*OneChange)&
        + (2*TwoChange)
enddo

end subroutine Transpose

!######################################################################
subroutine ConditionOnData(CurrentInd,Marker)
! Introduce the emission probabilities (the probability of observe a
! given genotype) into the forward variable of the HMM.
! It basically calculate the second term of the forward variables
use Global
use GlobalVariablesHmmMaCH
implicit none

integer, intent(in) :: CurrentInd, Marker

! Local variables
integer :: i, j, Index, genotype, nReads, RefAll, AltAll
double precision :: Factors(0:1), cond_probs(0:2)

! We treat missing genotypes as uninformative about the mosaic's
! underlying state. If we were to allow for deletions and the like,
! that may no longer be true.
! NOTE: gentoype can be refer either to genotypes or reads if working with sequence data (NGS)

if (HMMOption==RUN_HMM_NGS) then
    nReads = AlterAllele(CurrentInd,Marker)*MAX_READS_COUNT+ReferAllele(CurrentInd,Marker)
    RefAll = ReferAllele(CurrentInd,Marker)
    AltAll = AlterAllele(CurrentInd,Marker)
else
    genotype = GenosHmmMaCH(CurrentInd,Marker)
endif



if (genotype==MISSING) then
    return
else
    ! Index keeps track of the states already visited. The total number
    ! of states in this chunk of code is (nHapInSubH x (nHapInSubH-1)/2)
    Index=0
    if (HMMOption==RUN_HMM_NGS) then
        do i=0,2
            ! cond_probs(i)=Penetrance(Marker,i,0)*shotgunErrorMatrix(0,nReads)&
            !              +Penetrance(Marker,i,1)*shotgunErrorMatrix(1,nReads)&
            !              +Penetrance(Marker,i,2)*shotgunErrorMatrix(2,nReads)
            cond_probs(i)=Penetrance(Marker,i,0)*shotgunErrorMatrix(0,RefAll,AltAll)&
                         +Penetrance(Marker,i,1)*shotgunErrorMatrix(1,RefAll,AltAll)&
                         +Penetrance(Marker,i,2)*shotgunErrorMatrix(2,RefAll,AltAll)
        enddo
    endif

    do i=1,nHapInSubH
        if (HMMOption /= RUN_HMM_NGS) then
            ! Probability to observe genotype SubH(i) being the true
            ! genotype GenosHmmMaCH in locus Marker
            Factors(0) = Penetrance(Marker,SubH(i,Marker),genotype)
            ! Probability to observe genotype SubH(i)+1 being the true
            ! genotype GenosHmmMaCH in locus Marker
            Factors(1) = Penetrance(Marker,SubH(i,Marker)+1,genotype)
        else
            ! Probability to observe genotype SubH(i) being the true
            ! genotype GenosHmmMaCH in locus Marker
            Factors(0) = cond_probs(SubH(i,Marker))
            ! Probability to observe genotype SubH(i)+1 being the true
            ! genotype GenosHmmMaCH in locus Marker
            Factors(1) = cond_probs(SubH(i,Marker)+1)
        endif
        do j=1,i
            Index=Index+1
            ForwardProbs(Index,Marker)=&
                ForwardProbs(Index,Marker)*Factors(SubH(j,Marker))
        enddo
    enddo
endif

end subroutine ConditionOnData

!######################################################################
subroutine CalcPenetrance
! Initialize the Penetration matrix of the HMM as the emission
! probabilities matrix given in Appendix of Li et al. (2010)
use GlobalVariablesHmmMaCH
use omp_lib

implicit none

integer :: j, nprocs

! Penetrance(j,i,k) = P(P_j|S_j)
! i = G_j = {0,1,2}
! k = T(S_j) = T(x_j) + T(y_j) = {0,1,2}
! allocate(Penetrance(nSnpHmm,0:2,0:2))

nprocs = OMP_get_num_procs()
call OMP_set_num_threads(nprocs)

!$OMP PARALLEL DO DEFAULT(SHARED)
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
!$OMP END PARALLEL DO

end subroutine CalcPenetrance

!######################################################################
subroutine SetUpPrior
! Set up de initial state distribution, that is, the probability that
! the sequence starts with the state Sj:
!   PIj = P(t1=Sj) = ForwardProb(j,1)

use GlobalVariablesHmmMaCH
use omp_lib

implicit none
integer :: i,j,state
double precision :: prior

 prior=1.0/(nHapInSubH*nHapInSubH)

! Initially, every state is equally possible
!ForwardProbs(:,1)=1.0/(nHapInSubH*nHapInSubH)

! WARNING: I think this variable is treated in a wrong way. In MaCH
!          code, FordwardProbs is the array variable leftMatrices.
!          In there, every element of the array is another array with
!          s(s+1)/2 elements, as S(i,j)=S(j,i); where s=nHapInSubH.
!          So the probability of each state is two times 1 divided by
!          the total number of possible states.
!          So the code should be:

state=0
do i=1,nHapInSubH
   do j=1,i-1
      state=state+1
      ForwardProbs(state,1)=2.0*prior
   enddo
   state=state+1
   ForwardProbs(state,1)=prior
enddo

end subroutine SetUpPrior

!######################################################################
subroutine SetUpEquationsGenotypesHaploid
! Initialize the variables and parameters of the HMM model described in
! Li et al. 2010, Appendix

use Global
use GlobalVariablesHmmMaCH
use random
implicit none

integer :: i,j,p

!Initialise FullH

! If the number of phased gametes from AlphaImpute is above a threshold, then
! haploytpes produced from AlphaImpute are used in the model (FullH)
if (nGametesPhased/float(2*nAnisP)>phasedThreshold/100.0) then
    ! print *, nGametesPhased/float(2*nAnisP), phasedThreshold/100.0
    do i=1,nIndHmmMaCH
        FullH(i,:,:)=PhaseHmmMaCH(i,:,:)

        ! Missing alleles (and individuals that are not phased) are called at random
        do j=1,nSnpHmm
            if (PhaseHmmMaCH(i,j,1)==ALLELE_MISSING) then
                if (ran1(idum)>=0.5) then
                    FullH(i,j,1)=0
                else
                    FullH(i,j,1)=1
                endif
            endif

            if (PhaseHmmMaCH(i,j,2)==ALLELE_MISSING) then
                if (ran1(idum)>=0.5) then
                    FullH(i,j,2)=0
                else
                    FullH(i,j,2)=1
                endif
            endif

        enddo
    enddo

else
    ! If the number of phased gametes from AlphaImpute is below a threshold, then
    ! haplotypes of phased animals from AlphaImpute are used in the model, and
    ! haplotypes of genotyped animals are used otherwise.

    ! Initialise FullH with the Genotype information
    call SetUpEquationsGenotypesDiploid

    ! Overwrite haplotypes to use phased data in case phased haplotypes from
    ! AlphaImpute are available
    do i=1,nIndHmmMaCH      ! For every Individual in the Genotype file
        if (GlobalHmmPhasedInd(i,1)==.TRUE.) then
            FullH(i,:,1)=PhaseHmmMaCH(i,:,1)
            !  If there is missing information in the phased data, called allele at random
            do j=1,nSnpHmm
                if (PhaseHmmMaCH(i,j,1)==ALLELE_MISSING) then
                    if (ran1(idum)>=0.5) then
                        FullH(i,j,1)=0
                    else
                        FullH(i,j,1)=1
                    endif
                endif
            enddo
        endif

        if (GlobalHmmPhasedInd(i,2)==.TRUE.) then
            FullH(i,:,2)=PhaseHmmMaCH(i,:,2)
            !  If there is missing information in the phased data, called allele at random
            do j=1,nSnpHmm

                if (PhaseHmmMaCH(i,j,2)==ALLELE_MISSING) then
                    if (ran1(idum)>=0.5) then
                        FullH(i,j,2)=0
                    else
                        FullH(i,j,2)=1
                    endif
                endif
            enddo
        endif

    enddo
endif

ErrorUncertainty(:)=0
ErrorMatches(:)=0
ErrorMismatches(:)=0
Crossovers(:)=0

end subroutine SetUpEquationsGenotypesHaploid

!######################################################################
subroutine SetUpEquationsGenotypesDiploid
! Initialize the variables and parameters of the HMM model described in
! Li et al. 2010, Appendix

use Global
use GlobalVariablesHmmMaCH
use random
implicit none

integer :: i,j,p

!Initialise FullH
do i=1,nIndHmmMaCH      ! For every Individual in the Genotype file
    do j=1,nSnpHmm      ! For each SNP

        ! Phase homozygose locus
        if (GenosHmmMaCH(i,j)==0) then
            FullH(i,j,:)=0
        endif

        if (GenosHmmMaCH(i,j)==2) then
            FullH(i,j,:)=1
        endif

        ! Phase heterozygose case at random
        if (GenosHmmMaCH(i,j)==1) then
            if (ran1(idum)>=0.5) then
                FullH(i,j,1)=0
                FullH(i,j,2)=1
            else
                FullH(i,j,1)=1
                FullH(i,j,2)=0
            endif
        endif

        ! If locus is not genotyped, phase each haplotype at random
        if (GenosHmmMaCH(i,j)==MISSING) then
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

ErrorUncertainty(:)=0
ErrorMatches(:)=0
ErrorMismatches(:)=0
Crossovers(:)=0

end subroutine SetUpEquationsGenotypesDiploid


!######################################################################
subroutine SetUpEquationsReads
! Initialize the variables and parameters of the HMM model described in
! Li et al. 2010, Appendix

use Global
use GlobalVariablesHmmMaCH
use random
implicit none

integer :: i,j,p, alleles, mac, readObs, RefAll, AltAll
double precision :: prior_11, prior_12, prior_22
double precision :: posterior_11, posterior_12, posterior_22
double precision :: r, frequency, summ

! ! Initialization of HMM parameters
! Epsilon=0.00000001
! Thetas=0.01

!Initialise FullH
do j=1,nSnpHmm      ! For each SNP
    ! alleles = 0
    ! mac = 0
    ! do i=1,nIndHmmMaCH      ! For every Genotyped Individual
    !     readObs = GenosHmmMaCH(i,j)
    !     alleles = alleles + mod(readObs,MAX_READS_COUNT)
    !     mac = mac + readObs/MAX_READS_COUNT
    ! enddo
    ! alleles = alleles + mac

    ! if (alleles==0) then
    !     do i=1,nIndHmmMaCH
    !         FullH(i,j,:)=0
    !     enddo
    ! endif

    ! frequency = DBLE(mac) / DBLE(alleles)
    

    readObs = 0
    alleles = 0

    do i=1,nIndHmmMaCH
        readObs = readObs + GenosHmmMaCH(i,j)
        alleles = alleles + AlterAllele(i,j)
    enddo
    frequency =  dble(alleles) / dble(readObs)

    prior_11 = (1.0 - frequency)**2
    prior_12 = 2.0 * (1.0 - frequency) * frequency
    prior_22 = frequency**2
    ! print *, readObs, alleles, frequency, prior_11, prior_12, prior_22

    do i=1,nIndHmmMaCH
         ! readObs = GenosHmmMaCH(i,j)
         readObs = AlterAllele(i,j)*MAX_READS_COUNT+ReferAllele(i,j)
        RefAll = ReferAllele(i,j)
        AltAll = AlterAllele(i,j)         

         if (readObs==0) exit

         ! posterior_11 = prior_11 * ShotgunErrorMatrix(0,readObs)
         ! posterior_12 = prior_12 * ShotgunErrorMatrix(1,readObs)
         ! posterior_22 = prior_22 * ShotgunErrorMatrix(2,readObs)
         posterior_11 = prior_11 * ShotgunErrorMatrix(0,RefAll,AltAll)
         posterior_12 = prior_12 * ShotgunErrorMatrix(1,RefAll,AltAll)
         posterior_22 = prior_22 * ShotgunErrorMatrix(2,RefAll,AltAll)
         summ = posterior_11 + posterior_12 + posterior_22

         ! print *, summ, posterior_11, posterior_12, posterior_22

         if (summ==0) write(0,*) 'There is a problem here!'

         posterior_11 = posterior_11 / summ
         posterior_12 = posterior_12 / summ

         r = ran1(idum)

         if (r < posterior_11) then
            FullH(i,j,:) = 0
         elseif (r < posterior_11 + posterior_12) then
            if (ran1(idum)<0.5) then
                FullH(i,j,1) = 0
                FullH(i,j,2) = 1
            else
                FullH(i,j,1) = 1
                FullH(i,j,2) = 0
            endif
        else
            FullH(i,j,:)=1
        endif
    enddo
enddo

! call CalcPenetrance

ErrorUncertainty(:)=0
ErrorMatches(:)=0
ErrorMismatches(:)=0
Crossovers(:)=0
!print*, Crossovers(:)

end subroutine SetUpEquationsReads

!######################################################################
subroutine UpdateThetas

use GlobalVariablesHmmMaCH
implicit none

integer :: i, BaseCount=1, BaseIntervals=0
double precision :: BaseRates, BaseCrossovers=1

double precision :: Scale

Scale=1.0/(nIndHmmMaCH*2)
BaseCount=1
BaseIntervals=0

BaseCrossovers=1.0

! First we estimate a base line rate to be applied to intervals with
! 0 or 1 observed "crossovers"
do i=1,nSnpHmm-1
    if (Crossovers(i)<=1) then
        BaseCount=BaseCount+Crossovers(i)
        BaseIntervals=BaseIntervals+1
    endif
enddo

if (BaseIntervals==0) then
    BaseRates=BaseCount*Scale
else
    BaseRates=BaseCount*Scale/BaseIntervals
endif

! Then we update the rate for each interval using either the number
! of observed crossovers (if > 1) or the baseline rate
do i=1,nSnpHmm-1
    if (Crossovers(i)>1) then
        Thetas(i)=Crossovers(i)*Scale
    else
        Thetas(i)=BaseRates
    endif
enddo
!print*,Thetas(nSnpHmm-1)

end subroutine UpdateThetas

!######################################################################
subroutine UpdateErrorRate(rate)
! Group markers into those with low error rates, which are estimated
! as a group, and those with high error rates, which are estimated 
! individually

use GlobalVariablesHmmMaCH
implicit none

double precision,intent(out) :: rate

! Local variables
integer :: i,matches=0,mismatches=0,uncertain=0

rate=0.0
do i=1,nSnpHmm
    if (ErrorMismatches(i)<=2) then
        matches=matches+ErrorMatches(i)
        mismatches=mismatches+ErrorMismatches(i)
        uncertain=uncertain+ErrorUncertainty(i)
    else
        call UpdateError(ErrorMatches(i), ErrorMismatches(i), ErrorUncertainty(i), rate)
        call SetPenetrance(i,rate)
    endif
enddo

call UpdateError(matches, mismatches, uncertain, rate)

do i=1,nSnpHmm
    if (ErrorMismatches(i)<=2) call SetPenetrance(i, rate)
enddo

end subroutine UpdateErrorRate

!######################################################################
subroutine SetPenetrance(marker, Err)

use GlobalVariablesHmmMaCH
implicit none

integer,intent(in) :: marker
double precision, intent(in) :: Err

Epsilon(marker) = Err

Penetrance(marker,0,0)=(1.0-Err)**2
Penetrance(marker,0,1)=2.0*(1.0-Err)*Err
Penetrance(marker,0,2)=Err**2
Penetrance(marker,1,0)=(1.0-Err)*Err
Penetrance(marker,1,1)=((1.0-Err)**2)+(Err**2)
Penetrance(marker,1,2)=(1.0-Err)*Err
Penetrance(marker,2,0)=Err**2   
Penetrance(marker,2,1)=2.0*(1.0-Err)*Err
Penetrance(marker,2,2)=(1.0-Err)**2   

end subroutine SetPenetrance


!######################################################################
subroutine TotalCrossovers(Total)

use GlobalVariablesHmmMaCH
implicit none

integer,intent(out) :: Total

! Local variables
integer :: i

Total=0
do i=1,nSnpHmm-1
    Total=Total+Crossovers(i)
enddo

end subroutine TotalCrossovers

!######################################################################
subroutine ResetCrossovers

use GlobalVariablesHmmMaCH
implicit none

Crossovers(:)=0
call ResetErrors

end subroutine ResetCrossovers

!######################################################################
subroutine UpdateError(matches, mismatches, uncertain, rate)

use GlobalVariablesHmmMaCH
implicit none

integer,intent(in) :: matches,mismatches,uncertain
double precision, intent(out) :: rate

! Local variables
double precision :: previous=0.0, ratio

rate=0.0    ! Just in case...
if (matches+mismatches>0) then
    rate=mismatches/dble(matches+mismatches)
    if (uncertain>0) then
        do while((rate>1e-10).and.(abs(rate-previous) > rate*1e-4))
            ratio=rate*rate/(rate*rate+(1.0-rate)*(1.0-rate))
            previous=rate
            rate=(mismatches+ratio*uncertain*2.0)&
                 /(matches+mismatches+uncertain*2)
        enddo
    endif
else if (uncertain>0) then
    rate=0.0
endif

end subroutine UpdateError

!######################################################################
subroutine ResetErrors

use GlobalVariablesHmmMaCH
implicit none

ErrorUncertainty(:)=0
ErrorMatches(:)=0
ErrorMismatches(:)=0

end subroutine ResetErrors

!######################################################################
subroutine GetErrorRatebyMarker(marker, Err)
use GlobalVariablesHmmMaCH

implicit none
integer, intent(in) :: marker
double precision, intent(out) :: Err

Err = 0.0
Err = Epsilon(marker)

end subroutine GetErrorRatebyMarker

!######################################################################
subroutine GetErrorRate(mean)
use GlobalVariablesHmmMaCH

implicit none
double precision, intent(out) :: mean

double precision :: ErrorRate
integer :: i

mean = 0.0
do i=1,nSnpHmm
    call GetErrorRatebyMarker(i, ErrorRate)
    mean = mean + ErrorRate
enddo
mean = mean / nSnpHmm

end subroutine GetErrorRate

!######################################################################
subroutine ExtractTemplateHaps(forWhom,Shuffle1,Shuffle2)
! Set the Template of Haplotypes used in the HMM model.
! It takes the haplotypes from HD animals
! It differentiates between paternal (even) and maternal (odd) haps

use GlobalVariablesHmmMaCH

integer, intent(in) :: Shuffle1(nIndHmmMaCH), Shuffle2(nIndHmmMaCH)

! Local variables
integer :: HapCount, ShuffleInd1, ShuffleInd2

HapCount=0
ShuffleInd1=0
ShuffleInd2=0

#if DEBUG.EQ.1
    write(0,*) 'DEBUG: Create Haplotypes Template [ExtractTemplateHaps]'
#endif
! While the maximum number of haps in the template haplotypes set,
! H, is not reached...
do while (HapCount<nHapInSubH)
    if (mod(HapCount,2)==0) then
        ShuffleInd1=ShuffleInd1+1

        ! Select the paternal haplotype if the individual it belongs
        ! to is genotyped and it is not the current individual
        if ((Shuffle1(ShuffleInd1)/=CurrentInd)&
                .and.(GlobalHmmHDInd(ShuffleInd1)==1)) then
            HapCount=HapCount+1
            SubH(HapCount,:)=FullH(Shuffle1(ShuffleInd1),:,1)
        endif
    else
        ShuffleInd2=ShuffleInd2+1

        ! Select the maternal haplotype if the individual it belongs
        ! too is genotyped and it is not the current individual
        if ((Shuffle2(ShuffleInd2)/=CurrentInd)&
                .and.(GlobalHmmHDInd(ShuffleInd2)==1)) then
            HapCount=HapCount+1
            SubH(HapCount,:)=FullH(Shuffle2(ShuffleInd2),:,2)
        endif
    endif
enddo

end subroutine ExtractTemplateHaps

!######################################################################
subroutine ExtractTemplateByHaps(forWhom,Shuffle1,Shuffle2)
! Set the Template of Haplotypes used in the HMM model.
! It takes the haplotypes produced by AlphaImpute
! It differentiates between paternal (even) and maternal (odd) haps

use GlobalVariablesHmmMaCH

integer, intent(in) :: Shuffle1(nIndHmmMaCH), Shuffle2(nIndHmmMaCH)

! Local variables
integer :: HapCount, ShuffleInd1, ShuffleInd2

HapCount=0
ShuffleInd1=0
ShuffleInd2=0

#if DEBUG.EQ.1
    write(0,*) 'DEBUG: Create Haplotypes Template [ExtractTemplateByHaps]'
#endif
! While the maximum number of haps in the template haplotypes set,
! H, is not reached...
do while (HapCount<nHapInSubH)
    if (mod(HapCount,2)==0) then
        ShuffleInd1=ShuffleInd1+1

        ! Select the paternal haplotype if the individual it belongs
        ! to is genotyped and it is not the current individual
        if ((Shuffle1(ShuffleInd1)/=CurrentInd).and.&
                (GlobalHmmPhasedInd(ShuffleInd1,1)==.TRUE.)) then

            HapCount=HapCount+1
            SubH(HapCount,:)=FullH(Shuffle1(ShuffleInd1),:,1)
        endif
    else
        ShuffleInd2=ShuffleInd2+1

        ! Select the maternal haplotype if the individual it belongs
        ! too is genotyped and it is not the current individual
        if ((Shuffle2(ShuffleInd2)/=CurrentInd).and.&
                (GlobalHmmPhasedInd(ShuffleInd2,2)==.TRUE.)) then
            HapCount=HapCount+1
            SubH(HapCount,:)=FullH(Shuffle2(ShuffleInd2),:,2)
        endif
    endif
enddo

end subroutine ExtractTemplateByHaps

!######################################################################
subroutine ExtractTemplateHapsByAnimals(forWhom,Shuffle)
! Extract the Template of Haplotypes used in the HMM model.
! It consideres the two haps of a given individual.

use GlobalVariablesHmmMaCH

integer, intent(in) :: Shuffle(nIndHmmMaCH)

! Local variables
integer :: HapCount, ShuffleInd

HapCount=1
ShuffleInd=0

#if DEBUG.EQ.1
    write(0,*) 'DEBUG: Create Haplotypes Template [ExtractTemplateHapsByAnimals]'
#endif
! While the maximum number of haps in the template haplotypes set,
! H, is not reached...
do while (HapCount<nHapInSubH)
    ShuffleInd=ShuffleInd+1

    ! Select the paternal and maternal haplotypes from one individual
    ! if this individual is not being phased/imputed
    if ((Shuffle(ShuffleInd)/=CurrentInd)&
            .and.(GlobalHmmHDInd(ShuffleInd)==1)) then
        SubH(HapCount,:)=FullH(Shuffle(ShuffleInd),:,1)
        SubH(HapCount+1,:)=FullH(Shuffle(ShuffleInd),:,2)
        HapCount=HapCount+2
    endif
enddo

end subroutine ExtractTemplateHapsByAnimals

!######################################################################
subroutine SetShotgunError(ErrorRate)
use Global
use GlobalVariablesHmmMaCH

implicit none

double precision, intent(in) :: ErrorRate

! Local variables
double precision :: DFactorialInLog, ProdFactTmp
integer :: i,k,MaxReadCounts
! integer :: binomial(33,33)


do k=0,MAX_READS_COUNT-1
    do i=0,MAX_READS_COUNT-1
    ! do i=1,MAX_READS_COUNT
        ProdFactTmp=DFactorialInLog(k+i)-(DFactorialInLog(i)+DFactorialInLog(k))

        ! ShotgunErrorMatrix(0,k*MAX_READS_COUNT+i)=exp(ProdFactTmp+(dfloat(i)*log(ErrorRate))+(dfloat(k)*log(1.0-ErrorRate)))
        ! ShotgunErrorMatrix(1,k*MAX_READS_COUNT+i)=exp(ProdFactTmp+(dfloat(k+i)*log(0.5)))
        ! ShotgunErrorMatrix(2,k*MAX_READS_COUNT+i)=exp(ProdFactTmp+(dfloat(i)*log(1.0-ErrorRate))+(dfloat(k)*log(ErrorRate)))
        ShotgunErrorMatrix(0,k,i)=exp(ProdFactTmp+(dfloat(i)*log(ErrorRate))+(dfloat(k)*log(1.0-ErrorRate)))
        ShotgunErrorMatrix(1,k,i)=exp(ProdFactTmp+(dfloat(k+i)*log(0.5)))
        ShotgunErrorMatrix(2,k,i)=exp(ProdFactTmp+(dfloat(i)*log(1.0-ErrorRate))+(dfloat(k)*log(ErrorRate)))
    enddo
enddo

! binomial(1,1) = 1
! binomial(2,1) = 1
! binomial(2,2) = 1

! do i=3,33
!     binomial(i,1) = 1
!     binomial(i,i) = 1
!     do j=2,i 
!         binomial(i,j) = binomial(i-1,j) + binomial(i-1,j-1)
!     enddo
! enddo

! do i=1,16
!     do j=1,16
!         if (rate==0) then
!             if (j==1) then
!                 ShotgunErrorMatrix(0,j*16+i) = 1.0
!             else
!                 ShotgunErrorMatrix(0,j*16+i) = 0.0
!             endif
!             ShotgunErrorMatrix(1,j*16+i) = (0.5)**(i+j) * binomial(i+j,i)
!             if (i==1) then
!                 ShotgunErrorMatrix(2,j*16+i) = 1.0
!             else
!                 ShotgunErrorMatrix(2,j*16+i) = 0.0
!             endif
!         else
!             ShotgunErrorMatrix(0,j*16+i) = (1-rate)**i * (rate)**j * binomial(i+j,i)
!             ShotgunErrorMatrix(1,j*16+i) = (0.5)**(i+j) * binomial(i+j,i)
!             ShotgunErrorMatrix(2,j*16+i) = (rate)**i * (1-rate)**j * binomial(i+j,i)
!         endif
!     enddo
! enddo

end subroutine SetShotgunError

!######################################################################
subroutine ImputeAllelesNGS(CurrentInd,CurrentMarker,State1,State2)
use Global
use GlobalVariablesHmmMaCH
use omp_lib
use Par_Zig_mod

implicit none

integer, intent(in) :: CurrentInd, CurrentMarker, State1, State2

! Local variables
integer :: copied1, copied2, nReads, Thread, bin, imputed1, imputed2, Differences, RefAll, AltAll
double precision :: posterior_11, posterior_12, posterior_22, summ, random, rate

Thread = omp_get_thread_num()

copied1 = SubH(State1,CurrentMarker)
copied2 = SubH(State2,CurrentMarker)

! nReads = GenosHmmMaCH(CurrentInd,CurrentMarker)
nReads = AlterAllele(CurrentInd,CurrentMarker)*MAX_READS_COUNT+ReferAllele(CurrentInd,CurrentMarker)
RefAll = ReferAllele(CurrentInd,CurrentMarker)
AltAll = AlterAllele(CurrentInd,CurrentMarker)

! posterior_11 = Penetrance(CurrentMarker,copied1+copied2,0)*shotgunErrorMatrix(0,nReads)
! posterior_12 = Penetrance(CurrentMarker,copied1+copied2,1)*shotgunErrorMatrix(1,nReads)
! posterior_22 = Penetrance(CurrentMarker,copied1+copied2,2)*shotgunErrorMatrix(2,nReads)
posterior_11 = Penetrance(CurrentMarker,copied1+copied2,0)*shotgunErrorMatrix(0,RefAll,AltAll)
posterior_12 = Penetrance(CurrentMarker,copied1+copied2,1)*shotgunErrorMatrix(1,RefAll,AltAll)
posterior_22 = Penetrance(CurrentMarker,copied1+copied2,2)*shotgunErrorMatrix(2,RefAll,AltAll)

summ = posterior_11 + posterior_12 + posterior_22

posterior_11 = posterior_11 / summ
posterior_22 = posterior_22 / summ

random = par_uni(Thread)

! if (RefAll+AltAll==0) then
!     print *, "\n", summ, posterior_11 + posterior_22
!     print *, posterior_11, shotgunErrorMatrix(0,RefAll,AltAll)!, Penetrance(CurrentMarker,copied1+copied2,0)
!     print *, posterior_12/summ, shotgunErrorMatrix(1,RefAll,AltAll)!, Penetrance(CurrentMarker,copied1+copied2,1)
!     print *, posterior_22, shotgunErrorMatrix(2,RefAll,AltAll)!, Penetrance(CurrentMarker,copied1+copied2,2)
! endif

if (random < posterior_11) then
    ! print*, "posterior_11", random, posterior_11, AlterAllele(CurrentInd,CurrentMarker), ReferAllele(CurrentInd,CurrentMarker), nReads
    FullH(CurrentInd,CurrentMarker,1) = 0
    FullH(CurrentInd,CurrentMarker,2) = 0

elseif (random < posterior_11 + posterior_22) then
    ! print*, "posterior_11+posterior_22", random, posterior_11 + posterior_22, AlterAllele(CurrentInd,CurrentMarker), ReferAllele(CurrentInd,CurrentMarker), nReads
    FullH(CurrentInd,CurrentMarker,1) = 1
    FullH(CurrentInd,CurrentMarker,2) = 1

elseif (copied1 /= copied2) then
    call GetErrorRatebyMarker(CurrentMarker, rate)

    if (par_uni(Thread) < rate*rate / ((rate*rate) + (1-rate)*(1-rate))) then
        if (copied1 == 1) then 
            copied1 = 0
        else
            copied1 = 1
        endif

        if (copied2 == 1) then 
            copied2 = 0
        else
            copied2 = 1
        endif
    endif

    FullH(CurrentInd,CurrentMarker,1) = copied1
    FullH(CurrentInd,CurrentMarker,2) = copied2
else
    if (par_uni(Thread)<0.5) then
        FullH(CurrentInd,CurrentMarker,1) = 0
        FullH(CurrentInd,CurrentMarker,2) = 1
    else 
        FullH(CurrentInd,CurrentMarker,1) = 1
        FullH(CurrentInd,CurrentMarker,2) = 0
    endif
endif

imputed1 = FullH(CurrentInd,CurrentMarker,1)
imputed2 = FullH(CurrentInd,CurrentMarker,2)

! if (CurrentInd==1) then
!     print *, FullH(CurrentInd,1:10,1)
!     print *, FullH(CurrentInd,1:10,2)
! endif


Differences = abs(copied1 - imputed1) + abs(copied2 - imputed2)
! count the number of alleles matching
!$OMP ATOMIC
ErrorMatches(CurrentMarker)=ErrorMatches(CurrentMarker)+(2-Differences)
! count the number of mismatching alleles
!$OMP ATOMIC
ErrorMismatches(CurrentMarker)=ErrorMismatches(CurrentMarker)+Differences

end subroutine ImputeAllelesNGS

!######################################################################################################################################################

real(kind=8) function DFactorialInLog(n)

implicit none
integer,intent(in) :: n

! Local variables
integer :: i
real(kind=8) :: Ans

Ans=0.0
if (n==0) then
    Ans=1.0
else
    do i=1,n
        Ans=Ans+log(dfloat(i))
    enddo
endif

DFactorialInLog=Ans

end function DFactorialInLog

!######################################################################################################################################################
