!######################################################################

module GlobalVariablesHmmMaCH
implicit none

character(len=300) :: GenotypeFileName,CheckPhaseFileName,CheckGenoFileName
integer :: nIndHmmMaCH,GlobalRoundHmm,nSnpHmm
integer :: nHapInSubH,idum,nRoundsHmm,HmmBurnInRound
double precision :: Theta
integer(kind=1),allocatable,dimension(:,:) :: GenosHmmMaCH,SubH
integer(kind=1),allocatable,dimension(:,:,:) :: FullH
integer,allocatable,dimension(:) :: ErrorUncertainty,ErrorMatches,ErrorMismatches,Crossovers,GlobalHmmHDInd
double precision,allocatable,dimension(:) :: Thetas,Epsilon
double precision,allocatable,dimension(:,:) :: ForwardProbs
double precision,allocatable,dimension(:,:,:) :: Penetrance
real,allocatable,dimension (:,:) :: ProbImputeGenosHmm
!$omp threadprivate(ForwardProbs, SubH)
end module GlobalVariablesHmmMaCH

!######################################################################
subroutine MaCHController
use GlobalVariablesHmmMaCH
use omp_lib

implicit none
integer :: i, nprocs, nthreads, useProcs=1

call ParseMaCHData
call SetUpEquations

open (unit=6,form='formatted')

nprocs = OMP_get_num_procs()
call OMP_set_num_threads(useProcs)
nthreads = OMP_get_num_threads()

print*, ""
print*, " Impute genotypes by HMM"
print*, "    Using", useProcs, "processors of", nprocs

do GlobalRoundHmm=1,nRoundsHmm
    !write(6, 100) "   HMM Round   ",GlobalRoundHmm
    !100 format ('+', a17,i10)
    write(6, 100) char(13),"   HMM Round   ",GlobalRoundHmm
    100 format (a1, a17, i10)

    call ResetCrossovers

!    print*, "Allocate SubH"
!    allocate(SubH(nHapInSubH,nSnpHmm))
!    print*, "Allocate ForwardProbs"
!    allocate(ForwardProbs(nHapInSubH*(nHapInSubH+1)/2,nSnpHmm))

    !$OMP PARALLEL DO DEFAULT(shared)
    !$!OMP DO 
    do i=1,nIndHmmMaCH
        !print*, "Animal:", i, "Number of thread:", OMP_get_thread_num()
        !print*, "#. Thread = ", OMP_get_thread_num()
        !allocate(ForwardProbs(nHapInSubH*(nHapInSubH+1)/2,nSnpHmm))
        !allocate(SubH(nHapInSubH,nSnpHmm))
        !print*, 'Number of thread:', OMP_get_thread_num(), allocated(SubH), allocated(ForwardProbs)
        !print*, 'call MaCHForInd'
        call MaCHForInd(i)
        !deallocate(ForwardProbs)
        !deallocate(SubH)
    enddo
    !$!OMP END DO
    !$OMP END PARALLEL DO
    !deallocate(SubH)
    !deallocate(ForwardProbs)

    Theta = 0.01
    call UpdateThetas
    call UpdateErrorRate(Theta)
enddo


!close (6)

! Average genotype probability of the different hmm processes
ProbImputeGenosHmm=ProbImputeGenosHmm/(nRoundsHmm-HmmBurnInRound)

end subroutine MaCHController

!######################################################################
subroutine ParseMaCHData
use Global
use GlobalVariablesHmmMaCH

implicit none
integer :: i,j,k

! Number of SNPs and genotyped animals for the HMM algorithm
nSnpHmm=nSnp
nIndHmmMaCH=nAnisG

! ALLOCATE MEMORY

! Test allocation ForwardProb variable for HMM parallelisation
!allocate(ForwardProbs(nHapInSubH*nHapInSubH,nSnpHmm))

! Allocate a matrix to store the diploids of every Animal
! Template Diploids Library
allocate(GenosHmmMaCH(nIndHmmMaCH,nSnp))

! Allocate memory to store Animals contributing to the Template
! Haplotype Library
allocate(GlobalHmmID(nIndHmmMaCH))

! Allocate memory to store Animals Highly Dense Genotyped
allocate(GlobalHmmHDInd(nIndHmmMaCH))

! Allocate a matrix to store probabilities of the genotype of every
! Animal and inititalize to 0
allocate(ProbImputeGenosHmm(nIndHmmMaCH,nSnp))
ProbImputeGenosHmm=0.0

! Any animal hasn't been HD genotyped YET
! WARNING: If this variable only stores 1 and 0, then its type should
!          logical: GlobalHmmHDInd=.false.
GlobalHmmHDInd=0

k=0
do i=1,nAnisP
    ! Check if individual is genotype
    if (IndivIsGenotyped(i)==1) then
        k=k+1
        ! Add animal's diploid to the Diploids Library
        GenosHmmMaCH(k,:)=ImputeGenos(i,:)
        GlobalHmmID(k)=i
        ! Check if this animal is Highly Dense genotyped
        if ((float(count(GenosHmmMaCH(k,:)==9))/nSnp)<0.10) then
            ! WARNING: If this variable only stores 1 and 0, then its
            !          type should logical: GlobalHmmHDInd=.true.
            GlobalHmmHDInd(k)=1
        endif

        ! WARNING: This should have been previously done for ImputeGenos variable
        do j=1,nSnp
            if ((GenosHmmMaCH(k,j)<0).or.(GenosHmmMaCH(k,j)>2)) GenosHmmMaCH(k,j)=3
        enddo
    endif
enddo

! Check if the number of genotyped animals is correct
if (k/=nAnisG) then
    print*, "Error in ParseMaCHData"
    stop
endif

! Check if the number of Haplotypes the user has considered in the
! Spec file, Sub H (MaCH paper: Li et al. 2010), is reached.
if (nHapInSubH>2*sum(GlobalHmmHDInd(:))) then
    print*, "Data set is too small for the number of Haplotypes in Sub H specified"
    stop
endif

end subroutine ParseMaCHData

!######################################################################
subroutine MaCHForInd(CurrentInd)
! Create a Template Haplotype Library, H, and create HMM for each
! individual

use GlobalVariablesHmmMaCH
use omp_lib

implicit none

integer, intent(in) :: CurrentInd

! Local variables
integer :: HapCount,ShuffleInd1,ShuffleInd2, states
integer :: Shuffle1(nIndHmmMaCH),Shuffle2(nIndHmmMaCH)

!print*, "Allocate SubH"
!print*, "Allocate ForwardProbs"
allocate(ForwardProbs(nHapInSubH*(nHapInSubH+1)/2,nSnpHmm))
allocate(SubH(nHapInSubH,nSnpHmm))
!print*, 'Animal: ', CurrentInd, 'Number of thread:', OMP_get_thread_num(), allocated(SubH), allocated(ForwardProbs)

! Create vectors of random indexes
call RandomOrder(Shuffle1,nIndHmmMaCH,idum)
call RandomOrder(Shuffle2,nIndHmmMaCH,idum)
!print*, 'Inside MaCHForInd'

HapCount=0
ShuffleInd1=0
ShuffleInd2=0

! EXTRACT SUBH
! While the maximum number of haps in the template haplotypes set, H,
! is not reached...
! WARNING: This should be an independent subroutine (as in MaCH code)
do while (HapCount<nHapInSubH)
    ! Differentiate between paternal (even) and maternal (odd) haps
    if (mod(HapCount,2)==0) then
        ShuffleInd1=ShuffleInd1+1

        ! Select the paternal haplotype if the individual it belongs
        ! too is genotyped and it is not the current individual
        if ((Shuffle1(ShuffleInd1)/=CurrentInd).and.(GlobalHmmHDInd(ShuffleInd1)==1)) then
            HapCount=HapCount+1
            SubH(HapCount,:)=FullH(Shuffle1(ShuffleInd1),:,1)
        endif
    else
        ShuffleInd2=ShuffleInd2+1

        ! Select the maternal haplotype if the individual it belongs
        ! too is genotyped and it is not the current individual
        if ((Shuffle2(ShuffleInd2)/=CurrentInd).and.(GlobalHmmHDInd(ShuffleInd2)==1)) then
            HapCount=HapCount+1
            SubH(HapCount,:)=FullH(Shuffle2(ShuffleInd2),:,2)
        endif
    endif
enddo

! The number of parameters of the HMM, N and M (Rabiner (1989) notation)
! calculated as in Li et al. (2010) are:
!   nHapInSubH = Number of haplotypes in the template haplotype set, H
!   nHapInSubH*nHapInSubH = Number of states, N = H^2
!   nSnpHmm = Number of Observations, M.
!
! ForwardProbs are the accumulated probabilities(??), and for
! ForwardPrbos(:,1) this array are the prior probabilities
!
! Allocate all possible state sequencies
!write(*,*) "Antes allocate"
!allocate(ForwardProbs(nHapInSubH*nHapInSubH,nSnpHmm))

!write(*,*) "Despues allocate"

! WARNING: I think this variable is treated in a wrong way. In MaCH
!          code, FordwardProbs is the array variable leftMatrices.
!          In there, every element of the array is another array with
!          s(s+1)/2 elements, as S(i,j)=S(j,i); where s=nHapInSubH.
!          The code should then be:
!
! states = nHapInSubH*(nHapInSubH+1)/2
! allocate(ForwardProbs(states,nSnpHmm))
!allocate(SubH(nHapInSubH,nSnpHmm))

call ForwardAlgorithm(CurrentInd)
call SampleChromosomes(CurrentInd)

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
!
! TODO: Check where to implement UpdateThetas, UpdateErrorRate and
!       TotalCrossovers subroutines. According to MaCH code, they
!       should go outside this subroutine and inside MaCHController.

! Cumulative genotype probability of through hmm processes
if (GlobalRoundHmm>HmmBurnInRound) then
    ProbImputeGenosHmm(CurrentInd,:)=ProbImputeGenosHmm(CurrentInd,:)&
        +FullH(CurrentInd,:,1)+FullH(CurrentInd,:,2)
endif

deallocate(ForwardProbs)
deallocate(SubH)
!print*, 'Number of thread:', OMP_get_thread_num(), allocated(SubH), allocated(ForwardProbs)


end subroutine MaCHForInd

!######################################################################
subroutine SampleChromosomes(CurrentInd)
use GlobalVariablesHmmMaCH
implicit none

integer,intent(in) :: CurrentInd

! Local variables
integer :: i,j,k,l,SuperJ,Index,OffOn,State1,State2,TmpJ,TopBot,FirstState,SecondState,Tmp
double precision :: Summer,ran1,Choice,Sum00,Sum01,Sum10,Sum11
double precision :: Probs(nHapInSubH*(nHapInSubH+1)/2)

Summer=0.0
Index=0
Probs = ForwardProbs(:,nSnpHmm)

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
Choice=ran1(idum)*Summer
Summer=0.0
Index=0
OffOn=0
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
    call ImputeAlleles(CurrentInd,SuperJ+1,State1,State2)
    TmpJ=SuperJ

    ! Cumulative recombination fraction allows us to skip over
    ! uninformative positions: Alleles with missing genotype are skipped
    ! but the recombination information (Thetas(SuperJ) is accumulated
    ! and used in the next location.
    Theta=Thetas(SuperJ)
    do while ((GenosHmmMaCH(CurrentInd,SuperJ)==3).and.SuperJ>1)
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
    Choice=ran1(idum)*Summer

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
    if (ran1(idum)>0.5) then
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

call ImputeAlleles(CurrentInd,1,State1,State2)

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
implicit none

integer,intent(in) :: CurrentInd,TopBot,FromMarker,ToMarker,FromState,ToState

! Local variables
integer :: i,State,FromMarkerLocal
double precision :: Recomb,Theta1,ran1

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
    Recomb=ran1(idum)*Theta

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
        State=int(ran1(idum)*nHapInSubH)+1
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
implicit none

integer,intent(in) :: CurrentInd,CurrentMarker,State1,State2

! Local variables
integer :: Imputed1,Imputed2,Genotype,Differences
double precision :: ran1

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
    if (ran1(idum)>=0.5) then
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
implicit none

integer, intent(in) :: CurrentInd

! Local variables
integer :: i,j,PrecedingMarker

! Setup the initial state distributions
call SetUpPrior

j=1
! CurrentInd is the individual being studied and it is necessary to
! obtain the genotype of this individual in ConditionaOnData subroutine
! For j=1, ConditionaOnData will initialize the variable
! ForwardProbs(:,1) with the Prior Probabilities
call ConditionOnData(CurrentInd,j)

! WARNING: This variable, Theta, should be considered as local as is
!          global through out the HMM code for different purposes.
!          Look at subroutines Transpose, SampleChromosomes and
!          SamplePath
Theta=0.0

PrecedingMarker=1
do j=2,nSnpHmm
    ! Cumulative recombination fraction allows us to skip uninformative positions
    Theta=Theta+Thetas(j-1)-Theta*Thetas(j-1)
    ! Skip over uninformative positions to save time
    if ((GenosHmmMaCH(CurrentInd,j)/=3).or.(j==nSnpHmm)) then
        call Transpose(j,PrecedingMarker)
        call ConditionOnData(CurrentInd,j)
        PrecedingMarker=j
        Theta=0.0
    endif
enddo

end subroutine ForwardAlgorithm

!######################################################################
subroutine Transpose(CurrentMarker,PrecedingMarker)
! Calculates the probability of get a particular state at CurrentMarker
! from any other state at PrecedingMarker using the transition probabilities.
! It basically calculates the first term of the forward variable.
use GlobalVariablesHmmMaCH
implicit none

integer,intent(in) :: CurrentMarker,PrecedingMarker

!Local variables
integer :: i,j,Index
double precision :: Summer,Marginals(nHapInSubH),NoChange,OneChange,TwoChange

if (Theta==0.0) then
    ForwardProbs(:,CurrentMarker)=ForwardProbs(:,PrecedingMarker)
    return
endif

Summer=0.0
Index=0
Marginals(:)=0.0

! Calculate the sum of all forward variables and the marginal of a
! given state.

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
if (Summer < 1e-15) then
    NoChange=NoChange*1e30
    OneChange=OneChange*1e30
    TwoChange=TwoChange*1e30
endif

! This final loop actually transposes the probabilities for each state,
! that is, calculates the final probabilities of getting a particular state.
Index=0
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
use GlobalVariablesHmmMaCH
implicit none

integer, intent(in) :: CurrentInd, Marker

! Local variables
integer :: i,j,Index
double precision :: Factors(0:1)

! We treat missing genotypes as uninformative about the mosaic's
! underlying state. If we were to allow for deletions and the like,
! that may no longer be true.
if (GenosHmmMaCH(CurrentInd,Marker)==3) then
    return
else
    ! Index keeps track of the states already visited. The total number
    ! of states in this chunk of code is (nHapInSubH x (nHapInSubH-1)/2)
    Index=0
    do i=1,nHapInSubH
        ! Probability to observe genotype SubH(i) being the true
        ! genotype GenosHmmMaCH in locus Marker
        Factors(0)=Penetrance(Marker,SubH(i,Marker),&
                    GenosHmmMaCH(CurrentInd,Marker))

        ! Probability to observe genotype SubH(i)+1 being the true
        ! genotype GenosHmmMaCH in locus Marker
        Factors(1)=Penetrance(Marker,SubH(i,Marker)+1,&
                    GenosHmmMaCH(CurrentInd,Marker))
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
allocate(Penetrance(nSnpHmm,0:2,0:2))

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
implicit none
integer :: i,j,state=0
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
subroutine SetUpEquations
! Initialize the variables and parameters of the HMM model described in
! Li et al. 2010, Appendix

use GlobalVariablesHmmMaCH
implicit none

integer :: i,j
double precision :: ran1

! ALLOCATE MEMORY
! The full Template Haplotype Library, H (Li et al. 2010, Appendix)
allocate(FullH(nIndHmmMaCH,nSnpHmm,2))

! Subset of H, Sub H (setup by the user)
! This is a way to restrict the size of the template library, as the
! complexity of the algorith increases cubically with the sample size
! (the cost of each update increases quadratically and the number of
!  updates increases linearly with sample size)
!allocate(SubH(nHapInSubH,nSnpHmm))

! HMM PARAMETERS
! nSnpHmm is the number of states in the HMM

! Vector of Combination of genotyping error (Li et al. 2010, Appendix)
! Epsilon is related with the Penetrance Matrix of the HMM which gives
! the emision probabilities for each state/marker/snp.
allocate(Epsilon(nSnpHmm))

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

! Initialization of HMM parameters
Epsilon=0.00000001
Thetas=0.01

!Initialise FullH
do i=1,nIndHmmMaCH      ! For every Genotyped Individual
    do j=1,nSnpHmm      ! For each SNP

        ! Phase homozygose locus
        if (GenosHmmMaCH(i,j)==0) then
            FullH(i,j,:)=0
        endif

        if (GenosHmmMaCH(i,j)==2) then
            FullH(i,j,:)=1
        endif

        ! WARNING: ran1 is a function that given a negative integer as
        !          parameter, returns a random number between [0.0,1.0]
        !          idum variable is set by the user in the Spec hmm param.
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
        if (GenosHmmMaCH(i,j)==3) then
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

call CalcPenetrance

ErrorUncertainty(:)=0
ErrorMatches(:)=0
ErrorMismatches(:)=0
Crossovers(:)=0
!print*, Crossovers(:)

end subroutine SetUpEquations

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
print*,Thetas(nSnpHmm-1)

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
        call SetPenetrance(i)
    endif
enddo

call UpdateError(matches, mismatches, uncertain, rate)

do i=1,nSnpHmm
    if (ErrorMismatches(i)<=2) call SetPenetrance(i)
enddo

end subroutine UpdateErrorRate

!######################################################################
subroutine SetPenetrance(marker)

use GlobalVariablesHmmMaCH
implicit none

integer,intent(in) :: marker

! Local variables
double precision :: Err

Err=Epsilon(marker)

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











































