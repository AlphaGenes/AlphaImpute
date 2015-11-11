subroutine ForwardAlgorithmForHaplotype(CurrentInd, hap)
use Global
use GlobalVariablesHmmMaCH
use Par_Zig_mod
use omp_lib

implicit none
integer, intent(in) :: CurrentInd, hap
! double precision, intent(IN), allocatable :: Thetas(:)

! Local variables
integer :: marker

#if DEBUG.EQ.1
    write(0,*) 'DEBUG: [ForwardAlgorithmForHaplotype]'
#endif

marker = 1
call SetUpPriorHaplotype
call ConditionHaplotypeOnData(CurrentInd, marker, PhaseHmmMaCH(CurrentInd,marker,hap))

do marker=2,nSnpHmm
    ! call GetSmallMemoryBlock
    call TransposeHaplotype(marker-1, marker, Thetas(marker-1))
    call ConditionHaplotypeOnData(CurrentInd, marker, PhaseHmmMaCH(CurrentInd,marker,hap))
enddo

end subroutine ForwardAlgorithmForHaplotype

!######################################################################
subroutine SampleHaplotypeSource(CurrentInd,hap)
use Global
use GlobalVariablesHmmMaCH
use Par_Zig_mod
use omp_lib

implicit none
integer,intent(IN) :: CurrentInd, hap

! Local variables
integer :: i, state, marker, Thread

! double precision :: Probs(nHapInSubH*(nHapInSubH+1)/2)
double precision :: Probs(nHapInSubH)
double precision :: Summer, Choice, Theta, cross, nocross

#if DEBUG.EQ.1
    write(0,*) 'DEBUG: [SampleHaplotypeSource]'
#endif

Thread = omp_get_thread_num()
Summer=0.0
Probs = ForwardProbs(:,nSnpHmm)

! Calculate sum over all states
do state=1,nHapInSubH
    Summer = Summer + Probs(state)
enddo

! Sample number and select state
Choice = par_uni(Thread)*Summer
Summer=0.0

do i=1,nHapInSubH
    Summer = Summer + Probs(i)
    if (Summer >= Choice) exit
enddo

do marker=nSnpHmm-1,1,-1
    ! Track whether imputed state matches observed allele
    if (SubH(i,marker)==PhaseHmmMaCH(CurrentInd,marker,hap)) then
        ErrorMatches(marker)=ErrorMatches(marker)+1
    else
        ErrorMismatches(marker)=ErrorMismatches(marker)+1
    endif

    ! Impute if allele is missing
    if (PhaseHmmMaCH(CurrentInd,marker,hap)==3) then
        FullH(CurrentInd,marker,hap) = SubH(i,marker)
        ! call ImputeAllele(CurrentInd,marker,i,hap)
    endif

    Theta = Thetas(marker)
    Probs = ForwardProbs(:,marker)

    nocross = Probs(i) * (1.0 - Theta)
    Summer = 0.0

    do i=1,nHapInSubH
        Summer = Summer + Probs(i)
    enddo

    cross = Summer * Theta / nHapInSubH

    ! Sample number and decide how many state changes occurred between the
    ! two positions
    Choice = par_uni(Thread)*(nocross+cross)

    ! The most likely outcome is that no changes occur ...
    if (Choice <= nocross) continue

    ! TODO: Look what crossovers are
    crossovers(i)= crossovers(i)+1

    ! If a crossover occured, we need to sample a state according to probability
    Choice = par_uni(Thread)*(Summer)

    Summer = 0.0
    do i=1,nHapInSubH
        Summer = Summer + Probs(i)
        if (Summer >= Choice) exit
    enddo
enddo

! Track whether imputed state matches observed allele
if (SubH(i,1)==PhaseHmmMaCH(CurrentInd,1,hap)) then
    ErrorMatches(1)=ErrorMatches(1)+1
else
    ErrorMismatches(1)=ErrorMismatches(1)+1
endif

! Impute if allele is missing
if (PhaseHmmMaCH(CurrentInd,1,hap)==3) then
    FullH(CurrentInd,1,hap) = SubH(i,1)
    ! ImputeAllele(CurrentInd,1,i,hap)
endif

end subroutine SampleHaplotypeSource

!######################################################################
subroutine TransposeHaplotype(PrecedingMarker, CurrentMarker, Theta)
! Calculates the probability of get a particular state at CurrentMarker
! from any other state at PrecedingMarker using the transition probabilities.
use GlobalVariablesHmmMaCH

integer, intent(IN) :: PrecedingMarker, CurrentMarker
double precision, intent(IN) :: Theta

! Local variables
integer :: i
double precision :: Summer, NoChange, OneChange

if (Theta==0.0) then
    do i=1,nHapInSubH
        ForwardProbs(i,CurrentMarker)=ForwardProbs(i,PrecedingMarker)
    enddo
    return
endif

! Summer=0.0
Summer=sum(ForwardProbs(:,PrecedingMarker))
NoChange=1.0-Theta
OneChange=Summer*Theta/nHapInSubH

! Automatically rescale likelihoods when they get too small
if (Summer < 1e-15) then
    NoChange=NoChange*1e30
    OneChange=OneChange*1e30
endif

! This final loop actually transposes the probabilities for each state
do i=1,nHapInSubH
    ForwardProbs(i,CurrentMarker)=ForwardProbs(i,PrecedingMarker)*NoChange+OneChange
enddo

end subroutine TransposeHaplotype

!######################################################################
subroutine ConditionHaplotypeOnData(Marker, allele)
use Global
use GlobalVariablesHmmMaCH
use Par_Zig_mod
use omp_lib

implicit none
integer, intent(IN) :: Marker, allele

! Local variables
integer :: i
double precision :: factors(0:1), ErrorRate


call GetErrorRatebyMarker(Marker, ErrorRate)
factors(0) = ErrorRate
factors(1) = 1.0 - factors(0)

do i=1,nHapInSubH
    if (allele==SubH(i,Marker)) then
        ForwardProbs(i,Marker) = ForwardProbs(i,Marker) * factors(1)
    else
        ForwardProbs(i,Marker) = ForwardProbs(i,Marker) * factors(0)
    endif
enddo

end subroutine ConditionHaplotypeOnData

!######################################################################
subroutine SetUpPriorHaplotype
! Set up de initial state distribution, that is, the probability that
! the sequence starts with the state Sj:
!   PIj = P(t1=Sj) = ForwardProb(j,1)

use GlobalVariablesHmmMaCH
use omp_lib

implicit none
integer :: i
double precision :: prior

! Initially, every state is equally possible
prior=1.0/nHapInSubH
do i=1,nHapInSubH
   ForwardProbs(i,1)=prior
enddo

end subroutine SetUpPriorHaplotype
