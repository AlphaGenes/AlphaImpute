subroutine ForwardAlgorithmForHaplotype(CurrentInd)
use Global
use GlobalVariablesHmmMaCH
use Par_Zig_mod
use omp_lib

implicit none
integer, intent(in) :: CurrentInd

! Local variables
integer :: i

i = 1
call SetUpPriorHaplotype
call ConditionHaplotypeOnData(CurrentInd, i, SubH(i,Marker))

do i=2,nSnpHmm
	call GetSmallMemoryBlock
	call ConditionHaplotypeOnData(CurrentInd, i, SubH(i,Marker))
enddo

end subroutine ForwardAlgorithmForHaplotype

!######################################################################
subroutine SampleHaplotypeSource(CurrentInd)
use Global
use GlobalVariablesHmmMaCH
use Par_Zig_mod
use omp_lib

implicit none
integer,intent(IN) :: CurrentInd

! Local variables
integer :: state, marker, Thread

! double precision :: Probs(nHapInSubH*(nHapInSubH+1)/2)
double precision :: Probs(nHapInSubH)
double precision :: Summer, Choice, Theta, cross, nocross


Thread = omp_get_thread_num()
Summer=0.0
Probs = ForwardProbs(:,nSnpHmm)

do state=1,nHapInSubH
	Summer = Summer + Probs(state)
enddo

Summer=0.0
Choice = par_uni(Thread)*Summer

do haplotype=1,nHapInSubH
	Summer = Summer + Probs(haplotype)
	if (Summer >= Choice) exit
enddo

do marker=nSnpHmm-1,1,-1
	if (SubH(haplotype,marker)==SubH(nHapInSubH,marker)) then
		ErrorMatches(marker)=ErrorMatches(marker)+1
	else
		ErrorMismatches(marker)=ErrorMismatches(marker)+1
	endif

	Theta = Thetas(marker)
	Probs = ForwardProbs(:,marker)

	nocross = Probs(haplotype) * (1.0 - Theta)
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
	do haplotype=1,nHapInSubH
		Summer = Summer + Probs(haplotype)
		if (Summer >= Choice) exit
	enddo
enddo

if (SubH(haplotype,1)==SubH(haplotype,1)) then
	ErrorMatches(1)=ErrorMatches(marker)+1
else
	ErrorMismatches(1)=ErrorMismatches(marker)+1
endif

end subroutine SampleHaplotypeSource

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
	if (allele==factors(1)) then
		SubH(i,Marker) = factors(1)
	else
		SubH(i,Marker) = factors(0)
	endif
	ForwardProb(i,Marker) = ForwardProb(i,Marker) * SubH(i,Marker)
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
