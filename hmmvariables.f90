module GlobalVariablesHmmMaCH
implicit none

integer, parameter :: GENOTYPE_MISSING=3
integer, parameter :: ALLELE_MISSING=3
integer, parameter :: READ_MISSING=0
integer            :: MISSING=3

double precision, parameter :: SEQUENCING_ERROR=0.01
double precision, parameter :: EPSILON_ERROR=0.00000001

integer, parameter :: NUM_SEGMENTS=10

character(len=300) :: GenotypeFileName,CheckPhaseFileName,CheckGenoFileName
integer :: nIndHmmMaCH,GlobalRoundHmm,nSnpHmm,nGametesPhased,nAnimPhased
integer :: nHapInSubH,useProcs,nRoundsHmm,HmmBurnInRound,phasedThreshold,idum,windowLength
integer,allocatable,dimension(:,:) :: GenosHmmMaCH,SubH
integer(kind=1),allocatable,dimension(:,:,:) :: PhaseHmmMaCH,FullH
integer,allocatable,dimension(:) :: ErrorUncertainty,ErrorMatches,ErrorMismatches,Crossovers
integer,allocatable,dimension(:) :: GlobalHmmHDInd
logical,allocatable,dimension(:,:) :: GlobalHmmPhasedInd
double precision,allocatable,dimension(:) :: Thetas,Epsilon
double precision,allocatable,dimension(:,:) :: ForwardProbs
double precision,allocatable,dimension(:,:,:) :: Penetrance, ShotgunErrorMatrix
real,allocatable,dimension (:,:) :: ProbImputeGenosHmm
real,allocatable,dimension (:,:,:) :: ProbImputePhaseHmm
integer, allocatable :: GenosCounts(:,:,:)
!$omp threadprivate(ForwardProbs, SubH)

end module GlobalVariablesHmmMaCH
