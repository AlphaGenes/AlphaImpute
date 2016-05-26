module GlobalVariablesHmmMaCH
implicit none

integer, parameter :: GENOTYPE_MISSING=3
integer, parameter :: ALLELE_MISSING=3
integer, parameter :: READ_MISSING=0
integer            :: MISSING=3

double precision, parameter :: SEQUENCING_ERROR=0.01
double precision, parameter :: EPSILON_ERROR=0.01

integer, parameter :: NUM_SEGMENTS=1

character(len=300) :: GenotypeFileName,CheckPhaseFileName,CheckGenoFileName
character*(20), allocatable, dimension(:) :: AnimalsInbred,GlobalHmmMachID
integer :: nIndHmmMaCH,GlobalRoundHmm,nSnpHmm,nGametesPhased,nAnimPhased,nAnisInbred
integer :: nHapInSubH,useProcs,nRoundsHmm,HmmBurnInRound,idum,windowLength
real 	:: phasedThreshold,imputedThreshold
logical :: segmentOverlap
integer,allocatable,dimension(:,:) :: GenosHmmMaCH,SubH
integer(kind=1),allocatable,dimension(:,:,:) :: PhaseHmmMaCH,FullH
integer,allocatable,dimension(:) :: ErrorUncertainty,ErrorMatches,ErrorMismatches,Crossovers
integer,allocatable,dimension(:) :: GlobalHmmHDInd
logical,allocatable,dimension(:) :: GlobalInbredInd
logical,allocatable,dimension(:,:) :: GlobalHmmPhasedInd
double precision,allocatable,dimension(:) :: Thetas,Epsilon
double precision,allocatable,dimension(:,:) :: ForwardProbs
double precision,allocatable,dimension(:,:,:) :: Penetrance, ShotgunErrorMatrix
real,allocatable,dimension (:,:) :: ProbImputeGenosHmm
real,allocatable,dimension (:,:,:) :: ProbImputePhaseHmm
integer, allocatable :: GenosCounts(:,:,:)

!$omp threadprivate(ForwardProbs, SubH)

end module GlobalVariablesHmmMaCH
