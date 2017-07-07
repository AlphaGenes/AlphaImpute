
module PARAMETERS
    ! TODO params only used in read in params, but where to put them? 
    integer, parameter :: OPT_RESTART_ALL=0
    integer, parameter :: OPT_RESTART_GENEPROB=1
    integer, parameter :: OPT_RESTART_PHASING=2
    integer, parameter :: OPT_RESTART_IMPUTATION=3

    ! TODO - add output so that final step can be run seperately
    integer, parameter :: OPT_RESTART_RECOMB=4

    integer, parameter :: RUN_HMM_NULL=0
    integer, parameter :: RUN_HMM_NO=1
    integer, parameter :: RUN_HMM_YES=2
    integer, parameter :: RUN_HMM_ONLY=3
    integer, parameter :: RUN_HMM_PREPHASE=4
    integer, parameter :: RUN_HMM_NGS=5

    integer, parameter :: MAX_READS_COUNT=100 ! Maximum number of reads for reference and alternative alleles
    integer,parameter :: TestVersion=0      !If 1 then this is a development version with intermediate checking, if 0 it is not

    logical,parameter :: PicVersion=.FALSE. !If 1 then this is a PIC version with suitability for their system, if 0 it is not

    integer,parameter :: SleepParameter=1!00

    integer,parameter :: lengan=20,OffspringFillMin=10
    integer,parameter :: ImputeFromHDLibraryCountThresh=1,ImputeFromHDPhaseThresh=1
    integer,parameter :: ImputeFromParentCountThresh=1,ImputeFromGrandParentCountThresh=1
    real,parameter :: DisagreeThreshold=0.05,GeneProbThresh=0.99

end module PARAMETERS
!#############################################################################################################################################################################################################################

module Global
    use PARAMETERS
    use iso_fortran_env
    use PedigreeModule
    use AlphaPhaseResultsModule
    implicit none



    integer :: CountRawGenos
    integer :: MaxLeftRightSwitch,MinSpan
    integer :: nSnpIterate,AlphaPhasePresent,GeneProbPresent

    integer(kind=1),allocatable,dimension (:) :: SnpIncluded
    integer(kind=1),allocatable,dimension (:,:) :: MSTermInfo
    integer(kind=1),allocatable,dimension (:,:,:) :: GlobalWorkPhase
    integer,allocatable :: Setter(:),GpIndex(:,:),GlobalTmpCountInf(:,:)
    real(kind=real64),allocatable,dimension (:) :: Maf
    real,allocatable,dimension (:,:) :: ProbImputeGenos
    real,allocatable,dimension (:,:,:) :: ProbImputePhase

    ! character*(lengan),allocatable :: GenotypeId(:),GenderId(:)

    integer, allocatable :: animChip(:)
    type(PedigreeHolder) :: ped !TODO move out of global

    type(AlphaPhaseResultsContainer) :: apResults
    real(kind=real64), allocatable :: GenosProbs(:,:,:) !< output of geneprob


    ! HMM OUTPUT STUFF

        real,allocatable,dimension (:,:) :: ProbImputeGenosHmm
        real,allocatable,dimension (:,:,:) :: ProbImputePhaseHmm
        integer, allocatable :: GenosCounts(:,:,:)
        integer :: GlobalRoundHmm
        integer(kind=1),allocatable,dimension(:,:,:) :: FullH !< full phase and imputation information from HMM (nanis, snps, haps[2])


end module Global

