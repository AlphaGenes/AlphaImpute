!#############################################################################################################################################################################################################################

module Global
implicit none


integer,parameter :: WindowsLinux=0     !If 1 then compile for Windows / If 0 then compile for Linux

integer,parameter :: TestVersion=0  !If 1 then this is a development version with intermediate checking, if 0 it is not

integer,parameter :: PicVersion=0   !If 1 then this is a PIC version with suitability for their system, if 0 it is not

integer,parameter :: SleepParameter=1!00

integer,parameter :: lengan=20,MissingGenotypeCode=9,OffspringFillMin=10
integer,parameter :: ImputeFromHDLibraryCountThresh=1,ImputeFromHDPhaseThresh=1
integer,parameter :: ImputeFromParentCountThresh=1,ImputeFromGrandParentCountThresh=1
real,parameter :: DisagreeThreshold=0.05,GeneProbThresh=0.99

character (len=300) :: PedigreeFile,PhasePath,TrueGenosFile,GenotypeFile,GenderFile

integer :: nAnisG,nAnisRawPedigree,nSnp,nAnisP,IntEditStat,nPhaseInternal,nPhaseExternal,OutOpt,SexOpt,HetGameticStatus,HomGameticStatus
integer :: nProcessors,nProcessGeneProb,nProcessAlphaPhase,ManagePhaseOn1Off0,CountRawGenos,InternalIterations,nAnisInGenderFile
integer :: nAgreeImputeHDLib,nAgreeParentPhaseElim,nAgreeInternalHapLibElim,MaxLeftRightSwitch,MinSpan,ConservativeHapLibImputation
integer :: TrueGenos1None0,nSnpRaw,nObsDataRaw,nAgreePhaseElim,nAgreeGrandParentPhaseElim,PreProcess,UseGP,BypassGeneProb,HMMOption
integer :: nSnpIterate,NoPhasing,AlphaPhasePresent,GeneProbPresent,PrePhased,UserDefinedHD,PedFreePhasing,PhaseTheDataOnly,RestartOption

real :: PercGenoForHD,PercSnpMiss,SecondPercGenoForHD,GenotypeErrorPhase,WellPhasedThresh

integer(kind=1),allocatable,dimension (:) :: SnpIncluded,RecIdHDIndex,GenderRaw,RecGender,IndivIsGenotyped
integer(kind=1),allocatable,dimension (:,:) :: Genos,TempGenos,TmpGenos,MSTermInfo
integer(kind=1),allocatable,dimension (:,:) :: ImputeGenos,SireDam
integer(kind=1),allocatable,dimension (:,:,:) :: ImputePhase,TmpPhase,GlobalWorkPhase
integer,allocatable :: RecPed(:,:),Setter(:),CoreAndTailLengths(:),CoreLengths(:),GpIndex(:,:),BaseAnimals(:),GlobalTmpCountInf(:,:)
integer,allocatable :: GlobalHmmID(:)
real,allocatable,dimension (:) :: Maf
real,allocatable,dimension (:,:) :: ProbImputeGenos
real,allocatable,dimension (:,:,:) :: ProbImputePhase
character*(lengan),allocatable :: GenotypeId(:),GenderId(:)

end module Global

!#############################################################################################################################################################################################################################


module GlobalPedigree
use Global
implicit none

real(kind=4),allocatable :: xnumrelmatHold(:)
integer :: NRMmem, shell, shellmax, shellWarning
integer,allocatable:: seqid(:),seqsire(:),seqdam(:),RecodeGenotypeId(:),passedorder(:)
character*(lengan),allocatable :: ped(:,:),Id(:),sire(:),dam(:)
integer :: GlobalExtraAnimals       !Change John Hickey

end module GlobalPedigree
