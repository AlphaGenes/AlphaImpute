!-----------------------------------------------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-----------------------------------------------------------------------------------------------------------------------
!
! MODULE: AlphaImputeSpecFileModule
!
!> @file        AlphaImputeSpecFileModule.f90
!
! DESCRIPTION:
!> @brief       Module holding input parameters
!>
!> @details     This MODULE contains a class which contains all input parameters read in from a spec file.
!> It also contains the default container object for the spec file, defaultInput.
!
!> @author      David Wilson, david.wilson@roslin.ed.ac.uk
!
!> @date        Nov 07, 2016
!
!> @version     0.0.1 (alpha)
!
! REVISION HISTORY:
! 2016.11.07  DWilson - Initial Version
!
!-----------------------------------------------------------------------------------------------------------------------




module AlphaImputeSpecFileModule
	use iso_fortran_env
	use ConstantModule
	use baseSpecFileModule

	implicit none

	type, extends(baseSpecFile) ::  AlphaImputeInput
	! box 1
	character(len=300)::TrueGenotypeFile="None",GenderFile="None",InbredAnimalsFile="None", HapListFile="None",animalPhaseFile="None"
	character(len=300) :: PriorAllFreqsFile="None"
	integer(kind=1) :: TrueGenos1None0

	! box 3
	integer(kind=int32) :: MultiHD
	integer(kind=int32), allocatable :: nSnpByChip(:)
	real(real32) :: PercGenoForHD

	! box 4
	integer(kind=1) :: IntEditStat,OutOpt
	real(real32) :: PercSnpMiss,SecondPercGenoForHD

	! box 5
	integer(kind=1) :: NoPhasing,ManagePhaseOn1Off0,PedFreePhasing
	integer(kind=int32) ::  nPhaseExternal,nPhaseInternal
	character(len=:), allocatable :: PhasePath
	integer(kind=int32), allocatable :: CoreAndTailLengths(:),CoreLengths(:)
	real(kind=real32) :: GenotypeErrorPhase
	logical :: largeDatasets
	integer(kind=int32) :: PhaseSubsetSize, PhaseNIterations
	character(len=20) :: iterateMethod
	integer :: minoverlaphaplotype

	logical :: outputonlygenotypedanimals
	! box 6
	integer(kind=int32) :: InternalIterations
	integer(kind=1) :: ConservativeHapLibImputation
	real(kind=real32) :: WellPhasedThresh

	! box 7
	! idum is seed
	integer(kind=int32) :: idum
	integer(kind=int32) :: nHapInSubH,nRoundsHmm,HmmBurnInRound
	real(kind=real32) ::  imputedThreshold
	logical :: HapList=.FALSE.
	logical :: PriorAllFreqs=.FALSE.
	integer(kind=1) :: HMMOption
	integer :: HapListUnit
	integer :: PriorAllFreqsUnit
	real(kind=real32) :: phasedThreshold                   !< threshold of phase information accept


	! box 8

	logical :: PreProcess
	integer(kind=1) :: PhaseTheDataOnly

	integer(kind=1) :: UserDefinedHD,PrePhased,RestartOption
	integer :: prePhasedFileUnit
	! other
	integer(kind=int32) :: nSnpRaw,nAgreeInternalHapLibElim
	integer(kind=int32) :: useProcs
	logical :: cluster

	integer :: alphaphaseoutput

	logical :: useFerdosi

	logical :: modelrecomb

	logical :: phaseRunsSet !< have phasing runs been set?

	contains
		procedure :: ReadInParameterFile
		final :: destroyAlphaImputeInput
		! procedure :: validateAlphaImputeSpec
		procedure   :: validate => validateAlphaImputeSpec
		procedure :: writeSpec => writeOutSpecOptions
		procedure :: copy => copyAlphaImputeSpecFile
	end type AlphaImputeInput

	type(AlphaImputeInput),pointer :: defaultInput

	contains


		subroutine destroyAlphaImputeInput(in)

			type(AlphaImputeInput), intent(inout) :: in

			if (allocated(in%PhasePath)) then
				deallocate(in%PhasePath)
			endif

			if (allocated(in%coreLengths)) then
				deallocate(in%CoreAndTailLengths)
				deallocate(in%CoreLengths)
			endif

		end subroutine destroyAlphaImputeInput


		subroutine copyAlphaImputeSpecFile(old, new)

			class(AlphaImputeInput),intent(in) ::  old
			class(baseSpecFile),allocatable,intent(out) ::  new



			allocate(AlphaImputeInput :: new)

			select type(new)

			type is(AlphaImputeInput)
			new%resultFolderPath = old%resultFolderPath
			new%plinkinputfile = old%plinkinputfile
			new%PedigreeFile = old%PedigreeFile
			new%GenotypeFile = old%GenotypeFile
			new%plinkBinary = old%plinkBinary
			new%nsnp = old%nsnp
			new%SexOpt = old%SexOpt
			new%HetGameticStatus = old%HetGameticStatus
			new%HomGameticStatus = old%HomGameticStatus
			new%plinkOutput = old%plinkOutput
			new%stopAfterPlink = old%stopAfterPlink

			if (allocated(old%useChroms)) then
				new%useChroms = old%useChroms
			endif
			new%programName = old%programName
			new%version = old%version
			new%GenderFile = old%GenderFile
			new%InbredAnimalsFile = old%InbredAnimalsFile
			new%HapListFile = old%HapListFile
			new%animalPhaseFile = old%animalPhaseFile
			new%TrueGenotypeFile = old%TrueGenotypeFile
			new%PriorAllFreqsFile = old%PriorAllFreqsFile
			new%TrueGenos1None0 = old%TrueGenos1None0
			new%MultiHD = old%MultiHD

			if (allocated(old%nSnpByChip)) then
				new%nSnpByChip = old%nSnpByChip
			endif
			new%PercGenoForHD = old%PercGenoForHD
			new%IntEditStat = old%IntEditStat
			new%OutOpt = old%OutOpt
			new%PercSnpMiss = old%PercSnpMiss
			new%SecondPercGenoForHD = old%SecondPercGenoForHD
			new%NoPhasing = old%NoPhasing
			new%ManagePhaseOn1Off0 = old%ManagePhaseOn1Off0
			new%PedFreePhasing = old%PedFreePhasing
			new%nPhaseExternal = old%nPhaseExternal
			new%nPhaseInternal = old%nPhaseInternal

			if (allocated(old%phasePath)) then
				new%PhasePath = old%PhasePath
			endif

			if (allocated(old%CoreAndTailLengths)) then
				new%CoreAndTailLengths = old%CoreAndTailLengths
			endif

			if (allocated(old%coreLengths)) then
				new%CoreLengths = old%CoreLengths
			endif
			new%GenotypeErrorPhase = old%GenotypeErrorPhase
			new%largeDatasets = old%largeDatasets
			new%PhaseSubsetSize = old%PhaseSubsetSize
			new%PhaseNIterations = old%PhaseNIterations
			new%iterateMethod = old%iterateMethod
			new%minoverlaphaplotype = old%minoverlaphaplotype
			new%outputonlygenotypedanimals = old%outputonlygenotypedanimals
			new%InternalIterations = old%InternalIterations
			new%ConservativeHapLibImputation = old%ConservativeHapLibImputation
			new%WellPhasedThresh = old%WellPhasedThresh
			new%idum = old%idum
			new%nHapInSubH = old%nHapInSubH
			new%nRoundsHmm = old%nRoundsHmm
			new%HmmBurnInRound = old%HmmBurnInRound
			new%imputedThreshold = old%imputedThreshold
			new%HapList = old%HapList
			new%PriorAllFreqs = old%PriorAllFreqs
			new%HMMOption = old%HMMOption
			new%HapListUnit = old%HapListUnit
			new%PriorAllFreqsUnit = old%PriorAllFreqsUnit
			new%phasedThresholD = old%phasedThresholD
			new%PreProcess = old%PreProcess
			new%PhaseTheDataOnly = old%PhaseTheDataOnly
			new%UserDefinedHD = old%UserDefinedHD
			new%PrePhased = old%PrePhased
			new%RestartOption = old%RestartOption
			new%prePhasedFileUnit = old%prePhasedFileUnit
			new%nSnpRaw = old%nSnpRaw
			new%nAgreeInternalHapLibElim = old%nAgreeInternalHapLibElim
			new%useProcs = old%useProcs
			new%cluster = old%cluster
			new%alphaphaseoutput = old%alphaphaseoutput
			new%useFerdosi = old%useFerdosi
			new%modelrecomb = old%modelrecomb
			new%phaseRunsSet = old%phaseRunsSet



			class default
			write(error_unit, *) "ERROR: AlphaImpute given incorrect object type as input"
			call abort()

		end select

	end subroutine copyAlphaImputeSpecFile
	!---------------------------------------------------------------------------
	! DESCRIPTION:
	!> @brief      Constructor for AlphaImputeInput object
	!
	!> @details    Initialise new Bit Haplotype
	!
	!> @author     David Wilson, david.wilson@roslin.ed.ac.uk
	!
	!> @date       Nov 07, 2016
	!
	! PARAMETERS:
	!> @param[in]  specfile - The path of the specfile
	!---------------------------------------------------------------------------

	subroutine ReadInParameterFile(this,SpecFile)
		use AlphaHouseMod, only: parseToFirstWhitespace,splitLineIntoTwoParts,toLower
		use PARAMETERS
		use omp_lib

		class(AlphaImputeInput), intent(inout),target :: this
		integer :: unit,IOStatus,MultipleHDpanels,i,stat
		character(len=*), intent(in) :: SpecFile

		character(len=300) :: first, line
		character(len=:), allocatable::tag, tmptag
		character(len=300),dimension(:),allocatable :: second
		LOGICAL :: exists


		this%programName = "AlphaImpute"
		this%alphaphaseoutput = 0
		this%useFerdosi = .false.
		this%MultiHD = 0
		! this%nsnp= LARGESNUMBER
		this%minoverlaphaplotype = 0
		this%PreProcess = .false.
		this%cluster = .false.
		this%iterateMethod = "Off"
		this%PhaseNIterations = 1
		this%resultFolderPath = "Results"
		this%modelrecomb = .true.
		this%outputonlygenotypedanimals = .false.
		this%TrueGenos1None0=0
		this%hmmoption=RUN_HMM_NO
		this%nHapInSubH = 0
		this%HmmBurnInRound = 0
		this%nRoundsHMM = 1
		this%idum = -123456789
		this%phasedThreshold = 90.0
		this%imputedThreshold = 50.0
		this%plinkinputfile = ""
		this%nsnp = 0
		this%PercGenoForHD=90.0
		this%inteditstat = 0
		this%InternalIterations = 5
		this%outopt = 1
		this%RestartOption = 0
		this%managephaseon1off0 = 1
		this%PedFreePhasing = 0
		this%PhaseTheDataOnly=0
		this%noPhasing = 1
		this%useProcs = OMP_get_num_procs()
		this%phaseRunsSet = .false.

		this%nPhaseExternal = this%useProcs/2
		this%nPhaseInternal = this%useProcs
		this%plinkoutput = .false.
		this%stopAfterPlink = .false.

		inquire(FILE=SpecFile, EXIST=exists)


		if (exists) then

			open(newunit=unit, file=SpecFile, action="read", status="old")
			IOStatus = 0

			READFILE: do
				read(unit,"(A)", IOStat=IOStatus)  line


				if (IOStatus /= 0) exit


				if (len_trim(line)==0) then
					CYCLE
				end if

				call splitLineIntoTwoParts(trim(line), first, second)
				tag = parseToFirstWhitespace(first)
				tmptag = trim(tag)
				if (first(1:1)=="=" .or. first(1:1) == DEFAULTCOMMENT .or. len(trim(line))==0) then
					cycle
				else
					select case(tmptag)

						!runs full chromosome
					case("plinkinputfile")
						if (.not. allocated(second)) then
							write(error_unit, "(A)") "error, Plinkinputfile allocated incorrectly"
						else
							if (size(second) < 2) then
								write(error_unit, "(A)") "error, Plinkinputfile allocated incorrectly"
							else
								if (tolower(second(1)) == "binary") then
									this%plinkBinary = .true.
								else
									this%plinkBinary = .false.
								endif

								write(this%plinkinputfile, "(A)") second(2)
							endif

						end if
					case("plinkoutput")
						if (toLower(trim(second(1)))=="yes") then
							this%plinkoutput = .true.
						else if (toLower(trim(second(1)))=="no") then
							this%plinkoutput = .false.
						else
							write(error_unit, *) "WARNING - plinkoutput set incorrectly"
						endif

						! box 1 inputs
					case("pedigreefile")
						if (.not. allocated(second)) then
							write(*, "(A,A)") "No pedigree file specified. Using default filename: ", this%PedigreeFile
						else
							write(this%PedigreeFile, "(A)") second(1)
						end if
					case("genotypefile")
						if (.not. allocated(second)) then
							write(*, "(A,A)") "No genotype file specified. Using default filename: ", this%Genotypefile
						else
							write(this%Genotypefile, "(A)") second(1)
						endif
					case("truegenotypefile")
						if (.not. allocated(second)) then
							write(*, "(A,A)") "No true genotype file specified. Using default filename: ", this%TrueGenotypeFile
						else
							write(this%TrueGenotypeFile, "(A)") second(1)
							if (trim(toLower(this%TrueGenotypeFile))=="none") then
								this%TrueGenos1None0=0
							else
								this%TrueGenos1None0=1
							endif
						endif

						! box 2 inputs
					case("sexchrom")
						this%SexOpt=9
						this%HetGameticStatus=9
						this%HomGameticStatus=9
						if (toLower(trim(second(1)))=="yes") then
							this%genderFile = second(2)
							this%HetGameticStatus=9
							if (toLower(trim(second(3)))=="male") then        ! Species  with heterogametic males
								this%HetGameticStatus=1                              ! My father is heterogametic
								this%HomGameticStatus=2                              ! My mother is homogametic
							endif
							if (toLower(trim(second(3)))=="female") then      ! Species with heterogametic females
								this%HetGameticStatus=2                              ! My mother is heterogametic
								this%HomGameticStatus=1                              ! My father is homogametic
							endif
							if (this%HetGameticStatus==9) then
								print*, "Warning - heterogametic status is misspecified"
								stop
							endif
							this%SexOpt=1
						endif

						! Not sex chrom
						if (tolower(trim(second(1)))=="no") then
							this%SexOpt=0
						endif
						if (this%SexOpt==9) then
							print*, "Warning - Sex chromosome status is misspecified"
							stop
						endif

						! box 3 inputs
					case("nsnp")
						read(second(1),*,iostat=stat) this%nsnp
						if (stat /= 0) then

							print*, "Error: nsnp specified incorrectly"
							stop 3001
						endif

					case("numbersnp")
						read(second(1),*,iostat=stat) this%nsnp
						if (stat /= 0) then

							print*, "Error: numbersnp specified incorrectly"
							stop 3001
						endif

					case("multiplehdpanels")
						! Get the information of Multiple HD chips
						! MultipleHDpanels
						if (trim(toLower(second(1))) /= "no") then
							read(second(1),*,iostat=stat) MultipleHDpanels

							if (stat /= 0) then
								print*, "Error: hdanimalsthreshold specified incorrectly"
								stop
							endif
							if (MultipleHDpanels/=0) this%MultiHD=MultipleHDpanels

							if (allocated(this%nSnpByChip)) then
								deallocate(this%nSnpByChip)
							endif
							allocate(this%nSnpByChip(this%MultiHD))
						endif
					case("numbersnpxchip")
						do i=1,this%MultiHD
							read(second(i),*,iostat=stat) this%nSnpByChip(i)

							if (stat /= 0) then
								print*, "Error: hdanimalsthreshold specified incorrectly"
								stop
							endif
						enddo

					case("hdanimalsthreshold")
						read(second(1),*,iostat=stat) this%PercGenoForHD
						if (stat /= 0) then
							print*, "Error: hdanimalsthreshold specified incorrectly"
							stop
						endif

						! box 4
					case("internaledit")
						if (ToLower(trim(second(1))) == "yes") then
							this%inteditstat = 1
						else if (ToLower(trim(second(1))) == "no") then
							this%inteditstat = 0
						else
							write(error_unit,*) "error: internaledit is incorrectly defined in the spec file"
							stop
						endif
						if (this%IntEditStat==1 .AND. this%MultiHD/=0) then
							write(error_unit,*) "IntEditStat and MultipleHDpanels are incompatible,"
							write(error_unit,*) "Please, considere to use only one HD panel or to disable internal editing"
							stop 4001
						endif

					case("editingparameters")
						this%outopt=9
						if (this%IntEditStat==1) then
							if (size(second)<4) then
								goto 4000
							endif
							read(second(1),*) this%PercGenoForHD
							read(second(2),*) this%PercSnpMiss
							read(second(3),*) this%SecondPercGenoForHD

							if (toLower(second(4))=="allsnpout") this%outopt=1
							if (toLower(second(4))=="editedsnpout") this%outopt=0
							if (this%outopt==9) then
								goto 4000
							endif
						else
							! In case no editing is set and there is a single HD panel, a threshold to determine HD individuals is needed
							if (size(second) == 4) then

								read(second(1),*) this%PercGenoForHD
								read(second(2),*) this%PercSnpMiss
								read(second(3),*) this%SecondPercGenoForHD



								if (toLower(second(4))=="allsnpout") this%outopt=1
								if (toLower(second(4))=="editedsnpout") this%outopt=0
							else
								if (this%MultiHD==0) this%PercGenoForHD=90.0
							endif
							this%outopt=1
						endif

						cycle
						4000 print*, "Output options incorrectly specified"
						stop

					case("outputonlygenotypedanimals")
						if (ToLower(trim(second(1))) == "yes") then
							this%outputonlygenotypedanimals = .true.
						else if (ToLower(trim(second(1))) == "no") then
							this%outputonlygenotypedanimals = .false.
						endif

						! box 5
					case("numberphasingruns")
						this%noPhasing = 1
						this%phaseRunsSet = .true.
						if (ToLower(trim(second(1))) == "phasedone") then  !phasedone,path,nphaseruns
							if (size(second) /=3) then
								goto 4051
							endif
							this%managephaseon1off0=0
							if (allocated(this%phasePath)) then
								deallocate(this%phasePath)
							endif
							allocate(character(len(second(2))) :: this%phasePath)

							this%phasePath = second(2)
							read(second(3),*) this%nPhaseInternal
							this%nPhaseExternal = this%nPhaseInternal/2
						else if(ToLower(trim(second(1))) == "nophase") then
							this%noPhasing = 0
							this%managephaseon1off0 = 0
						else
							this%managephaseon1off0 = 1
							read(second(1),*) this%nPhaseExternal
							this%nPhaseInternal = 2*this%nPhaseExternal

							if(this%nPhaseExternal <1 ) then
								write(error_unit,*) "Error: Too few phasing runs requested. The minimum this program supports is 1, 4 are reccomended."
								stop 40512
							endif
						endif


						if (allocated(this%CoreAndTailLengths)) then
							deallocate(this%CoreAndTailLengths)
						endif
						if (allocated(this%CoreLengths)) then
							deallocate(this%CoreLengths)
						endif



						cycle
						4051 write(error_unit,*) "NumberPhasingRuns has been set incorrectly"
						stop 4051

					case("coreandtaillengths")
						if (this%managephaseon1off0 == 0) then
							cycle
						endif
						if (this%phaseRunsSet == .false.) then

							this%nPhaseExternal = size(second)
							this%phaseRunsSet = .true.
						endif

						allocate(this%CoreAndTailLengths(this%nPhaseExternal))


						do i=1,this%nPhaseExternal
							read(second(i), *) this%CoreAndTailLengths(i)

							if (this%nsnp /= 0 ) then
								if (this%CoreAndTailLengths(i) > this%nsnp) then

									write(error_unit, *) "Error: core and Tail lengths is given a number than largest number of snps specified in nsnps"
									write(error_unit, *) this%CoreAndTailLengths(i) ," vs ", this%nsnp
									stop 40523
								endif
							endif
						enddo

					case("corelengths")

						if (this%managephaseon1off0 == 0) then
							cycle
						endif

						if (this%phaseRunsSet == .false.) then
							this%nPhaseExternal = size(second)
							this%phaseRunsSet = .true.
						endif
						allocate(this%coreLengths(this%nPhaseExternal))
						if (size(second) /= size(this%corelengths)) then
							write(error_unit,*) "Error: numberofphasingruns is set to a different number of parameters than what is specified here, core number of core and tail lengths are different \n Please set this to the same number of parameters that are given for CoreAndTailLengths and CoreLengths"
							stop 40532
						endif
						do i=1,size(second)
							read(second(i),*,iostat=stat) this%corelengths(i)
							if (stat /= 0) then
								print*, "Error: hdanimalsthreshold specified incorrectly"
								stop
							endif
						enddo

					case("pedigreefreephasing")
						if (this%nPhaseExternal /= 0) then
							if (ToLower(trim(second(1))) == "no") then
								this%PedFreePhasing= 0
							elseif (ToLower(trim(second(1))) == "yes") then
								this%PedFreePhasing = 1
							else
								write(error_unit,*) "Error: PedFreePhasing has been set incorrectly."
								stop 4054
							endif
						endif

					case("genotypeerror")
						if (this%nPhaseExternal /= 0) then
							read(second(1), *, iostat=stat) this%GenotypeErrorPhase
							if (stat /= 0) then
								print*, "Error: hdanimalsthreshold specified incorrectly"
								stop
							endif
						endif

					case("numberofprocessorsavailable")
						read(second(1),*) this%useProcs
						write(error_unit,*) "WARNING: numberofprocessorsavailable is legacy and will be removed in future versions. Please use option ParallelProcessors instead"
						if (this%useProcs > OMP_get_num_procs()) then
							write(error_unit,*) "WARNING - more processors than are available are specified under numberofprocessorsavailable"
							write(error_unit,*) this%useProcs, " specified, ", OMP_get_num_procs(), " available."
						endif

					case("largedatasets")
						if (ToLower(trim(second(1)))== "yes") then
							this%largedatasets=.true.
							read(second(2),*) this%PhaseSubsetSize
							read(second(3),*) this%PhaseNIterations
							if (size(second) < 4 ) then
								this%iterateMethod  = "RandomOrder"
							else
								if (ToLower(trim(second(4)))== "off") then
									this%iterateMethod  = "Off"
								else if (ToLower(trim(second(4)))== "randomorder") then
									this%iterateMethod  = "RandomOrder"
								else if (ToLower(trim(second(4)))== "inputorder") then
									this%iterateMethod  = "InputOrder"
								else
									this%iterateMethod= "Off"

								endif
							endif

						else
							this%largedatasets=.false.

						endif

					case("minoverlaphaplotype")
						read(second(1),*) this%minoverlaphaplotype
						if (this%minoverlaphaplotype < 0) then
							write(error_unit,*) "ERROR: Min minoverlap haplotype size is set incorrectly!"
						endif


					case("alphaphaseoutput")
						if (ToLower(trim(second(1))) == "no") then
							this%alphaphaseoutput= 0
						elseif (ToLower(trim(second(1))) == "yes") then
							this%alphaphaseoutput = 1
						elseif (ToLower(trim(second(1))) == "binary") then
							this%alphaphaseoutput = 2
						elseif (ToLower(trim(second(1))) == "verbose") then
							this%alphaphaseoutput = 3
						else
							write(error_unit,*) "Error: alphaphaseoutput has been set incorrectly."
							stop 4054
						endif

						! box 6
					case("internaliterations")
						read(second(1), *) this%InternalIterations

					case("conservativehaplotypelibraryuse")
						if(ToLower(trim(second(1))) == "no") then
							this%ConservativeHapLibImputation = 0
						else if (ToLower(trim(second(1))) == "yes") then
							this%ConservativeHapLibImputation = 1
						else
							write (error_unit, *) "ConservativeHaplotypeLibraryUse not correctly set"
							stop
						endif

					case("wellphasedthreshold")
						read(second(1),*) this%WellPhasedThresh

						! box 7
					case("hmmoption")
						this%hmmoption=RUN_HMM_NULL
						if (toLower(trim(second(1)))=='no') this%hmmoption=RUN_HMM_NO
						if (toLower(trim(second(1)))=='yes') this%hmmoption=RUN_HMM_YES
						if (toLower(trim(second(1)))=='only') this%hmmoption=RUN_HMM_ONLY
						if (toLower(trim(second(1)))=='prephase') this%hmmoption=RUN_HMM_PREPHASE
						if (toLower(trim(second(1)))=="ngs") this%hmmoption=RUN_HMM_NGS
						if (this%hmmoption==RUN_HMM_NULL) then
							write(error_unit,*), "this%hmmoption not correctly specified"
							stop
						endif

						! case("hmmparameters")
						! 	if (.not. allocated(second)) then
						! 		write(error_unit,*) "templatehaplotypes not set correctly"
						! 		stop
						! 	endif

						! 	if (size(second) /= 5) then
						! 		write(error_unit,*) "hmmparameters not set correctly"
						! 		stop
						! 	endif

						! 	read(second(1), *) this%nHapInSubH
						! 	read(second(2), *) this%HmmBurnInRound
						! 	read(second(3), *) this%nRoundsHMM
						! 	read(second(4), *) this%useProcs
						! 	if (this%useProcs > OMP_get_num_procs()) then
						! 		write(error_unit,*) "WARNING - more processors than are available are specified under parallelprocessors"
						! 		write(error_unit,*) this%useProcs, " set vs ",OMP_get_num_procs()," available"
						! 	endif
						! 	read(second(5), * ) this%idum

					case("templatehaplotypes")
						if (.not. allocated(second)) then
							write(error_unit,*) "templatehaplotypes not set correctly"
							stop
						endif
						read(second(1), *,iostat=stat) this%nHapInSubH
						if (stat /= 0) then
							print *,"templatehaplotypes set incorrectly"
							stop
						endif

					case("burninrounds")
						read(second(1), *,iostat=stat) this%HmmBurnInRound
						if (stat /= 0) then
							print *,"burninrounds set incorrectly"
							stop
						endif
					case("rounds")
						read(second(1), *,iostat=stat)this%nRoundsHMM
						if (stat /= 0) then
							print *,"rounds set incorrectly"
							stop
						endif

					case("seed")
						read(second(1), *,iostat=stat)this%idum
						if (stat /= 0) then
							print *,"seed set incorrectly"
							stop
						endif

					case("thresholdformissingalleles")
						read(second(1), *,iostat=stat) this%phasedThreshold
						if (stat /= 0) then
							write(error_unit,*) "ERROR: thresholdformissingalleles set incorrectly"
							stop
						endif

					case("phasedanimalsthreshold")
						read(second(1), *,iostat=stat) this%phasedThreshold
						if (stat /= 0) then
							write(error_unit,*) "ERROR: phasedanimalsthreshold set incorrectly"
							stop
						endif

					case("thresholdimputed")
						read(second(1), *,iostat=stat) this%imputedThreshold
						if (stat /= 0) then
							write(error_unit,*) "ERROR: ThresholdImputed set incorrectly"
							stop
						endif

					case("wellimputedthreshold")
						read(second(1), *,iostat=stat) this%imputedThreshold
						if (stat /= 0) then
							write(error_unit,*) "wellimputedthreshold set incorrectly"
							stop
						endif
					case("haplotypeslist")
						if (.not. allocated(second)) then
							write(*, "(A,A)") "WARNING: No list of haploytpes specified"
						else
							if (trim(second(1)) /= "None") then
								this%HapList = .TRUE.
								this%HapListFile = trim(second(1))
							endif
						endif

					case("priorallelefrequencies")
						if (.not. allocated(second)) then
							write(*, "(A,A)") "WARNING: No list of haploytpes specified"
						else
							if (trim(second(1)) /= "None") then
								this%PriorAllFreqs = .TRUE.
								this%PriorAllFreqsFile = trim(second(1))
							endif
						endif

						!  box 8
					case("preprocessdataonly")
						if (ToLower(second(1))=="no") then
							this%PreProcess=.FALSE.
						else
							if (ToLower(second(1))=="yes") then
								this%PreProcess=.TRUE.
							else
								write(error_unit,*) "Stop - Preprocess of data option incorrectly specified"
								stop
							endif
						endif

					case("phasingonly")
						if (toLower(trim(second(1)))=="no") then
							this%PhaseTheDataOnly=0
						else
							if (toLower(trim(second(1)))=="yes") then
								this%PhaseTheDataOnly=1
							else
								write(error_unit,*) "Stop - Phasing only option incorrectly specified"
								stop
							endif
						endif

					case("userdefinedalphaphaseanimalsfile")
						this%UserDefinedHD=0
						if (tolower(trim(second(1)))/="none") then
							this%UserDefinedHD=1
							this%animalPhaseFile = second(1)
						endif

					case("prephasedfile")
						if (toLower(trim(second(1)))=="none") then
							this%PrePhased=0
						else
							this%PrePhased=1
							this%InbredAnimalsFile = second(1)
						endif
					case("bypassgeneprob")
						write(error_unit,*) "The Geneprob has been moved to legacy and is no longer in use"
					case("restartoption")
						read(second(1),*) this%restartOption

						if (this%restartoption >2) then
							write(error_unit,*) "Error - RestartOption can only be the following:"
							write(error_unit,*) "0 - run all"
							write(error_unit,*) "1 - make dirs run phasing"
							write(error_unit,*) "2 - run imputation/hmm"
							stop
						endif




					case("cluster")
						if (toLower(trim(second(1)))=="no") then
							this%cluster=.false.
						else
							if (toLower(trim(second(1)))=="yes") then
								this%cluster=.true.
							else
								write(error_unit,*) "Error: Cluster incorrectly specified, please use yes or no"
							endif
						endif
					case("parallelprocessors")
						read(second(1), *,iostat=stat) this%useProcs
						if (stat /=0) then
							write(error_unit,*) "ERROR: Parallel Processors set Incorrectly"
						endif
						if (this%useProcs > OMP_get_num_procs()) then
							write(error_unit,*) "WARNING - more processors than are available are specified under parallelprocessors"
							write(error_unit,*) this%useProcs, " set vs ",OMP_get_num_procs()," available"
						endif
					case("resultfolderpath")
						read(second(1), *) this%resultFolderPath
						case default
						write(*,"(A,A)") trim(tag), " is not valid for the AlphaImpute Spec File."
						cycle

					case("useferdosi")
						if(ToLower(trim(second(1))) == "no") then
							this%useFerdosi = .false.
						else if (ToLower(trim(second(1))) == "yes") then
							this%useFerdosi = .true.
						endif
					case("modelrecomb")
						if(ToLower(trim(second(1))) == "no") then
							this%modelrecomb = .false.
						else if (ToLower(trim(second(1))) == "yes") then
							this%modelrecomb = .true.
						endif

					case("stopafterplink")
						if(ToLower(trim(second(1))) == "no") then
							this%stopAfterPlink = .false.
						else if (ToLower(trim(second(1))) == "yes") then
							this%stopAfterPlink = .true.
						endif
					case("usechroms")
						allocate(this%usechroms(size(second)))
						do i= 1, size(second)
							read(second(i), *) this%useChroms(i)

						enddo

					end select
				end if
			end do READFILE
			deallocate(tag)
			deallocate(tmptag)
		endif

		! Set parameters for parallelisation
		if (this%nPhaseInternal==2) then
			this%nAgreeInternalHapLibElim=1
		endif
		if (this%nPhaseInternal==4) then
			this%nAgreeInternalHapLibElim=2
		endif
		if (this%nPhaseInternal==6) then
			this%nAgreeInternalHapLibElim=3
		endif
		if (this%nPhaseInternal>6) then
			this%nAgreeInternalHapLibElim=4
		endif



		this%PercGenoForHD=this%PercGenoForHD/100
		this%PercSnpMiss=this%PercSnpMiss/100
		this%SecondPercGenoForHD=this%SecondPercGenoForHD/100




		if (.not. allocated(this%CoreAndTailLengths) .or. .not. allocated(this%CoreLengths) .and. (this%ManagePhaseOn1Off0 /= 0)) then
			write(error_unit,*) "warning - CoreLengths or CoreAndTailLengths have not been specified, will be calculated instead"
		endif

		CALL OMP_SET_NUM_THREADS(this%useProcs)
	end subroutine ReadInParameterFile



	function validateAlphaImputeSpec(params) result(res)
		class(AlphaImputeInput)  :: params
		LOGICAL :: res
		integer :: i

		if (size(params%coreandtaillengths) /= size(params%coreLengths)) then

			write(error_unit,*) "ERROR: core and coreandtaillengths do not match"
			res = .false.
			return
		endif
		do i=1, size(params%CoreAndTailLengths)

			if (params%CoreAndTailLengths(i) > params%nsnp) then
				write(error_unit,*) "ERROR: core and tail lengths are more than NSNPS"
				res = .false.
				return
			endif

			if (params%coreLengths(i) > params%nsnp) then
				write(error_unit,*) "ERROR: core and tail lengths are more than NSNPS"
				res = .false.
				return
			endif
		enddo


		res = .true.
	end function validateAlphaImputeSpec

	!---------------------------------------------------------------------------
	!< @brief Generates core and tail length.
	!< @author  David Wilson david.wilson@roslin.ed.ac.uk
	!< @date    November 12 2017
	!---------------------------------------------------------------------------
	subroutine calculateCoresAndTails(nsnps, cores, tails, size)

		use omp_lib
		use PARAMETERS

		integer, intent(in) :: nsnps
		integer(kind=int32), allocatable, intent(out) :: cores(:),tails(:) !returned cores and tails
		integer, intent(inout) :: size !< if 0, will use number of processors
		integer :: tailSize,i

		if (size == 0)then
			size = OMP_get_num_procs()
		endif

		allocate(cores(size))
		allocate(tails(size))

		! set tails to be the following
		tailSize = (nsnps / size) / 4


		do i=1, size

			cores(i) = (nsnps / size)*i - tailSize
			tails(i) = (nsnps / size)*i
		enddo

		write(error_unit,*) "core lengths used:", cores

		write(error_unit,*) "tails lengths used:", tails
	end subroutine calculateCoresAndTails



	!---------------------------------------------------------------------------
	!< @brief Writes out the spec file options used to file AlphaImputeSpecFileUsed.txt
	!< @author  David Wilson david.wilson@roslin.ed.ac.uk
	!< @date    October 26, 2017
	!---------------------------------------------------------------------------
	subroutine writeOutSpecOptions(params, path)

		use Global

		class(AlphaImputeInput), intent(in) :: params
		integer :: unit,i
		character(len=:), allocatable :: tmpString
		character(len=*), optional, intent(in):: path


		if (present(path)) then
			open(newunit=unit, file=path, status='unknown')
		else
			open(newunit=unit, file="AlphaImputeSpecFileUsed.txt", status='unknown')
		endif

		write(unit, *) "PedigreeFile,",trim(params%PedigreeFile)
		write(unit, *) "GenotypeFile,",trim(params%GenotypeFile)
		write(unit, *) "TrueGenotypeFile,",trim(params%TrueGenotypeFile)

		if (params%plinkOutput) then
			write(unit, *) "plinkOutput,yes"
		else
			write(unit, *) "plinkOutput,no"
		endif

		if (params%plinkBinary) then
			write(unit, *) "plinkinputfile,binary,",trim(params%plinkinputfile)
		else
			write(unit, *) "plinkinputfile,text,",trim(params%plinkinputfile)
		endif
		if (params%sexOpt == 1) then
			write(unit, '(a,a,a)', advance="no") "sexchrom,","yes,",trim(params%genderFile)

			if (params%HetGameticStatus == 1) then
				write(unit, *) ",male"
			else if (params%HetGameticStatus == 2) then
				write(unit, *) ",female"
			endif
		else
			write(unit, *) "sexchrom,","no"
		endif
		write(unit, *) "nsnp,",params%nsnp
		write(unit, *) "multiplehdpanels,",params%MultiHd
		! TODO numbsnpxchip
		write(unit, *) "hdanimalsthreshold,",params%PercGenoForHD

		if (params%inteditstat == 1) then
			write(unit, *) "internaledit,","yes"
		else
			write(unit, *) "internaledit,","no"
		endif

		if (params%outopt == 1) then
			tmpString = "allsnpout"
		else
			tmpString = "editedsnpout"
		endif
		write(unit, *) "editingparameters,",params%PercGenoForHD*100,",",params%PercSnpMiss*100,",", params%SecondPercGenoForHD*100,",",tmpString
		deallocate(tmpString)

		if (params%outputonlygenotypedanimals) then
			write(unit, *) "outputonlygenotypedanimals,","yes"
		else
			write(unit, *) "outputonlygenotypedanimals,","no"
		endif

		if (params%noPhasing == 0) then
			write(unit, *) "numberphasingruns,","nophase"
		else if ( params%managephaseon1off0 == 0) then
			write(unit, *) "numberphasingruns,","phasedone",trim(params%phasePath),",", params%nPhaseInternal
		else
			write(unit, *) "numberphasingruns,", params%nPhaseExternal
		endif

		write(unit, "(a,*(a,I5,:))") "coreandtaillengths",(',',params%coreandtaillengths(i),i=1, size(params%coreandtaillengths))

		write(unit, "(a,*(a,I5,:))") "corelengths", (',',params%corelengths(i),i=1, size(params%coreLengths))

		if (params%PedFreePhasing == 1) then
			write(unit, *) "pedigreefreephasing,","yes"
		else
			write(unit, *) "pedigreefreephasing,","no"
		endif

		write(unit, *) "genotypeerror,", params%GenotypeErrorPhase

		if (params%largedatasets) then
			write(unit, *) "largedatasets,", "yes",params%PhaseSubsetSize,",",params%PhaseNIterations,",",trim(params%iterateMethod)
		else
			write(unit, *) "largedatasets,", "no"
		endif

		write(unit, *) "minoverlaphaplotype,", params%minoverlaphaplotype

		if (params%alphaphaseoutput== 0) then
			write(unit,*) "alphaphaseoutput,","no"
		elseif (params%alphaphaseoutput== 1) then
			write(unit,*) "alphaphaseoutput,","yes"
		elseif (params%alphaphaseoutput == 2) then
			write(unit,*) "alphaphaseoutput,","binary"
		elseif (params%alphaphaseoutput == 3) then
			write(unit,*) "alphaphaseoutput,","verbose"
		endif

		write(unit,*) "internaliterations,", params%internaliterations


		if (params%ConservativeHapLibImputation == 1) then
			write(unit, *) "conservativehaplotypelibraryuse,","yes"
		else
			write(unit, *) "conservativehaplotypelibraryuse,","no"
		endif

		write(unit,*) "wellphasedthreshold,", params%WellPhasedThresh


		if (params%hmmoption== RUN_HMM_NO) then
			write(unit,*) "hmmoption,","no"
		elseif (params%hmmoption== RUN_HMM_YES) then
			write(unit,*) "hmmoption,","yes"
		elseif (params%hmmoption == RUN_HMM_ONLY) then
			write(unit,*) "hmmoption,","only"
		elseif (params%hmmoption == RUN_HMM_PREPHASE) then
			write(unit,*) "hmmoption,","prephase"
		elseif (params%hmmoption == RUN_HMM_NGS) then
			write(unit,*) "hmmoption,","ngs"
		endif

		! write(unit,*) "hmmparameters,", params%nHapInSubH,params%HmmBurnInRound, params%nRoundsHmm, params%nRoundsHmm, params%useProcs, params%idum

		write(unit,*) "templatehaplotypes,",params%nHapInSubH

		write(unit,*) "burninrounds,",params%HmmBurnInRound
		write(unit,*) "rounds,",params%nRoundsHMM

		write(unit,*) "seed,",params%idum
		write(unit,*) "phasedanimalsthreshold,",params%phasedThreshold
		write(unit,*) "thresholdimputed,",params%imputedThreshold
		write(unit,*) "haplotypeslist,",trim(params%HapListFile)


		if (params%PreProcess) then
			write(unit, *) "preprocessdataonly,","yes"
		else
			write(unit, *) "preprocessdataonly,","no"
		endif

		if (params%PhaseTheDataOnly == 1) then
			write(unit, *) "phasingonly,","yes"
		else
			write(unit, *) "phasingonly,","no"
		endif

		write(unit, *) "userdefinedalphaphaseanimalsfile,",trim(params%animalPhaseFile)

		write(unit, *) "prephasedfile,",trim(params%InbredAnimalsFile)

		if (params%cluster) then
			write(unit, *) "cluster,","yes"
		else
			write(unit, *) "cluster,","no"
		endif

		write(unit,*) "restartoption,",params%restartoption

		if (params%cluster) then
			write(unit, *) "cluster,","yes"
		else
			write(unit, *) "cluster,","no"
		endif

		write(unit,*) "parallelprocessors,",params%useProcs
		write(unit,*) "resultfolderpath,",trim(params%resultfolderpath)

		if (params%useferdosi) then
			write(unit, *) "useferdosi,","yes"
		else
			write(unit, *) "useferdosi,","no"
		endif

		if (params%modelrecomb) then
			write(unit, *) "modelrecomb,","yes"
		else
			write(unit, *) "modelrecomb,","no"
		endif

		if (params%stopAfterPlink) then
			write(unit, *) "stopafterplink,","yes"
		else
			write(unit, *) "stopafterplink,","no"
		endif

		if (allocated(params%usechroms)) then
			write(unit, "(a,*(a,I5,:))") "usechroms",(',',params%usechroms(i),i=1, size(params%usechroms))
		else
			write(unit, *) "usechroms,0"
		endif

		! do i=1,size(params%usechroms)
		! 	write(unit, '(I3,a)', advance=no) params%usechroms(i),","
		! enddo




		close(unit)
	end subroutine writeOutSpecOptions



end module AlphaImputeSpecFileModule




