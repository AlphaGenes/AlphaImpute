!-----------------------------------------------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-----------------------------------------------------------------------------------------------------------------------
!
! MODULE: AlphaImputeInMod
!
!> @file        AlphaImputeInMod.f90
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




module AlphaImputeInMod
    use iso_fortran_env

    type AlphaImputeInput
        ! box 1
        character(len=300):: PedigreeFile = "Pedigree.txt",GenotypeFile="Genotypes.txt",TrueGenotypeFile="TrueGenotypes.txt",GenderFile="None",InbredAnimalsFile="None", HapListFile="None"
        integer(kind=1) :: TrueGenos1None0
        logical :: PlinkFormat, VCFFormat

        ! box 2
        integer(kind=1) :: SexOpt,HetGameticStatus,HomGameticStatus

        ! box 3
        integer(kind=int32) :: nSnp,MultiHD
        integer(kind=int32), allocatable :: nSnpByChip(:)
        real(real32) :: PercGenoForHD !TODO don't see why this is a real

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
        integer(kind=int32) :: nProcessors,nProcessGeneProb,nProcessAlphaPhase
        character(len=10) :: iterateMethod
        integer :: minoverlaphaplotype
        ! box 6
        integer(kind=int32) :: InternalIterations
        integer(kind=1) :: ConservativeHapLibImputation
        real(kind=real32) :: WellPhasedThresh

        ! box 7
        ! idum is seed
        integer(kind=int32) :: idum
        integer(kind=int32) :: nHapInSubH,useProcs,nRoundsHmm,HmmBurnInRound
        real(kind=real32) :: phasedThreshold,imputedThreshold
        logical :: HapList=.FALSE.
        integer(kind=1) :: HMMOption

        ! box 8

        logical :: PreProcess
        integer(kind=1) :: PhaseTheDataOnly

        integer(kind=1) :: UserDefinedHD,PrePhased,BypassGeneProb,RestartOption

        integer :: AnimalFileUnit, prePhasedFileUnit, pedigreeFileUnit,genotypeFileUnit,GenderFileUnit,HapListUnit

        ! other
        integer(kind=int32) :: nSnpRaw,nAgreeImputeHDLib,nAgreeParentPhaseElim,nAgreeGrandParentPhaseElim,nAgreePhaseElim,nAgreeInternalHapLibElim
    contains
        procedure :: ReadInParameterFile
    end type AlphaImputeInput

    type(AlphaImputeInput),target, allocatable :: defaultInput

contains

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

        class(AlphaImputeInput), intent(inout),target :: this
        integer :: unit,IOStatus,MultipleHDpanels,i
        character(len=*), intent(in) :: SpecFile

        character(len=300) :: first, line
        character(len=:), allocatable::tag
        character(len=300),dimension(:),allocatable :: second

        this%MultiHD = 0
        this%minoverlaphaplotype = 0

        open(newunit=unit, file=SpecFile, action="read", status="old")
        IOStatus = 0
        
        READFILE: do while (IOStatus==0)
            read(unit,"(A)", IOStat=IOStatus)  line
            if (len_trim(line)==0) then
                CYCLE
            end if

            call splitLineIntoTwoParts(trim(line), first, second)
            tag = parseToFirstWhitespace(first)
            if (first(1:1)=="=" .or. len(trim(line))==0) then
                cycle
            else
                select case(trim(tag))

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
                    this%PlinkFormat = .FALSE.
                    this%VCFFormat = .FALSE.
                    if (size(second) >1) then
                        select case(trim(toLower(second(2))))
                        case('plink')
                            this%PlinkFormat = .TRUE.
                        case('vcf')
                            this%VCFFormat = .TRUE.
                        end select
                    endif
                endif
            case("truegenotypefile")
                if (.not. allocated(second)) then
                    write(*, "(A,A)") "No true genotype file specified. Using default filename: ", this%TrueGenotypeFile
                else
                    write(this%TrueGenotypeFile, "(A)") second(1)
                    if (this%TrueGenotypeFile=="None") then
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
                if (trim(second(1))=="Yes") then
                    this%genderFile = second(2)
                    this%HetGameticStatus=9
                    if (trim(second(3))=="Male") then        ! Species  with heterogametic males
                        this%HetGameticStatus=1                              ! My father is heterogametic
                        this%HomGameticStatus=2                              ! My mother is homogametic
                    endif
                    if (trim(second(3))=="Female") then      ! Species with heterogametic females
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
                if (trim(second(1))=="No") then
                    this%SexOpt=0
                endif
                if (this%SexOpt==9) then
                    print*, "Warning - Sex chromosome status is misspecified"
                    stop
                endif

                ! box 3 inputs

            case("numbersnp")
                read(second(1),*) this%nsnp
                if (this%nsnp>240000) then
                    print*, "Contact John Hickey if you want to do more than 240,000 SNP"
                    stop 3001
                endif

            case("multiplehdpanels")
                ! Get the information of Multiple HD chips
                ! MultipleHDpanels
                if (trim(toLower(second(1))) /= "no") then
                    read(second(1),*) MultipleHDpanels
                    if (MultipleHDpanels/=0) this%MultiHD=MultipleHDpanels
                    ! if ((trim(MultipleHDpanels)/='Yes').and.(trim(MultipleHDpanels)/='No')) then
                    !     write (*,*) "Please, provide a valid option,"
                    !     write (*,*) "MultipleHDpanels only acepts 'No' or 'Yes'"
                    !     stop
                    ! endif
                    ! Snps of the multiple HD panels
                    if (allocated(this%nSnpByChip)) then
                        deallocate(this%nSnpByChip)
                    endif
                    allocate(this%nSnpByChip(this%MultiHD))
                endif
            case("numbersnpxchip")
                do i=1,this%MultiHD
                    read(second(i),*) this%nSnpByChip(i)
                enddo
            case("hdanimalsthreshold")
                read(second(1),*) this%PercGenoForHD

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
                    if (this%MultiHD==0) this%PercGenoForHD=90.0
                    this%outopt=1
                endif
                this%PercGenoForHD=this%PercGenoForHD/100
                this%PercSnpMiss=this%PercSnpMiss/100
                this%SecondPercGenoForHD=this%SecondPercGenoForHD/100
                cycle
                4000 print*, "Output options incorrectly specified"
                print*, "Beware!!!!! AlphaImpute is case sensitive"
                stop


                ! box 5
            case("numberphasingruns")
                this%noPhasing = 1

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
                else if(ToLower(trim(second(1))) == "nophase") then
                    this%noPhasing = 1
                    this%managephaseon1off0 = 0
                else
                    this%managephaseon1off0 = 1
                    read(second(1),*) this%nPhaseExternal

                    if (this%nPhaseExternal > 40) then
                        write(error_unit,*) "Error: Too many phasing runs required. The most that this program supports is 40."
                        stop 40511
                    else if(this%nPhaseExternal < 2) then
                        write(error_unit,*) "Error: Too few phasing runs requested. The minimum this program supports is 2, 10 is reccomended."
                        stop 40512
                    endif

                    this%nPhaseInternal = 2*this%nPhaseExternal
                    if (allocated(this%CoreAndTailLengths)) then
                        deallocate(this%CoreAndTailLengths)
                    endif
                    if (allocated(this%CoreLengths)) then
                        deallocate(this%CoreLengths)
                    endif
                    allocate(this%CoreAndTailLengths(this%nPhaseExternal))
                    allocate(this%CoreLengths(this%nPhaseExternal))

                endif
                cycle
                4051 write(error_unit,*) "NumberPhasingRuns has been set incorrectly"
                stop 4051

            case("coreandtaillengths")
                if (.not. allocated(this%CoreAndTailLengths)) then
                    write(error_unit,*) "Error: numberofphasingruns is not defined. Please define this before CoreAndTailLengths and CoreLengths"
                    stop 40521
                endif
                if (size(second) /= this%nPhaseExternal) then
                    write(error_unit,*) "Error: numberofphasingruns is set to a different number of parameters than what is specified here. \n Please set this to the same number of parameters that are given for CoreAndTailLengths and CoreLengths"
                    stop 40522
                endif
                do i=1,size(second)
                    read(second(i), *) this%CoreAndTailLengths(i)
                enddo

            case("corelengths")
                if (.not. allocated(this%corelengths)) then
                    write(error_unit,*) "Error: numberofphasingruns is not defined. Please define this before CoreAndTailLengths and CoreLengths"
                    stop 40531
                endif
                if (size(second) /= size(this%corelengths)) then
                    write(error_unit,*) "Error: numberofphasingruns is set to a different number of parameters than what is specified here. \n Please set this to the same number of parameters that are given for CoreAndTailLengths and CoreLengths"
                    stop 40532
                endif
                do i=1,size(second)
                    read(second(i),*) this%corelengths(i)
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
                    read(second(1), *) this%GenotypeErrorPhase
                endif

            case("numberofprocessorsavailable")
                read(second(1),*) this%nProcessors

            case("largedatasets")
                if (ToLower(trim(second(1)))== "yes") then
                    this%largedatasets=.true.
                    read(second(2),*) this%PhaseSubsetSize
                    read(second(3),*) this%PhaseNIterations
                else
                    this%largedatasets=.false.

                endif

            case("iteratemethod")
                 if (ToLower(trim(second(1)))== "off") then
                    this%iterateMethod  = "Off"
                 else if (ToLower(trim(second(1)))== "randomorder") then
                    this%iterateMethod  = "RandomOrder"
                 else if (ToLower(trim(second(1)))== "inputorder") then
                    this%iterateMethod  = "InputOrder"
                else
                    this%iterateMethod= "Off"

                endif

            
            case("minoverlaphaplotype")
                read(second(1),*) this%minoverlaphaplotype
                if (this%minoverlaphaplotype < 0) then
                    write(error_unit,*) "ERROR: Min minoverlap haplotype size is set incorrectly!"
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
                if (toLower(trim(second(1)))=='Prephase') this%hmmoption=RUN_HMM_PREPHASE
                if (toLower(trim(second(1)))=="ngs") this%hmmoption=RUN_HMM_NGS
                if (this%hmmoption==RUN_HMM_NULL) then
                    write(error_unit,*), "this%hmmoption not correctly specified"
                    stop
                endif

            case("templatehaplotypes")
                if (.not. allocated(second)) then
                    write(error_unit,*) "templatehaplotypes not set correctly"
                    stop
                endif
                read(second(1), *) this%nHapInSubH

            case("burninrounds")
                read(second(1), *) this%HmmBurnInRound
            case("rounds")
                read(second(1), *)this%nRoundsHMM
            case("parallelprocessors")
                read(second(1), *) this%useProcs
            case("seed")
                read(second(1), *)this%idum

            case("thresholdformissingalleles")
                read(second(1), *) this%phasedThreshold

            case("phasedanimalsthreshold")
                read(second(1), *) this%phasedThreshold
                
            case("thresholdimputed")
                read(second(1), *) this%imputedThreshold
                  
            case("wellimputedthreshold")
                read(second(1), *) this%imputedThreshold
            case("haplotypeslist")
                if (.not. allocated(second)) then
                    write(*, "(A,A)") "No list of haploytpes specified"
                else
                    if (trim(second(1)) /= "None") then
                        this%HapList = .TRUE.
                        this%HapListFile = trim(second(1))
                        open (newunit=this%HapListUnit, file=this%HapListFile, status='old')
                    endif
                endif

            !  box 8
            case("preprocessdataonly")
                if (second(1)=="No") then
                    this%PreProcess=.FALSE.
                else
                    if (second(1)=="Yes") then
                        this%PreProcess=.TRUE.
                    else
                        write(error_unit,*) "Stop - Preprocess of data option incorrectly specified"
                        stop
                    endif
                endif

            case("phasingonly")
                if (second(1)=="No") then
                    this%PhaseTheDataOnly=0
                else
                    if (second(1)=="Yes") then
                        this%PhaseTheDataOnly=1
                    else
                        write(error_unit,*) "Stop - Phasing only option incorrectly specified"
                        stop
                    endif
                endif

            case("userdefinedalphaphaseanimalsfile")
                this%UserDefinedHD=0
                if (second(1)/="None") then
                    this%UserDefinedHD=1
                    open (newunit=this%AnimalFileUnit,file=trim(second(1)),status="old")
                endif

            case("prephasedfile")
                if (second(1)=="None") then
                    this%PrePhased=0
                else
                    this%PrePhased=1
                    this%InbredAnimalsFile = second(1)
                    open (newunit=this%prePhasedFileUnit,file=trim(second(1)),status="old")
                endif
            case("bypassgeneprob")
                if (trim(second(1))=="No") then
                    this%BypassGeneProb=0
                else if (trim(second(1))=="Yes") then
                    this%BypassGeneProb=1
                else if (trim(second(1))=="Probabilities") then
                    this%bypassgeneprob=2
                else
                    write(error_unit,*) "bypassgeneprob not correctly specified"
                    stop
                endif
            case("restartoption")
                read(second(1),*) this%restartOption

            case default
                write(*,"(A,A)") trim(tag), " is not valid for the AlphaImpute Spec File."
                cycle
            end select
        end if
    end do READFILE
    deallocate(tag)
    open (newUnit=this%pedigreeFileUnit,file=trim(this%PedigreeFile),status="old")
    open (newUnit=this%genotypeFileUnit,file=trim(this%GenotypeFile),status="old")
    if (this%SexOpt==1) open (newUnit=this%genderFileUnit,file=trim(this%GenderFile),status="old")

    this%nProcessAlphaPhase=this%nProcessors-this%nProcessGeneProb ! Never used!

    ! Set parameters for parallelisation
    if (this%nPhaseInternal==2) then
        this%nAgreeImputeHDLib=1
        this%nAgreeParentPhaseElim=1
        this%nAgreePhaseElim=1
        this%nAgreeInternalHapLibElim=1
    endif
    if (this%nPhaseInternal==4) then
        this%nAgreeImputeHDLib=2
        this%nAgreeParentPhaseElim=2
        this%nAgreePhaseElim=2
        this%nAgreeInternalHapLibElim=2
    endif
    if (this%nPhaseInternal==6) then
        this%nAgreeImputeHDLib=3
        this%nAgreeParentPhaseElim=3
        this%nAgreePhaseElim=3
        this%nAgreeInternalHapLibElim=3
    endif
    if (this%nPhaseInternal>6) then
        this%nAgreeImputeHDLib=4
        this%nAgreeParentPhaseElim=4
        this%nAgreeGrandParentPhaseElim=4
        this%nAgreePhaseElim=4
        this%nAgreeInternalHapLibElim=4
    endif

    this%nSnpRaw = this%nsnp
    
    !$  CALL OMP_SET_NUM_THREADS(this%nProcessors)
end subroutine ReadInParameterFile

end module AlphaImputeInMod
