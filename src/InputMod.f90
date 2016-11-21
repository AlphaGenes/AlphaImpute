!-----------------------------------------------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-----------------------------------------------------------------------------------------------------------------------
!
! MODULE: InputMod
!
!> @file        InputMod.f90
!
! DESCRIPTION:
!> @brief       Module holding reading subroutines
!>
!> @details     This MODULE contains a class which contains all subroutines to read in a file.
!
!> @author      Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
!
!> @date        Nov 10, 2016
!
!> @version     0.0.1 (alpha)
!
! REVISION HISTORY:
! 2016.11.10  Rantolin - Initial Version
!
!-----------------------------------------------------------------------------------------------------------------------

module InputMod
  use iso_fortran_env

  private
  public:: CountInData, ReadInData, ReadSeq, ReadGenos


  contains

!---------------------------------------------------------------------------
! DESCRIPTION:
!> @brief      Count the number of individuals
!
!> @details    This subroutine encapsulates the two different the subroutines
!>             that count the number of animals and markers
!
!> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
!
!> @date       Nov 16, 2016
!
! PARAMETERS:
!> @param[in]
!---------------------------------------------------------------------------
  subroutine CountInData
  use Global
  use AlphaImputeInMod
  implicit none

  type(AlphaImputeInput), pointer :: inputParams
  inputParams => defaultInput

  if (inputParams%hmmoption == RUN_HMM_NGS) then
    call CountInSequenceData
  else
    call CountInGenotypeData
  end if
  end subroutine CountInData

!---------------------------------------------------------------------------
! DESCRIPTION:
!> @brief      Count the number of individuals
!
!> @details    This subroutine counts the number of individuals genotyped and
!>             the number of individuals in the pedigree.
!>             If no pedigree information is available, then these to counts
!>             are equal to the number of individuals genotyped.
!
!> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
!
!> @date       Nov 10, 2016
!
! PARAMETERS:
!> @param[in]
!---------------------------------------------------------------------------
  subroutine CountInGenotypeData

  use Global
  use GlobalVariablesHmmMaCH
  use AlphaImputeInMod
  implicit none

  integer :: k
  character (len=300) :: dumC
  type(AlphaImputeInput), pointer :: inputParams

  inputParams => defaultInput
  do
    read (inputParams%pedigreeFileUnit,*,iostat=k) dumC
    nAnisRawPedigree=nAnisRawPedigree+1
    if (k/=0) then
      nAnisRawPedigree=nAnisRawPedigree-1
      exit            ! This forces to exit if an error is found
    endif
  enddo
  rewind(inputParams%genotypeFileUnit)

  print*, " ",nAnisRawPedigree," individuals in the pedigree file"
  nObsDataRaw=nAnisRawPedigree

  do
    read (inputParams%genotypeFileUnit,*,iostat=k) dumC
    nAnisG=nAnisG+1
    if (k/=0) then
      nAnisG=nAnisG-1
      exit
    endif
  enddo

  rewind(inputParams%genotypeFileUnit)

  if (inputParams%hmmoption == RUN_HMM_NGS) then
    if(mod(nAnisG,2)==0) then
      nAnisG=nAnisG/2
    else
      write(0,*) "Error: The number of lines in the file of reads is not even. Is the file corrupt?"
      write(0,*) "The program will now stop"
      stop
    endif
  endif

  print*, " ",nAnisG," individuals in the genotype file"

  ! This is incoherent with functions ReadInData and ReadInParameterFile
  if (trim(inputParams%PedigreeFile)=="NoPedigree") nAnisRawPedigree=nAnisG

  end subroutine CountInGenotypeData

!---------------------------------------------------------------------------
! DESCRIPTION:
!> @brief      Count the number of individuals
!
!> @details    This subroutine counts the number of individuals sequenced and
!>             the number of individuals in the pedigree.
!>             If no pedigree information is available, then these to counts
!>             are equal to the number of individuals genotyped.
!
!> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
!
!> @date       Nov 10, 2016
!
! PARAMETERS:
!> @param[in]
!---------------------------------------------------------------------------
  subroutine CountInSequenceData
  use Global
  use GlobalVariablesHmmMaCH
  use AlphaImputeInMod
  implicit none

  integer :: k, markers
  character (len=1024*1024) :: dumC
  type(AlphaImputeInput), pointer :: inputParams

  inputParams => defaultInput
  do
    read (inputParams%pedigreeFileUnit,*,iostat=k) dumC
    nAnisRawPedigree=nAnisRawPedigree+1
    if (k/=0) then
      nAnisRawPedigree=nAnisRawPedigree-1
      exit            ! This forces to exit if an error is found
    endif
  enddo
  rewind(inputParams%genotypeFileUnit)

  print*, " ",nAnisRawPedigree," individuals in the pedigree file"
  nObsDataRaw=nAnisRawPedigree


  if (inputParams%VCFFormat) then
    read(inputParams%genotypeFileUnit,'(A)',iostat=k) dumC
    nAnisG = nWords(trim(dumC)) - 5
    print *, nAnisG

    markers = 0
    do
      read (inputParams%genotypeFileUnit,*,iostat=k) dumC
      markers = markers + 1
      if (k/=0) then
        markers = markers - 1
        exit
      endif
    enddo

    if (markers /= inputParams%nsnp) then
      write(0,*) "Error: The number of markers in the file of sequence data does not match the"
      write(0,*) "       number of markers provider by the user is not even. Is the file corrupt?"
      stop
    end if

  else

    do
      read (inputParams%genotypeFileUnit,*,iostat=k) dumC
      nAnisG = nAnisG + 1
      if (k/=0) then
        nAnisG = nAnisG - 1
        exit
      endif
    enddo

    rewind(inputParams%genotypeFileUnit)

    if (inputParams%hmmoption == RUN_HMM_NGS) then
      if(mod(nAnisG,2)==0) then
        nAnisG=nAnisG/2
      else
        write(0,*) "Error: The number of lines in the file of reads is not even. Is the file corrupt?"
        write(0,*) "The program will now stop"
        stop
      endif
    endif
  end if

  print*, " ",nAnisG," individuals in the genotype file"

  ! This is incoherent with functions ReadInData and ReadInParameterFile
  if (trim(inputParams%PedigreeFile)=="NoPedigree") nAnisRawPedigree=nAnisG

  end subroutine CountInSequenceData

!---------------------------------------------------------------------------
! DESCRIPTION:
!> @brief      Read files
!
!> @details    This subroutine reads the pedigree data and the gender file
!>             in case of Sex Chromosome
!
!> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
!
!> @date       Nov 10, 2016
!
! PARAMETERS:
!> @param[in]
!---------------------------------------------------------------------------
  subroutine ReadInData
  use GlobalPedigree
  use Global
  use AlphaImputeInMod
  implicit none

  integer :: i,j,k,CountLinesGender,GenCode,AnimalPresent
  character(len=300) :: dumC
  type(AlphaImputeInput), pointer :: inputParams
  integer, allocatable, dimension(:) :: Temp
  inputParams=> defaultInput

  allocate(Temp(inputParams%nsnp))
  allocate(GenotypeId(nAnisG))
  allocate(Genos(0:nAnisG,inputParams%nsnp))
  allocate(GenderId(nAnisRawPedigree))
  allocate(GenderRaw(nAnisRawPedigree))

  Genos(0,:)=9

  ! Read the pedigree information
  
  ped = initPedigree(inputParams%pedigreefile, nAnisRawPedigree)
  ! Read the genotype file
  if (inputParams%hmmoption /= RUN_HMM_NGS) then
      rewind(inputParams%genotypeFileUnit)
      if (.not. inputParams%PlinkFormat) then
        call ReadGenos(inputParams%genotypeFileUnit)
      else
        call ReadPlink(inputParams%genotypeFileUnit)
      end if
  endif
  close(inputParams%pedigreeFileUnit)
  close(inputParams%genotypeFileUnit)

  ! Read the gender file if imputing the sex chromosome
  if (inputParams%SexOpt==1) then
    ped = initPedigree(inputParams%pedigreefile, nAnisRawPedigree,inputParams%genderFile)
  else 
    ped = initPedigree(inputParams%pedigreefile, nAnisRawPedigree)
  endif
  deallocate(temp)
  end subroutine ReadInData

  !#############################################################################################################################################################################################################################
  subroutine ReadSeq(ReadsFileUnit)
  use GlobalPedigree
  use Global
  use alphaimputeinmod
  implicit none

  integer, intent(inout) :: ReadsFileUnit
  integer :: i,j

  type(AlphaImputeInput), pointer :: inputParams
  logical :: opened, named
  character(len=300) :: ReadsFile
  integer, allocatable,dimension (:) :: ReferAlleleLine, AlterAlleleLine

  inputParams => defaultInput
  allocate(ReferAllele(0:nAnisG,inputParams%nsnp))
  allocate(AlterAllele(0:nAnisG,inputParams%nsnp))
  allocate(ReferAlleleLine(inputParams%nsnp))
  allocate(AlterAlleleLine(inputParams%nsnp))

  Reads=0

#ifdef DEBUG
  write(0,*) "DEBUG: [ReadSeq] Reads size=", size(Reads,1)
#endif

  ! open (unit=3,file=trim(readsFileName),status="old")
  inquire(unit=ReadsFileUnit, opened=opened, named=named, name=ReadsFile)
  if (.NOT. opened .and. named) then
    open(unit=ReadsFileUnit, file=ReadsFile, status='unknown')
  else if (.NOT. named) then
    ! write(0, *) "ERROR - Something went wrong when trying to read the file of pre-phased data"
    open(newunit=ReadsFileUnit, file=inputParams%GenotypeFile)
  end if

#ifdef DEBUG
  write(0,*) "DEBUG: [ReadSeq] Reading sequence data..."
#endif

  if (inputParams%VCFFormat) then
    call readVCF(ReadsFileUnit, GenotypeId, ReferAllele, AlterAllele, inputParams%nsnp, inputParams%nsnp, 1, inputParams%nsnp, nAnisG)
  else
    do i=1,nAnisG
      read (ReadsFileUnit,*) GenotypeId(i), ReferAllele(i,:)
      read (ReadsFileUnit,*) GenotypeId(i), AlterAllele(i,:)
    end do
  end if

  do i=1,nAnisG
    do j=1,inputParams%nsnp
      if (ReferAllele(i,j)>=MAX_READS_COUNT) ReferAllele(i,j)=MAX_READS_COUNT-1
      if (AlterAllele(i,j)>=MAX_READS_COUNT) AlterAllele(i,j)=MAX_READS_COUNT-1
      Reads(i,j)=AlterAllele(i,j)+ReferAllele(i,j)
    enddo
  enddo

  open(unit=999, file='borrar.txt',status='unknown')
  do i =1, nAnisG
    write(999,'(a20,100i2)') ReferAllele(i,:)
    write(999,'(a20,100i2)') AlterAllele(i,:)
  end do
  close(999)

#ifdef DEBUG
  write(0,*) "DEBUG: [ReadSeq] Sequence data read"
#endif

  close(ReadsFileUnit)

  end subroutine ReadSeq

!#############################################################################################################################################################################################################################
  subroutine readVCF(ReadsFileUnit, Ids, RefAll, AltAll, nSnpIn, SnpUsed, StartSnp, EndSnp, nIndivIn)
    ! subroutine readRogerData(filename, Ids, position, quality, SequenceData, nSnpIn, SnpUsed, StartSnp, EndSnp, nIndivIn)
    use GlobalPedigree
    use Global
    use AlphaImputeInMod
    use omp_lib
    implicit none
    !filename has the name of the file to be read
    integer, intent(inout) :: ReadsFileUnit
    character(len=100), allocatable, dimension(:), intent(out):: Ids
    !info holds the quality and the poisition (

    character(len=100), dimension(:), allocatable::dumE, dumC
    ! real(real64), allocatable, dimension(:), intent(out):: quality
    ! integer(int32), dimension(:), allocatable, intent(out):: position
    integer(int32), dimension(:,:), allocatable, intent(out) :: RefAll, AltAll
    ! integer(int32), dimension(:,:,:), allocatable, intent(out) :: SequenceData

    real(real64), allocatable, dimension(:) :: quality
    integer(int32), dimension(:), allocatable :: position
    integer(int32), optional :: SnpUsed,nSnpIn,StartSnp,EndSnp,nIndivIn
    integer(int32) :: nSnp,pos, nIndiv
    integer(int32) :: i,j, k


    if (present(nSnpIn)) then
      nSnp = nSnpIn
    else
      write(*,*) "nSnp required at the moment."
      stop
    end if

    if (present(nIndivIn)) then
      nIndiv = nIndivIn
    else
      !Work out nIndiv
      write(*,*) "nIndiv required at the moment."
      stop
    end if

    ! allocate(SequenceData(nIndiv, SnpUsed, 2))
    allocate(RefAll(nIndiv, SnpUsed))
    allocate(AltAll(nIndiv, SnpUsed))
    allocate(position(SnpUsed))
    allocate(quality(SnpUsed))
    allocate(Ids(nIndiv))
    allocate(dumE(5+2*nIndiv))
    allocate(dumC(5+nIndiv))

    read(ReadsFileUnit, *) dumC
    ! write(*,"(5A)") "STUFF", trim(dumC(1)), trim(dumC(2)), trim(dumC(3)), "ENDSTUFF"
    do i =1, nIndiv
      write(Ids(i), *) dumC(i+5)
    end do

    pos=1
    do j = 1, nSnp
      read(ReadsFileUnit, *) dumE
      if ((j >= StartSnp).and.(j <= EndSnp)) then
        read(dumE(2), *) position(pos)
        read(dumE(5), *) quality(pos)

        !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED (nIndiv,j,pos,dumE,RefAll,AltAll,position)
        do i = 1, nIndiv
          k = i*2+4
          !write(*,'(4i5,1x,1a5)'),i,j,position(j),k,dumE(k)
          ! read(dumE(k), *) SequenceData(i,pos, 1)
          ! read(dumE(k+1), *) SequenceData(i, pos, 2)
          read(dumE(k), *) RefAll(i, pos)
          read(dumE(k+1), *) AltAll(i, pos)
          !write(*,'(5i20,2i10)'),i,j,pos,position(pos),k,SequenceData(i, pos,:)
        end do
        !$OMP END PARALLEL DO
        pos=pos+1
      endif
    end do
  end subroutine readVCF


  !#############################################################################################################################################################################################################################
  subroutine ReadGenos(GenoFileUnit)
  use GlobalPedigree
  use Global
  use alphaimputeinmod
  implicit none

  integer, intent(inout) :: GenoFileUnit
  integer :: i,j
  integer,allocatable,dimension(:) :: temp
  type(AlphaImputeInput), pointer :: inputParams
  logical :: opened, named
  character(len=300) :: GenoFile

  inputParams => defaultInput

  allocate(temp(inputParams%nSnp))

  if (allocated(Genos)) then
    deallocate(Genos)
  endif
  allocate(Genos(0:nAnisG,inputParams%nsnp))
  Genos(0,:)=9

  inquire(unit=GenoFileUnit, opened=opened, named=named, name=GenoFile)
  if (.NOT. opened .and. named) then
    open(unit=GenoFileUnit, file=GenoFile, status='unknown')
  else if (.NOT. named) then
    ! write(0, *) "ERROR - Something went wrong when trying to read the file of pre-phased data"
    open(newunit=GenoFileUnit, file=inputParams%GenotypeFile)
  end if

  do i=1,nAnisG
    read (GenoFileUnit,*) GenotypeId(i),Temp(:)
    do j=1,inputParams%nsnp
      if ((Temp(j)<0).or.(Temp(j)>2)) Temp(j)=9
    enddo
    Genos(i,:)=Temp(:)
  enddo
  close(GenoFileUnit)
  deallocate(temp)
  end subroutine ReadGenos

  !#############################################################################################################################################################################################################################
  subroutine ReadPlink(GenoFileUnit)
  use GlobalPedigree
  use Global
  use alphaimputeinmod
  implicit none

  integer, intent(inout) :: GenoFileUnit
  integer :: i,j
  character(len=2),allocatable,dimension(:) :: temp
  type(AlphaImputeInput), pointer :: inputParams
  logical :: opened, named

  character(len=300) :: dumC
  character(len=300) :: GenoFile

  inputParams => defaultInput

  allocate(temp(inputParams%nSnp))

  if (allocated(Genos)) then
    deallocate(Genos)
  endif
  allocate(Genos(0:nAnisG,inputParams%nsnp))
  Genos(0,:)=9

  inquire(unit=GenoFileUnit, opened=opened, named=named, name=GenoFile)
  if (.NOT. opened .and. named) then
    open(unit=GenoFileUnit, file=GenoFile, status='unknown')
  else if (.NOT. named) then
    open(newunit=GenoFileUnit, file=inputParams%GenotypeFile)
  end if

  ! Skip the first line as is the head
  read (GenoFileUnit,*) dumC
  do i=1,nAnisG
    read (GenoFileUnit,*) dumC,GenotypeId(i),dumC,dumC,dumC,dumC,Temp(:)
    do j=1,inputParams%nsnp
      if (Temp(j)=='NA') then
        Temp(j) = '9'
      else if (Temp(j)/='NA') then
      end if
    enddo
    read(Temp(:),*) Genos(i,:)
  enddo
  close(GenoFileUnit)
  deallocate(temp)

  end subroutine ReadPlink

!---------------------------------------------------------------------------
! DESCRIPTION:
!> @brief      Count the number of words in a line
!
!> @details    Count the number of words in a line separated by spaces
!
!> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
!
!> @date       Nov 16, 2016
!
! PARAMETERS:
!> @param[in] Line
!? @return    Number of words in the line
!---------------------------------------------------------------------------
  integer function nWords(line)
  implicit none
  character :: line*(*)
  logical :: back
  integer :: length, k

  back = .true.
  length = len_trim(line)

  k = index(line(1:length), ' ', back)
  if (k == 0) then
      nWords = 0
      return
  end if

  nWords = 1
  do
      ! starting with the right most blank space,
      ! look for the next non-space character down
      ! indicating there is another item in the line
      do
          if (k <= 0) exit
          if (line(k:k) == ' ') then
              k = k - 1
              cycle
          else
              nWords = nWords + 1
              exit
          end if
      end do

      ! once a non-space character is found,
      ! skip all adjacent non-space character
      do
          if ( k<=0 ) exit
          if (line(k:k) /= ' ') then
              k = k - 1
              cycle
          end if
          exit
      end do
      if (k <= 0) exit
  end do
  end function nWords


end module InputMod