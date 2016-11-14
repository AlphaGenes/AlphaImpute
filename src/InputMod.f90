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
  subroutine CountInData

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
      ! print*, 'boh'
      nAnisG=nAnisG-1
      exit
    endif
  enddo
  ! print *, 'hola', nAnisG

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

  end subroutine CountInData

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
  allocate(Ped(nAnisRawPedigree,3))
  allocate(Genos(0:nAnisG,inputParams%nsnp))
  allocate(GenderId(nAnisRawPedigree))
  allocate(GenderRaw(nAnisRawPedigree))

  Genos(0,:)=9

  ! Read the pedigree information
  rewind(inputParams%pedigreeFileUnit)
  do i=1,nAnisRawPedigree
      read(inputParams%pedigreeFileUnit,*) ped(i,:)
  enddo

  ! Read the genotype file
  if (inputParams%hmmoption /= RUN_HMM_NGS) then
      rewind(inputParams%genotypeFileUnit)
      if (inputParams%PlinkFormat) then
        call ReadGenos(inputParams%genotypeFileUnit)
      else
        call ReadPlink(inputParams%genotypeFileUnit)
      end if
  endif
  close(inputParams%pedigreeFileUnit)
  close(inputParams%genotypeFileUnit)

  ! Read the gender file if imputing the sex chromosome
  GenderRaw=9
  if (inputParams%SexOpt==1) then
      CountLinesGender=0
      do
          read (inputParams%GenderFileUnit,*,iostat=k) dumC
          CountLinesGender=CountLinesGender+1
          if (k/=0) then
              CountLinesGender=CountLinesGender-1
              exit
          endif
      enddo
      rewind(inputParams%GenderFileUnit)
      nAnisInGenderFile=CountLinesGender
      if (CountLinesGender/=nAnisRawPedigree) then
          print*, "Warning - number of lines in Gender file not the same as in pedigree file"
          stop
      endif
      do j=1,nAnisRawPedigree                 ! For each individual in the file
          read (inputParams%GenderFileUnit,*) dumC,GenCode
          if ((GenCode/=1).and.(GenCode/=2)) then
              print*, "Warning - Gender code incorrect for at least one animal"
              stop
          endif
          AnimalPresent=0
          do i=1,nAnisRawPedigree             ! For each individual in the pedigree
              if (trim(dumC)==trim(ped(i,1))) then
                  GenderId(i)=dumC
                  GenderRaw(i)=GenCode
                  AnimalPresent=1
                  exit
              endif
          enddo
          if (AnimalPresent==0) then
              print*, "Warning - Animal missing in gender file"
              stop
          endif
      enddo
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

  do i=1,nAnisG
    read (ReadsFileUnit,*) GenotypeId(i), ReferAlleleLine(:)
    read (ReadsFileUnit,*) GenotypeId(i), AlterAlleleLine(:)
    ReferAllele(i,:) = ReferAlleleLine
    AlterAllele(i,:) = AlterAlleleLine
    do j=1,inputParams%nsnp
      if (ReferAllele(i,j)>=MAX_READS_COUNT) ReferAllele(i,j)=MAX_READS_COUNT-1
      if (AlterAllele(i,j)>=MAX_READS_COUNT) AlterAllele(i,j)=MAX_READS_COUNT-1
      Reads(i,j)=AlterAllele(i,j)+ReferAllele(i,j)
    enddo
  enddo

#ifdef DEBUG
  write(0,*) "DEBUG: [ReadSeq] Sequence data read"
#endif

  close(ReadsFileUnit)

  end subroutine ReadSeq

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
    ! write(0, *) "ERROR - Something went wrong when trying to read the file of pre-phased data"
    open(newunit=GenoFileUnit, file=inputParams%GenotypeFile)
  end if

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


end module InputMod