!-----------------------------------------------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-----------------------------------------------------------------------------------------------------------------------
!
! MODULE: PhaseRounds
!
!> @file        PhaseRounds.f90
!
! DESCRIPTION: 
!> @brief       Module to define and read phasing rounds and their cores
!>
!> @details     This MODULE includes routines to read the information contain in the different phasing rounds located
!>              in the folder 'Phasing'
!
!> @author      Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
!
!> @date        July 15, 2016
!
!> @version     0.0.1 (alpha)
!
! REVISION HISTORY:
! 2016.07.15  RAntolin - Initial Version
!
!-----------------------------------------------------------------------------------------------------------------------
MODULE PhaseRounds
  use ISO_Fortran_Env
  implicit none
  PRIVATE

  ! private CountLines, getCoreIndexes_Path, getCoreIndexes_NoPath, ReadCores

  ! PUBLIC getCoreIndexes
  ! INTERFACE getCoreIndexes
  !   SUBROUTINE getFileNameCoreIndex_Path(PhasePath, phaseInternal) result(FileName)
  !   integer, intent(in) :: phaseInternal
  !   character (len=300), intent(in) :: PhasePath
  !   character(len=1000) :: FileName
  !   END SUBROUTINE getFileNameCoreIndex_Path

  !   SUBROUTINE getFileNameCoreIndex_NoPath(nPhaseInternal) result(FileName)
  !   integer, intent(in) :: phaseInternal
  !   character(len=1000) :: FileName
  !   END SUBROUTINE getFileNameCoreIndex_NoPath
  ! END interface getCoreIndexes

  PUBLIC getFileNameCoreIndex, getFileNameFinalPhase, getFileNameHapLib
  PUBLIC ReadCores, ReadPhased
  ! PUBLIC destroy

  INTERFACE getFileNameCoreIndex
    MODULE PROCEDURE getFileNameCoreIndex_NoPath, getFileNameCoreIndex_Path
  END INTERFACE getFileNameCoreIndex

  INTERFACE getFileNameFinalPhase
    MODULE PROCEDURE getFileNameFinalPhase_NoPath, getFileNameFinalPhase_Path
  END INTERFACE getFileNameFinalPhase

  INTERFACE getFileNameHapLib
    MODULE PROCEDURE getFileNameHapLib_NoPath, getFileNameHapLib_Path
  END INTERFACE getFileNameHapLib

  TYPE, PUBLIC :: CoreIndex
    ! PUBLIC
    integer(kind = 2), allocatable, dimension(:) :: StartSnp
    integer(kind = 2), allocatable, dimension(:) :: EndSnp
    integer(kind = 2)                            :: nCores
  CONTAINS
    ! PRIVATE
    final :: destroy_CoreIndex
  END TYPE CoreIndex

CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: 
!> @brief      Initialise new core index
!
!> @details    Initialise new core index
!
!> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
!
!> @date       Julyy 25, 2016
!
! PARAMETERS:
!> @param[inout]  CoreI  CoreIndex
!---------------------------------------------------------------------------  
  FUNCTION newCoreIndex(nCores) result(CoreI) 
    integer, intent(in) :: nCores
    type(CoreIndex)     :: CoreI

    integer :: UInputs

    CoreI%nCores = nCores
    allocate(CoreI%StartSnp(CoreI%nCores))
    allocate(CoreI%EndSnp(CoreI%nCores))

  END FUNCTION newCoreIndex

!---------------------------------------------------------------------------  
! DESCRIPTION: 
!> @brief      Deallocate core information
!
!> @details    Deallocate core information
!
!> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
!
!> @date       July 15, 2016
!
! PARAMETERS:
!> @param[out] definition Core
!---------------------------------------------------------------------------  
  SUBROUTINE destroy_CoreIndex(this)
    type(CoreIndex) :: this

    if (allocated(this%StartSnp)) then
      deallocate(this%StartSnp)
    end if
    if (allocated(this%EndSnp)) then
      deallocate(this%EndSnp)
    end if
  end SUBROUTINE destroy_CoreIndex

!---------------------------------------------------------------------------  
! DESCRIPTION: 
!> @brief      Read core information
!
!> @details    Read the start and end snp for each core in the phasing round
!
!> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
!
!> @date       Julyy 15, 2016
!
! PARAMETERS:
!> @param[in]  FileName  File containing the cores information
!> @param[out] CoreI     Core information (start snp, end snp)
!---------------------------------------------------------------------------  
  FUNCTION ReadCores(FileName) result(CoreI)
    use Utils

    character(len=1000), intent(in) :: FileName
    type(CoreIndex)                 :: CoreI

    integer :: i
    integer :: UInputs
    character(len=20) :: dum

    CoreI = newCoreIndex(CountLines(FileName))

    ! Get core information from file
    open (unit=UInputs,file=trim(FileName),status="old")
    do i=1,CoreI%nCores
        read (UInputs,*) dum, CoreI%StartSnp(i), CoreI%EndSnp(i)
    end do
    close(UInputs)
  END FUNCTION ReadCores

!---------------------------------------------------------------------------  
! DESCRIPTION: 
!> @brief      Return file of Core indexes
!
!> @details    Return the correct path of the file containing the core indexes 
!>             of a given intenal phase if a folder with phased data is provided
!
!> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
!
!> @date       July 25, 2016
!
! PARAMETERS:
!> @param[in] PhasePath       Path of phased data
!> @param[in] phaseInternal   Number of a internal phase
!> @return    File name
!---------------------------------------------------------------------------  
  FUNCTION getFileNameCoreIndex_Path(PhasePath, phaseInternal) result(FileName)
    integer, intent(in)             :: phaseInternal
    character(len=300), intent(in)  :: PhasePath
    character(len=1000)             :: FileName

#ifdef OS_UNIX
      write (FileName,'(a,"Phase",i0,"/PhasingResults/CoreIndex.txt")') trim(PhasePath),phaseInternal
#else
      write (FileName,'(a,"Phase",i0,"\PhasingResults\CoreIndex.txt")') trim(PhasePath),phaseInternal
#endif
  END FUNCTION getFileNameCoreIndex_Path

 !---------------------------------------------------------------------------  
 ! DESCRIPTION: 
 !> @brief      Return file of Core indexes
 !
 !> @details    Return the correct path of the file containing the core indexes 
 !>             of a given intenal phase
 !
 !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
 !
 !> @date       July 25, 2016
 !
 ! PARAMETERS:
 !> @param[in] phaseInternal   Number of a internal phase
 !> @return    File name
 !---------------------------------------------------------------------------  
  FUNCTION getFileNameCoreIndex_NoPath(phaseInternal) result(FileName)
    integer, intent(in) :: phaseInternal
    character(len=1000) :: FileName

#ifdef OS_UNIX
      write (FileName,'("./Phasing/Phase",i0,"/PhasingResults/CoreIndex.txt")')phaseInternal
#else
      write (FileName,'(".\Phasing\Phase",i0,"\PhasingResults\CoreIndex.txt")')phaseInternal
#endif
  END FUNCTION getFileNameCoreIndex_NoPath

!---------------------------------------------------------------------------  
! DESCRIPTION: 
!> @brief      Return file of the final phase of a internal phasing
!
!> @details    Return the correct path of the file containing the final phase 
!>             of a given intenal phase if a folder with phased data is provided
!
!> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
!
!> @date       July 25, 2016
!
! PARAMETERS:
!> @param[in] PhasePath       Path of phased data
!> @param[in] phaseInternal   Number of a internal phase
!> @return    File name
!---------------------------------------------------------------------------  
  FUNCTION getFileNameFinalPhase_Path(PhasePath, phaseInternal) result(FileName)
    integer, intent(in)             :: phaseInternal
    character(len=300), intent(in)  :: PhasePath
    character(len=1000)             :: FileName

#ifdef OS_UNIX
      write (FileName,'(a,"Phase",i0,"/PhasingResults/FinalPhase.txt")') trim(PhasePath),phaseInternal
#else
      write (FileName,'(a,"Phase",i0,"\PhasingResults\FinalPhase.txt")') trim(PhasePath),phaseInternal
#endif
  END FUNCTION getFileNameFinalPhase_Path

!---------------------------------------------------------------------------
! DESCRIPTION:
!> @brief      Return file of the final phase of a internal phasing
!
!> @details    Return the correct path of the file containing the final phase 
 !>             of a given intenal phase
!
!> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
!
!> @date       July 25, 2016
!
! PARAMETERS:
!> @param[in] PhasePath       Path of phased data
!> @param[in] phaseInternal   Number of a internal phase
!> @return    File name
!---------------------------------------------------------------------------
  FUNCTION getFileNameFinalPhase_NoPath(phaseInternal) result(FileName)
    integer, intent(in)            :: phaseInternal
    character(len=1000)            :: FileName

#ifdef OS_UNIX
      write (FileName,'("./Phasing/Phase",i0,"/PhasingResults/FinalPhase.txt")') phaseInternal
#else
      write (FileName,'(".\Phasing\Phase",i0,"\PhasingResults\FinalPhase.txt")')phaseInternal
#endif
  END FUNCTION getFileNameFinalPhase_NoPath

!---------------------------------------------------------------------------
! DESCRIPTION:
!> @brief      Return file of the haplotype library of a internal phasing
!
!> @details    Return the correct path of the file containing the haplotype library
!>             of a given intenal phase and core
!
!> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
!
!> @date       July 25, 2016
!
! PARAMETERS:
!> @param[in] PhasePath       Path of phased data
!> @param[in] phaseInternal   Number of a internal phase
!> @return    File name
!---------------------------------------------------------------------------
  FUNCTION getFileNameHapLib_Path(PhasePath, phaseInternal, core) result(FileName)
    integer, intent(in)             :: phaseInternal, core
    character(len=300), intent(in)  :: PhasePath
    character(len=1000)             :: FileName

#ifdef OS_UNIX
    write (FileName,'(a,"Phase",i0,"/PhasingResults/HaplotypeLibrary/HapLib",i0,".bin")') trim(PhasePath), phaseInternal, core
#else
    write (FileName,'(a,"Phase",i0,"\PhasingResults\HaplotypeLibrary\HapLib",i0,".bin")') trim(PhasePath), phaseInternal, core
#endif
  END FUNCTION getFileNameHapLib_Path

!---------------------------------------------------------------------------
! DESCRIPTION:
!> @brief      Return file of the haplotype library of a internal phasing
!
!> @details    Return the correct path of the file containing the haplotype library
!>             of a given intenal phase and core
!
!> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
!
!> @date       July 25, 2016
!
! PARAMETERS:
!> @param[in] PhasePath       Path of phased data
!> @param[in] phaseInternal   Number of a internal phase
!> @return    File name
!---------------------------------------------------------------------------
  FUNCTION getFileNameHapLib_NoPath(phaseInternal, core) result(FileName)
    integer, intent(in)            :: phaseInternal, core
    character(len=1000)            :: FileName

#ifdef OS_UNIX
    write (FileName,'("./Phasing/Phase",i0,"/PhasingResults/HaplotypeLibrary/HapLib",i0,".bin")') phaseInternal, core
#else
    write (FileName,'(".\Phasing\Phase",i0,"\PhasingResults\HaplotypeLibrary\HapLib",i0,".bin")') phaseInternal, core
#endif
  END FUNCTION getFileNameHapLib_NoPath


  SUBROUTINE ReadPhased(nAnis, nAnisPed, FileName, Ids, PhaseHD, PosHD)
    integer, intent(in)                     :: nAnis
    integer, intent(in)                     :: nAnisPed
    character(len=*), intent(in)            :: FileName
    character*(20), intent(in)              :: Ids(nAnisPed)
    integer, dimension (:), intent(out)     :: PosHD
    integer, dimension (:,:,:), intent(out) :: PhaseHD

    integer :: i, j
    integer :: UPhased
    character(len=300) :: dumC

  ! Get phase information from file
    open (unit=UPhased,file=trim(FileName),status="old")
    do i=1,nAnis
        read (UPhased,*) dumC,PhaseHD(i,:,1)
        read (UPhased,*) dumC,PhaseHD(i,:,2)
        ! Match HD phase information with individuals
        do j=1,nAnisPed
            if (trim(dumC)==trim(Ids(j))) then
                PosHD(j)=i
                exit
            endif
        enddo
    enddo
    close(UPhased)
  END SUBROUTINE ReadPhased

END MODULE PhaseRounds

