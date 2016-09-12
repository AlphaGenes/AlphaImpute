MODULE Utils
IMPLICIT NONE

CONTAINS

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief      Return the number of lines in a file
  !
  !> @details    Return the number of lines in a file
  !
  !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
  !
  !> @date       July 25, 2016
  !
  ! PARAMETERS:
  !> @param[in]  FileName  File name
  !> @return     Number of lines in the file
  !---------------------------------------------------------------------------
  FUNCTION CountLines(FileName) result(nLines)
    use ISO_Fortran_Env
    implicit none

    character(len=*), intent(in) :: FileName
    integer                         :: nLines

    integer :: f, UInputs
    character(len=300) :: dumC

    nLines=0
    UInputs = 111
    open (unit=UInputs,file=trim(FileName),status="old")
    do
      read (UInputs,*,iostat=f) dumC
      nLines=nLines+1
      if (f/=0) then
        nLines=nLines-1
        exit
      endif
    enddo
    close(UInputs)
  END FUNCTION CountLines

!######################################################################
SUBROUTINE RemoveInformationGameteSegment(gamete,nStart,nStop)
! Remove the allele information of a gamete segment

USE GlobalVariablesHmmMaCH

INTEGER,INTENT(IN) :: nStart, nStop
INTEGER(KIND=1),INTENT(OUT),DIMENSION(:) :: gamete

gamete(nStart:nStop) = ALLELE_MISSING

END SUBROUTINE RemoveInformationGameteSegment

!######################################################################
SUBROUTINE RemoveInformationGamete(gamete)
! Remove the allele information of a gamete

USE GlobalVariablesHmmMaCH

INTEGER(KIND=1),INTENT(OUT),DIMENSION(:) :: gamete

gamete(:) = ALLELE_MISSING

END SUBROUTINE RemoveInformationGamete

!######################################################################
SUBROUTINE RemoveAlleleInformationIndividual(ToWhom)
! Remove the genotype information of an individual

USE GlobalVariablesHmmMaCH

INTEGER, INTENT(IN) :: ToWhom

call RemoveAlleleInformationSegment(ToWhom, 1, nSnpHmm)

END SUBROUTINE RemoveAlleleInformationIndividual

!######################################################################
SUBROUTINE RemoveAlleleInformationSegment(ToWhom,nStart,nStop)
! Remove the genotype information of an individual

USE GlobalVariablesHmmMaCH

INTEGER, INTENT(IN) :: ToWhom, nStart, nStop

call RemoveInformationGamete(PhaseHmmMaCH(ToWhom,nStart:nStop,1))
call RemoveInformationGamete(PhaseHmmMaCH(ToWhom,nStart:nStop,2))

END SUBROUTINE RemoveAlleleInformationSegment

!######################################################################
SUBROUTINE RemoveGenotypeInformationIndividual(ToWhom)
! Remove the genotype information of an individual

USE GlobalVariablesHmmMaCH

INTEGER, INTENT(IN) :: ToWhom

call RemoveGenotypeInformationIndividualSegment(ToWhom,1,nSnpHmm)

END SUBROUTINE RemoveGenotypeInformationIndividual

!######################################################################
SUBROUTINE RemoveGenotypeInformationIndividualSegment(ToWhom,nStart,nStop)
! Remove the genotype information of an individual

USE GlobalVariablesHmmMaCH

INTEGER, INTENT(IN) :: ToWhom,nStart,nStop

GenosHmmMaCH(ToWhom,nStart:nStop) = MISSING

END SUBROUTINE RemoveGenotypeInformationIndividualSegment

!######################################################################
FUNCTION CountPhasedGametes RESULT( gametesPhased )

USE Global
USE GlobalVariablesHmmMaCH

INTEGER :: gametesPhased, ind
INTEGER :: CountPhasedAlleles
INTEGER(KIND=1), ALLOCATABLE :: gamete(:)

allocate(gamete(nSnpHmm))
gametesPhased=0
do ind=1,nAnisP
    gamete=ImputePhase(ind,:,1)
    if (float(count(gamete(:)==1 .OR. gamete(:)==0))/nSnpHmm >= imputedThreshold/100.0) then
        gametesPhased=gametesPhased+1
    endif

    gamete=ImputePhase(ind,:,2)
    if (float(count(gamete(:)==1 .OR. gamete(:)==0))/nSnpHmm >= imputedThreshold/100.0) then
        gametesPhased=gametesPhased+1
    endif
enddo
RETURN
END FUNCTION CountPhasedGametes

!######################################################################
FUNCTION CountPhasedAlleles(gamete) RESULT( allelesPhased )

USE Global
USE GlobalVariablesHmmMaCH

INTEGER :: allelesPhased
INTEGER(KIND=1),INTENT(IN) :: gamete(:)

allelesPhased = count(gamete(:)==1 .OR. gamete(:)==0)


RETURN
END FUNCTION CountPhasedAlleles

!######################################################################
FUNCTION CountMissingAlleles(gamete) RESULT( allelesMissing )

USE Global
USE GlobalVariablesHmmMaCH

INTEGER :: allelesMissing
INTEGER(KIND=1),INTENT(IN) :: gamete(:)

allelesMissing = nSnpHmm - CountPhasedAlleles(gamete)

RETURN
END FUNCTION CountMissingAlleles

!######################################################################
FUNCTION CountMissingAllelesByGametes(gamete1,gamete2) RESULT( allelesMissing )
! TODO: THIS SUBROUTINE IS NOT WORKING PROPERLY. ALLELESMISSING GETS SHIT

USE Global
USE GlobalVariablesHmmMaCH

INTEGER :: allelesMissing
INTEGER(KIND=1),INTENT(IN) :: gamete1(:),gamete2(:)

allelesMissing = nSnpHmm - CountGenotypedAllelesByGametes(gamete1,gamete2)

RETURN
END FUNCTION CountMissingAllelesByGametes

!######################################################################
FUNCTION CountGenotypedAllelesByGametes(gamete1, gamete2) RESULT( allelesGenotyped )

USE Global
USE GlobalVariablesHmmMaCH

INTEGER(KIND=1),INTENT(IN) :: gamete1(:), gamete2(:)

! Local variables
INTEGER :: allelesGenotyped
INTEGER :: i


allelesGenotyped = count((gamete1(:)==1 .OR. gamete1(:)==0) .AND. (gamete2(:)==1 .OR. gamete2(:)==0))

! do i=1,size(gamete1)
!     if ((gamete1(i)==1 .OR. gamete1(i)==0) .AND. (gamete2(i)==1 .OR. gamete2(i)==0)) then
!         allelesGenotyped = allelesGenotyped + 1
!     endif
! enddo

RETURN
END FUNCTION CountGenotypedAllelesByGametes

!######################################################################
FUNCTION CountGenotypedGenotypesByChromosome(chromosome) RESULT( genotypesMissing )

USE Global
USE GlobalVariablesHmmMaCH

INTEGER,INTENT(IN) :: chromosome(:)

! Local variables
INTEGER :: genotypesMissing

genotypesMissing = count(chromosome(:)==0 .or. chromosome(:)==2 .or. chromosome(:)==1)

RETURN
END FUNCTION CountGenotypedGenotypesByChromosome

!######################################################################
FUNCTION CountMissingdGenotypesByChromosome(chromosome) RESULT( genotypesMissing )

USE Global
USE GlobalVariablesHmmMaCH

INTEGER,INTENT(IN) :: chromosome(:)

! Local variables
INTEGER :: genotypesMissing

genotypesMissing = nSnpHmm - CountGenotypedGenotypesByChromosome(chromosome)

RETURN
END FUNCTION CountMissingdGenotypesByChromosome

!######################################################################

END MODULE
