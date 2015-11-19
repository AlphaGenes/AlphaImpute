MODULE Utils
IMPLICIT NONE

CONTAINS

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

GenosHmmMaCH(ToWhom,:) = MISSING

END SUBROUTINE RemoveGenotypeInformationIndividual

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
    if (float(count(gamete(:)==1 .OR. gamete(:)==0))/nSnpHmm >= WellPhasedThresh/100.0) then
        gametesPhased=gametesPhased+1
    endif

    gamete=ImputePhase(ind,:,2)
    if (float(count(gamete(:)==1 .OR. gamete(:)==0))/nSnpHmm >= WellPhasedThresh/100.0) then
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
FUNCTION CountGenotypedAllelesByGametes(gamete1, gamete2) RESULT( allelesGenotyped )

USE Global
USE GlobalVariablesHmmMaCH

INTEGER(KIND=1),INTENT(IN) :: gamete1(:), gamete2(:)

! Local variables
INTEGER :: allelesGenotyped
INTEGER :: i

do i=1,size(gamete1)
    if ((gamete1(i)==1 .OR. gamete1(i)==0) .AND. (gamete2(i)==1 .OR. gamete2(i)==0)) then
        allelesGenotyped = allelesGenotyped + 1
    endif
enddo

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
