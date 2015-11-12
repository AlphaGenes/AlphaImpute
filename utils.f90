!######################################################################
subroutine RemoveInformationGameteSegment(gamete,nStart,nStop)
! Remove the allele information of a gamete segment

use GlobalVariablesHmmMaCH

INTEGER,INTENT(IN) :: nStart, nStop
INTEGER,INTENT(OUT),ALLOCATABLE,DIMENSION(:) :: gamete

gamete(nStart:nStop) = ALLELE_MISSING

end subroutine RemoveInformationGameteSegment

!######################################################################
subroutine RemoveInformationGamete(gamete)
! Remove the allele information of a gamete

use GlobalVariablesHmmMaCH

INTEGER,INTENT(OUT),ALLOCATABLE,DIMENSION(:) :: gamete

gamete(:) = ALLELE_MISSING

end subroutine RemoveInformationGamete

!######################################################################
subroutine RemoveAlleleInformationIndividual(ToWhom)
! Remove the genotype information of an individual

use GlobalVariablesHmmMaCH

INTEGER, INTENT(IN) :: ToWhom

call RemoveAlleleInformationSegment(ToWhom, 1, nSnpHmm)

end subroutine RemoveAlleleInformationIndividual

!######################################################################
subroutine RemoveAlleleInformationSegment(ToWhom,nStart,nStop)
! Remove the genotype information of an individual

use GlobalVariablesHmmMaCH

INTEGER, INTENT(IN) :: ToWhom, nStart, nStop

call RemoveInformationGamete(PhaseHmmMaCH(ToWhom,:,1), nStart, nStop)
call RemoveInformationGamete(PhaseHmmMaCH(ToWhom,:,2), nStart, nStop)

end subroutine RemoveAlleleInformationSegment

!######################################################################
subroutine RemoveGenotypeInformationIndividual(ToWhom)
! Remove the genotype information of an individual

use GlobalVariablesHmmMaCH

INTEGER, INTENT(IN) :: ToWhom

GenosHmmMaCH(ToWhom,:) = MISSING

end subroutine RemoveGenotypeInformationIndividual

!######################################################################
FUNCTION CountPhasedGametes RESULT( gametesPhased )

use Global
use GlobalVariablesHmmMaCH

INTEGER :: gametesPhased, ind
INTEGER :: CountPhasedAlleles
INTEGER, ALLOCATABLE :: gamete(:)

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
end FUNCTION CountPhasedGametes

!######################################################################
FUNCTION CountPhasedAlleles(gamete) RESULT( allelesPhased )

use Global
use GlobalVariablesHmmMaCH

INTEGER :: allelesPhased
INTEGER,ALLOCATABLE :: gamete(:)

allelesPhased = count(gamete(:)==1 .OR. gamete(:)==0)
RETURN
end FUNCTION CountPhasedAlleles
