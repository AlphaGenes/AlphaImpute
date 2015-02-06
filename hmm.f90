subroutine UpdateThetas

use GlobalVariablesHmmMaCH
implicit none

integer :: i, BaseCount=1, BaseIntervals=0
double precision :: BaseRates, BaseCrossovers=1, BaseLength=0

double precision :: Scale

Scale=1.0/(individuals*2)

! First we estimate a base line rate to be applied to intervals with
! 0 or 1 observed "crossovers"
do i=1,nSnpHmm-1
    if (Crossovers(i)<=1) then
        BaseCount=BaseCount+Crossovers(i)
        BaseIntervals=BaseIntervals+1
    endif
enddo

if (BaseIntervals==0) then
    BaseRates=BaseCount*Scale
else
    BaseRates=BaseCount*Scale/BaseIntervals
endif

! Then we update the rate for each interval using either the number
! of observed crossovers (if > 1) or the baseline rate
do i=1,nSnpHmm
    if (Crossovers(i)>1) then
        Thetas(i)=Crossovers(i)*Scale
    else
        Thetas(i)=BaseRates
    endif
enddo

end subroutine UpdateThetas

!######################################################################
subroutine UpdateErrorRate(rate)
! Group markers into those with low error rates, which are estimated
! as a group, and those with high error rates, which are estimated 
! individually

use GlobalVariablesHmmMaCH
implicit none

double precision,intent(out) :: rate

! Local variables
integer :: i,matches=0,mismatches=0,uncertain=0

rate=0.0
do i=1,nSnpHmm
    if (ErrorMismatches(i)<=2) then
        matches=matches+ErrorMatches(i)
        mismatches=mismatches+ErrorMismatches(i)
        uncertain=uncertain+ErrorUncertainty(i)
    else
        call SetPenetrance(i)
    endif
enddo

call UpdateError(matches, mismatches, uncertain, rate)

do i=1,nSnpHmm
    if (ErrorMismatches(i)<=2) call SetPenetrance(rate)
enddo

end subroutine UpdateErrorRate

!######################################################################
subroutine SetPenetrance(marker)

use GlobalVariablesHmmMaCH
implicit none

integer,intent(in) :: marker

! Local variables
double precision :: Err

Err=Epsilon(marker)

Penetrance(marker,0,0)=(1.0-Err)**2
Penetrance(marker,0,1)=2.0*(1.0-Err)*Err
Penetrance(marker,0,2)=Err**2
Penetrance(marker,1,0)=(1.0-Err)*Err
Penetrance(marker,1,1)=((1.0-Err)**2)+(Err**2)
Penetrance(marker,1,2)=(1.0-Err)*Err
Penetrance(marker,2,0)=Err**2   
Penetrance(marker,2,1)=2.0*(1.0-Err)*Err
Penetrance(marker,2,2)=(1.0-Err)**2   

end subroutine SetPenetrance


!######################################################################
subroutine TotalCrossovers(Total)

use GlobalVariablesHmmMaCH
implicit none

integer,intent(out) :: Total

! Local variables
integer :: i

Total=0
do i=1,nSnpHmm-1
    Total=Total+Crossovers(i)
enddo

end subroutine TotalCrossovers

!######################################################################
subroutine ResetCrossovers

use GlobalVariablesHmmMaCH
implicit none

Crossovers(:)=0
call ResetErrors

enddo

end subroutine ResetCrossovers

!######################################################################
subroutine UpdateError(matches, mismatches, uncertain, rate)

use GlobalVariablesHmmMaCH
implicit none

integer,intent(in) :: matches,mismatches,uncertain
double precision, intent(out) :: rate

! Local variables
integer :: matches, mismatches, uncertain
double precision :: previous=0.0, ratio

rate=0.0    ! Just in case...
if (matches+mismatches>0) then
    rate=mismatches/dble(matches+mismatches)
    if (uncertain>0) then
        do while((rate>1e-10)).and.(abs(rate-previous) > rate*1e-4))
            ratio=rate*rate/(rate*rate+(1.0-rate)*(1.0-rate))
            previous=rate
            rate=(mismatches+ratio*uncertain*2.0)&
                 /(matches+mismatches+uncertain*2)
        enddo
    endif
else if (uncertain>0)
    rate=0.0
endif

end subroutine UpdateErrors

!######################################################################
subroutine ResetErrors

use GlobalVariablesHmmMaCH
implicit none

ErrorUncertainty(:)=0
ErrorMatches(:)=0
ErrorMismatches(:)=0

end subroutine ResetErrors











































