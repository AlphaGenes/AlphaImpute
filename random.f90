module random

use Par_Zig_mod

implicit none

public :: RandomOrderPar

contains

!#############################################################################################################################################################################################################################

subroutine RandomOrderPar(order,n,thread)
use Par_Zig_mod

 implicit none

!     Generate a random ordering of the integers 1 ... n.

integer, INTENT(IN)  :: n, thread
integer, INTENT(OUT) :: order(n)
!double precision par_uni

!     Local variables

integer :: i, j, k
double precision    :: wk

do i = 1, n
  order(i) = i
end do

!     Starting at the end, swap the current last indicator with one
!     randomly chosen from those preceeding it.

do i = n, 2, -1
  if (i<0) print*, i
  wk=par_uni(thread)
  j = 1 + i * wk
  if (j < i) then
    k = order(i)
    order(i) = order(j)
    order(j) = k
  end if
end do

RETURN
end subroutine RandomOrderPar

end module random
