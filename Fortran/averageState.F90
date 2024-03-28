subroutine averageState(vL, vR, vAvg)

#include "definition.h"

   implicit none
   real, dimension(NUMB_VAR), intent(IN) :: vL, vR !prim vars
   real, dimension(NUMB_VAR), intent(OUT) :: vAvg  !average state

   ! STUDENTS: PLEASE FINISH THIS SIMPLE AVERAGING SCHEME
   vAvg = (1.0d0/2.0d0)*(vL + vR)

   return
end subroutine averageState
