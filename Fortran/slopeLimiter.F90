module slopeLimiter
  
#include "definition.h"
  
  !use grid_data


contains

  subroutine minmod(a,b,delta)
    implicit none
    real, intent(IN) :: a, b
    real, intent(OUT) :: delta

    delta = 0.5 * (sign(1.0,a) + sign(1.0,b))*min(abs(a),abs(b))

    return
  end subroutine minmod

  
  subroutine mc(a,b,delta)
    implicit none
    real, intent(IN) :: a, b
    real, intent(OUT) :: delta

    ! STUDENTS: PLEASE FINISH THIS MC LIMITER
    delta = (sign(1.0, a) + sign(1.0, b)) * min(min(abs(a), abs(b)),0.25d0*abs(a + b))
    
    return
  end subroutine mc

  
  subroutine vanLeer(a,b,delta)
    implicit none
    real, intent(IN) :: a, b
    real, intent(OUT) :: delta

    ! STUDENTS: PLEASE FINISH THIS VAN LEER'S LIMITER
    if (a * b <= 0) then
      delta = 0.0d0
    else
      delta = (2.0d0*a*b)/(a + b)
    end if
    
    return
  end subroutine vanLeer
  
end module slopeLimiter
