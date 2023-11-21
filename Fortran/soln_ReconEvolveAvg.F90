subroutine soln_ReconEvolveAvg(dt)

#include "definition.h"  

  use grid_data
  use sim_data

  implicit none
  real, intent(IN) :: dt

  ! conservative left and right states
  real, dimension(NSYS_VAR,gr_imax) :: uL, uR


  call soln_reconstruct(dt)
  call soln_getFlux()

  return
end subroutine soln_ReconEvolveAvg
