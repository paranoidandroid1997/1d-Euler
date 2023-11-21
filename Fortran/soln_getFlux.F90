subroutine soln_getFlux()

#include "definition.h"  

  use grid_data
  use sim_data

  implicit none
  integer :: i

  if (sim_riemann == 'hll') then
     do i = gr_ibeg, gr_iend+1
        call hll(gr_vR(DENS_VAR:GAME_VAR,i-1),&
                 gr_vL(DENS_VAR:GAME_VAR,i  ),&
                 gr_flux(DENS_VAR:ENER_VAR,i))
     end do

  elseif (sim_riemann == 'roe') then
     do i = gr_ibeg, gr_iend+1
        call roe(gr_vR(DENS_VAR:GAME_VAR,i-1),&
                 gr_vL(DENS_VAR:GAME_VAR,i  ),&
                 gr_flux(DENS_VAR:ENER_VAR,i))
     end do
  endif

  return
end subroutine soln_getFlux
