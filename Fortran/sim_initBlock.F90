subroutine sim_initBlock()

#include "definition.h"

    use sim_data
    use grid_data, only: gr_V, gr_U, gr_i0, gr_imax, gr_xCoord, gr_dx, gr_ngc, gr_xbeg
    use primconsflux, only: prim2cons

    implicit none

    integer :: i
    real :: ekin, eint

    ! generate x-coordinate
    do i = gr_i0, gr_imax
        gr_xCoord(i) = gr_xbeg + (real(i - gr_ngc) - 0.5)*gr_dx
    end do

    if (IC_type == "sim_blast2") then
        do i = gr_i0, gr_imax
            if (gr_xCoord(i) < sim_shockLoc) then
                gr_V(DENS_VAR, i) = sim_densL
                gr_V(VELX_VAR, i) = sim_velxL
                gr_V(PRES_VAR, i) = sim_presL
                print *, "iiiiiiiiiiiiiiiiiiiiiiiiiiii"
                print *, sim_shockLoc
                print *, sim_shockLoc2
            elseif (gr_xCoord(i) <= sim_shockLoc2) then
                gr_V(DENS_VAR, i) = sim_densM
                gr_V(VELX_VAR, i) = sim_velxM
                gr_V(PRES_VAR, i) = sim_presM
                print *, "jjjjjjjjjjjjjjjjjjjjjjjjjjj"
            else
                gr_V(DENS_VAR, i) = sim_densR
                gr_V(VELX_VAR, i) = sim_velxR
                gr_V(PRES_VAR, i) = sim_presR
                print *, "kkkkkkkkkkkkkkkkkkkkkkkkkkk"
            end if

            gr_V(GAMC_VAR, i) = sim_gamma
            gr_V(GAME_VAR, i) = sim_gamma
            gr_V(EINT_VAR, i) = gr_V(PRES_VAR, i)/(gr_V(GAME_VAR, i) - 1.)/gr_V(DENS_VAR, i)
        end do
    elseif (IC_type == "user") then
        do i = gr_i0, gr_imax
            if (gr_xCoord(i) < sim_shockLoc) then
                gr_V(DENS_VAR, i) = sim_densL
                gr_V(VELX_VAR, i) = sim_velxL
                gr_V(PRES_VAR, i) = sim_presL
            else
                gr_V(DENS_VAR, i) = sim_densR + 0.20d0*sin(5.0*gr_xCoord(i))
                gr_V(VELX_VAR, i) = sim_velxR
                gr_V(PRES_VAR, i) = sim_presR
            end if

            gr_V(GAMC_VAR, i) = sim_gamma
            gr_V(GAME_VAR, i) = sim_gamma
            gr_V(EINT_VAR, i) = gr_V(PRES_VAR, i)/(gr_V(GAME_VAR, i) - 1.)/gr_V(DENS_VAR, i)
        end do
    else
        do i = gr_i0, gr_imax
            if (gr_xCoord(i) < sim_shockLoc) then
                gr_V(DENS_VAR, i) = sim_densL
                gr_V(VELX_VAR, i) = sim_velxL
                gr_V(PRES_VAR, i) = sim_presL
            else
                gr_V(DENS_VAR, i) = sim_densR
                gr_V(VELX_VAR, i) = sim_velxR
                gr_V(PRES_VAR, i) = sim_presR
            end if

            gr_V(GAMC_VAR, i) = sim_gamma
            gr_V(GAME_VAR, i) = sim_gamma
            gr_V(EINT_VAR, i) = gr_V(PRES_VAR, i)/(gr_V(GAME_VAR, i) - 1.)/gr_V(DENS_VAR, i)
        end do
    end if

    ! also initialize conservative vars
    do i = gr_i0, gr_imax
        call prim2cons(gr_V(:, i), gr_U(DENS_VAR:ENER_VAR, i))
    end do

end subroutine sim_initBlock
