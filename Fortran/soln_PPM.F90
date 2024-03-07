subroutine soln_PPM(dt)
#include "definition.h"

    use grid_data
    use sim_data
    use slopeLimiter
    use eigensystem

    implicit none
    real, intent(IN) :: dt

    ! STUDENTS: IMPLEMENT THE THIRD-ORDER PPM SCHEME

    !! My code

    ! Matrix to hold polynomial coefficients for 3 reconstructions
    real, dimension(DENS_VAR:PRES_VAR, 0:2) :: C_MAT

    ! Vectors to hold A plus and min
    real, dimension(DENS_VAR:PRES_VAR) :: A_PLUS
    real, dimension(DENS_VAR:PRES_VAR) :: A_MIN

    real, dimension(NUMB_WAVE) :: lambda
    real, dimension(NSYS_VAR, NUMB_WAVE) :: reig, leig
    logical :: conservative
    real, dimension(NSYS_VAR) :: vecL, vecL2, vecR, vecR2, sigL, sigR
    integer :: kWaveNum
    real :: lambdaDtDx
    real, dimension(NUMB_VAR)  :: delL, delR
    real, dimension(NUMB_VAR, gr_i0:gr_imax) :: delVChar, delVPrim
    real, dimension(NUMB_WAVE) :: delW, delC1, delC2
    integer :: nVar
    integer :: i

    ! we need conservative eigenvectors
    conservative = .false.

    do i = gr_ibeg - 2, gr_iend + 2 ! 3 ghost cells !!!
        call eigenvalues(gr_V(DENS_VAR:GAME_VAR, i), lambda)
        call left_eigenvectors(gr_V(DENS_VAR:GAME_VAR, i), conservative, leig)
        call right_eigenvectors(gr_V(DENS_VAR:GAME_VAR, i), conservative, reig)

        ! Primitive Limiting
        if (.not. sim_charLimiting) then
            stop
            ! Char Limiting
        elseif (sim_charLimiting) then
            do kWaveNum = 1, NUMB_WAVE
                delL(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR, i) - gr_V(DENS_VAR:PRES_VAR, i - 1)
                delR(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR, i + 1) - gr_V(DENS_VAR:PRES_VAR, i)
                if (sim_limiter == 'minmod') then
                    call minmod(dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delL(DENS_VAR:PRES_VAR)), &
                              & dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delR(DENS_VAR:PRES_VAR)), &
                              & delW(kWaveNum))
                elseif (sim_limiter == 'vanLeer') then
                    call vanLeer(dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delL(DENS_VAR:PRES_VAR)), &
                               & dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delR(DENS_VAR:PRES_VAR)), &
                               & delW(kWaveNum))
                elseif (sim_limiter == 'mc') then
                    call mc(dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delL(DENS_VAR:PRES_VAR)), &
                    &dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delR(DENS_VAR:PRES_VAR)), delW(kWaveNum))
                end if
            end do

            do nvar = DENS_VAR, PRES_VAR
                delVChar(nvar, i) = dot_product(reig(nvar, 1:NUMB_WAVE), delW(1:NUMB_WAVE))
            end do
        end if
    end do

    do i = gr_ibeg - 1, gr_iend + 1

        call eigenvalues(gr_V(DENS_VAR:GAME_VAR, i), lambda)
        call left_eigenvectors(gr_V(DENS_VAR:GAME_VAR, i), conservative, leig)
        call right_eigenvectors(gr_V(DENS_VAR:GAME_VAR, i), conservative, reig)

        do nVar = DENS_VAR, PRES_VAR
            if (.not. sim_charLimiting) then
                stop ! not implemented yet
            else
        A_PLUS(nvar) = (1.0d0/2.0d0)*(gr_V(nVar, i) + gr_V(nVar, i + 1)) - (1.0d0/6.0d0)*(delVChar(nVar, i + 1) - delVChar(nVar, i))
         A_MIN(nvar) = (1.0d0/2.0d0)*(gr_V(nVar, i - 1) + gr_V(nVar, i)) - (1.0d0/6.0d0)*(delVChar(nVar, i) - delVChar(nVar, i - 1))
            end if

            ! Checks
            ! cond 1
            !if ((A_PLUS(nVar) - gr_V(nVar, i)) * (gr_V(nVar, i) - A_MIN(nVar)) <= 0.0d0) then
            if ((A_PLUS(nVar) - gr_V(nVar, i))*(gr_V(nVar, i) - A_MIN(nVar)) <= 0.0d0) then
                A_PLUS(nVar) = gr_V(nVar, i)
                A_MIN(nVar) = gr_V(nVar, i)

                C_MAT(nVar, 2) = 0.0d0
                C_MAT(nVar, 1) = 0.0d0
                C_MAT(nVar, 0) = gr_V(nvar, i)
            else
                !cond 2
                if (-(A_PLUS(nVar) - A_MIN(nVar))**2.0d0 > 6.0d0*(A_PLUS(nVar) - A_MIN(nVar))*(gr_V(nVar, i) - ((A_PLUS(nVar) + A_MIN(nVar))/2.0d0))) then
                    A_PLUS(nVar) = 3.0d0*gr_V(nVar, i) - 2.0d0*A_MIN(nVar)
                    ! cond 3
                elseif ((A_PLUS(nVar) - A_MIN(nVar))**2.0d0 < 6.0d0*(A_PLUS(nVar) - A_MIN(nVar))*(gr_V(nVar, i) - ((A_PLUS(nVar) + A_MIN(nVar))/2.0d0))) then
                    A_MIN(nVar) = 3.0d0*gr_V(nVar, i) - 2.0d0*A_PLUS(nVar)
                end if
                C_MAT(nVar, 2) = (6.0d0/(gr_dx**2.0d0))*((A_PLUS(nVar) + A_MIN(nvar))/(2.0d0) - gr_V(nvar, i))
                C_MAT(nVar, 1) = (1.0d0/gr_dx)*(A_PLUS(nVar) - A_MIN(nVar))
                C_MAT(nVar, 0) = gr_V(nvar, i) - (C_MAT(nvar, 2)/12.0d0)*gr_dx**2.0d0
            end if

        end do

        ! set the initial sum to be zero
        sigL(DENS_VAR:PRES_VAR) = 0.0d0
        sigR(DENS_VAR:PRES_VAR) = 0.0d0
        vecL(DENS_VAR:PRES_VAR) = 0.0d0
        vecL2(DENS_VAR:PRES_VAR) = 0.0d0
        vecR(DENS_VAR:PRES_VAR) = 0.0d0
        vecR2(DENS_VAR:PRES_VAR) = 0.0d0

        do kWaveNum = 1, NUMB_WAVE
            delC1(kWaveNum) = dot_product(leig(:, kWaveNum), C_MAT(:, 1))*gr_dx
            delC2(kWaveNum) = dot_product(leig(:, kWaveNum), C_MAT(:, 2))*gr_dx**2.0d0

            lambdaDtDx = lambda(kWaveNum)*dt/gr_dx

            if (sim_riemann == 'roe') then
                ! STUDENTS: PLEASE FINISH THIS ROE SOLVER CASE
                ! THIS SHOULDN'T BE LONGER THAN THE HLL CASE
                if (lambda(kWaveNum) > 0.0d0) then
                    vecR(DENS_VAR:PRES_VAR) = 0.5*(1.0 - lambdaDtDx)*reig(DENS_VAR:PRES_VAR, kWaveNum)*delC1(kWaveNum)
                    vecR2(DENS_VAR:PRES_VAR) = (1.0d0/4.0d0)*(1.0d0 - 2.0d0*lambdaDtDx + (4.0d0/3.0d0)*(lambdaDtDx)**2.0d0)*reig(DENS_VAR:PRES_VAR, kWaveNum) * delC2(kWaveNum)
                    sigR(DENS_VAR:PRES_VAR) = sigR(DENS_VAR:PRES_VAR) + vecR(DENS_VAR:PRES_VAR) + vecR2(DENS_VAR:PRES_VAR)
                else if (lambda(kWaveNum) < 0.0d0) then
                    vecL(DENS_VAR:PRES_VAR) = 0.5*(-1.0 - lambdaDtDx)*reig(DENS_VAR:PRES_VAR, kWaveNum)*delC1(kWaveNum)
                    vecL2(DENS_VAR:PRES_VAR) = (1.0d0/4.0d0)*(1.0d0 + 2.0d0*lambdaDtDx + (4.0d0/3.0d0)*(lambdaDtDx)**2.0d0)*reig(DENS_VAR:PRES_VAR, kWaveNum) * delC2(kWaveNum)
                    sigL(DENS_VAR:PRES_VAR) = sigL(DENS_VAR:PRES_VAR) + vecL(DENS_VAR:PRES_VAR) + vecL2(DENS_VAR:PRES_VAR)
                end if

            elseif (sim_riemann == 'hll') then
                vecR(DENS_VAR:PRES_VAR) = 0.5*(1.0 - lambdaDtDx)*reig(DENS_VAR:PRES_VAR, kWaveNum)*delC1(kWaveNum)
                vecR2(DENS_VAR:PRES_VAR) = (1.0d0/4.0d0)*(1.0d0 - 2.0d0*lambdaDtDx + (4.0d0/3.0d0)*(lambdaDtDx)**2.0d0)*reig(DENS_VAR:PRES_VAR, kWaveNum) * delC2(kWaveNum)
                sigR(DENS_VAR:PRES_VAR) = sigR(DENS_VAR:PRES_VAR) + vecR(DENS_VAR:PRES_VAR) + vecR2(DENS_VAR:PRES_VAR)

                vecL(DENS_VAR:PRES_VAR) = 0.5*(-1.0 - lambdaDtDx)*reig(DENS_VAR:PRES_VAR, kWaveNum)*delC1(kWaveNum)
                vecL2(DENS_VAR:PRES_VAR) = (1.0d0/4.0d0)*(1.0d0 + 2.0d0*lambdaDtDx + (4.0d0/3.0d0)*(lambdaDtDx)**2.0d0)*reig(DENS_VAR:PRES_VAR, kWaveNum) * delC2(kWaveNum)
                sigL(DENS_VAR:PRES_VAR) = sigL(DENS_VAR:PRES_VAR) + vecL(DENS_VAR:PRES_VAR) + vecL2(DENS_VAR:PRES_VAR)
            end if
        end do

        ! Let's make sure we copy all the cell-centered values to left and right states
        ! this will be just FOG
        gr_vL(DENS_VAR:NUMB_VAR, i) = gr_V(DENS_VAR:NUMB_VAR, i)
        gr_vR(DENS_VAR:NUMB_VAR, i) = gr_V(DENS_VAR:NUMB_VAR, i)

        !Now PLM reconstruction for dens, velx, and pres
        gr_vL(DENS_VAR:PRES_VAR, i) = C_MAT(DENS_VAR:PRES_VAR, 0) + sigL(DENS_VAR:PRES_VAR)
        gr_vR(DENS_VAR:PRES_VAR, i) = C_MAT(DENS_VAR:PRES_VAR, 0) + sigR(DENS_VAR:PRES_VAR)

    end do

    return
end subroutine soln_PPM
