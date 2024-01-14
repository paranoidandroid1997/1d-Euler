subroutine soln_NN(dt)

#include "definition.h"  
    
    use grid_data
    use sim_data
    use slopeLimiter
    use eigensystem
    use NN
    use ftorch

    implicit none
    real, intent(IN) :: dt
    integer :: i

    real, dimension(NUMB_WAVE) :: lambda
    real, dimension(NSYS_VAR,NUMB_WAVE) :: reig, leig
    logical :: conservative
    real, dimension(NSYS_VAR) :: vecL,vecR,sigL,sigR
    integer :: kWaveNum
    real :: lambdaDtDx
    real, dimension(NUMB_VAR)  :: delV,delL,delR
    real, dimension(NUMB_WAVE) :: delW
    integer :: nVar

    integer :: classification

    !real(kind = 4), dimension(20, 1) :: input
    real(kind = 8), dimension(gr_imax) :: scaledPres
    real(kind = 8), dimension(gr_imax) :: div
    real(kind = 8) :: minPres, maxPres


    type(torch_module) :: model
    integer, parameter :: n_inputs = 1
    type(torch_tensor), dimension(n_inputs) :: model_input_arr
    type(torch_tensor) :: model_output
    real(kind = 4), dimension(20,1), target  :: input
    real(kind = 4), dimension(1, 2), target   :: output

    ! Set up number of dimensions of input tensor and axis order
    integer, parameter :: in_dims = 2
    integer :: in_layout(in_dims) = [1,2]
    integer, parameter :: out_dims = 2
    integer :: out_layout(out_dims) = [1, 2]

    model = torch_module_load("./models/nn-02/script-nn-02.pt")

    minPres = minval(gr_V(PRES_VAR, :))
    maxPres = maxval(gr_V(PRES_VAR, :))
    scaledPres = -1.0 + ( scaledPres -  minPres) * 2.0 / (maxPres - minPres)


    div = 0
    div(2:(gr_imax - 1)) = (gr_V(VELX_VAR, 3:gr_imax) - gr_V(VELX_VAR, 1:(gr_imax - 2)))/(2.0 * gr_dx)
    
    do i = gr_ibeg-1, gr_iend+1

        ! input(1:7,1) = scaledPres((i - 3):(i + 3 + 1))
        ! input(8, 1) = gr_dx
        ! call classify(input, classification)
        ! predictions(i) = classification

        input(1:5,1) = gr_V(DENS_VAR, (i - 2):(i + 2))
        input(6:10,1) = gr_V(VELX_VAR, (i - 2):(i + 2))
        input(11:15,1) = gr_V(PRES_VAR, (i - 2):(i + 2))
        input(16:20,1) = div((i - 2):(i + 2))
        ! call classify(input, classification)
        !call classify_ftorch(input, classification)

        model_input_arr(1) = torch_tensor_from_array(transpose(input), in_layout, torch_kCPU)
        model_output = torch_tensor_from_array(output, out_layout, torch_kCPU)

        call torch_module_forward(model, model_input_arr, n_inputs, model_output)

        if (output(1, 1) >= output(1, 2)) then
            classification = 0
        else
            classification = 1
        end if

        call torch_tensor_delete(model_input_arr(1))
        call torch_tensor_delete(model_output)

        predictions(i) = classification

        if (classification == 1 ) then
            !print *, "shock"
            gr_vL(DENS_VAR:GAME_VAR,i) = gr_V(DENS_VAR:GAME_VAR,i)
            gr_vR(DENS_VAR:GAME_VAR,i) = gr_V(DENS_VAR:GAME_VAR,i)

        else
            ! we need conservative eigenvectors
            conservative = .false.

            call eigenvalues(gr_V(DENS_VAR:GAME_VAR,i),lambda)
            call left_eigenvectors (gr_V(DENS_VAR:GAME_VAR,i),conservative,leig)
            call right_eigenvectors(gr_V(DENS_VAR:GAME_VAR,i),conservative,reig)

            ! primitive limiting
            if (.not. sim_charLimiting) then
            do kWaveNum = 1, NUMB_WAVE
                ! slope limiting
                ! deltas in primitive vars
                delL(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR,i  )-gr_V(DENS_VAR:PRES_VAR,i-1)
                delR(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR,i+1)-gr_V(DENS_VAR:PRES_VAR,i  )
                do nVar = DENS_VAR,PRES_VAR
                    if (sim_limiter == 'minmod') then
                        call minmod(delL(nVar),delR(nVar),delV(nVar))
                    elseif (sim_limiter == 'vanLeer') then
                        call vanLeer(delL(nVar),delR(nVar),delV(nVar))
                    elseif (sim_limiter == 'mc') then
                        call mc(delL(nVar),delR(nVar),delV(nVar))
                    endif
                enddo
                ! project primitive delta to characteristic vars
                delW(kWaveNum) = dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),delV(DENS_VAR:PRES_VAR))
            enddo
            elseif (sim_charLimiting) then
            !STUDENTS: PLEASE FINISH THIS CHARACTERISTIC LIMITING
            !(THE IMPLEMENTATION SHOULD NOT BE LONGER THAN THE PRIMITIVE LIMITING CASE)
            ! 9.30 in lecture notes
            do kWaveNum = 1, NUMB_WAVE
                delL(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR,i  )-gr_V(DENS_VAR:PRES_VAR,i-1)
                delR(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR,i+1)-gr_V(DENS_VAR:PRES_VAR,i) 
                if (sim_limiter == 'minmod') then
                call minmod(dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delL(DENS_VAR:PRES_VAR)), &
                &dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delR(DENS_VAR:PRES_VAR)), delW(kWaveNum))
                elseif (sim_limiter == 'vanLeer') then
                call vanLeer(dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delL(DENS_VAR:PRES_VAR)), &
                &dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delR(DENS_VAR:PRES_VAR)), delW(kWaveNum))
                elseif (sim_limiter == 'mc') then
                call mc(dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delL(DENS_VAR:PRES_VAR)), &
                &dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delR(DENS_VAR:PRES_VAR)), delW(kWaveNum))
                endif
            end do
            endif

            ! set the initial sum to be zero
            sigL(DENS_VAR:ENER_VAR) = 0.
            sigR(DENS_VAR:ENER_VAR) = 0.
            vecL(DENS_VAR:ENER_VAR) = 0.
            vecR(DENS_VAR:ENER_VAR) = 0.
            
            do kWaveNum = 1, NUMB_WAVE
            ! lambdaDtDx = lambda*dt/dx
            lambdaDtDx = lambda(kWaveNum)*dt/gr_dx

            
            if (sim_riemann == 'roe') then
                ! STUDENTS: PLEASE FINISH THIS ROE SOLVER CASE
                ! THIS SHOULDN'T BE LONGER THAN THE HLL CASE
                
                !! My Code
                if (lambda(kWaveNum) > 0.0d0) then
                vecR(DENS_VAR:PRES_VAR) = 0.5*(1.0 - lambdaDtDx)*reig(DENS_VAR:PRES_VAR,kWaveNum)*delW(kWaveNum)
                sigR(DENS_VAR:PRES_VAR) = sigR(DENS_VAR:PRES_VAR) + vecR(DENS_VAR:PRES_VAR) 
                else if (lambda(kWaveNum) < 0.0d0) then
                vecL(DENS_VAR:PRES_VAR) = 0.5*(-1.0 - lambdaDtDx)*reig(DENS_VAR:PRES_VAR,kWaveNum)*delW(kWaveNum)
                sigL(DENS_VAR:PRES_VAR) = sigL(DENS_VAR:PRES_VAR) + vecL(DENS_VAR:PRES_VAR) 
                end if
                
            elseif (sim_riemann == 'hll') then
                    vecR(DENS_VAR:PRES_VAR) = 0.5*(1.0 - lambdaDtDx)*reig(DENS_VAR:PRES_VAR,kWaveNum)*delW(kWaveNum)
                    sigR(DENS_VAR:PRES_VAR) = sigR(DENS_VAR:PRES_VAR) + vecR(DENS_VAR:PRES_VAR)

                    vecL(DENS_VAR:PRES_VAR) = 0.5*(-1.0 - lambdaDtDx)*reig(DENS_VAR:PRES_VAR,kWaveNum)*delW(kWaveNum)
                    sigL(DENS_VAR:PRES_VAR) = sigL(DENS_VAR:PRES_VAR) + vecL(DENS_VAR:PRES_VAR)
            endif

            ! Let's make sure we copy all the cell-centered values to left and right states
            ! this will be just FOG
            gr_vL(DENS_VAR:NUMB_VAR,i) = gr_V(DENS_VAR:NUMB_VAR,i)
            gr_vR(DENS_VAR:NUMB_VAR,i) = gr_V(DENS_VAR:NUMB_VAR,i)
            
            ! Now PLM reconstruction for dens, velx, and pres
            gr_vL(DENS_VAR:PRES_VAR,i) = gr_V(DENS_VAR:PRES_VAR,i) + sigL(DENS_VAR:PRES_VAR)
            gr_vR(DENS_VAR:PRES_VAR,i) = gr_V(DENS_VAR:PRES_VAR,i) + sigR(DENS_VAR:PRES_VAR)
            
            end do

        end if
    end do

    call torch_module_delete(model)
    
    
    return
end subroutine soln_NN
    