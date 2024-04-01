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
   integer, save :: count = 0

   real, dimension(NUMB_WAVE) :: lambda
   real, dimension(NSYS_VAR, NUMB_WAVE) :: reig, leig
   logical :: conservative
   real, dimension(NSYS_VAR) :: vecL, vecR, sigL, sigR
   integer :: kWaveNum
   real :: lambdaDtDx
   real, dimension(NUMB_VAR)  :: delV, delL, delR, delC
   real, dimension(NUMB_WAVE) :: delW
   integer :: nVar

   integer :: classification

   !real(kind = 4), dimension(20, 1) :: input
   real(kind=8), dimension(gr_imax) :: scaledPres
   real(kind=8), dimension(gr_imax) :: div
   real(kind=8) :: minPres, maxPres
   !real(kind=8) :: meanDens, meanVel, meanPres, meanDiv, stdDens, stdVel, stdPres, stdDiv
   real(kind=8), save :: meanDens, meanVel, meanPres, meanDiv, stdDens, stdVel, stdPres, stdDiv

   type(torch_module) :: model
   integer, parameter :: n_inputs = 1
   type(torch_tensor), dimension(n_inputs) :: model_input_arr
   type(torch_tensor) :: model_output
   !real(kind=4), dimension(20, 1), target  :: input
   real(kind=4), dimension(1, 2), target   :: output

   real(kind=8), dimension(20, 1) :: input

   ! Set up number of dimensions of input tensor and axis order
   integer, parameter :: in_dims = 2
   integer :: in_layout(in_dims) = [1, 2]
   integer, parameter :: out_dims = 2
   integer :: out_layout(out_dims) = [1, 2]

   integer :: l

   ! Line to load model if using Ftorch
   !model = torch_module_load("./models/nn-02/script-nn-02.pt")

   ! Global normalization vars
   div = 0
   div(2:(gr_imax - 1)) = (gr_V(VELX_VAR, 3:gr_imax) - gr_V(VELX_VAR, 1:(gr_imax - 2)))/(2.0*gr_dx)
   meanDiv = sum(div)/size(div)
   stdDiv = sqrt((1.0/size(div))*sum((div - meanDiv)**2))

   if (count == 0) then
      meanDens = sum(gr_V(DENS_VAR, :))/size(gr_V(DENS_VAR, :))
      meanPres = sum(gr_V(PRES_VAR, :))/size(gr_V(PRES_VAR, :))
      meanVel = sum(gr_V(VELX_VAR, :))/size(gr_V(VELX_VAR, :))
      stdDens = sqrt((1.0/size(gr_V(DENS_VAR, :)))*sum((gr_V(DENS_VAR, :) - meanDens)**2))
      stdPres = sqrt((1.0/size(gr_V(PRES_VAR, :)))*sum((gr_V(PRES_VAR, :) - meanPres)**2))
      stdVel = sqrt((1.0/size(gr_V(VELX_VAR, :)))*sum((gr_V(VELX_VAR, :) - meanVel)**2))
   end if
   
   do i = gr_ibeg - 1, gr_iend + 1
      !if (.False.) then
      if (div(i) >= 0) then
         classification = 2
      else if (.False.) then
         ! Local normalization
         input(1:5, 1) = gr_V(DENS_VAR, (i - 2):(i + 2))
         input(1:5, 1) = (input(1:5, 1) - meanDens)/stdDens
         input(6:10, 1) = gr_V(VELX_VAR, (i - 2):(i + 2))
         input(6:10, 1) = (input(6:10, 1) - meanVel)/stdVel
         input(11:15, 1) = gr_V(PRES_VAR, (i - 2):(i + 2))
         input(11:15, 1) = (input(11:15, 1) - meanPres)/stdPres
         input(16:20, 1) = div((i - 2):(i + 2))
         input(16:20, 1) = (input(16:20, 1) - meanDiv)/stdDiv

         model_input_arr(1) = torch_tensor_from_array(transpose(input), in_layout, torch_kCPU)
         model_output = torch_tensor_from_array(output, out_layout, torch_kCPU)

         call torch_module_forward(model, model_input_arr, n_inputs, model_output)

         !if (output(1, 1) >= output(1, 2)) then
         output(1, :) = output(1, :) - maxval(output(1, :))
         if (exp(output(1, 2))/(exp(output(1, 1)) + exp(output(1, 2))) >= 0.99) then
            classification = 1.0
         else
            !print *, exp(output(1,2))/(exp(output(1,1)) + exp(output(1,2)))
            !print *, output
            classification = 0.0
         end if

         call torch_tensor_delete(model_input_arr(1))
         call torch_tensor_delete(model_output)
      else
         ! Local normalization
         input(1:5, 1) = gr_V(DENS_VAR, (i - 2):(i + 2))
         input(1:5, 1) = (input(1:5, 1) - meanDens)/stdDens
         input(6:10, 1) = gr_V(VELX_VAR, (i - 2):(i + 2))
         input(6:10, 1) = (input(6:10, 1) - meanVel)/stdVel
         input(11:15, 1) = gr_V(PRES_VAR, (i - 2):(i + 2))
         input(11:15, 1) = (input(11:15, 1) - meanPres)/stdPres
         input(16:20, 1) = div((i - 2):(i + 2))
         input(16:20, 1) = (input(16:20, 1) - meanDiv)/stdDiv

         call classify(input, classification)
      end if
      predictions(i) = classification

      !end do

      if (classification == 1) then
         !print *, "shock"
         gr_vL(DENS_VAR:GAME_VAR, i) = gr_V(DENS_VAR:GAME_VAR, i)
         gr_vR(DENS_VAR:GAME_VAR, i) = gr_V(DENS_VAR:GAME_VAR, i)

      else
         ! we need conservative eigenvectors
         conservative = .false.

         call eigenvalues(gr_V(DENS_VAR:GAME_VAR, i), lambda)
         call left_eigenvectors(gr_V(DENS_VAR:GAME_VAR, i), conservative, leig)
         call right_eigenvectors(gr_V(DENS_VAR:GAME_VAR, i), conservative, reig)

         ! primitive limiting
         if (.not. sim_charLimiting) then
         do kWaveNum = 1, NUMB_WAVE
            ! slope limiting
            ! deltas in primitive vars
            delL(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR, i) - gr_V(DENS_VAR:PRES_VAR, i - 1)
            delR(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR, i + 1) - gr_V(DENS_VAR:PRES_VAR, i)
            do nVar = DENS_VAR, PRES_VAR
               if (sim_limiter == 'minmod') then
                  call minmod(delL(nVar), delR(nVar), delV(nVar))
               elseif (sim_limiter == 'vanLeer') then
                  call vanLeer(delL(nVar), delR(nVar), delV(nVar))
               elseif (sim_limiter == 'mc') then
                  call mc(delL(nVar), delR(nVar), delV(nVar))
               elseif (sim_limiter == "center") then
                  delV(nVar) = gr_V(nVar, i + 1) - gr_V(nVar, i - 1)
               end if
            end do
            ! project primitive delta to characteristic vars
            delW(kWaveNum) = dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delV(DENS_VAR:PRES_VAR))
         end do
         elseif (sim_charLimiting) then
         !STUDENTS: PLEASE FINISH THIS CHARACTERISTIC LIMITING
         !(THE IMPLEMENTATION SHOULD NOT BE LONGER THAN THE PRIMITIVE LIMITING CASE)
         ! 9.30 in lecture notes
         do kWaveNum = 1, NUMB_WAVE
            delL(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR, i) - gr_V(DENS_VAR:PRES_VAR, i - 1)
            delR(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR, i + 1) - gr_V(DENS_VAR:PRES_VAR, i)
            delC(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR, i + 1) - gr_V(DENS_VAR:PRES_VAR, i - 1)
            if (sim_limiter == 'minmod') then
               call minmod(dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delL(DENS_VAR:PRES_VAR)), &
               &dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delR(DENS_VAR:PRES_VAR)), delW(kWaveNum))
            elseif (sim_limiter == 'vanLeer') then
               call vanLeer(dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delL(DENS_VAR:PRES_VAR)), &
               &dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delR(DENS_VAR:PRES_VAR)), delW(kWaveNum))
            elseif (sim_limiter == 'mc') then
               call mc(dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delL(DENS_VAR:PRES_VAR)), &
               &dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delR(DENS_VAR:PRES_VAR)), delW(kWaveNum))
            elseif (sim_limiter == "center") then
               delW(kWaveNum) = 0.5*dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), delC(DENS_VAR:PRES_VAR))
            end if
         end do
         end if

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
                  vecR(DENS_VAR:PRES_VAR) = 0.5*(1.0 - lambdaDtDx)*reig(DENS_VAR:PRES_VAR, kWaveNum)*delW(kWaveNum)
                  sigR(DENS_VAR:PRES_VAR) = sigR(DENS_VAR:PRES_VAR) + vecR(DENS_VAR:PRES_VAR)
               else if (lambda(kWaveNum) < 0.0d0) then
                  vecL(DENS_VAR:PRES_VAR) = 0.5*(-1.0 - lambdaDtDx)*reig(DENS_VAR:PRES_VAR, kWaveNum)*delW(kWaveNum)
                  sigL(DENS_VAR:PRES_VAR) = sigL(DENS_VAR:PRES_VAR) + vecL(DENS_VAR:PRES_VAR)
               end if

            elseif (sim_riemann == 'hll') then
               vecR(DENS_VAR:PRES_VAR) = 0.5*(1.0 - lambdaDtDx)*reig(DENS_VAR:PRES_VAR, kWaveNum)*delW(kWaveNum)
               sigR(DENS_VAR:PRES_VAR) = sigR(DENS_VAR:PRES_VAR) + vecR(DENS_VAR:PRES_VAR)

               vecL(DENS_VAR:PRES_VAR) = 0.5*(-1.0 - lambdaDtDx)*reig(DENS_VAR:PRES_VAR, kWaveNum)*delW(kWaveNum)
               sigL(DENS_VAR:PRES_VAR) = sigL(DENS_VAR:PRES_VAR) + vecL(DENS_VAR:PRES_VAR)
            end if

            ! Let's make sure we copy all the cell-centered values to left and right states
            ! this will be just FOG
            gr_vL(DENS_VAR:NUMB_VAR, i) = gr_V(DENS_VAR:NUMB_VAR, i)
            gr_vR(DENS_VAR:NUMB_VAR, i) = gr_V(DENS_VAR:NUMB_VAR, i)

            ! Now PLM reconstruction for dens, velx, and pres
            gr_vL(DENS_VAR:PRES_VAR, i) = gr_V(DENS_VAR:PRES_VAR, i) + sigL(DENS_VAR:PRES_VAR)
            gr_vR(DENS_VAR:PRES_VAR, i) = gr_V(DENS_VAR:PRES_VAR, i) + sigR(DENS_VAR:PRES_VAR)

         end do

      end if
   end do

   ! Uncomment if using Ftorch
   !call torch_module_delete(model)

   return
end subroutine soln_NN
