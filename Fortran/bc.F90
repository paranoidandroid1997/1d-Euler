module bc

#include "definition.h"

   use grid_data
   use sim_data, only: sim_bcType
   implicit none

contains

   subroutine bc_apply()
      implicit none
      if (sim_bcType == 'outflow') then
         call bc_outflow(gr_V)
      elseif (sim_bcType == 'reflect') then
         call bc_reflect(gr_V)
      elseif (sim_bcType == 'periodic') then
      elseif (sim_bcType == 'user') then
         ! STUDENTS: PLEASE IMPLEMENT THIS FOR THE SHU-OSHER PROBLEM
         call bc_user(gr_V)
      end if

   end subroutine bc_apply

   subroutine bc_outflow(V)
      implicit none
      real, dimension(NUMB_VAR, gr_imax), intent(INOUT) :: V
      integer :: i

      do i = 1, gr_ngc
         ! on the left GC
         V(1:NUMB_VAR, i) = V(1:NUMB_VAR, gr_ibeg)

         ! on the right GC
         V(1:NUMB_VAR, gr_iend + i) = V(1:NUMB_VAR, gr_iend)
      end do

      return
   end subroutine bc_outflow

   subroutine bc_reflect(V)
      implicit none
      real, dimension(NUMB_VAR, gr_imax), intent(INOUT) :: V
      integer :: i, k0, k1, k2, thing

      do i = 1, gr_ngc
         k0 = 2*gr_ngc + 1
         k1 = gr_iend - gr_ngc

         ! on the left GC
         V(:, i) = V(:, k0 - i)
         V(VELX_VAR, i) = -V(VELX_VAR, k0 - i)

         ! on the right GC
         V(:, k1 + k0 - i) = V(:, k1 + i)
         V(VELX_VAR, k1 + k0 - i) = -V(VELX_VAR, k1 + i)

         !  k2 = gr_ngc + 2 * (k1) + 1
         !  thing = k1 + gr_ngc + i
         !  V(       :,k1 + gr_ngc + i) = V(         :,k2 - thing)
         !  V(VELX_VAR,k1 + gr_ngc + i) =-V(  VELX_VAR,k2 - thing)

      end do

      return
   end subroutine bc_reflect

   subroutine bc_user(V)
      implicit none
      real, dimension(NUMB_VAR, gr_imax), intent(INOUT) :: V
      ! STUDENTS: PLEASE IMPLEMENT THIS FOR THE SHU-OSHER PROBLEM

    !! My code

      integer :: i

      do i = gr_i0, gr_imax
         if (gr_xCoord(i) < -4.0d0) then
            V(DENS_VAR, i) = 3.857143d0
            V(VELX_VAR, i) = 2.629369
            V(PRES_VAR, i) = 10.33333
         elseif (gr_xCoord(i) > 4.0d0) then
            V(DENS_VAR, i) = 1.0d0 + 0.20d0*sin(5.0*gr_xCoord(i))
            V(VELX_VAR, i) = 0.0d0
            V(PRES_VAR, i) = 1.0d0
         end if
      end do

      return
   end subroutine bc_user

end module bc
