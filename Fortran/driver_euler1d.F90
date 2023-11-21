program driver_euler1d

#include "definition.h"
  
  use sim_data
  use grid_data
  use io, only : io_writeOutput
  use bc
  use eos, only : eos_all

  implicit none

  real :: t,dt
  integer :: nStep,ioCounter,ioTimeFreqCounter
  real :: ioCheckTime
  
  t = 0.
  nStep = 0
  ioCounter = 0
  ioTimeFreqCounter = 0

  
  ! grid_init should be called first before sim_init
  call grid_init()
  call sim_init()

  write(*,*)''
  write(*,*)'================================================='
  write(*,*)'     AM 260 - CFD, 1D Euler FVM Code             '
  write(*,*)'      Written by Prof. Dongwook Lee              '
  write(*,*)'================================================='
  write(*,*)''

  
  ! write the initial condition
  write(*,*)''
  write(*,*)'       Initial condition was written!            '
  write(*,*)'================================================='
  write(*,*)'   Steps      Time              dt               '
  write(*,*)'================================================='
  write(*,*)''
  call io_writeOutput(t, nStep,ioCounter)

  do while ( (t < sim_tmax) .and. (nStep .le. sim_nStep))
     call cfl(dt)

     ! reduce dt if the next time step is larger than tmax
     if (t+dt >= sim_tmax) then
        dt = sim_tmax - t
     endif
     
     call soln_ReconEvolveAvg(dt)
     call soln_update(dt)


     ! call BC on primitive vars
     call bc_apply()

     ! call eos to make sure all through GC regions
     ! call eos_all !! Let's not call this here as it messes up bc
     
     ! write outputs every ioNfreq cycle or ioTfreq cycle
     
     ioCheckTime = sim_ioTfreq*real(ioTimeFreqCounter+1)
     if (t-dt < ioCheckTime .and. t>ioCheckTime) then
        write(*,*)''
        write(*,*)' Output no.',ioCounter+1, 'has been written      '
        write(*,*)'================================================='
        write(*,*)'   Steps      Time              dt               '
        write(*,*)'================================================='
        write(*,*)''
        ioCounter = ioCounter + 1
        ioTimeFreqCounter = ioTimeFreqCounter + 1
        call io_writeOutput(t, nStep,ioCounter)
     endif

     if (sim_ioNfreq > 0) then
     if (mod(nStep, sim_ioNfreq) == 0) then
        write(*,*)''
        write(*,*)' Output no.',ioCounter+1, 'has been written      '
        write(*,*)'================================================='
        write(*,*)'   Steps      Time              dt               '
        write(*,*)'================================================='
        write(*,*)''
        ioCounter = ioCounter + 1
        call io_writeOutput(t, nStep,ioCounter)
     endif
     endif
     
     ! update your time and step count
     t = t + dt
     nStep = nStep + 1

     ! Stop broken cases
     if (dt <= 0.0d0) then
      stop
     end if

     write(*,900)nstep,t,dt
  enddo


  !! Let's write the final result before exiting
  write(*,*)''
  write(*,*)' Final output no.',ioCounter+1, 'has been written'
  write(*,*)'================================================='
  write(*,*)'        The final tmax has reached, bye!         '
  write(*,*)'================================================='
  write(*,*)''
  call io_writeOutput(t, nStep,ioCounter+1)

  !! finalize and deallocate memories
  call grid_finalize()

900 format(1x,i5,f16.8,1x,f16.8)
  
  return
end program driver_euler1d
