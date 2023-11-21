subroutine sim_init()

#include "definition.h"

  use sim_data
  use read_initFile

  
  implicit none


  call read_initFileInt ('slug.init','sim_order',  sim_order)
  call read_initFileInt ('slug.init','sim_nStep',  sim_nStep)
  call read_initFileReal('slug.init','sim_cfl',    sim_cfl)
  call read_initFileReal('slug.init','sim_tmax',   sim_tmax)
  call read_initFileReal('slug.init','sim_outputIntervalTime',sim_outputIntervalTime)
  call read_initFileChar('slug.init','sim_riemann',sim_riemann)
  call read_initFileChar('slug.init','sim_limiter',sim_limiter)
  call read_initFileChar('slug.init','sim_name',sim_name)
  call read_initFileBool('slug.init','sim_charLimiting',sim_charLimiting)

  call read_initFileReal('slug.init','sim_densL',    sim_densL)
  call read_initFileReal('slug.init','sim_velxL',    sim_velxL)
  call read_initFileReal('slug.init','sim_presL',    sim_presL)
  call read_initFileReal('slug.init','sim_densR',    sim_densR)
  call read_initFileReal('slug.init','sim_velxR',    sim_velxR)
  call read_initFileReal('slug.init','sim_presR',    sim_presR)
  call read_initFileReal('slug.init','sim_gamma',    sim_gamma)
  call read_initFileReal('slug.init','sim_shockLoc1', sim_shockLoc)
  call read_initFileReal('slug.init','sim_smallPres', sim_smallPres)

  call read_initFileChar('slug.init','sim_bcType',sim_bcType)

  call read_initFileReal('slug.init','sim_ioTfreq',  sim_ioTfreq)
  call read_initFileInt ('slug.init','sim_ioNfreq',  sim_ioNfreq)

  ! For Blast 2 IC
  call read_initFileChar("slug.init", "IC_type", IC_type)
  
  if (IC_type == "sim_blast2") then
    print *, "fdlkajfdkla"
    call read_initFileReal('slug.init','sim_densM',    sim_densM)
    call read_initFileReal('slug.init','sim_velxM',    sim_velxM)
    call read_initFileReal('slug.init','sim_presM',    sim_presM)
    call read_initFileReal('slug.init','sim_shockLoc2', sim_shockLoc2)
  end if

  call sim_initBlock()



end subroutine sim_init
